import concurrent.futures
import logging
import subprocess as sp
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Any

# Assuming these modules are available in the project structure
import bakta.features.annotation as ba
import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


HIT_EVALUE = 1E-4


log = logging.getLogger(__name__)


def run_cmscan_on_chunk(chunk_path: Path, output_path: Path, db_path: Path, z_value: float, env: dict):
    """
    Runs cmscan on a single chunk of sequences for ncRNA gene prediction.
    """
    cmd = [
        'cmscan',
        '--noali',
        '--cut_tc',
        '-g',
        '--nohmmonly',
        '--rfam',
        '--cpu', "1",  # Each chunk runs on a single CPU core
        '--tblout', str(output_path)
    ]
    if z_value > 0:
        cmd.extend(['-Z', str(z_value)])

    # Point to the ncRNA-genes database
    cmd.append(str(db_path.joinpath('ncRNA-genes')))
    cmd.append(str(chunk_path))
    log.debug('cmd=%s', cmd)
    proc = sp.run(
        cmd,
        env=env,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True
    )
    if proc.returncode != 0:
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('ncRNA prediction failed for chunk %s! cmscan-error-code=%d', chunk_path.name, proc.returncode)
        raise Exception(f'cmscan error for chunk {chunk_path.name}! error code: {proc.returncode}')


def determine_class(description: str) -> so.SO:
    """
    Determines the Sequence Ontology class for an ncRNA gene based on its description.
    """
    description = description.lower()
    if 'ribozyme' in description:
        return so.SO_NCRNA_GENE_RIBOZYME
    elif 'rnase p' in description:
        return so.SO_NCRNA_GENE_RNASEP
    elif 'antisense' in description:
        return so.SO_NCRNA_GENE_ANTISENSE
    else:
        return None


def predict_nc_rnas(data: Dict[str, Any], fasta_chunk_paths: List[Path]) -> List[Dict]:
    """
    Search for non-coding RNA genes using a parallelized, chunk-based approach.
    """
    final_output_path = cfg.tmp_path.joinpath('ncrna-genes.tsv')
    chunk_output_dir = cfg.tmp_path.joinpath('ncrna_genes_chunks_out')
    chunk_output_dir.mkdir(parents=True, exist_ok=True)

    z_value = (2 * data['stats']['size'] // 1000000) if data['stats']['size'] >= 1000000 else 0

    # Use tmp db path, if supplied
    db_path = cfg.tmp_db_path if cfg.tmp_db_path else cfg.db_path

    chunk_output_paths = [chunk_output_dir.joinpath(f'{chunk_path.name}.tblout') for chunk_path in fasta_chunk_paths]

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_cmscan_on_chunk, chunk_path, chunk_output_path, db_path, z_value, cfg.env)
            for chunk_path, chunk_output_path in zip(fasta_chunk_paths, chunk_output_paths)
        ]

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error('A cmscan run for ncRNA genes failed: %s', e)
                executor.shutdown(wait=False, cancel_futures=True)
                raise

    # Concatenate results
    with final_output_path.open('w') as outfile:
        header_written = False
        for chunk_path in chunk_output_paths:
            if not chunk_path.exists():
                continue
            with chunk_path.open() as infile:
                for line in infile:
                    if not header_written and line.startswith('#'):
                        outfile.write(line)
                    elif not line.startswith('#'):
                        outfile.write(line)
                if not header_written and chunk_path.stat().st_size > 0:
                    header_written = True

    # Clean up output chunks and directory
    for path in chunk_output_paths:
        if path.exists():
            path.unlink()
    chunk_output_dir.rmdir()

    log.info('ncRNA gene prediction completed successfully.')

    # Load Rfam to GO mapping
    rfam2go = {}
    rfam2go_path = cfg.db_path.joinpath('rfam-go.tsv')
    with rfam2go_path.open() as fh:
        for line in fh:
            if '\t' in line:
                (rfam, go) = line.strip().split('\t')
                if rfam in rfam2go:
                    rfam2go[rfam].append(go)
                else:
                    rfam2go[rfam] = [go]

    ncrnas = []
    
    sequences = {s['id']: s for s in data['sequences']}
    with final_output_path.open() as fh:
        for line in fh:
            if not line.startswith('#'):
                (
                    subject, accession, sequence_id, _, _, _, _,
                    start, stop, strand, trunc, _, _, _, score, evalue,
                    _, description
                ) = bc.RE_MULTIWHITESPACE.split(line.strip(), maxsplit=17)

                if strand == '-':
                    (start, stop) = (stop, start)
                (start, stop) = (int(start), int(stop))
                evalue = float(evalue)
                score = float(score)
                length = abs(stop - start) + 1
                
                truncated = None
                if trunc == "5'":
                    truncated = bc.FEATURE_END_5_PRIME
                elif trunc == "3'":
                    truncated = bc.FEATURE_END_3_PRIME

                if evalue > HIT_EVALUE:
                    log.debug(
                        'discard low E-value: seq=%s, start=%i, stop=%i, strand=%s, gene=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        sequence_id, start, stop, strand, subject, length, truncated, score, evalue
                    )
                    continue

                rfam_id = f'{bc.DB_XREF_RFAM}:{accession.split(".")[0]}'
                db_xrefs = [rfam_id]
                if rfam_id in rfam2go:
                    db_xrefs.extend(rfam2go[rfam_id])

                ncrna = OrderedDict()
                ncrna['type'] = bc.FEATURE_NC_RNA
                ncrna['class'] = determine_class(description)
                ncrna['sequence'] = sequence_id
                ncrna['start'] = start
                ncrna['stop'] = stop
                ncrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                
                gene = subject
                if ba.RE_PROTEIN_SYMBOL.fullmatch(gene):
                    gene = gene[0].lower() + gene[1:]
                    log.debug('fix gene: lowercase first char. new=%s, old=%s', gene, subject)
                ncrna['gene'] = gene
                ncrna['product'] = description

                if ncrna['class'] is not None:
                    db_xrefs.append(ncrna['class'].id)
                else:
                    db_xrefs.append(so.SO_NCRNA_GENE.id)

                if truncated:
                    ncrna['truncated'] = truncated

                ncrna['score'] = score
                ncrna['evalue'] = evalue
                ncrna['db_xrefs'] = db_xrefs

                nt = bu.extract_feature_sequence(ncrna, sequences[sequence_id])
                ncrna['nt'] = nt

                ncrnas.append(ncrna)
                log.info(
                    'seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                    ncrna['sequence'], ncrna['start'], ncrna['stop'], ncrna['strand'], ncrna['gene'], ncrna['product'], length, truncated, ncrna['score'], ncrna['evalue']
                )

    log.info('predicted=%i', len(ncrnas))
    return ncrnas
