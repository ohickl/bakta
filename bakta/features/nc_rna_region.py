import concurrent.futures
import logging
import subprocess as sp
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Any

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


HIT_EVALUE = 1E-4


log = logging.getLogger(__name__)


def run_cmscan_on_chunk(chunk_path: Path, output_path: Path, db_path: Path, z_value: float, env: dict):
    """
    Runs cmscan on a single chunk of sequences.
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
    # The -Z parameter is only added if it's greater than 0
    if z_value > 0:
        cmd.extend(['-Z', str(z_value)])

    cmd.append(str(db_path.joinpath('ncRNA-regions')))
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
        log.warning('ncRNA regions failed for chunk %s! cmscan-error-code=%d', chunk_path.name, proc.returncode)
        raise Exception(f'cmscan error for chunk {chunk_path.name}! error code: {proc.returncode}')


def determine_class(description: str) -> so.SO:
    """
    Determines the Sequence Ontology class based on the ncRNA description.
    """
    description = description.lower()
    if 'leader' in description:
        return so.SO_CIS_REG_ATTENUATOR
    elif 'ribosomal frameshifting' in description:
        return so.SO_CIS_REG_FRAMESHIFT
    elif 'insertion sequence' in description:
        return so.SO_CIS_REG_RECODING_STIMULATION_REGION
    elif 'riboswitch' in description or 'sensor' in description:
        return so.SO_CIS_REG_RIBOSWITCH
    elif 'thermoregulator' in description or 'thermometer' in description or 'rose' in description:
        return so.SO_CIS_REG_THERMOMETER
    elif 'ribosome binding site' in description:
        return so.SO_CIS_REG_RIBOSOME_BINDING_SITE
    else:
        return None


def predict_nc_rnas(data: Dict[str, Any], fasta_chunk_paths: List[Path]) -> List[Dict]:
    """
    Search for non-coding RNA regions using a parallelized, chunk-based approach.
    """
    final_output_path = cfg.tmp_path.joinpath('ncrna-regions.tsv')
    chunk_output_dir = cfg.tmp_path.joinpath('ncrna_chunks_out')
    chunk_output_dir.mkdir(parents=True, exist_ok=True)

    z_value = (2 * data['stats']['size'] // 1000000) if data['stats']['size'] >= 1000000 else 0

    chunk_output_paths = [chunk_output_dir.joinpath(f'{chunk_path.name}.tblout') for chunk_path in fasta_chunk_paths]

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_cmscan_on_chunk, chunk_path, chunk_output_path, cfg.db_path, z_value, cfg.env)
            for chunk_path, chunk_output_path in zip(fasta_chunk_paths, chunk_output_paths)
        ]

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error('A cmscan run failed: %s', e)
                executor.shutdown(wait=False, cancel_futures=True)
                raise

    # Concatenate results
    with final_output_path.open('w') as outfile:
        header_written = False
        for chunk_path in chunk_output_paths:
            if not chunk_path.exists():
                continue
            with chunk_path.open() as infile:
                is_first_line = True
                for line in infile:
                    # Write header from the first file that has one
                    if not header_written and line.startswith('#'):
                        outfile.write(line)
                    # Write content lines (non-header, non-footer)
                    elif not line.startswith('#'):
                        outfile.write(line)
                # After writing the content of the first file, mark header as written
                if not header_written:
                    header_written = True

    # Clean up output chunks and directory
    for path in chunk_output_paths:
        if path.exists():
            path.unlink()
    chunk_output_dir.rmdir()

    log.info('ncRNA regions prediction completed successfully.')

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
                (subject, accession, sequence_id, _, _, _, _,
                 start, stop, strand, trunc, _, _, _, score, evalue,
                 _, description) = bc.RE_MULTIWHITESPACE.split(line.strip(), maxsplit=17)

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

                ncrna_region = OrderedDict()
                ncrna_region['type'] = bc.FEATURE_NC_RNA_REGION
                ncrna_region['class'] = determine_class(description)
                ncrna_region['sequence'] = sequence_id
                ncrna_region['start'] = start
                ncrna_region['stop'] = stop
                ncrna_region['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                ncrna_region['label'] = subject
                ncrna_region['product'] = description

                if ncrna_region['class'] is not None:
                    db_xrefs.append(ncrna_region['class'].id)
                else:
                    db_xrefs.append(so.SO_REGULATORY_REGION.id)

                if truncated:
                    ncrna_region['truncated'] = truncated

                ncrna_region['score'] = score
                ncrna_region['evalue'] = evalue
                ncrna_region['db_xrefs'] = db_xrefs

                nt = bu.extract_feature_sequence(ncrna_region, sequences[sequence_id])
                ncrna_region['nt'] = nt

                ncrnas.append(ncrna_region)
                log.info(
                    'seq=%s, start=%i, stop=%i, strand=%s, label=%s, product=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                    ncrna_region['sequence'], ncrna_region['start'], ncrna_region['stop'], ncrna_region['strand'], ncrna_region['label'], ncrna_region['product'], length, truncated, ncrna_region['score'], ncrna_region['evalue']
                )

    log.info('predicted=%i', len(ncrnas))
    return ncrnas
