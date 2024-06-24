import concurrent.futures
import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.features.annotation as ba
import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


HIT_EVALUE = 1E-4


log = logging.getLogger('NC_RNA')


def run_cmscan_on_chunk(chunk_path: Path, output_path: Path, db_path: Path, z_value: float, env: dict):
    cmd = [
        'cmscan',
        '--noali',
        '--cut_tc',
        '-g',
        '--nohmmonly',
        '--rfam',
        '--cpu', "1",
        '--tblout', str(output_path),
        '-Z', str(z_value)
    ]
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
        log.warning('ncRNAs failed! cmscan-error-code=%d', proc.returncode)
        raise Exception(f'cmscan error! error code: {proc.returncode}')


def predict_nc_rnas(genome: dict, chunk_paths: Path):
    """Search for non-coding RNA genes."""
    output_path = cfg.tmp_path.joinpath('ncrna-genes.tsv')
    chunk_dir = cfg.tmp_path.joinpath('chunks_ncrna')
    chunk_dir.mkdir(parents=True, exist_ok=True)

    # Calculate the -Z parameter
    z_value = 2 * genome['size'] / 1000000

    chunk_output_paths = [chunk_dir.joinpath(f'chunk_{i}.tblout') for i in range(len(chunk_paths))]

    # Use tmp db path, if supplied
    if cfg.tmp_db_path:
        db_path = cfg.tmp_db_path
    else:
        db_path = cfg.db_path

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_cmscan_on_chunk, chunk_path, chunk_output_path, db_path, z_value, cfg.env)
            for chunk_path, chunk_output_path in zip(chunk_paths, chunk_output_paths)
        ]

        # Collect results in the order of submission
        for future in futures:
            try:
                future.result()
            except Exception as e:
                log.error('A cmscan run failed: %s', e)
                raise

    # Concatenate results (since the first three lines are headers, we take the header from the first file and skip the rest)
    final_lines = []
    header_seen = False
    with output_path.open('w') as outfile:
        for chunk_output_path in chunk_output_paths:
            with chunk_output_path.open() as infile:
                lines = infile.readlines()
                if not header_seen:
                    # Check if the file is empty
                    if not lines:
                        continue
                    outfile.writelines(lines[:-10])  # Write the header from the first file
                    final_lines.append(lines[-10:])  # Add footer to final lines
                    header_seen = True
                else:
                    outfile.writelines(lines[2:-10])  # Skip the two first lines (header) and last 10 lines (footer) for subsequent files
                    final_lines.append(lines[-10:])  # Add footer to final lines

    # Write the final lines to output
    final_lines_flat = [line for line_list in final_lines for line in line_list]
    with output_path.open('a') as outfile:
        for line in final_lines_flat:
            outfile.write(line)

    # Clean up output chunks
    for chunk_output_path in chunk_output_paths:
        chunk_output_path.unlink()
    # Remove chunk directory
    chunk_dir.rmdir()

    log.info('ncRNA prediction completed successfully.')

    rfam2go = {}
    rfam2go_path = cfg.db_path.joinpath('rfam-go.tsv')
    with rfam2go_path.open() as fh:
        for line in fh:
            (rfam, go) = line.split('\t')
            if(rfam in rfam2go):
                rfam2go[rfam].append(go)
            else:
                rfam2go[rfam] = [go]

    ncrnas = []
    contigs = {c['id']: c for c in genome['contigs']}
    with output_path.open() as fh:
        for line in fh:
            if(line[0] != '#'):
                (
                    subject, accession, contig_id, contig_acc, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description
                ) = bc.RE_MULTIWHITESPACE.split(line.strip(), maxsplit=17)

                if(strand == '-'):
                    (start, stop) = (stop, start)
                (start, stop) = (int(start), int(stop))
                evalue = float(evalue)
                score = float(score)
                length = stop - start + 1
                if(trunc == "5'"):
                    truncated = bc.FEATURE_END_5_PRIME
                elif(trunc == "3'"):
                    truncated = bc.FEATURE_END_3_PRIME
                else:
                    truncated = None

                if(evalue > HIT_EVALUE):
                    log.debug(
                        'discard low E value: contig=%s, start=%i, stop=%i, strand=%s, gene=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        contig_id, start, stop, strand, subject, length, truncated, score, evalue
                    )
                else:
                    rfam_id = f'{bc.DB_XREF_RFAM}:{accession}'
                    db_xrefs = [rfam_id]
                    if(rfam_id in rfam2go):
                        db_xrefs += rfam2go[rfam_id]

                    ncrna = OrderedDict()
                    ncrna['type'] = bc.FEATURE_NC_RNA
                    ncrna['class'] = determine_class(description)
                    ncrna['contig'] = contig_id
                    ncrna['start'] = start
                    ncrna['stop'] = stop
                    ncrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    
                    gene = subject
                    if(ba.RE_PROTEIN_SYMBOL.fullmatch(gene)):
                        gene = gene[0].lower() + gene[1:]
                        log.debug('fix gene: lowercase first char. new=%s, old=%s', gene, subject)
                    ncrna['gene'] = gene

                    if(truncated is None):
                        ncrna['product'] = description
                    elif(truncated == bc.FEATURE_END_UNKNOWN):
                        ncrna['product'] = f'(partial) {description}'
                    elif(truncated == bc.FEATURE_END_5_PRIME):
                        ncrna['product'] = f"(5' truncated) {description}"
                    elif(truncated == bc.FEATURE_END_3_PRIME):
                        ncrna['product'] = f"(3' truncated) {description}"

                    if(ncrna['class'] is not None):
                        db_xrefs.append(ncrna['class'].id)
                    else:
                        db_xrefs.append(so.SO_NCRNA_GENE.id)

                    if(truncated):
                        ncrna['truncated'] = truncated

                    ncrna['score'] = score
                    ncrna['evalue'] = evalue
                    ncrna['db_xrefs'] = db_xrefs

                    nt = bu.extract_feature_sequence(ncrna, contigs[contig_id])  # extract nt sequences
                    ncrna['nt'] = nt

                    ncrnas.append(ncrna)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e, nt=[%s..%s]',
                        ncrna['contig'], ncrna['start'], ncrna['stop'], ncrna['strand'], ncrna['gene'], ncrna['product'], length, truncated, ncrna['score'], ncrna['evalue'], nt[:10], nt[-10:]
                    )
    log.info('predicted=%i', len(ncrnas))
    return ncrnas


def determine_class(description: str) -> str:
    description = description.lower()
    if('ribozyme' in description):
        return so.SO_NCRNA_GENE_RIBOZYME
    elif('rnase p' in description):
        return so.SO_NCRNA_GENE_RNASEP
    elif('antisense' in description):
        return so.SO_NCRNA_GENE_ANTISENSE
    else:
        None
