import concurrent.futures
import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


HIT_COVERAGE = 0.3
HIT_COVERAGE_TRUNCATED = 0.8


log = logging.getLogger('R_RNA')


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
    cmd.append(str(db_path.joinpath('rRNA')))
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
        log.warning('rRNAs failed! cmscan-error-code=%d', proc.returncode)
        raise Exception(f'cmscan error! error code: {proc.returncode}')


def predict_r_rnas(genome: dict, contigs_path: Path):
    """Search for ribosomal RNA sequences."""

    output_path = cfg.tmp_path.joinpath('rrna.tsv')
    chunk_dir = cfg.tmp_path.joinpath('chunks')
    chunk_dir.mkdir(parents=True, exist_ok=True)

    # Calculate the -Z parameter
    z_value = 2 * genome['size'] / 1000000

    # Split the fasta file
    split_cmd = [
        'seqkit', 'split2',
        '-p', str(cfg.threads),
        '-O', str(chunk_dir),
        str(contigs_path)
    ]
    sp.run(split_cmd, check=True)

    # Determine the extension of the input fasta file
    contig_fasta_ext = contigs_path.suffix
    chunk_paths = list(chunk_dir.glob(f'*{contig_fasta_ext}'))
    chunk_output_paths = [chunk_dir.joinpath(f'chunk_{i}.tblout') for i in range(len(chunk_paths))]

    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_cmscan_on_chunk, chunk_path, chunk_output_path, cfg.db_path, z_value, cfg.env)
            for chunk_path, chunk_output_path in zip(chunk_paths, chunk_output_paths)
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error('A cmscan run failed: %s', e)
                raise

    # Concatenate results
    with output_path.open('w') as outfile:
        for chunk_output_path in chunk_output_paths:
            with chunk_output_path.open() as infile:
                outfile.write(infile.read())

    # Clean up chunks
    for chunk_path in chunk_paths:
        chunk_path.unlink()
    for chunk_output_path in chunk_output_paths:
        chunk_output_path.unlink()

    log.info('rRNA prediction completed successfully.')

    rrnas = []
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

                db_xrefs = [f'{bc.DB_XREF_GO}:0005840', f'{bc.DB_XREF_GO}:0003735']
                if(accession == 'RF00001'):
                    rrna_tag = '5S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF00001', f'{bc.DB_XREF_KOFAM}:K01985', so.SO_RRNA_5S.id]
                    consensus_length = 119
                elif(accession == 'RF00177'):
                    rrna_tag = '16S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF00177', f'{bc.DB_XREF_KOFAM}:K01977', so.SO_RRNA_16S.id]
                    consensus_length = 1533
                elif(accession == 'RF02541'):
                    rrna_tag = '23S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF02541', f'{bc.DB_XREF_KOFAM}:K01980', so.SO_RRNA_23S.id]
                    consensus_length = 2925
                else:
                    log.warning(
                        'unknown rRNA detected! accession=%s, contig=%s, start=%i, stop=%i, strand=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        accession, contig_id, start, stop, strand, length, truncated, score, evalue
                    )
                    continue

                coverage = length / consensus_length
                if(coverage < HIT_COVERAGE_TRUNCATED):
                    truncated = bc.FEATURE_END_UNKNOWN

                if(coverage < HIT_COVERAGE):
                    log.debug(
                        'discard low coverage: contig=%s, rRNA=%s, start=%i, stop=%i, strand=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e',
                        contig_id, rrna_tag, start, stop, strand, length, coverage, truncated, score, evalue
                    )
                else:
                    rrna = OrderedDict()
                    rrna['type'] = bc.FEATURE_R_RNA
                    rrna['contig'] = contig_id
                    rrna['start'] = start
                    rrna['stop'] = stop
                    rrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    if(accession == 'RF00001'):
                        rrna['gene'] = 'rrf'
                    elif(accession == 'RF00177'):
                        rrna['gene'] = 'rrs'
                    elif(accession == 'RF02541'):
                        rrna['gene'] = 'rrl'

                    if(truncated is None):
                        rrna['product'] = f'{rrna_tag} ribosomal RNA'
                    elif(truncated == bc.FEATURE_END_UNKNOWN):
                        rrna['product'] = f'(partial) {rrna_tag} ribosomal RNA'
                    elif(truncated == bc.FEATURE_END_5_PRIME):
                        rrna['product'] = f"(5' truncated) {rrna_tag} ribosomal RNA"
                    elif(truncated == bc.FEATURE_END_3_PRIME):
                        rrna['product'] = f"(3' truncated) {rrna_tag} ribosomal RNA"

                    if(truncated):
                        rrna['truncated'] = truncated

                    rrna['coverage'] = coverage
                    rrna['score'] = score
                    rrna['evalue'] = evalue
                    rrna['db_xrefs'] = db_xrefs

                    nt = bu.extract_feature_sequence(rrna, contigs[contig_id])  # extract nt sequences
                    rrna['nt'] = nt

                    rrnas.append(rrna)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e, nt=[%s..%s]',
                        rrna['contig'], rrna['start'], rrna['stop'], rrna['strand'], rrna['gene'], rrna['product'], length, coverage, truncated, score, evalue, nt[:10], nt[-10:]
                    )

    log.info('predicted=%i', len(rrnas))
    return rrnas
