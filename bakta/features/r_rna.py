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


HIT_COVERAGE = 0.3
HIT_COVERAGE_TRUNCATED = 0.8


log = logging.getLogger(__name__)


def run_cmscan_on_chunk(chunk_path: Path, output_path: Path, db_path: Path, z_value: float, env: dict):
    """
    Runs cmscan on a single chunk of sequences for rRNA prediction.
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

    # Point to the rRNA database
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
        log.warning('rRNA prediction failed for chunk %s! cmscan-error-code=%d', chunk_path.name, proc.returncode)
        raise Exception(f'cmscan error for chunk {chunk_path.name}! error code: {proc.returncode}')


def predict_r_rnas(data: Dict[str, Any], fasta_chunk_paths: List[Path]) -> List[Dict]:
    """
    Search for ribosomal RNA genes using a parallelized, chunk-based approach.
    """
    final_output_path = cfg.tmp_path.joinpath('rrna.tsv')
    chunk_output_dir = cfg.tmp_path.joinpath('rrna_chunks_out')
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
                log.error('A cmscan run for rRNAs failed: %s', e)
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

    log.info('rRNA prediction completed successfully.')

    rrnas = []
    
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

                db_xrefs = [f'{bc.DB_XREF_GO}:0005840', f'{bc.DB_XREF_GO}:0003735']
                accession = accession.split('.')[0]
                if accession == 'RF00001':
                    rrna_tag = '5S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF00001', f'{bc.DB_XREF_KOFAM}:K01985', so.SO_RRNA_5S.id]
                    consensus_length = 119
                elif accession == 'RF00177':
                    rrna_tag = '16S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF00177', f'{bc.DB_XREF_KOFAM}:K01977', so.SO_RRNA_16S.id]
                    consensus_length = 1533
                elif accession == 'RF02541':
                    rrna_tag = '23S'
                    db_xrefs += [f'{bc.DB_XREF_RFAM}:RF02541', f'{bc.DB_XREF_KOFAM}:K01980', so.SO_RRNA_23S.id]
                    consensus_length = 2925
                else:
                    log.warning(
                        'unknown rRNA detected! accession=%s, seq=%s, start=%i, stop=%i, strand=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        accession, sequence_id, start, stop, strand, length, truncated, score, evalue
                    )
                    continue

                coverage = length / consensus_length
                if coverage < HIT_COVERAGE_TRUNCATED and truncated is not None:
                    truncated = bc.FEATURE_END_UNKNOWN

                if coverage < HIT_COVERAGE:
                    log.debug(
                        'discard low coverage: seq=%s, rRNA=%s, start=%i, stop=%i, strand=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e',
                        sequence_id, rrna_tag, start, stop, strand, length, coverage, truncated, score, evalue
                    )
                    continue
                
                rrna = OrderedDict()
                rrna['type'] = bc.FEATURE_R_RNA
                rrna['sequence'] = sequence_id
                rrna['start'] = start
                rrna['stop'] = stop
                rrna['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                if accession == 'RF00001':
                    rrna['gene'] = 'rrf'
                elif accession == 'RF00177':
                    rrna['gene'] = 'rrs'
                elif accession == 'RF02541':
                    rrna['gene'] = 'rrl'
                
                rrna['product'] = f'{rrna_tag} ribosomal RNA'

                if truncated:
                    rrna['truncated'] = truncated
                    # Update product description for truncated features
                    if truncated == bc.FEATURE_END_UNKNOWN:
                        rrna['product'] = f'(partial) {rrna_tag} ribosomal RNA'
                    elif truncated == bc.FEATURE_END_5_PRIME:
                        rrna['product'] = f"(5' truncated) {rrna_tag} ribosomal RNA"
                    elif truncated == bc.FEATURE_END_3_PRIME:
                        rrna['product'] = f"(3' truncated) {rrna_tag} ribosomal RNA"

                rrna['coverage'] = coverage
                rrna['score'] = score
                rrna['evalue'] = evalue
                rrna['db_xrefs'] = db_xrefs

                nt = bu.extract_feature_sequence(rrna, sequences[sequence_id])
                rrna['nt'] = nt

                rrnas.append(rrna)
                log.info(
                    'seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, length=%i, coverage=%0.3f, truncated=%s, score=%1.1f, evalue=%1.1e',
                    rrna['sequence'], rrna['start'], rrna['stop'], rrna['strand'], rrna['gene'], rrna['product'], length, coverage, truncated, score, evalue
                )

    log.info('predicted=%i', len(rrnas))
    return rrnas
