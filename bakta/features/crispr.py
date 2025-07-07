import concurrent.futures
import logging
import re
import subprocess as sp
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Any

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu

RE_CRISPR = re.compile(r'(\d{1,8})\s+(\d{2})\s+(\d{1,3}\.\d)\s+(?:(\d{1,2})\s+)?([ATGCN]+)?\s+([ATGCN\.-]+)\s*(?:([ATGCN]+))?')

log = logging.getLogger(__name__)


def run_pilercr_on_chunk(chunk_path: Path, output_path: Path, env: dict, is_retry=False):
    """
    Runs PILER-CR on a single chunk of sequences. Includes a retry mechanism
    that splits the chunk into smaller pieces if the initial run fails.
    """
    try:
        cmd = [
            'pilercr',
            '-in', str(chunk_path),
            '-out', str(output_path),
            '-noinfo',
            '-quiet'
        ]
        log.debug('cmd=%s', cmd)
        proc = sp.run(
            cmd,
            env=env,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True
        )
        if proc.returncode != 0:
            log.debug(f'stdout=\'{proc.stdout}\', stderr=\'{proc.stderr}\'')
            log.warning(f'CRISPRs failed for chunk {chunk_path}! pilercr-error-code={proc.returncode}')
            raise Exception(f'PILER-CR error for chunk {chunk_path}! error code: {proc.returncode}, stderr: {proc.stderr}, command: {cmd}')
    except Exception as e:
        if not is_retry:
            log.info(f'Attempting to split chunk {chunk_path} into smaller pieces due to failure.')
            # Create a temporary directory for the smaller chunks
            tmp_dir = chunk_path.parent.joinpath(f'{chunk_path.stem}_split')
            tmp_dir.mkdir(parents=True, exist_ok=True)
            
            # Split the chunk into smaller pieces (e.g., 10 sequences per piece)
            seqs_per_chunk = 10
            split_cmd = [
                'seqkit', 'split',
                '--quiet',
                '-s', str(seqs_per_chunk),
                '-O', str(tmp_dir),
                str(chunk_path)
            ]
            sp.run(split_cmd, check=True)
            
            # Process each smaller chunk recursively
            small_chunk_paths = list(tmp_dir.glob(f'{chunk_path.name}.*'))
            small_chunk_output_paths = [p.with_suffix('.txt') for p in small_chunk_paths]
            
            for small_chunk, small_output in zip(small_chunk_paths, small_chunk_output_paths):
                try:
                    # Call recursively with is_retry=True to prevent infinite loops
                    run_pilercr_on_chunk(small_chunk, small_output, env, is_retry=True)
                except Exception as inner_e:
                    log.error(f'Failed to process smaller chunk {small_chunk}: {inner_e}')
            
            # Merge the outputs of the successfully processed smaller chunks
            concatenate_pilercr_outputs(small_chunk_output_paths, output_path)
            
            # Clean up temporary split files and directory
            for p in small_chunk_paths: p.unlink()
            for p in small_chunk_output_paths: p.unlink()
            tmp_dir.rmdir()
        else:
            # If it's already a retry, re-raise the exception to avoid infinite recursion
            raise


def concatenate_pilercr_outputs(chunk_output_paths: list, final_output_path: Path):
    """
    Concatenates PILER-CR outputs from multiple chunks into a single, correctly
    formatted output file, re-numbering array IDs to be unique.
    """
    header = []
    detail_reports, similarity_summaries, position_summaries = [], [], []
    crispr_array_det_rep_counter = 0
    crispr_array_sum_by_sim_counter = 0
    crispr_array_sum_by_pos_counter = 0

    report_line_seen = False
    summary_by_sim_line_seen = False
    summary_by_pos_line_seen = False

    summary_by_sim_line_header = 'Array          Sequence    Position      Length  # Copies  Repeat  Spacer  +  Consensus\n' \
                                 '=====  ================  ==========  ==========  ========  ======  ======  =  =========\n'

    for chunk_path in chunk_output_paths:
        if not chunk_path.exists() or chunk_path.stat().st_size == 0:
            log.warning(f'Chunk output file {chunk_path} is missing or empty, skipping.')
            continue
        
        with chunk_path.open() as infile:
            lines = infile.readlines()

            if not header and len(lines) >= 7:
                header = lines[:7]
            
            detail_report_started = False
            similarity_summary_started = False
            position_summary_started = False

            for line in lines[6:]:
                if 'DETAIL REPORT' in line:
                    detail_report_started, similarity_summary_started, position_summary_started = True, False, False
                    if not report_line_seen:
                        detail_reports.append(line)
                        report_line_seen = True
                elif 'SUMMARY BY SIMILARITY' in line:
                    detail_report_started, similarity_summary_started, position_summary_started = False, True, False
                    if not summary_by_sim_line_seen:
                        similarity_summaries.append(line)
                        similarity_summaries.append(summary_by_sim_line_header)
                        summary_by_sim_line_seen = True
                elif 'SUMMARY BY POSITION' in line:
                    detail_report_started, similarity_summary_started, position_summary_started = False, False, True
                    if not summary_by_pos_line_seen:
                        position_summaries.append(line)
                        summary_by_pos_line_seen = True
                else:
                    if line == '\n':
                        continue
                    elif detail_report_started:
                        if line.startswith('Array '):
                            crispr_array_det_rep_counter += 1
                            line = f'Array {crispr_array_det_rep_counter}\n'
                        detail_reports.append(line)
                    elif similarity_summary_started:
                        if line.startswith('Array ') or line.startswith('===='):
                            continue
                        if re.match(r'^\s*\d+\s+', line):
                            crispr_array_sum_by_sim_counter += 1
                            line = re.sub(r'^\s*\d+\s+', f'{crispr_array_sum_by_sim_counter: >5}  ', line)
                        similarity_summaries.append(line)
                    elif position_summary_started:
                        if re.match(r'^\s*\d+\s+', line):
                            crispr_array_sum_by_pos_counter += 1
                            line = re.sub(r'^\s*\d+\s+', f'{crispr_array_sum_by_pos_counter: >5}  ', line)
                        position_summaries.append(line)

    if crispr_array_det_rep_counter > 0:
        assert crispr_array_det_rep_counter == crispr_array_sum_by_sim_counter == crispr_array_sum_by_pos_counter

    if header:
        header[3] = f'In total {crispr_array_det_rep_counter} putative CRISPR arrays found.\n'

    with final_output_path.open('w') as outfile:
        if header:
            outfile.writelines(header)
        outfile.writelines(detail_reports)
        outfile.writelines(similarity_summaries)
        outfile.writelines(position_summaries)


def predict_crispr(data: Dict[str, Any], fasta_chunk_paths: List[Path]) -> List[Dict]:
    """
    Predict CRISPR arrays with PILER-CR using a chunking approach.
    """
    final_output_path = cfg.tmp_path.joinpath('crispr.txt')
    # Create a dedicated directory for pilercr outputs to avoid name clashes
    chunk_output_dir = cfg.tmp_path.joinpath('pilercr_chunks_out')
    chunk_output_dir.mkdir(parents=True, exist_ok=True)

    # Generate output paths based on the input chunk paths
    chunk_output_paths = [chunk_output_dir.joinpath(f'{chunk_path.name}.txt') for chunk_path in fasta_chunk_paths]

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_pilercr_on_chunk, chunk_path, chunk_output_path, cfg.env)
            for chunk_path, chunk_output_path in zip(fasta_chunk_paths, chunk_output_paths)
        ]

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error('A PILER-CR chunk run failed: %s', e)
                executor.shutdown(wait=False, cancel_futures=True)
                raise

    # Concatenate results
    concatenate_pilercr_outputs(chunk_output_paths, final_output_path)

    # Clean up output chunks and directory
    for path in chunk_output_paths:
        if path.exists():
            path.unlink()
    chunk_output_dir.rmdir()

    log.info('CRISPR prediction completed successfully.')

    # Parse crispr arrays
    crispr_arrays = {}
    
    sequences = {s['id']: s for s in data['sequences']}
    
    with final_output_path.open() as fh:
        output_section = None
        sequence_id = None
        array_id = None
        skip_lines = True
        crispr_array = None
        gap_count = 0
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line == 'DETAIL REPORT':
                output_section, skip_lines = 'DETAIL', False
            elif line == 'SUMMARY BY POSITION':
                output_section, skip_lines = 'POSITION', False
            elif line == 'SUMMARY BY SIMILARITY':
                output_section, skip_lines = 'SIMILARITY', False
            elif not skip_lines:
                if output_section == 'DETAIL':
                    if line.startswith('Array'):
                        gap_count = 0
                        array_id = line.split()[1]
                        crispr_array = OrderedDict([
                            ('type', bc.FEATURE_CRISPR),
                            ('strand', bc.STRAND_UNKNOWN),
                            ('repeats', []),
                            ('spacers', [])
                        ])
                        crispr_arrays[array_id] = crispr_array
                    elif line.startswith('>'):
                        sequence_id = line[1:]
                        crispr_array['sequence'] = sequence_id
                    elif not line.startswith('='):
                        m = RE_CRISPR.fullmatch(line)
                        if m:
                            position, repeat_length = int(m.group(1)), int(m.group(2))
                            repeat_seq, spacer_seq = m.group(6), m.group(7)
                            
                            crispr_repeat = OrderedDict()
                            crispr_repeat['strand'] = bc.STRAND_UNKNOWN
                            crispr_repeat['start'] = position - gap_count
                            crispr_repeat['stop'] = position + repeat_length - 1 - gap_count
                            crispr_array['repeats'].append(crispr_repeat)
                            gap_count += repeat_seq.count('-')

                            if spacer_seq:
                                spacer_seq = spacer_seq.upper()
                                crispr_spacer = OrderedDict()
                                crispr_spacer['strand'] = bc.STRAND_UNKNOWN
                                crispr_spacer['start'] = position + repeat_length - gap_count
                                crispr_spacer['stop'] = position + repeat_length + spacer_length - 1 - gap_count
                                crispr_spacer['sequence'] = spacer_seq
                                crispr_array['spacers'].append(crispr_spacer)
                                spacer_genome_seq = bu.extract_feature_sequence(crispr_spacer, sequences[sequence_id])
                                assert spacer_seq == spacer_genome_seq, f'spacer_seq: {spacer_seq}\nspacer_genome_seq: {spacer_genome_seq}'

                elif output_section == 'POSITION':
                    if line.startswith('>'):
                        sequence_id = line[1:]
                    elif not line.startswith(('A', '=')):
                        cols = line.split()
                        if len(cols) == 8:
                            (array_id, _, _, _, copies, repeat_length, spacer_length, repeat_consensus) = cols
                        else:
                            (array_id, _, _, _, copies, repeat_length, spacer_length, _, repeat_consensus) = cols
                        
                        crispr_array = crispr_arrays[array_id]
                        all_pos = [p[k] for p in crispr_array['repeats'] + crispr_array['spacers'] for k in ('start', 'stop')]
                        crispr_array['start'] = min(all_pos)
                        crispr_array['stop'] = max(all_pos)
                        crispr_array['product'] = f'CRISPR array with {copies} repeats of length {repeat_length}, consensus sequence {repeat_consensus} and spacer length {spacer_length}'
                        crispr_array['spacer_length'] = int(spacer_length)
                        crispr_array['repeat_length'] = int(repeat_length)
                        assert (len(crispr_array['repeats']) - int(copies)) <= 1, print(f"len(reps)={len(crispr_array['repeats'])}, int(copies)={int(copies)}")
                        crispr_array['repeat_consensus'] = repeat_consensus
                        crispr_array['db_xrefs'] = [so.SO_CRISPR.id]

                        nt = bu.extract_feature_sequence(crispr_array, sequences[sequence_id])
                        crispr_array['nt'] = nt
                        log.info(
                            'sequence=%s, start=%i, stop=%i, spacer-length=%i, repeat-length=%i, # repeats=%i, repeat-consensus=%s, nt=[%s..%s]',
                            crispr_array['sequence'], crispr_array['start'], crispr_array['stop'], crispr_array['spacer_length'], crispr_array['repeat_length'], len(crispr_array['repeats']), crispr_array['repeat_consensus'], nt[:10], nt[-10:]
                        )

    result = list(crispr_arrays.values())
    log.info('predicted=%i', len(result))
    return result
