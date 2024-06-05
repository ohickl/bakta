import concurrent.futures
import logging
import re
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


RE_CRISPR = re.compile(r'(\d{1,8})\s+(\d{2})\s+(\d{1,3}\.\d)\s+(?:(\d{2})\s+)?([ATGCN]+)?\s+([ATGCN\.-]+)\s*(?:([ATGCN]+))?')


log = logging.getLogger('CRISPR')


def run_pilercr_on_chunk(chunk_path: Path, output_path: Path, env: dict):
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
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        log.warning('CRISPRs failed! pilercr-error-code=%d', proc.returncode)
        raise Exception(f'PILER-CR error! error code: {proc.returncode}')


def predict_crispr(genome: dict, contigs_path: Path):
    """Predict CRISPR arrays with PILER-CR."""

    output_path = cfg.tmp_path.joinpath('crispr.txt')
    chunk_dir = cfg.tmp_path.joinpath('chunks_crispr')
    chunk_dir.mkdir(parents=True, exist_ok=True)

    # Split the fasta file
    split_cmd = [
        'seqkit', 'split',
        '--quiet',
        '-p', str(cfg.threads),
        '-O', str(chunk_dir),
        str(contigs_path)
    ]
    sp.run(split_cmd, check=True)

    # Determine the extension of the input fasta file
    contig_fasta_ext = contigs_path.suffix
    chunk_paths = sorted(chunk_dir.glob(f'*{contig_fasta_ext}'))  # Ensure the chunk paths are ordered
    chunk_output_paths = [chunk_dir.joinpath(f'chunk_{i}.txt') for i in range(len(chunk_paths))]

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_pilercr_on_chunk, chunk_path, chunk_output_path, cfg.env)
            for chunk_path, chunk_output_path in zip(chunk_paths, chunk_output_paths)
        ]

        # Collect results in the order of submission
        for future in futures:
            try:
                future.result()
            except Exception as e:
                log.error('A PILER-CR run failed: %s', e)
                raise

    # Concatenate results
    header = []
    detail_reports = []
    similarity_summaries = []
    position_summaries = []
    crispr_array_det_rep_counter = 0
    crispr_array_sum_by_sim_counter = 0
    crispr_array_sum_by_pos_counter = 0

    report_line_seen = False
    summary_by_sim_line_seen = False
    summary_by_pos_line_seen = False

    summary_by_sim_line_header = 'Array          Sequence    Position      Length  # Copies  Repeat  Spacer  +  Consensus\n' \
                                 '=====  ================  ==========  ==========  ========  ======  ======  =  =========\n'
    # summary_by_pos_line_header = 'Array          Sequence    Position      Length  # Copies  Repeat  Spacer    Distance  Consensus\n' \
    #                              '=====  ================  ==========  ==========  ========  ======  ======  ==========  =========\n'

    
    for chunk_output_path in chunk_output_paths:

        report_empty_skips = 0
        summary_by_pos_empty_skips = 0
        summary_by_sim_empty_skips = 0

        with chunk_output_path.open() as infile:
            lines = infile.readlines()

            if not header:
                header = lines[:7]
            
            detail_report_started = False
            similarity_summary_started = False
            position_summary_started = False

            for line in lines[6:]:
                if 'DETAIL REPORT' in line:
                    detail_report_started = True
                    similarity_summary_started = False
                    position_summary_started = False
                    if not report_line_seen:
                        detail_reports.append(line)
                        report_line_seen = True
                elif 'SUMMARY BY SIMILARITY' in line:
                    detail_report_started = False
                    similarity_summary_started = True
                    position_summary_started = False
                    if not summary_by_sim_line_seen:
                        similarity_summaries.append(line)
                        similarity_summaries.append(summary_by_sim_line_header)
                        summary_by_sim_line_seen = True
                elif 'SUMMARY BY POSITION' in line:
                    detail_report_started = False
                    similarity_summary_started = False
                    position_summary_started = True
                    if not summary_by_pos_line_seen:
                        position_summaries.append(line)
                        summary_by_pos_line_seen = True
                else:
                    if line == '\n':
                        continue
                    elif detail_report_started:
                        # Update counter of Array
                        if line.startswith('Array '):
                            crispr_array_det_rep_counter += 1
                            line = f'Array {crispr_array_det_rep_counter}\n'
                        detail_reports.append(line)
                    elif similarity_summary_started:
                        if line.startswith('Array ') or line.startswith('===='):
                            continue
                        # Update counter of Array
                        # We need to use regex, since the array number will be in a line that doesnt start with 'Array ` or `====` or is empty
                        # it will be a number potentially preceded by a spaces and followed by spaces
                        if re.match(r'^\s*\d+\s+', line):
                            crispr_array_sum_by_sim_counter += 1
                            # Replace chunk number with actual array number and format number to be preceded by spaces up to lenght 5
                            line = re.sub(r'^\s*\d+\s+', f'{crispr_array_sum_by_sim_counter: >5}  ', line)
                        similarity_summaries.append(line)
                    elif position_summary_started:
                        # Update counter of Array same as with similarity summary
                        if re.match(r'^\s*\d+\s+', line):
                            crispr_array_sum_by_pos_counter += 1
                            line = re.sub(r'^\s*\d+\s+', f'{crispr_array_sum_by_pos_counter: >5}  ', line)
                        position_summaries.append(line)

    # Assert all counters are equal
    assert crispr_array_det_rep_counter == crispr_array_sum_by_sim_counter == crispr_array_sum_by_pos_counter

    # Modify header to catch final count of arrays found
    header[3] = f'{contigs_path}: {crispr_array_det_rep_counter} putative CRISPR arrays found.\n'

    # Write the final output
    with output_path.open('w') as outfile:
        outfile.writelines(header)
        outfile.writelines(detail_reports)
        outfile.writelines(similarity_summaries)
        outfile.writelines(position_summaries)

    # Clean up chunks
    for chunk_path in chunk_paths:
        chunk_path.unlink()
    for chunk_output_path in chunk_output_paths:
        chunk_output_path.unlink()

    log.info('CRISPR prediction completed successfully.')

    # parse crispr arrays
    crispr_arrays = {}
    contigs = {c['id']: c for c in genome['contigs']}
    with output_path.open() as fh:
        output_section = None
        contig_id = None
        array_id = None
        skip_lines = True
        crispr_array = None
        gap_count = 0
        for line in fh:
            line = line.strip()
            if(line == ''):
                continue
            if(line == 'DETAIL REPORT'):
                output_section = 'DETAIL'
                skip_lines = False
            elif(line == 'SUMMARY BY POSITION'):
                output_section = 'POSITION'
                skip_lines = False
            elif(line == 'SUMMARY BY SIMILARITY'):
                output_section = 'SIMILARITY'
                skip_lines = False
            elif(skip_lines is False):
                if(output_section == 'DETAIL'):
                    if(line[0:5] == 'Array'):
                        gap_count = 0
                        array_id = line.split()[1]
                        crispr_array = OrderedDict()
                        crispr_array['type'] = bc.FEATURE_CRISPR
                        crispr_array['strand'] = bc.STRAND_UNKNOWN
                        crispr_array['repeats'] = []
                        crispr_array['spacers'] = []
                        crispr_arrays[array_id] = crispr_array
                    elif(line[0] == '>'):
                        contig_id = line[1:]
                        crispr_array['contig'] = contig_id
                    elif(line[0] != '='):
                        m = RE_CRISPR.fullmatch(line)
                        if(m is not None):
                            position = int(m.group(1))
                            repeat_length = int(m.group(2))
                            repeat_seq = m.group(6)
                            spacer_seq = m.group(7)
                            crispr_repeat = OrderedDict()
                            crispr_repeat['strand'] = bc.STRAND_UNKNOWN
                            crispr_repeat['start'] = position - gap_count
                            crispr_repeat['stop'] = position + repeat_length - 1 - gap_count
                            crispr_array['repeats'].append(crispr_repeat)
                            log.debug('repeat: array-id=%s, start=%i, stop=%i', array_id, crispr_repeat['start'], crispr_repeat['stop'])
                            gap_count += repeat_seq.count('-')  # correct wrong PILER-CR detail positions by gaps
                            if(spacer_seq is not None):
                                spacer_seq = spacer_seq.upper()
                                spacer_length = len(spacer_seq)
                                crispr_spacer = OrderedDict()
                                crispr_spacer['strand'] = bc.STRAND_UNKNOWN
                                crispr_spacer['start'] = position + repeat_length  - gap_count
                                crispr_spacer['stop'] = position + repeat_length + spacer_length - 1 - gap_count
                                crispr_spacer['sequence'] = spacer_seq
                                crispr_array['spacers'].append(crispr_spacer)
                                spacer_genome_seq = bu.extract_feature_sequence(crispr_spacer, contigs[contig_id])
                                log.debug('spacer: array-id=%s, start=%i, stop=%i, genome-seq=%s, spacer-seq=%s', array_id, crispr_spacer['start'], crispr_spacer['stop'], spacer_genome_seq, spacer_seq)
                                assert spacer_seq == spacer_genome_seq  # assure PILER-CR provided sequence equals sequence extracted from genome
                elif(output_section == 'POSITION'):
                    if(line[0] == '>'):
                        contig_id = line[1:]
                    elif(line[0] != 'A' and line[0] != '='):
                        cols = line.split()
                        if(len(cols) == 8):
                            (array_id, contig, position, length, copies, repeat_length, spacer_length, repeat_consensus) = cols
                        else:
                            (array_id, contig, position, length, copies, repeat_length, spacer_length, distance, repeat_consensus) = cols
                        crispr_array = crispr_arrays[array_id]
                        positions = [c['start'] for c in crispr_array['spacers']] + [c['stop'] for c in crispr_array['spacers']] + [c['start'] for c in crispr_array['repeats']] + [c['stop'] for c in crispr_array['repeats']]
                        crispr_array['start'] = min(positions)
                        crispr_array['stop'] = max(positions)
                        crispr_array['product'] = f'CRISPR array with {copies} repeats of length {repeat_length}, consensus sequence {repeat_consensus} and spacer length {spacer_length}'
                        crispr_array['spacer_length'] = int(spacer_length)
                        crispr_array['repeat_length'] = int(repeat_length)
                        assert (len(crispr_array['repeats']) - int(copies)) <= 1, print(f"len(reps)={len(crispr_array['repeats'])}, int(copies)={int(copies)}")
                        crispr_array['repeat_consensus'] = repeat_consensus
                        crispr_array['db_xrefs'] = [so.SO_CRISPR.id]

                        nt = bu.extract_feature_sequence(crispr_array, contigs[contig_id])  # extract nt sequences
                        crispr_array['nt'] = nt
                        log.info(
                            'contig=%s, start=%i, stop=%i, spacer-length=%i, repeat-length=%i, # repeats=%i, repeat-consensus=%s, nt=[%s..%s]',
                            crispr_array['contig'], crispr_array['start'], crispr_array['stop'], crispr_array['spacer_length'], crispr_array['repeat_length'], len(crispr_array['repeats']), crispr_array['repeat_consensus'], nt[:10], nt[-10:]
                        )
    crispr_arrays = crispr_arrays.values()                        
    log.info('predicted=%i', len(crispr_arrays))
    return crispr_arrays

# Test case with argparse
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Predict CRISPR arrays with PILER-CR.')
    parser.add_argument('contigs', type=Path, help='Path to the contigs file in FASTA format.')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use.')
    args = parser.parse_args()
    genome = {'contigs': [{'id': 'contig1'}, {'id': 'contig2'}]}

    # Fake config
    cfg.tmp_path = Path('.')
    cfg.threads = args.threads

    predict_crispr(genome, args.contigs)
