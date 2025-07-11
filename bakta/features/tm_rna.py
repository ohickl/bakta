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


log = logging.getLogger(__name__)


def run_aragorn_on_chunk(chunk_path: Path, txt_output_path: Path, translation_table: int, complete: bool, env: dict):
    """
    Runs aragorn on a single chunk of sequences to find tmRNAs.
    """
    cmd = [
        'aragorn',
        '-m',  # detect tmRNAs
        f'-gc{translation_table}',
        '-w',  # batch mode
        '-o', str(txt_output_path),
        str(chunk_path)
    ]
    if complete:
        cmd.append('-c')  # complete circular sequence(s)
    else:
        cmd.append('-l')  # linear sequence(s)

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
        log.warning('tmRNA prediction failed for chunk %s! aragorn-error-code=%d', chunk_path.name, proc.returncode)
        raise Exception(f'aragorn error for chunk {chunk_path.name}! error code: {proc.returncode}')


def predict_tm_rnas(data: Dict[str, Any], fasta_chunk_paths: List[Path]) -> List[Dict]:
    """
    Search for tmRNA genes using a parallelized, chunk-based approach.
    """
    final_txt_output_path = cfg.tmp_path.joinpath('tmrna.tsv')
    chunk_output_dir = cfg.tmp_path.joinpath('tmrna_chunks_out')
    chunk_output_dir.mkdir(parents=True, exist_ok=True)

    chunk_txt_paths = [chunk_output_dir.joinpath(f'{p.name}.tsv') for p in fasta_chunk_paths]

    translation_table = data['genome']['translation_table']
    is_complete = data['genome']['complete']

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_aragorn_on_chunk, chunk_path, txt_path, translation_table, is_complete, cfg.env)
            for chunk_path, txt_path in zip(fasta_chunk_paths, chunk_txt_paths)
        ]

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error('An aragorn run failed: %s', e)
                executor.shutdown(wait=False, cancel_futures=True)
                raise

    # Concatenate results and recalculate summary statistics
    total_sequences = 0
    total_tmrna_genes = 0
    total_no_hit_sequences = 0
    with final_txt_output_path.open('w') as outfile:
        for path in chunk_txt_paths:
            if not path.exists() or path.stat().st_size == 0:
                continue

            with path.open() as infile:
                lines = [line for line in infile if line.strip()]
                if not lines:
                    continue

                summary_line = None
                # Find and remove the summary line, which is usually near the end.
                for i in range(len(lines) - 1, -1, -1):
                    if "sequences" in lines[i] and "tmRNA" in lines[i]:
                        summary_line = lines.pop(i)
                        break
                
                # Write all other lines to the combined output, ignoring the end marker.
                for line in lines:
                    if not line.startswith('>end'):
                        outfile.write(line)

                # Safely parse the summary line if it was found
                if summary_line:
                    parts = summary_line.split()
                    try:
                        # Corresponds to: {N} sequences {N} tmRNA genes, nothing found in {N} sequences
                        if len(parts) >= 9:
                            total_sequences += int(parts[0])
                            total_tmrna_genes += int(parts[2])
                            total_no_hit_sequences += int(parts[8])
                    except (ValueError, IndexError) as e:
                        log.warning(f"Could not parse summary line from {path.name}: '{summary_line.strip()}'. Error: {e}")

    # Write the final recalculated summary line
    sensitivity = (total_tmrna_genes / total_sequences * 100) if total_sequences > 0 else 0
    with final_txt_output_path.open('a') as outfile:
         outfile.write(f"{total_sequences} sequences {total_tmrna_genes} tmRNA genes, nothing found in {total_no_hit_sequences} sequences, ({sensitivity:.2f}% sensitivity)\n")

    # Clean up output chunks and directory
    for path in chunk_txt_paths:
        if path.exists():
            path.unlink()
    chunk_output_dir.rmdir()
    
    log.info('tmRNA prediction completed successfully.')

    tmrnas = []
    
    sequences = {s['id']: s for s in data['sequences']}
    with final_txt_output_path.open() as fh:
        sequence_id = None
        for line in fh:
            line = line.strip()
            if not line:
                continue

            cols = line.split()
            if line.startswith('>'):
                # Skip summary lines which can start with a number, e.g. "5 sequences..."
                # A real sequence header from aragorn is just ">contig_name"
                if line.startswith('>end') or cols[0].isnumeric():
                    continue
                sequence_id = cols[0][1:]
            elif len(cols) == 5:
                (nr, tm_type, location, tag_location, tag_aa) = cols
                if tm_type != "tmRNA":
                    continue
                
                strand = bc.STRAND_FORWARD
                if location.startswith('c'):
                    strand = bc.STRAND_REVERSE
                    location = location[1:]
                
                (start_str, stop_str) = location[1:-1].split(',')
                start, stop = int(start_str), int(stop_str)
                
                tag_start_str, tag_stop_str = tag_location.split(',')
                tag_start, tag_stop = int(tag_start_str), int(tag_stop_str)

                if start > 0 and stop > 0:
                    tmrna = OrderedDict()
                    tmrna['type'] = bc.FEATURE_TM_RNA
                    tmrna['sequence'] = sequence_id
                    tmrna['start'] = start
                    tmrna['stop'] = stop
                    tmrna['strand'] = strand
                    tmrna['gene'] = 'ssrA'
                    tmrna['product'] = 'transfer-messenger RNA, SsrA'
                    tmrna['db_xrefs'] = [so.SO_TMRNA.id]
                    
                    if strand == bc.STRAND_FORWARD:
                         tmrna['tag'] = {
                            'start': start + tag_start - 1,
                            'stop': start + tag_stop - 1,
                            'aa': tag_aa.replace('*', '')
                        }
                    else: # reverse strand
                        tmrna['tag'] = {
                            'start': stop - tag_stop + 1,
                            'stop': stop - tag_start + 1,
                            'aa': tag_aa.replace('*', '')
                        }

                    nt = bu.extract_feature_sequence(tmrna, sequences[sequence_id])
                    tmrna['nt'] = nt

                    tag = tmrna['tag']
                    tag_nt = bu.extract_feature_sequence({'start': tag['start'], 'stop': tag['stop'], 'strand': strand}, sequences[sequence_id])
                    tag['nt'] = tag_nt

                    if start > stop:
                        tmrna['edge'] = True

                    tmrnas.append(tmrna)
                    log.info(
                        'seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s',
                        tmrna['sequence'], tmrna['start'], tmrna['stop'], tmrna['strand'], tmrna['gene'], tmrna['product']
                    )

    log.info('predicted=%i', len(tmrnas))
    return tmrnas
