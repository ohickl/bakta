import concurrent.futures
import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger('TM_RNA')


def run_aragorn_on_chunk(chunk_path: Path, txt_output_path: Path, translation_table: int, complete: bool, env: dict):
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
        log.warning('tmRNAs failed! aragorn-error-code=%d', proc.returncode)
        raise Exception(f'aragorn error! error code: {proc.returncode}')


def predict_tm_rnas(genome: dict, contigs_path: Path):
    """Search for tmRNA sequences."""

    txt_output_path = cfg.tmp_path.joinpath('tmrna.tsv')
    chunk_dir = cfg.tmp_path.joinpath('chunks_tmrna')
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
    chunk_txt_output_paths = [chunk_dir.joinpath(f'chunk_{i}.tsv') for i in range(len(chunk_paths))]

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_aragorn_on_chunk, chunk_path, chunk_txt_output_path, cfg.translation_table, cfg.complete, cfg.env)
            for chunk_path, chunk_txt_output_path in zip(chunk_paths, chunk_txt_output_paths)
        ]

        # Collect results in the order of submission
        for future in futures:
            try:
                future.result()
            except Exception as e:
                log.error('An aragorn run failed: %s', e)
                raise

    # Concatenate results
    final_line_list = []
    with txt_output_path.open('w') as outfile:
        for chunk_txt_output_path in chunk_txt_output_paths:
            with chunk_txt_output_path.open() as infile:
                lines = infile.readlines()
                final_line_list.append(lines[-1].split('\t')[-1])  # Append the final line
                outfile.writelines(lines[:-1])  # Write all but the final line

    # Recalculate final line stats with all parts
    # Formata is: `x sequences y tmRNA genes, nothing found in z sequences, (n.nn% sensitivity)`
    sequences = 0
    tmrna_genes = 0
    no_hit_sequences = 0

    for line in final_line_list:
        cols = line.split()
        sequences += int(cols[0])
        tmrna_genes += int(cols[2])
        no_hit_sequences += int(cols[8])

    sensitivity = tmrna_genes / sequences * 100 if sequences > 0 else 0

    # Write final line
    with txt_output_path.open('a') as outfile:
        outfile.write(f'>end\t{sequences} sequences {tmrna_genes} tmRNA genes, nothing found in {no_hit_sequences} sequences, ({sensitivity:.2f}% sensitivity)\n')

    # Clean up chunks
    for chunk_path in chunk_paths:
        chunk_path.unlink()
    for chunk_txt_output_path in chunk_txt_output_paths:
        chunk_txt_output_path.unlink()

    log.info('tmRNA prediction completed successfully.')

    tmrnas = []
    contigs = {c['id']: c for c in genome['contigs']}
    with txt_output_path.open() as fh:
        contig_id = None
        for line in fh:
            line = line.strip()
            cols = line.split()
            if(line[0] == '>'):
                contig_id = cols[0][1:]
            elif(len(cols) == 5):
                (nr, type, location, tag_location, tag_aa) = line.split()
                strand = bc.STRAND_FORWARD
                if(location[0] == 'c'):
                    strand = bc.STRAND_REVERSE
                    location = location[1:]
                (start, stop) = location[1:-1].split(',')
                start = int(start)
                stop = int(stop)

                if(start > 0 and stop > 0):  # prevent edge tmRNA on linear sequences
                    tmrna = OrderedDict()
                    tmrna['type'] = bc.FEATURE_TM_RNA
                    tmrna['contig'] = contig_id
                    tmrna['start'] = start
                    tmrna['stop'] = stop
                    tmrna['strand'] = strand
                    tmrna['gene'] = 'ssrA'
                    tmrna['product'] = 'transfer-messenger RNA, SsrA'
                    tmrna['db_xrefs'] = [so.SO_TMRNA.id]

                    nt = bu.extract_feature_sequence(tmrna, contigs[contig_id])  # extract nt sequences
                    tmrna['nt'] = nt

                    if(start > stop):
                        tmrna['edge'] = True  # mark tmRNA as edge feature

                    tmrnas.append(tmrna)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, nt=[%s..%s]',
                        tmrna['contig'], tmrna['start'], tmrna['stop'], tmrna['strand'], tmrna['gene'], tmrna['product'], nt[:10], nt[-10:]
                    )
    log.info('predicted=%i', len(tmrnas))
    return tmrnas
