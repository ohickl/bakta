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


def predict_tm_rnas(genome: dict, contigs_path: Path, n_threads: int, cfg: object):
    """Search for tmRNA sequences."""

    txt_output_path = cfg.tmp_path.joinpath('tmrna.tsv')
    chunk_dir = cfg.tmp_path.joinpath('chunks')
    chunk_dir.mkdir(parents=True, exist_ok=True)

    # Split the fasta file
    split_cmd = [
        'seqkit', 'split2',
        '-p', str(n_threads),
        '-O', str(chunk_dir),
        str(contigs_path)
    ]
    sp.run(split_cmd, check=True)

    # Determine the extension of the input fasta file
    contig_fasta_ext = contigs_path.suffix
    chunk_paths = list(chunk_dir.glob(f'*{contig_fasta_ext}'))
    chunk_txt_output_paths = [chunk_dir.joinpath(f'chunk_{i}.tsv') for i in range(len(chunk_paths))]

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_threads) as executor:
        futures = [
            executor.submit(run_aragorn_on_chunk, chunk_path, chunk_txt_output_path, cfg.translation_table, cfg.complete, cfg.env)
            for chunk_path, chunk_txt_output_path in zip(chunk_paths, chunk_txt_output_paths)
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error('An aragorn run failed: %s', e)
                raise

    # Concatenate results
    with txt_output_path.open('w') as outfile:
        for chunk_txt_output_path in chunk_txt_output_paths:
            with chunk_txt_output_path.open() as infile:
                outfile.write(infile.read())

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
