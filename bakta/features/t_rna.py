import concurrent.futures
import logging
import subprocess as sp
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Any

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger(__name__)


AMINO_ACID_DICT = {
    'ala': ('A', so.SO_TRNA_ALA),
    'gln': ('Q', so.SO_TRNA_GLN),
    'glu': ('E', so.SO_TRNA_GLU),
    'gly': ('G', so.SO_TRNA_GLY),
    'pro': ('P', so.SO_TRNA_PRO),
    'met': ('M', so.SO_TRNA_MET),
    'fmet':('fM', so.SO_TRNA_MET),
    'asp': ('D', so.SO_TRNA_ASP),
    'thr': ('T', so.SO_TRNA_THR),
    'val': ('V', so.SO_TRNA_VAL),
    'tyr': ('Y', so.SO_TRNA_TYR),
    'cys': ('C', so.SO_TRNA_CYS),
    'ile': ('I', so.SO_TRNA_ILE),
    'ile2':('I', so.SO_TRNA_ILE),
    'ser': ('S', so.SO_TRNA_SER),
    'leu': ('L', so.SO_TRNA_LEU),
    'trp': ('W', so.SO_TRNA_TRP),
    'lys': ('K', so.SO_TRNA_LYS),
    'asn': ('N', so.SO_TRNA_ASN),
    'arg': ('R', so.SO_TRNA_ARG),
    'his': ('H', so.SO_TRNA_HIS),
    'phe': ('F', so.SO_TRNA_PHE),
    'sec': ('U', so.SO_TRNA_SELCYS)
}

def run_trnascan_on_chunk(chunk_path: Path, txt_output_path: Path, fasta_output_path: Path, env: dict, threads: int = 1):
    """
    Runs tRNAscan-SE on a single chunk of sequences.
    """
    cmd = [
        'tRNAscan-SE',
        '-G',
        '--output', str(txt_output_path),
        '--fasta', str(fasta_output_path),
        '--thread', str(threads),
        str(chunk_path)
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
        log.warning('tRNA prediction failed for chunk %s! tRNAscan-SE-error-code=%d', chunk_path.name, proc.returncode)
        raise Exception(f'tRNAscan-SE error for chunk {chunk_path.name}! error code: {proc.returncode}')


def predict_t_rnas(data: Dict[str, Any], fasta_chunk_paths: List[Path]) -> List[Dict]:
    """
    Search for tRNA genes using a parallelized, chunk-based approach.
    """
    final_txt_output_path = cfg.tmp_path.joinpath('trna.tsv')
    final_fasta_output_path = cfg.tmp_path.joinpath('trna.fasta')
    chunk_output_dir = cfg.tmp_path.joinpath('trna_chunks_out')
    chunk_output_dir.mkdir(parents=True, exist_ok=True)

    chunk_txt_paths = [chunk_output_dir.joinpath(f'{p.name}.tsv') for p in fasta_chunk_paths]
    chunk_fasta_paths = [chunk_output_dir.joinpath(f'{p.name}.fasta') for p in fasta_chunk_paths]

    # tRNAscan-SE does not benefit from more than a few threads per process.
    # We parallelize by running on multiple chunks simultaneously.
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_trnascan_on_chunk, chunk_path, txt_path, fasta_path, cfg.env, 1)
            for chunk_path, txt_path, fasta_path in zip(fasta_chunk_paths, chunk_txt_paths, chunk_fasta_paths)
        ]

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error('A tRNAscan-SE run failed: %s', e)
                executor.shutdown(wait=False, cancel_futures=True)
                raise

    # Concatenate results
    with final_txt_output_path.open('w') as outfile:
        header_written = False
        for path in chunk_txt_paths:
            if path.exists() and path.stat().st_size > 0:
                with path.open() as infile:
                    lines = infile.readlines()
                    if not header_written:
                        outfile.writelines(lines)
                        header_written = True
                    else:
                        outfile.writelines(lines[3:]) # Skip header lines

    with final_fasta_output_path.open('w') as outfile:
        for path in chunk_fasta_paths:
            if path.exists():
                with path.open() as infile:
                    outfile.write(infile.read())

    # Clean up output chunks and directory
    for path in chunk_txt_paths + chunk_fasta_paths:
        if path.exists():
            path.unlink()
    chunk_output_dir.rmdir()

    log.info('tRNA prediction completed successfully.')

    trnas = {}
    
    sequences = {s['id']: s for s in data['sequences']}
    with final_txt_output_path.open() as fh:
        for line in fh.readlines()[3:]:
            (sequence_id, trna_id, start, stop, trna_type, anti_codon, _, _, score, note) = line.split('\t')

            start, stop, strand = int(start), int(stop), bc.STRAND_FORWARD
            if start > stop:
                start, stop = stop, start
                strand = bc.STRAND_REVERSE
            sequence_id = sequence_id.strip()

            trna = OrderedDict()
            trna['type'] = bc.FEATURE_T_RNA
            trna['sequence'] = sequence_id
            trna['start'] = start
            trna['stop'] = stop
            trna['strand'] = strand
            trna['gene'] = None
            trna['product'] = 'tRNA-Xxx'
            
            if trna_type != 'Undet' and trna_type != 'Sup':
                aa_code, so_term = AMINO_ACID_DICT.get(trna_type.lower(), ('', None))
                trna['gene'] = f'trn{aa_code}'
                trna['product'] = f'tRNA-{trna_type}({anti_codon.lower()})'
                trna['amino_acid'] = trna_type
                trna['anti_codon'] = anti_codon.lower()
                trna['db_xrefs'] = [so_term.id] if so_term else []
            else:
                trna['db_xrefs'] = []

            if 'pseudo' in note:
                trna[bc.PSEUDOGENE] = True

            trna['score'] = float(score)

            nt = bu.extract_feature_sequence(trna, sequences[sequence_id])
            trna['nt'] = nt

            key = f'{sequence_id}.trna{trna_id}'
            trnas[key] = trna
            log.info(
                'seq=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, score=%1.1f',
                trna['sequence'], trna['start'], trna['stop'], trna['strand'], trna.get('gene', ''), trna['product'], trna['score']
            )

    with final_fasta_output_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            if record.id in trnas:
                trna = trnas[record.id]
                if 'anti_codon' in trna and trna['amino_acid'].lower() not in ['fmet', 'ile2', 'sec', 'sup']:
                    anticodon_pos = trna['nt'].lower().find(trna['anti_codon'])
                    if anticodon_pos > -1:
                        if trna['strand'] == bc.STRAND_FORWARD:
                            ac_start = trna['start'] + anticodon_pos
                            ac_stop = ac_start + 2
                        else:
                            ac_stop = trna['stop'] - anticodon_pos
                            ac_start = ac_stop - 2
                        trna['anti_codon_pos'] = (ac_start, ac_stop)

    trnas_list = list(trnas.values())
    log.info('predicted=%i', len(trnas_list))
    return trnas_list
