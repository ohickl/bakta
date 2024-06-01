import concurrent.futures
import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

from Bio import SeqIO

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


log = logging.getLogger('T_RNA')


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


def run_trnascan_on_chunk(chunk_path: Path, txt_output_path: Path, fasta_output_path: Path, env: dict):
    cmd = [
        'tRNAscan-SE',
        '-G',
        '--output', str(txt_output_path),
        '--fasta', str(fasta_output_path),
        '--thread', "1",
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
        log.warning('tRNAs failed! tRNAscan-SE-error-code=%d', proc.returncode)
        raise Exception(f'tRNAscan-SE error! error code: {proc.returncode}')


def predict_t_rnas(genome: dict, contigs_path: Path):
    """Search for tRNA sequences."""

    txt_output_path = cfg.tmp_path.joinpath('trna.tsv')
    fasta_output_path = cfg.tmp_path.joinpath('trna.fasta')
    chunk_dir = cfg.tmp_path.joinpath('chunks_trna')
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
    chunk_fasta_output_paths = [chunk_dir.joinpath(f'chunk_{i}.fasta') for i in range(len(chunk_paths))]

    # Submit tasks to the executor
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg.threads) as executor:
        futures = [
            executor.submit(run_trnascan_on_chunk, chunk_path, chunk_txt_output_path, chunk_fasta_output_path, cfg.env)
            for chunk_path, chunk_txt_output_path, chunk_fasta_output_path in zip(chunk_paths, chunk_txt_output_paths, chunk_fasta_output_paths)
        ]

        # Collect results in the order of submission
        for future in futures:
            try:
                future.result()
            except Exception as e:
                log.error('A tRNAscan-SE run failed: %s', e)
                raise


    # Concatenate results (since the first three lines are headers, we take the header from the first file and skip the rest)
    header_seen = False
    with txt_output_path.open('w') as outfile:
        for chunk_txt_output_path in chunk_txt_output_paths:
            with chunk_txt_output_path.open() as infile:
                lines = infile.readlines()
                if not header_seen:
                    # Check if the file is empty
                    if not lines:
                        continue
                    outfile.writelines(lines)  # Write the header from the first file
                    header_seen = True
                else:
                    outfile.writelines(lines[3:])  # Skip the first three lines (header) for subsequent files

    with fasta_output_path.open('w') as outfile:
        for chunk_fasta_output_path in chunk_fasta_output_paths:
            with chunk_fasta_output_path.open() as infile:
                outfile.writelines(infile.readlines())

    # Clean up chunks
    for chunk_path in chunk_paths:
        chunk_path.unlink()
    for chunk_txt_output_path in chunk_txt_output_paths:
        chunk_txt_output_path.unlink()
    for chunk_fasta_output_path in chunk_fasta_output_paths:
        chunk_fasta_output_path.unlink()

    log.info('tRNA prediction completed successfully.')

    trnas = {}
    contigs = {c['id']: c for c in genome['contigs']}
    with txt_output_path.open() as fh:
        for line in fh.readlines()[3:]:  # skip first 3 lines
            (contig_id, trna_id, start, stop, trna_type, anti_codon, intron_begin, bounds_end, score, note) = line.split('\t')

            start, stop, strand = int(start), int(stop), bc.STRAND_FORWARD
            if(start > stop):  # reverse
                start, stop = stop, start
                strand = bc.STRAND_REVERSE
            contig_id = contig_id.strip()  # bugfix for extra single whitespace in tRNAscan-SE output

            trna = OrderedDict()
            trna['type'] = bc.FEATURE_T_RNA
            trna['contig'] = contig_id
            # Fix negative start position
            if int(start) < 0:
                log.warning(f'contig={contig_id}, trna_id={trna_id}, negative start position={start}, setting to 0.')
                start = '0'
            trna['start'] = start
            # Fix negative stop position
            if int(stop) < 0:
                log.warning(f'contig={contig_id}, trna_id={trna_id}, negative stop position={stop}, setting to 0.')
                stop = '0'
            trna['strand'] = strand
            trna['gene'] = None
            trna['product'] = 'tRNA-Xxx'
            if(trna_type != 'Undet' and trna_type != 'Sup'):
                aa_code = AMINO_ACID_DICT.get(trna_type.lower(), ('', None))[0]
                trna['gene'] = f'trn{aa_code}'
                trna['product'] = f'tRNA-{trna_type}({anti_codon.lower()})'
                trna['amino_acid'] = trna_type
                trna['anti_codon'] = anti_codon.lower()

            if('pseudo' in note):
                trna['pseudo'] = True

            trna['score'] = float(score)

            nt = bu.extract_feature_sequence(trna, contigs[contig_id])  # extract nt sequences
            trna['nt'] = nt

            trna['db_xrefs'] = []
            so_term = AMINO_ACID_DICT.get(trna_type.lower(), ('', None))[1]
            if(so_term):
                trna['db_xrefs'].append(so_term.id)

            key = f'{contig_id}.trna{trna_id}'
            trnas[key] = trna
            log.info(
                'contig=%s, start=%i, stop=%i, strand=%s, gene=%s, product=%s, score=%1.1f, nt=[%s..%s]',
                trna['contig'], trna['start'], trna['stop'], trna['strand'], trna.get('gene', ''), trna['product'], trna['score'], nt[:10], nt[-10:]
            )

    with fasta_output_path.open() as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            trna = trnas[record.id]
            nt = str(record.seq).upper()
            if('anti_codon' in trna and trna['amino_acid'].lower() not in ['fmet', 'ile2', 'sec', 'sup']):  # exclude fMet, Ile2 and Sec (INSDC wrong anticodon issue)
                anticodon_pos = trna['nt'].lower().find(trna['anti_codon'])
                if(anticodon_pos > -1):
                    if(trna['strand'] == bc.STRAND_FORWARD):
                        start = trna['start'] + anticodon_pos
                        stop = start + 2
                    else:
                        stop = trna['stop'] - anticodon_pos
                        start = stop - 2
                    trna['anti_codon_pos'] = (start, stop)
    trnas = list(trnas.values())
    log.info('predicted=%i', len(trnas))
    return trnas
