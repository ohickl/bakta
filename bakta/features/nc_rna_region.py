import concurrent.futures
import logging
import subprocess as sp

from collections import OrderedDict
from pathlib import Path

import bakta.config as cfg
import bakta.constants as bc
import bakta.so as so
import bakta.utils as bu


HIT_EVALUE = 1E-4


log = logging.getLogger('NC_RNA_REGION')


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
    cmd.append(str(db_path.joinpath('ncRNA-regions')))
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
        log.warning('ncRNA regions failed! cmscan-error-code=%d', proc.returncode)
        raise Exception(f'cmscan error! error code: {proc.returncode}')


def predict_nc_rna_regions(genome: dict, contigs_path: Path):
    """Search for non-coding RNA regions."""

    output_path = cfg.tmp_path.joinpath('ncrna-regions.tsv')
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

    log.info('ncRNA regions prediction completed successfully.')

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
                (subject, accession, contig_id, contig_acc, mdl, mdl_from, mdl_to,
                    start, stop, strand, trunc, passed, gc, bias, score, evalue,
                    inc, description) = bc.RE_MULTIWHITESPACE.split(line.strip(), maxsplit=17)

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

                    ncrna_region = OrderedDict()
                    ncrna_region['type'] = bc.FEATURE_NC_RNA_REGION
                    ncrna_region['class'] = determine_class(description)
                    ncrna_region['contig'] = contig_id
                    ncrna_region['start'] = start
                    ncrna_region['stop'] = stop
                    ncrna_region['strand'] = bc.STRAND_FORWARD if strand == '+' else bc.STRAND_REVERSE
                    ncrna_region['label'] = subject

                    if(truncated is None):
                        ncrna_region['product'] = description
                    elif(truncated == bc.FEATURE_END_UNKNOWN):
                        ncrna_region['product'] = f'(partial) {description}'
                    elif(truncated == bc.FEATURE_END_5_PRIME):
                        ncrna_region['product'] = f"(5' truncated) {description}"
                    elif(truncated == bc.FEATURE_END_3_PRIME):
                        ncrna_region['product'] = f"(3' truncated) {description}"

                    if(ncrna_region['class'] is not None):
                        db_xrefs.append(ncrna_region['class'].id)
                    else:
                        db_xrefs.append(so.SO_REGULATORY_REGION.id)

                    if(truncated):
                        ncrna_region['truncated'] = truncated

                    ncrna_region['score'] = score
                    ncrna_region['evalue'] = evalue
                    ncrna_region['db_xrefs'] = db_xrefs

                    nt = bu.extract_feature_sequence(ncrna_region, contigs[contig_id])  # extract nt sequences
                    ncrna_region['nt'] = nt

                    ncrnas.append(ncrna_region)
                    log.info(
                        'contig=%s, start=%i, stop=%i, strand=%s, label=%s, product=%s, length=%i, truncated=%s, score=%1.1f, evalue=%1.1e',
                        ncrna_region['contig'], ncrna_region['start'], ncrna_region['stop'], ncrna_region['strand'], ncrna_region['label'], ncrna_region['product'], length, truncated, ncrna_region['score'], ncrna_region['evalue']
                    )
    log.info('predicted=%i', len(ncrnas))
    return ncrnas


def determine_class(description: str) -> str:
    description = description.lower()
    if('leader' in description):
        return so.SO_CIS_REG_ATTENUATOR
    elif('ribosomal frameshifting' in description):
        return so.SO_CIS_REG_FRAMESHIFT
    elif('insertion sequence' in description):
        return so.SO_CIS_REG_RECODING_STIMULATION_REGION
    elif('riboswitch' in description or 'sensor' in description):
        return so.SO_CIS_REG_RIBOSWITCH
    elif('thermoregulator' in description or 'thermometer' in description or 'rose' in description):
        return so.SO_CIS_REG_THERMOMETER
    elif('ribosome binding site' in description):
        return so.SO_CIS_REG_RIBOSOME_BINDING_SITE
    else:
        None
