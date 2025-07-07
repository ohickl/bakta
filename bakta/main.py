import atexit
import logging
import sys
import time

from datetime import datetime
from pathlib import Path
import subprocess as sp

import bakta
import bakta.constants as bc
import bakta.config as cfg
import bakta.io.fasta as fasta
import bakta.io.json as json
import bakta.io.tsv as tsv
import bakta.io.gff as gff
import bakta.io.insdc as insdc
import bakta.expert.amrfinder as exp_amr
import bakta.expert.protein_sequences as exp_aa_seq
import bakta.expert.protein_hmms as exp_aa_hmms
import bakta.features.annotation as anno
import bakta.features.t_rna as t_rna
import bakta.features.tm_rna as tm_rna
import bakta.features.r_rna as r_rna
import bakta.features.nc_rna as nc_rna
import bakta.features.nc_rna_region as nc_rna_region
import bakta.features.crispr as crispr
import bakta.features.orf as orf
import bakta.features.cds as feat_cds
import bakta.features.s_orf as s_orf
import bakta.features.gaps as gaps
import bakta.features.ori as ori
import bakta.db as db
import bakta.utils as bu
import bakta.ups as ups
import bakta.ips as ips
import bakta.psc as psc
import bakta.pscc as pscc
import bakta.plot as plot


def main():
    args = bu.parse_arguments()  # parse arguments

    ############################################################################
    # Setup logging
    ############################################################################
    cfg.prefix = args.prefix if args.prefix else Path(args.genome).stem
    output_path = cfg.check_output_path(args.output, args.force)
    bu.setup_logger(output_path, cfg.prefix, args)
    log = logging.getLogger('MAIN')

    ############################################################################
    # Checks and configurations
    # - check parameters and setup global configuration
    # - test database
    # - test binary dependencies
    ############################################################################
    cfg.setup(args)  # check parameters and prepare global configuration
    cfg.db_info = db.check(cfg.db_path)
    bu.test_dependencies()
    if(cfg.verbose):
        print(f'Bakta v{cfg.version}')
        print('Options and arguments:')
        print(f'\tinput: {cfg.genome_path}')
        print(f"\tdb: {cfg.db_path}, version {cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}")
        if(cfg.replicons): print(f'\treplicon table: {cfg.replicons}')
        if(cfg.prodigal_tf): print(f'\tprodigal training file: {cfg.prodigal_tf}')
        if(cfg.regions): print(f'\tregion table: {cfg.regions}')
        if(cfg.user_proteins): print(f'\tuser proteins: {cfg.user_proteins}')
        if(cfg.user_hmms): print(f'\tuser hmms: {cfg.user_hmms}')
        print(f'\ttranslation table: {cfg.translation_table}')
        if(cfg.taxon): print(f'\ttaxon: {cfg.taxon}')
        if(cfg.plasmid): print(f'\tplasmid: {cfg.plasmid}')
        if(cfg.gram != '?'): print(f'\tgram: {cfg.gram}')
        if(cfg.locus): print(f'\tlocus prefix: {cfg.locus}')
        if(cfg.locus_tag): print(f'\tlocus tag prefix: {cfg.locus_tag}')
        if(cfg.meta): print(f'\tmeta mode: {cfg.meta}')
        if(cfg.complete): print(f'\tcomplete replicons: {cfg.complete}')
        print(f'\toutput: {cfg.output_path}')
        if(cfg.force): print(f'\tforce: {cfg.force}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        if(cfg.compliant): print(f'\tINSDC compliant: {cfg.compliant}')
        if(cfg.keep_sequence_headers): print(f'\tkeep/sequence headers: {cfg.keep_sequence_headers}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\tthreads: {cfg.threads}')
        if(cfg.debug): print(f'\tdebug: {cfg.debug}')
        if(cfg.skip_trna): print(f'\tskip tRNA: {cfg.skip_trna}')
        if(cfg.skip_tmrna): print(f'\tskip tmRNA: {cfg.skip_tmrna}')
        if(cfg.skip_rrna): print(f'\tskip rRNA: {cfg.skip_rrna}')
        if(cfg.skip_ncrna): print(f'\tskip ncRNA: {cfg.skip_ncrna}')
        if(cfg.skip_ncrna_region): print(f'\tskip ncRNA region: {cfg.skip_ncrna_region}')
        if(cfg.skip_crispr): print(f'\tskip CRISPR: {cfg.skip_crispr}')
        if(cfg.skip_cds): print(f'\tskip CDS: {cfg.skip_cds}')
        if(cfg.skip_sorf): print(f'\tskip sORF: {cfg.skip_sorf}')
        if(cfg.skip_gap): print(f'\tskip gap: {cfg.skip_gap}')
        if(cfg.skip_ori): print(f'\tskip oriC/V/T: {cfg.skip_ori}')
        if(cfg.skip_filter): print(f'\tskip feature overlap filters: {cfg.skip_filter}')
        if(cfg.skip_plot): print(f'\tskip plot: {cfg.skip_plot}')
        if(cfg.skip_write_genbank_embl): print(f'\tskip write GenBank/EMBL: {cfg.skip_write_genbank_embl}')
    
    if(cfg.debug):
        print(f"Bakta runs in DEBUG mode! Temporary data will not be destroyed at: {cfg.tmp_path}\n")
    else:
        atexit.register(bu.cleanup, log, cfg.tmp_path)  # register cleanup exit hook

    ############################################################################
    # Import genome
    # - parse sequences in Fasta file
    # - apply sequence length filter
    # - rename sequences
    ############################################################################
    print('Parse genome sequences...')
    try:
        sequences = fasta.import_sequences(cfg.genome_path)
        log.info('imported sequences=%i', len(sequences))
        print(f'\timported: {len(sequences)}')
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
    replicons = bu.parse_replicon_table(cfg.replicons) if cfg.replicons else None
    sequences, complete_genome = bu.qc_sequences(sequences, replicons)
    print(f'\tfiltered & revised: {len(sequences)}')
    no_chromosomes = len([seq for seq in sequences if seq['type'] == bc.REPLICON_CHROMOSOME])
    if(no_chromosomes > 0):
        print(f"\tchromosomes: {no_chromosomes}")
    no_plasmids = len([seq for seq in sequences if seq['type'] == bc.REPLICON_PLASMID])
    if(no_plasmids > 0):
        print(f"\tplasmids: {no_plasmids}")
    no_contigs = len([seq for seq in sequences if seq['type'] == bc.REPLICON_CONTIG])
    if(no_contigs > 0):
        print(f"\tcontigs: {no_contigs}")
    if(len(sequences) == 0):
        log.warning('no valid sequences!')
        sys.exit('Error: input file contains no valid sequences.')
    sequences_path = cfg.tmp_path.joinpath('sequences.fna')
    fasta.export_sequences(sequences, sequences_path)
    data = {
        'genome': {
            'genus': cfg.genus,
            'species': cfg.species,
            'strain': cfg.strain,
            'taxon': cfg.taxon,
            'complete': cfg.complete or complete_genome,
            'gram': cfg.gram,
            'translation_table': cfg.translation_table
        },
        'stats': {
            'size': sum([seq['length'] for seq in sequences])
        },
        'features': [],
        'sequences': sequences
    }
    if(cfg.plasmid):
        data['genome']['plasmid'] = cfg.plasmid
    print('\nStart annotation...')


    # Split assembly into chuncks for parallel processing
    print(f'split assembly into {cfg.threads} chunks...')

    # Change fasta chunks path for more speed tmp_db_path is set
    if cfg.tmp_db_path:
        assembly_chunks_parent = cfg.tmp_db_path
    else:
        assembly_chunks_parent = cfg.tmp_path

    fasta_chunk_dir = assembly_chunks_parent.joinpath('assembly_chunks')
    fasta_chunk_dir.mkdir(parents=True, exist_ok=True)

    # Split the fasta file
    print(f'\tsplitting assembly into {cfg.threads} chunks...')
    start_time = time.perf_counter()
    split_cmd = [
        'seqkit', 'split',
        '--quiet',
        '-p', str(cfg.threads),
        '-O', str(fasta_chunk_dir),
        str(sequences_path)
    ]
    sp.run(split_cmd, check=True)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\ttime: {time_min:.2f} min')

    # Determine the extension of the input fasta file and generate sorted chunk paths
    contig_fasta_ext = sequences_path.suffix
    fasta_chunk_paths = sorted(fasta_chunk_dir.glob(f'*{contig_fasta_ext}'))  # Ensure the chunk paths are ordered


    # Split assembly into chuncks for parallel processing
    print(f'split assembly into {cfg.threads} chunks...')

    # Change fasta chunks path for more speed tmp_db_path is set
    if cfg.tmp_db_path:
        assembly_chunks_parent = cfg.tmp_db_path
    else:
        assembly_chunks_parent = cfg.tmp_path

    fasta_chunk_dir = assembly_chunks_parent.joinpath('assembly_chunks')
    fasta_chunk_dir.mkdir(parents=True, exist_ok=True)

    # Split the fasta file
    print(f'\tsplitting assembly into {cfg.threads} chunks...')
    start_time = time.perf_counter()
    split_cmd = [
        'seqkit', 'split',
        '--quiet',
        '-p', str(cfg.threads),
        '-O', str(fasta_chunk_dir),
        str(sequences_path)
    ]
    sp.run(split_cmd, check=True)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\ttime: {time_min:.2f} min')

    # Determine the extension of the input fasta file and generate sorted chunk paths
    contig_fasta_ext = sequences_path.suffix
    fasta_chunk_paths = sorted(fasta_chunk_dir.glob(f'*{contig_fasta_ext}'))  # Ensure the chunk paths are ordered

    ############################################################################
    # tRNA prediction
    ############################################################################
    if(cfg.skip_trna):
        print('skip tRNA prediction...')
    else:
        print('predict tRNAs...')
        start_time = time.perf_counter()
        log.debug('start tRNA prediction')
        trnas = t_rna.predict_t_rnas(data, fasta_chunk_paths)
        data['features'].extend(trnas)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(trnas)} in {time_min:.2f} min")

    ############################################################################
    # tmRNA prediction
    ############################################################################
    if(cfg.skip_tmrna):
        print('skip tmRNA prediction...')
    else:
        print('predict tmRNAs...')
        start_time = time.perf_counter()
        log.debug('start tmRNA prediction')
        tmrnas = tm_rna.predict_tm_rnas(data, fasta_chunk_paths)
        data['features'].extend(tmrnas)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(tmrnas)} in {time_min:.2f} min")

    ############################################################################
    # rRNA prediction
    ############################################################################
    if(cfg.skip_rrna):
        print('skip rRNA prediction...')
    else:
        # Copy rRNA database to tmp db directory, if set
        if cfg.tmp_db_path:
            bu.rsync_copy(f'{cfg.db_path}/rRNA*', cfg.tmp_db_path)

        print('predict rRNAs...')
        start_time = time.perf_counter()
        log.debug('start rRNA prediction')
        rrnas = r_rna.predict_r_rnas(data, fasta_chunk_paths)
        data['features'].extend(rrnas)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(rrnas)} in {time_min:.2f} min")

        # Clean up tmp rRNA database, if set
        if cfg.tmp_db_path:
            for file in Path(cfg.tmp_db_path).iterdir():
                file.unlink()

    ############################################################################
    # ncRNA gene prediction
    ############################################################################
    if(cfg.skip_ncrna):
        print('skip ncRNA prediction...')
    else:
        # Copy ncRNA database to tmp db directory, if set
        if cfg.tmp_db_path:
            bu.rsync_copy(f'{cfg.db_path}/ncRNA-genes*', cfg.tmp_db_path)
            
        print('predict ncRNAs...')
        start_time = time.perf_counter()
        log.debug('start ncRNA prediction')
        ncrnas = nc_rna.predict_nc_rnas(data, fasta_chunk_paths)
        data['features'].extend(ncrnas)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(ncrnas)} in {time_min:.2f} min")

        # Clean up tmp ncRNA database, if set
        if cfg.tmp_db_path:
            for file in Path(cfg.tmp_db_path).iterdir():
                file.unlink()

    ############################################################################
    # ncRNA region prediction
    ############################################################################
    if(cfg.skip_ncrna_region):
        print('skip ncRNA region prediction...')
    else:
        # Copy ncRNA region database to tmp db directory, if set
        if cfg.tmp_db_path:
            bu.rsync_copy(f'{cfg.db_path}/ncRNA-regions*', cfg.tmp_db_path)
        
        print('predict ncRNA regions...')
        start_time = time.perf_counter()
        log.debug('start ncRNA region prediction')
        ncrnas = nc_rna.predict_nc_rnas(data, fasta_chunk_paths)
        data['features'].extend(ncrnas)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(ncrnas)} in {time_min:.2f} min")

        # Clean up tmp ncRNA database, if set
        if cfg.tmp_db_path:
            for file in Path(cfg.tmp_db_path).iterdir():
                file.unlink()

    ############################################################################
    # CRISPR prediction
    ############################################################################
    if(cfg.skip_crispr):
        print('skip CRISPR array prediction...')
    else:
        print('predict CRISPR arrays...')
        start_time = time.perf_counter()
        log.debug('start CRISPR prediction')
        crisprs = crispr.predict_crispr(data, fasta_chunk_paths)
        data['features'].extend(crisprs)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(crisprs)} in {time_min:.2f} min")

    # Clean up assembly chuncks and chunk dir
    print('remove assembly chunks and dirs...')
    start_time = time.perf_counter()
    for chunk_path in fasta_chunk_paths:
        chunk_path.unlink()
    fasta_chunk_dir.rmdir()
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\ttime: {time_min:.2f} min')

    ############################################################################
    # CDS prediction
    # - Prodigal prediction
    # - lookup UPS matches
    # - lookup IPS matches
    # - search PSC for unannotated CDSs
    # - conduct expert systems analysis
    # - lookup & combine annotations
    # - analyze hypotheticals
    ############################################################################
    if(cfg.skip_cds):
        print('skip CDS prediction...')
    else:
        print('predict & annotate CDSs...')
        start_time = time.perf_counter()
        log.debug('predict CDS')
        cdss = feat_cds.predict(data)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tpredicted: {len(cdss)} in {time_min:.2f} min")

        if(len(cdss) > 0):
            log.debug('detect spurious CDS')
            start_time = time.perf_counter()
            discarded_cdss = orf.detect_spurious(cdss)
            cdss = [cds for cds in cdss if 'discarded' not in cds]
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\tdiscarded spurious: {len(discarded_cdss)} in {time_min:.2f} min')
        
        if(len(cdss) > 0):
            log.debug('revise translational exceptions')
            start_time = time.perf_counter()
            no_revised = feat_cds.revise_translational_exceptions(data, cdss)
            cdss = [cds for cds in cdss if 'discarded' not in cds]
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\trevised translational exceptions: {no_revised} in {time_min:.2f} min')
        
        if(cfg.regions):
            log.debug('import user-provided CDS regions')
            start_time = time.perf_counter()
            imported_cdss = feat_cds.import_user_cdss(data, cfg.regions)
            cdss.extend(imported_cdss)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\timported CDS regions: {len(imported_cdss)} in {time_min:.2f} min')

        if(len(cdss) > 0):
            if(cfg.db_info['type'] == 'full'):
                log.debug('lookup CDS UPS/IPS')
            start_time = time.perf_counter()

            # Copy `bakta.db` here, if tmp_db_path is set
            if cfg.tmp_db_path:
                bu.rsync_copy(f'{cfg.db_path}/bakta.db', cfg.tmp_db_path)

                cdss_ups, cdss_not_found_ups = ups.lookup(cdss)
                cdss_ips, cdss_not_found_ips = ips.lookup(cdss_ups)

            # Immediately clean up the tmp db, if set to maximize space
            if cfg.tmp_db_path:
                Path(cfg.tmp_db_path.joinpath('bakta.db')).unlink()

                cdss_not_found = cdss_not_found_ups + cdss_not_found_ips
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\tdetected IPSs: {len(cdss_ips)} in {time_min:.2f} min')
            else:
                cdss_not_found = [*cdss]
                print(f'\tskip UPS/IPS detection with light db version')

            if(len(cdss_not_found) > 0):
                if(cfg.db_info['type'] == 'full'):
                    log.debug('search CDS PSC')
                    start_time = time.perf_counter()
                    cdss_psc, cdss_pscc, cdss_not_found = psc.search(cdss_not_found)
                    end_time = time.perf_counter()
                    time_min = (end_time - start_time) / 60
                    print(f'\tfound PSCs: {len(cdss_psc)} in {time_min:.2f} min')
                    print(f'\tfound PSCCs: {len(cdss_pscc)} in {time_min:.2f} min')
                else:
                    log.debug('search CDS PSCC')
                    start_time = time.perf_counter()
                    cdss_pscc, cdss_not_found = pscc.search(cdss_not_found)
                    end_time = time.perf_counter()
                    time_min = (end_time - start_time) / 60
                    print(f'\tfound PSCCs: {len(cdss_pscc)} in {time_min:.2f} min')
            
            print('\tlookup annotations...')
            start_time = time.perf_counter()
            log.debug('lookup CDS PSCs')

            # Copy `bakta.db` here, if tmp_db_path is set
            if cfg.tmp_db_path:
                bu.rsync_copy(f'{cfg.db_path}/bakta.db', cfg.tmp_db_path)

            psc.lookup(cdss)  # lookup PSC info
            pscc.lookup(cdss)  # lookup PSCC info

            # Immediately clean up the tmp db, if set to maximize space
            if cfg.tmp_db_path:
                Path(cfg.tmp_db_path.joinpath('bakta.db')).unlink()

            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

            print('\tconduct expert systems...')  # conduct expert systems annotation
            start_time = time.perf_counter()
            cds_aa_path = cfg.tmp_path.joinpath('cds.expert.faa')
            orf.write_internal_faa(cdss, cds_aa_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

            log.debug('conduct expert system: amrfinder')
            start_time = time.perf_counter()
            expert_amr_found = exp_amr.search(cdss, cds_aa_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\t\tamrfinder: {len(expert_amr_found)} in {time_min:.2f} min')
            log.debug('conduct expert system: aa seqs')

            # If tmp db dir is set, copy dimanod db to tmp db dir
            if cfg.tmp_db_path:
                bu.rsync_copy(f'{cfg.db_path}/expert-protein-sequences.dmnd', cfg.tmp_db_path)
                diamond_db_path = cfg.tmp_db_path.joinpath('expert-protein-sequences.dmnd')
            else:
                diamond_db_path = cfg.db_path.joinpath('expert-protein-sequences.dmnd')

            start_time = time.perf_counter()
            expert_aa_found = exp_aa_seq.search(cdss, cds_aa_path, 'expert_proteins', diamond_db_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\t\tprotein sequences: {len(expert_aa_found)} in {time_min:.2f} min')

            # Cleanup tmp diamond db
            if cfg.tmp_db_path:
                diamond_db_path.unlink()

            if(cfg.user_proteins):
                log.debug('conduct expert system: user aa seqs')
                start_time = time.perf_counter()
                user_aa_path = cfg.tmp_path.joinpath('user-proteins.faa')
                exp_aa_seq.write_user_protein_sequences(user_aa_path)
                user_aa_found = exp_aa_seq.search(cdss, cds_aa_path, 'user_proteins', user_aa_path)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\t\tuser protein sequences: {len(user_aa_found)} in {time_min:.2f} min')

            if(cfg.user_hmms):
                log.debug('conduct expert system: user HMM')
                user_hmm_found = exp_aa_hmms.search(cdss, cfg.user_hmms)
                print(f'\t\tuser HMM sequences: {len(user_hmm_found)}')

            print('\tcombine annotations and mark hypotheticals...')
            log.debug('combine CDS annotations')
            start_time = time.perf_counter()
            for cds in cdss:
                anno.combine_annotation(cds)  # combine IPS & PSC annotations and mark hypotheticals
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

            hypotheticals = [cds for cds in cdss if 'hypothetical' in cds and 'edge' not in cds and cds.get('start_type', 'Edge') != 'Edge']
            if(len(hypotheticals) > 0  and  not cfg.skip_pseudo):
                if(cfg.db_info['type'] == 'full'):
                    print('\tdetect pseudogenes...')
                    start_time = time.perf_counter()
                    log.debug('search pseudogene candidates')
                    pseudo_candidates = feat_cds.predict_pseudo_candidates(hypotheticals)
                    end_time = time.perf_counter()
                    time_min = (end_time - start_time) / 60
                    print(f'\t\tcandidates: {len(pseudo_candidates)} in {time_min:.2f} min')

                    start_time = time.perf_counter()
                    pseudogenes = feat_cds.detect_pseudogenes(pseudo_candidates, cdss, data) if len(pseudo_candidates) > 0 else []

                # Copy `bakta.db` here, if tmp_db_path is set
                if cfg.tmp_db_path:
                    bu.rsync_copy(f'{cfg.db_path}/bakta.db', cfg.tmp_db_path)

                    psc.lookup(pseudogenes, pseudo=True)
                    pscc.lookup(pseudogenes, pseudo=True)

                # Immediately clean up the tmp db, if set to maximize space
                if cfg.tmp_db_path:
                    Path(cfg.tmp_db_path.joinpath('bakta.db')).unlink()

                    for pseudogene in pseudogenes:
                        anno.combine_annotation(pseudogene)
                    end_time = time.perf_counter()
                    time_min = (end_time - start_time) / 60
                    print(f'\t\tverified: {len(pseudogenes)} in {time_min:.2f} min')
            
                else:
                    print(f'\tskip pseudogene detection with light db version')
            hypotheticals = [cds for cds in cdss if 'hypothetical' in cds]
            if(len(hypotheticals) > 0):
                log.debug('analyze hypotheticals')
                start_time = time.perf_counter()
                print(f'\tanalyze hypothetical proteins: {len(hypotheticals)}')
                pfam_hits = feat_cds.predict_pfam(hypotheticals)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f"\t\tdetected Pfam hits: {len(pfam_hits)} in {time_min:.2f} min")
                
                start_time = time.perf_counter()
                feat_cds.analyze_proteins(hypotheticals)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\t\tcalculated proteins statistics in {time_min:.2f} min')
            
            print('\trevise special cases...')
            start_time = time.perf_counter()
            feat_cds.revise_special_cases_annotated(data, cdss)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

        data['features'].extend(cdss)

    ############################################################################
    # sORF prediction
    # - in-mem sORF extraction
    # - overlap filtering (tRNA, tmRNA, rRNA, CDS)
    # - lookup UPS matches
    # - lookup IPS matches
    # - filter sORFs w/o IPS match
    ############################################################################
    if(cfg.skip_sorf):
        print('skip sORF prediction...')
    else:
        print('detect & annotate sORF...')
        log.debug('extract sORF')
        start_time = time.perf_counter()
        sorfs = s_orf.extract(data)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tdetected: {len(sorfs)} in {time_min:.2f} min')

        log.debug('apply sORF overlap filter')
        start_time = time.perf_counter()
        sorfs, discarded_sorfs = s_orf.overlap_filter(data, sorfs)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tdiscarded due to overlaps: {len(discarded_sorfs)} in {time_min:.2f} min')

        if(len(sorfs) > 0):
            log.debug('detect spurious sORF')
            start_time = time.perf_counter()
            discarded_sorfs = orf.detect_spurious(sorfs)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\tdiscarded spurious: {len(discarded_sorfs)} in {time_min:.2f} min')
            sorfs = [sorf for sorf in sorfs if 'discarded' not in sorf]

        log.debug('lookup sORF UPS/IPS')
        start_time = time.perf_counter()

        # Copy `bakta.db` here, if tmp_db_path is set
        if cfg.tmp_db_path:
            bu.rsync_copy(f'{cfg.db_path}/bakta.db', cfg.tmp_db_path)

        sorf_upss, sorfs_not_found = ups.lookup(sorfs)
        sorf_ipss, tmp = ips.lookup(sorf_upss)

        # Immediately clean up the tmp db, if set to maximize space
        if cfg.tmp_db_path:
            Path(cfg.tmp_db_path.joinpath('bakta.db')).unlink()

        sorfs_not_found.extend(tmp)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tdetected IPSs: {len(sorf_ipss)} in {time_min:.2f} min')

        sorf_pscs_psccs = []
        if(len(sorfs_not_found) > 0):
            if(cfg.db_info['type'] == 'full'):
                log.debug('search sORF PSC')
                start_time = time.perf_counter()
                cdss_not_found_tmp, sorfs_not_found = s_orf.search_pscs(sorfs_not_found)
                sorf_pscs_psccs.extend(cdss_not_found_tmp)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\tfound PSCs: {len(sorf_pscs_psccs)} in {time_min:.2f} min')
            else:
                log.debug('search sORF PSCC')
                start_time = time.perf_counter()
                sorf_psccs, sorfs_not_found = s_orf.search_psccs(sorfs_not_found)
                sorf_pscs_psccs.extend(sorf_psccs)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\tfound PSCCs: {len(sorf_pscs_psccs)} in {time_min:.2f} min')

        print("\tlookup annotations...")

        # Copy `bakta.db` here, if tmp_db_path is set
        if cfg.tmp_db_path:
            bu.rsync_copy(f'{cfg.db_path}/bakta.db', cfg.tmp_db_path)

        log.debug('lookup sORF PSCs')
        start_time = time.perf_counter()
        sorf_pscs_psccs.extend(sorf_ipss)
        psc.lookup(sorf_pscs_psccs)  # lookup PSC info
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tlooked up PSCs in {time_min:.2f} min')

        log.debug('lookup sORF PSCCs')
        start_time = time.perf_counter()
        pscc.lookup(sorf_pscs_psccs)  # lookup PSC info
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tlooked up PSCCs in {time_min:.2f} min')

        # Immediately clean up the tmp db, if set to maximize space
        if cfg.tmp_db_path:
            Path(cfg.tmp_db_path.joinpath('bakta.db')).unlink()

        print('\tfilter and combine annotations...')
        log.debug('filter sORF by annotations')
        start_time = time.perf_counter()
        sorfs_filtered = s_orf.annotation_filter(sorfs)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tfiltered sORFs in {time_min:.2f} min')

        log.debug('combine sORF annotations')
        start_time = time.perf_counter()
        for feat in sorfs_filtered:
            anno.combine_annotation(feat)  # combine IPS and PSC annotations
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tcombined annotations in {time_min:.2f} min')

        data['features'].extend(sorfs_filtered)
        print(f'\tfiltered sORFs: {len(sorfs_filtered)}')

    ############################################################################
    # gap annotation
    # - in-mem gap detection
    # - gap annotation
    ############################################################################
    if(cfg.skip_gap):
        print('skip gap annotation...')
    else:
        print('detect gaps...')
        log.debug('detect gaps')
        start_time = time.perf_counter()
        assembly_gaps = gaps.detect_assembly_gaps(data)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        data['features'].extend(assembly_gaps)
        print(f'\tfound: {len(assembly_gaps)} in {time_min:.2f} min')

    ############################################################################
    # oriC/T prediction
    ############################################################################
    if(cfg.skip_ori):
        print('skip oriC/T annotation...')
    else:
        print('detect oriCs/oriVs...')
        log.debug('detect oriC/V')
        start_time = time.perf_counter()
        oriCs = ori.predict_oris(data, sequences_path, bc.FEATURE_ORIC)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        data['features'].extend(oriCs)
        print(f'\tfound: {len(oriCs)} in {time_min:.2f} min')

        print('detect oriTs...')
        log.debug('detect oriT')
        start_time = time.perf_counter()
        oriTs = ori.predict_oris(data, sequences_path, bc.FEATURE_ORIT)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        data['features'].extend(oriTs)
        print(f'\tfound: {len(oriTs)} in {time_min:.2f} min')

    ############################################################################
    # Filter overlapping features
    ############################################################################
    if(cfg.skip_filter):
        print('skip feature overlap filters...')
    else:
        print('apply feature overlap filters...')
        start_time = time.perf_counter()
        anno.detect_feature_overlaps(data)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'Feature overlap filtering completed in {time_min:.2f} min')

    ############################################################################
    # Create annotations
    # - filter features based on precedence and overlaps
    # - sort features
    # - create locus tags for features
    ############################################################################
    print('select features and create locus tags...')
    log.debug('start feature selection and creation of locus tags')
    start_time = time.perf_counter()
    features_by_sequence = {seq['id']: [] for seq in data['sequences']}
    feature_id = 1
    feature_id_prefix = bu.create_locus_tag_prefix(sequences, length=10)
    for feature in data['features']:
        if('discarded' not in feature):
            feature['id'] = f'{feature_id_prefix}_{feature_id}'
            feature_id += 1
            seq_features = features_by_sequence.get(feature['sequence'])
            seq_features.append(feature)
    features = []
    for seq in data['sequences']:
        seq_features = features_by_sequence[seq['id']]
        seq_features.sort(key=lambda k: k['start'])
        features.extend(seq_features)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    data['features'] = features  # overwrite feature list by final sorted feature list
    log.info('selected features=%i', len(features))
    print(f'\tselected: {len(features)} in {time_min:.2f} min')

    locus_tag_prefix = cfg.locus_tag if cfg.locus_tag else bu.create_locus_tag_prefix(sequences)
    log.info('locus tag prefix=%s', locus_tag_prefix)
    start_time = time.perf_counter()
    locus_tag_nr = cfg.locus_tag_increment
    for feature in features:
        locus_tag = f'{locus_tag_prefix}_{locus_tag_nr:0{len(str(cfg.locus_tag_increment*len(features)))+1}}'
        if(feature['type'] in [bc.FEATURE_T_RNA, bc.FEATURE_TM_RNA, bc.FEATURE_R_RNA, bc.FEATURE_NC_RNA, bc.FEATURE_CDS, bc.FEATURE_SORF]):
            feature['locus'] = locus_tag
            locus_tag_nr += cfg.locus_tag_increment
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'Locus tags created in {time_min:.2f} min')

    ############################################################################
    # Improve annotations
    # - select CDS/sORF gene symbols based on adjacent genes
    ############################################################################
    print('improve annotations...')
    start_time = time.perf_counter()
    genes_with_improved_symbols = anno.select_gene_symbols([feature for feature in features if feature['type'] in [bc.FEATURE_CDS, bc.FEATURE_SORF]])
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\trevised gene symbols: {len(genes_with_improved_symbols)} in {time_min:.2f} min')

    ############################################################################
    # Print summary
    # - genome stats
    # - annotation stats
    ############################################################################
    start_time = time.perf_counter()
    bu.calc_genome_stats(data)
    print('\nGenome statistics:')
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'Genome statistics calculated in {time_min:.2f} min')

    print(f"\tGenome size: {data['stats']['size']:,} bp")
    print(f"\tContigs/replicons: {len(data['sequences'])}")
    print(f"\tGC: {100 * data['stats']['gc']:.1f} %")
    print(f"\tN50: {data['stats']['n50']:,}")
    print(f"\tN90: {data['stats']['n90']:,}")
    print(f"\tN ratio: {100 * data['stats']['n_ratio']:.1f} %")
    print(f"\tcoding density: {100 * data['stats']['coding_ratio']:.1f} %")
    print('\nannotation summary:')
    print(f"\ttRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_T_RNA])}")
    print(f"\ttmRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_TM_RNA])}")
    print(f"\trRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_R_RNA])}")
    print(f"\tncRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA])}")
    print(f"\tncRNA regions: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA_REGION])}")
    print(f"\tCRISPR arrays: {len([feat for feat in features if feat['type'] == bc.FEATURE_CRISPR])}")
    cdss = [feat for feat in features if feat['type'] == bc.FEATURE_CDS]
    print(f"\tCDSs: {len(cdss)}")
    print(f"\t\thypotheticals: {len([cds for cds in cdss if 'hypothetical' in cds])}")
    print(f"\t\tpseudogenes: {len([cds for cds in cdss if 'pseudogene' in cds])}")
    print(f"\tsORFs: {len([feat for feat in features if feat['type'] == bc.FEATURE_SORF])}")
    print(f"\tgaps: {len([feat for feat in features if feat['type'] == bc.FEATURE_GAP])}")
    print(f"\toriCs/oriVs: {len([feat for feat in features if (feat['type'] == bc.FEATURE_ORIC or feat['type'] == bc.FEATURE_ORIV)])}")
    print(f"\toriTs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIT])}")

    ############################################################################
    # Write output files
    # - measure runtime
    # - write optional output files in GFF3/GenBank/EMBL formats
    # - write comprehensive annotation results as JSON
    # - remove temp directory
    ############################################################################
    cfg.run_end = datetime.now()  # measure runtime

    print(f'\nExport annotation results to: {cfg.output_path}')
    print('\thuman readable TSV...')
    start_time = time.perf_counter()
    tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.tsv')
    tsv.write_features(data['sequences'], features_by_sequence, tsv_path)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tTSV written in {time_min:.2f} min')

    print('\tGFF3...')
    start_time = time.perf_counter()
    gff3_path = cfg.output_path.joinpath(f'{cfg.prefix}.gff3')
    gff.write_features(data, features_by_sequence, gff3_path)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tGFF3 written in {time_min:.2f} min')

    if(cfg.skip_write_genbank_embl):
        print('\tskip generation of GenBank & EMBL format files...')
    else:
        print('\tINSDC GenBank & EMBL...')
        start_time = time.perf_counter()
        genbank_path = cfg.output_path.joinpath(f'{cfg.prefix}.gbff')
        embl_path = cfg.output_path.joinpath(f'{cfg.prefix}.embl')
        insdc.write_features(data, features, genbank_path, embl_path)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tGenBank & EMBL written in {time_min:.2f} min')

    print('\tgenome sequences...')
    start_time = time.perf_counter()
    fna_path = cfg.output_path.joinpath(f'{cfg.prefix}.fna')
    fasta.export_sequences(data['sequences'], fna_path, description=True, wrap=True)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tGenome sequences written in {time_min:.2f} min')

    print('\tfeature nucleotide sequences...')
    start_time = time.perf_counter()
    ffn_path = cfg.output_path.joinpath(f'{cfg.prefix}.ffn')
    fasta.write_ffn(features, ffn_path)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tFeature nucleotide sequences written in {time_min:.2f} min')

    print('\ttranslated CDS sequences...')
    start_time = time.perf_counter()
    faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.faa')
    fasta.write_faa(features, faa_path)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tTranslated CDS sequences written in {time_min:.2f} min')

    print('\tfeature inferences...')
    tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.inference.tsv')
    tsv.write_feature_inferences(data['sequences'], features_by_sequence, tsv_path)

    if(cfg.skip_plot  or  cfg.meta):
        print('\tskip generation of circular genome plot...')
    else:
        print('\tcircular genome plot...')
        start_time = time.perf_counter()
        plot.write(data, features, cfg.output_path)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tCircular genome plot generated in {time_min:.2f} min')

    if(cfg.skip_cds is False):
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]
        print('\thypothetical TSV...')
        start_time = time.perf_counter()
        tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.tsv')
        tsv.write_hypotheticals(hypotheticals, tsv_path)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tHypothetical TSV written in {time_min:.2f} min')

        print('\ttranslated hypothetical CDS sequences...')
        start_time = time.perf_counter()
        faa_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.faa')
        fasta.write_faa(hypotheticals, faa_path)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tTranslated hypothetical CDS sequences written in {time_min:.2f} min')

    # calc & store runtime
    run_duration = (cfg.run_end - cfg.run_start).total_seconds()
    data['run'] = {
        'start': cfg.run_start.strftime('%Y-%m-%d %H:%M:%S'),
        'end': cfg.run_end.strftime('%Y-%m-%d %H:%M:%S'),
        'duration': f'{(run_duration / 60):.2f} min'
    }

    print('\tmachine readable JSON...')
    start_time = time.perf_counter()
    json_path = cfg.output_path.joinpath(f'{cfg.prefix}.json')
    json.write_json(data, features, json_path)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tJSON written in {time_min:.2f} min')

    print('\tGenome and annotation summary...')
    start_time = time.perf_counter()
    summary_path = cfg.output_path.joinpath(f'{cfg.prefix}.txt')
    with summary_path.open('w') as fh_out:
        fh_out.write('Sequence(s):\n')
        fh_out.write(f"Length: {data['stats']['size']:}\n")
        fh_out.write(f"Count: {len(data['sequences'])}\n")
        fh_out.write(f"GC: {100 * data['stats']['gc']:.1f}\n")
        fh_out.write(f"N50: {data['stats']['n50']:}\n")
        fh_out.write(f"N90: {data['stats']['n90']:}\n")
        fh_out.write(f"N ratio: {100 * data['stats']['n_ratio']:.1f}\n")
        fh_out.write(f"coding density: {100 * data['stats']['coding_ratio']:.1f}\n")
        fh_out.write('\nAnnotation:\n')
        fh_out.write(f"tRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_T_RNA])}\n")
        fh_out.write(f"tmRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_TM_RNA])}\n")
        fh_out.write(f"rRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_R_RNA])}\n")
        fh_out.write(f"ncRNAs: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA])}\n")
        fh_out.write(f"ncRNA regions: {len([feat for feat in features if feat['type'] == bc.FEATURE_NC_RNA_REGION])}\n")
        fh_out.write(f"CRISPR arrays: {len([feat for feat in features if feat['type'] == bc.FEATURE_CRISPR])}\n")
        fh_out.write(f"CDSs: {len(cdss)}\n")
        fh_out.write(f"pseudogenes: {len([cds for cds in cdss if 'pseudogene' in cds])}\n")
        fh_out.write(f"hypotheticals: {len([cds for cds in cdss if 'hypothetical' in cds])}\n")
        fh_out.write(f"sORFs: {len([feat for feat in features if feat['type'] == bc.FEATURE_SORF])}\n")
        fh_out.write(f"gaps: {len([feat for feat in features if feat['type'] == bc.FEATURE_GAP])}\n")
        fh_out.write(f"oriCs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIC])}\n")
        fh_out.write(f"oriVs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIV])}\n")
        fh_out.write(f"oriTs: {len([feat for feat in features if feat['type'] == bc.FEATURE_ORIT])}\n")
        fh_out.write('\nBakta:\n')
        fh_out.write(f'Software: v{cfg.version}\n')
        fh_out.write(f"Database: v{cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}\n")
        fh_out.write('DOI: 10.1099/mgen.0.000685\n')
        fh_out.write('URL: github.com/oschwengers/bakta\n')
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tSummary written in {time_min:.2f} min')

    print(f'\nIf you use these results please cite Bakta: https://doi.org/{bc.BAKTA_DOI}')
    print(f'Annotation successfully finished in {int(run_duration / 60):01}:{int(run_duration % 60):02} [mm:ss].')


if __name__ == '__main__':
    main()
