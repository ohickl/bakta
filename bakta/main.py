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
import bakta.features.annotation as anno
import bakta.features.t_rna as t_rna
import bakta.features.tm_rna as tm_rna
import bakta.features.r_rna as r_rna
import bakta.features.nc_rna as nc_rna
import bakta.features.nc_rna_region as nc_rna_region
import bakta.features.crispr as crispr
import bakta.features.orf as orf
import bakta.features.cds as feat_cds
import bakta.features.signal_peptides as sig_peptides
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
        print(f'Bakta v{bakta.__version__}')
        print('Options and arguments:')
        print(f'\tinput: {cfg.genome_path}')
        print(f"\tdb: {cfg.db_path}, version {cfg.db_info['major']}.{cfg.db_info['minor']}, {cfg.db_info['type']}")
        if(cfg.user_proteins): print(f'\tuser proteins: {cfg.user_proteins}')
        if(cfg.replicons): print(f'\treplicon table: {cfg.replicons}')
        if(cfg.prodigal_tf): print(f'\tprodigal training file: {cfg.prodigal_tf}')
        print(f'\toutput: {cfg.output_path}')
        if(cfg.force): print(f'\tforce: {cfg.force}')
        print(f'\ttmp directory: {cfg.tmp_path}')
        print(f'\tprefix: {cfg.prefix}')
        print(f'\tthreads: {cfg.threads}')
        if(cfg.debug): print(f'\tdebug: {cfg.debug}')
        if(cfg.meta): print(f'\tmeta mode: {cfg.meta}')
        print(f'\ttranslation table: {cfg.translation_table}')
        if(cfg.taxon): print(f'\ttaxon: {cfg.taxon}')
        if(cfg.plasmid): print(f'\tplasmid: {cfg.plasmid}')
        if(cfg.gram != '?'): print(f'\tgram: {cfg.gram}')
        if(cfg.locus): print(f'\tlocus prefix: {cfg.locus}')
        if(cfg.locus_tag): print(f'\tlocus tag prefix: {cfg.locus_tag}')
        if(cfg.complete): print(f'\tcomplete replicons: {cfg.complete}')
        if(cfg.compliant): print(f'\tINSDC compliant: {cfg.compliant}')
        if(cfg.keep_contig_headers): print(f'\tkeep contig headers: {cfg.keep_contig_headers}')
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
        if(cfg.skip_plot): print(f'\tskip plot: {cfg.skip_plot}')
        if(cfg.skip_write_genbank_embl): print(f'\tskip write GenBank/EMBL: {cfg.skip_write_genbank_embl}')
    
    if(cfg.debug):
        print(f"\nBakta runs in DEBUG mode! Temporary data will not be destroyed at: {cfg.tmp_path}")
    else:
        atexit.register(bu.cleanup, log, cfg.tmp_path)  # register cleanup exit hook

    ############################################################################
    # Import genome
    # - parse contigs in Fasta file
    # - apply contig length filter
    # - rename contigs
    ############################################################################
    print('\nparse genome sequences...')
    try:
        contigs = fasta.import_contigs(cfg.genome_path)
        log.info('imported sequences=%i', len(contigs))
        print(f'\timported: {len(contigs)}')
    except:
        log.error('wrong genome file format!', exc_info=True)
        sys.exit('ERROR: wrong genome file format!')
    replicons = bu.parse_replicon_table(cfg.replicons) if cfg.replicons else None
    contigs, complete_genome = bu.qc_contigs(contigs, replicons)
    print(f'\tfiltered & revised: {len(contigs)}')
    no_chromosomes = len([c for c in contigs if c['type'] == bc.REPLICON_CHROMOSOME])
    if(no_chromosomes > 0):
        print(f"\tchromosomes: {no_chromosomes}")
    no_plasmids = len([c for c in contigs if c['type'] == bc.REPLICON_PLASMID])
    if(no_plasmids > 0):
        print(f"\tplasmids: {no_plasmids}")
    no_contigs = len([c for c in contigs if c['type'] == bc.REPLICON_CONTIG])
    if(no_contigs > 0):
        print(f"\tcontigs: {no_contigs}")
    if(len(contigs) == 0):
        log.warning('no valid contigs!')
        sys.exit('Error: input file contains no valid contigs.')
    contigs_path = cfg.tmp_path.joinpath('contigs.fna')
    fasta.export_contigs(contigs, contigs_path)
    genome = {
        'genus': cfg.genus,
        'species': cfg.species,
        'strain': cfg.strain,
        'taxon': cfg.taxon,
        'gram': cfg.gram,
        'translation_table': cfg.translation_table,
        'size': sum([c['length'] for c in contigs]),
        'complete': cfg.complete or complete_genome,
        'features': {},
        'contigs': contigs
    }
    if(cfg.plasmid):
        genome['plasmid'] = cfg.plasmid
    print('\nstart annotation...')


    # Split assembly into chuncks for parallel processing
    print(f'split assembly into {cfg.threads} chunks...')
    fasta_chunk_dir = cfg.tmp_path.joinpath('assembly_chunks')
    fasta_chunk_dir.mkdir(parents=True, exist_ok=True)

    # Split the fasta file
    print(f'\tsplitting assembly into {cfg.threads} chunks...')
    start_time = time.perf_counter()
    split_cmd = [
        'seqkit', 'split',
        '--quiet',
        '-p', str(cfg.threads),
        '-O', str(fasta_chunk_dir),
        str(contigs_path)
    ]
    sp.run(split_cmd, check=True)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\ttime: {time_min:.2f} min')

    # Determine the extension of the input fasta file and generate sorted chunk paths
    contig_fasta_ext = contigs_path.suffix
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
        genome['features'][bc.FEATURE_T_RNA] = t_rna.predict_t_rnas(genome, fasta_chunk_paths)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(genome['features'][bc.FEATURE_T_RNA])} in {time_min:.2f} min")

    ############################################################################
    # tmRNA prediction
    ############################################################################
    if(cfg.skip_tmrna):
        print('skip tmRNA prediction...')
    else:
        print('predict tmRNAs...')
        start_time = time.perf_counter()
        log.debug('start tmRNA prediction')
        genome['features'][bc.FEATURE_TM_RNA] = tm_rna.predict_tm_rnas(genome, fasta_chunk_paths)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(genome['features'][bc.FEATURE_TM_RNA])} in {time_min:.2f} min")

    ############################################################################
    # rRNA prediction
    ############################################################################
    if(cfg.skip_rrna):
        print('skip rRNA prediction...')
    else:
        # Copy rRNA database to tmp db directory, if set
        if cfg.tmp_db_path:
            print('copy rRNA database to tmp db directory...')
            start_time = time.perf_counter()
            bu.rsync_copy(f'{cfg.db_path}/rRNA*', cfg.tmp_db_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

        print('predict rRNAs...')
        start_time = time.perf_counter()
        log.debug('start rRNA prediction')
        genome['features'][bc.FEATURE_R_RNA] = r_rna.predict_r_rnas(genome, fasta_chunk_paths)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(genome['features'][bc.FEATURE_R_RNA])} in {time_min:.2f} min")

        # Clean up tmp rRNA database, if set
        if cfg.tmp_db_path:
            print('remove rRNA database from tmp db directory...')
            start_time = time.perf_counter()
            for file in Path(cfg.tmp_db_path).iterdir():
                file.unlink()
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

    ############################################################################
    # ncRNA gene prediction
    ############################################################################
    if(cfg.skip_ncrna):
        print('skip ncRNA prediction...')
    else:
        # Copy ncRNA database to tmp db directory, if set
        if cfg.tmp_db_path:
            print('copy ncRNA database to tmp db directory...')
            start_time = time.perf_counter()
            bu.rsync_copy(f'{cfg.db_path}/ncRNA-genes*', cfg.tmp_db_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')
            
        print('predict ncRNAs...')
        start_time = time.perf_counter()
        log.debug('start ncRNA prediction')
        genome['features'][bc.FEATURE_NC_RNA] = nc_rna.predict_nc_rnas(genome, fasta_chunk_paths)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(genome['features'][bc.FEATURE_NC_RNA])} in {time_min:.2f} min")

        # Clean up tmp ncRNA database, if set
        if cfg.tmp_db_path:
            print('remove ncRNA database from tmp db directory...')
            start_time = time.perf_counter()
            for file in Path(cfg.tmp_db_path).iterdir():
                file.unlink()
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

    ############################################################################
    # ncRNA region prediction
    ############################################################################
    if(cfg.skip_ncrna_region):
        print('skip ncRNA region prediction...')
    else:
        # Copy ncRNA region database to tmp db directory, if set
        if cfg.tmp_db_path:
            print('copy ncRNA region database to tmp db directory...')
            start_time = time.perf_counter()
            bu.rsync_copy(f'{cfg.db_path}/ncRNA-regions*', cfg.tmp_db_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')
        
        print('predict ncRNA regions...')
        start_time = time.perf_counter()
        log.debug('start ncRNA region prediction')
        genome['features'][bc.FEATURE_NC_RNA_REGION] = nc_rna_region.predict_nc_rna_regions(genome, fasta_chunk_paths)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(genome['features'][bc.FEATURE_NC_RNA_REGION])} in {time_min:.2f} min")

        # Clean up tmp ncRNA region database, if set
        if cfg.tmp_db_path:
            print('remove ncRNA region database from tmp db directory...')
            start_time = time.perf_counter()
            for file in Path(cfg.tmp_db_path).iterdir():
                file.unlink()
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

    ############################################################################
    # CRISPR prediction
    ############################################################################
    if(cfg.skip_crispr):
        print('skip CRISPR array prediction...')
    else:
        print('predict CRISPR arrays...')
        start_time = time.perf_counter()
        log.debug('start CRISPR prediction')
        genome['features'][bc.FEATURE_CRISPR] = crispr.predict_crispr(genome, fasta_chunk_paths)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tfound: {len(genome['features'][bc.FEATURE_CRISPR])} in {time_min:.2f} min")

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
        cdss = feat_cds.predict(genome)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f"\tpredicted: {len(cdss)} in {time_min:.2f} min")

        if(len(cdss) > 0):
            log.debug('detect spurious CDS')
            start_time = time.perf_counter()
            discarded_cdss = orf.detect_spurious(cdss)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\tdiscarded spurious: {len(discarded_cdss)} in {time_min:.2f} min')
            cdss = [cds for cds in cdss if 'discarded' not in cds]
        
        if(len(cdss) > 0):
            log.debug('revise translational exceptions')
            start_time = time.perf_counter()
            no_revised = feat_cds.revise_translational_exceptions(genome, cdss)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\trevised translational exceptions: {no_revised} in {time_min:.2f} min')
            cdss = [cds for cds in cdss if 'discarded' not in cds]
        
        if(cfg.regions):
            log.debug('import user-provided CDS regions')
            start_time = time.perf_counter()
            imported_cdss = feat_cds.import_user_cdss(genome, cfg.regions)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\timported CDS regions: {len(imported_cdss)} in {time_min:.2f} min')
            cdss.extend(imported_cdss)

        if(len(cdss) > 0):
            log.debug('lookup CDS UPS/IPS')
            start_time = time.perf_counter()
            cdss_ups, cdss_not_found = ups.lookup(cdss)
            cdss_ips, sorf_pscs = ips.lookup(cdss_ups)
            cdss_not_found.extend(sorf_pscs)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\tdetected IPSs: {len(cdss_ips)} in {time_min:.2f} min')

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
            psc.lookup(cdss)  # lookup PSC info
            pscc.lookup(cdss)  # lookup PSCC info
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
                print('copy expert protein sequences diamond db to tmp db directory...')
                start_time = time.perf_counter()
                bu.rsync_copy(f'{cfg.db_path}/expert-protein-sequences.dmnd', cfg.tmp_db_path)
                diamond_db_path = cfg.tmp_db_path.joinpath('expert-protein-sequences.dmnd')
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\ttime: {time_min:.2f} min')
            else:
                diamond_db_path = cfg.db_path.joinpath('expert-protein-sequences.dmnd')

            start_time = time.perf_counter()
            expert_aa_found = exp_aa_seq.search(cdss, cds_aa_path, 'expert_proteins', diamond_db_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\t\tprotein sequences: {len(expert_aa_found)} in {time_min:.2f} min')

            # Cleanup tmp diamond db
            if cfg.tmp_db_path:
                print('remove expert protein sequences diamond db from tmp db directory...')
                start_time = time.perf_counter()
                diamond_db_path.unlink()
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\ttime: {time_min:.2f} min')

            if(cfg.user_proteins):
                log.debug('conduct expert system: user aa seqs')
                start_time = time.perf_counter()
                user_aa_path = cfg.tmp_path.joinpath('user-proteins.faa')
                exp_aa_seq.write_user_protein_sequences(user_aa_path)
                user_aa_found = exp_aa_seq.search(cdss, cds_aa_path, 'user_proteins', user_aa_path)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\t\tuser protein sequences: {len(user_aa_found)} in {time_min:.2f} min')

            if(cfg.gram != bc.GRAM_UNKNOWN):
                start_time = time.perf_counter()
                sig_peptides_found = sig_peptides.search(cdss, cds_aa_path)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\tsignal peptides: {len(sig_peptides_found)} in {time_min:.2f} min')

            print('\tcombine annotations and mark hypotheticals...')
            log.debug('combine CDS annotations')
            start_time = time.perf_counter()
            for cds in cdss:
                anno.combine_annotation(cds)  # combine IPS & PSC annotations and mark hypotheticals
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

            hypotheticals = [cds for cds in cdss if 'hypothetical' in cds and 'edge' not in cds and cds.get('start_type', 'Edge') != 'Edge']
            if(len(hypotheticals) > 0  and  not cfg.skip_pseudo  and  cfg.db_info['type'] == 'full'):
                print('\tdetect pseudogenes...')
                start_time = time.perf_counter()
                log.debug('search pseudogene candidates')
                pseudo_candidates = feat_cds.predict_pseudo_candidates(hypotheticals)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\t\tpseudogene candidates: {len(pseudo_candidates)} in {time_min:.2f} min')

                start_time = time.perf_counter()
                pseudogenes = feat_cds.detect_pseudogenes(pseudo_candidates, cdss, genome) if len(pseudo_candidates) > 0 else []
                psc.lookup(pseudogenes, pseudo=True)
                pscc.lookup(pseudogenes, pseudo=True)
                for pseudogene in pseudogenes:
                    anno.combine_annotation(pseudogene)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f'\t\tfound pseudogenes: {len(pseudogenes)} in {time_min:.2f} min')
            
            hypotheticals = [cds for cds in cdss if 'hypothetical' in cds]
            if(len(hypotheticals) > 0):
                log.debug('analyze hypotheticals')
                start_time = time.perf_counter()
                print(f'analyze hypothetical proteins: {len(hypotheticals)}')
                pfam_hits = feat_cds.predict_pfam(hypotheticals)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print(f"\tdetected Pfam hits: {len(pfam_hits)} in {time_min:.2f} min")
                
                start_time = time.perf_counter()
                feat_cds.analyze_proteins(hypotheticals)
                end_time = time.perf_counter()
                time_min = (end_time - start_time) / 60
                print('\tcalculated proteins statistics in {time_min:.2f} min')
            
            print('\trevise special cases...')
            start_time = time.perf_counter()
            feat_cds.revise_special_cases_annotated(genome, cdss)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f'\ttime: {time_min:.2f} min')

        genome['features'][bc.FEATURE_CDS] = cdss

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
        print('extract sORF...')
        log.debug('predict sORF')
        start_time = time.perf_counter()
        sorfs = s_orf.extract(genome)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tpotential: {len(sorfs)} in {time_min:.2f} min')

        log.debug('apply sORF overlap filter')
        start_time = time.perf_counter()
        sorfs, discarded_sorfs = s_orf.overlap_filter(genome, sorfs)
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
        sorf_upss, sorfs_not_found = ups.lookup(sorfs)
        sorf_ipss, tmp = ips.lookup(sorf_upss)
        sorfs_not_found.extend(tmp)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tdetected IPSs: {len(sorf_ipss)} in {time_min:.2f} min')

        sorf_pscs_psccs = []
        if(len(sorfs_not_found) > 0):
            if(cfg.db_info['type'] == 'full'):
                log.debug('search sORF PSC')
                start_time = time.perf_counter()
                sorf_pscs, sorfs_not_found = s_orf.search_pscs(sorfs_not_found)
                sorf_pscs_psccs.extend(sorf_pscs)
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

        genome['features'][bc.FEATURE_SORF] = sorfs_filtered
        print(f'\tfiltered sORFs: {len(sorfs_filtered)}')
        
        if(cfg.gram != bc.GRAM_UNKNOWN  and  len(sorfs_filtered) > 0):
            sorf_aa_path = cfg.tmp_path.joinpath('sorfs.faa')
            with sorf_aa_path.open(mode='wt') as fh:
                for sorf in sorfs_filtered:
                    fh.write(f">{sorf['aa_hexdigest']}-{sorf['contig']}-{sorf['start']}\n{sorf['aa']}\n")
            start_time = time.perf_counter()
            sig_peptides_found = sig_peptides.search(sorfs_filtered, sorf_aa_path)
            end_time = time.perf_counter()
            time_min = (end_time - start_time) / 60
            print(f"\tsignal peptides: {len(sig_peptides_found)} in {time_min:.2f} min")

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
        assembly_gaps = gaps.detect_assembly_gaps(genome)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        genome['features'][bc.FEATURE_GAP] = assembly_gaps
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
        oriCs = ori.predict_oris(genome, contigs_path, bc.FEATURE_ORIC)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        genome['features'][bc.FEATURE_ORIC] = oriCs
        print(f'\tfound: {len(oriCs)} in {time_min:.2f} min')

        print('detect oriTs...')
        log.debug('detect oriT')
        start_time = time.perf_counter()
        oriTs = ori.predict_oris(genome, contigs_path, bc.FEATURE_ORIT)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        genome['features'][bc.FEATURE_ORIT] = oriTs
        print(f'\tfound: {len(oriTs)} in {time_min:.2f} min')

    ############################################################################
    # Filter overlapping features
    ############################################################################
    print('apply feature overlap filters...')
    start_time = time.perf_counter()
    anno.detect_feature_overlaps(genome)
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
    features_by_contig = {k['id']: [] for k in genome['contigs']}
    feature_id = 1
    feature_id_prefix = bu.create_locus_tag_prefix(contigs, length=10)
    for feature_type in [
            bc.FEATURE_T_RNA,
            bc.FEATURE_TM_RNA,
            bc.FEATURE_R_RNA,
            bc.FEATURE_NC_RNA,
            bc.FEATURE_NC_RNA_REGION,
            bc.FEATURE_CRISPR,
            bc.FEATURE_CDS,
            bc.FEATURE_SORF,
            bc.FEATURE_GAP,
            bc.FEATURE_ORIC,
            bc.FEATURE_ORIV,
            bc.FEATURE_ORIT
        ]:
        feature_list = genome['features'].get(feature_type, [])
        for feature in feature_list:
            if('discarded' not in feature):
                feature['id'] = f'{feature_id_prefix}_{feature_id}'
                feature_id += 1
                contig_features = features_by_contig.get(feature['contig'])
                contig_features.append(feature)
    features = []
    for contig in genome['contigs']:
        contig_features = features_by_contig[contig['id']]
        contig_features.sort(key=lambda k: k['start'])
        features.extend(contig_features)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    log.info('selected features=%i', len(features))
    print(f'selected: {len(features)} in {time_min:.2f} min')

    locus_tag_nr = 5
    locus_tag_prefix = cfg.locus_tag if cfg.locus_tag else bu.create_locus_tag_prefix(contigs)
    log.info('locus tag prefix=%s', locus_tag_prefix)
    start_time = time.perf_counter()
    for feature in features:
        locus_tag = f'{locus_tag_prefix}_{locus_tag_nr:05}'
        if(feature['type'] in [bc.FEATURE_T_RNA, bc.FEATURE_TM_RNA, bc.FEATURE_R_RNA, bc.FEATURE_NC_RNA, bc.FEATURE_CDS, bc.FEATURE_SORF]):
            feature['locus'] = locus_tag
            locus_tag_nr += 5
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
    print('\ngenome statistics:')
    start_time = time.perf_counter()
    genome_stats = bu.calc_genome_stats(genome, features)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'Genome statistics calculated in {time_min:.2f} min')

    print(f"\tGenome size: {genome['size']:,} bp")
    print(f"\tContigs/replicons: {len(genome['contigs'])}")
    print(f"\tGC: {100 * genome_stats['gc']:.1f} %")
    print(f"\tN50: {genome_stats['n50']:,}")
    print(f"\tN ratio: {100 * genome_stats['n_ratio']:.1f} %")
    print(f"\tcoding density: {100 * genome_stats['coding_ratio']:.1f} %")

    print('\nannotation summary:')
    print(f"\ttRNAs: {len([f for f in features if f['type'] == bc.FEATURE_T_RNA])}")
    print(f"\ttmRNAs: {len([f for f in features if f['type'] == bc.FEATURE_TM_RNA])}")
    print(f"\trRNAs: {len([f for f in features if f['type'] == bc.FEATURE_R_RNA])}")
    print(f"\tncRNAs: {len([f for f in features if f['type'] == bc.FEATURE_NC_RNA])}")
    print(f"\tncRNA regions: {len([f for f in features if f['type'] == bc.FEATURE_NC_RNA_REGION])}")
    print(f"\tCRISPR arrays: {len([f for f in features if f['type'] == bc.FEATURE_CRISPR])}")
    cdss = [f for f in features if f['type'] == bc.FEATURE_CDS]
    print(f"\tCDSs: {len(cdss)}")
    print(f"\t\thypotheticals: {len([cds for cds in cdss if 'hypothetical' in cds])}")
    print(f"\t\tpseudogenes: {len([cds for cds in cdss if 'pseudogene' in cds])}")
    print(f"\t\tsignal peptides: {len([cds for cds in cdss if bc.FEATURE_SIGNAL_PEPTIDE in cds])}")
    print(f"\tsORFs: {len([f for f in features if f['type'] == bc.FEATURE_SORF])}")
    print(f"\tgaps: {len([f for f in features if f['type'] == bc.FEATURE_GAP])}")
    print(f"\toriCs/oriVs: {len([f for f in features if (f['type'] == bc.FEATURE_ORIC or f['type'] == bc.FEATURE_ORIV)])}")
    print(f"\toriTs: {len([f for f in features if f['type'] == bc.FEATURE_ORIT])}")

    ############################################################################
    # Write output files
    # - write optional output files in GFF3/GenBank/EMBL formats
    # - measure runtime
    # - write comprehensive annotation results as JSON
    # - remove temp directory
    ############################################################################
    print(f'\nexport annotation results to: {cfg.output_path}')

    print('\thuman readable TSV...')
    start_time = time.perf_counter()
    tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.tsv')
    tsv.write_tsv(genome['contigs'], features_by_contig, tsv_path)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tTSV written in {time_min:.2f} min')

    print('\tGFF3...')
    start_time = time.perf_counter()
    gff3_path = cfg.output_path.joinpath(f'{cfg.prefix}.gff3')
    gff.write_gff3(genome, features_by_contig, gff3_path)
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
        insdc.write_insdc(genome, features, genbank_path, embl_path)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tGenBank & EMBL written in {time_min:.2f} min')

    print('\tgenome sequences...')
    start_time = time.perf_counter()
    fna_path = cfg.output_path.joinpath(f'{cfg.prefix}.fna')
    fasta.export_contigs(genome['contigs'], fna_path, description=True, wrap=True)
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

    if(cfg.skip_plot  or  cfg.meta):
        print('\tskip generation of circular genome plot...')
    else:
        print('\tcircular genome plot...')
        start_time = time.perf_counter()
        plot.write_plot(features, contigs, cfg.output_path)
        end_time = time.perf_counter()
        time_min = (end_time - start_time) / 60
        print(f'\tCircular genome plot generated in {time_min:.2f} min')

    if(cfg.skip_cds is False):
        hypotheticals = [feat for feat in features if feat['type'] == bc.FEATURE_CDS and 'hypothetical' in feat]
        print('\thypothetical TSV...')
        start_time = time.perf_counter()
        tsv_path = cfg.output_path.joinpath(f'{cfg.prefix}.hypotheticals.tsv')
        tsv.write_hypotheticals_tsv(hypotheticals, tsv_path)
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

    # measure runtime at the latest possible
    cfg.run_end = datetime.now()
    run_duration = (cfg.run_end - cfg.run_start).total_seconds()
    genome['run'] = {
        'start': cfg.run_start.strftime('%Y-%m-%d %H:%M:%S'),
        'end': cfg.run_end.strftime('%Y-%m-%d %H:%M:%S'),
        'duration': f'{(run_duration / 60):.2f} min'
    }

    print('\tmachine readable JSON...')
    start_time = time.perf_counter()
    json_path = cfg.output_path.joinpath(f'{cfg.prefix}.json')
    json.write_json(genome, features, json_path)
    end_time = time.perf_counter()
    time_min = (end_time - start_time) / 60
    print(f'\tJSON written in {time_min:.2f} min')

    print('\tgenome and annotation summary...')
    start_time = time.perf_counter()
    summary_path = cfg.output_path.joinpath(f'{cfg.prefix}.txt')
    with summary_path.open('w') as fh_out:
        fh_out.write('Sequence(s):\n')
        fh_out.write(f"Length: {genome['size']:}\n")
        fh_out.write(f"Count: {len(genome['contigs'])}\n")
        fh_out.write(f"GC: {100 * genome_stats['gc']:.1f}\n")
        fh_out.write(f"N50: {genome_stats['n50']:}\n")
        fh_out.write(f"N ratio: {100 * genome_stats['n_ratio']:.1f}\n")
        fh_out.write(f"coding density: {100 * genome_stats['coding_ratio']:.1f}\n")
        fh_out.write('\nAnnotation:\n')
        fh_out.write(f"tRNAs: {len([f for f in features if f['type'] == bc.FEATURE_T_RNA])}\n")
        fh_out.write(f"tmRNAs: {len([f for f in features if f['type'] == bc.FEATURE_TM_RNA])}\n")
        fh_out.write(f"rRNAs: {len([f for f in features if f['type'] == bc.FEATURE_R_RNA])}\n")
        fh_out.write(f"ncRNAs: {len([f for f in features if f['type'] == bc.FEATURE_NC_RNA])}\n")
        fh_out.write(f"ncRNA regions: {len([f for f in features if f['type'] == bc.FEATURE_NC_RNA_REGION])}\n")
        fh_out.write(f"CRISPR arrays: {len([f for f in features if f['type'] == bc.FEATURE_CRISPR])}\n")
        fh_out.write(f"CDSs: {len(cdss)}\n")
        fh_out.write(f"pseudogenes: {len([cds for cds in cdss if 'pseudogene' in cds])}\n")
        fh_out.write(f"hypotheticals: {len([cds for cds in cdss if 'hypothetical' in cds])}\n")
        fh_out.write(f"signal peptides: {len([cds for cds in cdss if bc.FEATURE_SIGNAL_PEPTIDE in cds])}\n")
        fh_out.write(f"sORFs: {len([f for f in features if f['type'] == bc.FEATURE_SORF])}\n")
        fh_out.write(f"gaps: {len([f for f in features if f['type'] == bc.FEATURE_GAP])}\n")
        fh_out.write(f"oriCs: {len([f for f in features if f['type'] == bc.FEATURE_ORIC])}\n")
        fh_out.write(f"oriVs: {len([f for f in features if f['type'] == bc.FEATURE_ORIV])}\n")
        fh_out.write(f"oriTs: {len([f for f in features if f['type'] == bc.FEATURE_ORIT])}\n")
        fh_out.write('\nBakta:\n')
        fh_out.write(f'Software: v{bakta.__version__}\n')
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
