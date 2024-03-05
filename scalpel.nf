#!/usr/bin/env nextflow
/*
 * Author: PLASS lab - Franz AKE
 * Barcelona, SPAIN
 */


/* Main scalpel script for characterization of alternative polyadenylation */
nextflow.enable.dsl=2

/* Define default params variables */
/* ******************************* */
params.barcodes = null
params.clusters = null
params.transcriptome = null
params.reads = null
params.dt_threshold = 1200
params.dt_exon_end_threshold = 30
params.isoform_end_ip_threshold = 60
params.gene_fraction = "98%"
params.binsize = 20
params.cpus = 10
params.threads = 10
params.cellranger_repo = true



log.info """\
    ===============================
    SCALPEL - N F   P I P E L I N E
    ===============================
    Author: PLASS Lab
    *****************
    P-CMRC - Barcelona, SPAIN

    input files:
    - Annotation required files(required):
        - transcriptome reference [--transcriptome]: ${params.transcriptome}
        - annotation GTF reference [--gtf]: ${params.gtf}
        - internal priming annotation [--ipdb]: ${params.ipdb}
       

    - Reads processing files (required):
        - samples files [--samples]: ${params.samples}
        - fastqs files [--reads]: ${params.reads}

    - Params:
        Required:
        - sequencing type (required): ${params.sequencing}

        Optional:
        - barcodes [--barcodes] (optional): ${params.barcodes}
        - clusters of cells [--clusters] (optional): ${params.clusters}
        - transcriptomic distance threshold [--dt_threshold] (optional, default 600bp): ${params.dt_threshold}
        - transcriptomic end distance threhsold [--dt_exon_end_threshold] (optional, default 30bp): ${params.dt_exon_end_threshold}
        - minimal distance of Ip from isoform 3'ends (optional, 60): ${params.isoform_end_ip_threshold}
        - params.threads [--threads] (default 10): ${params.threads}
        - params.cpus [--cpus] (default 10): ${params.cpus}

    
    - Results:
        - "./results"

    """.stripIndent()


/* Chech required args */
/* ******************* */
if ( params.transcriptome == null | params.gtf == null | params.ipdb == null)
	error "Enter annotation reference files ! [--transcriptome / --gtf / --ipdb]"

if ( params.sequencing == null )
	error "Enter sequencing argument (chromium or dropseq) !  [--sequencing]"

if ( params.samples == null | params.reads == null)
	error "Enter sample files ! [--samples / --reads]"



/* Include subworkflow modules */
/* *************************** */
include { get_10X_sampleIDs; read_10Xrepo } from './src/loading_10X_files.nf'
include { salmon_indexing; salmon_bulk_quantification; tpm_bulk_average; gtf_splitting_bedfile_conversion; isoform_selection_weighting } from './src/annotation_processing.nf'
include { bam_splitting; bedfile_conversion; reads_mapping_filtering } from './src/reads_processing.nf'
include { ip_splitting; ip_filtering } from './src/internalp_filtering.nf'
include { probability_distribution; fragment_probabilities; cells_splitting; em_algorithm; cells_merging; dge_generation } from './src/isoform_quantification.nf'
include { merge_extracted_readsids; extraction_readids; bam_filtering; merging_results;  } from './src/publishing_results.nf'
include { differential_isoform_usage } from './src/downstream_analysis.nf'



/* Preprocessing of annotation files */
/* ********************************* */
workflow annotation_preprocessing {

    take:
        transcriptome_file
        gtf_file
        read_pairs

    main:
        salmon_indexing(transcriptome_file)
        salmon_bulk_quantification(salmon_indexing.out, gtf_file, read_pairs)
        tpm_bulk_average(salmon_bulk_quantification.out.collect(), gtf_file)
        gtf_splitting_bedfile_conversion(gtf_file)
        isoform_selection_weighting(tpm_bulk_average.out, gtf_splitting_bedfile_conversion.out.flatten())
    
    emit:
        selected_isoforms = isoform_selection_weighting.out
}


/* Preprocessing of annotation files */
/* ********************************* */
workflow chromium_repo_processing {

    take:
        repository_path

    main:
        get_10X_sampleIDs(repository_path)

        /* extract sampleIDs and associated paths */
        get_10X_sampleIDs.out.flatten().combine(Channel.from(repository_path)).set{ sample_ids }
        read_10Xrepo(sample_ids)

        /* formatting */
        read_10Xrepo.out.map{ it = [it[0], [it[1], it[2], it[3], it[4]]] }.set{ sample_ids_ch }

    emit:
        sample_links_extracted = sample_ids_ch
}

/* Processing of reads */
/* ******************* */
workflow reads_processing {

    take:
        selected_isoforms
        sample_links

    main:
        bam_splitting(selected_isoforms.flatMap { it = it[0] }.combine(sample_links.map{ sample_id, paths -> tuple( sample_id, paths[0,1,2]) }))
        bedfile_conversion(bam_splitting.out)

        /*format isoform channel*/
        selected_isoforms = sample_links.flatMap{it = it[0]}.combine(selected_isoforms)

        /** merging and mapping*/
        reads_mapping_filtering(bedfile_conversion.out.join(selected_isoforms, by: [0,1]))
    
    emit:
        splitted_bams = bam_splitting.out
        mappeds_reads = reads_mapping_filtering.out
}


/* Internal priming filtering */
/* ************************** */
workflow internalpriming_filtering {

    take:
        mappeds_reads

    main:
        iptargets = mappeds_reads.flatMap{ it = [it[1]]}.unique().combine(Channel.fromPath(params.ipdb))
        iptargets = mappeds_reads.flatMap{ it = [it[0]]}.unique().combine(ip_splitting(iptargets))        
        ip_filtering(mappeds_reads.join(iptargets, by:[0,1]))

    emit:
        filtered_reads = ip_filtering.out
}

/* Isoform quantification */
/* ************************** */
workflow isoform_quantification {

    take:
        filtered_reads
        raw_dge

    main:
        /* Calculate fragments transcriptomic porbabilities */
        filtered_reads.flatMap{ it = [it[0,1,3]]}.set{ unique_reads }
        unique_reads.map{ sample_id, chr, reads -> tuple( sample_id, [chr, reads]) }.groupTuple(by: 0).map{ sample_id, files -> tuple( sample_id, files.flatten() )}.set{ unique_reads }

        /* calculate probabilities for each sample */
        probability_distribution(unique_reads)
        probs = probability_distribution.out.map{ sample_id, prob_count, prob_figure -> tuple ( sample_id, prob_count) }

        /* calculate all probabililities */
        filtered_reads.flatMap{ it = [it[0,1,5]]}.join(filtered_reads.flatMap{ it = it[1]}.unique().combine(probs).flatMap{ it = [it[1,0,2]]}, by:[0,1]).set{ all_reads }
        fragment_probabilities(all_reads)
        cells_splitting(fragment_probabilities.out.groupTuple(by: 0))
        // cells_splitting.out.transpose().view()

        /* perform em algorithm */
        em_algorithm(cells_splitting.out.transpose())
        cells_merging(em_algorithm.out.groupTuple(by: 0))

        /* DGE generation */
        dge_generation(cells_merging.out.join(raw_dge))

    emit:
        dges = dge_generation.out
}


/* Results */
/* ******* */
workflow results {

    take:
        splitted_bams
        filtered_reads
    main:
        /* joining */
        filtered_reads.flatMap{ it = [it[0,1,2,4]]}.set{ selected_reads }
        splitted_bams.join(selected_reads, by: [0,1]).set{ selected_reads }
        filtered_reads.flatMap{ it = [it[0,1,5]]}.set{ all_filtered_reads }

        /* Extract readID/gene_name/transcript_name */
        extraction_readids(all_filtered_reads)
        extraction_readids.out.groupTuple(by: 0).map{ sample_id, files -> tuple(sample_id,files)}.set{filtered_reads2}
        merge_extracted_readsids(filtered_reads2)

        /* Filter bamfile */
        bam_filtering(selected_reads)

        /* merge and write */
        merging_results(bam_filtering.out.groupTuple(by: 0))
}


/* Analysis */
/* ******** */
workflow downstream_analysis {

    take:
        dges
    main:
        /* merging seurat objects */
        dges.flatMap{ it = [it[2]] }.collect().set{ seurat_objs }
        differential_isoform_usage(seurat_objs)

}


/* Main Entrypoint */
/* *************** */

workflow {

    /* Input files */
    /* Input fastq reads */
    read_pairs_ch = Channel.fromFilePairs( "${params.reads}/*{1,2}_001.fastq.gz", checkIfExists: true )

    /* annotation preprocessing (1) */
    annotation_preprocessing(params.transcriptome, params.gtf, read_pairs_ch)

    /* reads_processing - reading (1) */
    /* Chromium */
    if ( params.sequencing == "chromium" && params.cellranger_repo == true ) {

        /* Parsing of 10X repository */
        chromium_repo_processing("${params.samples}")
        chromium_repo_processing.out.sample_links_extracted.set{ sample_links_ch }

    } else {

        /* It is important that the samples extension get settled as *.bam / *.bai / *.barcodes / *.counts */
        Channel.fromFilePairs( "${params.samples}/*{.bam,.bam.bai,.barcodes.txt,.counts.txt}", size: 4, checkIfExists: true ).set{ sample_links_ch }

    }
    
    /* reads_processing - extracting (2) */   
    reads_processing(annotation_preprocessing.out.selected_isoforms, sample_links_ch)

    /* internalp filtering processing (3) */
    internalpriming_filtering(reads_processing.out.mappeds_reads)

    /* isoform quantification (4) */
    isoform_quantification(internalpriming_filtering.out, sample_links_ch.map{ sample_id, paths -> tuple( sample_id, paths[3])})

    /* results (5) */
    results(reads_processing.out.splitted_bams, internalpriming_filtering.out.filtered_reads)

    /* analysis (6) */
    if ( params.clusters != null )
        downstream_analysis(isoform_quantification.out)

}

