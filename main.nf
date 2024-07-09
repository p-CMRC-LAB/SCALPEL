#!/usr/bin/env nextflow

/*
=================================================================================================
Author: PLASS lab - Franz AKE
Main SCALPEL script for characterization of alternative polyadenylation at single-cell resolution
Barcelona, SPAIN
=================================================================================================
*/

/* - Define SCALPEL default params variables
=============================================================================
*/
params.dt_threshold = 600
params.de_threshold = 30
params.ip_threshold = 60
params.gene_fraction = "98%"
params.binsize = 20
params.output = "./results"
params.subsample = 1

/*optional args*/
params.barcodes = null
params.clusters = null
params.help = null


/* initialize require args to null */
params.transcriptome = null
params.gtf = null
params.ipdb = null
params.samplesheet = null
params.sequencing = null

/* - Import Functions / Modules / Workflows / Subworkflows
=================================================================================================
*/
include { salmon_transcriptome_indexing; salmon_bulk_quantification; tpm_counts_average; isoform_selection_weighting } from './workflows/annotation_preprocessing.nf'
include { samples_loading; bedfile_conversion; reads_mapping_and_filtering; ip_splitting; ip_filtering } from './workflows/reads_processing.nf'
include { probability_distribution; fragment_probabilities; cells_splitting; em_algorithm; cells_merging; dge_generation } from './workflows/isoform_quantification.nf'
include { differential_isoform_usage; generation_filtered_bams } from './workflows/apa_characterization.nf'


/* In case of Help args */
/* ==================== */

if( params.help != null )
    error( """\
    ===============================
    SCALPEL - N F   P I P E L I N E
    ===============================
    Author: PLASS Lab ; Franz AKE
    *****************
    P-CMRC - Barcelona, SPAIN

    input files:
    - Annotation required files(required):
        - transcriptome reference [--transcriptome]: ${params.transcriptome}
        - annotation GTF reference [--gtf]: ${params.gtf}
        - internal priming annotation [--ipdb]: ${params.ipdb}


    - Reads processing files (required):
        - samplesheet [--samplesheet]: ${params.samplesheet}

    - Params:
        Required:
        - sequencing type (required): ${params.sequencing}

        Optional:
        - barcodes whitelist [--barcodes] (optional): ${params.barcodes}
        - cell clusters annotation [--clusters] (optional): ${params.clusters}
        - transcriptomic distance threshold [--dt_threshold] (optional, default 600bp): ${params.dt_threshold}
        - transcriptomic end distance threhsold [--de_threshold] (optional, default 30bp): ${params.de_threshold}
        - minimal distance of internal priming sites (IP) from isoform 3'ends [--ip_threshold] (optional, 60nuc): ${params.ip_threshold}
        - gene fraction abundance threshold [--gene_fraction] (optional, default '98%'): ${params.gene_fraction}
        - binsize threshold for transcriptomic distance based probability [--binsize] (optional, default '20): ${params.binsize}
        - reads subsampling threshold [--subsample] (optional, default 1): ${params.subsample}

    """.stripIndent())


/* Check required args */
/* =================== */
if( params.samplesheet==null )
    error( "Provide samplesheet path  [--samplesheet]")

if( params.transcriptome==null)
    error( "Provide transcriptome FASTA path  [--transcriptome]" )

if( params.gtf==null )
    error( "Provide GTF path  [--gtf]")

if( params.ipdb==null )
    error( "Provide internal priming annotation file  [--ipdb]" )

if ( params.sequencing==null )
    error(" Provide sequencing type (dropseq / chromium)  [--sequencing]" )



/* - Print input Params & Files information to STDOUT
=============================================================================
*/
log.info """\
    ===============================
    SCALPEL - N F   P I P E L I N E
    ===============================
    Author: PLASS Lab ; Franz AKE
    *****************
    P-CMRC - Barcelona, SPAIN

    input files:
    - Annotation required files(required):
        - transcriptome reference [--transcriptome]: ${params.transcriptome}
        - annotation GTF reference [--gtf]: ${params.gtf}
        - internal priming annotation [--ipdb]: ${params.ipdb}


    - Reads processing files (required):
        - samplesheet [--samplesheet]: ${params.samplesheet}

    - Params:
        Required:
        - sequencing type (required): ${params.sequencing}

        Optional:
        - barcodes whitelist [--barcodes] (optional): ${params.barcodes}
        - cell clusters annotation [--clusters] (optional): ${params.clusters}
        - transcriptomic distance threshold [--dt_threshold] (optional, default 600bp): ${params.dt_threshold}
        - transcriptomic end distance threhsold [--de_threshold] (optional, default 30bp): ${params.de_threshold}
        - minimal distance of internal priming sites (IP) from isoform 3'ends [--ip_threshold] (optional, 60nuc): ${params.ip_threshold}
        - gene fraction abundance threshold [--gene_fraction] (optional, default '98%'): ${params.gene_fraction}
        - binsize threshold for transcriptomic distance based probability [--binsize] (optional, default '20): ${params.binsize}
        - reads subsampling threshold [--subsample] (optional, default 1): ${params.subsample}

""".stripIndent()


    

/* - Workflows
=============================================================================
*/

workflow annotation_preprocessing {
    /* workflow for loading and processing of annotation input files */
    take:
        genome_gtf
        genome_fasta
        samples_paths

    main:
        /* Salmon transcriptome indexing */
        salmon_transcriptome_indexing(file(genome_fasta))

        /* Salmon Bulk Quantification */
        salmon_bulk_quantification(salmon_transcriptome_indexing.out, file(genome_gtf), samples_paths.map{ it= tuple(it[0], file(it[1]), file(it[2])) })

        /* Averaging of isoforms pseudobulk counts between samples */
        tpm_counts_average(salmon_bulk_quantification.out.collect(), file(genome_gtf))

        /* extract bulk quantification and chromosome gtfs */
        tpm_counts_average.out.flatMap { it = it[0] }.set { bulk_quants }
        tpm_counts_average.out.flatMap { it = it[1] }.set { gtfs }

        /* Selection of isoforms */
        isoform_selection_weighting(bulk_quants.combine(gtfs), "${params.dt_threshold}", "${params.de_threshold}")

    emit:
        selected_isoforms = isoform_selection_weighting.out

}


workflow reads_processing {
    /* workflow for loading and processing of samples input files */
    take:
        samples_paths
	selected_isoforms
        ip_annots

    main:
        /* - Loading of Samples */
        samples_loading(samples_paths, selected_isoforms)

        /* - Conversion of BAM to BED file */
        bedfile_conversion(samples_loading.out.selected_bams)

        /*format isoform channel*/
        selected_isoforms = samples_paths.flatMap{it = it[0]}.combine(selected_isoforms)

        /* merging and mapping*/
        reads_mapping_and_filtering(bedfile_conversion.out.join(selected_isoforms, by: [0,1]))

        /* Internal priming filtering of reads */
        reads_mapping_and_filtering.out.set{ mappeds_reads }
        iptargets = mappeds_reads.flatMap{ it = [it[1]]}.unique().combine(Channel.fromPath(ip_annots))
        iptargets = mappeds_reads.flatMap{ it = [it[0]]}.unique().combine(ip_splitting(iptargets))
        ip_filtering(mappeds_reads.join(iptargets, by:[0,1]), "${params.ip_threshold}")

    emit:
        splitted_bams = samples_loading.out.selected_bams.map{ it=tuple(it[0], it[2]) }
        filtered_reads = ip_filtering.out
        sample_dge = samples_loading.out.sample_files.map{ it = tuple(it[1], it[4]) }
}



workflow isoform_quantification {
    /* workflow for isoform quantification at single-cell resolution */
    take:
        filtered_reads
        samples_paths

    main:
        /* Calculate fragments transcriptomic probabilities */
        filtered_reads.flatMap{ it = [it[0,1,3]]}.set{ unique_reads }
        unique_reads.map{ sample_id, chr, reads -> tuple( sample_id, [chr, reads]) }.groupTuple(by: 0).map{ sample_id, files -> tuple( sample_id, files.flatten() )}.set{ unique_reads }

        /* calculate probabilities for each sample */
        probability_distribution(unique_reads, "${params.gene_fraction}", "${params.binsize}")
        probs = probability_distribution.out.map{ sample_id, prob_count, prob_figure -> tuple ( sample_id, prob_count) }

        /* calculate all probabililities */
        filtered_reads.flatMap{ it = [it[0,1,5]]}.join(filtered_reads.flatMap{ it = it[1]}.unique().combine(probs).flatMap{ it = [it[1,0,2]]}, by:[0,1]).set{ all_reads }
        fragment_probabilities(all_reads)
        cells_splitting(fragment_probabilities.out.groupTuple(by: 0))

        /* perform em algorithm */
        em_algorithm(cells_splitting.out.transpose())
        cells_merging(em_algorithm.out.groupTuple(by: 0))

        /* DGE generation */
        dge_generation(cells_merging.out.join(samples_paths))

    emit:
        dges = dge_generation.out
}



workflow apa_characterization {
    /* workflow for differential isoform usage analysis */
    take:
        seurat_objs
        bams

    main:
        /* seurat objects merging */
        differential_isoform_usage( seurat_objs.collect() )

        /* Merge the filtered BAM files */
        generation_filtered_bams( bams )

    emit:
        dius = differential_isoform_usage.out
}


/* - MAIN Workflow entrypoint
=============================================================================
*/

workflow {

    /* - Process samplesheet input */
    ( Channel.fromPath(params.samplesheet) | splitCsv(header:false) ).set{ samples_paths }

    /* - Annotation Preprocessing  (A)
    ============================= */
    annotation_preprocessing( "${params.gtf}", "${params.transcriptome}", samples_paths )

    /* - Reads Preprocessing (B)
    ============================= */
    reads_processing(samples_paths, annotation_preprocessing.out, "${params.ipdb}")

    /* isoform quantification (C)
    ============================= */
    isoform_quantification(reads_processing.out.filtered_reads, reads_processing.out.sample_dge)

    /* APA characterization (D)
    =========================== */
    apa_characterization( isoform_quantification.out.flatMap{ it=it[2] }, reads_processing.out.splitted_bams.groupTuple(by: 0))

}
