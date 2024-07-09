#!/usr/bin/env nextflow

/*
Loading of SAMPLES data & preprocessing according to the sequencing type
=====================================================================
*/

/* loading of subworkflows */
include { read_10Xrepo } from '../subworkflows/loading_10X_files.nf'
include { bam_splitting } from '../subworkflows/bam_splitting.nf'


/*
- In Case of processing of 10X based samples with a CellRanger repository -> Proceeding to an extraction of the required sample files
=====================================================================================================================================
*/

workflow samples_loading {
    /* Loading of samples */
    take:
        samples_paths
        selected_isoforms

    main:

        /* - In case of 10X file type, extract required files... */
        if ( "${params.sequencing}" == "chromium" ) {

            /* extract sampleIDs and associated paths */
            read_10Xrepo(samples_paths.map{ it=tuple(it[0], it[3]) })

            /* formatting */
            read_10Xrepo.out.map{ it = tuple( it[0], file(it[1]), file(it[2]), file(it[4]) ) }.set{ samples_selects }

        } else if ( "${params.sequencing}" == "dropseq") {

            samples_paths.map{ it = tuple( it[0], file(it[3]), file(it[4]), file(it[5]) ) }.set{ samples_selects }

        } else
            error( "Incoherency in [--sequencing] args !!" )

        if (params.barcodes != null) {
            /* parse barcodes file */
            ( Channel.fromPath(params.barcodes) | splitCsv(header:false) ).set{ barcodes_paths }
            (samples_selects.join(barcodes_paths, by:[0])).set{ samples_selects }
            selected_isoforms.flatMap { it = it[0] }.combine(samples_selects).set{ samples_selects }

        } else {

            selected_isoforms.flatMap { it = it[0] }.combine(samples_selects.map{ it = tuple(it[0], it[1], it[2], it[3], null) }).set{ samples_selects }

        }

        /* processing of input BAM file... */
        bam_splitting( samples_selects )

    emit:
        selected_bams = bam_splitting.out
        sample_files = samples_selects
}


process bedfile_conversion{
    tag "${sample_id}, ${chr}"
    publishDir "./results/reads_processing/bedfile_conversion/${sample_id}"
    cache true
    label "small_mem"

    input:
        tuple val(sample_id), val(chr), path(bam)
    output:
        tuple val(sample_id), val(chr), path("${chr}.bed")
    script:
    """
        #Convertion of bam files to bed files
        bam2bed --all-reads --split --do-not-sort < ${bam} | gawk -v OFS="\\t" '{print \$1,\$2,\$3,\$6,\$4,\$14"::"\$15}' > ${chr}.bed
    """
}


process reads_mapping_and_filtering {
    tag "${sample_id}, ${chr}, ${bed}, ${exons}"
    publishDir "./results/reads_processing/mapping_filtering/${sample_id}"
    cache true
    label "big_mem"

    input:
        tuple val(sample_id), val(chr), path(bed), path(exons), path(exons_unique), path(collapsed)
    output:
        tuple val(sample_id), val(chr), path("${chr}.reads")
    script:
    """
        Rscript ${baseDir}/src/mapping_filtering.R ${bed} ${exons} ${chr}.reads
    """
}


process ip_splitting {
    tag "${chr}"
    publishDir "./results/internalp_filtering/ipdb_splitted"
    cache true
    label "small_mem"

    input:
        tuple val(chr), path(ipref)
    output:
        tuple val(chr), path("${chr}_sp.ipdb")
    script:
    """
    #select chromosome
    gawk '{ if(\$1=="${chr}") print }' ${ipref} > ${chr}_sp.ipdb
    """
}


process ip_filtering {
    tag "${sample_id}, ${chr}, ${ipdb}"
    publishDir "./results/internalp_filtering/internalp_filtered/${sample_id}"
    cache true
    label "big_mem"

    input:
        tuple val(sample_id), val(chr), path(reads), path(ipdb)
        val(ip_thr)
    output:
        tuple val(sample_id), val(chr), path("${chr}_mapped.ipdb"), path("${chr}_unique.reads"), path("${chr}.readid"), path("${chr}_ipf.reads")
    script:
        """
        #filtering
        Rscript ${baseDir}/src/ip_filtering.R ${reads} ${ipdb} ${ip_thr} ${chr}.readid ${chr}_unique.reads ${chr}_mapped.ipdb ${chr}_ipf.reads
        """
}
