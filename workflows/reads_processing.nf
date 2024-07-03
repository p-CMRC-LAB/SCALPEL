#!/usr/bin/env nextflow

/* 
Loading of SAMPLES data & preprocessing according to the sequencing type
=====================================================================
*/

/* loading of subworkflows */
include { get_10X_sampleIDs; read_10Xrepo } from '../subworkflows/loading_10X_files.nf'
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
        if ( "${params.sequencing}" == "10x" ){
            /*get CR path*/
            get_10X_sampleIDs(samples_paths.flatMap{ it=it[3] })

            /* extract sampleIDs and associated paths */
            get_10X_sampleIDs.out.flatten().combine(samples_paths.flatMap{ it=it[3] }).set{ sample_ids }
            read_10Xrepo(sample_ids)

            /* formatting */
            read_10Xrepo.out.map{ it = [it[0], [it[1], it[2], it[3], it[4]]] }.set{ sample_ids_ch }
        } else {
            samples_paths.map{ it = tuple( it[0], file(it[3]), file(it[4]), file(it[5]), file(it[6]) ) }.set{ sample_ids_ch }
	}

        /* - processing of input BAM file... */
        bam_splitting( selected_isoforms.flatMap { it = it[0] }.combine(sample_ids_ch) )

    emit:
        selected_bams = bam_splitting.out
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
