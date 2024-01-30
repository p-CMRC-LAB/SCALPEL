

/* Annotation_processing:
    set of process for preprocessing of Scalpel annotation files
*/

process salmon_indexing{
    publishDir "./results/annotation_processing/salmon_indexing", overwrite: true
    cpus params.cpus
    cache true

    input:
        path transcriptome_reference

    output:
        path 'transcriptome_index'

    script:
        """
        salmon index --threads ${params.threads} -t ${transcriptome_reference} -i transcriptome_index --gencode
        """
}

process salmon_bulk_quantification{
    tag "${pair_id}"
    publishDir "./results/annotation_processing/salmon_bulk_quantification", overwrite: true
    cpus params.cpus
    cache true

    input:
        path transcriptome_index
        path gtf_annotation_reference
        tuple val(pair_id), path(reads)
        cpus params.threads

    output:
        path "${pair_id}.sf"

    script:
        """
        salmon quant -i ${transcriptome_index} -l A -g ${gtf_annotation_reference} -1 ${reads[0]} -2 ${reads[1]} -o ${pair_id} --threads ${params.threads}
        mv ${pair_id}/quant.sf tmp.sf
        awk -v OFS="\t" 'NR > 1 {print \$0}' tmp.sf | sed "s/\$/\t${pair_id}/" > ${pair_id}.sf
        """
}

process tpm_bulk_average{
    tag "${bulk_quants}"
	publishDir "./results/annotation_processing/salmon_bulk_quantification", overwrite: true, mode: 'copy'
    cpus params.cpus
    cache true

	input:
    	path bulk_quants
        path gtf_annotation_reference

	output:
	    file "merge_quants.sf"

	script:
        """
        cat *.sf > all.sf
        #QUANTS=\$(echo ${bulk_quants} | sed 's/ /,/g')
        Rscript ${baseDir}/src/quantification_processing.R all.sf ${gtf_annotation_reference} merge_quants.sf
        """
}

process gtf_splitting_bedfile_conversion{
    publishDir "./results/annotation_processing/isoform_processing", overwrite: true
    cache true
    cpus params.cpus

	input:
        path gtf_annotation_reference

	output:
    	path "*.gtf"
	
    script:
        """
        Rscript ${baseDir}/src/gtf_split.R ${gtf_annotation_reference}
        """
}

process isoform_selection_weighting{
    tag "${gtf.baseName}"
    publishDir "./results/annotation_processing/isoform_processing", overwrite: true, mode: 'copy'
    cpus params.cpus
    cache true

	input:
        file merged_quants
        path gtf

	output:
        tuple val("${gtf.baseName}"), path("${gtf.baseName}.exons"), path("${gtf.baseName}.exons_unique")

	script:
        """
        Rscript ${baseDir}/src/gtf_processing2.R ${gtf} ${merged_quants} ${params.dt_threshold} ${params.dt_exon_end_threshold}  ${gtf.baseName}.exons ${gtf.baseName}.exons_unique
        """
}





