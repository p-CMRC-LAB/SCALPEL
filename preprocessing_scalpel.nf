#!/usr/bin/env groovy

//Author:   Franz AKE
//Date:     10/2020
//PLASS_LAB - CMR_c Barcelona


/************************************************************************************************/
//NEXTFLOW PIPELINE FOR SINGLE-CELL ALTERNATIVE POLYADENYLATION CHARACTERIZATION SC-APA (SCALPEL)
/*************************************************************************************************/

params.publish_rep = "$PWD/preprocessing/"

/*Required files*/
params.python_bin_path = 'python3'
//SALMON
params.salmon_index = null
params.salmon_path_bin = 'salmon'
params.salmon_quant_library_type = 'A'
params.salmon_quant_threads = 10
params.cpu_defined = 24

/*Get list of samples to preprocess*/
if (params.sample_names != null) {
	samples_ch = Channel.from(params.sample_names).splitCsv().flatten()
}else{
	error "provide sample names for the scalpel analysis ! [--sample_names]"
}

if (params.folder_in != null) {
	// Get fastq files into a channel
	fastq_files_ch = Channel.fromList(file("${params.folder_in}/*.fastq.gz"))
}else{
	error "provide fasta folder path ! [--folder_in]"
}

if (params.reference_fasta_transcript == null)
	error "provide [--reference_fasta_transcript]"



// ************************************************************************
// Nice Start Message priting for the user (Resuming all the inputed files)
// ************************************************************************

println """\

		 ===============================
		 SCALPEL - N F   P I P E L I N E
		 ===============================
		 Author: PLASS Lab
		 *****************
		 Last update : February 2022
		 P-CMRC - Barcelona, SPAIN

		 PREPROCESSING...

		 sample_name: ${params.sample_names}
		 folder_in: ${params.folder_in}
		 salmon_index: ${params.salmon_index}
		 Publishing repository: ${params.publish_rep}

		 """
.stripIndent()




/*SALMON inputs*/
/*Initiate salmon indexing if no salmon index provided by the user*/
if (params.salmon_index == null) {
	process salmon_indexing{
		tag "${reference_ft}"
		input:
		val salmon_path from params.salmon_path_bin
		path reference_ft from params.reference_fasta_transcript
		output:
		file "${reference_ft.baseName}.salmon_index" into salmon_index_ch
		script:
		"""
		${salmon_path} index -t ${reference_ft} -i ${reference_ft.baseName}.salmon_index
		"""
	}

	/*Quantification of transcripts for samples provided*/
	process salmon_quantification1{
		tag "${sample}, ${salmon_idx}"
		publishDir "${params.publish_rep}/salmon_quant/", overwrite: true
		input:
		val sample from samples_ch
		val salmon_path from params.salmon_path_bin
		file fastqs from fastq_files_ch.collect()
		file salmon_idx from salmon_index_ch
		val library_type from params.salmon_quant_library_type
		val threads from params.salmon_quant_threads
		output:
		file "quants/${sample}_quant" into quantification_folder
		file "${sample}.sf" into quant_files_ch
		script:
		"""
		${salmon_path} quant -i ${salmon_idx} -l ${library_type} -1 ${sample}*R1*gz -2 ${sample}*R2*gz --validateMappings -o quants/${sample}_quant
		cp quants/${sample}_quant/quant.sf ${sample}.sf1
		awk -v FS='\t' -v OFS='\t' '{print \$0,"${sample}"}' ${sample}.sf1 > ${sample}.sf
		"""
	}

}else{

	/*Quantification of transcripts for samples provided*/
	process salmon_quantification2{
		tag "${sample}, ${salmon_idx}"
		publishDir "${params.publish_rep}/salmon_quant/", overwrite: true
		input:
		val sample from samples_ch
		val salmon_path from params.salmon_path_bin
		file fastqs from fastq_files_ch.collect()
		path salmon_idx from params.salmon_index
		val library_type from params.salmon_quant_library_type
		val threads from params.salmon_quant_threads
		output:
		file "quants/${sample}_quant" into quantification_folder
		file "${sample}.sf" into quant_files_ch
		script:
		"""
		${salmon_path} quant -i ${salmon_idx} -l ${library_type} -1 ${sample}*R1*gz -2 ${sample}*R2*gz --validateMappings -o quants/${sample}_quant
		cp quants/${sample}_quant/quant.sf ${sample}.sf1
		awk -v FS='\t' -v OFS='\t' '{print \$0,"${sample}"}' ${sample}.sf1 > ${sample}.sf
		"""
	}
}

/*process salmon quantification files*/
process quantification_processing{
	tag "${quant} / ${sps}"
	publishDir "${params.publish_rep}/${sps}/", overwrite: true, mode: 'copy'
	input:
	file sps from params.sample_names
	file quant from quant_files_ch.collect()
	output:
	file "quant.filtered" into quant_filtered_ch
	script:
	"""
	python3 ${baseDir}/src/quantification_processing.py -salmon_quant "*.sf" -output_file 'quant.filtered'
	"""
}
