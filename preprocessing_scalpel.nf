#!/usr/bin/env groovy

//Author:   Franz AKE
//Date:     10/2020
//PLASS_LAB - CMR_c Barcelona


/************************************************************************************************/
//NEXTFLOW PIPELINE FOR SINGLE-CELL ALTERNATIVE POLYADENYLATION CHARACTERIZATION SC-APA (SCALPEL)
/*************************************************************************************************/

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
	error "provide sample names for the scalpel analysis"
}

if (params.folder_in != null) {
	// Get fastq files into a channel
	fastq_files_ch = Channel.fromList(file("${params.folder_in}/*.fastq.gz"))
}else{
	error "provide --folder_in path"
}

if (params.reference_fasta_transcript == null)
	error "provide --reference_fasta_transcript"



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
}else{
	salmon_index_ch = params.salmon_path
}


/*Quantification of transcripts for samples provided*/
process salmon_quantification{
	tag "${sample}"
	publishDir "$PWD/SALMON_QUANTS", overwrite: true, mode: 'copy'
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
	${salmon_path} quant -i ${salmon_idx} -l ${library_type} -1 ${sample}*read1.fastq* -2 ${sample}*2.fastq* --validateMappings -o quants/${sample}_quant
	cp quants/${sample}_quant/quant.sf ${sample}.sf1
	awk -v FS='\t' -v OFS='\t' '{print \$0,"${sample}"}' ${sample}.sf1 > ${sample}.sf
	"""
}

/*process salmon quantification files*/
process quantification_processing{
	tag "${quant}"
	publishDir "$PWD/SALMON_QUANTS", overwrite: true, mode: 'copy'
	input:
	file quant from quant_files_ch.collect()
	output:
	file "quant.filtered" into quant_filtered_ch
	script:
	"""
	python3 ${baseDir}/quantification_processing.py -salmon_quant "*.sf" -output_file 'quant.filtered'
	"""
}
