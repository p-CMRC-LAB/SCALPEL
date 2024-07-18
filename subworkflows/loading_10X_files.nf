
/* 10X_processing:
    set of process for parsing 10X repository
*/

process read_10Xrepo{
	tag "${sample_id}"
	publishDir "${params.output}/sample_files/", overwrite: true
	cache true
        label 'small_mem'

	input:
		tuple val(sample_id), val(sample_repo)
	output:
		tuple val("${sample_id}"), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.barcodes.txt"), path("${sample_id}.counts.txt")
	script:
	"""
	#Copy into current directory the input files
	ln -s ${sample_repo}/outs/possorted_genome_bam.bam ${sample_id}.bam
	ln -s ${sample_repo}/outs/possorted_genome_bam.bam.bai ${sample_id}.bam.bai

	#Generate Barcodes and counts files
	Rscript ${baseDir}/src/read_10X.R ${sample_repo}
	"""
}
