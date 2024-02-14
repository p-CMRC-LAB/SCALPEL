
/* 10X_processing:
    set of process for parsing 10X repository
*/

process get_10X_sampleIDs{
	tag "${sample_repo.baseName}"
    cpus 2
    cache true

    input:
        path sample_repo

    output:
    	path("*")

    script:
        """
        for temp in \$(ls -d ${sample_repo}/*/)
        do
        	sample_name=\$(basename \${temp::-1})
	        #Create a sampleID
	        touch \$sample_name
        done
        """
}


process read_10Xrepo{
	tag "${sample_id}"
	publishDir "./results/sample_files/", overwrite: true
	input:
		tuple path(sample_id), val(sample_repo)
	output:
		tuple val("${sample_id}"), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.barcodes.txt"), path("${sample_id}.counts.txt")
	script:
		"""
		#Copy into current directory the input files
		ln -s ${sample_repo}/${sample_id}/outs/possorted_genome_bam.bam ${sample_id.baseName}.bam
		ln -s ${sample_repo}/${sample_id}/outs/possorted_genome_bam.bam.bai ${sample_id.baseName}.bam.bai

		#Generate Barcodes and counts files
		Rscript ${baseDir}/src/read_10X.R ${sample_repo}/${sample_id}
		"""
}
