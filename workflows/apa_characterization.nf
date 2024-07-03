

/* 
	Downstream analysis of samples (default params)
*/


process differential_isoform_usage{
	tag "${objs}"
	publishDir "./results/", overwrite: true, mode: 'copy'
        label 'big_mem'

	input:
		file(objs)

	output:
		tuple path("DIU_table.csv"), file("iDGE_seurat.RDS")
	
	script:
		if( params.clusters == null )
			"""
			cp ${baseDir}/src/scalpel_library.R .
			Rscript ${baseDir}/src/analysis_samples.R ./ NULL DIU_table.csv
			"""

		else if( params.clusters != null )
			"""
			cp ${baseDir}/src/scalpel_library.R .
			Rscript ${baseDir}/src/analysis_samples.R ./ ${params.clusters} DIU_table.csv
			"""
}
