

/* 
	Downstream analysis of samples (default params)
*/


process differential_isoform_usage{
	tag "${objs}"
	publishDir "./results/final_results/", overwrite: true, mode: 'copy'

	input:
		file(objs)

	output:
		tuple path("differential_isoforms_usage_table.csv"), file("final_seurat_obj.RDS")
	
	script:
		if( params.clusters == null )
			"""
			cp ${baseDir}/src/scalpel_library.R .
			Rscript ${baseDir}/src/analysis_samples.R ./ NULL differential_isoforms_usage_table.csv
			"""

		else if( params.clusters != null )
			"""
			cp ${baseDir}/src/scalpel_library.R .
			Rscript ${baseDir}/src/analysis_samples.R ./ ${params.clusters} differential_isoforms_usage_table.csv
			"""
}