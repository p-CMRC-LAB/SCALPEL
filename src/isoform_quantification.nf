



process probability_distribution {
	tag "${sample_id}"
	publishDir "./results/isoform_quantification/fragment_probabilities/${sample_id}", overwrite:'true', mode: 'copy'
	cache true

	input:
		tuple val(sample_id), file(reads)

	output:
		tuple val(sample_id), path("${sample_id}_probabilities.txt"), file("${sample_id}_probabilities.pdf")

	script:
        """
        cat *_unique.reads > all_unique.reads
        Rscript ${baseDir}/src/compute_prob.R all_unique.reads ${params.gene_fraction} ${params.binsize} ${sample_id}_probabilities.txt ${sample_id}_probabilities.pdf

        #delete temp files
        rm all_unique.reads
        """
}



process fragment_probabilities{
	tag "${sample_id},${chr}"
	publishDir "./results/isoform_quantification/fragment_probabilities/${sample_id}", overwrite: true, mode: 'copy'
	cache true

	input:
		tuple val(sample_id), val(chr), path(reads), path(probabilities)

	output:
    	tuple val(sample_id), path("${chr}_frag.reads")

	script:
        """
        #Write cell Files
		Rscript ${baseDir}/src/fragment_probabilities.R ${reads} ${probabilities} ${chr}_frag.reads
        """
}



process cells_splitting{
	tag "${sample_id}"
	cache true
	
	input:
		tuple val(sample_id), file(reads)

	output:
		tuple val(sample_id), path("*.cell")
	
	script:
		"""
		Rscript ${baseDir}/src/merge_and_writecells.R "."
		"""
}


process em_algorithm{
	tag "${sample_id}, ${cell.baseName}"
	cache true

	input:
		tuple val(sample_id), path(cell)
	
	output:
		tuple val(sample_id), path("${cell.baseName}.pred")
	
	script:
		"""
		Rscript ${baseDir}/src/em_algorithm.R ${cell} ${cell.baseName}.pred
		"""
}


process cells_merging{
	tag "${sample_id}"
	publishDir "./results/isoform_quantification/fragment_probabilities/${sample_id}", overwrite: true
	cache true

	input:
		tuple val(sample_id), file(pred_cells)

	output:
		tuple val(sample_id), path("${sample_id}_isoforms_quantified.txt")

	script:
		"""
		Rscript ${baseDir}/src/merge_pred_cells.R \$PWD ${sample_id}_isoforms_quantified.txt ${params.threads}
		"""
}


process dge_generation{
	tag "${sample_id}, ${isoforms}, ${raw_dge}"
	publishDir "./results/final_results/${sample_id}", overwrite: true, mode: 'copy'
	cache true

	input:
		tuple val(sample_id), path(isoforms), path(raw_dge)

	output:
		tuple val(sample_id), path("${sample_id}_APADGE.txt"), path("${sample_id}_seurat.RDS")
	
	script:
		"""
		#Merge ALL cells prediction with DGE table
		python3 ${baseDir}/src/merge_DGE_PRED.py ${isoforms} ${raw_dge} ${sample_id}_APADGE.txt

		#merge it into seurat object
		Rscript ${baseDir}/src/APAtoseurat.R ${sample_id}_APADGE.txt ${sample_id} ${sample_id}_seurat.RDS
		"""
}



