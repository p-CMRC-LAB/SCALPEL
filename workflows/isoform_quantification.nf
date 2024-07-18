



process probability_distribution {
	tag "${sample_id}"
	publishDir "${params.output}/isoform_quantification/fragment_probabilities/${sample_id}", overwrite:'true'
	cache true
        label "big_mem"

	input:
	tuple val(sample_id), file(reads)
        val(gene_frac)
        val(binsz)

	output:
	tuple val(sample_id), path("${sample_id}_probabilities.txt"), file("${sample_id}_probabilities.pdf")

	script:
        """
        cat *_unique.reads > all_unique.reads
        Rscript ${baseDir}/src/compute_prob.R all_unique.reads ${gene_frac} ${binsz} ${sample_id}_probabilities.txt ${sample_id}_probabilities.pdf
        """
}



process fragment_probabilities{
	tag "${sample_id},${chr}"
	publishDir "${params.output}/isoform_quantification/fragment_probabilities/${sample_id}", overwrite: true, mode: 'copy'
	cache true
        label "big_mem"

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
        label "big_mem"
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
        label "small_mem"

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
	publishDir "${params.output}/isoform_quantification/fragment_probabilities/${sample_id}", overwrite: true
	cache true
        label "big_mem"

	input:
	tuple val(sample_id), file(pred_cells)

	output:
	tuple val(sample_id), path("${sample_id}_isoforms_quantified.txt")

	script:
	"""
        cat *.pred >> ${sample_id}_isoforms_quantified.txt
	"""
}


process dge_generation{
	tag "${sample_id}, ${isoforms}, ${raw_dge}"
	publishDir "${params.output}/final_results/${sample_id}", overwrite: true
	cache true
        label "big_mem"

	input:
	tuple val(sample_id), path(isoforms), path(raw_dge)

	output:
	tuple val(sample_id), path("${sample_id}_APADGE.txt"), path("${sample_id}_seurat.RDS")
	
	script:
	"""
	#Merge ALL cells prediction with DGE table
	Rscript ${baseDir}/src/merge_DGE_PRED.R ${raw_dge} ${isoforms} ${sample_id}_APADGE.txt

	#merge it into seurat object
	Rscript ${baseDir}/src/APAtoseurat.R ${sample_id}_APADGE.txt ${sample_id} ${sample_id}_seurat.RDS
	"""
}



