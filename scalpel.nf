#!/usr/bin/env groovy

//Author:   Franz AKE
//Date:     10/2020
//PLASS_LAB - CMR_c Barcelona


//*************************************************************************************
//NEXTFLOW PIPELINE FOR SINGLE-CELL ALTERNATIVE POLYADENYLATION CHARACTERIZATION SC-APA (SCALPEL)
//*************************************************************************************

// Input arguments...
// ******************

/*
	Exemple of launching command:
	#****************************
	- DROP_SEQ:
	#*********
	nextflow run -resume ~/15K/APA_TOOL/main.nf --sequencing 'DROP_SEQ' --folder_in {} --annot {} --ipdb {}

	- 10X:
	#*****
	nextflow run -resume ~/15K/APA_TOOL/main.nf --sequencing '10X' --folder_in {} --annot {} --ipdb {}
*/



// **************************************
/* Check presence of required arguments */
// **************************************


/* DEFAULT PARAMETERS */
//*********************

/*publish repository path*/
params.publish_rep = "$PWD/scalpel_results/"

/*progs*/
params.samtools_bin = 'samtools'
params.rbin = 'Rscript'
params.python_bin = 'python3'
params.bedmap_bin = 'bedmap'
params.chr_concordance = ''

/*Some params initilialization*/
params.dt_threshold = 600
params.dt_exon_end_threshold = 20
params.cpu_defined = 50
params.mt_remove = 1
params.fraction_read_overlapping = 2
params.subsampling = 1
params.gene_fraction = "90%"
params.binsize = 15


if ( params.help )
	error """"
	===============================
	SCALPEL - N F   P I P E L I N E
	===============================

	Execution:
	- In case of providing 10X cell ranger folder:
	usage: nextflow run -resume scalpel.nf --sequencing <chromium> --folder_in <10X_folder> --annot <genome_annotation_reference> --ipdb <internal_priming_ref_file> --quant_file <salmon_preprocessed_file>

	- If providing Dropseq files or Others:
	usage: nextflow run -resume scalpel.nf --sequencing <dropseq> --bam <BAM> --bai <BAI> --dge_matrix <DGE> --barcodes <barcodes> --annot <genome_annotation_reference> --ipdb <internal_priming_ref_file> --quant_file <salmon_preprocessed_file>

	Output options:
	--folder_in,						Path to 10X Cellranger results folder [required if 10X file analysis]
	--bam,							Path to indexed BAM file [required]
	--bai,							Path to BAM index file	[required]
	--dge_matrix,						Path to DGE count matrix file [required]
	--quant_file,						Path to salmon quantification file from preprocessing [required]
	--ipdb, 						Path to internal priming reference annotation file [required]
	--barcodes,						Path to file containing valid barcodes [required]
	--annot,						Path to genomic annotation reference file [required]
	--sequencing,						Sequencing type [chromium,dropseq]

	[--dt_threshold] (optional),				Transcriptomic distance threshold
	[--dt_exon_end_threshold] (optional)			Transcriptomic end distance threhsold
	[--cpu_defined] (optional)				Max cpus (default, 50)
	[--subsampling]						BAM file subsampling threshold (default 1, select all reads)
	[--gene_fraction]					theshold fraction gene
	[--binsize]						binsize fragment probability
	[--publish_rep] (optional)				Publishing repository
	[--chr_concordance]					Character at add in order to match chromosome name in BAM file and the genome reference annotation file
	"""


/* Some execution prechecks*/
if ( params.sequencing == null )
	error "Enter sequencing argument (chromium or dropseq) !  [--sequencing]"

if ( params.annot == null)
	error "Enter Genomic annotation GTF file ! [--annot]"

if ( params.ipdb == null )
	error "Enter Internal Priming reference file ! [--ipdb]"

if ( params.quant_file == null )
	error "Provide salmon quantification preprocessing file ! [--quant_file]"

if ( params.sequencing == 'dropseq' ) {
	if (params.folder_in == null) {
		println """ No dropseq FOLDER LOCATION SPECIFIED...."""
		if( params.bam == null || params.bai == null || params.dge_matrix == null)
			error "No BAM, the indexed BAI file, or DGE_matrix file inputed"
	}
} else {
	if (params.folder_in == null) {
		println """ No chromium FOLDER LOCATION SPECIFIED...."""
		if( params.bam == null || params.bai == null || params.dge_matrix == null)
			error "No BAM, the indexed BAI file, or DGE_matrix file inputed"
	} else {
		println """chromium FOLDER LOCATION SPECIFIED...."""
		//define required paths
		params.bam = "${params.folder_in}/outs/possorted_genome_bam.bam"
		params.bai = "${params.folder_in}/outs/possorted_genome_bam.bam.bai"
		params.dge_matrix = "${params.folder_in}/outs/filtered_feature_bc_matrix.h5"
		params.barcodes = "${params.folder_in}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
	}
}

if ( params.barcodes == null )
	error "Provide valid barcodes file ! [--barcodes]"



// ************************************************************************
// Nice Start Message priting for the user (Resuming all the inputed files)
// ************************************************************************

println """\

		 ===============================
		 SCALPEL - N F   P I P E L I N E
		 ===============================
		 Author: PLASS Lab
		 *****************
		 Last update : October 2022
		 P-CMRC - Barcelona, SPAIN

		 folder_in: ${params.folder_in}
		 BAM file (required): ${params.bam}
		 BAI file (required): ${params.bai}
		 DGE file (required): ${params.dge_matrix}
		 Quantification Salmon file (required): ${params.quant_file}
		 Internal priming reference file (required): ${params.ipdb}
		 Barcodes file [--barcodes]: ${params.barcodes}

		 Annotation reference file (required):  ${params.annot}
		 Sequencing type (required): ${params.sequencing}
		 Transcriptomic distance threshold [--dt_threshold] (optional): ${params.dt_threshold}
		 Transcriptomic end distance threhsold [--dt_exon_end_threshold] (optional): ${params.dt_exon_end_threshold}
		 Max cpus [--cpu_defined] (default 50): ${params.cpu_defined}
		 subsampling [--subsampling]: ${params.subsampling}
		 mapQ threshold [--mapq]: ${params.mapq}
		 theshold fraction gene [--gene_fraction]: ${params.gene_fraction}
		 binsize fragment probability [--binsize]: ${params.binsize}
		 chromosome character paste [--chr_concordance]: ${params.chr_concordance}

		 Publishing repository [--publish_rep] (optional): ${params.publish_rep}

		 """
.stripIndent()


/*****************/
/* PREPROCESSING */
/*****************/

/* Computing genomic Annotation file GTF & Internal Priming annotation(1) */
// ********************************************************************

/*S1*/
process Gtf_annotation_splitting{
	tag "${gtf_file}"
	publishDir "${params.publish_rep}/gtfs/", overwrite: true
	maxForks params.cpu_defined
	input:
	val pythonbin				from params.python_bin
	file gtf_file				from file(params.annot)
	file quantf					from file(params.quant_file)
	output:
	file "*"				into gtf_splitted_ch1, gtf_splitted_ch2, gtf_splitted_ch3 mode flatten
	script:
	"""
	${pythonbin} ${baseDir}/src/gtf_splitting.py ${gtf_file} ${quantf}
	rm *MT* || echo 'ok'
	rm *M* || echo 'ok'
	rm *G* || echo 'ok'
	rm *J* || echo 'ok'
	"""
}



/*S2*/
process Chromosome_processing{
	tag "${chr_file.baseName}"
	publishDir "${params.publish_rep}/exons/", overwrite: true, mode: 'copy'
	maxForks params.cpu_defined
	input:
	val pythonbin								from params.python_bin
	file chr_file								from gtf_splitted_ch1
	file quantf									from file(params.quant_file)
	val trs_distance 							from params.dt_threshold
	val trs_end_distance 						from params.dt_exon_end_threshold
	output:
	file "${chr_file.baseName}.exons"			into exons_ch1, exons_ch2, exons_ch3, exons_ch4
	file "${chr_file.baseName}.exons_unique"	into exons_unique_ch
	file "${chr_file.baseName}.exon_bedmap"		into exons_bedmap
	script:
	"""
	${pythonbin} ${baseDir}/src/exon_processing.py ${chr_file}  ${trs_distance} ${trs_end_distance}  ${chr_file.baseName}.exons ${chr_file.baseName}.exons_unique ${chr_file.baseName}.exon_bedmap
	"""
}



/*S3*/
process Bam_splitting{
	tag "${chr_id.baseName}"
	publishDir "${params.publish_rep}/reads/bams/", overwrite: true, mode: 'copy'
	maxForks params.cpu_defined
	input:
	val chrc				from params.chr_concordance
	val subsamp				from params.subsampling
	val seqtype				from params.sequencing
	val samtoolbin			from params.samtools_bin
	val pythonbin 			from params.python_bin
	file chr_id				from gtf_splitted_ch2
	file bam_file 			from file(params.bam)
	file bai_file 			from file(params.bai)
	file barcodes_file 		from file(params.barcodes)
	output:
	file "${chr_id}.ebam"	into exonic_bams_ch, exonic_bams_ch2

	script:
	if( params.sequencing == "dropseq" )
		"""
		#Dropseq
		#Filter reads and split by chromosome
		zcat -f ${barcodes_file} > bc.txt
		${samtoolbin} view -b --subsample ${subsamp} ${bam_file} ${chrc}${chr_id.baseName} -D XC:bc.txt --keep-tag "XC,XM" > ${chr_id}.ebam
		"""
	else
		"""
		#10X SEQ
		zcat -f ${barcodes_file} > bc.txt
		${samtoolbin} view -b --subsample ${subsamp} ${bam_file} ${chrc}${chr_id.baseName} -D CB:bc.txt --keep-tag "CB,UB" > ${chr_id}.ebam
		"""
}



/*S4*/
process Bed_formatting{
	tag "${bamf.baseName}"
	publishDir "${params.publish_rep}/reads/beds/", overwrite: true
	maxForks params.cpu_defined
	input:
	file bamf					from exonic_bams_ch
	val pythonbin				from params.python_bin
	output:
	file "${bamf.baseName}.bed" 		into bed_ch
	script:
	"""
	#Convertion of bam files to bed files (1)
	bam2bed --all-reads --split --do-not-sort < ${bamf} > temp_file

	#Extract relevant metadata infos and sort for overlapping purposes (2)
	${pythonbin} ${baseDir}/src/subset_metatada.py temp_file ${bamf.baseName}.bed
	"""
}


/*S5*/
process Bed_exons_overlapping{
	tag "${bed.baseName} /${fraction_rd}"
	publishDir "${params.publish_rep}/reads/overlap/", overwrite: true
	maxForks params.cpu_defined
	input:
	val bedmapbin 					from params.bedmap_bin
	file bed						from bed_ch
	file exon_file					from exons_bedmap.collect()
	val fraction_rd 				from params.fraction_read_overlapping
	output:
	file "${bed.baseName}.bed_gtf"	into bed_gtf_ch
	script:
	"""
	#mapped reads on exons bed (1)
	${bedmapbin} --echo --echo-map-id --bp-ovr ${fraction_rd} --delim '\t' ${bed} ${bed.baseName}.exon_bedmap > ${bed.baseName}.bed_gtf
	"""
}


/*S6*/
process Bed_exons_filtering{
	tag "${ebed.baseName}"
	publishDir "${params.publish_rep}/reads/overlap/", overwrite: true
	maxForks params.cpu_defined
	input:
	val pythonbin								from params.python_bin
	file ebed 									from bed_gtf_ch
	file exon_file 								from exons_ch3.collect()
	val trs_end_distance 						from params.dt_exon_end_threshold
	output:
	file "${ebed.baseName}.bed_gtf_filtered"	into bed_gtf_filtered
	script:
	"""
	${pythonbin} ${baseDir}/src/filtering.py ${ebed} ${ebed.baseName}.exons ${trs_end_distance} ${ebed.baseName}.bed_gtf_filtered
	"""
}


/*S7*/
process Split_ipdb_database{
	tag "${ipdbf}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/ipdb_splitted", overwrite: true
	input:
	file ipdbf	from file(params.ipdb)
	output:
	file "*"	into ipdb_splitted mode flatten
	script:
	"""
	#Optionally remove the chr character in the file (1)
	#sed 's/chr//' ${ipdbf} > nochr_ipdb

	#select columns
	awk -v OFS="\\t" '{print \$1,\$2,\$3,\$6}' ${ipdbf}  > nochr_ipdb2

	#split by chromosome
	gawk '{print > \$1".ipdb"}' nochr_ipdb2
	rm nochr_ipdb2
	#rm nochr_ipdb
	"""
}


/*S8*/
process Overlapping_ip_exons{
	tag "${ebed.baseName}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/reads/ip_filtered/", overwrite: true
	input:
	val pythonbin								from params.python_bin
	val bedmapbin 								from params.bedmap_bin
	file all_ips 								from ipdb_splitted.collect()
	file ebed									from bed_gtf_filtered
	file uniqs 									from exons_unique_ch.collect()
	output:
	file "${ebed.baseName}.ip_filtered" 					into exip_mapped, exip_mapped2
	file "${ebed.baseName}.fragment_filtered_unique" 		into fragment_filtered_uniq_ch
	file "${ebed.baseName}.ipp"								into ip_files_ch
	file "${ebed.baseName}.rid"								into rid_files_ch
	script:
	"""
	#select exons columns in exon file
	awk -v OFS="\\t" '{print \$1,\$10,\$11,\$8,\$5,\$15,\$14}' ${ebed} | sort -u > selected_exons.txt

	#add index number to the ip file and arrange columns
	nl ${ebed.baseName}.ipdb | awk -v OFS="\t" '{print \$2,\$3,\$4,\$1}' > selected_ipdb.txt

	#merge exon_file
	${bedmapbin} --echo --echo-map-id --fraction-map 0.9 --delim '\t' selected_exons.txt selected_ipdb.txt > ${ebed.baseName}.exip_mapped

	#filter internal priming associated reads
	${pythonbin} ${baseDir}/src/merge_ip_to_frags.py ${ebed} ${ebed.baseName}.exip_mapped ${ebed.baseName}.ipdb ${ebed.baseName}.exons_unique ${ebed.baseName}.ip_filtered  ${ebed.baseName}.ipp ${ebed.baseName}.fragment_filtered_unique
	awk '{print \$4}' ${ebed.baseName}.ip_filtered | sort -u > ${ebed.baseName}.rid
	"""
}


/*S10*/
process Probability_distribution{
	publishDir "${params.publish_rep}/reads/probability/", overwrite: true, mode: 'copy'
	input:
	val r_bin				from params.rbin
	val bins				from params.binsize
	val gfrac				from params.gene_fraction
	file unique_files		from fragment_filtered_uniq_ch.collect()
	output:
	file "BINS_PROB.txt"	into probability_ch
	file "BINS_PROB.pdf"	into probability_pdf_file
	script:
	"""
	${r_bin} ${baseDir}/src/compute_prob_distribution.R "\$PWD/" ${gfrac} ${bins} BINS_PROB.txt BINS_PROB.pdf
	"""
}



/*S13*/
process Fragment_probabilities{
	tag "${ebed.baseName}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/reads/cells/chr_cells/", overwrite: true
	maxForks params.cpu_defined
	input:
	file ebed			from exip_mapped
	file prob_file 		from probability_ch
	output:
	file "*.frag_prob"								into frags_channel
	script:
	"""
	#Write cell Files
	python3 ${baseDir}/src/fragment_probabilities.py ${ebed} ${prob_file} ${ebed.baseName}.frag_prob
	"""
}


/*S14*/
process Cells_joining{
	tag "${frag.baseName}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/reads/cells/cells/", overwrite: true
	maxForks params.cpu_defined
	input:
	file frag				from frags_channel.collect()
	output:
	file "*.cell"			into cells_channel mode flatten
	script:
	"""
	# merge all the frag files
	# ************************
	cat *.frag_prob > all_fragments.txt

	# split the files
	# ***************
	awk -v OFS="\\t" '{print > \$1".cell"}' all_fragments.txt
	"""
}


/*S15*/
process EM_algorithm{
	tag "${cell_file.baseName}"
	publishDir "${params.publish_rep}/reads/prediction/", overwrite: true
	maxForks params.cpu_defined
	input:
	file cell_file	from cells_channel
	output:
	file "${cell_file.baseName}.pred"	into prediction_ch
	script:
	"""
	python3 ${baseDir}/src/em_algorithm.py ${cell_file} ${cell_file.baseName}.pred
	"""
}


/*S16*/
process Merge_predicted_cells{
	tag "${all_preds}"
	maxForks params.cpu_defined
	input:
	file all_preds					from prediction_ch.collect()
	output:
	file "predicted_cells.txt"		into merged_pred_cells_ch
	script:
	"""
	Rscript ${baseDir}/src/merge_pred_cells.R \$PWD predicted_cells.txt
	"""
}


/*S17*/
process DGE_generation{
	tag ""
	publishDir "${params.publish_rep}/reads/apa_dge/", overwrite: true, mode: 'copy'
	input:
	file all_cells						from merged_pred_cells_ch
	file dge_mat						from file(params.dge_matrix)
	output:
	file "APADGE.txt"                   into apa_dge_file
	script:
	if (params.sequencing == "dropseq")
		"""
		#Merge ALL cells prediction with DGE table
		python3 ${baseDir}/src/merge_DGE_PRED.py ${all_cells} ${dge_mat} APADGE.txt
		"""
	else
		"""
		#preprocess h5 file
		Rscript ${baseDir}/src/processh5.R ${dge_mat} DGE.txt
		#Merge ALL cells prediction with DGE table
		python3 ${baseDir}/src/merge_DGE_PRED.py ${all_cells} DGE.txt APADGE.txt
		"""
}


/*S18*/
process Filter_BAMS{
	tag "${read_file.baseName}"
	maxForks params.cpu_defined
	input:
	file all_bams			from exonic_bams_ch2.collect()
	file read_file 			from exip_mapped2
	val samtoolbin			from params.samtools_bin
	output:
	file "${read_file.baseName}.ebamf"		into filtered_bams_ch
	script:
	"""
	#get reads_id from ip filtered reads
	python3 ${baseDir}/src/get_readid.py ${read_file} ${read_file.baseName}.rid
	#filter bam file and indexing
	${samtoolbin} view -b -N ${read_file.baseName}.rid ${read_file.baseName}.ebam > ${read_file.baseName}.ebamf
	${samtoolbin} index ${read_file.baseName}.ebamf
	"""
}


/*S19*/
process Merge_BAMS{
	publishDir "${params.publish_rep}/reads/filtered_bam/", overwrite: true, mode: 'copy'
	input:
	val samtoolbin			from params.samtools_bin
	file all_bams			from filtered_bams_ch.collect()
	output:
	file "final.bam"
	file "final.bam.bai"
	script:
	"""
	#merge all the bams files
	${samtoolbin} merge -f -o final1.bam ${all_bams}
	${samtoolbin} sort final1.bam > final.bam
	${samtoolbin} index final.bam
	"""
}







