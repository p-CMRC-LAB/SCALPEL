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
params.publish_rep = "$PWD/RESULTS/"

/*progs*/
params.samtools_bin = 'samtools'
params.rbin = 'Rscript'
params.python_bin = 'python3'
params.bedmap_bin = 'bedmap'
params.mapq = 10

/*Some params initilialization*/
params.dt_threshold = 500
params.dt_exon_end_threshold = 30
params.cpu_defined = 50
params.mt_remove = 1
params.fraction_read_overlapping = 2
params.barcodes = 'null'
params.quant_file = 'null'

/* Input params provided by the user*/
// BAM file
// BAI file
// DGE matrix file
// Reference genome file
// Sequencing format
// Sample_name
// (Optional) Barcodes file
// (Optional) Alevin quantification file


/* Some execution prechecks*/
if (params.sequencing == null)
	error "Enter sequencing argument (chromium or dropseq) !  [--sequencing]"
if (params.sample_name == null)
	error "Enter sample name (Name of the BAM/BAI file) [--sample_name]"

if (params.sequencing == 'dropseq') {
	if (params.folder_in == null) {
		println """ No dropseq FOLDER LOCATION SPECIFIED...."""
		if( params.bam == null || params.bai == null || params.dge_matrix == null)
			error "No BAM, the indexed BAI file, or DGE_matrix file inputed"
	} else {
		println """dropseq FOLDER LOCATION SPECIFIED...."""
		//define required paths
		params.bam = "${params.folder_in}/${params.sample_name}.bam"
		params.bai = "${params.folder_in}/${params.sample_name}.bam.bai"
		params.dge_matrix = "${params.folder_in}/DGE.txt"
	}
} else {
	if (params.folder_in == null) {
		// TO DEFINE LATER
	} else {
		println """chromium FOLDER LOCATION SPECIFIED...."""
		//define required paths
		params.bam = "${params.folder_in}/${params.sample_name}.bam"
		params.bai = "${params.folder_in}/${params.sample_name}.bam.bai"
		params.dge_matrix = "${params.folder_in}/filtered_feature_bc_matrix.h5"
		params.barcodes = "${params.folder_in}/filtered_feature_bc_matrix/barcodes.tsv.gz"
	}
}
if ( params.annot == null)
	error "Enter Genomic annotation GTF file ! [--annot]"
if ( params.ipdb == null )
	error "Enter Internal Priming reference file ! [--ipdb]"
if (params.quant_file == null)
	println """ No salmon alevin quantification specified """







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

		 sample_name: ${params.sample_name}
		 folder_in: ${params.folder_in}
		 BAM file (required): ${params.bam}
		 BAI file (required): ${params.bai}
		 DGE file (required): ${params.dge_matrix}
		 Barcodes file (optional): ${params.barcodes}
		 Quantification alevin file (optional): ${params.quant_file}

		 Annotation reference file (required):  ${params.annot}
		 Sequencing type (required): ${params.sequencing}
		 Transcriptomic distance threshold (optional): ${params.dt_threshold}
		 Transcriptomic end distance threhsold (optional): ${params.dt_exon_end_threshold}
		 Max cpus (optional): ${params.cpu_defined}

		 Internal priming reference file: ${params.ipdb}

		 Publishing repository: ${params.publish_rep}

		 """
.stripIndent()


/*****************/
/* PREPROCESSING */
/*****************/

/* Computing genomic Annotation file GTF & Internal Priming annotation(1) */
// ********************************************************************

/*S1*/
process gtf_annotation_splitting{
	tag "${gtf_file}"
	publishDir "${params.publish_rep}/gtfs/", overwrite: true
	maxForks params.cpu_defined
	input:
	val pythonbin				from params.python_bin
	file gtf_file				from file(params.annot)
	file quantf					from file(params.quant_file)
	output:
	file "*"				into gtf_splitted_ch1, gtf_splitted_ch2 mode flatten
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
process chromosome_processing{
	tag "${chr_file.baseName}"
	publishDir "${params.publish_rep}/exons/", overwrite: true
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
	script:
	"""
	${pythonbin} ${baseDir}/src/exon_processing.py ${chr_file}  ${trs_distance} ${trs_end_distance}  ${chr_file.baseName}.exons ${chr_file.baseName}.exons_unique
	"""
}


/*S3*/
process internal_priming_splitting{
	tag "${ipdb_DB}"
	publishDir "${params.publish_rep}/ipdb_splitted", overwrite: true
	maxForks params.cpu_defined
	input:
	file ipdb_DB	from file(params.ipdb)
	output:
	file "*"		into ipdb_files mode flatten
	script:
	"""
	#IPDB_splitting
	awk '{print>\$1".ip_extracted"}' ${ipdb_DB}
	#rename -d 'chr' *
	"""
}



/*S4*/
process bam_splitting{
	tag "${chr_id.baseName}"
	publishDir "${params.publish_rep}/reads/bams/", overwrite: true
	maxForks params.cpu_defined
	input:
	val seqtype				from params.sequencing
	val samtoolbin			from params.samtools_bin
	val pythonbin 			from params.python_bin
	val mq					from params.mapq
	file chr_id				from gtf_splitted_ch2
	file bam_file 			from file(params.bam)
	file bai_file 			from file(params.bai)
	file barcodes_file 		from file(params.barcodes)
	output:
	file "${chr_id}.bamf"			into exonic_bams_ch, exonic_bams_ch2

	script:
	if(params.barcodes == 'null')
		if (params.sequencing == "dropseq")
			"""
			#Filter reads based on mapping quality (20) and split by chromosome
			${samtoolbin} view -H ${bam_file} > ${chr_id}.header
			${samtoolbin} view ${bam_file} ${chr_id} -e '[gf]=~"(CODING|UTR)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.esam
			${samtoolbin} view ${bam_file} ${chr_id} -e '[gf]=~"(INTERGENIC|INTRONIC)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.isam
			#Filter intronic associated reads
			${pythonbin} ${baseDir}/src/bam_filtering.py ${chr_id}.esam ${chr_id}.isam ${seqtype} >> ${chr_id}.header
			#Reconvert into Bam file
			${samtoolbin} view -b ${chr_id}.header > ${chr_id}.bamf
			"""
		else
			"""
			#10X SEQ
			#Filter reads based on mapping quality (20) and split by chromosome
			${samtoolbin} view -h -b  ${bam_file} ${chr_id.baseName} --keep-tag "GX,GN,CB,UB" > ${chr_id}.bam
			"""
	else if(params.barcodes != "NONE")
		if(params.sequencing == "dropseq")
			"""
			#(special for this run)
			#Filter reads based on mapping quality (20) and split by chromosome
			${samtoolbin} view -H ${bam_file} > ${chr_id}.header
			${samtoolbin} view ${bam_file} ${chr_id} -D XC:${barcodes_file} -e '[gf]=~"(CODING|UTR)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.esam
			${samtoolbin} view ${bam_file} ${chr_id} -D XC:${barcodes_file} -e '[gf]=~"(INTERGENIC|INTRONIC)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.isam
			#Filter intronic associated reads
			${pythonbin} ${baseDir}/src/bam_filtering.py ${chr_id}.esam ${chr_id}.isam ${seqtype} tempf
			cat tempf >> ${chr_id}.header
			#Reconvert into Bam file
			${samtoolbin} view -b ${chr_id}.header > ${chr_id}.bamf
			"""
		else
			"""
			#10X SEQ
			#Filter reads based on mapping quality (20) and split by chromosome
			${samtoolbin} view -H ${bam_file} > ${chr_id}.header
			${samtoolbin} view -D CB:${barcodes_file} ${bam_file} ${chr_id.baseName} -e '[RE]=~"E"' --keep-tag "GX,GN,CB,UB"  > ${chr_id}.esam
			${samtoolbin} view -D CB:${barcodes_file} ${bam_file} ${chr_id.baseName} -e '[RE]=~"I"' --keep-tag "GX,GN,CB,UB" > ${chr_id}.isam
			#Filter intronic associated reads
			${pythonbin} ${baseDir}/src/bam_filtering.py ${chr_id}.esam ${chr_id}.isam >> ${chr_id}.header
			#Reconvert into Bam file
			${samtoolbin} view -b ${chr_id}.header > ${chr_id}.bamf
			"""
}


/*S5*/
process bed_formatting{
	tag "${bamf.baseName}"
	publishDir "${params.publish_rep}/reads/beds/", overwrite: true
	maxForks params.cpu_defined
	input:
	file bamf					from exonic_bams_ch
	val pythonbin				from params.python_bin
	output:
	file "${bamf.baseName}.bed" 	into bed_ch
	file "${bamf.baseName}.bedsb"	into bed_coords
	script:
	"""
	#Convertion of bam files to bed files (1) and preprocessing of spliced reads for late filtering ops (2)
	bam2bed --all-reads --split --do-not-sort < ${bamf} | ${pythonbin} ${baseDir}/src/splicing_filtering.py - ${bamf.baseName}.bed

	#Extract coords and sort for overlapping purposes (bedmap - 3)
	awk -v OFS='\t' '{print \$1,\$2,\$3,\$4}' ${bamf.baseName}.bed | sort -k1,1 -k2,2n -k3,3n > ${bamf.baseName}.bedsb
	"""
}


/*S6*/
process bed_exons_overlapping{
	tag "${bedcds.baseName} /${fraction_rd}"
	publishDir "${params.publish_rep}/reads/overlap/", overwrite: true
	maxForks params.cpu_defined
	input:
	val bedmapbin 					from params.bedmap_bin
	file bed						from bed_ch
	file exon_file					from exons_ch2.collect()
	file bedcds						from bed_coords
	val fraction_rd 				from params.fraction_read_overlapping
	output:
	file "${bedcds.baseName}.bed_gtf"	into bed_gtf_ch
	file "${bed}"						into bed_ch2
	script:
	"""
	#Extract exon coords and sort (1)
	awk -v OFS='\t' '{if (NR!=1) {print \$1,\$2,\$3,\$6,\$9,\$11}}' ${bedcds.baseName}.exons | sort -k1,1 -k2,2n -k3,3n > exonsb

	#mapped reads on exons bed (2)
	${bedmapbin} --echo --echo-map --bp-ovr ${fraction_rd} --delim "****" ${bedcds} exonsb > ${bedcds.baseName}.bed_gtf
	"""
}


/*S7*/
process bed_exons_filtering{
	tag "${ebed.baseName}"	publishDir "${params.publish_rep}/reads/overlap/", overwrite: true
	maxForks params.cpu_defined
	input:
	val pythonbin								from params.python_bin
	file ebed 									from bed_gtf_ch
	file bed 									from bed_ch2
	file exon_file 								from exons_ch3.collect()
	output:
	file "${ebed.baseName}.bed_gtf_filtered"	into bed_gtf_filtered
	script:
	"""
	${pythonbin} ${baseDir}/src/overlapping_processing.py ${bed} ${ebed.baseName}.exons ${ebed} ${ebed.baseName}.bed_gtf_filtered
	"""
}


/*S8*/
process bed_splicing_filtering{
	tag "${ebed.baseName}"
	publishDir "${params.publish_rep}/reads/fragments/", overwrite: true
	maxForks params.cpu_defined
	input:
	val pythonbin									from params.python_bin
	file ebed										from bed_gtf_filtered
	val trs_distance								from params.dt_threshold
	output:
	file "${ebed.baseName}.bed_spliced_filtered"	into spliced_filtered_ch
	script:
	"""
	${pythonbin} ${baseDir}/src/splicing_filtering2.py ${ebed} ${trs_distance} ${ebed.baseName}.bed_spliced_filtered
	"""
}


/*S9*/
process Internal_priming_filtering{
	tag "${frag_file.baseName}"
	publishDir "${params.publish_rep}/reads/fragments/", overwrite: true
	maxForks params.cpu_defined
	input:
	val bedmapbin 											from params.bedmap_bin
	val pythonbin											from params.python_bin
	file frag_file											from spliced_filtered_ch
	file exon_file											from exons_ch4.collect()
	file ipn_file 											from ipdb_files.collect()
	file gtf_uniq_file										from exons_unique_ch.collect()
	output:
	file "${frag_file.baseName}.fragment_filtered"			into fragment_filtered_ch1
	file "${frag_file.baseName}.fragment_filtered_unique" 	into fragment_filtered_uniq_ch
	file "${frag_file.baseName}.readid" 					into readid_filtered_ch
	script:
	"""
	#Extract exon coords and sort (1)
	awk -v OFS='\t' '{if (NR!=1) {print \$1,\$2,\$3,\$6,\$9,\$11}}' ${frag_file.baseName}.exons | sort -k1,1 -k2,2n -k3,3n > exonsb
	sort -k1,1 -k2,2n -k3,3n ${frag_file.baseName}.ip_extracted > iptb
	#in case we have a presence of chr character
	#sed -i 's/chr//g' iptb
	${bedmapbin} --echo --echo-map --fraction-ref 1 --delim "****" iptb exonsb > ipsb
	${pythonbin} ${baseDir}/src/internalp_filtering.py ${frag_file} ${frag_file.baseName}.exons ipsb ${frag_file.baseName}.exons_unique ${frag_file.baseName}.fragment_filtered  ${frag_file.baseName}.fragment_filtered_unique ${frag_file.baseName}.readid
	"""
}


/*S10*/
process Probability_distribution{
	publishDir "${params.publish_rep}/reads/probability/", overwrite: true, mode: 'copy'
	input:
	val r_bin				from params.rbin
	file unique_files		from fragment_filtered_uniq_ch.collect()
	output:
	file "BINS_PROB.txt"	into probability_ch
	file "BINS_PROB.pdf"	into probability_pdf_file optional true
	script:
	"""
	${r_bin} ${baseDir}/src/compute_prob_distribution.R "\$PWD/" BINS_PROB.txt BINS_PROB.pdf
	"""
}


/*S13*/
process Fragment_probabilities{
	tag "${fragment_file.baseName}"
	publishDir "${params.publish_rep}/reads/cells/", overwrite: true
	input:
	file fragment_file	from fragment_filtered_ch1.collect()
	file prob_file 		from probability_ch
	output:
	file "*.cell"	into cells_ch mode flatten
	script:
	"""
	mkdir files
	cp *.fragment_filtered  files/.
	python3 ${baseDir}/src/fragment_probabilities.py "\$PWD/files/" ${prob_file}
	"""
}


/*S14*/
process EM_algorithm{
	tag "${cell_file.baseName}"
	publishDir "${params.publish_rep}/reads/prediction/", overwrite: true, mode:'copy'
	maxForks 35
	input:
	file cell_file	from cells_ch
	output:
	file "${cell_file.baseName}.pred"	into prediction_ch
	script:
	"""
	python3 ${baseDir}/src/em_algorithm.py ${cell_file} ${cell_file.baseName}.pred
	"""
}


/*S15*/
process DGE_generation{
	tag "${all_cells}"
	publishDir "${params.publish_rep}/reads/apa_dge/", overwrite: true, mode: 'copy'

	input:
	file all_cells						from prediction_ch.collect()
	file dge_mat						from file(params.dge_matrix)
	output:
	file "APADGE.txt"                   into apa_dge_file
	file "ALL_predicted_cells"			into pred_file

	script:
	"""
	#1- merge all the predicted cell files in an single One
	#******************************************************
	Rscript ${baseDir}/src/merge_pred_cells.R \$PWD ALL_predicted_cells

	#Merge ALL cells prediction with DGE table
	#*****************************************
	python3 ${baseDir}/src/merge_DGE_PRED.py ALL_predicted_cells ${dge_mat} APADGE.txt
	"""
}


/*S16*/
process outputSAM_filtering{
	tag "${bam}"
	input:
	val samtoolbin		from params.samtools_bin
	file fragment 		from readid_filtered_ch
	file bam 	  		from exonic_bams_ch2.collect()
	output:
	file "${fragment.baseName}.bam" into filtered_bam_ch

	script:
	"""
	# 1- Extract header from BAM file and save
	${samtoolbin} view -H ${fragment.baseName}.bamf > ${fragment.baseName}.samh

	# 2- Convert BAM file into SAM (no header included)
	${samtoolbin} view ${fragment.baseName}.bamf > ${fragment.baseName}.sam

	# 3- Filter reads present in the fragment file into the sam file
	python3 ${baseDir}/src/sam_filtering.py ${fragment.baseName}.sam ${fragment} ${fragment.baseName}.sam_reads

	# 3b - Merge header and reads
	cat ${fragment.baseName}.sam_reads >> ${fragment.baseName}.samh

	# 4- Transform into a bam file
	${samtoolbin} view -b ${fragment.baseName}.samh > ${fragment.baseName}.bam
	"""
}


/*S17*/
process INPUT_BAM_merging{
	tag "${filtered_bams}"
	publishDir "${params.publish_rep}/reads/filtered_bam", overwrite: true, mode: 'copy'

	input:
	val samtoolbin		from params.samtools_bin
	file filtered_bams from filtered_bam_ch.collect()
	output:
	file "final.bam"

	"""
	#Merge filtered BAM Files
	#************************
	${samtoolbin} merge -o FINAL.bam1 ${filtered_bams} -@ 20
	${samtoolbin} sort -@ 20 FINAL.bam1 -O BAM -o final.bam
	"""
}
