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
<<<<<<< HEAD
params.mapq = 0
=======
params.mapq = 10
>>>>>>> 3a0e40369ab788e3b2e9b4fa9e8a61b986a28aff

/*Some params initilialization*/
params.dt_threshold = 750
params.dt_exon_end_threshold = 30
params.cpu_defined = 50
params.mt_remove = 1
params.fraction_read_overlapping = 2
params.subsampling = 1
params.gene_fraction = "90%"
params.binsize = 25


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
<<<<<<< HEAD
	[--mapq]						have mapping quality >= INT (default, 0)
=======
	[--mapq]						have mapping quality >= INT
>>>>>>> 3a0e40369ab788e3b2e9b4fa9e8a61b986a28aff
	[--gene_fraction]					theshold fraction gene
	[--binsize]						binsize fragment probability
	[--publish_rep] (optional)				Publishing repository
	[--chr_concordance]					Charachter at add in order to match chromosome name in BAM file and the genome reference annotation file
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
process gtf_annotation_splitting{
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
	file "${chr_file.baseName}.exon_bedmap"		into exons_bedmap
	file "${chr_file.baseName}.exon_bedmap2"		into exons_bedmap2
	script:
	"""
	${pythonbin} ${baseDir}/src/exon_processing.py ${chr_file}  ${trs_distance} ${trs_end_distance}  ${chr_file.baseName}.exons ${chr_file.baseName}.exons_unique ${chr_file.baseName}.exon_bedmap ${chr_file.baseName}.exon_bedmap2
	"""
}


/*S3*/
process internal_priming_splitting{
	tag "${ipdb_DB}"
	publishDir "${params.publish_rep}/ipdb_splitted", overwrite: true
	maxForks params.cpu_defined
	input:
	val pythonbin	from params.python_bin
	file ipdb_DB	from file(params.ipdb)
	file chr_id		from gtf_splitted_ch3
	output:
	file "*"		into ipdb_files mode flatten
	script:
	"""
	${pythonbin} ${baseDir}/src/parse_and_split_ipdb.py ${ipdb_DB} ${chr_id} ${chr_id}.ip_extracted
	"""
}


/*S4*/
process bam_splitting{
	tag "${chr_id.baseName}"
	publishDir "${params.publish_rep}/reads/bams/", overwrite: true
	maxForks params.cpu_defined
	input:
	val chrc				from params.chr_concordance
	val subsamp				from params.subsampling
	val seqtype				from params.sequencing
	val samtoolbin			from params.samtools_bin
	val pythonbin 			from params.python_bin
	val mq					from params.mapq
	file chr_id				from gtf_splitted_ch2
	file bam_file 			from file(params.bam)
	file bai_file 			from file(params.bai)
	file barcodes_file 		from file(params.barcodes)
	output:
	file "${chr_id}.bamf"			into exonic_bams_ch, exonic_bams_ch2 optional true

	script:
	if(params.barcodes == 'null')
		if (params.sequencing == "dropseq")
			"""
			#Filter reads based on mapping quality (20) and split by chromosome
			${samtoolbin} view -H ${bam_file} > ${chr_id}.header
			${samtoolbin} view --subsample ${subsamp} ${bam_file} ${chrc}${chr_id} -e '[gf]=~"(CODING|UTR)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.esam
			${samtoolbin} view --subsample ${subsamp} ${bam_file} ${chrc}${chr_id} -e '[gf]=~"(INTERGENIC|INTRONIC)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.isam
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
			${samtoolbin} view -q ${mq} --subsample ${subsamp} ${bam_file} ${chrc}${chr_id.baseName} -e '[RE]=~"E"' --keep-tag "CB,UB"  > ${chr_id}.esam
			${samtoolbin} view -q ${mq} --subsample ${subsamp} ${bam_file} ${chrc}${chr_id.baseName} -e '[RE]=~"I"' --keep-tag "CB,UB"  > ${chr_id}.isam
			#Filter intronic associated reads
			${pythonbin} ${baseDir}/src/bam_filtering.py ${chr_id}.esam ${chr_id}.isam ${seqtype} tempf
			cat tempf >> ${chr_id}.header
			#Reconvert into Bam file
			${samtoolbin} view -b ${chr_id}.header > ${chr_id}.bamf
			"""
	else if(params.barcodes != "NONE")
		if(params.sequencing == "dropseq")
			"""
			#(special for this run)
			#Filter reads based on mapping quality (20) and split by chromosome
			${samtoolbin} view -H ${bam_file} > ${chr_id}.header
			${samtoolbin} view -q ${mq} --subsample ${subsamp} ${bam_file} ${chrc}${chr_id} -D XC:${barcodes_file} -e '[gf]=~"(CODING|UTR)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.esam
			${samtoolbin} view -q ${mq} --subsample ${subsamp} ${bam_file} ${chrc}${chr_id} -D XC:${barcodes_file} -e '[gf]=~"(INTERGENIC|INTRONIC)"' --keep-tag "XC,XM,gf,gs"  > ${chr_id}.isam
			#Filter intronic associated reads
			${pythonbin} ${baseDir}/src/bam_filtering.py ${chr_id}.esam ${chr_id}.isam ${seqtype} tempf
			cat tempf >> ${chr_id}.header
			#Reconvert into Bam file
			${samtoolbin} view -b ${chr_id}.header > ${chr_id}.bamf
			"""
		else
			"""
			#10X SEQ
			zcat ${barcodes_file} > bc.txt
			#Filter reads based on mapping quality (20) and split by chromosome
			${samtoolbin} view -H ${bam_file} > ${chr_id}.header
			${samtoolbin} view -q ${mq} --subsample ${subsamp} -D CB:bc.txt ${bam_file} ${chrc}${chr_id.baseName}  --keep-tag "CB,UB"  > ${chr_id}.esam
			${samtoolbin} view -q ${mq} --subsample ${subsamp} -D CB:bc.txt ${bam_file} ${chrc}${chr_id.baseName} -e '[RE]=~"I"' --keep-tag "CB,UB"  > ${chr_id}.isam
			if [ -s ${chr_id}.esam ];then
				if [ -s ${chr_id}.isam ];then
					#Filter intronic associated reads
					${pythonbin} ${baseDir}/src/bam_filtering.py ${chr_id}.esam ${chr_id}.isam ${seqtype} tempf
					cat tempf >> ${chr_id}.header
					#Reconvert into Bam file
					${samtoolbin} view -b ${chr_id}.header > ${chr_id}.bamf
				else
					#if presence of no intronic reads...
					cat ${chr_id}.esam >> ${chr_id}.header
					#Reconvert into Bam file
					${samtoolbin} view -b ${chr_id}.header > ${chr_id}.bamf
				fi
			else
				echo "Empty file"
			fi
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
	file "${bamf.baseName}.bed" 		into bed_ch
	file "${bamf.baseName}.bed_bedmap"	into bed_coords
	script:
	"""
	#Convertion of bam files to bed files (1) and preprocessing of spliced reads for late filtering ops (2)
	bam2bed --all-reads --split --do-not-sort < ${bamf} | ${pythonbin} ${baseDir}/src/splicing_filtering.py - ${bamf.baseName}.bed ${bamf.baseName}.bed_bedmap

	#Extract coords and sort for overlapping purposes (bedmap - 3)
	#awk -v OFS='\t' '{print \$1,\$2,\$3,\$4}' ${bamf.baseName}.bed | sort -k1,1 -k2,2n -k3,3n > ${bamf.baseName}.bedsb
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
	file exon_file					from exons_bedmap.collect()
	file bedcds						from bed_coords
	val fraction_rd 				from params.fraction_read_overlapping
	output:
	file "${bedcds.baseName}.bed_gtf"	into bed_gtf_ch
	file "${bed}"						into bed_ch2
	script:
	"""
	#Extract exon coords and sort (1)
	#awk -v OFS='\t' '{if (NR!=1) {print \$1,\$2,\$3,\$6,\$9,\$11}}' ${bedcds.baseName}.exons | sort -k1,1 -k2,2n -k3,3n > exonsb

	#mapped reads on exons bed (2)
	${bedmapbin} --echo --echo-map --bp-ovr ${fraction_rd} --delim "****" ${bedcds} ${bedcds.baseName}.exon_bedmap > ${bedcds.baseName}.bed_gtf
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
	file "${ebed.baseName}.bed_gtf_filtered"	into bed_gtf_filtered optional true
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
	file exon_bd2											from exons_bedmap2.collect()
	file exon_file											from exons_ch4.collect()
	file ipn_file 											from ipdb_files.collect()
	file gtf_uniq_file										from exons_unique_ch.collect()
	output:
	file "${frag_file.baseName}.fragment_filtered"			into fragment_filtered_ch1
	file "${frag_file.baseName}.fragment_filtered_unique" 	into fragment_filtered_uniq_ch
	file "${frag_file.baseName}.readid" 					into readid_filtered_ch
	script:
	"""
	#extract exon coords and sort (1)
	#awk -v OFS='\t' '{if (NR!=1) {print \$1,\$2,\$3,\$6,\$9,\$11}}' ${frag_file.baseName}.exons | sort -k1,1 -k2,2n -k3,3n | sed 's/chr//g' > exonsb
	#sort -k1,1 -k2,2n -k3,3n ${frag_file.baseName}.ip_extracted > iptb
	#in case we have a presence of chr character
	#sed -i 's/chr//g' iptb

	${bedmapbin} --echo --echo-map --fraction-ref 1 --delim "****" ${frag_file.baseName}.ip_extracted ${frag_file.baseName}.exon_bedmap2 > ipsb
	${pythonbin} ${baseDir}/src/internalp_filtering.py ${frag_file} ${frag_file.baseName}.exons ipsb ${frag_file.baseName}.exons_unique ${frag_file.baseName}.fragment_filtered  ${frag_file.baseName}.fragment_filtered_unique ${frag_file.baseName}.readid
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
	file "BINS_PROB.pdf"	into probability_pdf_file optional true
	script:
	"""
	${r_bin} ${baseDir}/src/compute_prob_distribution.R "\$PWD/" ${gfrac} ${bins} BINS_PROB.txt BINS_PROB.pdf
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
	if (params.sequencing == "chromium")
		"""
		#1- merge all the predicted cell files in an single One
		#******************************************************
		Rscript ${baseDir}/src/merge_pred_cells.R \$PWD ALL_predicted_cells

		#Merge ALL cells prediction with DGE table
		#*****************************************
		python3 ${baseDir}/src/merge_DGE_PRED.py ALL_predicted_cells ${dge_mat} APADGE.txt
		"""
	else
		"""
		#preprocess h5 file
		Rscript ${baseDir}/src/processh5.R ${dge_mat} DGE.txt

		#1- merge all the predicted cell files in an single One
		#******************************************************
		Rscript ${baseDir}/src/merge_pred_cells.R \$PWD ALL_predicted_cells

		#Merge ALL cells prediction with DGE table
		#*****************************************
		python3 ${baseDir}/src/merge_DGE_PRED.py ALL_predicted_cells DGE.txt APADGE.txt
		"""
}


/*S16*/
process outputSAM_filtering{
	tag "${bam}"
	input:
	val seqtype			from params.sequencing
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
	python3 ${baseDir}/src/sam_filtering.py ${fragment.baseName}.sam ${fragment} ${seqtype} ${fragment.baseName}.sam_reads

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
