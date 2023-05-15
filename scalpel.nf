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

/*Some params initilialization*/
params.dt_threshold = 1000
params.dt_exon_end_threshold = 25
params.isoform_end_ip_threshold = 250

params.cpu_defined = 50
params.mt_remove = 1
params.fraction_read_overlapping = 0.1
params.subsampling = 1
params.gene_fraction = "80%"
params.binsize = 30


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
	--ipdb,							Path to internal priming reference annotation file [required]
	--barcodes,						Path to file containing valid barcodes [required]
	--annot,						Path to genomic annotation reference file [required]
	--sequencing,						Sequencing type [chromium,dropseq]

	[--dt_threshold] (optional),				Transcriptomic distance threshold (default, 1000)
	[--dt_exon_end_threshold] (optional)			Transcriptomic end distance between exons threhsold (default, 20)
	[--isoform_end_ip_threshold] (optional)			Minimal distance of the internal priming position from the isoform 3'end (default, 100)
	[--cpu_defined] (optional)				Max cpus (default, 50)
	[--subsampling] (optional)						BAM file subsampling threshold (default 1, select all reads)
	[--gene_fraction] (optional)					theshold fraction gene based on on expression abundance for probabilities estimation (default, 90%)
	[--binsize] (optional)						binsize on transcriptomic space for fragment probabilities estimation (default, 15)
	[--publish_rep] (optional)				Publishing repository (default, <scalpel_results> )
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
		if( params.bam == null || params.bai == null || params.dge_matrix == null)
			params.bam = "${params.folder_in}/outs/possorted_genome_bam.bam"
			params.bai = "${params.folder_in}/outs/possorted_genome_bam.bam.bai"
			params.dge_matrix = "${params.folder_in}/outs/filtered_feature_bc_matrix.h5"
			params.barcodes = "${params.folder_in}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
	}
}

if ( params.barcodes == null )
	error "Provide valid barcodes file ! [--barcodes]"




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
		 Minimal distance of Ip from isoform 3'ends (optional): ${params.isoform_end_ip_threshold}
		 Max cpus [--cpu_defined] (default 50): ${params.cpu_defined}
		 subsampling [--subsampling]: ${params.subsampling}
		 theshold fraction gene [--gene_fraction]: ${params.gene_fraction}
		 binsize fragment probability [--binsize]: ${params.binsize}

		 Publishing repository [--publish_rep] (optional): ${params.publish_rep}

		 """
.stripIndent()






/*****************/
/* PREPROCESSING */
/*****************/

/* Computing genomic Annotation file GTF & Internal Priming annotation(1) */
// ********************************************************************

/*S0*/
/* Process salmon quantification files */
process Quantification_processing{
	tag "${quant}"
	publishDir "${params.publish_rep}/", overwrite: true, mode: 'copy'
	input:
	val quant              from params.quant_file
	output:
	file "quant.filtered"   into quant_filtered_ch
	script:
	"""
	python3 ${baseDir}/src/quantification_processing.py -salmon_quant ${quant} -output_file 'quant.filtered'
	"""
}



/*S1*/
process Gtf_annotation_splitting{
	tag "${gtf_file}"
	publishDir "${params.publish_rep}/gtfs/", overwrite: true
	maxForks params.cpu_defined
	input:
	val pythonbin				from params.python_bin
	file gtf_file				from file(params.annot)
	file quantf					from quant_filtered_ch
	output:
	file "*"				into gtf_splitted_ch1, gtf_splitted_ch2, gtf_splitted_ch3 mode flatten
	script:
	"""
	${pythonbin} ${baseDir}/src/gtf_split.py ${gtf_file} ${quantf}
	"""
}



/*S2*/
process Chromosome_processing{
	tag "${chr_file.baseName}"
	publishDir "${params.publish_rep}/exons/", overwrite: true, mode: 'copy'
	maxForks params.cpu_defined
	input:
	val pythonbin						from params.python_bin
	file chr_file						from gtf_splitted_ch1
	val trs_distance					from params.dt_threshold
	val trs_end_distance					from params.dt_exon_end_threshold
	output:
	file "${chr_file.baseName}.exons"			into exons_ch1, exons_ch2, exons_ch3, exons_ch4
	file "${chr_file.baseName}.exons_unique"		into exons_unique_ch
	script:
	"""
	${pythonbin} ${baseDir}/src/gtf_process.py ${chr_file} ${trs_distance} ${trs_end_distance}  ${chr_file.baseName}.exons ${chr_file.baseName}.exons_unique
	"""
}



/*S3*/
process Split_ipdb{
	tag "${ipdbf}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/ipdb_splitted", overwrite: true
	input:
	file ipdbf		from file(params.ipdb)
	output:
	file "*"	into ipdb_splitted mode flatten
	script:
	"""
	#split by chromosome
	gawk '{print > \$1".ipdb_unsorted"}' ${ipdbf}
	"""
}


process Internal_priming_annotation{
	tag "${ipdbf}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/ipdb_splitted", overwrite: true, mode: 'copy'
	input:
	val pythonbin		from params.python_bin
	val dthr			from params.isoform_end_ip_threshold
	file ipdbf			from ipdb_splitted.collect()
	val bedmap_bin		from params.bedmap_bin
	file exfile			from exons_ch1
	output:
	file "*.ipdb_annotated"		into ipdb_annotated
	script:
	"""
	#select columns and ordering ipdb file
	gawk -v OFS="\\t" '{print \$1,\$2,\$3,\$6}' ${exfile.baseName}".ipdb_unsorted"  | sort -k1,1 -k2,2n -k3,3n > ${exfile.baseName}".ipdb_sorted"
	#select columns and ordering exfile
	gawk -v OFS="\\t" '{print \$1,\$2,\$3,\$12}' ${exfile} | sort -k1,1 -k2,2n -k3,3n > ${exfile}"_sorted"

	#Mapping
	${bedmap_bin} --echo --echo-map-id --fraction-ref 0.8 --delim "\t" ${exfile.baseName}".ipdb_sorted" ${exfile}"_sorted" > ${exfile.baseName}".ipdb"

	#Annotate ip positions
	${pythonbin} ${baseDir}/src/annotate_ips.py ${exfile.baseName}".ipdb" ${exfile} ${dthr} ${exfile.baseName}".ipdb_annotated"
	"""
}



/*S4*/
process Bam_splitting{
	tag "${chr_id.baseName}"
	publishDir "${params.publish_rep}/reads/bams/", overwrite: true, mode: 'copy'
	maxForks params.cpu_defined
	input:
	val ncores			from params.cpu_defined
	val subsamp			from params.subsampling
	val seqtype			from params.sequencing
	val samtoolbin			from params.samtools_bin
	val pythonbin			from params.python_bin
	file chr_id			from gtf_splitted_ch2
	file bam_file			from file(params.bam)
	file bai_file			from file(params.bai)
	file barcodes_file		from file(params.barcodes)
	output:
	file "${chr_id}.ebam"	into exonic_bams_ch, exonic_bams_ch2

	script:
	if( params.sequencing == "dropseq" )
		"""
		#Dropseq
		#Filter reads and split by chromosome
		zcat -f ${barcodes_file} > bc.txt
		${samtoolbin} view -@ ${ncores} -b --subsample ${subsamp} ${bam_file} ${chr_id.baseName} -D XC:bc.txt --keep-tag "XC,XM" > ${chr_id}.ebam
		"""
	else
		"""
		#10X SEQ
		zcat -f ${barcodes_file} | sort -u > bc.txt
		${samtoolbin} view -@ ${ncores} -b --subsample ${subsamp} ${bam_file} ${chr_id.baseName} -D CB:bc.txt --keep-tag "CB,UB" > ${chr_id}.ebam
		"""
}



/*S5*/
process Bed_formatting{
	tag "${bamf.baseName}"
	publishDir "${params.publish_rep}/reads/beds/", overwrite: true
	maxForks params.cpu_defined
	input:
	file bamf				from exonic_bams_ch
	val pythonbin				from params.python_bin
	output:
	file "${bamf.baseName}.bed"		into bed_ch
	script:
	"""
	#Convertion of bam files to bed files (1)
	bam2bed --all-reads --split --do-not-sort < ${bamf} > temp_file

	#Extract relevant metadata infos and sort for overlapping purposes (2)
	${pythonbin} ${baseDir}/src/subset_metatada.py temp_file ${bamf.baseName}.bed
	"""
}



/*S6*/
process Bed_exons_overlapping{
	tag "${bed.baseName} /${fraction_rd}"
	publishDir "${params.publish_rep}/reads/overlap/", overwrite: true
	maxForks params.cpu_defined
	input:
	val bedmapbin			from params.bedmap_bin
	file bed			from bed_ch
	file exfile			from exons_ch2.collect()
	val fraction_rd			from params.fraction_read_overlapping
	output:
	file "${bed.baseName}.bed_gtf"	into bed_gtf_ch
	script:
	"""
	#mapped reads on exons bed (1)
	gawk -v OFS="\\t" '{if(NR>1) print \$1,\$2,\$3,\$12}' ${exfile} | sort -k1,1 -k2,2n -k3,3n > ${bed.baseName}".exsorted"
	gawk -v OFS="\\t" '{if(NR>1) print}' | sort -k1,1 -k2,2n -k3,3n ${bed} > ${bed.baseName}"_readsorted"
	${bedmapbin} --echo --echo-map-id --fraction-ref ${fraction_rd} --delim '\\t' ${bed.baseName}"_readsorted" ${bed.baseName}".exsorted" > ${bed.baseName}.bed_gtf
	"""
}


/*S7*/
process Bed_exons_filtering{
        tag "${ebed.baseName}"
        publishDir "${params.publish_rep}/reads/overlap/", overwrite: true
        maxForks params.cpu_defined
        input:
        val pythonbin					from params.python_bin
        file ebed					from bed_gtf_ch
        file exon_file					from exons_ch3.collect()
        val trs_distance				from params.dt_threshold
        output:
        file "${ebed.baseName}.bed_gtf_filtered"        into bed_gtf_filtered
        script:
        """
        ${pythonbin} ${baseDir}/src/filtering.py ${ebed} ${ebed.baseName}.exons ${trs_distance} ${ebed.baseName}.bed_gtf_filtered
        """
}



/*S8*/
process internalp_filtering{
	tag "${ebed.baseName}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/reads/ip_filtered/", overwrite: true, mode: 'copy'
	input:
	val pythonbin								from params.python_bin
	file all_ips 								from ipdb_annotated.collect()
	file ebed									from bed_gtf_filtered
	output:
	file "${ebed.baseName}.ip_filtered" 					into exip_mapped, exip_mapped2
	file "${ebed.baseName}.ip_filtered_unique" 				into fragment_filtered_uniq_ch
	file "${ebed.baseName}.rid"								into rid_files_ch
	script:
	"""
	#filter internal priming associated reads
	${pythonbin} ${baseDir}/src/ip_filtering.py ${ebed} ${ebed.baseName}.ipdb_annotated ${ebed.baseName}.ip_filtered_tmp ${ebed.baseName}.ip_filtered_unique_tmp ${ebed.baseName}.rid_tmp
	#sorting and drop duplicates
	sort -u ${ebed.baseName}.ip_filtered_tmp > ${ebed.baseName}.ip_filtered
	sort -u ${ebed.baseName}.ip_filtered_unique_tmp > ${ebed.baseName}.ip_filtered_unique
	sort -u ${ebed.baseName}.rid_tmp > ${ebed.baseName}.rid
	"""
}


/*S9*/
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
	${r_bin} ${baseDir}/src/compute_prob.R "\$PWD/" ${gfrac} ${bins} BINS_PROB.txt BINS_PROB.pdf
	"""
}



/*S10*/
process Fragment_probabilities{
	tag "${ebed.baseName}"
	maxForks params.cpu_defined
	publishDir "${params.publish_rep}/reads/cells/chr_cells/", overwrite: true
	maxForks params.cpu_defined
	input:
	file ebed			from exip_mapped
	file prob_file 		from probability_ch
	output:
	file "*.frag_prob"	into frags_channel
	script:
	"""
	#Write cell Files
	python3 ${baseDir}/src/fragment_probabilities.py ${ebed} ${prob_file} ${ebed.baseName}.frag_prob
	"""
}


/*S11*/
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
	Rscript ${baseDir}/src/merge_and_writecells.R "."
	"""
}


/*S12*/
process EM_algorithm{
	tag "${cell_file}"
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



/*S13*/
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


/*S14*/
process DGE_generation{
	tag "${dge_mat}"
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



/*S15*/
process Filter_BAMS{
	tag "${chr.baseName}"
	maxForks params.cpu_defined
	input:
	file all_bams			from exonic_bams_ch2.collect()
	file rid				from rid_files_ch.collect()
	file chr				from gtf_splitted_ch3
	val samtoolbin			from params.samtools_bin
	output:
	file "${chr.baseName}.ebamf"		into filtered_bams_ch
	script:
	"""
	#get reads_id from ip filtered reads
	#filter bam file and indexing
	gawk '{print \$1}' ${chr.baseName}.rid | sort -u > read.ids
	${samtoolbin} view -b -N read.ids ${chr.baseName}.ebam > ${chr.baseName}.ebamf
	${samtoolbin} index ${chr.baseName}.ebamf
	"""
}



/*S16*/
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
	${samtoolbin} merge -@ 10 -f -o final1.bam ${all_bams}
	${samtoolbin} sort -@ 10 final1.bam > final.bam
	${samtoolbin} index final.bam
	"""
}



