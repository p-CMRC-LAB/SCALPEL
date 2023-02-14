SCALPEL: Characterization of Alternative polyadenylation
========================================================
\F. AKE, M. Plass, IDIBELL, Barcelona, SPAIN


Abstract
--------

This vignette provides an example on how to use SCALPEL to analyze 10X
or Dropseq single-cell RNA sequencing data.

Installation
------------

SCALPEL can be installed from github:

1. Clone the repo and enter in the folder

::

   git clone https://github.com/p-CMRC-LAB/SCALPEL.git

2. Enter into the scalpel folder and Install the required packages into
your base Conda environment or a new Conda environment dedicated to
SCALPEL

::

   conda install -c bioconda -c defaults -c conda-forge --file requirements.txt


Another solution (if conda installation takes long) can be to create a Conda environment, install Mamba (faster implementation of Conda) and install the packages using mamba:

::

   conda create --name scalpel_env python=3.9
   conda activate scalpel_env
   conda install -c conda-forge mamba
   mamba install -c bioconda -c defaults -c conda-forge --file SCALPEL/requirements.txt
  

Usage
-----

Prerequites
~~~~~~~~~~

Demo data download
^^^^^^^^^^^^^^^^^^

For the need of the analysis in this vignette, the data used is a 10X
dataset from the study from `Winterpacht A,
Lukassen <https://pubmed.ncbi.nlm.nih.gov/30204153/>`__ on Mouse. The
demo data is 10X processed folder containing a BAM file with aligned
reads from the published data
`GSE104556 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556>`__
from GEO. The FASTQ files related to the dataset are required too in the
SCALPEL analysis

`10X_FOLDER <https://drive.bio.idibell.cat/index.php/s/tQnxyecoiCB5qH2>`__

`FASTQ_FILES <https://drive.bio.idibell.cat/index.php/s/painctRPE5jsd6Q>`__

Internal priming reference file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the usage of SCALPEL, it is required to provide IP_REFERENCE_file
which reference all the internal priming positions.

`IP_REFERENCE_FILE <https://drive.bio.idibell.cat/index.php/s/EBMmiBGCEWBdmE7>`__
(HUMAN)

`IP_REFERENCE_FILE <https://drive.bio.idibell.cat/index.php/s/JaaYDaffZHWbiWn>`__
(MOUSE)

Reference genome files
^^^^^^^^^^^^^^^^^^^^^^

The reference genome file annotation (GTF) and the reference transcript
sequence (FASTA) can be downloaded on the
`GENCODE <https://www.gencodegenes.org/mouse/release_M10.html>`__
website.

`REF_GTF <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz>`__
(MOUSE)

`REF_FASTA <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.transcripts.fa.gz>`__
(MOUSE)

Preprocessing
~~~~~~~~~~~~

For running, SCALPEL requires a file quantifying the transcript
abundances in bulk. In order to obtain this file we can run the
preprocessing_scalpel.nf script into SCALPEL.

We can print the help documentation for running the SCALPEL
preprocessing

::

   N E X T F L O W  ~  version 22.10.1
   Launching `SCALPEL/preprocessing_scalpel.nf` [crazy_faggin] DSL1 - revision: 7df54d07e6
	===============================
	SCALPEL - N F   P I P E L I N E
	===============================

	Execution:
	- In case of providing 10X cell ranger folder:
	usage: nextflow run -resume SCALPEL/preprocessing_scalpel.nf --sample_names <SAMPLE1,SAMPLE2,...>  --folder_in <FASTQ_FOLDER_PATH> --reference_fasta_transcript <REF_FASTA>

	Output options:
	--sample_names,						Name of the samples to process (same as the FASTQ file names) [required]
	--folder_in,						Path to FASTQ files folder [required]
	--reference_fasta_transcript				Reference FASTA transcript file [required]
	--salmon_index,						Path of salmon index (optional) -- will skip the salmon index processing task

	[--python_bin_path] (optional)				Path to Python bin (default: python3)
	[--salmon_path_bin] (optional)				Path to Salmon bin (default: salmon)
	[--publish_rep] (optional)				Publishing repository (default: preprocessing)
	[--salmon_quant_library_type] (optional)		(default: A)
	[--salmon_quant_threads] (optional)			(default: 10)
	[--cpu_defined] (optional)				(default: 24)
	--tagR1							(default: R1)
	--tagR2							(default: R2)

Be careful that the sample name provided in this command match the sample name of the FASTQ file (Ex: SRR6129050)

::

   > tar -xf FASTQs.tar.gz
   > nextflow run -resume SCALPEL/preprocessing_scalpel.nf --sample_names SRR6129050 --folder_in <EXTRACTED_FASTQ_FOLDER_PATH> -- reference_fasta_transcript <REF_FASTA_PATH>

A file named *quant.filtered* is generated (by default into a **preprocessing** folder) and will be used later by Scalpel.

Scalpel execution
~~~~~~~~~~~~~~~~

You can print the Help documentation for running SCALPEL with the
command

::

   > nextflow run -resume scalpel.nf --help

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
	[--mapq]						have mapping quality >= INT (default, 0)
	[--gene_fraction]					theshold fraction gene
	[--binsize]						binsize fragment probability
	[--publish_rep] (optional)				Publishing repository
	[--chr_concordance]					Character at add in order to match chromosome name in BAM file and the genome reference annotation file

The 10X_folder dataset, and the others reference data files are
extracted, and SCALPEL can be run in this way:

::

   nextflow run -resume scalpel.nf --sequencing chromium --folder_in <10X_FOLDER_PATH> --annot <REG_GTF_PATH> --ipdb <IP_REFERENCE_FILE_PATH> --quant_file preprocessing/quant.filtered

the –-chr_concordance option is specified in the case than the REF_GTF file and the BAM file contain different chromosome names (chr,…),
and the --subsampling option enable to subsample only a fraction of the reads (default: 1 ~ all reads).

A **scalpel_results** folder containing intermediate and final result files is generated during the execution.

**Be careful to delete the work directory containing nextflow temporary files** when scalpel runs all its processs sucessfully and you don't plan to relaunch scalpel with modified parameters. (This folder can fill an high memory physical space depending of the size of input files analyzed)

.. image:: _static/scalpel_run.png
  :width: 1200
  :alt: scalpel_run.png


Results
-------

During the Nextflow execution or at the end, an image file (BINS_PROB.jpeg) showing the distribution of the fragments in the transcriptomic space is generated in the **scalpel_results/reads/probability**. Depending of the experiment, the **[–-gene_fraction]** and   **[–-dt_threshold]** can be modified in order to get a good fit between the fragment counts distribution and the empiric distribution (reads counts by intervals).
This file is located in **scalpel_results/reads/probability/BINS_PROB.jpeg**.

.. image:: _static/reads_distribution.jpeg
  :width: 1200
  :alt: reads_distribution.jpeg

Single-cell Analysis of quantified Isoforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We used `Seurat <https://satijalab.org/seurat/>`__ tool for the single cell analysis and APA characterization. We gonna perform a single-cell analysis of the quantified isoform along a classical single cell analysis using gene expression.

Processing of the SCALPEL DGE count files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: r

	library(Seurat)
	library(dplyr)
	library(data.table)
	library(clustree)
	library(stringr)
	library(patchwork)

	#Paths
	GENE_COUNTS_MATRIX_10X_FOLDER_PATH <- "10X_FOLDER/outs/filtered_feature_bc_matrix/"
	TRANSCRIPT_COUNTS_MATRIX_SCALPEL_FOLDER_PATH <- "scalpel_results/reads/apa_dge/APADGE.txt"

	# Create Seurat object of Gene counts matrix
	GENE_COUNTS_MATRIX <-  Seurat::Read10X(GENE_COUNTS_MATRIX_10X_FOLDER_PATH)
	s.obj <- Seurat::CreateSeuratObject(counts = GENE_COUNTS_MATRIX, project = "SRR6129050_GENE_EXP")

	#Create Seurat object of transcript counts matrix
	TRANSCRIPT_COUNTS_MATRIX <- read.table(file = TRANSCRIPT_COUNTS_MATRIX_SCALPEL_FOLDER_PATH, sep = "\t", header = T, row.names = 1)
	colnames(TRANSCRIPT_COUNTS_MATRIX) = stringr::str_replace(colnames(TRANSCRIPT_COUNTS_MATRIX), "\\.","\\-")
	sc.obj = Seurat::CreateSeuratObject(counts = TRANSCRIPT_COUNTS_MATRIX, project = 'SRR6129050_TRANSCRIPT_EXP')

::

	> s.obj
	# SEURAT_GENE_COUNT_OBJ
	An object of class Seurat 
	32285 features across 1300 samples within 1 assay 
	Active assay: RNA (32285 features, 0 variable features)
	
	> sc.obj
	# SEURAT_TRANSCRIPT_COUNT_OBJ
	An object of class Seurat 
	54938 features across 1300 samples within 1 assay 
	Active assay: RNA (54938 features, 0 variable features)


Quality filtering
'''''''''''''''''

Different approaches for the filtering of low quality cells can be realized at this step. A first approach can be to remove the bad quality cells into the Gene expression Assay resulting from the 10X analysis and then remove the same cell barcodes into the Transcript expression Assay from Scalpel.
A second approach could be to simply remove the bad quality cells directly into the Transcript expression Assay object.

.. code:: r

	#Let's filter our data using the 1st approach
	#Visualization of UMI and Genes counts in the Seurat object unfiltered
	Seurat::FeatureScatter(s.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	#Visualization of UMI and Genes counts in the Seurat object Filtered
	Seurat::FeatureScatter(subset(s.obj, nCount_RNA < 100e3), feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	#Filtering of 10X gene object
	s.obj1 = subset(s.obj, nCount_RNA < 100e3)

	#Let's subset the same cells barcodes into the Seurat object unfiltered from Scalpel
	#Filtering of 10X gene object
	sc.obj1 = subset(s.obj, cells = colnames(s.obj1))
	Seurat::FeatureScatter(subset(sc.obj1, nCount_RNA < 100e3), feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  + ggtitle("UMI_count x Transcript_count - Filtered")


.. image:: _static/QC_plots.png
  :width: 1200
  :alt: QC_plots.png

Normalization and data reduction
''''''''''''''''''''''''''''''''

Once the removing of low quality cells is effective into the Seurat object of Transcript expression from Scalpel, then, we can perform the Normalization, Data reduction and visualization of cells isoform expression.

.. code:: r
	
	#Normalization and Data reduction
	sc.obj1 = NormalizeData(sc.obj1)
	sc.obj1 = FindVariableFeatures(sc.obj1)
	sc.obj1 = ScaleData(sc.obj1)
	sc.obj1 = RunPCA(sc.obj1)
	ElbowPlot(sc.obj1, 50)
	pc_choice = 11
	
	#let's add metadata information about cell types and cell barcodes
	PATH_TO_BARCODES_CTYPES_FILE = "CEPH/benchmark/GSE104556/srr6129050_analysis/srr6129050_raw_10x_files_4/ctypes_barcodes.tsv"
	ctypes_barcodes = fread(PATH_TO_BARCODES_CTYPES_FILE, header = F, col.names = c("barcodes", "ctypes"))
	#let's perform a join to add the ctypes information into the meta.data slot
	sc.obj1$ctypes = left_join(sc.obj1@meta.data, ctypes_barcodes)$ctypes
	
	#let's visualize TSNE
	sc.obj1 = RunTSNE(sc.obj1, dims = 1:pc_choice)
	DimPlot(sc.obj1, reduction = "tsne", label = T, label.size = 7, pt.size = 0.7, group.by = "ctypes") + theme_classic(base_size = 14)

.. image:: _static/UMAP_ctypes.png
  :width: 1200
  :alt: UMAP_ctypes.png


Isoform quantification by clusters
''''''''''''''''''''''''''''''''''

.. code:: r

	#Get Transcripts quantification by clusters
	#get genes
	all_genes = rownames(sc.obj1)
	genes_tr_tab = (all_genes %>% str_split_fixed(pattern = "\\*\\*\\*", n = 2)) %>% data.table()
	genes_tr_tab$gene_tr = rownames(sc.obj1)
	colnames(genes_tr_tab) = c("gene", "transcript", "gene_transcript")
	#genes_tr_tab
	#filter out genes with only one Isoform present
	counts_genes_tab = genes_tr_tab$gene %>% table() %>% data.table() %>% filter(N > 1)
	genes_tr_tab_filtered1 = genes_tr_tab%>% filter(gene %in% counts_genes_tab$.) %>% arrange(gene_transcript)
	genes_tr_tab_filtered1
	#Matrix of counts
	sc.obj1 = ScaleData(sc.obj1, features = genes_tr_tab_filtered1$gene_transcript, do.center = F)

	##Scaling data matrix
	ALL_expression = AggregateExpression(sc.obj1, features = genes_tr_tab_filtered1$gene_transcript,
				     assays = 'RNA', group.by = 'seurat_clusters', verbose = T, slot = 'scale')$RNA %>% data.frame()
	ALL_expression$only_gene = (rownames(ALL_expression) %>% str_split_fixed(pattern = "\\*\\*\\*", n = 2))[,1]
	ALL_expression$gene_tr = rownames(ALL_expression)
	#Split all the table by genes
	ALL_expression_by_GENE = split(ALL_expression, ALL_expression$only_gene)



Let's visualize an example of a gene with its different isoform expression (**Eif4e Gene**)

.. image:: _static/EIF4E_table.png
  :width: 1200
  :alt: EIF4E_table.png


.. code:: r

   Reduce(`+`, lapply(ALL_expression_by_GENE$Eif4e$gene_tr, function(x){
     FeaturePlot(scalpel.seurat.filtered, features = x, pt.size = 0.1, order = T) + theme_classic(base_size = 6)
   })) + plot_layout(ncol = 3)


.. image:: _static/EIF14E_featurePlot.png
  :width: 1200
  :alt: EIF14E_featurePlot.png



Differential Isoform characterization
'''''''''''''''''''''''''''''''''''''

.. code:: r

	#let's look APA differences between the cell types CS1, CS2, ES1, ES2, SC1, SC2, RS1, RS2 and "SPG, Sertoli, Leydig"

	RES = lapply(names(ALL_expression_by_GENE), function(x){
	# print(x)
	a = ALL_expression_by_GENE[[x]]
	a = a[,c(1:9)]
	b = apply(a,2, function(x) x/sum(x))
	c = a[names(which(rowSums(b > 0.2) >= 1)),]

	if(nrow(c) > 1){
	d = suppressWarnings(chisq.test(c))
	c$gene = x
	c$p_value = d$p.value
	return(list(c, d))
	}else{
	return(NULL)
	}
	})
	#delete NULL occurences
	RES = RES[!sapply(RES,is.null)]
	#get tables extraction
	RES_TAB = lapply(RES, function(x) x[[1]])
	RES_TAB = do.call(rbind, RES_TAB)
	#adjust _pvalue
	RES_TAB$p_value.adjusted = p.adjust(RES_TAB$p_value,method = 'fdr')
	#filter
	RES_TAB_SIGNIF = RES_TAB %>% filter(p_value.adjusted < 0.01)
	RES_TAB_SIGNIF$gene_tr = rownames(RES_TAB_SIGNIF)
	RES_TAB_SIGNIF$transcript = str_split_fixed(RES_TAB_SIGNIF$gene_tr,pattern = "\\*\\*\\*",n=2)[,2]
	RES_TAB_SIGNIF = RES_TAB_SIGNIF %>% arrange(p_value.adjusted,gene)
	


Let's visualize some tops differential isoforms expression in the clusters

::

	> head(RES_TAB_SIGNIF, 20)


.. image:: _static/top_transcripts.png
  :width: 1200
  :alt: top_transcripts.png


Scalpel identied an differential expression of these isoform in the cell types analyzed. We can vizualize the Expression of these isoforms by using the **FeaturePlot** function of Seurat.

::

	> Reduce(`+`, lapply(ALL_expression_by_GENE$Cep57l1$gene_tr, function(x){
		Nebulosa::plot_density(sc.obj1, features = x) + theme_classic(base_size = 6)
	})) + plot_layout(ncol = 3)


.. image:: _static/Cep57L1_featplot.png
  :width: 1200
  :alt: Cep57L1_featplot.png


Now, let's look for the mapping of the reads in the input BAM file to see if they are in accordance with the quantification performed by scalpel.

.. code:: r

	library(Gviz)
	library(GenomicRanges)
	library(GenomicFeatures)

	output_path = "/CEPH/users/fake/test/"
	samtoolsbin = "/home/fake/.conda/envs/scalpel_env/bin/samtools"
	bam_file = paste0(output_path,"scalpel_results/reads/filtered_bam/final.bam")

	#Use the source code script in SCALPEL
	source("SCALPEL/src/coverage_visualization.R")
	#import GTF annotation file
	genome_gr = rtracklayer::import('~/CEPH/datas/mm10/gencode.vM10.annotation.gtf')
	genome_gr$transcript_id = paste0(genome_gr$transcript_name,'-',genome_gr$transcript_id)

	#Process metadata and filter clusters
	cell.info = sc.obj1@meta.data
	cell.info$cells = rownames(cell.info)
	cell.info$ctypes = as.character(cell.info$ctypes)
	cell.info$ctypes = str_replace(cell.info$ctypes, "SPG, Sertoli, Leydig", "SPG_Sertoli_Leydig")
	
	#split the bam files by barcodes
	#get barcodes cell associated to cluster1 and Bam file associated
	lapply(unique(cell.info$ctypes), function(x){
	    print(x)
	    fwrite(data.table((cell.info %>% filter(ctypes==x))$cells), file = paste0(output_path,x,".barcodes"), col.names = F,row.names = F)
	    system(paste0(samtoolsbin, " view -b -D CB:", output_path, x,".barcodes ", bam_file, " > ", output_path, x, ".bam"))
	    system(paste0(samtoolsbin, " index ", output_path, x, ".bam"))
	})

Once the sample BAM file splitted in accordance with the different cell types, we can visualize the coverage of the reads for each cell type and BAM file

.. code:: r
	
	#Enter the path of splitted bam file
	bamfiles = c("CS1.bam","CS2.bam","ES1.bam", "ES2.bam", "RS1.bam", "RS2.bam", "SC1.bam", "SC2.bam")
	#Attributes names
	bamnames = c("CS1","CS2","ES1","ES2","RS1","RS2","SC1","SC2")

	gene_in = "Eif4e"
	target = RES_TAB_SIGNIF %>% filter(gene == gene_in)

	#Coverage Plot
	genome_cover(genome_gr = genome_gr[genome_gr$transcript_name %in% target$transcript], bamfiles = bamfiles, bamnames = bamnames, gene_in = gene_in, sample_sizes = table(cell.info$clusters))


.. image:: _static/Coverage_2.png
  :width: 1200
  :alt: Coverage_2.png
