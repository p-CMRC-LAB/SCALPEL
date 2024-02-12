SCALPEL , a NextFlow based tool for the characterization of Alternative polyadenylation  at single-cell resolution
======================================================================================


<!-- PROJECT LOGO -->
<!-- <br />
<div align="center">
  <a href="">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/Dessin_scalpel.svg/1200px-Dessin_scalpel.svg.png" alt="SCALPEL" width="300" height="300">
  </a>
</div>-->

<div align="right">
  <a href="">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Scalpel_workflow2.svg" alt="SCALPEL" >
  </a>
</div>




<!-- ABOUT THE PROJECT -->
## About The Project

**SCALPEL** is a toolkit for the characterization of Alternative ployadenylation (APA) and the Isoform quantification at single-cell resolution.

Use the `BLANK_README.md` to get started.


<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

- [Conda](https://www.anaconda.com/) package tool or [Mamba](https://github.com/mamba-org/mamba) (Fast reiplementation of conda)
- [NextFlow package >= 22.10](https://www.nextflow.io/)
- [Salmon package >= 1.9.0](https://salmon.readthedocs.io/en/latest/index.html)

### Installation

1. Clone the repo and enter in the folder
```sh
> git clone https://github.com/p-CMRC-LAB/SCALPEL.git
```
2. Install the required packages using the requirement.txt file in the SCALPEL folder
```sh
> conda env create -f SCALPEL/requirements.yml
```
   
   Another solution (if conda installation takes long) can be to create a Conda environment, install Mamba (faster implementation of Conda) and install the packages using mamba:
```sh
> mamba env create --file SCALPEL/requirements.yml
```

## Usage

### Analysis on Lukassen et al dataset

#### Input files

For the need of the analysis in this vignette, the data used is a 10X dataset from the study from [Winterpacht A, Lukassen](https://pubmed.ncbi.nlm.nih.gov/30204153/) on Mouse. The demo data is 10X processed folder containing a BAM file with aligned reads from the published data [GSE104556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556) 
from GEO.

##### Sample

Scalpel can be applied on the 10X result folder, containing the default files (_possorted_genome_bam.bam, possorted_genome_bam.bam.bai, filtered_feature_bc_matrix.h5, barcodes.tsv.gz, ..._)

[10X_FOLDER](https://drive.bio.idibell.cat/index.php/s/eoYbCKA48eZXMDK)

##### Internal priming files which reference all the internal priming positions

[(Human) Internal priming annotation - GRCh38](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/GRCh38_2020_A_polyA.track.tar.gz)

[(Mouse) Internal priming annotation - mm10](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/mm10_polya.track.tar.gz)

##### Reference genome annotation files

The reference genome file annotation (GTF) and the reference transcript sequence (FASTA) can be downloaded on the [GENCODE](https://www.gencodegenes.org/mouse/release_M10.html) website associated to the organism.

[(Mouse) GTF_annotation_file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz)

[(Mouse) FASTA_reference_file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.transcripts.fa.gz)


#### Execution

For running, SCALPEL requires to provide specific input files path & parameters:
 - SAMPLE files folder path (\*.bam/\*.bai/\*.barcodes/\*.counts.txt)  **[\-\-samples]**
 - FASTQs files folder path (\*.fastq.gz) **[\-\-reads]**
 - FASTA transcriptome reference path **[--transcriptome]**
 - GTF annotation file path **[\-\-gtf]**
 - Internal priming annotation file **[\-\-ipdb]**
 - Sequencing type (chromium or dropseq) **[\-\-sequencing]**

The script can be applied into a folder containing several samples to be analyzed. All the samples files with the defined extensions mentionned above (\*.bam/\*.bai/\*.barcodes/\*.counts.txt) will be processed according each sample name.

**Differential isoform usage analysis on clusters:**
Differential analysis usage is done in the downstream analysis step following the isoform quantification by default considering the different samples provided. To perform a differential analysis on defined cells cluster, use the optional argument **[\-\-clusters]**

#### Processing
You can print the Scalpel help documentation by running the following command

```sh
> nextflow run -resume SCALPEL/scalpel.nf --help
```
  ===============================
	SCALPEL - NF  P I P E L I N E
	===============================

	Execution:
	Ex: nextflow run -resume scalpel.nf --sequencing <Sequencing type>
	 --samples <SAMPLE files folder path>
	 --reads <FASTQs files folder path> --transcriptome <FASTA transcriptome reference path>
	 --annot <GTF annotation file path> --ipdb <Internal priming annotation file> 
	
	Input files:
    - Annotation required files(required):
        - transcriptome reference [--transcriptome]
        - annotation GTF reference [--gtf]
        - internal priming annotation [--ipdb]
      
    - Reads processing files (required):
        - samples files [--samples]
        - fastqs files [--reads]
    
    - Params:
        Required:
        - sequencing type (required): ${params.sequencing}

        Optional:
        - transcriptomic distance threshold [--dt_threshold] (optional, default 600bp)
        - transcriptomic end distance threhsold [--dt_exon_end_threshold] (optional, default 30bp)
        - minimal distance of Ip from isoform 3'ends (optional, default 60bp)
        - params.threads [--threads] (default 30)
        - params.cpus [--cpus] (default 30)

**Important: The chromosome names must to be consistent between the BAM file and the annotation files. If the BAM file contains the '_chr_' character, the GTF and FASTA annotation files, and the internal priming reference annotation should contains the '_chr_' character, and inversely !**

**The internal priming reference files provided contains by default the '_chr_' character in the chromosome names !!**

**Be careful to delete the work directory containing nextflow temporary files** when scalpel runs all its processs sucessfully and you don’t plan to relaunch scalpel with modified parameters. (This folder can fill an high memory physical space depending of the size of input files analyzed)


### Results - Downstream analysis

#### Fragments distribution on transcriptomic space
During the Nextflow execution or at the end, an image file (BINS_PROB.jpeg) showing the distribution of the fragments (reads associated to the same Barcode and UMI tag) in the transcriptomic space is generated in the  **results/reads/probability/**. Depending of the experiment, the  **[–-gene_fraction]** , **[–-dt_threshold]**  or **[–-binsize]** can be 
modified in order to get a good fit between the fragment counts distribution and the empiric distribution (reads counts by intervals). This file is located in  **results/reads/probability/BINS_PROB.pdf**.

<br />
<div align="left">
  <a href="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/frag_coverage.png">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/frag_coverage.png" alt="SCALPEL" >
  </a>
</div>

#### Single-cell analysis of quantified isoforms

We used [Seurat](https://satijalab.org/seurat/)  for the single cell analysis and APA characterization. We perform a single-cell analysis of the quantified isoform matrix. 

1. Processing of the Scalpel DGE counts
```r
library(Seurat)
library(data.table)
library(dplyr)
library(Gviz)
#source the script scalpelin into the scalpel folder
source("SCALPEL/src/scalpelib_library.R")

#Paths
TRANSCRIPT_COUNTS_MATRIX_SCALPEL_FOLDER_PATH  <- "results/reads/apa_dge/APADGE.txt"

#Create a Seurat object of the estimated transcript counts matrix
TRANSCRIPT_COUNTS_MATRIX <- read.table(file=TRANSCRIPT_COUNTS_MATRIX_SCALPEL_FOLDER_PATH, sep="\t", header=T, row.names=1)
colnames(TRANSCRIPT_COUNTS_MATRIX) <- stringr::str_replace(colnames(TRANSCRIPT_COUNTS_MATRIX), "\\.","\\")
sc.obj <- Seurat::CreateSeuratObject(counts=TRANSCRIPT_COUNTS_MATRIX,  project='SRR6129050_TRANSCRIPT_EXP')

sc.obj
An object of class Seurat 
35328 features across 1300 samples within 1 assay 
Active assay: RNA (35328 features, 0 variable features)
```

2. Quality filtering

Different approaches for the filtering of low quality cells can be realized at this step. We will subset the Seurat object based on a set on curated cell barcodes provided into the scalpel folder, which are annotated with corresponding cell types.
```r
#A file of curated barcode cells and cell types associated from Lukassen et al paper can be found in SCALPEL/docs
CELL_BARCODES_PATH = "SCALPEL/docs/Lukassen_curated_cells.csv"

#Filter with list of cells in annotations
ctypes_barcodes = fread(CELL_BARCODES_PATH)
sc.obj = subset(sc.obj, cells = ctypes_barcodes$Barcode)

sc.obj
An object of class Seurat 
35328 features across 1237 samples within 1 assay 
Active assay: RNA (35328 features, 0 variable features)
```

3. Normalization and data reduction

Once the removing of low quality cells is effective into the Seurat object of Transcript expression from Scalpel, then, we can perform the Normalization, Data reduction and visualization of cells isoform expression.
```r
sc.obj = NormalizeData(sc.obj)
sc.obj = FindVariableFeatures(sc.obj, nfeatures = 2000)
sc.obj = ScaleData(sc.obj)
sc.obj = RunPCA(sc.obj)
pc_choice = 11

#UMAP
sc.obj = RunUMAP(sc.obj,dims = 1:pc_choice)

#Annotate cells with Lukassen et al cell types annotation
sc.obj$Barcode = colnames(sc.obj)
sc.obj$CellType = left_join(sc.obj@meta.data, ctypes_barcodes)$CellType
DimPlot(sc.obj, group.by = "CellType", label = T, label.size = 5)
```

<br />
<div align="center">
  <a href="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/UMAP_celltypes.png">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/UMAP_celltypes.png" alt="SCALPEL" >
  </a>
</div>


4. Differential isoform characterization

Let’s visualize an example of a gene with its different isoform expression (**Cnbp Gene**)

```r
#let's look for differential expression across all the conditions or in some specific condition
# For example, in the CS v RS1 cell types...
isoforms_markers = Find_isoforms(sc.obj %>% subset(CellType %in% c("CS","RS1")), condition = "CellType")
View(isoforms_markers)

# let's visualize one gene with differentially expressed isoforms expression
gene_in = "Cnbp"
gene_tab = isoforms_markers %>% filter(gene == gene_in)
gene_tab
			    CS     ES    RS1    RS2    SC1    SC2 gene p_value p_value.adjusted         gene_tr transcript
Cnbp***Cnbp-004 263.42 265.52  65.51 174.62  28.84  91.70 Cnbp       0                0 Cnbp***Cnbp-004   Cnbp-004
Cnbp***Cnbp-005  28.75  22.20 300.16  35.68 341.28 974.28 Cnbp       0                0 Cnbp***Cnbp-005   Cnbp-005

Nebulosa::plot_density(sc.obj, features=gene_tab$gene_tr)
FeaturePlot(sc.obj, features=gene_tab$gene_tr, order=T, pt.size=0.9, label=T, label.size=5, ncol=3)

```

<br />
<div align="center">
  <a href="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnbp_expression.png">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnbp_expression.png" alt="Cnbp expression" >
  </a>
</div>

<br />
<div align="center">
  <a href="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnbp_expression2.png">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnbp_expression2.png" alt="Cnpbp expression" >
  </a>
</div>


We can look for the Isoform coverage expression between the specific clusters. This visualization require to split the filtered BAM file produced after the scalpel execution according to the barcodes attached to selected clusters.

```r

#Let's split the input BAM file use for the analysis 
samtools_path = "/Users/franz/opt/anaconda3/envs/scalpel_env/bin/samtools"
BCTAG = "CB:"
#indicate the filtered BAMFILE preseent in the SCALPEL results folder
BAMFILE_PATH = "results/reads/filtered_bam/final.bam"
GTF_PATH = "~/CEPH/users/fake/SHARED_files/gencode.vM10.annotation.gtf"

# let's split the BAM file according to the clusters and filter the associated reads
lapply(c("CS1", "CS2", "SC1", "SC2", "RS1", "RS2"), function(CONDTYPE){
  print(CONDTYPE)

  #1)let's get the barcodes associated to each condition and write
  BC_FILE = paste0(CONDTYPE, ".txt")
  (sc.obj@meta.data %>% filter(Cluster == CONDTYPE)) %>% select(Barcode) %>% fwrite(BC_FILE)

  #2)subset bam file according to barcode cells
  FILTER_EXP = paste0(samtools_path, " view -b -D ", BCTAG, BC_FILE, " ", BAMFILE_PATH, " > ", CONDTYPE, 
  ".bam")
  system(FILTER_EXP)

  #3)Index bam file
  INDEX_EXP = paste0(samtools_path, " index ", CONDTYPE, ".bam")
})

```

Once the sample BAM file splitted in accordance with the different cell types, we can visualize the coverage of the reads associated to the differentially expressed isoform  on the filtered BAM of each cell type.

```r
# Look coverage on generated BAM files
genome_tb = rtracklayer::import(GTF_PATH)
bamfiles = c("CS.bam","ES.bam","RS2.bam","RS1.bam","SC2.bam","SC1.bam")
bamnames = c("CS","ES","RS2","RS1","SC2","SC1")

#Coverage Vizualisation 1 (all the isoforms of the gene processed by scalpel)
#NB: scalpel colapse in its execution all the isoforms with similar exons coordinates in the 3'end
CoveragePlot(genome_gr = genome_tb,
             gene = gene_in,
             bamfiles = bamfiles,
             bamnames = bamnames,
             seuratobj = sc.obj %>% subset(CellType %in% bamnames),
             condition = "CellType", genome = "mm10")
```

<br />
<div align="center">
  <a href="">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnbp_all_isoforms.png" alt="Cnbp isoforms expression" >
  </a>
</div>


We can subset the visualization to the only differential isoform usage detected by Scalpel
```r
#Coverage Vizualisation 2 (Subset only the isoforms with differential expression in cell types)
CoveragePlot(genome_gr = genome_tb,
             gene = gene_in,
             bamfiles = bamfiles,
             bamnames = bamnames,
             seuratobj = sc.obj %>% subset(CellType %in% bamnames),
             condition = "CellType", genome = "mm10",
             transcripts_tofilter = gene_tab$transcript)
```

<br />
<div align="center">
  <a href="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnpb_isoforms.png">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnpb_isoforms.png" alt="Cnbp isoform expression" >
  </a>
</div>


We can visualize the relative expression of the differentially expressed isoform in each cell type by using the function _plot_relativeExp()_
```r
#Let's observe Cnbp diff expressed isoforms in the CellTypes
plot_relativeExp(sc.obj, features = gene_tab$gene_tr, group.var = "CellType")
```

<br />
<div align="center">
  <a href="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnbp_relative_expression.png">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/Cnbp_relative_expression.png" alt="Cnbp isoform expression" >
  </a>
</div>



## Contact

Your Name - [@twitter](https://twitter.com/aerodx5) - fake@idibell.cat

Project Link: [SCALPEL Github](https://github.com/p-CMRC-LAB/SCALPEL)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Paper

Paper writing ongoing...
