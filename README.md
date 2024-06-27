SCALPEL , a NextFlow based tool for the quantification of transcript isoforms at single-cell resolution
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
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/SPERMATOGENESIS/SCALPEL_pipeline.png" alt="SCALPEL" >
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

### Installation

1. Clone the repo and enter in the folder
```sh
> git clone https://github.com/p-CMRC-LAB/SCALPEL.git
```
2. Install the required packages using the requirement.txt file in the SCALPEL folder
```sh
> conda env create -f SCALPEL/requirements.yml
> conda activate scalpel_env
```
3. Within the CONDA environnement, launch R and install the following R packages
```
> install.packages(c("stringi", "Seurat"))
```
   
Another solution (if conda installation takes long) can be to create a Conda environment, install Mamba (faster implementation of Conda) and install the packages using mamba:
```sh
> mamba env create --file SCALPEL/requirements.yml
```

## SCALPEL usage

### EX: Analysis on Lukassen et al dataset

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


### Downstream analysis

See [SCALPEL Wiki](https://github.com/p-CMRC-LAB/SCALPEL/wiki)


## Contact

Your Name - [@twitter](https://twitter.com/aerodx5) - fake@idibell.cat

Project Link: [SCALPEL Github](https://github.com/p-CMRC-LAB/SCALPEL)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Paper

Paper writing ongoing...
