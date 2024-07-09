# SCALPEL , a nextflow based tool for the quantification of isoforms at single-cell resolution

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

**SCALPEL** is a toolkit for the characterization of Alternative ployadenylation (APA) and the Isoform quantification using 3' tag-based scRNA-seq data.

Use the `BLANK_README.md` to get started.


<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

## Prerequisites

- [Conda](https://www.anaconda.com/) package tool or [Mamba](https://github.com/mamba-org/mamba) (Fast reiplementation of conda)

## Installation

1. Clone the repo and enter in the folder
```sh
> git clone https://github.com/p-CMRC-LAB/SCALPEL.git
```
2. Install the required packages using the requirement.txt file in the SCALPEL folder
```sh
> conda env create -f SCALPEL/requirements.yml
> conda activate scalpel_conda
```
3. Within the CONDA environnement, install the R package Seurat v5
```
> Rscript -e 'remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)'
```
   
Another solution (if conda installation takes long) can be to create a Conda environment, install Mamba (faster implementation of Conda) and install the packages using mamba:
```sh
> mamba env create --file SCALPEL/requirements.yml
> mamba activate scalpel_conda
> Rscript -e 'remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)'
```

## SCALPEL Usage
[](https://github.com/p-CMRC-LAB/SCALPEL/edit/dev/README.md#usage)

### Input files
For running, SCALPEL requires to provide specific input files path & parameters:

-   _SAMPLE path_  **[--samplesheet]**
-   _FASTA transcriptome reference path_  **[--transcriptome]**
-   _GTF annotation file path_  **[--gtf]**
-   _Internal priming annotation file for the organism_  **[--ipdb]**  (download link below)
-   _Sequencing type (chromium or dropseq)_  **[--sequencing]**

  1.  Required, **[\-\-samplesheet]**: Provide within a  **CSV**  (ex: samplesheet.csv) file the following paths : 
    -(In case of 10X based scRNA-seq sample [--sequencing  **chromium**] or DropSeq based scRNA-seq sample [--sequencing  **dropseq**]):    
```sh
<SAMPLE_NAME>,<FASTQ1_FILE_PATH>,<FASTQ2_FILE_PATH>,<10X_CELLRANGER_REPOSITORY_PATH>
```

```sh
<SAMPLE_NAME>,<FASTQ1_PATH>,<FASTQ2_PATH>,<BAM_PATH>,<BAM_INDEX_PATH>,<DGE_PATH>
```

  2. Optional, **[\-\-barcodes]**: Provide within a **CSV** (ex: barcodes_whitelist.csv) for each input sample, a barcode whitelist file path:
```sh
<SAMPLE_NAME>,<BARCODE_WHITELIST_FILE_PATH>
```

  3. (Optional, **[\-\-clusters]**: Provide within a **CSV** (ex: clusters.csv) for each input sample, a barcode whitelist file path:
```sh
<SAMPLE_NAME>,<BARCODE_CLUSTERS_FILE_PATH>
```

### Annotation files

- **Genome reference files**
Following the organism for the study, the reference genome annotation files (**[--gtf]**, GTF) and reference transcript sequence (**[--transcriptome]**, FASTA) can be downloaded on [GENCODE](https://www.gencodegenes.org/) repository.

- **Internal priming files which reference all the internal priming positions**
- [(Human) Internal priming annotation - GRCh38](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/GRCh38_2020_A_polyA.track.tar.gz) 
- [(Mouse) Internal priming annotation - mm10](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/mm10_polya.track.tar.gz) 


### SCALPEL execution
After activating the _scalpel_conda_ CONDA environment, SCALPEL can be executed using Nextflow: 
- All the computational ressource required for the execution by Nextflow can be defined within the _**SCALPEL/nextflow.config**_ file:
```
/* Define Nextflow configuration settings for SCALPEL pipeline execution */
/* ===================================================================== */

/* -> Processes
================*/
/* Enter here the desired Nextflow parameters for execution (see https://www.nextflow.io/docs/latest/index.html)*/

/* global parameters for execution */

executor {
    name = 'local' /* Adjust as needed (local/slurm/sge/etc...) */
    cpus = 60      /* Adjust as needed */
}

/* Parameters by process */

process {
    time = '120 m'
        withLabel: big_mem {
            cpus = 4
        }
        withLabel: small_mem {
            cpus = 2
        }
}
```

- **Execution**

1. Configurate sample file **[--samplesheet]**:
```sh
> cat samplesheet.csv
SRR6129050,SRR6129050_S1_L001_R1_001.fastq.gz,SRR6129050_S1_L001_R2_001.fastq.gz,/data/fake_data/DATAS/GSE104556/SAMPLES/SRR6129050/
SRR6129051,SRR6129051_S1_L001_R1_001.fastq.gz,SRR6129051_S1_L001_R2_001.fastq.gz,/data/fake_data/DATAS/GSE104556/SAMPLES/SRR6129051/
```

2. (**Optional**), Configurate barcode whitelist file **[--barcodes]**:
```sh
> cat barcodes_whitelist.csv
SRR6129050,/data/fake_data/NEW/SCALPEL/10X/SRR6129050_curatedBarcodes.txt
SRR6129051,/data/fake_data/NEW/SCALPEL/10X/SRR6129051_curatedBarcodes.txt

> head /data/fake_data/NEW/SCALPEL/10X/SRR6129050_curatedBarcodes.txt
AAACCTGAGCTTATCG-1
AAACCTGGTTGAGTTC-1
AAACCTGTCAACGAAA-1
AAACGGGCACAGGTTT-1
AAACGGGTCATTTGGG-1
...
```

3. **Running**
```sh
> nextflow run -resume SCALPEL/main.nf \
    --sequencing chromium \
    --samplesheet samplesheet.csv \
    --transcriptome gencode.vM10.transcripts.fa \
    --gtf gencode.vM21.annotation.gtf \
    --ipdb mm10.polyA.track \
    --barcodes barcodes_whitelist.csv \ 
```

### Results

During its execution, **SCALPEL** shows interactively in the Console all the information about the executed tasks:
Then, a _**./results**_ folder is generated encompassing all the final and intermediated files generated by **SCALPEL** that can be used for further downstream analysis:

_**./results/final_results**_ contains:
  - the final differential isoform usage table between the input samples
  - the Seurat object encompassing the iDGE of the input samples
  - the iDGE table of each input sample

```sh
> cat ./results/final_results
DIU_table.csv  iDGE_seurat.RDS.  SRR6129050/SRR6129050_APA_DGE.txt   SRR6129051/SRR6129051_APA_DGE.txt
```

Information about the collapsed isoforms by SCALPEL (see paper) can be found in _**./results/annotation_processing/isoform_processing**_

```sh
> cat ./results/annotation_processing/isoform_processing
chr10_collapsed_isoforms.txt  chr13_collapsed_isoforms.txt  chr16_collapsed_isoforms.txt ...
```

**NB:**
_**Important: The chromosome names must to be consistent between the BAM file and the annotation files. If the BAM file contains the '_chr_' character, the GTF and FASTA annotation files, and the internal priming reference annotation should contains the '_chr_' character, and inversely !**_

_**The internal priming reference files provided contains by default the '_chr_' character in the chromosome names !!**_

_**Be careful to delete the _work_ directory containing nextflow temporary files** when scalpel runs all its processs sucessfully and you don’t plan to relaunch scalpel with modified parameters. (This folder can fill an high memory physical space depending of the size of input files analyzed)_


## More about SCALPEL Usage

- [Example of SCALPEL application on 10X scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-10X-scRNA%E2%80%90seq)
- [Example of SCALPEL application on DropSeq scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-DropSeq-scRNA%E2%80%90seq)


## Downstream analysis

See [SCALPEL Wiki](https://github.com/p-CMRC-LAB/SCALPEL/wiki)


# Contact

Your Name - [@twitter](https://twitter.com/aerodx5) - fake@idibell.cat

Project Link: [SCALPEL Github](https://github.com/p-CMRC-LAB/SCALPEL)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


# Paper
[Access SCALPEL_PAPER](https://www.biorxiv.org/content/10.1101/2024.06.21.600022v1)

Quantification of transcript isoforms at the single-cell level using SCALPEL \
Franz Ake, Sandra M. Fernández-Moya, Marcel Schilling, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass \
bioRxiv 2024.06.21.600022; doi: https://doi.org/10.1101/2024.06.21.600022
