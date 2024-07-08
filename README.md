# SCALPEL , a Nextflow based tool for the quantification of transcript isoforms at single-cell resolution

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

### - Input files
For running, SCALPEL requires to provide specific input files path & parameters:

-   _SAMPLE path_  **[--samplesheet]**
-   _FASTA transcriptome reference path_  **[--transcriptome]**
-   _GTF annotation file path_  **[--gtf]**
-   _Internal priming annotation file for the organism_  **[--ipdb]**  (download link below)
-   _Sequencing type (chromium or dropseq)_  **[--sequencing]**

1.  Required, **[\-\-samplesheet]**: Provide within a  **CSV**  (ex: samplesheet.csv) file the following paths : 
    -   (In case of Dropseq based scRNA-seq sample [--sequencing  **dropseq**])  
        >_<SAMPLE_NAME>_,  _<FASTQ1_FILE_PATH>_,  _<FASTQ2_FILE_PATH>_,_<CELLRANGER_REPOSITORY_PATH>_
    -   (In case of 10X based scRNA-seq sample [--sequencing  **chromium**])  
        >_<SAMPLE_NAME>_,  _<FASTQ1_PATH>_,  _<FASTQ2_PATH>_,_<BAM_PATH>_, _<BAM_INDEX_PATH>_,_<DGE_PATH>_

2. Optional, **[\-\-barcodes]**: Provide within a **CSV** (ex: barcodes_whitelist.csv) for each input sample, a barcode whitelist file path:
   > _<SAMPLE_NAME>_,_<BARCODE_WHITELIST_FILE_PATH>_

3. (Optional, **[\-\-clusters]**: Provide within a **CSV** (ex: clusters.csv) for each input sample, a barcode whitelist file path:
    > _<SAMPLE_NAME>_,_<BARCODE_CLUSTERS_FILE_PATH>_

### - SCALPEL Execution
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

1. Configurate sample file:
```
> cat samplesheet.csv
SRR6129050,SRR6129050_S1_L001_R1_001.fastq.gz,SRR6129050_S1_L001_R2_001.fastq.gz,/data/fake_data/DATAS/GSE104556/SAMPLES/SRR6129050/
SRR6129051,SRR6129051_S1_L001_R1_001.fastq.gz,SRR6129051_S1_L001_R2_001.fastq.gz,/data/fake_data/DATAS/GSE104556/SAMPLES/SRR6129051/
```

2. (**Optional**), Configurate barcode whitelist file:
```
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

3. Running
```
> nextflow run -resume SCALPEL/main.nf \
    --sequencing chromium \
    --samplesheet samplesheet.csv \
    --transcriptome gencode.vM10.transcripts.fa \
    --gtf gencode.vM21.annotation.gtf \
    --ipdb mm10.polyA.track \
    --barcodes barcodes_whitelist.csv \ 
```

- **Results**

During its execution **SCALPEL** shows interactively in the Console all the information about the executed tasks:
```
===============================
SCALPEL - N F P I P E L I N E
===============================
Author: PLASS Lab ; Franz AKE
*****************
P-CMRC - Barcelona, SPAIN

input files:
- Annotation required files(required):
- transcriptome reference [--transcriptome]: gencode.vM10.transcripts.fa
- annotation GTF reference [--gtf]: gencode.vM10.annotation.gtf
- internal priming annotation [--ipdb]: mm10.polyA.track

- Reads processing files (required):
- samplesheet [--samplesheet]: samplesheet.txt

- Params:

Required:
- sequencing type (required): chromium

Optional:
- barcodes whitelist [--barcodes] (optional): barcodes_whitelist.csv
- cell clusters annotation [--clusters] (optional): null
- transcriptomic distance threshold [--dt_threshold] (optional, default 600bp): 600
- transcriptomic end distance threhsold [--de_threshold] (optional, default 30bp): 30
- minimal distance of internal priming sites (IP) from isoform 3'ends [--ip_threshold] (optional, 60nuc): 60
- gene fraction abundance threshold [--gene_fraction] (optional, default '98%'): 98%
- binsize threshold for transcriptomic distance based probability [--binsize] (optional, default '20): 20
- reads subsampling threshold [--subsample] (optional, default 1): 1

executor >  local (2600)
[9c/c2560c] process > annotation_preprocessing:salmon_transcriptome_indexing [100%] 1 of 1, cached: 1 ✔
[00/5af3e8] process > annotation_preprocessing:salmon_bulk_quantification (SRR6129050, SRR6129050_S1_L001_R1_001.fastq.gz, SRR6129050_S1_L001_R2_001.fastq.gz) [100%] 2 of 2, cached: 2 ✔
[2c/96a4af] process > annotation_preprocessing:tpm_counts_average ([SRR6129051.sf, SRR6129050.sf])  [100%] 1 of 1, cached: 1 ✔
[f7/9f45d1] process > annotation_preprocessing:isoform_selection_weighting (chr8, merge_quants.txt)  [100%] 22 of 22, cached: 22 ✔
[2c/35fe7c] process > reads_processing:samples_loading:read_10Xrepo (SRR6129051)  [100%] 2 of 2, cached: 2 ✔
[c5/86bdf8] process > reads_processing:samples_loading:bam_splitting (SRR6129051, chr10, /data/fake_data/NEW/SCALPEL/10X/SRR6129051_curatedBarcodes.txt)  [100%] 44 of 44, cached: 44 ✔
[fc/5d3dd0] process > reads_processing:bedfile_conversion (SRR6129051, chr11)  [100%] 44 of 44, cached: 44 ✔
[31/00f4b7] process > reads_processing:reads_mapping_and_filtering (SRR6129050, chr10, chr10.bed, chr10.exons)  [100%] 44 of 44, cached: 44 ✔
[da/cb39e5] process > reads_processing:ip_splitting (chr10)  [100%] 22 of 22, cached: 22 ✔
[c1/afa324] process > reads_processing:ip_filtering (SRR6129051, chr10, chr10_sp.ipdb)  [100%] 44 of 44, cached: 44 ✔
[c7/40f884] process > isoform_quantification:probability_distribution (SRR6129051)  [100%] 2 of 2, cached: 2 ✔
[57/7ae88f] process > isoform_quantification:fragment_probabilities (SRR6129051,chr4)  [100%] 44 of 44, cached: 44 ✔
[fb/adf53b] process > isoform_quantification:cells_splitting (SRR6129051)  [100%] 2 of 2, cached: 2 ✔
[69/ce54d8] process > isoform_quantification:em_algorithm (SRR6129051, TTTGTCATCTGTCCGT-1)  [100%] 2042 of 2042, cached: 2042 ✔
[6c/441287] process > isoform_quantification:cells_merging (SRR6129050)  [100%] 2 of 2, cached: 2 ✔
[c3/a2f9bf] process > isoform_quantification:dge_generation (SRR6129050, SRR6129050_isoforms_quantified.txt, SRR6129050.counts.txt)  [100%] 2 of 2, cached: 2 ✔
[03/00f705] process > apa_characterization:differential_isoform_usage ([SRR6129050_seurat.RDS, SRR6129051_seurat.RDS])  [100%] 1 of 1, cached: 1 ✔

**Completed at: 08-Jun-2024 10:33:52**
**Duration  : 39m 5s**
**CPU hours : 54.4 (41.9% cached)**
**Succeeded : 2'600**
**Cached  : 248**
```
Following the execution, a _**./results**_ folder is generated encompassing all the final and intermediated files generated by **SCALPEL**:

 - _**./results/final_results**_ contains
 -- the final differential isoform usage table between the input samples
 -- the Seurat object encompassing the iDGE of the input samples 
 -- the iDGE table of each input sample
```sh
DIU_table.csv  iDGE_seurat.RDS.  SRR6129050/SRR6129050_APA_DGE.txt   SRR6129051/SRR6129051_APA_DGE.txt
```

- Information about the collapsed isoforms by SCALPEL (see paper) can be found in _**./results/annotation_processing/isoform_processing**_
```sh
chr10_collapsed_isoforms.txt  chr13_collapsed_isoforms.txt  chr16_collapsed_isoforms.txt ...
```

## More about SCALPEL Usage

- Example of SCALPEL application on 10X scRNA-seq
- Example of SCALPEL application on DropSeq scRNA-seq



**Important: The chromosome names must to be consistent between the BAM file and the annotation files. If the BAM file contains the '_chr_' character, the GTF and FASTA annotation files, and the internal priming reference annotation should contains the '_chr_' character, and inversely !**

**The internal priming reference files provided contains by default the '_chr_' character in the chromosome names !!**

**Be careful to delete the _work_ directory containing nextflow temporary files** when scalpel runs all its processs sucessfully and you don’t plan to relaunch scalpel with modified parameters. (This folder can fill an high memory physical space depending of the size of input files analyzed)


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
