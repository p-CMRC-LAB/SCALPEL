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

### Usage

#### ex: SCALPEL execution on 10X based scRNA-seq dataset ([Lukassen et al](https://pubmed.ncbi.nlm.nih.gov/30204153/))

For the need of the analysis in this vignette, the data used is a 10X dataset from the study from [Winterpacht A, Lukassen](https://pubmed.ncbi.nlm.nih.gov/30204153/) on Mouse. The demo data is 10X processed folder containing a BAM file with aligned reads from the published data [GSE104556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556) 
from GEO.

- **Input files**

SCALPEL can be applied on a classic 10X result folder, containing the default files (_possorted_genome_bam.bam, possorted_genome_bam.bam.bai, filtered_feature_bc_matrix.h5, barcodes.tsv.gz, ..._).

For running, SCALPEL requires to provide specific input files path & parameters:
 - _SAMPLE path_ **[\-\-samplesheet]**
 - _FASTA transcriptome reference path_ **[--transcriptome]**
 - _GTF annotation file path_ **[\-\-gtf]**
 - _Internal priming annotation file_ **[\-\-ipdb]**
 - _Sequencing type (chromium or dropseq)_ **[\-\-sequencing]**

Provide within the **_samplesheet.csv_** file in SCALPEL repository, the following input file paths : \
_<SAMPLE_NAME>_, _<FASTQ1_PATH>_, _<FASTQ2_PATH>_,_<CELLRANGER_REPOSITORY_PATH>_

For this analysis on 10X scRNA-seq, the related sample files can be downloaded here:
 - [FASTQ_FILES]()
 - [CELLRANGER_REPOSITORY](https://drive.bio.idibell.cat/index.php/s/eoYbCKA48eZXMDK)

Following the organism for the study, the reference genome annotation files (GTF) and reference transcript sequence (FASTA) can be downloaded on [GENCODE](https://www.gencodegenes.org/) repository. 

For this analysis on Mouse, we used the release M10 on GENCODE:
 - [(Mouse) GTF_annotation_file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz) \
 - [(Mouse) FASTA_reference_file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.transcripts.fa.gz) \

##### Internal priming files which reference all the internal priming positions
- [(Human) Internal priming annotation - GRCh38](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/GRCh38_2020_A_polyA.track.tar.gz) \
- [(Mouse) Internal priming annotation - mm10](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/mm10_polya.track.tar.gz) \




Differential analysis usage is done in the downstream analysis step following the isoform quantification by default considering the different samples provided. \
_To perform a differential analysis on defined cell clusters, provide an annotation file using the optional argument_: **[\-\-clusters]** \
_To perform SCALPEL on a specific set of barcoded cells, provide a cell whitelist using the optional argument_: **[\-\-barcodes]**

#### Execution
You can print the Scalpel help documentation by running the following command

- **Help**

```sh
> nextflow run -resume SCALPEL/scalpel.nf --help
```

    SCALPEL - N F   P I P E L I N E
    ===============================
    Author: PLASS Lab ; Franz AKE
    *****************
    P-CMRC - Barcelona, SPAIN

    input files:
    - Annotation required files(required):
        - transcriptome reference [--transcriptome]: ${params.transcriptome}
        - annotation GTF reference [--gtf]: ${params.gtf}
        - internal priming annotation [--ipdb]: ${params.ipdb}


    - Reads processing files (required):
        - samplesheet [--samplesheet]: ${params.samplesheet}

    - Params:
        Required:
        - sequencing type (required): ${params.sequencing}

        Optional:
        - barcodes whitelist [--barcodes] (optional): ${params.barcodes}
        - cell clusters annotation [--clusters] (optional): ${params.clusters}
        - transcriptomic distance threshold [--dt_threshold] (optional, default 600bp): ${params.dt_threshold}
        - transcriptomic end distance threhsold [--de_threshold] (optional, default 30bp): ${params.de_threshold}
        - minimal distance of internal priming sites (IP) from isoform 3'ends [--ip_threshold] (optional, 60nuc): ${params.ip_threshold}

 - **Configuration of the input samplesheet file**

Provide within _samplesheet.csv_ the sample paths:
```
SRR6129050,/home/fake/CEPH/DATAS/GSE104556/FASTQs/SRR6129050_S1_L001_R1_001.fastq.gz,/home/fake/CEPH/DATAS/GSE104556/FASTQs/SRR6129050_S1_L001_R2_001.fastq.gz,/home/fake/CEPH/DATAS/GSE104556/SAMPLES/SRR6129050
SRR6129051,/home/fake/CEPH/DATAS/GSE104556/FASTQs/SRR6129051_S1_L001_R1_001.fastq.gz,/home/fake/CEPH/DATAS/GSE104556/FASTQs/SRR6129051_S1_L001_R2_001.fastq.gz,/home/fake/CEPH/DATAS/GSE104556/SAMPLES/SRR6129051
```

- **Configuration of parameters for the Nextflow execution**

Within the SCALPEL/nextflow.config, you can modify the existing parameters to suit the desired peformances for SCALPEL execution
```
/* Define Nextflow configuration settings for SCALPEL pipeline execution */
/* =====================================================================

Here, Please define the desired settings for SCALPEL execution !!

/* -> Processes
================
*/

process {
  cpus = 30
  executor = 'slurm'
  withLabel: big_mem {
    memory = 90.GB
  }
  withLabel: small_mem {
    memory = 30.GB
  }
}
```

- **SCALPEL execution**

Execute SCALPEL within the activated CONDA environment:
```
nextflow run -resume SCALPEL/main.nf -resume --samplesheet samplesheet.csv --sequencing chromium --transcriptome /home/fake/gencode.vM10.transcripts.fa --gtf /home/fake/gencode.vM10.annotation.gtf --ipdb /home/fake/mm10_polya.track
```

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
[Access SCALPEL_PAPER](https://www.biorxiv.org/content/10.1101/2024.06.21.600022v1)

Quantification of transcript isoforms at the single-cell level using SCALPEL \
Franz Ake, Sandra M. Fernández-Moya, Marcel Schilling, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass \
bioRxiv 2024.06.21.600022; doi: https://doi.org/10.1101/2024.06.21.600022
