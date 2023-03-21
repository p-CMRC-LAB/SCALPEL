

SCALPEL for characterization of Alternative polyadenylation  at single-cell resolution
======================================================================================


<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/Dessin_scalpel.svg/1200px-Dessin_scalpel.svg.png" alt="SCALPEL" width="300" height="300">
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

[Conda](https://www.anaconda.com/) package tool or [Mamba](https://github.com/mamba-org/mamba) (Fast reiplementation of conda)

### Installation

1. Clone the repo and enter in the folder
   ```sh
   git clone https://github.com/p-CMRC-LAB/SCALPEL.git
   ```
2. Install the required packages using the requirement.txt file in the SCALPEL folder
   ```sh
   conda install -c bioconda -c defaults -c conda-forge --file SCALPEL/requirements.txt
   ```
   
   Another solution (if conda installation takes long) can be to create a Conda environment, install Mamba (faster implementation of Conda) and install the packages using mamba:
   ```sh
   conda create --name scalpel_env -c conda-forge -c bioconda mamba python=3.9
   conda activate scalpel_env
   mamba install -c bioconda -c defaults -c conda-forge --file SCALPEL/requirements.txt
   ```

### Usage


See [SCALPEL usage documentation](https://readthedoctest-franz.readthedocs.io/en/latest/scRNAseq_10X.html) for detailed usage of the tool



<!-- CONTACT -->
## Contact

Your Name - [@twitter](https://twitter.com/aerodx5) - fake@idibell.cat

Project Link: [SCALPEL Github](https://github.com/p-CMRC-LAB/SCALPEL)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

