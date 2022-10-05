|PyPI| |Bioconda| |Downloads| |CI| |Notebooks| |Docs| |Codecov|

SCALPEL for characterization of Alternative polyadenylation  at single-cell resolution
======================================================================================


**SCALPEL** is a toolkit for the characterization of Alternative ployadenylation (APA) and the Isoform quantification at single-cell resolution.

Manuscript
^^^^^^^^^^



Installation
^^^^^^^^^^^^
Install SCALPEL by running::
	
	# Clone SCALPEL into the workspace
	git clone https://github.com/p-CMRC-LAB/SCALPEL.git

	# Create SCALPEL environment and install packages required
	conda create -n scalpel_env -c conda-forge -c bioconda python=3.9 mamba nextflow
	conda activate scalpel_env
	mamba install --file SCALPEL/requirements.txt


Running
^^^^^^^
Exemple of launching command::

DROP_SEQ::

	nextflow run -resume SCALPEL/main.nf --sequencing 'DROP_SEQ' --folder_in {} --annot {} --ipdb {}

10X::

	nextflow run -resume SCALPEL/main.nf --sequencing '10X' --folder_in {} --annot {} --ipdb {}
	


The parameter '--folder_in' point the path of a folder containing:

-required

BAM file (final.bam) / index BAM file (final.bai) / DGE (Digital Gene Expression) file (DGE.txt)

-optional

Barcodes file (bc_list.txt)

The parameter '--annot' point the genomic annotation file (.GTF) and the parameter '--annot' point the internal priming reference file (.txt)

