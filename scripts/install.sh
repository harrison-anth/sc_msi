#!/bin/bash

# Script to install dependencies for SINGLE-MSI

## This script is designed to be run from the scripts/ directory

## Download and install Conda and Snakemake
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

### Note: Must initalize new conda installation before running conda commands
### This can be accomplished by sourcing the user bashrc
source ~/.bashrc

### Use conda to create Snakemake environment
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake

## Create conda environments
conda env create -f ../conda_envs/atomic.yml
conda env create -f ../conda_envs/seurat.yml

## Install scATOMIC
conda activate atomic

### NB: While I have tried lifting over install instructions from the scATOMIC GitHub repository, this part can require some troubleshooting as there are many dependencies for scATOMIC
### If prompted, DO NOT update any of the R packages installed alongside scATOMIC as this can lead to version issues and will overwrite those in the atomic conda environment
Rscript -e 'devtools::install_version("dlm", version = "1.1.5", repos = "http://cran.us.r-project.org")
devtools::install_version("Rmagic", version = "2.0.3", repos = "http://cran.us.r-project.org")
if(!require(devtools)) install.packages("devtools")
if(!require(cutoff.scATOMIC)) devtools::install_github("inofechm/cutoff.scATOMIC", force = T)
if(!require(scATOMIC)) devtools::install_github("abelson-lab/scATOMIC")'


