# SC-MSI
## Computational pipeline to identify MSI-H cells and measure heterogeneity in the biomarker

### Information about the author(s) of this code:
Name(s): Harrison Anthony 
contact information: h dot anthony1 at universityofgalway dot ie
(Alternate contacts can be found on my github profile)

### License information
MIT; see LICENSE file for more information.

## Repository information

This repository functions as a distribution of the SC-MSI Snakemake pipeline used in our recept manuscript. The results and raw code used in the manuscript can be found in the legacy version of this
repository (https://github.com/harrison-anth/sc_msi_legacy)

### Before use

This pipeline has been tested on Ubuntu 20.04 with Snakemake 8.27.1 and Conda 24.1.2. 

Dependencies: 
Conda/Mamba (https://docs.conda.io/en/latest/)

Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

pandas (https://pandas.pydata.org/)

Cellranger (https://www.10xgenomics.com/support/software/cell-ranger/downloads)

Reference transcriptome

That's it!

Note: Snakemake is a great workflow manager because it is transferable across platforms, scalabale with HPC's, and reproducible. However, Snakemake 
(and other bioinformatics workflow tools) can become complicated very quickly when using many different wilcards 
and when being integrated into an HPC framework like SLURM. We have provided small BAM	and MTX files that serve as examples on how to build the
manifest, config, and profile for SC-MSI. An experienced Snakemake user would be needed to include changes to the pipeline (new rules, file formats, etc.,) 
and could very well affect the functionality of the pipeline. 

### SC-MSI quick use guide

If you have FASTQ files:

1.) Supply a manifest file that specifies the location of your FASTQ files and a key that has information on which samples should be integrated together (see the manifest directory for examples).

2.)  







### Conda environments
The conda_envs folder contains all conda environments used for data handling and all MSI tools.

sensor.yml was used to run MSIsensor MSIsensor2 and MSIsensor-pro

msings.yml was used to run mSINGS

MSIR.yml was used to run R code

mantis.yml was used to run MANTIS

vcf2maf.yml was used to convert from vcf file format to maf file format (note even with this environment, vcf2maf is notoriously finicky and will require 
a working local installation of VEP)

### Publication results
The TCGA results and scripts required to generate them are included in the tcga/ folder alongside a separate README.md
 detailing the pipeline. 

The results of all non-TCGA datasets are in the non_tcga/ folder. Each separate dataset has its own subfolder and README.md. 

### Graphs

All the code necessary to generate the graphs used in the paper are stored in the manuscript_graphs/ directory. The Rmarkdown files
used cater to my local paths and settings but can be adapated by any experienced R user. 
All graphs were generated using R version 4.1.2

### Baselines and reference files

The reference files are available from the GDC (https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)

The baselines created for each tool and used in this study are included in the baselines/ directory. 

### Test example

To ensure reproducibility of the results and quick implementation of our pipelines for other researchers who want to use these tools, we have included
a test example that allows for the user to quickly use each MSI tool on small BAM and gene count matrices that have been subset down to only the microsatellites found on 
chromosome 7. See the test/ directory for a separate README.md related to their use. 



Please feel free to reach out with any questions if this README has not answered your questions. 





# single_msi
Single cell rna-seq MSI tests to explore heterogeneity in MSI status as biomarker.

#data sources

6 scRNA samples from PRJNA932556

4 scRNA samples from PRJNA796219

8 scRNA samples from PRJNA650549

41 scRNA samples from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163558

40 samples from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205506
