# SC-MSI
## Computational pipeline to assess intratumoral heterogeneity of microsatellite instability status in single-cell data

### Submitted to *Cancer Research*

### Information about the author(s) of this code:
Name(s): Harrison Anthony 
contact information: h dot anthony1 at universityofgalway dot ie
(Alternate contacts can be found on my github profile)

## Repository information

This repository functions as a distribution of the SC-MSI Snakemake pipeline used in our recept manuscript. The results and raw code used in the manuscript can be found in the legacy version of this
repository (https://github.com/harrison-anth/sc_msi_legacy)

### Before use

The SC-MSI pipeline is verified to work on Ubuntu 20.04 with Snakemake 8.27.1 and Conda 24.1.2. Other Snakemake versions will work with this pipeline, but the variables in the example
config file can change between versions. It is also possible to use Snakemake/Conda on other operating systems, but we have not (and do not plan to) initialize the pipline outside of Linux.

To install Snakemake follow this guide (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)




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
