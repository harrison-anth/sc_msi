# SC-MSI
## Computational pipeline to identify MSI-H cells and measure heterogeneity in the biomarker

### Information about the author(s) of this code:
Name(s): Harrison Anthony 
contact information: h dot anthony1 at universityofgalway dot ie
(Alternate contacts can be found on my github profile)

### License information
MIT; see LICENSE file for more information.

### Repository information

This repository functions as a distribution of the SC-MSI Snakemake pipeline used in our recept manuscript. The results and raw code used in the manuscript can be found in the legacy version of this
repository (https://github.com/harrison-anth/sc_msi_legacy).

This pipeline has been tested on Ubuntu 20.04 with Snakemake 8.27.1 and Conda 24.1.2.

### Planned pipeline improvements

* Create separate files for each rule to help users incorporate multi-threading

* Combine FASTQ/MTX file workflows with automatic file detection

* Include small reference transcriptome to verify installation

* Containerize workflow to optionally not install each software

### Installation instructions

Download and install the following dependencies: 

1.) Conda/Mamba (https://docs.conda.io/en/latest/)

2.) Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

3.) Cellranger (https://www.10xgenomics.com/support/software/cell-ranger/downloads)

4.) Reference transcriptome (https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)

5.) Create the conda environments stored in the conda_envs/ directory:

```conda env create -f atomic.yml```

```conda env create -f seurat.yml```

5.) Activate the atomic Conda environment and install scATOMIC

https://github.com/abelson-lab/scATOMIC

That's it!

**Note:** Snakemake is a great workflow manager because it is transferable across platforms, scalabale with HPC's, and reproducible. However, Snakemake 
(and other bioinformatics workflow tools) can become complicated very quickly when using many different wilcards 
and when being integrated into an HPC framework like SLURM. We have provided small FASTQ and MTX files that serve as examples on how to build the
manifest, config, and profile for SC-MSI. An experienced Snakemake user would be needed to include changes to the pipeline (new rules, file formats, etc.,) 
and could very well affect the functionality of the pipeline. 

### SC-MSI quick use guide

1.) Verify installation by running the test sample. 

``` snakemake -s handle_fastq.snake ```

2.) Edit the settings in the config file to replace the default test settings

FASTQ files -- handle_fastq.config

MTX files -- handle_mtx.config

Note: There is a separate README in the manifests/ directory that describes each config file and how to create them. 

3.) Run the pipeline

FASTQ files -- ``` snakemake -s handle_fastq.snake --cores 1 --use-conda ```

MTX files -- ``` snakemake -s handle_mtx.snake ```

These are the very basic possible commands with Snakemake. It is recommended to take time to create a custom Snakemake profile to 
store user settings that enable multi-threading/multi-core processing. There are also many RAM intense applications (Cellranger and InferCNV) 
that will require tweaking the memory settings for each rule to fit the system settings.
