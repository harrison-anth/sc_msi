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

This pipeline has been tested on Ubuntu 20.04 with Snakemake 8.27.1, Conda 24.1.2, and Cell Ranger 7.2.0.

### Planned pipeline improvements

* Create separate files for each rule to help users incorporate multi-threading

* Combine FASTQ/MTX file workflows with automatic file detection

* Include small reference transcriptome to verify installation of FASTQ Snake file

* Containerize workflow to optionally not install each software

### Installation instructions

1.) Download and install the following dependencies: 

* Conda/Mamba (https://docs.conda.io/en/latest/)

* Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

* Cellranger (https://www.10xgenomics.com/support/software/cell-ranger/downloads)

* Reference transcriptome (https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)

2.) Create the conda environments stored in the conda_envs/ directory:

```conda env create -f atomic.yml```

```conda env create -f seurat.yml```

3.) Activate the atomic Conda environment and install scATOMIC

https://github.com/abelson-lab/scATOMIC

That's it!

**Note:** Snakemake is a great workflow manager because it is transferable across platforms, scalabale with HPC's, and reproducible. However, Snakemake 
(and other bioinformatics workflow tools) can become complicated very quickly when using many different wilcards 
and when being integrated into an HPC framework like SLURM. We have provided small FASTQ and MTX files that serve as examples on how to build the
manifest, config, and profile for SC-MSI. An experienced Snakemake user would be needed to include changes to the pipeline (new rules, file formats, etc.,) 
and could very well affect the functionality of the pipeline. 

### SC-MSI quick use guide

1.) Verify installation by running the test sample. 

``` snakemake -s handle_fastq.snake --cores 1 --use-conda ```

``` snakemake -s handle_mtx.snake --cores 1 --use-conda ```

2.) Move your FASTQ or MTX files into the raw_data/ directory

3.) Create manifest files that link individual ID's to sample ID's (see manifests/ directory for details)

4.) Edit the settings in the config file to replace the default test settings

FASTQ files -- handle_fastq.config

MTX files -- handle_mtx.config

5.) Re-run the pipeline with updated config files and additional settings (some examples below but reference the above link to Snakemake for more information)

**Run FASTQ Snake file with 1 core (single-threaded);
use Conda environments; 
don't stop if error encountered; 
increase time to detect output file; 
re-run previously incomplete rules/files triggered by modification time only**

``` snakemake -s handle_fastq.snake --cores 1 --use-conda --keep-going --latency-wait 120 --rerun-incomplete --rerun-triggers mtime```

** Run FASTQ Snake file with SLURM executor profile using 30 cores (parallel processing); use Conda environments**

 ``` snakemake -s handle_mtx.snake -p ../conda_envs/slurm_executor/ --cores 30 --use-conda --```

These are the very basic possible commands with Snakemake. It is recommended to take time to create a custom Snakemake profile to 
store user settings that enable multi-threading/multi-core processing. There are also many RAM intense applications (Cellranger and InferCNV) 
that will require tweaking the memory settings for each rule to fit the system settings.
