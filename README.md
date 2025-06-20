# SC-MSI
## Computational pipeline to identify MSI-H cells and measure heterogeneity in the biomarker

### Information about the author(s) of this code:
Name(s): Harrison Anthony 
contact information: h dot anthony1 at universityofgalway dot ie
(Alternate contacts can be found on my github profile)

### License information
MIT; see LICENSE file for more information.

### Repository information

This repository functions as a distribution of the SC-MSI Snakemake pipeline used in our recept manuscript. 
The results and raw code submitted as part of the publication can be found in the legacy version of this
repository (https://github.com/harrison-anth/sc_msi_legacy).

This pipeline has been tested on Ubuntu 20.04 with Snakemake 8.27.1, Conda 24.1.2, and Cell Ranger 7.2.0.

### Planned pipeline features to be added

* Create separate files for each rule to help users incorporate multi-threading and custom filter/cluster settings

* Combine FASTQ/MTX file workflows with automatic file detection

* Include small reference transcriptome to verify installation of FASTQ Snake file

* Containerize workflow to optionally not install each software

* Include additional/supplementary plots that were shown in the manuscript

**Feel free to message me or open an issue with any other suggestions you might have**

### Installation instructions

1.) Download and install the following dependencies: 

* Conda/Mamba (https://docs.conda.io/en/latest/)

* Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

**Note: Cell Ranger and reference transcriptome are only necessary for 
the FASTQ pipeline**
 
* Cellranger ```wget -O cellranger-7.2.0.tar.gz 
"https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz?Expires=1749374999&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=Oi4tu~mQ8B0T3hOUjGgm7Xy8f-Z92tdauD~0D~9M3GbwQKJW7txWVFWr6aQpmfI7mrYTM8965onbhcKNWBI2S9Uu8OYl7c8kWMLCzP1rs8rAeitL~eB~2kqc29NZrT5pnDjIkPTzh1huBj53qXmo7oDopIxdl2~llGwZ3oVd8dJSQfZtPe08NwPnOQgLa0IhD8TplsSf7uHGvyfR9lSBvvICSlkqzkuvKfE56Qg-fV5lRad7jLOjFTZNhPxZ4Mh4FpPuGxr-s2KWBa~OS0~4w3nP3PLUUVl5EjjHfW0XPusv0m-vT1jGskKWjN9Mzx0ttv8ddjNfiU04imz7o8BBfQ__"```

* Reference transcriptome (https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)

* Change the reference transcriptome path in the config file

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

**Run MTX Snake file with SLURM executor profile using 30 cores 
(parallel processing); use Conda environments**

 ``` snakemake -s handle_mtx.snake -p ../conda_envs/slurm_executor/ --cores 30 --use-conda --```

These are the very basic possible commands with Snakemake. It is recommended to take time to create a custom Snakemake profile to 
store user settings that enable multi-threading/multi-core processing. There are also many RAM intense applications (Cellranger and InferCNV) 
that will require tweaking the memory settings for each rule to fit the system settings.

6.) Go through the output files

* In the reports/ directory there will be patient and sample reports based on the information provided in the manifest files 

* There will be plots for each individual and each sample in the infer_cnv_results/ directory (see https://github.com/broadinstitute/infercnv/wiki/Running-InferCNV for information on how to interpret them)

* All individuals will also have summary stats written to a tsv file in the summary_stats/ directory

### Troubleshooting tips for Snakemake

1.) Re-run the pipeline with the --keep-going parameter to see which samples are not running correctly.

2.) Carefully go through each Snakemake log file to determine which step is throwing an ERROR

3.) If you are running with intense settings (multi-threaded/HPC) consider starting with a manifest file containing only 1 individual and 1 sample to check if 
your specifications are too high or insufficient. 

"Progress through trial and error depends not only on making trials, but on recognizing errors." ~ Virginia Postrel
