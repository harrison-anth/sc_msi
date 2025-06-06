# Information on required files and how they are structured

Note: downsampled fastq and matrix test files are included to verify the correct installation of the pipeline and provide an example. 

# Additional information: 

## Patient ID file
TXT file with one patient ID per line.

## Sample names file
TXT file with one sample name per line.

## Key file
TSV file that links patient ID's with the associated files for that individual. 

ID | site | filename | msi_status | msi_test

Where ID is a column of patient ID's, site specifies location of sample (tumor1, lymph, normal, etc.), filename is the pathable name of the sample in the fastq 
or matrix directory, msi_status is the PCR/IHC test result for the sample (can be left as unknown or NA), and msi_test is the type of test used (e.g. PCR or IHC).

 
