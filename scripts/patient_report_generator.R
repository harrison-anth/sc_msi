#### Create nice markdown patient reports ####

# Load libraries
library(R.utils)
library(rmarkdown)


# Receive arguments from the command line
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])


# Call Rmarkdown

rmarkdown::render('parallelized_reporter.rmd',output_file=paste0('../reports/',sample_name,'.html'),params=list(sample_id=sample_name))
