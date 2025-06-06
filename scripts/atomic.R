#### Run scATOMIC ####

# Quickly adjust the python paths (as R magic has trouble with this
Sys.setenv(RETICULATE_PYTHON = Sys.which("python3"))


# Load libraries
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(ggplot2)
library(R.utils)
library(foreach)
library(parallel)
library(cutoff.scATOMIC)
library(scATOMIC)
library(devtools)


# Receive arguments from command line
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
cores <- as.numeric(argus[2])
seed <- as.numeric(argus[3])

# Set seed
set.seed(seed = seed)

# Read in sample
rds1 <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))

# Create sparse matrix
sparse_matrix <- rds1@assays$RNA$counts

# Get cell type predictions
cell_predictions <- run_scATOMIC(sparse_matrix,mc.cores = cores)
results_temp <- create_summary_matrix(prediction_list = cell_predictions, 
use_CNVs = F, 
modify_results = T, 
mc.cores = cores, 
raw_counts = sparse_matrix, 
min_prop = 0.5 )

# Attach results
rds2<- AddMetaData(rds1, results_temp)

# Write out new Seurat object
saveRDS(rds2,paste0('../atomic/',sample_name,'.rds'))

# Write out cancer only barcodes
barcodes <- data.frame(barcode=names(rds2$pan_cancer_cluster),call=rds2$pan_cancer_cluster)
fwrite(barcodes,sep='\t',paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_cancer_barcodes.tsv'),row.names=FALSE,col.names=FALSE)


