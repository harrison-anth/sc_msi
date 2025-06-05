#### Annotate filtered h5 files ####

# Load libraries
library(Seurat)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(R.utils)
library(data.table)

# Receive arguments from command line
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
key <- fread(as.character(argus[2]))
msi_cutoff <- as.numeric(argus[3])

indv_key <- filter(key, filename == sample_name)

# Functions to add metadata to seurat object
add_sensor_rna <- function(s_obj,sample_name,threshold){
sensor_rna_results <- fread(paste0('../sensor_rna_results/',sample_name,'.txt'))
sensor_rna_results$msi_status[sensor_rna_results$`probability_of_MSI-H` >= msi_cutoff ] <- 'MSI-H'
sensor_rna_results$msi_status[sensor_rna_results$`probability_of_MSI-H` < msi_cutoff ] <- 'MSS'

temp <- t(sensor_rna_results)
colnames(temp) <- temp[1,]
temp <- temp[2:3,]
s_obj <- AddMetaData(s_obj,temp[1,],col.name='sensor_rna_status')
s_obj <- AddMetaData(s_obj,temp[2,],col.name='sensor_rna_prob')
return(s_obj)
}

new_sample <- readRDS(paste0('../filtered_h5/',sample_name,'.rds'))

new_sample <- add_sensor_rna(new_sample,sample_name)

# Add scAOMIC and other metadata

new_sample[['sample_name']] <- sample_name

new_sample[['tissue']] <- indv_key$site
atomic_cells <- readRDS(paste0('../atomic/',sample_name,'.rds'))

new_sample$scATOMIC_pred <- atomic_cells$scATOMIC_pred
new_sample$classification_confidence <- atomic_cells$classification_confidence
if('pan_cancer_cluster' %in% colnames(atomic_cells@meta.data)){
new_sample$pan_cancer_cluster <- atomic_cells$pan_cancer_cluster
}

saveRDS(new_sample, paste0('../annotated_h5/',sample_name,'.rds'))





