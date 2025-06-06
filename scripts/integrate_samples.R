#### Integrate samples by matching patient ID in key #### from common patient ID's

# Load libraries
library(scATOMIC)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(cutoff.scATOMIC)
library(copykat)
library(R.utils)
library(scATOMIC)
library(tidyverse)

# Receive command line arguments
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
key <- fread(as.character(argus[2]))
seed <- as.numeric(argus[3])

# Set seed
set.seed(seed = seed)

# Define primary function

integrate_data <- function(sample_name){

indv_key <- filter(key, patient_id == sample_name)
  print(paste("Begin processing sample",sample_name,sep=" "))
  print(paste("Integrating", nrow(indv_key),"samples", sep=" "))

#check to make sure integration is needed
num_tumor_samps <- nrow(filter(indv_key, site != "normal" & site != "Normal"))
num_tot_samps <- nrow(indv_key)

if(num_tot_samps <2 | num_tumor_samps < 2){
print('Fewer than 2 tumor samples, no need to integrate.')

indv_key <- filter(indv_key, site != "normal" & site != "Normal")

integrated <- readRDS(paste0('../annotated_h5/',sample_name,'.rds'))

saveRDS(integrated, paste0('../integrated_samples/',indv_key$filename[1],'.rds'))

} else{
all_objs <- list()

for(i in 1:nrow(indv_key)) {

assign(x = paste0('s_obj'),value = readRDS(paste0('../annotated_h5/',indv_key$filename[i],'.rds')))
s_obj$sensor_rna_prob <- as.numeric(s_obj$sensor_rna_prob)

assign(x = paste0(sample_name,'_', i),value = s_obj) %>% append(all_objs) -> all_objs
}

merged <- merge(all_objs[[1]],
	y=c(all_objs[2:length(all_objs)]),
	 project = sample_name)
pre <- merged %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA(.) %>%
FindNeighbors(.,dims=1:10, reduction = "pca") %>% FindClusters(., resolution = c(0.5)) %>%
RunUMAP(.,dims=1:10, reduction = "pca")

integrated <- IntegrateLayers(object = pre, method = CCAIntegration, orig.reduction = "pca", 
new.reduction="integrated.cca", verbose = TRUE,k.weight=50)

integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

integrated <- integrated %>% FindNeighbors(., reduction = "integrated.cca", dims = 1:10) %>% 
FindClusters(.,res=c(0.5), cluster.name = "integrated.clusters") %>% 
RunUMAP(., dims=1:10, reduction= "integrated.cca", reduction.name = "integrated.umap")




saveRDS(integrated, paste0('../integrated_samples/',sample_name,'.rds'))

}

# Prepare MSIsensor-RNA matrix for integrated level.


clustys <- AggregateExpression(integrated,return.seurat=TRUE,group.by=c('seurat_clusters'))
ct_mat <- t(as.matrix(clustys@assays$RNA$counts))

ct_mat <- as.data.frame(ct_mat) %>% rownames_to_column('SampleID')

# Create custom training data columns for pseudobulks

pre_model <- fread ('../references/sensor_rna.csv')

common_names <- intersect(colnames(ct_mat),colnames(pre_model))

new_model <- pre_model %>% select(all_of(c('SampleID','msi',common_names)))

fwrite(new_model,paste0('../temp/',sample_name,'_training_model.csv'),sep=',')

new_data <- ct_mat %>% select(all_of(c('SampleID',common_names)))

new_data2 <- new_data[ rowSums(new_data[,2:ncol(new_data)]) >0, ]

fwrite(new_data2,paste0('../temp/',sample_name,'.csv'),sep=',')

# Now just the cancer
cancer_int <- subset(integrated,subset=pan_cancer_cluster == "Cancer")
#recluster
DefaultAssay(cancer_int) <- "RNA"
cancer_fin <- cancer_int %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 2000) %>%
#ScaleData(., vars.to.regress = c("percent.mt")) %>%
ScaleData(.) %>%
RunPCA(.,npcs=10) %>% RunUMAP(., dims = 1:10,n.neighbors=10) %>% FindNeighbors(., dims = 1:10) %>%
FindClusters(., resolution = 0.8)

saveRDS(cancer_fin, paste0('../integrated_samples/',sample_name,'_cancer.rds'))

if(length(levels(droplevels(cancer_fin$seurat_clusters)))<2){
print('Only one cluster of cancer cells; not aggregating expression')
clustys <- cancer_fin
} else{
clustys <- AggregateExpression(cancer_fin,return.seurat=TRUE,group.by=c('seurat_clusters'))
}



# Create pseudobulk clusters for MSIsensor-RNA

ct_mat <- t(as.matrix(clustys@assays$RNA$counts))

ct_mat <- as.data.frame(ct_mat) %>% rownames_to_column('SampleID')

#create custom training data columns for pseudobulks

pre_model <- fread ('../references/sensor_rna.csv')

common_names <- intersect(colnames(ct_mat),colnames(pre_model))

new_model <- pre_model %>% select(all_of(c('SampleID','msi',common_names)))

fwrite(new_model,paste0('../temp/',sample_name,'_cancer_training_model.csv'),sep=',')

new_data <- ct_mat %>% select(all_of(c('SampleID',common_names)))

new_data2 <- new_data[ rowSums(new_data[,2:ncol(new_data)]) >0, ]

fwrite(new_data2,paste0('../temp/',sample_name,'_cancer.csv'),sep=',')

}


integrate_data(sample_name)

