##### Create filtered Seurat objects and cell barcode files for pseudobulk analyses ####

# Load libraries
library(Seurat)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(R.utils)

# Receive command line arguments 
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
key <- as.character(argus[2])
mtx <- as.character(argus[3])
seed <- as.numeric(argus[4])

# Set seed
set.seed(seed = seed)


process_data <- function(sample_name,mtx){
  print(paste("Begin processing sample",sample_name,sep=" "))
  print(paste("Read in file",sample_name,sep=" "))


if(mtx=="Y"){
indv_key <- filter(key, sample_name == filename)
temp <- ReadMtx(
  mtx=paste0('../raw_data/',sample_name,'.matrix.mtx.gz'),
  cells=paste0('../raw_data/',sample_name,'.barcodes.tsv.gz'),
  features=paste0('../raw_data/',sample_name,'.features.tsv.gz')
)

} else if(mtx=="N"){
temp <- Read10X_h5(filename = paste0('../cell_ranger_output/',sample_name,
'/outs/filtered_feature_bc_matrix.h5'))
 } else{print('SAMPLE NAME ERROR; CHECK KEY')
stop()}
temp <- CreateSeuratObject(counts=temp,
                           project="MSI",
                           min.cells=3,
                           min.features=100)


  #add percent mitochondrial
  temp[["percent.mt"]] <- PercentageFeatureSet(temp, pattern = "^MT-")
  temp <- subset(temp,subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 35)
  temp <- NormalizeData(temp,normalization.method = "LogNormalize",scale.factor = 10000)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(temp)
  temp <- ScaleData(temp,feature=all.genes)
  temp <- RunPCA(temp, features=VariableFeatures(object=temp))
  
  
  print(paste("Identifying clusters"))
  
  temp <- FindNeighbors(temp,dims=1:15)
  temp <- FindClusters(temp,resolution = 0.5)
  temp <- RunUMAP(temp, dims = 1:15)
  #find marker genes that are differentially expressed in clusters
  
  
saveRDS(temp, paste0('../filtered_h5/',sample_name,'.rds'))
 
return(temp)
 
  
}

# write out processed h5 as a new file for local review

gen_barcodes <- function(sample){

# Supply an object returned by the process_data() function

barcodes <- names(sample$nCount_RNA)
clusters <- sample$seurat_clusters
num_clust <- length(levels(clusters)) -1
tempy <- data.frame(barcodes=barcodes,cluster=as.character(clusters))
fwrite(tempy,file=paste0('../pseudobulk_barcodes/',sample_name,'/',sample_name,'_all_cell_barcodes.tsv'),
sep='\t',col.names=FALSE,row.names=FALSE)

  for(z in 0:num_clust){
    temp <- as.matrix(t(clusters))
    barcodes2 <- names(temp[,temp[] == z])
    fwrite(data.table(barcodes2),file = paste0('../pseudobulk_barcodes/',sample_name,
	'/',sample_name,'_cluster_',z,'.tsv'),sep='\t',col.names = FALSE)
    }
}


s_obj <- process_data(sample_name,mtx)

gen_barcodes(s_obj)
