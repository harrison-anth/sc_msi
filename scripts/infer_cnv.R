#### Run InferCNV to identify subclones ####

# Load libraries
library(data.table)
library(Matrix)
library(Seurat)
library(infercnv)
library(tidyverse)
library(R.utils)
library(stringr)

# Receive arguments from command line
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])
outdir <- as.character(argus[2])
seed <- as.numeric(argus[3])
cores <- as.numeric(argus[4])

# Set seed
set.seed(seed = 152727)

s_obj <- readRDS(paste0('../annotated_h5/',sample_name,'.rds'))

# Create annotations file
anno <- data.frame(cell_name=colnames(s_obj),cancer=s_obj$pan_cancer_cluster,
type=s_obj$scATOMIC_pred,msi=s_obj$sensor_rna_status)


#verify there are cancer cells in this sample
if(nrow(filter(anno, cancer == "Cancer"))){

print(paste0('ERROR: No cancer cells identified in ', sample_name, '. Re-assess scATOMIC settings or remove sample from manifest files.'))
stop()
}



# Verify there are MSI-H and MSS cells
if( nrow(filter(anno, msi == "MSI-H" & cancer == "Cancer")) <2 | nrow(filter(anno, msi == "MSS" & cancer == "Cancer")) < 2){
print('Too few MSI-H or MSS classifications; continuing infercnv without MSI labels')
final_anno <- data.frame(cell_name=anno$cell_name,type=anno$cancer)
fwrite(final_anno,file=paste0('../temp/',sample_name,'_anno.tsv'),sep='\t',col.names=FALSE)

} else{

anno$cancer2 <- paste0(anno$cancer,'_',anno$msi)
final_anno <- data.frame(cell_name=anno$cell_name,type=ifelse(anno$cancer =="Normal",yes=anno$cancer,no=anno$cancer2))

#filter NA's 
final_anno <- filter(final_anno, type != "Cancer_NA")

fwrite(final_anno,file=paste0('../temp/',sample_name,'_anno.tsv'),sep='\t',col.names=FALSE)
}

#create infercnv object
cnv_obj <- CreateInfercnvObject(raw_counts_matrix=s_obj@assays$RNA$counts,
annotations_file=paste0('../temp/',sample_name,'_anno.tsv'),
delim='\t',gene_order_file='../references/hg38_gencode_v27.txt',
ref_group_names=c("Normal"))

#fix scipen for infercnv run
options(scipen = 100)

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(cnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=paste0(outdir),
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             num_threads=cores
                             )

