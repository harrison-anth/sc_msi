#### Create plots for each sample prior to integration ####

# Load libraries
library(Seurat)
library(gridExtra)
library(patchwork)
library(tidyverse)
library(scAnnotatR)
library(data.table)
library(cowplot)
library(Nebulosa)
library(ggpubr)
library(R.utils)

# Receive arguments from command line
argus <- (commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1])
sample_name <- as.character(argus[1])

plotter <- function(sample_name){

s_obj <- readRDS(paste0('../annotated_h5/',sample_name,'.rds'))
atomic_cells <- readRDS(paste0('../atomic/',sample_name,'.rds'))

s_obj$scATOMIC_pred <- atomic_cells$scATOMIC_pred
s_obj$pan_cancer_cluster <- atomic_cells$pan_cancer_cluster
s_obj$sensor_rna_prob <- as.numeric(s_obj$sensor_rna_prob)

plots <- ggarrange( 
DimPlot(s_obj)+ggtitle(paste0('Clusters')),
DimPlot(s_obj,group.by='pan_cancer_cluster',cols=c('red','blue','grey'))+ggtitle('scATOMIC'),
DimPlot(s_obj,group.by='sensor_rna_status',cols=c('red','blue'))+ggtitle('MSIsensor-RNA status'),
FeaturePlot(s_obj,features=c('sensor_rna_prob'),cols=c('blue','red'))+ggtitle('MSIsensor-RNA score'),
ncol=2,nrow=2)

plots_final <- annotate_figure(plots,
top=text_grob(paste0('Sample: ', sample_name),
color='black',face='bold',
size=20))

ggsave(paste0('../reports/',sample_name,'_sample_report.pdf'),plots_final)

}


plotter(sample_name)
