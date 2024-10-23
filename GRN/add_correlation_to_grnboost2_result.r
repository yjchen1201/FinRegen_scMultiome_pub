#!/usr/bin/Rscript
#setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step4_new_TF_target_correlation/")
# Rscript gene_gene_correlation.r celltype stage
library(dplyr)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
celltype = args[1]
stage = args[2]

network =read.table(paste0(celltype,"_",stage,"_grnboost2_network.tsv"),sep="\t")
Corr_res = readRDS(paste0(celltype,"_",stage,"_Corr_res.rds"))

colnames(network) <- c("TF","Gene","Score")
network$index = paste(network$TF,network$Gene,sep="-->")
k1 = which(Corr_res$Var1 %in% network$TF == T & Corr_res$Var2 %in% network$Gene == T)
Corr_res_cl = Corr_res[k1,]
Corr_res_index = paste(Corr_res_cl$Var1,Corr_res_cl$Var2,sep="-->")
m = match(network$index,Corr_res_index)
network$Corr = Corr_res_cl$value[m]
saveRDS(network,file=paste0(celltype,"_",stage, "_grnboost2_network_wCorr.rds"))

print(paste0(celltype,"_",stage, "_grnboost2_network_wCorr saved!"))
