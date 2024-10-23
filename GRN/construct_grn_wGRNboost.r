#!/usr/bin/Rscript
# Rscript construct_grn_wGRNboost.r celltype_name
args = commandArgs(trailingOnly=TRUE)
celltype = args[1]

# if no spearman correlation added to grnboost2 result, run script commented below
# setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step4_new_TF_target_correlation/")
# network =read.table(paste0(celltype,"_grnboost2_network.tsv"),sep="\t")
# load(paste0(celltype,"_Corr_res"))

# colnames(network) <- c("TF","Gene","Score")
# network$index = paste(network$TF,network$Gene,sep="-->")
# k1 = which(Corr_res$Var1 %in% network$TF == T & Corr_res$Var2 %in% network$Gene == T)
# Corr_res_cl = Corr_res[k1,]
# Corr_res_index = paste(Corr_res_cl$Var1,Corr_res_cl$Var2,sep="-->")
# m = match(network$index,Corr_res_index)
# network$Corr = Corr_res_cl$value[m]
# saveRDS(network,file=paste0(celltype,"_grnboost2_network_wCorr.rds"))


#### load footprint, motif activity, peak-gene, and gene-gene correlations
#load("_TF_activity") ### from Step1 ####
setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step1_motif_TF_correlation/")
TF_activity <-readRDS(paste0(celltype,"_Corr_TFgene_Motif_Res.rds"))

#load("_PtoG_Ann") ### from Step2: Signac link peaks + peak annotation####
setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step2_p2g_links")
PtoG_Ann <-readRDS(paste0(celltype,"_PtoG_Ann.rds"))

# load("Globle_motif_cleaned") ### from Step3 ####
setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step3_new_Predicting_TF_Binding_Sites")
Globle_motif_cleaned_0dpa = readRDS(paste0(celltype,'_Globle_motif_cleaned_0dpa.rds'))
Globle_motif_cleaned_1dpa = readRDS(paste0(celltype,'_Globle_motif_cleaned_1dpa.rds'))
Globle_motif_cleaned_2dpa = readRDS(paste0(celltype,'_Globle_motif_cleaned_2dpa.rds'))
Globle_motif_cleaned_4dpa = readRDS(paste0(celltype,'_Globle_motif_cleaned_4dpa.rds'))
Globle_motif_cleaned_6dpa = readRDS(paste0(celltype,'_Globle_motif_cleaned_6dpa.rds'))

#load("_network") ### from Step4_new ####
setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step4_new_TF_target_correlation/")
network = readRDS(paste0(celltype,"_grnboost2_network_wCorr.rds"))

# start to construct grn table
setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step5_new_construct_TF_peak_target_links")
source("step5_construct_TF_peak_target_links.r")

GRNs_0dpa = Combined_all_things(GG_network=network,PG_cor=PtoG_Ann,TP_cor = TF_activity,Foot=Globle_motif_cleaned_0dpa)
saveRDS(GRNs_0dpa,file=paste0(celltype,"_GRNs_0dpa.rds"))
write.csv(GRNs_0dpa,file=paste0(celltype,"_GRNs_0dpa.csv"),row.names=F)

GRNs_1dpa = Combined_all_things(GG_network=network,PG_cor=PtoG_Ann,TP_cor = TF_activity,Foot=Globle_motif_cleaned_1dpa)
saveRDS(GRNs_1dpa,file=paste0(celltype,"_GRNs_1dpa.rds"))
write.csv(GRNs_1dpa,file=paste0(celltype,"_GRNs_1dpa.csv"),row.names=F)

GRNs_2dpa = Combined_all_things(GG_network=network,PG_cor=PtoG_Ann,TP_cor = TF_activity,Foot=Globle_motif_cleaned_2dpa)
saveRDS(GRNs_2dpa,file=paste0(celltype,"_GRNs_2dpa.rds"))
write.csv(GRNs_2dpa,file=paste0(celltype,"_GRNs_2dpa.csv"),row.names=F)

GRNs_4dpa = Combined_all_things(GG_network=network,PG_cor=PtoG_Ann,TP_cor = TF_activity,Foot=Globle_motif_cleaned_4dpa)
saveRDS(GRNs_4dpa,file=paste0(celltype,"_GRNs_4dpa.rds"))
write.csv(GRNs_4dpa,file=paste0(celltype,"_GRNs_4dpa.csv"),row.names=F)

GRNs_6dpa = Combined_all_things(GG_network=network,PG_cor=PtoG_Ann,TP_cor = TF_activity,Foot=Globle_motif_cleaned_6dpa)
saveRDS(GRNs_6dpa,file=paste0(celltype,"_GRNs_6dpa.rds"))
write.csv(GRNs_6dpa,file=paste0(celltype,"_GRNs_6dpa.csv"),row.names=F)


# STEP6: Identification of enriched gene regulatory sub-networks
## Filter by stage DEGs/DARs

setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step5_new_construct_TF_peak_target_links")
#### load the stage GRNs (total GRNs)
# GRNs_0dpa =readRDS(paste0(celltype,"_GRNs_0dpa.rds"))
# GRNs_1dpa =readRDS(paste0(celltype,"_GRNs_1dpa.rds"))
# GRNs_2dpa =readRDS(paste0(celltype,"_GRNs_2dpa.rds"))
# GRNs_4dpa =readRDS(paste0(celltype,"_GRNs_4dpa.rds"))
# GRNs_6dpa =readRDS(paste0(celltype,"_GRNs_6dpa.rds"))


# genes need to be expressed in this cell type and in certain stage
rna.exp = read.csv("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/0_integrated_wnn_merged_celltype_rna_expression_bycelltype_bystage.csv",row.names=1)
colnames(rna.exp) = gsub("mfin","0dpa",colnames(rna.exp))

Genes_need_0dpa = rownames(rna.exp)[which(rna.exp[,paste0(celltype,"_0dpa")]>0)]
Genes_need_1dpa = rownames(rna.exp)[which(rna.exp[,paste0(celltype,"_1dpa")]>0)]
Genes_need_2dpa = rownames(rna.exp)[which(rna.exp[,paste0(celltype,"_2dpa")]>0)]
Genes_need_4dpa = rownames(rna.exp)[which(rna.exp[,paste0(celltype,"_4dpa")]>0)]
Genes_need_6dpa = rownames(rna.exp)[which(rna.exp[,paste0(celltype,"_6dpa")]>0)]

#### Next convert to triple pairs: TF peak target ######
#### only kept the genes expressed in MG group cells ####
setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step6_new_identify_enriched_GRN_sub-networks")
source("step6_new_filter_subGRNs.r")
GRNs_triple_0dpa = Convert_to_triple_and_filterbygenes(GRNs_0dpa, Genes_need_0dpa)
GRNs_triple_1dpa = Convert_to_triple_and_filterbygenes(GRNs_1dpa, Genes_need_1dpa)
GRNs_triple_2dpa = Convert_to_triple_and_filterbygenes(GRNs_2dpa, Genes_need_2dpa)
GRNs_triple_4dpa = Convert_to_triple_and_filterbygenes(GRNs_4dpa, Genes_need_4dpa)
GRNs_triple_6dpa = Convert_to_triple_and_filterbygenes(GRNs_6dpa, Genes_need_6dpa)
write.csv(GRNs_triple_0dpa,file=paste0(celltype,"_GRNs_triple_0dpa.csv"),row.names=F)
write.csv(GRNs_triple_1dpa,file=paste0(celltype,"_GRNs_triple_1dpa.csv"),row.names=F)
write.csv(GRNs_triple_2dpa,file=paste0(celltype,"_GRNs_triple_2dpa.csv"),row.names=F)
write.csv(GRNs_triple_4dpa,file=paste0(celltype,"_GRNs_triple_4dpa.csv"),row.names=F)
write.csv(GRNs_triple_6dpa,file=paste0(celltype,"_GRNs_triple_6dpa.csv"),row.names=F)

#### Next convert to double pairs ######
#### only kept the genes expressed in MG group cells ####
GRNs_double_0dpa = Convert_to_double_and_filterbygenes(GRNs_0dpa, Genes_need_0dpa)
GRNs_double_1dpa = Convert_to_double_and_filterbygenes(GRNs_1dpa, Genes_need_1dpa)
GRNs_double_2dpa = Convert_to_double_and_filterbygenes(GRNs_2dpa, Genes_need_2dpa)
GRNs_double_4dpa = Convert_to_double_and_filterbygenes(GRNs_4dpa, Genes_need_4dpa)
GRNs_double_6dpa = Convert_to_double_and_filterbygenes(GRNs_6dpa, Genes_need_6dpa)

write.csv(GRNs_double_0dpa,file=paste0(celltype,"_GRNs_double_0dpa.csv"),row.names=F)
write.csv(GRNs_double_1dpa,file=paste0(celltype,"_GRNs_double_1dpa.csv"),row.names=F)
write.csv(GRNs_double_2dpa,file=paste0(celltype,"_GRNs_double_2dpa.csv"),row.names=F)
write.csv(GRNs_double_4dpa,file=paste0(celltype,"_GRNs_double_4dpa.csv"),row.names=F)
write.csv(GRNs_double_6dpa,file=paste0(celltype,"_GRNs_double_6dpa.csv"),row.names=F)




setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step6_new_identify_enriched_GRN_sub-networks")
# GRNs_triple_0dpa=read.csv(file=paste0(celltype,"_GRNs_triple_0dpa.csv"))
# GRNs_triple_1dpa=read.csv(file=paste0(celltype,"_GRNs_triple_1dpa.csv"))
# GRNs_triple_2dpa=read.csv(file=paste0(celltype,"_GRNs_triple_2dpa.csv"))
# GRNs_triple_4dpa=read.csv(file=paste0(celltype,"_GRNs_triple_4dpa.csv"))
# GRNs_triple_6dpa=read.csv(file=paste0(celltype,"_GRNs_triple_6dpa.csv"))

source("step6_new_filter_subGRNs.r")

##### Next load the log2 fold change files between stages #######

##### Differential DARs 
rrr = readRDS(paste0("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/plot/for_pub/Fig3/RRRs/",celltype,"_RRRs.rds"))
rrr = unique(as.character(rrr))
peak_avg_log2FC = readRDS("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/plot/for_pub/Fig3/peak_avg_log2FC/major_celltype_all_stage_0vsx_avglog2FC_peak_expression.rds")
peak_avg_log2FC =-peak_avg_log2FC  #change comparison direction, so log2FC>0 means higher expression in postinjury
rownames(peak_avg_log2FC) = sub("-",":",rownames(peak_avg_log2FC))
columns = c(paste0(celltype,"_0vs1_log2FC"),paste0(celltype,"_0vs2_log2FC"),paste0(celltype,"_0vs4_log2FC"),paste0(celltype,"_0vs6_log2FC"))
ct_rrr_avg_log2FC = as.data.frame(peak_avg_log2FC[rrr,columns])
ct_rrr_avg_log2FC$peaks =as.character(rownames(ct_rrr_avg_log2FC))
all_peak_log2fc = as.data.frame(peak_avg_log2FC[,columns])
all_peak_log2fc$peaks =as.character(rownames(all_peak_log2fc))

##### Differential DEGs 
gene_avg_log2FC = read.csv("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/plot/for_pub/Fig2/gene_avg_log2FC/major_celltype_all_stage_0vsx_avglog2FC_gene_expression.csv",row.names=1)
gene_avg_log2FC = -gene_avg_log2FC #change comparison direction, so log2FC>0 means higher expression in postinjury
rrg.list = readRDS("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/plot/for_pub/Fig2/major_celltype_allStage_DEGs_heatmap_modules/all_celltype_all_bystage_DEGs.rds")
rrg = rrg.list[[celltype]]
columns = c(paste0(celltype,"_0vs1_log2FC"),paste0(celltype,"_0vs2_log2FC"),paste0(celltype,"_0vs4_log2FC"),paste0(celltype,"_0vs6_log2FC"))
ct_rrg_avg_log2FC = as.data.frame(gene_avg_log2FC[rrg,columns])
ct_rrg_avg_log2FC$gene_short =as.character(rownames(ct_rrg_avg_log2FC))
all_gene_log2fc = as.data.frame(gene_avg_log2FC[,columns])
all_gene_log2fc$gene =rownames(all_gene_log2fc)


##### call specific GRNs by the fold change of Genes and Peaks ######
## in our case:
###1. have to be stage DARs in cell type
###2. log2FC relative to 0dpa need to be >0 

GRNs_triple_specific_1dpa <- Filter_triple_network_according_to_foldchange(GRNs_MG_triple = GRNs_triple_1dpa,DEGs = ct_rrg_avg_log2FC[,c(1,5)],DARs=ct_rrr_avg_log2FC[,c(1,5)],all_peak_log2fc = all_peak_log2fc[,c(1,5)])
GRNs_triple_specific_2dpa <- Filter_triple_network_according_to_foldchange(GRNs_MG_triple = GRNs_triple_2dpa,DEGs = ct_rrg_avg_log2FC[,c(2,5)],DARs=ct_rrr_avg_log2FC[,c(2,5)],all_peak_log2fc = all_peak_log2fc[,c(2,5)])
GRNs_triple_specific_4dpa <- Filter_triple_network_according_to_foldchange(GRNs_MG_triple = GRNs_triple_4dpa,DEGs = ct_rrg_avg_log2FC[,c(3,5)],DARs=ct_rrr_avg_log2FC[,c(3,5)],all_peak_log2fc = all_peak_log2fc[,c(3,5)])
GRNs_triple_specific_6dpa <- Filter_triple_network_according_to_foldchange(GRNs_MG_triple = GRNs_triple_6dpa,DEGs = ct_rrg_avg_log2FC[,c(4,5)],DARs=ct_rrr_avg_log2FC[,c(4,5)],all_peak_log2fc = all_peak_log2fc[,c(4,5)])
write.csv(GRNs_triple_specific_1dpa,file=paste0(celltype,"_GRNs_triple_specific_1dpa.csv"),row.names=F)
write.csv(GRNs_triple_specific_2dpa,file=paste0(celltype,"_GRNs_triple_specific_2dpa.csv"),row.names=F)
write.csv(GRNs_triple_specific_4dpa,file=paste0(celltype,"_GRNs_triple_specific_4dpa.csv"),row.names=F)
write.csv(GRNs_triple_specific_6dpa,file=paste0(celltype,"_GRNs_triple_specific_6dpa.csv"),row.names=F)




library(dplyr)
library(tibble)
library(ggplot2)
library(ggnetwork)
library(network)
library(sna)

setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step6_new_identify_enriched_GRN_sub-networks")


# GRNs_triple_specific_1dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_1dpa.csv"))
# GRNs_triple_specific_2dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_2dpa.csv"))
# GRNs_triple_specific_4dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_4dpa.csv"))
# GRNs_triple_specific_6dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_6dpa.csv"))

head(GRNs_triple_specific_1dpa,1)
#       TF    Motif                   Peak Target TF_Motif_Cor PtoG_Class
# 1 bach1b MA1633.2 chr5:18894307-18894901   aacs    0.1738943       Body
#     GG_Corr                               Triple type    TF_fold Target_fold
# 1 0.0774867 bach1b->chr5:18894307-18894901->aacs  Act 0.04541422  0.06823267
#   Peak_fold Peak_fold_2
# 1   

library(dplyr)
library(tibble)
library(ggplot2)
library(ggnetwork)
library(network)
library(sna)



setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step6_new_identify_enriched_GRN_sub-networks")

# GRNs_triple_specific_1dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_1dpa.csv"))
# GRNs_triple_specific_2dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_2dpa.csv"))
# GRNs_triple_specific_4dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_4dpa.csv"))
# GRNs_triple_specific_6dpa=read.csv(file=paste0(celltype,"_GRNs_triple_specific_6dpa.csv"))
GRNs_triple_specific_1dpa$stage = "1dpa"
GRNs_triple_specific_2dpa$stage = "2dpa"
GRNs_triple_specific_4dpa$stage = "4dpa"
GRNs_triple_specific_6dpa$stage = "6dpa"

GRN.all = rbind(GRNs_triple_specific_1dpa,GRNs_triple_specific_2dpa,GRNs_triple_specific_4dpa,GRNs_triple_specific_6dpa)
GRN.all = GRN.all %>% filter(abs(GG_Corr)>=0.2)

#myCols <- circlize::colorRamp2(seq(-0.5,0,1),colors = BuenColors::jdb_palette("solar_flare"))


library(ComplexHeatmap)
myCols= circlize::colorRamp2(c(-0.5, 0, 0.5), c("#2E7341", "white", "#54008B"))



GRNs = GRNs_triple_specific_1dpa %>% filter(abs(GG_Corr)>=0.2)
net.d <- GRNs %>% reshape2::dcast(Target ~ TF,value.var ="GG_Corr",fun.aggregate = mean) %>%
tibble::column_to_rownames("Target") %>% as.matrix()
net.d[is.nan(net.d)] <- 0
net.d=t(net.d)
p1 = Heatmap(net.d,
                  col=myCols,
                  clustering_distance_rows = "pearson",
                  clustering_distance_columns = "pearson",
                  name="Score",border = TRUE,
                  column_names_gp = gpar(fontsize=5,fontface="italic"),
                  row_names_gp = gpar(fontsize=5,fontface="italic"))





GRNs = GRNs_triple_specific_2dpa%>% filter(abs(GG_Corr)>=0.2)
net.d <- GRNs %>% reshape2::dcast(Target ~ TF,value.var ="GG_Corr",fun.aggregate = mean) %>%
tibble::column_to_rownames("Target") %>% as.matrix()
net.d[is.nan(net.d)] <- 0
net.d=t(net.d)
p2 = Heatmap(net.d,
                  col=myCols,
                  clustering_distance_rows = "pearson",
                  clustering_distance_columns = "pearson",
                  name="Score",border = TRUE,
                  column_names_gp = gpar(fontsize=5,fontface="italic"),
                  row_names_gp = gpar(fontsize=5,fontface="italic"))


GRNs = GRNs_triple_specific_4dpa%>% filter(abs(GG_Corr)>=0.2)
net.d <- GRNs %>% reshape2::dcast(Target ~ TF,value.var ="GG_Corr",fun.aggregate = mean) %>%
tibble::column_to_rownames("Target") %>% as.matrix()
net.d[is.nan(net.d)] <- 0
net.d=t(net.d)
p3 = Heatmap(net.d,
                  col=myCols,
                  clustering_distance_rows = "pearson",
                  clustering_distance_columns = "pearson",
                  name="Score",border = TRUE,
                  column_names_gp = gpar(fontsize=5,fontface="italic"),
                  row_names_gp = gpar(fontsize=5,fontface="italic"))


GRNs = GRNs_triple_specific_6dpa %>% filter(abs(GG_Corr)>=0.2)
net.d <- GRNs %>% reshape2::dcast(Target ~ TF,value.var ="GG_Corr",fun.aggregate = mean) %>%
tibble::column_to_rownames("Target") %>% as.matrix()
net.d[is.nan(net.d)] <- 0
net.d=t(net.d)
p4 = Heatmap(net.d,
                  col=myCols,
                  clustering_distance_rows = "pearson",
                  clustering_distance_columns = "pearson",
                  name="Score",border = TRUE,
                  column_names_gp = gpar(fontsize=5,fontface="italic"),
                  row_names_gp = gpar(fontsize=5,fontface="italic"))

pdf(paste0(celltype,"_TFtarget_heatmap_Cor0.2.pdf"))
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()


####### if plot key TF dotplot
gene_avg_log2FC = read.csv("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/plot/for_pub/Fig2/gene_avg_log2FC/major_celltype_all_stage_0vsx_avglog2FC_gene_expression.csv",row.names=1)
gene_avg_log2FC = -gene_avg_log2FC #change comparison direction, so log2FC>0 means higher expression in postinjury
colnames(rna.exp) = gsub("mfin","0dpa",colnames(rna.exp))
columns = c(paste0(celltype,"_0vs1_log2FC"),paste0(celltype,"_0vs2_log2FC"),paste0(celltype,"_0vs4_log2FC"),paste0(celltype,"_0vs6_log2FC"))


GRN.all = rbind(GRNs_triple_specific_1dpa,GRNs_triple_specific_2dpa,GRNs_triple_specific_4dpa,GRNs_triple_specific_6dpa) %>% filter(abs(GG_Corr)>=0.2)
ct.rna.exp = gene_avg_log2FC[unique(GRN.all$TF),columns]
colnames(ct.rna.exp)=c("1dpa","2dpa","4dpa","6dpa")
# GRN.all = GRN.all %>% filter(abs(Peak_fold_2)>0.5) 

dist_matrix <- dist(ct.rna.exp)
hc <- hclust(dist_matrix, method = "complete")
num_clusters <- 4  # or another number that makes sense for your data
clusters <- cutree(hc, k = num_clusters)
gene_order <- rownames(ct.rna.exp)[hc$order]


library(dplyr)
library(tidyverse)
count = GRN.all %>% group_by(stage,TF,type) %>% tally(name = "DuplicateCount")


ct.rna.exp <- tibble::rownames_to_column(ct.rna.exp, var = "gene")

# Use pivot_longer to convert to long format
ct.rna.exp.long <- ct.rna.exp %>% 
  pivot_longer(
    cols = -gene,  # Select all columns except 'gene' to pivot
    names_to = "stage", 
    values_to = "log2fc"
  )

colnames(ct.rna.exp.long)[1] <-  "TF" 
count <- count %>% left_join(ct.rna.exp.long, by = c("TF", "stage")) %>% arrange(stage,-DuplicateCount)
colnames(count)[4]="TargetCount"
write.csv(count,file=paste0(celltype,"_bystage_TF_targetCount_Cor0.2.csv"),row.names=F)

ordered_data <- count %>%arrange(type, desc(TargetCount))
ordered_data$TF <- factor(ordered_data$TF, levels = gene_order)

library(scales) 
pdf(paste0(celltype,"_TFtarget_count_bystage_dot_plot_Cor0.2.pdf"), width = 5, height = 7)
print(
    ggplot(ordered_data, aes(x = stage, y = TF, size = TargetCount, color = log2fc)) +
  geom_point() +scale_size(range = c(0.5, 4)) + 
  scale_color_gradient2(low = "#1B5C1D", high = "#5C2E64", mid = "white", midpoint = 0, limits = c(-1, 1), oob = squish) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,color="black"),
  axis.text.y = element_text(size=6,face="italic",color="black")) +
  labs(size = "Target Count", color = "Log2 Fold Change") +
  ggtitle(paste0(celltype," TF by Stage"))
  )
dev.off()

