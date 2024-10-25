
## GRN construction
```r
# reference: https://github.com/Pinlyu3/Zebrafish-retina-GRNs
suppressPackageStartupMessages({
        library(Signac)
        library(Seurat) 
        library(GenomeInfoDb)
        library(plyr)
        library(dplyr)
        library(tidyverse)
        library(rtracklayer)
        library(GenomicRanges)
        library(BSgenome.Drerio.UCSC.danRer11)
        library(future)
        set.seed(1234)
        plan("multicore", workers = 10)
        options(future.globals.maxSize = 400000000 * 1024^2)
})

# Please check GRN folder for rest functions
############ step1. TF gene expression and TF motif activity spearman correlation ############

source("source_function.r")
dir.create("step1_motif_TF_correlation")
setwd("step1_motif_TF_correlation")

require(TFBSTools)
# Add motif
JASPAR2022=readRDS("0_JASPAR_2022_vetebrate_nonredundant_pfm_matrix.rds")    #this can be load from JASPAR2022 package

# Add chromvar motif assay
# DefaultAssay(obj)="peaks"
# obj <- AddMotifs(
#   object = obj,
#   genome = BSgenome.Drerio.UCSC.danRer11,
#   pfm = JASPAR2022
# )
# # Run chromvar
# library(BiocParallel)
# register(SerialParam()) #have to add this, otherwise excessed mem...
# obj <- RunChromVAR(
#   object = obj,
#   genome = BSgenome.Drerio.UCSC.danRer11
# )
# #Extract motif matrix
# motif_matrix = obj@assays$chromvar@data
# rownames(motif_matrix) <- sub("_.+$", "", rownames(motif_matrix))
# #Extract gene expression matrix
# DefaultAssay(obj)="RNA"
# obj=NormalizeData(obj)
# gene_matrix = obj@assays$RNA@data

#Import motif gene table
Fish_All_motifs_table_short=read.csv("JASPAR2022_motif_gene_table_cleaned.csv")
colnames(Fish_All_motifs_table_short) = c("Motif","Gene")

# Get TFs gene motif correlation, gene gene correlation, gene matrix
# obj.list is a list of all major cell type seurat subset object
for (i in seq_along(obj.list)) {
    stage = names(obj.list)[i]
    stage.obj = obj.list[[i]]
    motif_matrix = stage.obj@assays$chromvar@data
    rownames(motif_matrix) <- sub("_.+$", "", rownames(motif_matrix))
    gene_matrix = stage.obj@assays$RNA@data
    peak_matrix = stage.obj@assays$peaks@data
    # save gene expression matrix for later grnboost2(TF-target inferrence) use
    #write.csv(motif_matrix,file=paste(celltype,stage,"motif_matrix.csv",sep="_"),row.names=T)
    #write.csv(peak_matrix,file=paste(celltype,stage,"peak_matrix.csv",sep="_"),row.names=T)

    # calculate motif-tfexp correaction in spearman

    gene_matrix = t(gene_matrix)
    k = which(colSums(gene_matrix) == 0)
    k
    gene_matrix = gene_matrix[,-k,drop = FALSE]
    dim(gene_matrix)
    write.csv(gene_matrix,file=paste(celltype,stage,"gene_exp_matrix_transposed.csv",sep="_"),row.names=F)

    Corr_RNA_Motif_Res = Get_corr_gene_motif(motif_matrix, gene_matrix,Fish_All_motifs_table_short)
    Corr_RNA_Motif_Res$TAG = paste(celltype,stage,sep="_")
    saveRDS(Corr_RNA_Motif_Res,file=paste(celltype,stage,'Corr_TFgene_Motif_Res.rds',sep="_"))
    colnames(Corr_RNA_Motif_Res)[c(3)] <- paste("Gene_Motif",colnames(Corr_RNA_Motif_Res)[c(3)],sep='_') 
    library(openxlsx)
    write.xlsx(Corr_RNA_Motif_Res, file = paste(celltype,stage,"TFs_gene_motif_corr.xlsx",sep="_"))

    Corr_res = sparse.cor3(gene_matrix)
    Corr_res = melt(Corr_res)
    saveRDS(Corr_res,file=paste0(celltype,"_",stage,"_Corr_res.rds"))
    # gene gene correlation 

    print(paste(i,celltype, stage,"TFs gene motif correlation, gene gene correlation, gene matrix saved!",sep=" "))
    }


###### grnboost2 to infer TF-target pairs
# Example:
python grnboost2.py Mesenchymal

###### Combined cor and grnboost2 result 
#!/usr/bin/Rscript
# Rscript combined celltype stage
library(dplyr)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

celltype = args[1]

for (stage in c("0dpa", "1dpa", "2dpa", "4dpa","6dpa")){
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
}


#################### Step2 p2g link, Signac 100kb distance, default setting for other parameter ########
# Example
Rscript p2g_preinjury_and_regeneration.r -o Mesenchymal_subset_object.rds  -p "Mesenchymal" 

## Then add peak annotation
require(GenomicRanges)
setwd("step2_p2g_links/")
source("add_P2G_peakAnnotation.r")

#Gtfs_protein_trans_GR_Body, it's a list of gene body info
load("step2_p2g_links/Gtfs_protein_trans_GR_Body")
#Gtfs_protein_trans_GR_TSS  it's a list of gene tss info
load("step2_p2g_links/Gtfs_protein_trans_GR_TSS")
All_Peaks_GRNs = readRDS("atac_all_peaks.rds")
GTF_TSS_table_res2 = Get_TSS_table(Gtfs_protein_trans_GR_TSS,All_Peaks_GRNs,extend=1000) #Promoter region: 1000bp up/down from TSS
GTF_Body_table_res = Get_Body_table(Gtfs_protein_trans_GR_Body,All_Peaks_GRNs)

celltypes= c("Superficial_Epithelial",
      "Intermediate_Epithelial",  
      "Basal_Epithelial",
      "Epidermal_Mucous",
      "Hematopoietic",
      "Mesenchymal")

for (celltype in celltypes){
    PG_cor = readRDS(paste0(celltype,"_regenerating_peak_gene_links_100kb.rds"))
    PG_cor$peak <- sub("-", ":", PG_cor$peak)
    PtoG_Ann = Get_all_peak_Gene_tables(GTF_TSS_table_res,GTF_Body_table_res,PtoG=PG_cor,Gtfs_protein_trans_GR_TSS)
    saveRDS(PtoG_Ann,file=paste0(celltype,"_regenerating_PtoG_Ann.rds"))

    PG_cor = readRDS(paste0(celltype,"_preinjury_peak_gene_links_100kb.rds"))
    PG_cor$peak <- sub("-", ":", PG_cor$peak)
    PtoG_Ann = Get_all_peak_Gene_tables(GTF_TSS_table_res,GTF_Body_table_res,PtoG=PG_cor,Gtfs_protein_trans_GR_TSS)
    saveRDS(PtoG_Ann,file=paste0(celltype,"_preinjury_PtoG_Ann.rds"))
}


###################### Step3 Combine p2g motif footprint together ###################
# filter tf motif cor by 0.05, then add tag
Corr_RNA_Motif_Res =readRDS(paste0(celltype,'_Corr_TFgene_Motif_Res.rds'))
Motif_with_tag = Get_all_need_Motif_tag(Corr_RNA_Motif_Res,0.05)


#JASPAR2022<-readRDS("0_JASPAR_2022_vetebrate_nonredundant_pfm_matrix.rds")
Fish_combined_motifs = JASPAR2022
motif_footprint <- Match_motifs(Peak_with_tag,Motif_with_tag,Fish_combined) #scan peaks for all motif region
#save(motif_footprint,file='motif_footprint')

#### select the activator and repressor motifs which passed the cutoff ####
#### cutoff = 0.05 #####
dir.create("step3_Predicting_TF_Binding_Sites/")
setwd("step3_Predicting_TF_Binding_Sites/")
source("step3_Predicting_TF_Binding_Sites.r")

celltypes= c("Superficial_Epithelial",
      "Intermediate_Epithelial",  
      "Basal_Epithelial",
      "Epidermal_Mucous",
      "Hematopoietic",
      "Mesenchymal")
stages=c("0dpa","1dpa","2dpa","4dpa","6dpa")

for (celltype in celltypes){
    for (stage in stages){
        TF_activity = readRDS(paste0("step1_motif_TF_correlation/",celltype,"_",stage,"_Corr_TFgene_Motif_Res.rds"))
        colnames(TF_activity)[3] ="Gene_Motif_Cor"
        colnames(TF_activity)[1] ="ID"
        Motif_with_tag = Get_all_need_Motif_tag(TF_activity, 0.05) #filter out TF with low correlation with motif activity

        ######## merge all peak regions to match TF binding motifs ##########
        if (stage =="0dpa"){
            PtoG_Anno = readRDS(paste0("step2_p2g_links/",celltype,"_preinjury_PtoG_Ann.rds"))
            } else {
                PtoG_Anno = readRDS(paste0("step2_p2g_links/",celltype,"_regenerating_PtoG_Ann.rds"))
                }

        Peak_with_tag = Get_all_need_Peak_tag(PtoG_Anno)

        library(TFBSTools)
        JASPAR2022<-readRDS("JASPAR_2022_vetebrate_nonredundant_pfm_matrix.rds")
        Fish_combined = JASPAR2022
        #all_peak_motif_match =readRDS("Total_footprint_Motif_in_WNNcelltype_calledpeaks.rds")
        motif_footprint <- Match_motifs(Peak_with_tag,Motif_with_tag,Fish_combined)
        saveRDS(motif_footprint,file=paste0(celltype,"_",stage,'_motif_footprint.rds'))
        print(paste(celltype,stage,"footprint at annotated p2g saved!",sep=" "))
        }
    }


#### STEP4 (merge to step3): Predicting cell-type specific TFs binding in cis-regulatory elements
##### after Add peak annotation to p2g, then add footprint signal 

#### loading corrected signal by rtracklayer #####
setwd("step3_Predicting_TF_Binding_Sites/")
source("step3_Predicting_TF_Binding_Sites.r")
# to generate TOBIAS corrected bigwig file
# Example:
# TOBIAS ATACorrect --cores 4 --read_shift 4 -5 --bam Mesenchymal_1dpa_atac_sorted.bam --outdir Mesenchymal_1dpa --genome genome.fa --peaks called_peaks.bed --blacklist Blacklist_danRer10_to_danRer11_YueLab_srt.bed > Mesenchymal_1dpa_bam2bw_TOBIAS.log 2>&1 &
celltypes= c("Superficial_Epithelial",
      "Intermediate_Epithelial",  
      "Basal_Epithelial",
      "Epidermal_Mucous",
      "Hematopoietic",
      "Mesenchymal")
stages=c("0dpa","1dpa","2dpa","4dpa","6dpa")

for (celltype in celltypes){
    for (stage in stages){
        motif_footprint = readRDS(paste0(celltype,'_',stage,'_motif_footprint.rds'))
        if (stage=="0dpa"){
            footprint_stage <- rtracklayer::import.bw(paste0("step4-3_TOBIAS_ATACcorrect/",celltype,"_mfin/",celltype,"_mfin_atac_sorted_corrected.bw"))
            } else {
                footprint_stage <- rtracklayer::import.bw(paste0("step4-3_TOBIAS_ATACcorrect/",celltype,"_",stage,"/",celltype,"_",stage,"_atac_sorted_corrected.bw"))
                }
                
            motif_cleaned= Get_footprint_Score_all_para(motif_footprint,footprint_signal=footprint_stage)
            saveRDS(motif_cleaned,file=paste0(celltype,"_",stage,'_motif_footprint_cleaned.rds'))
            print(paste(celltype,stage,"motif footprint cleaned!",sep=" "))
            }
        }



#### Step5_construct_TF_peak_target_links
##You need::
##TF RNA-motif correlation (spearman correlation):/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step1_motif_TF_correlation
##p2g /scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step2_p2g_links
##footprint and TF (TOBIAS) /scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step4-3_TOBIAS_ATACcorrect
##TF gene-target gene (grnboost2) and gene gene correlation /scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step4_new_TF_target_correlation

celltypes= c("Superficial_Epithelial",
      "Intermediate_Epithelial",  
      "Basal_Epithelial",
      "Epidermal_Mucous",
      "Hematopoietic",
      "Mesenchymal")
stages=c("0dpa","1dpa","2dpa","4dpa","6dpa")

for (celltype in celltypes){
    for (stage in stages){
#load("_TF_activity") ### from Step1 ####
    setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step1_motif_TF_correlation/")
    TF_activity <-readRDS(paste0(celltype,"_",stage,"_Corr_TFgene_Motif_Res.rds"))

    #load("_PtoG_Ann") ### from Step2 ####
    setwd("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step2_p2g_links")
    if (stage=="0dpa"){
        PtoG_Ann <-readRDS(paste0(celltype,"_preinjury_PtoG_Ann.rds"))
        } else {
            PtoG_Ann <-readRDS(paste0(celltype,"_regenerating_PtoG_Ann.rds"))
            }

    # load("Globle_motif_cleaned") ### from Step3 ####
    setwd("../step3_Predicting_TF_Binding_Sites")
    motif_cleaned = readRDS(paste0(celltype,"_",stage,'_motif_footprint_cleaned.rds'))

    #load("_network") ### from Step4_new ####
    setwd("../step4_new_TF_target_correlation/")
    network = readRDS(paste0(celltype,"_",stage,"_grnboost2_network_wCorr.rds"))

    # start to construct grn table
    setwd("../step5_new_construct_TF_peak_target_links")
    source("..source_function.r")

    GRNs_stage = Combined_all_things(GG_network=network,PG_cor=PtoG_Ann,TP_cor = TF_activity,Foot=motif_cleaned)
    saveRDS(GRNs_stage,file=paste0(celltype,"_",stage,"_GRNs.rds"))
    write.csv(GRNs_stage,file=paste0(celltype,"_",stage,"_GRNs.csv"),row.names=F)
    }
    }


#### Add RNA expression to GRN table 
require(openxlsx)
require(GenomicRanges)
require(dplyr)
require(tidyverse)
# gene log2FC compare to 0 dpa
gene_avg_log2FC = read.csv("major_celltype_all_stage_0vsx_avglog2FC_gene_expression.csv",row.names=1)
gene_avg_log2FC = -gene_avg_log2FC #change comparison direction, so log2FC>0 means higher expression in postinjury
colnames(gene_avg_log2FC) = gsub("mfin","0dpa",colnames(gene_avg_log2FC))
#columns = c(paste0(celltype,"_0vs1_log2FC"),paste0(celltype,"_0vs2_log2FC"),paste0(celltype,"_0vs4_log2FC"),paste0(celltype,"_0vs6_log2FC"))

#gene average expression 
gene.exp = read.csv("0_integrated_wnn_merged_celltype_rna_expression_bycelltype_bystage_avgCount_log1p.csv",row.names=1)
colnames(gene.exp) = gsub("mfin","0dpa",colnames(gene.exp))

#rrg annotation
rrg.list = readRDS("all_celltype_all_bystage_DEGs.rds")

# peak log2FC compare to 0 dpa
peak_avg_log2FC = readRDS("major_celltype_all_stage_0vsx_avglog2FC_peak_expression.rds")
peak_avg_log2FC =-peak_avg_log2FC  #change comparison direction, so log2FC>0 means higher expression in postinjury
rownames(peak_avg_log2FC) = sub("-",":",rownames(peak_avg_log2FC))
# columns = c(paste0(celltype,"_0vs1_log2FC"),paste0(celltype,"_0vs2_log2FC"),paste0(celltype,"_0vs4_log2FC"),paste0(celltype,"_0vs6_log2FC"))

# peak average expression 
peak.exp = read.csv("0_integrated_wnn_merged_celltype_peak_expression_bycelltype_bystage.csv",row.names=1)
rownames(peak.exp) = sub("-",":",rownames(peak.exp))
colnames(peak.exp) = gsub("mfin","0dpa",colnames(peak.exp))


setwd("step5_construct_TF_peak_target_links/GRN_wAllAnno")

celltypes= c("Superficial_Epithelial",
      "Intermediate_Epithelial",  
      "Basal_Epithelial",
      "Epidermal_Mucous",
      "Hematopoietic",
      "Mesenchymal")

stages<-c("0dpa","1dpa","2dpa","4dpa","6dpa")

for (celltype in celltypes){
  for (stage in stages){
    grn = readRDS(paste0("../",celltype,"_",stage,"_GRNs.rds"))
    
    if (stage!="0dpa"){
        # gene expressions and annotation
        gene.log2fc.column = paste0(celltype,"_0vs",sub("dpa","",stage),"_log2FC") 
        gene.avgexpLog1P.column = paste0(celltype, "_", stage)
        rrg = rrg.list[[celltype]]
        rrg = data.frame(gene = rrg, RRG_anno = paste0(celltype,"_RRG"))

        #peak expression and annotation
        peak.log2fc.column = paste0(celltype,"_0vs",sub("dpa","",stage),"_log2FC") 
        peak.avgexp.column = paste0(celltype, "_", stage)
        rrr = readRDS(paste0(celltype,"_RRRs.rds"))
        rrr$peak = as.character(rrr)
        rrr = data.frame(peak = rrr$peak, RRR_anno = rrr$name)

        #Add gene.log2fc, gene.avgexpLog1P.column, rrg, peak.log2fc.column,peak.avgexp.column
        #gene_avg_log2FC,gene.exp,rrg,peak_avg_log2FC,peak.exp,rrr
        grn$TF_log2fc_vs0dpa = gene_avg_log2FC[grn$TF, gene.log2fc.column]
        grn$Target_log2fc_vs0dpa = gene_avg_log2FC[grn$Target, gene.log2fc.column]

        grn$TF_AvgExp_log1p = gene.exp[grn$TF, gene.avgexpLog1P.column]
        grn$Target_AvgExp_log1p = gene.exp[grn$Target, gene.avgexpLog1P.column]

        grn$TF_RRG_anno = paste0(celltype,"_NonRRG") 
        grn$Target_RRG_anno = paste0(celltype,"_NonRRG") 
        grn$TF_RRG_anno[grn$TF %in% rrg$gene] = paste0(celltype,"_RRG")
        grn$Target_RRG_anno[grn$Target %in% rrg$gene] = paste0(celltype,"_RRG")

        grn$peak_log2fc_vs0dpa = peak_avg_log2FC[grn$Peak,peak.log2fc.column]
        grn$peak_AvgExp = peak.exp[grn$Peak,peak.avgexp.column]

        grn <- merge(grn, rrr[, c("peak", "RRR_anno")], by.x = "Peak", by.y="peak", all.x = TRUE)
        grn$RRR_anno[is.na(grn$RRR_anno)] =paste0(celltype,"_NonRRR")         
        } else {
            #gene.log2fc.column = paste0(celltype,"_0vs",sub("dpa","",stage),"_log2FC") 
            gene.avgexpLog1P.column = paste0(celltype, "_", stage)
            rrg = rrg.list[[celltype]]
            rrg = data.frame(gene = rrg, RRG_anno = paste0(celltype,"_RRG"))

            #peak expression and annotation
            #peak.log2fc.column = paste0(celltype,"_0vs",sub("dpa","",stage),"_log2FC") 
            peak.avgexp.column = paste0(celltype, "_", stage)
            rrr = readRDS(paste0("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/plot/for_pub/Fig3/RRRs/",celltype,"_RRRs.rds"))
            rrr$peak = as.character(rrr)
            rrr = data.frame(peak = rrr$peak, RRR_anno = rrr$name)

            #Add gene.log2fc, gene.avgexpLog1P.column, rrg, peak.log2fc.column,peak.avgexp.column
            #gene_avg_log2FC,gene.exp,rrg,peak_avg_log2FC,peak.exp,rrr
            grn$TF_log2fc_vs0dpa = 0
            grn$Target_log2fc_vs0dpa = 0

            grn$TF_AvgExp_log1p = gene.exp[grn$TF, gene.avgexpLog1P.column]
            grn$Target_AvgExp_log1p = gene.exp[grn$Target, gene.avgexpLog1P.column]

            grn$TF_RRG_anno = paste0(celltype,"_NonRRG") 
            grn$Target_RRG_anno = paste0(celltype,"_NonRRG") 
            grn$TF_RRG_anno[match(rrg$gene, grn$TF)] = paste0(celltype,"_RRG")
            grn$Target_RRG_anno[match(rrg$gene, grn$Target)] = paste0(celltype,"_RRG")

            grn$peak_log2fc_vs0dpa = 0
            grn$peak_AvgExp = peak.exp[grn$Peak,peak.avgexp.column]
            # grn$RRR_anno = paste0(celltype,"_NonRRR") 
            # grn$RRR_anno[which(grn$Peak %in% rrr$peak)] = rrr$RRR_anno[which(rrr$peak %in% grn$Peak)]   
            grn$RRR_anno = rrr$RRR_anno[match(grn$Peak,rrr$peak)]  
            grn$RRR_anno[is.na(grn$RRR_anno)] =paste0(celltype,"_NonRRR") 
            }

        grn$ct_stage=paste(celltype,stage,sep="_")
        saveRDS(grn,paste(celltype,stage,"GRN_wAllAnno.rds", sep="_"))
        write.csv(grn,paste(celltype,stage,"GRN_wAllAnno.csv", sep="_"),row.names=F)
        print(paste(celltype, stage,"GRN with annotation done!",sep=" "))
    }
    }




## Dot Plot to show top TFs 
### Normalize centrality degree by the maximum degree centrality observed in that stage.
# Set working directory and load required libraries
setwd("../plot")
require(igraph)
require(dplyr)
require(tidyverse)
require(ggplot2)
require(scales)

# Read RNA expression data and adjust column names
rna.exp <- read.csv("0_integrated_wnn_merged_celltype_rna_expression_bycelltype_bystage_avgCount_log1p.csv", row.names = 1)
colnames(rna.exp) <- gsub("mfin", "0dpa", colnames(rna.exp))

# Define stages and cell types
stages <- c("0dpa", "1dpa", "2dpa", "4dpa", "6dpa")
celltypes <- c("Superficial_Epithelial", "Intermediate_Epithelial", "Basal_Epithelial", "Epidermal_Mucous", "Hematopoietic", "Mesenchymal")

for (celltype in celltypes) {
  ct_tf_list <- list()
  for (stage in stages) {
    grn <- readRDS(paste0("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step5_new_construct_TF_peak_target_links/GRN_wAllAnno/", celltype, "_", stage, "_GRN_wAllAnno.rds"))
    
    # Apply filters based on stage
    if (stage == "0dpa") {
      grn.filter <- grn %>%
        filter(TF_AvgExp_log1p > 0) %>%
        filter(Target_AvgExp_log1p > 0) %>%
        filter(peak_AvgExp > 0)
    } else {
      grn.filter <- grn %>%
        filter(TF_AvgExp_log1p > 0) %>%
        filter(abs(TF_log2fc_vs0dpa) >= 0.25) %>%
        filter(Target_AvgExp_log1p > 0) %>%
        filter(peak_AvgExp > 0) %>%
        filter(peak_log2fc_vs0dpa > 0) %>%
        filter(abs(GG_Corr) >= 0.1)
    }

    df <- grn.filter %>% select(TF, Target)
    g <- graph_from_data_frame(df)
    degree_centrality <- degree(g, mode = "out")
    max_degree <- max(degree_centrality, na.rm = TRUE)
    normalized_centrality <- degree_centrality / max_degree #Normalize centrality degree by the maximum degree centrality observed in that stage.
    
    centrality_df <- data.frame(TF = names(normalized_centrality), centrality_degree = normalized_centrality) %>%
                     filter(centrality_degree > 0) %>%
                     arrange(desc(centrality_degree)) %>%
                     head(30) %>%
                     mutate(stage = stage)
    
    tf.exp <- grn %>% select(TF, TF_log2fc_vs0dpa, TF_AvgExp_log1p) %>% unique()
    centrality_df$log2FC <- tf.exp$TF_log2fc_vs0dpa[match(centrality_df$TF, tf.exp$TF)]
    ct.stage.exp <- grn.filter %>% select(TF, TF_AvgExp_log1p) %>% unique()
    centrality_df$avg_exp <- ct.stage.exp$TF_AvgExp_log1p[match(centrality_df$TF, ct.stage.exp$TF)]
    
    ct_tf_list[[stage]] <- centrality_df
  }
  
  # Combine and prepare data for plotting
  ordered_data <- do.call(rbind, ct_tf_list)
  ordered_data$TF <- factor(ordered_data$TF, levels = rev(unique(ordered_data$TF)))

  # Generate dot plot for Log2 Fold Change
  pdf(paste0(celltype, "_TFtarget_count_bystage_dot_plot_top30_by_centrality_out_log2FC_GGcor0.1_maxNorm.pdf"), width = 4, height = 8)
  print(
    ggplot(ordered_data, aes(x = stage, y = TF, size = centrality_degree, fill = log2FC)) +
      geom_point(shape = 21, stroke = 0.5) +
      scale_size(range = c(1, 5)) +
      scale_fill_gradient2(low = "#1B5C1D", high = "#5C2E64", mid = "white", midpoint = 0, limits = c(-1, 1), oob = squish) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, color = "black"),
            axis.text.y = element_text(size = 10, face = "italic", color = "black")) +
      labs(size = "Centrality degree", fill = "Log2 Fold Change") +
      ggtitle(paste0(celltype, " TF by Stage"))
  )
  dev.off()

  # Generate dot plot for Average Expression
  pdf(paste0(celltype, "_TFtarget_count_bystage_dot_plot_top30_by_centrality_out_avgExp_log1P_GGcor0.1_maxNorm.pdf"), width = 4, height = 8)
  print(
    ggplot(ordered_data, aes(x = stage, y = TF, size = centrality_degree, fill = avg_exp)) +
      geom_point(shape = 21, stroke = 0.5) +
      scale_size(range = c(1, 5)) +
      scale_fill_gradient2(high = "#E7D39A", low = "#781399", mid = "#F06D70", midpoint = 1.5, limits = c(0, round(max(ordered_data$avg_exp), 0)), oob = squish) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, color = "black"),
            axis.text.y = element_text(size = 10, face = "italic", color = "black")) +
      labs(size = "Centrality degree", fill = "Log1p Expression") +
      ggtitle(paste0(celltype, " TF by Stage"))
  )
  dev.off()
  
  print(paste(celltype, "plot done!"))
}
```

## Fig5A
```r
stages<-c("0dpa", "1dpa","2dpa","4dpa","6dpa")
celltype = "Mesenchymal"
  grn.list = list()
  for (stage in stages){
    #stage="6dpa"
    grn = readRDS(paste0("step5_new_construct_TF_peak_target_links/GRN_wAllAnno/",celltype,"_",stage,"_GRN_wAllAnno.rds"))
    grn.filter = grn %>% 
    #filter(Target_AvgExp_log1p>0) %>% 
    filter(peak_AvgExp>0) #%>% 
    #filter(peak_log2fc_vs0dpa >0) #%>%
    grn.list[[stage]]=grn
    }

grn.all = do.call("rbind",grn.list)
write.csv(grn.all,file=paste0(celltype,"_GRN_wAllAnno_allstage.csv"))

chromvar.mk=read.csv(paste0(celltype,"_chromvar_findAllMarkers_with_Avg_gene_expression_table.csv"),row.names=1)
avg_chromvar <- read.csv("Mesenchymal_chromvar_by_stage_avgexp.csv",row.names=1)

# sort the heatmap prior to plotting find column index for maxium value for each TF activity
avg_chromvar=avg_chromvar[unique(chromvar.mk$motifID_gene),]

avg_chromvar = as.matrix(t(scale(t(avg_chromvar))))
# Convert to dataframe for manipulation
avg_chromvar_df <- as.data.frame(avg_chromvar)
# Add 'max' column
avg_chromvar_df$max <- max.col(avg_chromvar_df)
# Sort by 'max' column
avg_chromvar_df <- avg_chromvar_df[order(avg_chromvar_df$max), ]
# Remove 'max' column
avg_chromvar_df <- dplyr::select(avg_chromvar_df, -max)
avg_chromvar <- as.matrix(avg_chromvar_df)
colnames(avg_chromvar) <- c("0dpa","1dpa","2dpa","4dpa","6dpa")

selected_motifID_gene = c('MA0489.2_Jun','MA0476.1_FOS' ,'MA0477.2_FOSL1','MA1141.1_FOS::JUND','MA1142.1_FOSL1::JUND','MA1137.1_FOSL1::JUNB','MA0490.2_JUNB','MA1128.1_FOSL1::JUN',
'MA0511.2_RUNX2','MA1990.1_Gli1','MA0734.3_Gli2',
'MA0899.1_HOXA10','MA0651.2_HOXC11','MA0905.1_HOXC10','MA0901.2_HOXB13','MA1506.1_HOXD10',"MA0873.1_HOXD12",
'MA0703.2_LMX1B','MA0700.2_LHX2','MA0666.2_MSX1','MA0709.1_Msx3','MA0708.2_MSX2',
'MA0887.1_EVX1','MA0888.1_EVX2',
'MA0716.1_PRRX1','MA0497.1_MEF2C',
'MA1628.1_Zic1::Zic2','MA1629.1_Zic2')

require(ComplexHeatmap)
position_to_lable1 = which(rownames(avg_chromvar) %in% selected_motifID_gene)
RowAnnotation1 = rowAnnotation(genes = anno_mark(at = position_to_lable1, labels = rownames(avg_chromvar)[position_to_lable1], labels_gp = gpar(fontsize = 8)))

require(colorRamp2)
#instead of col_name_order
ht1 = Heatmap(avg_chromvar, name = "motif activity", 
        row_title = "TF", column_title = "motif activity",cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE,show_column_names = TRUE,
        column_names_rot = 0,
        col =colorRamp2(c(-1.5, 0, 1.5), c( "#3F3869", "#FDF4EF","#BC7A3C")),#c("#BC7A3C", "#FDF4EF", "#3F3869")),
        right_annotation = RowAnnotation1)
pdf("Mesenchymal_chromvar_findallmark_stage_mk_heatmap_wTFlabels2.pdf",width=5,height=4)
ht1
dev.off()
```

## Fig5B
```r
library(igraph)
labeled_genes <-c("msx1b","msx3","prrx1a","twist1a","twist2","prdm1a","runx2a","lmx1ba","lhx2b","zic2a","zic2b")
stages<-c("1dpa","2dpa","4dpa","6dpa")
celltypes= c("Superficial_Epithelial",
      "Intermediate_Epithelial",  
      "Basal_Epithelial",
      "Epidermal_Mucous",
      "Hematopoietic",
      "Mesenchymal")
celltype="Mesenchymal"
ct_tf_list = list()
for (stage in stages){
    grn = readRDS(paste0("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step5_new_construct_TF_peak_target_links/GRN_wAllAnno/",celltype,"_",stage,"_GRN_wAllAnno.rds"))
    if (stage =="0dpa"){
        grn.filter = grn %>% 
        filter (TF_AvgExp_log1p> 0) %>%
        filter(Target_AvgExp_log1p>0) %>% 
        filter(peak_AvgExp>0)
        }else{
            grn.filter = grn %>% 
            filter (TF_AvgExp_log1p> 0) %>%
            filter (abs(TF_log2fc_vs0dpa)>= 0.25 ) %>% 
            filter(Target_AvgExp_log1p>0) %>% 
            filter(peak_AvgExp>0) %>% 
            filter(peak_log2fc_vs0dpa >0) 
    }

    df = grn.filter %>% select(TF,Target)
    g <- graph_from_data_frame(df) 
    degree_centrality <- degree(g,mode = "out")
    ct_tf_list[[stage]] = degree_centrality[degree_centrality > 0] %>%sort(decreasing= T) %>% data.frame()
    ct_tf_list[[stage]] = data.frame("TF"=rownames(ct_tf_list[[stage]]), "centrality_degree"=ct_tf_list[[stage]][,1])
    ct_tf_list[[stage]]$stage = paste0(stage)
    tf.exp = grn %>% select(TF,TF_log2fc_vs0dpa,TF_AvgExp_log1p) %>% unique()

    ct_tf_list[[stage]]$log2FC = tf.exp$TF_log2fc_vs0dpa[match(ct_tf_list[[stage]]$TF,tf.exp$TF)]

    ct.stage.exp = grn.filter %>%select(TF,TF_AvgExp_log1p) %>% unique()
    ct_tf_list[[stage]]$avg_exp = ct.stage.exp$TF_AvgExp_log1p[match(ct_tf_list[[stage]]$TF,ct.stage.exp$TF)]    
    }
    # generate dotplot
    ordered_data <-  do.call(rbind, ct_tf_list)
    ordered_data = ordered_data[ordered_data$TF %in% labeled_genes,]
    ordered_data$TF <- factor(ordered_data$TF, levels = rev(unique(ordered_data$TF)))   
    library(scales) 
    library(ggplot2)  
    pdf(paste0(celltype, "_selected_TFtarget_count_bystage_dot_plot_top30_by_centrality_out_avgExp_log1P.pdf"), width = 4, height = 4)
    print(
    ggplot(ordered_data, aes(x = stage, y = TF, size = centrality_degree, fill = avg_exp)) +
    geom_point(shape = 21, stroke = 0.5) + 
    scale_size(range = c(1, 5)) + 
    scale_fill_gradient2(high = "#E7D39A", low = "#781399", mid = "#F06D70", midpoint = 1.5, limits = c(0, round(max(ordered_data$avg_exp),0)), oob = squish) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 10, face = "italic", color = "black")) +
    labs(size = "Centrality degree", fill = "Log1p Expression") +
    ggtitle(paste0(celltype, " TF by Stage"))
    )
    dev.off()
```

## Fig5C
```R
require(ggplot2)
require(ggrepel)
stages<-c("0dpa", "1dpa","2dpa","4dpa","6dpa")

# Assuming your dataframe is named grn.filter
#for (celltype in celltypes){
celltype = "Mesenchymal"
  ct_tf_list = list()
  for (stage in stages){
    #stage="6dpa"
    grn = readRDS(paste0("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step5_new_construct_TF_peak_target_links/GRN_wAllAnno/",celltype,"_",stage,"_GRN_wAllAnno.rds"))
    grn.filter = grn %>% 
    filter(Target_AvgExp_log1p>0) %>% 
    filter(peak_AvgExp>0) %>% 
    filter(peak_log2fc_vs0dpa >0)

# Filter for runx2a transcription factor
    tf = "runx2a"
    tf_data <- subset(grn.filter, TF == {tf}) %>% select(TF, Target, TF_Gene,GG_Score, GG_Corr,TF_log2fc_vs0dpa,Target_log2fc_vs0dpa,TF_AvgExp_log1p,Target_AvgExp_log1p) %>% unique()
    pdf(paste(celltype,stage,tf,"target_gene_scatter_plot_GGscore.pdf",sep="_"),width=5,height=4)
 
    top_tf_targets <- tf_data[order(tf_data$GG_Score, tf_data$Target_log2fc_vs0dpa, decreasing = TRUE), ][1:30,]

    print(
      ggplot(tf_data, aes(y =GG_Score, x = Target_log2fc_vs0dpa)) +
      geom_point() +
      ggrepel::geom_label_repel(data = top_tf_targets, 
                       aes(label = Target), 
                       max.overlaps = Inf,
                       box.padding = 0.35, 
                       point.padding = 0.5,
                       segment.color = 'grey50', color = "darkred") +
                       theme_minimal() +
                       theme(axis.line = element_line(color = "black")) +
                       labs(
                        title = paste0(tf," targeted Genes in ",celltype," ",stage),
                        y = "GG_Score(TF-Target)",
                        x = "Target_log2fc_vs_0dpa")
        )
    dev.off()
}

```

## Fig5D
```markdown
Screenshot from WashU Epigenome Browser
```
