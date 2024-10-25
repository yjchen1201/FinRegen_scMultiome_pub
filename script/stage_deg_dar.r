#!/usr/bin/env Rscript
# example: Rscript stage_deg_dar.r seurat_object.rds celltype_name
args = commandArgs(trailingOnly=TRUE)
obj<-readRDS(paste(args[1]))
celltype = (paste(args[2]))
sample="WNN_celltype"

# load packages
suppressPackageStartupMessages({
library(Signac)
library(Seurat)
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Drerio.UCSC.danRer11)
library(paletteer)
set.seed(1234)
})

obj=readRDS("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/WNN_celltype_subset_object/BE_subset_object.rds")

#DEG
# functions
GetDEGs<-function(celltype,object,assay,ident1,ident2) {
  #sample=sample
  obj = object
  assay= assay
  compare <-paste0(ident1,"_vs_",ident2)
  deg <- FindMarkers(obj, 
                    assay=assay,
                    group.by = 'Stage',
                     #subset.ident=celltype, 
                     ident.1 = ident1,
                     ident.2 = ident2,    
                     test.use = "wilcox", 
                     logfc.threshold = 0.25,
                     min.pct = 0.1) %>% dplyr::filter(p_val_adj < 0.05)%>% tibble::rownames_to_column(var="gene")

  deg$compare <- paste(celltype,compare,sep="_")
  deg$category <-"DEG"
  deg.up<- deg %>% dplyr::filter(avg_log2FC < -0.25) %>% dplyr::arrange(avg_log2FC) 
  deg.down<- deg %>% dplyr::filter(avg_log2FC > 0.25) %>% dplyr::arrange(-avg_log2FC) 
  return(list(deg_all=deg,deg_up=deg.up,deg_down=deg.down))
}


# DEGs
assay <-"RNA"
DefaultAssay(obj)<-assay
obj<-NormalizeData(obj)
celltype= unique(obj$celltype_wnn)


DEG_assay="RNA"
print(paste("Find DEGs for ",celltype,"...",sep=""))
### prevent writing results from last loop to current loop
deg_0vs1<-list()
deg_0vs2<-list()
deg_0vs4<-list()
deg_0vs6<-list()
deg_1vs2<-list()
deg_1vs4<-list()
deg_1vs6<-list()
deg_2vs4<-list()
deg_2vs6<-list()
deg_4vs6<-list()


print(paste("Find DEGs for ",celltype," 0vs1...",sep=""))
tryCatch(deg_0vs1<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="mfin",ident2="1dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 0vs2...",sep=""))
tryCatch(deg_0vs2<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="mfin",ident2="2dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 0vs4...",sep=""))
tryCatch(deg_0vs4<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="mfin",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 0vs6...",sep=""))
tryCatch(deg_0vs6<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="mfin",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 1vs2...",sep=""))
tryCatch(deg_1vs2<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="1dpa",ident2="2dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 1vs4...",sep=""))
tryCatch(deg_1vs4<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="1dpa",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 1vs6...",sep=""))
tryCatch(deg_1vs6<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="1dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 2vs4...",sep=""))
tryCatch(deg_2vs4<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="2dpa",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 2vs6...",sep=""))
tryCatch(deg_2vs6<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="2dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DEGs for ",celltype," 4vs6...",sep=""))
tryCatch(deg_4vs6<-GetDEGs(celltype,object=obj,assay=DEG_assay,ident1="4dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


deg.list<-list()
deg.list<-list(
    deg_0vs1=deg_0vs1$deg_all, deg_0vs1_up=deg_0vs1$deg_up, deg_0vs1_down=deg_0vs1$deg_down,
    deg_0vs2=deg_0vs2$deg_all, deg_0vs2_up=deg_0vs2$deg_up, deg_0vs2_down=deg_0vs2$deg_down,
    deg_0vs4=deg_0vs4$deg_all, deg_0vs4_up=deg_0vs4$deg_up, deg_0vs4_down=deg_0vs4$deg_down,
    deg_0vs6=deg_0vs6$deg_all, deg_0vs6_up=deg_0vs6$deg_up, deg_0vs6_down=deg_0vs6$deg_down,
    deg_1vs2=deg_1vs2$deg_all, deg_1vs2_up=deg_1vs2$deg_up, deg_1vs2_down=deg_1vs2$deg_down, 
    deg_1vs4=deg_1vs4$deg_all, deg_1vs4_up=deg_1vs4$deg_up, deg_1vs4_down=deg_1vs4$deg_down, 
    deg_1vs6=deg_1vs6$deg_all, deg_1vs6_up=deg_1vs6$deg_up, deg_1vs6_down=deg_1vs6$deg_down, 
    deg_2vs4=deg_2vs4$deg_all, deg_2vs4_up=deg_2vs4$deg_up, deg_2vs4_down=deg_2vs4$deg_down, 
    deg_2vs6=deg_2vs6$deg_all, deg_2vs6_up=deg_2vs6$deg_up, deg_2vs6_down=deg_2vs6$deg_down, 
    deg_4vs6=deg_4vs6$deg_all, deg_4vs6_up=deg_4vs6$deg_up, deg_4vs6_down=deg_4vs6$deg_down
    )

require(openxlsx)
print(paste("write DEGs to xlsx ",celltype,"...",sep=""))
write.xlsx(deg.list, file = paste(sample,celltype,"DEGs_stage_comparison_usingRNAassay_renormalized.xlsx",sep="_"), sheetName = names(deg.list), rowNames = F,overwrite=T)
stat <-lapply(deg.list, function(x) {print(nrow(x))}) %>% unlist() %>% data.frame() 
write.table(stat,file=paste(sample,celltype,"DEGs_stage_comparison_usingRNAassay_renormalized_stats.txt",sep="_"),quote=F,sep="\t",col.names=F)



# DARs
GetDARs<-function(celltype,object,assay,ident1,ident2) {
  obj = object
  assay= assay
  compare <-paste0(ident1,"_vs_",ident2)
  dar <- FindMarkers(obj, 
                    assay=assay,
                    group.by = 'Stage',
                     #subset.ident=celltype, 
                     ident.1 = ident1,
                     ident.2 = ident2,    
                     test.use = 'LR', 
                     latent.vars = paste("nCount",assay,sep="_"),
                     logfc.threshold = 0.25,
                     min.pct = 0.1) %>% dplyr::filter(p_val_adj < 0.05)

  cf <- ClosestFeature(object=obj,regions=rownames(dar), annotation=Annotation(obj), sep=c(':','-')) #make sure that you have annotaion added to seurat object
  dar <- cbind(dar, gene=cf$gene_name, distance=cf$distance) 
  dar$compare <- paste(celltype,compare,sep="_")
  dar$category <-"DAR"
  dar <- tibble::rownames_to_column(dar,var="region")
  dar.open<- dar %>% dplyr::filter(avg_log2FC < -0.25) %>% dplyr::arrange(avg_log2FC) 
  dar.close<- dar %>% dplyr::filter(avg_log2FC > 0.25) %>% dplyr::arrange(-avg_log2FC) 
  
  return(list(dar_all=dar,dar_open=dar.open,dar_close=dar.close))
}


print("Get DARs")
sample="WNNcelltype"
assay="peaks"
DefaultAssay(obj) <- assay


DefaultAssay(obj) <- "peaks"
obj <- FindTopFeatures(obj, min.cutoff = 10)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

DAR_assay=assay
print(paste("Find DARs for ",celltype,"...",sep=""))

### prevent writing results from last loop to current loop
dar_0vs1<-list()
dar_0vs2<-list()
dar_0vs4<-list()
dar_0vs6<-list()
dar_1vs2<-list()
dar_1vs4<-list()
dar_1vs6<-list()
dar_2vs4<-list()
dar_2vs6<-list()
dar_4vs6<-list()


print(paste("Find DARs for ",celltype," 0vs1...",sep=""))
tryCatch(dar_0vs1<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="mfin",ident2="1dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 0vs2...",sep=""))
tryCatch(dar_0vs2<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="mfin",ident2="2dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 0vs4...",sep=""))
tryCatch(dar_0vs4<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="mfin",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 0vs6...",sep=""))
tryCatch(dar_0vs6<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="mfin",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 1vs2...",sep=""))
tryCatch(dar_1vs2<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="1dpa",ident2="2dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 1vs4...",sep=""))
tryCatch(dar_1vs4<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="1dpa",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 1vs6...",sep=""))
tryCatch(dar_1vs6<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="1dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 2vs4...",sep=""))
tryCatch(dar_2vs4<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="2dpa",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 2vs6...",sep=""))
tryCatch(dar_2vs6<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="2dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",celltype," 4vs6...",sep=""))
tryCatch(dar_4vs6<-GetDARs(celltype,object=obj,assay=DAR_assay,ident1="4dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

dar.list<-list()
dar.list<-list(
    dar_0vs1=dar_0vs1$dar_all, dar_0vs1_open=dar_0vs1$dar_open, dar_0vs1_close=dar_0vs1$dar_close,
    dar_0vs2=dar_0vs2$dar_all, dar_0vs2_open=dar_0vs2$dar_open, dar_0vs2_close=dar_0vs2$dar_close,
    dar_0vs4=dar_0vs4$dar_all, dar_0vs4_open=dar_0vs4$dar_open, dar_0vs4_close=dar_0vs4$dar_close,
    dar_0vs6=dar_0vs6$dar_all, dar_0vs6_open=dar_0vs6$dar_open, dar_0vs6_close=dar_0vs6$dar_close,
    dar_1vs2=dar_1vs2$dar_all, dar_1vs2_open=dar_1vs2$dar_open, dar_1vs2_close=dar_1vs2$dar_close, 
    dar_1vs4=dar_1vs4$dar_all, dar_1vs4_open=dar_1vs4$dar_open, dar_1vs4_close=dar_1vs4$dar_close, 
    dar_1vs6=dar_1vs6$dar_all, dar_1vs6_open=dar_1vs6$dar_open, dar_1vs6_close=dar_1vs6$dar_close, 
    dar_2vs4=dar_2vs4$dar_all, dar_2vs4_open=dar_2vs4$dar_open, dar_2vs4_close=dar_2vs4$dar_close, 
    dar_2vs6=dar_2vs6$dar_all, dar_2vs6_open=dar_2vs6$dar_open, dar_2vs6_close=dar_2vs6$dar_close, 
    dar_4vs6=dar_4vs6$dar_all, dar_4vs6_open=dar_4vs6$dar_open, dar_4vs6_close=dar_4vs6$dar_close
    )

require(openxlsx)
print(paste("write DARs to xlsx ",celltype,"...",sep=""))
write.xlsx(dar.list, file = paste(sample,celltype,"DARs_stage_comparison_renormalized.xlsx",sep="_"), sheetName = names(dar.list), rowNames = F,overwrite=T)
stat <-lapply(dar.list, function(x) {print(nrow(x))}) %>% unlist() %>% data.frame() 
write.table(stat,file=paste(sample,celltype,"DARs_stage_comparison_renormalized_stats.txt",sep="_"),quote=F,sep="\t", col.names=F)
