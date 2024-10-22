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

obj=readRDS("0_integrated_wCelltypePeaks.rds")
Idents(obj)=obj$wnn_clusters

Idents(obj)="wnn_clusters"
obj <- RenameIdents(
   obj, 
    '1' = 'Superficial_Epithelial',
    '10' = 'Superficial_Epithelial',
    '0' = 'Intermediate_Epithelial',
    '3' = 'Intermediate_Epithelial',
    '4' = 'Intermediate_Epithelial',
    '7' = 'Intermediate_Epithelial',
    '2' = 'Basal_Epithelial',
    '5' = 'Basal_Epithelial',
    '11' = 'Epidermal_Mucous_1',
    '14' = 'Epidermal_Mucous_2',
    '9' = 'Hematopoietic_1',
    '13' = 'Hematopoietic_2',
    '6' = 'Mesenchymal_1',
    '8' = 'Mesenchymal_2',
    '12' = 'Unknown_1',
    '15' = 'Metaphocyte',
    '16' = 'Pigment',
    '17' = 'Endothelial',
    '18' = 'Unknown_2'
    ) 
obj$celltype_final = Idents(obj)
obj$celltype_final = factor(
    obj$celltype_final,
    levels = c('Superficial_Epithelial',
    "Intermediate_Epithelial", 
    "Basal_Epithelial",
    "Epidermal_Mucous_1",
    "Epidermal_Mucous_2",
    'Metaphocyte',
    "Hematopoietic_1",
    "Hematopoietic_2",
    "Mesenchymal_1",
    "Mesenchymal_2",
    'Endothelial',
    "Pigment",
    'Unknown_1',
    'Unknown_2'
    ))

Idents(obj) = "celltype_final"
celltype_final<-levels(obj$celltype_final)

celltype_color_final<-c(
    "#8B2DB2", #SE
    "#CE6DBD", #IE
    "#9C9EDE", #BE
    "#F28E2BFF", #EM1
    "#FFBE7DFF" ,#EM2
    "#E15759FF" ,#Metaphocyte
    "#59A14FFF" ,#HEM1
    "#8CD17DFF", #HEM2
    "#B6992DFF" ,#MES1
    "#F1CE63FF" ,#MES2
    "#33608C",   #Endo
    "#86BCB6FF", #Pigment
    "#FF9D9AFF", #C12 Tuftlike  unknown_1
    "#79706EFF"#C18 neuromast like unknown_2 
    ) 

names(celltype_color_final)<-celltype_final


# Fig 1B
pdf(paste0("WNN_integrated_wnnUMAP_plot_byCelltype.pdf"),width = 12, height = 9)
print(DimPlot(obj, cols = celltype_color_final,reduction = "wnn.umap", label = F, label.size = 12, pt.size =0.1, repel = TRUE) + ggtitle("Colored_by_celltype") +
theme(text = element_text(size = 20),axis.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 8))))
dev.off()

obj$Stage = gsub("mfin","0dpa",obj$Stage)
obj$Stage<-factor(obj$Stage,levels=c("0dpa","1dpa","2dpa","4dpa","6dpa"))
pdf(paste0("WNN_integrated_wnnUMAP_plot_byCelltypeUMAP_split_byStage.pdf"),width = 8, height = 20)
print(DimPlot(obj, cols = celltype_color_final,reduction = "wnn.umap", split.by = "Stage", ncol=1, label = F, label.size = 10, pt.size =0.1, repel = TRUE)+
theme(text = element_text(size = 20),axis.text = element_text(size = 20))+
guides(color = guide_legend(override.aes = list(size = 8))))
dev.off()

# Fig 1C
DefaultAssay(obj)="RNA"
genes=c(
    "krt1-19d","cldn23a", #superficial epithelial
    "cldna","tp63",#"krt94",#Intermediate epithelial
    "cldni","fras1",#Basal epithelial
    "cldnh",#mucous and metaphocyte
    "agr2",#agr2EM1
    "pvalb8",#pvalb8 EM2
    "grn2", "cdh28",#"mpeg1.1", #metaphocyte
    "ptprc", #immune cells
    "mpeg1.1", #"spi1b", #macrophage,"lgmn","cxcr3.2",
    "tnfrsf9b", #activated T cells
    "msx1b", "tnfaip6",#"twist1a", #Mesenchymal "msx3",
    "prrx1b", #fibroblast mesenchymal
    "cdh11", #osteoblast mesenchymal
    "kdrl","fli1a",#,#Endothelial,"plvapb"
    "trpm1b","gpr143", #pan pigment
    "si:ch211-213d14.1", "gng13a", #Tuft-like
    "fndc7a","stm" #Neuromast-like    
    )

pdf(paste0("Allcelltype_marker_gene_Vlnplot.pdf"),width=12,height=12)
VlnPlot(obj, fill.by = "ident", cols = celltype_color_final, features = genes, stack = TRUE, flip = TRUE)
dev.off()

# Fig 1D
# ATAC coverage plot in maker gene promoter region
library(ggplot2)
library(cowplot)
library(patchwork)

DefaultAssay(obj)="peaks"
mk.genes = c(
  "krt1-19d", "cldna", "cldni", "agr2", "pvalb8", "grn2",
  "ptprc", "tnfaip6", "kdrl", "trpm1b", "gng13a", "fndc7a", "cldnh")

tss =GetTSSPositions(Annotation(obj), biotypes = "protein_coding")
mk.tss = tss[which(tss$gene_name %in% mk.genes)]

mk.tss.flank = Extend(mk.tss, upstream = 1000, downstream = 1000, from.midpoint = T)
mk.tss.flank = mk.tss.flank[match(mk.genes,mk.tss.flank$gene_name)]
regions = paste(
  data.frame(mk.tss.flank)$seqnames,
  data.frame(mk.tss.flank)$start,
  data.frame(mk.tss.flank)$end,
  sep="-")


p=list()
for (i in 1:length(regions)){
    p[[i]] <-CoveragePlot(
      object = obj,
      assay = "peaks",
      region = regions[i],
      annotation = T,
      extend.upstream = 0,
      extend.downstream =0,
      peaks = FALSE
    ) & scale_fill_manual(values = celltype_color_final) & theme(text = element_text(size = 20))
}
names(p)=mk.tss.flank$gene_name

p.select<-p
library(cowplot)
p.all.select <- do.call(plot_grid, c(p.select, list(ncol = length(p.select))))


pdf("Allcelltype_markergene_coveragePlot_1kb.pdf",width=60, height=10)
p.all.select
dev.off()
