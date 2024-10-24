## Fig6A,6B
```r
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(paletteer)
library(tidyverse)

obj=readRDS("MES_seurat_obj_SCT_RNA_integrated_only.rds")
my_cols=paletteer_d("ggthemes::Classic_Cyclic") #MES WNN

pdf(paste0("MES_UMAP_plot.pdf"),width = 15, height = 15)
DimPlot(obj, cols = my_cols,reduction = "umap", label = F, label.size = 10, pt.size =3, repel = TRUE) + 
ggtitle("RNA")+
theme(text = element_text(size = 40),axis.text = element_text(size = 30))+
guides(color = guide_legend(override.aes = list(size = 15)))
dev.off()

pdf(paste0("MES_UMAP_plot_SplitbyStage.pdf"),width = 40, height = 13)
DimPlot(obj, cols = my_cols,reduction = "umap", split.by="Stage", label = F, label.size = 10, pt.size =3, repel = TRUE) + 
ggtitle("RNA")+
theme(text = element_text(size = 40),axis.text = element_text(size = 30))+
guides(color = guide_legend(override.aes = list(size = 15)))
dev.off()
```


## Fig6C
```r
genes_fig6 = c(
    "cdh11","runx2a","mmp9","col10a1a","sp7","ifitm5","bglap","lmx1ba","pthlha","evx1",
    "and2","hapln1a","tph1b","robo3","mmp13a.1","loxa","fhl1a","rgs5a","mfap5","abi3bpb")

library(scCustomize)
pdf("Fig6_marker_gene_feature_plot_selected.pdf",width=22,height=15)
FeaturePlot_scCustom(
    seurat_object = obj, reduction="umap",combine = TRUE,
    colors_use = viridis_light_high,na_color = "#440154",
    features = genes_fig6) +patchwork::plot_layout(ncol = 5, nrow = 4) & 
theme(plot.title = element_text(size = 20))& NoAxes() & theme(plot.title = element_text(vjust = -10))
dev.off()
```

## Fig6D
```r
links = readRDS("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/WNN_celltype_subset_object/p2g/Mesenchymal_peak_gene_links_100kb.rds")


obj=readRDS("MES_C3rm_Dim1-15_res0.8_cellreordered_seurat_obj.rds")
DefaultAssay(obj)="peaks"
Links(obj)=links

pdf("MES_lineage_marker_gene_peak_links.pdf")
print(
    CoveragePlot(
  object = obj,
  region = "cdh11",
  features = gene,
  expression.assay = "RNA",
  #idents = idents.plot,
  extend.upstream = 50000,
  extend.downstream = 50000,
  heights = 1
) & scale_fill_manual(values = my_cols)& theme(text = element_text(size = 12))
)
dev.off()
```

## Fig6E
```r
dir.create("centrality_degree_comparison")
for (i in c(1:13)){
    for (j in (c(1:13))){
    cp1 = paste0("MES_",i-1)
    cp2 = paste0("MES_",j-1)
    if (i!=j){
        reshaped_df.filter = reshaped_df[,c("X",cp1,cp2)]
        top_tfs1 = reshaped_df.filter[order(reshaped_df.filter[,cp1], decreasing = TRUE),][c(1:10),]
        top_tfs2 = reshaped_df.filter[order(reshaped_df.filter[,cp2], decreasing = TRUE),][c(1:10),]
        top_tfs = unique(rbind(top_tfs1,top_tfs2))
        pdf(paste("centrality_degree_comparison/MES_top10_eigenvector_centrality",cp1,"vs",cp2,"scatterplot.pdf",sep="_"),width=4,height=4)
        print(
        ggplot(reshaped_df.filter, aes_string(x = cp1, y = cp2)) +
        geom_point() + 
        theme_bw() + 
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        ggrepel::geom_label_repel(data = top_tfs, 
                        aes(label = X), 
                        max.overlaps = Inf,
                        box.padding = 0.35, 
                        point.padding = 0.5,
                        segment.color = 'grey50', color = "darkred") +
                        labs(
                         title = paste0(""),
                         y = paste(cp2,"eigenvector centrality",sep=" "),
                         x = paste(cp1,"eigenvector centrality",sep=" "))
         )
        dev.off()
        }
    }
}
```

## Fig6F
```r
genes = c("mef2cb","lmx1ba")
pdf("Fig6_mef2cb_lmx1ba_feature_plot.pdf",width=8,height=4)
FeaturePlot_scCustom(
    seurat_object = obj, reduction="umap",combine = TRUE,
    colors_use = viridis_light_high,na_color = "#440154",
    features = genes)+patchwork::plot_layout(ncol = 2) & 
theme(plot.title = element_text(size = 20))& NoAxes() & theme(plot.title = element_text(vjust = -10))
dev.off()
```
