## Figure 3A
```r
# Assign colors to cell types
library(ggplot2)
library(reshape2)

celltype_color <- c("#8B2DB2","#CE6DBD","#9C9EDE","#F28E2BFF", "#59A14FFF","#B6992DFF", "#FF9D9AFF",
  "#86BCB6FF","#E15759FF")
celltypes <- c("Superficial_Epithelial","Intermediate_Epithelial","Basal_Epithelial","Epidermal_Mucous","Hematopoietic",
  "Mesenchymal","Pigment","Endothelial", "Metaphocyte")
 names(celltype_color) <- celltypes

# create your comparison order file
data_order <- read.table("DAR_order.txt", header=FALSE) ## generated with stage_deg_dar.r
data_order.df <- as.data.frame(data_order, row.names =rownames(data_order))

for (i in celltypes){
    fileNM <- paste0("WNN_celltype_",i,"_DARs_stage_comparison_stats.txt")
    data <- read.table(fileNM)
    rownames(data)<-data$V1
    data_order.df[i] <- data[data_order$V1,"V2"]
  }

long <- melt(data_order.df, id.vars = c("V1"), variable.name = "CellType")
names(long) <- c("Stage_Comparison", "CellType", "Number")
long[is.na(long)] <- 0
long = long %>% separate(Stage_Comparison, c('Category', 'Comparison','Type'))
long<-dcast(long, ...~Type)
write.csv(long, "WNNcelltype_DARs_stage_comparison_summary_for_back2back_barplot.csv", row.names=FALSE) 
df <-read.csv("WNNcelltype_DARs_stage_comparison_summary_for_back2back_barplot.csv") 

# Create the plot
# Set the width of the bars
bar_width <- 1
long %>% mutate(CellType = factor(CellType, levels = celltypes)) %>%
  ggplot(aes(x = CellType)) +
  geom_bar(aes(y = open, fill = CellType),colour="black", stat = "identity", width = bar_width) +
  geom_bar(aes(y = -close, fill = CellType),colour="black", stat = "identity", width = bar_width) +
  scale_y_continuous(name="DAR counts", limits=c(-4000 , 4000),labels = abs) +
  #scale_fill_manual(name = "celltype", values = c("Up-regulated" = up_color, "Down-regulated" = down_color)) +
  scale_fill_manual(name = "celltype", values = celltype_color) +
  scale_color_manual(name = "celltype", values = celltype_color) +
  #coord_flip() +
  geom_vline(xintercept = 0) +
  labs(title = "DAR by stage comparison", y = "Counts") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line.x = element_line(colour = "black"))+
  theme(panel.spacing = unit(.05, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 1), strip.background = element_rect(color = "black",size = 1),strip.text.x = element_text(colour = "black",face="bold"))+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ facet_grid(cols = vars(Comparison)) +
  theme(axis.text = element_text(face="bold",color="black"),axis.title = element_text(face="bold",color="black"))
ggsave("DAR_back2back_barplot.pdf",width=10,height=4)
```

## Figure 3B
```r
library(reshape2)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(paletteer)
library(UpSetR)
library(ComplexHeatmap)
set.seed(1234)

celltype_color<-c("#8B2DB2","#CE6DBD","#9C9EDE","#F28E2BFF", "#59A14FFF","#B6992DFF")
celltypes <- c("Superficial_Epithelial","Intermediate_Epithelial","Basal_Epithelial","Epidermal_Mucous","Hematopoietic","Mesenchymal")
names(celltype_color)<-celltypes

combined_rrrs_list = list()
combined_rrrs_list = lapply(celltypes,function(X) readRDS(paste0(X,file="_RRRs.rds")))
names(combined_rrrs_list)=celltypes
peak_str.list = lapply(combined_rrrs_list,function(X) {
  peaks=as.character(X)
  peaks =gsub(":","-",peaks)
  return(unique(peaks))
  })

lt = peak_str.list
lt <- lapply(lt, function(x) sort(x)) #have to sort before making matrix

m = make_comb_mat(lt)
m=m[comb_size(m) >= 20]
m_list[[celltype]]=m

pdf("major_celltype_RRRs_upsetPlot.pdf",height=4,width=15)
ss = set_size(m)
cs = comb_size(m)
ht = UpSet(m, 
    set_order = names(ss),
    comb_order = order(comb_degree(m), -cs),
    top_annotation = HeatmapAnnotation(
        "Intersections" = anno_barplot(cs, 
            ylim = c(0, max(cs)*1.1),
            border = FALSE, 
            gp = gpar(fill = "black",fontface="bold",fontsize = 8), 
            height = unit(4, "cm")
        ), 
        annotation_name_side = "left", 
        annotation_name_rot = 90),
    left_annotation = rowAnnotation(
        "Gene Counts" = anno_barplot(-ss, 
            baseline = 0,
            axis_param = list(
                at = c(0, -2500, -4000, -6500,-8000),
                labels = c(0, 2000, 4000, 6000, 8000),
                labels_rot = 45),
            border = F, 
            gp = gpar(fill = celltype_color,fontsize = 10), 
            width = unit(3.5, "cm")
        ),
        set_name = anno_text(set_name(m), 
            location = 0.5, 
            just = "center",
            width = max_text_width(set_name(m)) + unit(4, "mm"))
    ), 
    right_annotation = NULL,
    show_row_names = FALSE)
  ht = draw(ht)
  od = column_order(ht)
  decorate_annotation("Intersections", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
        default.units = "native", just = c("left", "bottom"), 
        gp = gpar(fontsize = 10, col = "#404040", hjust = 0.5), rot = 45)
        
})
dev.off()
```


## Figure 3C
```bash
#!/bin/bash
## 1. deeptools: RRR region signal across all stages in each major cell type
# Example:
## ${celltype}_atac_bigwig_list.txt saves paths to all bigwig files
## ${celltype}_openRRRs.bed is the opening RRRs of ${celltype}
celltype="Mesenchymal"
computeMatrix reference-point \
--referencePoint center \
-p 8 -bs 10 \
-S `cat ${celltype}_atac_bigwig_list.txt` \
-R ${celltype}_openRRRs.bed \
-b 5000 \
-a 5000 \
--skipZeros -o ${celltype}_openRRRs_bigwig_matrix.mat.gz \
--outFileNameMatrix ${celltype}_openRRRs_bigwig_matrix_scaled.tab \
--outFileSortedRegions ${celltype}_openRRRs_bigwig_matrix_sortregion.bed \
> computeMatrix_${celltype}_openRRRs_bigwig_matrix.log 2>&1 &

plotHeatmap \
-m ${celltype}_openRRRs_bigwig_matrix.mat.gz \
--plotTitle "${celltype} openRRRs" \
--regionsLabel "${celltype} openRRRs" \
--samplesLabel "0 dpa" "1 dpa" "2 dpa" "4 dpa" "6 dpa" \
--colorMap "PRGn" "BrBG" \
--perGroup \
--legendLocation upper-right \
--heatmapHeight 10 \
--heatmapWidth 5 \
--xAxisLabel "distance" \
--startLabel "start" \
--endLabel "end" \
--plotFileFormat pdf \
--verbose \
-out ${celltype}_openRRRs_byStage_signal.pdf
```

```r
## 2. RRR nearest gene expression change box plot
library(tidyverse)
library(ggprism)
library(rstatix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(data.table)

rna.log2fc =readRDS("major_celltype_all_stage_0vsx_avglog2FC_gene_expression.rds")
atac.log2fc =readRDS("major_celltype_all_stage_0vsx_avglog2FC_peak_expression.rds")
peakAnnoList =readRDS("All_major_celltype_RRRs_annotation_chipseeker.rds")

# change the log2 fold change direction, above 0 means higher feature expression after injury
rna.log2fc = -rna.log2fc 
atac.log2fc = -atac.log2fc

celltypes<- c("Superficial_Epithelial","Intermediate_Epithelial","Basal_Epithelial","Epidermal_Mucous","Hematopoietic","Mesenchymal")
plot.list=list()
for (ct in celltypes){
# subset rna.log2fc, atac.log2fc
columns = paste0(ct,c("_0vs1_log2FC","_0vs2_log2FC","_0vs4_log2FC","_0vs6_log2FC"))

rna.log2fc.ct = rna.log2fc[,columns]
atac.log2fc.ct = atac.log2fc[,columns]

colnames(rna.log2fc.ct) = paste(colnames(rna.log2fc.ct),"RNA",sep="_")
colnames(atac.log2fc.ct) = paste(colnames(atac.log2fc.ct),"ATAC",sep="_")

#ct RRR
rrr.gr =readRDS(paste0(ct,"_RRRs.rds"))

# ct RRR annotation (got from ChIPseeker)
peakAnno = peakAnnoList[[ct]]@anno #with gene symbol with RRR category already
peakAnno$peak = paste(data.frame(peakAnno)[,1],data.frame(peakAnno)[,2],data.frame(peakAnno)[,3],sep="-")
peakAnno.df = data.frame(peakAnno)

#merge atac peak log2fc
peakAnno.df = merge(x= peakAnno.df, y = atac.log2fc.ct,by.x = "peak", by.y = "row.names", all.x = TRUE)
#merge rna gene log2fc
peakAnno.df = merge(x= peakAnno.df, y = rna.log2fc.ct,by.x = "geneSymbol", by.y = "row.names", all.x = TRUE)

# separate RRRs to promoter(+/- 1kb TSS), proximal(+/- 10kb TSS, except promoter) and distal (+/- 100kb TSS, except promoter and proximal)
# change complex intron annotation to "intron" only
peakAnno.df$annotation[grepl("^Intron", peakAnno.df$annotation)] ="Intron" 
peakAnno.df$annotation[grepl("^Downstream", peakAnno.df$annotation)] ="Intergenic" 

# add promoter(+/-1kb), proximal intronic/intergenic(+/-10kb), distal annotation
peakAnno.df$peak_category <- ifelse(
  peakAnno.df$annotation == "Promoter", "Promoter", #I set promoter as +/-1kb in ChiPseeker analysis
  ifelse(
    (peakAnno.df$annotation %in% c("Distal Intergenic", "Intergenic", "Intron")) & abs(peakAnno.df$distanceToTSS) <= 10000, "Proximal",
    ifelse(
      (peakAnno.df$annotation %in% c("Distal Intergenic","Intergenic", "Intron")) & abs(peakAnno.df$distanceToTSS) > 10000 & abs(peakAnno.df$distanceToTSS) <= 100000, "Distal",
      NA
    )
  )
)

# filter out the openRRR from 0vsX dpa comparisons, distanceToTSS 100kb, remove NA
peakAnno.df.openrrr = peakAnno.df %>% filter(name ==paste(ct,"openRRR",sep="_")) 
peakAnno.df.openrrr = peakAnno.df.openrrr %>% filter(abs(distanceToTSS)<=100000) #~300 filtered out in mes

peakAnno.df.openrrr$peak_category = factor(peakAnno.df.openrrr$peak_category,levels = c("Promoter","Proximal","Distal"))

#extract promoter, proximal 
nearest.gene.exp = peakAnno.df.openrrr[,c(1,23:27)]
nearest.gene.exp.promoter = nearest.gene.exp %>% filter(peak_category =="Promoter") %>% unique() 
promoter.gene.name = unique(nearest.gene.exp.promoter$geneSymbol)
nearest.gene.exp.proximal = nearest.gene.exp %>% filter(peak_category =="Proximal") %>% filter(!geneSymbol %in% promoter.gene.name) %>% unique() 
nearest.gene.exp.distal = nearest.gene.exp %>% filter(peak_category =="Distal") %>% filter(!geneSymbol %in% promoter.gene.name) %>% unique() 

nearest.gene.exp.promoter = nearest.gene.exp.promoter %>% na.omit()
nearest.gene.exp.proximal = nearest.gene.exp.proximal %>% na.omit()
nearest.gene.exp.distal = nearest.gene.exp.distal %>% na.omit()

df = rbind(nearest.gene.exp.promoter, nearest.gene.exp.proximal,nearest.gene.exp.distal)

df$peak_category = factor(df$peak_category, levels = c("Promoter", "Proximal", "Distal"))

colnames(df) = c("geneSymbol","1dpa", "2dpa", "4dpa", "6dpa","peak_category")

df_long <- as.data.frame(df) %>%
  pivot_longer(cols = c("1dpa", "2dpa", "4dpa", "6dpa"), names_to = "Time", values_to = "Log2FC_Expression") %>%
  mutate(Time = factor(Time, levels = c("1dpa", "2dpa", "4dpa", "6dpa")))

stat.test <- df_long %>%
  group_by(peak_category) %>%
  pairwise_wilcox_test(Log2FC_Expression ~ Time, p.adjust.method = "BH") %>%
  add_xy_position()

bxp <- ggboxplot(
    df_long, x = "Time", y = "Log2FC_Expression",facet.by = "peak_category",outlier.shape = NA,
    fill = "Time",
    palette = c("#8F429B", "#C2A5C9", "#F3F2F3", "#8CC88D"))+
    geom_hline(yintercept = 0, linetype="dashed", color = "darkgrey") + 
    labs(title = paste0(ct))+theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(limits = c(-2, 4)) 
 
stat.test = stat.test %>% mutate(cp = paste(group1,group2,sep="vs")) 
stat.test.filtered = stat.test %>% filter(cp %in% c("1dpavs2dpa","2dpavs4dpa", "4dpavs6dpa")) %>% add_xy_position() %>% mutate(y.position = y.position-4)
#stat.test.filtered = stat.test %>% filter(cp %in% c("1dpavs2dpa","2dpavs4dpa", "4dpavs6dpa","1dpavs4dpa","1dpavs6dpa")) %>% add_xy_position()

bxp = bxp + stat_pvalue_manual(stat.test.filtered, label = "p.adj", tip.length = 0.01) 
plot.list[[ct]] = bxp

}

combined_plot = cowplot::plot_grid(plotlist = plot.list, ncol = 3)
pdf("RRR_nearest_gene_expresson_change_boxplot.pdf", width = 15, height = 8)
print(combined_plot)
dev.off()
```

## Figure 3D,3E
```bash
# RRRs GREAT analysis 
#!/bin/bash
module load kentUCSC

for i in Superficial_Epithelial Intermediate_Epithelial Basal_Epithelial Epidermal_Mucous Hematopoietic Mesenchymal; \
do file_name=$i"_RRRs.bed" new_file=$i"_RRRs_danRer7.bed" umapped=$i"_RRRs_danRer11to7_unmapped.bed"; echo $file_name; \
liftOver $file_name danRer11ToDanRer7.over.chain.gz $new_file $umapped; done
# Output regions in danRer7 version are used as the input for GREAT analysis
```
