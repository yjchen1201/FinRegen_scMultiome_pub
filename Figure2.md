## Fig2A
```r
library(reshape2)
library(ggplot2)
celltype_color <- c("#8B2DB2","#CE6DBD","#9C9EDE","#F28E2BFF", "#59A14FFF","#B6992DFF", "#FF9D9AFF",
  "#86BCB6FF","#E15759FF")
celltypes <- c("Superficial_Epithelial","Intermediate_Epithelial","Basal_Epithelial","Epidermal_Mucous","Hematopoietic",
  "Mesenchymal","Pigment","Endothelial", "Metaphocyte")
 names(celltype_color) <- celltypes

# create your comparison order file
data_order <- read.table("DEG_order.txt", header=FALSE) ## Before this, I wrote an order file with Sublime Text
data_order.df <- as.data.frame(data_order, row.names =rownames(data_order))

# Load DEG stats list generated using stage_deg_dar.r
for (i in celltypes){
    fileNM <- paste0("WNN_celltype_",i,"_DEGs_stage_comparison_usingRNAassay_stats.txt")
    data <- read.table(fileNM)
    rownames(data)<-data$V1
    data_order.df[i] <- data[data_order$V1,"V2"]
  }


long <- melt(data_order.df, id.vars = c("V1"), variable.name = "CellType")
names(long) <- c("Stage_Comparison", "CellType", "Number")
long[is.na(long)] <- 0
long = long %>% separate(Stage_Comparison, c('Category', 'Comparison','Type'))
long<-dcast(long, ...~Type)
write.csv(long, "WNNcelltype_DEGs_stage_comparison_summary_for_back2back_barplot.csv", row.names=FALSE) 
df <-read.csv("WNNcelltype_DEGs_stage_comparison_summary_for_back2back_barplot.csv") 

# Create the plot
# Set the width of the bars
bar_width <- 1
long %>% mutate(CellType = factor(CellType, levels = celltypes)) %>%
  ggplot(aes(x = CellType)) +
  geom_bar(aes(y = up, fill = CellType),colour="black", stat = "identity", width = bar_width) +
  geom_bar(aes(y = -down, fill = CellType),colour="black", stat = "identity", width = bar_width) +
  scale_y_continuous(name="DEG counts", limits=c(-4000 , 4000),labels = abs) +
  scale_fill_manual(name = "celltype", values = celltype_color) +
  scale_color_manual(name = "celltype", values = celltype_color) +
  geom_vline(xintercept = 0) +
  labs(title = "DEG by stage comparison", y = "Counts") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line.x = element_line(colour = "black"))+
  theme(panel.spacing = unit(.05, "lines"),panel.border = element_rect(color = "black", fill = NA, size = 1), strip.background = element_rect(color = "black",size = 1),strip.text.x = element_text(colour = "black",face="bold"))+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ facet_grid(cols = vars(Comparison)) +
  theme(axis.text = element_text(face="bold",color="black"),axis.title = element_text(face="bold",color="black"))
ggsave("DEG_back2back_barplot.pdf",width=10,height=4)
```

## Fig 2B
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

ct_deg_list =list()
for (ct in celltypes){
    print(paste0("deg for ",ct))

    fileNM <- paste("WNN_celltype_",ct,"_DEGs_stage_comparison_usingRNAassay.xlsx", sep="") #generated using stage_deg_dar.r
   
    deg_list<-list()
    deg_sheet_names = c("deg_0vs1","deg_0vs2","deg_0vs4","deg_0vs6","deg_1vs2","deg_1vs4","deg_1vs6","deg_2vs4","deg_2vs6","deg_4vs6")
    tryCatch(deg_list <- lapply(deg_sheet_names, function(X) read.xlsx(fileNM, sheet = X)), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    names(deg_list) = deg_sheet_names
    deg_list <-lapply(deg_list, function(X) X=X$gene)
    degs = unique(unlist(deg_list))
    ct_deg_list[[ct]]=degs
    rm(deg_list)
    
}

lt = ct_deg_list
lt <- lapply(lt, function(x) sort(x)) #have to sort before making matrix

m = make_comb_mat(lt)
m=m[comb_size(m) >= 10]

pdf("major_celltype_bystage_DEGs_upsetPlot.pdf",height=4,width=15)
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
                at = c(0, -1000, -2000, -3000, -4000,-5000),
                labels = c(0, 1000, 2000, 3000, 4000,5000),
                labels_rot = 0),
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
        gp = gpar(fontsize = 10, col = "#404040", hjust = 0.5), rot = 0)
        
})
dev.off()
```

## Fig2C
```r
rna=readRDS("/integrated_dietseurat_rna_only.rds")
Idents(rna)="celltype_wnn_merge"
celltypes<-c("Superficial_Epithelial",
    "Intermediate_Epithelial",  
    "Basal_Epithelial",
    "Epidermal_Mucous",
    "Hematopoietic",
    "Mesenchymal")

rna.subset = subset(rna, idents = celltypes)
rna.sub = rna
rna.subset$ct_stage = paste(rna.subset$celltype_wnn_merge,rna.subset$Stage,sep="_")

celltype_colors<-c("#8B2DB2","#CE6DBD","#9C9EDE","#F28E2B", "#59A14F","#B6992D")

stages<-c("0dpa", "1dpa","2dpa","4dpa","6dpa")
combinations_df <- expand.grid(stages,celltypes)
ct_stage_combinations <- with(combinations_df, paste(Var2, Var1, sep = "_"))

rna.subset$Stage = gsub("mfin","0dpa",rna.subset$Stage)
rna.subset$Stage = factor(rna.subset$Stage,levels=stages)
rna.subset$ct_stage = gsub("mfin","0dpa",rna.subset$ct_stage)
rna.subset$ct_stage = factor(rna.subset$ct_stage,levels=ct_stage_combinations)

DefaultAssay(rna)="RNA"
rna=NormalizeData(rna)
celltype_wnn = levels(Idents(rna))

fc_all_ct_all_stage.list =list()
for (ct in levels(Idents(rna))){
  print(paste(ct))
  fc_1dpa = FoldChange(
  rna,
  ident.1 = "mfin",
  ident.2 = "1dpa",
  group.by ='Stage' ,
  subset.ident = ct,
  assay = "RNA",
  slot = "data",
  #reduction = "wnn.umap",
  features = rownames(rna),
  #pseudocount.use = NULL,
  #mean.fxn = NULL,
  base = 2,
  fc.name = NULL
)

fc_2dpa = FoldChange(
  rna,
  ident.1 = "mfin",
  ident.2 = "2dpa",
  group.by ='Stage' ,
  subset.ident = ct,
  assay = "RNA",
  slot = "data",
  features = rownames(rna),
  base = 2,
  fc.name = NULL
)

fc_4dpa = FoldChange(
  rna,
  ident.1 = "mfin",
  ident.2 = "4dpa",
  group.by ='Stage' ,
  subset.ident = ct,
  assay = "RNA",
  slot = "data",
  features = rownames(rna),
  base = 2,
  fc.name = NULL
)

fc_6dpa = FoldChange(
  rna,
  ident.1 = "mfin",
  ident.2 = "6dpa",
  group.by ='Stage' ,
  subset.ident = ct,
  assay = "RNA",
  slot = "data",
  features = rownames(rna),
  base = 2,
  fc.name = NULL
)

fc_all_stage = cbind(fc_1dpa[,1],fc_2dpa[,1],fc_4dpa[,1],fc_6dpa[,1])
rownames(fc_all_stage)=rownames(rna)
colnames(fc_all_stage)=c("0vs1_log2FC", "0vs2_log2FC", "0vs4_log2FC", "0vs6_log2FC")
colnames(fc_all_stage) = paste(ct,colnames(fc_all_stage),sep="_")
fc_all_ct_all_stage.list[[ct]] = fc_all_stage
rm(fc_all_stage,fc_1dpa,fc_2dpa,fc_4dpa,fc_6dpa)
}

combined_major_ct_df <-do.call(cbind, fc_all_ct_all_stage.list[1:6])
write.csv(combined_major_ct_df,file="/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/DEGs_DARs/DEGs_stage_comparison_RNAassay_logFC0.25_minpct0.1_padj0.05/all_major_cellctype_all_stage_0vsx_log2FC_gene_expression.csv",row.names=T)


#Average Expression
rna$ct_stage=paste(Idents(rna),rna$Stage,sep="_")
gsub("mfin","0dpa",rna$ct_stage) ->rna$ct_stage

stages<-c("0dpa", "1dpa","2dpa","4dpa","6dpa")
celltype_wnn = levels(Idents(rna))
combinations_df <- expand.grid(stages,celltype_wnn)
combinations <- with(combinations_df, paste(Var2, Var1, sep = "_"))

rna$ct_stage = factor(rna$ct_stage,levels=combinations)
Idents(rna)="ct_stage"

DefaultAssay(rna)="RNA"
rna=NormalizeData(rna)
avg.exp <- AverageExpression(rna, assays = "RNA" , slot = "data")
avg.exp = log1p(avg.exp$RNA)

library(UpSetR)
library(ComplexHeatmap)
set.seed(1234)

celltypes<- c("Superficial_Epithelial","Intermediate_Epithelial","Basal_Epithelial","Epidermal_Mucous","Hematopoietic","Mesenchymal")
all.ct.deg = readRDS("major_celltype_all_bystage_DEGs.rds")


lt = all.ct.deg[celltypes]
lt <- lapply(lt, function(x) sort(x)) #have to sort before making matrix

m = make_comb_mat(lt)
m=m[comb_size(m) >= 10]

all.shared.deg=extract_comb(m,"111111");write.csv(all.shared.deg,file="major_celltype_intersectedDEGs.csv")


celltype_colors<-c("#8B2DB2","#CE6DBD","#9C9EDE","#F28E2BFF", "#59A14FFF","#B6992DFF")#, "#FF9D9AFF","#86BCB6FF","#E15759FF")
celltypes <- c("Superficial_Epithelial","Intermediate_Epithelial","Basal_Epithelial","Epidermal_Mucous","Hematopoietic","Mesenchymal")#,"Pigment","Endothelial", "Metaphocyte")
names(celltype_colors)<-celltypes

for (ct in celltypes){

    print(paste0("deg for ",ct))  
    rm(all_degs)
    all_degs = readRDS(paste0("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/DEGs_DARs/DEGs_stage_comparison_RNAassay_logFC0.25_minpct0.1_padj0.05/",ct,"_all_stage_comparison_DEG_genename_list.rds"))
    all_degs.unlist = unlist(all_degs) %>% unique() #5107

    fileNM <- paste("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/DEGs_DARs/DEGs_stage_comparison_RNAassay_logFC0.25_minpct0.1_padj0.05/WNN_celltype_",ct,"_DEGs_stage_comparison_usingRNAassay.xlsx", sep="")
   
    deg_list<-list()
    deg_sheet_names = c("deg_0vs1","deg_0vs2","deg_0vs4","deg_0vs6")#,"deg_1vs2","deg_1vs4","deg_1vs6","deg_2vs4","deg_2vs6","deg_4vs6")
    tryCatch(deg_list <- lapply(deg_sheet_names, function(X) read.xlsx(fileNM, sheet = X)), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    names(deg_list) = deg_sheet_names
    rm(df1,df2,df3,df4)
    df1 = deg_list[[1]][which(deg_list[[1]]$gene %in% all_degs.unlist),] %>% dplyr::select(gene,avg_log2FC) %>% dplyr::rename(`0vs1_log2FC` = `avg_log2FC`)
    df2 = deg_list[[2]][which(deg_list[[2]]$gene %in% all_degs.unlist),] %>% dplyr::select(gene,avg_log2FC) %>% dplyr::rename(`0vs2_log2FC` = `avg_log2FC`)
    df3 = deg_list[[3]][which(deg_list[[3]]$gene %in% all_degs.unlist),] %>% dplyr::select(gene,avg_log2FC) %>% dplyr::rename(`0vs4_log2FC` = `avg_log2FC`)
    df4 = deg_list[[4]][which(deg_list[[4]]$gene %in% all_degs.unlist),] %>% dplyr::select(gene,avg_log2FC) %>% dplyr::rename(`0vs6_log2FC` = `avg_log2FC`)
    df4 = deg_list[[4]][which(deg_list[[4]]$gene %in% all_degs.unlist),] %>% dplyr::select(gene,avg_log2FC) %>% dplyr::rename(`0vs6_log2FC` = `avg_log2FC`)

    rm(dfs); dfs <- list(df1, df2, df3, df4)

    # Use Reduce to merge all data frames
    rm(df_all); df_all <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), dfs) 
    df_all = df_all %>% column_to_rownames(var = "gene")
    df_all[is.na(df_all)] = 0  #change all NA to 0
    df_all = round(df_all, 3)*-1 #only preserve 3 decimals, and change the FC value by *-1
    colnames(df_all) = paste(ct,colnames(df_all),sep="_")
    write.csv(df_all,file=paste0(ct,"_all_degs_log2FC_compareto0dpa.csv"))
}


df.list=list()
df_combined =data.frame()
for (ct in celltypes){
    rm(df.ct); df.ct = read.csv(paste0("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/DEGs_DARs/DEGs_stage_comparison_RNAassay_logFC0.25_minpct0.1_padj0.05/",ct,"_all_degs_log2FC_compareto0dpa.csv"),row.names=1)
    #df.list[[ct]] = df.ct
    df_combined = merge(df_combined,df.ct,by="row.names", all = TRUE)
    df_combined  = df_combined  %>%  column_to_rownames(var = "Row.names")
    }

df_combined[is.na(df_combined)] = 0  #change all NA to 0

row_reordered.df = read.csv(file=paste0("major_celltype_intersectedDEGs_all_DEGs_log2fc_clustering_heatmap_genes_km",km,".csv"),row.names=1)
df_combined.intersected =-df_combined[row_reordered.df$gene,]

library(colorRamp2)
require(circlize)
require(ComplexHeatmap)
set.seed(1234)

km=2
col_fun = colorRamp2(c(-1, 0, 1), c("#009B9F", "#EFECEE", "#C75DAA")) #pink and green
ht <- Heatmap(as.matrix(df_combined.intersected), col = col_fun,name = "log2FC", cluster_rows=T,cluster_columns=T,show_row_names = FALSE,row_km = km,row_km_repeats = 100,row_gap = unit(2, "mm"))
pdf(paste0("allcell_intersected_degs_test_km",km,".pdf"),width=7,height=12)
ht = draw(ht)
row_order <- row_order(ht)
dev.off()

row_reordered <- row_order[as.character(c(1:km))]
row_reordered.df = stack(row_reordered)
row_reordered.df$module = paste("Module",row_reordered.df$ind,sep=" ")
row_reordered.df$gene = rownames(df_combined.intersected[row_reordered.df$value,])

row_split = rep(unique(row_reordered.df$module),unlist(lapply(row_reordered , length))) #generate label for rows to split each modules
row_split = factor(row_split,levels=unique(row_split))

require(circlize)
require(ComplexHeatmap)
col_fun = colorRamp2(c(-1, 0, 1), c("#009B9F", "#EFECEE", "#C75DAA")) #pink and green

ht_opt(
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
    legend_labels_gp = gpar(fontsize = 8, fontface = "bold"), 
    heatmap_column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    heatmap_column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    heatmap_row_title_gp = gpar(fontsize = 8, fontface = "bold")
)

celltype_col <- data.frame(celltype = rep(c("Superficial_Epithelial","Intermediate_Epithelial", "Basal_Epithelial", "Epidermal_Mucous","Hematopoietic","Mesenchymal"), c(4,4,4,4,4,4)))
celltype_col$celltype = factor(celltype_col$celltype,levels= unique(celltype_col$celltype))

celltype_colors<-c("#8B2DB2","#CE6DBD","#9C9EDE","#F28E2BFF", "#59A14FFF","#B6992DFF")#, "#FF9D9AFF","#86BCB6FF","#E15759FF")
celltypes <- c("Superficial_Epithelial","Intermediate_Epithelial","Basal_Epithelial","Epidermal_Mucous","Hematopoietic","Mesenchymal")#,"Pigment","Endothelial", "Metaphocyte")
names(celltype_colors)<-celltypes


ha1 = HeatmapAnnotation(
    `Cell types` = celltype_col$celltype,  # assuming celltype_col$celltype is a factor or character vector of cell types
    Stage = rep(c("1dpa", "2dpa", "4dpa", "6dpa"), 6),  # stage
    col = list(
        `Cell types` = celltype_colors,  # named vector for cell types
        Stage = c("1dpa" = "#8F7EE5", "2dpa" = "#6551CC", "4dpa" = "#422CB2", "6dpa" = "#260F99")
    ),
    na_col = "black", 
    annotation_name_side = "left",
    gp = gpar(col = "black")
)

genes_to_label <- c(
  "btg2","socs3a", "klf4", "stat3", #downregu M1
  "mdka","nrp1a","anxa1c", "inhbaa"#upregu regeneationGO
  )  # Replace with your genes
position_to_lable = match(genes_to_label,row_reordered.df$gene)
RowAnnotation = rowAnnotation(genes = anno_mark(at = position_to_lable, labels = genes_to_label, labels_gp = gpar(fontsize = 16)))

# plot with new row order and column clustering
pdf(paste0("major_celltype_intersectedDEGs_0vsXDEGs_log2fc_heatmap_km",km,"_new.pdf"),width=7,height=12)
ht = Heatmap(
  df_combined.intersected[row_reordered.df$gene,],
  col = col_fun,
  row_split =row_reordered.df$module, 
  row_gap = unit(3, "mm"), # between modules
  row_title_rot = 0,
  column_names_rot = 65,
  #column_names_gp = gpar(fontsize = 10),
  name = "avg_log2FC", 
  column_title = "Intersected regeneration reposonsive DEGs",
  cluster_rows=F,
  cluster_columns=F,
  show_row_names = FALSE,
  top_annotation = ha1,
  #right_annotation = ra,
  rect_gp = gpar(col = "white", lwd = 0),
  heatmap_legend_param = list(legend_height = unit(2, "cm")),
  row_names_gp = gpar(fontsize = 20),
  right_annotation = RowAnnotation) #for label specific genes

draw(
  ht, 
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
  )  

dev.off()
```

## Fig2D
```markdown
#### Extract shared gene list, use Metascape for GO term analysis
```

## Fig2E
```r
# check by stage avg expression scale normalization
rna.exp =read.csv("0_integrated_wnn_merged_celltype_rna_expression_bycelltype_bystage.csv",row.names=1)
all.ct.deg =readRDS("all_celltype_all_bystage_DEGs.rds")

colnames(rna.exp)=gsub("mfin","0dpa",colnames(rna.exp))

celltypes<-c("Superficial_Epithelial",
    "Intermediate_Epithelial",  
    "Basal_Epithelial",
    "Epidermal_Mucous",
    "Hematopoietic",
    "Mesenchymal")
stages <- c("0dpa","1dpa","2dpa","4dpa","6dpa")

#rm(combinations,combined_strings)
combinations <- expand.grid(stages = stages,celltypes = celltypes)
combined_strings <- paste(combinations$celltypes, combinations$stages, sep = "_")
rna.exp = rna.exp[,combined_strings]

library(ComplexHeatmap)
library(ggplot2)
set.seed(1234)
for (celltype in celltypes){
 
  rrg = all.ct.deg[[celltype]] 
  rna.exp.filtered = rna.exp[rrg, grep(paste0("^",celltype), colnames(rna.exp))]
  rna.exp.filtered[is.na(rna.exp.filtered)] = 0  #change all NA to 0
  rna.exp.filtered = t(scale(t(rna.exp.filtered)))
# #elbow plot
# set.seed(123)
# k.max <- 20
# data <- rna.exp.filtered
# wss <- sapply(1:k.max, 
#               function(k){kmeans(data, k, nstart=20,iter.max = 15 )$tot.withinss})
# wss
# plot(1:k.max, wss,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")
  km=3
  #pdf(paste0(celltype,"_0vsX_allDEGs_heatmap_km",km,".pdf"),width=7,height=7)
  ht <-Heatmap(rna.exp.filtered, 
          name = "Z-score",
          cluster_rows = TRUE, 
          cluster_columns = F,
          clustering_distance_rows = "euclidean",
          clustering_distance_columns = "euclidean",
          clustering_method_rows = "ward.D2",
          clustering_method_columns = "ward.D2",
          row_km = km,row_km_repeats = 100,row_gap = unit(2, "mm"),
          show_row_names = FALSE,
          show_column_names = TRUE)
  ht = draw(ht)
  row_order <- row_order(ht)
  dev.off()

  row_reordered <- row_order[as.character(c(1:km))]
  row_reordered.df = stack(row_reordered)
  row_reordered.df$module = paste("Module",row_reordered.df$ind,sep=" ")
  row_reordered.df$gene = rownames(rna.exp.filtered[row_reordered.df$value,])

  write.csv(row_reordered.df,file=paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,".csv"))

  row_split = rep(unique(row_reordered.df$module),unlist(lapply(row_reordered , length))) #generate label for rows to split each modules
  row_split = factor(row_split,levels=unique(row_split))

  library(circlize)
  colnames(rna.exp.filtered) <- c("0 dpa","1 dpa","2 dpa","4 dpa","6 dpa")
  pdf(paste0(celltype,"_allStage_allDEGs_heatmap_km",km,".pdf"),width=4,height=4)
  ht = Heatmap(
    rna.exp.filtered[row_reordered.df$value,],
    #col = colorRamp2(c(-2, 0, 2), c("#798233", "#EFECEE", "#D66982")),
    row_split =row_split, 
    row_gap = unit(1, "mm"), # between modules
    column_names_rot = 0,
    show_column_names = TRUE,
    column_names_centered = TRUE,
    #column_names_gp = gpar(fontsize = 10),
    name = "Z-score", 
    column_title = celltype,
    row_title =NULL,
    show_row_names = FALSE,
    cluster_rows=F,
    cluster_columns=F,
    #top_annotation = ,
    left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill="white"),
          labels = c("Module 1", "Module 2", "Module 3"),#,"Module 4","Module 5"), 
          labels_gp = gpar(col = "black", fontsize = 10))))

  draw(ht)  
  dev.off()

  df_long = read.csv(paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,".csv"),row.names=1)
  require(clusterProfiler)
  require(org.Dr.eg.db)
  gene.df <- bitr(df_long$gene, fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Dr.eg.db) #  6.77% fail to map...

  df_long$ENTREZID = gene.df$ENTREZID[match(df_long$gene, gene.df$SYMBOL)]
  # df_long = df_long %>% filter(ENTREZID!="NA")


  df_long = df_long %>% filter(ENTREZID!="NA")
  formula_res <- compareCluster(ENTREZID~module, data=df_long, fun="enrichGO",ont= "BP", OrgDb='org.Dr.eg.db', pvalueCutoff=0.05)
  #formula_res <- compareCluster(ENTREZID~module, data=df_long, fun="enrichKEGG", organism="dre",pvalueCutoff=0.05)

  #enrichKEGG
  simp_ego <- simplify(formula_res, cutoff=0.7, by="p.adjust", select_fun=min)
  simp_ego <- setReadable(simp_ego, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
  simp_ego <- simp_ego %>% filter(Count>=10)
  saveRDS(simp_ego,file=paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,"_moduleGOterm.rds"))

  pdf(paste0(celltype,"_0vsX_allDEGs_heatmap_genes_km",km,"_moduleGOenrich.pdf"), width=10,height=6)
  print(
    dotplot(simp_ego,showCategory=10)+
              theme(axis.text.x = element_text(angle =45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 10),
              title = element_text(size = 10))
  #             +theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
  #             axis.text.y = element_text(size = 14),
  #             axis.title = element_text(size = 16),
  #             title = element_text(size = 18))
  )
  dev.off()


dotplot(go,showCategory=5)+
              theme(axis.text.x = element_text(angle =45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 10),
              title = element_text(size = 10))
  #             +theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
  #             axis.text.y = element_text(size = 14),
  #             axis.title = element_text(size = 16),
  #             title = element_text(size = 18))
  
  dev.off()

  #write.csv(row_reordered.df,file=paste0("major_celltype_intersectedDEGs_log2fc_clustering_heatmap_genes_km",km,".csv"))
  df_long = read.csv(paste0("../",celltype,"_allStage_allDEGs_heatmap_genes_km",km,".csv"),row.names=1)
  require(clusterProfiler)
  require(org.Dr.eg.db)
  gene.df <- bitr(df_long$gene, fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org.Dr.eg.db) #  6.77% fail to map...

  df_long$ENTREZID = gene.df$ENTREZID[match(df_long$gene, gene.df$SYMBOL)]
  # df_long = df_long %>% filter(ENTREZID!="NA")


  df_long = df_long %>% filter(ENTREZID!="NA")
  formula_res <- compareCluster(ENTREZID~module, data=df_long, fun="enrichGO",ont= "BP", OrgDb='org.Dr.eg.db', pvalueCutoff=0.05)
  #formula_res <- compareCluster(ENTREZID~module, data=df_long, fun="enrichKEGG", organism="dre",pvalueCutoff=0.05)

  #enrichKEGG
  simp_ego <- simplify(formula_res, cutoff=0.7, by="p.adjust", select_fun=min)
  simp_ego <- setReadable(simp_ego, OrgDb = org.Dr.eg.db, keyType="ENTREZID")
  simp_ego <- simp_ego %>% filter(Count>=10)
  saveRDS(simp_ego,file=paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,"_moduleGOterm.rds"))

simp_ego = readRDS(paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,"_moduleGOterm.rds"))
write.csv(data.frame(simp_ego),file=paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,"_moduleGOterm.csv"),row.names=F)



simp_ego@compareClusterResult<- simp_ego@compareClusterResult %>%
  group_by(module) %>%
  nest() %>%
  mutate(data = map(data, ~arrange(.x, desc(Count)))) %>%unnest(data)
# for(celltype in celltypes){
# simp_ego = readRDS(paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,"_moduleGOterm.rds"))
  pdf(paste0(celltype,"_allStage_allDEGs_heatmap_genes_km",km,"_moduleGOenrich.pdf"), width=6,height=6)
  print(
    dotplot(simp_ego,showCategory=10)+
              theme(axis.text.x = element_text(angle =45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 8),
              axis.title = element_text(size = 8),
              title = element_text(size = 8))
  )
  dev.off()
}
```

## Fig 2F, 2G
```r
library(clusterProfiler)
library(org.Dr.eg.db)
library(ggplot2)
require(viridis)
library(scCustomize)
library(data.table)

mes =readRDS("MES_subset_seurat_object.rds")
celltype="Mesenchymal"
DefaultAssay(mes)="RNA"

go = readRDS(paste0("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/plot/for_pub/Fig2/major_celltype_allStage_DEGs_heatmap_modules/",celltype,"_allStage_allDEGs_heatmap_genes_km3_moduleGOterm.rds"))
km=3
pdf(paste0(celltype,"_allstageDEGs_heatmap_genes_km",km,"_moduleGOenrich.pdf"), width=6,height=6)
  print(
    dotplot(go,showCategory=10)+
              theme(axis.text.x = element_text(angle =45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title = element_text(size = 10),
              title = element_text(size = 10))
  )
  dev.off()

key.word = "extracellular matrix organization" 
df = data.frame(go@compareClusterResult)
# temp.df=df[df$Description %like% key.word,] #data.table
# temp.df <- temp.df %>%separate_rows(geneID, sep = "/") 

genes.select  = c(
    "hmcn1","adamts2","col17a1a","col4a2","col8a2","postna","lum","fbln2","scara3",
    "col1a1a","col11a2","col11a1b","col9a2","col10a1a","col10a1b","col27a1b","col4a6","col5a3b","adamts3","adamts17")

pdf(paste0(celltype,"_allstage_DEGs_3module_GOterm_",key.word,"_gene_Clustered_Dotplot_selected.pdf"),width=8,height=3)
Clustered_DotPlot(
    seurat_object = mes, features = genes.select,k=2,group.by="Stage",
    cluster_ident = F,flip=T)
dev.off()

key.word = "skeletal system development" 
df = data.frame(go@compareClusterResult)
# temp.df=df[df$Description %like% key.word,] #data.table
# temp.df <- temp.df %>%separate_rows(geneID, sep = "/") 

genes.select  = c(
    "ednraa","loxl2b","ucmaa","skia","smoc2","furinb","lpar1","rbms3","dlx3b",
    "and1","and2","hapln1b","mmp9","dcn","nog1","sp7","runx2a","hoxa9a","prrx1b")

pdf(paste0(celltype,"_allstage_DEGs_3module_GOterm_",key.word,"_gene_Clustered_Dotplot_selected.pdf"),width=8,height=3)
Clustered_DotPlot(
    seurat_object = mes, features = genes.select,k=2,group.by="Stage",
    cluster_ident = F,flip=T)
dev.off()

key.word = "plasma membrane bounded cell projection morphogenesis" 
df = data.frame(go@compareClusterResult)
# temp.df=df[df$Description %like% key.word,] #data.table
# temp.df <- temp.df %>%separate_rows(geneID, sep = "/") 

genes.select  = c(
    "plxna1a","sema6bb","ephb3a","ephb2b","sema4ba","sema3ab","sema3fa","sema5ba","plxnb2a",
    "sema3e","sema6ba","plxna3","sema6e","epha3","plxnb3","efnb1","ephb6","sema5a","efnb3a")

pdf(paste0(celltype,"_allstage_DEGs_3module_GOterm_",key.word,"_gene_Clustered_Dotplot_selected.pdf"),width=8,height=3)
Clustered_DotPlot(
    seurat_object = mes, features = genes.select,k=2,group.by="Stage",
    cluster_ident = F,flip=T)
dev.off()
```
