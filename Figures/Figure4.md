## Fig4A
```r
library(RColorBrewer)
library(colorRamp2)
library(ComplexHeatmap)

celltypes<-c("Superficial_Epithelial",
    "Intermediate_Epithelial",  
    "Basal_Epithelial",
    "Epidermal_Mucous",
    "Hematopoietic",
    "Mesenchymal")

p2g.list<-list()
for (ct in celltypes){
  path= paste0("p2g_links_100kb/",ct,"_peak_gene_links_100kb.csv")
  p2g = read.csv(path)
  p2g$celltype = ct #ct name
  p2g$p2g = paste(p2g$peak, p2g$gene, sep="@") #p2g string used for filtering later
  p2g.list[[ct]] = p2g
  #p2g = merge(p2g, peakAnno.df[,c("peak","annotation","nearestTSSGene","nearestTSSGene_dist")], by = "peak", all.x = TRUE) # add chipseeker annotation info 
  }

gene.list=list()
for (ct in celltypes){
  p2g = readRDS(paste0("p2g_links_100kb/",ct,"_peak_gene_links_100kb.rds"))
  numlink = data.frame(p2g) %>% group_by(gene) %>% tally() %>% arrange(desc(n))
  gene.list[[ct]] = numlink
  gene.list[[ct]]$celltype = ct
}

for (ct in celltypes){
  gene = gene.list[[ct]]$gene[gene.list[[ct]]$n>=5]
  p2g.list[[ct]] = p2g.list[[ct]][which(p2g.list[[ct]]$gene %in% gene),]
  }

p2g.mm <- do.call('rbind',p2g.list)
p2g.pair =p2g.mm %>% select(peak, gene) %>% unique() #49763
peak.exp = read.csv("0_integrated_wnn_merged_celltype_peak_expression_bycelltype_bystage.csv",row.name=1)
peak.exp [is.na(peak.exp)] <- 0
colnames(peak.exp)=gsub("mfin","0dpa",colnames(peak.exp))

gene.exp = read.csv("0_integrated_wnn_merged_celltype_rna_expression_bycelltype_bystage.csv",row.name=1)
gene.exp [is.na(gene.exp)] <- 0
colnames(gene.exp)=gsub("mfin","0dpa",colnames(gene.exp))

# order gene and peak matrix rows by paired p2g link list
peak.exp = peak.exp[p2g.pair$peak, ]  
gene.exp = gene.exp[p2g.pair$gene, ] 

# get major cell type columns only
stages <- c("0dpa","1dpa","2dpa","4dpa","6dpa")
combinations <- expand.grid(stages = stages,celltypes = celltypes)
combined_strings <- paste(combinations$celltypes, combinations$stages, sep = "_")

peak.exp = peak.exp[,combined_strings]
gene.exp = gene.exp[,combined_strings]

# scale matrix
mat_peak <- t(scale(t(as.matrix(peak.exp))))
mat_gene <- t(scale(t(as.matrix(gene.exp))))

# get row order
set.seed(123)
km.row <- kmeans(mat_peak, 6, nstart = 25)
row_order <- order(km.row$cluster)

# get column order
set.seed(123)
km.column <- kmeans(t(mat_peak), 6, nstart = 25)
column_order <- order(km.column$cluster)

######## test and check the orders #########
ht_peak = Heatmap(mat_peak[row_order,column_order], name = "Peak", 
                  col = peak_col_fun,
                  cluster_rows=F,cluster_columns=F,
                  show_row_names = FALSE,
                  show_column_names = T,
                  column_title = "Peak",
                  top_annotation = stage_annotation)

# For the paired gene heatmap
ht_gene = Heatmap(mat_gene[row_order,column_order], name = "Gene", 
                  col = gene_col_fun,
                  cluster_rows=F,cluster_columns=F,
                  show_row_names = FALSE,
                  show_column_names = T,
                  row_names_gp = grid::gpar(fontsize = 4),
                  column_title = "Gene",
                  top_annotation = stage_annotation)

# Combine the two heatmaps
ht_list = ht_peak + ht_gene
draw(ht_list, #row_split = celltype_bar_anno,
column_title = "Peak to gene", 
column_title_gp = gpar(fontsize = 12), 
merge_legends = TRUE, heatmap_legend_side = "bottom")   
dev.off()

# reorder columns 
celltypes<-c("Basal_Epithelial",
    "Mesenchymal",
    "Intermediate_Epithelial", 
    "Superficial_Epithelial",
    "Epidermal_Mucous",
    "Hematopoietic"
    )

stages <- c("0dpa","1dpa","2dpa","4dpa","6dpa")
combinations <- expand.grid(stages = stages,celltypes = celltypes)
column_order_final <- paste(combinations$celltypes, combinations$stages, sep = "_")

# get cluster split vector
cluster_vector <- km.row$cluster[row_order]

require(colorRamp2)
peak_col_fun = colorRamp2(c(-2, 0, 2), c("#481C66", "#009F94", "#FDE333"))
gene_col_fun = colorRamp2(c(-2, 0, 2), c("#001889", "#C73D7B", "#FFFE9E"))

pdf("all_p2g_peak_heatmap_p2g100kb_w5links_column_reordered.pdf",width=5,height=7)
stages <- c("0dpa","1dpa","2dpa","4dpa","6dpa")
stages_column = list(stage = c(rep(stages, 6)))
stage_annotation <- HeatmapAnnotation(df = stages_column, 
                                       col = list(stage = c("0dpa" = "#d1d0ff", "1dpa" = "#b19fed", 
                                                           "2dpa" = "#906ddb", "4dpa" =  "#703cc9", 
                                                           "6dpa" = "#4f0ab7")), 
                                                           show_legend = TRUE, height = unit(0.2, "cm"))


ht_peak = Heatmap(mat_peak[row_order,column_order_final], name = "Peak", 
                  col = peak_col_fun,
                  cluster_rows=F,cluster_columns=F,
                  row_split = cluster_vector,   # Apply the cluster vector here
                  show_row_names = FALSE,
                  show_column_names = T,
                  column_title = "Peak",
                  top_annotation = stage_annotation)

# For the paired gene heatmap
ht_gene = Heatmap(mat_gene[row_order,column_order_final], name = "Gene", 
                  col = gene_col_fun,
                  cluster_rows=F,cluster_columns=F,
                  row_split = cluster_vector,   # Apply the cluster vector here too
                  show_row_names = FALSE,
                  show_column_names = T,
                  row_names_gp = grid::gpar(fontsize = 4),
                  column_title = "Gene",
                  top_annotation = stage_annotation)

# Combine the two heatmaps
ht_list = ht_peak + ht_gene
draw(ht_list, #row_split = celltype_bar_anno,
column_title = "Peak to gene", 
column_title_gp = gpar(fontsize = 12), 
merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

# extract genes from each cluster
cluster_vector <- km.row$cluster[row_order]
gene_order = p2g.pair$gene[row_order]
clustered_genes <- split(gene_order, cluster_vector)
clustered_genes <- lapply(clustered_genes, function(x) unique(sort(x))) #have to sort before making matrix

for (i in 1:length(clustered_genes)){
    cluster_genes <- clustered_genes[[i]]
    print(length(cluster_genes))
    write.table(cluster_genes, paste0("all_p2g_peak_heatmap_p2g100kb_w5links_cluster",i,"_genes.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
}
```

## Fig4B
```r
# Cluster GO terms bar plot
library(openxlsx)
library(ggplot2)
library(dplyr)
go <-read.xlsx("selected_GOs_p2g_cluster.xlsx",sheet="5peaks")
df = go %>% select(Description, LogP,cluster, celltype) 
# Create an interaction term and arrange
df <- df %>%
  arrange(cluster, -LogP) %>%
  mutate(cluster_Description = interaction(cluster, Description, lex.order = TRUE))

# Reverse the order of levels for the factor
df$cluster_Description <- with(df, factor(cluster_Description, levels = rev(unique(cluster_Description))))

# Plot
p <- ggplot(df, aes(x = cluster_Description, y = -LogP, fill = -LogP)) + 
  geom_bar(stat = "identity") + 
  scale_fill_gradientn(colors=rev( paletteer::paletteer_c("grDevices::Purple-Blue", 30)[1:15] )) + # Change color scale
  coord_flip() + # flip coordinates for better visualization
  labs(title = "", 
       y = "-log10P", x = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank()) 

pdf("gene_cluster_GO_terms.pdf")
print(p)
dev.off()
```


## Figure4C-E
```markdown
#### genome browser plots are generated using WashU Epigenome Browser
```
```bash
## Example of generating browser bigwig track file ##
mkdir "$atac_output_dir"
cd "$atac_output_dir"
#reference: https://timoast.github.io/sinto/basic_usage.html#filter-cell-barcodes-from-bam-file
sinto filterbarcodes -b "$atac_bam" -c "$barcode_file" -p "$thread" > "$sample"_sinto_filterbarcodes.log 2>&1 &
sinto filterbarcodes -b "$atac_bam" -c "$combined_barcode_file" -p "$thread" > "$sample"_sinto_combined_HEM_MES_EM_filterbarcodes.log 2>&1 &
# Merge two replicates for each time point
module load python3 samtools 
samtools merge -@ 10 ${celltype}"_"${stage_rep1}".bam" ${celltype}"_"${stage_rep2}".bam" ${celltype}"_"${stage}".bam"
samtools sort -@ 10 ${celltype}"_"${stage}".bam" -o ${celltype}"_"${stage}".sorted.bam" &
samtools index -@ 12 ${celltype}"_"${stage}".sorted.bam"
# Convert bam file to bigwig file. bigwig file is the input file of genome browser to generate signal track.
bamCoverage  -of bigwig --normalizeUsing CPM --binSize 10 --effectiveGenomeSize 1400000000 \
-b "RNA_"$stage"_${celltype}_sorted.bam" \
-o "RNA_"$stage"_${celltype}_sorted.bam"  \
-p 6  > "RNA_"$stage"_${celltype}_sorted.bamtobigw" 2>&1 &
```
