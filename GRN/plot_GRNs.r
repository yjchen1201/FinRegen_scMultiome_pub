
library(dplyr)
library(tibble)
library(ggplot2)
library(ggnetwork)
library(network)
library(sna)
library(ggnewscale)

setwd("step6_new_identify_enriched_GRN_sub-networks")
celltypes= c("Superficial_Epithelial",
      "Intermediate_Epithelial",  
      "Basal_Epithelial",
      "Epidermal_Mucous",
      "Hematopoietic",
      "Mesenchymal")
stages=c("1dpa","2dpa","4dpa","6dpa")
types=c("Act","Rep")
for (celltype in celltypes){
  print(paste(celltype))
  for (stage in stages){
    print(paste(celltype,stage,sep=" "))
    grn = read.csv(file=paste0(celltype,"_GRNs_triple_specific_",stage,".csv"))
    for (type in types) {
      grn.filter <- grn %>% 
      filter(Peak_fold_2 >0) %>% 
      filter(type == {{type}}) %>% 
      arrange(desc(abs(GG_Corr))) %>% 
      select(TF, Target, GG_Corr) %>% 
      unique() %>%
      top_n(200)

      # Load required libraries

      # Create a unique list of all entities (TFs and Targets)
      nodes <- unique(c(grn.filter$TF, grn.filter$Target))

      # Create a mapping (e.g., using an index)
      node_mapping <- setNames(seq_along(nodes), nodes)

# Replace TF and Target in the data frame with their unique identifiers
grn.filter$TF <- node_mapping[grn.filter$TF]
grn.filter$Target <- node_mapping[grn.filter$Target]

# Aggregate the edges by taking the mean (or sum) of the scores
grn.filter.agg <- grn.filter %>% group_by(TF, Target) #%>% summarize(GG_Corr = mean(GG_Corr), .groups = 'drop')

# Now create the network with the aggregated data
net <- network(grn.filter.agg, directed = TRUE)

# network.vertex.names(net) <- grn.filter$TF

# # Create a layout for your network
set.seed(123)
net_layout <- ggnetwork(net, layout = 'fruchtermanreingold')

# Add a 'type' column to distinguish between TF and Target
net_layout$group <- ifelse(net_layout$vertex.names %in% grn.filter$TF, 'TF', 'Target')
net_layout$gene_name = names(node_mapping)[match(net_layout$vertex.names,node_mapping)]
net_layout$size <- ifelse(net_layout$group == 'TF', 3, 1)  # Adjust the numbers as needed

################# use log2FC #####################

gene_avg_log2FC = read.csv("major_celltype_all_stage_0vsx_avglog2FC_gene_expression.csv",row.names=1)
gene_avg_log2FC = -gene_avg_log2FC #change comparison direction, so log2FC>0 means higher expression in postinjury
rrg.list = readRDS("all_celltype_all_bystage_DEGs.rds")
rrg = rrg.list[[celltype]]
columns = c(paste0(celltype,"_0vs1_log2FC"),paste0(celltype,"_0vs2_log2FC"),paste0(celltype,"_0vs4_log2FC"),paste0(celltype,"_0vs6_log2FC"))
ct_rrg_avg_log2FC = as.data.frame(gene_avg_log2FC[rrg,columns])
ct_rrg_avg_log2FC$gene_short =as.character(rownames(ct_rrg_avg_log2FC))
all_gene_log2fc = as.data.frame(gene_avg_log2FC[,columns])
colnames(all_gene_log2fc) = c("1dpa","2dpa","4dpa","6dpa")
exp_column = stage
expression_data <- data.frame(
  gene_name = rownames(all_gene_log2fc),
  exp.log2FC = all_gene_log2fc[, exp_column]
)

net_layout$exp.log2FC <- expression_data$exp.log2FC[match(net_layout$gene_name, expression_data$gene_name)]
# Plot the network with different gradient colors for TFs and Targets
#library(ggnewscale)
# Create a column for absolute GG_Corr values
net_layout$GG_Corr_abs <- abs(net_layout$GG_Corr)
# Load the ggnewscale package for multiple color scales
library(ggnewscale)
pdf(paste(celltype, stage, type, "top200_GRN.pdf", sep = "_"), width = 5, height = 5)
print(
  ggplot(net_layout, aes(x = x, y = y, xend = xend, yend = yend)) +
  # Edges with gradient color based on absolute GG_Corr
  geom_edges(aes(color = GG_Corr_abs)) +
  scale_color_gradient(low = "#CFCFCF", high = "#363636", name = "Correlation Strength") +
  
  # New color scale for TF nodes
  new_scale_color() +
  geom_nodes(data = subset(net_layout, group == "TF"),
             aes(color = exp.log2FC, size = size)) +
  scale_color_gradient(low = "lightblue", high = "#626496", name = "TF Expression") +
  
  # New color scale for Target nodes
  new_scale_color() +
  geom_nodes(data = subset(net_layout, group == "Target"),
             aes(color = exp.log2FC, size = size)) +
  scale_color_gradient(low = "pink", high = "red", name = "Target Expression") +
  
  # Node text layers
  geom_nodetext(data = subset(net_layout, group == "TF"),
                aes(label = gene_name), vjust = 0.5, check_overlap = TRUE, size = 2) +
  geom_nodetext(data = subset(net_layout, group == "Target"),
                aes(label = gene_name), vjust = 1.5, check_overlap = TRUE, size = 1.5) +
  
  theme_void() +
  guides(size = "none", color = guide_legend(title = "Expression")) +
  ggtitle(paste(celltype, stage, type, "GRN", sep = " "))
  
)
dev.off()

}
}
}
