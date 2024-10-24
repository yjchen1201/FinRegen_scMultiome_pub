#!/usr/bin/env Rscript
# Example: 
## module load R/4.0.3
## Rscript p2g_preinjury_and_regeneration.r -o Epidermal_Mucous_subset_object.rds  --p "Epidermal_Mucous" 
library("optparse")

# option_list = list(
#   make_option(c("-o", "--obj_file_path"), type="character", default=NULL, 
#               help="seurat object", metavar="character"),
#     make_option(c("-p", "--prefix"), type="character", default=NULL, 
#               help="output file name [default= %default]", metavar="character")
# ); 

suppressPackageStartupMessages({
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(plyr)
library(dplyr)
library(ggplot2)
library(BSgenome.Drerio.UCSC.danRer11)
#set multicores for mtx count
library(future)
set.seed(1234)
plan("multicore", workers = 6)
options(future.globals.maxSize = 90000 * 1024^2)
})
# opt_parser = OptionParser(option_list=option_list); opt = parse_args(opt_parser);


args = commandArgs(trailingOnly=TRUE)

obj_file_path = args[1]
prefix = args[2]


obj =readRDS(obj_file_path)
Idents(obj)="Stage"

regen.obj = subset(x = obj, idents = c("mfin"), invert = TRUE)
preinjury.obj = subset(x = obj, idents = c("mfin"))


prefix=prefix

DefaultAssay(regen.obj) <- "RNA"
regen.obj <- SCTransform(regen.obj)

DefaultAssay(regen.obj) <- "peaks"
regen.obj <- FindTopFeatures(regen.obj, min.cutoff = 5)
regen.obj<- RunTFIDF(regen.obj)
# first compute the GC content for each peak
regen.obj <- RegionStats(regen.obj, genome = BSgenome.Drerio.UCSC.danRer11)

# link peaks to genes
regen.obj <- LinkPeaks(
  object = regen.obj,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 100000, #50000
  genes.use = rownames(regen.obj@assays$SCT@data),
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  min.cells = 10,
  n_sample = 200
)
regen.links=Links(regen.obj)
saveRDS(regen.links,file=paste0(prefix,"_regenerating_peak_gene_links_100kb.rds"))



DefaultAssay(preinjury.obj) <- "RNA"
preinjury.obj <- SCTransform(preinjury.obj)

DefaultAssay(preinjury.obj) <- "peaks"
preinjury.obj <- FindTopFeatures(preinjury.obj, min.cutoff = 5)
preinjury.obj<- RunTFIDF(preinjury.obj)
# first compute the GC content for each peak
preinjury.obj <- RegionStats(preinjury.obj, genome = BSgenome.Drerio.UCSC.danRer11)

# link peaks to genes
preinjury.obj <- LinkPeaks(
  object = preinjury.obj,
  peak.assay = "peaks",
  expression.assay = "SCT",
  distance = 100000, #50000
  genes.use = rownames(preinjury.obj@assays$SCT@data),
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  min.cells = 10,
  n_sample = 200
)
preinjury.links=Links(preinjury.obj)
saveRDS(preinjury.links,file=paste0(prefix,"_preinjury_peak_gene_links_100kb.rds"))

