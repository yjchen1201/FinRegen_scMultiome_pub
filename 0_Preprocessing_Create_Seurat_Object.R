#!/usr/bin/env Rscript
# Retrieve system arguments
args = commandArgs(trailingOnly=TRUE)
# test if there 4 arguments: if fewer, return an error
if (length(args) < 4) {
  stop("Please provide absolute full path for working directory, sample name, stage name, cellrangerarc output directory.n", call.=FALSE)
} else if (length(args) > 4) {
  # return an error when there are more 4 errors
  stop("Please remove unused arguments.", call.=FALSE)
}

working_path <- args[1]  #working_path<-"/scratch/ichen/FinRegene_10xMultiome/Signac_processed"
sample <- args[2] #sample<-"6dpa1"
stage<-args[3] #"6dpa"
cellrangerARC_output_path<-args[4] #"/scratch/ichen/FinRegene_10xMultiome/cellrangerARC_processed/6dpa1/outs"


count_mtx_path<-paste(cellrangerARC_output_path,"filtered_feature_bc_matrix.h5",sep="/")
frag_path <- paste(cellrangerARC_output_path,"atac_fragments.tsv.gz",sep="/")
metadata.path<-paste(cellrangerARC_output_path,"per_barcode_metrics.csv",sep="/")
stats_file<-paste(sample,"stats_record.txt",sep="_")


# Load packages 
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Drerio.UCSC.danRer11)
#set multicores for mtx count
library(future)
set.seed(1234)
plan("multicore", workers = 12)
options(future.globals.maxSize = 120000 * 1024^2)


# get gene annotations for danRer11
annotation<-readRDS("/scratch/ichen/FinRegen_10xMultiome/annotations/Danio_rerio_GRCz11.104_Ensembl_gtf_annotation_onlyChr_wEGFP.rds")
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "danRer11"

sample_dir <- paste(working_path,sample,sep="/")
dir.create(sample_dir) #create sample directory
setwd(paste(working_path,sample,sep="/"))
file.create(stats_file) 
step <-"0"

counts <- Read10X_h5(count_mtx_path)

# create a Seurat object containing the RNA adata

obj<-CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  project=sample,
  min.cells = 5
)
# Add sample name to metadata orig.ident column and stage to Stage column
obj$orig.ident<-sample
obj<-SetIdent(obj, value="orig.ident")
obj@meta.data$Stage <- stage

# create ATAC assay and add it to the object
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = paste(frag_path),
  annotation = annotation,
  min.cells = 5
)

#read raw per barcode metadata 
cellrangerARC.metadata<- read.csv(
  file = paste(metadata.path),
  header = TRUE,
  row.names = 1
)

####################### 0. Quality Control Plot #######################
step<-"0"
DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-") 

DefaultAssay(obj) <- "ATAC"
obj<- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

# optional: save the original object
saveRDS(obj,file=paste(step,sample,"cellrangerARC_seurat_obj.rds",sep="_"))


# pdf(paste(step,sample,"QC_plot_wDots.pdf",sep="_"),width = 23, height = 5)
# VlnPlot(
#    object = obj,
#    features = c("nCount_RNA","nFeature_RNA","log10GenesPerUMI","nCount_ATAC", "nFeature_ATAC","TSS.enrichment", "nucleosome_signal","percent.mt"),
#    ncol = 8,
#    pt.size = 0.1
# )
# dev.off()

# pdf(paste(step,sample,"QC_plot.pdf",sep="_"),width = 23, height = 5)
# VlnPlot(
#    object = obj,
#    features = c("nCount_RNA","nFeature_RNA","log10GenesPerUMI","nCount_ATAC", "nFeature_ATAC","TSS.enrichment", "nucleosome_signal","percent.mt"),
#    ncol = 8,
#    pt.size = 0
# )
# dev.off()

# #### Histogram showing distribution of nFeature and nCount
# pdf(paste(step,sample,"nCounts_nFeature_distribution_histPlot.pdf",sep="_"), width=18,height=3)
# par(mfrow=c(1,7))
# hist(log10(obj$nCount_RNA),n=1000,col="darkgrey",border="darkgrey",main="nCount_RNA")
# hist(log10(obj$nFeature_RNA),n=1000,col="darkgrey",border="darkgrey",main="nFeature_RNA")
# hist(log10(obj$nCount_ATAC),n=1000,col="darkgrey",border="darkgrey",main="nCount_ATAC")
# hist(log10(obj$nFeature_ATAC),n=1000,col="darkgrey",border="darkgrey",main="nFeature_ATAC")
# hist(obj$nucleosome_signal,n=1000,col="darkgrey",border="darkgrey",main="Nucleosome_signal")
# hist(obj$TSS.enrichment,n=1000,col="darkgrey",border="darkgrey",main="TSS_enrichment")
# hist(obj$percent.mt,n=1000,col="darkgrey",border="darkgrey",main="Pct_mt")
# dev.off()


# obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
# pdf(paste(step, sample,"fragment_distribution.pdf",sep="_"))
# FragmentHistogram(object = obj, group.by = 'nucleosome_group')
# dev.off()

# subset object based on QC
obj_subset <- subset(
  x = obj,
  subset = nCount_ATAC < 2e5 & #exception: based on distribution, 6dpa2 use 5e5 here
    nCount_ATAC > 1000 &
    nCount_RNA < 5e4 & #exception: based on distribution, 6dpa1 use 6e5 here
    nCount_RNA > 1000 & 
    percent.mt < 10 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
)

cat("Cells after basic QC filtering:", length(colnames(obj_subset)), "\n", file = stats_file, append = TRUE)


####### 1. Macs2 peak calling using filtered cells #######
# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
cat("step1_ATAC_callpeak:","\n", file = stats_file, append = TRUE)

step<-"1"
assay="ATAC"
DefaultAssay(obj_subset) <- assay

macs2_path <- paste(working_path,sample,"macs2_callpeak",sep="/")
dir.create(macs2_path)

macs2peaks<-CallPeaks(
  object=obj_subset,
  macs2.path = "/opt/apps/python2/bin/macs2",
  outdir = macs2_path,
  extsize = 200,
  shift = -100,
  effective.genome.size = 1.4e+09,
  additional.args = "--bdg --SPMR -f BEDPE",
  name = sample,
  cleanup = FALSE,
  verbose = TRUE,
)
saveRDS(macs2peaks, file=paste(step,sample,"macs2_peaks.rds",sep="_"))
length(macs2peaks)
write.table(as.data.frame(macs2peaks),file=paste(step,sample,"macs2_peaks_wBEDPE.txt",sep="_"),sep="\t",quote=FALSE,row.names=FALSE)

cat("macs2 called peak counts:", length(macs2peaks), "\n", file = stats_file, append = TRUE)

# remove peaks on nonstandard chromosomes (we don't have, since we removed the nonstandard chr during cellranger-arc count)
macs2peaks <- keepStandardChromosomes(macs2peaks, pruning.mode = "coarse")

# remove peaks in genomic blacklist regions
blacklist_danRer11 <-rtracklayer::import("/scratch/ichen/zebrafish_blacklist/Blacklist_danRer10_to_danRer11_YueLab_srt.bed")
macs2peaks <- subsetByOverlaps(x = macs2peaks, ranges = blacklist_danRer11, invert = TRUE)

cat("macs2 called peaks after removing blacklist region peaks:", length(macs2peaks), "\n", file = stats_file, append = TRUE)

saveRDS(macs2peaks, file=paste(step,sample,"macs2_peaks_afterfiltering.rds",sep="_"))

# create macs2peak-cell matrix
obj_subset_macs2counts <- FeatureMatrix(
  fragments = Fragments(obj_subset),
  features = macs2peaks,
  cells = colnames(obj_subset),
  process_n = 2000,
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
obj_subset[["ATAC_macs2"]] <- CreateChromatinAssay(
  counts = obj_subset_macs2counts,
  fragments = frag_path,
  annotation = annotation,
  min.cells = 5
)

cat("step1 macs2peak assay feature-cell:", length(rownames(obj_subset[["ATAC_macs2"]])),length(colnames(obj_subset[["ATAC_macs2"]])), "\n", file = stats_file, append = TRUE)

# Calculate reads in peak ratio
obj_subset <- FRiP(object = obj_subset,
               assay = 'ATAC_macs2',
               total.fragments = 'nCount_ATACmacs2')

# add number of peak region fragments for Signac called peaks (equals the value in the "nCount_peaks" value in the object)
obj_subset$macs2peak_region_fragments <- obj_subset$FRiP * obj_subset$nCount_ATACmacs2

#  Add pct_reads_in_peaks for Signac called peaks to seurat object
obj_subset$pct_reads_in_macs2peaks <- obj_subset$macs2peak_region_fragments / obj_subset$nCount_ATACmacs2 * 100


# QC violin plot
pdf(paste(step,sample, "macs2peaks_QC_plot_wDots.pdf",sep="_"), width = 20, height = 10) 
VlnPlot(
  object = obj_subset,
  group.by= "orig.ident",
  features = c("nCount_RNA","nFeature_RNA","nCount_ATAC", "nCount_ATAC_macs2", "nFeature_ATAC", "nFeature_ATAC_macs2","TSS.enrichment", "nucleosome_signal","percent.mt",
  'FRiP','pct_reads_in_macs2peaks', 'macs2peak_region_fragments'),
  pt.size = 0.1,
  ncol = 6,
  col="pink"
)
dev.off()

pdf(paste(step,sample, "macs2peaks_QC_plot.pdf",sep="_"), width = 20, height = 10) 
VlnPlot(
  object = obj_subset,
  group.by= "orig.ident",
  features = c("nCount_RNA","nFeature_RNA","nCount_ATAC", "nCount_ATAC_macs2", "nFeature_ATAC", "nFeature_ATAC_macs2","TSS.enrichment", "nucleosome_signal","percent.mt",
  'FRiP','pct_reads_in_macs2peaks', 'macs2peak_region_fragments'),
  pt.size = 0,
  ncol = 6,
  col="pink"
)
dev.off()

saveRDS(obj_subset,file=paste(step,sample,"subset_seurat_object.rds",sep="_"))


################# diet seurat object to RNA and ATAC ASSAY ##############
# RNA assay
DefaultAssay(obj_subset)<-"RNA"
obj_rna<-DietSeurat(
  object =obj_subset,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  assays = "RNA"
)

####### SCT Normalization #########
# Import cell cycle genes and calculate cell cycle scores after normalization
fishCCgenes <- readLines(con = "/bar/yhou/SingleCell/FinReg10xSCRNAseurat_Rproj/AGGallSamp/cyclingCells/cell_cycle_vignette_files/regev_lab_cell_cycle_genes_asFish.txt")
s.genes <- fishCCgenes[1:42]
g2m.genes <- fishCCgenes[43:96]

# SCTransform and cell cycle regression 
obj_rna <- SCTransform(obj_rna, assay = "RNA",new.assay.name = "SCT",vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
obj_rna <- CellCycleScoring(obj_rna, s.features = s.genes, g2m.features = g2m.genes)
obj_rna$Phase <- factor(obj_rna$Phase, levels = c('G1', 'S', 'G2M'))
obj_rna$CC.Difference <- obj_rna$S.Score - obj_rna$G2M.Score
obj_rna <- ScaleData(obj_rna, vars.to.regress = "CC.Difference", features = rownames(obj_rna))

DefaultAssay(obj_rna)<-"RNA"
obj_rna <- NormalizeData(obj_rna)

DefaultAssay(obj_rna)<-"SCT"
saveRDS(obj_rna,file=paste(step,sample,"subset_seurat_RNA_object.rds",sep="_"))


# 
DefaultAssay(obj_subset)<-"ATAC_macs2"
obj_atac<-DietSeurat(
  object= obj_subset,
  counts = TRUE,
  data = TRUE,
  assays = "ATAC_macs2"
)

obj_atac<- RunTFIDF(obj_atac)
obj_atac <- FindTopFeatures(obj_atac, min.cutoff = 'q5')
obj_atac <- RunSVD(obj_atac)

saveRDS(obj_atac,file=paste(step,sample,"subset_seurat_ATAC_object.rds",sep="_"))

