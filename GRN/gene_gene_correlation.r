#!/usr/bin/Rscript
# Rscript gene_gene_correlation.r celltype stage
args = commandArgs(trailingOnly=TRUE)
celltype = args[1]
stage = args[2]

set.seed(123)
library(dplyr)
library(doParallel)
library(reshape2)
# celltype="Basal_Epithelial"


sparse.cor3 <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  cSums <- colSums(x)
  # Calculate the population covariance matrix.
  # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
  # The code is optimized to minize use of memory and expensive operations
  covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
  crossp <- as.matrix(crossprod(x))
  covmat <- covmat+crossp
  sdvec <- sqrt(diag(covmat)) # standard deviations of columns
  covmat/crossprod(t(sdvec)) # correlation matrix
}

mat = read.csv(paste0(celltype,"_",stage,"_gene_exp_matrix.csv"),row.names=1) 
mat =t(mat) #gene as column
mat = as.matrix(mat)

k = which(colSums(mat) == 0)
k
mat_cl = mat[,-k,drop = FALSE]
dim(mat_cl)
Corr_res = sparse.cor3(mat_cl)
Corr_res = melt(Corr_res)
saveRDS(Corr_res,file=paste0(celltype,"_",stage,"_Corr_res.rds"))

