# TF activity TF gene expression correlation 
Get_corr_gene_motif <- function(motif_matrix,gene_matrix,Fish_All_motifs_table_short){

  ########
  k = which(Matrix::rowSums(gene_matrix) ==0) # gene in row, cell in column
  if(length(k) > 1){
    gene_matrix = gene_matrix[-k,] #gene 24380->20267 in mes
  }
  ######## clean cols #####
  cells = colnames(gene_matrix)[which(colnames(gene_matrix) %in% colnames(motif_matrix) == T)]
  gene_matrix_cl = gene_matrix[,which(colnames(gene_matrix) %in% cells == T)]
  motif_matrix_cl = motif_matrix[,which(colnames(motif_matrix) %in% cells == T)]
  #######
  gM = match(cells,colnames(gene_matrix_cl))
  gene_matrix_cl = gene_matrix_cl[,gM]
  mM = match(cells,colnames(motif_matrix_cl))
  motif_matrix_cl = motif_matrix_cl[,mM]
  ########
  all.equal(colnames(gene_matrix_cl),colnames(motif_matrix_cl))
  ########
  Fish_All_motifs_table_short$Cor = 0
  ########
  Fish_All_motifs_table_short_cl = Fish_All_motifs_table_short[which(Fish_All_motifs_table_short$Gene %in% rownames(gene_matrix_cl) == T),]
  Fish_All_motifs_table_short_cl$Cor = 0
  ########
  ########
  library(parallel)
  cl <- makeCluster(6)
  cor_res = parApply(cl,Fish_All_motifs_table_short_cl,1,cor_functon,motif_matrix_cl,gene_matrix_cl)
  stopCluster(cl)
  #########
  Fish_All_motifs_table_short_cl$Cor = cor_res
  return(Fish_All_motifs_table_short_cl)
}


cor_functon <- function(x,motif_matrix_cl,gene_matrix_cl){
  tmp1 = x[1]
  tmp2 = x[2]
  #######
  vector1 = which(rownames(motif_matrix_cl) == tmp1)
  vector2 = which(rownames(gene_matrix_cl) == tmp2)
  if(length(vector1) == 1 & length(vector2) == 1){
    v1 = motif_matrix_cl[vector1,]
    v2 = gene_matrix_cl[vector2,]
    ### all.equal(names(v1),names(v2))
    Cor = cor(v1,v2,method = c("spearman"))
    ###
    return(Cor)
  }else{
    Cor = 0
    return(Cor)
  }
}


# gene-gene correlation 
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

###################  add peak annotatation to peak-to-gene link ##############
Get_TSS_table <- function(GR_list,GR_peak,extend=1000){
  ##########
  outtable_list <- list()
  for(i in 1:length(GR_list)){
    tmp_GR = GR_list[[i]]
    ######
    start(tmp_GR) = start(tmp_GR) - extend
    end(tmp_GR) = end(tmp_GR) + extend
    ######
    k = which(countOverlaps(GR_peak,tmp_GR) > 0)
    if(length(k) > 0){
      GR_peak_sub = GR_peak[k]
      #####
      outtable_tmp = data.frame(Gene=tmp_GR$genes[1],Peak=as.character(GR_peak_sub),Class="TSS")
      outtable_list = c(outtable_list,list(outtable_tmp))
    }
    if(i %% 1000 == 0){
      print(i)
    }
  }
  outtable_list_merge = do.call("rbind",outtable_list)
  ##########
  return(outtable_list_merge)
}

Get_Body_table <- function(GR_list,GR_peak){
  ##########
  outtable_list <- list()
  for(i in 1:length(GR_list)){
    tmp_GR = GR_list[[i]]
    ######
    ######
    k = which(countOverlaps(GR_peak,tmp_GR) > 0)
    if(length(k) > 0){
      GR_peak_sub = GR_peak[k]
      #####
      outtable_tmp = data.frame(Gene=tmp_GR$genes[1],Peak=as.character(GR_peak_sub),Class="Body")
      outtable_list = c(outtable_list,list(outtable_tmp))
    }
    if(i %% 1000 == 0){
      print(i)
    }
  }
  outtable_list_merge = do.call("rbind",outtable_list)
  ##########
  return(outtable_list_merge)
}


Get_all_peak_Gene_tables <- function(GTF_TSS_table_res,GTF_Body_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
  ####
  GTF_TSS_table_res_1 = Process_TSS_peaks(GTF_TSS_table_res,PtoG=PtoG,Gtfs_protein_trans_GR_TSS)
  GTF_Body_table_res_1 = Process_Body_peaks(GTF_Body_table_res,PtoG=PtoG,Gtfs_protein_trans_GR_TSS)
  ####
  GTF_TSSBody_table_res = Merge_TSS_and_Body(GTF_TSS_table_res_1,GTF_Body_table_res_1)
  ####
  GTF_Inter_table_res_1 = Process_Inter_peaks(GTF_TSSBody_table_res,PtoG=PtoG,Gtfs_protein_trans_GR_TSS)
  ####
  ALL_peak_table = rbind(GTF_TSS_table_res_1,GTF_Body_table_res_1,GTF_Inter_table_res_1)
  ####
  ALL_peak_table_anno = Process_peak_table_to_annotation(ALL_peak_table)
  ####
  return(ALL_peak_table_anno)
}



Process_Body_peaks <- function(GTF_Body_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
  ###### First add PtoG link correaltions !!! ###########
  GTF_Body_table_res$Index = paste0(GTF_Body_table_res$Gene,"#",GTF_Body_table_res$Peak)
  #GTF_Body_table_res$Index= gsub(":","-",GTF_Body_table_res$Index) #to match Signac peak name style
  ######
  GTF_Body_table_res$PtoG_cor = "ND"
  GTF_Body_table_res$PtoG_pvalue = "ND"
  ######
  PtoG_index = paste0(PtoG$gene,"#",PtoG$peak)
  ######
  m = match(GTF_Body_table_res$Index,PtoG_index)
  GTF_Body_table_res$PtoG_cor = PtoG$score[m]
  GTF_Body_table_res$PtoG_pvalue = PtoG$pvalue[m]
  #######
  k = which(is.na(GTF_Body_table_res$PtoG_pvalue) == T)
  GTF_Body_table_res_cl = GTF_Body_table_res[-k,]
  #######
  ####### Next see the distance to the nearnest TSS gene ########
  #######
  GTF_Body_table_res_cl_GR = GRanges(GTF_Body_table_res_cl$Peak)
  #######
  dis_all = c()
  #######
  for(i in 1:length(GTF_Body_table_res_cl_GR)){
    tmp = GTF_Body_table_res_cl_GR[i]
    tmp_gene = GTF_Body_table_res_cl$Gene[i]
    ######
    k = which(names(Gtfs_protein_trans_GR_TSS) == tmp_gene)
    tmp_gene_tss = Gtfs_protein_trans_GR_TSS[[k]]
    ######
    tmp_dis = distance(tmp,tmp_gene_tss)
    out = min(na.omit(tmp_dis))
    dis_all = c(dis_all,out)
    if(i %% 1000 == 0){
      print(i)
    }
  }
  #######
  #######
  GTF_Body_table_res_cl$DistanceTSS = dis_all
  #######
  return(GTF_Body_table_res_cl)
}

Process_TSS_peaks <- function(GTF_TSS_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
  GTF_TSS_table_res$Index = paste0(GTF_TSS_table_res$Gene,"#",GTF_TSS_table_res$Peak)
  #GTF_TSS_table_res$Index= gsub(":","-",GTF_TSS_table_res$Index) #to match Signac peak name style

  ###### First add PtoG link correaltions !!! ###########
  ######
  GTF_TSS_table_res$PtoG_cor = "ND"
  GTF_TSS_table_res$PtoG_pvalue = "ND"
  ######
  PtoG_index = paste0(PtoG$gene,"#",PtoG$peak)
  ######
  m = match(GTF_TSS_table_res$Index,PtoG_index)
  GTF_TSS_table_res$PtoG_cor = PtoG$score[m]
  GTF_TSS_table_res$PtoG_pvalue = PtoG$pvalue[m]
  #######
  k = which(is.na(GTF_TSS_table_res$PtoG_pvalue) == T)
  GTF_TSS_table_res$PtoG_pvalue[k] = "ND"
  GTF_TSS_table_res$PtoG_cor[k] = "ND"
  GTF_TSS_table_res_cl = GTF_TSS_table_res
  #######
  ####### Next see the distance to the nearnest TSS gene ########
  #######
  #######
  #######
  GTF_TSS_table_res_cl$DistanceTSS = 0
  #######
  return(GTF_TSS_table_res_cl)
}

Process_Inter_peaks <- function(GTF_TSSBody_table_res,PtoG,Gtfs_protein_trans_GR_TSS){
  New_table = data.frame(Gene = PtoG$gene, Peak = PtoG$peak, Class="Inter",Index=paste0(PtoG$gene,"#",PtoG$peak),PtoG_cor=PtoG$score,PtoG_pvalue=PtoG$pvalue)
  #New_table$Index = gsub("#","-",New_table$Index)
  ###### 
  ###### next remove peaks from GTF_TSSBody ######
  k = which(New_table$Peak %in% GTF_TSSBody_table_res$Peak == T)
  ######
  New_table_cl = New_table[-k,]
  ######
  ###### First add PtoG link correaltions !!! ###########
  ######
  #######
  ####### Next see the distance to the nearnest TSS gene ########
  #######
  New_table_cl_GR = GRanges(New_table_cl$Peak)
  #######
  dis_all = c()
  #######
  for(i in 1:length(New_table_cl_GR)){
    tmp = New_table_cl_GR[i]
    tmp_gene = New_table_cl$Gene[i]
    ######
    #print(tmp_gene)
    k = which(names(Gtfs_protein_trans_GR_TSS) == tmp_gene)
    if(length(k) > 0){
      tmp_gene_tss = Gtfs_protein_trans_GR_TSS[[k]]
      ######
      tmp_dis = distance(tmp,tmp_gene_tss)
      out = min(na.omit(tmp_dis))
      dis_all = c(dis_all,out)
    }
    if(length(k) == 0){
      dis_all = c(dis_all,"ND")
    }
    if(i %% 1000 == 0){
      print(i)
    }
  }
  #######
  #######
  New_table_cl$DistanceTSS = dis_all
  #######
  k = which(New_table_cl$DistanceTSS == "ND")
  #######
  New_table_cl$DistanceTSS[k] = 10000000
  return(New_table_cl)
}

Merge_TSS_and_Body <- function(GTF_TSS_table_res,GTF_Body_table_res){
  ########
  GTF_TSS_table_res$Index = paste0(GTF_TSS_table_res$Gene,"#",GTF_TSS_table_res$Peak)
  GTF_Body_table_res$Index = paste0(GTF_Body_table_res$Gene,"#",GTF_Body_table_res$Peak)
  ########
  ######## rm Body with overlapped with tss ########
  ########
  k = which(GTF_Body_table_res$Index %in% GTF_TSS_table_res$Index == T)
  GTF_Body_table_res_cl = GTF_Body_table_res[-k,]
  ########
  GTF_TSSBody_table_res = rbind(GTF_TSS_table_res,GTF_Body_table_res_cl)
  ########
  return(GTF_TSSBody_table_res)
}

# This original function will cause error when index Ks are empty
# Process_peak_table_to_annotation <- function(ALL_peak_table){
#   #######
#   k1 = which(ALL_peak_table$Class %in% c("Body","Inter") == T & ALL_peak_table$DistanceTSS > 100000)
#   ALL_peak_table_cl = ALL_peak_table[-k1,]
#   #######
#   k2 = which(ALL_peak_table_cl$Class == 'TSS' & ALL_peak_table_cl$PtoG_cor < 0)
#   #######
#   ALL_peak_table_cl2 = ALL_peak_table_cl[-k2,]
#   ######
#   k3 = which(ALL_peak_table_cl2$Class %in% c("Body","Inter") == T & ALL_peak_table_cl2$PtoG_cor < 0)
#   ALL_peak_table_cl3 = ALL_peak_table_cl2[-k3,]
#   ######
#   ######
#   ALL_peak_table_cl3$Index_Peak = paste0(ALL_peak_table_cl3$Gene,"::",ALL_peak_table_cl3$Class)
#   ##### ALL_peak_table_cl3[which(ALL_peak_table_cl3$Gene == 'mmp9'),]
#   #####
#   ##### add tags for the peak #####
#   ALL_peak_table_cl4 = split(ALL_peak_table_cl3,ALL_peak_table_cl3$Index_Peak)
#   #####
#   add_index <- function(x){
#     x$Index_Peak2 = paste0(x$Index_Peak,"_",1:dim(x)[1])
#     return(x)
#   }
#   #####
#   ALL_peak_table_cl4 = lapply(ALL_peak_table_cl4,add_index)
#   ALL_peak_table_cl4 = do.call("rbind",ALL_peak_table_cl4)
#   #####
#   return(ALL_peak_table_cl4)
# }


Process_peak_table_to_annotation <- function(ALL_peak_table){
  #######
  ALL_peak_table_cl <- subset(ALL_peak_table, !(Class %in% c("Body","Inter") & DistanceTSS > 100000))

  ALL_peak_table_cl2 <- subset(ALL_peak_table_cl, !(Class == 'TSS' & PtoG_cor < 0))
  ######
  ALL_peak_table_cl3 = subset(ALL_peak_table_cl2, !(Class %in% c("Body","Inter") & PtoG_cor < 0))
  ######
  ######
  ALL_peak_table_cl3$Index_Peak = paste0(ALL_peak_table_cl3$Gene,"::",ALL_peak_table_cl3$Class)
  ##### ALL_peak_table_cl3[which(ALL_peak_table_cl3$Gene == 'mmp9'),]
  #####
  ##### add tags for the peak #####
  ALL_peak_table_cl4 = split(ALL_peak_table_cl3,ALL_peak_table_cl3$Index_Peak)
  #####
  add_index <- function(x){
    x$Index_Peak2 = paste0(x$Index_Peak,"_",1:dim(x)[1])
    return(x)
  }
  #####
  ALL_peak_table_cl4 = lapply(ALL_peak_table_cl4,add_index)
  ALL_peak_table_cl4 = do.call("rbind",ALL_peak_table_cl4)
  #####
  return(ALL_peak_table_cl4)
}
