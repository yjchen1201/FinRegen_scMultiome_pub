# TF activity TF gene expression correlation 
Get_corr_gene_motif <- function(motif_matrix,gene_matrix,Fish_All_motifs_table_short){

  ########
  k = which(Matrix::rowSums(gene_matrix) < 10)
  if(length(k) > 1){
    gene_matrix = gene_matrix[-k,] 
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


Get_all_need_Motif_tag <- function(x,cutoff=0.05){
  ########
  k1 = which(x$Gene_Motif_Cor > cutoff | x$Gene_Motif_Cor < -cutoff)
  #########
  x_cl = x[k1,]
  #########
  k3 = which(x_cl$Gene_Motif_Cor > cutoff)
  k4 = which(x_cl$Gene_Motif_Cor < -cutoff)
  x_cl$class = 'Unknown'
  x_cl$class[k3] = 'pos'
  x_cl$class[k4] = 'neg'
  #########
  print(length(x_cl$ID[!duplicated(x_cl$ID)]))
  x_cl[which(x_cl$Gene == 'ascl1a'),]
  return(x_cl)
}

Get_all_need_Peak_tag <- function(x1,x2,x3,x4){
  all_tab = rbind(x1,x2,x3,x4)
  ######
  all_peak = all_tab$Peak
  all_peak = all_peak[!duplicated(all_peak)]
  ######
  return(all_peak)
}


Match_motifs <- function(Peak_with_tag,Motif_with_tag,Fish_combined){
  ############
  Peak_with_tag_1 = Merge_peaks(Peak_with_tag,extend=150)
  ############
  Peak_with_tag_2 = Merge_GRlist(Peak_with_tag_1)
  ############ filter the Fish_combined #####
  k = which(names(Fish_combined) %in% Motif_with_tag$ID == T)
  Fish_combined_cl = Fish_combined[k]
  ############
  library(parallel)
  cl <- makeCluster(35)
  ############
  footprint_Motif_list = parLapply(cl,Fish_combined_cl,Find_Motif,Peak_with_tag_2=Peak_with_tag_2)
  stopCluster(cl)
  ############
  head(names(footprint_Motif_list))
  ###########
  return(footprint_Motif_list)
  ###########
}

Find_Motif <- function(x,Peak_with_tag_2){
  library(motifmatchr)
  library("BSgenome.Drerio.UCSC.danRer11")
  footprint = matchMotifs(x,Peak_with_tag_2,genome = BSgenome.Drerio.UCSC.danRer11,out='positions',p.cutoff = 5e-05)
  ######
  footprint_tab = data.frame(footprint[[1]])
  ######
  return(footprint_tab)
}




Combined_all_things <- function(GG_network,PG_cor,TP_cor,Foot){
  library(GenomicRanges)
  ##### first filtr the GG_network #######
  ##### keep 95% intervals ###############
  GG_network = GG_network[order(GG_network$Score,decreasing=T),]
  GG_cutoff = quantile(GG_network$Score,0.9)
  #####
  GG_network_cl = GG_network[which(GG_network$Score > GG_cutoff),]
  #####
  # filter abs(correlation) below 0.03
  k1 = which(GG_network_cl$Corr > 0.03)
  k2 = which(GG_network_cl$Corr < -0.03)
  #####
  GG_network_cl$Class = 'unknown'
  GG_network_cl$Class[k1] = 'pos'
  GG_network_cl$Class[k2] = 'neg'
  #####
  GG_network_clcl = GG_network_cl[which(GG_network_cl$Class != 'unknown'),]
  GG_network_clcl$TF_Gene = paste0(GG_network_clcl$TF,'->',GG_network_clcl$Gene)
  GG_network_clclcl = GG_network_clcl[,c('TF_Gene','Score','Corr')]
  colnames(GG_network_clclcl) = c('TF_Gene','GG_Score','GG_Corr')
  #####
  ## head(GG_network_clcl[which(GG_network_clcl$TF == 'foxj1a'),])
  #####
  ##### Next filter the footprint #####
  #####
  k_foot <- which(c(Foot$left_delta + Foot$right_delta) > 0.1)
  ######
  Foot_cl = Foot[k_foot]
  ######
  ###### Next link footprint to Peaks !!! ###########
  ######
  Footprint_Peak_table = All_footprint_to_peaks(Foot_cl,PG_cor)
  ###### Next merge the footpint peak to motif gene list !!! ##########
  ######
  colnames(TP_cor) = c("Motif","TF","TF_Motif_Cor","TAG") #modified from orignal code
  TP_cor_cl = TP_cor[which(TP_cor$TF_Motif_Cor > 0.05 | TP_cor$TF_Motif_Cor < -0.05),]
  TP_cor_cl = TP_cor_cl[,c("TF","Motif","TF_Motif_Cor")]
  ###### Next we will merge big table ########
  library(data.table)
  library(dplyr)
  Footprint_Peak_table = data.table(Footprint_Peak_table)
  TP_cor_cl = data.table(TP_cor_cl)
  merged_table <- inner_join(Footprint_Peak_table, TP_cor_cl, by = "Motif",relationship = "many-to-many")
  ###### test merge ##### correct !! #
  #Footprint_Peak_table = data.frame(Footprint_Peak_table)
  #TP_cor_cl = data.frame(TP_cor_cl)
  #merged_table_test <- merge(Footprint_Peak_table, TP_cor_cl)
  #dim(merged_table_test)
  ####### then next we will merge TF_gene corrections ######
  merged_table$TF_Gene = paste0(merged_table$TF,"->",merged_table$Target)
  ######
  GG_network_clclcl <- data.table(GG_network_clclcl)
  ######
  merged_table2 = inner_join(merged_table, GG_network_clclcl, by = "TF_Gene")
  ######
  k_filter1 = which(merged_table2$TF_Motif_Cor > 0 & merged_table2$GG_Corr > 0)
  k_filter2 = which(merged_table2$TF_Motif_Cor < 0 & merged_table2$GG_Corr < 0)
  #######
  merged_table3 = merged_table2[c(k_filter1,k_filter2),]
  #######
  print(dim(merged_table3))
  return(merged_table3)
}

All_footprint_to_peaks <- function(footprint_GR,peak_table){
  
  peak_table_GR = GRanges(peak_table$Peak)
  start(peak_table_GR) = start(peak_table_GR)-0  #original script is to extend 150 up and down,so total peak size will be 801bp for all peak, ours is different, since we didnt define peak region based on summit
  end(peak_table_GR) = end(peak_table_GR)+0
  #######
  #######
  ####### Next we extend the peak_table #######
  ####### find the overlaps with peaks and footprint ####
  ####### left is peak right is the footprint ###
  countIndex = findOverlaps(peak_table_GR,footprint_GR)
  ####### 
  left_table = peak_table[queryHits(countIndex),]
  #######
  right_table = data.frame(footprint_GR)[subjectHits(countIndex),]
  ####### Clean_left_table ########
  left_table_cl = left_table[,c("Index","Peak","Class","PtoG_cor","Gene")]
  colnames(left_table_cl) <- c("PtoG_Index","Peak","PtoG_Class","PtoG_cor","Target")
  ### 
  k = which(left_table_cl$PtoG_cor == 'ND')
  left_table_cl$PtoG_cor[k] = 0
  left_table_cl$PtoG_cor = round(as.numeric(left_table_cl$PtoG_cor),3)
  ####### head(right_table)
  right_table$footprint = paste0(right_table$seqnames,':',right_table$start,'-',right_table$end)
  right_table_cl = right_table[,c('motif','footprint','mid_score','left_delta','right_delta')]
  colnames(right_table_cl) <- c("Motif",'Footprint_region','Footprint_Mid',"Footprint_Left_del","Footprint_Right_del")
  #######
  merge_table = cbind(left_table_cl,right_table_cl)
  #######
  return(merge_table)
}


