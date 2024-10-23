Convert_to_triple_and_filterbygenes <- function(GRNs,Genes_need){
  ###########
  GRNs = GRNs[which(GRNs$PtoG_cor!=0),] 
  GRNs_cl = GRNs[,c('TF','Motif','Peak','Target','TF_Motif_Cor','PtoG_Class','GG_Corr')]
  ###########
  GRNs_cl_Triple = paste0(GRNs_cl$TF,"->",GRNs_cl$Peak,"->",GRNs_cl$Target)
  ###########
  GRNs_cl$Triple = GRNs_cl_Triple
  GRNs_cl_Triple_cl = GRNs_cl[!duplicated(GRNs_cl$Triple),]
  ###########
  k = which(GRNs_cl_Triple_cl$TF %in% Genes_need == T & GRNs_cl_Triple_cl$Target %in% Genes_need== T)
  ###########
  GRNs_cl_Triple_clcl = GRNs_cl_Triple_cl[k,]
  ############
  print(dim(GRNs_cl_Triple_clcl))
  return(GRNs_cl_Triple_clcl)
}


Convert_to_double_and_filterbygenes <- function(GRNs,Genes_need){
  ###########
  GRNs_cl = GRNs[,c('TF','Motif','Peak','Target','TF_Motif_Cor','PtoG_Class','GG_Corr')]
  ###########
  GRNs_cl_Double = paste0(GRNs_cl$TF,"->",GRNs_cl$Target)
  ###########
  GRNs_cl$Double = GRNs_cl_Double
  GRNs_cl_Double_cl = GRNs_cl[!duplicated(GRNs_cl$Double),]
  ###########
  k = which(GRNs_cl_Double_cl$TF %in% Genes_need == T & GRNs_cl_Double_cl$Target %in% Genes_need == T)
  ###########
  GRNs_cl_Double_clcl = GRNs_cl_Double_cl[k,]
  ############
  print(dim(GRNs_cl_Double_clcl))
  return(GRNs_cl_Double_clcl)
}


Filter_triple_network_according_to_foldchange <- function(GRNs_MG_triple,DEGs,DARs,all_peak_log2fc){
  colnames(DEGs) = c("avg_log2FC","gene_short")
  colnames(DARs) = c("Log2FC","peaks")


  ######## ad activator/repressor tag ######
  GRNs_MG_triple$type = 'Act'
  k = which(GRNs_MG_triple$GG_Corr < 0)
  GRNs_MG_triple$type[k] = 'Rep'
  ########
  ######## for Activator/repressor #######
  ########
  GRNs_MG_triple_Act = GRNs_MG_triple[which(GRNs_MG_triple$type == 'Act'),]
  GRNs_MG_triple_Rep = GRNs_MG_triple[which(GRNs_MG_triple$type == 'Rep'),]

  ######## for Activator :Activtor with positive TF-gene correlation,#######
  genes_1 = DEGs$gene_short[which(DEGs$avg_log2FC < 0)]#downregu in regeneration
  genes_2 = DEGs$gene_short[which(DEGs$avg_log2FC > 0)]#upregulated in regeneration
  peaks = DARs$peaks[which(DARs$Log2FC < 0)]
  # select upregu TF and target, also the peak is not in closing trend
  k1 = which(GRNs_MG_triple_Act$TF %in% genes_2 == T & GRNs_MG_triple_Act$Target %in% genes_2 == T & GRNs_MG_triple_Act$Peak %in% peaks == F)
  GRNs_MG_triple_Act_cl = GRNs_MG_triple_Act[k1,]
  
  
  ######## for Repressor: Repressor with negative TF-gene correlation #######
  genes_1 = DEGs$gene_short[which(DEGs$avg_log2FC < 0)] #downregu in regeneration
  genes_2 = DEGs$gene_short[which(DEGs$avg_log2FC > 0)] #upregulated in regeneration
  peaks = DARs$peaks[which(DARs$Log2FC < 0)]
  # upregu TF and downregu target, and peak is not in closing trend
  k2 = which(GRNs_MG_triple_Rep$TF %in% genes_1 == T & GRNs_MG_triple_Rep$Target %in% genes_2 == T & GRNs_MG_triple_Rep$Peak %in% peaks == F)
  GRNs_MG_triple_Rep_cl = GRNs_MG_triple_Rep[k2,]
  #######
  #######
  DEGs_pos = DEGs[which(DEGs$avg_log2FC > 0),]
  DEGs_neg = DEGs[which(DEGs$avg_log2FC < 0),]
  Peak_pos = DARs[which(DARs$Log2FC > 0),]
  ########
  ######## for GRNs_MG_triple_Act_cl #####
  GRNs_MG_triple_Act_cl$TF_fold = sapply(as.list(GRNs_MG_triple_Act_cl$TF),Get_highest_score_G,tab=DEGs_pos)
  GRNs_MG_triple_Act_cl$Target_fold = sapply(as.list(GRNs_MG_triple_Act_cl$Target),Get_highest_score_G,tab=DEGs_pos)
  GRNs_MG_triple_Act_cl$Peak_fold = sapply(as.list(GRNs_MG_triple_Act_cl$Peak),Get_highest_score_P,tab=Peak_pos)
  GRNs_MG_triple_Act_cl$Peak_fold_2 = all_peak_log2fc[,1][match(GRNs_MG_triple_Act_cl$Peak,all_peak_log2fc$peaks)]
  ########
  GRNs_MG_triple_Rep_cl$TF_fold = sapply(as.list(GRNs_MG_triple_Rep_cl$TF),Get_lowest_score_G,tab=DEGs_neg)
  GRNs_MG_triple_Rep_cl$Target_fold = sapply(as.list(GRNs_MG_triple_Rep_cl$Target),Get_highest_score_G,tab=DEGs_pos)
  GRNs_MG_triple_Rep_cl$Peak_fold = sapply(as.list(GRNs_MG_triple_Rep_cl$Peak),Get_highest_score_P,tab=Peak_pos)
  GRNs_MG_triple_Rep_cl$Peak_fold_2 = all_peak_log2fc[,1][match(GRNs_MG_triple_Rep_cl$Peak,all_peak_log2fc$peaks)]

  ########
  GRNs_MG_triple_cl_out = rbind(GRNs_MG_triple_Act_cl,GRNs_MG_triple_Rep_cl)
  ########
  return(GRNs_MG_triple_cl_out)
}

Get_highest_score_G <- function(x,tab){
  k = which(tab$gene_short == x)
  #####
  if(length(k) == 0){
    out = 0
    return(out)
  }
  if(length(k) > 0){
    out = max(tab$avg_log2FC[k])
    return(out)
  }
}


Get_highest_score_P <- function(x,tab){
  k = which(tab$peaks == x)
  #####
  if(length(k) == 0){
    out = 0
    return(out)
  }
  if(length(k) > 0){
    out = max(tab$Log2FC[k])
    return(out)
  }
}


Get_lowest_score_G <- function(x,tab){
  k = which(tab$gene_short == x)
  #####
  if(length(k) == 0){
    out = 0
    return(out)
  }
  if(length(k) > 0){
    out = min(tab$avg_log2FC[k])
    return(out)
  }
}