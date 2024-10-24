Get_all_need_Motif_tag <- function(x,cutoff=0.05){
  #x = rbind(x,y)
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
  #x_cl[which(x_cl$Gene == 'ascl1a'),]
  return(x_cl)
}


Get_all_need_Peak_tag <- function(x1){
  ######
  all_peak = x1$Peak
  all_peak = all_peak[!duplicated(all_peak)]
  ######
  return(all_peak)
}


# this is modified from original Match_motifs function
## I already did motif finding using whole peak set, so just extract motif position from this file 
Match_motifs <- function(Globle_need_Peak,Globle_need_Motif,Fish_combined){
  ############
  ############ filter the Fish_combined #####
  # for all called peaks, i did motif match already
  all_peak_motif_match =readRDS("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step4-5_All_peak_motif_binding_region/Total_footprint_Motif_in_WNNcelltype_calledpeaks.rds")
  
  k = which(names(Fish_combined) %in% Globle_need_Motif$ID == T)
  ## fiter by motifs
  footprint_Motif_list = all_peak_motif_match[names(Fish_combined)[k]]
  ## filter by cell type linked peaks
  Global_need_Peak_gr = Signac::StringToGRanges(Globle_need_Peak,sep = c(":", "-"))

  overlapping_motifs <- list()
  for (i in seq_along(footprint_Motif_list)) {
      # Find overlaps
      overlap_hits <- findOverlaps(footprint_Motif_list[[i]], Global_need_Peak_gr)
      # Check if there are any overlaps
      if (length(queryHits(overlap_hits)) > 0) {
          # If overlaps exist, extract them
          overlaps <- footprint_Motif_list[[i]][queryHits(overlap_hits)]
          overlaps = data.frame(sort(overlaps))
          overlapping_motifs[[names(footprint_Motif_list)[i]]] <-  overlaps

      }
      # If no overlaps, the loop continues to the next iteration automatically
  }
  overlapping_motifs ->footprint_Motif_list
  head(names(footprint_Motif_list))
  ###########
  return(footprint_Motif_list)
  ###########
}



Match_motifs <- function(Globle_need_Peak,Globle_need_Motif,Fish_combined){
  ############
  ############ filter the Fish_combined #####
  # for all called peaks, i did motif match already
  all_peak_motif_match =readRDS("/scratch/ichen/IReNA/IReNA-v2/use_Signac_p2g/use_Signac_peaklink_JASPAR_motif/step4-5_All_peak_motif_binding_region/Total_footprint_Motif_in_WNNcelltype_calledpeaks.rds")
  
  k = which(names(Fish_combined) %in% Globle_need_Motif$ID == T)
  ## fiter by motifs
  footprint_Motif_list = all_peak_motif_match[names(Fish_combined)[k]]
  ## filter by cell type linked peaks
  Global_need_Peak_gr = Signac::StringToGRanges(Globle_need_Peak,sep = c(":", "-"))

  overlapping_motifs <- list()
  for (i in seq_along(footprint_Motif_list)) {
      # Find overlaps
      overlap_hits <- findOverlaps(footprint_Motif_list[[i]], Global_need_Peak_gr)
      # Check if there are any overlaps
      if (length(queryHits(overlap_hits)) > 0) {
          # If overlaps exist, extract them
          overlaps <- footprint_Motif_list[[i]][queryHits(overlap_hits)]
          overlaps = data.frame(sort(overlaps))
          overlapping_motifs[[names(footprint_Motif_list)[i]]] <-  overlaps

      }
      # If no overlaps, the loop continues to the next iteration automatically
  }
  overlapping_motifs ->footprint_Motif_list
  head(names(footprint_Motif_list))
  ###########
  return(footprint_Motif_list)
  ###########
}



Get_footprint_Score_all_para <- function(Globle_motif_footprint,footprint_signal){
  library(GenomicRanges)
  chr = sapply(Globle_motif_footprint,function(x) x$seqnames)
  start = sapply(Globle_motif_footprint,function(x) x$start)
  end = sapply(Globle_motif_footprint,function(x) x$end)
  score = sapply(Globle_motif_footprint,function(x) x$score)
  #########
  len = sapply(chr,function(x) length(x))
  names_tag = rep(names(len),len)
  #########
  chr =do.call(c,chr)
  start=do.call(c,start)
  end=do.call(c,end)
  score=do.call(c,score)
  #########
  new_GR = GRanges(seqnames=chr,IRanges(start,end),score=score,motif=names_tag)
  seqlevels(new_GR) <- paste0("chr", seqlevels(new_GR))
  #########
  x_GRlist <- Get_left_mid_right_region_list(new_GR)
  x_left = cal_footprint_score_sub1(x_GRlist[[1]],footprint_signal)
  x_mid = cal_footprint_score_sub1(x_GRlist[[2]],footprint_signal)
  x_right = cal_footprint_score_sub1(x_GRlist[[3]],footprint_signal)
  #########
  x_res = Merge_footprint_score(x_left,x_mid,x_right)
  #########
  return(x_res)
}

cal_footprint_score_sub1 <- function(tmp_GR,tmp_F){
  library(GenomicRanges)
  ####
  res = findOverlaps(tmp_GR,tmp_F)
  ####
  score = tmp_F$score[subjectHits(res)]
  ####
  res_score_merge = tapply(score,queryHits(res),sum)
  ####
  tmp_GR$footprint_score = 0
  ####
  tmp_GR$footprint_score[as.numeric(names(res_score_merge))] = as.numeric(res_score_merge)
  ####
  tmp_GR$wid = width(tmp_GR)
  ####
  tmp_GR$footprint_score_norm = tmp_GR$footprint_score / tmp_GR$wid
  ####
  return(tmp_GR)
}

Get_left_mid_right_region_list <- function(x){
  library(GenomicRanges)
  #######
  x_L = x
  x_R = x
  x_M = x
  #######
  width = width(x)
  ######
  start(x_L) = start(x_L) - width*3
  end(x_L) = end(x_L) - width
  ########
  start(x_R) = start(x_R) + width
  end(x_R) = end(x_R) + width*3
  ########
  #######
  return(list(x_L,x,x_R))
}

Merge_footprint_score <- function(x_L,x_M,x_R){
  ###
  #############
  x = x_M
  #############
  x$left_score = 0
  x$mid_score = 0
  x$right_score = 0
  ##############
  x$left_delta = 0
  x$right_delta = 0
  #############
  x$left_score = x_L$footprint_score_norm
  x$mid_score = x_M$footprint_score_norm
  x$right_score = x_R$footprint_score_norm
  #############
  x$left_delta = x$left_score - x$mid_score
  x$right_delta = x$right_score - x$mid_score
  #############
  ## summary(x$mid_score)
  return(x)
  ###
}



