#' Color scales
#' @description
#' `borders_def` Determines the color scale for each modality.
#' @param node_list character vector indicating the complete set of nodes 
#' in the resulting network structure.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' @param omics_meth_original METH matrix containing original beta values
#' Each component of the list is a matrix with samples in rows and features 
#' in columns.
#' @importFrom RColorBrewer brewer.pal
#' @return List of 5 elements indicating the color scale for each modality
borders_def <- function(node_list, layers_def, omics, omics_meth_original)
{
  ge_cols <- brewer.pal(9, "Blues")
  ge_mat <- omics[[layers_def$omics[1]]]
  ge_common <- intersect(unique(node_list), colnames(ge_mat))
  omics_ge_gs <- as.matrix(ge_mat[,ge_common])
  colnames(omics_ge_gs) <- ge_common
  
  ge_mat_cut_med_low <- cut(ge_mat[ge_mat<=median(ge_mat)], 
    seq(from=min(ge_mat), to=median(ge_mat), length.out=5), 
    include.lowest = TRUE)
  ge_mat_cut_med_high <- cut(ge_mat[ge_mat>median(ge_mat)], 
    seq(from=median(ge_mat), to=max(ge_mat), length.out=6), 
    include.lowest = TRUE)
  
  borders_ge_b1 <- unlist(lapply(strsplit(levels(ge_mat_cut_med_low),","),
    FUN=function(l) l[1]))
  borders_ge_b1[1] <- sub("[","(",borders_ge_b1[1], fixed = TRUE)
  borders_ge_b1 <- as.numeric(sub("(","",borders_ge_b1, fixed = TRUE))
  borders_ge_b2 <- unlist(lapply(strsplit(levels(ge_mat_cut_med_high),",")
    ,FUN=function(l) l[1]))
  borders_ge_b2[1] <- sub("[","(",borders_ge_b2[1], fixed = TRUE)
  borders_ge_b2 <- as.numeric(sub("(","",borders_ge_b2, fixed = TRUE))
  borders_ge_b <- c(borders_ge_b1,borders_ge_b2)
  
  borders_ge_t1 <- as.numeric(sub("]","", 
    unlist(lapply(strsplit(levels(ge_mat_cut_med_low),","),
    FUN=function(l) l[2]))))
  borders_ge_t2 <- as.numeric(sub("]","", 
    unlist(lapply(strsplit(levels(ge_mat_cut_med_high),","),
    FUN=function(l) l[2]))))
  borders_ge_t <- c(borders_ge_t1,borders_ge_t2)
  borders <- sort(unique(c(borders_ge_b,borders_ge_t)))
  expr_group <- cut(colMeans(omics_ge_gs), breaks = borders,
                    include.lowest = TRUE, labels = FALSE)
  names(expr_group) <- colnames(omics_ge_gs)
  node_list <- matrix(data = c(node_list, 
    as.numeric(expr_group[match(node_list, names(expr_group))])), 
      nrow = length(node_list), 
      dimnames = list(c(), c("label", "color")))
  ind_cols <- paste(paste("(", paste(borders[-length(borders)],
    borders[-1]),sep=""),"]",sep="")
  borders_cnv <- NULL
  borders_meth <- NULL
  
  if(any(mapply(FUN=function(mod)
    any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))
  {
    cnv_cols <- brewer.pal(11, "PiYG")
    cnv_common <- intersect(node_list[,"label"][regexpr("eid",
      node_list[,"label"])>0],
      colnames(omics[[names(which(mapply(FUN=function(mod)
      any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]))
    omics_cnv_gs <- as.matrix(omics[["cnv"]][,cnv_common])
    cnv_mat_cut_neg <- cut(omics[["cnv"]][omics[["cnv"]]<=0], 
      seq(from=min(omics[["cnv"]], na.rm = TRUE), to=0, length.out=6), 
      include.lowest = TRUE)
    cnv_mat_cut_pos <- cut(omics[["cnv"]][omics[["cnv"]]>0], 
      seq(from=0, to=max(omics[["cnv"]], na.rm = TRUE), length.out=7), 
      include.lowest = TRUE)
  
    borders_cnv_b1 <- unlist(lapply(strsplit(levels(cnv_mat_cut_neg),","),
      FUN=function(l) l[1]))
    borders_cnv_b1[1] <- sub("[","(",borders_cnv_b1[1], fixed = TRUE)
    borders_cnv_b1 <- as.numeric(sub("(","",borders_cnv_b1, fixed = TRUE))
    borders_cnv_b2 <- unlist(lapply(strsplit(levels(cnv_mat_cut_pos),","),
      FUN=function(l) l[1]))
    borders_cnv_b2[1] <- sub("[","(",borders_cnv_b2[1], fixed = TRUE)
    borders_cnv_b2 <- as.numeric(sub("(","",borders_cnv_b2, fixed = TRUE))
    borders_cnv_b <- c(borders_cnv_b1,borders_cnv_b2)
    
    borders_cnv_t1 <- as.numeric(sub("]","", 
    unlist(lapply(strsplit(levels(cnv_mat_cut_neg),","),
      FUN=function(l) l[2]))))
    borders_cnv_t2 <- as.numeric(sub("]","", 
    unlist(lapply(strsplit(levels(cnv_mat_cut_pos),","),
      FUN=function(l) l[2]))))
    borders_cnv_t <- c(borders_cnv_t1,borders_cnv_t2)
    borders_cnv <- sort(unique(c(borders_cnv_b,borders_cnv_t)))
    
    cnv_group <- cut(colMeans(omics_cnv_gs, na.rm = TRUE), 
      breaks = borders_cnv, include.lowest = TRUE, labels = FALSE) + 
      length(ge_cols)
    names(cnv_group) <- colnames(omics_cnv_gs)
    node_list[regexpr("eid",node_list[,"label"])>0,"color"] <-
      as.numeric(cnv_group[match(node_list[regexpr("eid", 
      node_list[,"label"])>0,"label"], names(cnv_group))])
    
    ge_cols <- c(ge_cols, cnv_cols)
  } # end if(any(mapply(FUN=function(mod)...
  
  if(any(mapply(omics,
    FUN=function(list) any(regexpr("eid:",colnames(list), 
    ignore.case = TRUE)<0))))
  {
    meth_cols <- brewer.pal(9, "YlOrRd")
    meth_common <-
      intersect(node_list[,"label"],colnames(omics_meth_original))
    omics_meth_gs <- as.matrix(omics_meth_original[,meth_common])
    colnames(omics_meth_gs) <- meth_common
    
    meth_mat_cut_05_low <- cut(omics_meth_original[omics_meth_original<=0.5], 
      seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5, length.out=5), 
      include.lowest = TRUE)
    meth_mat_cut_05_high <- cut(omics_meth_original[omics_meth_original>0.5], 
      seq(from=0.5, to=max(omics_meth_original, na.rm = TRUE), length.out=6), 
      include.lowest = TRUE)
    
    borders_meth_b1 <- unlist(lapply(strsplit(levels(meth_mat_cut_05_low),","), 
      FUN=function(l) l[1]))
    borders_meth_b1[1] <- sub("[","(",borders_meth_b1[1], fixed = TRUE)
    borders_meth_b1 <- as.numeric(sub("(","",borders_meth_b1, fixed = TRUE))
    borders_meth_b2 <- unlist(lapply(strsplit(levels(meth_mat_cut_05_high),","),
      FUN=function(l) l[1]))
    borders_meth_b2[1] <- sub("[","(",borders_meth_b2[1], fixed = TRUE)
    borders_meth_b2 <- as.numeric(sub("(","",borders_meth_b2, fixed = TRUE))
    borders_meth_b <- c(borders_meth_b1,borders_meth_b2)
    borders_meth_t1 <- as.numeric(sub("]","", 
      unlist(lapply(strsplit(levels(meth_mat_cut_05_low),","),
      FUN=function(l) l[2]))))
    borders_meth_t2 <- as.numeric(sub("]","", 
      unlist(lapply(strsplit(levels(meth_mat_cut_05_high),","),
      FUN=function(l) l[2]))))
    borders_meth_t <- c(borders_meth_t1,borders_meth_t2)
    borders_meth <- sort(unique(c(borders_meth_b,borders_meth_t)))
    meth_group <- cut(colMeans(omics_meth_gs, na.rm = TRUE), 
      breaks = borders_meth, include.lowest = TRUE, labels = FALSE) + 
      length(ge_cols)
    names(meth_group) <- colnames(omics_meth_gs)
    node_list[!is.na(match(node_list[,"label"], 
      colnames(omics_meth_original))),"color"] <- 
      as.numeric(meth_group[match(node_list[,"label"][!is.na(match(node_list[,
      "label"], colnames(omics_meth_original)))], names(meth_group))])
    ge_cols <- c(ge_cols, meth_cols)
  } # end if(any(mapply(omics,FUN=function(list)...
  return(list(borders=borders,
              borders_meth=borders_meth,
              borders_cnv=borders_cnv,
              node_list=node_list,
              ge_cols=ge_cols))
}
