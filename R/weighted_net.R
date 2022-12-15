#' Resulting network definition
#' @description
#' `weighted_net` Defines the resulting network structure and determines 
#' the color scale for each modality.
#' @param cpdag_weights data.frame output from the edge_weights function.
#' @param B_prior_mat_weighted matrix one of the outputs 
#' of the bn_module function.
#' @param gene_annot data.frame containing the entrez ID and corresponding gene
#' symbol for conversion.
#' @param PK data.frame with known interactions.
#' @param OMICS_mod_res list output from the omics_module function.
#' @param edge_weights character vector includes either "MCMC_freq" to reflect
#' the edge weights frequency over the final set of network structures or 
#' "empB" to reflect the empirical biological knowledge estimated by IntOMICS.
#' @param gene_ID character vector includes either "gene_symbol" or "entrezID"
#' to reflect gene identifiers used in the final figure.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @importFrom igraph degree
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph as_ids
#' @importFrom methods is
#' 
#' @examples
#' data(list=c("OMICS_mod_res", "BN_mod_res", "gene_annot", "TFtarg_mat", 
#' "PK"), package="IntOMICS")
#' res_weighted <- edge_weights(mcmc_res = BN_mod_res, burn_in = 10000, 
#'  thin = 500, edge_freq_thres = 0.3) 
#' weighted_net_res <- weighted_net(cpdag_weights = res_weighted, 
#'  gene_annot = gene_annot, PK = PK, OMICS_mod_res = OMICS_mod_res, 
#'  gene_ID = "gene_symbol", TFtargs = TFtarg_mat,
#'  B_prior_mat_weighted = B_prior_mat_weighted(BN_mod_res)) 
#'
#' @return List of 7 elements needed to plot the final regulatory network
#' @export
weighted_net <- function(cpdag_weights, gene_annot, PK=NULL, OMICS_mod_res, 
                         edge_weights = "MCMC_freq", gene_ID, TFtargs = NULL,
                         B_prior_mat_weighted)
{
  if(!is.data.frame(cpdag_weights) | !all(colnames(cpdag_weights) %in% 
                                          c("edge", "from.x", "to.x", 
                                          "strength.x", "direction.x", 
                                          "from.y", "to.y", "strength.y",
                                          "direction.y", "strength")))
  {
    message('Invalid input "cpdag_weights". Must be data.frame with colnames 
          c("edge", "from.x", "to.x", "strength.x", "direction.x", "from.y", 
          "to.y", "strength.y", "direction.y", "strength").')  
  }
  
  if(!is.data.frame(gene_annot) | !all(colnames(gene_annot) %in% c("entrezID","gene_symbol")) | 
     !all(regexpr("EID:",gene_annot$entrezID, fixed = TRUE)==1))
  {
    message('Invalid input "gene_annot". Must be data.frame with colnames 
            c("entrezID","gene_symbol") and entrezID must 
            be in EID:XXXX format indicating Entrez IDs.')
  }
  
  if(!is.data.frame(PK) | !all(colnames(PK) %in% c("src_entrez","dest_entrez",
                                                   "edge_type")) | 
     !all(regexpr("EID:",c(PK$src_entrez,PK$dest_entrez), fixed = TRUE)==1))
  {
    message('Invalid input "PK". Must be data.frame with colnames 
            c("src_entrez","dest_entrez","edge_type") and gene names must 
            be in EID:XXXX format indicating Entrez IDs.')
  }
  
  if(!is.list(OMICS_mod_res) | 
     !all(colnames(OMICS_mod_res) %in% c("pf_UB_BGe_pre", "B_prior_mat", 
                                         "annot", "omics", "layers_def", "omics_meth_original")))
  {
    message('Invalid input "OMICS_mod_res". Must be named list, 
            output from omics_module().')
  } 
  
  if(!(edge_weights %in% c("MCMC_freq","empB")))
  {
    message('Invalid input "edge_weights". Must be either "MCMC_freq" or "empB"')  
  }
  
  if(!gene_ID %in% c("entrezID","gene_symbol"))
  {
    message('Invalid input "gene_ID". Must be either "entrezID" or "gene_symbol".')
  }
  
  if(!is.matrix(TFtargs) | 
     !all(regexpr("EID:",unlist(dimnames(TFtargs)), fixed = TRUE)==1))
  {
    message('Invalid input "TFtargs". Must be matrix and dimnames must 
            be in EID:XXXX format indicating Entrez IDs.')
  }
  
  if(!is.matrix(B_prior_mat_weighted) | 
     !any(regexpr("EID:",unlist(dimnames(B_prior_mat_weighted)), fixed = TRUE)==1))
  {
    message('Invalid input "B_prior_mat_weighted". Must be matrix and dimnames must 
            be in correct format indicating Entrez IDs (and METH probes if available).')
  }
  
  if(!is.null(PK))
  {
    PK <- PK[PK$src_entrez %in% unlist(lapply(OMICS_mod_res$omics,colnames)),]
    PK <- PK[PK$dest_entrez %in% unlist(lapply(OMICS_mod_res$omics,colnames)),]
  }
  
  edge_list <- matrix(data = c(cpdag_weights$from.x, cpdag_weights$to.x,
                               cpdag_weights$strength, 
                               rep(NA,length(cpdag_weights$strength)),
                               rep(NA,length(cpdag_weights$strength))), 
                      nrow = length(cpdag_weights$strength), 
                      dimnames = list(c(),
                                      c("from", "to", "weight", "edge_type", "edge")))
  node_list <- unique(c(edge_list[,"from"], edge_list[,"to"]))
  edge_list[,"edge"] <- paste(edge_list[,"from"], 
                              edge_list[,"to"], sep="_")
  return_list <- edge_types(B_prior_mat_weighted = B_prior_mat_weighted, 
                            PK = PK, gene_annot = gene_annot, 
                            edge_list = edge_list, node_list = node_list, 
                            OMICS_mod_res = OMICS_mod_res, 
                            edge_weights = edge_weights, TFtargs = TFtargs)
  
  if(gene_ID=="gene_symbol")
  {
    # rename GE edges to gene symbols
    return_list$edge_list[return_list$edge_list[,"from"] %in% 
                            gene_annot$entrezID,"from"] <- 
      gene_annot$gene_symbol[match(return_list$edge_list[,"from"], 
                                   gene_annot$entrezID, nomatch = 0)]
    return_list$edge_list[return_list$edge_list[,"to"] %in% 
                            gene_annot$entrezID,"to"] <- 
      gene_annot$gene_symbol[match(return_list$edge_list[,"to"], 
                                   gene_annot$entrezID, nomatch = 0)]
    # rename CNV edges to gene symbols
    return_list$edge_list[toupper(return_list$edge_list[,"from"]) %in% 
                            gene_annot$entrezID,"from"] <- 
      tolower(gene_annot$gene_symbol[
        match(toupper(return_list$edge_list[,"from"]), gene_annot$entrezID, 
              nomatch = 0)])
    return_list$edge_list[toupper(return_list$edge_list[,"to"]) %in% 
                            gene_annot$entrezID,"to"] <- 
      tolower(gene_annot$gene_symbol[
        match(toupper(return_list$edge_list[,"to"]), gene_annot$entrezID, 
              nomatch = 0)])
    
    return_list$edge_list[,"edge"] <- paste(return_list$edge_list[,"from"], 
                                            return_list$edge_list[,"to"], 
                                            sep="_")
    # rename GE nodes to gene symbols
    return_list$node_list[return_list$node_list[,"label"] %in% 
                            gene_annot$entrezID,"label"] <- 
      gene_annot$gene_symbol[match(return_list$node_list[,"label"], 
                                   gene_annot$entrezID, nomatch = 0)]
    # rename CNV nodes to gene symbols
    return_list$node_list[toupper(return_list$node_list[,"label"]) %in% 
                            gene_annot$entrezID,"label"] <- 
      tolower(gene_annot$gene_symbol[
        match(toupper(return_list$node_list[,"label"]), gene_annot$entrezID, 
              nomatch = 0)])
  }
  
  net_weighted <-
    graph_from_edgelist(return_list$edge_list[,c("from","to")])
  igraph::V(net_weighted)$color <-
    return_list$node_list[match(as_ids(igraph::V(net_weighted)),
                                return_list$node_list[,"label"]),"color"]
  palette <- return_list$node_palette
  names(palette) <- seq_len(length(palette))
  palette <- palette[unique(igraph::V(net_weighted)$color)]
  igraph::V(net_weighted)$label <-
    return_list$node_list[match(as_ids(igraph::V(net_weighted)),
                                return_list$node_list[,"label"]),"label"]
  igraph::E(net_weighted)$edge <- return_list$edge_list[
    match(sub("|", "_", as_ids(igraph::E(net_weighted)), fixed = TRUE), 
          return_list$edge_list[,"edge"]),"edge_type"]
  igraph::E(net_weighted)$weight <- return_list$edge_list[
    match(sub("|", "_", as_ids(igraph::E(net_weighted)), fixed = TRUE),
          return_list$edge_list[,"edge"]),"weight"]
  
  igraph::V(net_weighted)$degree <- degree(net_weighted, mode = "in")
  igraph::V(net_weighted)$degree <- normalise(igraph::V(net_weighted)$degree,
                                              to = c(3, 11))
  
  return(list(edge_list = return_list$edge_list, node_palette = palette,
              node_list = return_list$node_list, net_weighted = net_weighted, 
              borders_GE = return_list$borders_GE,
              borders_CNV = return_list$borders_CNV, 
              borders_METH = return_list$borders_METH))
}
