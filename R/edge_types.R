#' Resulting edge types definition
#' @description
#' `edge_types` Defines the resulting network structure.
#' @param B_prior_mat_weighted matrix one of the outputs 
#' of the bn_module function.
#' @param PK data.frame with known interactions.
#' @param gene_annot data.frame containing the entrez ID and corresponding 
#' gene symbol for conversion.
#' @param edge_list matrix indicating the interaction between nodes, 
#' the first column indicates the source node, the second column indicates 
#' the target node.
#' @param node_list character vector indicating the complete set of nodes 
#' in the resulting network structure.
#' @param OMICS_mod_res list output from the omics_module function.
#' @param edge_weights character vector includes either "MCMC_freq" to reflect
#' the edge weights frequency over the final set of network structures or 
#' "empB" to reflect the empirical biological knowledge estimated by IntOMICSr.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @return List of 6 elements needed to plot the final regulatory network edges
edge_types <- function(B_prior_mat_weighted, PK = NULL, gene_annot, edge_list, 
                       node_list, OMICS_mod_res, edge_weights, TFtargs = NULL)
{
    omics <- OMICS_mod_res$omics
    layers_def <- OMICS_mod_res$layers_def
    omics_meth_original <- OMICS_mod_res$omics_meth_original
  
    if(!is.null(PK))
    {
      PK <- paste(PK$src_entrez, PK$dest_entrez, sep="_")
    } # end if(!is.null(PK))
    
    if(edge_weights=="empB")
    {
      edge_list[,"edge_type"] <- "empirical"
      
      if(!is.null(TFtargs))
      {
        TF_pk <- as.matrix(TFtargs[intersect(edge_list[,"to"], 
                                             rownames(TFtargs)), 
                                   intersect(edge_list[,"from"],
                                             colnames(TFtargs))])
        colnames(TF_pk) <- intersect(edge_list[,"from"], colnames(TFtargs))
        if(ncol(TF_pk)>=1)
        {
          TF_pk <- paste(colnames(TF_pk)[which(TF_pk==1, 
                                               arr.ind = TRUE)[,2]], 
                         rownames(TF_pk)[which(TF_pk==1, 
                                               arr.ind = TRUE)[,1]], sep="_")
          edge_list[match(intersect(edge_list[,"edge"],TF_pk),
                          edge_list[,"edge"]),"edge_type"] <- "TF"
        } # end if(ncol(TF_pk)>=1)
      } # end if(!is.null(TFtargs))
      
      if(!is.null(PK))
      {
        edge_list[match(intersect(edge_list[,"edge"],PK),
                        edge_list[,"edge"]), "edge_type"] <- "PK"
      } # end if(!is.null(PK))
      
      edge_list[,"weight"] <- 
        round(apply(edge_list, 1, 
                    FUN=function(row) B_prior_mat_weighted[row["from"],
                                                           row["to"]]),2)
    } else {
      if(!is.null(PK))
      {
        edge_list[match(setdiff(edge_list[,"edge"],PK),
                        edge_list[,"edge"]),"edge_type"] <- "new"
        edge_list[match(intersect(edge_list[,"edge"],PK),
                        edge_list[,"edge"]), "edge_type"] <- "PK"
      } # end if(!is.null(PK))
    } # end if else (edge_weights=="empB")
    borders_output <- borders_def(node_list=node_list, layers_def=layers_def, 
                                  omics=omics, 
                                  omics_meth_original = omics_meth_original)
    
    return(list(edge_list = edge_list, node_list = borders_output$node_list, 
        borders_GE = borders_output$borders, 
        borders_CNV = borders_output$borders_cnv, 
        borders_METH = borders_output$borders_meth, 
        node_palette = borders_output$ge_cols))
}
