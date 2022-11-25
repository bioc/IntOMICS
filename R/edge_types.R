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
#' "empB" to reflect the empirical biological knowledge estimated by IntOMICS.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @return List of 6 elements needed to plot the final regulatory network edges
#' @keywords internal
#' @export 
edge_types <- function(B_prior_mat_weighted, PK = NULL, gene_annot, edge_list, 
                       node_list, OMICS_mod_res, edge_weights, TFtargs = NULL)
{
    omics <- OMICS_mod_res$omics
    layers_def <- OMICS_mod_res$layers_def
    omics_meth_original <- OMICS_mod_res$omics_meth_original
  
    if(any(regexpr("EID:",node_list)>0))
    {
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
                rownames(TFtargs)), intersect(edge_list[,"from"],
                colnames(TFtargs))])
                colnames(TF_pk) <- intersect(edge_list[,"from"], 
                colnames(TFtargs))
                if(ncol(TF_pk)>=1)
                {
                    TF_pk <- paste(colnames(TF_pk)[which(TF_pk==1, 
                    arr.ind = TRUE)[,2]], rownames(TF_pk)[which(TF_pk==1, 
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
            round(as.numeric(unlist(lapply(seq_along(edge_list[,2]), 1, 
                FUN=function(row)
                B_prior_mat_weighted[edge_list[row,"from"], 
                edge_list[row,"to"]]))), 2)
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
        omics=omics, omics_meth_original = omics_meth_original)
    } else {
        if(!is.null(PK))
        {
            PK_src_dest <- as.character(gene_annot$gene_symbol[match(
            PK$src_entrez,gene_annot$entrezID)])
            PK_src_dest[regexpr("eid",PK_src_dest)>0] <- tolower(
            as.character(gene_annot$gene_symbol[match(toupper(PK_src_dest[
            regexpr("eid",PK_src_dest)>0]), gene_annot$entrezID)]))
            PK_src_dest[is.na(PK_src_dest)] <- 
            PK$src_entrez[is.na(PK_src_dest)]
            PK <- paste(PK_src_dest, as.character(gene_annot$gene_symbol[
                match(PK$dest_entrez,gene_annot$entrezID)]), sep="_")
        } # end if(!is.null(PK))
        
        if(edge_weights=="empB")
        {
            edge_list[,"edge_type"] <- "empirical"
            targs_eid <- gene_annot$entrezID[match(edge_list[,"to"],
                gene_annot$gene_symbol)]
            TFs_eid <- gene_annot$entrezID[match(edge_list[,"from"],
                gene_annot$gene_symbol, nomatch = 0)]
            if(!is.null(TFtargs))
            {
                TF_pk <- as.matrix(TFtargs[intersect(targs_eid, rownames(TFtargs)),
                    intersect(TFs_eid, colnames(TFtargs))])
                colnames(TF_pk) <- intersect(TFs_eid, colnames(TFtargs))
                if(ncol(TF_pk)>=1)
                {
                    TF_pk <- paste(gene_annot$gene_symbol[match(colnames(TF_pk)
                        [which(TF_pk==1, arr.ind = TRUE)[,2]],
                        gene_annot$entrezID)], gene_annot$gene_symbol[match(
                        rownames(TF_pk)[which(TF_pk==1, arr.ind = TRUE)[,1]],
                        gene_annot$entrezID)], sep="_")
                    edge_list[match(intersect(edge_list[,"edge"],TF_pk),
                    edge_list[,"edge"]),"edge_type"] <- "TF"
                } # end if(ncol(TF_pk)>=1)
            } # end if(!is.null(TFtargs))
            
          
            if(!is.null(PK))
            {
                edge_list[match(intersect(edge_list[,"edge"],PK),
                edge_list[,"edge"]), "edge_type"] <- "PK"
            } # end if(!is.null(PK))
            
            rownames(B_prior_mat_weighted)[!is.na(match(rownames(
            B_prior_mat_weighted), gene_annot$entrezID))] <-
            gene_annot$gene_symbol[match(rownames(
            B_prior_mat_weighted), gene_annot$entrezID, nomatch = 0)]
            rownames(B_prior_mat_weighted)[!is.na(match(toupper(
            rownames(B_prior_mat_weighted)), gene_annot$entrezID))] <-
            tolower(gene_annot$gene_symbol[match(toupper(rownames(
            B_prior_mat_weighted)), gene_annot$entrezID, nomatch = 0)])
            colnames(B_prior_mat_weighted) <-
            rownames(B_prior_mat_weighted)
            edge_list[,"weight"] <-
            round(as.numeric(unlist(lapply(seq_along(edge_list[,2]),
                FUN=function(row) B_prior_mat_weighted[
                edge_list[row,"from"],edge_list[row,"to"]]))),2)
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
        omics=omics, omics_meth_original = omics_meth_original)
    } # end if else(any(regexpr("EID:",node_list)>0))
    return(list(edge_list = edge_list, node_list = borders_output$node_list, 
        borders_GE = borders_output$borders, borders_CNV = borders_output$borders_cnv, 
        borders_METH = borders_output$borders_meth, node_palette = borders_output$ge_cols))
}