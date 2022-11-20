#' OMICS_module
#' @description
#' `OMICS_module` data preprocessing + B_prior_mat definition + 
#'  partition function upper bound estimation + 
#'  all possible parent sets per node definition + 
#'  BGe score computation for all possible parent sets
#'
#' @param omics MatchedAssayExperiment containing the gene expression 
#' (possibly copy number variation and methylation data).
#' @param PK data.frame with known interactions.
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @param gene_annot data.frame containing the entrez ID and corresponding 
#' gene symbol for conversion.
#' @param lm_METH logical asking whether to use linear regression to filter
#' methylation data (default=TRUE).
#' @param r_squared_thres numeric vector to define the R^2 used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.3).
#' @param p_val_thres numeric vector to define the p-value used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.05).
#' @param TFBS_belief numeric vector to define the belief concerning the TF 
#' and its target interaction (default=0.75).
#' @param nonGE_belief numeric vector to define the belief concerning
#' interactions of features except GE-GE (default=0.5).
#' @param woPKGE_belief numeric vector to define the belief concerning GE-GE
#' interactions without prior knowledge (default=0.5).
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", 
#' "gene_annot"), package="IntOMICS")
#' OMICS_mod_res <- OMICS_module(omics = omics, PK = PK, 
#'     layers_def = layers_def, TFtargs = TFtarg_mat, annot = annot, 
#'     gene_annot = gene_annot, r_squared_thres = 0.3, lm_METH = TRUE)
#'
#' @return List of 6 elements needed to init MCMC simulation          
#' @export
OMICS_module <- function(omics, PK=NULL, layers_def, TFtargs=NULL, annot=NULL, 
lm_METH = TRUE, r_squared_thres = 0.3, p_val_thres = 0.05, TFBS_belief = 0.75, 
nonGE_belief = 0.5, woPKGE_belief = 0.5, gene_annot)
{
    omics_list <- list()
    for(i in seq(1,nrow(layers_def)))
    {
        omics_list[[i]] <- t(SummarizedExperiment::assay(omics[[layers_def$omics[i]]]))
        if(!is.numeric(omics_list[[i]]))
        {
            omics_list[[i]] <- apply(omics_list[[i]],2,as.numeric)
        } # end if(!is.numeric(omics_list[[i]]))
    } # end for(i in seq(1,nrow(layers_def)))
    names(omics_list) <- layers_def$omics
    colnames(omics_list$ge) <- gene_annot$entrezID[match(colnames(omics_list$ge),gene_annot$gene_symbol)]
    colnames(omics_list$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics_list$cnv),gene_annot$gene_symbol)])
    names(annot) <- gene_annot$entrezID[match(names(annot),gene_annot$gene_symbol)]

    if(TFBS_belief==woPKGE_belief)
    {
        message("GE-GE interactions without PK have the same belief as TFs-targets interactions. TFs-targets interactions will be considered as GE-GE interactions without prior knowledge.")
    }
  
    layers_def <- layers_def[order(layers_def$layer, decreasing = TRUE),]
    omics <- omics_list[layers_def$omics[order(layers_def$layer, 
        decreasing = TRUE)]]
  
    B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
        annot = annot, lm_METH = lm_METH, r_squared_thres = r_squared_thres, 
        p_val_thres = p_val_thres, TFtargs = TFtargs, 
        TFBS_belief = TFBS_belief, nonGE_belief = nonGE_belief, 
        woPKGE_belief = woPKGE_belief)
    pf_UB_res <- pf_UB_est(omics = B$omics, layers_def = layers_def, 
        B_prior_mat = B$B_prior_mat, annot = B$annot)
    return(list(pf_UB_BGe_pre = pf_UB_res, B_prior_mat = B$B_prior_mat, 
        annot = B$annot, omics = B$omics, layers_def = layers_def,
        omics_meth_original = B$omics_meth_original))
}
