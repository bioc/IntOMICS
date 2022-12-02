#' Convert omics MultiAssayExperiment to list
#' @description
#' `omics_to_list` converts omics MultiAssayExperiment to list
#'
#' @param omics MultiAssayExperiment containing the gene expression 
#' (possibly copy number variation and methylation data).
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param gene_annot data.frame containing the entrez ID and corresponding 
#' gene symbol for conversion.
#' @importFrom SummarizedExperiment assay
#' @import HDF5Array
#' 
#' @examples
#' data(list=c("layers_def", "omics", "gene_annot"), package="IntOMICS")
#' omics_to_list(omics = omics,layers_def = layers_def,
#' gene_annot = gene_annot)
#'
#' @return List of omics modalities         
#' @export
omics_to_list <- function(omics, layers_def, gene_annot)
{
    omics_list <- list()
    for(i in seq(1,nrow(layers_def)))
    {
        omics_list[[i]] <- t(assay(omics[[layers_def$omics[i]]]))
        if(!is.numeric(omics_list[[i]]))
        {
          omics_list[[i]] <- apply(omics_list[[i]],2,as.numeric)
        } # end if(!is.numeric(omics_list[[i]]))
    } # end for(i in seq(1,nrow(layers_def)))
    names(omics_list) <- layers_def$omics
    colnames(omics_list$ge) <- gene_annot$entrezID[
      match(colnames(omics_list$ge),gene_annot$gene_symbol)]
    colnames(omics_list$cnv) <- tolower(gene_annot$entrezID[
      match(colnames(omics_list$cnv),gene_annot$gene_symbol)])
    omics <- omics_list[layers_def$omics[order(layers_def$layer, 
                                               decreasing = TRUE)]]
    return(omics)
}
