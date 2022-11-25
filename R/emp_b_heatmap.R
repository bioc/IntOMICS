#' Heatmap of empB - B
#' @description
#' `emp_b_heatmap` plot a heatmap with empB - B values (depicts the difference
#' between prior knowledge and the empirical knowledge)
#' @param mcmc_res list output from the bn_module function.
#' @param OMICS_mod_res list output from the omics_module function.
#' @param gene_annot data.frame containing the entrez ID and corresponding gene
#' symbol for conversion.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @importFrom gplots heatmap.2 bluered
#'
#' @examples
#' data(list=c("TFtarg_mat", "gene_annot", "OMICS_mod_res",
#' "BN_mod_res"), package="IntOMICS")
#' emp_b_heatmap(mcmc_res = BN_mod_res, OMICS_mod_res = OMICS_mod_res, 
#'     gene_annot = gene_annot, TFtargs = TFtarg_mat)
#'             
#' @return Figure heatmap
#' @export
emp_b_heatmap <- function(mcmc_res, OMICS_mod_res, gene_annot, TFtargs)
{
    mat <- mcmc_res$B_prior_mat_weighted - OMICS_mod_res$B_prior_mat
    mat <- mat[regexpr("EID:",rownames(mat))>0,
        regexpr("EID:",rownames(mat))>0]
    mat[which(OMICS_mod_res$B_prior_mat[regexpr("EID:",
        rownames(OMICS_mod_res$B_prior_mat))>0,regexpr("EID:",
        rownames(OMICS_mod_res$B_prior_mat))>0]==1)] <- NA
    TFtargs <- as.matrix(TFtargs[intersect(rownames(TFtargs),colnames(mat)),
        intersect(colnames(TFtargs),rownames(mat))])
    colnames(TFtargs) <- intersect(colnames(TFtargs),rownames(mat))
    if(ncol(TFtargs)>0)
    {
        inds <- which(TFtargs==1, arr.ind = TRUE)
        for(i in seq_len(nrow(inds)))
        {
            mat[colnames(TFtargs)[inds[i,2]],rownames(TFtargs)[inds[i,1]]] <- NA
        } # end for
    } # end if 
    rownames(mat) <- gene_annot$gene_symbol[match(rownames(mat),
        gene_annot$entrezID)]
    colnames(mat) <- rownames(mat)
    diag(mat) <- NA
    heatmap.2(mat, col=bluered, na.color="gray", Rowv = NA, 
        Colv = NA, trace = 'none', scale = 'none', main = "empB - B",
        dendrogram = "none")
}
