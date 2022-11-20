#' Node energy function
#' @description
#' `energy_function_node_specific`  For each node returns its energy over all
#' parent set configurations, the empty parent set is included.
#' @param all_parents_config matrix with all possible parent set configurations
#' (column indicates parents of given int_node).
#' @param B_prior_mat a biological prior matrix.
#' @param int_node character vector with given node name.
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", "gene_annot"),
#' package="IntOMICS")
#' colnames(omics$ge) <- gene_annot$entrezID[match(colnames(omics$ge),
#' gene_annot$gene_symbol)]
#' colnames(omics$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics$cnv),
#' gene_annot$gene_symbol)])
#' B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
#'    annot = annot, lm_METH = TRUE, r_squared_thres = 0.3, p_val_thres = 0.05,
#'    TFtargs = TFtarg_mat, TFBS_belief = 0.75, nonGE_belief = 0.5, 
#'    woPKGE_belief = 0.5)
#' all_parents_config <- matrix(c("EID:2535", "EID:2535", 
#' "EID:1857", "EID:2932"), byrow = TRUE, nrow=2)
#' energy_function_node_specific(all_parents_config = all_parents_config,
#' B_prior_mat = B$B_prior_mat, int_node = "EID:7482")
#'
#' @return Numeric vector of length 1            
#' @export
energy_function_node_specific <- function(all_parents_config, B_prior_mat,
int_node)
{
    epsilon <- apply(all_parents_config,2,FUN=function(z)
    if(is.na(z[1]))
    {
        sum(B_prior_mat[,int_node])
    } else {
        sum(1-B_prior_mat[z,int_node]) + 
            sum(B_prior_mat[-match(z,rownames(B_prior_mat)),int_node])
    })
    return(epsilon)
}
