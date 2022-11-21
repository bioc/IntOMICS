#' Number of reverse edge candidates
#' @description
#' `fan_in_reverse` Determine the number of edges that can be reversed using
#' the fan-in restriction in the largest layer.
#' @param positions character vector indicating the interaction between two
#' nodes (the first string indicates the source node, the second string
#' indicates the target node).
#' @param net_layer_max adjacency matrix of the network containing only GE
#' nodes.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#'   
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", "gene_annot"),
#' package="IntOMICS")
#' colnames(omics$ge) <- gene_annot$entrezID[match(colnames(omics$ge),
#' gene_annot$gene_symbol)]
#' colnames(omics$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics$cnv),
#' gene_annot$gene_symbol)])
#' B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
#'     annot = annot, lm_METH = TRUE, r_squared_thres = 0.3,
#'     p_val_thres = 0.05, TFtargs = TFtarg_mat, TFBS_belief = 0.75, 
#'     nonGE_belief = 0.5, woPKGE_belief = 0.5)
#' adjacency_matrix <- B$B_prior_mat
#' adjacency_matrix[,] <- 0
#' adjacency_matrix[1,2] <- 1
#' layer_max <- colnames(B$omics[[layers_def$omics[1]]])
#' fan_in_reverse(positions = c(row=1,col=2), 
#' net_layer_max = adjacency_matrix[layer_max,layer_max], 
#' layers_def = layers_def)
#'
#' @return Numeric vector of length 1: reverse edge candidates
fan_in_reverse <- function(positions, net_layer_max, layers_def)
{
    net_layer_max[positions["col"],positions["row"]] <- 1
    possible_rev_edges <- sum(net_layer_max[,positions["col"]]) <=
        layers_def$fan_in_ge[1]
    return(possible_rev_edges)
}