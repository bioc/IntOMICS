#' Neighborhood size
#' @description
#' `neighborhood_size` This function is determines number of network structures
#' that can be reached from the current network structure.
#' @param net adajcency matrix of given network.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param B_prior_mat a biological prior matrix.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", "gene_annot"),
#' package="IntOMICS")
#' omics <- omics_to_list(omics = omics, gene_annot = gene_annot, 
#'                        layers_def = layers_def)
#' B <- b_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
#'      annot = annot, lm_METH = TRUE, r_squared_thres = 0.3,
#'      p_val_thres = 0.05, TFtargs = TFtarg_mat, TFBS_belief = 0.75, 
#'      nonGE_belief = 0.5, woPKGE_belief = 0.5)
#' adjacency_matrix <- B$B_prior_mat
#' adjacency_matrix[,] <- 0
#' neighborhood_size(net = adjacency_matrix, layers_def = layers_def,
#' B_prior_mat = B$B_prior_mat, omics = B$omics)
#'
#' @return Numeric of length 1: neighborhood size
#' @keywords internal
#' @export
neighborhood_size <- function(net, layers_def, B_prior_mat, omics)
{
    remove.edge.size <- sum(net)
    layer_max <- colnames(omics[[layers_def$omics[1]]])
    net_layer_max <- net[layer_max,layer_max]
    reverse_edge_pos <- which(net_layer_max==1, arr.ind = TRUE)
    reverse.edge.size <- sum(apply(reverse_edge_pos,1,
        FUN=function(x) fan_in_reverse(positions = x, 
        net_layer_max = net_layer_max, layers_def = layers_def)))
    layer_lower <- unlist(lapply(omics[layers_def$omics[-1]],colnames))
    net_layer_lower <- net[layer_lower,layer_max]
    B_prior_mat_layer_lower <- B_prior_mat[layer_lower,layer_max]
    add.edge.size <- sum(net_layer_lower[B_prior_mat_layer_lower > 0]==0)
    add.edge.size <- add.edge.size + sum(net_layer_max[,
        colSums(net_layer_max) <= (layers_def$fan_in_ge[1]-1)]==0)
    nbhd_size <- remove.edge.size + reverse.edge.size + add.edge.size
    return(nbhd_size)
}
