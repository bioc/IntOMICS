#' BGe score  
#' @description
#' `BGe_score` Computes the BGe score of given network using precomputed sets
#' of possible parents.
#' @param adjacency_matrix adjacency matrix of given network.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' 
#' @examples
#' data("OMICS_mod_res", package="IntOMICS")
#' BGe_score(adjacency_matrix = OMICS_mod_res$B_prior_mat,
#'     layers_def = OMICS_mod_res$layers_def, omics = OMICS_mod_res$omics, 
#'     parent_set_combinations =
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node)
#'
#' @return Numeric vector of length 1: BGe score of given adjacency matrix
#' @export
BGe_score <- function(adjacency_matrix, omics, layers_def,
parent_set_combinations, BGe_score_all_configs_node)
{
    nodes_cand <- rownames(adjacency_matrix)[
        which(regexpr("EID",rownames(adjacency_matrix))>0)]
    score_nodes <- sum(unlist(lapply(nodes_cand,
        FUN=function(node) BGe_node(node = node, 
        adjacency_matrix = adjacency_matrix, 
        parent_set_combinations = parent_set_combinations,
        BGe_score_all_configs_node = BGe_score_all_configs_node))))
    return(score_nodes)
}
