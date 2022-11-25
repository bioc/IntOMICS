#' BGe score for specific node 
#' @description
#' `bge_node` Computes the BGe score of given node using precomputed sets 
#' of all possible parents.
#' @param node character vector with given node name.
#' @param adjacency_matrix adjacency matrix of given network.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @return Numeric vector of length 1: BGe score of given node
#' @keywords internal
#' @export
bge_node <- function(node, adjacency_matrix, parent_set_combinations,
BGe_score_all_configs_node)
{
    parents <- names(which(adjacency_matrix[,node]==1))
    if(length(parents)>0)
    {
        parents_ind <- 
        which(apply(parent_set_combinations[[node]][[length(parents)]], 
            2, FUN=function(column)
            length(intersect(column,parents))==length(parents)))
        score_node <- 
        BGe_score_all_configs_node[[node]][[length(parents)]][parents_ind]
    } else {
        score_node <- BGe_score_all_configs_node[[node]][[
            1]][is.na(parent_set_combinations[[node]][[1]])]
    } # end if(length(parents)>0)
    return(score_node)
}
