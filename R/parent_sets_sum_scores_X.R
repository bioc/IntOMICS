#' MBR sum of scores
#' @description
#' `parent_sets_sum_scores_X` This function determines the sum of BGe scores 
#' of given node's parents.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param selected_node character vector with given node name.
#' @param descendants character vector with descendants of selected_node 
#' in given network structure.
#' @param parent_set character vector with parents of selected_node in given
#' network structure.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @importFrom matrixStats logSumExp
#' @return List of 3 elements
parent_sets_sum_scores_X <- function(parent_set_combinations, 
selected_node, descendants, parent_set, BGe_score_all_configs_node)
{
    BGe_marked <- lapply(parent_set_combinations[[selected_node]], 
        FUN=function(list) apply(list, 2, FUN=function(column)
        length(intersect(column, descendants))>0 | length(intersect(column,
        parent_set))==length(parent_set)))
    names(BGe_marked) <-
        paste(as.character(seq_len(length(BGe_marked))),"_",sep="")
    BGe_marked_compressed <- lapply(BGe_marked,FUN=function(list)
        which(list==FALSE))
    possible_parent_sets_ind <- unlist(BGe_marked_compressed, use.names = TRUE)
    if(length(possible_parent_sets_ind)==0)
    {
        new_parent_set <- NA
        sum_score_unmarked <- 0
    } else if(length(possible_parent_sets_ind)==1)
    {
        score_unmarked <- unlist(Map(function(pos, scores) scores[!pos], 
            BGe_marked, BGe_score_all_configs_node[[selected_node]]))
        ind <- as.numeric(unlist(lapply(strsplit(names(
        possible_parent_sets_ind),"_"), FUN=function(list) list[1])))
        new_parent_set <- parent_set_combinations[[selected_node]][[
            ind]][,possible_parent_sets_ind]
        sum_score_unmarked <- score_unmarked
    } else {
        score_unmarked <- unlist(Map(function(pos, scores) scores[!pos], 
            BGe_marked, BGe_score_all_configs_node[[selected_node]]))
        new_parent_set_ind <- 
            sample(x = seq_len(length(possible_parent_sets_ind)),
            size = 1, prob = range_01(score_unmarked - sum(score_unmarked)))
        ind <- as.numeric(unlist(lapply(strsplit(names(
            possible_parent_sets_ind[new_parent_set_ind]),"_"),
            FUN=function(list) list[1])))
        new_parent_set <- parent_set_combinations[[selected_node]][[
            ind]][,possible_parent_sets_ind[new_parent_set_ind]]
        sum_score_unmarked <- matrixStats::logSumExp(score_unmarked)
    } # end if(length(possible_parent_sets_ind)==0)
    return(list(new_parent_set = new_parent_set, 
    sum_score_unmarked = sum_score_unmarked, BGe_marked = BGe_marked))
}
