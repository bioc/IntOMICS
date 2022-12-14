#' MBR sum of children scores
#' @description
#' `parent_sets_sum_scores_children_x` This function determines the sum of BGe
#' scores of given node's children.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param selected_node character vector with given node name.
#' @param children_selected_node character vector with children 
#' of selected_node in given network structure.
#' @param child_order numeric vector random order of children_selected_node.
#' @param dag_tmp_bn object of class bn reflecting given network structure.
#' @param new_parent_set logical asking whether to define new parent set 
#' for selected_node.
#' @param source_net_adjacency adajcency matrix of given network.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @importFrom bnlearn descendants
#' @importFrom matrixStats logSumExp
#' @importFrom bnlearn amat
#' @return List of 3 elements
parent_sets_sum_scores_children_x <- function(parent_set_combinations, 
selected_node, children_selected_node, child_order, dag_tmp_bn, 
new_parent_set, source_net_adjacency, BGe_score_all_configs_node)
{
    sum_score_unmarked <- c()
    if(new_parent_set)
    {
        for(j in child_order)
        {
            descendants <- descendants(x = dag_tmp_bn, 
                node = children_selected_node[j])
            BGe_marked <-
                lapply(parent_set_combinations[[children_selected_node
                [j]]], FUN=function(list) apply(list, 2, FUN=function(column)
                length(intersect(column, descendants))>0 | 
                !any(column==selected_node)))
            BGe_marked[[1]][is.na(parent_set_combinations[[
                children_selected_node[j]]][[1]])] <- TRUE
            names(BGe_marked) <-
            paste(as.character(seq(1,length(BGe_marked))),"_",sep="")
            BGe_marked_compressed <- lapply(BGe_marked,FUN=function(list)
                which(list==FALSE))
            possible_parent_sets_ind <- unlist(BGe_marked_compressed, 
                use.names = TRUE)
            if(length(possible_parent_sets_ind)==0)
            {
                new_parent_set <- NA
                sum_score_unmarked[j] <- NA
            } else if(length(possible_parent_sets_ind)==1)
            {
                sum_score_unmarked[j] <- unlist(Map(function(pos, scores)
                    scores[!pos], BGe_marked,
                    BGe_score_all_configs_node[[children_selected_node[j]]]))
                ind <- as.numeric(unlist(lapply(strsplit(
                    names(possible_parent_sets_ind),"_"),FUN=function(list)
                    list[1])))
                new_parent_set <- parent_set_combinations[[
                children_selected_node[j]]][[ind]][,possible_parent_sets_ind]
            } else {
                score_unmarked <- unlist(Map(function(pos, scores) 
                scores[!pos], BGe_marked, BGe_score_all_configs_node[[
                children_selected_node[j]]]))
                new_parent_set_ind <- 
                sample(x = seq(1,length(possible_parent_sets_ind)), size = 1, 
                    prob = range_01(score_unmarked - sum(score_unmarked)))
                ind <- as.numeric(unlist(lapply(strsplit(names(
                    possible_parent_sets_ind[new_parent_set_ind]),"_"),
                    FUN=function(list) list[1])))
                new_parent_set <-
                parent_set_combinations[[children_selected_node[
                j]]][[ind]][,possible_parent_sets_ind[new_parent_set_ind]]
                sum_score_unmarked[j] <- logSumExp(score_unmarked)
            } # end if(length(possible_parent_sets_ind)==0)
            bnlearn::amat(dag_tmp_bn)[new_parent_set,
                children_selected_node[j]] <- 1
        } # end for j
        return(list(new_parent_set = new_parent_set, 
            sum_score_unmarked = sum_score_unmarked, dag_tmp_bn = dag_tmp_bn,
            BGe_marked = BGe_marked))
    } else {
    for(j in child_order)
    {
        descendants <- descendants(x = dag_tmp_bn, 
            node = children_selected_node[j])
        BGe_marked <- lapply(parent_set_combinations[[
            children_selected_node[j]]], FUN=function(list) apply(list, 2, 
            FUN=function(column) length(intersect(column, descendants))>0 | 
            !any(column==selected_node)))
        BGe_marked[[1]][is.na(parent_set_combinations[[children_selected_node[
        j]]][[1]])] <- TRUE
        names(BGe_marked) <-
        paste(as.character(seq(1,length(BGe_marked))),"_",sep="")
        BGe_marked_compressed <- lapply(BGe_marked,FUN=function(list)
            which(list==FALSE))
        possible_parent_sets_ind <- unlist(BGe_marked_compressed, 
            use.names = TRUE)
        if(length(possible_parent_sets_ind)==0)
        {
            sum_score_unmarked[j] <- NA
        } else if(length(possible_parent_sets_ind)==1)
        {
            sum_score_unmarked[j] <- unlist(Map(function(pos, scores) 
                scores[!pos], BGe_marked, BGe_score_all_configs_node[[
                children_selected_node[j]]]))
        } else {
            sum_score_unmarked[j] <-
                logSumExp(unlist(Map(function(pos, scores)
                scores[!pos], BGe_marked, BGe_score_all_configs_node[[
                children_selected_node[j]]])))
        } # end if(length(possible_parent_sets_ind)==0)
        bnlearn::amat(dag_tmp_bn)[names(which(source_net_adjacency[,
        children_selected_node[j]]==1)), children_selected_node[j]] <- 1
    } # end for j
    return(list(sum_score_unmarked = sum_score_unmarked, 
        dag_tmp_bn = dag_tmp_bn, BGe_marked = BGe_marked))
  } # end if(new_parent_set)
}
