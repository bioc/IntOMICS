#' Markov Blanket Resampling
#' @description
#' `MBR` This function performs the markov blanket resampling method according
#' to Su and Borsuk, 2016.
#' @param source_net_adjacency adajcency matrix of given network.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @importFrom bnlearn descendants amat empty.graph
#'
#' @examples
#' data("OMICS_mod_res", package="IntOMICS")
#' MBR(source_net_adjacency = OMICS_mod_res$B_prior_mat, 
#'     layers_def = OMICS_mod_res$layers_def, omics = OMICS_mod_res$omics, 
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node,
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations)
#'
#' @return List of 10 elements needed to define adjacency matrix 
#' with markov blanket resampling
MBR <- function(source_net_adjacency, layers_def, omics,
BGe_score_all_configs_node, parent_set_combinations) 
{
    selected_node <- sample(colnames(omics[[layers_def$omics[1]]]),1)
    dag_tmp <- source_net_adjacency
    current_parent_set <- names(which(dag_tmp[,selected_node]==1))
    if(length(current_parent_set)==0)
    {
        current_parent_set <- NA
    } # end if(length(current_parent_set)==0)
    children_selected_node <- names(which(dag_tmp[selected_node,]==1))
    dag_tmp[,selected_node] <- 0
    dag_tmp[-which(rownames(dag_tmp)==selected_node),
    children_selected_node] <- 0
    dag_tmp_bn <- bnlearn::empty.graph(rownames(dag_tmp),1)
    bnlearn::amat(dag_tmp_bn) <- dag_tmp
    descendants_selected_node <- bnlearn::descendants(x = dag_tmp_bn, 
        node = selected_node)
    selected_node_parents_scores <- parent_sets_sum_scores_X(
        selected_node = selected_node, parent_set = current_parent_set,
        parent_set_combinations = parent_set_combinations, 
        descendants = descendants_selected_node, 
        BGe_score_all_configs_node = BGe_score_all_configs_node)
    if(!is.na(selected_node_parents_scores$new_parent_set[1]))
    {
        bnlearn::amat(dag_tmp_bn)[selected_node_parents_scores$new_parent_set,
        selected_node] <- 1
    }
    if(length(children_selected_node)>0)
    {
        child_order <- sample(seq_len(length(children_selected_node)),
            length(children_selected_node))
        node_child_pars_scores <- parent_sets_sum_scores_childrenX(
            parent_set_combinations = parent_set_combinations, 
            selected_node = selected_node,
            children_selected_node = children_selected_node, 
            child_order = child_order,
            dag_tmp_bn = dag_tmp_bn,
            new_parent_set = TRUE,
            source_net_adjacency = source_net_adjacency, 
            BGe_score_all_configs_node = BGe_score_all_configs_node)
        dag_tmp_bn <- node_child_pars_scores$dag_tmp_bn
    } # if(length(children_selected_node)>0)
    candidate_net_adjacency <- bnlearn::amat(dag_tmp_bn)
    bnlearn::amat(dag_tmp_bn)[,selected_node] <- 0
    bnlearn::amat(dag_tmp_bn)[-which(rownames(dag_tmp)==selected_node),
    children_selected_node] <- 0
    selected_node_parents_scores_candidate <- parent_sets_sum_scores_X(
        parent_set_combinations = parent_set_combinations, 
        selected_node = selected_node, 
        descendants = descendants_selected_node, 
        parent_set = selected_node_parents_scores$new_parent_set, 
        BGe_score_all_configs_node = BGe_score_all_configs_node)
    bnlearn::amat(dag_tmp_bn)[current_parent_set, selected_node] <- 1
    if(length(children_selected_node)>0)
    {
        node_child_pars_scores_candidate <- parent_sets_sum_scores_childrenX(
            parent_set_combinations = parent_set_combinations,
            selected_node = selected_node, dag_tmp_bn = dag_tmp_bn, 
            children_selected_node = children_selected_node, 
            child_order = child_order, new_parent_set = FALSE,
            source_net_adjacency = source_net_adjacency, 
            BGe_score_all_configs_node = BGe_score_all_configs_node)
        r_source_candidate <- (
        selected_node_parents_scores$sum_score_unmarked + 
            sum(node_child_pars_scores$sum_score_unmarked)) - 
            (selected_node_parents_scores_candidate$sum_score_unmarked + 
            sum(node_child_pars_scores_candidate$sum_score_unmarked))
    } else {
        r_source_candidate <- 
        (selected_node_parents_scores$sum_score_unmarked) - 
        (selected_node_parents_scores_candidate$sum_score_unmarked)
    }

    return(list(adjacency = candidate_net_adjacency, nbhd.size = c(), 
        proposal.distr = c(), energy = c(), prior = c(), BGe = c(), 
        likelihood_part = c(), likelihood = c(), 
        acceptance = r_source_candidate, edge_move = c()))
}
