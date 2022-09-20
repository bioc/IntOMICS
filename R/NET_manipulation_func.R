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
#' @return Numeric vector of length 1: reverse edge candidates
#' @export
fan_in_reverse <- function(positions, net_layer_max, layers_def)
{
    net_layer_max[positions["col"],positions["row"]] <- 1
    possible_rev_edges <- sum(net_layer_max[,positions["col"]]) <=
        layers_def$fan_in_ge[1]
    return(possible_rev_edges)
}

#' Random initial network
#' @description
#' `init.net.mcmc` This function is used to sample random initial network. 
#' The edges are sampled only between GE nodes.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param B_prior_mat a biological prior matrix.
#'
#' @return List of 2 elements: random adjacency network and empty network
#' @export
init.net.mcmc <- function(omics, layers_def, B_prior_mat)
{
    empty.net <- matrix(0, nrow = sum(mapply(ncol,omics)), 
        ncol = sum(mapply(ncol,omics)),
        dimnames = list(unlist(mapply(colnames,omics)),
        unlist(mapply(colnames,omics))))
    init.net <- sample.chain(empty_net = empty.net, 
        omics_ge = omics[[layers_def$omics[1]]])
    rownames(init.net@dag) <- rownames(empty.net)
    colnames(init.net@dag) <- rownames(empty.net)
    init.net@dag <- init.net@dag[rownames(B_prior_mat),rownames(B_prior_mat)]
    if(any(init.net@dag==1 & B_prior_mat==0))
    {
        init.net@dag[init.net@dag==1 & B_prior_mat==0] <- 0
    }
    while(!is.acyclic(init.net@dag) |
    any(colSums(init.net@dag[colnames(omics[[layers_def$omics[1]]]),
    colnames(omics[[layers_def$omics[1]]])]) > layers_def$fan_in_ge[1]))
    {
        init.net <- sample.chain(empty_net = empty.net, 
            omics_ge = omics[[layers_def$omics[1]]])
        rownames(init.net@dag) <- rownames(empty.net)
        colnames(init.net@dag) <- rownames(empty.net)
        init.net@dag <- init.net@dag[rownames(B_prior_mat),
            rownames(B_prior_mat)]
        if(any(init.net@dag==1 & B_prior_mat==0))
        {
            init.net@dag[init.net@dag==1 & B_prior_mat==0] <- 0
        }
    }
    source.net <- list(adjacency = init.net@dag, nbhd.size = c(), 
        proposal.distr = c(), energy = c(), prior = c(), BGe = c(), 
        likelihood_part = c(), likelihood = c(), acceptance = c(), 
        edge_move = c())
    return(list(source.net = source.net, empty.net = empty.net))
}

#' Acyclic network identification.
#' @description
#' `is.acyclic` This function is from bnstruct R package. Check if the directed
#' graph is acyclic.
#' @param g adajcency matrix of given network/graph.
#'
#' @return boolean of length 1
#' @export
is.acyclic <- function(g)
{
    rem <- rep(FALSE,nrow(g))
    while( !all(rem) ) # still some edges to remove
    {
        leaves <- (rowSums(g) == 0)
        if( !any(leaves & !rem) )
            return(FALSE)
        g[,leaves] <- 0L
        rem <- rem | leaves
    }
  return(TRUE)
}

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
#' @return Numeric of length 1: neighborhood size
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

#' Random initial network edge generation
#' @description
#' `sample.chain` This function is used to sample random initial network. 
#' The edges are sampled only between GE nodes.
#' @param empty_net adjacency matrix of an empty network/graph 
#' (all values are 0).
#' @param omics_ge matrix with gene expression data (samples in rows and
#' features in columns).
#' @importFrom bnstruct BNDataset BN learn.params dag
#'
#' @return BN object with conditional probabilities
#' @export
sample.chain <- function(empty_net, omics_ge)
{
    dataset_BND <- bnstruct::BNDataset(data = empty_net, 
        discreteness = rep('d',ncol(empty_net)),
        variables = c(colnames(empty_net)), node.sizes = rep(2,ncol(empty_net)),
        starts.from=0)
    net <- bnstruct::BN(dataset_BND)
    net.dag <- bnstruct::dag(net)
    n <- ncol(omics_ge)
    chain <- sample(n,n)
    for(i in seq(from=2, to=n))
    {
        net.dag[chain[i-1],chain[i]] <- 1
    }
    bnstruct::dag(net) <- net.dag
    return(bnstruct::learn.params(net,dataset_BND))
}

#' Source network for MCMC simulation
#' @description
#' `source_net_def` This function is used to create the initial network 
#' with its features necessary for MCMC simulation.
#' @param init.net.mcmc.output list output of the init.net.mcmc function.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param parent_set_combinations list of all possible parent set 
#' configuration for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @param B_prior_mat a biological prior matrix.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param energy_all_configs_node list of nodes energy for all possible parent
#' set configurations.
#' @param len numeric vector initial width of the sampling interval 
#' for hyperparameter beta.
#' @importFrom stats runif
#' @importFrom matrixStats logSumExp
#'       
#' @return List of 10 elements needed to define the initial adjacency matrix            
#' @export
source_net_def <- function(init.net.mcmc.output, parent_set_combinations,
omics, BGe_score_all_configs_node, B_prior_mat, layers_def,
energy_all_configs_node, len)
{
    beta_init <- stats::runif(1, min = 0, max = 10)
    beta.source <- list(value = beta_init, prior = c())
    source.net <- init.net.mcmc.output$source.net
    source.net$BGe <- BGe_score(adjacency_matrix = source.net$adjacency, 
        omics = omics, layers_def = layers_def, 
        parent_set_combinations = parent_set_combinations,
        BGe_score_all_configs_node = BGe_score_all_configs_node)
    source.net$nbhd.size <- neighborhood_size(net = source.net$adjacency, 
        layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
    source.net$energy <- sum(epsilon(net = source.net$adjacency, 
        B_prior_mat = B_prior_mat))
    partition_func_UB_beta_source <- sum(mapply(energy_all_configs_node,
        FUN=function(x) matrixStats::logSumExp(-beta.source$value*x)))
    source.net$prior <- (-beta.source$value*source.net$energy) -
        partition_func_UB_beta_source
    source.net$likelihood_part <- source.net$BGe + source.net$prior
    beta.source$prior <- source.net$prior
    beta.source$len <- len
    acceptance_saved <- vector("numeric")
    acceptance_beta_saved <- vector("numeric")
    method_choice_saved <- vector("numeric")
    nets <- list()
    nets[[1]] <- source.net
    betas <- list()
    betas[[1]] <- beta.source
    return(list(source.net=source.net, beta.source = beta.source,
        partition_func_UB_beta_source=partition_func_UB_beta_source, 
        acceptance_saved = acceptance_saved, B_prior_mat = B_prior_mat,
        acceptance_beta_saved = acceptance_beta_saved, betas = betas,
        method_choice_saved = method_choice_saved, nets = nets,
        energy_all_configs_node = energy_all_configs_node))
}

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
#' @return List of 10 elements needed to define adjacency matrix 
#' with markov blanket resampling
#' @export
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

#' MBR sum of children scores
#' @description
#' `parent_sets_sum_scores_childrenX` This function determines the sum of BGe
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
#' @return List of 3 elements
#' @export
parent_sets_sum_scores_childrenX <- function(parent_set_combinations, 
selected_node, children_selected_node, child_order, dag_tmp_bn, 
new_parent_set, source_net_adjacency, BGe_score_all_configs_node)
{
    sum_score_unmarked <- c()
    if(new_parent_set)
    {
        for(j in child_order)
        {
            descendants <- bnlearn::descendants(x = dag_tmp_bn, 
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
                sum_score_unmarked[j] <- matrixStats::logSumExp(score_unmarked)
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
        descendants <- bnlearn::descendants(x = dag_tmp_bn, 
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
                matrixStats::logSumExp(unlist(Map(function(pos, scores)
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
#' @export
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

#' Markov Chain conventional single edge proposal move
#' @description
#' `MC3` This function samples a conventional single edge proposal move.
#' @param source_net list with adjacency matrix and other parameters needed 
#' for MCMC simulation.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param B_prior_mat a biological prior matrix.
#' @param beta.source named list with hyperparameter beta value and other
#' parameters needed for MCMC simulation.
#' @param partition_func_UB_beta_source numeric vector the upper bound 
#' of the partition function needed to define the prior distribution of 
#' network structure.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#'
#' @return List of 10 elements needed to define adjacency matrix 
#" with conventional single edge move
#' @export
MC3 <- function(source_net, omics, layers_def, B_prior_mat, beta.source, 
partition_func_UB_beta_source, parent_set_combinations,
BGe_score_all_configs_node, annot)
{
    ge_nodes <- rownames(source_net$adjacency)[regexpr("EID",
        rownames(source_net$adjacency))>0]
    vec <- seq_len(length(c(source_net$adjacency)))
    vec <- vec[c(B_prior_mat)>0]
    edge_proposal_res <- edge_proposal(net = source_net$adjacency, 
        candidates = vec, layers_def = layers_def, ge_nodes = ge_nodes, 
        omics = omics, B_prior_mat = B_prior_mat)
    while(edge_proposal_res$no_action | !is.acyclic(edge_proposal_res$net))
    {
        vec <- vec[vec!=edge_proposal_res$edge]
        edge_proposal_res <- edge_proposal(net = source_net$adjacency, 
            candidates = vec, layers_def = layers_def, ge_nodes = ge_nodes, 
            omics = omics, B_prior_mat = B_prior_mat)
    } # end while(edge_proposal_res$no_action...
    if(edge_proposal_res$edge_move=="add")
    {
        parents_source_target <- names(which(source_net$adjacency[,
            edge_proposal_res$col]==1))
        lp <- length(parents_source_target)
        if(lp>0)
        {
            ind_BGe_source <- apply(parent_set_combinations[[colnames(
                source_net$adjacency)[edge_proposal_res$col]]][[lp]], 
                2, FUN=function(col) all(parents_source_target %in% col))
            BGe_node_source <-
            BGe_score_all_configs_node[[colnames(source_net$adjacency
                )[edge_proposal_res$col]]][[lp]][ind_BGe_source]
            ind_BGe_candidate <-
            apply(parent_set_combinations[[colnames(source_net$adjacency
                )[edge_proposal_res$col]]][[lp+1]], 2, FUN=function(col)
                length(intersect(c(parents_source_target,
                colnames(source_net$adjacency
                )[edge_proposal_res$row]),col))>lp)
            BGe_node_candidate <-
            BGe_score_all_configs_node[[colnames(source_net$adjacency
                )[edge_proposal_res$col]]][[lp+1]][ind_BGe_candidate]
        } else {
            ind_BGe_source <- is.na(c(parent_set_combinations[[
                colnames(source_net$adjacency)[edge_proposal_res$col]]][[1]]))
            BGe_node_source <-
            BGe_score_all_configs_node[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[1]][ind_BGe_source]
            ind_BGe_candidate <-
                apply(parent_set_combinations[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[1]], 2, FUN=function(col)
                (colnames(source_net$adjacency)[edge_proposal_res$row] %in%
                col))
            BGe_node_candidate <-
            BGe_score_all_configs_node[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[1]][ind_BGe_candidate]
        }# end if else (lp>0)
        BGe <- source_net$BGe - BGe_node_source + BGe_node_candidate
    } else if (edge_proposal_res$edge_move=="delete")
    {
        parents_source_target <-
        names(which(source_net$adjacency[,edge_proposal_res$col]==1))
        lp <- length(parents_source_target)
        if(lp>1)
        {
            ind_BGe_source <- apply(parent_set_combinations[[
                colnames(source_net$adjacency)[edge_proposal_res$col]]][[lp]], 
                2, FUN=function(col) all(parents_source_target %in% col))
            BGe_node_source <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$col]]][[lp
                ]][ind_BGe_source]
            ind_BGe_candidate <-
            apply(parent_set_combinations[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[lp-1]], 2, FUN=function(col)
                length(intersect(setdiff(parents_source_target, 
                colnames(source_net$adjacency)
                [edge_proposal_res$row]),col))==(lp-1))
            BGe_node_candidate <-
                BGe_score_all_configs_node[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[lp-1]][ind_BGe_candidate]
        } else {
            ind_BGe_source <-
            which(c(parent_set_combinations[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[1]])==parents_source_target)
            BGe_node_source <-
                BGe_score_all_configs_node[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[1]][ind_BGe_source]
            ind_BGe_candidate <-
            is.na(c(parent_set_combinations[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[1]]))
            BGe_node_candidate <-
            BGe_score_all_configs_node[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[1]][ind_BGe_candidate]
        } # end if else (lp>0)
        BGe <- source_net$BGe - BGe_node_source + BGe_node_candidate
    } else {
        parents_source_target <-
        names(which(source_net$adjacency[,edge_proposal_res$col]==1))
        lp <- length(parents_source_target)
        if(lp>1)
        {
            ind_BGe_source <-
            apply(parent_set_combinations[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[lp]], 2, FUN=function(col)
                all(parents_source_target %in% col))
            BGe_node_source <-
            BGe_score_all_configs_node[[colnames(source_net$adjacency)
                [edge_proposal_res$col]]][[lp]][ind_BGe_source]
            ind_BGe_candidate <- apply(parent_set_combinations[[
                colnames(source_net$adjacency)[edge_proposal_res$col
                ]]][[lp-1]], 2, FUN=function(col)
                length(intersect(setdiff(parents_source_target,
                colnames(source_net$adjacency)[edge_proposal_res$row]),
                col))==(lp-1))
            BGe_node_candidate <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$col]]][[
                lp-1]][ind_BGe_candidate]
        } else {
            ind_BGe_source <- which(c(parent_set_combinations[[
                colnames(source_net$adjacency)[edge_proposal_res$col
                ]]][[1]])==parents_source_target)
            BGe_node_source <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$col]]][[1
                ]][ind_BGe_source]
            ind_BGe_candidate <- is.na(c(parent_set_combinations[[
                colnames(source_net$adjacency)[edge_proposal_res$col]]][[1]]))
            BGe_node_candidate <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$col]]][[1
                ]][ind_BGe_candidate]
        } # end if else (lp>0)
        BGe <- source_net$BGe - BGe_node_source + BGe_node_candidate
        parents_candidate_target <-
        names(which(source_net$adjacency[,edge_proposal_res$row]==1))
        lp <- length(parents_candidate_target)
        if(lp>0)
        {
            ind_BGe_source <- apply(parent_set_combinations[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[lp]], 2, 
                FUN=function(col) all(parents_candidate_target %in% col))
            BGe_node_source <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[lp
                ]][ind_BGe_source]
            ind_BGe_candidate <- apply(parent_set_combinations[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[lp+1]], 2, 
                FUN=function(col) length(intersect(c(parents_candidate_target,
                colnames(source_net$adjacency)
                [edge_proposal_res$col]),col))>lp)
            BGe_node_candidate <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[lp+1
                ]][ind_BGe_candidate]
        } else {
            ind_BGe_source <- is.na(c(parent_set_combinations[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[1]]))
            BGe_node_source <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[1
                ]][ind_BGe_source]
            ind_BGe_candidate <- apply(parent_set_combinations[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[1]], 2, 
                FUN=function(col)
                length(intersect(colnames(source_net$adjacency
                )[edge_proposal_res$col],col))==1)
            BGe_node_candidate <- BGe_score_all_configs_node[[colnames(
                source_net$adjacency)[edge_proposal_res$row]]][[1
                ]][ind_BGe_candidate]
        }# end if else (lp>0)
        BGe <- BGe - BGe_node_source + BGe_node_candidate
    } # end if(edge_proposal_res$edge_move=="add") ...
  
    nbhd.size <- neighborhood_size(net = edge_proposal_res$net, 
        layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
  
    energy <- sum(epsilon(net = edge_proposal_res$net, 
        B_prior_mat = B_prior_mat))
    prior <- (-beta.source$value*energy) - partition_func_UB_beta_source
    likelihood_part <- BGe + prior
    return(list(adjacency = edge_proposal_res$net, nbhd.size = nbhd.size, 
        proposal.distr = c(), energy = energy, prior = prior, BGe = BGe, 
        likelihood_part = likelihood_part, likelihood = c(), acceptance = c(), 
        edge_move = edge_proposal_res$edge_move))
}

#' Markov Chain conventional single edge proposal move with BGe score fixed
#' @description
#' `MC3_constantBGe` This function samples a conventional single edge proposal
#' move with ficed BGe score.
#' @param source_net list with adjacency matrix and other parameters needed 
#' for MCMC simulation.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param B_prior_mat a biological prior matrix.
#' @param beta.source named list with hyperparameter beta value and other
#' parameters needed for MCMC simulation.
#' @param partition_func_UB_beta_source numeric vector the upper bound 
#' of the partition function needed to define the prior distribution 
#' of network structure.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#'
#' @return List of 10 elements needed to define adjacency matrix 
#' with conventional single edge move
#' @export
MC3_constantBGe <- function(source_net, omics, layers_def, B_prior_mat, 
beta.source, partition_func_UB_beta_source, parent_set_combinations, 
BGe_score_all_configs_node, annot)
{
    ge_nodes <- rownames(source_net$adjacency)[regexpr("EID",
        rownames(source_net$adjacency))>0]
    vec <- seq_len(length(c(source_net$adjacency)))
    vec <- vec[c(B_prior_mat)>0]
    edge_proposal_res <- edge_proposal(omics = omics,
        net = source_net$adjacency, candidates = vec, layers_def = layers_def,
        ge_nodes = ge_nodes, B_prior_mat = B_prior_mat)
    while(edge_proposal_res$no_action | !is.acyclic(edge_proposal_res$net))
    {
        vec <- vec[vec!=edge_proposal_res$edge]
        edge_proposal_res <- edge_proposal(net = source_net$adjacency,
            candidates = vec, layers_def = layers_def, ge_nodes = ge_nodes, 
            omics = omics, B_prior_mat = B_prior_mat)
    } # end while(edge_proposal_res$no_action...
    nbhd.size <- neighborhood_size(net = edge_proposal_res$net, 
        layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
    energy <- sum(epsilon(net = edge_proposal_res$net, 
        B_prior_mat = B_prior_mat))
    prior <- (-beta.source$value*energy) - partition_func_UB_beta_source
    likelihood_part <- source_net$BGe + prior
    return(list(adjacency = edge_proposal_res$net, nbhd.size = nbhd.size, 
        proposal.distr = c(), energy = energy, prior = prior, likelihood = c(),
        BGe = source_net$BGe, likelihood_part = likelihood_part,
        acceptance = c(), edge_move = edge_proposal_res$edge_move))
}

#' Markov Chain conventional single edge proposal move
#' @description
#' `edge_proposal` This function samples a conventional single edge proposal
#' moves (identify those edges that are possible to change in given network
#' structure)
#' @param net adajcency matrix of given network.
#' @param candidates numeric vector with IDs of potential edge to be changed.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param ge_nodes character vector with GE node names
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and features 
#' in columns.
#' @param B_prior_mat a biological prior matrix.
#'                      
#' @return List of 6 elements needed to define candidates for conventional
#' single edge proposal move            
#' @export
edge_proposal <- function(net, candidates, layers_def, ge_nodes, omics,
B_prior_mat) {
    edge <- sample(candidates,1); div <- nrow(net); row <- edge %% div
    if(row==0)
    {
        row <- div; col <- edge %/% div
    } else {
        col <- edge %/% div + 1
    } # end if else (row==0)
    no_action <- FALSE
    if(net[edge]==1)
    {
        if(B_prior_mat[col,row]>0 &
            sum(net[ge_nodes,row])<layers_def$fan_in_ge[
            which.max(layers_def$layer)])
        {
            edge_move <- sample(c("delete","reverse"),1)
            if(edge_move=="delete")
            {
                net[edge] <- 0
            } else {
                net[col,row] <- 1; net[row,col] <- 0
            } # end if else (edge_move=="delete")
        } else {
            edge_move <- "delete"; net[edge] <- 0
        } # end if else ()
    } else {
        edge_move <- "add"
        candidate_layer <- names(which(mapply(omics,FUN=function(mat) 
            length(intersect(colnames(mat),rownames(net)[col]))==1)==TRUE))
        if(candidate_layer==layers_def$omics[which.max(layers_def$layer)])
        {
            if(sum(net[ge_nodes,col])<layers_def$fan_in_ge[
                which.max(layers_def$layer)])
            {
                net[edge] <- 1
            } else {
                no_action <- TRUE
            } # end if(source_net_adjacency[edge]==0 ...
        } else {
            if(regexpr("eid",rownames(net))>0)
            {
                net[edge] <- 1
            } else {
                if(sum(net[colnames(omics[[candidate_layer]]),
                    col])<layers_def$fan_in_ge[layers_def$omics==
                    candidate_layer])
                {
                    net[edge] <- 1
                } else {
                    no_action <- TRUE
                }
            } # end if else (regexpr("eid",rownames(net))>0)
        } # end if else (candidate_layer==layers_def$omics...
    } # end if else (net[edge]==1)
    return(list(net = net, edge = edge, no_action = no_action, row = row, 
        col = col, edge_move = edge_move))
}
