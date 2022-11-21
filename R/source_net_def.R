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
#' @examples
#' data("OMICS_mod_res", package="IntOMICS")
#' init.net <- init.net.mcmc(omics = OMICS_mod_res$omics, 
#'     layers_def = OMICS_mod_res$layers_def, 
#'     B_prior_mat = OMICS_mod_res$B_prior_mat)
#' source_net_def(init.net.mcmc.output = init.net, 
#'     omics = OMICS_mod_res$omics, 
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations,
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     B_prior_mat = OMICS_mod_res$B_prior_mat, 
#'     layers_def = OMICS_mod_res$layers_def, len = 5,
#'     energy_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node)
#'              
#' @return List of 10 elements needed to define the initial adjacency matrix            
source_net_def <- function(init.net.mcmc.output, parent_set_combinations,
omics, BGe_score_all_configs_node, B_prior_mat, layers_def,
energy_all_configs_node, len)
{
    beta_init <- runif(1, min = 0, max = 10)
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
        FUN=function(x) logSumExp(-beta.source$value*x)))
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
