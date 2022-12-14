#' Sampling phase 
#' @description
#' `mcmc_simulation_sampling_phase` This function performs the final sampling
#' of network structures with estimated hyperparameters. 
#' It if part of sampling_phase function.
#' @param first numeric vector iteration to start.
#' @param last numeric vector iteration to stop.
#' @param sim_init list output from the source_net_def function or from two
#' independent simulations from the mcmc_simulation_sampling_phase function.
#' @param prob_mbr numeric vector probability of the MBR step.
#' @param B_prior_mat a biological prior matrix.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and features 
#' in columns.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param len numeric vector initial width of the sampling interval 
#' for hyperparameter beta.
#' @param thin numeric vector thinning frequency of the resulting 
#' MCMC simulation.
#' @param energy_all_configs_node list of nodes energy for all possible parent
#' set configurations.
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @importFrom bnlearn cpdag
#' @importFrom stats runif
#' @importFrom bnlearn amat
#' @importFrom bnlearn empty.graph
#' @return List of 1 element: sampling phase result before MCMC convergence
mcmc_simulation_sampling_phase <- function(first, last, sim_init, prob_mbr,
B_prior_mat, omics, parent_set_combinations, BGe_score_all_configs_node, 
layers_def, len, thin, energy_all_configs_node, annot)
{
    source.net <- sim_init$nets[[length(sim_init$nets)]]
    beta.source <- sim_init$betas[[length(sim_init$betas)]]
    start <- length(sim_init$nets)

    sim_init_forks <- list()
    for(j in seq_len(3))
    {
        sim_init_fork <- sim_init
        for (i in first:last)
        {
            sim_init$method_choice_saved[i] <- sample(x = c("MC3", "MBR"), 
                size = 1, prob = c(1-prob_mbr, prob_mbr))
            if(sim_init$method_choice_saved[i]=="MC3")
            {
                candidate.net <- mc3_constant_bge(source_net = source.net,
                    layers_def = layers_def, B_prior_mat = B_prior_mat, 
                    beta.source = beta.source, annot = annot, omics = omics, 
                    partition_func_UB_beta_source =
                    sim_init_fork$partition_func_UB_beta_source, 
                    parent_set_combinations = parent_set_combinations, 
                    BGe_score_all_configs_node = BGe_score_all_configs_node)
                candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
                source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
                candidate.net$likelihood <- candidate.net$likelihood_part +
                    source.net$proposal.distr
                source.net$likelihood <- source.net$likelihood_part +
                    candidate.net$proposal.distr
                sim_init_fork$acceptance_saved[i] <- candidate.net$likelihood -
                    source.net$likelihood
  
                u <- log(runif(1))
                if (u < sim_init_fork$acceptance_saved[i])
                {
                    source.net <- candidate.net
                } # end if (u < acceptance_saved[i])
                sim_init_fork$nets[[i]] <- source.net
            } else {
                candidate.net <- mbr(omics = omics, 
                    source_net_adjacency = source.net$adjacency, 
                    layers_def = layers_def, 
                    BGe_score_all_configs_node = BGe_score_all_configs_node, 
                    parent_set_combinations = parent_set_combinations)
                sim_init_fork$acceptance_saved[i] <- candidate.net$acceptance
                candidate.net$BGe <-bge_score(omics = omics, 
                    adjacency_matrix = candidate.net$adjacency, 
                    layers_def = layers_def, 
                    parent_set_combinations = parent_set_combinations, 
                    BGe_score_all_configs_node = BGe_score_all_configs_node)
                candidate.net$nbhd.size <- neighborhood_size(omics = omics,
                    net = candidate.net$adjacency, layers_def = layers_def, 
                    B_prior_mat = B_prior_mat)
                candidate.net$energy <- sum(epsilon(B_prior_mat = B_prior_mat,
                net = candidate.net$adjacency))
                candidate.net$prior <- 
                (-beta.source$value*candidate.net$energy) -
                sim_init_fork$partition_func_UB_beta_source
                candidate.net$likelihood_part <- 
                candidate.net$BGe + candidate.net$prior
        
                u <- log(runif(1))
                if (u < sim_init_fork$acceptance_saved[i])
                {
                    source.net <- candidate.net
                } # end if (u < acceptance_saved[i])
                sim_init_fork$nets[[i]] <- source.net
            } # end if else (sim_init_fork$method_choice_saved[i]=="MC3")
            if(i==last)
            {
                sim_init_fork$cpdags[[length(sim_init_fork$cpdags)+1]] <-
                empty.graph(rownames(sim_init_fork$nets[[i
                ]]$adjacency))
                bnlearn::amat(sim_init_fork$cpdags[[
                length(sim_init_fork$cpdags)]]) <-
                sim_init_fork$nets[[i]]$adjacency
                sim_init_fork$cpdags[[length(sim_init_fork$cpdags)]] <-
                    cpdag(sim_init_fork$cpdags[[
                    length(sim_init_fork$cpdags)]])
            } # end if(i==last)
        } # end for (i in first:last)
        sim_init_fork$nets[[i]]$BGe <- bge_score(omics = omics,
            adjacency_matrix =sim_init_fork$nets[[i]]$adjacency, 
            layers_def = layers_def, 
            parent_set_combinations = parent_set_combinations, 
            BGe_score_all_configs_node = BGe_score_all_configs_node)
        sim_init_fork$nets[[i]]$likelihood_part <-sim_init_fork$nets[[i]]$BGe +
            sim_init_fork$nets[[i]]$prior
        sim_init_forks[[j]] <- sim_init_fork
    } # end for j
    sim_init <-sim_init_forks[[which.max(mapply(sim_init_forks,
        FUN=function(list) list$nets[[length(list$nets)]]$likelihood_part))]]
    return(sim_init)
}
