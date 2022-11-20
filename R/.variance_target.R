#' Second adaption phase variance tuning
#' @description
#' `variance_target` This phase identifies the proposal distribution that 
#' has a similar covariance structure with the target distribution. 
#' This is part of second_adapt_phase.
#' @param transient.phase_net list output of the variance_target or
#' transient.phase function.
#' @param constant numeric vector used to multiply the beta_sd to determine 
#' the variance of the distribution of the hyperparameter beta.
#' @param fin numeric vector iteration to stop.
#' @param B_prior_mat a biological prior matrix.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param prob_mbr numeric vector probability of the MBR step.
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @importFrom utils tail
#' @importFrom stats rnorm sd
#'
#' @examples
#' data(list=c("OMICS_mod_res", "first.adapt.phase_net"),
#'     package="IntOMICS")
#' if(interactive()){transient.phase_net <- transient_phase(
#'     annot = OMICS_mod_res$annot, 
#'     first.adapt.phase_net = first.adapt.phase_net, 
#'     omics = OMICS_mod_res$omics, prob_mbr = 0.07, 
#'     B_prior_mat = OMICS_mod_res$B_prior_mat, 
#'     layers_def = OMICS_mod_res$layers_def, 
#'     energy_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations) 
#'     variance_target(annot = OMICS_mod_res$annot,
#'     constant = 1.586667, B_prior_mat = OMICS_mod_res$B_prior_mat,
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     layers_def = OMICS_mod_res$layers_def, omics = OMICS_mod_res$omics, 
#'     prob_mbr = 0.07, transient.phase_net = transient.phase_net, fin = 200)}
#'
#' @return Large List of 3 elements: second adaptive phase result 
#' with possible MCMC mixing; acceptance rate of hyperparameter beta; 
#' SD of hyperparameter beta
#' @export
variance_target <- function(transient.phase_net, constant, fin, B_prior_mat, 
omics, parent_set_combinations, BGe_score_all_configs_node, layers_def,
prob_mbr, annot)
{
    beta_sd <- stats::sd(mapply(utils::tail(transient.phase_net$betas,1000),
        FUN=function(list) list$value)[-1000])
    source.net <- transient.phase_net$nets[[length(transient.phase_net$nets)]]
    beta.source <-
    transient.phase_net$betas[[length(transient.phase_net$betas)]]
    start <- length(transient.phase_net$nets)

    for(i in (start+1):(start+fin))
    {
        transient.phase_net$method_choice_saved[i] <- 
        sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
        if(transient.phase_net$method_choice_saved[i]=="MC3")
        {
            candidate.net <- MC3(source_net = source.net, annot = annot,
                B_prior_mat = B_prior_mat, beta.source = beta.source, 
                partition_func_UB_beta_source =
                transient.phase_net$partition_func_UB_beta_source, 
                omics = omics, layers_def = layers_def,
                parent_set_combinations = parent_set_combinations,
                BGe_score_all_configs_node = BGe_score_all_configs_node)
            candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
            source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
            candidate.net$likelihood <- 
            candidate.net$likelihood_part + source.net$proposal.distr
            source.net$likelihood <- source.net$likelihood_part +
                candidate.net$proposal.distr
            transient.phase_net$acceptance_saved[i] <- 
            candidate.net$likelihood - source.net$likelihood
     
            u <- log(stats::runif(1))
            if (u < transient.phase_net$acceptance_saved[i])
            {
                source.net <- candidate.net
                beta.source$prior <- source.net$prior
            } # end if (u < transient.phase_net$acceptance_saved[i])
            transient.phase_net$nets[[i]] <- source.net
        } else {
            candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                layers_def = layers_def, omics = omics, 
                BGe_score_all_configs_node = BGe_score_all_configs_node, 
                parent_set_combinations = parent_set_combinations)
            transient.phase_net$acceptance_saved[i] <- candidate.net$acceptance
            candidate.net$BGe <- BGe_score(omics = omics,
                adjacency_matrix = candidate.net$adjacency, 
                layers_def = layers_def, 
                parent_set_combinations = parent_set_combinations, 
                BGe_score_all_configs_node = BGe_score_all_configs_node)
            candidate.net$nbhd.size <- neighborhood_size(omics = omics,
                net = candidate.net$adjacency, 
                B_prior_mat = B_prior_mat, layers_def = layers_def)
            candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency,
                B_prior_mat = B_prior_mat))
            candidate.net$prior <- (-beta.source$value*candidate.net$energy) -
                transient.phase_net$partition_func_UB_beta_source
            candidate.net$likelihood_part <- 
            candidate.net$BGe + candidate.net$prior
      
            u <- log(stats::runif(1))
            if (u < transient.phase_net$acceptance_saved[i])
            {
                source.net <- candidate.net
                beta.source$prior <- source.net$prior
            } # end if (u < acceptance_saved[i])
            transient.phase_net$nets[[i]] <- source.net
            partition_func_UB_beta_source <-
                sum(mapply(transient.phase_net$energy_all_configs_node,    
                FUN=function(x) 
                matrixStats::logSumExp(-beta.source$value*x)))
        } # end if(method.choice=="MC3")
    
        beta.candidate <- list(value = stats::rnorm(1, 
            mean = beta.source$value, sd = beta_sd*constant), prior = c(), 
            len = beta_sd*constant)
        if(beta.candidate$value < 0.5)
        {
            beta.candidate$value <- 0.5
        } # end if(beta.candidate$value < 0.5)
    
        partition_func_UB_beta_candidate <-
        sum(mapply(transient.phase_net$energy_all_configs_node,
            FUN=function(x) matrixStats::logSumExp(-beta.candidate$value*x)))
        beta.candidate$prior <- (-beta.candidate$value*source.net$energy) -
            partition_func_UB_beta_candidate
    
        transient.phase_net$acceptance_beta_saved[i] <- 
        beta.candidate$prior - beta.source$prior
        u_beta <- log(stats::runif(1))
        if (u_beta < transient.phase_net$acceptance_beta_saved[i])
        {
            beta.source <- beta.candidate
            transient.phase_net$partition_func_UB_beta_source <-
            partition_func_UB_beta_candidate
        } # end if (u_beta < transient.phase_net$acceptance_beta_saved[i])
        transient.phase_net$betas[[i]] <- beta.source
    } # end for(i in (start+1):(start+200))
  
    acceptance.trace_betas <-
    unlist(lapply(utils::tail(transient.phase_net$betas, 200),
        FUN=function(list) list$prior))
    acceptance.trace_betas <-
    c(1,acceptance.trace_betas[seq_len((length(acceptance.trace_betas)-1))] -
        acceptance.trace_betas[seq(from=2,to=length(acceptance.trace_betas))])
    acceptance.trace_betas[acceptance.trace_betas!=0] <- 1
    acceptance.rate_betas <- 
    sum(acceptance.trace_betas==1)/length(acceptance.trace_betas)
  
  return(list(variance.target_net = transient.phase_net, 
  acceptance.rate_betas = acceptance.rate_betas, beta_sd = beta_sd))
}
