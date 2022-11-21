#' Acceptance rate checking  
#' @description
#' `acceptance_check` This phase verify if the acceptance is in range 
#' of 0.28 and 0.6.
#' @param first.adapt.phase_net list output of the first.adapt.phase 
#' or source_net_def function.
#' @param round_check numeric vector after each round_check iterations 
#' for which we calculate the beta acceptance rate.
#' @param last_iter_check numeric vector number of the acceptance rate 
#' for the past last_iter_check iterations.
#' @param prob_mbr numeric vector probability of the MBR step.
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param parent_set_combinations list of all possible parent set 
#' configuration for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param annot named list containing the associated methylation 
#' probes of given gene.
#' @importFrom numbers mod
#' @importFrom stats runif rnorm
#' @importFrom matrixStats logSumExp
#' @importFrom utils tail
#' 
#' @examples
#' data(list=c("OMICS_mod_res", "first.adapt.phase_net"), package="IntOMICS")
#' if(interactive()){acceptance_check(round_check = 100, prob_mbr = 0.07,
#'     first.adapt.phase_net = first.adapt.phase_net, last_iter_check = 100,
#'     layers_def = OMICS_mod_res$layers_def, 
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations,
#'     omics = OMICS_mod_res$omics, annot = OMICS_mod_res$annot)}
#'
#' @return List of 1 element: first adaption phase result 
#' before given acceptance rate
acceptance_check <- function(first.adapt.phase_net, round_check, 
last_iter_check, prob_mbr, layers_def, parent_set_combinations,
BGe_score_all_configs_node, omics, annot)
{
    source.net <-
    first.adapt.phase_net$nets[[length(first.adapt.phase_net$nets)]]
    beta.source <-
    first.adapt.phase_net$betas[[length(first.adapt.phase_net$betas)]]
    i <- length(first.adapt.phase_net$nets)
    acceptance.rate_betas <- 0

    while(!(acceptance.rate_betas > 0.28 & acceptance.rate_betas < 0.60))
    {
        i <- i+1
        first.adapt.phase_net$method_choice_saved[i] <- 
        sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
        if(first.adapt.phase_net$method_choice_saved[i]=="MC3")
        {
            candidate.net <- MC3(source_net = source.net, omics = omics, 
                layers_def = layers_def, beta.source = beta.source, 
                B_prior_mat = first.adapt.phase_net$B_prior_mat, 
                partition_func_UB_beta_source = 
                first.adapt.phase_net$partition_func_UB_beta_source, 
                annot = annot, parent_set_combinations =
                parent_set_combinations, 
                BGe_score_all_configs_node = BGe_score_all_configs_node)
            candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
            source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
            candidate.net$likelihood <- candidate.net$likelihood_part +
                source.net$proposal.distr
            source.net$likelihood <- source.net$likelihood_part +
                candidate.net$proposal.distr
            first.adapt.phase_net$acceptance_saved[i] <-
            candidate.net$likelihood - source.net$likelihood
      
            u <- log(stats::runif(1))
            if (u < first.adapt.phase_net$acceptance_saved[i])
            {
                source.net <- candidate.net
                beta.source$prior <- source.net$prior
            }
            first.adapt.phase_net$nets[[i]] <- source.net
        } else {
            candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                layers_def = layers_def, omics = omics, 
                BGe_score_all_configs_node = BGe_score_all_configs_node, 
                parent_set_combinations = parent_set_combinations)
            first.adapt.phase_net$acceptance_saved[i] <-
            candidate.net$acceptance
            candidate.net$BGe <-BGe_score(omics = omics, 
                adjacency_matrix = candidate.net$adjacency, 
                layers_def = layers_def, 
                parent_set_combinations = parent_set_combinations, 
                BGe_score_all_configs_node = BGe_score_all_configs_node)
            candidate.net$nbhd.size <- neighborhood_size(omics = omics,
                net = candidate.net$adjacency, layers_def = layers_def, 
                B_prior_mat = first.adapt.phase_net$B_prior_mat)
            candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency, 
                B_prior_mat = first.adapt.phase_net$B_prior_mat))
            candidate.net$prior <- (-beta.source$value*candidate.net$energy) -
                first.adapt.phase_net$partition_func_UB_beta_source
            candidate.net$likelihood_part <- candidate.net$BGe +
                candidate.net$prior
      
            u <- log(stats::runif(1))
            if (u < first.adapt.phase_net$acceptance_saved[i])
            {
                source.net <- candidate.net
                beta.source$prior <- source.net$prior
            } # end if (u < acceptance_saved[i])
            first.adapt.phase_net$nets[[i]] <- source.net
        } # end if(method.choice=="MC3")
    
        beta.candidate <- list(value = stats::rnorm(n = 1, 
            mean = beta.source$value, 
            sd = beta.source$len), prior = c(), len = beta.source$len)
        if(beta.candidate$value < 0.5)
        {
            beta.candidate$value <- 0.5
        } # end if(beta.candidate$value < 0.5)
    
        partition_func_UB_beta_candidate <-
        sum(mapply(first.adapt.phase_net$energy_all_configs_node,
            FUN=function(x) matrixStats::logSumExp(-beta.candidate$value*x)))
        beta.candidate$prior <- (-beta.candidate$value*source.net$energy) -
            partition_func_UB_beta_candidate
    
        first.adapt.phase_net$acceptance_beta_saved[i] <- 
        beta.candidate$prior - beta.source$prior
        u_beta <- log(stats::runif(1))
    
        if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
        {
            beta.source <- beta.candidate
            first.adapt.phase_net$partition_func_UB_beta_source <-
            partition_func_UB_beta_candidate
        } # end if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
    
        if(numbers::mod(length(first.adapt.phase_net$nets), round_check)==0)
        {
            acceptance.trace_betas <-
            unlist(lapply(utils::tail(first.adapt.phase_net$betas,
                last_iter_check), FUN=function(list) list$prior))
            acceptance.trace_betas <-
            c(1,
            acceptance.trace_betas[seq_len((length(acceptance.trace_betas)-1))]
             - acceptance.trace_betas[seq(from=2,
             to=length(acceptance.trace_betas))])
            acceptance.trace_betas[acceptance.trace_betas!=0] <- 1
            acceptance.rate_betas <- sum(acceptance.trace_betas==1)/
                length(acceptance.trace_betas)
            # modify len if necessary
            if(acceptance.rate_betas > 0.44)
            {
                beta.source$len <- exp(log(beta.source$len) + 0.05)
            } else {
                beta.source$len <- exp(log(beta.source$len) - 0.05)
            } # end if else (acceptance.rate_betas > 0.44)
        } # end if(numbers::mod(i,round_check)==0)
        first.adapt.phase_net$betas[[i]] <- beta.source
    } # end while(!(acceptance.rate_betas > 0.28 & ...
    return(first.adapt.phase_net)
}
