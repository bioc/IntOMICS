#' transient phase 
#' @description
#' `transient_phase` This phase verify if the chain is moving towards 
#' the mode of target distribution.
#' @param first.adapt.phase_net list output of the first.adapt.phase 
#' function.
#' @param omics named list containing the gene expression 
#' (possibly copy number variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param B_prior_mat a biological prior matrix.
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param energy_all_configs_node list of nodes energy for all possible 
#' parent set configurations.
#' @param prob_mbr numeric vector probability of the MBR step.
#' @param BGe_score_all_configs_node list of nodes BGe score for all
#' possible parent set configurations.
#' @param parent_set_combinations list of all possible parent set 
#' configuration for all nodes available.
#' @param annot named list containing the associated methylation 
#' probes of given gene.
#' @importFrom stats runif lm rnorm
#' @importFrom utils tail
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", 
#'    "omics", "gene_annot", "OMICS_mod_res", "first.adapt.phase_net"), 
#'    package="IntOMICS")
#' if(interactive()){transient_phase(first.adapt.phase_net = 
#'    first.adapt.phase_net, omics = OMICS_mod_res$omics, 
#'    B_prior_mat = OMICS_mod_res$B_prior_mat, prob_mbr = 0.07, 
#'     layers_def = OMICS_mod_res$layers_def, annot = OMICS_mod_res$annot,
#'     energy_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations)}
#'
#' @return List of 1 element: first adaption phase and transient phase result
transient_phase <- function(first.adapt.phase_net, omics, B_prior_mat,
    layers_def, energy_all_configs_node, prob_mbr, BGe_score_all_configs_node, 
    parent_set_combinations, annot)
{
    beta_1st_adapt <- first.adapt.phase_net$betas[[length(
        first.adapt.phase_net$betas)]]$len
    source.net <- first.adapt.phase_net$nets[[length(
        first.adapt.phase_net$nets)]]
    beta.source <- first.adapt.phase_net$betas[[length(
        first.adapt.phase_net$betas)]]
    start <- length(first.adapt.phase_net$nets)
 
    for(i in (start+1):(start+1000))
    {
        first.adapt.phase_net$method_choice_saved[i] <- sample(x = c("MC3",
            "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
        if(first.adapt.phase_net$method_choice_saved[i]=="MC3")
        {
            candidate.net <- MC3(source_net = source.net, annot = annot,
                layers_def =  layers_def, B_prior_mat = B_prior_mat, 
                beta.source = beta.source, omics = omics, 
                partition_func_UB_beta_source = 
                first.adapt.phase_net$partition_func_UB_beta_source, 
                parent_set_combinations = parent_set_combinations, 
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
            } # end if (u < first.adapt.phase_net$acceptance_saved[i])
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
                B_prior_mat = B_prior_mat)
            candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency,
                B_prior_mat = B_prior_mat))
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
            partition_func_UB_beta_source <- sum(mapply(
                first.adapt.phase_net$energy_all_configs_node,
                FUN=function(x) matrixStats::logSumExp(-beta.source$value*x)))
        } # end if(method.choice=="MC3")
    
        beta.candidate <- list(value = stats::rnorm(1, 
            mean = beta.source$value, sd = beta_1st_adapt), 
            prior = c(), len = beta_1st_adapt)
        if(beta.candidate$value < 0.5)
        {
            beta.candidate$value <- 0.5
        } # end if(beta.candidate$value < 0.5)
        partition_func_UB_beta_candidate <- sum(mapply(
            first.adapt.phase_net$energy_all_configs_node, 
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
        first.adapt.phase_net$betas[[i]] <- beta.source
    } # end for(i in (start+1):(start+1000))
  
    beta_means <-
    colMeans(matrix(mapply(utils::tail(first.adapt.phase_net$betas,1000),
        FUN=function(list) list$value),nrow=200))
    reg_dat <- data.frame(beta_means = utils::tail(beta_means,5), 
        iter = seq_len(5))
    model <- stats::lm(beta_means ~ iter, data = reg_dat)
    p.val <- summary(model)$coefficients[1,4]
  
    while(p.val < 0.1)
    {
        for(i in (length(first.adapt.phase_net$nets)+1):(length(
        first.adapt.phase_net$nets)+200))
        {
            first.adapt.phase_net$method_choice_saved[i] <- 
            sample(x = c("MC3", "MBR"), size = 1, 
                prob = c(1-prob_mbr, prob_mbr))
            if(first.adapt.phase_net$method_choice_saved[i]=="MC3")
            {
                candidate.net <- MC3(source_net = source.net, annot = annot,
                    layers_def =  layers_def, B_prior_mat = B_prior_mat, 
                    beta.source = beta.source, omics = omics, 
                    partition_func_UB_beta_source = 
                    first.adapt.phase_net$partition_func_UB_beta_source,
                    parent_set_combinations = parent_set_combinations,
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
                } # end if (u < first.adapt.phase_net$acceptance_saved[i])
                first.adapt.phase_net$nets[[i]] <- source.net
            } else {
                candidate.net <- MBR(omics = omics, 
                    source_net_adjacency = source.net$adjacency, 
                    layers_def = layers_def, 
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
                    net = candidate.net$adjacency,
                    layers_def = layers_def, 
                    B_prior_mat = B_prior_mat)
            candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency,
                B_prior_mat = B_prior_mat))
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
            partition_func_UB_beta_source <-
            sum(mapply(first.adapt.phase_net$energy_all_configs_node, 
                FUN=function(x) matrixStats::logSumExp(-beta.source$value*x)))
        } # end if(method.choice=="MC3")
        beta.candidate <- list(value = stats::rnorm(1, sd = beta_1st_adapt,
            mean = beta.source$value), prior = c(), len = beta_1st_adapt)
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
        first.adapt.phase_net$betas[[i]] <- beta.source
      } # end for(i in (length(first.adapt.phase_net$nets)+1)...
    
      beta_means <-
      colMeans(matrix(mapply(utils::tail(first.adapt.phase_net$betas,1000),
           FUN=function(list) list$value),nrow=200))
      reg_dat <- data.frame(beta_means = beta_means, iter = seq_len(5))
      model <- stats::lm(beta_means ~ iter, data = reg_dat)
      p.val <- summary(model)$coefficients[1,4]
    } # end while(p.val < 0.1)
    return(first.adapt.phase_net)
}
