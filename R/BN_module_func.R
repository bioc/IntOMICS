#' #' BN module
#' @description
#' `BN_module` Performs automatically tuned MCMC sampling from posterior 
#' distribution together with conventional MCMC sampling using empirical
#' biological prior matrix to sample network structures from posterior
#' distribution.
#' @param burn_in numeric vector the minimal length of burn-in period 
#' of the MCMC simulation.
#' @param thin numeric vector thinning frequency of the resulting MCMC
#' simulation.
#' @param OMICS_mod_res list output from the OMICS_module function.
#' @param minseglen numeric vector minimal number of iterations 
#' with the c_rms value below the c_rms threshold.
#' @param len numeric vector initial width of the sampling interval 
#' for hyperparameter beta.
#' @param prob_mbr numeric vector probability of the MBR step.
#'
#' @examples
#' data("OMICS_mod_res", package="IntOMICS")
#' if(interactive()){BN_mod_res <- BN_module(burn_in = 100000, 
#'     thin = 500, OMICS_mod_res = OMICS_mod_res, 
#'     minseglen = 50000, len = 5, prob_mbr = 0.07)}
#' 
#' @return Large List of 3 elements: empirical biological matrix, 
#' sampling phase result and hyperparameter beta tuning trace
#' @export
BN_module <- function(burn_in, thin, OMICS_mod_res, minseglen, len = 5,
    prob_mbr = 0.07) {
    energy_all_configs_node <- 
    OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node
    BGe_score_all_configs_node <- 
    OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node
    parent_set_combinations <- 
    OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations
    annot <- OMICS_mod_res$annot
    ### 1st adaption phase ###
    first.adapt.phase_net <- first_adapt_phase(len = len,
        omics = OMICS_mod_res$omics, annot = annot,
        B_prior_mat = OMICS_mod_res$B_prior_mat, prob_mbr = prob_mbr, 
        energy_all_configs_node = energy_all_configs_node,
        layers_def = OMICS_mod_res$layers_def,
        BGe_score_all_configs_node = BGe_score_all_configs_node,
        parent_set_combinations = parent_set_combinations)
    ### transient phase ###
    transient.phase_net <- transient_phase(omics = OMICS_mod_res$omics,
        first.adapt.phase_net = first.adapt.phase_net, 
        B_prior_mat = OMICS_mod_res$B_prior_mat, prob_mbr = prob_mbr,
        layers_def = OMICS_mod_res$layers_def, annot = annot,
        energy_all_configs_node = energy_all_configs_node, 
        BGe_score_all_configs_node = BGe_score_all_configs_node,
        parent_set_combinations = parent_set_combinations)
    ### 2nd adaption phase ###
    second.adapt.phase_net <- second_adapt_phase(omics = OMICS_mod_res$omics,
        transient.phase_net = transient.phase_net, prob_mbr = prob_mbr,
        B_prior_mat = OMICS_mod_res$B_prior_mat, annot = annot,
        energy_all_configs_node = energy_all_configs_node, 
        layers_def = OMICS_mod_res$layers_def,
        BGe_score_all_configs_node = BGe_score_all_configs_node,
        parent_set_combinations = parent_set_combinations) 
    ### sampling phase ###
    sampling.phase_net <- sampling_phase(omics = OMICS_mod_res$omics,
        second.adapt.phase_net = second.adapt.phase_net, thin = thin,
        layers_def = OMICS_mod_res$layers_def, prob_mbr = prob_mbr,
        minseglen = minseglen, burn_in = burn_in, annot = annot)
    sampling.phase_net$mcmc_sim_part_res$seed1 <- 
    sampling.phase_net$mcmc_sim_part_res$seed1[c("betas","cpdags")]
    sampling.phase_net$mcmc_sim_part_res$seed2 <- 
    sampling.phase_net$mcmc_sim_part_res$seed2[c("betas","cpdags")]
    beta_tuning <- second.adapt.phase_net$betas
    return(list(sampling.phase_res = sampling.phase_net,
        B_prior_mat_weighted = second.adapt.phase_net$B_prior_mat_weighted, 
        beta_tuning = beta_tuning))
}

#' 1st adaption phase
#' @description
#' `first_adapt_phase` 1st adaption phase of the adaptive MCMC: 
#' the variance of the proposal distribution is changed to achieve 
#' the MC acceptance rate of 0.44.
#' @param omics named list containing the gene expression 
#' (possibly copy number variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param B_prior_mat a biological prior matrix.
#' @param energy_all_configs_node list of nodes energy for all possible 
#' parent set configurations.
#' @param len numeric vector initial width of the sampling interval 
#' for hyperparameter beta.
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param prob_mbr numeric vector probability of the MBR step.
#' @param BGe_score_all_configs_node list of nodes BGe score 
#' for all possible parent set configurations.
#' @param parent_set_combinations list of all possible parent set 
#' configuration for all nodes available.
#' @param annot named list containing the associated methylation 
#' probes of given gene.
#'
#' @examples
#' data("OMICS_module", package="IntOMICS")
#' if(interactive()){first_adapt_phase(omics = OMICS_mod_res$omics, 
#'    B_prior_mat = OMICS_mod_res$B_prior_mat, len = 5, 
#'    energy_all_configs_node = 
#'    OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
#'    layers_def = OMICS_mod_res$layers_def, prob_mbr = 0.07,
#'    BGe_score_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'    parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
#'    annot = OMICS_mod_res$annot)}
#'
#' @return List of 1 element: first adaption phase result
#' @export
first_adapt_phase <- function(omics, B_prior_mat, energy_all_configs_node, 
len, layers_def, prob_mbr, BGe_score_all_configs_node, 
parent_set_combinations, annot) {
    init.net1 <- init.net.mcmc(omics = omics, layers_def = layers_def,
    B_prior_mat = B_prior_mat)
    first.adapt.phase_net <- source_net_def(omics = omics, len = len,
        init.net.mcmc.output = init.net1, 
        parent_set_combinations = parent_set_combinations,
        BGe_score_all_configs_node = BGe_score_all_configs_node,
        B_prior_mat = B_prior_mat, layers_def = layers_def,
        energy_all_configs_node = energy_all_configs_node)
    first.adapt.phase_net <- acceptance_check(round_check = 100,
        first.adapt.phase_net = first.adapt.phase_net, 
        last_iter_check = 100, prob_mbr = prob_mbr, 
        layers_def = layers_def, omics = omics, annot = annot, 
        parent_set_combinations = parent_set_combinations, 
        BGe_score_all_configs_node = BGe_score_all_configs_node)
    first.adapt.phase_net <- acceptance_check(round_check = 100,
        first.adapt.phase_net = first.adapt.phase_net, 
        last_iter_check = 200, omics = omics, annot = annot, 
        prob_mbr = prob_mbr, layers_def = layers_def, 
        parent_set_combinations = parent_set_combinations, 
        BGe_score_all_configs_node = BGe_score_all_configs_node)
    first.adapt.phase_net <- acceptance_check(round_check = 200,
        first.adapt.phase_net = first.adapt.phase_net, 
        last_iter_check = 400, prob_mbr = prob_mbr, 
        layers_def = layers_def, omics = omics, annot = annot,
        parent_set_combinations = parent_set_combinations, 
        BGe_score_all_configs_node = BGe_score_all_configs_node)
    return(first.adapt.phase_net)
}

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
#' @export
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

#' Second adaption phase   
#' @description
#' `second_adapt_phase` This phase identifies the proposal distribution 
#' that has a similar covariance structure with the target distribution.
#' @param transient.phase_net list output of the transient.phase function.
#' @param omics named list containing the gene expression 
#' (possibly copy number variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param B_prior_mat a biological prior matrix.
#' @param energy_all_configs_node list of nodes energy for all possible 
#' parent set configurations.
#' @param prob_mbr numeric vector probability of the MBR step.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#'  parent set configurations.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param annot named list containing the associated methylation 
#' probes of given gene.
#' @param woPKGE_belief numeric vector to define the belief concerning 
#' GE-GE interactions without prior knowledge (default=0.5).
#' @importFrom utils tail
#' @importFrom stats lm 
#' 
#' @examples
#' data(list=c("OMICS_mod_res", "first.adapt.phase_net"), package="IntOMICS")
#' if(interactive()){transient.phase_net <- transient_phase(prob_mbr = 0.07,
#'     first.adapt.phase_net = first.adapt.phase_net, 
#'     omics = OMICS_mod_res$omics, B_prior_mat = OMICS_mod_res$B_prior_mat, 
#'     layers_def = OMICS_mod_res$layers_def, annot = OMICS_mod_res$annot,
#'     energy_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations) 
#'     second_adapt_phase(omics = OMICS_mod_res$omics,
#'     transient.phase_net = transient.phase_net, woPKGE_belief = 0.5,
#'     layers_def = OMICS_mod_res$layers_def, annot = OMICS_mod_res$annot,
#'     energy_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node, 
#'     prob_mbr = 0.07, B_prior_mat = OMICS_mod_res$B_prior_mat, 
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations)}
#'
#' @return List of 1 element: first adaption phase + transient phase + 
#' second adaption phase result
#' @export
second_adapt_phase <- function(transient.phase_net, omics, layers_def,
B_prior_mat, energy_all_configs_node, prob_mbr, BGe_score_all_configs_node,
parent_set_combinations, annot, woPKGE_belief = 0.5) 
{
    second.adapt.phase_net <- variance_target(fin = 200, 
        transient.phase_net = transient.phase_net, constant = 1.586667, 
        B_prior_mat = B_prior_mat, omics = omics, annot = annot,
        parent_set_combinations = parent_set_combinations, 
        BGe_score_all_configs_node = BGe_score_all_configs_node, 
        layers_def = layers_def, prob_mbr = prob_mbr)
    i <- 1
    while(second.adapt.phase_net$acceptance.rate_betas < 0.02)
    {
        i <- i+1
        constant <- 2.38/(1.5^i)
        second.adapt.phase_net <- variance_target(fin = 200,
            transient.phase_net = transient.phase_net, 
            constant = constant, B_prior_mat = B_prior_mat, omics = omics, 
            parent_set_combinations = parent_set_combinations, 
            BGe_score_all_configs_node = BGe_score_all_configs_node, 
            layers_def = layers_def, prob_mbr = prob_mbr, annot = annot)
    } # end while(second.adapt.phase_net$acceptance.rate_betas < 0.02)
    constant <- 2.38/(1.5^i)
    squared.jump_second.adapt.phase_net <-
        squared_jumping(omics = omics, B_prior_mat = B_prior_mat,
        second.adapt.phase_net = second.adapt.phase_net$variance.target_net, 
        constant = constant, prob_mbr = prob_mbr, annot = annot,
        beta_sd = second.adapt.phase_net$beta_sd,
        fin = (nrow(B_prior_mat)^2)*5, layers_def = layers_def,
        parent_set_combinations = parent_set_combinations, 
        BGe_score_all_configs_node = BGe_score_all_configs_node)

    betas_check <-
        mapply(utils::tail(squared.jump_second.adapt.phase_net$betas,
        1001), FUN=function(list) list$value)
    if(length(unique(betas_check))>1)
    {
        betas_check <- colMeans(matrix((betas_check[-1] -
            betas_check[-1001])^2,nrow=200))
        reg_dat <- data.frame(beta_means = betas_check, iter = seq_len(5))
        model <- stats::lm(beta_means ~ iter, data = reg_dat)
        squared.jump_second.adapt.phase_net$p.val <- 
        summary(model)$coefficients[2,4]
    
        while(squared.jump_second.adapt.phase_net$p.val < 0.1)
        {
            squared.jump_second.adapt.phase_net <- 
            squared_jumping(omics = omics, B_prior_mat = B_prior_mat,
            second.adapt.phase_net = squared.jump_second.adapt.phase_net, 
            constant = constant, beta_sd = second.adapt.phase_net$beta_sd,
            fin = 200, parent_set_combinations = parent_set_combinations,
            BGe_score_all_configs_node = BGe_score_all_configs_node, 
            layers_def = layers_def, prob_mbr = prob_mbr, annot = annot)
      
            betas_check <-
            mapply(utils::tail(squared.jump_second.adapt.phase_net$betas,1001),
                FUN=function(list) list$value)
            if(length(unique(betas_check))>1)
            {
                betas_check <- colMeans(matrix((betas_check[-1] -
                    betas_check[-1001])^2,nrow=200))
                reg_dat <- data.frame(beta_means = betas_check, 
                    iter = seq_len(5))
                model <- stats::lm(beta_means ~ iter, data = reg_dat)
                squared.jump_second.adapt.phase_net$p.val <- 
                summary(model)$coefficients[2,4]
            } else {
                squared.jump_second.adapt.phase_net$p.val <- 1
            }# end if else (length(unique(betas_check))>1)
        } # end while(squared.jump_second.adapt.phase_net$p.val < 0.1)
    } # end if(length(unique(betas_check))>1)
    squared.jump_second.adapt.phase_net$constant <- constant
    squared.jump_second.adapt.phase_net$beta_sd <-
    second.adapt.phase_net$beta_sd
    B_prior_mat_weighted <- c(B_prior_mat)
    conditions <- c(B_prior_mat)==woPKGE_belief &
    c(squared.jump_second.adapt.phase_net$iter_edges[,,1])>0
    B_prior_mat_weighted[conditions] <-
    c(squared.jump_second.adapt.phase_net$iter_edges[,,2])[conditions] /
        c(squared.jump_second.adapt.phase_net$iter_edges[,,1])[conditions]
    squared.jump_second.adapt.phase_net$B_prior_mat_weighted <-
    matrix(B_prior_mat_weighted, nrow=nrow(B_prior_mat), 
    dimnames = list(rownames(B_prior_mat), colnames(B_prior_mat)))
    squared.jump_second.adapt.phase_net$partition_func_UB <- 
    pf_UB_est(omics = omics, layers_def = layers_def, annot = annot,
        B_prior_mat = squared.jump_second.adapt.phase_net$B_prior_mat_weighted)
    return(squared.jump_second.adapt.phase_net)
}

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
#' @export
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

#' BGe score  
#' @description
#' `BGe_score` Computes the BGe score of given network using precomputed sets
#' of possible parents.
#' @param adjacency_matrix adjacency matrix of given network.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' 
#' @examples
#' data("OMICS_mod_res", package="IntOMICS")
#' BGe_score(adjacency_matrix = OMICS_mod_res$B_prior_mat,
#'     layers_def = OMICS_mod_res$layers_def, omics = OMICS_mod_res$omics, 
#'     parent_set_combinations =
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node)
#'
#' @return Numeric vector of length 1: BGe score of given adjacency matrix
#' @export
BGe_score <- function(adjacency_matrix, omics, layers_def,
parent_set_combinations, BGe_score_all_configs_node)
{
    nodes_cand <- rownames(adjacency_matrix)[
        which(regexpr("EID",rownames(adjacency_matrix))>0)]
    score_nodes <- sum(unlist(lapply(nodes_cand,
        FUN=function(node) BGe_node(node = node, 
        adjacency_matrix = adjacency_matrix, 
        parent_set_combinations = parent_set_combinations,
        BGe_score_all_configs_node = BGe_score_all_configs_node))))
    return(score_nodes)
}

#' BGe score for specific node 
#' @description
#' `BGe_node` Computes the BGe score of given node using precomputed sets 
#' of all possible parents.
#' @param node character vector with given node name.
#' @param adjacency_matrix adjacency matrix of given network.
#' @param parent_set_combinations list of all possible parent set configuration
#' for all nodes available.
#' @param BGe_score_all_configs_node list of nodes BGe score for all possible
#' parent set configurations.
#' 
#' @examples
#' data("OMICS_mod_res", package="IntOMICS")
#' BGe_node(node = "EID:2535", adjacency_matrix = OMICS_mod_res$B_prior_mat,
#'     parent_set_combinations =
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
#'     BGe_score_all_configs_node =
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node)
#'
#' @return Numeric vector of length 1: BGe score of given node
#' @export
BGe_node <- function(node, adjacency_matrix, parent_set_combinations,
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

#' Epsilon  
#' @description
#' `epsilon` This function returns the epsilon value for each variable/node 
#' of the network. 
#' The sum of the epsilons of all variables/nodes in the network gives us 
#' the energy of given network.
#' @param net adjacency matrix of given network.
#' @param B_prior_mat a biological prior matrix.
#'
#' @examples
#' data("OMICS_mod_res", package="IntOMICS")
#' adjacency_matrix <- OMICS_mod_res$B_prior_mat
#' adjacency_matrix[,] <- 0
#' epsilon(net = adjacency_matrix, B_prior_mat = OMICS_mod_res$B_prior_mat)
#'
#' @return Numeric vector of length 1: epsilon of given adjacency matrix
#' (needed to compute energy of given adjacency matrix)
#' @export
epsilon <- function(net, B_prior_mat)
{
    epsilon <- rep(NA,nrow(net))
    for(i in seq_len(nrow(net)))
    {
        iter_feature <- rownames(net)[i]
        parent <- names(net[,iter_feature])[net[,iter_feature]==1]
        noparent <- setdiff(rownames(net),parent)
        epsilon[i] <- sum(1-B_prior_mat[parent,iter_feature]) +
            sum(B_prior_mat[noparent,iter_feature])
    } # end for i
    return(epsilon)
}

#' Sampling phase 
#' @description
#' `mcmc.simulation_sampling.phase` This function performs the final sampling
#' of network structures with estimated hyperparameters. 
#' It if part of sampling_phase function.
#' @param first numeric vector iteration to start.
#' @param last numeric vector iteration to stop.
#' @param sim_init list output from the source_net_def function or from two
#' independent simulations from the mcmc.simulation_sampling.phase function.
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
#' importFrom bnlearn empty.graph amat cpdag
#' importFrom stats runif
#' @return List of 1 element: sampling phase result before MCMC convergence
#' @export
mcmc.simulation_sampling.phase <- function(first, last, sim_init, prob_mbr,
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
                candidate.net <- MC3_constantBGe(source_net = source.net,
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
  
                u <- log(stats::runif(1))
                if (u < sim_init_fork$acceptance_saved[i])
                {
                    source.net <- candidate.net
                } # end if (u < acceptance_saved[i])
                sim_init_fork$nets[[i]] <- source.net
            } else {
                candidate.net <- MBR(omics = omics, 
                    source_net_adjacency = source.net$adjacency, 
                    layers_def = layers_def, 
                    BGe_score_all_configs_node = BGe_score_all_configs_node, 
                    parent_set_combinations = parent_set_combinations)
                sim_init_fork$acceptance_saved[i] <- candidate.net$acceptance
                candidate.net$BGe <-BGe_score(omics = omics, 
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
        
                u <- log(stats::runif(1))
                if (u < sim_init_fork$acceptance_saved[i])
                {
                    source.net <- candidate.net
                } # end if (u < acceptance_saved[i])
                sim_init_fork$nets[[i]] <- source.net
            } # end if else (sim_init_fork$method_choice_saved[i]=="MC3")
            if(i==last)
            {
                sim_init_fork$cpdags[[length(sim_init_fork$cpdags)+1]] <-
                bnlearn::empty.graph(rownames(sim_init_fork$nets[[i
                ]]$adjacency))
                bnlearn::amat(sim_init_fork$cpdags[[
                length(sim_init_fork$cpdags)]]) <-
                sim_init_fork$nets[[i]]$adjacency
                sim_init_fork$cpdags[[length(sim_init_fork$cpdags)]] <-
                    bnlearn::cpdag(sim_init_fork$cpdags[[
                    length(sim_init_fork$cpdags)]])
            } # end if(i==last)
        } # end for (i in first:last)
        sim_init_fork$nets[[i]]$BGe <- BGe_score(omics = omics,
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

#' Sampling phase 
#' @description
#' `sampling_phase` Now we apply 2 MCMC simulations and check the RMS value. 
#' After the burn-in period, we discard the values from the first half 
#' of this phase.
#' @param second.adapt.phase_net list output of the second.adapt.phase
#' function.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param prob_mbr numeric vector probability of the MBR step.
#' @param thin numeric vector thinning frequency of the resulting MCMC
#' simulation.
#' @param minseglen numeric vector minimal number of iterations 
#' with the c_rms value below the c_rms threshold.
#' @param burn_in numeric vector the minimal length of burn-in period 
#' of the MCMC simulation.
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @importFrom bnlearn nodes
#' @importFrom utils tail
#' @importFrom bnlearn custom.strength
#' @importFrom stats quantile
#'
#' @examples
#' data(list=c("first.adapt.phase_net", "OMICS_mod_res"),
#'     package="IntOMICS")
#' if(interactive()){transient.phase_net <- transient_phase(prob_mbr = 0.07, 
#'     first.adapt.phase_net = first.adapt.phase_net, 
#'     omics = OMICS_mod_res$omics, B_prior_mat = OMICS_mod_res$B_prior_mat, 
#'     layers_def = OMICS_mod_res$layers_def, annot = OMICS_mod_res$annot,
#'     energy_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
#'     BGe_score_all_configs_node = 
#'     OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     parent_set_combinations = 
#'     OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations) 
#'     second.adapt.phase_net <- second_adapt_phase(prob_mbr = 0.07, 
#'     transient.phase_net = transient.phase_net, woPKGE_belief = 0.5, 
#'     omics = OMICS_mod_res$omics, B_prior_mat = OMICS_mod_res$B_prior_mat, 
#'     layers_def = OMICS_mod_res$layers_def, annot = OMICS_mod_res$annot,
#'     energy_all_configs_node =
#'     OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
#'     BGe_score_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#'     parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations) 
#'     sampling_phase(omics = OMICS_mod_res$omics, 
#'     second.adapt.phase_net = second.adapt.phase_net, 
#'     layers_def = OMICS_mod_res$layers_def, prob_mbr = 0.07, 
#'     thin = 500, minseglen = 50000, burn_in = 100000, 
#'     annot = OMICS_mod_res$annot)}
#'
#' @return List of 2 elements: sampling phase result; RMS used to evaluate 
#' MCMC convergence
#' @export
sampling_phase <- function(second.adapt.phase_net, omics, layers_def, prob_mbr,
thin, minseglen, burn_in, annot) 
{
    init.net_sampling <- init.net.mcmc(omics = omics, layers_def = layers_def,
        B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted)
    init.net_sampling <- source_net_def(omics = omics, 
        init.net.mcmc.output = init.net_sampling, 
        parent_set_combinations =
        second.adapt.phase_net$partition_func_UB$parents_set_combinations,
        BGe_score_all_configs_node =
        second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node,
        B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted,
        layers_def = layers_def,
        energy_all_configs_node =
        second.adapt.phase_net$partition_func_UB$energy_all_configs_node,
        len = utils::tail(second.adapt.phase_net$betas,1)[[1]][["len"]])
    rms <- c()
    seeds_res <- list(seed1=list(),seed2=list())
    seeds_res$seed1$nets <- utils::tail(second.adapt.phase_net$nets,1)
    seeds_res$seed1$betas <- utils::tail(second.adapt.phase_net$betas,1)
    seeds_res$seed1$partition_func_UB_beta_source <-
        second.adapt.phase_net$partition_func_UB_beta_source
    seeds_res$seed1$nets[[1]]$adjacency <-
        init.net_sampling$source.net$adjacency
    seeds_res$seed1$nets[[1]]$nbhd.size <- neighborhood_size(omics = omics, 
        net = seeds_res$seed1$nets[[1]]$adjacency, layers_def = layers_def, 
        B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted)
    seeds_res$seed1$nets[[1]]$energy <- 
        sum(epsilon(net = seeds_res$seed1$nets[[1]]$adjacency, 
        B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted))
    seeds_res$seed1$nets[[1]]$prior <- 
    (-seeds_res$seed1$betas[[1]]$value*seeds_res$seed1$nets[[1]]$energy) -
        seeds_res$seed1$partition_func_UB_beta_source
    seeds_res$seed1$nets[[1]]$BGe <- BGe_score(omics = omics, 
        adjacency_matrix = seeds_res$seed1$nets[[1]]$adjacency, 
        layers_def = layers_def, 
        parent_set_combinations =
        second.adapt.phase_net$partition_func_UB$parents_set_combinations,
        BGe_score_all_configs_node =
        second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node)
    seeds_res$seed1$nets[[1]]$likelihood_part <- 
    seeds_res$seed1$nets[[1]]$BGe + seeds_res$seed1$nets[[1]]$prior
    seeds_res$seed1$betas[[1]]$prior <- seeds_res$seed1$nets[[1]]$prior
    seeds_res$seed1$nets[[1]]$proposal.distr <- c()
    seeds_res$seed1$acceptance_saved <- vector("numeric")
    seeds_res$seed1$method_choice_saved <- vector("numeric")
    seeds_res$seed1$layers <- second.adapt.phase_net$layers
    seeds_res$seed1$cpdags <- list()
    seeds_res$seed2 <- seeds_res$seed1
    seeds_res$seed2$nets[[1]]$adjacency[
    seeds_res$seed2$nets[[1]]$adjacency==1] <- 0
    seeds_res$seed2$nets[[1]]$nbhd.size <- neighborhood_size(omics = omics,
        net = seeds_res$seed2$nets[[1]]$adjacency, layers_def = layers_def, 
        B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted)
    seeds_res$seed2$nets[[1]]$energy <- 
    sum(epsilon(net = seeds_res$seed2$nets[[1]]$adjacency, 
        B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted))
    seeds_res$seed2$nets[[1]]$prior <- 
    (-seeds_res$seed2$betas[[1]]$value*seeds_res$seed2$nets[[1]]$energy) -
        seeds_res$seed2$partition_func_UB_beta_source
    seeds_res$seed2$nets[[1]]$BGe <- BGe_score(omics = omics, 
        adjacency_matrix = seeds_res$seed2$nets[[1]]$adjacency, 
        layers_def = layers_def, 
        parent_set_combinations =
        second.adapt.phase_net$partition_func_UB$parents_set_combinations,
        BGe_score_all_configs_node =
        second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node)
    seeds_res$seed2$nets[[1]]$likelihood_part <- 
    seeds_res$seed2$nets[[1]]$BGe + seeds_res$seed2$nets[[1]]$prior
    seeds_res$seed2$betas[[1]]$prior <- seeds_res$seed2$nets[[1]]$prior
    mcmc_sim_part_res <- lapply(seeds_res, FUN=function(list_l)
    mcmc.simulation_sampling.phase(first = 1, last = thin, sim_init = list_l,
    prob_mbr = prob_mbr, omics = omics, annot = annot,
    B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted,
    parent_set_combinations =
    second.adapt.phase_net$partition_func_UB$parents_set_combinations,
    BGe_score_all_configs_node =
    second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node,
    layers_def = layers_def, len = seeds_res$seed1$betas[[1]]$len, thin = thin,
    energy_all_configs_node =
    second.adapt.phase_net$partition_func_UB$energy_all_configs_node))
    cpdags1 <- mcmc_sim_part_res$seed1$cpdags
    cpdags2 <- mcmc_sim_part_res$seed2$cpdags
    cpdag_weights1 <- bnlearn::custom.strength(cpdags1, 
        nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
    cpdag_weights2 <- bnlearn::custom.strength(cpdags2, 
        nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
    cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
    cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]
    total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
    N <- nrow(total)
    dist_i <- abs(total$strength.x - total$strength.y)^2 / 2
    rms <- c(rms,sqrt(1/N*sum(dist_i)))
 
    while(length(mcmc_sim_part_res$seed1$nets)<(2*burn_in))
    {
        mcmc_sim_part_res <- lapply(mcmc_sim_part_res, FUN=function(list_l) 
            mcmc.simulation_sampling.phase(first = length(list_l$nets)+1,
            last = length(list_l$nets)+thin, sim_init = list_l, annot = annot, 
            prob_mbr = prob_mbr, omics = omics, thin = thin, 
            B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted,
            parent_set_combinations =
            second.adapt.phase_net$partition_func_UB$parents_set_combinations,
            BGe_score_all_configs_node =
        second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node, 
            layers_def = layers_def, len = seeds_res$seed1$betas[[1]]$len, 
            energy_all_configs_node =
            second.adapt.phase_net$partition_func_UB$energy_all_configs_node))
        cpdags1 <- unique(mcmc_sim_part_res$seed1$cpdags)
        cpdags2 <- unique(mcmc_sim_part_res$seed2$cpdags)
        cpdag_weights1 <- custom.strength(cpdags1, 
            nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
        cpdag_weights2 <- custom.strength(cpdags2, 
            nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
        cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
        cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]
        total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
        N <- nrow(total)
        dist_i <- abs(total$strength.x - total$strength.y)^2 / 2
        rms <- c(rms,sqrt(1/N*sum(dist_i)))
    } # end while(length(mcmc_sim_part_res$seed1$nets)<(2*burn_in))
    rms_strength <- abs(diff(rms))
    strength_threshold <- stats::quantile(rms_strength, 0.75, na.rm = TRUE)
  
    while(any(utils::tail(rms_strength,minseglen/thin)>strength_threshold))
    {
        mcmc_sim_part_res <- lapply(mcmc_sim_part_res, FUN=function(list_l)
            mcmc.simulation_sampling.phase(first = length(list_l$nets)+1,
            last = length(list_l$nets)+thin, sim_init = list_l, omics = omics,
            prob_mbr = prob_mbr, B_prior_mat =
            second.adapt.phase_net$B_prior_mat_weighted, 
            parent_set_combinations =
            second.adapt.phase_net$partition_func_UB$parents_set_combinations,
            BGe_score_all_configs_node =
        second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node,
            layers_def = layers_def, len = seeds_res$seed1$betas[[1]]$len, 
            thin = thin, annot = annot,
            energy_all_configs_node =
            second.adapt.phase_net$partition_func_UB$energy_all_configs_node))
        cpdags1 <- unique(mcmc_sim_part_res$seed1$cpdags)
        cpdags2 <- unique(mcmc_sim_part_res$seed2$cpdags)
        cpdag_weights1 <- custom.strength(cpdags1, 
            nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
        cpdag_weights2 <- custom.strength(cpdags2, 
            nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
        cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
        cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]
        total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
        N <- nrow(total)
        dist_i <- abs(total$strength.x - total$strength.y)^2 / 2
        rms <- c(rms,sqrt(1/N*sum(dist_i)))
        rms_strength <- abs(diff(rms))
    } # end while(any(utils::tail(rms_strength,...
    return(list(mcmc_sim_part_res = mcmc_sim_part_res, rms = rms))
}

#' Squared jumping of adaptive MCMC algorithm
#' @description
#' `squared_jumping` Squared jumping of adaptive MCMC algorithm is used to tune
#' the variance of the beta parameter.
#' @param second.adapt.phase_net list output of the variance_target 
#' or squared_jumping function.
#' @param constant numeric vector used to multiply the beta_sd to determine 
#' the variance of the distribution of the hyperparameter beta.
#' @param beta_sd numeric vector used to determine the variance 
#' of the distribution of the hyperparameter beta.
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
#' @param annot named list containing the associated methylation probes of
#' given gene.
#' @importFrom stats rnorm
#' @importFrom matrixStats logSumExp
#'
#' @examples
#' data(list=c("OMICS_mod_res", "first.adapt.phase_net", "annot"), 
#' package="IntOMICS")
#' if(interactive()){transient.phase_net <- transient_phase(prob_mbr = 0.07,
#' first.adapt.phase_net = first.adapt.phase_net, 
#' omics = OMICS_mod_res$omics, B_prior_mat = OMICS_mod_res$B_prior_mat, 
#' layers_def = OMICS_mod_res$layers_def, annot = OMICS_mod_res$annot,
#' energy_all_configs_node = 
#' OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
#' BGe_score_all_configs_node = 
#' OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#' parent_set_combinations = 
#' OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations) 
#' second_adapt_phase_net <- second_adapt_phase(prob_mbr = 0.07, 
#' transient.phase_net = transient.phase_net, woPKGE_belief = 0.5,
#' layers_def = OMICS_mod_res$layers_def, annot = OMICS_mod_res$annot,
#' energy_all_configs_node = 
#' OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node, 
#' B_prior_mat = OMICS_mod_res$B_prior_mat, 
#' BGe_score_all_configs_node = 
#' OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
#' parent_set_combinations = 
#' OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations,
#' omics = OMICS_mod_res$omics)
#' squared_jumping(omics = OMICS_mod_res$omics,
#' second.adapt.phase_net = second.adapt.phase_net$variance.target_net, 
#' constant = 2.38/(1.5), fin = (nrow(OMICS_mod_res$B_prior_mat)^2)*5, 
#' beta_sd = second.adapt.phase_net$beta_sd, prob_mbr = 0.07,
#' B_prior_mat = OMICS_mod_res$B_prior_mat, 
#' layers_def = OMICS_mod_res$layers_def,
#' parent_set_combinations=OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations,
#' BGe_score_all_configs_node=OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node,
#' annot = annot)}
#'
#' @return List of 1 element: second adaptive phase result with stopped 
#' MCMC mixing
#' @export
squared_jumping <- function(second.adapt.phase_net, constant, fin, beta_sd,
B_prior_mat, omics, parent_set_combinations, BGe_score_all_configs_node,
layers_def, prob_mbr, annot)
{
    source.net <-
    second.adapt.phase_net$nets[[length(second.adapt.phase_net$nets)]]
    beta.source <-
    second.adapt.phase_net$betas[[length(second.adapt.phase_net$betas)]]
    start <- length(second.adapt.phase_net$nets)
    second.adapt.phase_net$iter_edges <- array(0, dim=c(dim(B_prior_mat),2), 
        dimnames = list(rownames(B_prior_mat), colnames(B_prior_mat),
        c("frequency", "acceptance")))

    for(i in (length(second.adapt.phase_net$nets)+1):
    (length(second.adapt.phase_net$nets)+fin))
    {
        second.adapt.phase_net$method_choice_saved[i] <- 
        sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
        if(second.adapt.phase_net$method_choice_saved[i]=="MC3")
        {
            candidate.net <- MC3(source_net = source.net, omics = omics, 
            layers_def =  layers_def, B_prior_mat = B_prior_mat, 
            beta.source = beta.source, annot = annot,
            partition_func_UB_beta_source =
            second.adapt.phase_net$partition_func_UB_beta_source,
            parent_set_combinations = parent_set_combinations,
            BGe_score_all_configs_node = BGe_score_all_configs_node)
            candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
            source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
            candidate.net$likelihood <- 
            candidate.net$likelihood_part + source.net$proposal.distr
            source.net$likelihood <- 
            source.net$likelihood_part + candidate.net$proposal.distr
            second.adapt.phase_net$acceptance_saved[i] <-
                candidate.net$likelihood - source.net$likelihood
            candidate_edge <- 
            which(candidate.net$adjacency!=source.net$adjacency, 
            arr.ind = TRUE)
            if(candidate.net$edge_move=="reverse")
            {
                second.adapt.phase_net$iter_edges[candidate_edge[1,"row"],
                    candidate_edge[1,"col"], 1] <-
                second.adapt.phase_net$iter_edges[candidate_edge[1,"row"],
                    candidate_edge[1,"col"], 1] + 1
                second.adapt.phase_net$iter_edges[candidate_edge[2,"row"],
                    candidate_edge[2,"col"], 1] <-
                second.adapt.phase_net$iter_edges[candidate_edge[2,"row"],
                    candidate_edge[2,"col"], 1] + 1
            } else {
                second.adapt.phase_net$iter_edges[candidate_edge[,"row"],
                    candidate_edge[,"col"], 1] <- 
                second.adapt.phase_net$iter_edges[candidate_edge[,"row"],
                    candidate_edge[,"col"], 1] + 1
            } # end if else (candidate.net$edge_move=="reverse")
      
            u <- log(stats::runif(1))
            if (u < second.adapt.phase_net$acceptance_saved[i])
            {
                if(candidate.net$edge_move=="add")
                {
                    second.adapt.phase_net$iter_edges[candidate_edge[,"row"],
                        candidate_edge[,"col"], 2] <-
                    second.adapt.phase_net$iter_edges[candidate_edge[,"row"],
                        candidate_edge[,"col"], 2] + 1
                } else if (candidate.net$edge_move=="reverse")
                { second.adapt.phase_net$iter_edges[
                candidate_edge[source.net$adjacency[candidate_edge]==0,"row"],
                candidate_edge[source.net$adjacency[candidate_edge]==0,"col"],
                2] <- 
                second.adapt.phase_net$iter_edges[
                candidate_edge[source.net$adjacency[candidate_edge]==0,"row"],
                candidate_edge[source.net$adjacency[candidate_edge]==0,"col"]
                ,2] + 1
                } # end if else if (candidate.net$edge_move=="add")
                source.net <- candidate.net
                beta.source$prior <- source.net$prior
            } else {
                if(candidate.net$edge_move=="reverse")
                    { second.adapt.phase_net$iter_edges[
                    candidate_edge[source.net$adjacency[candidate_edge]==1,
                    "row"],candidate_edge[source.net$adjacency[
                    candidate_edge]==1,"col"],2] <-
                    second.adapt.phase_net$iter_edges[candidate_edge[
                    source.net$adjacency[candidate_edge]==1,"row"],
                    candidate_edge[source.net$adjacency[candidate_edge]==1,
                    "col"], 2] + 1
                 } else if(candidate.net$edge_move=="delete")
                 { second.adapt.phase_net$iter_edges[candidate_edge[,"row"],
                     candidate_edge[,"col"], 2] <-
                     second.adapt.phase_net$iter_edges[candidate_edge[,"row"], 
                         candidate_edge[,"col"], 2] + 1
                 }# end if else if(candidate.net$edge_move=="reverse")
            } # end if else(u < second.adapt.phase_net$acceptance_saved[i])
            second.adapt.phase_net$nets[[i]] <- source.net
        } else {
            candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                layers_def = layers_def, omics = omics, 
                BGe_score_all_configs_node = BGe_score_all_configs_node, 
                parent_set_combinations = parent_set_combinations)
                second.adapt.phase_net$acceptance_saved[i] <-
                candidate.net$acceptance
                candidate.net$BGe <-
                BGe_score(adjacency_matrix = candidate.net$adjacency, 
                    omics = omics, layers_def = layers_def, 
                    parent_set_combinations = parent_set_combinations, 
                    BGe_score_all_configs_node = BGe_score_all_configs_node)
                candidate.net$nbhd.size <- neighborhood_size(omics = omics,
                    net = candidate.net$adjacency, layers_def = layers_def, 
                    B_prior_mat = B_prior_mat)
                candidate.net$energy <- sum(epsilon(B_prior_mat = B_prior_mat, 
                    net = candidate.net$adjacency))
                candidate.net$prior <- 
                (-beta.source$value*candidate.net$energy) -
                    second.adapt.phase_net$partition_func_UB_beta_source
                candidate.net$likelihood_part <- 
                candidate.net$BGe + candidate.net$prior
      
                u <- log(stats::runif(1))
                if (u < second.adapt.phase_net$acceptance_saved[i])
                {
                    source.net <- candidate.net
                    beta.source$prior <- source.net$prior
                } # end if (u < acceptance_saved[i])
                second.adapt.phase_net$nets[[i]] <- source.net
            } # end if(method.choice=="MC3")
    
            beta.candidate <- list(value = stats::rnorm(1, 
                mean = beta.source$value, sd = beta_sd*constant), prior = c(),
                len = beta_sd*constant)
            if(beta.candidate$value < 0.5)
            {
                beta.candidate$value <- 0.5
            } # end if(beta.candidate$value < 0.5)
            partition_func_UB_beta_candidate <-
            sum(mapply(second.adapt.phase_net$energy_all_configs_node,
            FUN=function(x) matrixStats::logSumExp(-beta.candidate$value*x)))
            beta.candidate$prior <- (-beta.candidate$value*source.net$energy) -
                partition_func_UB_beta_candidate
            second.adapt.phase_net$acceptance_beta_saved[i] <-
            beta.candidate$prior - beta.source$prior
            
            u_beta <- log(stats::runif(1))
            if (u_beta < second.adapt.phase_net$acceptance_beta_saved[i])
            {
                beta.source <- beta.candidate
                second.adapt.phase_net$partition_func_UB_beta_source <-
                partition_func_UB_beta_candidate
            } # end if (u_beta <...
            second.adapt.phase_net$betas[[i]] <- beta.source
        } # end for(i in (length(second.adapt.phase_net$nets)+1)...
    return(second.adapt.phase_net)
}

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
