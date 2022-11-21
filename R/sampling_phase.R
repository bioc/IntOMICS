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
#'     thin = 500, minseglen = 500, burn_in = 10000, 
#'     annot = OMICS_mod_res$annot)}
#'
#' @return List of 2 elements: sampling phase result; RMS used to evaluate 
#' MCMC convergence
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
