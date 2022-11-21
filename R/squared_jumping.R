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
      
            u <- log(runif(1))
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
      
                u <- log(runif(1))
                if (u < second.adapt.phase_net$acceptance_saved[i])
                {
                    source.net <- candidate.net
                    beta.source$prior <- source.net$prior
                } # end if (u < acceptance_saved[i])
                second.adapt.phase_net$nets[[i]] <- source.net
            } # end if(method.choice=="MC3")
    
            beta.candidate <- list(value = rnorm(1, 
                mean = beta.source$value, sd = beta_sd*constant), prior = c(),
                len = beta_sd*constant)
            if(beta.candidate$value < 0.5)
            {
                beta.candidate$value <- 0.5
            } # end if(beta.candidate$value < 0.5)
            partition_func_UB_beta_candidate <-
            sum(mapply(second.adapt.phase_net$energy_all_configs_node,
            FUN=function(x) logSumExp(-beta.candidate$value*x)))
            beta.candidate$prior <- (-beta.candidate$value*source.net$energy) -
                partition_func_UB_beta_candidate
            second.adapt.phase_net$acceptance_beta_saved[i] <-
            beta.candidate$prior - beta.source$prior
            
            u_beta <- log(runif(1))
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
