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
        mapply(tail(squared.jump_second.adapt.phase_net$betas,
        1001), FUN=function(list) list$value)
    if(length(unique(betas_check))>1)
    {
        betas_check <- colMeans(matrix((betas_check[-1] -
            betas_check[-1001])^2,nrow=200))
        reg_dat <- data.frame(beta_means = betas_check, iter = seq_len(5))
        model <- lm(beta_means ~ iter, data = reg_dat)
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
            mapply(tail(squared.jump_second.adapt.phase_net$betas,1001),
                FUN=function(list) list$value)
            if(length(unique(betas_check))>1)
            {
                betas_check <- colMeans(matrix((betas_check[-1] -
                    betas_check[-1001])^2,nrow=200))
                reg_dat <- data.frame(beta_means = betas_check, 
                    iter = seq_len(5))
                model <- lm(beta_means ~ iter, data = reg_dat)
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
