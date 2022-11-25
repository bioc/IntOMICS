#' Edge weights of MCMC simulation
#' @description
#' `edge_weights` Returns list of edges with corresponding posterior 
#' probabilities (possibly filtered low reliable edges). 
#' @param mcmc_res list output from the bn_module function.
#' @param burn_in numeric vector the minimal length of burn-in period 
#' of the MCMC simulation.
#' @param thin numeric vector thinning frequency of the resulting 
#' MCMC simulation.
#' @param edge_freq_thres numerical vector the quantile of all edge weights
#' used to filter the most reliable edges.
#' @importFrom bnlearn nodes
#' @importFrom graphics text abline
#' @importFrom stats quantile
#' @importFrom bnlearn custom.strength
#'
#' @examples
#' data("BN_mod_res", package="IntOMICS")
#' res_weighted <- trace_plots(mcmc_res = BN_mod_res, burn_in = 10000, 
#'        thin = 500, edge_freq_thres = 0.3) 
#'
#' @return data.frame with edges and corresponding edge weights;
#' edge_freq_thres used to filter relevant edges
#' @export
edge_weights <- function(mcmc_res, burn_in, thin, edge_freq_thres = NULL)
{
  if(!is.list(mcmc_res) | !all(names(mcmc_res) %in% 
                               c("sampling.phase_res","B_prior_mat_weighted",
                                 "beta_tuning")))
  {
    message('Invalid input "mcmc_res". Must be named list with names 
          c("sampling.phase_res","B_prior_mat_weighted","beta_tuning").')  
  }
  
  if(!is.numeric(burn_in) | !is.numeric(thin) | 
     length(burn_in)!=1 | length(thin)!=1)
  {
    message('Invalid input. "burn_in" or "thin" must be numeric of length 1.')  
  }
  
  if(!is.null(edge_freq_thres) & !is.numeric(edge_freq_thres) | 
     length(edge_freq_thres)>1)
  {
    message('Invalid input "edge_freq_thres". 
          Must be "NULL" or numeric of length 1.')
  }
  
  cpdag_f <- (burn_in/thin+1)
  cpdag_l <- length(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed1$cpdags)
  cpdags1 <- 
    unique(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed1$cpdags[
      seq(from = cpdag_f, to = cpdag_l)])
  cpdags2 <- 
    unique(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed2$cpdags[
      seq(from = cpdag_f, to = cpdag_l)])
  
  cpdag_weights1 <- custom.strength(cpdags1, 
                                    nodes = nodes(cpdags1[[1]]), weights = NULL)
  cpdag_weights2 <- custom.strength(cpdags2, 
                                    nodes = nodes(cpdags2[[1]]), weights = NULL)
  cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
  cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]
  
  cpdag_weights1$edge <- paste(cpdag_weights1$from, cpdag_weights1$to,
                               sep="_")
  cpdag_weights2$edge <- paste(cpdag_weights2$from, cpdag_weights2$to,
                               sep="_")
  cpdag_weights <- merge(cpdag_weights1, cpdag_weights2, by = "edge")
  cpdag_weights$strength <- round(rowMeans(cbind(cpdag_weights$strength.x,
                                                 cpdag_weights$strength.y)),2)
  if(!is.null(edge_freq_thres))
  {
    strength_quant <- quantile(x = cpdag_weights$strength, 
                               probs = edge_freq_thres)
    cpdag_weights <- cpdag_weights[cpdag_weights$strength >=
                                     strength_quant,]
  }
  return(cpdag_weights)
}
