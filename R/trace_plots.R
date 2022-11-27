#' Trace plots of MCMC simulation
#' @description
#' `trace_plots` Create trace plots of MCMC simulation and filter low reliable 
#' edges based on the edge_freq_thres parameter. 
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
#' @return MCMC simulation trace plots
#' @export
trace_plots <- function(mcmc_res, burn_in, thin, edge_freq_thres = NULL)
{
  if(!is(mcmc_res,'list') | !all(names(mcmc_res) %in% 
                               c("sampling.phase_res","B_prior_mat_weighted",
                                 "beta_tuning")))
  {
    message('Invalid input "mcmc_res". Must be named list with names 
          c("sampling.phase_res","B_prior_mat_weighted","beta_tuning").')  
  }
  
  if(!is(burn_in,'numeric') | !is(thin,'numeric') | 
     length(burn_in)>1 | length(thin)>1)
  {
    message('Invalid input. "burn_in" or "thin" must be numeric of length 1.')  
  }
  
  if(!is(edge_freq_thres,'NULL') & !is(edge_freq_thres,'numeric') | 
     length(edge_freq_thres)>1)
  {
    message('Invalid input "edge_freq_thres". 
          Must be "NULL" or numeric of length 1.')
  }
  
  df1 <- data.frame(beta = mapply(mcmc_res$beta_tuning,
                                  FUN=function(x) x$value), 
                    k=seq_len(length(mapply(mcmc_res$beta_tuning,
                                            FUN=function(x) x$value))), 
                    accept = 1)
  rms_strength <- abs(diff(mcmc_res$sampling.phase_res$rms))
  strength_threshold <- quantile(rms_strength, 0.75, na.rm = TRUE)
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
  total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
  
  plot(df1$beta ~ df1$k, type = "l", col= "darkblue", xlab = "iteration",
       ylab = "beta", main = "Beta values of adaptive MCMC")
  plot(total$strength.x ~ total$strength.y, xlab="MCMC run 2",
       ylab = "MCMC run 1", 
       main = "Consistency of edges posterior probabilities")
  abline(0,1, col="orange")
  plot(rms_strength, main="Convergence RMS strength (C.RMS.str)", pch = 18,
       col="gray30")
  abline(h=strength_threshold, col="#E69F00", lwd = 1.5)
  text(label = paste("3rd quartile of C.RMS.str = ",
                     round(strength_threshold,3),sep=""), x = 100, 
       y = strength_threshold+0.015, col="#E69F00")
}
