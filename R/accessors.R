#' Estimated beta accessor
#' @description
#' `estimated_beta` This is accessor function for MCMC_sapling_res-class.
#' @param x MCMC_sapling_res-class, output from the bn_module function
#'
#' @examples
#' if(interactive()){data("BN_mod_res", package="IntOMICS")
#' estimated_beta(BN_mod_res)}
#'
#' @return Numeric, trace of root mean square used for c_rms measure
#' @export
estimated_beta <- function(x) x@estimated_beta

#' Estimated len accessor
#' @description
#' `estimated_len` This is accessor function for MCMC_sapling_res-class.
#' @param x MCMC_sapling_res-class, output from the bn_module function
#'
#' @examples
#' if(interactive()){data("BN_mod_res", package="IntOMICS")
#' estimated_len(BN_mod_res)}
#'
#' @return Numeric, width of the sampling interval for hyperparameter beta
#' @export
estimated_len <- function(x) x@estimated_len

#' Empirical biological knowledge accessor
#' @description
#' `B_prior_mat_weighted` This is accessor function for MCMC_sapling_res-class.
#' @param x MCMC_sapling_res-class, output from the bn_module function
#'
#' @examples
#' if(interactive()){data("BN_mod_res", package="IntOMICS")
#' B_prior_mat_weighted(BN_mod_res)}
#'
#' @return Matrix, empirical biological knowledge
#' @export
B_prior_mat_weighted <- function(x) x@B_prior_mat_weighted

#' Beta tuning accessor
#' @description
#' `beta_tuning` This is accessor function for MCMC_sapling_res-class.
#' @param x MCMC_sapling_res-class, output from the bn_module function
#'
#' @examples
#' if(interactive()){data("BN_mod_res", package="IntOMICS")
#' beta_tuning(BN_mod_res)}
#'
#' @return Matrix, results from adaptive phases that contains hyperparameter
#' beta tuning
#' @export
beta_tuning <- function(x) x@beta_tuning

#' CPDAGs from the first simulation accessor
#' @description
#' `CPDAGs_sim1` This is accessor function for MCMC_sapling_res-class.
#' @param x MCMC_sapling_res-class, output from the bn_module function
#'
#' @examples
#' if(interactive()){data("BN_mod_res", package="IntOMICS")
#' CPDAGs_sim1(BN_mod_res)}
#'
#' @return List, CPDAGs from the first independent MCMC simulation
#' @export
CPDAGs_sim1 <- function(x) x@CPDAGs_sim1

#' CPDAGs from the second simulation accessor
#' @description
#' `CPDAGs_sim2` This is accessor function for MCMC_sapling_res-class.
#' @param x MCMC_sapling_res-class, output from the bn_module function
#'
#' @examples
#' if(interactive()){data("BN_mod_res", package="IntOMICS")
#' CPDAGs_sim2(BN_mod_res)}
#'
#' @return List, CPDAGs from the second independent MCMC simulation
#' @export
CPDAGs_sim2 <- function(x) x@CPDAGs_sim2

#' c_rms trace accessor
#' @description
#' `rms` This is accessor function for MCMC_sapling_res-class.
#' @param x MCMC_sapling_res-class, output from the bn_module function
#'
#' @examples
#' if(interactive()){data("BN_mod_res", package="IntOMICS")
#' rms(BN_mod_res)}
#'
#' @return Numeric, trace of root mean square used for c_rms measure
#' @export
rms <- function(x) x@rms
