#' @name MCMC_sapling_res
#' @aliases MCMC_sapling_res-class
#' @title The MCMC_sapling_res class
#' @description Container of an MCMC sampling phase results generated 
#' by the function \code{\link[IntOMICS]{bn_module}}.
#' 
#' @slot estimated_beta Numeric, estimated value of
#' hyperparameter beta
#' @slot estimated_len Numeric, estimated width of the sampling 
#' interval for hyperparameter beta
#' @slot B_prior_mat_weighted Empirical biological knowledge matrix,
#' interactions from the biological prior knowledge and TFs-target 
#' interactions are constant (unless if "TFBS_belief" is not equal
#' to "woPKGE_belief").
#' @slot CPDAGs_sim1 List of CPDAGs from the first independent MCMC 
#' simulation (thinned DAGs from the MCMC simulation converted 
#' into CPDAGs, ducplicated CPDAGs discarded)
#' @slot CPDAGs_sim2 List of CPDAGs from the second independent MCMC 
#' simulation (thinned DAGs from the MCMC simulation converted 
#' into CPDAGs, ducplicated CPDAGs discarded)
#' @slot beta_tuning Matrix of results from adaptive phases that
#' contains hyperparameter beta tuning
#' \describe{
#' \item{\code{value}}{trace of hyperparameter beta}
#' \item{len}{trace of width of the sampling interval for hyperparameter 
#' beta}}
#' @slot rms Numeric, trace of root mean square used for c_rms
#' measure to evaluate the convergence of MCMC simulation
#'
#' @examples
#'
#' # A MCMC_sapling_res object created by the bn_module function.
#' data("OMICS_mod_res", package="IntOMICS")
#' if(interactive()){BN_mod_res <- bn_module(burn_in = 500, 
#'     thin = 20, OMICS_mod_res = OMICS_mod_res, 
#'     minseglen = 5, len = 5, prob_mbr = 0.07)}
#' 
#' @export

MCMC_sapling_res <- setClass("MCMC_sapling_res", 
    slots = c(estimated_beta = "numeric", estimated_len = "numeric",
    B_prior_mat_weighted = "matrix", beta_tuning = "matrix",
    CPDAGs_sim1 = "list", CPDAGs_sim2 = "list", rms = "numeric"))
