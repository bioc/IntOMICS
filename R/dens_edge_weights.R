#' Density plot of edge weights inferred by IntOMICS
#' @description
#' `dens_edge_weights` Creates density plot of edge weights. 
#' @param weighted_net_res list output from the weighted_net function.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_density
#' @importFrom stats quantile
#'
#' @examples
#' data(list=c("OMICS_mod_res", "BN_mod_res", "gene_annot", "TFtarg_mat", 
#' "PK"), package="IntOMICS")
#' res_weighted <- edge_weights(mcmc_res = BN_mod_res, burn_in = 10000, 
#'  thin = 500, edge_freq_thres = 0.3) 
#' weighted_net_res <- weighted_net(cpdag_weights = res_weighted, 
#'  gene_annot = gene_annot, PK = PK, OMICS_mod_res = OMICS_mod_res, 
#'  gene_ID = "gene_symbol", TFtargs = TFtarg_mat,
#'  B_prior_mat_weighted = B_prior_mat_weighted) 
#' dens_edge_weights(weighted_net_res)
#'
#' @return density plot of edge weights
#' @export
dens_edge_weights <- function(weighted_net_res)
{
  df <- data.frame(edge_weight = as.numeric(weighted_net_res$edge_list[,"weight"]))
  q3 <- quantile(df$edge_weight, 0.75)
  p <- ggplot(df, aes(x=edge_weight, y = ..scaled..)) +
    geom_density(col="dodgerblue")
  p+geom_vline(xintercept=q3, size=0.5, color="black", linetype = "dashed")
}
