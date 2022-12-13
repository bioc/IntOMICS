#' Density plot of edge weights inferred by IntOMICS
#' @description
#' `dens_edge_weights` Creates density plot of edge weights. 
#' @param net list output from the weighted_net function.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 after_stat
#' @importFrom stats quantile
#' @importFrom rlang .data
#'
#' @examples
#' data(list=c("OMICS_mod_res", "BN_mod_res", "gene_annot", "TFtarg_mat", 
#' "PK"), package="IntOMICS")
#' res_weighted <- edge_weights(mcmc_res = BN_mod_res, burn_in = 10000, 
#'  thin = 500, edge_freq_thres = 0.3) 
#' weighted_net_res <- weighted_net(cpdag_weights = res_weighted, 
#'  gene_annot = gene_annot, PK = PK, OMICS_mod_res = OMICS_mod_res, 
#'  gene_ID = "gene_symbol", TFtargs = TFtarg_mat,
#'  B_prior_mat_weighted = BN_mod_res@B_prior_mat_weighted)
#' dens_edge_weights(weighted_net_res)
#'
#' @return density plot of edge weights
#' @export
dens_edge_weights <- function(net)
{
  if(!is(net,'list') | 
     is(names(net),'NULL'))
  {
    message('Invalid input "net". Must be named list, 
              output from weighted_net().')
  }
  
  df <- data.frame(edge_weight = as.numeric(net$edge_list[,"weight"]))
  q3 <- quantile(df$edge_weight, 0.75)
  p <- ggplot(df, aes(x=.data$edge_weight, y = after_stat(.data$scaled))) +
    geom_density(col="dodgerblue")
  p+geom_vline(xintercept=q3, size=0.5, color="black", linetype = "dashed")
}
