#' Regulatory network plot with edge labels
#' @description
#' `ggraph_weighted_net` Figure of the regulatory network.
#' @param net list output from the trace_plots function.
#' @param node_size numeric node size
#' @param node_label_size numeric node label size
#' @param edge_label_size numeric edge label size
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot geom_rect geom_text theme_minimal theme annotate
#' @importFrom utils head tail
#' @importFrom igraph edge
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom ggplot2 scale_colour_manual arrow ylim
#' @importFrom ggraph geom_node_text circle
#' @importFrom grid unit
#' @importFrom methods is
#' @importFrom cowplot align_plots plot_grid
#' @importFrom rlang .data
#'
#' @examples
#' if(interactive()){data(list=c("OMICS_mod_res", "BN_mod_res", "gene_annot", "TFtarg_mat", 
#' "PK"), package="IntOMICS")
#' res_weighted <- edge_weights(mcmc_res = BN_mod_res, burn_in = 10000, 
#'  thin = 500, edge_freq_thres = 0.3)
#' weighted_net_res <- weighted_net(cpdag_weights = res_weighted, 
#'  gene_annot = gene_annot, PK = PK, OMICS_mod_res = OMICS_mod_res, 
#'  gene_ID = "gene_symbol", TFtargs = TFtarg_mat,
#'  B_prior_mat_weighted = B_prior_mat_weighted(BN_mod_res))
#' library(ggraph)
#' ggraph_weighted_net(weighted_net_res)}
#'
#' @return Figure of weighted network
#' @export
ggraph_weighted_net <- function(net, node_size = 10, node_label_size = 4, 
                                edge_label_size = 4)
{
  if(!is(net, 'list') | 
     is(names(net),'NULL'))
  {
    message('Invalid input "net". Must be named list, 
            output from weighted_net().')
  }
  
  if(!is(node_size,'numeric') | !is(node_label_size,'numeric') | 
     !is(edge_label_size,'numeric') | length(node_size)>1 |
     length(node_label_size)>1 | length(edge_label_size)>1)
  {
    message('Invalid input. "node_size", "node_label_size", "edge_label_size",
            and must be numeric of length 1.')  
  }

  # regulatory network
  rn <- ggraph(net$net_weighted, layout = 'dh') + 
    geom_edge_link(aes(end_cap = circle(.data$node2.degree + 7, "pt"),
                     edge_color = edge, label = .data$weight), 
                 label_size = edge_label_size,
                 arrow = arrow(angle = 20, length = unit(0.1, "inches"),
                               ends = "last", type = "closed")) +
    geom_node_point(aes(color = factor(.data$color)), size = node_size) +
    scale_colour_manual(values = net$node_palette, guide = "none") +
    geom_node_text(aes(label = .data$label), size = node_label_size)
  
  leg <- legend_custom_ggplot(net = net)
  plot_grid(rn, leg, ncol = 1, rel_heights = c(3, 1))
}
