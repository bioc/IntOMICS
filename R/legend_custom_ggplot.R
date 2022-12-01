#' Node color legend
#' @description
#' `legend_custom_ggplot` The color scale for each modality.
#' @param net list output from the trace_plots function.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot geom_rect geom_text theme_minimal theme annotate
#' @importFrom ggplot2 ylim
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom ggplot2 element_blank
#' @importFrom methods is
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
#' legend_custom_ggplot(weighted_net_res)
#'
#' @return Figure with color key
#' @export
legend_custom_ggplot <- function(net)
{
    if(!is(net,'list') | 
      is(names(net),'NULL'))
    {
      message('Invalid input "net". Must be named list, 
              output from weighted_net().')
    }
  
    xl <- 1.2;yb <- 1;xr <- 1.3;yt <- 2
    df <- data.frame(x_min = head(seq(yb,yt,(yt-yb)/9),-1),
                     x_max = tail(seq(yb,yt,(yt-yb)/9),-1),
                     y_min = xl, y_max = xr,
                     col = factor(brewer.pal(9, "Blues")),
                     modality = "GE")
    df <- rbind(df,df[nrow(df),])
    df$val <- net$borders_GE
    df$x_lab <- df$x_min
    df$x_lab[nrow(df)] <- df$x_max[nrow(df)]
    
    if(!is.null(net$borders_CNV))
    {
      df2 <- data.frame(x_min = head(seq(yb, yt, (yt-yb)/11), -1),
                        x_max = tail(seq(yb ,yt , (yt-yb)/11), -1),
                        y_min = xl+0.15, y_max = xr+0.15,
                        col = factor(brewer.pal(11, "PiYG")),
                        modality = "CNV") 
      df2 <- rbind(df2,df2[nrow(df2),])
      df2$val <- round(net$borders_CNV,2)
      df2$x_lab <- df2$x_min
      df2$x_lab[nrow(df2)] <- df2$x_max[nrow(df2)]
      df <- rbind(df, df2)
    } # end if(!is.null(net$borders_CNV))
    
    if(!is.null(net$borders_METH))
    {
      df3 <- data.frame(x_min = head(seq(yb, yt, (yt-yb)/9), -1),
                        x_max = tail(seq(yb ,yt , (yt-yb)/9), -1),
                        y_min = xl+0.2, y_max = xr+0.2,
                        col = factor(brewer.pal(9, "YlOrRd")),
                        modality = "METH")
      df3 <- rbind(df3,df3[nrow(df3),])
      df3$val <- round(net$borders_METH,2)
      df3$x_lab <- df3$x_min
      df3$x_lab[nrow(df3)] <- df3$x_max[nrow(df3)]
      df <- rbind(df, df3)
    } # end if(!is.null(net$borders_METH))
    
    ggplot(df, aes(xmin = x_min, 
                   xmax = x_max, 
                   ymin = y_min, 
                   ymax = y_max)) +
      geom_rect(fill = df$col, colour = "grey50") +
      geom_text(aes(x = x_lab, y = y_min-0.01, label = val), size=4, data = df, 
                check_overlap = TRUE, angle = 45) +
      theme_minimal() + theme(panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), 
                              panel.background = element_blank(),
                              axis.title=element_blank(),
                              axis.text=element_blank(),
                              axis.ticks=element_blank()) +
      annotate("text", x = 0.94, y = unique(df$y_min+(xr-xl)/2), 
               label = unique(df$modality), size=5)+
      ylim(min(df$y_min)-0.05, max(df$y_max)+0.05)
      
}
