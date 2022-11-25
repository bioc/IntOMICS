#' Node color
#' @description
#' `legend_custom` Determines the color scale for each modality.
#' @param net list output from the trace_plots function.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom utils head
#' @importFrom utils tail
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
#' legend_custom(weighted_net_res)
#'
#' @return Figure with color key
#' @export
legend_custom <- function(net)
{
    xl <- 1.2;yb <- 1;xr <- 1.3;yt <- 2
    par(oma=c(0,0,0,0))
    plot(NA,type="n",ann=FALSE,xlim=c(0.94,2),ylim=c(1.2,1.71),xaxt="n",
        yaxt="n", bty="n")
    text(x = 0.95,y = 1.25,labels = "GE")
    rect(head(seq(yb,yt,(yt-yb)/9),-1), 
        xl, tail(seq(yb,yt,(yt-yb)/9),-1), xr,
        col=brewer.pal(9, "Blues"))
    text(x = unique(c(head(seq(yb,yt,(yt-yb)/9),-1), 
        tail(seq(yb,yt,(yt-yb)/9),-1))),y = 1.32,
        labels = round(net$borders_GE,2))
    if(!is.null(net$borders_CNV))
    {
        text(x = 0.94,y = 1.45,labels = "CNV")
        rect(head(seq(yb, yt, (yt-yb)/11), -1), xl+0.2,
            tail(seq(yb ,yt , (yt-yb)/11), -1), xr+0.2,
            col=brewer.pal(11, "PiYG"))
        text(x = unique(c(head(seq(yb,yt,(yt-yb)/11),-1), 
            tail(seq(yb,yt,(yt-yb)/11),-1))),y = 1.52,
            labels = round(net$borders_CNV,2))
    } # end if(!is.null(net$borders_CNV))
    if(!is.null(net$borders_METH))
    {
        rect(head(seq(yb, yt, (yt-yb)/9), -1), xl+0.4, 
            tail(seq(yb ,yt , (yt-yb)/9), -1), xr+0.4,
            col=brewer.pal(9, "YlOrRd"))
        text(x = unique(c(head(seq(yb,yt,(yt-yb)/9),-1), 
            tail(seq(yb,yt,(yt-yb)/9),-1))),y = 1.72,
            labels = round(net$borders_METH,2))
        text(x = 0.94,y = 1.65,labels = "METH")
    } # end if(!is.null(net$borders_METH))
}
