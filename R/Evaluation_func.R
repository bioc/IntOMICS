#' Trace plots of MCMC simulation & resulting network definition
#' @description
#' `trace_plots` Create trace plots of MCMC simulation and filter low reliable 
#' edges based on the edge_freq_thres parameter. 
#' Defines the resulting network structure and determines the color scale 
#' for each modality.
#' @param mcmc_res list output from the BN_module function.
#' @param burn_in numeric vector the minimal length of burn-in period 
#' of the MCMC simulation.
#' @param thin numeric vector thinning frequency of the resulting 
#' MCMC simulation.
#' @param figures_dir character vector the path of the folder to save figures.
#' @param gene_annot data.frame containing the entrez ID and corresponding gene
#' symbol for conversion.
#' @param PK data.frame with known interactions.
#' @param OMICS_mod_res list output from the OMICS_module function.
#' @param edge_weights character vector includes either "MCMC_freq" to reflect
#' the edge weights frequency over the final set of network structures or 
#' "empB" to reflect the empirical biological knowledge estimated by IntOMICS.
#' @param edge_freq_thres numerical vector the quantile of all edge weights
#' used to filter the most reliable edges.
#' @param gene_ID character vector includes either "gene_symbol" or "entrezID"
#' to reflect gene identifiers used in the final figure.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @importFrom bnlearn nodes custom.strength
#' @importFrom igraph degree
#' @importFrom graphics text abline
#' @importFrom igraph graph_from_edgelist V E as_ids
#' @importFrom stats quantile
#' @importFrom grDevices svg dev.off
#'
#' @examples
#' data(list=c("OMICS_mod_res", "BN_mod_res", "gene_annot", "TFtarg_mat", 
#' "PK"), package="IntOMICS")
#' res_weighted <- trace_plots(mcmc_res = BN_mod_res, figures_dir = "figures/", 
#'        burn_in = 10000, thin = 500, gene_annot = gene_annot, 
#'        PK = PK, OMICS_mod_res = OMICS_mod_res, gene_ID = "gene_symbol", 
#'        edge_freq_thres = 0.3, TFtargs = TFtarg_mat) 
#'
#' @return List of 7 elements needed to plot the final regulatory network
#' @export
trace_plots <- function(mcmc_res, burn_in, thin, figures_dir, gene_annot, PK=NULL,
OMICS_mod_res, edge_weights = "MCMC_freq", edge_freq_thres = NULL, gene_ID,
TFtargs = NULL)
{
    if(!(edge_weights %in% c("MCMC_freq","empB")))
    {
        message('edge_weights argument must be either "MCMC_freq" or "empB"')  
    }

    if(!dir.exists(figures_dir)){dir.create(figures_dir)}

    df1 <- data.frame(beta = mapply(mcmc_res$beta_tuning,
        FUN=function(x) x$value), k=seq_len(length(mapply(mcmc_res$beta_tuning,
        FUN=function(x) x$value))), accept = 1)
  
    rms_strength <- abs(diff(mcmc_res$sampling.phase_res$rms))
    strength_threshold <- stats::quantile(rms_strength, 0.75, na.rm = TRUE)
  
    cpdags1 <-
    unique(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed1$cpdags[seq(from 
    = (burn_in/thin+1), to =
    length(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed1$cpdags))])
    cpdags2 <-
    unique(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed2$cpdags[seq(from
    = (burn_in/thin+1), 
    to = length(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed2$cpdags))])
  
    cpdag_weights1 <- bnlearn::custom.strength(cpdags1, 
        nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
    cpdag_weights2 <- bnlearn::custom.strength(cpdags2, 
        nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
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
        strength_quant <- stats::quantile(x = cpdag_weights$strength, 
            probs = edge_freq_thres)
        cpdag_weights <- cpdag_weights[cpdag_weights$strength >=
            strength_quant,]
    }
    total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
  
    grDevices::jpeg(paste(figures_dir,"beta_values.jpeg",sep="/"))
    plot(df1$beta ~ df1$k, type = "l", col= "darkblue", xlab = "iteration",
        ylab = "beta", main = "Beta values of adaptive MCMC")
    grDevices::dev.off()
  
    grDevices::jpeg(paste(figures_dir,"post_prob_edges.jpeg",sep="/"))
    plot(total$strength.x ~ total$strength.y, xlab="MCMC run 2",
        ylab = "MCMC run 1", 
        main = "Consistency of edges posterior probabilities")
    graphics::abline(0,1, col="orange")
    grDevices::dev.off()
  
    grDevices::jpeg(paste(figures_dir,"convergence_RMS.jpeg",sep="/"))
    plot(rms_strength, main="Convergence RMS strength (C.RMS.str)", pch = 18,
        col="gray30")
    graphics::abline(h=strength_threshold, col="#E69F00", lwd = 1.5)
    graphics::text(label = paste("3rd quartile of C.RMS.str = ",
        round(strength_threshold,3),sep=""), x = 100, 
        y = strength_threshold+0.015, col="#E69F00")
    grDevices::dev.off()
  
    if(!is.null(PK))
    {
        PK <- PK[PK$src_entrez %in% unlist(lapply(OMICS_mod_res$omics,colnames)),]
        PK <- PK[PK$dest_entrez %in% unlist(lapply(OMICS_mod_res$omics,colnames)),]
    }

  
    if(gene_ID=="entrezID")
    {
        edge_list <- matrix(data = c(cpdag_weights$from.x, cpdag_weights$to.x,
            cpdag_weights$strength, rep(NA,length(cpdag_weights$strength)),
            rep(NA,length(cpdag_weights$strength))), 
            nrow = length(cpdag_weights$strength), dimnames = list(c(),
            c("from", "to", "weight", "edge_type", "edge")))
        node_list <- unique(c(edge_list[,"from"], edge_list[,"to"]))
        edge_list[,"edge"] <- paste(edge_list[,"from"], 
            edge_list[,"to"], sep="_")
        return_list <- edge_types(mcmc_res = mcmc_res, PK = PK, 
            gene_annot = gene_annot, edge_list = edge_list, 
            node_list = node_list, OMICS_mod_res = OMICS_mod_res, 
            edge_weights = edge_weights, TFtargs = TFtargs)
    } else if (gene_ID=="gene_symbol") {
        from <- as.character(gene_annot$gene_symbol[match(cpdag_weights$from.x,
            gene_annot$entrezID)])
        from[is.na(from)] <- cpdag_weights$from.x[is.na(from)]
        from[regexpr("eid",from)>0] <-
        tolower(as.character(gene_annot$gene_symbol[match(toupper(
            from[regexpr("eid", from)>0]), gene_annot$entrezID)]))
        to <- as.character(gene_annot$gene_symbol[match(cpdag_weights$to.x,
            gene_annot$entrezID)])
        edge_list <- matrix(data = c(from, to, cpdag_weights$strength,
            rep(NA,length(cpdag_weights$strength)),
            rep(NA,length(cpdag_weights$strength))),
            nrow = length(cpdag_weights$strength), dimnames = list(c(),
            c("from", "to", "weight", "edge_type", "edge")))
        node_list <- unique(c(edge_list[,"from"], edge_list[,"to"]))
        edge_list[,"edge"] <- paste(edge_list[,"from"], edge_list[,"to"],
            sep="_")
        return_list <- edge_types(mcmc_res = mcmc_res, PK = PK, 
            gene_annot = gene_annot, edge_list = edge_list, 
            node_list = node_list, OMICS_mod_res = OMICS_mod_res, 
            edge_weights = edge_weights, TFtargs = TFtargs)
    } else {
        message('gene_ID argument must be either "entrezID" or "gene_symbol"')
    } # end if else if else (gene_ID=="entrezID")
  
    net_weighted <-
    igraph::graph_from_edgelist(return_list$edge_list[,c("from","to")])
    igraph::V(net_weighted)$color <-
    return_list$node_list[match(igraph::as_ids(igraph::V(net_weighted)),
        return_list$node_list[,"label"]),"color"]
    palette <- return_list$node_palette
    names(palette) <- seq_len(length(palette))
    palette <- palette[unique(igraph::V(net_weighted)$color)]
    igraph::V(net_weighted)$label <-
    return_list$node_list[match(igraph::as_ids(igraph::V(net_weighted)),
        return_list$node_list[,"label"]),"label"]
    igraph::E(net_weighted)$edge <- return_list$edge_list[match(sub("|", "_",
        igraph::as_ids(igraph::E(net_weighted)), fixed = TRUE),
        return_list$edge_list[,"edge"]),"edge_type"]
    igraph::E(net_weighted)$weight <- return_list$edge_list[match(sub("|", "_",
        igraph::as_ids(igraph::E(net_weighted)), fixed = TRUE),
        return_list$edge_list[,"edge"]),"weight"]

    igraph::V(net_weighted)$degree <- igraph::degree(net_weighted, mode = "in")
    igraph::V(net_weighted)$degree <- normalise(igraph::V(net_weighted)$degree,
        to = c(3, 11))
  
    return(list(edge_list = return_list$edge_list, node_palette = palette,
        node_list = return_list$node_list, borders_GE = return_list$borders_GE,
        borders_CNV = return_list$borders_CNV, net_weighted = net_weighted,
        borders_METH = return_list$borders_METH))
}

#' Node color
#' @description
#' `legend_custom` Determines the color scale for each modality.
#' @param net list output from the trace_plots function.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par text rect
#' @importFrom utils head tail
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "gene_annot",
#' "OMICS_mod_res", "BN_mod_res"), package="IntOMICS")
#' res_weighted <- trace_plots(mcmc_res = BN_mod_res, figures_dir = "figures", 
#'      burn_in = 10000, thin = 500, gene_annot = gene_annot, PK = PK, 
#'      OMICS_mod_res = OMICS_mod_res, gene_ID = "gene_symbol", 
#'      edge_freq_thres = 0.3, TFtargs = TFtarg_mat) 
#' legend_custom(res_weighted)
#'
#' @return Figure with color key
#' @export
legend_custom <- function(net)
{
    xl <- 1.2;yb <- 1;xr <- 1.3;yt <- 2
    graphics::par(oma=c(0,0,0,0))
    plot(NA,type="n",ann=FALSE,xlim=c(0.94,2),ylim=c(1.2,1.71),xaxt="n",
        yaxt="n", bty="n")
    graphics::text(x = 0.95,y = 1.25,labels = "GE")
    graphics::rect(utils::head(seq(yb,yt,(yt-yb)/9),-1), 
        xl, utils::tail(seq(yb,yt,(yt-yb)/9),-1), xr,
        col=RColorBrewer::brewer.pal(9, "Blues"))
    graphics::text(x = unique(c(utils::head(seq(yb,yt,(yt-yb)/9),-1), 
        utils::tail(seq(yb,yt,(yt-yb)/9),-1))),y = 1.32,
        labels = round(net$borders_GE,2))
    if(!is.null(net$borders_CNV))
    {
        graphics::text(x = 0.94,y = 1.45,labels = "CNV")
        graphics::rect(utils::head(seq(yb, yt, (yt-yb)/11), -1), xl+0.2,
            utils::tail(seq(yb ,yt , (yt-yb)/11), -1), xr+0.2,
            col=RColorBrewer::brewer.pal(11, "PiYG"))
        graphics::text(x = unique(c(utils::head(seq(yb,yt,(yt-yb)/11),-1), 
            utils::tail(seq(yb,yt,(yt-yb)/11),-1))),y = 1.52,
            labels = round(net$borders_CNV,2))
    } # end if(!is.null(net$borders_CNV))
    if(!is.null(net$borders_METH))
    {
        graphics::rect(utils::head(seq(yb, yt, (yt-yb)/9), -1), xl+0.4, 
            utils::tail(seq(yb ,yt , (yt-yb)/9), -1), xr+0.4,
            col=RColorBrewer::brewer.pal(9, "YlOrRd"))
        graphics::text(x = unique(c(utils::head(seq(yb,yt,(yt-yb)/9),-1), 
            utils::tail(seq(yb,yt,(yt-yb)/9),-1))),y = 1.72,
            labels = round(net$borders_METH,2))
        graphics::text(x = 0.94,y = 1.65,labels = "METH")
    } # end if(!is.null(net$borders_METH))
}

#' Heatmap of empB - B
#' @description
#' `empB_heatmap` plot a heatmap with empB - B values (depicts the difference
#' between prior knowledge and the empirical knowledge)
#' @param mcmc_res list output from the BN_module function.
#' @param OMICS_mod_res list output from the OMICS_module function.
#' @param gene_annot data.frame containing the entrez ID and corresponding gene
#' symbol for conversion.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @importFrom gplots heatmap.2 bluered
#'
#' @examples
#' data(list=c("TFtarg_mat", "gene_annot", "OMICS_mod_res",
#' "BN_mod_res"), package="IntOMICS")
#' empB_heatmap(mcmc_res = BN_mod_res, OMICS_mod_res = OMICS_mod_res, 
#'     gene_annot = gene_annot, TFtargs = TFtarg_mat)
#'             
#' @return Figure heatmap
#' @export
empB_heatmap <- function(mcmc_res, OMICS_mod_res, gene_annot, TFtargs)
{
    mat <- mcmc_res$B_prior_mat_weighted - OMICS_mod_res$B_prior_mat
    mat <- mat[regexpr("EID:",rownames(mat))>0,
        regexpr("EID:",rownames(mat))>0]
    mat[which(OMICS_mod_res$B_prior_mat[regexpr("EID:",
        rownames(OMICS_mod_res$B_prior_mat))>0,regexpr("EID:",
        rownames(OMICS_mod_res$B_prior_mat))>0]==1)] <- NA
    TFtargs <- as.matrix(TFtargs[intersect(rownames(TFtargs),colnames(mat)),
        intersect(colnames(TFtargs),rownames(mat))])
    colnames(TFtargs) <- intersect(colnames(TFtargs),rownames(mat))
    if(ncol(TFtargs)>0)
    {
        inds <- which(TFtargs==1, arr.ind = TRUE)
        for(i in seq_len(nrow(inds)))
        {
            mat[colnames(TFtargs)[inds[i,2]],rownames(TFtargs)[inds[i,1]]] <- NA
        } # end for
    } # end if 
    rownames(mat) <- gene_annot$gene_symbol[match(rownames(mat),
        gene_annot$entrezID)]
    colnames(mat) <- rownames(mat)
    diag(mat) <- NA
    gplots::heatmap.2(mat, col=gplots::bluered, na.color="gray", Rowv = NA, 
        Colv = NA, trace = 'none', scale = 'none', main = "empB - B",
        dendrogram = "none")
}
