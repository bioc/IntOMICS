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
#'        burn_in = 100000, thin = 500, gene_annot = gene_annot, 
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

#' Arrow of directed edges tuning
#' @description
#' `normalise` This function is from the ambient package. 
#' It is used to determine the position of directed edge arrows.
#' @param x numeric vector to be modified.
#' @param from numeric vector range of x.
#' @param to numeric vector range of normalised x.
#' @importFrom stats median
#'
#' @examples
#' x <- seq(1,10)
#' normalise(x, from = range(x), to = c(0, 1))
#'
#' @return Numeric vector
#' @export
normalise <- function (x, from = range(x), to = c(0, 1)) {
    x <- (x - from[1])/(from[2] - from[1])
    if (!identical(to, c(0, 1))) {
        x <- x * (to[2] - to[1]) + to[1]
    }
    x
}

#' Resulting network definition
#' @description
#' `edge_types` Defines the resulting network structure and determines 
#' the color scale for each modality. This is part of trace_plots.
#' @param mcmc_res list output from the BN_module function.
#' @param PK data.frame with known interactions.
#' @param gene_annot data.frame containing the entrez ID and corresponding 
#' gene symbol for conversion.
#' @param edge_list matrix indicating the interaction between nodes, 
#' the first column indicates the source node, the second column indicates 
#' the target node.
#' @param node_list character vector indicating the complete set of nodes 
#' in the resulting network structure.
#' @param OMICS_mod_res list output from the OMICS_module function.
#' @param edge_weights character vector includes either "MCMC_freq" to reflect
#' the edge weights frequency over the final set of network structures or 
#' "empB" to reflect the empirical biological knowledge estimated by IntOMICS.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @importFrom RColorBrewer brewer.pal
#'
#' @return List of 6 elements needed to plot the final regulatory network edges
#' @export
edge_types <- function(mcmc_res, PK = NULL, gene_annot, edge_list, node_list,
OMICS_mod_res, edge_weights, TFtargs = NULL)
{
    omics <- OMICS_mod_res$omics
    layers_def <- OMICS_mod_res$layers_def
    omics_meth_original <- OMICS_mod_res$omics_meth_original
  
    if(any(regexpr("EID:",node_list)>0))
    {
        if(!is.null(PK))
        {
            PK <- paste(PK$src_entrez, PK$dest_entrez, sep="_")
        } # end if(!is.null(PK))
        
        if(edge_weights=="empB")
        {
            edge_list[,"edge_type"] <- "empirical"
            
            if(!is.null(TFtargs))
            {
                TF_pk <- as.matrix(TFtargs[intersect(edge_list[,"to"],
                rownames(TFtargs)), intersect(edge_list[,"from"],
                colnames(TFtargs))])
                colnames(TF_pk) <- intersect(edge_list[,"from"], 
                colnames(TFtargs))
                if(ncol(TF_pk)>=1)
                {
                    TF_pk <- paste(colnames(TF_pk)[which(TF_pk==1, 
                    arr.ind = TRUE)[,2]], rownames(TF_pk)[which(TF_pk==1, 
                    arr.ind = TRUE)[,1]], sep="_")
                    edge_list[match(intersect(edge_list[,"edge"],TF_pk),
                    edge_list[,"edge"]),"edge_type"] <- "TF"
                } # end if(ncol(TF_pk)>=1)
            } # end if(!is.null(TFtargs))
            
            if(!is.null(PK))
            {
                edge_list[match(intersect(edge_list[,"edge"],PK),
                edge_list[,"edge"]), "edge_type"] <- "PK"
            } # end if(!is.null(PK))
            
            edge_list[,"weight"] <-
            round(as.numeric(unlist(lapply(seq_along(edge_list[,2]), 1, 
                FUN=function(row)
                mcmc_res$B_prior_mat_weighted[edge_list[row,"from"], 
                edge_list[row,"to"]]))), 2)
        } else {
            if(!is.null(PK))
            {
                edge_list[match(setdiff(edge_list[,"edge"],PK),
                edge_list[,"edge"]),"edge_type"] <- "new"
                edge_list[match(intersect(edge_list[,"edge"],PK),
                edge_list[,"edge"]), "edge_type"] <- "PK"
            } # end if(!is.null(PK))
        
        } # end if else (edge_weights=="empB")
    
        ge_cols <- RColorBrewer::brewer.pal(9, "Blues")
        ge_common <- intersect(unique(node_list),
            colnames(omics[[layers_def$omics[1]]]))
        omics_ge_gs <- as.matrix(omics[[layers_def$omics[1]]][,ge_common])
        colnames(omics_ge_gs) <- ge_common
    
        borders_ge_b1 <- unlist(lapply(strsplit(levels(cut(omics[[
            layers_def$omics[1]]][omics[[layers_def$omics[1]]]<=
            stats::median(omics[[layers_def$omics[1]]])],
            seq(from=min(omics[[layers_def$omics[1]]]), 
            to=stats::median(omics[[layers_def$omics[1]]]), length.out=5), 
            include.lowest = TRUE)),","),FUN=function(l) l[1]))
        borders_ge_b1[1] <- sub("[","(",borders_ge_b1[1], fixed = TRUE)
        borders_ge_b1 <- as.numeric(sub("(","",borders_ge_b1, fixed = TRUE))
        borders_ge_b2 <-
        unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1
            ]]][omics[[layers_def$omics[1]]]>
            stats::median(omics[[layers_def$omics[1]]])],
            seq(from=stats::median(omics[[layers_def$omics[1]]]), 
            to=max(omics[[layers_def$omics[1]]]), length.out=6), 
            include.lowest = TRUE)),","),FUN=function(l) l[1]))
        borders_ge_b2[1] <- sub("[","(",borders_ge_b2[1], fixed = TRUE)
        borders_ge_b2 <- as.numeric(sub("(","",borders_ge_b2, fixed = TRUE))
        borders_ge_b <- c(borders_ge_b1,borders_ge_b2)
    
        borders_ge_t1 <- as.numeric(sub("]","",
            unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1
            ]]][omics[[layers_def$omics[1]]]<=
            stats::median(omics[[layers_def$omics[1]]])], 
            seq(from=min(omics[[layers_def$omics[1]]]), 
            to=stats::median(omics[[layers_def$omics[1]]]), length.out=5), 
            include.lowest = TRUE)),","),FUN=function(l) l[2]))))
        borders_ge_t2 <- as.numeric(sub("]","",
            unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1
            ]]][omics[[layers_def$omics[1]]]>
            stats::median(omics[[layers_def$omics[1]]])], 
            seq(from=stats::median(omics[[layers_def$omics[1]]]), 
            to=max(omics[[layers_def$omics[1]]]), length.out=6), 
            include.lowest = TRUE)),","),FUN=function(l) l[2]))))
        borders_ge_t <- c(borders_ge_t1,borders_ge_t2)
        borders <- sort(unique(c(borders_ge_b,borders_ge_t)))
        expr_group <- cut(colMeans(omics_ge_gs), breaks = borders,
            include.lowest = TRUE, labels = FALSE)
        names(expr_group) <- colnames(omics_ge_gs)
        node_list <- matrix(data = c(node_list,
            as.numeric(expr_group[match(node_list, names(expr_group))])),
            nrow = length(node_list), dimnames = list(c(), c("label", "color")))
        ind_cols <- paste(paste("(", paste(borders[-length(borders)],
            borders[-1]),sep=""),"]",sep="")
        borders_cnv <- NULL
        borders_meth <- NULL
    
        if(any(mapply(FUN=function(mod)
            any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))
        {
          cnv_cols <- RColorBrewer::brewer.pal(11, "PiYG")
          cnv_common <- intersect(node_list[,"label"][regexpr("eid",
              node_list[,"label"])>0],
          colnames(omics[[names(which(mapply(FUN=function(mod)
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]))
          omics_cnv_gs <- as.matrix(omics[[names(which(mapply(FUN=function(mod)
              any(regexpr("eid:",colnames(mod))>0), 
              omics)==TRUE))]][,cnv_common])
          borders_cnv_b1 <- unlist(lapply(strsplit(levels(cut(omics
              [[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), 
              omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]<=0], 
              seq(from=min(omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]],
              na.rm = TRUE), to=0, length.out=6), 
              include.lowest = TRUE)),","),FUN=function(l) l[1]))
          borders_cnv_b1[1] <- sub("[","(",borders_cnv_b1[1], fixed = TRUE)
          borders_cnv_b1 <- as.numeric(sub("(","",borders_cnv_b1, fixed = TRUE))
          borders_cnv_b2 <- unlist(lapply(strsplit(levels(cut(omics
              [[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), 
              omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]>0], 
              seq(from=0, to=max(omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]],
              na.rm = TRUE), length.out=7), 
              include.lowest = TRUE)),","),FUN=function(l) l[1]))
          borders_cnv_b2[1] <- sub("[","(",borders_cnv_b2[1], fixed = TRUE)
          borders_cnv_b2 <- as.numeric(sub("(","",borders_cnv_b2, fixed = TRUE))
          borders_cnv_b <- c(borders_cnv_b1,borders_cnv_b2)
      
          borders_cnv_t1 <- as.numeric(sub("]","",
              unlist(lapply(strsplit(levels(cut(omics
              [[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), 
              omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]<=0], 
              seq(from=min(omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]],
              na.rm = TRUE), to=0, length.out=6), 
              include.lowest = TRUE)),","),FUN=function(l) l[2]))))
          borders_cnv_t2 <- as.numeric(sub("]","",
              unlist(lapply(strsplit(levels(cut(omics
              [[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), 
              omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]>0], 
              seq(from=0, to=max(omics[[names(which(mapply(FUN=function(mod) 
              any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]],
              na.rm = TRUE), length.out=7), 
              include.lowest = TRUE)),","),FUN=function(l) l[2]))))
          borders_cnv_t <- c(borders_cnv_t1,borders_cnv_t2)
          borders_cnv <- sort(unique(c(borders_cnv_b,borders_cnv_t)))
      
          cnv_group <- cut(colMeans(omics_cnv_gs, na.rm = TRUE), 
              breaks = borders_cnv, include.lowest = TRUE, labels = FALSE) + 
              length(ge_cols)
          names(cnv_group) <- colnames(omics_cnv_gs)
          node_list[regexpr("eid",node_list[,"label"])>0,"color"] <-
              as.numeric(cnv_group[match(node_list[regexpr("eid",
              node_list[,"label"])>0,"label"], names(cnv_group))])
      
          ge_cols <- c(ge_cols, cnv_cols)
        } # end if(any(mapply(FUN=function(mod)...
    
        if(any(mapply(omics,
        FUN=function(list) any(regexpr("eid:",colnames(list), 
        ignore.case = TRUE)<0))))
        {
            meth_cols <- RColorBrewer::brewer.pal(9, "YlOrRd")
            meth_common <-
                intersect(node_list[,"label"],colnames(omics_meth_original))
            omics_meth_gs <- as.matrix(omics_meth_original[,meth_common])
            colnames(omics_meth_gs) <- meth_common
            borders_meth_b1 <- unlist(lapply(strsplit(levels(cut(
                omics_meth_original[omics_meth_original<=0.5],
                seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5,
                length.out=5), include.lowest = TRUE)),","),
                FUN=function(l) l[1]))
            borders_meth_b1[1] <- sub("[","(",borders_meth_b1[1], fixed = TRUE)
            borders_meth_b1 <- as.numeric(sub("(","",borders_meth_b1, 
            fixed = TRUE))
            borders_meth_b2 <- unlist(lapply(strsplit(levels(cut(
                omics_meth_original[omics_meth_original>0.5], seq(from=0.5,
                to=max(omics_meth_original, na.rm = TRUE), length.out=6), 
                include.lowest = TRUE)),","),FUN=function(l) l[1]))
            borders_meth_b2[1] <- sub("[","(",borders_meth_b2[1], fixed = TRUE)
            borders_meth_b2 <- as.numeric(sub("(","",borders_meth_b2, 
                fixed = TRUE))
            borders_meth_b <- c(borders_meth_b1,borders_meth_b2)
            borders_meth_t1 <- as.numeric(sub("]","",
                unlist(lapply(strsplit(levels(cut(omics_meth_original[
                omics_meth_original<=0.5], seq(from=min(omics_meth_original,
                na.rm = TRUE), to=0.5, length.out=5), 
                include.lowest = TRUE)),","),FUN=function(l) l[2]))))
            borders_meth_t2 <- as.numeric(sub("]","",
                unlist(lapply(strsplit(levels(cut(omics_meth_original[
                omics_meth_original>0.5], seq(from=0.5,
                to=max(omics_meth_original, na.rm = TRUE), length.out=6),
                include.lowest = TRUE)),","),FUN=function(l) l[2]))))
            borders_meth_t <- c(borders_meth_t1,borders_meth_t2)
            borders_meth <- sort(unique(c(borders_meth_b,borders_meth_t)))
            meth_group <- cut(colMeans(omics_meth_gs, na.rm = TRUE), 
                breaks = borders_meth, include.lowest = TRUE, labels = FALSE) +
                length(ge_cols)
            names(meth_group) <- colnames(omics_meth_gs)
            node_list[!is.na(match(node_list[,"label"],
            colnames(omics_meth_original))),"color"] <- 
            as.numeric(meth_group[match(node_list[,"label"][!is.na(
                match(node_list[,"label"], colnames(omics_meth_original)))],
                names(meth_group))])
            ge_cols <- c(ge_cols, meth_cols)
        } # end if(any(mapply(omics,FUN=function(list)...
    
    } else {
        if(!is.null(PK))
        {
            PK_src_dest <- as.character(gene_annot$gene_symbol[match(
            PK$src_entrez,gene_annot$entrezID)])
            PK_src_dest[regexpr("eid",PK_src_dest)>0] <- tolower(
            as.character(gene_annot$gene_symbol[match(toupper(PK_src_dest[
            regexpr("eid",PK_src_dest)>0]), gene_annot$entrezID)]))
            PK_src_dest[is.na(PK_src_dest)] <- 
            PK$src_entrez[is.na(PK_src_dest)]
            PK <- paste(PK_src_dest, as.character(gene_annot$gene_symbol[
                match(PK$dest_entrez,gene_annot$entrezID)]), sep="_")
        } # end if(!is.null(PK))
        
        if(edge_weights=="empB")
        {
            edge_list[,"edge_type"] <- "empirical"
            targs_eid <- gene_annot$entrezID[match(edge_list[,"to"],
                gene_annot$gene_symbol)]
            TFs_eid <- gene_annot$entrezID[match(edge_list[,"from"],
                gene_annot$gene_symbol, nomatch = 0)]
            if(!is.null(TFtargs))
            {
                TF_pk <- as.matrix(TFtargs[intersect(targs_eid, rownames(TFtargs)),
                    intersect(TFs_eid, colnames(TFtargs))])
                colnames(TF_pk) <- intersect(TFs_eid, colnames(TFtargs))
                if(ncol(TF_pk)>=1)
                {
                    TF_pk <- paste(gene_annot$gene_symbol[match(colnames(TF_pk)
                        [which(TF_pk==1, arr.ind = TRUE)[,2]],
                        gene_annot$entrezID)], gene_annot$gene_symbol[match(
                        rownames(TF_pk)[which(TF_pk==1, arr.ind = TRUE)[,1]],
                        gene_annot$entrezID)], sep="_")
                    edge_list[match(intersect(edge_list[,"edge"],TF_pk),
                    edge_list[,"edge"]),"edge_type"] <- "TF"
                } # end if(ncol(TF_pk)>=1)
            } # end if(!is.null(TFtargs))
            
          
            if(!is.null(PK))
            {
                edge_list[match(intersect(edge_list[,"edge"],PK),
                edge_list[,"edge"]), "edge_type"] <- "PK"
            } # end if(!is.null(PK))
            
            rownames(mcmc_res$B_prior_mat_weighted)[!is.na(match(rownames(
            mcmc_res$B_prior_mat_weighted), gene_annot$entrezID))] <-
            gene_annot$gene_symbol[match(rownames(
            mcmc_res$B_prior_mat_weighted), gene_annot$entrezID, nomatch = 0)]
            rownames(mcmc_res$B_prior_mat_weighted)[!is.na(match(toupper(
            rownames(mcmc_res$B_prior_mat_weighted)), gene_annot$entrezID))] <-
            tolower(gene_annot$gene_symbol[match(toupper(rownames(
            mcmc_res$B_prior_mat_weighted)), gene_annot$entrezID, nomatch = 0)])
            colnames(mcmc_res$B_prior_mat_weighted) <-
            rownames(mcmc_res$B_prior_mat_weighted)
            edge_list[,"weight"] <-
            round(as.numeric(unlist(lapply(seq_along(edge_list[,2]),
                FUN=function(row) mcmc_res$B_prior_mat_weighted[
                edge_list[row,"from"],edge_list[row,"to"]]))),2)
        } else {
            if(!is.null(PK))
            {
                edge_list[match(setdiff(edge_list[,"edge"],PK),
                edge_list[,"edge"]),"edge_type"] <- "new"
                edge_list[match(intersect(edge_list[,"edge"],PK),
                edge_list[,"edge"]), "edge_type"] <- "PK"
            } # end if(!is.null(PK))
        } # end if else (edge_weights=="empB")
    
        ge_cols <- RColorBrewer::brewer.pal(9, "Blues")
        ge_common <- intersect(gene_annot$entrezID[match(unique(node_list),
            gene_annot$gene_symbol)], colnames(omics[[layers_def$omics[1]]]))
        omics_ge_gs <- as.matrix(omics[[layers_def$omics[1]]][,ge_common])
        colnames(omics_ge_gs) <- gene_annot$gene_symbol[match(ge_common,
            gene_annot$entrezID)]
        borders_ge_b1 <- unlist(lapply(strsplit(levels(cut(
            omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]<=
            stats::median(omics[[layers_def$omics[1]]])],
            seq(from=min(omics[[layers_def$omics[1]]]), 
            to=stats::median(omics[[layers_def$omics[1]]]), length.out=5), 
            include.lowest = TRUE)),","),FUN=function(l) l[1]))
        borders_ge_b1[1] <- sub("[","(",borders_ge_b1[1], fixed = TRUE)
        borders_ge_b1 <- as.numeric(sub("(","",borders_ge_b1, fixed = TRUE))
        borders_ge_b2 <- unlist(lapply(strsplit(levels(cut(
            omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]>
            stats::median(omics[[layers_def$omics[1]]])], 
            seq(from=stats::median(omics[[layers_def$omics[1]]]), 
            to=max(omics[[layers_def$omics[1]]]), length.out=6), 
            include.lowest = TRUE)),","),FUN=function(l) l[1]))
        borders_ge_b2[1] <- sub("[","(",borders_ge_b2[1], fixed = TRUE)
        borders_ge_b2 <- as.numeric(sub("(","",borders_ge_b2, fixed = TRUE))
        borders_ge_b <- c(borders_ge_b1,borders_ge_b2)
        
        borders_ge_t1 <- as.numeric(sub("]","", unlist(lapply(strsplit(
        levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]<=
            stats::median(omics[[layers_def$omics[1]]])],
            seq(from=min(omics[[layers_def$omics[1]]]), 
            to=stats::median(omics[[layers_def$omics[1]]]), length.out=5), 
            include.lowest = TRUE)),","),FUN=function(l) l[2]))))
        borders_ge_t2 <- as.numeric(sub("]","", unlist(lapply(
            strsplit(levels(cut(omics[[layers_def$omics[1
            ]]][omics[[layers_def$omics[1]]]>
            stats::median(omics[[layers_def$omics[1]]])], 
            seq(from=stats::median(omics[[layers_def$omics[1]]]), 
            to=max(omics[[layers_def$omics[1]]]), length.out=6), 
            include.lowest = TRUE)),","),FUN=function(l) l[2]))))
        borders_ge_t <- c(borders_ge_t1,borders_ge_t2)
        borders <- sort(unique(c(borders_ge_b,borders_ge_t)))
        expr_group <- cut(colMeans(omics_ge_gs), breaks = borders,
            include.lowest = TRUE, labels = FALSE)
        names(expr_group) <- colnames(omics_ge_gs)
        node_list <- matrix(data = c(node_list,
            as.numeric(expr_group[match(node_list, names(expr_group))])),
            nrow = length(node_list), dimnames = list(c(), c("label", "color")))
        ind_cols <- paste(paste("(", paste(borders[-length(borders)],
            borders[-1]),sep=""),"]",sep="")
        borders_cnv <- NULL
        borders_meth <- NULL
    
        if(any(mapply(FUN=function(mod)
        any(regexpr("eid:",colnames(mod))>0), 
        omics)==TRUE))
        {
            cnv_cols <- RColorBrewer::brewer.pal(11, "PiYG")
            cnv_common <- intersect(tolower(gene_annot$entrezID[
                match(toupper(node_list[is.na(node_list[,"color"]),"label"]), 
                gene_annot$gene_symbol)]),
                colnames(omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:", colnames(mod))>0), omics)==TRUE))]]))
            omics_cnv_gs <-
                as.matrix(omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:",colnames(mod))>0), 
                omics)==TRUE))]][,cnv_common])
            colnames(omics_cnv_gs) <- 
                node_list[node_list[,"label"]==tolower(node_list[,"label"]),
                "label"][gene_annot$entrezID[match(toupper(node_list[node_list[,
                "label"]==tolower(node_list[,"label"]),"label"]), 
                gene_annot$gene_symbol)] %in% toupper(colnames(omics[[names(
                which(mapply(FUN=function(mod) any(regexpr("eid:",
                colnames(mod))>0), omics)==TRUE))]]))]
            borders_cnv_b1 <- unlist(lapply(strsplit(levels(cut(
                omics[[names(which(mapply(FUN=function(mod)
                any(regexpr("eid:",colnames(mod))>0), 
                omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:",colnames(mod))>0),
                omics)==TRUE))]]<=0],
                seq(from=min(omics[[names(which(mapply(FUN=function(mod)
                any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]],
                na.rm = TRUE), to=0, length.out=6), 
                include.lowest = TRUE)),","), FUN=function(l) l[1]))
            borders_cnv_b1[1] <- sub("[","(",borders_cnv_b1[1], fixed = TRUE)
            borders_cnv_b1 <- as.numeric(sub("(","",borders_cnv_b1, 
                fixed = TRUE))
            borders_cnv_b2 <- unlist(lapply(strsplit(levels(cut(
                omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:", colnames(mod))>0), 
                omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:", colnames(mod))>0),
                omics)==TRUE))]]>0], seq(from=0,
                to=max(omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:", colnames(mod))>0), omics)==TRUE))]],
                na.rm = TRUE), length.out=7), include.lowest = TRUE)),","),
                FUN=function(l) l[1]))
            borders_cnv_b2[1] <- sub("[","(",borders_cnv_b2[1], fixed = TRUE)
            borders_cnv_b2 <- as.numeric(sub("(","",borders_cnv_b2, 
                fixed = TRUE))
            borders_cnv_b <- c(borders_cnv_b1,borders_cnv_b2)
            borders_cnv_t1 <- as.numeric(sub("]","", unlist(lapply(
                strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:",colnames(mod))>0), 
                omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:",colnames(mod))>0),
                omics)==TRUE))]]<=0], seq(from=min(omics[[names(which(mapply(
                FUN=function(mod) any(regexpr("eid:",colnames(mod))>0),
                omics)==TRUE))]], na.rm = TRUE), to=0, length.out=6),
                include.lowest = TRUE)),","), FUN=function(l) l[2]))))
            borders_cnv_t2 <- as.numeric(sub("]","",
                unlist(lapply(strsplit(levels(cut(omics[[
                names(which(mapply(FUN=function(mod)
                any(regexpr("eid:",colnames(mod))>0), 
                omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]>0], 
                seq(from=0, to=max(omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]],
                na.rm = TRUE), length.out=7), include.lowest = TRUE)),","),
                FUN=function(l) l[2]))))
            borders_cnv_t <- c(borders_cnv_t1,borders_cnv_t2)
            borders_cnv <- sort(unique(c(borders_cnv_b,borders_cnv_t)))
            cnv_group <- cut(colMeans(omics_cnv_gs, na.rm = TRUE), 
            breaks = borders_cnv, include.lowest = TRUE, 
            labels = FALSE) + length(ge_cols)
            names(cnv_group) <- colnames(omics_cnv_gs)
            node_list[node_list[,"label"]==tolower(node_list[,"label"]),
            "color"][gene_annot$entrezID[match(toupper(node_list[node_list[,
            "label"]==tolower(node_list[,"label"]),"label"]), 
            gene_annot$gene_symbol)] %in% toupper(colnames(
            omics[[names(which(mapply(FUN=function(mod)
            any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]))] <-
                as.numeric(cnv_group[match(node_list[node_list[,
                "label"]==tolower(node_list[,"label"]),
                "label"][gene_annot$entrezID[match(toupper(
                node_list[node_list[,"label"]==tolower(node_list[,
                "label"]),"label"]), gene_annot$gene_symbol)] %in%
                toupper(colnames(omics[[names(which(mapply(FUN=function(mod) 
                any(regexpr("eid:",colnames(mod))>0), omics)==TRUE))]]))], 
                names(cnv_group))])
            ge_cols <- c(ge_cols, cnv_cols)
        } # end if(any(mapply(FUN=function(mod)...
    
        if(any(mapply(omics,FUN=function(list)
            any(regexpr("eid:",colnames(list), ignore.case = TRUE)<0))))
        {
            meth_cols <- RColorBrewer::brewer.pal(9, "YlOrRd")
            meth_common <-
                intersect(node_list[,"label"],colnames(omics_meth_original))
                omics_meth_gs <- as.matrix(omics_meth_original[,meth_common])
                colnames(omics_meth_gs) <- meth_common
            borders_meth_b1 <-
                unlist(lapply(strsplit(levels(cut(
                omics_meth_original[omics_meth_original<=0.5], 
                seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5,
                length.out=5), include.lowest = TRUE)),","), 
                FUN=function(l) l[1]))
            borders_meth_b1[1] <- sub("[","(",borders_meth_b1[1], fixed = TRUE)
            borders_meth_b1 <- as.numeric(sub("(","",borders_meth_b1, 
                fixed = TRUE))
            borders_meth_b2 <-  unlist(lapply(strsplit(levels(cut(
                omics_meth_original[omics_meth_original>0.5], 
                seq(from=0.5, to=max(omics_meth_original, na.rm = TRUE),
                length.out=6), include.lowest = TRUE)),","), 
                FUN=function(l) l[1]))
            borders_meth_b2[1] <- sub("[","(",borders_meth_b2[1], fixed = TRUE)
            borders_meth_b2 <- as.numeric(sub("(","",borders_meth_b2, 
                fixed = TRUE))
            borders_meth_b <- c(borders_meth_b1,borders_meth_b2)
            borders_meth_t1 <- as.numeric(sub("]","",
                unlist(lapply(strsplit(levels(cut(
                omics_meth_original[omics_meth_original<=0.5], 
                seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5,
                length.out=5), include.lowest = TRUE)),","), 
                FUN=function(l) l[2]))))
            borders_meth_t2 <- as.numeric(sub("]","",
                unlist(lapply(strsplit(levels(cut(
                omics_meth_original[omics_meth_original>0.5], 
                seq(from=0.5, to=max(omics_meth_original, na.rm = TRUE),
                length.out=6), include.lowest = TRUE)),","), 
                FUN=function(l) l[2]))))
            borders_meth_t <- c(borders_meth_t1, borders_meth_t2)
            borders_meth <- sort(unique(c(borders_meth_b,borders_meth_t)))
            meth_group <- cut(colMeans(omics_meth_gs, na.rm = TRUE), 
                breaks = borders_meth, include.lowest = TRUE, labels = FALSE) +
                length(ge_cols)
            names(meth_group) <- colnames(omics_meth_gs)

            node_list[!is.na(match(node_list[,"label"],
            colnames(omics_meth_original))),"color"] <- 
                as.numeric(meth_group[match(node_list[,"label"][
                !is.na(match(node_list[,"label"],
                colnames(omics_meth_original)))], names(meth_group))])
            ge_cols <- c(ge_cols, meth_cols)
        } # end if(any(mapply(omics,FUN=function(list)...
    } # end if else(any(regexpr("EID:",node_list)>0))
    return(list(edge_list = edge_list, node_list = node_list, 
        borders_GE = borders, borders_CNV = borders_cnv, 
        borders_METH = borders_meth, node_palette = ge_cols))
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
#'      burn_in = 100000, thin = 500, gene_annot = gene_annot, PK = PK, 
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
