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
