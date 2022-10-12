#' OMICS_module
#' @description
#' `OMICS_module` data preprocessing + B_prior_mat definition + 
#'  partition function upper bound estimation + 
#'  all possible parent sets per node definition + 
#'  BGe score computation for all possible parent sets
#'
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and 
#' features in columns.
#' @param PK data.frame with known interactions.
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @param gene_annot data.frame containing the entrez ID and corresponding 
#' gene symbol for conversion.
#' @param lm_METH logical asking whether to use linear regression to filter
#' methylation data (default=TRUE).
#' @param r_squared_thres numeric vector to define the R^2 used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.3).
#' @param p_val_thres numeric vector to define the p-value used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.05).
#' @param TFBS_belief numeric vector to define the belief concerning the TF 
#' and its target interaction (default=0.75).
#' @param nonGE_belief numeric vector to define the belief concerning
#' interactions of features except GE-GE (default=0.5).
#' @param woPKGE_belief numeric vector to define the belief concerning GE-GE
#' interactions without prior knowledge (default=0.5).
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", 
#' "gene_annot"), package="IntOMICS")
#' OMICS_mod_res <- OMICS_module(omics = omics, PK = PK, 
#'     layers_def = layers_def, TFtargs = TFtarg_mat, annot = annot, 
#'     gene_annot = gene_annot, r_squared_thres = 0.3, lm_METH = TRUE)
#'
#' @return List of 6 elements needed to init MCMC simulation          
#' @export
OMICS_module <- function(omics, PK=NULL, layers_def, TFtargs=NULL, annot=NULL, 
lm_METH = TRUE, r_squared_thres = 0.3, p_val_thres = 0.05, TFBS_belief = 0.75, 
nonGE_belief = 0.5, woPKGE_belief = 0.5, gene_annot)
{
    colnames(omics$ge) <- gene_annot$entrezID[match(colnames(omics$ge),gene_annot$gene_symbol)]
    colnames(omics$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics$cnv),gene_annot$gene_symbol)])
    names(annot) <- gene_annot$entrezID[match(names(annot),gene_annot$gene_symbol)]
    
    if(TFBS_belief==woPKGE_belief)
    {
        message("GE-GE interactions without PK have the same belief as TFs-targets interactions. TFs-targets interactions will be considered as GE-GE interactions without prior knowledge.")
    }
  
    layers_def <- layers_def[order(layers_def$layer, decreasing = TRUE),]
    omics <- omics[layers_def$omics[order(layers_def$layer, 
        decreasing = TRUE)]]
  
    B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
        annot = annot, lm_METH = lm_METH, r_squared_thres = r_squared_thres, 
        p_val_thres = p_val_thres, TFtargs = TFtargs, 
        TFBS_belief = TFBS_belief, nonGE_belief = nonGE_belief, 
        woPKGE_belief = woPKGE_belief)
    pf_UB_res <- pf_UB_est(omics = B$omics, layers_def = layers_def, 
        B_prior_mat = B$B_prior_mat, annot = B$annot)
    return(list(pf_UB_BGe_pre = pf_UB_res, B_prior_mat = B$B_prior_mat, 
        annot = B$annot, omics = B$omics, layers_def = layers_def,
        omics_meth_original = B$omics_meth_original))
}

#' B_prior_mat
#' @description 
#' 'B_prior_mat' creates the biological prior matrix.
#'
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). 
#' Each component of the list is a matrix with samples in rows and features 
#' in columns.
#' @param PK data.frame with known interactions.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @param lm_METH logical asking whether to use linear regression to filter
#' methylation data (default=TRUE).
#' @param r_squared_thres numeric vector to define the R^2 used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.3).
#' @param p_val_thres numeric vector to define the p-value used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.05).
#' @param TFBS_belief numeric vector to define the belief concerning the TF and
#' its target interaction (default=0.75).
#' @param nonGE_belief numeric vector to define the belief concerning
#' interactions of features except GE-GE (default=0.5).
#' @param woPKGE_belief numeric vector to define the belief concerning GE-GE
#' interactions without prior knowledge (default=0.5).
#' @importFrom bestNormalize orderNorm
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", "gene_annot"),
#' package="IntOMICS")
#' colnames(omics$ge) <- gene_annot$entrezID[match(colnames(omics$ge),
#' gene_annot$gene_symbol)]
#' colnames(omics$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics$cnv),
#' gene_annot$gene_symbol)])
#' B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
#'      annot = annot, lm_METH = TRUE, r_squared_thres = 0.3,
#'      p_val_thres = 0.05, TFtargs = TFtarg_mat, TFBS_belief = 0.75, 
#'      nonGE_belief = 0.5,  woPKGE_belief = 0.5)
#'
#' @return List of 4 elements: prior biological matrix and data
#' preprocessing           
#' @export
B_prior_mat <- function(omics, PK, layers_def, TFtargs, annot, lm_METH, 
r_squared_thres, p_val_thres, TFBS_belief, nonGE_belief, woPKGE_belief)
{
    if(any(regexpr("EID:",colnames(omics[[layers_def$omics[1]]]))<0))
    {
        stop("Gene names in GE matrix are not in the correct form. Please, use EID:XXXX.")
    } # end if(regexpr("EID:",colnames(omics[[layers_def$omics[1]]]))<0)
  
    if(!is.null(TFtargs))
    {
        if(any(regexpr("EID:",colnames(TFtargs))<0))
        {
            message("Gene names in TFtargs are not in the correct form. Please, use EID:XXXX.")
        } # end if(any(regexpr("EID:",colnames(TFtargs))<0))
    } # end if(!is.null(TFtargs))
  
    if(is.null(PK))
    {
        message("There in no PK. Please, consider adding PK to increase the IntOMICS prediction accuracy.")
    } else {
        if(is.list(omics))
        {
            if(length(intersect(c(PK$src_entrez,PK$dest_entrez),
                unlist(mapply(colnames,omics))))==0)
            {
                message("There are no valid PK interactions. Please, check if PK fits omics features.")
            } # end if(length(intersect(c(PK$src_entrez,PK$dest_entrez)...
        } else {
            if(length(intersect(c(PK$src_entrez,
                PK$dest_entrez),colnames(omics)))==0)
            {
                message("There are no valid PK interactions. Please, check if PK fits omics features.")
            } # end if(length(intersect(c(PK$src_entrez,PK$dest_entrez),...
        } # end if else (is.list(omics))
    } # end if else if(is.null(PK))

    if(!all(sort(names(omics))==sort(layers_def$omics)))
    {
        stop("Names of the list 'omics' does not match omics names in the 'layers_def'.")
    } # end if(!all(order(names(omics))==order(layers_def$omics)))
  
    pk_belief <- c(1,0)
    features <- colnames(omics[[layers_def$omics[1]]])
    B_layer_max <- matrix(woPKGE_belief, ncol=length(features), 
        nrow=length(features), dimnames=list(features,features))
    TFs <- intersect(colnames(TFtargs),rownames(B_layer_max))
    targets <- intersect(rownames(TFtargs),rownames(B_layer_max))
    if(length(TFs)>0 & length(targets)>0)
    {
        TFtargs_spec <- matrix(TFtargs[targets, TFs], nrow = length(targets), 
            dimnames=list(targets,TFs))
        for(j in seq_len(ncol(TFtargs_spec)))
        {
            B_layer_max[colnames(TFtargs_spec)[j],
            names(which(TFtargs_spec[,j]==1))] <- TFBS_belief
        } # end for j
    } # end if(length(TFs)>0 & length(targets)>0)
    PK_present <- PK[PK$edge_type=="present",]
    PK_absent <- PK[PK$edge_type=="absent",]
  
    for(i in seq_len(nrow(B_layer_max)))
    {
        if (sum(PK_present$src_entrez==rownames(B_layer_max)[i])!=0)
        {
            B_layer_max[rownames(B_layer_max)[i], 
            intersect(PK_present[PK_present$src_entrez==rownames(
            B_layer_max)[i], "dest_entrez"],rownames(B_layer_max))] <-
            pk_belief[1]
        } # end if sum...
    
        if (sum(PK_absent$src_entrez==rownames(B_layer_max)[i]) != 0)
        {
            B_layer_max[rownames(B_layer_max)[i], 
            intersect(PK_absent[PK_absent$src_entrez==rownames(B_layer_max)[i],
            "dest_entrez"],rownames(B_layer_max))] <- pk_belief[2]
        } # end if sum...
    
    } # end for i...
    B <- B_layer_max
    diag(B) <- 0
    new_annot <- list()  
    omics_meth_original <- matrix(nrow = 0, ncol = 0)
    if(length(layers_def$omics)>1)
    {
        for(j in seq(from=2,to=length(layers_def$omics)))
        {
            features_lower <- colnames(omics[[layers_def$omics[j]]])
            if(any(regexpr("eid:",
                colnames(omics[[layers_def$omics[j]]]))>0)==TRUE)
            {
                B_layer_lower <- matrix(0, ncol=ncol(B),
                    nrow=length(features_lower),
                    dimnames=list(features_lower[match(colnames(B), 
                    toupper(features_lower), nomatch=0)],colnames(B)))
                diag_sim <- match(toupper(rownames(B_layer_lower)),
                    colnames(B_layer_lower))
                B_layer_lower[cbind(seq_len(length(features_lower)), 
                diag_sim)] <- nonGE_belief
                B_layer1 <- matrix(0, ncol=length(features_lower), 
                    nrow=nrow(B)+length(features_lower), 
                    dimnames=list(c(rownames(B), features_lower), 
                    rownames(B_layer_lower)))
                B <- cbind(rbind(B,B_layer_lower),B_layer1)
            } else {
                annot <- annot[intersect(names(annot),
                    colnames(omics[[layers_def$omics[1]]]))]
                annot <- lapply(annot, FUN=function(s) intersect(s,
                    colnames(omics[[layers_def$omics[j]]])))
                omics_meth_original <- omics[[layers_def$omics[j]]]
                omics[[layers_def$omics[j]]][,!apply(omics[[
                layers_def$omics[j]]],2,FUN=function(col) all(is.na(col)))]
                omics[[layers_def$omics[j]]] <-
                apply(omics[[layers_def$omics[j]]], 2, FUN=function(column)
                    bestNormalize::orderNorm(column)$x.t)
        
                if(lm_METH)
                {
                    new_annot <- lapply(seq_along(annot), function(list)
                        lm_meth(omics[[layers_def$omics[1]]], 
                        omics[[layers_def$omics[j]]], names(annot)[[list]], 
                        annot[[list]], r_squared_thres = r_squared_thres, 
                        p_val_thres = p_val_thres))
                    names(new_annot) <- names(annot)
                } else {
                    new_annot <- annot
                } # end if(lm_METH)
        
                if(sum(!mapply(is.null,new_annot))!=0)
                {
                    new_annot <- new_annot[!mapply(is.null,new_annot)]
                    features_lower <- unlist(new_annot)
                    B_layer_lower <- matrix(0, ncol=ncol(B),
                        nrow=length(features_lower),
                        dimnames=list(features_lower,colnames(B)))
                    for(a in seq_len(length(new_annot)))
                    {
                        B_layer_lower[intersect(features_lower,
                        new_annot[[a]]),names(new_annot)[a]] <- nonGE_belief
                    } # end for a
                    B_layer1 <- matrix(0, ncol=length(features_lower),
                        nrow=nrow(B)+length(features_lower), 
                        dimnames=list(c(rownames(B), features_lower),
                        rownames(B_layer_lower)))
                    B <- cbind(rbind(B,B_layer_lower),B_layer1)
                    if(length(unlist(new_annot)==1))
                    {
                        omics[[layers_def$omics[j]]] <-
                        as.matrix(omics[[layers_def$omics[j]]][,
                            unlist(new_annot)])
                        colnames(omics[[layers_def$omics[j]]]) <-
                            unlist(new_annot)
                    } else {
                        omics[[layers_def$omics[j]]] <-
                            omics[[layers_def$omics[j]]][,unlist(new_annot)]
                    } # end if else(length(unlist(new_annot)==1))
                } else {
                    new_annot <- list()
                    omics[[layers_def$omics[j]]] <- matrix(nrow = 0, ncol = 0)
                } # end if else (sum(!mapply(is.null,new_annot))==0)
            } # end if else (any(regexpr("eid:",...
        } # end for j
    } # end if(length(layers_def$omics)>1)
    return(list(B_prior_mat = B, annot = new_annot, omics = omics,
        omics_meth_original = omics_meth_original))
}

#' BGe score
#' @description
#' `DAGcorescore` The log of the BGe score simplified as much as possible. 
#' This function is from BiDAG package.
#' @param j character vector a node to be scored
#' @param parentnodes character vector the parents of the node j
#' @param n numeric vector number of nodes in the netwrok
#' @param param an object of class scoreparameters, which includes all
#' necessary information for calculating the BDe/BGe score
#'
#' @return Numeric vector of length 1            
#' @export
DAGcorescore <- function(j,parentnodes,n,param) {
  
    TN <- param$TN
    awpN <- param$awpN
    scoreconstvec <- param$scoreconstvec
    lp <- length(parentnodes) #number of parents
    awpNd2 <- (awpN-n+lp+1)/2
    A <- TN[j,j]
    switch(as.character(lp),
            "0"={# just a single term if no parents
            corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
            },
           
            "1"={# no need for matrices
            D <- TN[parentnodes,parentnodes]
            logdetD <- log(D)
            B <- TN[j,parentnodes]
            logdetpart2 <- log(A-B^2/D)
            corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
            },
           
            "2"={
            D <- TN[parentnodes,parentnodes]
             
            dettwobytwo <- function(D) {
            D[1,1]*D[2,2]-D[1,2]*D[2,1]
            }
             
            detD <- dettwobytwo(D)
            logdetD <- log(detD)
            B <- TN[j,parentnodes]
            logdetpart2 <- log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
            corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
            },
           
            {# otherwise we use cholesky decomposition to perform both
            D<-as.matrix(TN[parentnodes,parentnodes])
            choltemp<-chol(D)
            logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
            B<-TN[j,parentnodes]
            logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
            corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
            })
    return(corescore)
}

#' Node energy function
#' @description
#' `energy_function_node_specific`  For each node returns its energy over all
#' parent set configurations, the empty parent set is included.
#' @param all_parents_config matrix with all possible parent set configurations
#' (column indicates parents of given int_node).
#' @param B_prior_mat a biological prior matrix.
#' @param int_node character vector with given node name.
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", "gene_annot"),
#' package="IntOMICS")
#' colnames(omics$ge) <- gene_annot$entrezID[match(colnames(omics$ge),
#' gene_annot$gene_symbol)]
#' colnames(omics$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics$cnv),
#' gene_annot$gene_symbol)])
#' B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
#'    annot = annot, lm_METH = TRUE, r_squared_thres = 0.3, p_val_thres = 0.05,
#'    TFtargs = TFtarg_mat, TFBS_belief = 0.75, nonGE_belief = 0.5, 
#'    woPKGE_belief = 0.5)
#' all_parents_config <- matrix(c("EID:2535", "EID:2535", 
#' "EID:1857", "EID:2932"), byrow = TRUE, nrow=2)
#' energy_function_node_specific(all_parents_config = all_parents_config,
#' B_prior_mat = B$B_prior_mat, int_node = "EID:7482")
#'
#' @return Numeric vector of length 1            
#' @export
energy_function_node_specific <- function(all_parents_config, B_prior_mat,
int_node)
{
    epsilon <- apply(all_parents_config,2,FUN=function(z)
    if(is.na(z[1]))
    {
        sum(B_prior_mat[,int_node])
    } else {
        sum(1-B_prior_mat[z,int_node]) + 
            sum(B_prior_mat[-match(z,rownames(B_prior_mat)),int_node])
    })
    return(epsilon)
}

#' Linear regression GE~METH
#' @description
#' `lm_meth` The linear regression model for a dependent variable GE and
#' explanatory variable METH. Returns METH with significant coefficient, 
#' R^2 > threshold and R~Gaussian residuals.
#' @param ge_mat matrix of gene expression with samples in rows and 
#' features in columns.
#' @param meth_mat matrix of DNA methylaton with samples in rows and 
#' features in columns.
#' @param gene character vector with given node name.
#' @param meth_probes character vector methylation probes associated 
#' with a gene.
#' @param r_squared_thres numeric vector to define the R^2 used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.3).
#' @param p_val_thres numeric vector to define the p-value used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.05).
#' @importFrom stats lm shapiro.test
#'
#' @examples
#' data(list=c("annot", "omics"), package="IntOMICS")
#' lm_meth(ge_mat = omics$ge, meth_mat = omics$meth, 
#'     gene = "WNT2B", meth_probes = annot[["WNT2B"]], 
#'     r_squared_thres = 0.3, p_val_thres = 0.05)
#' 
#' @return Character vector with methylation probes           
#' @export
lm_meth <- function(ge_mat, meth_mat, gene, meth_probes, r_squared_thres,
p_val_thres)
{
    meth_probes_sig <- c()
    if(length(meth_probes)>0)
    {
        for(f in seq_len(length(meth_probes)))
        {
            res <- stats::lm(ge_mat[,gene] ~ meth_mat[,meth_probes[f]])
            if(nrow(summary(res)$coefficients)>1)
            {
                cond1 <- summary(res)$coefficients[2,"Pr(>|t|)"] < p_val_thres
                cond2 <- summary(res)$r.squared > r_squared_thres
                cond3 <- stats::shapiro.test(summary(res)$resid)$p.value > 0.1
                if(cond1 & cond2 & cond3)
                {
                    meth_probes_sig <- c(meth_probes_sig, meth_probes[f])
                } # end if(cond1 & cond2 & cond3)
            } # end if(nrow(summary(res)$coefficients)>1)
        } # end for f
    } # end if(length(meth_probes)>0)
  return(meth_probes_sig)
}

#' Partition function upper bound
#' @description
#' `pf_UB_est` Partition function upper bound estimation with beta = 0. 
#' For each node returns energy over all possible parent set configurations 
#' and BGe score.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param B_prior_mat a biological prior matrix.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @importFrom utils combn
#' @importFrom stats na.omit
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", "gene_annot"),
#' package="IntOMICS")
#' colnames(omics$ge) <- gene_annot$entrezID[match(colnames(omics$ge),
#' gene_annot$gene_symbol)]
#' colnames(omics$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics$cnv),
#' gene_annot$gene_symbol)])
#' B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
#'      annot = annot, lm_METH = TRUE, r_squared_thres = 0.3,
#'      p_val_thres = 0.05, TFtargs = TFtarg_mat, TFBS_belief = 0.75, 
#'      nonGE_belief = 0.5, woPKGE_belief = 0.5)
#' pf_UB_est(omics = B$omics, B_prior_mat = B$B_prior_mat, 
#' layers_def = layers_def, annot = B$annot)
#' 
#' @return List of 4 elements needed to simulate MCMC sampling            
#' @export
pf_UB_est <- function(omics, B_prior_mat, layers_def, annot)
{
    comb_all <- list()
    for(i in seq_len(ncol(omics[[layers_def$omics[1]]])))
    {
        int_node <- colnames(omics[[layers_def$omics[1]]])[i]
        potentials_layer <- intersect(rownames(B_prior_mat)
            [B_prior_mat[,int_node]>0],
            colnames(omics[[layers_def$omics[1]]]))
        comb_some <- list()
    
        for(rep in seq_len(layers_def$fan_in_ge[1]))
        {
            if(length(potentials_layer)>=rep)
            {
                comb_some[[rep]] <- utils::combn(potentials_layer,rep)
            } # end if(length(potentials_layer)>=rep)
        } # end for rep
        comb_some[[length(comb_some)+1]] <- matrix(NA,1,1)
        if(length(layers_def$omics)>1)
        {
            modalities <- layers_def$omics[-1]
            if(any(mapply(FUN=function(mod)
                any(regexpr("eid:",colnames(mod))>0)==TRUE, 
                omics[modalities])) & tolower(int_node) %in%
                unlist(lapply(omics[modalities],colnames)))
            {
                comb_some[seq(length(comb_some)+1,
                length.out = length(comb_some))] <-
                lapply(comb_some,FUN=function(list)
                rbind(list,tolower(int_node)))
            } # end if(any(mapply(FUN=function(mod)...
      
            if(any(mapply(omics, FUN=function(list)
            any(regexpr("eid:",colnames(list), 
            ignore.case = TRUE)<0))) & length(annot[[int_node]])>0)
            {
                modality <- names(which(mapply(omics, FUN=function(list)
                    any(regexpr("eid:",colnames(list), 
                    ignore.case = TRUE)<0))==TRUE))
                max_fan_in <- max(layers_def$fan_in_ge[layers_def$omics==
                    modality],length(annot[[int_node]]), na.rm = TRUE)
                comb_some_meth <- list()
                for(rep in seq_len(max_fan_in))
                {
                    if(length(annot[[int_node]])>=rep)
                    {
                        comb_some_meth[[rep]] <- combn(annot[[int_node]],rep)
                    } # end if(length(annot[[int_node]]<=rep)
                } # end for rep
                comb_some_meth_add <- list()
                for(meth_pr in seq_len(length(comb_some_meth)))
                {
                    comb_some_meth_add_2 <- list()
                    for(a in seq_len(ncol(comb_some_meth[[meth_pr]])))
                    {
                        comb_some_meth_add_2[[a]] <- lapply(comb_some, 
                            FUN=function(par_def) apply(par_def,2,
                            FUN=function(par_def_col) c(par_def_col,
                            comb_some_meth[[meth_pr]][,a])))
                    } # end for a
                    comb_some_meth_add <- c(comb_some_meth_add,
                        comb_some_meth_add_2)
                } # end for meth_pr
                comb_some_meth_add <- unlist(comb_some_meth_add, 
                    recursive = FALSE)
                comb_some <- c(comb_some, comb_some_meth_add)
            } # if if(any(mapply(omics,FUN=function(list)...
        } # end if(length(layers_def$omics)>1)
    
        comb_some <- lapply(comb_some,stats::na.omit)
        comb_some[[1]] <- cbind(comb_some[[1]], NA)
        comb_some <- comb_some[mapply(comb_some,FUN=function(x) nrow(x))!=0]
        parents_config <- list()
        for(l in seq_len(max(mapply(comb_some,FUN=function(x) nrow(x)))))
        {
            parents_config[[l]] <- do.call(cbind, comb_some[mapply(comb_some,
                FUN=function(x) nrow(x))==l])
        } # end for l
        comb_all[[i]] <- parents_config
    } # end for i
    names(comb_all) <- colnames(omics[[layers_def$omics[1]]])
    if(length(layers_def$omics)>1)
    {
        comb_all_others <- vector(mode = "list", length = sum(mapply(ncol,
            omics[setdiff(layers_def$omics,layers_def$omics[1])])))
        comb_all_others <- lapply(comb_all_others, FUN=function(list) list <-
            matrix(NA))
        names(comb_all_others) <- unlist(mapply(colnames,
            omics[setdiff(layers_def$omics,layers_def$omics[1])]))
        comb_all <- c(comb_all, comb_all_others)
    } # end if(length(layers_def$omics)>1)
    energy_all_configs_node <- list()
    for(i in seq_len(ncol(omics[[layers_def$omics[1]]])))
    {
        energy_all_configs_node[[i]] <- unlist(lapply(comb_all[[i]],
            FUN=function(list) energy_function_node_specific(list, B_prior_mat,
            names(comb_all)[i])))
    } # end for i
    if(length(layers_def$omics)>1)
    {
        for(i in c((ncol(omics[[layers_def$omics[1]]])+1) :
        (ncol(omics[[layers_def$omics[1]]]) +
        sum(mapply(ncol,omics[setdiff(layers_def$omics,
        layers_def$omics[1])])))))
        {
            energy_all_configs_node[[i]] <- 
            sum(B_prior_mat[,names(comb_all)[i]])
        } # end for i
    } # end if(length(layers_def$omics)>1)
    partition_func_UB <- sum(log(mapply(energy_all_configs_node,FUN=function(x)
        sum(exp(-0*x)))))
  
    data <- do.call(cbind, omics[mapply(nrow,omics)>0])
    myScore <- scoreparameters_BiDAG_BGe(n = ncol(data), data = data)
    n <- ncol(myScore$data)
    BGe_score_list <- list()
    for(i in seq_len(ncol(omics[[layers_def$omics[1]]])))
    {
        BGe_score_list[[i]] <- lapply(comb_all[[i]], FUN=function(list)
            apply(list, 2, FUN=function(column) 
        if(is.na(column[1]))
        {
            DAGcorescore(names(comb_all)[i], integer(length = 0), 
                n = myScore$n, param = myScore)
        } else {
            DAGcorescore(names(comb_all)[i], column, n = myScore$n, 
                param = myScore)
        }))
    } # end for(i in seq_len(ncol(omics[[layers_def$omics[1]]])))
    if(length(layers_def$omics)>1)
    {
        for(i in c((ncol(omics[[layers_def$omics[1]]])+1) :
        (ncol(omics[[layers_def$omics[1]]]) +
        sum(mapply(ncol,omics[setdiff(layers_def$omics,
        layers_def$omics[1])])))))
        {
            BGe_score_list[[i]] <- matrix(DAGcorescore(names(comb_all)[i], 
                integer(length = 0), n = myScore$n, param = myScore))
        } # end for i
    } # end if(length(layers_def$omics)>1)
    names(BGe_score_list) <- names(comb_all)
    return(list(partition_func_UB = partition_func_UB, 
        parents_set_combinations = comb_all, 
        energy_all_configs_node = energy_all_configs_node, 
        BGe_score_all_configs_node = BGe_score_list))
}

#' Range between 0 and 1
#' @description
#' `range_01` This function re-scales a numeric vector so that it ranges
#' between 0 and 1.
#' @param x numeric vector.
#'
#' @examples
#' range_01(stats::rnorm(10))
#' 
#' @return Numeric vector with normalised values           
#' @export
range_01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' BGe score parameters
#' @description
#' `scoreparameters_BiDAG_BGe` Returns parameters needed for calculation 
#' of the BGe score. This function is from BiDAG package.
#' @param data matrix with features in columns and a number of rows equal 
#' to the number of samples.
#' @param n numeric number of columns in data matrix.
#' @param bgepar list which contains parameters for BGe score computation.
#' @importFrom stats cov
#' @return Object of class scoreparameters, which includes all necessary
#' information for calculating the BDe/BGe score
#' @export
scoreparameters_BiDAG_BGe <- function (n, data, 
bgepar = list(am = 1, aw = NULL))
{
    if (anyNA(data)) {
        message("Dataset contains missing data (covariance matrix computation: complete.obs parameter - missing values are handled by casewise deletion)")
    }

    if (all(is.character(colnames(data)))) {
        nodeslabels <- colnames(data)
    } else {
        nodeslabels <- unlist(lapply(seq_len(n), 
            function(x) paste("v", x, sep = "")))
    }
    colnames(data) <- nodeslabels
    initparam <- list()
    initparam$labels <- nodeslabels
    initparam$type <- "bge"
    initparam$DBN <- FALSE
    initparam$weightvector <- NULL
    initparam$data <- data

    initparam$bgnodes <- NULL
    initparam$static <- NULL
    initparam$mainnodes <- seq_len(n)
  
    initparam$bgn <- 0
    initparam$n <- n
    initparam$nsmall <- n
    initparam$labels.short <- initparam$labels
    initparam$logedgepmat <- NULL

    N <- nrow(data)
    if(N==1)
    {
        covmat <- matrix(0,nrow=n, ncol=n,
            dimnames=list(initparam$labels.short,initparam$labels.short))
    } else {
        covmat <- stats::cov(data, use = "complete.obs") * (N - 1)
    } # end if else (N==1)
    means <- colMeans(data, na.rm = TRUE)
    bgepar$aw <- n + bgepar$am + 1
  
    initparam$am <- bgepar$am
    initparam$aw <- bgepar$aw
    initparam$N <- N
    initparam$means <- means
    mu0 <- numeric(n)
    T0scale <- bgepar$am * (bgepar$aw - n - 1)/(bgepar$am + 1)
    T0 <- diag(T0scale, n, n)
    initparam$TN <- T0 + covmat + ((bgepar$am * N)/(bgepar$am + N)) * 
        (mu0 - means) %*% t(mu0 - means)
    initparam$awpN <- bgepar$aw + N
    constscorefact <- -(N/2) * log(pi) + (1/2) * log(bgepar$am/(bgepar$am +  N))
    initparam$muN <- (N * means + bgepar$am * mu0)/(N + bgepar$am)
    initparam$SigmaN <- initparam$TN/(initparam$awpN - n - 1)
    initparam$scoreconstvec <- numeric(n)
    for (j in seq_len(n)) {
        awp <- bgepar$aw - n + j
        initparam$scoreconstvec[j] <- constscorefact - lgamma(awp/2) + 
        lgamma((awp + N)/2) + ((awp + j - 1)/2) * log(T0scale)
    }
    attr(initparam, "class") <- "scoreparameters"
    initparam
}
