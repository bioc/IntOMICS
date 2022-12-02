#' biological prior matrix
#' @description 
#' 'b_prior_mat' creates the biological prior matrix.
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
#' omics <- omics_to_list(omics = omics, gene_annot = gene_annot, 
#'                        layers_def = layers_def)
#' B <- b_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
#'      annot = annot, lm_METH = TRUE, r_squared_thres = 0.3,
#'      p_val_thres = 0.05, TFtargs = TFtarg_mat, TFBS_belief = 0.75, 
#'      nonGE_belief = 0.5,  woPKGE_belief = 0.5)
#'
#' @return List of 4 elements: prior biological matrix and data
#' preprocessing  
#' @keywords internal
#' @export   
b_prior_mat <- function(omics, PK, layers_def, TFtargs, annot, lm_METH, 
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
        message('Please, consider adding PK to increase the IntOMICS prediction accuracy.')
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
                    orderNorm(column)$x.t)
        
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
