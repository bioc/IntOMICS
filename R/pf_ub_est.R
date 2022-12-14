#' Partition function upper bound
#' @description
#' `pf_ub_est` Partition function upper bound estimation with beta = 0. 
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
#' @return List of 4 elements needed to simulate MCMC sampling         
pf_ub_est <- function(omics, B_prior_mat, layers_def, annot)
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
                comb_some[[rep]] <- combn(potentials_layer,rep)
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
    
        comb_some <- lapply(comb_some,na.omit)
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
    myScore <- score_parameters_bidag_bge(n = ncol(data), data = data)
    n <- ncol(myScore$data)
    BGe_score_list <- list()
    for(i in seq_len(ncol(omics[[layers_def$omics[1]]])))
    {
        BGe_score_list[[i]] <- lapply(comb_all[[i]], FUN=function(list)
            apply(list, 2, FUN=function(column) 
        if(is.na(column[1]))
        {
            dag_core_score(names(comb_all)[i], integer(length = 0), 
                n = myScore$n, param = myScore)
        } else {
            dag_core_score(names(comb_all)[i], column, n = myScore$n, 
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
            BGe_score_list[[i]] <- matrix(dag_core_score(names(comb_all)[i], 
                integer(length = 0), n = myScore$n, param = myScore))
        } # end for i
    } # end if(length(layers_def$omics)>1)
    names(BGe_score_list) <- names(comb_all)
    return(list(partition_func_UB = partition_func_UB, 
        parents_set_combinations = comb_all, 
        energy_all_configs_node = energy_all_configs_node, 
        BGe_score_all_configs_node = BGe_score_list))
}
