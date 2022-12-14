#' Epsilon  
#' @description
#' `epsilon` This function returns the epsilon value for each variable/node 
#' of the network. 
#' The sum of the epsilons of all variables/nodes in the network gives us 
#' the energy of given network.
#' @param net adjacency matrix of given network.
#' @param B_prior_mat a biological prior matrix.
#' @return Numeric vector of length 1: epsilon of given adjacency matrix
#' (needed to compute energy of given adjacency matrix)
epsilon <- function(net, B_prior_mat)
{
    epsilon <- rep(NA,nrow(net))
    for(i in seq_len(nrow(net)))
    {
        iter_feature <- rownames(net)[i]
        parent <- names(net[,iter_feature])[net[,iter_feature]==1]
        noparent <- setdiff(rownames(net),parent)
        epsilon[i] <- sum(1-B_prior_mat[parent,iter_feature]) +
            sum(B_prior_mat[noparent,iter_feature])
    } # end for i
    return(epsilon)
}
