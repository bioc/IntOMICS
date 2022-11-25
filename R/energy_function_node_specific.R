#' Node energy function
#' @description
#' `energy_function_node_specific`  For each node returns its energy over all
#' parent set configurations, the empty parent set is included.
#' @param all_parents_config matrix with all possible parent set configurations
#' (column indicates parents of given int_node).
#' @param B_prior_mat a biological prior matrix.
#' @param int_node character vector with given node name.
#' @return Numeric vector of length 1    
#' @keywords internal
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
