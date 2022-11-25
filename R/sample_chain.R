#' Random initial network edge generation
#' @description
#' `sample_chain` This function is used to sample random initial network. 
#' The edges are sampled only between GE nodes.
#' @param empty_net adjacency matrix of an empty network/graph 
#' (all values are 0).
#' @param omics_ge matrix with gene expression data (samples in rows and
#' features in columns).
#' @importFrom bnstruct BNDataset
#' @importFrom bnstruct BN
#' @importFrom bnstruct dag
#' @importFrom bnstruct learn.params
#' @return BN object with conditional probabilities
#' @keywords internal
#' @export 
sample_chain <- function(empty_net, omics_ge)
{
    dataset_BND <- BNDataset(data = empty_net, 
        discreteness = rep('d',ncol(empty_net)),
        variables = c(colnames(empty_net)), node.sizes = rep(2,ncol(empty_net)),
        starts.from=0)
    net <- BN(dataset_BND)
    net.dag <- bnstruct::dag(net)
    n <- ncol(omics_ge)
    chain <- sample(n,n)
    for(i in seq(from=2, to=n))
    {
        net.dag[chain[i-1],chain[i]] <- 1
    }
    bnstruct::dag(net) <- net.dag
    return(learn.params(net,dataset_BND))
}
