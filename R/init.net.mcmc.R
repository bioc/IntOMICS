#' Random initial network
#' @description
#' `init.net.mcmc` This function is used to sample random initial network. 
#' The edges are sampled only between GE nodes.
#' @param omics named list containing the gene expression (possibly copy number
#' variation and methylation data). Each component of the list is a matrix 
#' with samples in rows and features in columns.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @param B_prior_mat a biological prior matrix.
#'
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", "gene_annot"),
#' package="IntOMICS")
#' colnames(omics$ge) <- gene_annot$entrezID[match(colnames(omics$ge),
#' gene_annot$gene_symbol)]
#' colnames(omics$cnv) <- tolower(gene_annot$entrezID[match(colnames(omics$cnv),
#' gene_annot$gene_symbol)])
#' B <- B_prior_mat(omics = omics, PK = PK, annot = annot, lm_METH = TRUE,
#'      layers_def = layers_def, r_squared_thres = 0.3, p_val_thres = 0.05,
#'      TFtargs = TFtarg_mat, TFBS_belief = 0.75, nonGE_belief = 0.5, 
#'      woPKGE_belief = 0.5)
#' init.net.mcmc(omics = B$omics, layers_def = layers_def, 
#'      B_prior_mat = B$B_prior_mat)
#'
#' @return List of 2 elements: random adjacency network and empty network
init.net.mcmc <- function(omics, layers_def, B_prior_mat)
{
    empty.net <- matrix(0, nrow = sum(mapply(ncol,omics)), 
        ncol = sum(mapply(ncol,omics)),
        dimnames = list(unlist(mapply(colnames,omics)),
        unlist(mapply(colnames,omics))))
    init.net <- sample.chain(empty_net = empty.net, 
        omics_ge = omics[[layers_def$omics[1]]])
    rownames(init.net@dag) <- rownames(empty.net)
    colnames(init.net@dag) <- rownames(empty.net)
    init.net@dag <- init.net@dag[rownames(B_prior_mat),rownames(B_prior_mat)]
    if(any(init.net@dag==1 & B_prior_mat==0))
    {
        init.net@dag[init.net@dag==1 & B_prior_mat==0] <- 0
    }
    while(!is.acyclic(init.net@dag) |
    any(colSums(init.net@dag[colnames(omics[[layers_def$omics[1]]]),
    colnames(omics[[layers_def$omics[1]]])]) > layers_def$fan_in_ge[1]))
    {
        init.net <- sample.chain(empty_net = empty.net, 
            omics_ge = omics[[layers_def$omics[1]]])
        rownames(init.net@dag) <- rownames(empty.net)
        colnames(init.net@dag) <- rownames(empty.net)
        init.net@dag <- init.net@dag[rownames(B_prior_mat),
            rownames(B_prior_mat)]
        if(any(init.net@dag==1 & B_prior_mat==0))
        {
            init.net@dag[init.net@dag==1 & B_prior_mat==0] <- 0
        }
    }
    source.net <- list(adjacency = init.net@dag, nbhd.size = c(), 
        proposal.distr = c(), energy = c(), prior = c(), BGe = c(), 
        likelihood_part = c(), likelihood = c(), acceptance = c(), 
        edge_move = c())
    return(list(source.net = source.net, empty.net = empty.net))
}
