#' Number of reverse edge candidates
#' @description
#' `fan_in_reverse` Determine the number of edges that can be reversed using
#' the fan-in restriction in the largest layer.
#' @param positions character vector indicating the interaction between two
#' nodes (the first string indicates the source node, the second string
#' indicates the target node).
#' @param net_layer_max adjacency matrix of the network containing only GE
#' nodes.
#' @param layers_def data.frame containing the modality ID, corresponding layer
#' in BN and maximal number of parents from given layer to GE nodes.
#' @return Numeric vector of length 1: reverse edge candidates
#' @keywords internal
#' @export 
fan_in_reverse <- function(positions, net_layer_max, layers_def)
{
    net_layer_max[positions["col"],positions["row"]] <- 1
    possible_rev_edges <- sum(net_layer_max[,positions["col"]]) <=
        layers_def$fan_in_ge[1]
    return(possible_rev_edges)
}
