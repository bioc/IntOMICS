#' Acyclic network identification.
#' @description
#' `is_acyclic` This function is from bnstruct R package. Check if the directed
#' graph is acyclic.
#' @param g adajcency matrix of given network/graph.
#'
#' @examples
#' is_acyclic(matrix(c(1,rep(0,20)), nrow=3))
#'
#' @return boolean of length 1
#' @keywords internal
#' @export
is_acyclic <- function(g)
{
    rem <- rep(FALSE,nrow(g))
    while( !all(rem) ) # still some edges to remove
    {
        leaves <- (rowSums(g) == 0)
        if( !any(leaves & !rem) )
            return(FALSE)
        g[,leaves] <- 0L
        rem <- rem | leaves
    }
  return(TRUE)
}
