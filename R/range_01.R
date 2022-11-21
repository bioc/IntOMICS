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
range_01 <- function(x){(x-min(x))/(max(x)-min(x))}
