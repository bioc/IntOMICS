#' Arrow of directed edges tuning
#' @description
#' `normalise` This function is from the ambient package. 
#' It is used to determine the position of directed edge arrows.
#' @param x numeric vector to be modified.
#' @param from numeric vector range of x.
#' @param to numeric vector range of normalised x.
#' @importFrom stats median
#'
#' @examples
#' x <- seq(1,10)
#' normalise(x, from = range(x), to = c(0, 1))
#'
#' @return Numeric vector
normalise <- function (x, from = range(x), to = c(0, 1)) {
    x <- (x - from[1])/(from[2] - from[1])
    if (!identical(to, c(0, 1))) {
        x <- x * (to[2] - to[1]) + to[1]
    }
    x
}
