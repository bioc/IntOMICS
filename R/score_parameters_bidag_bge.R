#' BGe score parameters
#' @description
#' `score_parameters_bidag_bge` Returns parameters needed for calculation 
#' of the BGe score. This function is from BiDAG package.
#' @param data matrix with features in columns and a number of rows equal 
#' to the number of samples.
#' @param n numeric number of columns in data matrix.
#' @param bgepar list which contains parameters for BGe score computation.
#' @importFrom stats cov
#' @return Object of class scoreparameters, which includes all necessary
#' information for calculating the BDe/BGe score
#' @keywords internal
#' @export
score_parameters_bidag_bge <- function (n, data, 
bgepar = list(am = 1, aw = NULL))
{
    if (anyNA(data)) {
        message("Dataset contains missing data (covariance matrix computation: complete.obs parameter - missing values are handled by casewise deletion)")
    }

    if (all(is.character(colnames(data)))) {
        nodeslabels <- colnames(data)
    } else {
        nodeslabels <- unlist(lapply(seq_len(n), 
            function(x) paste("v", x, sep = "")))
    }
    colnames(data) <- nodeslabels
    initparam <- list()
    initparam$labels <- nodeslabels
    initparam$type <- "bge"
    initparam$DBN <- FALSE
    initparam$weightvector <- NULL
    initparam$data <- data

    initparam$bgnodes <- NULL
    initparam$static <- NULL
    initparam$mainnodes <- seq_len(n)
  
    initparam$bgn <- 0
    initparam$n <- n
    initparam$nsmall <- n
    initparam$labels.short <- initparam$labels
    initparam$logedgepmat <- NULL

    N <- nrow(data)
    if(N==1)
    {
        covmat <- matrix(0,nrow=n, ncol=n,
            dimnames=list(initparam$labels.short,initparam$labels.short))
    } else {
        covmat <- cov(data, use = "complete.obs") * (N - 1)
    } # end if else (N==1)
    means <- colMeans(data, na.rm = TRUE)
    bgepar$aw <- n + bgepar$am + 1
  
    initparam$am <- bgepar$am
    initparam$aw <- bgepar$aw
    initparam$N <- N
    initparam$means <- means
    mu0 <- numeric(n)
    T0scale <- bgepar$am * (bgepar$aw - n - 1)/(bgepar$am + 1)
    T0 <- diag(T0scale, n, n)
    initparam$TN <- T0 + covmat + ((bgepar$am * N)/(bgepar$am + N)) * 
        (mu0 - means) %*% t(mu0 - means)
    initparam$awpN <- bgepar$aw + N
    constscorefact <- -(N/2) * log(pi) + (1/2) * log(bgepar$am/(bgepar$am +  N))
    initparam$muN <- (N * means + bgepar$am * mu0)/(N + bgepar$am)
    initparam$SigmaN <- initparam$TN/(initparam$awpN - n - 1)
    initparam$scoreconstvec <- numeric(n)
    for (j in seq_len(n)) {
        awp <- bgepar$aw - n + j
        initparam$scoreconstvec[j] <- constscorefact - lgamma(awp/2) + 
        lgamma((awp + N)/2) + ((awp + j - 1)/2) * log(T0scale)
    }
    attr(initparam, "class") <- "scoreparameters"
    initparam
}
