#' BGe score
#' @description
#' `DAGcorescore` The log of the BGe score simplified as much as possible. 
#' This function is from BiDAG package.
#' @param j character vector a node to be scored
#' @param parentnodes character vector the parents of the node j
#' @param n numeric vector number of nodes in the netwrok
#' @param param an object of class scoreparameters, which includes all
#' necessary information for calculating the BDe/BGe score
#'
#' @return Numeric vector of length 1            
#' @export
DAGcorescore <- function(j,parentnodes,n,param) {
  
    TN <- param$TN
    awpN <- param$awpN
    scoreconstvec <- param$scoreconstvec
    lp <- length(parentnodes) #number of parents
    awpNd2 <- (awpN-n+lp+1)/2
    A <- TN[j,j]
    switch(as.character(lp),
            "0"={# just a single term if no parents
            corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
            },
           
            "1"={# no need for matrices
            D <- TN[parentnodes,parentnodes]
            logdetD <- log(D)
            B <- TN[j,parentnodes]
            logdetpart2 <- log(A-B^2/D)
            corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
            },
           
            "2"={
            D <- TN[parentnodes,parentnodes]
             
            dettwobytwo <- function(D) {
            D[1,1]*D[2,2]-D[1,2]*D[2,1]
            }
             
            detD <- dettwobytwo(D)
            logdetD <- log(detD)
            B <- TN[j,parentnodes]
            logdetpart2 <- log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
            corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
            },
           
            {# otherwise we use cholesky decomposition to perform both
            D<-as.matrix(TN[parentnodes,parentnodes])
            choltemp<-chol(D)
            logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
            B<-TN[j,parentnodes]
            logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
            corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
            })
    return(corescore)
}
