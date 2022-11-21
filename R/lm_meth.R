#' Linear regression GE~METH
#' @description
#' `lm_meth` The linear regression model for a dependent variable GE and
#' explanatory variable METH. Returns METH with significant coefficient, 
#' R^2 > threshold and R~Gaussian residuals.
#' @param ge_mat matrix of gene expression with samples in rows and 
#' features in columns.
#' @param meth_mat matrix of DNA methylaton with samples in rows and 
#' features in columns.
#' @param gene character vector with given node name.
#' @param meth_probes character vector methylation probes associated 
#' with a gene.
#' @param r_squared_thres numeric vector to define the R^2 used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.3).
#' @param p_val_thres numeric vector to define the p-value used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.05).
#' @importFrom stats lm shapiro.test
#'
#' @examples
#' data(list=c("annot", "omics"), package="IntOMICS")
#' lm_meth(ge_mat = omics$ge, meth_mat = omics$meth, 
#'     gene = "WNT2B", meth_probes = annot[["WNT2B"]], 
#'     r_squared_thres = 0.3, p_val_thres = 0.05)
#' 
#' @return Character vector with methylation probes           
lm_meth <- function(ge_mat, meth_mat, gene, meth_probes, r_squared_thres,
p_val_thres)
{
    meth_probes_sig <- c()
    if(length(meth_probes)>0)
    {
        for(f in seq_len(length(meth_probes)))
        {
            res <- lm(ge_mat[,gene] ~ meth_mat[,meth_probes[f]])
            if(nrow(summary(res)$coefficients)>1)
            {
                cond1 <- summary(res)$coefficients[2,"Pr(>|t|)"] < p_val_thres
                cond2 <- summary(res)$r.squared > r_squared_thres
                cond3 <- shapiro.test(summary(res)$resid)$p.value > 0.1
                if(cond1 & cond2 & cond3)
                {
                    meth_probes_sig <- c(meth_probes_sig, meth_probes[f])
                } # end if(cond1 & cond2 & cond3)
            } # end if(nrow(summary(res)$coefficients)>1)
        } # end for f
    } # end if(length(meth_probes)>0)
  return(meth_probes_sig)
}
