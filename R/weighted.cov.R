#' @title weighted.cov
#' 
#' @description This function is to compute an unbiased sample weighted covariance. The function uses only pairwise complete observations.
#'
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#'
#' @param x An input vector to compute a covariance, cov(x,y)
#' @param y An input vector to compute a covariance, cov(x,y)
#' @param w A weight vector
#' 
#' @return Returns an unbiased sample weighted covariance
#'
#' @export
weighted.cov <- function(x,y,w){
  
  if (is.null(w)){
    res <- cov(x,y,use="pairwise.complete.obs")
  } else{
    
    ind <- !is.na(x) & !is.nan(x) & !is.na(y) & !is.nan(y) & !is.na(w) & !is.nan(w)
    x <- x[ind]
    y <- y[ind]
    w <- w[ind]

    v1 <- sum(w)
    v2 <- sum(w^2)
    
    res <- sum(w*(x-weighted.mean(x=x,w=w))*(y-weighted.mean(x=y,w=w))) / (v1-v2/v1)
  }
  
  return(res)
}