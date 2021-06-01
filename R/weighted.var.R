#' @title weighted.var
#' 
#' @description This function is to compute an unbiased sample weighted variance.
#'
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#'
#' @param x A vector to compute a variance, var(x)
#' @param w A weight vector
#' 
#' @return Returns an unbiased sample weighted variance
#' 
#' @examples 
#' ## If you do not specify weights, 
#' ## it returns the usual unweighted sample variance
#' weighted.var(x=c(1,3,5)) 
#' 
#' weighted.var(x=c(1,3,5),w=c(0.1,0.5,0.4))
#'
#' @export
weighted.var <- function(x,w=NULL){
  return(weighted.cov(x=x,y=x,w=w))
}