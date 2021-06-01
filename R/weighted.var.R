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
#' @export
weighted.var <- function(x,w){
  return(weighted.cov(x=x,y=x,w=w))
}