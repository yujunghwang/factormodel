#' @title cproxyme
#' @description This function estimates a linear factor model using continuous variables.
#' The linear factor model to estimate has the following form.
#' proxy = intercept + factorloading * (latent variable) + measurement error
#' The measurement error is assumed to follow a Normal distribution with a mean zero and a variance, which needs to be estimated.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references \describe{
#' \item{Cunha, F., Heckman, J. J., & Schennach, S. M. (2010)}{Estimating the technology of cognitive and noncognitive skill formation. Econometrica, 78(3), 883-931. <doi:10.3982/ECTA6551:>}
#' \item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper.}}
#' @importFrom utils install.packages
#' @importFrom gtools combinations
#' @import     stats
#'
#' @param dat A proxy variable data frame list.
#'
#' @param anchor This is a column index of an anchoring proxy variable. Default is 1. That is, the code will use the first column in dat data frame as an achoring variable.
#'
#' @param weights An optional weight vector
#'
#' @return Returns a list of 3 components : \describe{
#' \item{alpha0}{This is a vector of intercepts in a linear factor model. The k-th entry is the intercept of k-th proxy variable factor model.}
#'
#' \item{alpha1}{This is a vector of factor loadings. The k-th entry is the factor loading of k-th proxy variable. The factor loading of anchoring variable is normalized to 1.}
#'
#' \item{varnu}{This is a vector of variances of measurement errors in proxy variables. The k-th entry is the variance of k-th proxy measurement error. The measurement error is assumed to follow a Normal distribution with mean 0.}
#'
#' \item{mtheta}{This is a mean of the latent variable. It is equal to the mean of the anchoring proxy variable.}
#'
#' \item{vartheta}{This is a variance of the latent variable.}}
#'
#' @examples
#' dat1 <- data.frame(proxy1=c(1,2,3),proxy2=c(0.1,0.3,0.6),proxy3=c(2,3,5))
#' cproxyme(dat=dat1,anchor=1)
#' ## you can specify weights
#' cproxyme(dat=dat1,anchor=1,weights=c(0.1,0.5,0.4))
#'
#' @export
cproxyme <- function(dat,anchor=1,weights=NULL){

  # load libraries
  requireNamespace("stats")
  requireNamespace("utils")
  requireNamespace("gtools")

  #############
  # check if inputs are there in a correct form
  #############

  if (!is.matrix(dat) & !is.data.frame(dat)){
    stop("please provide proxy data in either matrix or data frame format.")
  }

  # check if the weights vector has a correct length
  if (!is.null(weights)){
    if (length(weights)!=dim(dat)[1]){
      stop("Incorrect length for the weights vector. The length must be equal to the number of rows of 'dat'.")
    }
  }

  # check if any weights vector includes NA or NaN or Inf
  if (sum(is.na(weights))>0|sum(is.nan(weights))>0|sum(is.infinite(weights))>0){
    stop("The weights vector can not include any NAs or NaNs or Infs.")
  }

  #############
  # prepare data in a right form
  #############

  # set the data frame size
  N  <- dim(dat)[1] # num of obs
  np <- dim(dat)[2] # num of proxy vars

  # vectors to save results
  alpha0 <- rep(NA,np)
  alpha1 <- rep(NA,np)
  varnu  <- rep(NA,np)

  # mean of latent variable
  if (is.null(weights)){
    mtheta <- mean(dat[,anchor],na.rm=TRUE)
  } else{
    mtheta <- weighted.mean(x=dat[,anchor],w=weights,na.rm=TRUE)
  }

  # anchoring variable coefficients
  alpha1[anchor] <- 1
  alpha0[anchor] <- 0

  # variance of latent variable
  # use all combos excluding the anchoring variable
  other <- c(1:np)[-anchor]
  combo <- combinations((np-1),2,v=other)
  nc <- dim(combo)[1]

  temp <-0
  for (k in 1:nc){
    temp <- temp + (weighted.cov(x=dat[,anchor],y=dat[,combo[k,1]],w=weights)*weighted.cov(x=dat[,anchor],y=dat[,combo[k,2]],w=weights)) / weighted.cov(x=dat[,combo[k,1]],y=dat[,combo[k,2]],w=weights)
  }
  vartheta <- temp/nc
  rm(temp)

  # factor loading
  for (k in 1:np){
    if (k!=anchor){
      alpha1[k] <- weighted.cov(x=dat[,anchor],y=dat[,k],w=weights)/vartheta
    }
  }

  # variance of measurement error
  for (k in 1:np){
    varnu[k] <- weighted.var(x=dat[,k],w=weights) - (alpha1[k]^2) * vartheta
  }

  # intercept
  for (k in 1:np){
    if (k!=anchor){

      if (is.null(weights)){
        alpha0[k] <- mean(dat[,k],na.rm=TRUE) - alpha1[k] * mtheta
      } else{
        alpha0[k] <- weighted.mean(x=dat[,k],w=weights,na.rm=TRUE) - alpha1[k] * mtheta
      }

    }
  }

  return(list(alpha0=alpha0,alpha1=alpha1,varnu=varnu,mtheta=mtheta,vartheta=vartheta))
}
