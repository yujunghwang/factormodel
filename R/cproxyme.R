#' @title cproxyme
#' @description This function estimates a linear factor model using continuous variables.
#' The linear factor model to estimate has the following form.
#' proxy = intercept + factorloading * (latent variable) + measurement error
#' The measurement error is assumed to follow a Normal distribution with a mean zero and a variance, which needs to be estimated.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references Cunha, F., Heckman, J. J., & Schennach, S. M. (2010). Estimating the technology of cognitive and noncognitive skill formation. Econometrica, 78(3), 883-931.
#'             Hwang, Yujung (2021). Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper.
#' @importFrom utils install.packages
#' @importFrom gtools combinations
#' @import     stats
#'
#' @param dat A proxy variable data frame list.
#'
#' @param anchor This is a column index of an anchoring proxy variable. Default is 1. That is, the code will use the first column in dat data frame as an achoring variable.
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
#' @export
cproxyme <- function(dat,anchor=1){

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
  mtheta <- mean(dat[,anchor],na.rm=TRUE)

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
    temp <- temp + (cov(dat[,anchor],dat[,combo[k,1]],use="pairwise.complete.obs")*cov(dat[,anchor],dat[,combo[k,2]],use="pairwise.complete.obs")) / cov(dat[,combo[k,1]],dat[,combo[k,2]],use="pairwise.complete.obs")
  }
  vartheta <- temp/nc
  rm(temp)

  # factor loading
  for (k in 1:np){
    if (k!=anchor){
      alpha1[k] <- cov(dat[,anchor],dat[,k],use="pairwise.complete.obs")/vartheta
    }
  }

  # variance of measurement error
  for (k in 1:np){
    varnu[k] <- var(dat[,k],na.rm=TRUE) - (alpha1[k]^2) * vartheta
  }

  # intercept
  for (k in 1:np){
    if (k!=anchor){
      alpha0[k] <- mean(dat[,k],na.rm=TRUE) - alpha1[k] * mtheta
    }
  }

  return(list(alpha0=alpha0,alpha1=alpha1,varnu=varnu,mtheta=mtheta,vartheta=vartheta))
}
