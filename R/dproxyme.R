#' @title dproxyme
#' @description This function estimates measurement stochastic matrices of discrete proxy variables.
#' @author Yujung Hwang, \email{yujungghwang@gmail.com}
#' @references \describe{
#'  \item{Dempster, Arthur P., Nan M. Laird, and Donald B. Rubin (1977)}{"Maximum likelihood from incomplete data via the EM algorithm." Journal of the Royal Statistical Society: Series B (Methodological) 39.1 : 1-22. \doi{10.1111/j.2517-6161.1977.tb01600.x}}
#'  \item{Hu, Yingyao (2008)}{Identification and estimation of nonlinear models with misclassification error using instrumental variables: A general solution. Journal of Econometrics, 144(1), 27-61. \doi{10.1016/j.jeconom.2007.12.001}}
#'  \item{Hu, Yingyao (2017)}{The econometrics of unobservables: Applications of measurement error models in empirical industrial organization and labor economics. Journal of Econometrics, 200(2), 154-168. \doi{10.1016/j.jeconom.2017.06.002}}
#'  \item{Hwang, Yujung (2021)}{Identification and Estimation of a Dynamic Discrete Choice Models with Endogenous Time-Varying Unobservable States Using Proxies. Working Paper. \doi{10.2139/ssrn.3535098}}
#'  \item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper. \doi{10.2139/ssrn.3866876}}}
#' @importFrom utils install.packages
#' @importFrom nnet multinom
#' @import     stats
#' @importFrom pracma eye
#'
#' @param dat A proxy variable data frame list.
#' @param sbar A number of discrete types. Default is 2.
#' @param initvar A column index of a proxy variable to initialize the EM algorithm. Default is 1. That is, the proxy variable in the first column of "dat" is used for initialization.
#' @param initvec This vector defines how to group the initvar to initialize the EM algorithm.
#' @param seed Seed. Default is 210313 (birthday of this package).
#' @param tol A tolerance for EM algorithm. Default is 0.005.
#' @param maxiter A maximum number of iterations for EM algorithm. Default is 200.
#' @param miniter A minimum number of iterations for EM algorithm. Default is 10.
#' @param minobs Compute likelihood of a proxy variable only if there are more than "minobs" observations. Default is 100.
#' @param maxiter2 Maximum number of iterations for "multinom". Default is 1000.
#' @param trace Whether to trace EM algorithm progress. Default is FALSE.
#' @param weights An optional weight vector
#'
#' @return Returns a list of 5 components : \describe{
#' \item{M_param}{This is a list of estimated measurement (stochastic) matrices.
#'               The k-th matrix is a measurement matrix of a proxy variable saved in the kth column of dat data frame (or matrix).
#'               The ij-th element in a measurement matrix is the conditional probability of observing j-th (largest) proxy response value conditional on that the latent type is i.}
#'
#' \item{M_param_col}{This is a list of column labels of 'M_param' matrices}
#'
#' \item{M_param_row}{This is a list of row labels of 'M_param' matrices. It is simply c(1:sbar).}
#'
#' \item{mparam}{This is a list of multinomial logit coefficients which were used to compute 'M_param' matrices. These coefficients are useful to compute the likelihood of proxy responses.}
#'
#' \item{typeprob}{This is a type probability matrix of size N-by-sbar. The ij-th entry of this matrix gives the probability of observation i to have type j.}}
#'
#' @examples
#' dat1 <- data.frame(proxy1=c(1,2,3),proxy2=c(2,3,4),proxy3=c(4,3,2))
#' ## default minimum num of obs to run an EM algorithm is 10
#' dproxyme(dat=dat1,sbar=2,initvar=1,minobs=3)
#' ## you can specify weights
#' dproxyme(dat=dat1,sbar=2,initvar=1,minobs=3,weights=c(0.1,0.5,0.4))
#'
#'
#' @export
dproxyme <- function(dat,sbar=2,initvar=1,initvec=NULL,seed=210313,tol=0.005,maxiter=200,miniter=10,minobs=100,maxiter2=1000,trace=FALSE,weights=NULL){

  set.seed(seed)

  # load libraries
  requireNamespace("nnet")
  requireNamespace("pracma")
  requireNamespace("stats")
  requireNamespace("utils")

  #############
  # check if inputs are there in a correct form
  #############

  if (!is.matrix(dat) & !is.data.frame(dat)){
    stop("please provide proxy data in either matrix or data frame format.")
  }

  if (sbar<2){
    stop("please enter the number of types greater than or equal to 2.")
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

  pdat <- list()
  for (k in 1:np){
    pdat[[k]] <- makeDummy(tZ=kronecker(rep(1,sbar),dat[,k]))
  }

  M_param_col <- list()
  M_param_row <- list()
  for (k in 1:np){
    # M_param matrices column names
    M_param_col[[k]] <- sort(unique(dat[,k]))
    # M_param matrices row names
    M_param_row[[k]] <- c(1:sbar)
  }

  # regressor for proxies
  ttype <- kronecker(c(1:sbar),rep(1,N))
  ittype <- as.integer(ttype==2)
  if (sbar>2){
    for (u1 in 3:sbar){
      ittype <- cbind(ittype,as.integer(ttype==u1))
    }
  }

  # initialize the type prob
  tinit <- kronecker(rep(1,sbar),dat[,initvar])

  # group values of dat[,initvar] if initvec is not provided
  if (is.null(initvec)){

    izval <- sort(unique(dat[,initvar]))
    inz <- length(izval)

    initvec <- list( izval[c(1:floor(inz/sbar))] )
    for (u1 in 2:sbar){
      initvec[[u1]] <- izval[c((floor((u1-1)*inz/sbar)+1):floor(u1*inz/sbar))]
    }

  }

  tqst <- rep(NA,N*sbar)
  for (u1 in 1:sbar){
    tqst[ttype==u1] <- (tinit[ttype==u1]%in%initvec[[u1]])+0.01*runif(sum(ttype==u1))
  }
  dim(tqst) <- c(N,sbar)

  aa <- apply(tqst,1,sum)
  for (u1 in 1:sbar){
    tqst[,u1]<- tqst[,u1]/aa
  }
  rm(aa)
  dim(tqst)<-c(N*sbar)

  # set prior
  prior <- matrix(1/sbar,nrow=N,ncol=sbar)

  tlm <- matrix(rep(1,N*sbar*np),ncol=np)
  tqst <- matrix(tqst,ncol=1) ### matrix to use multinom
  typemat <- cbind(rep(1,sbar),eye(sbar)[,2:sbar])


  if (!is.null(weights)){
    tweights <- kronecker(rep(1,sbar),weights)
    dim(tweights) <- c(N*sbar,1)
  }

  ####################
  # start EM
  ####################

  i    <-1
  dif  <-1

  avll <- rep(NA,maxiter)

  mparam   <- list()
  M_param  <- list()
  CM_param <- list()

  while ((i<=maxiter & dif>tol) | i <=miniter ){

    otqst <- tqst

    for (k in 1:np){
      if (sum(is.na(pdat[[k]])==0)>minobs){
        # run estimation if there are sufficiently many proxy obs.

        if (is.null(weights)){
          mm <- multinom(pdat[[k]]~ittype,weights=tqst,maxit=maxiter2,trace=FALSE)
        } else{
          mm <- multinom(pdat[[k]]~ittype,weights=(tqst*tweights),maxit=maxiter2,trace=FALSE)
        }

        mparam[[k]]   <- coef(mm)
        temp <- cbind(rep(0,dim(typemat)[2]),t(coef(mm)))
        numer <- exp(typemat%*%temp)
        denom <- matrix(rep(apply(numer,1,sum),length(M_param_col[[k]])),ncol=length(M_param_col[[k]]))
        M_param[[k]] <- numer/denom
        rm(temp)

        CM_param[[k]] <- t(apply(M_param[[k]],1,cumsum))

        nc <- dim(M_param[[k]])[2]

        for (v in 1:sbar){
          tlm[ttype==v,k] <- apply(t(matrix(rep(M_param[[k]][v,],sum(ttype==v)),nrow=nc)) * pdat[[k]][ttype==v,],1,sum)
        }
      }
    }

    tlm[is.na(tlm)] <-1
    tlm2 <- matrix(apply(tlm,1,prod),ncol=sbar)

    # update tqst
    num <- prior*tlm2
    den <- matrix(rep(apply(num,1,sum),sbar),ncol=sbar)
    tqst <- num/den

    # update prior
    pp <- apply(tqst,2,mean)
    prior <- rep(pp[1], N)
    for (k in 1:(sbar-1)){
      prior <- cbind(prior, rep(pp[k+1],N))
    }

    # reformat tqst
    dim(tqst) <- c(N*sbar,1)

    # compute dif
    dif <- max(abs(tqst-otqst))

    if (trace==TRUE){
      print(paste0("Iteration : ",i, " , difference : ",dif))
    }

    i<-i+1
  }

  # return type probability
  dim(tqst) <- c(N,sbar)

  return(list(M_param=M_param,M_param_col=M_param_col,M_param_row=M_param_row,mparam=mparam,typeprob=tqst))
}
