

# This function estimates measurement equation for discrete proxy variables
dproxyme <- function(dat,sbar,seed=12321){
  # input
  # dat : proxy variable data frame list
  # sbar : number of discrete types
  # seed : seed

  set.seed(seed)

  # install libraries if not installed
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("nnet")) install.packages("nnet")
  if (!require("stats")) install.packages("stats")
  if (!require("pracma")) install.packages("pracma")

  library(dplyr)
  library(nnet)
  library(stats)
  library(pracma)

  #############
  # check if inputs are there in correct form
  #############





  #############
  # prepare data in a right form
  #############

  pdat <- data.frame()
  pdat$id <- kronecker(rep(1,sbar),c(1:dim(dat)[1]))



  # set the data frame size
  N  <- dim(pdat)[1]
  np <- dim(pdat)[2]



  # set prior
  prior <- rep(1/sbar,N*T)
  for (k in 1:(sbar-1)){
    prior <- cbind(prior, rep(1/sbar,N*T))
  }

  tlm <- matrix(rep(1,N*T*sbar*np),ncol=np)

  tqst <- matrix(tqst,ncol=1) ### matrix to use multinom

  #typemat <- data.frame(ittype2=(c(0,1,0)),ittype3=(c(0,0,1)))
  typemat <- cbind(rep(1,sbar),eye(sbar)[,2:sbar])


  # start EM
  i    <-1
  dif  <-1
  tol  <-0.005
  niter<-200

  avll <- rep(NA,niter)

  mparam   <- list()
  M_param  <- list()
  CM_param <- list()

  while ((i<=niter & dif>tol) | i <=10 ){

    otqst <- tqst

    for (k in 1:np){
      if (sum(is.na(pdat[[k]])==0)>100){
        # run estimation if there are sufficiently many proxy obs.
        #mm <- multinom(pdat[[k]]~ittype2+ittype3,weights=tqst,maxit=1000,trace=FALSE)
        mm <- multinom(pdat[[k]]~ittype,weights=tqst,maxit=1000,trace=FALSE)
        mparam[[k]]   <- coef(mm)
        temp <- cbind(rep(0,dim(typemat)[2]),t(coef(mm)))
        numer <- exp(typemat%*%temp)
        denom <- matrix(rep(apply(numer,1,sum),7),ncol=7)
        M_param[[k]] <- numer/denom
        rm(temp)

        #M_param[[k]]  <- predict(mm,newdata=typemat,"probs")
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
    prior <- rep(pp[1], N*T)
    for (k in 1:(sbar-1)){
      prior <- cbind(prior, rep(pp[k+1],N*T))
    }

    # reformat tqst
    dim(tqst) <- c(N*T*sbar,1)

    # compute dif
    dif <- max(abs(tqst-otqst))

    #print(paste0("Iteration : ",i, " , difference : ",dif))

    i<-i+1
  }

  return(list(M_param,mparam,tlm,tqst,prior))
}
