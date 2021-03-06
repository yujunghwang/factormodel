% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dproxyme.R
\name{dproxyme}
\alias{dproxyme}
\title{dproxyme}
\usage{
dproxyme(
  dat,
  sbar = 2,
  initvar = 1,
  initvec = NULL,
  seed = 210313,
  tol = 0.005,
  maxiter = 200,
  miniter = 10,
  minobs = 100,
  maxiter2 = 1000,
  trace = FALSE,
  weights = NULL
)
}
\arguments{
\item{dat}{A proxy variable data frame list.}

\item{sbar}{A number of discrete types. Default is 2.}

\item{initvar}{A column index of a proxy variable to initialize the EM algorithm. Default is 1. That is, the proxy variable in the first column of "dat" is used for initialization.}

\item{initvec}{This vector defines how to group the initvar to initialize the EM algorithm.}

\item{seed}{Seed. Default is 210313 (birthday of this package).}

\item{tol}{A tolerance for EM algorithm. Default is 0.005.}

\item{maxiter}{A maximum number of iterations for EM algorithm. Default is 200.}

\item{miniter}{A minimum number of iterations for EM algorithm. Default is 10.}

\item{minobs}{Compute likelihood of a proxy variable only if there are more than "minobs" observations. Default is 100.}

\item{maxiter2}{Maximum number of iterations for "multinom". Default is 1000.}

\item{trace}{Whether to trace EM algorithm progress. Default is FALSE.}

\item{weights}{An optional weight vector}
}
\value{
Returns a list of 5 components : \describe{
\item{M_param}{This is a list of estimated measurement (stochastic) matrices.
The k-th matrix is a measurement matrix of a proxy variable saved in the kth column of dat data frame (or matrix).
The ij-th element in a measurement matrix is the conditional probability of observing j-th (largest) proxy response value conditional on that the latent type is i.}

\item{M_param_col}{This is a list of column labels of 'M_param' matrices}

\item{M_param_row}{This is a list of row labels of 'M_param' matrices. It is simply c(1:sbar).}

\item{mparam}{This is a list of multinomial logit coefficients which were used to compute 'M_param' matrices. These coefficients are useful to compute the likelihood of proxy responses.}

\item{typeprob}{This is a type probability matrix of size N-by-sbar. The ij-th entry of this matrix gives the probability of observation i to have type j.}}
}
\description{
This function estimates measurement stochastic matrices of discrete proxy variables.
}
\examples{
dat1 <- data.frame(proxy1=c(1,2,3),proxy2=c(2,3,4),proxy3=c(4,3,2))
## default minimum num of obs to run an EM algorithm is 10
dproxyme(dat=dat1,sbar=2,initvar=1,minobs=3)
## you can specify weights
dproxyme(dat=dat1,sbar=2,initvar=1,minobs=3,weights=c(0.1,0.5,0.4))


}
\references{
\describe{
\item{Dempster, Arthur P., Nan M. Laird, and Donald B. Rubin (1977)}{"Maximum likelihood from incomplete data via the EM algorithm." Journal of the Royal Statistical Society: Series B (Methodological) 39.1 : 1-22. \doi{10.1111/j.2517-6161.1977.tb01600.x}}
\item{Hu, Yingyao (2008)}{Identification and estimation of nonlinear models with misclassification error using instrumental variables: A general solution. Journal of Econometrics, 144(1), 27-61. \doi{10.1016/j.jeconom.2007.12.001}}
\item{Hu, Yingyao (2017)}{The econometrics of unobservables: Applications of measurement error models in empirical industrial organization and labor economics. Journal of Econometrics, 200(2), 154-168. \doi{10.1016/j.jeconom.2017.06.002}}
\item{Hwang, Yujung (2021)}{Identification and Estimation of a Dynamic Discrete Choice Models with Endogenous Time-Varying Unobservable States Using Proxies. Working Paper. \doi{10.2139/ssrn.3535098}}
\item{Hwang, Yujung (2021)}{Bounding Omitted Variable Bias Using Auxiliary Data. Working Paper. \doi{10.2139/ssrn.3866876}}}
}
\author{
Yujung Hwang, \email{yujungghwang@gmail.com}
}
