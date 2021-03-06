% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted.var.R
\name{weighted.var}
\alias{weighted.var}
\title{weighted.var}
\usage{
weighted.var(x, w = NULL)
}
\arguments{
\item{x}{A vector to compute a variance, var(x)}

\item{w}{A weight vector}
}
\value{
Returns an unbiased sample weighted variance
}
\description{
This function is to compute an unbiased sample weighted variance.
}
\examples{
## If you do not specify weights, 
## it returns the usual unweighted sample variance
weighted.var(x=c(1,3,5)) 

weighted.var(x=c(1,3,5),w=c(0.1,0.5,0.4))

}
\author{
Yujung Hwang, \email{yujungghwang@gmail.com}
}
