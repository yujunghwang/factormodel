factormodel
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Introduction

The R package **factormodel** provides functions to estimate a factor
model using either discrete or continuous proxy variables. Such model is
useful when proxy variables include measurement errors.

When proxy variables are discrete, you can use ‘dproxyme’ function. The
function estimates a finite-mixture model using an EM algorithm
(Dempster, Laird, Rubin, 1977).

A function ‘dproxyme’ returns a list of estimated measurement
(stochastic) matrices and a type probability matrix from discrete proxy
variables. The ij-th element in a measurement matrix is the conditional
probability of observing j-th (largest) proxy response value conditional
on that the latent type is i. The type probability matrix is of size N
(num of obs) by sbar (num of type). The ij-th element of the type
probability is the probability of observation i to belong to the type j.
For further explanation on identification of measurement stochastic
matrices, see Hu(2008) and Hu(2017).

When proxy variables are continuous, you can use ‘cproxyme’ function.
The function estimates a linear factor model assuming a continuous
latent variable.

A function ‘cproxyme’ returns a list of linear factor model coefficients
and the variance of measurement errors in each proxy variable. For
further explanation on identification of linear factor model, see Cunha,
Heckman, Schennach (2010).

## Installation

You can install a package **factormodel** using either CRAN or github.

``` r
install.packages("factormodel")
```

or

``` r
# install.packages("devtools")
devtools::install_github("yujunghwang/factormodel")
```

## Example 1 : dproxyme

Below example shows how to use ‘dproxyme’ function to estimate
measurement (stochastic) matrices of proxy variables and type
probabilities. The code first simulates fake data using a data
generating process provided below and then estimates the parameters
using ‘dproxyme’ function.

``` r
library(factormodel)
library(nnet)
library(pracma)
library(stats)
library(utils)

# DGP
# set parameters
nsam <- 5000

M1 <- rbind(c(0.8,0.1,0.1),c(0.1,0.2,0.7))
M2 <- rbind(c(0.7,0.2,0.1),c(0.2,0.2,0.6))
M3 <- rbind(c(0.9,0.05,0.05),c(0.1,0.1,0.8))

CM1 <- t(apply(M1,1,cumsum))
CM2 <- t(apply(M2,1,cumsum))
CM3 <- t(apply(M3,1,cumsum))

# 40% of sample is type 1, 60% is type 2
truetype <- as.integer(runif(nsam)<=0.4) +1

# generate fake data
dat <- data.frame(msr1=rep(NA,nsam),msr2=rep(NA,nsam),msr3=rep(NA,nsam))

for (k in 1:nsam){
  dat$msr1[k] <- which(runif(1)<=CM1[truetype[k],])[1]
  dat$msr2[k] <- which(runif(1)<=CM2[truetype[k],])[1]
  dat$msr3[k] <- which(runif(1)<=CM3[truetype[k],])[1]
}

# estimate using dproxyme
oout <- dproxyme(dat=dat,sbar=2,initvar=1,initvec=NULL,seed=210313,tol=0.005,maxiter=200,miniter=10,minobs=100,maxiter2=1000,trace=FALSE)

# check whether the estimated measurement stochastic matrices are same with the true # measurement stochastic matrices
print(oout$M_param)
#> [[1]]
#>                         2          3
#> [1,] 0.8109318 0.09434914 0.09471906
#> [2,] 0.1005454 0.21215160 0.68730298
#> 
#> [[2]]
#>                        2          3
#> [1,] 0.7004681 0.2075708 0.09196105
#> [2,] 0.1955446 0.2081434 0.59631203
#> 
#> [[3]]
#>                         2          3
#> [1,] 0.9038930 0.04134642 0.05476054
#> [2,] 0.1069748 0.10807967 0.78494554

# check type probability
print(head(oout$typeprob))
#>             [,1]        [,2]
#> [1,] 0.002224585 0.997775415
#> [2,] 0.990309611 0.009690389
#> [3,] 0.940488990 0.059511010
#> [4,] 0.990309611 0.009690389
#> [5,] 0.997283248 0.002716752
#> [6,] 0.997283248 0.002716752

# compare this to true type
print(head(truetype))
#> [1] 2 1 1 1 1 1
```

## Example 2 : cproxyme

Below example shows how to use ‘cproxyme’ function to estimate a linear
factor model. The code first simulates fake data using a data generating
process provided below and then estimates the parameters using
‘cproxyme’ function.

``` r
library(factormodel)
library(stats)
library(utils)
library(gtools)
#> 
#> Attaching package: 'gtools'
#> The following object is masked from 'package:pracma':
#> 
#>     logit

set.seed(seed=210315)

# DGP
# set parameters
nsam <- 5000 # number of observations
np <- 3 # number of proxies

true_mtheta <- 2
true_vartheta <- 1.5
true_theta <- rnorm(nsam, mean=true_mtheta, sd=sqrt(true_vartheta))

# first proxy variable is an anchoring variable
true_alpha0 <- c(0,2,5)
true_alpha1 <- c(1,0.5,2)
true_varnu  <- c(0.5,2,1)

# simulate fake data
dat <- matrix(NA,nrow=nsam,ncol=np)
for (k in 1:np){
  dat[,k] <- true_alpha0[k] + true_alpha1[k]*true_theta + rnorm(nsam,mean=0,sd=sqrt(true_varnu[k]))
}

# estimate parameters using cproxyme
oout <- cproxyme(dat=dat,anchor=1)

# print estimated parameters
print(oout$alpha0)
#> [1] 0.000000 2.032455 5.086428
print(oout$alpha1)
#> [1] 1.0000000 0.4913708 1.9630702
print(oout$varnu)
#> [1] 0.4827664 2.0430161 1.0605388
print(oout$mtheta)
#> [1] 1.990616
print(oout$vartheta)
#> [1] 1.586096
```

## Conclusion

This vignette showed how to use functions in \`factormodel’ R package.

## References

  - [Cunha, F., Heckman, J. J., & Schennach, S. M. (2010). Estimating
    the technology of cognitive and noncognitive skill formation.
    Econometrica, 78(3), 883-931.](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA6551?casa_token=MNFL_7OY05UAAAAA:zysMye4e8rFMMDnBzvu9D1LWJ1XFEa9nhQkI0jl6lXWlKsy4xskj6qmrUHFJgNaRxDS1YTlUR8LiOQ)

  - [Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum
    likelihood from incomplete data via the EM algorithm. Journal of the
    Royal Statistical Society: Series B
    (Methodological), 39(1), 1-22.](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1977.tb01600.x)

  - [Hu, Yingyao (2008). Identification and estimation of nonlinear
    models with misclassification error using instrumental variables: A
    general solution. Journal of
    Econometrics, 144(1), 27-61.](https://www.sciencedirect.com/science/article/pii/S0304407607002436?casa_token=b9ManDs-MlQAAAAA:E02Ae5SIzmrGbIbCAFeSk-BI2pR9ZcZMSc7q28S8VJVDzj0gl-sKOS9fWTklX7nydQixJxMr)

  - [Hu, Yingyao (2017). The econometrics of unobservables: Applications
    of measurement error models in empirical industrial organization and
    labor economics. Journal of
    econometrics, 200(2), 154-168.](https://www.sciencedirect.com/science/article/pii/S0304407617300830?casa_token=gGghKpWlo0kAAAAA:4DcT91SeVK56FF1XFsn34pMKnBLv46VBBM3zhBj7Mj4q2q94LSZCMY0LHRtxUHvCma1QVjDV)

  - [Hwang, Yujung (2021). Identification and Estimation of a Dynamic
    Discrete Choice Models with Endogenous Time-Varying Unobservable
    States Using Proxies. Working
    Paper.](https://sites.google.com/view/yujunghwang/research?authuser=0)

  - [Hwang, Yujung (2021). Bounding Omitted Variable Bias Using
    Auxiliary Data. Working
    Paper.](https://sites.google.com/view/yujunghwang/research?authuser=0)
