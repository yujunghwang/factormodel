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
(stochastic) matrices from discrete proxy variables. The ij-th element
in a measurement matrix is the conditional probability of observing j-th
(largest) proxy response value conditional on that the latent type is i.
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

``` r
library(factormodel)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
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
#> [1,] 0.8160435 0.09073877 0.09321774
#> [2,] 0.1053269 0.19131273 0.70336038
#> 
#> [[2]]
#>                        2          3
#> [1,] 0.7110982 0.1894228 0.09947892
#> [2,] 0.2058876 0.2022934 0.59181899
#> 
#> [[3]]
#>                         2          3
#> [1,] 0.9102475 0.04229629 0.04745621
#> [2,] 0.1295547 0.10203853 0.76840678
```

## Example 2 : cproxyme

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
    econometrics, 200(2), 154-168.](https://www.sciencedirect.com/science/article/pii/S0304407617300830?casa_token=gGghKpWlo0kAAAAA:4DcT91SeVK56FF1XFsn34pMKnBLv46VBBM3zhBj7Mj4q2q94LSZCMY0LHRtxUHvCma1QVjDV).
