# rqPen: Penalized quantile regression

## Overview

This R package provides tools for estimating a quantile regression model with a penalized objective function. Implements a variety of 
penalties, including group penalties. 

## Installation

For most up to date versions use the following code. However, be warned the github package is often in a state of testing and debugging.
``` r
devtools::install_github("bssherwood/rqpen")
```

## Example

``` r
library(rqPen)
n<- 200
p<- 30
x0<- matrix(rnorm(n*p),n,p)
X<- cbind(x0, x0^2, x0^3)[,order(rep(1:p,3))]
y<- -2+X[,1]+0.5*X[,2]-X[,3]-0.5*X[,7]+X[,8]-0.2*X[,9]+rt(n,2)
group<- rep(1:p, each=3)

# lasso estimation
# one tau
fit1 <- rq.pen(x,y)
# several values of tau
fit2 <- rq.pen(x,y,tau=c(.2,.5,.8))

# Group SCAD estimation
fit3 <- rq.group.pen(x,y,groups=group,penalty="gSCAD")

# cross validation
cv1 <- rq.pen.cv(x,y)
plot(cv1)

cv2 <- rq.pen.cv(x,y,tau=c(.2,.5,.8))
plot(cv2)

cv3 <- rq.group.pen(x,y,groups=group,penalty="gSCAD")
plot(cv3)

# BIC selection of tuning parameters
qs1 <- qic.select(fit1)
qs2 <- qic.select(fit2)
qs3 <- qic.select(fit3)
```

See, https://github.com/bssherwood/rqpen/blob/master/ignore/rqPenArticle.pdf, for a vignette. The Huber approach for rq.pen relies on the R package hqreg and work presented in "Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression". The Huber approach in rq.group.pen relies on R package hrqglas and work presented in An Efficient Approach to Feature Selection and Estimation for Quantile Regression with Grouped Variables


## References

Sherwood, B. and Li, S. (2022) An Efficient Approach to Feature Selection and Estimation for Quantile Regression with Grouped Variables, Statistics and computing, 75. 

Yi, C. and Huang, J. (2015) Semismooth Newton Coordinate Descent Algorithm for Elastic-Net Penalized Huber Loss Regression and Quantile Regression, Journal of Computational and Graphical Statistics, 26:3, 547-557. 
