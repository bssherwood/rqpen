\name{QICD.nonpen}
\alias{QICD.nonpen}
\title{Penalized Quantile Regression with some nonpenalized coefficients 
with QICD Algorithm}
\usage{
QICD.nonpen(y, x, z, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                 initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, 
				 coef.cutoff=1e-08, a=3.7, method="br", ...)
}
\arguments{
  \item{y}{ Vector of response values.}
  \item{x}{ n x p matrix of observations with each row corresponding to one observation. Penalties (and variable selection) will be applied to these coefficients. Do not include column of 1's; set intercept=TRUE, if intercept is desired.}
  \item{z}{ n x q matrix of observations with each row corresponding to one observation. Penalties (and variable selection) WILL NOT be applied to these coefficients. Do not include column of 1's; set intercept=TRUE, if intercept is desired.}
  \item{tau}{ Conditional quantile being modelled.}
  \item{lambda}{ Tuning parameter for LASSO, SCAD, and MCP penalties. Must be positive.}
  \item{intercept}{ If TRUE, an intercept is included in the model.  If FALSE, no intercept is included.}
  \item{penalty}{ Penalty function for the coefficients.  Either "SCAD", "MCP", or "LASSO". See details for description of penalties.}
  \item{initial_beta}{ Vector of initial values for QICD algorithm. The vector should contain the intercept first (if intercept=TRUE) and then the p coefficients correspinding to x.  Initial values for the coefficients corresponding to z can be passed after the the x coefficients, but will be ignored. If NULL, exact LASSO estimates will be computed and used as initial values.}
  \item{maxin}{ Maximum number of iterations on the minimization step of the QICD algorithm.}
  \item{maxout}{ Maximum number of iterations on the majorization step of the QICD algorithm.}
  \item{eps}{ Threshold for convergence of algorithm.}
  \item{coef.cutoff}{ Coefficients with magnitude less than this value are set to 0.}
  \item{a}{ Second tuning parameter for SCAD and MCP penalties.  Must be greater than 2 for SCAD and greater than 1 for MCP.}
  \item{method}{ Method used in QICD algorithm. Default is "br" method. When sample size is several thousand, method "fn" should be used for increased speed. See quantreg package for more details.}
  \item{...}{ The extra arguments will not be used.}
}
\value{
Returns a vector containing the intercept (if intercept=TRUE) and the estimated coefficients for each column in x and z.
}
\description{
QICD.nonpen produces penalized quantile regression estimates with some nonpenalized coefficients using the QICD algorithm.  If no initial values are given, LASSO estimates will be used.  This function can handle the LASSO, SCAD, and MCP penalties. This can be useful when you would like to perform variable selection only on some covariates and would like to guarantee that other covariates remain in the model.
}
\details{
The QICD algorithm should only be used for the LASSO penalty if initial_beta can be set to LASSO estimates with a similar lambda (similar to a "warm start").  Otherwise, exact LASSO estimates will be used as initial values for the QICD algorithm: this will cause unnecessary computations and could lead to less accurate estimates.
}
\examples{
library(splines)
n = 50
p = 100
x <- matrix(rnorm(n*p), nrow=n)
z1 <- runif(n)
z2 <- runif(n)
z <- cbind(bs(z1), bs(z2))
y <- 0 + x[,1] - 3*x[,5] + z1^3 + sin(2*pi*z2) + rnorm(n)
fit1 <- QICD.nonpen(y,x,z, tau=.5, lambda=1, intercept=TRUE, penalty="SCAD")
fit2 <- QICD.nonpen(y,x,z, tau=.7, lambda=1, intercept=TRUE, penalty="SCAD")
}
\references{
[1] Wang, L., Wu, Y. and Li, R. (2012). Quantile regression of analyzing heterogeneity in ultra-high dimension. \emph{J. Am. Statist. Ass}, \bold{107}, 214--222.

[2] Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.

[3] Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models. \emph{Ann. Statist.}, \bold{36}, 1509--1533.

[4] Peng, B. and Wang, L. (2015). An iterative coordinate-descent algorithm for high-dimensional nonconvex penalized quantile regression. \emph{J. Comp. Graph.}, \bold{24}, 676--694.
}
\author{Adam Maidman}