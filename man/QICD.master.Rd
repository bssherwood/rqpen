\name{QICD.master}
\alias{QICD.master}
\title{Master QICD Function for Regular QICD, group QICD, and Partially Penalized QICD with Multiple Lambdas}
\usage{
QICD.master(y, x, z=NULL, groups=NULL, tau=.5, lambda, intercept=TRUE, 
			penalty="SCAD", initial_beta, maxin=100, maxout=20, 
			eps = 1e-05, coef.cutoff=1e-08, a=3.7, ...)
}
\arguments{
  \item{y}{ Vector of response values.}
  \item{x}{ n x p matrix of observations with each row corresponding to one observation. Penalties (and variable selection) will be applied to these coefficients. Do not include column of 1's; set intercept=TRUE, if intercept is desired.}
  \item{z}{ n x q matrix of observations with each row corresponding to one observation. Penalties (and variable selection) WILL NOT be applied to these coefficients. Do not include column of 1's; set intercept=TRUE, if intercept is desired. Set to NULL if all coefficients (except for intercept) should be penalized. Currently no support for z and groups.}
  \item{groups}{ Vector of length p with the group number of the corresponding coefficient.  Coefficients in the same group will either all be 0, or none will be 0. Set to NULL if no groups. Currently no support for groups and z.}
  \item{tau}{ Conditional quantile being modelled.}
  \item{lambda}{ Tuning parameters for LASSO, SCAD, and MCP penalties. Must be positive.}
  \item{intercept}{ If TRUE, an intercept is included in the model.  If FALSE, no intercept is included.}
  \item{penalty}{ Penalty function for the coefficients.  Either "SCAD", "MCP", or "LASSO". See details for description of penalties.}
  \item{initial_beta}{ Vector of initial values for QICD algorithm. The vector should contain the intercept first (if intercept=TRUE) and then the p coefficients correspinding to x. Initial values for the coefficients corresponding to z can be passed after the the x coefficients, but will be ignored. These initial values will be used for all values of lambda. It is recommended to use LASSO estimates (with appropriately chosen lambda) as initial_beta.}
  \item{maxin}{ Maximum number of iterations on the minimization step of the QICD algorithm.}
  \item{maxout}{ Maximum number of iterations on the majorization step of the QICD algorithm.}
  \item{eps}{ Threshold for convergence of algorithm.}
  \item{coef.cutoff}{ Coefficients with magnitude less than this value are set to 0.}
  \item{a}{ Second tuning parameter for SCAD and MCP penalties.  Must be greater than 2 for SCAD and greater than 1 for MCP.}
  \item{...}{ For partially penalized penalties, the method can be changed. See QICD.nonpen function.}
}
\value{
Returns the following:
  \item{coefficients}{ Matrix of estimated coefficients corresponding to each value of lambda. The ith column corresponds to the ith lambda value in lambda.}
  \item{lambda}{ Unique values of lambda sorted into ascending order.}
}
\description{
QICD.master produces penalized quantile regression estimates for all three cases of the QICD algorithm (QICD, QICD.nonpen, QICD.group). This function will find estimates for multiple lambdas.
}
\examples{
n = 50
p = 100
x <- matrix(rnorm(n*p), nrow=n)
y <- 0 + x[,1] - 3*x[,5] + rnorm(n)
lambda <- exp(-20:0)
# fit1 <- QICD.master(y,x, tau=.5, lambda=lambda, intercept=TRUE, penalty="SCAD")
# fit2 <- QICD.master(y,x, tau=.7, lambda=lambda, intercept=TRUE, penalty="SCAD")
# head(fit1)
# fit2
}
\references{
[1] Wang, L., Wu, Y. and Li, R. (2012). Quantile regression of analyzing heterogeneity in ultra-high dimension. \emph{J. Am. Statist. Ass}, \bold{107}, 214--222.

[2] Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.

[3] Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models. \emph{Ann. Statist.}, \bold{36}, 1509--1533.

[4] Peng, B. and Wang, L. (2015). An iterative coordinate-descent algorithm for high-dimensional nonconvex penalized quantile regression. \emph{J. Comp. Graph.}, \bold{24}, 676--694.
}
\author{Adam Maidman}