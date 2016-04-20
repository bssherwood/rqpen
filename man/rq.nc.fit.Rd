\name{rq.nc.fit}
\alias{rq.nc.fit}
\title{Non-convex penalized quantile regression}
\usage{
rq.nc.fit(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
          penalty="SCAD",a=3.7,iterations=10,
          converge_criteria=1e-06,alg=ifelse(p<50,"LP","QICD"),penVars=NULL, ...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{tau}{ Conditional quantile being modeled.}
  \item{lambda}{ Tuning parameter.}  
  \item{weights}{ Weights for the objective function.}
  \item{intercept}{ Whether model should include an intercept. Constant does not need to be included in "x".}
  \item{penalty}{ MCP or SCAD.}
  \item{a}{ Second tuning parameter.} 
  \item{iterations}{ Number of iterations to be done for iterative LLA algorithm.} 
  \item{converge_criteria}{ Difference in betas from iteration process that would satisfy convergence.}
  \item{alg}{ Defaults for small p to linear programming (LP), see Wang, Wu and Li (2012) for details. Otherwise a coordinate descent algorithm is used (QICD), see Peng and Wang (2015) for details. Both methods rely on the One-step sparse estimates algorithm.}
  \item{penVars}{Variables that should be penalized. With default value of NULL all variables are penalized.}
  \item{...}{Additional items to be sent to rq.lasso.fit.}   
}
\value{
Returns the following:
\item{coefficients}{ Coefficients from the penalized model.} 
\item{PenRho}{ Penalized objective function value.} 
\item{residuals}{ Residuals from the model.}
\item{rho}{ Objective function evaluation without the penalty.}
\item{tau}{ Conditional quantile being modeled.}
\item{n}{ Sample size.}  
\item{penalty}{ Penalty used, SCAD or MCP.} 
}
\description{
Produces penalized quantile regression models for a range of lambdas and penalty of choice. If lambda is unselected than an iterative algorithm is used to find a maximum lambda such that the penalty is large enough to produce an intercept only model. Then range of lambdas goes from the maximum lambda found to "eps" on the log scale. Local linear approximation approach used by Wang, Wu and Li to extend LLA as proposed by Zou and Li (2008) to the quantile regression setting.  
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
scadModel <- rq.nc.fit(x,y,lambda=1)
}
\references{
[1] Wang, L., Wu, Y. and Li, R. (2012). Quantile regression of analyzing heterogeneity in ultra-high dimension. \emph{J. Am. Statist. Ass}, \bold{107}, 214--222.

[2] Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.

[3] Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models. \emph{Ann. Statist.}, \bold{36}, 1509--1533.

[4] Peng, B. and Wang, L. (2015). An iterative coordinate-descent algorithm for high-dimensional nonconvex penalized quantile regression. \emph{J. Comp. Graph.}, \bold{24}, 676--694.
}
\author{Ben Sherwood}