\name{cv.rq.pen}
\alias{cv.rq.pen}
\title{Cross Validated quantile regression}
\usage{
cv.rq.pen(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",
          intercept=TRUE,criteria="CV",cvFunc="check",nfolds=10,
          foldid=NULL,nlambda=100,eps=.0001,init.lambda=1, penVars=NULL,...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{tau}{ Conditional quantile being modelled.}            
  \item{lambda}{ Vector of lambdas. Default is for lambdas to be automatically generated.}
  \item{weights}{ Weights for the objective function.} 
  \item{penalty}{ Type of penalty: "LASSO", "SCAD" or "MCP".}
  \item{intercept}{ Whether model should include an intercept. Constant does not need to be included in "x".}
  \item{criteria}{ How models will be evaluated. Either cross-validation "CV", BIC "BIC" or large P BIC "PBIC".}
  \item{cvFunc}{ If cross-validation is used how errors are evaluated. Check function "check", "SqErr" (Squared Error) or "AE" (Absolute Value).}
  \item{nfolds}{ K for K-folds cross-validation.}
  \item{foldid}{ Group id for cross-validation. Function will randomly generate groups if not specified.}
  \item{nlambda}{ Number of lambdas for which models are fit.}
  \item{eps}{ Smallest lambda used.}
  \item{init.lambda}{ Initial lambda used to find the maximum lambda. Not needed if lambda values are set.}
  \item{penVars}{Variables that should be penalized. With default value of NULL all variables are penalized.}
  \item{...}{Additional arguments to be sent to rq.lasso.fit or rq.nc.fit.}   
}
\value{
Returns the following:
\item{models}{ List of penalized models fit. Number of models will match number of lambdas and correspond to cv$lambda.}
\item{cv}{ Data frame with "lambda" and second column is the evaluation based on the criteria selected.}
\item{lambda.min}{ Lambda which provides the smallest statistic for the selected criteria.}
\item{penalty}{ Penalty selected.} 
}
\description{
Produces penalized quantile regression models for a range of lambdas and penalty of choice. If lambda is unselected than an iterative algorithm is used to find a maximum lambda such that the penalty is large enough to produce an intercept only model. Then range of lambdas goes from the maximum lambda found to "eps" on the log scale. For non-convex penalties local linear approximation approach used by Wang, Wu and Li to extend LLA as proposed by Zou and Li (2008) to the quantile regression setting.  
}
\examples{
\dontrun{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- cv.rq.pen(x,y)
}
}
\references{
[1] Wang, L., Wu, Y. and Li, R. Quantile regression of analyzing heterogeneity in ultra-high dimension. \emph{J. Am. Statist. Ass}, \bold{107}, 214--222.

[2] Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.

[3] Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models. \emph{Ann. Statist.}, \bold{36}, 1509--1533.
}
\author{Ben Sherwood}