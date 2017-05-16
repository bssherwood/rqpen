\name{cv.rq.group.pen}
\alias{cv.rq.group.pen}
\title{Cross Validated quantile regression with group penalty}
\usage{
cv.rq.group.pen(x, y, groups, tau = 0.5, lambda = NULL, penalty = "LASSO", 
    intercept = TRUE, criteria = "CV", cvFunc = "check", nfolds = 10, 
    foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="QICD", 
	penGroups=NULL,
    ...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{groups}{ Vector assigning columns of x to groups.}
  \item{tau}{ Conditional quantile being modelled.}            
  \item{lambda}{ Vector of lambdas. Default is for lambdas to be automatically generated.} 
  \item{penalty}{ Type of penalty: "LASSO", "SCAD" or "MCP".}
  \item{intercept}{ Whether model should include an intercept. Constant does not need to be included in "x".}
  \item{criteria}{ How models will be evaluated. Either cross-validation "CV", BIC "BIC" or large P BIC "PBIC".}
  \item{cvFunc}{ If cross-validation is used how errors are evaluated. Check function "check", "SqErr" (Squared Error) or "AE" (Absolute Value).}
  \item{nfolds}{ K for K-folds cross-validation.}
  \item{foldid}{ Group id for cross-validation. Function will randomly generate groups if not specified.}
  \item{nlambda}{ Number of lambdas for which models are fit.}
  \item{eps}{ Smallest lambda used.}
  \item{init.lambda}{ Initial lambda used to find the maximum lambda. Not needed if lambda values are set.}
  \item{alg}{Algorithm used for fit. "QICD" or "LP".}
  \item{penGroups}{Specify which groups will be penalized. Default is to penalize all groups.}
  \item{...}{Additional arguments to be sent to rq.group.fit or groupQICDMultLambda.}   
}
\value{
Returns the following:          
\item{beta}{ Matrix of coefficients for different values of lambda}
\item{residuals}{ Matrix of residuals for different values of lambda.}
\item{rho}{Vector of rho, unpenalized portion of the objective function, for different values of lambda.}
\item{cv}{ Data frame with "lambda" and second column is the evaluation based on the criteria selected.}
\item{lambda.min}{ Lambda which provides the smallest statistic for the selected criteria.}
\item{penalty}{ Penalty selected.} 
\item{intercept}{Whether intercept was included in model.}
\item{groups}{Group structure for penalty function.}
}
\description{
Similar to cv.rq.pen function, but uses group penalty. Group penalties use the L1 norm instead of L2 for computational convenience. QICD is a group penalty extension of the algorithm presented by Peng and Wang (2015). LP does a linear programming version of the group penalty. 
}
\examples{                           
\dontrun{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- cv.rq.group.pen(x,y,groups=c(rep(1,4),rep(2,4)))
}
}
\references{
[1] Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.

[2] Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.

}
\author{Ben Sherwood}