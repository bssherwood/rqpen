\name{QICD}
\alias{QICD}
\title{QICD Algorithm for Nonconvex Penalized Quantile Regression}
\usage{
QICD(y, x, tau = 0.5, lambda = NULL, weights=NULL, beta_initial = NULL, 
    intercept = TRUE, penalty = "SCAD", a = 3.7, converge_criteria = 1e-06, 
    maxin = 100, maxout = 20, method = "br",...)
}
\arguments{
  \item{y}{Vector of response values.}
  \item{x}{Matrix of predictors.  Note: If intercept is desired, do not include a column of 1's.  If intercept is desired, then set "intercept=TRUE".}
  \item{tau}{Conditional quantile being modelled}
  \item{lambda}{Tuning parameter.}
  \item{weights}{Weights for the objective function.}
  \item{beta_initial}{Vector of values for coefficients (not intercept) to initialize algorithm.  If NULL, the LASSO model is fit and used as initializing values.}
  \item{intercept}{Whether model should include an intercept. Constant does not need to be included in "x".  Takes values "TRUE" or "FALSE".}
  \item{penalty}{MCP or SCAD.  Takes values "MCP" or "SCAD".}
  \item{a}{Second tuning parameter.}
  \item{converge_criteria}{Difference in betas from iteration process that would satisfy convergence.}
  \item{maxin}{Maximum number of iterations of the coordinate descent part of the QICD algorithm to be done.}
  \item{maxout}{Maximum number of iterations of the Minimize-Majorize (MM) part of the QICD algorithm to be done.}
  \item{method}{Use method "br" or "fn" as outlined in quantreg package. Method "fn" better for large n or p.}
  \item{...}{Additional arguments sent to rq.lasso.fit for initial estimates}
}
\value{
Returns the following:
\item{coefficients}{ Coefficients from the penalized model.} 
\item{PenRho}{ Penalized objective function value.} 
\item{residuals}{ Residuals from the model.}
\item{rho}{ Objective function evaluation without the penalty.}
\item{tau}{ Conditional quantile being modelled.}
\item{n}{ Sample size.}  
\item{penalty}{ Penalty used, SCAD or MCP.} 
}
\description{
Produces penalized quantile regression models for specified lambdas and penalty of choice. For now, lambda must be specified.  The model is fit using the QICD algorithm from Peng and Wang (2015).
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
scadModel <- QICD(y,x,lambda=1, intercept=TRUE, beta_initial=rep(0,8))
}
\references{
[1] Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for 
    High-Dimensional Nonconvex Penalized Quantile Regression. 
	\emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.
}
\author{Adam Maidman}