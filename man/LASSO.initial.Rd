\name{LASSO.initial}
\alias{LASSO.initial}
\title{LASSO Penalized Quantile Regression}
\usage{
LASSO.initial(y, x, tau, lambda, intercept, coef.cutoff)
}
\arguments{
  \item{y}{ Vector of responses.} 
  \item{x}{ n x p matrix of covariates.}
  \item{tau}{ Conditional quantile being modelled.}
  \item{lambda}{ Tuning parameter. Must be positive.}
  \item{intercept}{ If TRUE, an intercept is included in the model.  If FALSE, no intercept is included.}
  \item{coef.cutoff}{ Coefficients with magnitude less than this value are set to 0.}
}
\value{
Returns a vector containing the intercept (if intercept=TRUE) and the estimated coefficients for each column in x.
}
\description{
LASSO.initial obtains coefficient estimates for Lasso penalized quantile regression.  It is called by the QICD and QICD.group functions to obtain initial estimates when they are not provided.
}
\details{
This is a barebones function that only provides coefficient estimates. It will not provide any warnings or errors, so you need to check that inputs are accurate and appropriate.  The rq.lasso.fit function can be used to obtain more information from the Lasso fit.
}

\examples{
n = 50
p = 100
x <- matrix(rnorm(n*p),nrow=n)
y <- 0 + x[,1] - 3*x[,5] + rnorm(n)
fit1 <- LASSO.initial(y,x, tau=.5, lambda=1, intercept=TRUE, coef.cutoff=1e-08)
fit2 <- LASSO.initial(y,x, tau=.5, lambda=.1, intercept=TRUE, coef.cutoff=1e-08)
}
\references{
[1] Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
\emph{Journal of the Royal Statistical Society. Series B}, \bold{58}, 267--288.
}
\author{Adam Maidman}