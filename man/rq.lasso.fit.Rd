\name{rq.lasso.fit}
\alias{rq.lasso.fit}
\title{Quantile Regression with LASSO penalty}
\usage{
rq.lasso.fit(x,y,tau=.5,lambda=NULL,weights=NULL,
             intercept=TRUE,coef.cutoff=.00000001,
             method=NULL, penVars=NULL, ...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{tau}{ Conditional quantile being modelled.}   
  \item{lambda}{ Tuning parameter.}
  \item{weights}{ Weights for the objective function.}
  \item{intercept}{ Whether model should include an intercept. Constant does not need to be included in "x".}
  \item{coef.cutoff}{ Coefficients below this value will be set to zero.}
  \item{method}{Use method "br" or "fn" as outlined in quantreg package. Method "fn" better for large n or p. Default will choose "br" or "fn" based on n and p.} 
  \item{penVars}{Variables that should be penalized. With default value of NULL all variables are penalized.}
  \item{...}{Additional items to be sent to rq. Note this will have to be done carefully as rq is run on the augmented data to account for penalization and could provide strange results if this is not taken into account.}
}
\value{
Returns the following:
\item{coefficients}{ Coefficients from the penalized model.} 
\item{PenRho}{ Penalized objective function value.}
\item{residuals}{ Residuals from the model.}
\item{rho}{ Objective function evaluation without the penalty.}
\item{tau}{ Conditional quantile being modelled.}
\item{n}{ Sample size.}  
}
\description{
Fits a quantile regression model with the LASSO penalty. Algorithm is similar to LASSO code presented in Koenker and Mizera (2014).    
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
lassoModel <- rq.lasso.fit(x,y,lambda=1)
}
\references{
[1] Koenker, R. and Mizera, I. (2014). Convex optimization in R. \emph{Journal of Statistical Software}, \bold{60}, 1--23.  

[2] Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
\emph{Journal of the Royal Statistical Society. Series B}, \bold{58}, 267--288.

[3] Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica
Sinica}, \bold{19}, 801--817.   
}
\author{Ben Sherwood}