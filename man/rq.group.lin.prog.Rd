\name{rq.group.lin.prog}
\alias{rq.group.lin.prog}
\title{Quantile Regresion with Group Penalty using linear programming algorithm}
\usage{
rq.group.lin.prog(x,y,groups,tau,lambda,intercept=TRUE,eps=1e-05,
     penalty="SCAD", a=3.7, coef.cutoff=1e-08, initial_beta=NULL,
     iterations=10,converge_criteria=.0001,...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{groups}{ Vector assigning columns of x to groups.}
  \item{tau}{ Conditional quantile being modelled.}            
  \item{lambda}{ Vector of lambdas. Default is for lambdas to be automatically generated.}   
  \item{intercept}{ Whether model should include an intercept. Constant does not need to be included in "x".}
  \item{eps}{Multiplier for smallest lambda.}
  \item{penalty}{ Type of penalty: "LASSO", "SCAD" or "MCP".}  
  \item{a}{Additional parameter for non-convex penalties.}
  \item{coef.cutoff}{Estimates below this value are set to zero.}
  \item{initial_beta}{Initial beta estimate.}
  \item{iterations}{Maximum number of iterations.}
  \item{converge_criteria}{Convergence criteria}
  \item{...}{Additional arguments to be sent to rq.lasso.fit.}   
}
\value{
Returns the following:          
\item{coefficients}{Coefficients of the model.}
\item{residuals}{ Residuals from the fitted model.}
\item{rho}{Unpenalized portion of the objective function.}
\item{tau}{ Quantile being modelled.}
\item{n}{Sample size.}
\item{intercept}{Whether intercept was included in model.}
\item{penalty}{Penalty used for fitting the model.}
\item{class}{rqPen and rqNC}
}
\description{
Linear programming implementation of quantile regression with a group penalty.  
}
\examples{                           
\dontrun{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- rq.group.lin.prog(x,y,groups=c(rep(1,4),rep(2,4)), tau=.5, lambda=1)
}
}
\author{Ben Sherwood}