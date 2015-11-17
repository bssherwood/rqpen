\name{groupMultLambda}
\alias{groupMultLambda}
\title{Quantile Regression with Group Penalty for multiple lambdas}
\usage{
groupMultLambda(x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
	penalty="LASSO", alg="QICD", ...) 
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{groups}{ Vector assigning columns of x to groups.}
  \item{tau}{ Conditional quantile being modelled.}            
  \item{lambda}{ Vector of lambdas. Default is for lambdas to be automatically 
					generated.}   
  \item{intercept}{ Whether model should include an intercept. Constant does not 
					need to be included in "x".}
  \item{penalty}{ Type of penalty: "LASSO", "SCAD" or "MCP".}   
  \item{alg}{"QICD" for QICD implementation. Otherwise linear programming approach is implemented.}
  \item{...}{Additional parameters to be sent to rq.group.fit.}
}
\value{
Returns a list of rq.group.pen objects. Each element of the list is a fit for a 
			different value of lambda. 
}
\description{
Fit multiple models with L1 group penalty. QICD algorithm is using an adaptation of the algorithm 
		presented by Peng and Wang (2015). 
}
\examples{
\dontrun{
x <- matrix(rnorm(400),nrow=100)
y <- 1 + x[,1] - 3*x[,3] + rnorm(100)
cv_model <- groupMultLambda(x,y,groups=c(rep(1,2),rep(2,2)),lambda=seq(.1,.5,.1))
}
}
\references{
[1] Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with 
	grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.

[2] Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for 
	High-Dimensional Nonconvex Penalized Quantile Regression. 
	\emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.

}
\author{Ben Sherwood}