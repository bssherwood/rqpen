\name{groupQICD}
\alias{groupQICD}
\title{Quantile Regression with Group Penalty using QICD algorithm}
\usage{
groupQICD(x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                         max_iter = 100, eps = 1e-05, penalty = "SCAD", 
                         a = 3.7, coef.cutoff = 1e-08, 
                         initial_beta = NULL) 
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{groups}{ Vector assigning columns of x to groups.}
  \item{tau}{ Conditional quantile being modelled.}            
  \item{lambda}{ Vector of lambdas. Default is for lambdas to be 
				automatically generated.}   
  \item{intercept}{ Whether model should include an intercept. Constant 
					does not need to be included in "x".}
  \item{max_iter}{Maximum number of iterations in QICD algorithm. For both 
					the inner and outer loop of the algorithm}
  \item{eps}{Squared difference estimates below this value are considered 
					to have converged.}
  \item{penalty}{ Type of penalty: "LASSO", "SCAD" or "MCP".}   
  \item{a}{Additional tuning parameter for SCAD and MCP}  
  \item{coef.cutoff}{Coefficients below this value are set to zero.}
  \item{initial_beta}{Initial value of beta to start algorithm. Default is all zeros.}
}
\value{
Returns coefficients of the model. This is mostly used as a worker function for 
		other functions such as rq.group.pen. 
}
\description{
Fit the model with L1 group penalty using an adaptation of the QICD algorithm presented 
		by Peng and Wang (2015). 
}
\examples{
x <- matrix(rnorm(400),nrow=100)
y <- 1 + x[,1] - 3*x[,3] + rnorm(100)
cv_model <- groupQICD(x,y,groups=c(rep(1,2),rep(2,2)),lambda=.5)
}
\references{
[1] Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with 
	grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.

[2] Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for 
    High-Dimensional Nonconvex Penalized Quantile Regression. 
	\emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.

}
\author{Ben Sherwood}