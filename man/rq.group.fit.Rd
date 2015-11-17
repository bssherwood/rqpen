\name{rq.group.fit}
\alias{rq.group.fit}
\title{Quantile Regresion with Group Penalty}
\usage{
rq.group.fit(x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                penalty = "LASSO", alg = "QICD", ...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{groups}{ Vector assigning columns of x to groups.}
  \item{tau}{ Conditional quantile being modelled.}            
  \item{lambda}{ Vector of lambdas. Default is for lambdas to be automatically generated.}   
  \item{intercept}{ Whether model should include an intercept. Constant does not need to be included in "x".}
  \item{penalty}{ Type of penalty: "LASSO", "SCAD" or "MCP".}     
  \item{alg}{If set to QICD algorithm will use coordinate descent algorithm. Otherwise, will use a linear programming algorithm, which is an extension of algorithm used in rq.lasso.fit.}
  \item{...}{Additional arguments to be sent to rq.lasso.fit or groupQICD.}   
}
\value{
Returns the following:          
\item{coefficients}{Coefficients of the model.}
\item{residuals}{ Residuals from the fitted model.}
\item{rho}{Unpenalized portion of the objective function.}
\item{tau}{ Quantile being modelled.}
\item{n}{Sample size.}
\item{intercept}{Whether intercept was included in model.}
}
\description{
Similar to cv.rq.pen function, but uses group penalty. Group penalties use the L1 norm instead of L2 for computational convenience. We use a group penalty extension of the QICD algorithm presented by Peng and Wang (2015). 
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- rq.group.fit(x,y,groups=c(rep(1,4),rep(2,4)),lambda=1)
}
\references{
[1] Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.

[2] Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.

}
\author{Ben Sherwood}