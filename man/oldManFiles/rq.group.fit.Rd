\name{rq.group.fit}
\alias{rq.group.fit}
\title{Quantile Regresion with Group Penalty}
\usage{
rq.group.fit(x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                penalty = "SCAD", alg = "QICD", a=3.7, penGroups=NULL, ...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{groups}{ Vector assigning columns of x to groups.}
  \item{tau}{ Conditional quantile being modelled.}            
  \item{lambda}{ Vector of lambdas. Default is for lambdas to be automatically generated.}   
  \item{intercept}{ Whether model should include an intercept. Constant does not need to be included in "x".}
  \item{penalty}{ Type of penalty: "SCAD" or "MCP".}     
  \item{alg}{If set to QICD algorithm will use coordinate descent algorithm. Otherwise, will use a linear programming algorithm, which is an extension of algorithm used in rq.lasso.fit.}
  \item{a}{The additional tuning parameter for SCAD and MCP.}
  \item{penGroups}{Specify which groups will be penalized. Default is to penalize all groups.}
  \item{...}{Additional arguments to be sent to rq.lasso.fit or groupQICD.}   
}
\value{
Returns the following:          
\item{coefficients}{Coefficients of the model.}
\item{residuals}{ Residuals from the fitted model.}
\item{rho}{Unpenalized portion of the objective function.}
\item{tau}{ Quantile being modeled.}
\item{n}{Sample size.}
\item{intercept}{Whether intercept was included in model.}
}
\description{
Similar to cv.rq.pen function, but uses group penalty. Group penalties use the L1 norm instead of L2 for computational convenience. As a result of this the group lasso penalty is the same as the typical lasso penalty and thus you should only use a SCAD or MCP penalty. Only the SCAD and MCP penalties incorporate the group structure into the penalty. The group lasso penalty is implemented because it is needed for the SCAD and MCP algorithm. We use a group penalty extension of the QICD algorithm presented by Peng and Wang (2015). 
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- rq.group.fit(x,y,groups=c(rep(1,4),rep(2,4)),lambda=1,penalty="SCAD")
}
\references{
[1] Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.

[2] Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.

}
\author{Ben Sherwood; Adam Maidman}