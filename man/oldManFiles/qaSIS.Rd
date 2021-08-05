\name{qaSIS}
\alias{qaSIS}
\title{Quantile Adaptive Sure Independence Screening}
\usage{
qaSIS(x,y,tau,linear=FALSE,...)
}
\arguments{
  \item{x}{ Matrix of predictors.}
  \item{y}{ Vector of response values.}
  \item{tau}{ Conditional quantile being modelled.}   
  \item{linear}{ If true linear screening will be done, otherwise splines will be fit to each column vector.}
  \item{...}{Additional items to be sent to bs function from splines package.}
}
\value{
Returns the ranking of important predictors, from highest to lowest.
}
\description{
Implements quantile adaptive screening as outlined by He, Wang and Hong
}
\examples{
x <- pnorm(matrix(rnorm(800),nrow=100))
y <- 1 + 2*cos(2*pi*x[,1])+exp(2*x[,2]) + rnorm(100,sd=.1)
var_ranks <- qaSIS(x,y)
}
\references{
[1] He, X., Wang, L. and Hong, H. (2013). Quantile-Adaptive Model-free Variable Screening for High-dimensional Heterogeneous Data. \emph{The Annals of Statistics}, \bold{41}, 342--369.  
}
\author{Ben Sherwood}