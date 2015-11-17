\name{kernel_estimates}
\alias{kernel_estimates}
\title{Kernel based estimates of Y|X}
\usage{
kernel_estimates(x,y,h,...)
}
\arguments{
  \item{x}{matrix of predictors}
  \item{y}{Vector of response}
  \item{h}{Scalar bandwidth tuning parameter}
  \item{...}{Additional arguments to be sent to kernesti.regr from regrpo package.}
}
\value{
Estimates a conditional density. For future use of implementing inverse probability weights (IPW) to handle missing data. 
}
\description{
Provides fitted values of Y.  
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- rbinom(100,1,.5)
cond_fit <- kernel_estimates(x,y,1)
}
\author{Ben Sherwood}