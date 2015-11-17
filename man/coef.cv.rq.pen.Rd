\name{coef.cv.rq.pen}
\alias{coef.cv.rq.pen}
\title{Penalized Quantile Regression Coefficients}
\usage{
\method{coef}{cv.rq.pen}(object, lambda='min',...)
}
\arguments{
  \item{object}{cv.rq.pen object}
  \item{lambda}{Tuning parameter lambda. Default is to select minimum lambda from cross-validation method. User can also select a specific value of lambda, but it needs to be a lambda that was part of the fit of cv.rq.pen object.}
  \item{...}{Additional arguments, currently not used.}
}
\value{
Coefficients for selected value of lambda. 
}
\description{
Returns coefficients for a cv.rq.pen object with default values being coefficients associated with minimum cross-validation value. 
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- cv.rq.pen(x,y)
coefficients(cv_model)
}
\author{ Ben Sherwood }