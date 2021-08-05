\name{coef.cv.rq.group.pen}
\alias{coef.cv.rq.group.pen}
\title{Group Penalized Quantile Regression Coefficients}
\usage{
\method{coef}{cv.rq.group.pen}(object, lambda='min',...)
}
\arguments{
  \item{object}{cv.rq.group.pen object}
  \item{lambda}{Tuning parameter lambda. Default is to select minimum lambda from cross-validation method. User can also select a specific value of lambda, but it needs to be a lambda that was part of the fit of cv.rq.pen object.}
  \item{...}{Additional arguments, but currently not used.}
}
\value{
Coefficients for selected value of lambda. 
}
\description{
Returns coefficients for a cv.rq.pen object with default values being coefficients associated with minimum cross-validation value. 
}
\examples{
\dontrun{
x <- matrix(rnorm(400),nrow=100)
y <- 1 + x[,1] - 3*x[,3] + rnorm(100)
cv_model <- cv.rq.group.pen(x,y,groups=c(rep(1,2),rep(2,2)))
coefficients(cv_model)
}
}
\author{ Ben Sherwood }