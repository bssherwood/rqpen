\name{plot.cv.rq.group.pen}
\alias{plot.cv.rq.group.pen}
\title{Plot cv.rq.group.pen}
\usage{
\method{plot}{cv.rq.group.pen}(x,y=NULL,...) 
}
\arguments{
  \item{x}{cv.rq.group.pen object.}
  \item{y}{holder value to match up with default plot program}
  \item{...}{Additional values to function, but not currently used}
}
\value{
Plots the validation criteria against the lambda values.   
}
\description{
Plots the validation criteria against the lambda values.   
}
\examples{
\dontrun{
x <- matrix(rnorm(400),nrow=100)
y <- 1 + x[,1] - 3*x[,3] + rnorm(100)
cv_model <- cv.rq.group.pen(x,y,groups=c(rep(1,2),rep(2,2)))
plot(cv_model)
}
}
\author{Ben Sherwood}