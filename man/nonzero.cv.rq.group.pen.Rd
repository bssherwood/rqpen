\name{nonzero.cv.rq.group.pen}
\alias{nonzero.cv.rq.group.pen}
\title{Nonzero}
\usage{
\method{nonzero}{cv.rq.group.pen}(obj) 
}
\arguments{
  \item{obj}{cv.rq.group.pen object.}
}
\value{
Returns true if all elements in a group are non-zero and FALSE if they are not. Chooses the model associated with lambda.min.  
}
\description{
Calls nonzero function. 
}
\examples{
\dontrun{
x <- matrix(rnorm(400),nrow=100)
y <- 1 + x[,1] - 3*x[,3] + rnorm(100)
cv_model <- cv.rq.group.pen(x,y,groups=c(rep(1,2),rep(2,2)))
nonzero(cv_model)
}
}
\author{Ben Sherwood}