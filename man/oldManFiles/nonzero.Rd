\name{nonzero}
\alias{nonzero}
\title{Nonzero}
\usage{
nonzero(obj) 
}
\arguments{
  \item{obj}{ Model.}
}
\value{
Returns if coefficients or groups are nonzero or not. TRUE if they are not zero and FALSE if they are.  
}
\description{
Calls nonzero function based on object. 
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