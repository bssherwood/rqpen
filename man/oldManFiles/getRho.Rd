\name{getRho}
\alias{getRho}
\title{Objective Function Value}
\usage{
getRho(model)
}
\arguments{
  \item{model}{Object with rho as an attribute}
}
\value{
Rho is the value of the unpenalized portion of the objective function.
}
\description{
Returns unpenalized portion of the objective function for a penalized quantile regression object.   
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
lassoModel <- rq.lasso.fit(x,y,lambda=.5)
getRho(lassoModel)
}
\author{Ben Sherwood}