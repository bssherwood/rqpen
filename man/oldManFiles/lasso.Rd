\name{lasso}
\alias{lasso}
\title{Lasso}
\usage{
lasso(x,lambda)
}
\arguments{
  \item{x}{Number to be evaluated}
  \item{lambda}{Tuning parameter lambda}
}
\value{
lambda*abs(x)
}
\description{
LASSO penalty function.
}
\examples{
lasso(3,1)
lasso(-3,1)
lasso(-3,2)
}
\references{
[1] Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
\emph{Journal of the Royal Statistical Society. Series B}, \bold{58}, 267--288.
}
\author{ Ben Sherwood}