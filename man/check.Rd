\name{check}
\alias{check}
\title{Quantile check function}
\usage{
check(x, tau)
}
\arguments{
  \item{x}{Number to be evaluated.}
  \item{tau}{Number between 0 and 1 for quantile of interest.}
}
\value{
x*(tau-I(x < 0)
}
\description{
Evaluates the check function for quantile tau at value x. Check function is the objective function defined
in Koenker and Bassett (1978).
}
\examples{
check(2,.5)
check(-2,.5)
check(2,.2)
check(2,.8)
}
\references{
[1] Koenker, R. and Bassett, G. (1978). Regression Quantiles, \emph{Econometrica}, \bold{46}, 33--50. 
}
\author{Ben Sherwood}