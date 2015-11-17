\name{scad}
\alias{scad}
\title{scad}
\usage{
scad(x,lambda,a)
}
\arguments{
  \item{x}{Number to be evaluated}
  \item{lambda}{Tuning parameter lambda}
  \item{a}{Tuning parameter a}
}
\value{
SCAD penalty function with tuning parameters lambda and "a" evaluated at "x".
}
\description{
SCAD penalty function as described in Fan and Li (2001).
}
\examples{
scad(3,1)
scad(-3,1)
scad(.001,2)
}
\references{
[1] Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties, \emph{J. Am. Statist. Ass}, \bold{96}, 1348--1360.
}
\author{ Ben Sherwood }