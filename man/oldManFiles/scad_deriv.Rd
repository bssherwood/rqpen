\name{scad_deriv}
\alias{scad_deriv}
\title{SCAD Derivative}
\usage{
scad_deriv(x,lambda,a)
}
\arguments{
  \item{x}{Number to be evaluated}
  \item{lambda}{Tuning parameter lambda}
  \item{a}{Tuning parameter a. Default value of 3.7 as suggested in Fan and Li (2001)}
}
\value{
Derivative SCAD penalty function with tuning parameters lambda and "a" evaluated at "x".
}
\description{
Derivative of SCAD penalty function as described in Fan and Li (2001).
}
\examples{
scad_deriv(3,1)
scad_deriv(-3,1)
scad_deriv(.001,2)
}
\references{
[1] Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized likelihood and its oracle properties, \emph{J. Am. Statist. Ass}, \bold{96}, 1348--1360.
}
\author{ Ben Sherwood}