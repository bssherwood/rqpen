\name{mcp}
\alias{mcp}
\title{MCP}
\usage{
mcp(x,lambda,a)
}
\arguments{
  \item{x}{Number to be evaluated}
  \item{lambda}{Tuning parameter lambda}
  \item{a}{Tuning parameter a}
}
\value{
MCP function with tuning parameters lambda and "a" evaluated at "x".
}
\description{
MCP function as described in Fan and Li (2001).
}
\examples{
mcp(3,1)
mcp(-3,1)
mcp(.001,2)
}
\references{
[1] Zhang, C. Nearly unbiased variable selection under minimax concave penalty. (2010). \emph{Ann. Statist.}, \bold{38}, 894--942.
}
\author{ Ben Sherwood }