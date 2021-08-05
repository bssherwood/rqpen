\name{randomly_assign}
\alias{randomly_assign}
\title{Randomly Assign}
\usage{
randomly_assign(n,k)
}
\arguments{
  \item{n}{ Number of samples.}
  \item{k}{ Number of groups.}
}
\value{
A vector of length n having entries of 1,...k.
}
\description{
Randomly assign n samples into k groups
}
\examples{
randomly_assign(37,5)
randomly_assign(11,3)
}
\author{Ben Sherwood}