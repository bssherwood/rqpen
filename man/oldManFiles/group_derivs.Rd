\name{group_derivs}
\alias{group_derivs}
\title{Derivative of a group penalty}
\usage{
group_derivs(deriv_func,groups,coefs,lambda,a=3.7)
}
\arguments{
  \item{deriv_func}{Target derivative function.}
  \item{groups}{ Vector assigning columns of x to groups.}
  \item{coefs}{ Coefficients.}            
  \item{lambda}{ Lambda value for deriv_func.}   
  \item{a}{Additional parameter for deriv_func.}   
}
\value{
Returns the derivative of the L1 group penalty function.
}
\description{
Used to estimate non-convex group penalties. 
}
\author{Ben Sherwood}