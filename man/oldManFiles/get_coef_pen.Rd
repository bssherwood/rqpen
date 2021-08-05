\name{get_coef_pen}
\alias{get_coef_pen}
\title{Returns the coefficient part of the penalized objective function}
\usage{
get_coef_pen(coefs,lambda,intercept,penVars,penalty="LASSO")
}
\arguments{
  \item{coefs}{Coefficients to be transformed.}
  \item{lambda}{The penalty tuning parameter.}
  \item{intercept}{If the intercept is part of the coefficients.}
  \item{penVars}{Variables that were penalized.}
  \item{penalty}{The penalty function used.}
}
\description{
Returns the value for the penalized portion of the penalized objective function. 
}
\author{Ben Sherwood}