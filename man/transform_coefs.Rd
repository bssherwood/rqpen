\name{transform_coefs}
\alias{transform_coefs}
\title{Transform coefficients back to original scale}
\usage{
transform_coefs(coefs,mu_x,sigma_x)
}
\arguments{
  \item{coefs}{Coefficients to be transformed.}
  \item{mu_x}{Means of the original predictors.}
  \item{sigma_x}{Standard deviations of the original predictors.}
}
\description{
Takes the coefficients fit after scaling the predictors to have mean zero and a standard 
deviation of one and transforms them back to the original scale of the predictors. 
}
\author{Ben Sherwood}