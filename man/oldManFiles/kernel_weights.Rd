\name{kernel_weights}
\alias{kernel_weights}
\title{Nonparametric estimate of IPW weights}
\usage{
kernel_weights(obs_data,obs_ind,...)
}
\arguments{
  \item{obs_data}{Matrix of variables with complete observations}
  \item{obs_ind}{Vector of whether sample is observed or not (1-observed, 0-not)}
  \item{...}{Additional arguments to be sent to kernel_estimates.}
}
\value{
Estimates of weights. 
}
\description{
This is for downstream analysis for fitting models with missing data. Future work is to fully incorporate these into penalized models. Tuning parameter for condtional density is esimated using approach of Chen, Wan and Zhou (2015), which is a simplified approach of Sepanski et al. (1994)
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- rbinom(100,1,.5)
wt_est <- kernel_weights(x,y)
}
\references{
[1] Chen, X., Wan, A. and Zhou, Y. Efficient quantile regression analysis with missing observations. (2015).  \emph{J. Amer. Statist. Assoc.}, \bold{110}, 723--741.
[2] Sepanski, J., Knickerbocker, R. and Carroll, R. A semiparametric correction for attenuation. (1994). \emph{J. Amer. Statist. Assoc.}, \bold{89}, 1366--1373.
}
\author{Ben Sherwood}