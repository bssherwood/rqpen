\name{beta_plots}
\alias{beta_plots}
\title{Plots of Betas}
\usage{
beta_plots(model,voi=NULL,logLambda=TRUE,loi=NULL,...)
}
\arguments{
  \item{model}{"cv.rq.pen" object.}
  \item{voi}{Index of betas to include. Default is all of the lambdas from "cv.rq.pen" object.}
  \item{logLambda}{Plot of lambdas is on the log scale.}
  \item{loi}{Index of lambdas to include. Default is all of the lambdas from "cv.rq.pen" object.}
  \item{...}{Additional arguments to be sent to plot.}
}
\value{
Plot of how beta estimates change with lambda.
}
\description{
Plots how the beta estimates changes with the different values of lambda.  
}
\examples{
set.seed(1)
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
lassoModels <- cv.rq.pen(x,y)
b_plot <- beta_plots(lassoModels)
}
\author{Ben Sherwood}