\name{cv_plots}
\alias{cv_plots}
\title{Plots of Cross-validation results}
\usage{
cv_plots(model,logLambda=TRUE,loi=NULL,...)
}
\arguments{
  \item{model}{"cv.rq.pen" object.}
  \item{logLambda}{Plot of lambdas is on the log scale.}
  \item{loi}{Index of lambdas to be plotted. Default is all of the lambdas from "cv.rq.pen" object.}
  \item{...}{Additional items to be sent to plot function.}
}
\value{
Plot of how cross validation statistic changes with lambda.
}
\description{
Slightly misnamed as user could choose BIC as a criteria for "cv.rq.pen" object. Function is able to discern between the two 
types of evaluation criteria and provides appropriate labels for the plot.  
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
lassoModels <- cv.rq.pen(x,y)
cv_plot <- cv_plots(lassoModels)
}
\author{Ben Sherwood}