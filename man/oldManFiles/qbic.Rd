\name{qbic}
\alias{qbic}
\title{Quantile Regresion BIC}
\usage{
qbic(model, largeP=FALSE)
}
\arguments{
  \item{model}{Model of class "rqPen".}
  \item{largeP}{Large P version using an additional penalty factor of log(s) where "s" is the total number of covariates considered.}  
}
\value{
Numeric value representing BIC of selected model.  
}
\description{
Quantile regression BIC with large p alternative as described in Lee, Noh and Park (2013). 
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
l_model <- rq.lasso.fit(x,y, lambda=1)
nc_model <- rq.nc.fit(x,y, lambda=1)
qbic(l_model)
qbic(nc_model)
qbic(l_model, largeP=TRUE)
qbic(nc_model, largeP=TRUE)
}
\references{
[1] Lee, E., Noh, H. and Park, B. (2014). Model selection via Bayesian Information Criterion for quantile regression models., \emph{J. Am. Statist. Ass}, \bold{109}, 216--229.
}
\author{Ben Sherwood}