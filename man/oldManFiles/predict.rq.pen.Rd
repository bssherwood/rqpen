\name{predict.rq.pen}
\alias{predict.rq.pen}
\title{Prediction from a quantile regression penalized model}
\usage{
\method{predict}{rq.pen}(object, newx,...)
}
\arguments{
  \item{object}{ rq.pen object.}
  \item{newx}{ Matrix of covariates used for prediction.}
  \item{...}{Needed for consistency with generic predict.}
}
\value{
Returns predicted values from the model.   
}
\description{
Returns predicted values from "rq.pen" object for "newx" covariates.    
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
lassoModel <- rq.lasso.fit(x,y,lambda=1)
n2 <- matrix(rnorm(80),nrow=10)
preds <- predict(lassoModel,n2)
}
\author{Ben Sherwood}