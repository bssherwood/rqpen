\name{predict.cv.rq.pen}
\alias{predict.cv.rq.pen}
\title{Prediction from a cv quantile regression penalized model}
\usage{
\method{predict}{cv.rq.pen}(object, newx, lambda,...)
}
\arguments{
  \item{object}{ "cv.rq.pen" object.} 
  \item{newx}{ Matrix of covariates used for prediction.}
  \item{lambda}{ Lambda associated with the model from which predictions should be made. Default is to use the lambda that provides the minimum criteria (cross-validation or BIC) that was selected by cv.rq.pen.}
  \item{...}{Needed for consistency with generic predict.}
}
\value{
Returns predicted values from the model for the selected lambda.   
}
\description{
Returns predicted values from "rqPen" model associated with lambda for "newx" covariates.    
}
\examples{
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- cv.rq.pen(x,y)
n2 <- matrix(rnorm(80),nrow=10)
preds <- predict(cv_model,n2)
}
\author{Ben Sherwood}