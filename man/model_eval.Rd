\name{model_eval}
\alias{model_eval}
\title{Model Evaluation}
\usage{
model_eval(model, test_x, test_y, test_w=NULL, func="check",...)
}
\arguments{
  \item{model}{ Model of class "rqPen".}
  \item{test_x}{ Covariates used for prediction.}
  \item{test_y}{ Response to compare predictions against.}
  \item{test_w}{ Weights for a weighted mean, typically used if weights were used to fit the model.}
  \item{func}{ Function used for evaluation. Options: "check" (Quantile Check), "SqErr" (Squared Error), "AE" (Absolute Value)}
  \item{...}{Additional arguments to be sent to evaluation function. For instance check requires tau which is not part of the model_eval function.}
}
\value{
Mean of prediction errors using the chosen function. 
}
\description{
Used for cross-validation. For a model of class "rqPen" it provides the average prediction error given the evaluation function of choice.
}
\examples{
x <- matrix(rnorm(800),ncol=8)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
l_model <- rq.lasso.fit(x,y,lambda=1)
nc_model <- rq.nc.fit(x,y,lambda=1)
newx <- matrix(rnorm(16),ncol=8)
newy <- 1 + newx[,1] - 3*newx[,5] + rnorm(2)
model_eval(l_model, newx, newy)
model_eval(l_model, newx, newy, func="SqErr")
model_eval(nc_model, newx, newy)
}
\author{Ben Sherwood}