library(rqPen)
library(hqreg)

p <- 80
n <- 100

x <- matrix(rnorm(n*p),ncol=p)

y <- 1 + x[,1] + x[,3] - x[,60] + rnorm(n)

r1 <- rq.lasso(x,y)
h1 <- hqreg(x,y,method="quantile")
g1 <- glmnet(x,y)
r2 <- rq.lasso(x,y,tau=seq(.1,.9,.1))
r2 <- rq.lasso(x,y,tau=seq(.1,.9,.1),tau.pen)
