library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hqreg)
library(glmnet)

p <- 8
n <- 100

x <- matrix(rnorm(n*p),ncol=p)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(n,3)

# select debugging
r1 <- rq.lasso(x,y,alg="huber",tau=.475)
r1a <- rq.lasso(x,y,alg="br",tau=.475)
r2 <- qic.select(r1,method="PBIC")
r2a <- qic.select(r1a,method="BIC",septau=TRUE)






h1 <- hqreg(x,y,method="quantile")
g1 <- glmnet(x,y)

gfits <- cbind(1,x)%*%coefficients(g1)
gresids <- y-gfits

tLL <- g1$nulldev - deviance(g1)
k <- g1$df
n <- g1$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc

BIC<-log(n)*k - tLL
BIC


r2 <- rq.lasso(x,y,tau=seq(.1,.9,.1))
r2 <- rq.lasso(x,y,tau=seq(.1,.9,.1),tau.pen)
