library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)


library(hqreg)
library(glmnet)




p <- 8
n <- 100

x <- matrix(rnorm(n*p),ncol=p)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(n,353)

# select debugging
obj <- rq.enet(x,y,tau=.475)
obj2 <- rq.enet(x,y,tau=c(.1,.7))

obj3 <- rq.nc(x,y,penalty="SCAD",tau=.4)
obj4 <- rq.nc(x,y,penalty="SCAD",tau=c(.5,.9))
obj5 <- rq.nc(x,y,penalty="SCAD",alg="QICD")

obj6 <- rq.group.pen(x,y,group=c(1,1,1,1,2,2,3,3))


obj2 <- rq.lla(obj,x,y)
coefficients(qic.select(obj2))

obj  <- rq.lasso(x,y,alg="huber",tau=.475)
obj2 <- rq.lla(obj,x,y)
obj3 <- rq.lla(obj,x,y,penalty="MCP")


r1a <- rq.lasso(x,y,alg="br",tau=.475,penalty.factor=c(0,0,0,0,1,0,2,3))
r2 <- qic.select(r1,method="PBIC")
r2a <- qic.select(r1a,method="BIC",septau=TRUE)

x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
cv_model <- cv.rq.pen(x,y)
cv_model$models




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
