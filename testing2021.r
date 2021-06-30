library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas)


library(hqreg)
library(glmnet)

set.seed(1)
p <- 8
n <- 100

x <- matrix(rnorm(n*p,sd=10),ncol=p)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(n,3)
g <- c(1,1,1,1,2,2,3,3)

obj   <- rq.nc(x,y,tau=.25, penalty="aLasso")


obj9 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=.25,penalty="gAdLasso")


obj9 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=.25)
obj10 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=c(.25,.75))
#coefficients(obj9)
#coefficients(obj10)

obj9 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=.25,penalty="gAdLasso")
obj10 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=c(.25,.75), penalty="gAdLasso")
#coefficients(obj9)
#coefficients(obj10)

obj9 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=.25,penalty="gSCAD")
obj10 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=c(.25,.75), penalty="gSCAD")
#coefficients(obj9)
#coefficients(obj10)

obj9 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=.25,penalty="gMCP")
obj10 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=c(.25,.75), penalty="gMCP")
#coefficients(obj9)
#coefficients(obj10)



obj10 <- hrq_glasso(x,y,c(1,1,1,1,2,2,3,3),tau=.25,w.lambda=c(1,1,1))

obj11 <- rq.group.pen(x,y,groups=g,penalty="gSCAD")


obj   <- rq.nc(x,y,tau=.25, penalty="aLasso", alg="lp",scalex=TRUE)
objns <- rq.nc(x,y,tau=.25, penalty="aLasso", alg="lp",scalex=FALSE)


obj6 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,3,3,3), penalty="gSCAD", norm=1)
obj7 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,3,3,3), penalty="gAdLasso", norm=1)
obj8 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,3,3,3), penalty="gMCP", norm=1,tau=c(.25,.75))
obj9 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=.25)
obj9 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3),tau=c(.25,.8))


hrq_glasso(x,y,group.index=c(1,1,1,1,2,2,2,3),tau=c(.25,.75))


# select debugging
obj <- rq.enet(x,y,tau=.475)
obj2 <- rq.enet(x,y,tau=c(.1,.7))

obj3 <- rq.nc(x,y,penalty="SCAD",tau=.4)
obj4 <- rq.nc(x,y,penalty="SCAD",tau=c(.5,.9))
#obj5 <- rq.nc(x,y,penalty="SCAD",alg="QICD")

obj6 <- rq.group.pen(x,y,groups=c(1,1,1,1,2,2,3,3))

obj <- rq.lasso(x,y,tau=.25)


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
