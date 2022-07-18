rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas)


library(hqreg)
library(glmnet)

set.seed(1)

x <- matrix(rnorm(100*8,sd=10),ncol=8)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(100,3)
g <- c(1,1,1,1,2,2,3,3)
tvals <- c(.25,.75)

m1 <- rq.pen.cv(x,y)

coefficients(m1$fit)
coefficients(m1$fit,lambdaIndex=c(3,5,7))
coefficients(m1$fit,tau=.5)
coefficients(m1$fit,a=1)
coefficients(m1$fit,a=1,tau=.5)
coefficients(m1$fit,tau=.5,lambdaIndex=c(3,5,7))
coefficients(m1$fit,lambda=m1$fit$models[[1]]$lambda[c(3,5,23)])


m2 <- rq.pen.cv(x,y,alg="br")
m3 <- rq.pen.cv(x,y,penalty="Ridge")
m4 <- rq.pen.cv(x,y,penalty="ENet",a=.5)
m5 <- rq.pen.cv(x,y,penalty="aLASSO")
m6 <- rq.pen.cv(x,y,penalty="SCAD")
m7 <- rq.pen.cv(x,y,penalty="MCP")

m1 <- rq.pen.cv(x,y,tau=c(.1,.3,.7))
coefficients(m1$fit)
coefficients(m1$fit,lambdaIndex=c(3,5,7))
coefficients(m1$fit,tau=.3)
coefficients(m1$fit,a=1)
coefficients(m1$fit,a=1,tau=.3)
coefficients(m1$fit,tau=.3,lambdaIndex=c(3,5,7))
coefficients(m1$fit,lambda=m1$fit$models[[1]]$lambda[c(3,5,23)])


m2 <- rq.pen.cv(x,y,alg="br",tau=c(.1,.3,.7))
m3 <- rq.pen.cv(x,y,penalty="Ridge",tau=c(.1,.3,.7))
m4 <- rq.pen.cv(x,y,penalty="ENet",a=.5,tau=c(.1,.3,.7))
m5 <- rq.pen.cv(x,y,penalty="aLASSO",tau=c(.1,.3,.7))
m6 <- rq.pen.cv(x,y,penalty="SCAD",tau=c(.1,.3,.7))
m7 <- rq.pen.cv(x,y,penalty="MCP",tau=c(.1,.3,.7),a=c(3,4,5))

m1 <- rq.group.pen.cv(x,y,tau=c(.1,.3,.7),groups=g)
m2 <- rq.group.pen.cv(x,y,alg="br",tau=c(.1,.3,.7),groups=g)
m5 <- rq.group.pen.cv(x,y,penalty="gAdLASSO",tau=c(.1,.3,.7),a=c(1,2,3),groups=g)
m6 <- rq.group.pen.cv(x,y,penalty="gSCAD",tau=c(.1,.3,.7),a=c(3,4,5),groups=g)
m7 <- rq.group.pen.cv(x,y,penalty="gMCP",tau=c(.1,.3,.7),a=c(3,4,5),groups=g)
m8 <- rq.group.pen.cv(x,y,penalty="gMCP",tau=c(.1,.3,.7),a=c(3,4,5),groups=g,norm=1)
m9 <- rq.group.pen.cv(x,y,penalty="gMCP",tau=c(.1,.3,.7),a=c(3,4,5),groups=g,norm=1)








