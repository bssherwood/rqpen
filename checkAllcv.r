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
m2 <- rq.pen.cv(x,y,alg="br")
m3 <- rq.pen.cv(x,y,penalty="Ridge")
m4 <- rq.pen.cv(x,y,penalty="ENet",a=.5)
m5 <- rq.pen.cv(x,y,penalty="aLASSO")
m6 <- rq.pen.cv(x,y,penalty="SCAD")
m7 <- rq.pen.cv(x,y,penalty="MCP")

m1 <- rq.pen.cv(x,y,tau=c(.1,.3,.7))
m2 <- rq.pen.cv(x,y,alg="br",tau=c(.1,.3,.7))
m3 <- rq.pen.cv(x,y,penalty="Ridge",tau=c(.1,.3,.7))
m4 <- rq.pen.cv(x,y,penalty="ENet",a=.5,tau=c(.1,.3,.7))
m5 <- rq.pen.cv(x,y,penalty="aLASSO",tau=c(.1,.3,.7))
m6 <- rq.pen.cv(x,y,penalty="SCAD",tau=c(.1,.3,.7))
m7 <- rq.pen.cv(x,y,penalty="MCP",tau=c(.1,.3,.7))

m1 <- rq.group.pen.cv(x,y,tau=c(.1,.3,.7))
m2 <- rq.group.pen.cv(x,y,alg="br",tau=c(.1,.3,.7))
m3 <- rq.group.pen.cv(x,y,penalty="Ridge",tau=c(.1,.3,.7))
m4 <- rq.group.pen.cv(x,y,penalty="ENet",a=.5,tau=c(.1,.3,.7),a=seq(.1,.9,.1))
m5 <- rq.group.pen.cv(x,y,penalty="aLASSO",tau=c(.1,.3,.7),a=c(1,2,3))
m6 <- rq.group.pen.cv(x,y,penalty="SCAD",tau=c(.1,.3,.7),a=c(3,4,5))
m7 <- rq.group.pen.cv(x,y,penalty="MCP",tau=c(.1,.3,.7),a=c(3,4,5))








