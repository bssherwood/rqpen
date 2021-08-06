rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas) #blah


library(hqreg)
library(glmnet)

set.seed(1)

x <- matrix(rnorm(100*8,sd=10),ncol=8)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(100,3)
g <- c(1,1,1,1,2,2,3,3)
tvals <- c(.25,.75)



m1 <- rq.pen(x,y,alg="huber")
m2 <- rq.pen(x,y,alg="br")
m1$models
m2$models

m2 <- rq.pen(x,y,penalty="ENet",a=.5)
m2$models

m1 <- rq.pen(x,y,penalty="SCAD")
m2 <- rq.pen(x,y,penalty="SCAD",alg="br")
m3 <- rq.pen(x,y,penalty="SCAD",alg="QICD")
m1$models
m2$models
m3$models

m1 <- rq.pen(x,y,penalty="aLASSO")
m1$models

m1 <- rq.pen(x,y,penalty="MCP")
m2 <- rq.pen(x,y,alg="br", penalty="MCP")
m3 <- rq.pen(x,y,alg="QICD",penalty="MCP")
m1$models
m2$models
m3$models

m1 <- rq.pen(x,y,tau=tvals)
m2 <- rq.pen(x,y,alg="br",tau=tvals)
m1$models
m2$models

m2 <- rq.pen(x,y,penalty="ENet",tau=tvals,a=.5)
m2$models

m1 <- rq.pen(x,y,penalty="SCAD",tau=tvals)
m2 <- rq.pen(x,y,penalty="SCAD",alg="br",tau=tvals)
m3 <- rq.pen(x,y,penalty="SCAD",alg="QICD",tau=tvals)
m1$models
m2$models
m3$models

m1 <- rq.pen(x,y,tau=tvals,penalty="aLASSO")
m1$models

m1 <- rq.pen(x,y,tau=tvals,penalty="MCP")
m2 <- rq.pen(x,y,alg="br", penalty="MCP",tau=tvals)
m3 <- rq.pen(x,y,alg="QICD",penalty="MCP",tau=tvals)
m1$models
m2$models
m3$models

m1 <- rq.pen(x,y,a=c(.2,.5,.7),penalty="ENet")
m1$models

m1 <- rq.pen(x,y,a=c(3,4,5),penalty="SCAD")
m2 <- rq.pen(x,y,alg="br", a = c(3,4,5),penalty="SCAD")
m3 <- rq.pen(x,y,alg="QICD", a = c(3,4,5),penalty="SCAD")
m1$models
m2$models
m3$models 

m1 <- rq.pen(x,y,penalty="aLASSO", a=c(1,2,3))
m1$models

m1 <- rq.pen(x,y,penalty="MCP", a=c(3,4,5))
m2 <- rq.pen(x,y,alg="br", penalty="MCP",a=c(3,4,5))
m3 <- rq.pen(x,y,alg="QICD",penalty="MCP",a=c(3,4,5))
m1$models
m2$models
m3$models

m1 <- rq.pen(x,y,a=c(.2,.5,.7),tau=tvals,penalty="ENet")
m1$models

m1 <- rq.pen(x,y,tau=tvals,a=c(3,4,5),penalty="SCAD")
m2 <- rq.pen(x,y,alg="br",tau=tvals,a=c(3,4,5),penalty="SCAD")
m3 <- rq.pen(x,y,alg="QICD",tau=tvals,a=c(3,4,5),penalty="SCAD")
m1$models
m2$models
m3$models

m1 <- rq.pen(x,y,tau=tvals,penalty="aLASSO",a=c(1,2,3))
m1$models

m1 <- rq.pen(x,y,tau=tvals,penalty="MCP",a=c(3,4,5))
m2 <- rq.pen(x,y,alg="br", penalty="MCP",tau=tvals,a=c(3,4,5))
m3 <- rq.pen(x,y,alg="QICD",penalty="MCP",tau=tvals,a=c(3,4,5))
m1$models
m2$models
m3$models







