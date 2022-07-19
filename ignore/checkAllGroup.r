rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
3
library(rqPen)

set.seed(1)

x <- matrix(rnorm(25*30,sd=10),ncol=30)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(25,3)
g <- rep(seq(1:5),6)
tvals <- c(.25,.75)

#run examples
r1 <- rq.group.pen(x,y,groups=g)
r5 <- rq.group.pen(x,y,groups=g,tau=tvals)
#Linear programming approach with group SCAD penalty and L1-norm
m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gSCAD",norm=1,a=seq(3,4))
# No penalty for the first group
m3 <- rq.group.pen(x,y,groups=g,group.pen.factor=c(0,rep(1,4)))
# Smaller penalty for the median
m4 <- rq.group.pen(x,y,groups=g,tau=c(.25,.5,.75),tau.penalty.factor=c(1,.25,1))

h1 <- hrq_glasso(x,y,g)

r1 <- rq.group.pen(x,y,groups=g)
r2 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO")
r3 <- rq.group.pen(x,y,groups=g,penalty="gSCAD") #still a problem here with the updates of coefficients, I think. 
r4 <- rq.group.pen(x,y,groups=g,penalty="gMCP")
r5 <- rq.group.pen(x,y,groups=g,tau=tvals)

r1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals)
r2 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals)
r3 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals)

r1 <- rq.group.pen(x,y,groups=g,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gAdLASSO",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gSCAD",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gMCP",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",norm=1)

r1 <- rq.group.pen(x,y,groups=g,tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",tau=tvals,norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gAdLASSO",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",tau=tvals,norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gSCAD",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",tau=tvals,norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gMCP",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",tau=tvals,norm=1)

rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas)


library(hqreg)
library(glmnet)

set.seed(1)

x <- matrix(rnorm(25*30,sd=10),ncol=30)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(25,3)
g <- rep(seq(1:5),6)
tvals <- c(.25,.75)

rq.group.pen(x,y,groups=g,a=seq(3,5))

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,norm=1,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gAdLASSO",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gSCAD",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gMCP",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="br",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

