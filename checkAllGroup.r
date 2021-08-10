rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
3
library(rqPen)
library(hrqglas)


library(hqreg)
library(glmnet)

set.seed(1)

x <- matrix(rnorm(25*30,sd=10),ncol=30)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(25,3)
g <- rep(seq(1:5),6)
tvals <- c(.25,.75)

h1 <- hrq_glasso(x,y,g)

r1 <- rq.group.pen(x,y,groups=g)
r2 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO")
r3 <- rq.group.pen(x,y,groups=g,penalty="gSCAD")
r4 <- rq.group.pen(x,y,groups=g,penalty="gMCP")
r5 <- rq.group.pen(x,y,groups=g,tau=tvals)

rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals)
rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals)
rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals)

rq.group.pen(x,y,groups=g,norm=1)
rq.group.pen(x,y,groups=g,alg="lp",norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",norm=1)

rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",norm=1)

rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",norm=1)

rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",norm=1)

rq.group.pen(x,y,groups=g,tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="lp",tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",tau=tvals,norm=1)

rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",tau=tvals,norm=1)

rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",tau=tvals,norm=1)

rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",tau=tvals,norm=1)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",tau=tvals,norm=1)

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
m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",a=seq(3,5))
m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",a=seq(3,5))
m1$models

m1 <- rq.group.pen(x,y,groups=g,tau=tvals,a=seq(3,5))
m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,a=seq(3,5))
m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,a=seq(3,5))
m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,a=seq(3,5))
m1$models

m1 <- rq.group.pen(x,y,groups=g,norm=1,a=seq(3,5))
m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

