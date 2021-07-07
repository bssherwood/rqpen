library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas)


library(hqreg)
library(glmnet)

set.seed(1)
p <- 30
n <- 25

x <- matrix(rnorm(n*p,sd=10),ncol=p)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(n,3)
g <- c(1,1,1,1,2,2,3,3)
tvals <- c(.25,.75)

rq.group.pen(x,y,groups=g)
rq.group.pen(x,y,groups=g,alg="lp")
rq.group.pen(x,y,groups=g,alg="qicd")

rq.group.pen(x,y,groups=g,penalty="gAdLASSO")
rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO")
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO")

rq.group.pen(x,y,groups=g,penalty="gSCAD")
rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD")
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD")

rq.group.pen(x,y,groups=g,penalty="gMCP")
rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP")
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP")

rq.group.pen(x,y,groups=g,tau=tvals)
rq.group.pen(x,y,groups=g,alg="lp",tau=tvals)
rq.group.pen(x,y,groups=g,alg="qicd",tau=tvals)

rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",tau=tvals)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",tau=tvals)

rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",tau=tvals)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",tau=tvals)

rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals)
rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",tau=tvals)
rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",tau=tvals)


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

