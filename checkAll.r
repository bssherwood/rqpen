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

rq.lasso(x,y)
rq.lasso(x,y,alg="br")

rq.enet(x,y)
rq.enet(x,y,a=.5)

rq.nc(x,y)
rq.nc(x,y,alg="br")
rq.nc(x,y,alg="QICD")

rq.nc(x,y,penalty="aLASSO")
rq.nc(x,y,alg="br", penalty="aLASSO")
rq.nc(x,y,alg="QICD",penalty="aLASSO")

rq.nc(x,y,penalty="MCP")
rq.nc(x,y,alg="br", penalty="MCP")
rq.nc(x,y,alg="QICD",penalty="MCP")

rq.lasso(x,y,tau=tvals)
rq.lasso(x,y,alg="br",tau=tvals)

rq.enet(x,y,tau=tvals)
rq.enet(x,y,a=.5,tau=tvals)

rq.nc(x,y,tau=tvals)
rq.nc(x,y,alg="br",tau=tvals)
rq.nc(x,y,alg="QICD",tau=tvals)

rq.nc(x,y,tau=tvals,penalty="aLASSO")
rq.nc(x,y,alg="br", penalty="aLASSO",tau=tvals)
rq.nc(x,y,alg="QICD",penalty="aLASSO",tau=tvals)

rq.nc(x,y,tau=tvals,penalty="MCP")
rq.nc(x,y,alg="br", penalty="MCP",tau=tvals)
rq.nc(x,y,alg="QICD",penalty="MCP",tau=tvals)

rq.enet(x,y,a=c(.2,.5,.7))

rq.nc(x,y,a=c(3,4,5))
rq.nc(x,y,alg="br", a = c(3,4,5))
rq.nc(x,y,alg="QICD", a = c(3,4,5))

rq.nc(x,y,penalty="aLASSO", a=c(1,2,3))
rq.nc(x,y,alg="br", penalty="aLASSO", a=c(1,2,3))
rq.nc(x,y,alg="QICD",penalty="aLASSO", a=c(1,2,3))

rq.nc(x,y,penalty="MCP", a=c(3,4,5))
rq.nc(x,y,alg="br", penalty="MCP",a=c(3,4,5))
rq.nc(x,y,alg="QICD",penalty="MCP",a=c(3,4,5))

rq.enet(x,y,a=c(.2,.5,.7),tau=tvals)

rq.nc(x,y,tau=tvals,a=c(3,4,5))
rq.nc(x,y,alg="br",tau=tvals,a=c(3,4,5))
rq.nc(x,y,alg="QICD",tau=tvals,a=c(3,4,5))

rq.nc(x,y,tau=tvals,penalty="aLASSO",a=c(1,2,3))
rq.nc(x,y,alg="br", penalty="aLASSO",tau=tvals,a=c(1,2,3))
rq.nc(x,y,alg="QICD",penalty="aLASSO",tau=tvals,a=c(1,2,3))

rq.nc(x,y,tau=tvals,penalty="MCP",a=c(3,4,5))
rq.nc(x,y,alg="br", penalty="MCP",tau=tvals,a=c(3,4,5))
rq.nc(x,y,alg="QICD",penalty="MCP",tau=tvals,a=c(3,4,5))








