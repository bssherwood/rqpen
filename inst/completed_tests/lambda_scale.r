library(mvtnorm)
library(grpreg)
library(devtools)
install_github("bssherwood/rqpen",force=TRUE)
library(rqPen)
library(splines)


n <- 100
p <- 8
x <- matrix( rnorm(n*p), ncol=p)

y <- 1 + x[,1] + x[,3] - x[,8] + rnorm(n)

#test scale for non-group penalties
lp_scad <- cv.rq.pen(x,y,penalty="SCAD",alg="LP")
qicd_scad <- cv.rq.pen(x,y,penalty="SCAD",alg="QICD")

lp_gscad$cv
qicd_gscad$cv
#lambdas are on the same scale

#group penalty lambdas are on the same scale
g <- rep(seq(1,4),each=2)

lp_gscad <- cv.rq.group.pen(x,y,groups=g,penalty="SCAD",alg="LP")
qicd_gscad <- cv.rq.group.pen(x,y,groups=g,penalty="SCAD",alg="QICD")

lp_gscad$cv
qicd_gscad$cv
