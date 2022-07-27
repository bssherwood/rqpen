library(rqPen)
set.seed(1)

n <- 25
p <- 20

x <- matrix(rnorm(n*p),ncol=p)
y <- x[,1] + x[,3] - x[,8] + rnorm(n)

x <- cbind(x,x) # x is non-singular now and with more columns than rows



r1 <- cv.rq.pen(x,y)
r2 <- cv.rq.pen(x,y,penalty="SCAD",alg="LP")