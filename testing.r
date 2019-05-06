library(devtools)
install_github("bssherwood/rqpen", force=TRUE)

library(rqPen)

x <- matrix(rnorm(800),ncol=8)
x_diff <- x
x_diff[,8] <- x_diff[,8]/100

y <- 1 + x[,1] - x[,8] + rnorm(100)

groups <- rep(c(1,2),each=4)

a <- groupMultLambda(x,y,groups,lambda=c(.001,.01),tau=.5)
b <- groupMultLambda(x_diff,y,groups,lambda=c(.001,.01),tau=.5)
