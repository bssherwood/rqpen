library(devtools)
install_github("bssherwood/rqpen", force=TRUE)

x <- matrix(rnorm(800),ncol=8)
x_diff <- x_
x_diff[,8] <- x_diff[,8]/100

a <- groupMultLambda(x,y,groups,lambda=c(.001,.01),tau=.5)
b <- groupMultLambda(x_diff,y,groups,lambda=c(.001,.01),tau=.5)
