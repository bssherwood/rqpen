library(devtools)
install_github("bssherwood/rqpen")#, force=TRUE)
3
library(rqPen)
library(quantreg)
data(barro)
set.seed(1)
y <- barro$y.net
x <- as.matrix(barro[,-1])

r1 <- rq.gq.pen.cv(x,y, tau=c(.3,.5,.7))
coefficients(r1)
#Best model is totally sparse. Kind of surprising given that this is a motivating data set in the quantreg package. 

summary(y)
summary(predict(r1,newx=x))
summary(predict(r1,newx=x,lambdaIndex=2, useDefaults = FALSE))
summary(predict(r1,newx=x,lambdaIndex=50, useDefaults = FALSE))

x <- matrix(rnorm(8000),ncol=8)
y <- 1 + x[,1] + x[,3] + x[,8] + rnorm(1000,sd=.1)
y <- y/10

g1 <- rq.gq.pen.cv(x,y,tau=c(.1,.5,.9))