library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)

x <- matrix(rnorm(800),ncol=8)
y <- 1 + x[,1] + x[,8] + rnorm(100)

m1 <- rq.enet(x,y)
m2 <- rq.enet(x,y,a=.5)
m3 <- rq.enet(x,y,a=1)
m4 <- rq.lasso(x,y)

