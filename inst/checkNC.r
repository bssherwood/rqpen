library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)

x <- matrix(rnorm(800),ncol=8)
y <- 1 + x[,1] + x[,8] + rnorm(100)

m1 <- rq.nc(x,y,a=c(3,4))

