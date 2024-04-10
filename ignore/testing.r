library(devtools)
install_github("bssherwood/rqpen", force=TRUE)
3
library(rqPen)
n<- 200
p<- 20
X<- matrix(rnorm(n*p),n,p)
y<- -2+X[,1]+0.5*X[,2]-X[,3]-0.5*X[,7]+X[,8]-0.2*X[,9]+rt(n,2)
taus <- seq(0.1, 0.9, 0.2)
folds <- rqPen:::randomly_assign(n,10)

cvfit <- rq.gq.pen.cv(x=X, y=y, tau=taus, foldid=folds)
cvfit2 <- cv.hrq_tau_glasso(x=X, y=y, tau=taus, folds=folds)


library(rqPen)

x <- matrix(rnorm(800),ncol=80)
x_diff <- x
x_diff[,8] <- x_diff[,8]/100

y <- 1 + x[,1] - x[,8] + rnorm(100)

groups <- rep(c(1,2),each=4)

a <- groupMultLambda(x,y,groups,lambda=c(.001,.01),tau=.5)
coefficients(a)
b <- groupMultLambda(x_diff,y,groups,lambda=c(.001,.01),tau=.5)
coefficients(b)

a <- rq.group.fit(x,y,groups,lambda=.001,tau=.25,penalty="SCAD")
coefficients(a)
b <- rq.group.fit(x_diff,y,groups,lambda=.001,tau=.25,penalty="SCAD")
coefficients(b)

a <- rq.nc.fit(x,y,lambda=.001,tau=.25,penalty="SCAD")
coefficients(a)
b <- rq.nc.fit(x_diff,y,lambda=.001,tau=.25,penalty="SCAD")
coefficients(b)

folds <- randomly_assign(100,10)

a <- cv.rq.pen(x,y,foldid=folds)
coefficients(a)
b <- cv.rq.pen(x_diff,y,foldid=folds)
coefficients(b)

a <- cv.rq.pen(x,y,foldid=folds, penalty="SCAD")
coefficients(a)
b <- cv.rq.pen(x_diff,y,foldid=folds, penalty="SCAD")
coefficients(b)

a <- cv.rq.pen(x,y,foldid=folds, penalty="SCAD", alg="LP")
coefficients(a)
b <- cv.rq.pen(x_diff,y,foldid=folds, penalty="SCAD", alg="LP")
coefficients(b)

a <- cv.rq.group.pen(x,y,groups,foldid=folds)
coefficients(a)
b <- cv.rq.group.pen(x_diff,y,groups,foldid=folds)
coefficients(b)

a <- cv.rq.group.pen(x,y,groups,alg="LP",foldid=folds)
coefficients(a)
b <- cv.rq.group.pen(x_diff,y,groups,alg="LP",foldid=folds)
coefficients(b)

x <- matrix(rnorm(800),ncol=80)

y <- 1 + x[,1] - x[,8] + rnorm(10)

groups <- rep(seq(1,8),each=10)

test <- cv.rq.group.pen(x,y,groups,alg="QICD",criteria="BIC",penalty="SCAD")


