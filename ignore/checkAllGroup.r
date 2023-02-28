library(rqPen)

set.seed(1)

x <- matrix(rnorm(25*9,sd=10),ncol=9)
testx <- matrix(rnorm(81),ncol=9)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(25,3)
g <- rep(seq(1:3),3)
tvals <- c(.25,.75)

# lasso
m1 <- rq.group.pen(x,y,groups=g)

c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots1.pdf")
plot(m1)
dev.off()

q1 <- qic.select(m1)


c1 <- coefficients(q1)

p1 <- predict(q1,testx)

# adaptive lasso 
m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO")
m2 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1)
m3 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1,alg="br")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])

pdf("plots1.pdf")
plot(m1)
plot(m2)
plot(m3)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)

#group SCAD
m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD")
m2 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1)
m3 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1,alg="br")
m3 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1,alg="qicd")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx,lambda=m4$lambda[5])

pdf("plots1.pdf")
plot(m1)
plot(m2)
plot(m3)
plot(m4)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4,method="PBIC")

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)

#group MCP
m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP")
m2 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1)
m3 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1,alg="br")
m3 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1,alg="qicd")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx,lambda=m4$lambda[5])

pdf("plots1.pdf")
plot(m1)
plot(m2)
plot(m3)
plot(m4)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4,method="PBIC")

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)

r1 <- rq.group.pen(x,y,groups=g)
r2 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO")
r3 <- rq.group.pen(x,y,groups=g,penalty="gSCAD") #still a problem here with the updates of coefficients, I think. 
r4 <- rq.group.pen(x,y,groups=g,penalty="gMCP")
r5 <- rq.group.pen(x,y,groups=g,tau=tvals)

r1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals)
r2 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals)
r3 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals)

r1 <- rq.group.pen(x,y,groups=g,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",norm=1)

r1 <- rq.group.pen(x,y,groups=g,tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",tau=tvals,norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",tau=tvals,norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",tau=tvals,norm=1)

r1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,norm=1)
r2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",tau=tvals,norm=1)
r3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",tau=tvals,norm=1)

rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas)


library(hqreg)
library(glmnet)

set.seed(1)

x <- matrix(rnorm(25*30,sd=10),ncol=30)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(25,3)
g <- rep(seq(1:5),6)
tvals <- c(.25,.75)

rq.group.pen(x,y,groups=g,a=seq(3,5))

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,norm=1,a=seq(3,5))
#m1$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gAdLASSO",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gSCAD",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

m1 <- rq.group.pen(x,y,groups=g,penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m3 <- rq.group.pen(x,y,groups=g,alg="qicd",penalty="gMCP",tau=tvals,norm=1,a=seq(3,5))
m1$models
m2$models
m3$models

