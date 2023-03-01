library(rqPen)

set.seed(1)

x <- matrix(rnorm(25*9,sd=10),ncol=9)
testx <- matrix(rnorm(81),ncol=9)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(25,3)
g <- rep(seq(1:3),3)
tvals <- c(.25,.75)

# lasso
m1 <- rq.group.pen.cv(x,y,groups=g)

c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots1.pdf")
plot(m1)
dev.off()

q1 <- qic.select(m1)


c1 <- coefficients(q1)

p1 <- predict(q1,testx)

# adaptive lasso 
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO")
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",norm=1)
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",norm=1,alg="br")

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
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD")
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1)
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,alg="br")
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,alg="qicd")

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
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP")
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1)
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,alg="br")
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,alg="qicd")

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

# lasso mult tau 
m1 <- rq.group.pen.cv(x,y,groups=g,tau=tvals)

c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots1.pdf")
plot(m1,tau=.25)
dev.off()

q1 <- qic.select(m1)


c1 <- coefficients(q1)

p1 <- predict(q1,testx)

# adaptive lasso mult tau 
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",tau=tvals)
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",norm=1,tau=tvals)
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",norm=1,alg="br",tau=tvals)

c1 <- coefficients(m1)
c2 <- coefficients(m2,tau=.25)
c3 <- coefficients(m3)

p1 <- predict(m1,testx,tau=.75)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])

pdf("plots1.pdf")
plot(m1,tau=.25)
plot(m2,tau=.25)
plot(m3,tau=.25)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3,tau=.25)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx,tau=.75)
p3 <- predict(q3,testx)

#group SCAD
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",tau=tvals)
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,tau=tvals)
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,alg="br",tau=tvals)
m4 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,alg="qicd",tau=tvals)

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4,septau=FALSE)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx,lambda=m4$lambda[5])

pdf("plots1.pdf")
plot(m1,tau=.25)
plot(m2,tau=.75)
plot(m3,tau=.25)
plot(m4,tau=.25)
dev.off()

q1 <- qic.select(m1,septau=FALSE)
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
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",tau=tvals)
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,tau=tvals)
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,alg="br",tau=tvals)
m4 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,alg="qicd",tau=tvals)

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4,septau=FALSE)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx,lambda=m4$lambda[5])

pdf("plots1.pdf")
plot(m1,tau=.25)
plot(m2,tau=.75)
plot(m3,tau=.25)
plot(m4,tau=.25)
dev.off()

q1 <- qic.select(m1,septau=FALSE)
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


# adaptive lasso mult tau 
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",tau=tvals,a=c(3,4,5))
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",norm=1,tau=tvals,a=c(3,4,5))
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gAdLASSO",norm=1,alg="br",tau=tvals,a=c(3,4,5))

c1 <- coefficients(m1)
c2 <- coefficients(m2,tau=.25)
c3 <- coefficients(m3,septau=FALSE)

p1 <- predict(m1,testx,tau=.75)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])

pdf("plots1.pdf")
plot(m1,tau=.25,a=3)
plot(m2,tau=.25,a=4)
plot(m3,tau=.25,a=5)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3,tau=.25)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx,tau=.75)
p3 <- predict(q3,testx)

#group SCAD
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",tau=tvals,a=c(3,4,5))
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,tau=tvals,a=c(3,4,5))
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,alg="br",tau=tvals,a=c(3,4,5))
m4 <- rq.group.pen.cv(x,y,groups=g,penalty="gSCAD",norm=1,alg="qicd",tau=tvals,a=c(3,4,5))

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4,septau=FALSE)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx,lambda=m4$lambda[5])

pdf("plots1.pdf")
plot(m1,tau=.25,a=3)
plot(m2,tau=.75,a=4)
plot(m3,tau=.25,a=5)
plot(m4,tau=.25,a=3)
dev.off()

q1 <- qic.select(m1,septau=FALSE)
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
m1 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",tau=tvals,a=c(3,4,5))
m2 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,tau=tvals,a=c(3,4,5))
m3 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,alg="br",tau=tvals,a=c(3,4,5))
m4 <- rq.group.pen.cv(x,y,groups=g,penalty="gMCP",norm=1,alg="qicd",tau=tvals,a=c(3,4,5))

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4,septau=FALSE)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx,lambda=m4$lambda[5])

pdf("plots1.pdf")
plot(m1,tau=.25,a=3)
plot(m2,tau=.75,a=4)
plot(m3,tau=.25,a=5)
plot(m4,tau=.25,a=3)
dev.off()

q1 <- qic.select(m1,septau=FALSE)
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


