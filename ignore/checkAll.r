library(rqPen)

set.seed(1)

x <- matrix(rnorm(100*8,sd=10),ncol=8)
testx <- matrix(rnorm(800),ncol=8)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(100,3)
g <- c(1,1,1,1,2,2,3,3)
tvals <- c(.25,.75)

# lasso with one quantile
m1 <- rq.pen(x,y,alg="huber")
m2 <- rq.pen(x,y,alg="br")
m3 <- rq.pen(x,y,alg="fn")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m1$lambda[5])

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

# elastic net with one quantile

m2 <- rq.pen(x,y,penalty="ENet",a=.5)

c2 <- coefficients(m2)

p2 <- predict(m2,testx)

pdf("plots2.pdf")
plot(m2)
dev.off()

q2 <- qic.select(m2,method="AIC")
c2 <- coefficients(q2)
p2 <- predict(q2,testx)

#SCAD with one quantile

m1 <- rq.pen(x,y,penalty="SCAD")
m2 <- rq.pen(x,y,penalty="SCAD",alg="br")
m3 <- rq.pen(x,y,penalty="SCAD",alg="QICD")
m4 <- rq.pen(x,y,penalty="SCAD",alg="huber")
m5 <- rq.pen(x,y,penalty="SCAD",alg="fn")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)

# Adaptive lasso

m1 <- rq.pen(x,y,penalty="aLASSO")

c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots2.pdf")
plot(m1)
dev.off()

q2 <- qic.select(m1,method="AIC")
c2 <- coefficients(q2)
p2 <- predict(q2,testx)

# MCP with one quantile 
m1 <- rq.pen(x,y,penalty="MCP")
m2 <- rq.pen(x,y,alg="br", penalty="MCP")
m3 <- rq.pen(x,y,alg="QICD",penalty="MCP")
m4 <- rq.pen(x,y,alg="huber",penalty="MCP")
m5 <- rq.pen(x,y,alg="fn",penalty="MCP")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)

# Lasso with multiple quantiles 

m1 <- rq.pen(x,y,tau=tvals)
m2 <- rq.pen(x,y,alg="br",tau=tvals)
m3 <- rq.pen(x,y,alg="huber",tau=tvals)
m4 <- rq.pen(x,y,alg="fn",tau=tvals)

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3,lambda=m3$lambda[5])
c4 <- coefficients(m4,tau=.25)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)

pdf("plots1.pdf")
plot(m1)
plot(m2)
plot(m3)
plot(m4)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)

# Elastic net with multiple quantiles

m1 <- rq.pen(x,y,tau=tvals, penalty="ENet", a=.25)
m2 <- rq.pen(x,y,a=.5,tau=tvals, penalty="ENet")

c1 <- coefficients(m1)
c2 <- coefficients(m2)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)

pdf("plots1.pdf")
plot(m1)
plot(m2)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")


c1 <- coefficients(q1)
c2 <- coefficients(q2)


p1 <- predict(q1,testx)
p2 <- predict(q2,testx)


#SCAD with multiple quantiles
m1 <- rq.pen(x,y,tau=tvals, penalty="SCAD")
m2 <- rq.pen(x,y,alg="br",tau=tvals, penalty="SCAD")
m3 <- rq.pen(x,y,alg="QICD",tau=tvals, penalty="SCAD")
m4 <- rq.pen(x,y,alg="huber",tau=tvals, penalty="SCAD")
m5 <- rq.pen(x,y,alg="fn",tau=tvals, penalty="SCAD")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1, tau=.25)
plot(m2,tau=.25)
plot(m3,tau=.75)
plot(m4,tau=.25)
plot(m5,tau=.75)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)

#Adaptive Lasso
m1 <- rq.pen(x,y,tau=tvals,penalty="aLASSO")

c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots2.pdf")
plot(m1)
dev.off()

q2 <- qic.select(m1,method="AIC")
c2 <- coefficients(q2)
p2 <- predict(q2,testx)

#MCP with quantiles
m1 <- rq.pen(x,y,tau=tvals, penalty="MCP")
m2 <- rq.pen(x,y,alg="br",tau=tvals, penalty="MCP")
m3 <- rq.pen(x,y,alg="QICD",tau=tvals, penalty="MCP")
m4 <- rq.pen(x,y,alg="huber",tau=tvals, penalty="MCP")
m5 <- rq.pen(x,y,alg="fn",tau=tvals, penalty="MCP")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1, tau=.25)
plot(m2,tau=.25)
plot(m3,tau=.75)
plot(m4,tau=.25)
plot(m5,tau=.75)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)

# SCAD with multiple a
m1 <- rq.pen(x,y,a=c(3,4,5), penalty="SCAD")
m2 <- rq.pen(x,y,alg="huber", a = c(3,4,5), penalty="SCAD")
m3 <- rq.pen(x,y,alg="QICD", a = c(3,4,5), penalty="SCAD")
m4 <- rq.pen(x,y,alg="fn", a = c(3,4,5), penalty="SCAD")
m5 <- rq.pen(x,y,alg="br", a = c(3,4,5), penalty="SCAD")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1)
plot(m2)
plot(m3)
plot(m4)
plot(m5)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)

#Adaptive lasso with multiple a
m1 <- rq.pen(x,y,penalty="aLASSO", a=c(1,2,3))

c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots2.pdf")
plot(m1,a=2)
dev.off()

q2 <- qic.select(m1,method="AIC")
c2 <- coefficients(q2)
p2 <- predict(q2,testx)

#MCP with multiple a
m1 <- rq.pen(x,y,penalty="MCP", a=c(3,4,5))
m2 <- rq.pen(x,y,alg="huber", penalty="MCP",a=c(3,4,5))
m3 <- rq.pen(x,y,alg="QICD",penalty="MCP",a=c(3,4,5))
m4 <- rq.pen(x,y,alg="br", penalty="MCP",a=c(3,4,5))
m5 <- rq.pen(x,y,alg="fn", penalty="MCP",a=c(3,4,5))

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1,a=3)
plot(m2,a=4)
plot(m3,a=5)
plot(m4,a=3)
plot(m5,a=3)
dev.off()

q1 <- qic.select(m1)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)

# Enet with multiple a and tau 
m1 <- rq.pen(x,y,a=c(.2,.5,.7),tau=tvals, penalty="ENet")
c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots2.pdf")
plot(m1,a=.2,tau=.25)
dev.off()

q2 <- qic.select(m1,method="AIC")
c2 <- coefficients(q2)
p2 <- predict(q2,testx)

#Multiple a and tau with SCAD
m1 <- rq.pen(x,y,tau=tvals,a=c(3,4,5), penalty="SCAD")
m2 <- rq.pen(x,y,alg="br",tau=tvals,a=c(3,4,5), penalty="SCAD")
m3 <- rq.pen(x,y,alg="QICD",tau=tvals,a=c(3,4,5), penalty="SCAD")
m4 <- rq.pen(x,y,alg="huber",tau=tvals,a=c(3,4,5), penalty="SCAD")
m5 <- rq.pen(x,y,alg="fn",tau=tvals,a=c(3,4,5), penalty="SCAD")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1,a=3, tau=.25)
plot(m2,a=4,tau=.75)
plot(m3,a=5,tau=.25)
plot(m4,a=3,tau=.25)
plot(m5,a=3,tau=.25)
dev.off()

q1 <- qic.select(m1,septau=FALSE)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)

#Multiple quantiles and tau values for adaptive lasso
m1 <- rq.pen(x,y,tau=tvals,penalty="aLASSO",a=c(1,2,3))
c1 <- coefficients(m1)

p1 <- predict(m1,testx)

pdf("plots2.pdf")
plot(m1,a=2,tau=.25)
dev.off()

q2 <- qic.select(m1,method="AIC")
c2 <- coefficients(q2)
p2 <- predict(q2,testx)

#Multiple quantiles and tau values for MCP

m1 <- rq.pen(x,y,tau=tvals,penalty="MCP",a=c(3,4,5))
m2 <- rq.pen(x,y,alg="br", penalty="MCP",tau=tvals,a=c(3,4,5))
m3 <- rq.pen(x,y,alg="QICD",penalty="MCP",tau=tvals,a=c(3,4,5))
m4 <- rq.pen(x,y,alg="huber",penalty="MCP",tau=tvals,a=c(3,4,5))
m5 <- rq.pen(x,y,alg="fn",tau=tvals,a=c(3,4,5), penalty="MCP")

c1 <- coefficients(m1)
c2 <- coefficients(m2)
c3 <- coefficients(m3)
c4 <- coefficients(m4)
c5 <- coefficients(m5)

p1 <- predict(m1,testx)
p2 <- predict(m2,testx)
p3 <- predict(m3,testx,lambda=m3$lambda[5])
p4 <- predict(m4,testx)
p5 <- predict(m5,testx)

pdf("plots1.pdf")
plot(m1,a=3, tau=.25)
plot(m2,a=4,tau=.75)
plot(m3,a=5,tau=.25)
plot(m4,a=3,tau=.25)
plot(m5,a=3,tau=.25)
dev.off()

q1 <- qic.select(m1,septau=FALSE)
q2 <- qic.select(m2,method="AIC")
q3 <- qic.select(m3,method="PBIC")
q4 <- qic.select(m4)
q5 <- qic.select(m5)

c1 <- coefficients(q1)
c2 <- coefficients(q2)
c3 <- coefficients(q3)
c4 <- coefficients(q4)
c5 <- coefficients(q5)

p1 <- predict(q1,testx)
p2 <- predict(q2,testx)
p3 <- predict(q3,testx)
p4 <- predict(q4,testx)
p5 <- predict(q5,testx)