library(AmesHousing)
ames <- make_ames()
x_g <- cbind(model.matrix(~ Lot_Shape+Garage_Type+Full_Bath
                          +Fireplaces+Lot_Config - 1,ames))
y_g <- log(ames$Sale_Price)
g <-  c(rep(1,4),rep(2,6),3,4,rep(5,4))
gpf <- sqrt(xtabs(~g))
gpf[5] <- 0
taus <- seq(.1,.9,.1)
qAmes1 <- rq.group.pen(x_g,y_g,groups=g, group.pen.factor = gpf, tau=taus,
                      penalty="gSCAD",norm=1,tau.penalty.factor=taus*(1-taus))
qAmes2 <- rq.group.pen(x_g,y_g,groups=g, group.pen.factor = gpf, tau=taus,
                      penalty="gSCAD",tau.penalty.factor=taus*(1-taus))

# Potential cause could be the small number of FR3