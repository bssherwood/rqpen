library(rqPen)
source("my_rqPen.R")

set.seed(1)
x <- matrix(rnorm(800),nrow=100)
y <- 1 + x[,1] - 3*x[,5] + rnorm(100)

tau <- c(0.2, 0.5, 0.8)

myModel <- myrq.lasso.fit(x, y, lambda = 0.1, tau = tau)

# no discrepancies!
for (i in seq_along(tau)) {
    cat("Tau:", tau[i], fill = T)
    rqfit <- rq.lasso.fit(x, y, lambda = 0.1, tau = tau[i])
    cat("Max coef discrepancy:", max(abs(myModel$coefficients[i, ] - rqfit$coefficients)), fill = T)
    cat("Max residual discrepancy:", max(abs(myModel$residuals[, i] - rqfit$residuals)), fill = T)
    cat("PenRho discrepancy:", max(abs(myModel$PenRho[i] - rqfit$PenRho)), fill = T)
    cat("rho discrepancy:", max(abs(myModel$rho[i] - rqfit$rho)), fill = T)
}

set.seed(4)
n <- 1000; p <- 100
x <- matrix(rnorm(n * p), nrow = n)
y <- rowSums(x[, 1:10]) + rnorm(n)
tau <- seq(from = 0.1, to = 0.9, length.out = 20)
system.time({myModel <- myrq.lasso.fit(x, y, lambda = 0.1, tau = tau)})
system.time({
    sapply(seq_along(tau),
           function(i) rq.lasso.fit(x, y, lambda = 0.1, tau = tau[i]))
})
