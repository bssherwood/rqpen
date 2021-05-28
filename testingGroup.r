library(mvtnorm)
library(grpreg)
library(devtools)
library(rqPen)
library(splines)
library(parallel)

x_cor <- .5

#all functions are setup to have values between [-1,1]
f1 <- function(z){
 2*z^3-1
}
 
f2 <- function(z){
 sin(2*pi*z)
}

f3 <- function(z){
 (z-.5)^2*8-1
}
 
get_splines <- function(z,df){
 b_splines <- NULL
 spline_list <- list()
 for(i in 1:dim(z)[2]){
	spline_list[[i]] <- bs(z[,i],df)
 }
 spline_list
}
 
get_new_splines <- function(z,b_list){
 list_length <- length(b_list)
 new_b <- NULL
 for(i in 1:list_length){
	new_b <- cbind(new_b, predict(b_list[[i]],z[,i]))
 }
 new_b
}

val_n <- 1000

index <- 1
set.seed(index)
print(paste("working on index", index))

generateMvNorm <- function(mean,sd, rho=0, n,obs, structure){
   if(length(mean)==1){
     mean <- rep(mean, obs)
   }
   if(structure=="ar1"){
     Sigma <- diag(obs)
     Sigma <- rho^abs(row(Sigma) - col(Sigma))
     Sigma <- Sigma*sd
   }
   else if(structure=="equal"){
     Sigma <- matrix(rep(rho, obs^2), nrow=obs)
     diag(Sigma) <- 1
     Sigma <- Sigma*sd
   }
   else{
     print("Unrecognized type choose ar1 or equal")
   }
   rmvnorm(n, mean, Sigma)
}

n <- 500 #as.numeric(args[1])
b_df <- 3 #as.numeric(args[2])
p <- 300 #as.numeric(args[3]

x <- generateMvNorm(0,1,rho=x_cor,n,p,"ar1")
z <- pnorm(x)
#using t-distribution
y <- 1 + f1(z[,1]) + f2(z[,2]) + f3(z[,3]) + rt(n,3)#rnorm(n,1)#rt(n,3)
b_list <- get_splines(z,b_df)
b <- do.call(cbind,b_list)

groups <- rep(seq(1,p),each=b_df)
active_vars <- seq(1,3)#c(1,3,8)#seq(1,3)
inactive_vars <- seq(4,p)#c(2,4,5,6,7)#seq(4,p)
 
x_new <- generateMvNorm(0,1,x_cor, val_n, p, "ar1")
z_new <- pnorm(x_new)

y_new <- 1 + f1(z_new[,1]) + f2(z_new[,2]) + f3(z_new[,3]) + rt(n,3)#rnorm(n) #rt(n,3)
b_new <- cbind(1,get_new_splines(z_new, b_list)) #warnings expected here

time_4 <- system.time( rq_group_scad_lp  <- cv.rq.group.pen(b,y, groups, criteria="PBIC", penalty="SCAD", alg="LP",lambda=seq(.01,.03,.001)))

