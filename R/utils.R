######################################################
###  Functions to be called in the main algorithm  ###
######################################################

#' The quantile check function
#'
#' @param r A vector or scaler
#' @param tau Percentile
#'
#' @return Quantile loss
#' @noRd
rq.loss<- function(r, tau){
  (abs(r)+(2*tau-1)*r)/2
}

# Quantile check function for augmented data (length(tau)>=2)
rq.loss.aug<- function(r, tau, n){
  tau.vec <- rep(tau, each=n)
  (abs(r)+(2*tau.vec-1)*r)/2
}



##############################


## Huberized quantile loss (not used in this implementation)
rq.huber<- function(r, tau, gmma){
  r<- as.vector(r)
  (huber.loss(r,gmma)+(2*tau-1)*r)/2
}
rq.huber.aug<- function(r, tau, gmma, n){
  tau.vec <- rep(tau, each=n)
  r<- as.vector(r)
  (huber.loss(r,gmma)+(2*tau.vec-1)*r)/2
}
##############################


#' Scale back coefficients
#'
#' @param coefs Coefficient vector
#' @param mu_x 
#' @param sigma_x 
#' @param intercept 
#' @noRd
transform_coefs <- function(coefs,mu_x,sigma_x,intercept=TRUE){
  if(intercept){
    new_coefs<- coefs[-1]/sigma_x
    intercept<- coefs[1] - sum(coefs[-1]*mu_x/sigma_x)
    new_coefs<- c(intercept, new_coefs)
  } else{
    new_coefs<- coefs/sigma_x
  }
  new_coefs
}

transform_coefs.aug <- function(coefs,mu_x,sigma_x, ntau,intercept=TRUE){
  if(intercept){
    new_coefs<- coefs[-(1:ntau)]/sigma_x
    temp <- matrix(coefs[-(1:ntau)]*rep(mu_x/sigma_x, each=ntau), length(mu_x), ntau, byrow=TRUE)
    intercept<- coefs[1:ntau] - apply(temp, 2, sum)
    new_coefs<- c(intercept, new_coefs)
  } else{
    new_coefs<- coefs/rep(sigma_x, each=ntau)
  }
  new_coefs
}
############################


#  
#' First order derivative w.r.t. residual
#'
#' @param r residual
#' @param tau Percentile
#' @param gmma Huber parameter
#' @param n sample size
#'
#' @return First order derivative, a vector of n
#' @noRd
#'
rq.huber.deriv<- function(r, tau, gmma){
  r<- as.vector(r)
  le.ind<- which(abs(r) <= gmma)
  if(length(le.ind)!= 0){
    l.vec<- r
    l.vec[le.ind]<- (r[le.ind]/gmma+(2*tau-1))/2
    l.vec[-le.ind]<- (sign(r[-le.ind])+(2*tau-1))/2
  } else{
    l.vec <- (sign(r)+(2*tau-1))/2
  }
  return(l.vec)
} # end of function

rq.huber.deriv.aug<- function(r, tau, gmma, n){
  tau.vec <- rep(tau, each=n)
  r<- as.vector(r)
  le.ind<- which(abs(r) <= gmma)
  if(length(le.ind)!= 0){
    l.vec<- r
    l.vec[le.ind]<- (r[le.ind]/gmma+(2*tau.vec[le.ind]-1))/2
    l.vec[-le.ind]<- (sign(r[-le.ind])+(2*tau.vec[-le.ind]-1))/2
  } else{
    l.vec <- (sign(r)+(2*tau.vec-1))/2
  }
  return(l.vec)
} # end of function


neg.gradient.aug <- function(r,weights,tau,gmma,x,n,apprx){  # input r, weights, x needs to be the augmented version
  if(apprx=="huber"){
    wt_deriv <- as.vector(weights*rq.huber.deriv.aug(r, tau, gmma, n))
  }else{
    stop("huber approximation is the only allowed approach")
  }
  
  if(is.null(dim(x))){
    mean(x*wt_deriv)
  } else{
    apply(x*wt_deriv,2,sum)/n
  }
} # end of function

## l2 norm
l2norm<- function(x){
  sqrt(sum(x^2))
} # end of function


#' A function used for kkt condition check in order to verify strong rule
#'
#' @param r Residual
#' @param weights Observation weights
#' @param w Weight of shrinkage parameter. Default is square root of each group size.
#' @param gmma Huber parameter
#' @param tau Percentile
#' @param group.index A vector of group index, e.g., (1,1,1,2,2,2,3,3)
#' @param inactive.ind A vector of index for group with zero coefficients according to strong rule, e.g., (1,2,4) 
#' @param lambda Shrinkage parameter
#' @param x Reduced design matrix
#' @param n Sample size
#' @param apprx Approximation method
#'
#' @return True and False to indicate if KKT condition is met.
#' @noRd
kkt_check<- function(r, weights, w, gmma, tau, group.index, inactive.ind, lambda, x, n, apprx){
  grad<- -neg.gradient(r, weights, tau, gmma, x, apprx)
  grad.norm<- tapply(grad, group.index, l2norm)/w
  bad_spots <- grad.norm > lambda
  if(sum(bad_spots[inactive.ind])==0){
    list(kkt=TRUE,new_groups=NULL)
  } else{
    list(kkt=FALSE,new_groups=inactive.ind[which(bad_spots[inactive.ind])])
  }
} # end of function

kkt_check.aug<- function(r, weights, w, gmma, tau, group.index, inactive.ind, lambda, x, n, apprx){
  grad<- -neg.gradient.aug(r, weights, tau, gmma, x, n, apprx)
  grad.norm<- tapply(grad, group.index, l2norm)/w
  bad_spots <- grad.norm > lambda
  if(sum(bad_spots[inactive.ind])==0){
    list(kkt=TRUE,new_groups=NULL)
  } else{
    list(kkt=FALSE,new_groups=inactive.ind[which(bad_spots[inactive.ind])])
  }
} # end of function

## use Rcpp implemented neg_gradient_aug()
kkt_check_aug<- function(r, weights, w, gmma, tau, group.index, inactive.ind, lambda, x, ntau){
  grad<- -neg_gradient_aug(r, weights, tau, gmma, x, ntau) 
  grad.norm<- tapply(grad, group.index, l2norm)/w
  bad_spots <- grad.norm > lambda
  if(sum(bad_spots[inactive.ind])==0){
    list(kkt=TRUE,new_groups=NULL)
  } else{
    list(kkt=FALSE,new_groups=inactive.ind[which(bad_spots[inactive.ind])])
  }
} # end of function

getGQCoefs <- function(beta,taupos,p,k){
  beta[seq(taupos,k*(p-1)+taupos,by=k),]
}
