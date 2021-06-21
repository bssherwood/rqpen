kernel_weights <- function(obs_data,obs_ind,...){
   if(is.null(dim(obs_data))){
      d <- 1
      n <- length(obs_data)
   } else{
      d <- dim(obs_data)[2]
      n <- dim(obs_data)[1]
   }
   tune_h <- sd(obs_data)*n^{-1/(d+2)}
   kernel_est <- kernel_estimates(obs_data,obs_ind,tune_h,...)
   1/kernel_est
}

check <- function(x,tau=.5){
   x*(tau - (x<0))
}

pos_part <- function(x){
  ifelse( x < 0, x, 0 ) # min(x,0) # 
}

lasso <- function(x,lambda=1){
   lambda*abs(x)
}

scad <- function(x, lambda=1, a=3.7){
  absx <- abs(x)
  ifelse( absx < lambda, # Evaluate this
          lambda*absx, # If true (absx < lambda)
          ifelse( absx < a*lambda, # If false, evaluate this
                  ( (a^2-1)*lambda^2 - (absx-a*lambda)^2 ) / ( 2*(a-1) ), # If true (absx < a*lambda)
                  (a+1)*lambda^2 / 2 # If false (absx > a*lambda)
                ) 
        )
}


scad_deriv <- function(x, lambda=1,a=3.7){
	ll <- length(lambda)
	if(ll !=1 && ll != length(x)){
		stop("lambda must be of length 1 or same length as x")
	}
	absx <- u <- abs(x)
	u[] <- 0
	index <- absx < a*lambda & absx > 0
	if(ll==1){
		u[ index ] <- ifelse( absx[ index ] <= lambda, lambda, ( a*lambda - absx[ index ] )/( a-1 ) )
		u[index] <- u[index]*sign( x[index] )
		u[ x == 0 ] <- lambda # because we take derivative as x approaces 0 from above
	} else{
		u[index] <- ifelse( absx[ index ] <= lambda[index], lambda[index], ( a*lambda[index] - absx[ index ] )/( a-1 ) )
		u[index] <- u[index]*sign(x[index])
		u[ x==0] <- lambda[x==0]
	}
	u
}

#scad_1_deriv <- function(x,lambda=1,a=3.7){
#  sign(x)*lambda
#}
#
#scad_2_deriv <- function(x,lambda=1,a=3.7){
#  sign(x)*lambda*(1-( pos_part(a*lambda-abs(x)) / ( lambda*(a-1))))*(abs(x) > lambda)
#}

mcp <- function(x, lambda=1, a=3){
  absx <- abs(x)
  ifelse( absx < a*lambda, # Evaluate this
          lambda*(absx - absx^2/(2*a*lambda)), # If true (absx < a*lambda)
          a*lambda^2/2 # If false
        )
}


mcp_deriv <- function(x, lambda=1, a=3){
	ll <- length(lambda)
	if(ll !=1 && ll != length(x)){
		stop("lambda must be of length 1 or same length as x")
	}
	u <- x
	u[] <- 0
	index <- abs(x) < a*lambda
	if(ll == 1){
		u[ index ] <- ifelse( x[index] == 0,
							lambda,
							lambda*sign(x[index]) - x[index]/a
						  )
	} else{
		u[ index ] <- ifelse( x[index] == 0,
							lambda[index],
							lambda[index]*sign(x[index]) - x[index]/a
						  )
	}
	u
}

alasso_wt <- function(x,lambda=1,a=1){
	lambda*(1/abs(x))^a
}



square <- function(x){
  x^2
}

randomly_assign <- function(n,k){
#randomly assign n samples into k groups
   small_set <- floor(n/k)
   group_assign <- NULL
   if(n %% k == 0){
     group_assign <-  rep(seq(1,k),n/k)
   } else{
     remainder <- n %% k
     for(i in 1:remainder){
        group_assign <- c(group_assign, rep(i,small_set+1))
     }
     group_assign <- c(group_assign, rep(seq((i+1),k),small_set))
   }
   sample(group_assign)
}

rq.lasso.fit.mult <- function(x,y,tau_seq=c(.1,.3,.5,.7,.9),lambda=NULL,weights=NULL,intercept=TRUE,coef.cutoff=.00000001,...){
   warning("rq.lasso() fits for multiple tau and provides a sequence of lambda and faster algorithms")
   model_list <- list()
   iter <- 1
   for(tau in tau_seq){
      model_list[[iter]] <- rq.lasso.fit(x,y,tau,lambda,weights,intercept,coef.cutoff,...)
      iter <- iter+1
   }
   model_list
}

transform_coefs <- function(coefs,mu_x,sigma_x,intercept=TRUE){
  new_coefs <- coefs
  if(intercept){
	  intercept <- coefs[1]
	  for(j in 2:length(coefs)){
		new_coefs[j] <- coefs[j]/sigma_x[j-1]
		intercept <- intercept - coefs[j]*mu_x[j-1]/sigma_x[j-1]
	  }
	  new_coefs[1] <- intercept
  } else{
  	for(j in 1:length(coefs)){
		new_coefs[j] <- coefs[j]/sigma_x[j]
	}
  }
  new_coefs
}

get_coef_pen <- function(coefs,lambda,intercept,penVars,penalty="LASSO",a=NULL){
	if(intercept){
		coefs <- coefs[-1]
	}
	if(is.null(penVars)==FALSE){
		coefs <- coefs[penVars]
	}
	if( penalty == "LASSO" ){ ### PenRho for LASSO
      sum( abs( coefs )*lambda )
    } else if( penalty == "SCAD" ){ ### PenRho for SCAD
      sum( scad( coefs, lambda, a ))
    } else { ### PenRho for MCP
      sum( mcp( coefs, lambda, a ))
    }
}

# modified version (rescaled) in order to approximate LAD
huber.loss <- function(r,gamma){
  le.ind<- which(abs(r) <= gamma)
  ##If statement is a Ben addition. There was a problem if all the values of r had a magnitude larger than gamma. 
  if(length(le.ind)==0){
    return( abs(r)-gamma/2)
  } else{
    #Ben: I put in the below comment because huber.loss was not working well with cv function
    #l.vec<- as.vector(r)
    #l.vec[le.ind]<- (r[le.ind])^2/2/gamma
    #l.vec[-le.ind]<- abs(r[-le.ind])-gamma/2
    #return(l.vec)
    r[le.ind]<- (r[le.ind])^2/2/gamma
    r[-le.ind]<- abs(r[-le.ind])-gamma/2
    return(r)
  }
}
##############################
## compare loss functions
# r<- seq(-10,10,0.01)
# gamma<- 1
# plot(r, tanh.loss(r, gamma), type = "l", ylab = "loss")
# lines(r, huber.loss(r, gamma), col=2)
# lines(r, r^2, col=3)
# legend("bottomright", legend = c("tanh", "huber", "squared"), lty = 1, col=1:3)
# dev.copy2pdf(file="loss_compare.pdf")

rq.huber<- function(r, tau, gamma){
  r<- as.vector(r)
  (huber.loss(r,gamma)+(2*tau-1)*r)/2
}

#  returns a vector of n
rq.huber.deriv<- function(r, tau, gamma){
  r<- as.vector(r)
  le.ind<- which(abs(r) <= gamma)
  if(length(le.ind)!= 0){
    l.vec<- r
    l.vec[le.ind]<- (r[le.ind]/gamma+(2*tau-1))/2
    l.vec[-le.ind]<- (sign(r[-le.ind])+(2*tau-1))/2
  } else{
    l.vec <- (sign(r)+(2*tau-1))/2
  }
  return(l.vec)
} 

neg.gradient <- function(r,weights,tau,gamma,x,apprx){
  if(apprx=="huber"){
    wt_deriv <- as.vector(weights*rq.huber.deriv(r, tau, gamma))
  }else{
    wt_deriv <- as.vector(weights*rq.tanh.deriv(r, tau, gamma))
  }
  
  if(is.null(dim(x))){
    mean(x*wt_deriv)
  } else{
    apply(x*wt_deriv,2,mean)
  }
}

getLamMax <- function(x,y,tau=.5,gamma=.2,gamma.max=4,gamma.q=.1,penalty="lasso"){
	n <- length(y)
	returnVal <- 0
	
	#to-do: remove for loop
	for(tau_val in tau){
		inter <- quantile(y,tau_val)
		r <- y - inter
	
		gamma0<- min(gamma.max, max(gamma, quantile(abs(r), probs = gamma.q)))
		grad_k<- -neg.gradient(r, rep(1,n), tau_val, gamma=gamma0, x, apprx="huber")
		returnVal <- max(c(returnVal,abs(grad_k)))
	}
	returnVal
}

# If tau.pen is set to true then the reported lambdas are actually lambda*sqrt(tau*(1-tau))
rq.lasso <- function(x,y,tau=.5,lambda=NULL,nlambda=100,eps=.0001, penalty.factor = rep(1, ncol(x)),
						alg=ifelse(sum(dim(x))<200,"huber","br"),scalex=TRUE,tau.pen=FALSE,...){
	if(alg == "lp"){
	#use br as the default for linear programming 
		alg <- "br"
	}
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	lpf <- length(penalty.factor)
	pfmat <- FALSE
	
	if(sum(tau <= 0 | tau >=1)>0){
		stop("tau needs to be between 0 and 1")
	}
	if(sum(penalty.factor<0)>0){
		stop("penalty factors must be positive")
	}
	if(nt==1){
		if(lpf!=p){
			stop("penalty factor must be of length p")
		}
		if(tau.pen){
			penalty.factor <- penalty.factor*sqrt(tau*(1-tau))
			warning("tau.pen set to true for a single value of tau should return simliar answers as tau.pen set to false. Makes more sense to use it when fitting multiple taus at the same time")
		}
	} else{
		if(length(penalty.factor) == p*nt){
			pfmat <- TRUE
		}
		else if(length(penalty.factor) !=p){
			stop("penalty factor must be a length p vector or a nt by p matrix, where nt is the number of taus")
		}
		if(tau.pen){
			pfmat <- TRUE
			if(length(penalty.factor)==p){
				penalty.factor <- matrix(rep(penalty.factor,nt),nrow=nt,byrow=TRUE)
			} else{
				penalty.factor <- penalty.factor*sqrt(tau*(1-tau))
			}
		}
		
	}
	if(is.null(lambda)){
		lamMax <- getLamMax(x,y,tau)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	if(alg=="huber"){
		if(length(lambda)==1){
			stop("The Huber algorithm requires at least 2 values of lambda")
		}
		returnVal <- rq.lasso.huber(x,y,tau,lambda,penalty.factor,scalex,pfmat,...)
	} else{
		if(nt > 1){
			models <- list()
			for(i in 1:nt){
				coefs <- NULL
				for(lam in lambda){
					sublam <- lam*penalty.factor
					subm <- rq.lasso.fit(x,y,tau[i],lambda=sublam, method=alg,scalex=scalex, ...)
					coefs <- cbind(coefs,coefficients(subm))
				}
				models[[i]] <- rq.lasso.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor)
			}
		} else{
			coefs <- NULL
			for(lam in lambda){
				sublam <- lam*penalty.factor
				subm <- rq.lasso.fit(x,y,tau[i],lambda=sublam, method=alg,scalex=scalex, ...)
				coefs <- cbind(coefs,coefficients(subm))
			}
			models <- rq.lasso.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor)
		}
		returnVal <- list(models=models, n=n, p=p,alg=alg,tau=tau,lambda=lambda,penalty.factor=penalty.factor)
	}
	returnVal$penalty <- "lasso"
	class(returnVal) <- "rq.pen.seq"
	returnVal
}

rq.enet <- function(x,y,tau=.5,lambda=NULL,nlambda=100,eps=.0001, penalty.factor = rep(1, ncol(x)),scalex=TRUE,tau.pen=FALSE,alpha=0,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	lpf <- length(penalty.factor)
	pfmat <- FALSE
	
	if(sum(tau <= 0 | tau >=1)>0){
		stop("tau needs to be between 0 and 1")
	}
	if(sum(penalty.factor<0)>0){
		stop("penalty factors must be positive")
	}
	if(nt==1){
		if(lpf!=p){
			stop("penalty factor must be of length p")
		}
		if(tau.pen){
			penalty.factor <- penalty.factor*sqrt(tau*(1-tau))
			warning("tau.pen set to true for a single value of tau should return simliar answers as tau.pen set to false. Makes more sense to use it when fitting multiple taus at the same time")
		}
	} else{
		if(length(penalty.factor) == p*nt){
			pfmat <- TRUE
		}
		else if(length(penalty.factor) !=p){
			stop("penalty factor must be a length p vector or a nt by p matrix, where nt is the number of taus")
		}
		if(tau.pen){
			pfmat <- TRUE
			if(length(penalty.factor)==p){
				penalty.factor <- matrix(rep(penalty.factor,nt),nrow=nt,byrow=TRUE)
			} else{
				penalty.factor <- penalty.factor*sqrt(tau*(1-tau))
			}
		}
		
	}
	if(is.null(lambda)){
		lamMax <- getLamMax(x,y,tau)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	if(length(lambda)==1){
			stop("The Huber algorithm requires at least 2 values of lambda and elastic net only uses the Huber algorithm")
	}
	returnVal <- rq.lasso.huber(x,y,tau,lambda,penalty.factor,scalex,pfmat,alpha=alpha,...)
	if(alpha==0){
		returnVal$penalty <- "ridge"
	} else{
		returnVal$penalty <- "enet"
		returnVal$alpha <- alpha
	}
	class(returnVal) <- "rq.pen.seq"
	returnVal	
}


rq.lla <- function(obj,x,y,penalty="SCAD",a=ifelse(penalty=="SCAD",3.7,3),...){
	nt <- length(obj$tau)
	if(penalty=="SCAD"){
		derivf <- scad_deriv
	} else if(penalty=="MCP"){
		derivf <- mcp_deriv
	} else if(penalty=="aLasso"){
		derivf <- alasso_wt
	} else{
		stop("Penalty must be SCAD, MCP or aLasso")
	}
	lampen <- as.numeric(obj$penalty.factor %*% t(obj$lambda))
	ll <- length(obj$lambda)
	if(nt == 1){
		pfs <- matrix(derivf(as.numeric(abs(coefficients(obj$models)[-1,])),lampen,a=a),ncol=ll)
		for(i in 1:ll){
			if(obj$alg=="huber"){
				update_est <- coefficients(rq.lasso(x,y,obj$tau,lambda=c(2,1),penalty.factor=pfs[,i],alg=obj$alg,...)$models)[,2]

			} else{
				update_est <- coefficients(rq.lasso(x,y,obj$tau,lambda=1,penalty.factor=pfs[,i],alg=obj$alg,...)$models)
			}
			obj$models$coefficients[,i] <- update_est
		}
	} else{
		for(j in 1:nt){
			pfs <- matrix(derivf(as.numeric(abs(coefficients(obj$models[[j]])[-1,])),lampen,a=a),ncol=ll)
			for(i in 1:ll){
				if(obj$alg=="huber"){
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=c(2,1),penalty.factor=pfs[,i],alg=obj$alg,...)$models)[,2]

				} else{
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=1,penalty.factor=pfs[,i],alg=obj$alg,...)$models)
				}
				obj$models[[j]]$coefficients[,i] <- update_est
			}
		}
	}
	obj$penalty <- penalty
	obj
}

rq.group.lla <- function(obj,x,y,penalty="SCAD",a=ifelse(penalty=="SCAD",3.7,3),norm=2,...){
	nt <- length(obj$tau)
	if(penalty=="SCAD"){
		derivf <- scad_deriv
	} else if(penalty=="MCP"){
		derivf <- mcp_deriv
	} else if(penalty=="aLasso"){
		derivf <- alasso_wt
	} else{
		stop("Penalty must be SCAD, MCP or aLasso")
	}
	lampen <- as.numeric(obj$penalty.factor %*% t(obj$lambda))
	ll <- length(obj$lambda)
	if(nt == 1){
		pfs <- matrix(derivf(as.numeric(abs(coefficients(obj$models)[-1,])),lampen,a=a),ncol=ll)
		for(i in 1:ll){
			if(obj$alg=="huber"){
				update_est <- coefficients(rq.lasso(x,y,obj$tau,lambda=c(2,1),penalty.factor=pfs[,i],alg=obj$alg,...)$models)[,2]

			} else{
				update_est <- coefficients(rq.lasso(x,y,obj$tau,lambda=1,penalty.factor=pfs[,i],alg=obj$alg,...)$models)
			}
			obj$models$coefficients[,i] <- update_est
		}
	} else{
		for(j in 1:nt){
			pfs <- matrix(derivf(as.numeric(abs(coefficients(obj$models[[j]])[-1,])),lampen,a=a),ncol=ll)
			for(i in 1:ll){
				if(obj$alg=="huber"){
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=c(2,1),penalty.factor=pfs[,i],alg=obj$alg,...)$models)[,2]

				} else{
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=1,penalty.factor=pfs[,i],alg=obj$alg,...)$models)
				}
				obj$models[[j]]$coefficients[,i] <- update_est
			}
		}
	}
	obj$penalty <- penalty
	obj
}

rq.nc <- function(x, y, tau=.5,group=1:ncol(X),  penalty=c("gLasso","gAdLasso","gSCAD","gMCP"),a=NULL,lambda=NULL,nlambda=100,eps=.0001,alg="huber", ...) {
	#should look at how ncvreg generates the lambda sequence and combine that with the Huber based approach
	penalty <- match.arg(penalty)
	nt <- length(tau)
	dims <- dim(x)
	p <- dims[2]
	n <- dims[1] 
	if( max(tau) >= 1 | min(tau) <= 0){
		stop("Values for tau must be between 0 and 1") 
	}
	if(penalty == "aLasso" & alg != "huber"){
		alg <- "huber"
		warning("Algorithm switched to huber becaused that is the only one available for adaptive lasso.")
	}
	if(is.null(lambda)){
		lamMax <- getLamMax(x,y,tau,penalty)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	
	if(alg != "qicd" & alg != "QICD"){
		if(penalty=="aLasso"){
			if(is.null(a)){
				a <- 1
			}
			init.model <- rq.enet(x,y,tau,...)
		} else{
			if(is.null(a)){
				if(penalty == "SCAD"){
					a <- 3.7
				} else{
					a <- 3
				}
			}
			init.model <- rq.lasso(x,y,tau,alg=alg,lambda=lambda,...)
		}
		rq.lla(init.model,x,y,penalty,a,...)
	} else{
		if(nt > 1){
			models <- list()
			for(i in 1:nt){
				coefs <- NULL
				for(lam in lambda){
					sublam <- lam
					subm <- rq.nc.fit(x,y,tau[i],lambda=sublam, alg="QICD", ...)
					coefs <- cbind(coefs,coefficients(subm))
				}
				models[[i]] <- rq.lasso.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor=rep(1,p))
			}
		} else{
			coefs <- NULL
			for(lam in lambda){
				sublam <- lam
				subm <- rq.nc.fit(x,y,tau[i],lambda=sublam, alg="QICD",...)
				coefs <- cbind(coefs,coefficients(subm))
			}
			models <- rq.lasso.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor=rep(1,p))
		}
		returnVal <- list(models=models, n=n, p=p,alg=alg,tau=tau,lambda=lambda,penalty.factor=rep(1,p))
		returnVal
	}
}

rq.group.pen <- function(x,y, tau=.5,group=1:ncol(X), penalty=c("gLasso","gAdLasso","gSCAD","gMCP"),lambda=NULL,nlambda=100,eps=.0001,alg=c("huber","lp","qicd"), a=NULL, norm=2, group.pen.factor=rep(1,length(groups)), tau.pen=FALSE, ...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	g <- max(unique(group))
	nt <- length(tau)
	lpf <- length(group.pen.factor)
	pfmat <- FALSE
	if(max(group==1:ncol(X))==1){
		warning("p groups for p predictors, not really using a group penalty")
	}
	penalty <- match.arg(penalty)
	alg <- match.arg(alg)
	if(penalty=="gLasso" & norm==1){
		stop("Group Lasso with composite norm of 1 is the same as regular lasso, use norm = 2 if you want group lasso")
	}
	if(norm == 2 & alg != "huber"){
		alg <- "huber"
		warning("algorithm switched to huber, which is the only option for 2-norm")
	}
	if(penalty=="gAdLasso" & alg != "huber"){
		warning("huber algorithm used to derive ridge regression initial estimates for adaptive lasso. Second stage of algorithm used lp")
	}
	if(penalty=="gAdLasso" & alg == "qicd"){
		warning("No qicd algorithm for adaptive lasso, so switched to huber. If lp is used it will only be for the second stage.")
	}
	if(sum(tau <= 0 | tau >=1)>0){
		stop("tau needs to be between 0 and 1")
	}
	if(sum(penalty.factor<0)>0){
		stop("penalty factors must be positive")
	}
	if(nt==1){
		if(lpf!=g){
			stop("group penalty factor must be of length g")
		}
		if(tau.pen){
			penalty.factor <- penalty.factor*sqrt(tau*(1-tau))
			warning("tau.pen set to true for a single value of tau should return simliar answers as tau.pen set to false. Makes more sense to use it when fitting multiple taus at the same time")
		}
	} else{
		if(length(penalty.factor) == g*nt){
			pfmat <- TRUE
		}
		else if(max(penalty.factor) !=g){
			stop("penalty factor must be a vector with g groups or a nt by g matrix, where nt is the number of taus and g is the number of groups")
		}
		if(tau.pen){
			pfmat <- TRUE
			if(length(penalty.factor)==p){
				penalty.factor <- matrix(rep(penalty.factor,nt),nrow=nt,byrow=TRUE)
			} else{
				penalty.factor <- penalty.factor*sqrt(tau*(1-tau))
			}
		}
		
	}
	if(is.null(a))){
		if(penalty=="gLasso"){
			a <- 1
		} else if( penalty == "gAdLasso"){
			a <- 1
		} else if (penalty == "gSCAD"){ 
			a <- 3.7
		} else if (penalty == "gMCP"){
			a <- 3
		}
	}
	if(is.null(lambda)){
		lamMax <- getLamMax(x,y,tau,penalty)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	if(norm == 1){
		if(penalty == "gAdLasso"){
			init.model <- rq.enet(x,y,tau,...)
		} else{
			if(alg == "qicd"){
				init.alg <- "lp"
			} else{
				init.alg <- alg
			}
			init.model <- rq.lasso(x,y,tau,alg=init.alg,lambda=lambda,...)
		}
		#then figure out how to get group derivative for each coefficient. I might have some code that does that already. 
	}	
}

rq.lasso.modelreturn <- function(coefs,x,y,tau,lambda,penalty.factor){
	return_val <- NULL
	return_val$coefficients <- coefs
	return_val$lambda <- lambda
	return_val$penalty.factor <- penalty.factor
	if(is.null(colnames(x))){
		x_names <- paste("x",1:p,sep="")
	} else{
		x_names <- colnames(x)
	}
	x_names <- c("intercept",x_names)
	
	rownames(return_val$coefficients) <- x_names
	#need to think about how the fits will be, along with the rest. Maybe should I be transposing matrix. Maybe check code to see how other betas are done. 
	fits <- cbind(1,x)%*% return_val$coefficients
	return_val$fitted <- fits
	return_val$residuals <- y - fits
	return_val$rho <- apply(check(return_val$residuals,tau),2,mean)
	penalty.factor <- c(0,penalty.factor)
	return_val$PenRho <- return_val$rho + apply(penalty.factor*abs(return_val$coefficients),2,sum)
	return_val$tau <- tau
	return_val$df <- apply(return_val$coefficients!=0,2,sum)
	return_val
} 



rq.lasso.huber.onetau <- function(x,y,tau,lambda,penalty.factor=rep(1,ncol(x)),scalex=TRUE,...){
	dims <- dim(x)
	p <- dims[2]
	if(scalex){
		hqModel <- hqreg(x,y,method="quantile",tau=tau,lambda=lambda,penalty.factor=penalty.factor,...)
	} else{
		hqModel <- hqreg_raw(x,y,method="quantile",tau=tau,lambda=lambda,penalty.factor=penalty.factor,...)
	}
	rq.lasso.modelreturn(hqModel$beta,x,y,tau,lambda,penalty.factor)
}

rq.lasso.huber <- function(x,y,tau,lambda,penalty.factor=rep(1,ncol(x)),scalex=TRUE,pfmat=FALSE,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	
	if(length(tau)==1){		
		models <- rq.lasso.huber.onetau(x,y,tau=tau,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,...)
	} else{
		penf <- penalty.factor
		models <- list()
		for(i in 1:nt){
			subtau <- tau[i]
			if(pfmat){
				penf <- penalty.factor[i,]
			}
			models[[i]] <- rq.lasso.huber.onetau(x,y,tau=subtau,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,...)
		}
		attributes(models)$names <- paste0("tau",tau)
	}
	returnVal <- list(models=models, n=n, p=p,alg="huber",tau=tau,lambda=lambda,penalty.factor=penalty.factor)
	returnVal
}


print.rq.lasso <- function(x,...){
    if(x$nt==1){
		print(data.frame(df=x$models$df,lambda=x$models$lambda))
	} else{
		print(c("Quantile regression with lasso penalty for quantiles:",x$tau))
	}
}

getModelCoefs <- function(x,index){
	coefficients(x)[,index]
}

coef.rq.pen.seq <- function(x,index=NULL){
	lt <- length(x$tau)
	if(lt==1){
		if(is.null(index)){
			returnVal <- coefficients(x$models)
		} else{
			returnVal <- coefficients(x$models)[,index]
		}
	} else{
		if(is.null(index)){
			returnVal <- lapply(x$models,coef)
		} else if(length(index) == 1){
			returnVal <- sapply(x$models,getModelCoefs,index)
		} else if(lt != length(index)){
			stop("index must be one value or one value for each tau")
		} else{
			returnVal <- NULL
			for(i in 1:lt){
				returnVal <- cbind(returnVal,coefficients(x$models[[i]])[,index[i]]) 
			}
			colnames(returnVal) <- x$tau
		}
	}
	returnVal
}




rq.lasso.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                         coef.cutoff=1e-08, method="br",penVars=NULL,scalex=TRUE, ...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# lambda takes values of 1 or p
# coef.cutoff is a threshold to set to zero. 
# Choose the method used to estimate the coefficients ("br", "fn" or any other method used by quantreg)
### According to quantreg manual and my experience, "fn" is much faster for big n
### The "n" can grow rapidly using lin. prog. approach  
# penVars - variables to be penalized, doesn't work if lambda has multiple entries (Ben: I think it does though it is a little bit strange to do)
   if(is.null(dim(x))){
      stop('x needs to be a matrix with more than 1 column')
   }
   p <- dim(x)[2]
   if(p == 1){
	  stop('x needs to be a matrix with more than 1 column')
   }
   n <- dim(x)[1]
   if(n != length(y)){
      stop('length of y and rows of x do not match')
   }
   if(is.null(lambda)==TRUE | (length(lambda) != 1 & length(lambda) != dim(x)[2])){
      stop(paste('input of lambda must be of length 1 or', dim(x)[2]))
   }
   if( sum(lambda < 0) > 0){
      stop(paste('lambda must be positive and we have a lambda of ', lambda, sep=""))
   }
   if(scalex){
	  original_x <- x
	  x <- scale(x)
	  mu_x <- attributes(x)$`scaled:center`
	  sigma_x <- attributes(x)$`scaled:scale`
   }

   if(is.null(penVars) !=TRUE){# & length(lambda) == 1){
      if(length(lambda)==1){
		  mult_lambda <- rep(0,p)
		  mult_lambda[penVars] <- lambda
		  lambda <- mult_lambda
	  } else{
		lambda[-penVars] <- 0
	  }
   }
   lambda <- lambda*n # need this to account for the fact that rq does not normalize the objective function
   if(length(lambda)==1){
      pen_x <- rbind(diag(rep(lambda,p)),diag(rep(-lambda,p)))
   } else{
      pen_x <- rbind(diag(lambda), diag(-lambda))
      pen_x <- pen_x[rowSums(pen_x==0)!=dim(pen_x)[2],]#drop all zero rows
   }
   aug_n <- dim(pen_x)[1]
   aug_x <- rbind(x,pen_x)
   if(intercept){
      aug_x <- cbind(c(rep(1,n),rep(0,aug_n)), aug_x)
   }
   aug_y <- c(y, rep(0,aug_n))
   if(is.null(weights)){
     model <- rq(aug_y ~ aug_x+0, tau=tau, method=method)
   } else{
     if(length(weights) != n){
       stop("Length of weights does not match length of y")
     }
     orig_weights <- weights
     weights <- c(weights, rep(1,aug_n))
     model <- rq(aug_y ~ aug_x+0, tau=tau, weights=weights, method=method)
   }
   p_star <- p+intercept
   coefs <- coefficients(model)[1:p_star]
   return_val <- NULL
   return_val$coefficients <- coefs
   if(is.null(colnames(x))){
     x_names <- paste("x",1:p,sep="")
   } else{
     x_names <- colnames(x)
   }
   if(intercept){
     x_names <- c("intercept",x_names)
   }
   attributes(return_val$coefficients)$names <- x_names
   return_val$coefficients[abs(return_val$coefficients) < coef.cutoff] <- 0
   if(scalex){
   #need to update for penVars
	 return_val$coefficients <- transform_coefs(return_val$coefficients,mu_x,sigma_x,intercept)
	 if(intercept){
		fits <- cbind(1,original_x) %*% return_val$coefficients
	 } else{
		fits <- original_x %*% return_val$coefficients
	 }
	 return_val$residuals <- y - fits
	 return_val$PenRho <- sum(sapply(return_val$residuals,check,tau))+get_coef_pen(return_val$coefficients,lambda,intercept,penVars)	 
   } else{
	 return_val$PenRho <- model$rho
	 return_val$residuals <- model$residuals[1:n]   
   }
   if(is.null(weights)){   
     return_val$rho <- sum(sapply(return_val$residuals,check,tau))
   } else{
     return_val$rho <- sum(orig_weights*sapply(return_val$residuals,check,tau))
   }
   return_val$tau <- tau
   return_val$n <- n                  
   return_val$intercept <- intercept
   class(return_val) <- c("rq.pen", "rqLASSO")
   return_val
}

predict.rq.pen <- function(object, newx,...){
  coefs <- object$coefficients
  if(object$intercept){
     newx <- cbind(1,newx)
  }
  newx %*% coefs
}

predict.cv.rq.pen <- function(object, newx, lambda="lambda.min",...){
  if(lambda == "lambda.min"){
     target_pos <- which(object$cv$lambda == object$lambda.min)
  } else{
     target_pos <- which(object$cv$lambda == lambda)
  }
  predict(object$models[[target_pos]],newx,...)
}

getRho <- function(model){
    model$rho
}

coef.cv.rq.group.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  object$beta[,target_model]
}

group_derivs <- function(deriv_func,groups,coefs,lambda,a=3.7){
   if(length(lambda)==1){
      lambda <- rep(lambda,length(groups))
   }
   derivs <- NULL
   for(g in 1:length(unique(groups))){
      g_index <- which(groups==g)
      current_lambda <- lambda[g]
      coefs_l1 <- sum(abs(coefs[g_index]))
      derivs <- c(derivs, deriv_func(coefs_l1,current_lambda,a))
   }
   derivs
}

rq.group.lin.prog <- function(x,y,groups,tau,lambda,intercept=TRUE,eps=1e-05,penalty="SCAD", a=3.7, coef.cutoff=1e-08,
                                initial_beta=NULL,iterations=1,converge_criteria=.0001,penGroups=NULL,...){
    group_num <- length(unique(groups))
    if(length(lambda) == 1){
       lambda <- rep(lambda,group_num)
    }
    if (length(lambda) != group_num) {
        stop("lambdas do not match with group number")
    }
    if (sum(groups == 0) > 0) {
        stop("0 cannot be used as a group")
    }
    if (dim(x)[2] != length(groups)) {
        stop("length of groups must be equal to number of columns in x")
    }
    if (penalty == "SCAD") {
        deriv_func <- scad_deriv
    }
    if (penalty == "MCP") {
        deriv_func <- mcp_deriv
    } 
    
    new_lambda <- NULL
    group_count <- xtabs(~groups)
    for (g in 1:group_num) {
        new_lambda <- c(new_lambda, rep(lambda[g], each = group_count[g]))
    }
	if(is.null(penGroups)==FALSE){
		zero_pen_spots <- which(!groups %in% penGroups)
		new_lambda[zero_pen_spots] <- 0
	}
    if(is.null(initial_beta)){
       initial_beta <- rq.lasso.fit(x,y,tau,new_lambda, intercept=intercept, coef.cutoff=coef.cutoff,...)$coefficients
    }
    
    coef_by_group_deriv <- group_derivs(deriv_func, groups, initial_beta,lambda,a)
    lambda_update <- coef_by_group_deriv[groups]
    old_beta <- initial_beta
        
    iter_complete <- FALSE
    iter_num <- 0
    
    #pen_range <- (1+intercept):(dim(x)[2]+intercept)
    coef_range <- (1+intercept):(dim(x)[2]+intercept)
    
    while(!iter_complete){
      sub_fit <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda_update,intercept=intercept,...)
      coef_by_group_deriv <- group_derivs(deriv_func,groups,sub_fit$coefficients[coef_range],lambda,a)
      lambda_update <- coef_by_group_deriv[groups]
	  if(is.null(penGroups)==FALSE){
		lambda_update[zero_pen_spots] <- 0
	  }
      iter_num <- 1
      new_beta <- sub_fit$coefficients
      beta_diff <- sum( (old_beta - new_beta)^2)
      if(iter_num == iterations | beta_diff < converge_criteria){
        iter_complete <- TRUE
        if(iter_num == iterations & beta_diff > converge_criteria){
          warning(paste("did not converge after ", iterations, " iterations", sep=""))
        }
      } else{
        old_beta <- new_beta
      }
    }
    sub_fit$penalty <- penalty
    class(sub_fit) <-  c("rq.pen", "rqNC")
    sub_fit
}



groupMultLambda <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, penalty="LASSO", 
    #initial_beta = NULL,
    alg="QICD_warm",penGroups=NULL, ...) 
{
    if(alg != "QICD_warm"){
		#don't know how to do warm start with linear programming approach 
		return_val <- list()
		pos <- 1
		for (lam in lambda) {
			return_val[[pos]] <- rq.group.fit(x = x, y = y, groups = groups, 
				tau = tau, lambda = lam, intercept = intercept, penalty=penalty,alg=alg, penGroups=penGroups,
				...)
			#initial_beta <- return_val[[pos]]$coefficients
			pos <- pos + 1
		}
	} else{
		p <- dim(x)[2]
		pos <- 1
		alg = "QICD"
		
		return_val <- list()
		if(intercept){
			initial_beta <- c(quantile(y,tau), rep(0,p))
		} else{
			initial_beta <- rep(0,p)
		}
		
		for(lam in lambda){
			return_val[[pos]] <- rq.group.fit(x=x, y=y, groups=groups, tau=tau, lambda= lam, intercept=intercept, penalty="LASSO", alg=alg, initial_beta=initial_beta, penGroups=penGroups, ...)
			initial_beta <- coefficients(return_val[[pos]])
			pos <- pos + 1
		}
		
		#if penalty is not lasso then update those initial estimates
		if(penalty != "LASSO"){
			pos <- 1
			for(lam in lambda){
				initial_beta <- coefficients(return_val[[pos]]) #use lasso estimate as initial estimate
				return_val[[pos]] <- rq.group.fit(x=x, y=y, groups=groups, tau=tau, lambda= lam, intercept=intercept, penalty=penalty, alg=alg, initial_beta=initial_beta, penGroups=penGroups, ...)
				pos <- pos + 1
			}
		}
	}
    return_val
}

nonzero <- function (obj) 
{
    UseMethod("nonzero")
}

nonzero.cv.rq.group.pen <- function (obj) 
{
    coefs <- coefficients(obj)
    if (obj$intercept) {
		coefs <- coefs[-1]
    }
    tapply(coefs, obj$groups, sum) != 0
}

print.cv.rq.pen <- function(x,...){
   cat("\nCoefficients:\n")
   print(coefficients(x,...))
   cat("\nCross Validation (or BIC) Results\n")
   print(x$cv)
}

print.rq.pen <- function(x,...){
    cat("\nCoefficients:\n")
	print(coefficients(x,...))
}