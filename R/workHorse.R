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

lasso <- function(x,lambda=1,a=1){
   lambda*a*abs(x)
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

getLamMax <- function(x,y,tau=.5,gamma=.2,gamma.max=4,gamma.q=.1,penalty="LASSO",scalex=TRUE){
	n <- length(y)
	returnVal <- 0
	if(scalex){
		x <- scale(x)
	}
	
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

getLamMaxGroup <- function(x,y,group.index,tau=.5,group.pen.factor,gamma=.2,gamma.max=4,gamma.q=.1,penalty="gLASSO",scalex=TRUE){
	returnVal <- 0
	n <- length(y)
	if(scalex){
		x <- scale(x)
	}
	for(tau_val in tau){
		r <- y - quantile(y,tau_val)
		gamma0<- min(gamma.max, max(gamma, quantile(abs(r), probs = gamma.q)))

		grad_k<- -neg.gradient(r, rep(1,n), tau_val, gamma=gamma0, x, apprx="huber")
		grad_k.norm<- tapply(grad_k, group.index, hrqglas:::l2norm)
  
		lambda.max<- max(c(returnVal,grad_k.norm/group.pen.factor))
	}
	lambda.max
}

# If tau.pen is set to true then the reported lambdas are actually lambda*sqrt(tau*(1-tau))
rq.lasso <- function(x,y,tau=.5,lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001), penalty.factor = rep(1, ncol(x)),
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
		lamMax <- getLamMax(x,y,tau,scalex=scalex)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	if(alg=="huber"){
		if(length(lambda)==1){
			stop("The Huber algorithm requires at least 2 values of lambda")
		}
		returnVal <- rq.lasso.huber(x,y,tau,lambda,penalty.factor,scalex,pfmat,...)
	} else{
		models <- vector(mode="list",length=nt)
		modelnames <- NULL
		for(i in 1:nt){
			coefs <- NULL
			j <- 1
			for(lam in lambda){
				if(pfmat){
					sublam <- lam*penalty.factor[j,]
				} else{
					sublam <- lam*penalty.factor
				}
				subm <- rq.lasso.fit(x,y,tau[i],lambda=sublam, method=alg,scalex=scalex, ...)
				coefs <- cbind(coefs,coefficients(subm))
				j <- j + 1
			}
			models[[i]] <- rq.pen.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor,"LASSO",1)
			modelnames <- c(modelnames,paste0("tau",tau[i],"a",1))
		}
		names(models) <- modelnames
		returnVal <- list(models=models, n=n, p=p,alg=alg,tau=tau, a=1)
	}
	modelIndex <- 1:length(returnVal$models)
	avals <- sapply(returnVal$models,modelA)
	tauvals <- sapply(returnVal$models,modelTau)
	modelsInfo <- data.frame(modelIndex=modelIndex,a=avals,tau=tauvals)
	returnVal$modelsInfo <- modelsInfo
	returnVal$penalty <- "LASSO"
	class(returnVal) <- "rq.pen.seq"
	returnVal
}

rq.enet <- function(x,y,tau=.5,lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001), penalty.factor = rep(1, ncol(x)),scalex=TRUE,tau.pen=FALSE,a=0,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	lpf <- length(penalty.factor)
	pfmat <- FALSE
	
	if(sum(tau <= 0 | tau >=1)>0){
		stop("tau needs to be between 0 and 1")
	}
	if(sum(a < 0 | a > 1) > 0){
		stop("a needs to be >= 0 and <= 1")
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
		lamMax <- getLamMax(x,y,tau,scalex=scalex)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	if(length(lambda)==1){
			stop("The Huber algorithm requires at least 2 values of lambda and elastic net only uses the Huber algorithm")
	}
	returnVal <- rq.lasso.huber(x,y,tau,lambda,penalty.factor,scalex,pfmat,a=a,...)
	avals <- sapply(returnVal$models,modelA)
	tauvals <- sapply(returnVal$models,modelTau)
	modelsInfo <- data.frame(modelIndex=1:length(returnVal$models),a=avals,tau=tauvals)
	returnVal$modelsInfo <- modelsInfo
	returnVal$penalty <- "ENet"
	returnVal$a <- a
	
	class(returnVal) <- "rq.pen.seq"
	returnVal	
}


rq.lla <- function(obj,x,y,penalty="SCAD",a=ifelse(penalty=="SCAD",3.7,3),...){
	nt <- length(obj$tau)
	na <- length(a)
	if(penalty=="SCAD"){
		derivf <- scad_deriv
	} else if(penalty=="MCP"){
		derivf <- mcp_deriv
	} else if(penalty=="aLASSO"){
		derivf <- alasso_wt
	} else{
		stop("Penalty must be SCAD, MCP or aLASSO")
	}
	
	newModels <- vector(mode="list",length=nt*na)
	pos <- 1
	modelNames <- NULL
	for(j in 1:nt){
		lampen <- as.numeric(obj$models[[j]]$penalty.factor %*% t(obj$models[[j]]$lambda))
		ll <- length(obj$models[[j]]$lambda)
		for(k in 1:na){	
			pfs <- matrix(derivf(as.numeric(abs(coefficients(obj$models[[j]])[-1,])),lampen,a=a[k]),ncol=ll)
			newModels[[pos]] <- obj$models[[j]]
			for(i in 1:ll){
				if(obj$alg=="huber"){
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=c(2,1),penalty.factor=pfs[,i],alg=obj$alg,...)$models[[1]])[,2]

				} else{
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=1,penalty.factor=pfs[,i],alg=obj$alg,...)$models[[1]])
				}
				newModels[[pos]]$coefficients[,i] <- update_est
			}
			newModels[[pos]]$a <- a[k]
			modelNames <- c(modelNames,paste0("tau",obj$tau[j],"a",a[k]))
			pos <- pos+1
		}
	}
	names(newModels) <- modelNames
	obj$models <- newModels	
	if(penalty=="aLASSO"){
		obj$penalty.factor <- pfs
	}
	obj$a <- a
	obj$penalty <- penalty
	obj$modelsInfo <- createModelsInfo(obj$models)
	obj
}

clearModels <- function(model,pos){
	endpos <- pos - 1
	model$coefficients <- model$coefficients[,1:endpos]
	model$lambda <- model$lambda[1:endpos]
	model$rho <- model$rho[1:endpos]
	model$PenRho <- model$PenRho[1:endpos]
	model$nzero <- model$nzero[1:endpos]
	model
}

rq.group.lla <- function(obj,x,y,groups,penalty=c("gAdLASSO","gSCAD","gMCP"),a=NULL,norm=2, group.pen.factor,...){
	#for loop calculation of penalty factors that could maybe be removed
	nt <- length(obj$tau)
	p <- ncol(x)
	g <- max(groups)
	penalty <- match.arg(penalty)
	obj$penalty <- penalty
	derivf <- getDerivF(penalty)
	na <- length(a)
	if(norm !=1 & norm != 2){
		stop("Norm needs to be set to 1 or 2.")
	}
	if(norm ==2 & obj$alg != "huber"){
		stop("Huber approximation is the only available algorithm for 2-norm based penalties.")
	}
	gpfmat <- NULL
	newModels <- vector(mode="list",length=nt*na)
	pos <- 1
	for(j in 1:nt){
		lampen <- group.pen.factor %*% t(obj$models[[j]]$lambda)
		ll <- length(obj$models[[j]]$lambda)
		for(k in 1:na){
			newModels[[pos]] <- obj$models[[j]]	
			for(i in 1:ll){	
				coef_by_group_deriv <- group_derivs(derivf, groups, abs(coefficients(obj$models[[j]])[-1,i]),lampen[,i],a[k],norm=norm)
				if(sum(coef_by_group_deriv)==0){
					newModels[[pos]] <- clearModels(newModels[[pos]],i)
					break
				} else{
					if(penalty == "gAdLASSO"){
						gpfmat <- cbind(gpfmat,coef_by_group_deriv)
					}
					if(obj$alg=="huber"){
						if(norm == 1){
							penalty.factor <- mapvalues(groups,seq(1,g),coef_by_group_deriv)
							update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=c(2,1),penalty.factor=penalty.factor,alg=obj$alg,...)$models[[1]])[,2]
						} else{
							update_est <- coefficients(rq.group.pen(x,y,obj$tau[j],groups,lambda=1,group.pen.factor=coef_by_group_deriv, alg=obj$alg,...)$models[[1]])
						}
					} else{
						penalty.factor <- mapvalues(groups,seq(1,g),coef_by_group_deriv)
						update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=1,penalty.factor=penalty.factor,alg=obj$alg,...)$models[[1]])
					}
					newModels[[pos]]$coefficients[,i] <- update_est
				}
			}
			newModels[[pos]] <- rq.pen.modelreturn(newModels[[pos]]$coefficients,x,y,obj$tau[j],newModels[[pos]]$lambda,rep(1,p),penalty,a[k])	
			newModels[[pos]]$penalty.factor <- NULL			
			if(penalty == "gAdLASSO"){
				newModels[[pos]]$group.pen.factor <- gpfmat 
			} else{
				newModels[[pos]]$group.pen.factor <- group.pen.factor
			}
			dimnames(newModels[[pos]]$group.pen.factor) <- NULL			
			pos <- pos + 1
		}
	}
	obj$models <- newModels
	obj$a <- a
	obj  <- updateGroupPenRho(obj,norm,groups)
	obj$groups <- groups
	obj$penalty <- penalty
	#obj$class <- c(obj$class, "rq.group.pen.seq")
	obj
}


getA <- function(a,penalty){
	if(penalty=="aLASSO" | penalty == "gAdLASSO"){
		if(is.null(a)){
			a <- 1
		} else if(sum(a < 0)>0){
			stop("for adaptive lasso the tuning parameter \"a\" must be positive")
		}
	} else{
		if(is.null(a)){
			if(penalty == "SCAD" | penalty == "gSCAD"){
				a <- 3.7
				penalty <- "SCAD"
			} 
			if(penalty == "MCP" | penalty == "gMCP"){
				a <- 3
				penalty <- "MCP"
			}
		} else{
			if(penalty == "SCAD" & sum(a <= 2) > 0 ){
				stop("Tuning parameter \"a\" must be larger than 2 for SCAD")
			}
			if(penalty == "MCP" & sum(a <= 1) > 0){
				stop("Tuning parameter \"a\" must be larger than 1 for MCP")
			}
		}
	}
	a
}

createModelsInfo <- function(models){
	modelIndex <- 1:length(models)
	avals <- sapply(models,modelA)
	tauvals <- sapply(models,modelTau)
	modelsInfo <- data.frame(modelIndex=modelIndex,a=avals,tau=tauvals)
	modelsInfo
}

rq.nc <- function(x, y, tau=.5,  penalty=c("SCAD","aLASSO","MCP"),a=NULL,lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001),alg="huber",scalex=TRUE,...) {
	#should look at how ncvreg generates the lambda sequence and combine that with the Huber based approach
	penalty <- match.arg(penalty)
	nt <- length(tau)
	dims <- dim(x)
	p <- dims[2]
	n <- dims[1] 
	if( max(tau) >= 1 | min(tau) <= 0){
		stop("Values for tau must be between 0 and 1") 
	}
	if(penalty == "aLASSO" & alg != "huber"){
		stop("Only algorithm for adaptive lasso is huber.")
	}
	if(is.null(lambda)){
		lamMax <- getLamMax(x,y,tau,penalty=penalty,scalex=scalex)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	a <- getA(a,penalty)
	na <- length(a)
	
	if(alg != "qicd" & alg != "QICD"){
	#QICD implementation not provided in rq.lasso
		if(penalty=="aLASSO"){
			init.model <- rq.enet(x,y,tau,...)
		} else{
			init.model <- rq.lasso(x,y,tau,alg=alg,lambda=lambda,...)
		}
		rq.lla(init.model,x,y,penalty,a,...)
	} else{
		pos <- 1
		models <- vector(mode="list",length=nt*na)
		modelNames <- NULL
		for(i in 1:nt){
			for(j in 1:na){
				coefs <- NULL
				for(lam in lambda){
					sublam <- lam
					subm <- rq.nc.fit(x,y,tau[i],lambda=sublam, alg="QICD", a=a[j], ...)
					coefs <- cbind(coefs,coefficients(subm))
				}
				models[[pos]] <- rq.pen.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor=rep(1,p),penalty,a[j])
				modelNames <- c(modelNames,paste0("tau",tau[i],"a",a[j]))
				pos <- pos + 1
			}
		}
		names(models) <- modelNames
		returnVal <- list(models=models, n=n, p=p,alg=alg,tau=tau,lambda=lambda,penalty.factor=rep(1,p),penalty=penalty,a=a)		
		returnVal$modelsInfo <- createModelsInfo(models)
		class(returnVal) <- "rq.pen.seq"
		returnVal
	}
}

group_derivs <- function(deriv_func,groups,coefs,lambda,a=3.7,norm=1){
   if(length(lambda)==1){
      lambda <- rep(lambda,length(groups))
   }
   derivs <- NULL
   for(g in 1:length(unique(groups))){
      g_index <- which(groups==g)
      current_lambda <- lambda[g]
	  if(norm == 1){
		gnorm <- sum(abs(coefs[g_index]))
	  } else if(norm==2){
		gnorm <- sqrt(sum(coefs[g_index]^2))
	  } else{
		stop("Invalid norm, use 1 or 2")
	  }
      derivs <- c(derivs, deriv_func(gnorm,current_lambda,a))
   }
   derivs
}

getDerivF <- function(penalty){
	if(penalty=="SCAD" | penalty=="gSCAD"){
		derivf <- scad_deriv
	} else if(penalty=="MCP" | penalty == "gMCP"){
		derivf <- mcp_deriv
	} else if(penalty=="aLASSO" | penalty == "gAdLASSO"){
		derivf <- alasso_wt
	}
	derivf
}

rq.glasso <- function(x,y,tau,groups, lambda, group.pen.factor,pfmat,scalex,lambda.discard=FALSE,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	g <- length(unique(groups))
	nt <- length(tau)
	
	penf <- group.pen.factor
	models <- vector(mode="list",length=nt)
	
	for(i in 1:nt){
		subtau <- tau[i]
		if(pfmat){
			penf <- group.pen.factor[i,]
		}
		models[[i]] <- hrq_glasso(x,y,group.index=groups,tau=subtau,lambda=lambda,w.lambda=penf,scalex=scalex,lambda.discard=lambda.discard,...)
		models[[i]] <- rq.pen.modelreturn(models[[i]]$beta,x,y,subtau,models[[i]]$lambda,rep(1,p),"gLASSO",1)
		models[[i]]$penalty.factor <- NULL
		models[[i]]$group.pen.factor <- penf
	}
	attributes(models)$names <- paste0("tau",tau)
		
	returnVal <- list(models=models, n=n, p=p,alg="huber",tau=tau,penalty="gLASSO",a=1)
	returnVal <- updateGroupPenRho(returnVal,2,groups)
	returnVal$a <- 1
	returnVal$modelsInfo <- data.frame(modelIndex=1:nt,a=1,tau=tau)
	returnVal
}

rq.group.pen <- function(x,y, tau=.5,groups=1:ncol(x), penalty=c("gLASSO","gAdLASSO","gSCAD","gMCP"),lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001),alg=c("huber","lp","qicd"), a=NULL, norm=2, group.pen.factor=rep(1,length(unique(groups))),tau.pen=FALSE,scalex=TRUE, ...){
	penalty <- match.arg(penalty)
	alg <- match.arg(alg)
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	g <- length(unique(groups))
	nt <- length(tau)
	na <- length(a)
	lpf <- length(group.pen.factor)
	if(penalty != "gLASSO"){
		a <- getA(a,penalty)
	} else{
		if(is.null(a)==FALSE){
			warning("The tuning parameter a is not used for group lasso")
		}
	}
	pfmat <- FALSE
	if(g==p){
		warning("p groups for p predictors, not really using a group penalty")
	}
	penalty <- match.arg(penalty)
	alg <- match.arg(alg)
	if(norm != 1 & norm != 2){
		stop("norm must be 1 or 2")
	}
	if(penalty=="gLASSO" & norm==1){
		stop("Group Lasso with composite norm of 1 is the same as regular lasso, use norm = 2 if you want group lasso")
	}
	if(norm == 1 & penalty == "gAdLASSO"){
		warning("Group adapative lasso with 1 norm results in a lasso estimator where lambda weights are the same for each coefficient in a group. However, it does not force groupwise sparsity, there can be zero and non-zero coefficients within a group.")
	}
	if(norm == 2 & alg != "huber"){
		stop("If setting norm = 2 then algorithm must be huber")
	}
	if(penalty=="gAdLASSO" & alg != "huber"){
		warning("huber algorithm used to derive ridge regression initial estimates for adaptive lasso. Second stage of algorithm used lp")
	}
	if(penalty=="gAdLASSO" & alg == "qicd"){
		stop("No qicd algorithm for adaptive lasso.")
	}
	if(sum(tau <= 0 | tau >=1)>0){
		stop("tau needs to be between 0 and 1")
	}
	if(sum(group.pen.factor<0)>0){
		stop("penalty factors must be positive")
	}
	if(nt==1){
		if(lpf!=g){
			stop("group penalty factor must be of length g")
		}
		if(tau.pen){
			penalty.factor <- group.pen.factor*sqrt(tau*(1-tau))
			warning("tau.pen set to true for a single value of tau should return simliar answers as tau.pen set to false. Makes more sense to use it when fitting multiple taus at the same time")
		}
	} else{
		if(length(group.pen.factor) == g*nt){
			pfmat <- TRUE
		}
		else if(length(group.pen.factor) !=g){
			stop("penalty factor must be a vector with g groups or a nt by g matrix, where nt is the number of taus and g is the number of groups")
		}
		if(tau.pen){
			pfmat <- TRUE
			if(length(group.pen.factor)==p){
				group.pen.factor <- matrix(rep(group.pen.factor,nt),nrow=nt,byrow=TRUE)
			} else{
				group.pen.factor <- group.pen.factor*sqrt(tau*(1-tau))
			}
		}
		
	}
	if(is.null(lambda)){
		lamMax <- getLamMaxGroup(x,y,groups,tau,group.pen.factor,penalty=penalty,scalex)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
			
	if(pfmat){
		#maybe some applies mapvalues
	} else{
		penalty.factor <- mapvalues(groups,seq(1,g),group.pen.factor)
	}
	
	if(norm == 1){
		if(penalty == "gAdLASSO"){
			init.model <- rq.enet(x,y,tau,lambda=lambda,penalty.factor=penalty.factor,...)
		} else{
			if(alg == "qicd"){
				init.alg <- "lp"
			} else{
				init.alg <- alg
			}
			init.model <- rq.lasso(x,y,tau,alg=init.alg,lambda=lambda,tau.pen=FALSE,penalty.factor=penalty.factor,scalex=scalex,...)
		}
		return_val <- rq.group.lla(init.model,x,y,groups,penalty=penalty,a=a,norm=norm,group.pen.factor=group.pen.factor,...)
	} else{
		if(penalty == "gLASSO"){
			return_val <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,pfmat,scalex,...)
		} else{
			if(penalty == "gAdLASSO"){
				init.model <- rq.enet(x,y,tau,lambda=lambda,penalty.factor=penalty.factor,...)
			} else{
				init.model <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,pfmat,scalex,...)
			}
			return_val <- rq.group.lla(init.model,x,y,groups,penalty=penalty,a=a,norm=norm,group.pen.factor=group.pen.factor,...) 
		}
	}
	class(return_val) <- "rq.pen.seq"
	return_val
}



rq.pen.modelreturn <- function(coefs,x,y,tau,lambda,penalty.factor,penalty,a){
# for loop that could be removed
	penfunc <- getPenfunc(penalty)
	return_val <- NULL
	return_val$coefficients <- coefs
	return_val$lambda <- lambda
	return_val$penalty.factor <- penalty.factor
	p <- ncol(x)
	if(is.null(colnames(x))){
		x_names <- paste("x",1:p,sep="")
	} else{
		x_names <- colnames(x)
	}
	x_names <- c("intercept",x_names)
	
	#need to think about how the fits will be, along with the rest. Maybe should I be transposing matrix. Maybe check code to see how other betas are done. 
	fits <- cbind(1,x)%*% return_val$coefficients
	#return_val$fitted <- fits
	res <- y - fits
	if(is.null(dim(return_val$coefficients))==TRUE){
		return_val$rho <- mean(check(res,tau))	
		return_val$PenRho <- return_val$rho + sum(penfunc(return_val$coefficients[-1],lambda*return_val$penalty.factor,a))
		return_val$nzero <- sum(return_val$coefficients!=0)
	} else{
		return_val$rho <- apply(check(res,tau),2,mean)	
		rownames(return_val$coefficients) <- x_names
		for(i in 1:length(return_val$rho)){
			return_val$PenRho[i] <- return_val$rho[i] + sum(penfunc(return_val$coefficients[-1,i],lambda[i]*return_val$penalty.factor,a))
		}
		return_val$nzero <- apply(return_val$coefficients!=0,2,sum)
	}
	return_val$tau <- tau
	return_val$a <- a
	return_val
}

group_derivs <- function(deriv_func,groups,coefs,lambda,a=3.7,norm=1){
   if(length(lambda)==1){
      lambda <- rep(lambda,length(groups))
   }
   derivs <- NULL
   for(g in 1:length(unique(groups))){
      g_index <- which(groups==g)
      current_lambda <- lambda[g]
	  if(norm == 1){
		gnorm <- sum(abs(coefs[g_index]))
	  } else if(norm==2){
		gnorm <- sqrt(sum(coefs[g_index]^2))
	  } else{
		stop("Invalid norm, use 1 or 2")
	  }
      derivs <- c(derivs, deriv_func(gnorm,current_lambda,a))
   }
   derivs
}

getDerivF <- function(penalty){
	if(penalty=="SCAD" | penalty=="gSCAD"){
		derivf <- scad_deriv
	} else if(penalty=="MCP" | penalty == "gMCP"){
		derivf <- mcp_deriv
	} else if(penalty=="aLASSO" | penalty == "gAdLASSO"){
		derivf <- alasso_wt
	}
	derivf
}

ridge <- function(x,lambda,a=1){
	a*lambda*square(x)
}

elastic <- function(x,lambda,a){
	a*lasso(x,lambda)+(1-a)*ridge(x,lambda)
}


getPenfunc <- function(penalty){
	if(penalty == "LASSO" | penalty == "gLASSO" | penalty=="aLASSO" | penalty=="gAdLASSO"){
	#adaptive lasso is lasso with different weights, so i think this will work
		penfunc <- lasso
	}
	else if(penalty=="SCAD" | penalty=="gSCAD"){
		penfunc <- scad
	} else if(penalty=="MCP" | penalty == "gMCP"){
		penfunc <- mcp
	} else if(penalty=="ENet"){
		penfunc <- elastic
	}
	penfunc
}



updateGroupPenRho <- function(obj,norm,groups){
	
	for(j in 1:length(obj$models)){
		a <- obj$models[[j]]$a
		if(length(obj$models[[j]]$lambda)==1){
			obj$models[[j]]$PenRho <- obj$models[[j]]$rho + sum(getGroupPen(obj$models[[j]]$coefficients[-1],groups,obj$models[[j]]$lambda,obj$models[[j]]$group.pen.factor,obj$penalty,norm,a)) 
		} else{
			for(i in 1:length(obj$models[[j]]$lambda)){
				if(obj$penalty=="gAdLASSO"){
					obj$models[[j]]$PenRho[i] <- obj$models[[j]]$rho[i] + sum(getGroupPen(obj$models[[j]]$coefficients[-1,i],groups,obj$models[[j]]$lambda[i],obj$models[[j]]$group.pen.factor[,i],obj$penalty,norm,a))
				} else{			
					obj$models[[j]]$PenRho[i] <- obj$models[[j]]$rho[i] + sum(getGroupPen(obj$models[[j]]$coefficients[-1,i],groups,obj$models[[j]]$lambda[i],obj$models[[j]]$group.pen.factor,obj$penalty,norm,a))
				}
			}
		}
	}	
	obj
}


getGroupPen <- function(coefs,groups,lambda,group.pen.factor,penalty,norm,a){
   lambda <- lambda*group.pen.factor
   penfunc <- getPenfunc(penalty)
   pens <- NULL
   for(g in 1:length(unique(groups))){
      g_index <- which(groups==g)
      current_lambda <- lambda[g]
	  if(norm == 1){
		gnorm <- sum(abs(coefs[g_index]))
	  } else if(norm==2){
		gnorm <- sqrt(sum(coefs[g_index]^2))
	  } else{
		stop("Invalid norm, use 1 or 2")
	  }
      pens <- c(pens, penfunc(gnorm,current_lambda,a))
   }
   pens	
}

rq.lasso.huber.onetau <- function(x,y,tau,lambda,penalty.factor=rep(1,ncol(x)),scalex=TRUE,a=1,...){
	dims <- dim(x)
	p <- dims[2]
	if(scalex){
		hqModel <- hqreg(x,y,method="quantile",tau=tau,lambda=lambda,penalty.factor=penalty.factor,alpha=a,...)
	} else{
		hqModel <- hqreg_raw(x,y,method="quantile",tau=tau,lambda=lambda,penalty.factor=penalty.factor,alpha=a,...)
	}
	rq.pen.modelreturn(hqModel$beta,x,y,tau,lambda,penalty.factor,"LASSO",a)
}

rq.lasso.huber <- function(x,y,tau,lambda,penalty.factor=rep(1,ncol(x)),scalex=TRUE,pfmat=FALSE,a=1,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	na <- length(a)
	
	penf <- penalty.factor
	models <- vector(mode="list",length=nt*na)
	pos <- 1
	modelnames <- NULL
	for(i in 1:nt){
		for(j in 1:na){
			subtau <- tau[i]
			if(pfmat){
				penf <- penalty.factor[i,]
			}
			models[[pos]] <- rq.lasso.huber.onetau(x,y,tau=subtau,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,a=a[j],...)
			pos <- pos+1
			modelnames <- c(modelnames,paste0("tau",tau[i],"a",a[j]))
		}
	}
	attributes(models)$names <- modelnames
		
	returnVal <- list(models=models, n=n, p=p,alg="huber",tau=tau,a=a)
	returnVal
}


print.rq.pen.seq <- function(x,...){
	nt <- length(x$tau)
	na <- length(x$a)
    if(nt==1 & na==1){
		print(data.frame(nzero=x$models[[1]]$nzero,lambda=x$models[[1]]$lambda))
	} else if(nt > 1 & na > 1){
		print(paste(c(paste(c("Quantile regression with ", x$penalty, " penalty for quantiles:",x$tau), collapse=" ")," and tuning parameters a:", x$a),collapse=" "))
	} else if( na > 1){
		print(paste(c(paste(c("Quantile regression with ", x$penalty, " penalty for quantile:",x$tau), collapse=" ")," and tuning parameters a:", x$a),collapse=" "))
	} else{
			print(paste(c("Quantile regression with ", x$penalty, " penalty for quantiles:",x$tau), collapse=" "))
	}	
}

getModelCoefs <- function(x,index){
	coefficients(x)[,index]
}

closeEnough <- function(targets,original){
#bad for loop	
	returnVals <- rep(0,length(original))
	for(targ in targets){
		returnVals <- returnVals + elementwise.all.equal(original,targ)
	}
	returnVals > 0 
}

whichMatch <- function(targets,original){
#bad for loop	
	returnVals <- NULL
	for(targ in targets){
		returnVals <- c(returnVals, which(elementwise.all.equal(original,targ))) 
	}
	returnVals
}

elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})


#I should get these to return a frame with a tau, and model index so I could extract that information easily
coef.rq.pen.seq <- function(x,tau=NULL,a=NULL,lambda=NULL,modelIndex=NULL,lambdaIndex=NULL){
	if( (is.null(tau)==FALSE | is.null(a)==FALSE) & is.null(modelIndex) == FALSE){
		stop("Set tau and a or set modelIndex, not both")
	}
	if( (is.null(lambda)==FALSE & is.null(lambdaIndex)==FALSE)){
		stop("Use lambda or lambdaIndex, not both")
	}
	lt <- length(x$tau)
	na <- length(x$a)	
	if((is.null(tau) == FALSE & is.null(a)==FALSE)){
		modelIndex <- intersect(whichMatch(tau,x$modelsInfo$tau),whichMatch(a,x$modelsInfo$a))
	} else if(is.null(tau)==FALSE){
		modelIndex <- whichMatch(tau,x$modelsInfo$tau)
	} else if(is.null(a) == FALSE){
		modelIndex <- whichMatch(a,x$modelsInfo$a)
	}
	else if(is.null(modelIndex)){
		modelIndex <- 1:length(x$models)
	}
	if(length(modelIndex)==0){
		stop("Invalid tau or a provided")
	}
	if(is.null(lambda)==FALSE){
		lambdaIndex <- whichMatch(lambda,x$models[[1]]$lambda)
	} else if(is.null(lambdaIndex)){
		lambdaIndex <- 1:length(x$models[[1]]$lambda)
	}
	if(length(lambdaIndex)==0){
		stop("Invalid lambda provided")
	}
	targetModels <- x$models[modelIndex]
	lapply(targetModels,getModelCoefs,lambdaIndex)
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
	 res <- y - fits
	 return_val$PenRho <- sum(sapply(res,check,tau))+get_coef_pen(return_val$coefficients,lambda,intercept,penVars)	 
   } else{
	 return_val$PenRho <- model$rho
	 res <- model$residuals[1:n]   
   }
   if(is.null(weights)){   
     return_val$rho <- sum(sapply(res,check,tau))
   } else{
     return_val$rho <- sum(orig_weights*sapply(res,check,tau))
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



rq.group.lin.prog <- function(x,y,groups,tau,lambda,intercept=TRUE,eps=1e-05,penalty="SCAD", a=3.7, coef.cutoff=1e-08,
                                initial_beta=NULL,iterations=1,converge_criteria=.0001,penGroups=NULL,...){
	#warning("Warning rq.group.lin.prog() is old and potentially inefficient. Recommend using rq.group.pen() instead")
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