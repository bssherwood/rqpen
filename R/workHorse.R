kernel_estimates <- function(x,y,h,...){
  kernel_estimates <- NULL
  if(is.null(dim(x))){
    n <- length(x)
    d <- 1
  } else{
    n <- dim(x)[1]
    d <- dim(x)[2]
  }
  for(i in 1:n){
     if(d == 1){
        kernel_estimates <- c(kernel_estimates,kernesti.regr(x[i],x,y,h=h,...))
     } else{
        kernel_estimates <- c(kernel_estimates,kernesti.regr(x[i,],x,y,h=h,...))
     }
  }
  kernel_estimates
}

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

model_eval <- function(model, test_x, test_y, test_w=NULL, func="check",...){
#func: "check" (Quantile Check), "SqErr" (Squared Error), "AE" (Absolute Value)
  if(model$intercept){
    test_x <- cbind(1,test_x)
  }
  fits <- test_x %*% coefficients(model)
  eval_func <- switch(which(c("check","SqErr","AE")==func), check, square, abs)
  if(is.null(test_w)){
	  mean(eval_func(test_y-fits,...)) 
  } else{
	  weighted.mean(eval_func(test_y-fits,...), test_w)
  }
}

qbic <- function(model, largeP=FALSE){
  tau <- model$tau
  n <- model$n
  nzero <- sum(model$coefficients !=0)
  if(largeP){
    bic <- log(model$rho) + nzero*log(n)*log(length(model$coefficients))/(2*n)
  }else{
    bic <- log(model$rho) + nzero*log(n)/(2*n)
  }
  bic
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

#randomly assign n samples into k groups
randomly_assign <- function(n,k){
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

#Old function that is no longer exported. Kept for reproducible reasons. 
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
  #print(length(coefs))
  #print(length(sigma_x))
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

rq.huber<- function(r, tau, gamma){
  r<- as.vector(r)
  (huber.loss(r,gamma)+(2*tau-1)*r)/2
}

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

#returns a lambda value that (likely) produces a completely sparse model. 
#Only likely and not guranteed, because it relies on huber approximation to the quantile loss
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

# Finds lambda max for a group penalty. 
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
#
rq.lasso <- function(x,y,tau=.5,lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001), penalty.factor = rep(1, ncol(x)),
						alg=ifelse(sum(dim(x))<200,"huber","br"),scalex=TRUE,tau.penalty.factor=rep(1,length(tau)),
						coef.cutoff=1e-8,max.iter=10000,converge.eps=1e-7,gamma=IQR(y)/10,...){
	if(alg == "lp"){
	#use br as the default for linear programming 
		alg <- "br"
	}
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	lpf <- length(penalty.factor)
	ltpf <- length(tau.penalty.factor)
	
	if(sum(tau <= 0 | tau >=1)>0){
		stop("tau needs to be between 0 and 1")
	}
	if(sum(penalty.factor<0)>0){
		stop("penalty factors must be positive")
	}
	if(lpf !=p){
		stop("penalty factor must be of length p")
	}
	if(ltpf !=nt){
		stop("tau penalty factor must be of length tau")
	}
	
	if(is.null(lambda)){
		lamMax <- getLamMax(x,y,tau,scalex=scalex)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	if(alg=="huber"){
		if(length(lambda)==1){
			stop("The Huber algorithm requires at least 2 values of lambda")
		}
		returnVal <- rq.lasso.huber(x,y,tau,lambda,penalty.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma=gamma,...)
	} else{
		models <- vector(mode="list",length=nt)
		modelnames <- NULL
		for(i in 1:nt){
			coefs <- NULL
			j <- 1
			for(lam in lambda){
				sublam <- lam*penalty.factor*tau.penalty.factor[i]
				subm <- rq.lasso.fit(x,y,tau[i],lambda=sublam, method=alg,scalex=scalex, coef.cutoff=coef.cutoff,...)
				coefs <- cbind(coefs,coefficients(subm))
				j <- j + 1
			}
			models[[i]] <- rq.pen.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor*tau.penalty.factor[i],"LASSO",1)
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
	returnVal$lambda <- lambda
	returnVal$penalty <- "LASSO"
	class(returnVal) <- "rq.pen.seq"
	returnVal
}

#coef.cutoff is actually ignored, but used here so ... works correctly, possibly bad code. 
rq.enet <- function(x,y,tau=.5,lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001), 
			penalty.factor = rep(1, ncol(x)),scalex=TRUE,tau.penalty.factor=rep(1,length(tau)),
			a=0,max.iter=10000,converge.eps=1e-7,gamma=IQR(y)/10,coef.cutoff=NULL,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	lpf <- length(penalty.factor)
	ltpf <- length(tau.penalty.factor)	
		
	if(sum(tau <= 0 | tau >=1)>0){
		stop("tau needs to be between 0 and 1")
	}
	if(sum(a < 0 | a > 1) > 0){
		stop("a needs to be >= 0 and <= 1")
	}
	if(sum(penalty.factor<0)>0){
		stop("penalty factors must be positive")
	}
	if(lpf!=p){
		stop("penalty factor must be of length p")	
	}
	if(ltpf!=nt){
		stop("tau penalty factor must be of length tau")
	}
	
	if(is.null(lambda)){
		lamMax <- getLamMax(x,y,tau,scalex=scalex)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}
	if(length(lambda)==1){
			stop("The Huber algorithm requires at least 2 values of lambda and elastic net only uses the Huber algorithm")
	}
	returnVal <- rq.lasso.huber(x,y,tau,lambda,penalty.factor,scalex,tau.penalty.factor,a=a,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,...)
	avals <- sapply(returnVal$models,modelA)
	tauvals <- sapply(returnVal$models,modelTau)
	modelsInfo <- data.frame(modelIndex=1:length(returnVal$models),a=avals,tau=tauvals)
	returnVal$modelsInfo <- modelsInfo
	returnVal$penalty <- "ENet"
	returnVal$lambda <- lambda
	returnVal$a <- a
	
	class(returnVal) <- "rq.pen.seq"
	returnVal	
}

#kind of hacky code with respect to penalty calculation 
rq.lla <- function(obj,x,y,penalty,a=ifelse(penalty=="SCAD",3.7,3),penalty.factor,tau.penalty.factor,scalex=TRUE,
					coef.cutoff=1e-8,max.iter=10000,converge.eps=1e-7,gamma=IQR(y)/10,...){
	nt <- length(obj$tau)
	na <- length(a)
	if(penalty=="SCAD"){
		derivf <- scad_deriv
		penf <- scad
	} else if(penalty=="MCP"){
		derivf <- mcp
		penf <- mcp
	} else if(penalty=="aLASSO"){
		derivf <- alasso_wt
	} else{
		stop("Penalty must be SCAD, MCP or aLASSO")
	}
	
	newModels <- vector(mode="list",length=nt*na)
	pos <- 1
	modelNames <- NULL
	ll <- length(obj$lambda)

	for(j in 1:nt){
		lampen <- penalty.factor*tau.penalty.factor[j]#as.numeric(obj$models[[j]]$penalty.factor %*% t(obj$models[[j]]$lambda))
		for(k in 1:na){	
			#pfs <- matrix(derivf(as.numeric(abs(coefficients(obj$models[[j]])[-1,])),lampen,a=a[k]),ncol=ll)
			newModels[[pos]] <- obj$models[[j]]
			penSums <- NULL
			endHit <- FALSE
			for(i in 1:ll){
				llapenf <- derivf(as.numeric(abs(coefficients(obj$models[[j]]))[-1,i]),lampen,a=a[k])
				if(sum(llapenf)==0){
					if(!endHit){
						update_est <- coefficients(rq(y~x,tau=obj$tau[j]))
						endHit <- TRUE
					}
				}
				else if(obj$alg=="huber"){
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=c(2,1),penalty.factor=llapenf,scalex=scalex,alg=obj$alg,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,...)$models[[1]])[,2]

				} else{
					update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=1,penalty.factor=llapenf,alg=obj$alg,scalex=scalex,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,...)$models[[1]])
				}
				newModels[[pos]]$coefficients[,i] <- update_est
				if(penalty=="aLASSO"){
					subPenSum <- sum(llapenf*abs(update_est[-1]))
				} else{
					subPenSum <- sum(f(update_est[-1],lampen,a[k]))
				}
				penSums <- c(penSums, subPenSum)
			}
			newModels[[pos]]$a <- a[k]
			newModels[[pos]] <- rq.pen.modelreturn(newModels[[pos]]$coefficients,x,y,newModels[[pos]]$tau,obj$lambda,local.penalty.factor=penalty.factor*tau.penalty.factor[j],penalty,a[k])
			newModels[[pos]]$PenRho <- penSums + newModels[[pos]]$rho
			modelNames <- c(modelNames,paste0("tau",obj$tau[j],"a",a[k]))
			pos <- pos+1
		}
	}
	names(newModels) <- modelNames
	obj$models <- newModels	
	#if(penalty=="aLASSO"){
	##	obj$penalty.factor <- pfs
	#}
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


rq.group.lla <- function(obj,x,y,groups,penalty=c("gAdLASSO","gSCAD","gMCP"),a=NULL,norm=2, group.pen.factor,
						tau.penalty.factor,scalex,coef.cutoff,max.iter,converge.eps,gamma,...){
	#for loop calculation of penalty factors that could maybe be removed
	nt <- length(obj$tau)
	p <- ncol(x)
	n <- nrow(x)
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
	newModels <- vector(mode="list",length=nt*na)
	pos <- 1
	modelNames <- NULL
	for(j in 1:nt){
		lampen <- group.pen.factor %*% t(obj$lambda)*tau.penalty.factor[j]
		ll <- length(obj$lambda)
		for(k in 1:na){
			newModels[[pos]] <- obj$models[[j]]	
			endHit <- FALSE
			penVals <- NULL
			for(i in 1:ll){	
				coef_by_group_deriv <- group_derivs(derivf, groups, abs(coefficients(obj$models[[j]])[-1,i]),lampen[,i],a[k],norm=norm)
				if(sum(coef_by_group_deriv)==0){
					if(n > p + 1){
						update_est <- coefficients(rq(y~x,tau=obj$tau[j]))
						endHit <- TRUE
					} else{
						i <- i -1 
						break
					}
				} else{
					if(obj$alg=="huber"){
						if(norm == 1){
							penalty.factor <- mapvalues(groups,seq(1,g),coef_by_group_deriv)
							update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=c(2,1),penalty.factor=penalty.factor,alg=obj$alg,scalex=scalex,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,...)$models[[1]])[,2]
						} else{
							update_est <- coefficients(rq.group.pen(x,y,obj$tau[j],groups,lambda=1,group.pen.factor=coef_by_group_deriv, alg=obj$alg,scalex=scalex,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,...)$models[[1]])
						}
					} else{
					#norm must be 1 as huber is only valid approach for norm 2
						penalty.factor <- mapvalues(groups,seq(1,g),coef_by_group_deriv)
						update_est <- coefficients(rq.lasso(x,y,obj$tau[j],lambda=1,penalty.factor=penalty.factor,alg=obj$alg,scalex=scalex,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,...)$models[[1]])
					}
				}
				newModels[[pos]]$coefficients[,i] <- update_est
				if(penalty == "gAdLASSO"){
					penVals <- c(penVals,sum(getGroupPen(update_est,groups,1,coef_by_group_deriv,penalty,norm,1)))
				}
				if(endHit){
					break
				}
			}
			if(i < ll){
				newModels[[pos]] <- clearModels(newModels,i)
			}
			
			newModels[[pos]] <- rq.pen.modelreturn(newModels[[pos]]$coefficients,x,y,obj$tau[j],obj$lambda[1:i],rep(1,p),penalty,a[k])	
			newModels[[pos]]$penalty.factor <- NULL	
			if(penalty=="gAdLASSO"){
				newModels[[pos]]$PenRho <- newModels[[pos]]$rho + penVals
			}
			
			dimnames(newModels[[pos]]$group.pen.factor) <- NULL		
			modelNames <- c(modelNames,paste0("tau",obj$tau[j],"a",a[k]))			
			pos <- pos + 1
		}
	}
	names(newModels) <- modelNames
	obj$models <- newModels
	obj$a <- a
	obj$modelsInfo <- createModelsInfo(obj$models)
	if(penalty != "gAdLASSO"){
	#special case of adaptive lasso done in the above for loops
	#code improve: could make this less clunky. 
		obj  <- updateGroupPenRho(obj,norm,groups,group.pen.factor,tau.penalty.factor)
	}
	obj$groups <- groups
	obj$penalty <- penalty
	#obj$class <- c(obj$class, "rq.group.pen.seq")
	obj
}

updateGroupPenRho <- function(obj,norm,groups,group.pen.factor,tau.penalty.factor){
	
	for(j in 1:length(obj$models)){
		a <- obj$models[[j]]$a
		taupos <- which(obj$tau == obj$models[[j]]$tau)
		if(length(obj$lambda)==1){
			obj$models[[j]]$PenRho <- obj$models[[j]]$rho + sum(getGroupPen(obj$models[[j]]$coefficients[-1],groups,obj$lambda,obj$models[[j]]$group.pen.factor*tau.penalty.factor[taupos],obj$penalty,norm,a)) 
		} else{
			for(i in 1:length(obj$lambda)){
				obj$models[[j]]$PenRho[i] <- obj$models[[j]]$rho[i] + sum(getGroupPen(obj$models[[j]]$coefficients[-1,i],groups,obj$lambda[i],obj$models[[j]]$group.pen.factor*tau.penalty.factor[taupos],obj$penalty,norm,a))
			}
		}
	}	
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
	rownames(modelsInfo) <- paste0("tau",tauvals,"a",avals)
	modelsInfo
}

rq.nc <- function(x, y, tau=.5,  penalty=c("SCAD","aLASSO","MCP"),a=NULL,lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001),alg="huber",scalex=TRUE,
					penalty.factor = rep(1, ncol(x)),tau.penalty.factor=rep(1,length(tau)),coef.cutoff=1e-8,max.iter=10000,converge.eps=1e-7,gamma=IQR(y)/10,...) {
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
			init.model <- rq.enet(x,y,tau,lambda=lambda,scalex=scalex,penalty.factor=penalty.factor,
						tau.penalty.factor=tau.penalty.factor,coef.cutoff=coef.cutoff,max.iter=max.iter,
						converge.eps=converge.eps,gamma=gamma,...)
		} else{
			init.model <- rq.lasso(x,y,tau,alg=alg,lambda=lambda,scalex=scalex,penalty.factor=penalty.factor,
						tau.penalty.factor=tau.penalty.factor,coef.cutoff=coef.cutoff,max.iter=max.iter,
						converge.eps=converge.eps,gamma=gamma,...)
		}
		rq.lla(init.model,x,y,penalty,a,penalty.factor,tau.penalty.factor,scalex,coef.cutoff,max.iter,converge.eps,gamma,...)
	} else{
		if(length(unique(penalty.factor))!=1 & length(unique(tau.penalty.factor))!=1){
			warning("The QICD algorithm takes predictor and tau penalty factors and turns them into a zero or 1. Zero if the weight is zero and one otherwise. Other algorithms are better if you want a more nuanced approach.")
		}
		pos <- 1
		models <- vector(mode="list",length=nt*na)
		modelNames <- NULL
		for(i in 1:nt){
			for(j in 1:na){
				coefs <- NULL
				for(lam in lambda){
					sublam <- lam*penalty.factor*tau.penalty.factor[i]
					penvars <- which(sublam != 0)
					subm <- rq.nc.fit(x,y,tau[i],lambda=lam, alg="QICD", a=a[j], coef.cutoff=coef.cutoff,converge_criteria=converge.eps,penVars=penvars,internal=TRUE,...)
					coefs <- cbind(coefs,coefficients(subm))
				}
				models[[pos]] <- rq.pen.modelreturn(coefs,x,y,tau[i],lambda,local.penalty.factor=penalty.factor*tau.penalty.factor[j],penalty,a[j])
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

rq.glasso <- function(x,y,tau,groups, lambda, group.pen.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	g <- length(unique(groups))
	nt <- length(tau)
	
	models <- vector(mode="list",length=nt)
	for(i in 1:nt){
		subtau <- tau[i]
		penf <- group.pen.factor*tau.penalty.factor[i]
		models[[i]] <- hrq_glasso(x,y,group.index=groups,tau=subtau,lambda=lambda,w.lambda=penf,scalex=scalex,gamma=gamma,max_iter=max.iter,epsilon=converge.eps,...)
		models[[i]] <- rq.pen.modelreturn(models[[i]]$beta,x,y,subtau,models[[i]]$lambda,rep(1,p),"gLASSO",1)
	}
	attributes(models)$names <- paste0("tau",tau,"a",1)
		
	returnVal <- list(models=models, n=n, p=p,alg="huber",tau=tau,penalty="gLASSO",a=1,lambda=lambda)
	returnVal <- updateGroupPenRho(returnVal,2,groups,group.pen.factor,tau.penalty.factor)
	returnVal$a <- 1
	returnVal$modelsInfo <- data.frame(modelIndex=1:nt,a=1,tau=tau)
	rownames(returnVal$modelsInfo) <- paste0("tau",tau,"a",1)
	returnVal
}

qaSIS <- function(x,y,tau=.5,linear=FALSE,...){#n.cores=1,...){
	if(linear){
		eval_function<- function(x,y,tau){ 
							q1 <- rq(y ~ x, tau)
							sum((fitted(q1)-quantile(y,tau))^2)
						}
		eval_results <- apply(x,2,eval_function,y,tau,...)
	} else{
		eval_function2 <- function(x,y,tau,...){ 
							 b <- bs(x,...)
							 q1 <- rq(y ~ b, tau)
							 sum((fitted(q1)-quantile(y,tau))^2)
						 }
		eval_results <- apply(x,2,eval_function2,y,tau,...)
	}
	#if(n.cores==1){
		
	#} else{
	#	p <- dim(x)[2]
	#	mc_func <- function(idx,...){ eval_function(x[,idx],y,...)}
	#	mc_results <- mclapply(1:p, mc_func, mc.cores=n.cores, ...)
	#	eval_results <- do.call(c,mc_results)
	#}
	order( eval_results, decreasing=TRUE)
}


#need to update this and think about what penalty.factor needs to be used. 
rq.pen.modelreturn <- function(coefs,x,y,tau,lambda,local.penalty.factor,penalty,a){
# for loop that could be removed
	penfunc <- getPenfunc(penalty)
	return_val <- NULL
	return_val$coefficients <- coefs
	#return_val$penalty.factor <- penalty.factor
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
		return_val$PenRho <- return_val$rho + sum(penfunc(return_val$coefficients[-1],lambda*local.penalty.factor,a))
		return_val$nzero <- sum(return_val$coefficients!=0)
	} else{
		return_val$rho <- apply(check(res,tau),2,mean)	
		rownames(return_val$coefficients) <- x_names
		for(i in 1:length(return_val$rho)){
			return_val$PenRho[i] <- return_val$rho[i] + sum(penfunc(return_val$coefficients[-1,i],lambda[i]*local.penalty.factor,a))
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



updateGroupPenRho <- function(obj,norm,groups,group.pen.factor,tau.penalty.factor){
	
	for(j in 1:length(obj$models)){
		a <- obj$models[[j]]$a
		taupos <- which(obj$tau == obj$models[[j]]$tau)
		if(length(obj$lambda)==1){
			obj$models[[j]]$PenRho <- obj$models[[j]]$rho + sum(getGroupPen(obj$models[[j]]$coefficients[-1],groups,obj$lambda,obj$models[[j]]$group.pen.factor*tau.penalty.factor[taupos],obj$penalty,norm,a)) 
		} else{
			for(i in 1:length(obj$lambda)){
				obj$models[[j]]$PenRho[i] <- obj$models[[j]]$rho[i] + sum(getGroupPen(obj$models[[j]]$coefficients[-1,i],groups,obj$lambda[i],obj$models[[j]]$group.pen.factor*tau.penalty.factor[taupos],obj$penalty,norm,a))
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

rq.lasso.huber.onetau <- function(x,y,tau,lambda,penalty.factor=rep(1,ncol(x)),scalex=TRUE,max.iter,converge.eps,a=1,gamma,...){
	dims <- dim(x)
	p <- dims[2]
	if(scalex){
		hqModel <- hqreg(x,y,method="quantile",tau=tau,lambda=lambda,penalty.factor=penalty.factor,max.iter=max.iter,eps=converge.eps,alpha=a,gamma=gamma,...)
	} else{
		hqModel <- hqreg_raw(x,y,method="quantile",tau=tau,lambda=lambda,penalty.factor=penalty.factor,max.iter=max.iter,eps=converge.eps,alpha=a,gamma=gamma,...)
	}
	rq.pen.modelreturn(hqModel$beta,x,y,tau,lambda,penalty.factor,"LASSO",a)
}

rq.lasso.huber <- function(x,y,tau,lambda,penalty.factor=rep(1,ncol(x)),scalex=TRUE,tau.penalty.factor=rep(1,length(tau)),max.iter,converge.eps,a=1,gamma,...){
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	nt <- length(tau)
	na <- length(a)
	
	models <- vector(mode="list",length=nt*na)
	pos <- 1
	modelnames <- NULL
	for(i in 1:nt){
		subtau <- tau[i]
		penf <- penalty.factor*tau.penalty.factor[i]
		for(j in 1:na){			
			models[[pos]] <- rq.lasso.huber.onetau(x,y,tau=subtau,lambda=lambda,penalty.factor=penf,scalex=scalex,max.iter,converge.eps,a=a[j],gamma=gamma,...)
			pos <- pos+1
			modelnames <- c(modelnames,paste0("tau",tau[i],"a",a[j]))
		}
	}
	attributes(models)$names <- modelnames
		
	returnVal <- list(models=models, n=n, p=p,alg="huber",tau=tau,a=a)
	returnVal
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

getRho <- function(model){
    model$rho
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

subtract <- function(predicted,obs){
	obs-predicted
}

modelTau <- function(object){
	object$tau
}

modelA <- function(object){
	object$a	
}

modelNz <- function(object, index){
	object$nz[index]
}

byTauResults <- function(cvErr,tauvals,avals,models,se,lambda){
#for loops!
	mn <- length(tauvals)
	overallMin <- apply(cvErr,1,min)
	overallSpot <- apply(cvErr,1,which.min)
	index <- 1:mn
	
	btr <- data.table(tau=tauvals,minCv=overallMin,lambdaIndex=overallSpot,a=avals,modelsIndex=index)
	btr <- btr[, .SD[which.min(minCv)],by=tau]
	
	cvse <- lambda1se <- lambda1seIndex <- lambdaVals <- nz <-  NULL
	for(i in 1:nrow(btr)){
		subse <- se[btr[[5]][i],btr[[3]][i]] #5 is model index and 3 is lambda index
		cvse <- c(cvse,subse)
		se1Above <- btr[[2]][1] + subse
		subLambda <- lambda[btr[[3]][i]]
		lambdaVals <- c(lambdaVals, subLambda)
		subLambda1sePos <- which(cvErr[btr[[5]][1],] <= se1Above)[1]
		lambda1seIndex <- c(lambda1seIndex,subLambda1sePos)
		subLambda1se <- lambda[subLambda1sePos]
		lambda1se <- c(lambda1se,subLambda1se)
		nz <- c(nz, models[[btr[[5]][i]]]$nz[btr[[3]][i]])
	}
	
	btr <- cbind(btr, lambda=lambdaVals, cvse = cvse, lambda1se=lambda1se, lambda1seIndex=lambda1seIndex, nonzero=nz)
	btr <- setcolorder(btr, c(1,2,6,3,8,9,4,7,5,10))
	btr
}

groupTauResults <- function(cvErr, tauvals,a,avals,models,tauWeights,lambda){
# improve note: maybe code in for loop could be improved upon by checking at each iteration if the better cv value has been found or not
	nl <- length(lambda)
	na <- length(a)
	gcve <- matrix(rep(0,na*nl),ncol=nl)
	for(i in 1:na){
		subErr <- subset(cvErr, avals==a[i])
		gcve[i,] <- tauWeights %*% subErr
	}
	rownames(gcve) <- paste0("a",a)
	minIndex <- which(gcve==min(gcve),arr.ind=TRUE)
	returnA <- a[minIndex[1]]
	modelIndex <- which(avals==returnA)
	targetModels <- models[modelIndex]
	tauvals <- sapply(targetModels,modelTau)
	lambdavals <- lambda[minIndex[1,2]]
	nz <- sapply(targetModels, modelNz, minIndex[1,2])
	minCv <- cvErr[modelIndex,minIndex[1,2]]
	list(returnTable=data.table(tau=tauvals,lambda=lambdavals,a=returnA,minCv=minCv,lambdaIndex=minIndex[1,2],modelsIndex=modelIndex, nonzero=nz),gcve=gcve)
}

re_order_nonpen_coefs <- function(nonpen_coefs, penVars, intercept=TRUE){
	p <- length(nonpen_coefs)
	new_coefs <- rep(NA,p)
	if(intercept){
		penVars <- penVars+1
		pen_output <- 2:(length(penVars)+1)
	} else{
		pen_output <- 1:length(penVars)
	}
	new_coefs[penVars] <- nonpen_coefs[pen_output]
	new_coefs[-penVars] <- nonpen_coefs[-pen_output]
	new_coefs
}

getModels <- function(x,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL){
	if( (is.null(tau)==FALSE | is.null(a)==FALSE) & is.null(modelsIndex) == FALSE){
		stop("Set tau and a or set modelsIndex, not both")
	}
	if( (is.null(lambda)==FALSE & is.null(lambdaIndex)==FALSE)){
		stop("Use lambda or lambdaIndex, not both")
	}
	lt <- length(x$tau)
	na <- length(x$a)	
	if((is.null(tau) == FALSE & is.null(a)==FALSE)){
		modelsIndex <- intersect(whichMatch(tau,x$modelsInfo$tau),whichMatch(a,x$modelsInfo$a))
	} else if(is.null(tau)==FALSE){
		modelsIndex <- whichMatch(tau,x$modelsInfo$tau)
	} else if(is.null(a) == FALSE){
		modelsIndex <- whichMatch(a,x$modelsInfo$a)
	}
	else if(is.null(modelsIndex)){
		modelsIndex <- 1:length(x$models)
	}
	if(length(modelsIndex)==0){
		stop("Invalid tau or a provided")
	}
	if(is.null(lambda)==FALSE){
		lambdaIndex <- whichMatch(lambda,x$lambda)
	} else if(is.null(lambdaIndex)){
		lambdaIndex <- 1:length(x$lambda)
	}
	if(length(lambdaIndex)==0){
		stop("Invalid lambda provided")
	}
	targetModels <- x$models[modelsIndex]
	list(targetModels=targetModels,lambdaIndex=lambdaIndex,modelsIndex=modelsIndex)
}


plotgroup.rq.pen.seq.cv <- function(x,logLambda,main,...){
	a <- x$fit$a
	na <- length(a)
	# besta <- x$gtr$a[1]
	# bestaidx <- which(a==besta)
	# a <- a[-bestaidx]
	# a <- c(besta,a)
	if(logLambda){
		lambdas <- log(x$fit$lambda)
		xtext <- expression(Log(lambda))
	} else{
		lambdas <- x$fit$lambda
		xtext <- expression(lambda)
	}
	if(is.null(main)){
		main <- "Cross validation results summarized for all tau"
	}
	maxerr <- max(x$gcve)
	plot(lambdas, x$gcve[1,],ylab="Cross Validation Error",ylim=c(0,maxerr),xlab=xtext,type="n",main=main,...)
	for(i in 1:na){
		points(lambdas,x$gcve[i,],col=i,pch=1)
	}
	bestlamidx <- x$gtr$lambdaIndex[1]
	lines(rep(lambdas[bestlamidx],2),c(-5,maxerr+5),lty=2)
	if(na > 1){
		legend("topleft",paste("a=",a),col=1:na,pch=1)
	}
}




error.bars <- function (x, upper, lower, width = 0.02, ...) 
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

plotsep.rq.pen.seq.cv <- function(x,tau,logLambda,main,...){
	if(is.null(tau)){
		tau <- x$fit$tau
	}
	keepers <- which(closeEnough(tau,x$fit$modelsInfo$tau))
	minfo <- x$fit$modelsInfo[keepers,]
	avals <- unique(minfo$a)
	na <- length(avals)
	nt <- length(tau)
	if(! length(main) %in% c(0,1,nt)){
		stop("main needs to be null or length one or the length of tau")
	}
	if(nt > 1){
		par(ask=TRUE)
	}
	if(logLambda){
		lambdas <- log(x$fit$lambda)
		xtext <- expression(Log(lambda))
	} else{
		lambdas <- x$fit$lambda
		xtext <- expression(lambda)
	}
	for(i in 1:nt){
		if(is.null(main)){
			mainText <- paste("Cross validation results for ", expression(tau)," = ", tau[i])
		} else if(length(main)==1){
			mainText <- main
		} else{
			mainText <- main[i]
		}
		subkeepers <- which(closeEnough(tau[i],minfo$tau))
		subinfo <- minfo[subkeepers,]
		suberr <- x$cverr[subkeepers,]
		
		bestkeep <- which(closeEnough(tau[i],x$btr$tau))
		subbtr <- x$btr[bestkeep,]
		besterr <- x$cverr[subbtr$modelsIndex,]
		cvsd <- x$cvse[subbtr$modelsIndex,]
		
		plot(lambdas,suberr[1,],ylim=c(0,max(max(suberr),max(besterr+cvsd))),ylab="Cross Validation Error", xlab=xtext,main=mainText,type="n",...)
		for(j in 1:na){
			points(lambdas,suberr[j,],col=j)
		}
		
		#then get index and plot error bars for this. 
		#get the best and create the segments for that
		segments(lambdas,besterr-cvsd,lambdas,besterr+cvsd)
		segments(lambdas-.01,besterr-cvsd,lambdas+.01,besterr-cvsd)
		segments(lambdas-.01,besterr+cvsd,lambdas+.01,besterr+cvsd)
		if(na>1){
			legend("topleft",paste("a=",subinfo$a),col=1:na,pch=1)
		}
		lidx <- subbtr$lambdaIndex
		lidxse <- subbtr$lambda1seIndex
		
		lines(rep(lambdas[lidx],2),c(-5,max(besterr+cvsd)+1),lty=2)
		lines(rep(lambdas[lidxse],2),c(-5,max(besterr+cvsd)+1),lty=2)
	}
	if(nt > 1){
		par(ask=FALSE)
	}
}
