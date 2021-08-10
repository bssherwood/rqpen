qic <- function(model,n, method="BIC"){
	tau <- model$tau
	df <- sum(model$coefficients != 0)
	if(method=="PBIC"){
		log(model$rho*n) + df*log(n)*log(length(model$coefficients))/(2*n)
	} else if (method == "BIC"){		    
		log(model$rho*n) + df*log(n)/(2*n)
	} else if (method== "AIC"){
		log(model$rho*n) + df/n
	} else{
		stop("invalid method")
	}
}



qic.select <- function(obj, method="BIC",septau=FALSE,weights=NULL){
# code help: Maybe think about how the qic values are returned for the septau=TRUE case
	if(class(obj) == "rq.pen.seq.cv"){
		obj <- obj$fit
	} else if(class(obj) != "rq.pen.seq"){
		stop("obj must be of class rq.pen.seq or rq.pen.seq.cv")
	}
	if(is.null(weights)==FALSE & septau){
		warning("Weights are only used when septau is set to true.")
	}
	if(is.null(weights) & !septau){
		weights <- rep(1,length(obj$tau))
	}
	
	n <- obj$n
	if(septau & length(weights)==1){
		warning("septau set to TRUE, but only one quantile modeled")
	}
	nt <- length(obj$tau)
	na <- length(obj$a)
	nl <- length(obj$lambda)
	
	qic_vals <- sapply(obj$models,qic,n,method)
	if(septau){
		minQIC <- apply(qic_vals,2,min)
		lambdaIndex <- apply(qic_vals,2,which.min)
		modelsInfo <- cbind(obj$modelsInfo,minQIC,lambdaIndex)
		modelsInfo <- data.table(modelsInfo)
		modelsInfo <- modelsInfo[, .SD[which.min(minQIC)],by=tau]
		
		coefs <- vector(mode="list", length=nt)
		for(i in 1:nt){
			coefs[[i]] <- coef(obj$models[[modelsInfo$modelIndex[i]]])[,modelsInfo$lambdaIndex[i]]
		}
		gic <- NULL
	} else{
		gic <- matrix(rep(0,na*nl),ncol=nl)
		tqic_vals <- t(qic_vals)
		for(i in 1:na){
			subIC <- subset(tqic_vals, obj$modelsInfo$a==obj$a[i])
			gic[i,] <- weights %*% subIC
		}
		minIndex <- which(gic==min(gic),arr.ind=TRUE)
		returnA <- obj$a[minIndex[1]]
		modelsInfo <- subset(obj$modelsInfo, a==returnA)
		modelsInfo <- cbind(modelsInfo,minQIC=tqic_vals[modelsInfo$modelIndex,minIndex[2]],lambdaIndex=minIndex[2])
		coefs <- coef(obj, lambdaIndex=minIndex[2], modelsIndex=modelsInfo$modelIndex)
	}
	coefIndex <- 1:nt
	modelsInfo <- cbind(modelsInfo, coefIndex)
	
	return_val <- list(coefficients = coefs, ic=qic_vals,modelsInfo=modelsInfo, gic=gic)
	class(return_val) <- "qic.select"
	return_val
}


print.qic.select <- function(x,...){
   print(coefficients(x))
}

predict.qic.select <- function(x, newdata, ...){
	coefs <- do.call(cbind(coefficients(x)))
	cbind(1,newdata) %*% coefs
}


print.rq.pen.seq <- function(x,...){
  nt <- length(x$tau)
  na <- length(x$a)
  if(nt==1 & na==1){
    print(data.frame(nzero=x$models[[1]]$nzero,lambda=x$lambda))
  } else if(nt > 1 & na > 1){
    print(paste(c(paste(c("Quantile regression with ", x$penalty, " penalty for quantiles:",x$tau), collapse=" ")," and tuning parameters a:", x$a),collapse=" "))
  } else if( na > 1){
    print(paste(c(paste(c("Quantile regression with ", x$penalty, " penalty for quantile:",x$tau), collapse=" ")," and tuning parameters a:", x$a),collapse=" "))
  } else{
    print(paste(c("Quantile regression with ", x$penalty, " penalty for quantiles:",x$tau), collapse=" "))
  }	
}

coef.cv.rq.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  coefficients(object$models[[target_model]])
}


#' Title
#'
#' @param x matrix of predictors
#' @param y vector of responses
#' @param tau vector of quantiles
#' @param lambda vector of lambda, if not set will be generated automatically
#' @param penalty choice of penalty
#' @param a additional tuning parameter, not used for lasso or ridge penalties
#' @param nlambda number of lambda, ignored if lambda is set
#' @param eps If not pre-specified the lambda vector will be from lambda_max to lambda_max times eps
#' @param penalty.factor penalty factor for the predictors
#' @param alg algorithm used
#' @param scalex Whether x should be scaled before fitting the model. Coefficients are returned on the orginal scale. 
#' @param tau.penalty.factor A penalty factor for each quantile.
#' @param coef.cuttoff Some of the linear programs will provide very small, but not sparse solutions. Estimates below this number will be set to zero. This is ignored if a non-linear programming algorithm is used. 
#' @param max.iter Maximum number of iterations of non-linear programming algorithms.
#' @param converge.eps Convergence threshold for non-linear programming algorithms. 
#' @param gamma tuning parameter for Huber loss, not applicable for non-huber algorithms. 
#'
#' @return
#' @export
#'
#' @examples
rq.pen <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001), 
	penalty.factor = rep(1, ncol(x)),alg=ifelse(sum(dim(x))<200,"huber","br"),scalex=TRUE,tau.penalty.factor=rep(1,length(tau)),
	coef.cutoff=1e-8,max.iter=10000,converge.eps=1e-7,gamma=IQR(y)/10,lambda.discard=TRUE,...){
	penalty <- match.arg(penalty)
	if(min(penalty.factor) < 0 | min(tau.penalty.factor) < 0){
		stop("Penalty factors must be non-negative.")
	}
	if(sum(penalty.factor)==0 | sum(tau.penalty.factor)==0){
		stop("Cannot have zero for all entries of penalty factors. This would be an unpenalized model")
	}
	if(scalex){
		x <- scale(x)
	}
	if(penalty=="LASSO"){
		fit <- rq.lasso(x,y,tau,lambda,nlambda,eps,penalty.factor,alg,scalex=FALSE,tau.penalty.factor,coef.cutoff,max.iter,converge.eps,gamma,lambda.discard=lambda.discard,...)
	} else if(penalty=="Ridge"){
		if(alg != "huber"){
			stop("huber alg is only option for Ridge penalty")
		}
		fit <-  rq.enet(x,y,tau,lambda,nlambda,eps, penalty.factor,scalex=FALSE,tau.penalty.factor,a=0,max.iter,converge.eps,gamma,lambda.discard=lambda.discard,...)
	} else if(penalty == "ENet"){
		if(alg != "huber"){
			stop("huber alg is only option for ENet penalty")
		}
		if(is.null(a)){
			stop("Specify a value for a for Enet penalty")
		}
		fit <- rq.enet(x,y,tau,lambda,nlambda,eps, penalty.factor,scalex=FALSE,tau.penalty.factor,a,max.iter,converge.eps,gamma,lambda.discard=lambda.discard,...)
	} else if(penalty == "aLASSO" | penalty=="SCAD" | penalty == "MCP"){
		fit <- rq.nc(x,y,tau,penalty,a,lambda,nlambda=nlambda,eps=eps,penalty.factor=penalty.factor,alg=alg,scalex=FALSE,tau.penalty.factor=tau.penalty.factor,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma,lambda.discard=lambda.discard,...)
	}
	if(scalex){
		for(i in 1:length(fit$models)){
			if(!is.null(dim(fit$models[[i]]$coefficients))){
				fit$models[[i]]$coefficients <- apply(fit$models[[i]]$coefficients,2,transform_coefs,attributes(x)$`scaled:center`,attributes(x)$`scaled:scale`)
			} else{
				fit$models[[i]]$coefficients <- transform_coefs(fit$models[[i]]$coefficients,attributes(x)$`scaled:center`,attributes(x)$`scaled:scale`)
			}
		}
	}
	fit$call <- match.call()
	fit
}

updateCoefs <- function(obj,mu_x,sigma_x){
	transform_coefs(return_val$coefficients,mu_x,sigma_x,intercept)
}
predict.models <- function(object, newx){
	cbind(1,newx) %*% object$coefficients
}

predict.errors <- function(object, newx, newy){
	preds <- predict(object,newx)
	errors <- lapply(preds,subtract,newy)
}

check.errors <- function(object,newx,newy){
	lapply(object$models,check.errors.model,newx,newy)
}

check.errors.model <- function(object,newx,newy){
	preds <- predict.models(object,newx)
	errors <- newy- preds
	check(errors,object$tau)
}

predict.rq.pen.seq <- function(object, newx){
	lapply(object$models, predict.models,newx=newx)
}

rq.pen.cv <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,nlambda=100,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),printProgress=FALSE,...){
	n <- length(y)
	if(is.null(weights)==FALSE){
		stop("weights not currently implemented. Can use older function cv.rq.pen, but it supports fewer penalties and is slower.")
	}
	if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
    }
	if(is.null(a)==FALSE & ( penalty[1] == "LASSO" | penalty[1] == "Ridge")){
		warning("For Ridge and LASSO the tuning parameter a is ignored")
	}
	fit <- rq.pen(x,y,tau,lambda=lambda,penalty=penalty,a=a,...)
	nt <- length(tau)
	na <- length(fit$a)
	nl <- length(fit$lambda)
	if(!groupError){
		indErrors <- matrix(rep(0,nt*na*nl),nrow=nl)
	}
	foldErrors <- fe2ndMoment <- matrix(rep(0,nt*na*nl),ncol=nl)
    for(i in 1:nfolds){
		if(printProgress){
			print(paste("Working on fold",i))
		}
		train_x <- x[foldid!=i,]
		train_y <- y[foldid!=i]
		test_x <- x[foldid==i,,drop=FALSE]
		test_y <- y[foldid==i]
		trainModel <- rq.pen(train_x,train_y,tau,lambda=fit$lambda,penalty=penalty,a=fit$a,...)
		if(is.null(cvFunc)){
			testErrors <- check.errors(trainModel,train_x,train_y)
		} else{
			testErrors <- lapply(predict.errors(trainModel,test_x,test_y),cvFunc)
		}
		if(!groupError){
			indErrors <- indErrors + sapply(testErrors,apply,2,sum)
		}
		#code improve maybe should replace below with sapply, but then we get the transpose which effects downstream code
		foldMeans <- do.call(rbind, lapply(testErrors,apply,2,cvSummary))
		foldErrors <- foldErrors + foldMeans
		fe2ndMoment <- fe2ndMoment + foldMeans^2
    }
	fe2ndMoment <- fe2ndMoment/nfolds
	foldErrors <- foldErrors/nfolds
	stdErr <- sqrt( (nfolds/(nfolds-1))*(fe2ndMoment - foldErrors^2))
	tauvals <- sapply(fit$models,modelTau)
	avals <- sapply(fit$models,modelA)
	if(groupError){
		btr <- byTauResults(foldErrors,tauvals,avals,fit$models,stdErr,fit$lambda)
		gtr <- groupTauResults(foldErrors, tauvals,fit$a,avals,fit$models,tauWeights,fit$lambda)
	} else{
		indErrors <- t(indErrors)/n
		btr <- byTauResults(indErrors,tauvals,avals,fit$models,stdErr,fit$lambda)
		gtr <- groupTauResults(indErrors, tauvals,fit$a,avals,fit$models,tauWeights,fit$lambda)
	}

	returnVal <- list(cverr = foldErrors, cvse = stdErr, fit = fit, btr=btr, gtr=gtr$returnTable, gcve=gtr$gcve, call=match.call())
	class(returnVal) <- "rq.pen.seq.cv"
	returnVal
}

print.rq.pen.seq.cv <- function(x,...){
	if(length(x$fit$tau)==1){
		cat("\nCross validation tuning parameter choices\n")
		print(x$btr)
	} else{
		cat("\nCross validation tuning parameter optimized for each quantile\n")
		print(x$btr)
		cat("\nCross validation tuning parameter optimized across all quantiles\n")
		print(x$gtr)
	}
}

coef.rq.pen.seq.cv <- function(x,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,tau=NULL,...){
	if(!useDefaults){
		coefficients(x$models,tau=tau,...)
	} else{
		if(is.null(tau)){
			tau <- m1$fit$tau
		}
		if(septau){
			keepers <- which(closeEnough(tau,x$btr$tau))
			btr <- x$btr[keepers,]
			models <- x$fit$models[btr$modelsIndex]
			if(cvmin){
				lambdaIndex <- btr$lambdaIndex
			} else{
				lambdaIndex <- btr$lambda1seIndex
			}
			returnVal <- vector(mode="list", length=length(models))
			names(returnVal) <- names(models)
			for(i in 1:length(returnVal)){
				returnVal[[i]] <- coef(x$fit,modelsIndex=btr$modelsIndex[i],lambdaIndex=lambdaIndex[i])[[1]]
			}
			returnVal
		} else{
			if(!cvmin){
				stop("One standard error approach not implemented for group choice of tuning parameter")
			} else{
				keepers <- which(closeEnough(tau,x$gtr$tau))
				gtr <- x$gtr[keepers,]#subset(x$gtr, closeEnough(tau,x$gtr$tau))
				coef(x$fit,modelsIndex=gtr$modelsIndex,lambdaIndex=gtr$lambdaIndex[1])
			}
		}
	}
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

rq.group.pen.cv <- function(x,y,tau=.5,groups=1:ncol(x),lambda=NULL,a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),printProgress=FALSE,...){
	n <- length(y)
	if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
    }
	fit <- rq.group.pen(x,y,tau,lambda=lambda,groups=groups,a=a,...)
	nt <- length(tau)
	na <- length(fit$a)
	nl <- length(fit$lambda)
	if(!groupError){
		indErrors <- matrix(rep(0,nt*na*nl),nrow=nl)
	}
	foldErrors <- fe2ndMoment <- matrix(rep(0,nt*na*nl),ncol=nl)
    for(i in 1:nfolds){
		if(printProgress){
			print(paste("Working on fold",i))
		}
		train_x <- x[foldid!=i,]
		train_y <- y[foldid!=i]
		test_x <- x[foldid==i,,drop=FALSE]
		test_y <- y[foldid==i]
		trainModel <- rq.group.pen(train_x,train_y,tau,groups=groups,lambda=fit$lambda,a=fit$a,...)
		if(is.null(cvFunc)){
			testErrors <- check.errors(trainModel,train_x,train_y)
		} else{
			testErrors <- lapply(predict.errors(trainModel,test_x,test_y),cvFunc)
		}
		if(!groupError){
			indErrors <- indErrors + sapply(testErrors,apply,2,sum)
		}
		#code improve maybe should replace below with sapply, but then we get the transpose which effects downstream code
		foldMeans <- do.call(rbind, lapply(testErrors,apply,2,cvSummary))
		foldErrors <- foldErrors + foldMeans
		fe2ndMoment <- fe2ndMoment + foldMeans^2
    }
	fe2ndMoment <- fe2ndMoment/nfolds
	foldErrors <- foldErrors/nfolds
	stdErr <- sqrt( (nfolds/(nfolds-1))*(fe2ndMoment - foldErrors^2))
	tauvals <- sapply(fit$models,modelTau)
	avals <- sapply(fit$models,modelA)
	if(groupError){
		btr <- byTauResults(foldErrors,tauvals,avals,fit$models,stdErr)
		gtr <- groupTauResults(foldErrors, tauvals,fit$a,avals,fit$models,tauWeights)
	} else{
		indErrors <- t(indErrors)/n
		btr <- byTauResults(indErrors,tauvals,avals,fit$models,stdErr)
		gtr <- groupTauResults(indErrors, tauvals,fit$a,avals,fit$models,tauWeights)
	}

	returnVal <- list(cverr = foldErrors, cvse = stdErr, fit = fit, btr=btr, gtr=gtr$returnTable, gcve=gtr$gcve)
	class(returnVal) <- "rq.pen.seq.cv"
	returnVal
}


cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",criteria = "CV",intercept=TRUE,cvFunc="check",nfolds=10,foldid=NULL,nlambda=100,eps=.0001,init.lambda=1,penVars=NULL,alg=ifelse(ncol(x)<50,"LP","QICD"),internal=FALSE,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
# nfolds: number of folds for cross validation
# foldid: preset id of folds
# penVar: variables to be penalized, default is all non-intercept terms

  # Pre-algorithm setup/ get convenient values
  if(!internal){
	warning("Recommend using rq.pen.cv instead. This is an older function that is kept for reproducibality reasons, but tends to be much slower. It will not be exported to the namespace in future versions")
  }
  m.c <- match.call() # This stores all the arguments in the function call as a list

  p <- dim(x)[2]
  if(is.null(penVars)){
    penVars <- 1:p
  }
  p_range <- penVars + intercept
  n <- dim(x)[1]
  pen_func <- switch(which(c("LASSO","SCAD","MCP")==penalty), lasso, scad, mcp)
 
  ### QICD ###
  if( alg=="QICD" & penalty!="LASSO" ){
    if(criteria=="CV"){
      stop("CV criteria not implemented for QICD algorithm with nonconvex penalties. Please use BIC or PBIC instead")
    }
    m.c[["alg"]] <- "LP" #maybe this should be moved inside the is.null initial beta if statement. I don't think it matters, but might be cleaner code
    penname <- penalty

    if( !all(penVars==1:p) ){ # Some unpenalized coefficients
      z    <- as.matrix(x[,-penVars])
      xpen <- as.matrix(x[,penVars])
      QICD_func <- "QICD.nonpen"
      mapback <- order( c(penVars, (1:p)[-penVars]) ) # reorders the coefficients properly if some (non-intercept) coefficients are not penalized 
      if( intercept )
        mapback <- c(1, 1+mapback)
    } else { # All penalized coefficients
      z <- NULL
      xpen <- x
      QICD_func <- "QICD"
      mapback <- 1:p # no reordering necessary if all (non-intercept) coefficients are penalized
      if( intercept )
        mapback <- c(1, 1+mapback)
    }

    # The QICD algorithm needs good starting values, use LASSO solution 
    ## Speed things up using BIC, not k-fold, to select lambda for LASSO
    ## Can skip this part if starting values are provided
    if( is.null(m.c[["initial_beta"]]) ){
      m.c[["penalty"]] <- "LASSO"
      m.c[["criteria"]] <- "BIC"
	  if(is.null(m.c[["lambda"]])==FALSE){
		m.c[["lambda"]] <- NULL
	  }
      suppressWarnings(
        m.c[["initial_beta"]] <- coefficients( eval.parent(m.c) )
        # QICD.start <- coefficients( cv.rq.pen(x,y,tau=tau,lambda=lambda,penalty="LASSO",intercept=intercept,criteria="BIC",nlambda=nlambda,eps=eps,init.lambda=lambda,penVars=penVars,...) ) # Use the LASSO with BIC
      )
    }

    # Start in middle of lambda vector
    ## Increase lambda until intercept only model (or all penlized coefficients are zero)
    ## Decrease lambda until full model (Sparsity assumption => full model is bad)
    m.c[[1]] <- as.name(QICD_func)
    m.c[["penalty"]] <- penalty
    m.c[["x"]] <- xpen
    m.c$z <- z

    if( is.null(lambda) ){
      lambdas <- exp( seq(-7, 1, length=100) ) 
    } else {
      lambdas <- lambda
    }
    
    coefs <- matrix(NA, p+intercept, length(lambdas))

    ### Loop for increasing lambda
    for( i in floor(length(lambdas)/2):length(lambdas) ){
      m.c[["lambda"]] <- lambdas[i]
      coefs[,i] <- eval.parent( m.c )[mapback]
      # coefs[,i] <- QICD_func( y=y,x=x, z=z, tau=tau, lambda=lambdas[i], intercept=intercept, 
      #                       penalty=penalty, initial_beta=QICD.start, ... )[mapback]
      if( all(coefs[p_range,i] == 0) )
        break()
    }

    ### Loop for decreasing lambda
    for( i in floor(length(lambdas)/2):2 -1 ){
      m.c[["lambda"]] <- lambdas[i]
      coefs[,i] <- eval.parent( m.c )[mapback]
      # coefs[,i] <- QICD_func( y=y,x=x, z=z, tau=tau, lambda=lambdas[i], intercept=intercept, 
      #                       penalty=penalty, initial_beta=QICD.start, ... )[mapback]
      if( all(coefs[p_range,i] != 0) )
        break()
    }

    #### Remove the NA columns from coefs and corresponding lambdas
    lambdas.keep <- which( !is.na(coefs[1,]) )
    lambdas <- lambdas[lambdas.keep]
    coefs <- coefs[,lambdas.keep]
    rownames(coefs) <- names( m.c[["initial_beta"]] )
    XB <- x%*%coefs[p_range,]
    if( intercept )
      XB <- XB + matrix(coefs[1,], n, ncol(coefs), byrow=TRUE)

    residuals <- y - XB
    rho <- colSums( check(residuals, tau=tau) )
    if( is.null(m.c[["a"]]) )
      a <- 3.7
    PenRho <- rho + colSums(apply( rbind(lambdas, coefs), 2, 
                    function(xx) pen_func(xx[1+p_range], lambda=xx[1], a=a) ))

    cv <- data.frame(lambda=lambdas, cve=NA)
    
	if( criteria=="AIC" ){
	  cv$cve <- log(rho) + colSums(coefs!=0)/n
      names(cv)[2] <- "BIC"
	}
    if( criteria=="BIC" ){
      cv$cve <- log(rho) + colSums(coefs!=0)*log(n)/(2*n)
      names(cv)[2] <- "BIC"
    } else { # PBIC
      cv$cve <- log(rho) + colSums(coefs!=0)*log(n)*log(nrow(coefs))/(2*n)
      names(cv)[2] <- "PBIC"
    }

    lambda.min <- lambdas[which.min(cv[,2])]

    # Final cleanup for QICD
    ## First get models for each lambda
    models <- vector( "list", length(lambdas) )
    for( j in 1:length(models) ){
      models[[j]]$coefficients <- coefs[,j]
      models[[j]]$PenRho <- PenRho[j]
      models[[j]]$residuals <- residuals[,j]
      models[[j]]$rho <- rho[j]
      models[[j]]$tau <- tau
      models[[j]]$n <- n
      models[[j]]$intercept <- intercept
      models[[j]]$penalty <- penalty
      class(models[[j]]) <- c("rq.pen", "rqNC")
    }

    return_val <- list( models=models, cv=cv, lambda.min=lambda.min, penalty=penalty )
    class(return_val) <- "cv.rq.pen"
  } else{
  ############



  # If no lambdas provided, find reasonable lambdas to use
  if(is.null(lambda)){
  # find a lambda that sets all coefficients to zero. 
  # Strategy is to get \lambda \sum p_\lambda(|\beta_j}) >= \sum \rho_\tau(y-quantile(y,tau)
  # do this by fitting model with lambda = init.lambda and then set new lambda such that 
  # \lambda* = \sum \rho_\tau(y-quantile(y,tau)) / \sum p_\lambda(|beta_j|) repeat as needed
     sample_q <- quantile(y,tau)
     inter_only_rho <- sum(check(y-sample_q,tau))
     #lambda_star <- rep(0,p)
     #lambda_star[penVars] <- init.lambda
	 lambda_star <- init.lambda
     searching <- TRUE
     while(searching){
       if(penalty=="LASSO"){
         init_fit <- rq.lasso.fit(x,y,tau,lambda=lambda_star,weights,intercept,penVars=penVars,...)
       } else{
         init_fit <- rq.nc.fit(x,y,tau,lambda=lambda_star,weights,intercept,penVars=penVars,internal=TRUE,...)
       }
       if(sum(init_fit$coefficients[p_range])==0){
         searching <- FALSE     
       } else{
         lambda_star <- inter_only_rho / sum(sapply(init_fit$coefficients[p_range],pen_func,1)) 
         #1 used here because solving for lambda
       }
     }
     lambda_min <- eps*lambda_star
     lambda <- exp(seq(log(max(lambda_min)),log(max(lambda_star)),length.out=nlambda))#max is included to handle cases where
     # some variables are unpenalized and thus lambda is a multivalued vector with some zeros
  }
  # lambda is the vector of reasonable choices of lambda to use in the penalty

  models <- list()
  fit_models <- TRUE
  lam_pos <- 1
  if(penalty=="LASSO"){
   while(fit_models){
		if(fit_models){
			models[[lam_pos]] <- rq.lasso.fit(x,y,tau,lambda[lam_pos],weights,intercept,penVars=penVars,...)
      
		}
		if(sum(abs(coefficients(models[[lam_pos]])[p_range]))==0 || lam_pos==length(lambda)){
		#if we got a fully sparse model, no need to fit more sparse models
			fit_models <- FALSE
			lambda <- lambda[1:lam_pos]
		}
    lam_pos <- lam_pos + 1
	 }
  } else{
  	while(fit_models){
  		if(fit_models){
  			models[[lam_pos]] <- rq.nc.fit(x,y,tau,lambda[lam_pos],weights,intercept,penalty=penalty,penVars=penVars,...)
      }
  		if(sum(abs(coefficients(models[[lam_pos]])[p_range]))==0 || lam_pos==length(lambda)){
  			fit_models <- FALSE
  			lambda <- lambda[1:lam_pos]
  		}
      lam_pos <- lam_pos + 1
  	 }
     #models <- lapply(lambda,rq.nc.fit, x=x,y=y,tau=tau,weights=weights,intercept=intercept,penalty=penalty,
     #                                   penVars=penVars,...)
  }
  cv_results <- NULL
  if(criteria=="CV"){
    if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
    }
    for(i in 1:nfolds){
      train_x <- x[foldid!=i,]
      train_y <- y[foldid!=i]
      test_x <- x[foldid==i,,drop=FALSE]
      test_y <- y[foldid==i]
	  train_weights <- weights[foldid!=i] #not sure this line is needed
	  if(is.null(weights)){
		train_weights <- test_weights <- NULL
	  } else{
	    train_weights <- weights[foldid!=i]
		test_weights <- weights[foldid==i]
	  }
      if(penalty=="LASSO"){
         cv_models <- lapply(lambda,rq.lasso.fit, x=train_x,y=train_y,tau=tau,weights=train_weights,intercept=intercept,penVars=penVars,...)
      } else{
         cv_models <- lapply(lambda,rq.nc.fit, x=train_x,y=train_y,tau=tau,weights=train_weights,intercept=intercept,penalty=penalty,penVars=penVars,...)
      }
      if(cvFunc=="check"){
         cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, test_y, test_weights, tau=tau))
      } else{
         cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, test_y, test_weights, func=cvFunc))
      } 
    }
    cv_results <- apply(cv_results,1,mean)
  }
  if(criteria=="BIC"){
    cv_results <- sapply(models,qbic)
  }
  if(criteria=="PBIC"){
    cv_results <- sapply(models,qbic,largeP=TRUE)
  }
  lambda.min <- lambda[which.min(cv_results)]
  return_val <- NULL
  return_val$models <- models
  return_val$cv <- data.frame(lambda=lambda, cve=cv_results)
  colnames(return_val$cv)[2] <- criteria
  return_val$lambda.min <- lambda.min
  return_val$penalty <- penalty
  class(return_val) <- "cv.rq.pen"
  }

  return_val
}

rq.nc.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                      penalty="SCAD",a=3.7,iterations=1,converge_criteria=1e-06,
                      alg=ifelse(p<50,"LP","QICD"),penVars=NULL,internal=FALSE,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# lambda takes values of 1 or p
# penalty SCAD or MCP
# penVars - variables to be penalized, doesn't work if lambda has multiple entries
  if(!internal){
	warning("Recommend using rq.pen() instead. This is an older function and usually will run slower. It will not be exported to the namespace in future versions.")
  }
  p <- ncol(x)
  n <- nrow(x)

  if( alg=="QICD" ){
  ### QICD Algorithm ###
	coefnames <- paste("x",1:p, sep="") ### Coefficient names
    if( length(lambda) != 1 )
      stop( "QICD Algorithm only allows 1 lambda value")	
    ### Check if we are using QICD or QICD.nonpen
    if( is.null(penVars) | length(penVars) == p){ ### No unpenalized coefficients
      
      coefs <- QICD(y, x, tau, lambda, intercept, penalty, eps=converge_criteria, a=a, ...)
      penbeta <- intercept + 1:p ### Use later to calculate objective function
    } else { ### Some unpenalized coefficients
      z    <- as.matrix(x[,-penVars])
      xpen <- as.matrix(x[,penVars])
      #coefnames <- paste("x",1:ncol(xpen), sep="") ### Coefficient names
      #coefnames <- c( coefnames, paste("z",1:ncol(z), sep="") )
      penbeta <- intercept + penVars ### Use later to calculate objective function
      coefs <- QICD.nonpen(y, xpen, z, tau, lambda, intercept, penalty, eps=converge_criteria, a=a, ...)
	  coefs <- re_order_nonpen_coefs(coefs, penVars, intercept)
    }

    ### Add extra information to QICD output
    
    if( intercept ){ ### Residuals
      residuals <- c( y - x%*%(coefs[-1]) - coefs[1] )
	  names(coefs) <- c("intercept",coefnames)
    } else {
      residuals <- c( y - x%*%coefs )
	  names(coefs) <- coefnames
    }
    rho <- sum( check(residuals) )
	#1/n*sum( check(residuals) ) ### rho
    if( penalty == "LASSO" ){ ### PenRho for LASSO
      PenRho <- sum( abs( coefs[penbeta] )*lambda )
    } else if( penalty == "SCAD" ){ ### PenRho for SCAD
      PenRho <- sum( scad( coefs[penbeta], lambda, a ))
    } else { ### PenRho for MCP
      PenRho <- sum( mcp( coefs[penbeta], lambda, a ))
    }
    PenRho <- rho + PenRho

    sub_fit <- list( coefficients=coefs,  PenRho=PenRho, residuals=residuals,
                     rho=rho, tau=tau, n=n , intercept=intercept)
  ######################
  } else {  
  ###  LP Algorithm  ###
	   if(penalty=="SCAD"){
		 deriv_func <- scad_deriv
	   }
	   if(penalty=="MCP"){
		 deriv_func <- mcp_deriv
	   }
	   if(is.null(dim(x))){                                                                                    
		  stop('x needs to be a matrix with more than 1 column')
	   }
	   if(n != length(y)){
		  stop('length of y and rows of x do not match')
	   }
	   if(is.null(lambda)==TRUE | (length(lambda) != 1 & length(lambda) != dim(x)[2])){
		  stop(paste('input of lambda must be of length 1 or', dim(x)[2]))
	   }
	   if( sum(lambda < 0) > 0){
		  stop('lambda must be positive')
	   }
	   if(is.null(penVars) !=TRUE & length(lambda) == 1){
		  mult_lambda <- rep(0,p)
		  mult_lambda[penVars] <- lambda
		  lambda <- mult_lambda
	   }
	   if(length(lambda) != 1){
		  penVars <- (1:p)[lambda != 0]
		  pen_range <- intercept + penVars
	   } else{
		  pen_range <- intercept + 1:p
		  if(is.null(penVars)){
			 penVars <- 1:p
		  }
	   }
	   #lambda_update <- n*lambda
	   lambda_update <- lambda
	   iter_complete <- FALSE
	   iter_num <- 0
	   old_beta <- rep(0, p+intercept)
	   while(!iter_complete){
		  sub_fit <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda_update,weights=weights,intercept=intercept,...)
		  if(length(lambda)==1){
			  lambda_update <- sapply(abs(sub_fit$coefficients[pen_range]),deriv_func, lambda=lambda, a=a)
		  } else{
			  lambda_update[-penVars] <- 0
			  lambda_update[penVars] <- mapply(deriv_func, x=abs(sub_fit$coefficients[pen_range]),
															 lambda=lambda[penVars],
															 MoreArgs=list(a=a))
		  }
		  #lambda_update <- n*lambda_update
		  iter_num <- iter_num + 1
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
    ######################
    }
	 sub_fit$penalty <- penalty
	 class(sub_fit) <-  c("rq.pen", "rqNC")
	 return_val <- sub_fit
	return(return_val)
}

coef.rq.pen.seq <- function(x,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL){
	models <- getModels(x,tau,a,lambda,modelsIndex,lambdaIndex)
	lapply(models$targetModels,getModelCoefs,models$lambdaIndex)
}

plot.rq.pen.seq.cv <- function(x,septau=TRUE,tau=NULL,logLambda=FALSE,main=NULL,...){
	if(septau){
		plotsep.rq.pen.seq.cv(x,tau,logLambda,main,...)
	} else{
		if(is.null(tau)==FALSE){
			stop("Tau cannot be set if septau set to FALSE")
		}
		plotgroup.rq.pen.seq.cv(x,logLambda,main,...)
	}
}

plot.rq.pen.seq <- function(x,vars=NULL,logLambda=FALSE,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL,main=NULL, ...){
	models <- getModels(x,tau=tau,a=a,lambda=lambda,modelsIndex=modelsIndex,lambdaIndex=lambdaIndex)
	tm <- models$targetModels
	li <- models$lambdaIndex
	ml <- length(tm)
	if(!is.null(main) & ( length(main) != ml | length(main) !=1)){
		stop(paste("Main needs to be of length one or length ", ml, ", the number of models being plotted"))
	}
	if(ml > 1){
		par(ask=TRUE)
	}
	for(i in 1:ml){
		if(is.null(main)){
			mainText <- paste("Plot for tau = ", tm[[i]]$tau, " and a = ", tm[[i]]$a)
		} else if(length(main)==1){
			mainText <- main
		} else{
			main <- main[i]
		}
		if(logLambda){
			lambdas <- log(x$lambda)
			xtext <- expression(Log(lambda))
		} else{
			lambdas <- x$lambda
			xtext <- expression(lambda)
		}
		betas <- tm[[i]]$coefficients[-1,li]
		plot(lambdas, betas[1,], type="n",ylim=c(min(betas),max(betas)),ylab="Coefficient Value",xlab=xtext,main=mainText,...)
		for(i in 1:dim(betas)[1]){
			lines(lambdas, betas[i,],col=i)
		}
	}
	if(ml > 1){
		par(ask=FALSE)	
	}
}

beta_plots <- function(model,voi=NULL,logLambda=TRUE,loi=NULL,...){
#voi - index variables of interest
#logLambda - lambdas plotted on log scale
#loi - index of target lambdas
  if( "cv.rq.group.pen" %in% class(model)){
	betas <- t( model$beta)
	if(model$intercept){
		betas <- betas[,-1]
	}
  }
  else{
	betas <- t(sapply(model$models, coefficients))
	if(is.null(voi)==FALSE){
		betas <- betas[,voi]
	}
	if(colnames(betas)[1]=="intercept"){
		betas <- betas[,-1]
	}
  }
  if(logLambda){
	lambdas <- log(model$cv$lambda)
	xlabText <- "Log Lambda"
  } else{
	lambdas <- model$cv$lambda
	xlabText <- "Lambda"
  }
  
  if(is.null(loi)==FALSE){
     lambdas <- lambdas[loi]
  }
  plot(lambdas, betas[,1], type="n",ylim=c(min(betas),max(betas)),ylab="Coefficient Value",xlab=xlabText,...)
  for(i in 1:dim(betas)[2]){
    lines(lambdas, betas[,i],col=i)
  }  
}

bytau.plot <- function(x,...){
	UseMethod("bytau.plot")
} 

bytau.plot.rq.pen.seq <- function(x,a=NULL,lambda=NULL,lambdaIndex=NULL,...){
	if(is.null(a) & length(x$a)>1){
		stop("Need to specify value for a")
	} else if(is.null(a)){
		a <- x$a
	}
	if(is.null(lambdaIndex) & is.null(lambda)){
		stop("need to specify a value for lambda")
	}
	if(!is.null(lambdaIndex) & !is.null(lambda)){
		stop("only specify lambdaIndex or lambda")
	}
	if(is.null(lambdaIndex)){
		lambdaIndex <- which(x$lambda==lambda)
	}
	if(length(lambdaIndex)>1){
		stop("Function only supports a single value of lambda or lambdaIndex")
	}
	coefs <- do.call(rbind,coefficients(x,a=a,lambdaIndex=lambdaIndex))
	par(ask=TRUE)
	p <- ncol(coefs)
	tau <- x$tau
	for(i in 1:p){
		plot(tau,coefs[,i],xlab=expression(tau),ylab="Coefficient",main=colnames(coefs)[i],pch=16,...)
		lines(tau,coefs[,i])
	}
	par(ask=FALSE)
}

bytau.plot.rq.pen.seq.cv <- function(x,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,...){
	coefs <- do.call(rbind,coefficients(x,septau,cvmin,useDefaults,tau=x$fit$tau,...))
	if(nrow(coefs) != length(x$fit$tau)){
		stop("Too many coefficients returned, function only works with one lambda value")
	}
	par(ask=TRUE)
	p <- ncol(coefs)
	tau <- x$fit$tau
	for(i in 1:p){
		plot(tau,coefs[,i],xlab=expression(tau),ylab="Coefficient",main=colnames(coefs)[i],pch=16,...)
		lines(tau,coefs[,i])
	}
	par(ask=FALSE)
	
}

coef.rq.pen.seq.cv <- function(x,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,tau=NULL,...){
	if(!useDefaults){
		coefficients(x$models,tau=tau,...)
	} else{
		if(is.null(tau)){
			tau <- m1$fit$tau
		}
		if(septau){
			keepers <- which(closeEnough(tau,x$btr$tau))
			btr <- x$btr[keepers,]
			models <- x$fit$models[btr$modelsIndex]
			if(cvmin){
				lambdaIndex <- btr$lambdaIndex
			} else{
				lambdaIndex <- btr$lambda1seIndex
			}
			returnVal <- vector(mode="list", length=length(models))
			names(returnVal) <- names(models)
			for(i in 1:length(returnVal)){
				returnVal[[i]] <- coef(x$fit,modelsIndex=btr$modelsIndex[i],lambdaIndex=lambdaIndex[i])[[1]]
			}
			returnVal
		} else{
			if(!cvmin){
				stop("One standard error approach not implemented for group choice of tuning parameter")
			} else{
				keepers <- which(closeEnough(tau,x$gtr$tau))
				gtr <- x$gtr[keepers,]#subset(x$gtr, closeEnough(tau,x$gtr$tau))
				coef(x$fit,modelsIndex=gtr$modelsIndex,lambdaIndex=gtr$lambdaIndex[1])
			}
		}
	}
}

cv_plots <- function(model,logLambda=TRUE,loi=NULL,...){
#logLambda - lambdas plotted on log scale
#loi - index of target lambdas
  cv_data <- model$cv
  if(logLambda){
     cv_data$lambda <- log(cv_data$lambda)
     colnames(cv_data)[1] <- "logLambda"
  }
  if(is.null(loi)==FALSE){
     cv_data <- cv_data[loi,]                                        
  }                                    
  plot(cv_data[,1], cv_data[,2], ylab=colnames(cv_data)[2], xlab=colnames(cv_data)[1],...)
}

cv.rq.group.pen <- function (x, y, groups, tau = 0.5, lambda = NULL, penalty = "SCAD", 
    intercept = TRUE, criteria = "CV", cvFunc = "check", nfolds = 10, 
    foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="QICD",penGroups=NULL,
    ...) 
{
  warning("Recommend that you use rq.group.pen.cv() instead. This is an older and slower version that is only kept for reproducibality. It will not be exported to the namespace in future versions.")
  if(penalty=="LASSO"){
	warning("The Lasso group penalties use the L1 norm and thus the Lasso group penalty is the same as the standard Lasso penalty and therefore does not account for group structure. The group lasso method is only implemented because it is needed for the SCAD and MCP algorithms. Otherwise it should be avoided. ")
  }
	if(is.null(penGroups)){
		p_range <- 1:dim(x)[2] + intercept
    } else{
		p_range <- which(groups %in% penGroups) + intercept
	}
  n <- dim(x)[1]
  pen_func <- switch(which(c("LASSO", "SCAD", "MCP") == penalty), 
        lasso, scad, mcp)
  if(alg=="QICD" & penalty != "LASSO"){
    if(criteria== "CV"){
       stop("QICD algorithm wtih non-convex penalties setup only to use BIC or PBIC as the criteria")
    }
    #start with lasso fit
    lasso_fit <- cv.rq.group.pen.old(x,y,groups,tau,lambda,penalty="LASSO",intercept,criteria="BIC",cvFunc,nfolds,foldid,
                                   nlambda, eps, init.lambda,alg="LP",penGroups,...)
    #then iterate through lasso models to get new models
    model_coefs <- NULL
    lambda_vals <- lasso_fit$cv$lambda
    lasso_beta <- lasso_fit$beta
    for(model_num in 1:dim(lasso_beta)[2]){
       model_coefs <- cbind(model_coefs, QICD.group(y, x, groups, tau, lambda_vals[model_num], intercept, penalty, 
                 initial_beta=lasso_beta[,model_num], eps = eps,...))
    }
    return_val <- NULL
    return_val$beta <- model_coefs
    if(intercept){
       fits <- cbind(1,x) %*% model_coefs
    } else{
       fits <- x %*% model_coefs
    }
    return_val$residuals <- y - fits
    return_val$rho <- apply(check(return_val$residuals,tau),2,sum)
    non_zero_count <- apply(model_coefs!=0,2,sum)
    if( criteria=="BIC" ){
      cve <- log(return_val$rho) + non_zero_count*log(n)/(2*n)
    } else { # PBIC
      cve <- log(return_val$rho) + non_zero_count*log(n)*log(nrow(model_coefs))/(2*n)
    }
       
    return_val$cv <- data.frame(lambda = lambda_vals, cve = cve)
    colnames(return_val$cv)[2] <- criteria
    return_val$lambda.min <- lambda_vals[which.min(cve)]
    return_val$penalty <- penalty
    return_val$intercept <- intercept
    return_val$groups <- groups
    class(return_val) <- c("cv.rq.group.pen", "cv.rq.pen")    
  }      
  else{
    if (is.null(lambda)) {
        sample_q <- quantile(y, tau)
        inter_only_rho <- sum(check(y - sample_q, tau))
        lambda_star <- init.lambda
        searching <- TRUE
        search_num <- 1
        while(searching){
            if (search_num == 1) {
                init_fit <- rq.group.fit(x, y, groups, tau, lambda_star, 
                  intercept, penalty,alg,penGroups=penGroups, ...)
                search_num <- 2
            }
            else {
                init_fit <- rq.group.fit(x, y, groups, tau, lambda_star, 
                  intercept, penalty,alg, initial_beta = beta_update, penGroups=penGroups,
                  ...)
            }
            beta_update <- init_fit$coefficients
            if (sum(init_fit$coefficients[p_range]) == 0) {
                searching <- FALSE
            }
            else {
				option_1 <- (inter_only_rho-init_fit$rho)/ sum(sapply(init_fit$coefficients[p_range],pen_func, 1))
				option_2 <- lambda_star*search_num
				if(option_1 >= option_2){
					lambda_star <- option_1
				} else{
					lambda_star <- option_2
					search_num <- search_num + 1
				}
                # lambda_star = max((inter_only_rho-init_fit$rho)/ sum(sapply(init_fit$coefficients[p_range],pen_func, 1)),
                #                    lambda_star+1)
                 #this is sort of a hack, need to think of a better way to pick lambda_star
                #lambda_star = inter_only_rho/sum(sapply(init_fit$coefficients[p_range], 
                #  pen_func, 1))
                #weird behavior if we don't set penalty lambda to 1
                #in general this idea needs to be improved upon
            }
       }
          
       lambda_min <- eps * lambda_star
       lambda <- exp(seq(log(lambda_star), log(lambda_min),  length.out = nlambda))
       fine_tune <- TRUE
       fine_tune_pos <- length(lambda)-1
       while(fine_tune){
         test_fit <- rq.group.fit(x,y,groups,tau,lambda[fine_tune_pos],intercept,penalty,alg,penGroups=penGroups,...)
         if (sum(test_fit$coefficients[p_range]) != 0) {
         #want only one intercept only model
                  fine_tune <- FALSE
         } else{
           fine_tune_pos <- fine_tune_pos - 1
         }
       }
       if(fine_tune_pos != (length(lambda)-1)){
         lambda_star <- lambda[fine_tune_pos+1]
         lambda_min <- eps*lambda_star
         lambda <- exp(seq(log(lambda_star), log(lambda_min),  length.out = nlambda))
       }
    }
    
    models <- groupMultLambda(x = x, y = y, groups = groups, 
        tau = tau, lambda = lambda, intercept = intercept, penalty=penalty,alg=alg,penGroups=penGroups, ...)
    model_coefs  <- sapply(models,coefficients)
    model_resids <- sapply(models,residuals)
    model_rhos <- sapply(models, getRho)
    
    cv_results <- NULL
    if (criteria == "CV") {
        if (is.null(foldid)) {
            foldid <- randomly_assign(n, nfolds)
        }
        for (i in 1:nfolds) {
            train_x <- x[foldid != i, ]
            train_y <- y[foldid != i]
            #test_x <- x[foldid == i, ]
			test_x <- x[foldid==i,,drop=FALSE]
            test_y <- y[foldid == i]
            cv_models <- groupMultLambda(x = train_x, y = train_y, 
                groups = groups, tau = tau, lambda = lambda, 
                intercept = intercept,penalty=penalty,alg=alg,penGroups=penGroups, ...)
            if (cvFunc == "check") {
                cv_results <- cbind(cv_results, sapply(cv_models, 
                  model_eval, test_x, test_y, tau = tau))
            }
            else {
                cv_results <- cbind(cv_results, sapply(cv_models, 
                  model_eval, test_x, test_y, func = cvFunc))
            }
        }
        cv_results <- apply(cv_results, 1, mean)
    }
    if (criteria == "BIC") {
        cv_results <- sapply(models, qbic)
    }
    if (criteria == "PBIC") {
        cv_results <- sapply(models, qbic, largeP = TRUE)
    }
    lambda.min <- lambda[which.min(cv_results)]
    return_val <- NULL
    #return_val$models <- models
    return_val$beta <- model_coefs
    return_val$residuals <- model_resids
    return_val$rho <- model_rhos   
    return_val$cv <- data.frame(lambda = lambda, cve = cv_results)
    colnames(return_val$cv)[2] <- criteria
    return_val$lambda.min <- lambda.min
    return_val$penalty <- penalty
    return_val$intercept <- intercept
    return_val$groups <- groups
    class(return_val) <- c("cv.rq.group.pen", "cv.rq.pen")
  }
  return_val
}

rq.group.fit <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                penalty = "SCAD", alg="QICD", a=3.7,penGroups=NULL, ...) 
{
  ### Some cleaning/checking before getting to the algorithms
  p <- ncol(x)
  n <- nrow(x)
  #if(is.null(penGroups) & max(penGroups) > max(groups)){ stop("penalize groups not coefficients")}  
  if (!penalty %in% c("LASSO", "SCAD", "MCP")) {
      stop("Penalty must be LASSO, SCAD or MCP")
  }
  if(penalty=="LASSO"){
	warning("Group penalties use the L1 norm and the Lasso group penalty is the same as the standard Lasso penalty and therefore does not account for group structure. The group lasso method is only implemented because it is needed for the SCAD and MCP algorithms. Otherwise it should be avoided. ")
  }
  if(is.null(dim(x))){ stop("x must be matrix with at least 1 column") }
  if(length(groups)!=ncol(x)){
    stop("length(groups) must be equal to ncol(x)")
  }
  if( lambda <= 0 ){ stop("lambda must be positive")}
  
  if(penalty=="LASSO"){
	pen_func <- lasso
  }
  if(penalty=="SCAD"){
	pen_func <- scad
  }
  if(penalty=="MCP"){
	pen_func <- mcp
  }

  if (alg == "QICD") {
  ### QICD Algorithm ###
    if( length(lambda) != 1 )
      stop( "QICD Algorithm only allows 1 lambda value")

    coefs <- QICD.group(y, x, groups, tau, lambda, intercept, penalty,a=a, ...)

    ### Add extra information to QICD output
    coefnames <- paste("x",1:p, sep="") ### Coefficient names
    if(intercept)
      coefnames <- c("(Intercept)", coefnames)
    names(coefs) <- coefnames
    if( intercept ){ ### Residuals
      residuals <- c( y - x%*%(coefs[-1]) - coefs[1] )
	  pen_vars <- coefs[-1]
    } else {
      residuals <- c( y - x%*%coefs )
	  pen_vars <- coefs
    }
	if(penalty=="LASSO"){
		pen_val <- sum(pen_func(tapply(abs(pen_vars),groups,sum),lambda=lambda))
	} else{
		pen_val <- sum(pen_func(tapply(abs(pen_vars),groups,sum),lambda=lambda,a=a))
	}
	
    rho <- sum( check(residuals) ) # rho (similiar to quantreg)
		#1/n*sum( check(residuals) ) ### rho
	PenRho <- rho+pen_val

    return_val <- list( coefficients=coefs, PenRho=PenRho, residuals=residuals, rho=rho, tau=tau, n=n, intercept=intercept, penalty=penalty)
	class(return_val) <- c("rq.group.pen", "rq.pen")
  ######################
  } else {
        group_num <- length(unique(groups))
        if (length(lambda) == 1) {
            lambda <- rep(lambda, group_num)
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
          
        if (penalty == "LASSO") {			
			new_lambda <- NULL
			group_count <- xtabs(~groups)
			for (g in 1:group_num) {
				 new_lambda <- c(new_lambda, rep(lambda[g], each = group_count[g]))
			}
			if(is.null(penGroups)==FALSE){
				new_lambda[-which(groups %in% penGroups)] <- 0
			}
            return_val <- rq.lasso.fit(x, y, tau, new_lambda, 
                intercept = intercept, ...)
            class(return_val) <- c("rq.group.pen", "rq.pen", 
                "rqLASSO")
        }
        else {
            return_val <- rq.group.lin.prog(x,y,groups,tau,lambda,intercept=intercept,penalty=penalty,penGroups=penGroups,a=a,...)
            class(return_val) <- c("rq.group.pen", "rq.pen")
        }
    }
    return_val
}

plot.cv.rq.group.pen <- function (x,y=NULL,...) 
{
    plot(x$cv[, 1], x$cv[, 2])
}

rq.lasso.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                         coef.cutoff=1e-08, method="br",penVars=NULL,scalex=TRUE,lambda.discard=TRUE, ...){
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


coef.cv.rq.group.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  object$beta[,target_model]
}


rq.group.pen <- function(x,y, tau=.5,groups=1:ncol(x), penalty=c("gLASSO","gAdLASSO","gSCAD","gMCP"),
						lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.01,.0001),alg=c("huber","lp","qicd"), 
						a=NULL, norm=2, group.pen.factor=rep(1,length(unique(groups))),tau.penalty.factor=rep(1,length(tau)),
						scalex=TRUE,coef.cutoff=1e-8,max.iter=10000,converge.eps=1e-7,gamma=IQR(y)/10, ...){
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
	
	if(lpf!=g){
		stop("group penalty factor must be of length g")
	}	
	
	if(is.null(lambda)){
		lamMax <- getLamMaxGroup(x,y,groups,tau,group.pen.factor,penalty=penalty,scalex)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}

	penalty.factor <- mapvalues(groups,seq(1,g),group.pen.factor)
	
	if(penalty == "gLASSO"){
		return_val <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma,...)
	} else{
		if(penalty == "gAdLASSO"){
			init.model <- rq.enet(x,y,tau,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,tau.penalty.factor=tau.penalty.factor,
									a=0,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,...)
		} else{
			if(norm == 1){
				if(alg == "qicd"){
					init.alg <- "lp"
				} else{
					init.alg <- alg
				}
				init.model <- rq.lasso(x,y,tau,alg=init.alg,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,
							tau.penalty.factor=tau.penalty.factor,max.iter=max.iter,coef.cutoff=coef.cutoff,converge.eps=converge.eps,
							gamma=gamma,...)
			} else{
				init.model <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma,...)
			}
		}
		return_val <- rq.group.lla(init.model,x,y,groups,penalty=penalty,a=a,norm=norm,group.pen.factor=group.pen.factor,tau.penalty.factor=tau.penalty.factor,scalex=scalex,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,...)
	}
	class(return_val) <- "rq.pen.seq"
	return_val$call <- match.call()	
	return_val$lambda <- lambda
	return_val
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


