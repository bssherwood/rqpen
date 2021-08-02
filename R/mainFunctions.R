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
	if(class(obj) == "cv.rq.pen.seq"){
		obj <- obj$fit
	} else if(class(obj) != "rq.pen.seq"){
		stop("obj must be of class rq.pen.seq or cv.rq.pen.seq")
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
	nl <- length(obj$models[[1]]$lambda)
	
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

# print.rq.pen.seq <- function(x,...){
	# if(length(x$tau)==1){
		# cat("\n df by lambda:\n")
		# print(data.frame(lambda=x$lambda,df=x$models$df))
	# } else{
		# cat("\n df by lambda:\n")
		# printdata <- NULL
		# for(i in 1:length(x$tau)){
			# printdata <- cbind(printdata,x$models[[i]]$df) 
		# }
		# printdata <- cbind(x$lambda,printdata)
		# printdata <- data.frame(printdata)
		# colnames(printdata) <- c("lambda",x$tau)
	# }
# }


coef.cv.rq.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  coefficients(object$models[[target_model]])
}


rq.pen <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,...){
	penalty <- match.arg(penalty)
	if(penalty=="LASSO"){
		fit <- rq.lasso(x,y,tau,lambda,...)
	} else if(penalty=="Ridge"){
		fit <- rq.enet(x,y,tau,lambda,a=0,...)
	} else if(penalty == "ENet"){
		fit <- rq.enet(x,y,tau,lambda,a=a,...)
	} else if(penalty == "aLASSO" | penalty=="SCAD" | penalty == "MCP"){
		fit <- rq.nc(x,y,tau,penalty,a,lambda,...)
	}
	fit
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

subtract <- function(predicted,obs){
	obs-predicted
}

modelTau <- function(object){
	object$tau
}

modelA <- function(object){
	object$a	
}

modelLambda <- function(object, index){
	object$lambda[index]
}

modelNz <- function(object, index){
	object$nz[index]
}

byTauResults <- function(cvErr,tauvals,avals,models,se){
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
		subLambda <- models[[btr[[5]][i]]]$lambda[btr[[3]][i]]
		lambdaVals <- c(lambdaVals, subLambda)
		subLambda1sePos <- which(cvErr[btr[[5]][1],] <= se1Above)[1]
		lambda1seIndex <- c(lambda1seIndex,subLambda1sePos)
		subLambda1se <- models[[btr[[5]][i]]]$lambda[subLambda1sePos]
		lambda1se <- c(lambda1se,subLambda1se)
		nz <- c(nz, models[[btr[[5]][i]]]$nz[btr[[3]][i]])
	}
	
	btr <- cbind(btr, lambda=lambdaVals, cvse = cvse, lambda1se=lambda1se, lambda1seIndex=lambda1seIndex, nonzero=nz)
	btr <- setcolorder(btr, c(1,2,6,3,8,9,4,7,5,10))
	btr
}

groupTauResults <- function(cvErr, tauvals,a,avals,models,tauWeights){
# improve note: maybe code in for loop could be improved upon by checking at each iteration if the better cv value has been found or not
	nl <- length(models[[1]]$lambda)
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
	lambdavals <- sapply(targetModels,modelLambda,minIndex[1,2])
	nz <- sapply(targetModels, modelNz, minIndex[1,2])
	minCv <- cvErr[modelIndex,minIndex[1,2]]
	list(returnTable=data.table(tau=tauvals,lambda=lambdavals,a=returnA,minCv=minCv,lambdaIndex=minIndex[1,2],modelsIndex=modelIndex, nonzero=nz),gcve=gcve)
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
	nl <- length(fit$models[[1]]$lambda)
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
		btr <- byTauResults(foldErrors,tauvals,avals,fit$models,stdErr)
		gtr <- groupTauResults(foldErrors, tauvals,fit$a,avals,fit$models,tauWeights)
	} else{
		indErrors <- t(indErrors)/n
		btr <- byTauResults(indErrors,tauvals,avals,fit$models,stdErr)
		gtr <- groupTauResults(indErrors, tauvals,fit$a,avals,fit$models,tauWeights)
	}

	returnVal <- list(cverr = foldErrors, cvse = stdErr, fit = fit, btr=btr, gtr=gtr$returnTable, gcve=gtr$gcve)
	class(returnVal) <- "cv.rq.pen.seq"
	returnVal
}

print.cv.rq.pen.seq <- function(x,...){
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

coef.cv.rq.pen.seq <- function(x,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,tau=NULL,...){
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

# coef.cv.rq.pen.seq <- function(x,tau=NULL,lambda=NULL,a=NULL,tauType=c("indTau","groupTau"),cvCrit=c("min","1se")){
	# allNull <- is.null(tau) & is.null(lambda) & is.null(a)
	# lt <- length(tau)
	# la <- length(a)
	# ll <- length(lambda)
	# if(is.null(tau) | is.null(lambda) | is.null(a) & !allNull){
		# stop("All of tau, lambda and a need to be set if defaults are not used")
	# }
	# if( ll != 1 & la !=1 & ll != lt & la != lt){
		# stop("length of lambda and a must be one or the same as tau")
	# }
	# if(allNull){
		# tauType <- match.arg(default)
		# cvCrit <- match.arg(cvCrit)
		# if(tauType=="indTau"){
			
		# } else{
			# if(cvCrit=="1se"){
				# stop("One standard error approach not implemented for group choice of tuning parameter")
			# }
		# }
	# }
# }


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
	nl <- length(fit$models[[1]]$lambda)
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
	class(returnVal) <- "cv.rq.pen.seq"
	returnVal
}


cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",criteria = "CV",intercept=TRUE,cvFunc="check",nfolds=10,foldid=NULL,nlambda=100,eps=.0001,init.lambda=1,penVars=NULL,alg=ifelse(ncol(x)<50,"LP","QICD"),...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
# nfolds: number of folds for cross validation
# foldid: preset id of folds
# penVar: variables to be penalized, default is all non-intercept terms

  # Pre-algorithm setup/ get convenient values
  warning("Recommend using rq.pen.cv instead. This is an older function that is kept for reproducibality reasons, but tends to be much slower. It will not be exported to the namespace in future versions")
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
         init_fit <- rq.nc.fit(x,y,tau,lambda=lambda_star,weights,intercept,penVars=penVars,...)
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


rq.nc.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                      penalty="SCAD",a=3.7,iterations=1,converge_criteria=1e-06,
                      alg=ifelse(p<50,"LP","QICD"),penVars=NULL,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# lambda takes values of 1 or p
# penalty SCAD or MCP
# penVars - variables to be penalized, doesn't work if lambda has multiple entries
  warning("Recommend using rq.pen() instead. This is an older function and usually will run slower. It will not be exported to the namespace in future versions.")
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
		lambdaIndex <- whichMatch(lambda,x$models[[1]]$lambda)
	} else if(is.null(lambdaIndex)){
		lambdaIndex <- 1:length(x$models[[1]]$lambda)
	}
	if(length(lambdaIndex)==0){
		stop("Invalid lambda provided")
	}
	targetModels <- x$models[modelsIndex]
	list(targetModels=targetModels,lambdaIndex=lambdaIndex,modelsIndex=modelsIndex)
}

coef.rq.pen.seq <- function(x,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL){
	models <- getModels(x,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL)
	lapply(models$targetModels,getModelCoefs,models$lambdaIndex)
}

plot.cv.rq.pen.seq <- function(x,septau=TRUE,tau=NULL,logLambda=FALSE,main=NULL,...){
	if(septau){
		plotsep.cv.rq.pen.seq(x,tau,logLambda,main,...)
	} else{
		if(is.null(tau)==FALSE){
			stop("Tau cannot be set if septau set to FALSE")
		}
		plotgroup.cv.rq.pen.seq(x,logLambda,main,...)
	}
}

plotgroup.cv.rq.pen.seq <- function(x,logLambda,main,...){
	a <- x$fit$a
	na <- length(a)
	# besta <- x$gtr$a[1]
	# bestaidx <- which(a==besta)
	# a <- a[-bestaidx]
	# a <- c(besta,a)
	if(logLambda){
		lambdas <- log(x$fit$models[[1]]$lambda)
		xtext <- expression(Log(lambda))
	} else{
		lambdas <- x$fit$models[[1]]$lambda
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

# plotgroup.cv.rq.pen.seq <- function(x,a,logLambda,main,...){
# #code challenge implicitly assumes lambda is the same for all models. 
	# if(is.null(a)){
		# a <- x$fit$a
	# }
	# na <- length(a)
	# if(na > 1){
		# par(ask=TRUE)
	# }
	# for(i in 1:na){
		# aindex <- which(x$fit$a == a[i])
		# if(is.null(main)){
			# mainText <- paste("Group cross validation results for a = ",a[i])
		# } else if(length(main)==1){
			# mainText <- main
		# } else{
			# main <- main[i]
		# }
		# if(logLambda){
			# lambdas <- log(x$fit$models[[1]]$lambda)
			# xtext <- expression(Log(lambda))
		# } else{
			# lambdas <- x$fit$models[[1]]$lambda
			# xtext <- expression(lambda)
		# }
		# plot(lambdas, m1$gcve[aindex,],ylab="Cross Validation Error",xlab=xtext,main=mainText,col="red",pch=16,...)
		# bestIndex <- subset(x$gtr, a==a[i])$lambdaIndex[1]
		# lines(rep(lambdas[bestIndex],2),c(-5,max(m1$gcve[aindex,])+1),lty=2)
	# }
	# if(na > 1){
		# par(ask=FALSE)
	# }
# }



error.bars <- function (x, upper, lower, width = 0.02, ...) 
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

plotsep.cv.rq.pen.seq <- function(x,tau,logLambda,main,...){
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
		lambdas <- log(x$fit$models[[1]]$lambda)
		xtext <- expression(Log(lambda))
	} else{
		lambdas <- x$fit$models[[1]]$lambda
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
		cvsd <- x$cvse[bestidx,]
		
		plot(lambdas,suberr[1,],ylim=c(0,max(max(suberr),max(besterr+cvsd))),ylab="Cross Validation Error", xlab=xtext,main=mainText,type="n",...)
		for(j in 1:na){
			points(lambdas,suberr[j,],col=j,...)
		}
		
		#then get index and plot error bars for this. 
		#get the best and create the segments for that
		segments(lambdas,besterr-cvsd,lambdas,besterr+cvsd)
		segments(lambdas-.01,besterr-cvsd,lambdas+.01,besterr-cvsd)
		segments(lambdas-.01,besterr+cvsd,lambdas+.01,besterr+cvsd)
		if(na>1){
			legend("topleft",paste("a=",subinfo$a),col=1:na)
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

# plotsep.cv.rq.pen.seq <- function(x,tau=NULL,a=NULL,modelsIndex=NULL,logLambda=FALSE,main=NULL,...){
	# models <- getModels(x$fit,tau=tau,a=a,modelsIndex=modelsIndex)
	# tm <- models$targetModels
	# li <- models$lambdaIndex
	# mi <- models$modelsIndex
	# ml <- length(tm)
	# if(ml > 1){
		# par(ask=TRUE)
	# }
	# for(i in 1:ml){
		# if(is.null(main)){
			# mainText <- paste("Cross validation results for ", expression(tau)," = ", tm[[i]]$tau, " and a = ", tm[[i]]$a, "\n")
		# } else if(length(main)==1){
			# mainText <- paste(main, "\n")
		# } else{
			# main <- paste(main[i], "\n")
		# }
		# if(logLambda){
			# lambdas <- log(tm[[i]]$lambda[li])
			# xtext <- expression(Log(lambda))
		# } else{
			# lambdas <- tm[[i]]$lambda[li]
			# xtext <- expression(lambda)
		# }
		# err <- x$cverr[mi[i],]
		# cvsd <- x$cvse[mi[i],]
		# plot(lambdas, err, ylim=c(0,max(err+cvsd)),ylab="Cross Validation Error",xlab=xtext,main=mainText,col="red",pch=16,...)
		# segments(lambdas,err-cvsd,lambdas,err+cvsd)
		# segments(lambdas-.01,err-cvsd,lambdas+.01,err-cvsd)
		# segments(lambdas-.01,err+cvsd,lambdas+.01,err+cvsd)
		# mInfo <- subset(x$btr, modelsIndex == mi[i])
		# lines(rep(lambdas[mInfo$lambdaIndex],2),c(-5,max(err+cvsd)+1),lty=2)
		# lines(rep(lambdas[mInfo$lambda1seIndex],2),c(-5,max(err+cvsd)+1),lty=2)
		# axis(side=3, at = lambdas, labels=paste(tm[[i]]$nzero),tick=FALSE,line=0)
	# }
	# if(ml > 1){
		# par(ask=FALSE)	
	# }
# }

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
			lambdas <- log(tm[[i]]$lambda[li])
			xtext <- expression(Log(lambda))
		} else{
			lambdas <- tm[[i]]$lambda[li]
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

getRho <- function(model){
    model$rho
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


