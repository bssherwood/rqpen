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

qic.select <- function(obj, method="BIC",septau=FALSE,weights=rep(1,length(obj$tau))){
	n <- obj$n
	if(length(weights)==1){
		if(septau){
			warning("septau set to TRUE, but only one quantile modeled")
		}
		qic_vals <- qic(obj$models,n,method)
		min_spot <- which.min(qic_vals)
		coefs <- coefficients(obj)[,min_spot]
	} else{
		qic_vals <- sapply(obj$models,qic,n,method)
		if(septau){
			min_vals <- apply(qic_vals,2,which.min)
			coefs <- coefficients(obj,min_vals)
		} else{
			qic_vals <- apply(qic_vals %*% diag(weights),1,sum)
			coefs <- coefficients(obj,which.min(qic_vals))
		}
	}
	return_val <- list(coefficients = coefs, ic=qic_vals, lambda=obj$lambda, penalty.factor=obj$penalty.factor)
	class(return_val) <- "qic.select"
	return_val
}


print.qic.select <- function(x,...){
   cat("\n IC values by Lambda:\n")
   print(data.frame(lambda=x$lambda,IC=x$ic))
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



# coef.rq.pen.seq <- function(object, lambda=NULL, tau=NULL){
	# nt <- length(ojb$tau)
	# if(nt == 1){	
		# coefficients(object$models)
	# } else{
		# lapply(object$models, coefficients)
	# }
# }

rq.pen <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","enet","aLASSO","SCAD","MCP"),a=NULL,...){
	penalty <- match.arg(penalty)
	if(penalty=="LASSO"){
		fit <- rq.lasso(x,y,tau,lambda,...)
	} else if(penalty == "enet"){
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

cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty=c("LASSO","ridge","enet","aLASSO","SCAD","MCP"),a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,nlambda=100,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),...){
#need to think about how to handle this for multi vs one tau. Also multi-a vs single a. Do the four types or something like that and then run the code
	n <- length(y)
	if(is.null(weights)==FALSE){
		stop("weights not currently implemented. Can use cv.rq.pen.old, but it supports fewer penalties and is slower.")
	}
	if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
    }
	fit <- rq.pen(x,y,tau,lambda=lambda,penalty=penalty,a=a,...)
	if(!groupError){
		indErrors <- vector(mode="list",length=nfolds)
	}
	nt <- length(tau)
	na <- length(fit$a)
	nl <- length(fit$models[[1]]$lambda)
	foldErrors <- fe2ndMoment <- matrix(rep(0,nt*na*nl),ncol=nl)
    for(i in 1:nfolds){
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
			indErrors[[i]] <- testErrors # will this be a problem? A list of a list? 
		}
		foldMeans <- do.call(rbind, lapply(testErrors,apply,2,cvSummary))
		foldErrors <- foldErrors + foldMeans
		fe2ndMoment <- fe2ndMoment + foldMeans^2
    }
	fe2ndMoment <- fe2ndMoment/nfolds
	foldErrors <- foldErrors/nfolds
	stdErr <- sqrt( (nfold/(nfold-1))*(fe2ndMoment - foldErrors^2))
    cv_results <- apply(cv_results,1,mean)
	tauvals <- sapply(fit$models,modelTau)
	avals <- sapply(fit$models,modelA)
	#Then get the min by this index or something....
}



cv.rq.pen.old <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",criteria = "CV",intercept=TRUE,cvFunc="check",nfolds=10,foldid=NULL,nlambda=100,eps=.0001,init.lambda=1,penVars=NULL,alg=ifelse(ncol(x)<50,"LP","QICD"),...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
# nfolds: number of folds for cross validation
# foldid: preset id of folds
# penVar: variables to be penalized, default is all non-intercept terms

  # Pre-algorithm setup/ get convenient values
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




# cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",
                      # criteria = "CV",
                      # intercept=TRUE,cvFunc="check",nfolds=10,foldid=NULL,
                      # nlambda=100,eps=.0001,init.lambda=1,penVars=NULL,
                      # alg=ifelse(ncol(x)<50,"LP","QICD"),...){
# # x is a n x p matrix without the intercept term
# # y is a n x 1 vector
# # criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
# # nfolds: number of folds for cross validation
# # foldid: preset id of folds
# # penVar: variables to be penalized, default is all non-intercept terms

  # # Pre-algorithm setup/ get convenient values
  # m.c <- match.call() # This stores all the arguments in the function call as a list

  # p <- dim(x)[2]
  # if(is.null(penVars)){
    # penVars <- 1:p
  # }
  # p_range <- penVars + intercept
  # n <- dim(x)[1]
  # pen_func <- switch(which(c("LASSO","SCAD","MCP")==penalty), lasso, scad, mcp)
 
  # ### QICD ###
  # if( alg=="QICD" & penalty!="LASSO" ){
    # if(criteria=="CV"){
      # stop("CV criteria not implemented for QICD algorithm with nonconvex penalties. Please use BIC or PBIC instead")
    # }
    # m.c[["alg"]] <- "LP" #maybe this should be moved inside the is.null initial beta if statement. I don't think it matters, but might be cleaner code
    # penname <- penalty

    # if( !all(penVars==1:p) ){ # Some unpenalized coefficients
      # z    <- as.matrix(x[,-penVars])
      # xpen <- as.matrix(x[,penVars])
      # QICD_func <- "QICD.nonpen"
      # mapback <- order( c(penVars, (1:p)[-penVars]) ) # reorders the coefficients properly if some (non-intercept) coefficients are not penalized 
      # if( intercept )
        # mapback <- c(1, 1+mapback)
    # } else { # All penalized coefficients
      # z <- NULL
      # xpen <- x
      # QICD_func <- "QICD"
      # mapback <- 1:p # no reordering necessary if all (non-intercept) coefficients are penalized
      # if( intercept )
        # mapback <- c(1, 1+mapback)
    # }

    # # The QICD algorithm needs good starting values, use LASSO solution 
    # ## Speed things up using BIC, not k-fold, to select lambda for LASSO
    # ## Can skip this part if starting values are provided
    # if( is.null(m.c[["initial_beta"]]) ){
      # m.c[["penalty"]] <- "LASSO"
      # m.c[["criteria"]] <- "BIC"
	  # if(is.null(m.c[["lambda"]])==FALSE){
		# m.c[["lambda"]] <- NULL
	  # }
      # suppressWarnings(
        # m.c[["initial_beta"]] <- coefficients( eval.parent(m.c) )
        # # QICD.start <- coefficients( cv.rq.pen(x,y,tau=tau,lambda=lambda,penalty="LASSO",intercept=intercept,criteria="BIC",nlambda=nlambda,eps=eps,init.lambda=lambda,penVars=penVars,...) ) # Use the LASSO with BIC
      # )
    # }

    # # Start in middle of lambda vector
    # ## Increase lambda until intercept only model (or all penlized coefficients are zero)
    # ## Decrease lambda until full model (Sparsity assumption => full model is bad)
    # m.c[[1]] <- as.name(QICD_func)
    # m.c[["penalty"]] <- penalty
    # m.c[["x"]] <- xpen
    # m.c$z <- z

    # if( is.null(lambda) ){
      # lambdas <- exp( seq(-7, 1, length=100) ) 
    # } else {
      # lambdas <- lambda
    # }
    
    # coefs <- matrix(NA, p+intercept, length(lambdas))

    # ### Loop for increasing lambda
    # for( i in floor(length(lambdas)/2):length(lambdas) ){
      # m.c[["lambda"]] <- lambdas[i]
      # coefs[,i] <- eval.parent( m.c )[mapback]
      # # coefs[,i] <- QICD_func( y=y,x=x, z=z, tau=tau, lambda=lambdas[i], intercept=intercept, 
      # #                       penalty=penalty, initial_beta=QICD.start, ... )[mapback]
      # if( all(coefs[p_range,i] == 0) )
        # break()
    # }

    # ### Loop for decreasing lambda
    # for( i in floor(length(lambdas)/2):2 -1 ){
      # m.c[["lambda"]] <- lambdas[i]
      # coefs[,i] <- eval.parent( m.c )[mapback]
      # # coefs[,i] <- QICD_func( y=y,x=x, z=z, tau=tau, lambda=lambdas[i], intercept=intercept, 
      # #                       penalty=penalty, initial_beta=QICD.start, ... )[mapback]
      # if( all(coefs[p_range,i] != 0) )
        # break()
    # }

    # #### Remove the NA columns from coefs and corresponding lambdas
    # lambdas.keep <- which( !is.na(coefs[1,]) )
    # lambdas <- lambdas[lambdas.keep]
    # coefs <- coefs[,lambdas.keep]
    # rownames(coefs) <- names( m.c[["initial_beta"]] )
    # XB <- x%*%coefs[p_range,]
    # if( intercept )
      # XB <- XB + matrix(coefs[1,], n, ncol(coefs), byrow=TRUE)

    # residuals <- y - XB
    # rho <- colSums( check(residuals, tau=tau) )
    # if( is.null(m.c[["a"]]) )
      # a <- 3.7
    # PenRho <- rho + colSums(apply( rbind(lambdas, coefs), 2, 
                    # function(xx) pen_func(xx[1+p_range], lambda=xx[1], a=a) ))

    # cv <- data.frame(lambda=lambdas, cve=NA)
    
    # if( criteria=="BIC" ){
      # cv$cve <- log(rho) + colSums(coefs!=0)*log(n)/(2*n)
      # names(cv)[2] <- "BIC"
    # } else { # PBIC
      # cv$cve <- log(rho) + colSums(coefs!=0)*log(n)*log(nrow(coefs))/(2*n)
      # names(cv)[2] <- "PBIC"
    # }

    # lambda.min <- lambdas[which.min(cv[,2])]

    # # Final cleanup for QICD
    # ## First get models for each lambda
    # models <- vector( "list", length(lambdas) )
    # for( j in 1:length(models) ){
      # models[[j]]$coefficients <- coefs[,j]
      # models[[j]]$PenRho <- PenRho[j]
      # models[[j]]$residuals <- residuals[,j]
      # models[[j]]$rho <- rho[j]
      # models[[j]]$tau <- tau
      # models[[j]]$n <- n
      # models[[j]]$intercept <- intercept
      # models[[j]]$penalty <- penalty
      # class(models[[j]]) <- c("rq.pen", "rqNC")
    # }

    # return_val <- list( models=models, cv=cv, lambda.min=lambda.min, penalty=penalty )
    # class(return_val) <- "cv.rq.pen"
  # } else{
  # ############



  # # If no lambdas provided, find reasonable lambdas to use
  # if(is.null(lambda)){
  # # find a lambda that sets all coefficients to zero. 
  # # Strategy is to get \lambda \sum p_\lambda(|\beta_j}) >= \sum \rho_\tau(y-quantile(y,tau)
  # # do this by fitting model with lambda = init.lambda and then set new lambda such that 
  # # \lambda* = \sum \rho_\tau(y-quantile(y,tau)) / \sum p_\lambda(|beta_j|) repeat as needed
     # sample_q <- quantile(y,tau)
     # inter_only_rho <- sum(check(y-sample_q,tau))
     # #lambda_star <- rep(0,p)
     # #lambda_star[penVars] <- init.lambda
	   # lambda_star <- init.lambda
     # searching <- TRUE
     # while(searching){
       # if(penalty=="LASSO"){
         # init_fit <- rq.lasso.fit(x,y,tau,lambda=lambda_star,weights,intercept,penVars=penVars,...)
       # } else{
         # init_fit <- rq.nc.fit(x,y,tau,lambda=lambda_star,weights,intercept,penVars=penVars,...)
       # }
       # if(sum(init_fit$coefficients[p_range])==0){
         # searching <- FALSE     
       # } else{
         # lambda_star <- inter_only_rho / sum(sapply(init_fit$coefficients[p_range],pen_func,1)) 
         # #1 used here because solving for lambda
       # }
     # }
     # lambda_min <- eps*lambda_star
     # lambda <- exp(seq(log(max(lambda_min)),log(max(lambda_star)),length.out=nlambda))#max is included to handle cases where
     # # some variables are unpenalized and thus lambda is a multivalued vector with some zeros
  # }
  # # lambda is the vector of reasonable choices of lambda to use in the penalty
  # models <- list()
  # fit_models <- TRUE
  # lam_pos = floor(length(lambda)/2)
  # max_lambda_pos <- length(lambda)
  # min_lambda_pos <- 1
  # if(penalty=="LASSO"){
    # while (fit_models) {
      # if (fit_models) {
        # models[[lam_pos]] = rq.lasso.fit(x,y,tau,lambda[lam_pos],
                                         # weights,intercept,
                                         # penVars=penVars,...)
      # }
      # if (sum(abs(coefficients(models[[lam_pos]])[p_range])) != 0) {
        # if(lam_pos > min_lambda_pos) {
          # min_lambda_pos = lam_pos
        # }
        # lam_pos = lam_pos + ceiling((max_lambda_pos - lam_pos)/2)
      # } else {
        # if(lam_pos == max_lambda_pos) {
          # end_pos = lam_pos
          # fit_models = FALSE
        # }
        # if(lam_pos < max_lambda_pos) {
          # max_lambda_pos = lam_pos
        # }
        # lam_pos = lam_pos - floor((lam_pos - min_lambda_pos)/2)
      # }
    # }
    
    # lambda = lambda[1:lam_pos]
  # } else{
  	# while(fit_models){
  		# if(fit_models){
  			# models[[lam_pos]] <- rq.nc.fit(x,y,tau,lambda[lam_pos],weights,intercept,penalty=penalty,penVars=penVars,...)
      # }
  		# if(sum(abs(coefficients(models[[lam_pos]])[p_range]))==0 || lam_pos==length(lambda)){
  			# fit_models <- FALSE
  			# lambda <- lambda[1:lam_pos]
  		# }
      # lam_pos <- lam_pos + 1
  	 # }
     # #models <- lapply(lambda,rq.nc.fit, x=x,y=y,tau=tau,weights=weights,intercept=intercept,penalty=penalty,
     # #                                   penVars=penVars,...)
  # }
  # cv_results <- NULL
  # if(criteria=="CV"){
    # if(is.null(foldid)){
      # foldid <- randomly_assign(n,nfolds)
    # }
    # for(i in 1:nfolds){
      # train_x <- x[foldid!=i,]
      # train_y <- y[foldid!=i]
      # test_x <- x[foldid==i,,drop=FALSE]
      # test_y <- y[foldid==i]
	  # train_weights <- weights[foldid!=i] #not sure this line is needed
	  # if(is.null(weights)){
		# train_weights <- test_weights <- NULL
	  # } else{
	    # train_weights <- weights[foldid!=i]
		# test_weights <- weights[foldid==i]
	  # }
      # if(penalty=="LASSO"){
         # cv_models <- lapply(lambda,rq.lasso.fit, x=train_x,y=train_y,tau=tau,
                             # weights=train_weights,intercept=intercept,penVars=penVars,...)
      # } else{
         # cv_models <- lapply(lambda,rq.nc.fit, x=train_x,y=train_y,tau=tau,
                             # weights=train_weights,intercept=intercept,
                             # penalty=penalty,penVars=penVars,...)
      # }
      # if(cvFunc=="check"){
         # cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, 
                                                # test_y, test_weights, tau=tau))
      # } else{
         # cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, 
                                                # test_y, test_weights, func=cvFunc))
      # } 
    # }
    # cv_results <- apply(cv_results,1,mean)
  # }
  # if(criteria=="BIC"){
    # cv_results <- sapply(models,qbic)
  # }
  # if(criteria=="PBIC"){
    # cv_results <- sapply(models,qbic,largeP=TRUE)
  # }
  # lambda.min <- lambda[which.min(cv_results)]
  # return_val <- NULL
  # return_val$models <- models
  # return_val$cv <- data.frame(lambda=lambda, cve=cv_results)
  # colnames(return_val$cv)[2] <- criteria
  # return_val$lambda.min <- lambda.min
  # return_val$penalty <- penalty
  # class(return_val) <- "cv.rq.pen"
  # }

  # return_val
# }


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
  } else{
	lambdas <- model$cv$lambda
  }
  
  if(is.null(loi)==FALSE){
     lambdas <- lambdas[loi]
  }                                    
  plot(lambdas, betas[,1], type="n",ylim=c(min(betas),max(betas)),ylab="Coefficient Value",xlab="Log Lambda",...)
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



cv.rq.group.pen.old <- function (x, y, groups, tau = 0.5, lambda = NULL, penalty = "SCAD", 
    intercept = TRUE, criteria = "CV", cvFunc = "check", nfolds = 10, 
    foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="QICD",penGroups=NULL,
    ...) 
{
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


