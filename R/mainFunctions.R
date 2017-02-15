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



model_eval <- function(model, test_x, test_y, func="check",...){
#func: "check" (Quantile Check), "SqErr" (Squared Error), "AE" (Absolute Value)
  if(model$intercept){
    test_x <- cbind(1,test_x)
  }
  fits <- test_x %*% coefficients(model)
  eval_func <- switch(which(c("check","SqErr","AE")==func), check, square, abs)
  mean(eval_func(test_y-fits,...)) 
}


qbic <- function(model, largeP=FALSE){
  tau <- model$tau
  n <- model$n
  df <- sum(model$coefficients !=0)
  if(largeP){
    bic <- log(model$rho) + df*log(n)*log(length(model$coefficients))/(2*n)
  }else{
    bic <- log(model$rho) + df*log(n)/(2*n)
  }
  bic
}

coef.cv.rq.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  coefficients(object$models[[target_model]])
}



cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",intercept=TRUE,criteria="CV",cvFunc="check",nfolds=10,foldid=NULL,nlambda=100,eps=.0001,init.lambda=1,penVars=NULL,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
# nfolds: number of folds for cross validation
# foldid: preset id of folds
# penVar: variables to be penalized, default is all non-intercept terms
  p <- dim(x)[2]
  if(is.null(penVars)){
    penVars <- 1:p
  }
  p_range <- penVars + intercept
  n <- dim(x)[1]
  pen_func <- switch(which(c("LASSO","SCAD","MCP")==penalty), lasso, scad, mcp)
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
      test_x <- x[foldid==i,]
      test_y <- y[foldid==i]
      if(penalty=="LASSO"){
         cv_models <- lapply(lambda,rq.lasso.fit, x=train_x,y=train_y,tau=tau,weights=weights,intercept=intercept,penVars=penVars,...)
      } else{
         cv_models <- lapply(lambda,rq.nc.fit, x=train_x,y=train_y,tau=tau,weights=weights,intercept=intercept,penalty=penalty,penVars=penVars,...)
      }
      if(cvFunc=="check"){
         cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, test_y, tau=tau))
      } else{
         cv_results <- cbind(cv_results, sapply(cv_models,model_eval, test_x, test_y, func=cvFunc))
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
                      penalty="SCAD",a=3.7,iterations=10,converge_criteria=1e-06,
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
  betas <- t(sapply(model$models, coefficients))
  if(is.null(voi)==FALSE){
    betas <- betas[,voi]
  }
  if(colnames(betas)[1]=="intercept"){
    betas <- betas[,-1]
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

cv.rq.group.pen <- function (x, y, groups, tau = 0.5, lambda = NULL, penalty = "LASSO", 
    intercept = TRUE, criteria = "CV", cvFunc = "check", nfolds = 10, 
    foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="QICD",penGroups=NULL,
    ...) 
{
	if(is.null(penGroups)){
		p_range <- 1:dim(x)[2] + intercept
    } else{
		p_range <- which(groups %in% penGroups) + intercept
	}
    n <- dim(x)[1]
    pen_func <- switch(which(c("LASSO", "SCAD", "MCP") == penalty), 
        lasso, scad, mcp)
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
            test_x <- x[foldid == i, ]
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
    return_val
}

rq.group.fit <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                penalty = "LASSO", alg="QICD", a=3.7,penGroups=NULL, ...) 
{
  ### Some cleaning/checking before getting to the algorithms
  p <- ncol(x)
  n <- nrow(x)
  #if(is.null(penGroups) & max(penGroups) > max(groups)){ stop("penalize groups not coefficients")}  
  if (!penalty %in% c("LASSO", "SCAD", "MCP")) {
      stop("Penalty must be LASSO, SCAD or MCP")
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

    coefs <- QICD.group(y, x, groups, tau, lambda, intercept, penalty, ...)

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
            return_val <- rq.group.lin.prog(x,y,groups,tau,lambda,intercept=intercept,penalty=penalty,penGroups=penGroups,...)
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
	} else{
		eval_function <- function(x,y,tau,...){ 
							 b <- bs(x,...)
							 q1 <- rq(y ~ b, tau)
							 sum((fitted(q1)-quantile(y,tau))^2)
						 }
	}
	#if(n.cores==1){
		eval_results <- apply(x,2,eval_function,y,tau,...)
	#} else{
	#	p <- dim(x)[2]
	#	mc_func <- function(idx,...){ eval_function(x[,idx],y,...)}
	#	mc_results <- mclapply(1:p, mc_func, mc.cores=n.cores, ...)
	#	eval_results <- do.call(c,mc_results)
	#}
	order( eval_results, decreasing=TRUE)
}


