kernel_estimates <- function(x,y,h,...){
  kernel_estimates <- NULL
  if(is.null(dim(x))){
    n <- length(x)
  } else{
    n <- dim(x)[1]
  }
  for(i in 1:n){
     kernel_estimates <- c(kernel_estimates,kernesti.regr(x[i,],x,y,h=h,...))
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
  mean(eval_func(fits-test_y,...)) 
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

cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",intercept=TRUE,criteria="CV",cvFunc="check",nfolds=10,foldid=NULL,nlambda=100,eps=.0001,init.lambda=1,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
# nfolds: number of folds for cross validation
# foldid: preset id of folds
  p_range <- 1:dim(x)[2] + intercept
  n <- dim(x)[1]
  pen_func <- switch(which(c("LASSO","SCAD","MCP")==penalty), lasso, scad, mcp)
  if(is.null(lambda)){
  # find a lambda that sets all coefficients to zero. 
  # Strategy is to get \lambda \sum p_\lambda(|\beta_j}) >= \sum \rho_\tau(y-quantile(y,tau)
  # do this by fitting model with lambda = init.lambda and then set new lambda such that 
  # \lambda* = \sum \rho_\tau(y-quantile(y,tau)) / \sum p_\lambda(|beta_j|) repeat as needed
     sample_q <- quantile(y,tau)
     inter_only_rho <- sum(check(y-sample_q,tau))
     lambda_star <- init.lambda
     searching <- TRUE
     while(searching){
       if(penalty=="LASSO"){
         init_fit <- rq.lasso.fit(x,y,tau,lambda=lambda_star,weights,intercept,...)
       } else{
         init_fit <- rq.nc.fit(x,y,tau,lambda=lambda_star,weights,intercept,...)
       }
       if(sum(init_fit$coefficients[p_range])==0){
         searching <- FALSE     
       } else{
         lambda_star = inter_only_rho / sum(sapply(init_fit$coefficients[p_range],pen_func,1)) 
         #1 used here because solving for lambda
       }
     }
     lambda_min <- eps*lambda_star
     lambda <- exp(seq(log(lambda_min),log(lambda_star),length.out=nlambda))
  }
  if(penalty=="LASSO"){
     models <- lapply(lambda,rq.lasso.fit, x=x,y=y,tau=tau,weights=weights,intercept=intercept,...)
  } else{
     models <- lapply(lambda,rq.nc.fit, x=x,y=y,tau=tau,weights=weights,intercept=intercept,penalty=penalty,...)
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
         cv_models <- lapply(lambda,rq.lasso.fit, x=x,y=y,tau=tau,weights=weights,intercept=intercept,...)
      } else{
         cv_models <- lapply(lambda,rq.nc.fit, x=x,y=y,tau=tau,weights=weights,intercept=intercept,penalty=penalty,...)
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



rq.nc.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,penalty="SCAD",a=3.7,iterations=10,converge_criteria=.0001,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# lambda takes values of 1 or p
# penalty SCAD or MCP
   if(penalty=="SCAD"){
     deriv_func <- scad_deriv
   }
   if(penalty=="MCP"){
     deriv_func <- mcp_deriv
   }
   if(is.null(dim(x))){                                                                                    
      stop('x needs to be a matrix with more than 1 column')
   }
   p <- dim(x)[2]
   n <- dim(x)[1]
   if(n != length(y)){
      stop('length of y and rows of x do not match')
   }
   if(is.null(lambda)==TRUE | (length(lambda) != 1 & length(lambda) != dim(x)[2])){
      stop(paste('input of lambda must be of length 1 or', dim(x)[2]))
   }
   if( sum(lambda <= 0) > 0){
      stop('lambda must be positive')
   }
   if(length(lambda) != 1){
      pen_vars <- 1:p[lambda != 0]
      pen_range <- intercept + pen_vars
   } else{
      pen_range <- intercept + 1:p
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
          lambda_update[!pen_range] <- 0
          lambda_update[pen_range] <- mapply(deriv_func, beta=abs(sub_fit$coefficients[pen_range]),
                                                         lambda=lambda,
                                                         MoreArgs=list(a=a))
      }
      #lambda_update <- n*lambda_update
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
    foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="QICD", 
    ...) 
{
    p_range <- 1:dim(x)[2] + intercept
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
                  intercept, penalty,alg, ...)
                search_num <- 2
            }
            else {
                init_fit <- rq.group.fit(x, y, groups, tau, lambda_star, 
                  intercept, penalty,alg, initial_beta = beta_update, 
                  ...)
            }
            beta_update <- init_fit$coefficients
            if (sum(init_fit$coefficients[p_range]) == 0) {
                searching <- FALSE
            }
            else {
                 lambda_star = max((inter_only_rho-init_fit$rho)/ sum(sapply(init_fit$coefficients[p_range],pen_func, 1)),
                                    lambda_star+1)
                 #this is sort of a hack, need to think of a better way to pick lambda_star
                #lambda_star = inter_only_rho/sum(sapply(init_fit$coefficients[p_range], 
                #  pen_func, 1))
                #weird behavior if we don't set penalty lambda to 1
                #in general this idea needs to be improved upon
            }
       }
          
       lambda_min <- eps * lambda_star
       lambda <- exp(seq(log(lambda_min), log(lambda_star), length.out = nlambda))
       fine_tune <- TRUE
       fine_tune_pos <- length(lambda)-1
       while(fine_tune){
         test_fit <- rq.group.fit(x,y,groups,tau,lambda[fine_tune_pos],intercept,penalty,alg,...)
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
         lambda <- exp(seq(log(lambda_min), log(lambda_star), length.out = nlambda))
       }
    }
    
    models <- groupMultLambda(x = x, y = y, groups = groups, 
        tau = tau, lambda = lambda, intercept = intercept, penalty=penalty,alg=alg, ...)
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
                intercept = intercept,penalty=penalty,alg=alg, ...)
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
                penalty = "LASSO",alg="QICD", ...) 
{
    p <- dim(x)[2]
    if (!penalty %in% c("LASSO", "SCAD", "MCP")) {
        stop("Penalty must be LASSO, SCAD or MCP")
    }
    if (alg == "QICD") {
        coefs <- groupQICD(x = x, y = y, groups = groups, lambda = lambda, 
            tau = tau, intercept = intercept, penalty = penalty, 
            ...)
        return_val <- list()
        return_val$coefficients <- coefs
        if (is.null(colnames(x))) {
            x_names <- paste("x", 1:p, sep = "")
        }
        else {
            x_names <- colnames(x)
        }
        if (intercept) {
            x_names <- c("intercept", x_names)
        }
        attributes(return_val$coefficients)$names <- x_names
        if (intercept) {
            x <- cbind(1, x)
        }
        return_val$residuals <- y - x %*% coefs
        return_val$rho <- sum(check(return_val$residuals))
        return_val$tau <- tau
        return_val$n <- dim(x)[1]
        return_val$intercept <- intercept
        class(return_val) <- c("rq.group.pen", "rq.pen", "rqNC")
    }
    else {
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
        new_lambda <- NULL
        group_count <- xtabs(~groups)
        for (g in 1:group_num) {
             new_lambda <- c(new_lambda, rep(lambda[g], each = group_count[g]))
        }
          
        if (penalty == "LASSO") {
            return_val <- rq.lasso.fit(x, y, tau, new_lambda, 
                intercept = intercept, ...)
            class(return_val) <- c("rq.group.pen", "rq.pen", 
                "rqLASSO")
        }
        else {
            return_val <- rq.group.lin.prog(x,y,groups,tau,lambda,intercept=intercept,penalty=penalty,...)
            class(return_val) <- c("rq.group.pen", "rq.pen")
        }
    }
    return_val
}

plot.cv.rq.group.pen <- function (x,y=NULL,...) 
{
    plot(x$cv[, 1], x$cv[, 2])
}




