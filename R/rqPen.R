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
  min(x,0)
}

lasso <- function(x,lambda=1){
   lambda*abs(x)
}

scad <- function(x, lambda=1, a=3.7){
  if(abs(x) < lambda){
    lambda*abs(x)
  }
  else if( abs(x) >= a*lambda){
  #third case but easier to program
    (a+1)*lambda^2 / 2
  }
  else{
    ( (a^2-1)*lambda^2 - (abs(x)-a*lambda)^2) / ( 2*(a-1)) 
  }
}

scad_deriv <- function(x, lambda=1,a=3.7){
  if(abs(x) <= lambda){
  #really undefined but should be penalized as lambda using the taylor expansion technique
    return_val <- lambda
  } 
  if(lambda == 0){
    return_val <- 0
  }
  if(abs(x) > lambda){
     return_val <- max(a*lambda-abs(x),0)/(a-1)
  }
  if(return_val == 0){
    0
  } else if(x == 0){
    lambda
  }else{
    sign(x)*return_val
  }
}

#scad_1_deriv <- function(x,lambda=1,a=3.7){
#  sign(x)*lambda
#}
#
#scad_2_deriv <- function(x,lambda=1,a=3.7){
#  sign(x)*lambda*(1-( pos_part(a*lambda-abs(x)) / ( lambda*(a-1))))*(abs(x) > lambda)
#}

mcp <- function(x, lambda=1, a=3){
    if(abs(x) < a*lambda){
      lambda*(abs(x)-x^2/(2*a*lambda))
    } else{
      a*lambda^2 / 2
    }
}

mcp_deriv <- function(x, lambda=1, a=3){
  if(x == 0){
     lambda
  } else if(abs(x) < a*lambda){
    lambda*sign(x)- x*(1/a)
  } else{
    0
  }
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

rq.lasso.fit.mult <- function(x,y,tau_seq=c(.1,.3,.5,.7,.9),lambda=NULL,weights=NULL,intercept=TRUE,coef.cutoff=.00000001,...){
   model_list <- list()
   iter <- 1
   for(tau in tau_seq){
      model_list[[iter]] <- rq.lasso.fit(x,y,tau,lambda,weights,intercept,coef.cutoff,...)
      iter <- iter+1
   }
   model_list
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



rq.lasso.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,coef.cutoff=.00000001,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# lambda takes values of 1 or p
# coef.cutoff is a threshold to set to zero. 
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
   if( sum(lambda < 0) > 0){
      stop(paste('lambda must be positive and we have a lambda of ', lambda, sep=""))
   }
   lambda <- lambda*n # need this to account for the fact that rq does not normalize the objective function
   if(length(lambda)==1){
      pen_x <- rbind(diag(rep(lambda,p)),diag(rep(-lambda,p)))
   } else{
      pen_x <- rbind(diag(lambda), diag(-lambda))
      pen_x <- pen_x[rowSums(pen_x==0)!=dim(pen_x)[2],]
   }
   aug_n <- dim(pen_x)[1]
   aug_x <- rbind(x,pen_x)
   if(intercept){
      aug_x <- cbind(c(rep(1,n),rep(0,aug_n)), aug_x)
   }
   aug_y <- c(y, rep(0,aug_n))
   if(is.null(weights)){
     model <- rq(aug_y ~ aug_x+0, tau=tau, ...)
   } else{
     if(length(weights) != n){
       stop("Length of weights does not match length of y")
     }
     weights <- c(weights, rep(1,aug_n))
     model <- rq(aug_y ~ aug_x+0, tau=tau, weights=weights,...)
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
   return_val$PenRho <- model$rho
   return_val$residuals <- model$residuals[1:n]      
   return_val$rho <- sum(sapply(return_val$residuals,check,tau))
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
      sub_fit <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda_update,weights=weights,intercept=intercept)
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

coef.cv.rq.group.pen <- function(object, lambda='min',...){
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  object$beta[,target_model]
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
                                initial_beta=NULL,iterations=10,converge_criteria=.0001,...){
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


groupQICD <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
    max_iter = 100, eps = 1e-05, penalty = "SCAD", a = 3.7, coef.cutoff = 1e-08, 
    initial_beta = NULL) 
{
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
    if (intercept) {
        x <- cbind(1, x)
        groups <- c(0, groups)
    }
    if (penalty == "SCAD") {
        deriv_func <- scad_deriv
    }
    if (penalty == "MCP") {
        deriv_func <- mcp_deriv
    }
    if (is.null(initial_beta)) {
        new_outer_parameters <- rep(0, dim(x)[2])
    }
    else {
        new_outer_parameters <- initial_beta
    }
    n <- dim(x)[1]
    lambda <- n*lambda #again need to normalize for fact that rq minimizes unnormalized penalty function
    outer_loop_done <- FALSE
    outer_count <- 0
    while (!outer_loop_done) {
        inner_loop_done <- FALSE
        inner_count <- 0
        old_outer_parameters <- new_outer_parameters
        new_inner_parameters <- old_outer_parameters
        if (penalty == "LASSO") {
            loop_lambdas <- lambda
        }
        if (penalty == "SCAD" | penalty == "MCP") {
            loop_lambdas <- NULL
            for (g in unique(groups)) {
                if (g != 0) {
                  if (lambda[g] != 0) {
                    active_vars <- which(groups == g)
                    loop_lambdas <- c(loop_lambdas, deriv_func(sum(abs(new_outer_parameters[active_vars])), 
                      lambda[g], a))
                  }
                }
            }
        }
        while (!inner_loop_done) {
            old_inner_parameters <- new_inner_parameters
            for (g in unique(groups)) {
                active_vars <- which(groups == g)
                y_star <- y - x[, -active_vars] %*% new_inner_parameters[-active_vars]
                if (g == 0 || loop_lambdas[g] == 0) {
                  q1 <- rq(y_star ~ x[, active_vars] + 0, tau = tau)
                  new_inner_parameters[active_vars] <- q1$coefficients
                }
                else {
                  y_aug <- c(y_star, rep(0, 2 * length(active_vars)))
                  x_aug <- x
                  lam_val <- loop_lambdas[g]
                  for (spots in active_vars) {
                    row_1_aug <- row_2_aug <- rep(0, dim(x)[2])
                    row_1_aug[spots] <- lam_val
                    row_2_aug[spots] <- -lam_val
                    x_aug <- rbind(x_aug, row_1_aug, row_2_aug)
                  }
                  q1 <- rq(y_aug ~ x_aug[, active_vars] + 0, 
                    tau = tau)
                  new_inner_parameters[active_vars] <- q1$coefficients
                }
            }
            if (sum((new_inner_parameters - old_inner_parameters)^2) < 
                eps | inner_count == max_iter) {
                inner_loop_done <- TRUE
            }
            else {
                inner_count <- inner_count + 1
            }
        }
        new_outer_parameters <- new_inner_parameters
        if (sum((new_outer_parameters - old_outer_parameters)^2) < 
            eps | outer_count == max_iter) {
            outer_loop_done <- TRUE
        }
        else {
            outer_count <- outer_count + 1
        }
    }
    new_outer_parameters[abs(new_outer_parameters) < coef.cutoff] <- 0
    new_outer_parameters
}

groupMultLambda <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, penalty="LASSO", 
    #initial_beta = NULL,
    alg="QICD", ...) 
{
    return_val <- list()
    #if (is.null(initial_beta)) {
    #    initial_beta <- rep(0, dim(x)[2] + intercept)
    #}
    pos <- 1
    for (lam in lambda) {
        return_val[[pos]] <- rq.group.fit(x = x, y = y, groups = groups, 
            tau = tau, lambda = lam, intercept = intercept, penalty=penalty,alg=alg, 
            ...)
        #initial_beta <- return_val[[pos]]$coefficients
        pos <- pos + 1
    }
    return_val
}

nonzero <- function (obj) 
{
    UseMethod("nonzero")
}

plot.cv.rq.group.pen <- function (x,y=NULL,...) 
{
    plot(x$cv[, 1], x$cv[, 2])
}

nonzero.cv.rq.group.pen <- function (obj) 
{
    coefs <- coefficients(obj)
    if (obj$intercept) {
        coefs <- coefs[-1]
    }
    tapply(coefs, obj$groups, sum) != 0
}









