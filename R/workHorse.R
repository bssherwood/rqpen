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
  absx <- u <- abs(x)
  u[] <- 0
  index <- absx < a*lambda & absx > 0
  u[ index ] <-
       ifelse( absx[ index ] <= lambda, 
               lambda,
               ( a*lambda - absx[ index ] )/( a-1 ) 
             )
  u[index] <- u[index]*sign( x[index] )
  u[ x == 0 ] <- lambda # because we take derivative as x approaces 0 from above

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
  u <- x
  u[] <- 0
  index <- abs(x) < a*lambda
  u[ index ] <- ifelse( x[index] == 0,
                        lambda,
                        lambda*sign(x[index]) - x[index]/a
                      )

  u
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
   model_list <- list()
   iter <- 1
   for(tau in tau_seq){
      model_list[[iter]] <- rq.lasso.fit(x,y,tau,lambda,weights,intercept,coef.cutoff,...)
      iter <- iter+1
   }
   model_list
}

rq.lasso.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                         coef.cutoff=1e-08, method="br",penVars=NULL, ...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# lambda takes values of 1 or p
# coef.cutoff is a threshold to set to zero. 
# Choose the method used to estimate the coefficients ("br" or "fn")
### According to quantreg manual and my experience, "fn" is much faster for big n
### The "n" can grow rapidly using lin. prog. approach  
# penVars - variables to be penalized, doesn't work if lambda has multiple entries

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

   # ##############################################################################################
   # ### This uses the LASSO.fit or LASSO.fit.nonpen functions to obtain coefficient estimates. ###
   # ### Note that the coefficients might need to be reordered to match x.                      ###
   # ##############################################################################################
   # if( !is.null(weights) & length(weights) != n )
   #    stop("Length of weights does not match length of y")

   # ### Find indices of penalized and nonpenalized coefficients
   # nonpenVars <-     ### indices of nonpenalized coefficients (lambdas of 0 or not included in penVars), NULL if no nonpenalized oefficients
   # penVars    <-     ### indices of penalized coefficients

   # xnew <- as.matrix( x[,penVars] )
   # if( is.null(nonpenVars) ){
   #   coefs <- LASSO.fit(y, xnew, tau, intercept, coef.cutoff, weights)
   # } else {
   #   znew <- as.matrix( x[,nonpenVars] )
   #   coefs <- LASSO.fit.nonpen(y, xnew, znew, tau, intercept, coef.cutoff, weights)
   #   coefs[ intercept + 1:p ] <- coefs[ intercept + order(c(penVars, nonpenVars)) ]
   # }

   # ### coefs are the coefficients in the correct order with the intercept first
   # ##############################################################################################
   # ##############################################################################################
   # ##############################################################################################

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
   return_val$PenRho <- model$rho
   return_val$residuals <- model$residuals[1:n]   
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
       initial_beta <- rq.lasso.fit(x,y,tau,new_lambda, intercept=intercept, coef.cutoff=coef.cutoff,method="br",...)$coefficients
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