#' Calculate information criterion for penalized quantile regression models
#'
#' @param model model from a rq.pen.seq() object
#' @param n Sample size
#' @param method Choice of BIC, AIC or PBIC, a large p BIC. 
#'
#' @return 
#' Let \eqn{\hat{\beta}} be the coefficient vectors for the estimated model. Function returns the value 
#' \deqn{\log(\sum_{i=1}^n \rho_\tau(y_i-x_i^\top\hat{\beta})) + d*b/(2n),} where d is the number of nonzero coefficients and b depends on the method used. For AIC \eqn{b=2},
#' for BIC \eqn{b=log(n)} and for PBIC \eqn{d=log(n)*log(p)} where p is the dimension of \eqn{\hat{\beta}}. Returns this value for each coefficient vector in the model, so one
#' for every value of \eqn{\lambda}. 
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(runif(800),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' m1 <- rq.pen(x,y,tau=c(.25,.75))
#' # returns the IC values for tau=.25
#' qic(m1$models[[1]],m1$n) 
#' # returns the IC values for tau=.75
#' qic(m1$models[[2]],m1$n) 
#' @references 
#' \insertRef{qrbic}{rqPen}
#'@author Ben Sherwood, \email{ben.sherwood@ku.edu}
qic <- function(model,n, method=c("BIC","AIC","PBIC")){
  method <- match.arg(method)
	tau <- model$debug
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



#' Selects tuning parameter \eqn{\lambda} and a according to information criterion of choice. For a given \eqn{\hat{\beta}} the information criterion is calculated
#' as
#' \deqn{\log(\sum_{i=1}^n \rho_\tau(y_i-x_i^\top\hat{\beta})) + d*b/(2n),} where d is the number of nonzero coefficients and b depends on the method used. For AIC \eqn{b=2},
#' for BIC \eqn{b=log(n)} and for PBIC \eqn{d=log(n)*log(p)} where p is the dimension of \eqn{\hat{\beta}}.
#' If septau set to FALSE then calculations are made across the quantiles. Let \eqn{\hat{\beta}^q} be the coefficient vector for the qth quantile of Q quantiles. In addition let \eqn{d_q} and \eqn{b_q} 
#' be d and b values from the qth quantile model. Note, for all of these we are assuming eqn and a are the same. Then the summary across all quantiles is 
#' \deqn{\sum_{q=1}^Q w_q[ \log(\sum_{i=1}^n  \rho_\tau(y_i-x_i^\top\hat{\beta}^q)) + d_q*b_q/(2n)],}
#' where \eqn{w_q} is the weight assigned for the qth quantile model. 
#'
#' @param obj A rq.pen.seq or rq.pen.seq.cv object. 
#' @param method Choice of BIC, AIC or PBIC, a large p BIC.
#' @param septau If optimal values of \eqn{\lambda} and a can vary with \eqn{\tau}. Default is TRUE. 
#' @param weights Weights for each quantile. Useful if you set septau to FALSE but want different weights for the different quantiles. If not specified default is to have \eqn{w_q=1} for all quantiles.
#'
#' @return 
#' \itemize{
#' \item{coefficients}{Coefficients of the selected models.}
#' \item{ic}{Information criterion values for all considered models.}
#' \item{modelsInfo}{Model info for the selected models related to the original object obj.}
#' \item{gic}{Information criterion summarized across all quantiles. Only returned if septau set to FALSE}
#' }
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(runif(800),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' m1 <- rq.pen(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.75))
#' qic.select(m1)
#' @references 
#' \insertRef{qrbic}{rqPen}
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
qic.select <- function(obj, method=c("BIC","AIC","PBIC"),septau=TRUE,weights=NULL){
# code help: Maybe think about how the qic values are returned for the septau=TRUE case. Also, potential issue with different values of lambda
	method <- match.arg(method)
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
		modelsInfo <- data.table(obj$modelsInfo,minQIC,lambdaIndex,lambda=obj$lambda[lambdaIndex])
		#modelsInfo <- data.table(modelsInfo)
		tau <- modelsInfo$tau
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
		a <- obj$a
		modelsInfo <- subset(obj$modelsInfo, a==returnA)
		modelsInfo <- cbind(modelsInfo,minQIC=tqic_vals[modelsInfo$modelIndex,minIndex[2]],lambdaIndex=minIndex[2],lambda=obj$lambda[minIndex[2]])
		coefs <- coef(obj, lambdaIndex=minIndex[2], modelsIndex=modelsInfo$modelIndex)
	}
	#coefIndex <- 1:nt
	#modelsInfo <- cbind(modelsInfo, coefIndex)
	coefs <- do.call(cbind,coefs)
	colnames(coefs) <- paste0("tau=",obj$tau)
	
	return_val <- list(coefficients = coefs, ic=qic_vals,modelsInfo=modelsInfo, gic=gic)
	class(return_val) <- "qic.select"
	return_val
}


#' Print a qic.select object
#'
#' @param x qic.select object
#' @param ... optional arguments
#'
#' @return Prints the coefficients of the qic.select object
#' @export
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
print.qic.select <- function(x,...){
   print(coefficients(x))
}

#' Predictions from a qic.select object
#'
#' @param object qic.select object
#' @param newdata Data matrix to make predictions from. 
#' @param ... optional arguments
#'
#' @return A matrix of predicted values.
#' @export
#'
#' @examples
#' x <- matrix(runif(800),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' m1 <- rq.pen(x,y,tau=c(.25,.75))
#' q1 <- qic.select(m1)
#' newx <- matrix(runif(80),ncol=8)
#' preds <- predict(q1,newx)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
predict.qic.select <- function(object, newdata, ...){
	#coefs <- do.call(cbind,coefficients(object))
	cbind(1,newdata) %*% coefficients(object)
}


#' Print a rq.pen.seq object
#'
#' @param x rq.pen.seq object
#' @param ... optional arguments
#'
#' @return If only one model, prints a data.frame of the number of nonzero coefficients and lambda. Otherwise prints information about the quantiles being modeled and choices for a.
#' @export
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
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


#' Returns Coefficients of a cv.rq.pen object
#' 
#' Warning: this function will be deprecated and not exported in future versions of rqPen, due to the switch from cv.rq.pen() to rq.pen.cv().
#'
#' @param object cv.rq.pen object 
#' @param lambda Value of lambda, default is to use the minimum value. 
#' @param ... Additional parameters.
#'
#' @return Coefficients for a given lambda, or the lambda associated with the minimum cv value. 
#' @export
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
coef.cv.rq.pen <- function(object, lambda='min',...){
  deprecate_soft("3.0","coef.cv.rq.pen()","coef.rq.pen.seq.cv()")
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  coefficients(object$models[[target_model]])
}


#' Fit a quantile regression model using a penalized quantile loss function. 
#'
#' @param x matrix of predictors
#' @param y vector of responses
#' @param tau vector of quantiles
#' @param lambda vector of lambda, if not set will be generated automatically
#' @param penalty choice of penalty
#' @param a Additional tuning parameter, not used for lasso or ridge penalties. However, will be set to the elastic net values of 1 and 0 respectively. Defaults are ENet(0), aLASSO(1), SCAD(3.7) and MCP(3).
#' @param nlambda number of lambda, ignored if lambda is set
#' @param eps If not pre-specified the lambda vector will be from lambda_max to lambda_max times eps
#' @param penalty.factor penalty factor for the predictors
#' @param alg Algorithm used.
#' @param scalex Whether x should be scaled before fitting the model. Coefficients are returned on the original scale. 
#' @param tau.penalty.factor A penalty factor for each quantile.
#' @param coef.cutoff Some of the linear programs will provide very small, but not sparse solutions. Estimates below this number will be set to zero. This is ignored if a non-linear programming algorithm is used. 
#' @param max.iter Maximum number of iterations of non-linear programming algorithms.
#' @param converge.eps Convergence threshold for non-linear programming algorithms. 
#' @param lambda.discard Algorithm may stop for small values of lambda if the coefficient estimates are not changing drastically. One example of this is it is possible for the LLA weights of the non-convex functions to all become zero and smaller values of lambda are extremely likely to produce the same zero weights. 
#' @param ... Extra parameters. 
#' 
#' @description  
#' Let q index the Q quantiles of interest. Let \eqn{\rho_\tau(a) = a[\tau-I(a<0)]}. Fits quantile regression models by minimizing the penalized objective function of
#' \deqn{\frac{1}{n} \sum_{q=1}^Q \sum_{i=1}^n \rho_\tau(y_i-x_i^\beta^q) + \sum_{q=1}^Q  \sum_{j=1}^p P(\beta^q_p,w_q*v_j*\lambda,a).}
#' Where \eqn{w_q} and \eqn{v_j} are designated by penalty.factor and tau.penalty.factor respectively. Value of P() depends on the penalty. Briefly, but see references or vignette for more details,
#' \itemize{
#' \item{LASSO:}{ \eqn{P(\beta,\lambda,a)=\lambda|\beta|}}
#' \item{SCAD:}{ \eqn{P(\beta,\lambda,a)=SCAD(\beta,\lambda,a)}}
#' \item{MCP:}{ \eqn{P(\beta,\lambda,a)=MCP(\beta,\lambda,a)}}
#' \item{Ridge:}{ \eqn{P(\beta,\lambda,a)=\lambda\beta^2}}
#' \item{Elastic Net:}{ \eqn{P(\beta,\lambda,a)=a*\lambda|\beta|+(1-a)*\lambda*\beta^2}}
#' \item{Adaptive LASSO:}{ \eqn{P(\beta,\lambda,a)=\frac{\lambda |\beta|}{|\beta_0|^a}}}
#' }
#' For Adaptive LASSO the values of \eqn{\beta_0} come from a Ridge solution with the same value of \eqn{\lambda}. 
#' @return An rq.pen.seq object. 
#' \itemize{
#' \item{models: }{ A list of each model fit for each tau and a combination.}
#' \item{n:}{ Sample size.}
#' \item{p:}{ Number of predictors.}
#' \item{alg:}{ Algorithm used. Options are "huber", "qicd" or any method implemented in rq(), such as "br". }
#' \item{tau:}{ Quantiles modeled.}
#' \item{a:}{ Tuning parameters a used.}
#' \item{modelsInfo:}{ Information about the quantile and a value for each model.}
#' \item{lambda:}{ Lambda values used for all models. If a model has fewer coefficients than lambda, say k. Then it used the first k values of lambda. Setting lambda.discard to TRUE will gurantee all values use the same lambdas, but may increase computational time noticeably and for little gain.}
#' \item{penalty:}{ Penalty used.}
#' \item{call:}{ Original call.}
#' }
#' Each model in the models list has the following values. 
#' \itemize{
#' \item{coefficients:}{ Coefficients for each value of lambda.}
#' \item{rho:}{ The unpenalized objective function for each value of lambda.}
#' \item{PenRho:}{ The penalized objective function for each value of lambda.}
#' \item{nzero:}{ The number of nonzero coefficients for each value of lambda.}
#' \item{tau:}{ Quantile of the model.}
#' \item{a:}{ Value of a for the penalized loss function.}
#' }
#' 
#' If the Huber algorithm is used than \eqn{\rho_\tau(y_i-x_i^\top\beta)} is replaced by a Huber-type approximation. Specifically, it is replaced by \eqn{h^\tau_\gamma(y_i-x_i^\top\beta)/2} where 
#' \deqn{h^\tau_\gamma(a) = a^2/(2\gamma)I(|a| \leq \gamma) + (|a|-\gamma/2)I(|a|>\gamma)+(2\tau-1)a.}
#' Where if \eqn{\tau=.5}, we get the usual Huber loss function. The Huber implementation calls the package hqreg which implements the methods of Yi and Huang (2017) 
#' for Huber loss with elastic net penalties. For non-elastic net penalties the LLA algorithm of Zou and Li (2008) is used to approximate those loss functions
#' with a lasso penalty with different weights for each predictor. 
#' @export
#'
#' @examples
#' n <- 200
#' p <- 8
#' x <- matrix(runif(n*p),ncol=p)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' r1 <- rq.pen(x,y) #Lasso fit for median
#' # Lasso for multiple quantiles
#' r2 <- rq.pen(x,y,tau=c(.25,.5,.75))
#' # Elastic net fit for multiple quantiles, which must use Huber algorithm
#' r3 <- rq.pen(x,y,penalty="ENet",a=c(0,.5,1),alg="huber")
#' # First variable is not penalized
#' r4 <- rq.pen(x,y,penalty.factor=c(0,rep(1,7)))
#' tvals <- c(.1,.2,.3,.4,.5)
#' #Similar to penalty proposed by Belloni and Chernouzhukov. 
#' #To be exact you would divide the tau.penalty.factor by n. 
#' r5 <- rq.pen(x,y,tau=tvals, tau.penalty.factor=sqrt(tvals*(1-tvals)))
#' @references 
#' \insertRef{lla}{rqPen}
#' 
#' \insertRef{huber_cd}{rqPen}
#' 
#' \insertRef{qr_lasso}{rqPen}
#' 
#' \insertRef{qr_cd}{rqPen}
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} and Adam Maidman
rq.pen <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.05,.01), 
	penalty.factor = rep(1, ncol(x)),alg=ifelse(sum(dim(x))<200 & penalty %in% c("LASSO","aLASSO","SCAD","MCP"),"br","huber"),scalex=TRUE,tau.penalty.factor=rep(1,length(tau)),
	coef.cutoff=1e-8,max.iter=10000,converge.eps=1e-7,lambda.discard=TRUE,...){
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
	if(alg=="lasso" || alg=="scad"){
	  stop("Choice of lasso or scad for algorithm is invalid, use ``huber'' or a non-lasso, non-scad method from quantreg::rq()")
	}
	if(penalty=="LASSO"){
		fit <- rq.lasso(x,y,tau,lambda,nlambda,eps,penalty.factor,alg,scalex=FALSE,tau.penalty.factor,coef.cutoff,max.iter,converge.eps,lambda.discard=lambda.discard,...)
	} else if(penalty=="Ridge"){
		if(alg != "huber"){
			stop("huber alg is only option for Ridge penalty")
		}
		fit <-  rq.enet(x,y,tau,lambda,nlambda,eps, penalty.factor,scalex=FALSE,tau.penalty.factor,a=0,max.iter,converge.eps,lambda.discard=lambda.discard,...)
	} else if(penalty == "ENet"){
		if(alg != "huber"){
			stop("huber alg is only option for ENet penalty")
		}
		if(is.null(a)){
			stop("Specify a value for a for ENet penalty")
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
	if(lambda.discard){
		lmax <- max(sapply(fit$models,lambdanum))
		fit$lambda <- fit$lambda[1:lmax]
	}
	fit$call <- match.call()
	fit
}

#' Returns coefficients of a rq.pen.seq object
#'
#' @param object rq.pen.seq object
#' @param tau Quantile of interest. Default is NULL, which will return all quantiles. Should not be specified if modelsIndex is used.  
#' @param a Tuning parameter of a. Default is NULL, which returns coefficients for all values of a. Should not be specified if modelsIndex is used. 
#' @param lambda Tuning parameter of \eqn{\lambda}. Default is NULL, which returns coefficients for all values of \eqn{\lambda}.
#' @param modelsIndex Index of the models for which coefficients should be returned. Does not need to be specified if tau or a are specified. 
#' @param lambdaIndex Index of the lambda values for which coefficients should be returned. Does not need to be specified if lambda is specified. 
#' @param ... Additional parameters. 
#'
#' @return A list of a matrix of coefficients for each tau and a combination
#' @export
#'
#' @examples
#' x <- matrix(runif(800),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' m1 <- rq.pen(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.75),lambda=c(.1,.05,.01))
#' allCoefs <- coef(m1)
#' targetCoefs <- coef(m1,tau=.25,a=.5,lambda=.1)
#' idxApproach <- coef(m1,modelsIndex=2)
#' bothIdxApproach <- coef(m1,modelsIndex=2,lambdaIndex=1)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
coef.rq.pen.seq <- function(object,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL,...){
  models <- getModels(object,tau,a,lambda,modelsIndex,lambdaIndex)
  lapply(models$targetModels,getModelCoefs,models$lambdaIndex)
}


#' Predictions from rq.pen.seq.cv object
#'
#' @param object rq.pen.seq.cv object
#' @param newx Matrix of predictors 
#' @param tau Quantile of interest. Default is NULL, which will return all quantiles. Should not be specified if modelsIndex is used.  
#' @param septau Whether tuning parameter should be optimized separately for each quantile. 
#' @param cvmin If TRUE then minimum error is used, if FALSE then one standard error rule is used. 
#' @param useDefaults Whether the default results are used. Set to FALSE if you you want to specify specific models and lambda values. 
#' @param ... Additional parameters sent to coef.rq.pen.seq.cv(). 
#'
#' @return A list of a matrix of predictions for each tau and a combination
#' @export
#'
#' @examples
#' x <- matrix(runif(1600),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(200)
#' m1 <- rq.pen.cv(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.75),lambda=c(.1,.05,.01))
#' newx <- matrix(runif(80),ncol=8)
#' cvpreds <- predict(m1,newx)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
predict.rq.pen.seq.cv <- function(object, newx,tau=NULL,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,...){
  coefs <- coefficients(object,septau=septau,cvmin=cvmin,useDefaults=useDefaults,tau=tau,...)
	lapply(coefs, quick.predict,newx=newx)
}

#' Predictions from rq.pen.seq.cv object
#'
#' @param object rq.pen.seq.cv object
#' @param newx Matrix of predictors 
#' @param tau Quantile of interest. Default is NULL, which will return all quantiles. Should not be specified if modelsIndex is used.  
#' @param a Tuning parameter of a. Default is NULL, which returns coefficients for all values of a. Should not be specified if modelsIndex is used. 
#' @param lambda Tuning parameter of \eqn{\lambda}. Default is NULL, which returns coefficients for all values of \eqn{\lambda}.
#' @param modelsIndex Index of the models for which coefficients should be returned. Does not need to be specified if tau or a are specified. 
#' @param lambdaIndex Index of the lambda values for which coefficients should be returned. Does not need to be specified if lambda is specified. 
#' @param ... Additional parameters. 
#'
#' @return A list of a matrix of predictions for each tau and a combination
#' @export
#'
#' @examples
#' x <- matrix(runif(800),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' m1 <- rq.pen(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.75),lambda=c(.1,.05,.01))
#' newx <- matrix(runif(80),ncol=8)
#' allCoefs <- predict(m1,newx)
#' targetCoefs <- predict(m1,newx,tau=.25,a=.5,lambda=.1)
#' idxApproach <- predict(m1,newx,modelsIndex=2)
#' bothIdxApproach <- predict(m1,newx,modelsIndex=2,lambdaIndex=1)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
predict.rq.pen.seq <- function(object, newx,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL,...){
  coefs <- coefficients(object,tau,a,lambda,modelsIndex,lambdaIndex)
  lapply(coefs, quick.predict,newx=newx)
}

#' Does k-folds cross validation for rq.pen. If multiple values of a are specified then does a grid based search for best value of \eqn{\lambda} and a.
#'
#' @param x Matrix of predictors.
#' @param y Vector of responses.
#' @param tau Quantiles to be modeled. 
#' @param lambda Values of \eqn{\lambda}. Default will automatically select the \eqn{\lambda} values. 
#' @param penalty Choice of penalty between LASSO, Ridge, Elastic Net (ENet), Adaptive Lasso (aLASSO), SCAD and MCP.
#' @param a Tuning parameter of a. LASSO and Ridge has no second tuning parameter, but for notation is set to 1 or 0 respectively, the values for elastic net. Defaults are Ridge ()
#' @param cvFunc Loss function for cross-validation. Defaults to quantile loss, but user can specify their own function.
#' @param nfolds Number of folds.
#' @param foldid Ids for folds. If set will override nfolds.
#' @param nlambda Number of lambda, ignored if lambda is set.
#' @param groupError If set to false then reported error is the sum of all errors, not the sum of error for each fold.  
#' @param cvSummary Function to summarize the errors across the folds, default is mean. User can specify another function, such as median.
#' @param tauWeights Weights for the different tau models. 
#' @param printProgress If set to TRUE prints which partition is being worked on. 
#' @param ... Additional arguments passed to rq.pen()
#' 
#' @details 
#' Two cross validation results are returned. One that considers the best combination of a and lambda for each quantile. The second considers the best combination of the tuning 
#' parameters for all quantiles. Let \eqn{y_{b,i}} and \eqn{x_{b,i}} index the observations in 
#' fold b. Let \eqn{\hat{\beta}_{\tau,a,\lambda}^{-b}} be the estimator for a given quantile and tuning parameters that did not use the bth fold. Let \eqn{n_b} be the number of observations in fold
#' b. Then the cross validation error for fold b is 
#' \deqn{\mbox{CV}(b,\tau) = \frac{1}{n_b} \sum_{i=1}^{n_b} \rho_\tau(y_{b,i}-x_{b,i}^\top\hat{\beta}_{\tau,a,\lambda}^{-b}).}
#' Note that \eqn{\rho_\tau()} can be replaced by a different function by setting the cvFunc parameter. The function returns two different cross-validation summaries. The first is btr, by tau results. 
#' It provides the values of \code{lambda} and \code{a} that minimize the average, or whatever function is used for \code{cvSummary}, of \eqn{\mbox{CV}(b)}. In addition it provides the 
#' sparsest solution that is within one standard error of the minimum results. 
#' 
#' The other approach is the group tau results, gtr. Consider the case of estimating Q quantiles of \eqn{\tau_1,\ldots,\tau_Q} It returns the values of \code{lambda} and \code{a} that minimizes the average, or again whatever function is used for \code{cvSummary}, of 
#' \deqn{\sum_{q=1}^Q\mbox{CV}(b,\tau_q).} If only one quantile is modeled then the gtr results can be ignored as they provide the same minimum solution as btr. I THINK WRITING THIS WAY GIVES 
#' ME AN IDEA ON HOW TO DO STANDARD ERROR FOR THIS SETTTING. 
#' 
#' @return
#' \itemize{
#' \item{cverr:}{ Matrix of cvSummary function, default is average, cross-validation error for each model, tau and a combination, and lambda.}
#' \item{cvse:}{ Matrix of the standard error of cverr foreach model, tau and a combination, and lambda.}
#' \item{fit:}{ The rq.pen.seq object fit to the full data.}
#' \item{btr:}{ A data.table of the values of a and lambda that are best as determined by the minimum cross validation error and the one standard error rule, which fixes a. In btr the values of lambda and a are selected seperately for each quantile.}
#' \item{gtr:}{ A data.table for the combination of a and lambda that minimize the cross validation error across all tau.}
#' \item{gcve:}{ Group, across all quantiles, cross-validation error results for each value of a and lambda.}
#' \item{call:}{ Original call to the function.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' x <- matrix(runif(800),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' r1 <- rq.pen.cv(x,y) #lasso fit for median
#' # Elastic net fit for multiple values of a and tau
#' r2 <- rq.pen.cv(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.5,.75)) 
#' #same as above but more weight given to median when calculating group cross validation error. 
#' r3 <- rq.pen.cv(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.5,.75),tauWeights=c(.25,.5,.25))
#' # uses median cross-validation error instead of mean.
#' r4 <- rq.pen.cv(x,y,cvSummary=median)  
#'#Cross-validation with no penalty on the first variable.
#' r5 <- rq.pen.cv(x,y,penalty.factor=c(1,rep(0,7)))
#' }
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
rq.pen.cv <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,nlambda=100,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),printProgress=FALSE,...){
	n <- length(y)
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
	
	min_nl <- min(sapply(fit$models,lambdanum))
	if(min_nl != nl){
		warning("Different models had a different number of lambdas. To avoid this set lambda.discard=FALSE. Results presented are for the shortest lambda sequence")
		#code improvement, bad for loop
		for(i in 1:length(fit$models)){
			fit$models[[i]] <- clearModels(fit$models[[i]],min_nl)
		}
		fit$lambda <- fit$lambda[1:min_nl]
		nl <- min_nl
	}		
	
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
  		trainModel <- rq.pen(train_x,train_y,tau,lambda=fit$lambda,penalty=penalty,a=fit$a,lambda.discard=FALSE,alg=fit$alg,...)
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

#' Prints a rq.pen.seq.cv object
#'
#' @param x A req.pen.seq.cv object. 
#' @param ... Additional arguments. 
#'
#' @return Print of btr and gtr from a rq.pen.seq.cv object. If only one quantile is modeled then only btr is returned. 
#' @export
#'
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


#' Returns coefficients from a rq.pen.seq.cv object. 
#'
#' @param object An rq.pen.seq.cv object.
#' @param septau Whether tuning parameter should be optimized separately for each quantile. 
#' @param cvmin If TRUE then minimum error is used, if FALSE then one standard error rule is used. 
#' @param useDefaults Whether the default results are used. Set to FALSE if you you want to specify specific models and lambda values. 
#' @param tau Quantiles of interest. 
#' @param ... Additional parameters sent to coef.rq.pen.seq()
#'
#' @return Returns coefficients
#' @export
#'
#' @examples
#'  set.seed(1)
#'  x <- matrix(rnorm(800),nrow=100)
#'  y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#'  lassoModels <- rq.pen.cv(x,y,tau=seq(.1,.9,.1))
#'  coefficients(lassoModels,septau=FALSE)
#'  coefficients(lassoModels,cvmin=FALSE)
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
coef.rq.pen.seq.cv <- function(object,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,tau=NULL,...){
  if(!useDefaults){
    coefficients(object$models,tau=tau,...)
  } else{
    if(is.null(tau)){
      tau <- object$fit$tau
    }
    if(septau){
      keepers <- which(closeEnough(tau,object$btr$tau))
      btr <- object$btr[keepers,]
      models <- object$fit$models[btr$modelsIndex]
      if(cvmin){
        lambdaIndex <- btr$lambdaIndex
      } else{
        lambdaIndex <- btr$lambda1seIndex
      }
      returnVal <- vector(mode="list", length=length(models))
      names(returnVal) <- names(models)
      for(i in 1:length(returnVal)){
        returnVal[[i]] <- coef(object$fit,modelsIndex=btr$modelsIndex[i],lambdaIndex=lambdaIndex[i])[[1]]
      }
      returnVal
    } else{
      if(!cvmin){
        stop("One standard error approach not implemented for group choice of tuning parameter")
      } else{
        keepers <- which(closeEnough(tau,object$gtr$tau))
        gtr <- object$gtr[keepers,]
        coef(object$fit,modelsIndex=gtr$modelsIndex,lambdaIndex=gtr$lambdaIndex[1])
      }
    }
  }
}

# 
# 
# coef.rq.pen.seq.cv <- function(x,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,tau=NULL,...){
# 	if(!useDefaults){
# 		coefficients(x$models,tau=tau,...)
# 	} else{
# 		if(is.null(tau)){
# 			tau <- m1$fit$tau
# 		}
# 		if(septau){
# 			keepers <- which(closeEnough(tau,x$btr$tau))
# 			btr <- x$btr[keepers,]
# 			models <- x$fit$models[btr$modelsIndex]
# 			if(cvmin){
# 				lambdaIndex <- btr$lambdaIndex
# 			} else{
# 				lambdaIndex <- btr$lambda1seIndex
# 			}
# 			returnVal <- vector(mode="list", length=length(models))
# 			names(returnVal) <- names(models)
# 			for(i in 1:length(returnVal)){
# 				returnVal[[i]] <- coef(x$fit,modelsIndex=btr$modelsIndex[i],lambdaIndex=lambdaIndex[i])[[1]]
# 			}
# 			returnVal
# 		} else{
# 			if(!cvmin){
# 				stop("One standard error approach not implemented for group choice of tuning parameter")
# 			} else{
# 				keepers <- which(closeEnough(tau,x$gtr$tau))
# 				gtr <- x$gtr[keepers,]#subset(x$gtr, closeEnough(tau,x$gtr$tau))
# 				coef(x$fit,modelsIndex=gtr$modelsIndex,lambdaIndex=gtr$lambdaIndex[1])
# 			}
# 		}
# 	}
# }

#' Prints a cv.rq.pen object.  
#'
#' @param x cv.rq.pen object
#' @param ... Optional arguments, not used. 
#' 
#' @details Warning this function is deprecated and will not be exported in future releases. 
#'
#' @return Prints cross validation or information criterion values by lambda. 
#' @export
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
print.cv.rq.pen <- function(x,...){
   deprecate_soft("3.0","print.cv.rq.pen()","print.cv.rq.pen.seq()")
   cat("\nCoefficients:\n")
   print(coefficients(x,...))
   cat("\nCross Validation (or BIC) Results\n")
   print(x$cv)
}


#' Prints an rq.pen object
#' 
#' @description Warning this function is deprecated and will not be exported in future releases. 
#'
#' @param x The rq.pen object
#' @param ... Additional parameters sent to function
#'
#' @return Prints the coefficients of the object.
#' @export
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
print.rq.pen <- function(x,...){
    cat("\nCoefficients:\n")
	print(coefficients(x,...))
}

#' Performs cross validation for a group penalty. 
#'#'
#' @param x Matrix of predictors. 
#' @param y Vector of responses. 
#' @param tau Vector of quantiles. 
#' @param groups Vector of group assignments for the predictors.
#' @param lambda Vector of lambda values, if set to NULL they will be generated automatically.
#' @param a Vector of the other tuning parameter values. 
#' @param cvFunc Function used for cross-validation error, default is quantile loss. 
#' @param nfolds Number of folds used for cross validation. 
#' @param foldid Fold assignments, if not set this will be randomly created. 
#' @param groupError If errors are to be reported as a group or as the average for each fold. 
#' @param cvSummary The 
#' @param tauWeights Weights for the tau penalty. 
#' @param printProgress If set to TRUE will print which fold the process is working on. 
#' @param ... Additional parameters that will be sent to rq.group.pen().
#'
#' @return
#' \itemize{
#' \item{cverr}{Matrix of cvSummary function, default is average, cross-validation error for each model, tau and a combination, and lambda.}
#' \item{cvse}{Matrix of the standard error of cverr foreach model, tau and a combination, and lambda.}
#' \item{fit}{The rq.pen.seq object fit to the full data.}
#' \item{btr}{A data.table of the values of a and lambda that are best as determined by the minimum cross validation error and the one standard error rule, which fixes a. In btr the values of lambda and a are selected seperately for each quantile.}
#' \item{gtr}{A data.table for the combination of a and lambda that minimize the cross validation error across all tau.}
#' \item{gcve}{Group, across all quantiles, cross-validation error results for each value of a and lambda.}
#' \item{call}{Original call to the function.}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100*8,sd=1),ncol=8)
#' y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(100,3)
#' g <- c(1,1,1,1,2,2,3,3)
#' tvals <- c(.25,.75)
#' \dontrun{
#' m1 <- rq.group.pen.cv(x,y,tau=c(.1,.3,.7),groups=g)
#' m2 <- rq.group.pen.cv(x,y,penalty="gAdLASSO",tau=c(.1,.3,.7),groups=g)
#' m3 <- rq.group.pen.cv(x,y,penalty="gSCAD",tau=c(.1,.3,.7),a=c(3,4,5),groups=g)
#' m4 <- rq.group.pen.cv(x,y,penalty="gMCP",tau=c(.1,.3,.7),a=c(3,4,5),groups=g)
#' }
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} and Shaobo Li \email{shaobo.li@ku.edu}
rq.group.pen.cv <- function(x,y,tau=.5,groups=1:ncol(x),lambda=NULL,a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),printProgress=FALSE,...){
	n <- length(y)
	if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
    }
	fit <- rq.group.pen(x,y,tau,lambda=lambda,groups=groups,a=a,...)
	nt <- length(tau)
	na <- length(fit$a)
	nl <- length(fit$lambda)
	min_nl <- min(sapply(fit$models,lambdanum))
	if(min_nl != nl){
		warning("Different models had a different number of lambdas. To avoid this set lambda.discard=FALSE. Results presented are for the shortest lambda sequence")
		#code improvement, bad for loop
		for(i in 1:length(fit$models)){
			fit$models[[i]] <- clearModels(fit$models[[i]],min_nl)
		}
		fit$lambda <- fit$lambda[1:min_nl]
		nl <- min_nl
	}	
	
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
		trainModel <- rq.group.pen(train_x,train_y,tau,groups=groups,lambda=fit$lambda,a=fit$a,alg=fit$alg,...)
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

	returnVal <- list(cverr = foldErrors, cvse = stdErr, fit = fit, btr=btr, gtr=gtr$returnTable, gcve=gtr$gcve,call=match.call())
	class(returnVal) <- "rq.pen.seq.cv"
	returnVal
}
 


#' Cross Validated quantile regression
#'
#' @param x Matrix of predictors.
#' @param y Vector of response values.
#' @param tau  Conditional quantile being modelled.
#' @param lambda  Vector of lambdas. Default is for lambdas to be automatically generated.
#' @param weights Weights for the objective function.
#' @param penalty Type of penalty: "LASSO", "SCAD" or "MCP".
#' @param intercept Whether model should include an intercept. Constant does not need to be included in "x".
#' @param criteria How models will be evaluated. Either cross-validation "CV", BIC "BIC" or large P BIC "PBIC".
#' @param cvFunc If cross-validation is used how errors are evaluated. Check function "check", "SqErr" (Squared Error) or "AE" (Absolute Value).
#' @param nfolds  K for K-folds cross-validation.
#' @param foldid  Group id for cross-validation. Function will randomly generate groups if not specified.
#' @param nlambda Number of lambdas for which models are fit.
#' @param eps Smallest lambda used.
#' @param init.lambda Initial lambda used to find the maximum lambda. Not needed if lambda values are set.
#' @param penVars Variables that should be penalized. With default value of NULL all variables are penalized.
#' @param alg Algorithm that will be used, either linear programming (LP) or coordinate descent (QICD) algorithm from Peng and Wang (2015).
#' @param internal If this is an internal call to this function. 
#' @param ... Additional arguments to be sent to rq.lasso.fit or rq.nc.fit.
#' 
#' @description 
#' Warning: this function is depracated and will not be exported in future rqPen releases. Produces penalized quantile regression models for a range of lambdas and penalty of choice. 
#' If lambda is unselected than an iterative algorithm is used to find a maximum lambda such  that the penalty is large enough to produce an intercept only model. Then range of lambdas 
#' goes from the maximum lambda found to "eps" on the log scale. For non-convex penalties local linear approximation approach used by Wang, Wu and Li to extend LLA as proposed by Zou and Li (2008) to the quantile regression setting.  
#'
#' @return Returns the following:
#' \itemize{
#' \item{models}{List of penalized models fit. Number of models will match number of lambdas and correspond to cv$lambda.}
#' \item{cv}{Data frame with "lambda" and second column is the evaluation based on the criteria selected.}
#' \item{lambda.min}{Lambda which provides the smallest statistic for the selected criteria.}
#' \item{penalty}{Penalty selected.}
#' }
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(800),nrow=100)
#' y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#' cv_model <- cv.rq.pen(x,y)
#' }
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty="LASSO",intercept=TRUE,criteria = "CV",cvFunc="check",nfolds=10,foldid=NULL,nlambda=100,eps=.0001,init.lambda=1,penVars=NULL,alg=ifelse(ncol(x)<50,"LP","QICD"),internal=FALSE,...){
# x is a n x p matrix without the intercept term
# y is a n x 1 vector
# criteria used to select lambda is cross-validation (CV), BIC, or PBIC (large P)
# nfolds: number of folds for cross validation
# foldid: preset id of folds
# penVar: variables to be penalized, default is all non-intercept terms

  # Pre-algorithm setup/ get convenient values
  deprecate_soft("3.0","cv.rq.pen()","rq.pen.cv()")
  if(length(tau)>1){
      stop("cv.rq.pen() only allows for a single value of tau. The new and improved rq.pen.cv() allows for multiple")
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

#' Non-convex penalized quantile regression
#' 
#' @param x Matrix of predictors.
#' @param y Vector of response values.
#' @param tau  Conditional quantile being modelled.
#' @param lambda  Vector of lambdas. Default is for lambdas to be automatically generated.
#' @param weights Weights for the objective function.
#' @param intercept Whether model should include an intercept. Constant does not need to be included in "x".
#' @param penalty Type of penalty: "LASSO", "SCAD" or "MCP".
#' @param a Additional tuning parameter for SCAD and MCP
#' @param iterations Number of iterations to be done for iterative LLA algorithm.
#' @param converge_criteria Difference in betas from iteration process that would satisfy convergence.
#' @param alg Defaults for small p to linear programming (LP), see Wang, Wu and Li (2012) for details. Otherwise a coordinate descent algorithm is used (QICD), see Peng and Wang (2015) for details. Both methods rely on the One-step sparse estimates algorithm.
#' @param penVars Variables that should be penalized. With default value of NULL all variables are penalized.
#' @param internal Whether call to this function has been made internally or not. 
#' @param ... Additional items to be sent to rq.lasso.fit.
#' 
#' @description 
#' Warning: this function is deprecated and will not be exported in future releases. Produces penalized quantile regression models for a range of lambdas and penalty of choice. If lambda is unselected than an iterative algorithm is used to 
#' find a maximum lambda such that the penalty is large enough to produce an intercept only model. Then range of lambdas goes from the maximum lambda found to "eps" on the 
#' log scale. Local linear approximation approach used by Wang, Wu and Li to extend LLA as proposed by Zou and Li (2008) to the quantile regression setting.  
#'
#' @return Returns the following:
#' \itemize{
#' \item{coefficients}{Coefficients from the penalized model.}
#' \item{PenRho}{Penalized objective function value.}
#' \item{residuals}{ Residuals from the model.}
#' \item{rho}{ Objective function evaluation without the penalty.}
#' \item{coefficients}{ Coefficients from the penalized model.} 
#' \item{tau}{ Conditional quantile being modeled.}
#' \item{n}{ Sample size.}  
#' \item{penalty}{ Penalty used, SCAD or MCP.} 
#' \item{penalty}{Penalty selected.}
#' }
#' @export
#'
#' @examples
#' x <- matrix(rnorm(800),nrow=100)
#' y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#' scadModel <- rq.nc.fit(x,y,lambda=1)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} and Adam Maidman. 
#' @references 
#' \itemize{
#' \item Wang, L., Wu, Y. and Li, R. (2012). Quantile regression of analyzing heterogeneity in ultra-high dimension. \emph{J. Am. Statist. Ass}, \bold{107}, 214--222.
#' \item Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.
#' \item Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models. \emph{Ann. Statist.}, \bold{36}, 1509--1533.
#' \item Peng, B. and Wang, L. (2015). An iterative coordinate-descent algorithm for high-dimensional nonconvex penalized quantile regression. \emph{J. Comp. Graph.}, \bold{24}, 676--694.
#' }
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



#' Plots cross validation results from a rq.pen.seq.cv object
#'
#' @param x The rq.pen.seq.cv object
#' @param septau If set to true then optimal tuning parameters are selected seperately for each quantile and there will be a different plot for each quanitle. 
#' @param tau Quantiles of interest.
#' @param logLambda Whether log(lambda) is used for the x-axis
#' @param main Title to the plot
#' @param ... Additional parameters sent to the plot function. 
#' 
#' @description Provides plots of cross-validation results by lambda. If septau is set to TRUE then plots the cross-validation results for each quantile. If septau is set to FALSE
#' then provides one plot for cross-validation results across all quantiles. 
#'
#' @return Plots of the cross validation results by lambda. 
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100*8,sd=1),ncol=8)
#' y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(100,3)
#' m1 <- rq.pen.cv(x,y,tau=c(.1,.3,.7))
#' plot(m1)
#' plot(m1,septau=FALSE)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
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

#' Plot of coefficients of rq.pen.seq object as a function of lambda
#'
#' @param x rq.pen.seq object
#' @param vars Variables of interest
#' @param logLambda Whether lambda should be reported on the log scale
#' @param tau Quantiles of interest
#' @param a Tuning parameter a values of interest.
#' @param lambda Values of lambda of interest. 
#' @param modelsIndex Specific models of interest.
#' @param lambdaIndex Specific lambda values of interest. 
#' @param main Title of the plots. Can be a vector of multiple titles if multiple plots are created. 
#' @param ... Additional arguments sent to plot
#'
#' @return Returns plot(s) of coefficients as they change with lambda.  
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100*8,sd=10),ncol=8)
#' y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(100,3)
#' m1 <- rq.pen(x,y,tau=c(.1,.5,.7),penalty="SCAD",a=c(3,4))
#' plot(m1,a=3,tau=.7)
#' plot(m1)
#' mlist <- list()
#' for(i in 1:6){
#' mlist[[i]] <- paste("Plot",i)
#' }
#' plot(m1,main=mlist)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
plot.rq.pen.seq <- function(x,vars=NULL,logLambda=FALSE,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL,main=NULL, ...){
	models <- getModels(x,tau=tau,a=a,lambda=lambda,modelsIndex=modelsIndex,lambdaIndex=lambdaIndex)
	tm <- models$targetModels
	li <- models$lambdaIndex
	ml <- length(tm)
	if(!is.null(main) & length(main) != ml & length(main) !=1){
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
		  mainText <- main[i]
		}
		if(logLambda){
			lambdas <- log(x$lambda)
			xtext <- expression(Log(lambda))
		} else{
			lambdas <- x$lambda
			xtext <- expression(lambda)
		}
	  maxli <- ncol(tm[[i]]$coefficients)
	  subli <- li[which(li<=maxli)]
		betas <- tm[[i]]$coefficients[-1,subli]
		plot(lambdas[subli], betas[1,], type="n",ylim=c(min(betas),max(betas)),ylab="Coefficient Value",xlab=xtext,main=mainText,...)
		for(i in 1:dim(betas)[1]){
			lines(lambdas[subli], betas[i,],col=i)
		}
	}
	if(ml > 1){
		par(ask=FALSE)	
	}
}

#' Plots of coefficients by lambda for cv.rq.group.pen and cv.rq.pen
#'
#' @param model cv.rq.pen or cv.rq.group.pen object
#' @param voi Index of betas to include. Default is all of them.
#' @param logLambda Plot of lambdas is on the log scale.
#' @param loi Index of lambdas to use, default is all of them. 
#' @param ... Additional arguments to be sent to plot()
#' 
#' @description Warning: this function is deprecated and will not be exported in future versions. 
#'
#' @return Plot of how beta estimates change with lambda.
#' @export
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   x <- matrix(rnorm(800),nrow=100)
#'   y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#'   lassoModels <- cv.rq.pen(x,y)
#'   b_plot <- beta_plots(lassoModels)
#' }
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
beta_plots <- function(model,voi=NULL,logLambda=TRUE,loi=NULL,...){
#voi - index variables of interest
#logLambda - lambdas plotted on log scale
#loi - index of target lambdas
  deprecate_soft("3.0","beta_plots()","print.rq.pen.seq()")
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

#' Plot of how coefficients change with tau
#'
#' @param x A rq.pen.seq or rq.pen.seq.cv object. 
#' @param ... Additional arguments see bytau.plot.rq.pen.seq() or bytau.plot.rq.pen.seq.cv() for more information. 
#'
#' @return Returns the plot of how coefficients change with tau. 
#' @export
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
bytau.plot <- function(x,...){
	UseMethod("bytau.plot")
} 

#' Plot of how coefficients change with tau. 
#'
#' @param x An rq.pen.seq object
#' @param a The tuning parameter a of interest
#' @param lambda The lambda value of interest. 
#' @param lambdaIndex The lambda index of interest. Only specify lambdaIndex or lambda, not both. 
#' @param ... Additional parameters sent to plot()
#'
#' @return A plot of coefficient values by tau. 
#' @export
#'
#' @examples
#'   set.seed(1)
#'   x <- matrix(rnorm(800),nrow=100)
#'   y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#'   lassoModels <- rq.pen(x,y,tau=seq(.1,.9,.1))
#'   bytau.plot(lassoModels,lambda=lassoModels$lambda[5])
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
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

#' Plot of coefficients varying by quantiles for rq.pen.seq.cv object
#'
#' @param x An rq.pen.seq.cv object
#' @param septau Whether optimal tuning parameters are estimated separately for each quantile.
#' @param cvmin Whether the minimum cv error should be used or the one standard error rule. 
#' @param useDefaults Set to FALSE if you want to use something besides minimum cv or 1se. 
#' @param ... Additional parameters sent to plot() 
#' 
#' @description Produces plots of how coefficient estimates vary by quantile for models selected by using cross validation.
#'
#' @return Returns plots of coefficient estimates varying by quantile. 
#' @export
#'
#' @examples
#'  set.seed(1)
#'  x <- matrix(runif(800),nrow=100)
#'  y <- 1 + x[,1] - 3*x[,5] + (1+x[,4])*rnorm(100)
#'  lmcv <- rq.pen.cv(x,y,tau=seq(.1,.9,.1))
#'  bytau.plot(lmcv)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
bytau.plot.rq.pen.seq.cv <- function(x,septau=TRUE,cvmin=TRUE,useDefaults=TRUE,...){
	coefs <- do.call(rbind,coefficients(x,septau,cvmin,TRUE,tau=x$fit$tau))
	if(ncol(coefs) != length(x$fit$tau)){
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


#' Plots of cross validation results as a function of lambda. 
#'
#' @param model A cv.rq.pen() object.
#' @param logLambda Whether lambda values should be logged or not. 
#' @param loi Lambda indexes of interest, if null all lambda values will be used. 
#' @param ... Additional parameters sent to plot function.
#'
#' @return returns a cross validation plot
#' @export
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
cv_plots <- function(model,logLambda=TRUE,loi=NULL,...){
#logLambda - lambdas plotted on log scale
#loi - index of target lambdas
  deprecate_soft("3.0","cv_plots()", "plot.rq.pen.seq.cv()")
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

#' Old cross validation function for group penalty
#' 
#' 
#'
#' @param x Matrix of predictors. 
#' @param y Vector of responses. 
#' @param groups Vector of groups.
#' @param tau Quantile being modeled.
#' @param lambda Vector of lambdas. Default is for lambdas to be automatically generated.
#' @param penalty Type of penalty: "LASSO", "SCAD" or "MCP".
#' @param intercept Whether model should include an intercept. Constant does not need to be included in "x".
#' @param criteria How models will be evaluated. Either cross-validation "CV", BIC "BIC" or large P BIC "PBIC".
#' @param cvFunc If cross-validation is used how errors are evaluated. Check function "check", "SqErr" (Squared Error) or "AE" (Absolute Value).
#' @param nfolds K for K-folds cross-validation.
#' @param foldid Group id for cross-validation. Function will randomly generate groups if not specified.
#' @param nlambda Number of lambdas for which models are fit.
#' @param eps Multiple of lambda max for Smallest lambda used.
#' @param init.lambda Initial lambda used to find the maximum lambda. Not needed if lambda values are set.
#' @param alg Algorithm used for fit. "QICD" or "LP".
#' @param penGroups Specify which groups will be penalized. Default is to penalize all groups.
#' @param ... Additional arguments to be sent to rq.group.fit or groupQICDMultLambda.   
#' 
#' @description 
#' This function is deprecated. Recommend using rq.group.pen.cv() instead. 
#'
#' @return
#' Returns the following: 
#' \itemize{         
#' \item{beta}{ Matrix of coefficients for different values of lambda}
#' \item{residuals}{ Matrix of residuals for different values of lambda.}
#' \item{rho}{Vector of rho, unpenalized portion of the objective function, for different values of lambda.}
#' \item{cv}{ Data frame with "lambda" and second column is the evaluation based on the criteria selected.}
#' \item{lambda.min}{ Lambda which provides the smallest statistic for the selected criteria.}
#' \item{penalty}{ Penalty selected.} 
#' \item{intercept}{Whether intercept was included in model.}
#' \item{groups}{Group structure for penalty function.}
#' }
#' 
#' @examples 
#' \dontrun{
#' x <- matrix(rnorm(800),nrow=100)
#' y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#' cv_model <- cv.rq.group.pen(x,y,groups=c(rep(1,4),rep(2,4)),criteria="BIC")
#' }
#' @references 
#' \itemize{
#' \item Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.
#' \item Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.
#' }
#' @export
#'
cv.rq.group.pen <- function (x, y, groups, tau = 0.5, lambda = NULL, penalty = "SCAD", 
    intercept = TRUE, criteria = "CV", cvFunc = "check", nfolds = 10, 
    foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="QICD",penGroups=NULL,
    ...) 
{
  deprecate_soft("3.0","cv.rq.group.pen()","rq.group.pen.cv()")
  #warning("Recommend that you use rq.group.pen.cv() instead. This is an older and slower version that is only kept for reproducibality. It will not be exported to the namespace in future versions.")
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
    lasso_fit <- cv.rq.group.pen(x,y,groups,tau,lambda,penalty="LASSO",intercept,criteria="BIC",cvFunc,nfolds,foldid,
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

#' Estimates a quantile regression model with a group penalized objective function.
#'
#' @param x Matrix of predictors.
#' @param y Vector of responses.
#' @param groups Vector of group assignments.
#' @param tau Single quantile to be modeled.
#' @param lambda Single value or seperate value for each group.
#' @param intercept Whether intercept should be included in the model or not.
#' @param penalty Type of penalty used: SCAD, MCP or LASSO. 
#' @param alg Type of algorithm used: QICD or LP. 
#' @param a Additional tuning parameter for SCAD and MCP. 
#' @param penGroups Vector of TRUE and FALSE entries for each group determing if they should be penalized. Default is TRUE for all groups.
#' @param ... Additional arguments sent to rq.group.lin.prog()
#'
#' @return Returns the following:      
#' \itemize{    
#' \item{coefficients}{Coefficients of the model.}
#' \item{residuals}{ Residuals from the fitted model.}
#' \item{rho}{Unpenalized portion of the objective function.}
#' \item{tau}{ Quantile being modeled.}
#' \item{n}{Sample size.}
#' \item{intercept}{Whether intercept was included in model.}
#' }
#' 
#' @description Warning: function is deprecated and will not be exported in future R packages. Recommend using rq.group.pen() instead. 
#' Similar to cv.rq.pen function, but uses group penalty. Group penalties use the L1 norm instead of L2 for computational convenience. 
#' As a result of this the group lasso penalty is the same as the typical lasso penalty and thus you should only use a SCAD or MCP penalty. 
#' Only the SCAD and MCP penalties incorporate the group structure into the penalty. The group lasso penalty is implemented because it is 
#' needed for the SCAD and MCP algorithm. We use a group penalty extension of the QICD algorithm presented by Peng and Wang (2015). 
#' @export
#' 
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} and Adam Maidman
#' 
#' @references 
#' \itemize{
#' \item Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.
#' \item Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.
#' }
rq.group.fit <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                penalty = "SCAD", alg="QICD", a=3.7,penGroups=NULL, ...) 
{
  deprecate_soft("3.0","rq.group.fit()","rq.group.pen()")
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

#' Cross validation plot for cv.rq.group.pen object
#'
#' @param x A cv.rq.group.pen object
#' @param ... Additional parameters for plot function.
#'
#' @return A cross validation plot. 
#' @export
#'
plot.cv.rq.group.pen <- function (x,...) 
{
    deprecate_soft("3.0","plot.cv.rq.group.pen()","plot.rq.group.pen.seq()")
    plot(x$cv[, 1], x$cv[, 2],...)
}

#' Estimates a quantile regression model with a lasso penalized quanitle loss function. 
#'
#' @param x Matrix of predictors. 
#' @param y Vector of responses.
#' @param tau Quantile of interest. 
#' @param lambda Tuning parameter. 
#' @param weights Weights for the objective function.
#' @param intercept Whether model should include an intercept. Constant does not need to be included in "x".
#' @param coef.cutoff Coefficients below this value will be set to zero.
#' @param method Use method "br" or "fn" as outlined in quantreg package. We have found "br" to be more stable for penalized regression problems.
#' @param penVars Variables that should be penalized. With default value of NULL all variables are penalized.
#' @param scalex If set to true the predictors will be scaled to have mean zero and standard deviation of one before fitting the model. The output returned will be on the original scale of the data.
#' @param lambda.discard If TRUE lambda sequence will stop early if for small values of lambda the estimates do not change much. 
#' @param ... Additional items to be sent to rq. Note this will have to be done carefully as rq is run on the augmented data to account for penalization and could provide strange results if this is not taken into account.
#'
#' @return
#' Returns the following:
#' \itemize{
#' \item{coefficients}{ Coefficients from the penalized model.} 
#' \item{PenRho}{ Penalized objective function value.}
#' \item{residuals}{ Residuals from the model.}
#' \item{rho}{ Objective function evaluation without the penalty.}
#' \item{tau}{ Conditional quantile being modeled.}
#' \item{n}{ Sample size.}  
#' }
#' 
#' @description Fits a quantile regression model with the LASSO penalty. Uses the augmented data approach similar to the proposal in Sherwood and Wang (2016).   
#' 
#' @export
#'
#' @examples
#' x <- matrix(rnorm(800),nrow=100)
#' y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#' lassoModel <- rq.lasso.fit(x,y,lambda=.1)
#' 
#' @references 
#' \itemize{
#' \item Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. \emph{Journal of the Royal Statistical Society. Series B}, \bold{58}, 267--288.
#' \item Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.  
#' \item Sherwood, B. and Wang, L. (2016) Partially linear additive quantile regression in ultra-high dimension. \emph{Annals of Statistics} \bold{44}, 288--317. 
#' }
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
  deprecate_soft("3.0","rq.lasso.fit()","rq.pen()")
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

#' Prediction for a rq.pen object
#'
#' @param object An rq.pen object. 
#' @param newx Matrix of new data to make predictions with. 
#' @param ... Additional parameters that are currenlty ignored
#'
#' @return A vector of predictions. 
#' 
#' @description This function is deprecated and will not be exported in future versions. 
#' 
#' @export
#'
predict.rq.pen <- function(object, newx,...){
  deprecate_soft("3.0","predict.rq.pen()","predict.rq.pen.seq()")
  coefs <- object$coefficients
  if(object$intercept){
     newx <- cbind(1,newx)
  }
  newx %*% coefs
}

#' Prediction for a cv.rq.pen object
#'
#' @param object A cv.rq.pen object. 
#' @param newx Matrix of new data to make predictions with. 
#' @param lambda Lambda value used, default is the value associated with the minimum cross validation result. 
#' @param ... Additional parameters that are currenlty ignored
#'
#' @return A vector of predictions. 
#' 
#' @description This function is deprecated and will not be exported in future versions. 
#' 
#' @export
#'
predict.cv.rq.pen <- function(object, newx, lambda="lambda.min",...){
  deprecate_soft("3.0","predict.cv.rq.pen()","predict.rq.pen.seq.cv()")
  if(lambda == "lambda.min"){
     target_pos <- which(object$cv$lambda == object$lambda.min)
  } else{
     target_pos <- which(object$cv$lambda == lambda)
  }
  predict(object$models[[target_pos]],newx,...)
}


#' Coefficients from a cv.rq.group.pen object
#'
#' @param object A cv.rq.group.pen object.
#' @param lambda The lambda value, default is to use the one associated with the minimum cv error. 
#' @param ... Additional parameters. 
#'
#' @return Vector of coefficients. 
#' @export
#'
coef.cv.rq.group.pen <- function(object, lambda='min',...){
  deprecate_soft("3.0","coef.cv.rq.group.pen()","coef.rq.pen.seq.cv()")
  if(lambda=='min'){
     lambda <- object$lambda.min
  }
  target_model <- which(object$cv[,1] == lambda)
  object$beta[,target_model]
}


#' Fits quantile regression models using a group penalized objective function.
#'
#' @param x Matrix of predictors. 
#' @param y Vector of responses.
#' @param tau Vector of quantiles. 
#' @param groups Vector of group assignments for predictors. 
#' @param penalty Penalty used, choices are group lasso ("gLASSO"), group adaptive lasso ("gAdLASSO"), group SCAD ("gSCAD") and group MCP ("gMCP")
#' @param lambda Vector of lambda tuning parameters. Will be autmoatically generated if it is not set. 
#' @param nlambda The number of lambda tuning parameters.
#' @param eps The value to be multiplied by the largest lambda value to determine the smallest lambda value. 
#' @param alg Algorithm used. Choices are Huber approximation ("huber"), linear programming ("lp") or quantile iterative coordinate descent ("qicd").
#' @param a The additional tuning parameter for adaptive lasso, SCAD and MCP. 
#' @param norm Whether a L1 or L2 norm is used for the grouped coefficients. 
#' @param group.pen.factor Penalty factor for each group.
#' @param tau.penalty.factor Penalty factor for each quantile.
#' @param scalex Whether X should be centered and scaled so that the columns have mean zero and standard deviation of one. If set to TRUE, the coefficients will be returned to the original scale of the data.
#' @param coef.cutoff Coefficient cutoff where any value below this number is set to zero. Useful for the lp algorithm, which are prone to finding almost, but not quite, sparse solutions. 
#' @param max.iter The maximum number of iterations for the algorithm. 
#' @param converge.eps The convergence criteria for the algorithms. 
#' @param gamma The tuning parameter for the Huber loss. 
#' @param lambda.discard Whether lambdas should be discarded if for small values of lambda there is very little change in the solutions. 
#' @param ... Additional parameters 
#' 
#' @description  
#' Let the predictors be divided into G groups with G corresponding vectors of coefficients, \eqn{\beta_1,\ldots,\beta_G}. 
#' Let \eqn{\rho_\tau(a) = a[\tau-I(a<0)]}. Fits quantile regression models for Q quantiles by minimizing the penalized objective function of
#' \deqn{\sum_{q=1}^Q \frac{1}{n} \sum_{i=1}^n \rho_\tau(y_i-x_i^T\beta^q) + \sum_{q=1}^Q  \sum_{g=1}^G P(||\beta^q_g||_k,w_q*v_j*\lambda,a).}
#' Where \eqn{w_q} and \eqn{v_j} are designated by penalty.factor and tau.penalty.factor respectively. The value of \eqn{k} is chosen by \code{norm}.
#' Value of P() depends on the penalty. Briefly, but see references or vignette for more details,
#' \itemize{
#' \item{Group LASSO (gLASSO)}{\eqn{P(||\beta||_k,\lambda,a)=\lambda||\beta||_k}}
#' \item{Group SCAD}{\eqn{P(||\beta||_k,\lambda,a)=SCAD(||\beta||_k,\lambda,a)}}
#' \item{Group MCP}{\eqn{P(||\beta||_k,\lambda,a)=MCP(||\beta||_k,\lambda,a)}}
#' \item{Group Adaptive LASSO}{\eqn{P(||\beta||_k,\lambda,a)=\frac{\lambda ||\beta||_k}{|\beta_0|^a}}}
#' }
#' Note if \eqn{k=1} and the group lasso penalty is used then this is identical to the regular lasso and thus function will stop and
#' suggest that you use rq.pen() instead. For Adaptive LASSO the values of \eqn{\beta_0} come from a Ridge solution with the same value of \eqn{\lambda}.
#'  If the Huber algorithm is used than \eqn{\rho_\tau(y_i-x_i^\top\beta)} is replaced by a Huber-type approximation. Specifically, it is replaced by \eqn{h^\tau_\gamma(y_i-x_i^\top\beta)/2} where 
#' \deqn{h^\tau_\gamma(a) = a^2/(2\gamma)I(|a| \leq \gamma) + (|a|-\gamma/2)I(|a|>\gamma)+(2\tau-1)a.}
#' Where if \eqn{\tau=.5}, we get the usual Huber loss function. 
#' 
#' @return An rq.pen.seq object. 
#' \itemize{
#' \item{models}{A list of each model fit for each tau and a combination.}
#' \item{n}{Sample size.}
#' \item{p}{Number of predictors.}
#' \item{alg}{Algorithm used.}
#' \item{tau}{Quantiles modeled.}
#' \item{penalty}{Penalty used.}
#' \item{a}{Tuning parameters a used.}
#' \item{lambda}{Lambda values used for all models. If a model has fewer coefficients than lambda, say k. Then it used the first k values of lambda. Setting lambda.discard to TRUE will gurantee all values use the same lambdas, but may increase computational time noticeably and for little gain.}
#' \item{modelsInfo}{Information about the quantile and a value for each model.}
#' \item{call}{Original call.}
#' }
#' Each model in the models list has the following values. 
#' \itemize{
#' \item{coefficients}{Coefficients for each value of lambda.}
#' \item{rho}{The unpenalized objective function for each value of lambda.}
#' \item{PenRho}{The penalized objective function for each value of lambda.}
#' \item{nzero}{The number of nonzero coefficients for each value of lambda.}
#' \item{tau}{Quantile of the model.}
#' \item{a}{Value of a for the penalized loss function.}
#' }
#' 

#'
#' @return
#' @export
#'
#' @examples
#'  
#' set.seed(1)
#' x <- matrix(rnorm(200*8,sd=1),ncol=8)
#' y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(200,3)
#' g <- c(1,1,1,2,2,2,3,3)
#' tvals <- c(.25,.75)
#' r1 <- rq.group.pen(x,y,groups=g)
#' r5 <- rq.group.pen(x,y,groups=g,tau=tvals)
#' #Linear programming approach with group SCAD penalty and L1-norm
#' m2 <- rq.group.pen(x,y,groups=g,alg="lp",penalty="gSCAD",norm=1,a=seq(3,4))
#' # No penalty for the first group
#' m3 <- rq.group.pen(x,y,groups=g,group.pen.factor=c(0,rep(1,2)))
#' # Smaller penalty for the median
#' m4 <- rq.group.pen(x,y,groups=g,tau=c(.25,.5,.75),tau.penalty.factor=c(1,.25,1))
#' 
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}, Shaobo Li \email{shaobo.li@ku.edu} and Adam Maidman
#' @references 
#' \insertRef{qr_cd}{rqPen}
rq.group.pen <- function(x,y, tau=.5,groups=1:ncol(x), penalty=c("gLASSO","gAdLASSO","gSCAD","gMCP"),
						lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.05,.01),alg=c("huber","lp","qicd"), 
						a=NULL, norm=2, group.pen.factor=rep(1,length(unique(groups))),tau.penalty.factor=rep(1,length(tau)),
						scalex=TRUE,coef.cutoff=1e-8,max.iter=10000,converge.eps=1e-7,gamma=IQR(y)/10, lambda.discard=TRUE, ...){
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
			if(a[1] != 1 | length(a) !=1){
				stop("The tuning parameter a is not used for group lasso. Leave it set to NULL or set to 1")
			}
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
	if(!is.matrix(x)){
	  stop("x must be a matrix")
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
	if(any(tau.penalty.factor <= 0) | any(group.pen.factor < 0)){
	  stop("group penalty factors must be positive and tau penalty factors must be non-negative")
	}
	if(sum(group.pen.factor)==0){
	  stop("Some group penalty factors must be non-zero")
	}
	if(lpf!=g){
		stop("group penalty factor must be of length g")
	}	
	if(is.null(lambda)){
		lamMax <- getLamMaxGroup(x,y,groups,tau,group.pen.factor,penalty=penalty,scalex=scalex,tau.penalty.factor=tau.penalty.factor,norm=norm,gamma=gamma)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}

	penalty.factor <- mapvalues(groups,seq(1,g),group.pen.factor)
	
	if(penalty == "gLASSO"){
		return_val <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma,lambda.discard=lambda.discard,...)
	} else{
		if(penalty == "gAdLASSO"){
			init.model <- rq.enet(x,y,tau,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,tau.penalty.factor=tau.penalty.factor,
									a=0,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,lambda.discard=lambda.discard,...)
		} else{
			if(norm == 1){
				if(alg == "qicd"){
					init.alg <- "lp"
				} else{
					init.alg <- alg
				}
				init.model <- rq.lasso(x,y,tau,alg=init.alg,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,
							tau.penalty.factor=tau.penalty.factor,max.iter=max.iter,coef.cutoff=coef.cutoff,converge.eps=converge.eps,
							gamma=gamma,lambda.discard=lambda.discard,...)
			} else{
				init.model <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma,lambda.discard=lambda.discard,...)
			}
		}
		return_val <- rq.group.lla(init.model,x,y,groups,penalty=penalty,a=a,norm=norm,group.pen.factor=group.pen.factor,tau.penalty.factor=tau.penalty.factor,scalex=scalex,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,lambda.discard=lambda.discard,...)
	}
	class(return_val) <- "rq.pen.seq"
	return_val$call <- match.call()	
	return_val$lambda <- lambda
	if(lambda.discard){
		lmax <- max(sapply(return_val$models,lambdanum))
		return_val$lambda <- return_val$lambda[1:lmax]
	}
	return_val
}

#' Prints a cv.rq.pen object
#'
#'
#' @param x A cv.rq.pen object
#' @param ... Additional arguments
#' 
#' @description Warning: this function is deprecated and will not be exported in future releases. 
#'
#' @return Prints coefficients and cross validation results. 
#' @export
#'
print.cv.rq.pen <- function(x,...){
   deprecate_soft("3.0","print.cv.rq.pen()","print.rq.pen.seq.cv()")
   cat("\nCoefficients:\n")
   print(coefficients(x,...))
   cat("\nCross Validation (or BIC) Results\n")
   print(x$cv)
}

#' Prints a rq.pen object
#'
#' @param x A rq.pen object
#' @param ... Additional arguments
#'
#' @return
#' @export
#'
print.rq.pen <- function(x,...){
  deprecate_soft("3.0","print.rq.pen()","print.rq.pen.seq()")
  cat("\nCoefficients:\n")
	print(coefficients(x,...))
}


