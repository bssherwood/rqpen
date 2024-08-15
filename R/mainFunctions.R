#' Calculate information criterion for penalized quantile regression models. Currently not exported. 
#'
#' @param model model from a rq.pen.seq() object
#' @param n Sample size
#' @param method Choice of BIC, AIC or PBIC, a large p BIC. 
#'
#' @return 
#' Let \eqn{\hat{\beta}} be the coefficient vectors for the estimated model. Function returns the value 
#' \deqn{\log(\sum_{i=1}^n w_i \rho_\tau(y_i-x_i^\top\hat{\beta})) + d*b/(2n),} where d is the number of nonzero coefficients and b depends on the method used. For AIC \eqn{b=2},
#' for BIC \eqn{b=log(n)} and for PBIC \eqn{d=log(n)*log(p)} where p is the dimension of \eqn{\hat{\beta}}. The values of w_i default to one and are set using weights when fitting the models. Returns this value for each coefficient vector in the model, so one
#' for every value of \eqn{\lambda}. 
#' @keywords internal
#' @examples \dontrun{
#' set.seed(1)
#' x <- matrix(runif(800),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(100)
#' m1 <- rq.pen(x,y,tau=c(.25,.75))
#' # returns the IC values for tau=.25
#' qic(m1$models[[1]],m1$n) 
#' # returns the IC values for tau=.75
#' qic(m1$models[[2]],m1$n)
#' } 
#' @references 
#' \insertRef{qrbic}{rqPen}
#'@author Ben Sherwood, \email{ben.sherwood@ku.edu}
qic <- function(model,n, method=c("BIC","AIC","PBIC")){
  method <- match.arg(method)
	tau <- model$debug
	df <- model$nzero
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


#' Select tuning parameters using IC 
#' 
#' Selects tuning parameter \eqn{\lambda} and a according to information criterion of choice. For a given \eqn{\hat{\beta}} the information criterion is calculated
#' as
#' \deqn{\log(\sum_{i=1}^n w_i \rho_\tau(y_i-x_i^\top\hat{\beta})) + d*b/(2n),} where d is the number of nonzero coefficients and b depends on the method used. For AIC \eqn{b=2},
#' for BIC \eqn{b=log(n)} and for PBIC \eqn{d=log(n)*log(p)} where p is the dimension of \eqn{\hat{\beta}}.
#' If septau set to FALSE then calculations are made across the quantiles. Let \eqn{\hat{\beta}^q} be the coefficient vector for the qth quantile of Q quantiles. In addition let \eqn{d_q} and \eqn{b_q} 
#' be d and b values from the qth quantile model. Note, for all of these we are assuming eqn and a are the same. Then the summary across all quantiles is 
#' \deqn{\sum_{q=1}^Q w_q[ \log(\sum_{i=1}^n m_i \rho_\tau(y_i-x_i^\top\hat{\beta}^q)) + d_q*b_q/(2n)],}
#' where \eqn{w_q} is the weight assigned for the qth quantile model. 
#'
#'
#' @param obj A rq.pen.seq or rq.pen.seq.cv object. 
#' @param ... Additional arguments see qic.select.rq.pen.seq() or qic.select.rq.pen.seq.cv() for more information. 
#'
#' @return Returns a qic.select object. 
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
qic.select <- function(obj,...){
  UseMethod("qic.select")
} 

#' Select tuning parameters using IC
#' 
#' Selects tuning parameter \eqn{\lambda} and a according to information criterion of choice. For a given \eqn{\hat{\beta}} the information criterion is calculated
#' as
#' \deqn{\log(\sum_{i=1}^n w_i \rho_\tau(y_i-x_i^\top\hat{\beta})) + d*b/(2n),} where d is the number of nonzero coefficients and b depends on the method used. For AIC \eqn{b=2},
#' for BIC \eqn{b=log(n)} and for PBIC \eqn{d=log(n)*log(p)} where p is the dimension of \eqn{\hat{\beta}}.
#' If septau set to FALSE then calculations are made across the quantiles. Let \eqn{\hat{\beta}^q} be the coefficient vector for the qth quantile of Q quantiles. In addition let \eqn{d_q} and \eqn{b_q} 
#' be d and b values from the qth quantile model. Note, for all of these we are assuming eqn and a are the same. Then the summary across all quantiles is 
#' \deqn{\sum_{q=1}^Q w_q[ \log(\sum_{i=1}^n m_i \rho_\tau(y_i-x_i^\top\hat{\beta}^q)) + d_q*b_q/(2n)],}
#' where \eqn{w_q} is the weight assigned for the qth quantile model. 
#'
#' @param obj A rq.pen.seq or rq.pen.seq.cv object. 
#' @param method Choice of BIC, AIC or PBIC, a large p BIC.
#' @param septau If optimal values of \eqn{\lambda} and a can vary with \eqn{\tau}. Default is TRUE. 
#' @param tauWeights Weights for each quantile. Useful if you set septau to FALSE but want different weights for the different quantiles. If not specified default is to have \eqn{w_q=1} for all quantiles.
#' @param ... Additional arguments.
#'
#' @return 
#' \describe{
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
qic.select.rq.pen.seq <- function(obj, method=c("BIC","AIC","PBIC"),septau=ifelse(obj$penalty!="gq",TRUE,FALSE),tauWeights=NULL,...){
# code help: Maybe think about how the qic values are returned for the septau=TRUE case. Also, potential issue with different values of lambda
	method <- match.arg(method)
	if(is.null(tauWeights)==FALSE & septau){
		warning("Weights are only used when septau is set to true.")
	}
	if(obj$penalty=="gq" & septau){
	  septau = FALSE
	  warning("septau set to false because group quantile penalty was used, which is a joint optimization across all quantiles")
	}
	if(is.null(tauWeights) & !septau){
		tauWeights <- rep(1,length(obj$tau))
	}
	
	n <- obj$n
	if(septau & length(tauWeights)==1){
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
		
		coefs <- matrix(nrow=nrow(coef(obj)), ncol=nt)#vector(mode="list", length=nt)
		for(i in 1:nt){
			coefs[,i] <- coef(obj$models[[modelsInfo$modelIndex[i]]])[,modelsInfo$lambdaIndex[i]]
		}
		rownames(coefs) <- names(coef(obj$models[[modelsInfo$modelIndex[i]]])[,modelsInfo$lambdaIndex[i]])
		gic <- NULL
	} else{
		gic <- matrix(rep(0,na*nl),ncol=nl)
		tqic_vals <- t(qic_vals)
		for(i in 1:na){
			subIC <- subset(tqic_vals, obj$modelsInfo$a==obj$a[i])
			gic[i,] <- tauWeights %*% subIC
		}#
		minIndex <- which(gic==min(gic),arr.ind=TRUE)
		returnA <- obj$a[minIndex[1]]
		a <- obj$a
		modelsInfo <- subset(obj$modelsInfo, a==returnA)
		modelsInfo <- cbind(modelsInfo,minQIC=tqic_vals[modelsInfo$modelIndex,minIndex[2]],lambdaIndex=minIndex[2],lambda=obj$lambda[minIndex[2]])
		coefs <- coef(obj, lambdaIndex=minIndex[2], modelsIndex=modelsInfo$modelIndex)
	}
	#coefIndex <- 1:nt
	#modelsInfo <- cbind(modelsInfo, coefIndex)
	#coefs <- do.call(cbind,coefs)
	if(!is.null(ncol(coefs))){
	  colnames(coefs) <- paste0("tau=",obj$tau)
	}
	
	return_val <- list(coefficients = coefs, ic=qic_vals,modelsInfo=modelsInfo, gic=gic)
	class(return_val) <- "qic.select"
	return_val
}

#' Select tuning parameters using IC
#' 
#' Selects tuning parameter \eqn{\lambda} and a according to information criterion of choice. For a given \eqn{\hat{\beta}} the information criterion is calculated
#' as
#' \deqn{\log(\sum_{i=1}^n w_i \rho_\tau(y_i-x_i^\top\hat{\beta})) + d*b/(2n),} where d is the number of nonzero coefficients and b depends on the method used. For AIC \eqn{b=2},
#' for BIC \eqn{b=log(n)} and for PBIC \eqn{d=log(n)*log(p)} where p is the dimension of \eqn{\hat{\beta}}.
#' If septau set to FALSE then calculations are made across the quantiles. Let \eqn{\hat{\beta}^q} be the coefficient vector for the qth quantile of Q quantiles. In addition let \eqn{d_q} and \eqn{b_q} 
#' be d and b values from the qth quantile model. Note, for all of these we are assuming eqn and a are the same. Then the summary across all quantiles is 
#' \deqn{\sum_{q=1}^Q w_q[ \log(\sum_{i=1}^n  \rho_\tau(y_i-x_i^\top\hat{\beta}^q)) + d_q*b_q/(2n)],}
#' where \eqn{w_q} is the weight assigned for the qth quantile model. 
#'
#' @param obj A rq.pen.seq.cv object. 
#' @param method Choice of BIC, AIC or PBIC, a large p BIC.
#' @param septau If optimal values of \eqn{\lambda} and a can vary with \eqn{\tau}. Default is TRUE. 
#' @param weights Weights for each quantile. Useful if you set septau to FALSE but want different weights for the different quantiles. If not specified default is to have \eqn{w_q=1} for all quantiles.
#' @param ... Additional arguments. 
#'
#' @return 
#' \describe{
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
#' m1 <- rq.pen.cv(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.75))
#' qic.select(m1)
#' @references 
#' \insertRef{qrbic}{rqPen}
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
qic.select.rq.pen.seq.cv <- function(obj, method=c("BIC","AIC","PBIC"),septau=ifelse(obj$fit$penalty!="gq",TRUE,FALSE),weights=NULL,...){
  # code help: Maybe think about how the qic values are returned for the septau=TRUE case. Also, potential issue with different values of lambda
  #print("septau",septau)
  septau = septau
  method <- match.arg(method)
  #print("septau",septau)
  obj <- obj$fit
  #print("septau",septau)
  if(is.null(weights)==FALSE & septau){
    warning("Weights are only used when septau is set to true.")
  }
  if(obj$penalty=="gq" & septau){
    septau = FALSE
    warning("septau set to false because group quantile penalty was used, which is a joint optimization across all quantiles")
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
    
    coefs <- matrix(nrow=nrow(coef(obj)), ncol=nt)#vector(mode="list", length=nt)
    for(i in 1:nt){
      coefs[,i] <- coef(obj$models[[modelsInfo$modelIndex[i]]])[,modelsInfo$lambdaIndex[i]]
    }
	rownames(coefs) <- names(coef(obj$models[[modelsInfo$modelIndex[i]]])[,modelsInfo$lambdaIndex[i]])
    gic <- NULL
  } else{
    gic <- matrix(rep(0,na*nl),ncol=nl)
    tqic_vals <- t(qic_vals)
    for(i in 1:na){
      subIC <- subset(tqic_vals, obj$modelsInfo$a==obj$a[i])
      gic[i,] <- weights %*% subIC
    }#
    minIndex <- which(gic==min(gic),arr.ind=TRUE)
    returnA <- obj$a[minIndex[1]]
    a <- obj$a
    modelsInfo <- subset(obj$modelsInfo, a==returnA)
    modelsInfo <- cbind(modelsInfo,minQIC=tqic_vals[modelsInfo$modelIndex,minIndex[2]],lambdaIndex=minIndex[2],lambda=obj$lambda[minIndex[2]])
    coefs <- coef(obj, lambdaIndex=minIndex[2], modelsIndex=modelsInfo$modelIndex)
  }
  #coefIndex <- 1:nt
  #modelsInfo <- cbind(modelsInfo, coefIndex)
  #coefs <- do.call(cbind,coefs)
  if(!is.null(ncol(coefs))){
    colnames(coefs) <- paste0("tau=",obj$tau)
  }
  
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
#' @method print qic.select
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
print.qic.select <- function(x,...){
   print(coefficients(x))
}

#' Predictions from a qic.select object
#'
#' @param object qic.select object
#' @param newx Data matrix to make predictions from. 
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
predict.qic.select <- function(object, newx, ...){
	if(is.null(dim(newx))){
	  c(1,newx) %*% coefficients(object)
	} else{
	  cbind(1,newx) %*% coefficients(object)
	}
}


#' Print a rq.pen.seq object
#'
#' @param x rq.pen.seq object
#' @param ... optional arguments
#'
#' @return If only one model, prints a data.frame of the number of nonzero coefficients and lambda. Otherwise prints information about the quantiles being modeled and choices for a.
#' @export
#'
#' @method print rq.pen.seq
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
print.rq.pen.seq <- function(x,...){
  nt <- length(x$tau)
  na <- length(x$a)
  if(nt==1 & na==1){
    print(data.frame(nzero=x$models[[1]]$nzero,lambda=x$lambda))
  } else{
	cat("\n Number of nonzero coefficients by lambda for all models\n")
	print(data.frame(lambda=x$lambda,sapply(x$models,getNZero)))
  }
  # }else if(nt > 1 & na > 1){
    # print(paste(c(paste(c("Quantile regression with ", x$penalty, " penalty for quantiles:",x$tau), collapse=" ")," and tuning parameters a:", x$a),collapse=" "))
  # } else if( na > 1){
    # print(paste(c(paste(c("Quantile regression with ", x$penalty, " penalty for quantile:",x$tau), collapse=" ")," and tuning parameters a:", x$a),collapse=" "))
  # } else{
    # print(paste(c("Quantile regression with ", x$penalty, " penalty for quantiles:",x$tau), collapse=" "))
  # }	
}


#' Returns Coefficients of a cv.rq.pen object
#' 
#' Warning: this function is no longer exported, due to the switch from cv.rq.pen() to rq.pen.cv().
#'
#' @param object cv.rq.pen object 
#' @param lambda Value of lambda, default is to use the minimum value. 
#' @param ... Additional parameters.
#'
#' @return Coefficients for a given lambda, or the lambda associated with the minimum cv value. 
#' 
#' @keywords internal
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
coef.cv.rq.pen <- function(object, lambda='min',...){
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
#' @param weights Weights for the quantile objective function.  
#' @param ... Extra parameters. 
#' 
#' @description  
#' Let q index the Q quantiles of interest. Let \eqn{\rho_\tau(a) = a[\tau-I(a<0)]}. Fits quantile regression models by minimizing the penalized objective function of
#' \deqn{\frac{1}{n} \sum_{q=1}^Q \sum_{i=1}^n m_i \rho_\tau(y_i-x_i^\top\beta^q) + \sum_{q=1}^Q  \sum_{j=1}^p P(\beta^q_p,w_q*v_j*\lambda,a).}
#' Where \eqn{w_q} and \eqn{v_j} are designated by penalty.factor and tau.penalty.factor respectively, and \eqn{m_i} is designated by weights. Value of \eqn{P()} depends on the penalty. See references or vignette for more details,
#' \describe{
#' \item{LASSO:}{ \eqn{P(\beta,\lambda,a)=\lambda|\beta|}}
#' \item{SCAD:}{ \eqn{P(\beta,\lambda,a)=SCAD(\beta,\lambda,a)}}
#' \item{MCP:}{ \eqn{P(\beta,\lambda,a)=MCP(\beta,\lambda,a)}}
#' \item{Ridge:}{ \eqn{P(\beta,\lambda,a)=\lambda\beta^2}}
#' \item{Elastic Net:}{ \eqn{P(\beta,\lambda,a)=a*\lambda|\beta|+(1-a)*\lambda*\beta^2}}
#' \item{Adaptive LASSO:}{ \eqn{P(\beta,\lambda,a)=\frac{\lambda |\beta|}{|\beta_0|^a}}}
#' }
#' For Adaptive LASSO the values of \eqn{\beta_0} come from a Ridge solution with the same value of \eqn{\lambda}. Three different algorithms are implemented
#' \describe{
#' \item{huber:}{ Uses a Huber approximation of the quantile loss function. See Yi and Huang 2017 for more details.}
#' \item{br:}{ Solution is found by re-formulating the problem so it can be solved with the rq() function from quantreg with the br algorithm.} 
#' }
#' The huber algorithm offers substantial speed advantages without much, if any, loss in performance. However, it should be noted that it solves an approximation of the quantile loss function.
#' @return An rq.pen.seq object. 
#' \describe{
#' \item{models: }{ A list of each model fit for each tau and a combination.}
#' \item{n:}{ Sample size.}
#' \item{p:}{ Number of predictors.}
#' \item{alg:}{ Algorithm used. Options are "huber" or any method implemented in rq(), such as "br". }
#' \item{tau:}{ Quantiles modeled.}
#' \item{a:}{ Tuning parameters a used.}
#' \item{modelsInfo:}{ Information about the quantile and a value for each model.}
#' \item{lambda:}{ Lambda values used for all models. If a model has fewer coefficients than lambda, say k. Then it used the first k values of lambda. Setting lambda.discard to TRUE will gurantee all values use the same lambdas, but may increase computational time noticeably and for little gain.}
#' \item{penalty:}{ Penalty used.}
#' \item{call:}{ Original call.}
#' }
#' Each model in the models list has the following values. 
#' \describe{
#' \item{coefficients:}{ Coefficients for each value of lambda.}
#' \item{rho:}{ The unpenalized objective function for each value of lambda.}
#' \item{PenRho:}{ The penalized objective function for each value of lambda.}
#' \item{nzero:}{ The number of nonzero coefficients for each value of lambda.}
#' \item{tau:}{ Quantile of the model.}
#' \item{a:}{ Value of a for the penalized loss function.}
#' }
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
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(n)
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
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}, Shaobo Li, and Adam Maidman
rq.pen <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.05,.01), 
	penalty.factor = rep(1, ncol(x)),alg=c("huber","br","fn"),scalex=TRUE,tau.penalty.factor=rep(1,length(tau)),
	coef.cutoff=1e-8,max.iter=5000,converge.eps=1e-4,lambda.discard=TRUE,weights=NULL,...){
	penalty <- match.arg(penalty)
	alg <- match.arg(alg)
	if(length(y)!=nrow(x)){
	  stop("length of x and number of rows in x are not the same")
	}
	if(min(apply(x,2,sd))==0){
		stop("At least one of the x predictors has a standard deviation of zero")
	}
	if(is.null(weights)==FALSE){
	  if( penalty=="ENet" | penalty=="Ridge"){
		  stop("Cannot use weights with elastic net or ridge penalty. Can use it with lasso, though may be much slower than unweighted version.")
	  }
	  if(penalty=="aLASSO"){
	    warning("Weights are ignored when getting initial (Ridge) estimates for adaptive Lasso")
	  }
	  if(length(weights)!=length(y)){
	    stop("number of weights does not match number of responses")
	  }
	  if(sum(weights<=0)>0){
	    stop("all weights most be positive")
	  }
	}
	if(is.matrix(y)==TRUE){
		y <- as.numeric(y)
	}
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
		fit <- rq.lasso(x,y,tau,lambda,nlambda,eps,penalty.factor,alg,scalex=FALSE,tau.penalty.factor,coef.cutoff,max.iter,converge.eps,lambda.discard=lambda.discard,weights=weights,...)
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
		fit <- rq.nc(x,y,tau,penalty,a,lambda,nlambda=nlambda,eps=eps,penalty.factor=penalty.factor,alg=alg,scalex=FALSE,tau.penalty.factor=tau.penalty.factor,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,lambda.discard=lambda.discard,weights=weights,...)
	} 
	#else if(penalty == "gQuant"){
		# fit <- rq.pen.gq1y(x,y,tau,lambda,nlambda=100,eps,scalex=FALSE,penalty.factor=penalty.factor, 
							# tau.penalty.factor=tau.penalty.factor,
							# max.iter=max.iter,converge.eps=converge.eps,lambda.discard=lambda.discard,
							# gamma=gamma)
	# }
	#TODO: this should be done in modelreturn but updating all that code is complicated
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
	#If lambda.discard is used we need to make sure the sequence is the same across all quantiles. 
	#A to-do thing would be to consider how to allow downstream functions such as coef() to allow for different sequence lengths,
	# but that seems unnecessarliy complicated right now. 
		lmin <- min(sapply(fit$models,lambdanum))
		fit$lambda <- fit$lambda[1:lmin]
		for(j in 1:length(fit$models)){
			fit$models[[j]]$coefficients <- fit$models[[j]]$coefficients[,1:lmin]
			fit$models[[j]]$rho <- fit$models[[j]]$rho[1:lmin]
			fit$models[[j]]$PenRho <- fit$models[[j]]$PenRho[1:lmin]
			fit$models[[j]]$nzero <- fit$models[[j]]$nzero[1:lmin]
		}
	}
	fit$weights <- weights
	fit$call <- match.call()
	fit
}


rq.pen.new <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.05,.01), 
                   penalty.factor = rep(1, ncol(x)),alg=c("huber","br","fn"),scalex=TRUE,tau.penalty.factor=rep(1,length(tau)),
                   coef.cutoff=1e-8,max.iter=5000,converge.eps=1e-4,lambda.discard=TRUE,weights=NULL,...){
  penalty <- match.arg(penalty)
  alg <- match.arg(alg)
  if(length(y)!=nrow(x)){
    stop("length of x and number of rows in x are not the same")
  }
  if(min(apply(x,2,sd))==0){
    stop("At least one of the x predictors has a standard deviation of zero")
  }
  if(is.null(weights)==FALSE){
    if( penalty=="ENet" | penalty=="Ridge"){
      stop("Cannot use weights with elastic net or ridge penalty. Can use it with lasso, though may be much slower than unweighted version.")
    }
    if(penalty=="aLASSO"){
      warning("Weights are ignored when getting initial (Ridge) estimates for adaptive Lasso")
    }
    if(length(weights)!=length(y)){
      stop("number of weights does not match number of responses")
    }
    if(sum(weights<=0)>0){
      stop("all weights most be positive")
    }
  }
  if(is.matrix(y)==TRUE){
    y <- as.numeric(y)
  }
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
    fit <- rq.lasso(x,y,tau,lambda,nlambda,eps,penalty.factor,alg,scalex=FALSE,tau.penalty.factor,coef.cutoff,max.iter,converge.eps,lambda.discard=lambda.discard,weights=weights,...)
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
    fit <- rq.nc(x,y,tau,penalty,a,lambda,nlambda=nlambda,eps=eps,penalty.factor=penalty.factor,alg=alg,scalex=FALSE,tau.penalty.factor=tau.penalty.factor,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,lambda.discard=lambda.discard,weights=weights,...)
  } 
  #else if(penalty == "gQuant"){
  # fit <- rq.pen.gq1y(x,y,tau,lambda,nlambda=100,eps,scalex=FALSE,penalty.factor=penalty.factor, 
  # tau.penalty.factor=tau.penalty.factor,
  # max.iter=max.iter,converge.eps=converge.eps,lambda.discard=lambda.discard,
  # gamma=gamma)
  # }
  #TODO: this should be done in modelreturn but updating all that code is complicated
  if(scalex){
    fit$models <- updateModelReturn(x,y,models,penalty,a,tau,tau.penalty.factor,penalty.factor)
  }
  if(lambda.discard){
    #If lambda.discard is used we need to make sure the sequence is the same across all quantiles. 
    #A to-do thing would be to consider how to allow downstream functions such as coef() to allow for different sequence lengths,
    # but that seems unnecessarliy complicated right now. 
    lmin <- min(sapply(fit$models,lambdanum))
    fit$lambda <- fit$lambda[1:lmin]
    for(j in 1:length(fit$models)){
      fit$models[[j]]$coefficients <- fit$models[[j]]$coefficients[,1:lmin]
      fit$models[[j]]$rho <- fit$models[[j]]$rho[1:lmin]
      fit$models[[j]]$PenRho <- fit$models[[j]]$PenRho[1:lmin]
      fit$models[[j]]$nzero <- fit$models[[j]]$nzero[1:lmin]
    }
  }
  fit$weights <- weights
  fit$call <- match.call()
  fit
}

updateModelReturn <- function(x,y,models,penalty,a,tau,tau.penalty.factor,penalty.factor){
  n <- length(y)
  if(is.null(weights)){
    weights <- rep(1,n)
  }
  penfunc <- getPenfunc(penalty)
  x.orig = t(apply(x, 1, function(r)r*attr(x,'scaled:scale') + attr(x, 'scaled:center')))
  
  for(i in 1:length(models)){
    taupos <- which(tau==models[[i]]$tau)
    local.penalty.factor <- penalty.factor*tau.penalty.factor[taupos]
    if(!is.null(dim(models[[i]]$coefficients))){
      models[[i]]$coefficients <- apply(models[[i]]$coefficients,2,transform_coefs,attributes(x)$`scaled:center`,attributes(x)$`scaled:scale`)
      res <- y-cbind(1,x.orig)%*%models[[i]]$coefficients
      models[[i]]$rho <- apply(check(res,models[[i]]$tau)*weights,2,mean)
      for(i in 1:length(models[[i]]$rho)){
        models[[i]]$PenRho[i] <- models[[i]]$rho[i] + sum(penfunc(models[[i]]$coefficients[-1,i],lambda[i]*local.penalty.factor,a))
      }
      
    } else{
      models[[i]]$coefficients <- transform_coefs(models[[i]]$coefficients,attributes(x)$`scaled:center`,attributes(x)$`scaled:scale`)
      res <- y-cbind(1,x.orig)%*%models[[i]]$coefficients
      models[[i]]$rho <- mean(check(res,tau)*weights)	
      models[[i]]$PenRho <-  models[[i]]$rho + sum(penfunc( models[[i]]$coefficients[-1],lambda*local.penalty.factor,a))
    }
  }  
}

rq.pen.modelreturn <- function(coefs,x,y,tau,lambda,local.penalty.factor,penalty,a,weights=NULL){
  
  
  fits <- cbind(1,x)%*% return_val$coefficients
  res <- y - fits
  if(is.null(dim(return_val$coefficients))==TRUE){
    return_val$rho <- mean(check(res,tau)*weights)	
    return_val$PenRho <- return_val$rho + sum(penfunc(return_val$coefficients[-1],lambda*local.penalty.factor,a))
    return_val$nzero <- sum(return_val$coefficients!=0)
  } else{
    return_val$rho <- apply(check(res,tau)*weights,2,mean)	
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
  modelsCombined <- lapply(models$targetModels,getModelCoefs,models$lambdaIndex)
  if(is.null(colnames(modelsCombined[[1]]))){
  #just one lambda for each tau-a combo
    modelNames <- names(modelsCombined)
  } else{
    modelNames <- as.vector( t(outer(names(modelsCombined),colnames(modelsCombined[[1]]),paste)))
  }
  coefReturn <- do.call(cbind,modelsCombined)
  colnames(coefReturn) <- modelNames
  coefReturn
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
#' @return A matrix of predictions for each tau and a combination
#' @export
#'
#' @examples
#' x <- matrix(runif(1600),ncol=8)
#' y <- 1 + x[,1] + x[,8] + (1+.5*x[,3])*rnorm(200)
#' m1 <- rq.pen.cv(x,y,penalty="ENet",a=c(0,.5,1),tau=c(.25,.75),lambda=c(.1,.05,.01))
#' newx <- matrix(runif(80),ncol=8)
#' cvpreds <- predict(m1,newx)
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
predict.rq.pen.seq.cv <- function(object, newx,tau=NULL,septau=ifelse(object$fit$penalty!="gq",TRUE,FALSE),cvmin=TRUE,useDefaults=TRUE,...){
  if(object$fit$penalty=="gq" & septau){
    septau = FALSE
    warning("septau set to false because group quantile penalty was used, which is a joint optimization across all quantiles")
  }
  coefs <- coefficients(object,septau=septau,cvmin=cvmin,useDefaults=useDefaults,tau=tau,...)
  if(is.null(dim(newx))){
    c(1,newx) %*% coefs  
  } else{
    cbind(1,newx) %*% coefs
  }
}

#' Predictions from rq.pen.seq object
#'
#' @param object rq.pen.seq object
#' @param newx Matrix of predictors 
#' @param tau Quantile of interest. Default is NULL, which will return all quantiles. Should not be specified if modelsIndex is used.  
#' @param a Tuning parameter of a. Default is NULL, which returns coefficients for all values of a. Should not be specified if modelsIndex is used. 
#' @param lambda Tuning parameter of \eqn{\lambda}. Default is NULL, which returns coefficients for all values of \eqn{\lambda}.
#' @param modelsIndex Index of the models for which coefficients should be returned. Does not need to be specified if tau or a are specified. 
#' @param lambdaIndex Index of the lambda values for which coefficients should be returned. Does not need to be specified if lambda is specified. 
#' @param ... Additional parameters passed to coef.rq.pen.seq() 
#'
#' @return A matrix of predictions for each tau and a combination
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
  coefs <- coefficients(object,tau,a,lambda,modelsIndex,lambdaIndex,...)
  #lapply(coefs, quick.predict,newx=newx)
  if(is.null(dim(newx))){
    c(1,newx) %*% coefs
  } else{
    cbind(1,newx) %*% coefs
  }
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
#' @param tauWeights Weights for the different tau models. Only used in group tau results (gtr). 
#' @param printProgress If set to TRUE prints which partition is being worked on. 
#' @param weights Weights for the quantile loss objective function.
#' @param ... Additional arguments passed to rq.pen()
#' 
#' @details 
#' Two cross validation results are returned. One that considers the best combination of a and lambda for each quantile. The second considers the best combination of the tuning 
#' parameters for all quantiles. Let \eqn{y_{b,i}}, \eqn{x_{b,i}}, and \eqn{m_{b,i}} index the response, predictors, and weights of observations in 
#' fold b. Let \eqn{\hat{\beta}_{\tau,a,\lambda}^{-b}} be the estimator for a given quantile and tuning parameters that did not use the bth fold. Let \eqn{n_b} be the number of observations in fold
#' b. Then the cross validation error for fold b is 
#' \deqn{\mbox{CV}(b,\tau) = \frac{1}{n_b} \sum_{i=1}^{n_b} m_{b,i} \rho_\tau(y_{b,i}-x_{b,i}^\top\hat{\beta}_{\tau,a,\lambda}^{-b}).}
#' Note that \eqn{\rho_\tau()} can be replaced by a different function by setting the cvFunc parameter. The function returns two different cross-validation summaries. The first is btr, by tau results. 
#' It provides the values of \code{lambda} and \code{a} that minimize the average, or whatever function is used for \code{cvSummary}, of \eqn{\mbox{CV}(b)}. In addition it provides the 
#' sparsest solution that is within one standard error of the minimum results. 
#' 
#' The other approach is the group tau results, gtr. Consider the case of estimating Q quantiles of \eqn{\tau_1,\ldots,\tau_Q} with quantile (tauWeights) of \eqn{v_q}. The gtr returns the values of \code{lambda} and \code{a} that minimizes the average, or again whatever function is used for \code{cvSummary}, of 
#' \deqn{\sum_{q=1}^Q v_q\mbox{CV}(b,\tau_q).} If only one quantile is modeled then the gtr results can be ignored as they provide the same minimum solution as btr. 
#' 
#' @return
#' An rq.pen.seq.cv object. 
#' \describe{
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
#' r5 <- rq.pen.cv(x,y,penalty.factor=c(0,rep(1,7)))
#' }
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}
rq.pen.cv <- function(x,y,tau=.5,lambda=NULL,penalty=c("LASSO","Ridge","ENet","aLASSO","SCAD","MCP"),a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,nlambda=100,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),printProgress=FALSE,weights=NULL,...){
	n <- length(y)
	if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
	} else{
    nfolds <- length(unique(foldid))
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
		trainidx <- foldid!=i
		testidx <- foldid==i
		
  		train_x <- x[trainidx,]
  		train_y <- y[trainidx]
		train_wts <- weights[trainidx]
		
  		test_x <- x[testidx,,drop=FALSE]
  		test_y <- y[testidx]
		test_wts <- weights[testidx]
  		
		trainModel <- rq.pen(train_x,train_y,tau,lambda=fit$lambda,penalty=penalty,a=fit$a,lambda.discard=FALSE,weights=train_wts,...)
  		if(is.null(cvFunc)){
  			testErrors <- check.errors(trainModel,test_x,test_y)
  		} else{
  			testErrors <- lapply(predErrors(trainModel,test_x,test_y),cvFunc)
  		}
		if(is.null(weights)==FALSE){
			testErrors <- lapply(testErrors,"*",test_wts)
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
		gtr <- groupTauResults(foldErrors, tauvals,fit$a,avals,fit$models,tauWeights,fit$lambda,stdErr)
	} else{
		indErrors <- t(indErrors)/n
		btr <- byTauResults(indErrors,tauvals,avals,fit$models,stdErr,fit$lambda)
		gtr <- groupTauResults(indErrors, tauvals,fit$a,avals,fit$models,tauWeights,fit$lambda,stdErr)
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
#' @method print rq.pen.seq.cv
print.rq.pen.seq.cv <- function(x,...){
	if(length(x$fit$tau)==1){
		cat("\nCross validation tuning parameter choices\n")
		print(x$btr)
	} else{
	  if(x$fit$penalty != "gq"){
  		cat("\nCross validation tuning parameter optimized for each quantile\n")
  		print(x$btr)
	  }
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
#'  \dontrun{
#'  set.seed(1)
#'  x <- matrix(rnorm(800),nrow=100)
#'  y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#'  lassoModels <- rq.pen.cv(x,y,tau=seq(.1,.9,.1))
#'  coefficients(lassoModels,septau=FALSE)
#'  coefficients(lassoModels,cvmin=FALSE)
#' }
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
coef.rq.pen.seq.cv <- function(object,septau=ifelse(object$fit$penalty!="gq",TRUE,FALSE),cvmin=TRUE,useDefaults=TRUE,tau=NULL,...){
  if(object$fit$penalty=="gq" & septau){
    septau <- FALSE
    warning("septau set to false because group quantile penalty was used, which is a joint optimization across all quantiles")
  }
  if(!useDefaults){
    coefficients(object$fit,tau=tau,...)
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
      nm <- length(models)
      returnVal <- matrix(0,nrow=nrow(coef(object$fit)),ncol=nm)  #vector(mode="list", length=length(models))
      rownames(returnVal) <- rownames(coef(object$fit))
      colnames(returnVal) <- names(models)
      for(i in 1:nm){
        returnVal[,i] <- coef(object$fit,modelsIndex=btr$modelsIndex[i],lambdaIndex=lambdaIndex[i])
      }
      returnVal
    } else{
	  keepers <- which(closeEnough(tau,object$gtr$tau))
	  gtr <- object$gtr[keepers,]
	  if(!cvmin){
        li <- gtr$lambda1seIndex[1]
      } else{
        li <- gtr$lambdaIndex[1]
      }
	  coef(object$fit,modelsIndex=gtr$modelsIndex,lambdaIndex=li)
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


#' Performs cross validation for a group penalty. 
#'
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
#' @param tauWeights Weights for the tau penalty only used in group tau results (gtr). 
#' @param printProgress If set to TRUE will print which fold the process is working on. 
#' @param weights Weights for the quantile loss function. Used in both model fitting and cross-validation. 
#' @param ... Additional parameters that will be sent to rq.group.pen().
#'
#' @return
#' An rq.pen.seq.cv object. 
#' \describe{
#' \item{cverr}{Matrix of cvSummary function, default is average, cross-validation error for each model, tau and a combination, and lambda.}
#' \item{cvse}{Matrix of the standard error of cverr foreach model, tau and a combination, and lambda.}
#' \item{fit}{The rq.pen.seq object fit to the full data.}
#' \item{btr}{A data.table of the values of a and lambda that are best as determined by the minimum cross validation error and the one standard error rule, which fixes a. In btr the values of lambda and a are selected seperately for each quantile.}
#' \item{gtr}{A data.table for the combination of a and lambda that minimize the cross validation error across all tau.}
#' \item{gcve}{Group, across all quantiles, cross-validation error results for each value of a and lambda.}
#' \item{call}{Original call to the function.}
#' }
#' 
#' @details 
#' Two cross validation results are returned. One that considers the best combination of a and lambda for each quantile. The second considers the best combination of the tuning 
#' parameters for all quantiles. Let \eqn{y_{b,i}}, \eqn{x_{b,i}}, and \eqn{m_{b,i}} index the response, predictors, and weights of observations in 
#' fold b. Let \eqn{\hat{\beta}_{\tau,a,\lambda}^{-b}} be the estimator for a given quantile and tuning parameters that did not use the bth fold. Let \eqn{n_b} be the number of observations in fold
#' b. Then the cross validation error for fold b is 
#' \deqn{\mbox{CV}(b,\tau) = \frac{1}{n_b} \sum_{i=1}^{n_b} m_{b,i} \rho_\tau(y_{b,i}-x_{b,i}^\top\hat{\beta}_{\tau,a,\lambda}^{-b}).}
#' Note that \eqn{\rho_\tau()} can be replaced by a different function by setting the cvFunc parameter. The function returns two different cross-validation summaries. The first is btr, by tau results. 
#' It provides the values of \code{lambda} and \code{a} that minimize the average, or whatever function is used for \code{cvSummary}, of \eqn{\mbox{CV}(b)}. In addition it provides the 
#' sparsest solution that is within one standard error of the minimum results. 
#' 
#' The other approach is the group tau results, gtr. Consider the case of estimating Q quantiles of \eqn{\tau_1,\ldots,\tau_Q} with quantile (tauWeights) of \eqn{v_q}. The gtr returns the values of \code{lambda} and \code{a} that minimizes the average, or again whatever function is used for \code{cvSummary}, of 
#' \deqn{\sum_{q=1}^Q v_q\mbox{CV}(b,\tau_q).} If only one quantile is modeled then the gtr results can be ignored as they provide the same minimum solution as btr. 
#' 
#' 
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
#' 
rq.group.pen.cv <- function(x,y,tau=.5,groups=1:ncol(x),lambda=NULL,a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),printProgress=FALSE,weights=NULL,...){
	n <- length(y)
	if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
	} else{
	  nfolds <- length(unique(foldid))    
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
		trainIdx <- foldid!=i
		testIdx <- foldid==i
		
		train_x <- x[trainIdx,]
		train_y <- y[trainIdx]
		train_wts <- weights[trainIdx]
		
		test_x <- x[testIdx,,drop=FALSE]
		test_y <- y[testIdx]
		test_wts <- weights[testIdx]
		
		trainModel <- rq.group.pen(train_x,train_y,tau,groups=groups,lambda=fit$lambda,a=fit$a,lambda.discard=FALSE,weights=train_wts,...)
		if(is.null(cvFunc)){
			testErrors <- check.errors(trainModel,test_x,test_y)
		} else{
			testErrors <- lapply(predErrors(trainModel,test_x,test_y),cvFunc)
		}
		if(is.null(weights)==FALSE){
			testErrors <- lapply(testErrors,"*",test_wts)
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
		gtr <- groupTauResults(foldErrors, tauvals,fit$a,avals,fit$models,tauWeights,fit$lambda,stdErr)
	} else{
		indErrors <- t(indErrors)/n
		btr <- byTauResults(indErrors,tauvals,avals,fit$models,stdErr,fit$lambda)
		gtr <- groupTauResults(indErrors, tauvals,fit$a,avals,fit$models,tauWeights,fit$lambda,stdErr)
	}

	returnVal <- list(cverr = foldErrors, cvse = stdErr, fit = fit, btr=btr, gtr=gtr$returnTable, gcve=gtr$gcve,call=match.call())
	class(returnVal) <- "rq.pen.seq.cv"
	returnVal
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
#' @param alg Defaults for small p to linear programming (LP), see Wang, Wu and Li (2012) for details. QICD is no longer available. 
#' @param penVars Variables that should be penalized. With default value of NULL all variables are penalized.
#' @param internal Whether call to this function has been made internally or not. 
#' @param ... Additional items to be sent to rq.lasso.fit.
#' 
#' @description 
#' Warning: this function is no longer exported. Produces penalized quantile regression models for a range of lambdas and penalty of choice. If lambda is unselected than an iterative algorithm is used to 
#' find a maximum lambda such that the penalty is large enough to produce an intercept only model. Then range of lambdas goes from the maximum lambda found to "eps" on the 
#' log scale. Local linear approximation approach used by Wang, Wu and Li to extend LLA as proposed by Zou and Li (2008) to the quantile regression setting.  
#'
#' @return Returns the following:
#' \describe{
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
#' 
#' @keywords internal
#'
#' @examples \dontrun{
#' x <- matrix(rnorm(800),nrow=100)
#' y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#' scadModel <- rq.nc.fit(x,y,lambda=1)
#' }
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
                      alg="LP",penVars=NULL,internal=FALSE,...){
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
plot.rq.pen.seq.cv <- function(x,septau=ifelse(x$fit$penalty!="gq",TRUE,FALSE),tau=NULL,logLambda=TRUE,main=NULL,...){
  if(septau & x$fit$penalty=="gq"){
    warning("septau set to false because group quantile penalty was used, which is a joint optimization across all quantiles")
    septau <- FALSE
  }
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
plot.rq.pen.seq <- function(x,vars=NULL,logLambda=TRUE,tau=NULL,a=NULL,lambda=NULL,modelsIndex=NULL,lambdaIndex=NULL,main=NULL, ...){
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
		if(is.null(vars)){
			for(i in 1:dim(betas)[1]){
				lines(lambdas[subli], betas[i,],col=i)
			}
		}else{
			for(i in vars){
				lines(lambdas[subli], betas[i,],col=i)
			}
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
#' @description Warning: this function is no longer exported.  
#'
#' @return Plot of how beta estimates change with lambda.
#' 
#' @keywords internal
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
#' @param vars Index of the variables to plot with 1 being the intercept, 2 being the first predictor, etc. Default is to include all variables. 
#' @param ... Additional parameters sent to coef()
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
bytau.plot.rq.pen.seq <- function(x,a=NULL,lambda=NULL,lambdaIndex=NULL,vars=NULL,...){
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
	lli <- length(lambdaIndex)
	if(lli>1){
		stop("Function only supports a single value of lambda or lambdaIndex")
	}
	if(lli==0){
		stop("Lambda value must be one used when fitting the models")
	}
	coefs <- coefficients(x,a=a,lambdaIndex=lambdaIndex,...)
	if(is.null(vars)){
	  pindex <- 1:nrow(coefs)
	} else{
	  pindex <- vars
	}
	lp <- length(pindex)
	if(lp > 1){	
		opar <- par(ask=TRUE) 
		on.exit(par(opar))
	}
	tau <- x$tau
	for(i in pindex){
		plot(tau,coefs[i,],xlab=expression(tau),ylab="Coefficient",main=rownames(coefs)[i],pch=16)
		lines(tau,coefs[i,])
	}
	#if(lp > 1){	par(ask=FALSE) }
}

#' Plot of coefficients varying by quantiles for rq.pen.seq.cv object
#'
#' @param x An rq.pen.seq.cv object
#' @param septau Whether optimal tuning parameters are estimated separately for each quantile.
#' @param cvmin Whether the minimum cv error should be used or the one standard error rule. 
#' @param useDefaults Set to FALSE if you want to use something besides minimum cv or 1se. 
#' @param vars Index of the variables to plot with 1 being the intercept, 2 being the first predictor, etc. Default is to include all variables. 
#' @param ... Additional parameters sent to coef() 
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
bytau.plot.rq.pen.seq.cv <- function(x,septau=ifelse(x$fit$penalty!="gq",TRUE,FALSE),cvmin=TRUE,useDefaults=TRUE,vars=NULL,...){
	coefs <- coefficients(x,septau,cvmin,useDefaults,tau=x$fit$tau,...)
	if(ncol(coefs) != length(x$fit$tau)){
		stop("Too many coefficients returned, function only works with one lambda value")
	}
	if(is.null(vars)){
	  pindex <- 1:nrow(coefs)
	} else{
	  pindex <- vars
	}
	lp <- length(pindex)
	if(lp > 1){	
		opar <- par(ask=TRUE) 
		on.exit(par(opar))
	}
	tau <- x$fit$tau
	for(i in pindex){
		plot(tau,coefs[i,],xlab=expression(tau),ylab=paste("Coefficient estimate"),main=rownames(coefs)[i],pch=16)
		lines(tau,coefs[i,])
	}
	#if(lp > 1){	par(ask=FALSE) }
}




#' Plots of cross validation results as a function of lambda. 
#'
#' @param model A cv.rq.pen() object.
#' @param logLambda Whether lambda values should be logged or not. 
#' @param loi Lambda indexes of interest, if null all lambda values will be used. 
#' @param ... Additional parameters sent to plot function.
#'
#' @return returns a cross validation plot
#'
#' @keywords internal
#'
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} 
cv_plots <- function(model,logLambda=TRUE,loi=NULL,...){
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
#' @param alg Algorithm used for fit. Only "LP", "QICD" is no longer available. 
#' @param penGroups Specify which groups will be penalized. Default is to penalize all groups.
#' @param ... Additional arguments to be sent to rq.group.fit  
#' 
#' @description 
#' This function is no longer exported. Recommend using rq.group.pen.cv() instead. 
#'
#' @return
#' Returns the following: 
#' \describe{
#' \item{beta}{ Matrix of coefficients for different values of lambda}
#' \item{residuals}{ Matrix of residuals for different values of lambda.}
#' \item{rho}{Vector of rho, unpenalized portion of the objective function, for different values of lambda.}
#' \item{cv}{ Data frame with "lambda" and second column is the evaluation based on the criteria selected.}
#' \item{lambda.min}{ Lambda which provides the smallest statistic for the selected criteria.}
#' \item{penalty}{ Penalty selected.} 
#' \item{intercept}{Whether intercept was included in model.}
#' \item{groups}{Group structure for penalty function.}
#' }
#' @keywords internal
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
cv.rq.group.pen <- function (x, y, groups, tau = 0.5, lambda = NULL, penalty = "SCAD", 
    intercept = TRUE, criteria = "CV", cvFunc = "check", nfolds = 10, 
    foldid = NULL, nlambda = 100, eps = 1e-04, init.lambda = 1,alg="huber",penGroups=NULL,
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
#' @param alg Only LP, QICD no longer available
#' @param a Additional tuning parameter for SCAD and MCP. 
#' @param penGroups Vector of TRUE and FALSE entries for each group determing if they should be penalized. Default is TRUE for all groups.
#' @param ... Additional arguments sent to rq.group.lin.prog()
#'
#' @return Returns the following:    
#' \describe{  
#' \item{coefficients}{Coefficients of the model.}
#' \item{residuals}{ Residuals from the fitted model.}
#' \item{rho}{Unpenalized portion of the objective function.}
#' \item{tau}{ Quantile being modeled.}
#' \item{n}{Sample size.}
#' \item{intercept}{Whether intercept was included in model.}
#' }
#' @description Warning: function is no longer exported. Recommend using rq.group.pen() instead. 
#' Similar to cv.rq.pen function, but uses group penalty. Group penalties use the L1 norm instead of L2 for computational convenience. 
#' As a result of this the group lasso penalty is the same as the typical lasso penalty and thus you should only use a SCAD or MCP penalty. 
#' Only the SCAD and MCP penalties incorporate the group structure into the penalty. The group lasso penalty is implemented because it is 
#' needed for the SCAD and MCP algorithm.
#' 
#' @keywords internal
#' 
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu} and Adam Maidman
#' 
#' @references 
#' \itemize{
#' \item Yuan, M. and Lin, Y. (2006). Model selection and estimation in regression with grouped variables. \emph{J. R. Statist. Soc. B}, \bold{68}, 49-67.
#' \item Peng, B. and Wang, L. (2015). An Iterative Coordinate Descent Algorithm for High-Dimensional Nonconvex Penalized Quantile Regression. \emph{Journal of Computational and Graphical Statistics}, \bold{24}, 676-694.
#' }
rq.group.fit <- function (x, y, groups, tau = 0.5, lambda, intercept = TRUE, 
                penalty = "SCAD", alg="LP", a=3.7,penGroups=NULL, ...) 
{
  #warning("recommend using rq.group.pen() instead, this is an older function with fewer options and does not include faster algorithms")
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
    
    return_val
}

#' Cross validation plot for cv.rq.group.pen object
#'
#' @param x A cv.rq.group.pen object
#' @param ... Additional parameters for plot function.
#'
#' @keywords internal
#' 
#' @return A cross validation plot. 
#'
plot.cv.rq.group.pen <- function (x,...) 
{
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
#' \describe{
#' \item{coefficients}{ Coefficients from the penalized model.} 
#' \item{PenRho}{ Penalized objective function value.}
#' \item{residuals}{ Residuals from the model.}
#' \item{rho}{ Objective function evaluation without the penalty.}
#' \item{tau}{ Conditional quantile being modeled.}
#' \item{n}{ Sample size.}  
#' }
#' @description Fits a quantile regression model with the LASSO penalty. Uses the augmented data approach similar to the proposal in Sherwood and Wang (2016).   
#' 
#' @keywords internal
#'
#' @examples \dontrun{
#' x <- matrix(rnorm(800),nrow=100)
#' y <- 1 + x[,1] - 3*x[,5] + rnorm(100)
#' lassoModel <- rq.lasso.fit(x,y,lambda=.1)
#' }
#' @references 
#' \itemize{
#' \item Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. \emph{Journal of the Royal Statistical Society. Series B}, \bold{58}, 267--288.
#' \item Wu, Y. and Liu, Y. (2009). Variable selection in quantile regression. \emph{Statistica Sinica}, \bold{19}, 801--817.  
#' \item Sherwood, B. and Wang, L. (2016) Partially linear additive quantile regression in ultra-high dimension. \emph{Annals of Statistics} \bold{44}, 288--317. 
#' }
rq.lasso.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                         coef.cutoff=1e-08, method="br",penVars=NULL,scalex=TRUE,lambda.discard=TRUE, ...){
  #warning("Recommend using rq.pen() instead, this is an older functions with fewer options and does not provide access to faster algorithms.")
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
#' @keywords internal
#' 
#' @description This function is no longer exported. 
#' 
#'
predict.rq.pen <- function(object, newx,...){
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
#' @keywords internal
#' 
#' @description This function is no longer exported. 
#' 
#'
predict.cv.rq.pen <- function(object, newx, lambda="lambda.min",...){
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
#' @keywords internal
#'
#'
#' @return Vector of coefficients. 
#'
coef.cv.rq.group.pen <- function(object, lambda='min',...){
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
#' @param alg Algorithm used. Choices are Huber approximation ("huber") or linear programming ("lp").
#' @param a The additional tuning parameter for adaptive lasso, SCAD and MCP. 
#' @param norm Whether a L1 or L2 norm is used for the grouped coefficients. 
#' @param group.pen.factor Penalty factor for each group. Default is 1 for all groups if norm=1 and square root of group size if norm=2. 
#' @param tau.penalty.factor Penalty factor for each quantile.
#' @param scalex Whether X should be centered and scaled so that the columns have mean zero and standard deviation of one. If set to TRUE, the coefficients will be returned to the original scale of the data.
#' @param coef.cutoff Coefficient cutoff where any value below this number is set to zero. Useful for the lp algorithm, which are prone to finding almost, but not quite, sparse solutions. 
#' @param max.iter The maximum number of iterations for the algorithm. 
#' @param converge.eps The convergence criteria for the algorithms. 
#' @param gamma The tuning parameter for the Huber loss. 
#' @param lambda.discard Whether lambdas should be discarded if for small values of lambda there is very little change in the solutions. 
#' @param weights Weights used in the quanitle loss objective function. 
#' @param ... Additional parameters 
#' 
#' @description  
#' Let the predictors be divided into G groups with G corresponding vectors of coefficients, \eqn{\beta_1,\ldots,\beta_G}. 
#' Let \eqn{\rho_\tau(a) = a[\tau-I(a<0)]}. Fits quantile regression models for Q quantiles by minimizing the penalized objective function of
#' \deqn{\sum_{q=1}^Q \frac{1}{n} \sum_{i=1}^n m_i \rho_\tau(y_i-x_i^\top\beta^q) + \sum_{q=1}^Q  \sum_{g=1}^G P(||\beta^q_g||_k,w_q*v_j*\lambda,a).}
#' Where \eqn{w_q} and \eqn{v_j} are designated by penalty.factor and tau.penalty.factor respectively and \eqn{m_i} can be set by weights. The value of \eqn{k} is chosen by \code{norm}.
#' Value of P() depends on the penalty. Briefly, but see references or vignette for more details,
#' \describe{
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
#' \describe{
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
#' \describe{
#' \item{coefficients}{Coefficients for each value of lambda.}
#' \item{rho}{The unpenalized objective function for each value of lambda.}
#' \item{PenRho}{The penalized objective function for each value of lambda.}
#' \item{nzero}{The number of nonzero coefficients for each value of lambda.}
#' \item{tau}{Quantile of the model.}
#' \item{a}{Value of a for the penalized loss function.}
#' }
#' @export
#'
#' @examples
#' \dontrun{ 
#' set.seed(1)
#' x <- matrix(rnorm(200*8,sd=1),ncol=8)
#' y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(200,3)
#' g <- c(1,1,1,2,2,2,3,3)
#' tvals <- c(.25,.75)
#' r1 <- rq.group.pen(x,y,groups=g)
#' r5 <- rq.group.pen(x,y,groups=g,tau=tvals)
#' #Linear programming approach with group SCAD penalty and L1-norm
#' m2 <- rq.group.pen(x,y,groups=g,alg="br",penalty="gSCAD",norm=1,a=seq(3,4))
#' # No penalty for the first group
#' m3 <- rq.group.pen(x,y,groups=g,group.pen.factor=c(0,rep(1,2)))
#' # Smaller penalty for the median
#' m4 <- rq.group.pen(x,y,groups=g,tau=c(.25,.5,.75),tau.penalty.factor=c(1,.25,1))
#' }
#' @author Ben Sherwood, \email{ben.sherwood@ku.edu}, Shaobo Li \email{shaobo.li@ku.edu} and Adam Maidman
#' @references 
#' \insertRef{qr_cd}{rqPen}
rq.group.pen <- function(x,y, tau=.5,groups=1:ncol(x), penalty=c("gLASSO","gAdLASSO","gSCAD","gMCP"),
						lambda=NULL,nlambda=100,eps=ifelse(nrow(x)<ncol(x),.05,.01),alg=c("huber","br"), 
						a=NULL, norm=2, group.pen.factor=NULL,tau.penalty.factor=rep(1,length(tau)),
						scalex=TRUE,coef.cutoff=1e-8,max.iter=5000,converge.eps=1e-4,gamma=IQR(y)/10, lambda.discard=TRUE,weights=NULL, ...){
	penalty <- match.arg(penalty)
	alg <- match.arg(alg)
	
	g <- length(unique(groups))
	if(min(apply(x,2,sd))==0){
		stop("At least one of the x predictors has a standard deviation of zero")
	}
	if(is.null(group.pen.factor)){
	  if(norm == 2){
	    group.pen.factor <- sqrt(xtabs(~groups))
	  } else{
	    group.pen.factor <- rep(1,g)
	  }
	}
	if(norm != 1 & norm != 2){
		stop("norm must be 1 or 2")
	}
	if(length(y)!=nrow(x)){
	  stop("length of x and number of rows in x are not the same")
	}
	if(is.null(weights)==FALSE){
	  if( penalty=="ENet" | penalty=="Ridge"){
	    stop("Cannot use weights with elastic net or ridge penalty. Can use it with lasso, though may be much slower than unweighted version.")
	  }
	  if(penalty=="aLASSO"){
	    warning("Weights are ignored when getting initial (Ridge) estimates for adaptive Lasso")
	  }
	  if(length(weights)!=length(y)){
	    stop("number of weights does not match number of responses")
	  }
	  if(sum(weights<=0)>0){
	    stop("all weights most be positive")
	  }
	}
	if(is.matrix(y)==TRUE){
	  y <- as.numeric(y)
	}
	if(min(group.pen.factor) < 0 | min(tau.penalty.factor) < 0){
	  stop("Penalty factors must be non-negative.")
	}
	if(sum(group.pen.factor)==0 | sum(tau.penalty.factor)==0){
	  stop("Cannot have zero for all entries of penalty factors. This would be an unpenalized model")
	}
	
	dims <- dim(x)
	n <- dims[1]
	p <- dims[2]
	if(length(groups)!=p){
	  stop("length of groups is not equal to number of columns in x")
	}	
	if(is.null(weights) & penalty=="gAdLASSO"){
		warning("Initial estimate for group adaptive lasso ignores the weights")
	}
	
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
	if(is.matrix(y)==TRUE){
		y <- as.numeric(y)
	}
	#penalty <- match.arg(penalty)
	#alg <- match.arg(alg)
	
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
		lamMax <- getLamMaxGroup(x,y,groups,tau,group.pen.factor,penalty=penalty,scalex=scalex,tau.penalty.factor=tau.penalty.factor,norm=norm,gamma=gamma,weights=weights)
		lambda <- exp(seq(log(lamMax),log(eps*lamMax),length.out=nlambda))
	}

	penalty.factor <- mapvalues(groups,seq(1,g),group.pen.factor)
	
	if(penalty == "gLASSO"){
		return_val <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma,lambda.discard=lambda.discard,weights=weights,...)
	} else{
		if(penalty == "gAdLASSO"){
			init.model <- rq.enet(x,y,tau,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,tau.penalty.factor=tau.penalty.factor,
									a=0,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,lambda.discard=lambda.discard,...)
		} else{
			if(norm == 1){
			  init.alg <- alg
				init.model <- rq.lasso(x,y,tau,alg=init.alg,lambda=lambda,penalty.factor=penalty.factor,scalex=scalex,
							tau.penalty.factor=tau.penalty.factor,max.iter=max.iter,coef.cutoff=coef.cutoff,converge.eps=converge.eps,
							gamma=gamma,lambda.discard=lambda.discard,weights,...)
			} else{
				init.model <- rq.glasso(x,y,tau,groups, lambda, group.pen.factor,scalex,tau.penalty.factor,max.iter,converge.eps,gamma,lambda.discard=lambda.discard,weights=weights,...)
			}
		}
		return_val <- rq.group.lla(init.model,x,y,groups,penalty=penalty,a=a,norm=norm,group.pen.factor=group.pen.factor,tau.penalty.factor=tau.penalty.factor,scalex=scalex,coef.cutoff=coef.cutoff,max.iter=max.iter,converge.eps=converge.eps,gamma=gamma,lambda.discard=lambda.discard,weights=weights,...)
	}
	class(return_val) <- "rq.pen.seq"
	return_val$call <- match.call()	
	return_val$lambda <- lambda
	if(lambda.discard){
		lmin <- min(sapply(return_val$models,lambdanum))
		return_val$lambda <- return_val$lambda[1:lmin]
	}
	return_val$weights <- weights
	return_val
}


