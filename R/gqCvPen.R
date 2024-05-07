
#' Title Cross validation for consistent variable selection across multiple quantiles. 
#'
#' @param x covariate matrix. Not needed if \code{model_obj} is supplied.
#' @param y univariate response. Not needed if \code{model_obj} is supplied.
#' @param tau a sequence of tau to be modeled, must be at least of length 3. 
#' @param lambda Values of \eqn{\lambda}. Default will automatically select the \eqn{\lambda} values.
#' @param nfolds number of folds
#' @param cvFunc loss function to be evaluated for cross-validation. Supported loss functions include quantile ("rq") and squared loss("se"). Default is the quantile loss.
#' @param tauWeights weights for different quantiles in calculating the cv error. Default is equal weight.
#' @param foldid indices of pre-split testing obervations 
#' @param printProgress If set to TRUE prints which partition is being worked on. 
#' @param weights Weights for the quantile loss objective function.
#' @param ... other arguments for \code{rq.gq.pen.cv} sent to \code{rq.gq.pen}
#'
#' @return
#' An rq.pen.seq.cv object. 
#' \describe{
#' \item{cverr:}{ Matrix of cvSummary function, default is average, cross-validation error for each model, tau and a combination, and lambda.}
#' \item{cvse:}{ Matrix of the standard error of cverr foreach model, tau and a combination, and lambda.}
#' \item{fit:}{ The rq.pen.seq object fit to the full data.}
#' \item{btr:}{ Let blank, unlike rq.pen.seq.cv() or rq.group.pen.cv(), because optmizes the quantiles individually does not make sense with this penalty.}
#' \item{gtr:}{ A data.table for the combination of a and lambda that minimize the cross validation error across all tau.}
#' \item{gcve:}{ Group, across all quantiles, cross-validation error results for each value of a and lambda.}
#' \item{call:}{ Original call to the function.}
#' }
#' 
#' @details 
#' Let \eqn{y_{b,i}} and \eqn{x_{b,i}} index the observations in 
#' fold b. Let \eqn{\hat{\beta}_{\tau,a,\lambda}^{-b}} be the estimator for a given quantile and tuning parameters that did not use the bth fold. Let \eqn{n_b} be the number of observations in fold
#' b. Then the cross validation error for fold b is 
#' \deqn{\mbox{CV}(b,\tau) = \sum_{q=1}^Q \frac{1}{n_b} \sum_{i=1}^{n_b} m_{b,i}v_q \rho_\tau(y_{b,i}-x_{b,i}^\top\hat{\beta}_{\tau_q,a,\lambda}^{-b}).}
#' Where, \eqn{m_{b,i}} is the weight for the ith observation in fold b and \eqn{v_q} is a quantile specific weight. Note that \eqn{\rho_\tau()} can be replaced squared error loss. Provides results about how the average of the cross-validation error changes with \eqn{\lambda}. Uses a
#' Huber approximation in the fitting of model, as presented in Sherwood and Li (2022).
#'  
#' @export
#'
#' @examples
#' \dontrun{ 
#' n<- 200
#' p<- 10
#' X<- matrix(rnorm(n*p),n,p)
#' y<- -2+X[,1]+0.5*X[,2]-X[,3]-0.5*X[,7]+X[,8]-0.2*X[,9]+rt(n,2)
#' taus <- seq(0.1, 0.9, 0.2)
#' cvfit<- rq.gq.pen.cv(x=X, y=y, tau=taus)
#' cvCoefs <- coefficients(cvfit)
#' }
#' @references 
#' \insertRef{heteroIdQR}{rqPen}
#' 
#' \insertRef{huberGroup}{rqPen}
#' 
#' @author Shaobo Li and Ben Sherwood, \email{ben.sherwood@ku.edu} 
rq.gq.pen.cv <- function(x=NULL, y=NULL, tau=NULL, lambda=NULL, nfolds=10, cvFunc=c("rq","se"), tauWeights=NULL,  foldid=NULL, printProgress=FALSE, weights=NULL, ...){
  cvFunc <- match.arg(cvFunc)
  ## two ways to call this function
  #if(!is.null(model_obj)){
   # y<- model_obj$y
  #  x<- model_obj$x
   # tau<- model_obj$tau
  #  ntau <- length(tau)
  #  lambda<- model_obj$lambda
  #  fullmodel<- model_obj
  #}else{
  
  
  if(is.null(weights)){
    weights <- rep(1,n)
  }
  
  fullmodel<- rq.gq.pen(x=x, y=y, tau=tau,lambda=lambda, weights=weights, ...)
  lambda<- fullmodel$lambda
  tau<- fullmodel$tau
  ntau <- length(tau)
  #}
  
  if(is.null(tauWeights)){
    tauWeights <- rep(1,ntau)
  }
  # Code statement: probably a good idea, but removing for consistency with other functions
  # if(is.null(tauWeights)){
  #   tauWeights <- rep(1, ntau)/ntau
  # }else{
  #   tauWeights <- tauWeights/sum(tauWeights)
  # } 
  
  n <- length(y)
  if(is.null(foldid)){
    foldid <- sample(rep(1:nfolds, length=n))
  } else{
    nfolds <- max(foldid)
  }
  
  nlambda<- length(lambda)
  
  mse<- mqe<- matrix(NA, nlambda*ntau, nfolds)
  rownames(mse) <- rownames(mqe) <- paste("lam", rep(1:nlambda, each=ntau), "_", rep(paste("tau", tau, sep = ""), nlambda), sep = "")
  colnames(mse) <- colnames(mqe) <- paste("fold", 1:nfolds, sep = "")
  row_label <- expand.grid(tau, lambda)
  for(i in 1:nfolds){
    ind<- which(foldid==i)
    train_x<- x[-ind,]
    train_y<- y[-ind]
    test_x<- x[ind,]
    test_y<- y[ind]
    train_wts <- weights[-ind]
    test_wts <- weights[ind]
    
    train_model<- rq.gq.pen(x=train_x, y=train_y, tau=tau, lambda=lambda, lambda.discard=FALSE, weights=train_wts, ...) #,...
    pred<- predict(train_model, newx = test_x)
    #print(pred[,c(5,10,15)])
    if(cvFunc == "se"){
      se<- (test_y-pred)^2*test_wts
      mse[,i]<- apply(se,2,mean)#as.vector(do.call(c, lapply(se, apply, 2,mean))) 
    } 
    if(cvFunc == "rq"){
	    #double for loop that could be removed
      #hacky code
  	  test_err <- test_y-pred
  	  tpos <- 1
  	  for(tauval in tau){
  	    posseq <- seq(tpos,(nlambda-1)*ntau+tpos,ntau)
  	    err_seq <- seq((tpos-1)*nlambda+1,nlambda*tpos)
    		for(k in 1:nlambda){
    			mqe[posseq[k],i] <- mean(rq.loss(test_wts*test_err[,err_seq[k]],tauval))
    		}
  	    tpos <- tpos+1
  	  }
        #eq<- sapply(1:nlambda, function(xx) rq.loss.aug(rep(test_y,ntau)-as.vector(pred[,seq(xx,xx+nlambda*(ntau-1),nlambda)]), tau, n=nrow(test_x)))
        #mqe[,i]<- as.vector(apply(eq, 2,mean))
    } 
  }
  
  #med.ind <- which.min(abs(tau-0.5))
  if(cvFunc == "rq"){
    cv.mqe<- matrix(apply(mqe, 1, mean), ntau, nlambda)
    cv.mqe.wt<- apply(sapply(1:nfolds, function(xx) tapply(mqe[,xx]*rep((tauWeights), nlambda), rep(1:nlambda, each=ntau), sum)), 1, mean)
    gcv <- cv.mqe.wt
    colnames(cv.mqe) <- paste("lam", 1:nlambda, sep = "")
    rownames(cv.mqe) <- paste("tau", tau, sep = "")
    
    cv.mqe.1se<- matrix(apply(mqe, 1, sd), ntau, nlambda)
    cv.mqe.wt.1se<- apply(sapply(1:nfolds, function(xx) tapply(mqe[,xx]*rep((tauWeights), nlambda), rep(1:nlambda, each=ntau), sum)), 1, sd)
    colnames(cv.mqe.1se) <- paste("lam", 1:nlambda, sep = "")
    rownames(cv.mqe.1se) <- paste("tau", tau, sep = "")
    
    ind.eachtau.min <- apply(cv.mqe, 1, which.min) ## choose lambda based on each tau
    lambda.min <- lambda[ind.eachtau.min]
    lambda.1se <- sapply(1:ntau, function(xx) lambda[cv.mqe[xx,] <= min(cv.mqe[xx,])+cv.mqe.1se[xx, ind.eachtau.min[xx]]][1])
    cv.min <- cv.mqe[,ind.eachtau.min]
    ind.eachtau.1se <- sapply(1:length(lambda.1se), function(xx) which(lambda==lambda.1se[xx]))
    cv.1se <- cv.mqe[,ind.eachtau.1se]
    
    ## weighted cv score
    ind.lambda.min.wt <- which.min(cv.mqe.wt)
    lambda.min.wt <- lambda[ind.lambda.min.wt]
    ind.lambda.1se.wt <- which(cv.mqe.wt <= min(cv.mqe.wt)+cv.mqe.wt.1se[ind.lambda.min.wt])[1]
    lambda.1se.wt <- lambda[ind.lambda.1se.wt]
    cv.min.wt <- cv.mqe.wt[ind.lambda.min.wt]
    cv.1se.wt <- cv.mqe.wt[ind.lambda.1se.wt]
    
    output<- list(beta=fullmodel$beta, lambda=lambda, cv_all= cv.mqe, err_all=mqe, se=cv.mqe.wt.1se,
                  lambda_min=lambda.min.wt, lambda_1se=lambda.1se.wt, cv_min=cv.min.wt, cv_1se=cv.1se.wt,
                  cvup=cv.mqe.wt+cv.mqe.wt.1se, cvlo=cv.mqe.wt-cv.mqe.wt.1se, 
                  eachtau=cbind(optimalmodel=ind.eachtau.min, lambda=lambda.min, optimalmodel_1se=ind.eachtau.1se, lambda_1se=lambda.1se),
                  foldid=foldid, ntau=ntau, n.nonzero.beta=fullmodel$n.nonzero.beta, p=ncol(test_x), tau=fullmodel$tau, X=fullmodel$X,
                  me=mqe)
  }
  else{
    if(cvFunc=="se"){
      cv.mse<- matrix(apply(mse, 1, mean), ntau, nlambda)
      cv.mse.wt<- apply(sapply(1:nfolds, function(xx) tapply(mse[,xx]*rep((tauWeights), nlambda), rep(1:nlambda, each=ntau), sum)), 1, mean)
      gcv <- cv.mse.wt
      colnames(cv.mse) <- paste("lam", 1:nlambda, sep = "")
      rownames(cv.mse) <- paste("tau", tau, sep = "")
      
      cv.mse.1se<- matrix(apply(mse, 1, sd), ntau, nlambda)
      cv.mse.wt.1se<- apply(sapply(1:nfolds, function(xx) tapply(mse[,xx]*rep((tauWeights), nlambda), rep(1:nlambda, each=ntau), sum)), 1, sd)
      colnames(cv.mse.1se) <- paste("lam", 1:nlambda, sep = "")
      rownames(cv.mse.1se) <- paste("tau", tau, sep = "")
      
      ind.eachtau.min <- apply(cv.mse, 1, which.min) ## choose lambda based on each tau
      lambda.min <- lambda[ind.eachtau.min]
      lambda.1se <- sapply(1:ntau, function(xx){
        ind.eachtau.1se <- which(cv.mse[xx,] <= min(cv.mse[xx,])+cv.mse.1se[xx, ind.eachtau.min[xx]])
        max(lambda[ind.eachtau.1se])
      })
      cv.min <- cv.mse[,ind.eachtau.min]
      ind.eachtau.1se <- sapply(1:length(lambda.1se), function(xx) which(lambda==lambda.1se[xx]))
      cv.1se <- cv.mse[,ind.eachtau.1se]
      
      ## weighted cv score
      ind.lambda.min.wt <- which.min(cv.mse.wt)
      lambda.min.wt <- lambda[ind.lambda.min.wt]
      ind.lambda.1se.wt <- which(cv.mse.wt <= min(cv.mse.wt)+cv.mse.wt.1se[ind.lambda.min.wt])[1]
      lambda.1se.wt <- lambda[ind.lambda.1se.wt]
      cv.min.wt <- cv.mse.wt[ind.lambda.min.wt]
      cv.1se.wt <- cv.mse.wt[ind.lambda.1se.wt]
      
      output<- list(beta=fullmodel$beta, lambda=lambda, cv_all= cv.mse, err_all=mse, se=cv.mse.wt.1se,
                    lambda_min=lambda.min.wt, lambda_1se=lambda.1se.wt, cv_min=cv.min.wt, cv_1se=cv.1se.wt,
                    cvup=cv.mse.wt+cv.mse.wt.1se, cvlo=cv.mse.wt-cv.mse.wt.1se,
                    eachtau=cbind(optimalmodel=ind.eachtau.min, lambda=lambda.min, optimalmodel_1se=ind.eachtau.1se, lambda_1se=lambda.1se),
                    foldid=foldid, n.nonzero.beta=fullmodel$n.nonzero.beta, ntau=ntau, p=ncol(test_x), tau=fullmodel$tau, X=fullmodel$X,
                    me=mse)
    }
  }
  #output
  nz <- sapply(fullmodel$models, getNZero)
  gtr <- data.table(tau=tau, minCv = output$cv_min, lambda=output$lambda_min, lambdaIndex=ind.lambda.min.wt,
                    lambda1se=output$lambda_1se, lambda1seIndex=ind.lambda.1se.wt, a=1, cvse=output$se[ind.lambda.min.wt], 
                    modelsIndex=1:ntau, nonzero=nz[ind.lambda.min.wt], nzse=nz[ind.lambda.1se.wt])
  returnVal <- list(me=output$me,err=output$err_all,cverr=output$cv_all, cvse=output$se, fit=fullmodel, btr=NULL, gtr=gtr, gcve=matrix(gcv,ncol=nlambda), call=match.call())
  class(returnVal) <- "rq.pen.seq.cv"
  returnVal
}# end of function
############################
#' 
#' 
#' #' Title Print a summary from \code{cv.hrq_tau_glasso} 
#' #'
#' #' @param cv.fit The CV object from \code{cv.hrq_tau_glasso}
#' #'
#' #' @export
#' #'
#' print.cv.hrq_tau_glasso <- function(cv.fit){
#'   weighted <- c(which(cv.fit$lambda==cv.fit$lambda_min), cv.fit$lambda_min, 
#'                 which(cv.fit$lambda==cv.fit$lambda_1se), cv.fit$lambda_1se)
#'   out <- rbind(cv.fit$eachtau, weighted)
#'   rownames(out)[nrow(cv.fit$eachtau)+1] <- "weighted"
#'   print(out)
#' }
#' 
#' ## coefficient
#' #' Title Getting the coefficient estimates from \code{cv.hrq_tau_glasso}.
#' #'
#' #' @param cv.fit The CV object from \code{cv.hrq_tau_glasso}
#' #' @param s lambda value, or can be character either "lambda.min" or "lambda.1se". If not specified, "lambda.min" is used.
#' #'
#' #' @return coefficient estimates corresponding to the specified lambda values. 
#' #' @export
#' #'
#' coef.cv.hrq_tau_glasso<- function(cv.fit, s){
#'   lambda<- cv.fit$lambda
#'   if(missing(s)){
#'     sval<- cv.fit$lambda_min
#'   }else{
#'     if(is.numeric(s)){
#'       sval <- s
#'     }else{
#'       if(s == "lambda.min") sval<- cv.fit$lambda_min
#'       if(s == "lambda.1se") sval<- cv.fit$lambda_1se
#'     } 
#'   } 
#'   
#'   if(length(sval)==1){
#'     ind<- which.min(abs(sval-lambda))[1]
#'     beta<- matrix(cv.fit$beta[,ind], cv.fit$p+1, cv.fit$ntau, byrow = T)
#'     rownames(beta)<- c("Intercept", unique(rownames(cv.fit$beta)[-c(1:cv.fit$ntau)]))
#'     colnames(beta) <- paste("tau", cv.fit$tau, sep = "")
#'     
#'   }else{
#'     ind <- sapply(1:length(sval), function(xx) which(lambda==sval[xx]))
#'     beta <- list()
#'     for(i in 1:length(ind)){
#'       beta[[i]] <- matrix(cv.fit$beta[,ind[i]], cv.fit$p+1, cv.fit$ntau, byrow = T)
#'       rownames(beta[[i]])<- c("Intercept", unique(rownames(cv.fit$beta)[-c(1:cv.fit$ntau)]))
#'       colnames(beta[[i]]) <- paste("tau", cv.fit$tau, sep = "")
#'     }
#'     names(beta) <- paste("best_lam_tau", tau, sep = "")
#'   }
#'   return(beta)
#' }

