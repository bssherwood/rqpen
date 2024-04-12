
#' Title Cross validation for consistent variable selection across multiple quantiles
#'
#' @param x covariate matrix. Not needed if \code{model_obj} is supplied.
#' @param y univariate response. Not needed if \code{model_obj} is supplied.
#' @param tau a sequence of tau to be modeled
#' @param lambda The sequence of lambdas.
#' @param nfolds number of folds
#' @param loss loss function to be evaluated. Supported loss functions include quantile ("rq") and squared loss("se"). Default is the quantile loss.
#' @param wt_tau_loss weights for different quantiles in calculating the cv error. Default is equal weight.
#' @param foldid indices of pre-split testing obervations 
#' @param ... other arguments for \code{gq.cv.pen}
#'
#' @return The full solution path is returned. It also returns the vector of CV score 
#' as well as the optimal lambda values in terms of min and 1se of the CV error. 
#' \item{beta}{The estimated coefficients for all lambdas, stored in sparse matrix format, where each column corresponds to a lambda.}
#' \item{lambda}{The sequence of lambdas.}
#' \item{cv_all}{An ntau*nlambda matrix of all values of CV error for each tau and each lambda.}
#' \item{err_all}{Errors for all \code{nfolds} folds evaluation. Has row of (ntau*nlambda) and column of nfolds.}
#' \item{lambda_min}{The optimal lambda that minimizes the CV error}
#' \item{lambda_1se}{The largest lambda such that CV error is within 1 standard error of the minimum CV error.}
#' \item{cv_min}{The value of CV error corresponding to \code{lambda_min}.}
#' \item{cv_1se}{The value of CV error corresponding to \code{lambda_1se}.}
#' \item{foldid}{The vector of indices for k folds split.}
#' \item{cvup}{CV error + 1 standard error}
#' \item{cvlo}{CV error - 1 standard error}
#' \item{n.nonzero.beta}{The number of selected covariates for each lambda.}
#' \item{eachtau}{A summary table of optimal lambdas with respect to each individual tau.}
#'  
#' @export
#'
#' @examples
#' \dontrun{ 
#' n<- 200
#' p<- 20
#' X<- matrix(rnorm(n*p),n,p)
#' y<- -2+X[,1]+0.5*X[,2]-X[,3]-0.5*X[,7]+X[,8]-0.2*X[,9]+rt(n,2)
#' taus <- seq(0.1, 0.9, 0.2)
#' cvfit<- rq.gq.cv.pen(x=X, y=y, tau=taus)
#' }
#' 
rq.gq.pen.cv <- function(x=NULL, y=NULL, tau=NULL, lambda=NULL, nfolds=10, loss=c("rq","se"), wt_tau_loss=NULL,  foldid=NULL, ...){
  loss <- match.arg(loss)
  ## two ways to call this function
  #if(!is.null(model_obj)){
   # y<- model_obj$y
  #  x<- model_obj$x
   # tau<- model_obj$tau
  #  ntau <- length(tau)
  #  lambda<- model_obj$lambda
  #  fullmodel<- model_obj
  #}else{
  fullmodel<- rq.gq.pen(x=x, y=y, tau=tau,lambda=lambda, ...)
  lambda<- fullmodel$lambda
  tau<- fullmodel$tau
  ntau <- length(tau)
  #}
  
  if(is.null(wt_tau_loss)){
    wt_tau_loss <- rep(1, ntau)/ntau
  }else{
    wt_tau_loss <- wt_tau_loss/sum(wt_tau_loss)
  } 
  
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
    
    train_model<- rq.gq.pen(x=train_x, y=train_y, tau=tau, lambda=lambda, lambda.discard=FALSE,...) #,...
    pred<- predict(train_model, newx = test_x)
    
    if(loss == "se"){
      se<- (test_y-pred)^2
      mse[,i]<- apply(se,2,mean)#as.vector(do.call(c, lapply(se, apply, 2,mean))) 
    } 
    if(loss == "rq"){
	  #double for loop that could be removed
	  test_err <- test_y-pred
	  pos <- 1
	  for(tauval in tau){
		for(k in 1:nlambda){
			mqe[pos,i] <- mean(rq.loss(test_err[,pos],tauval))
			pos <- pos + 1
		}
	  }
      #eq<- sapply(1:nlambda, function(xx) rq.loss.aug(rep(test_y,ntau)-as.vector(pred[,seq(xx,xx+nlambda*(ntau-1),nlambda)]), tau, n=nrow(test_x)))
      #mqe[,i]<- as.vector(apply(eq, 2,mean))
    } 
  }
  
  #med.ind <- which.min(abs(tau-0.5))
  if(loss == "rq"){
    cv.mqe<- matrix(apply(mqe, 1, mean), ntau, nlambda)
    cv.mqe.wt<- apply(sapply(1:nfolds, function(xx) tapply(mqe[,xx]*rep((wt_tau_loss), nlambda), rep(1:nlambda, each=ntau), sum)), 1, mean)
    gcv <- cv.mqe.wt
    colnames(cv.mqe) <- paste("lam", 1:nlambda, sep = "")
    rownames(cv.mqe) <- paste("tau", tau, sep = "")
    
    cv.mqe.1se<- matrix(apply(mqe, 1, sd), ntau, nlambda)
    cv.mqe.wt.1se<- apply(sapply(1:nfolds, function(xx) tapply(mqe[,xx]*rep((wt_tau_loss), nlambda), rep(1:nlambda, each=ntau), sum)), 1, sd)
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
    if(loss=="se"){
      cv.mse<- matrix(apply(mse, 1, mean), ntau, nlambda)
      cv.mse.wt<- apply(sapply(1:nfolds, function(xx) tapply(mse[,xx]*rep((wt_tau_loss), nlambda), rep(1:nlambda, each=ntau), sum)), 1, mean)
      gcv <- cv.mse.wt
      colnames(cv.mse) <- paste("lam", 1:nlambda, sep = "")
      rownames(cv.mse) <- paste("tau", tau, sep = "")
      
      cv.mse.1se<- matrix(apply(mse, 1, sd), ntau, nlambda)
      cv.mse.wt.1se<- apply(sapply(1:nfolds, function(xx) tapply(mse[,xx]*rep((wt_tau_loss), nlambda), rep(1:nlambda, each=ntau), sum)), 1, sd)
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

