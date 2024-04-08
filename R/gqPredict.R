
############################################################
#' Title Prediction for quantile regression and consistent variable selection
#'
#' @param fit model object from \code{hrq_tau_glasso}
#' @param newX new covariate matrix. If not supplied, fitted value is returned.
#' @param s optimal value of lambda. If not supplied, a list of prediction is returned, where each element of the list corresponds to a lambda.
#'
#' @return predicted response for each tau
#' @export
#'
predict.hrq_tau_glasso<- function(fit, newX=NULL, s=NULL){
  if(is.null(newX)) newX <- fit$X
  lambda<- fit$lambda
  tau <- fit$tau
  ntau <- length(tau)
  p <- ncol(newX)
  n <- nrow(newX)
  
  Xaug <- kronecker(diag(ntau), cbind(1,newX))
  groupOrder <- order(rep(1:(p+1), ntau)) 
  Xaug <- Matrix(Xaug[,groupOrder])
  
  if(length(lambda)==1){
    beta<- as.vector(fit$beta)
    if(length(beta)==(p+1)*ntau){
      pred<- Xaug%*%beta  # vector of n*ntau
      # pred_tau <- cbind(rep(tau, each=nrow(newX)), pred)
      # colnames(pred_tau) <- c("tau", "pred")
      pred_tau <- matrix(pred, n, ntau)
      colnames(pred_tau) <- paste("tau", tau, sep="")
    }else{
      stop("newX has wrong dimension!")
    }
    
  }else{
    if(is.null(s)){
      beta<- fit$beta
      if(nrow(beta)==(p+1)*ntau){
        pred<- Xaug%*%beta  # matrix of n*ntau by nlambda
        # pred_tau <- cbind(rep(tau, each=nrow(newX)), pred)
        # colnames(pred_tau) <- c("tau", paste("pred_lam", 1:length(lambda), sep = ""))
        pred_tau <- list()
        for(i in 1:length(lambda)){
          pred_lambda_i <- matrix(pred[,i], n, ntau)
          colnames(pred_lambda_i) <- paste("tau", tau, sep="")
          pred_tau[[paste("lambda",i, sep = "")]] <- pred_lambda_i
        }
      }else{
        stop("new X has wrong dimension!")
      }
    }else{
      ind <- which.min(abs(lambda-s))
      beta<- as.vector(fit$beta[,ind])
      if(length(beta)==(p+1)*ntau){
        pred<- Xaug%*%beta  # vector of n*ntau
        #pred_tau <- cbind(rep(tau, each=n), pred)
        #colnames(pred_tau) <- c("tau", "pred")
        pred_tau <- matrix(pred, n, ntau)
        colnames(pred_tau) <- paste("tau", tau, sep="")
      }else{
        stop("new X has wrong dimension!")
      }
    }
  }
  
  return(pred_tau)
}


#' Title Prediction (from CV object) for quantile regression and consistent variable selection
#'
#' @param fit model object from \code{hrq_tau_glasso}
#' @param newX new covariate matrix. If not supplied, fitted value is returned.
#' @param s optimal value of lambda. If not supplied, a list of prediction is returned, where each element of the list corresponds to a lambda.
#'
#' @return predicted response for each tau
#' @export
predict.cv.hrq_tau_glasso<- function(fit, newX=NULL, s=NULL){
  
  if(is.null(newX)){
    newX <- fit$X
  }
  
  lambda<- fit$lambda
  tau <- fit$tau
  ntau <- length(tau)
  p <- ncol(newX)
  n <- nrow(newX)
  
  Xaug <- kronecker(diag(ntau), cbind(1,newX))
  groupOrder <- order(rep(1:(p+1), ntau)) 
  Xaug <- Matrix(Xaug[,groupOrder])
  
  if(length(lambda)==1){
    beta<- as.vector(fit$beta)
    if(length(beta)==(p+1)*ntau){
      pred<- Xaug%*%beta  # vector of n*ntau
      # pred_tau <- cbind(rep(tau, each=n), pred)
      # colnames(pred_tau) <- c("tau", "pred")
      pred_tau <- matrix(pred, n, ntau)
      colnames(pred_tau) <- paste("tau", tau, sep="")
    }else{
      stop("newX has wrong dimension!")
    }
    
  }else{
    if(is.null(s)){
      beta<- fit$beta
      if(nrow(beta)==(p+1)*ntau){
        pred<- Xaug%*%beta  # matrix of n*ntau by nlambda
        # pred_tau <- cbind(rep(tau, each=n), pred)
        # colnames(pred_tau) <- c("tau", paste("pred_lam", 1:length(lambda), sep = ""))
        pred_tau <- list()
        for(i in 1:length(lambda)){
          pred_lambda_i <- matrix(pred[,i], n, ntau)
          colnames(pred_lambda_i) <- paste("tau", tau, sep="")
          pred_tau[[paste("lambda",i, sep = "")]] <- pred_lambda_i
        }
        
      }else{
        stop("new X has wrong dimension!")
      }
    }else{
      ind <- which.min(abs(lambda-s))
      beta<- as.vector(fit$beta[,ind])
      if(length(beta)==(p+1)*ntau){
        pred<- Xaug%*%beta  # vector of n*ntau
        #pred_tau <- cbind(rep(tau, each=n), pred)
        #colnames(pred_tau) <- c("tau", "pred")
        pred_tau <- matrix(pred, n, ntau)
        colnames(pred_tau) <- paste("tau", tau, sep="")
      }else{
        stop("new X has wrong dimension!")
      }
    }
  }
  
  return(pred_tau)
}


