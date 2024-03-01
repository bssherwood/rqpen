
############################################################
#' Title Prediction for quantile regression and consistent variable selection
#'
#' @param fit model object from \code{hrq_tau_glasso}
#' @param newX new covariate matrix
#' @param s optimal value of lambda
#'
#' @return predicted response for each tau
#' @export
#'
predict.hrq_tau_glasso<- function(fit, newX, s=NULL){
  lambda<- fit$lambda
  ntau <- length(fit$tau)
  p <- ncol(newX)
  
  Xaug <- kronecker(diag(ntau), cbind(1,newX))
  groupOrder <- order(rep(1:(p+1), ntau)) 
  Xaug <- Matrix(Xaug[,groupOrder])
  
  if(length(lambda)==1){
    beta<- as.vector(fit$beta)
    if(length(beta)==(p+1)*ntau){
      pred<- Xaug%*%beta  # vector of n*ntau
      pred_tau <- cbind(rep(tau, each=nrow(newX)), pred)
      colnames(pred_tau) <- c("tau", "pred")
    }else{
      stop("newX has wrong dimension!")
    }
    
  }else{
    if(is.null(s)){
      beta<- fit$beta
      if(nrow(beta)==(p+1)*ntau){
        pred<- Xaug%*%beta  # matrix of n*ntau by nlambda
        pred_tau <- cbind(rep(tau, each=nrow(newX)), pred)
        colnames(pred_tau) <- c("tau", paste("pred_lam", 1:length(lambda), sep = ""))
      }else{
        stop("new X has wrong dimension!")
      }
    }else{
      ind <- which.min(abs(lambda-s))
      beta<- as.vector(fit$beta[,ind])
      if(length(beta)==(p+1)*ntau){
        pred<- Xaug%*%beta  # vector of n*ntau
        pred_tau <- cbind(rep(tau, each=nrow(newX)), pred)
        colnames(pred_tau) <- c("tau", "pred")
      }else{
        stop("new X has wrong dimension!")
      }
    }
  }
  
  return(pred_tau)
}


