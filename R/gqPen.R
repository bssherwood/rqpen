#### core function solvebeta is based on Rcpp implementation

#' Title Quantile regression estimation and consistent variable selection across multiple quantiles
#'
#' @param x covariate matrix
#' @param y a univariate response variable
#' @param tau a sequence of quantiles to be modeled
#' @param lambda shrinkage parameter. Default is NULL, and the algorithm provides a solution path.
#' @param nlambda Number of lambda values to be considered. 
#' @param weights observation weights. Default is NULL, which means equal weights.
#' @param penalty.factor weights for the shrinkage parameter for each covariate. Default is equal weight.
#' @param tau.penalty.factor weights for different quantiles. Default is equal weight.
#' @param gmma tuning parameter for the Huber loss
#' @param maxIter maximum number of iteration. Default is 200.
#' @param lambda.discard Default is TRUE, meaning that the solution path stops if the relative deviance changes sufficiently small. It usually happens near the end of solution path. However, the program returns at least 70 models along the solution path. 
#' @param epsilon The epsilon level convergence. Default is 1e-4.
#' @param beta0 Initial estimates. Default is NULL, and the algorithm starts with the intercepts being the quantiles of response variable and other coefficients being zeros.
#'
#' @description  
#' Uses the group lasso penalty across the quantiles to provide consistent selection across all, Q, modeled quantiles. Let \eqn{\beta^q}
#' be the coefficients for the $q$th quantiles, \eqn{\beta_j} be the Q-dimensional vector of the jth coefficient for each quantile, and
#' \eqn{\rho_\tau(u)} is the quantile loss function. The method minimizes
#' \deqn{\sum_{q=1}^Q \frac{1}{n} \sum_{i=1}^n \rho_\tau(y_i-x_i^\top\beta^q) + \lambda \sum_{j=1}^p\sqrt{Q} ||\beta_j||_2  .}
#'
#' @return returns a matrix of estimated coefficients in the sparse matrix format. 
#' Each column corresponds to a lambda value, and is of length (ntau*(p+1)) 
#' which can be restructured to a coefficient matrix for each tau and each covariate. 
#' Returned values also include the sequence of lambda, the null deviance, 
#' values of penalized loss, and unpenalized loss across the sequence of lambda. 
#' \item{beta}{The estimated coefficients for all taus and a sequence of lambdas, stored in sparse matrix format, where each column corresponds to a lambda.}
#' \item{lambda}{The sequence of lambdas.}
#' \item{null.dev}{The null deviance.}
#' \item{pen.loss}{The value of penalized loss for each lambda.}
#' \item{loss}{The value of unpenalized loss for each lambda.}
#' \item{index.nonzero.beta}{Indices stored in a p*nlambda matrix, indicating if the coefficient is zero.}
#' \item{n.nonzero.beta}{The number of nonzero coefficients for each lambda.}
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
#' fit<- rq.gq.pen(X, y, taus)
#' matrix(fit$beta[,13], p+1, length(taus), byrow=TRUE)
#' }
#' 
rq.gq.pen <- function(x, y, tau, lambda=NULL, nlambda=100, weights=NULL, penalty.factor=NULL, tau.penalty.factor=NULL, gmma=0.2, 
                          maxIter=200, lambda.discard=TRUE, epsilon=1e-4, beta0=NULL){
  
  ## basic info about dimensions
  ntau <- length(tau)
  np<- dim(x)
  n<- np[1]; p<- np[2]
  #ng<- p
  nng<- rep(ntau, p)
  
  ## some initial checks
  if(ntau < 3){
    stop("please provide at least three tau values!")
  }
  if(is.null(penalty.factor)) penalty.factor<- sqrt(nng)
  if(is.null(weights)) weights<- rep(1, n)
  if(is.null(tau.penalty.factor)) tau.penalty.factor<- rep(1, ntau)
  
  # ## standardize X
  # if(scalex){
  #   x <- scale(x)
  #   mu_x <- attributes(x)$`scaled:center`
  #   sigma_x <- attributes(x)$`scaled:scale`
  # }
  
  ## augmenting data
  Xaug <- kronecker(diag(ntau), x)
  groupIndex <- rep(1:p, ntau)
  ## reorder the group for Xaug, put the same group together, e.g. (1,1,2,2,2,3,3)
  groupOrder <- order(groupIndex) 
  Xaug <- Matrix(Xaug[,groupOrder])
  groupIndex <- sort(groupIndex)
  
  Yaug <- rep(y, ntau)
  weights_aug <- rep(weights, ntau)
  tau_aug <- rep(tau, each=n)
  
  ## initial value
  gmma0<- gmma
  if(is.null(beta0)){
    # intercept 
    b.int<- quantile(y, probs = tau)
    r<- Yaug-rep(b.int, each=n)
    # null devience
    dev0<- sum(weights_aug*rq.loss.aug(r,tau,n))
    gmma.max<- 4; gmma<- min(gmma.max, max(gmma0, quantile(abs(r), probs = 0.1)))
    beta0<- c(t(rbind(b.int, matrix(0,p, ntau)))) 
  } else{
    r <- Yaug - c(cbind(1, x)%*%beta0)   # beta0 must be of the dimension (p+1)*ntau
    beta0 <- c(t(beta0))   
    dev0<- sum(rq.loss.aug(r,tau,n))
    gmma.max<- 4; gmma<- min(gmma.max, max(gmma0, quantile(abs(r), probs = 0.1)))
  }
  r0<- r
  
  ## get sequence of lambda if not supplied
  # l2norm of gradient for each group
  #grad_kR<- -neg.gradient.aug(r=r0, weights=weights_aug, tau=tau, gmma=gmma, x=Xaug, n=n, apprx=apprx)
  grad_k<- -neg_gradient_aug(r0, weights_aug, tau_aug, gmma, Xaug, ntau) ## Rcpp
  #grad_k.norm<- tapply(grad_k, groupIndex, l2norm)
  grad_k.norm<- tapply(grad_k, groupIndex, weighted_norm, normweights=tau.penalty.factor)
  
  lambda.max<- max(grad_k.norm/penalty.factor)
  lambda.flag<- 0
  if(is.null(lambda)){
    lambda.min<- ifelse(n>p, lambda.max*0.001, lambda.max*0.01)
    #lambda<- seq(lambda.max, lambda.min, length.out = 100)
    lambda<- exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda+1))
  }else{
    # user supplied lambda
    lambda.discard<- FALSE
	#kind of a hack here, but not sure why we have the following code. 
    if(lambda.max/max(lambda) > 1.01){
      lambda<- exp(c(log(lambda.max), log(sort(lambda, decreasing = TRUE))))
    }else{
      if(length(lambda)>1 & min(lambda)<lambda.max){
        lambda.flag<- 1
        lambda.user<- lambda
        lambda<- exp(c(log(lambda.max), log(sort(lambda[lambda<lambda.max], decreasing = TRUE))))
      }else{
        #warning("lambda is too large, all coefficients are shrunk to 0!")
        betaEst <- rbind(b.int, matrix(0,p, ntau))
        rownames(betaEst) <- c("Intercept", paste("V", 1:p, sep = ""))
        colnames(betaEst) <- paste("tau_", tau, sep = "")
        #return(list(beta=betaEst))
      }
    }
  }
  
  ## QM condition in Yang and Zou, Lemma 1 (2) -- PD matrix H
  # H<- 2*t(Xaug)%*%diag(weights_aug)%*%Xaug/(n*gmma)
  # # get eigen values of sub matrices for each group
  # eigen.sub.H<- rep(0, ng)
  # for(k in 1:ng){
  #   ind<- which(group.index==k)
  #   sub.H<- H[ind, ind]
  #   eigen.sub.H[k]<- max(eigen(sub.H)$value)+1e-6
  # }
  
  ## obtain maximum eigenvalues of submatrix for each group.
  ##  we do not need to compute eigenvalues because for each group, X^TWX is diagonal and all elements are the same.
  eigen.sub.H <- colSums(weights*x^2)/(n*gmma)
  
  
  ### for loop over lambda's
  if(length(lambda)==2){
    lambda<- lambda[2]
    update<- solvebetaCpp(Xaug, Yaug, n, tau_aug, gmma, weights_aug, groupIndex, lambda, penalty.factor, 
                          tau.penalty.factor, eigen.sub.H, beta0, maxIter, epsilon, ntau)
    beta1<- update$beta_update
    iter_num <- update$n_iter
    converge_status <- update$converge
    update.r <- as.vector(update$resid)
    
    # if(scalex){
    #   beta1 <- transform_coefs.aug(beta1, mu_x, sigma_x, ntau)
    #   X <- x*matrix(sigma_x,n,p,byrow = T)+matrix(mu_x,n,p,byrow = T)
    # }else{
    #   beta1 <- beta1
    #   X <- x
    # }
    
    dev1<- sum(weights_aug*rq.loss.aug(update.r, tau, n))
    
    pen.loss<- dev1/n+lambda*sum(eigen.sub.H*sapply(1:p, function(xx) l2norm(beta1[-(1:ntau)][groupIndex==xx])))
    group.index.out<- unique(groupIndex[beta1[-(1:ntau)]!=0])
    betaEst <- matrix(beta1, (p+1), ntau, byrow = T)
    
    colnames(betaEst) <- paste("tau_", tau, sep = "")
    if(is.null(colnames(x))){
      rownames(betaEst)<- c("Intercept", paste("V", 1:p, sep = ""))
    }else{
      rownames(betaEst)<- c("Intercept", colnames(x))
    }
    
    output<- list(beta=betaEst, lambda=lambda, null.dev=dev0, pen.loss=pen.loss, loss=dev1/n, tau=tau, 
                  n.nonzero.beta=length(group.index.out), index.nonzero.beta=group.index.out, X=x, y=y)
    #output.hide<- list(converge=update$converge, iter=update$iter, rel_dev=dev1/dev0)
    #class(output)<- "hrq_tau_glasso"
    #return(output)
    
  }else{
    
    group.index.out<- matrix(0, p, length(lambda))
    n.grp<- rep(0, length(lambda)); n.grp[1]<- 0
    beta.all<- matrix(0, (p+1)*ntau, length(lambda))
    beta.all[,1]<- beta0
    kkt_seq<- rep(NA, length(lambda)); kkt_seq[1]<- TRUE  # for internal use
    converge<- rep(0, length(lambda)); converge[1]<- TRUE # for internal use
    iter<- rep(0, length(lambda)); iter[1]<- 1
    outer_iter_count <- rep(0,length(lambda)); outer_iter_count[1] <- 0 # for internal use
    rel_dev<- rep(0, length(lambda)); rel_dev[1]<- 1
    loss<- rep(0, length(lambda)); loss[1]<- dev0/n
    pen.loss<- rep(0, length(lambda)); pen.loss[1]<- dev0/n
    gmma.seq<- rep(0, length(lambda)); gmma.seq[1]<- gmma
    active.ind<- NULL
    for(j in 2:length(lambda)){ #
      #print(j)
      if(length(active.ind)<p){
        ## use strong rule to determine active group at (i+1) (pre-screening)
        grad_k<- -neg_gradient_aug(r0, weights_aug, tau_aug, gmma, Xaug, ntau)
        grad_k.norm<- tapply(grad_k, groupIndex, l2norm)
        active.ind<- which(grad_k.norm>=penalty.factor*(2*lambda[j]-lambda[j-1])) 
        n.active.ind<- length(active.ind)
        
        if(length(active.ind)==p){ # first time strong rule suggests length(active.ind)=p
          # next
          outer_iter<- 0
          kkt_met<- NA
          max_iter<- 50

          # update beta and residuals
          update<- solvebetaCpp(Xaug, Yaug, n, tau_aug, gmma, weights_aug, groupIndex, 
                                lambda[j], penalty.factor, tau.penalty.factor, eigen.sub.H, beta0, maxIter, epsilon, ntau)
          beta0<- update$beta_update
          update.iter <- update$n_iter
          update.converge <- update$converge
          update.r <- as.vector(update$resid)

        }
        
        if(length(active.ind)==0){
          inactive.ind<- 1:p
          outer_iter<- 0
          kkt_met<- NA
          update.iter<- 0
          update.converge<- NA
          update.r<- r0
        }
        else{
          #if(length(active.ind)>ng/2) max_iter<- 50
          
          inactive.ind<- (1:p)[-active.ind]
          
          ## outer loop to update beta and check KKT after each update
          kkt_met <- FALSE
          outer_iter <- 0 
          
          while(!kkt_met & length(inactive.ind)>0){
            outer_iter<- outer_iter+1
            
            # reduced data
            reduced.ind<- which(groupIndex %in% active.ind)
            reduced.group.index<- groupIndex[reduced.ind]
            u.reduced.group.index <- unique(reduced.group.index)
            reduced.ng<- length(u.reduced.group.index)
            x.sub<- Xaug[,reduced.ind]
            reduced.eigen <- eigen.sub.H[u.reduced.group.index]
            
            # update beta and residuals
            update<- solvebetaCpp(x.sub, Yaug, n, tau_aug, gmma, weights_aug, reduced.group.index, 
                                  lambda[j], penalty.factor, tau.penalty.factor, eigen.sub.H, beta0[c(1:ntau, reduced.ind+ntau)],
                                  maxIter, epsilon, ntau)
            beta.update <- update$beta_update
            update.r <- as.vector(update$resid)
            update.converge <- update$converge
            update.iter <- update$n_iter
            update.grad<- update$gradient
            beta0[c(1:ntau, reduced.ind+ntau)]<- beta.update
            
            ## check inactive set by KKT condition
            kkt_results <- kkt_check_aug(r=update.r,weights=weights_aug,w=penalty.factor,gmma=gmma,tau=tau_aug,group.index=groupIndex,
                                         inactive.ind=inactive.ind,lambda=lambda[j],x=Xaug,ntau=ntau)
            if(kkt_results$kkt==TRUE){
              kkt_met <- TRUE
            } else{
              active.ind <- c(active.ind, kkt_results$new_groups)
              inactive.ind<- (1:p)[-active.ind]
              #break
            }
          }
        }
      }
      else{ # nonsparse estimates, length(active_ind)=ng
        #print("nonsparse solution")
        outer_iter<- 0
        kkt_met<- NA
        max_iter<- 50
        
        # update beta and residuals
        update<- solvebetaCpp(x=Xaug, y=Yaug, n=n, tau=tau_aug, gmma=gmma, weights=weights_aug, groupIndex=groupIndex, 
                              lambdaj=lambda[j], wlambda=penalty.factor, wtau=tau.penalty.factor, eigenval=eigen.sub.H, betaini=beta0, 
                              maxIter=maxIter, epsilon=epsilon, ntau=ntau)
        beta0 <- update$beta_update
        update.r <- as.vector(update$resid)
        update.converge <- update$converge
        update.iter <- update$n_iter
        update.grad<- update$grad

      }
      
      outer_iter_count[j] <- outer_iter
      kkt_seq[j]<- kkt_met
      converge[j]<- update.converge
      iter[j]<- update.iter
      gmma.seq[j]<- gmma
      beta.all[,j]<- beta0
      r0<- update.r 
      
      gmma.max<- 4; gmma<- min(gmma.max, max(gmma0, quantile(abs(r0), probs = 0.1)))
      
      
      # group index and number of groups
      grp.ind<- unique(groupIndex[beta0[-(1:ntau)]!=0])
      group.index.out[grp.ind,j]<- grp.ind
      n.grp[j]<- length(grp.ind)
      
      # alternative deviance
      dev1<- sum(weights_aug*rq.loss.aug(r0, tau, n))
      loss[j]<- dev1/n
      pen.loss[j]<- dev1/n+lambda[j]*sum(eigen.sub.H*sapply(1:p, function(xx) l2norm(beta0[-(1:ntau)][groupIndex==xx])))
      rel_dev[j]<- dev1/dev0
      rel_dev_change<- rel_dev[j]-rel_dev[j-1]
      if(abs(rel_dev_change)<1e-3 & j>70 & lambda.discard) break
      
    } # end of for loop of lambda
    
    
    stop.ind<- which(rel_dev!=0)
    if(length(stop.ind)==0) stop.ind<- 1:length(lambda)
    
    if(lambda.flag==0){
      stop.ind<- stop.ind[-1]
      # if(scalex){
      #   beta.final<- apply(beta.all[,stop.ind], 2, transform_coefs.aug, mu_x,sigma_x, ntau)
      #   X <- x*matrix(sigma_x,n,p,byrow = T)+matrix(mu_x,n,p,byrow = T)
      # } else{
      #   beta.final <- beta.all[,stop.ind]
      #   X <- x
      # }
      beta.final <- beta.all[,stop.ind]
      
      if(is.null(colnames(x))){
        rownames(beta.final)<- c(paste("Intercept(tau=", tau,")", sep = ""), rep(paste("V", 1:p, sep = ""), each=ntau))
      }else{
        rownames(beta.final)<- c(paste("Intercept(tau=", tau,")", sep = ""), rep(colnames(x), each=ntau))
      }
      
      output<- list(beta=Matrix(beta.final), lambda=lambda[stop.ind], null.dev=dev0, pen.loss=pen.loss[stop.ind], 
                    loss=loss[stop.ind], n.nonzero.beta=n.grp[stop.ind], index.nonzero.beta=Matrix(group.index.out[,stop.ind]),
                    tau=tau, X=x, y=y) 
      
      #output_hide<- list(iter=iter[stop.ind], rel_dev=rel_dev[stop.ind], outer.iter=outer_iter_count[stop.ind], 
      #                   kkt=kkt_seq[stop.ind], gmma=gmma.seq[stop.ind])
      #output1<- c(output, output_hide)
      #class(output) <- "hrq_tau_glasso"
      #return(output)
    }
    
    ####################
    if(lambda.flag==1){
      length.diff<- length(lambda.user)-length(lambda)+1
      beta.all<- cbind(matrix(beta.all[,1], nrow = ntau*(p+1), ncol = length.diff),beta.all[,-1])
      # if(scalex){
      #   beta.final<- apply(beta.all, 2, transform_coefs.aug, mu_x,sigma_x,ntau)
      #   X <- x*matrix(sigma_x,n,p,byrow = T)+matrix(mu_x,n,p,byrow = T)
      # } else{
      #   beta.final <- beta.all
      #   X <- x
      # }
      beta.final <- beta.all
      
      if(is.null(colnames(x))){
        rownames(beta.final)<- c(paste("Intercept(tau=", tau,")", sep = ""), rep(paste("V", 1:p, sep = ""), each=ntau))
      }else{
        rownames(beta.final)<- c(paste("Intercept(tau=", tau,")", sep = ""), rep(colnames(x), each=ntau))
      }
      
      output<- list(beta=Matrix(beta.final), lambda=lambda.user, null.dev=dev0, pen.loss=c(rep(pen.loss[1], length.diff), pen.loss[-1]), 
                    loss=c(rep(loss[1], length.diff), loss[-1]), n.nonzero.beta=c(rep(n.grp[1], length.diff), n.grp[-1]), 
                    index.nonzero.beta=Matrix(cbind(matrix(group.index.out[,1], nrow = nrow(group.index.out), ncol = length.diff),group.index.out[,-1])), 
                    tau=tau, X=x, y=y)
      #output_hide<- list(iter=c(rep(iter[1], length.diff), iter), rel_dev=c(rep(rel_dev[1], length.diff), rel_dev), outer.iter=outer_iter_count[stop.ind], 
                       #  kkt=c(rep(kkt_seq[1], length.diff), kkt_seq[-1]), gmma=c(rep(gmma.seq[1], length.diff), gmma.seq[-1]))
      
      #class(output) <- "hrq_tau_glasso"
      warning(paste("first ", length.diff, " lambdas results in pure sparse estimates!", sep = ""))
      #return(output)
      
    }
    
  } # end of else condition
  models <- list()
  modelsInfo <- modelnames <-  NULL
  for(i in 1:ntau){
    coefs <- getGQCoefs(output$beta,i,p+1,ntau)
    models[[i]] <- rq.pen.modelreturn(coefs,x,y,tau[i],lambda,penalty.factor*tau.penalty.factor[i],"gq",1,weights=weights)
    modelnames <- c(modelnames,paste0("tau",tau[i],"a",1))
    modelsInfo <- rbind(modelsInfo,c(i,1,tau[i]))
  }
  rownames(modelsInfo) <- modelnames
  modelsInfo <- data.frame(modelsInfo)
  colnames(modelsInfo) <- c("modelIndex","a","tau")
  names(models) <- modelnames
  returnVal <- list(models=models, n=n, p=p,alg="huber",tau=tau, a=1,modelsInfo=modelsInfo,lambda=lambda[1:ncol(coefs)],penalty="gq",call=match.call())
  class(returnVal) <- "rq.pen.seq"
  returnVal
} # end of function
########################################


