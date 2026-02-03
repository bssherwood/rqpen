#### core function solvebeta is based on Rcpp implementation

#' Title Quantile regression estimation and consistent variable selection across multiple quantiles
#'
#' @param x covariate matrix
#' @param y a univariate response variable
#' @param tau a sequence of quantiles to be modeled, must be of at least length 3. 
#' @param lambda shrinkage parameter. Default is NULL, and the algorithm provides a solution path.
#' @param nlambda Number of lambda values to be considered. 
#' @param eps If not pre-specified the lambda vector will be from lambda_max to lambda_max times eps
#' @param weights observation weights. Default is NULL, which means equal weights.
#' @param penalty.factor weights for the shrinkage parameter for each covariate. Default is equal weight.
#' @param scalex Whether x should be scaled before fitting the model. Coefficients are returned on the original scale. 
#' @param tau.penalty.factor weights for different quantiles. Default is equal weight.
#' @param gmma tuning parameter for the Huber loss
#' @param max.iter maximum number of iteration. Default is 200.
#' @param lambda.discard Default is TRUE, meaning that the solution path stops if the relative deviance changes sufficiently small. It usually happens near the end of solution path. However, the program returns at least 70 models along the solution path. 
#' @param converge.eps The epsilon level convergence. Default is 1e-4.
#' @param beta0 Initial estimates. Default is NULL, and the algorithm starts with the intercepts being the quantiles of response variable and other coefficients being zeros.
#'
#'
#' @description  
#' Uses the group lasso penalty across the quantiles to provide consistent selection across all, K, modeled quantiles. Let \eqn{\beta^q}
#' be the coefficients for the kth quantiles, \eqn{\beta_j} be the Q-dimensional vector of the jth coefficient for each quantile, and
#' \eqn{\rho_\tau(u)} is the quantile loss function. The method minimizes
#' \deqn{\sum_{q=1}^Q \frac{1}{n} \sum_{i=1}^n m_i \rho_\tau(y_i-x_i^\top\beta^q) + \lambda \sum_{j=1}^p ||\beta_j||_{2,w}  .}
#' Uses a Huber approximation in the fitting of model, as presented in Sherwood and Li (2022). Where,
#' \deqn{||\beta_j||_{2,w} = \sqrt{\sum_{k=1}^K w_kv_j\beta_{kj}^2},} where \eqn{w_k} is a quantile weight 
#' that can be specified by \code{tau.penalty.factor}, \eqn{v_j} is a predictor weight that can be assigned by \code{penalty.factor}, 
#' and \eqn{m_i} is an observation weight that can be set by \code{weights}. 
#'
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
#' fit<- rq.gq.pen(X, y, taus)
#' #use IC to select best model, see rq.gq.pen.cv() for a cross-validation approach
#' qfit <- qic.select(fit)
#' }
#' @references 
#' \insertRef{heteroIdQR}{rqPen}
#' 
#' \insertRef{huberGroup}{rqPen}
#' 
#' @author Shaobo Li and Ben Sherwood, \email{ben.sherwood@ku.edu} 
rq.gq.pen <- function(x, y, tau, lambda=NULL, nlambda=100,  eps = ifelse(nrow(x) < ncol(x), 0.01, 0.001),
                          weights=NULL, penalty.factor=NULL, scalex=TRUE, tau.penalty.factor=NULL, gmma=0.2, 
                          max.iter=200, lambda.discard=TRUE, converge.eps=1e-4, beta0=NULL){
  ## basic info about dimensions
  ntau <- length(tau)
  np<- dim(x)
  n<- np[1]; p<- np[2]
  #ng<- p
  nng<- rep(ntau, p)
  if(is.null(penalty.factor)) penalty.factor<- rep(1,p)
  if(is.null(weights)) weights<- rep(1, n)
  if(is.null(tau.penalty.factor)) tau.penalty.factor<- rep(1, ntau)
  
  ## some initial checks
  if(ntau < 3){
    stop("please provide at least three tau values!")
  }
  if(min(apply(x,2,sd))==0){
		stop("At least one of the x predictors has a standard deviation of zero")
  }
  if(length(tau)>length(unique(tau))){
    stop("All entries of tau should be unique")
  }
  if(is.unsorted(tau)){
    stop("Quantile values should be entered in ascending order")
  }
  
  if(length(y)!=nrow(x)){
    stop("length of y and number of rows in x are not the same")
  }
  if(is.null(weights)==FALSE){
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
  if(length(penalty.factor)!=p){
    stop("penalty.factor should be a length p vector")
  }
  if(length(tau.penalty.factor) != ntau){
    stop("length of tau.penalty.factor should be the number of quantiles being modeled")
  }
  if(sum(penalty.factor==0) >= n){
    stop("The number of unpenalized variables must be smaller than n")
  }
  ## standardize X
  if(scalex){
    x <- scale(x)
    mu_x <- attributes(x)$`scaled:center`
    sigma_x <- attributes(x)$`scaled:scale`
  }
  
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
  if(is.null(beta0) & min(penalty.factor)>0){
      # intercept 
      b.int<- quantile(y, probs = tau)
      r<- Yaug-rep(b.int, each=n)
      # null devience
      dev0<- sum(weights_aug*rq.loss.aug(r,tau,n))
      gmma.max<- 4; gmma<- min(gmma.max, max(gmma0, quantile(abs(r), probs = 0.1)))
      beta0<- c(t(rbind(b.int, matrix(0,p, ntau))))
  } else{
    if(is.null(beta0)){
      #issue is that some vars are not penalized
      beta0 <- matrix(0,nrow=p+1,ncol=ntau)
      keepvars <- which(penalty.factor == 0)
      initmodel <- rq(y ~ x[,keepvars], tau=tau)
      beta0[c(1,keepvars+1),] <- coefficients(initmodel)
    }
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
  penvars <- which(penalty.factor!=0)
  
  lambda.max<- max(grad_k.norm[penvars]/penalty.factor[penvars])
  lambda.preset <- FALSE
  if(is.null(lambda)){
    lambda.min<- lambda.max*eps #ifelse(n>p, lambda.max*0.001, lambda.max*0.01)
    #lambda<- seq(lambda.max, lambda.min, length.out = 100)
    lambda<- exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda+1))
  } else{
    lambda.preset <- TRUE
    #a hack to have this case work with the rest of the code, which assumes first lambda provides a sparse solution.
    diff <- lambda[1] - lambda[2]
    lambda <- c(lambda[1] + diff, lambda)
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
    for(j in 1:length(lambda)){ #
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
                                lambda[j], penalty.factor, tau.penalty.factor, eigen.sub.H, beta0, max.iter, converge.eps, ntau)
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
                                  max.iter, converge.eps, ntau)
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
                              maxIter=max.iter, epsilon=converge.eps, ntau=ntau)
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
	    if(j > 1){
		    rel_dev_change<- rel_dev[j]-rel_dev[j-1]
		    if(abs(rel_dev_change)<1e-3 & j>70 & lambda.discard) break
      }
    } # end of for loop of lambda
    
    
    stop.ind<- which(rel_dev!=0)
    if(lambda.preset){
      stop.ind <- stop.ind[-1]#removing the first initial beta for this group because the initial lambda was not actually part of the preset lambda values
                              #see code around setting lambda values for more details on why this is done. 
    }
    if(length(stop.ind)==0) stop.ind<- 1:length(lambda)
    
      stop.ind<- stop.ind#[-1]
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

    
  models <- list()
  modelsInfo <- modelnames <-  NULL
  for(i in 1:ntau){
    coefs <- getGQCoefs(output$beta,i,p+1,ntau)
    if(scalex){
      coefs <- apply(coefs, 2, transform_coefs,mu_x,sigma_x)
    }
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


