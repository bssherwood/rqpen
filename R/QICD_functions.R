### This is a stripped down version of rq that
### only estimates betas without bells and whistles uses "br" method
#### Useful for problems with sample size up to several thousand
shortrq.fit.br <- function (x, y, tau = 0.5)
{
    tol <- .Machine$double.eps^(2/3)
    eps <- tol
    big <- .Machine$double.xmax
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    ny <- NCOL(y)
    nsol <- 2
    ndsol <- 2
    # # Check for Singularity of X since br fortran isn't very reliable about this
    # if (qr(x)$rank < p)
    #     stop("Singular design matrix")
    lci1 <- FALSE
    qn <- rep(0, p)
    cutoff <- 0
 
    z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n +
        5), as.integer(p + 3), as.integer(p + 4), as.double(x),
        as.double(y), as.double(tau), as.double(tol), flag = as.integer(1),
        coef = double(p), resid = double(n), integer(n), double((n +
            5) * (p + 4)), double(n), as.integer(nsol), as.integer(ndsol),
        sol = double((p + 3) * nsol), dsol = double(n * ndsol),
        lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn),
        cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 *
            p), as.double(big), as.logical(lci1))
    # if (z$flag != 0)
    #     warning(switch(z$flag, "Solution may be nonunique", "Premature end - possible conditioning problem in x"))
    coef <- z$coef
    coef
}

### This is a stripped down version of rq that
### only estimates betas without bells and whistles uses "fn" method
#### Useful for problems with huge sample size
shortrq.fit.fnb <- function (x, y, tau = 0.5, beta = 0.99995, eps = 1e-06)
{
    n <- length(y)
    p <- ncol(x)
    rhs <- (1 - tau) * apply(x, 2, sum)
    d   <- rep(1,n)
    u   <- rep(1,n)
    wn <- rep(0,10*n)
    wn[1:n] <- (1-tau) #initial value of dual solution
    z <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
        c = as.double(-y), rhs = as.double(rhs), d = as.double(d),as.double(u),
        beta = as.double(beta), eps = as.double(eps),
        wn = as.double(wn), wp = double((p + 3) * p),
        it.count = integer(3), info = integer(1))
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    return ( coefficients )
}




QICD <- function(y, x, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                 initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                 a=3.7, scalex=TRUE, ...)
#y: response variable, length n vector
#x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
#tau is the quantile value
#lambda is the tuning parameter (numeric value > 0)
#intercept is a logical value,should intercept be fitted (default=TRUE) or set to zero (FALSE)
#maxin: maximum number of iterations for inside coordinate descent,default value is 100
#maxout: maximum number of iterations for outside MM step, default value is 20
#initial_beta: initial value for x-covariates, the default value is NULL (lasso estimates will be used)
#eps is the convergence threshold for coordinate descent and majorization minimization step
#penalty is the name of the penalty function ("SCAD", "MCP", "LASSO")
#a is scale parameter, the default value is 3.7 (>2 for SCAD, >1 for MCP, not used in LASSO)
#coef.cutoff is threshold for determining nonzero coefficients
#initial_beta is vector containing initial values for intercept (if included) and coefficients.
### Should be in the form (intercept, coefficients) intercept should be left out if intercept=FALSE
{
  cleanInputs(y, x, lambda, initial_beta, intercept, penalty, a)
  
  if(scalex){
	x <- scale(x)
	mu_x <- attributes(x)$`scaled:center`
	sigma_x <- attributes(x)$`scaled:scale`
  }

  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
      pentype <- as.integer(1)
  } else{
      pentype <- as.integer(2)
  }

  if( is.null(initial_beta) ){
    # initial_beta <- LASSO.fit(y, x, tau, lambda, intercept, coef.cutoff)
    initial_beta <- coefficients( cv.rq.pen.old(x, y, tau=tau, intercept=intercept, 
                            penalty="LASSO", criteria="BIC") )
  }

  if( intercept ){
    beta <- initial_beta[-1]
    intval <- initial_beta[1]
  } else{
    beta <- initial_beta
    intval <- 0    
  }


  n         <- length(y)
  y         <- as.double(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps    <- as.double(eps)
  maxin     <- as.integer(maxin)

  if( length(lambda) == 1 ){
  	lambda <- rep( lambda, p)
  } else if ( length(lambda) != p ){
  	stop("lambda must be of length 1 or p")
  }
  lambda    <- as.double(lambda)
  
  i=0
  distance <- eps+1
  groupl1 <- rep(0, p)
  beta0 <- beta
  intval0 <- intval
  residuals <- as.double(y - x%*%beta - intval)

  while( (i < maxout) & (distance >= eps) ){

	#lambda <- lambda*n
    out <- .C("penderiv", as.double(beta), p, a, lambda, pentype)
	#penweight <- as.double(out[[1]])
    penweight <- as.double( n*out[[1]] )

    out <- .C("QCD", xdoub, as.double(beta), as.double(intval), penweight, residuals,
                         nint, p, int, tau, eps, maxin)
    beta <- out[[2]]
    intval <- out[[3]]
    residuals <- as.double( out[[5]] )
    i <- i+1
    distance <- sqrt( sum((beta - beta0)^2) + (intval0 - intval)^2 )
    beta0 <- beta
    intval0 <- intval
  }

  if(i == maxout & distance > eps){
      warning(paste("did not converge after ", maxout, " iterations", sep=""))
  }    

  beta[ abs(beta) < coef.cutoff ] <- 0
  coefficients <- beta
  if(intercept){
    coefficients <- c(intval, beta)
  }
  if(scalex){
	coefficients <- transform_coefs(coefficients,mu_x,sigma_x,intercept)
  }
  return( coefficients )
}



QICD.nonpen <- function(y, x, z, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                 initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                 a=3.7, method="br", scalex=TRUE, ...)
#y: response variable, length n vector
#x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector.
#z is nxq matrix of bases; the coefficients for these columns will be unpenalized
#tau is the quantile value
#lambda is the tuning parameter (numeric value > 0)
#intercept is a logical value,should intercept be fitted (default=TRUE) (intercept should be included when using splines)
#maxin: maximum number of iterations for inside coordinate descent,default value is 100
#maxout: maximum number of iterations for outside MM step, default value is 20
#eps is the convergence threshold for coordinate descent and majorization minimization step
#penalty is the name of the penalty function ("SCAD", "MCP", "LASSO")
#a is scale parameter, the default value is 3.7 (>2 for SCAD, >1 for MCP, not used in LASSO)
#coef.cutoff is threshold for determining nonzero coefficients
#initial_beta is vector containing initial values for intercept (if included) and x coefficients.
### Should be in the form (intercept, coefficients) intercept should be left out if intercept=FALSE.
### The intercept should be included to be consistent with other methods, but intercept and z coefficients
### will be initialized to rq( y-x%*%beta ~ z ).
#method for quantile regression: can be "br" or "fn", see top for description.
{
  cleanInputs(y, x, lambda, initial_beta[ 1:(intercept+ncol(x)) ], intercept, penalty, a)
  if(scalex){
	x <- scale(x)
	mu_x <- attributes(x)$`scaled:center`
	sigma_x <- attributes(x)$`scaled:scale`
  }

  if(is(z,"matrix") == FALSE ){                                                                                    
    stop('z needs to be a matrix')
  }

  if( nrow(z) != length(y) ){
    stop('length of y and rows of z do not match')
  }

  if( method == "br"){
    zreg <- shortrq.fit.br
  } else if ( method == "fn" ){
    zreg <- shortrq.fit.fnb 
  } else {
    stop("Incorrect method.  Choose br or fn")
  }

  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
      pentype <- as.integer(1)
  } else{
      pentype <- as.integer(2)
  }

  zmat <- z
  if( intercept )
    zmat <- cbind(1,z)

  # Check for Singularity of X since br fortran isn't very reliable about this
  if (qr(zmat)$rank < ncol(zmat))
      stop("z is a singular matrix (make sure intercept column is not included)")  

  if( is.null(initial_beta) ){
    initial_beta <- LASSO.fit.nonpen(y, x, z, tau, lambda, intercept, coef.cutoff) 
    initial_beta <- initial_beta[ 1:(intercept+ncol(x)) ] ### Only keep the coefficients for x 
  } else {
    initial_beta <- initial_beta[ 1:(intercept+ncol(x)) ] ### Only keep the coefficients for x 
  }

 if( intercept ){
    beta <- initial_beta[-1]
    beta <- beta
  } else{
    beta <- initial_beta
  }



  n         <- length(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps       <- as.double(eps)
  maxin     <- as.integer(maxin)

  if( length(lambda) == 1 ){
  	lambda <- rep( lambda, p)
  } else if ( length(lambda) != p ){
  	stop("lambda must be of length 1 or p")
  }
  lambda    <- as.double(lambda)

  i=0
  distance <- distance.inner <-  eps+1
  beta0 <- beta00 <- beta
  
  xb <- x%*%beta
  xresiduals <- y - xb
  zbeta <- zreg(zmat, xresiduals, tau = tau)
  zbeta0 <- zbeta00 <- zbeta
  zb <- zmat%*%zbeta
  residuals <- as.double(y - xb - zb)

  while( (i < maxout) & (distance >= eps) ){

    ii=0
    out <- .C("penderiv", as.double(beta), p, a, lambda, pentype)
    penweight <- as.double( n*out[[1]] )

    while( (ii < maxin) & distance.inner >= eps ){


      out <- .C("QCD", xdoub, as.double(beta), as.double(0), penweight, residuals,
                       nint, p, as.integer(0), tau, eps, as.integer(1))
      beta <- out[[2]]

      xb <- x%*%beta
      xresiduals <- y - xb
      zbeta <- zreg(zmat, xresiduals, tau = tau)
      zb <- zmat%*%zbeta
      residuals <- as.double(y - xb - zb)

      ii <- ii+1
      distance.inner <- sqrt( sum((beta - beta00)^2) + sum((zbeta - zbeta00)^2) )
      beta00 <- beta
      zbeta00 <- zbeta
    }

    i <- i+1
    distance <- sqrt( sum((beta - beta0)^2) + sum((zbeta - zbeta0)^2) )
    beta0 <- beta
    zbeta0 <- zbeta
  }

  if(i == maxout & distance > eps){
      warning(paste("did not converge after ", maxout, " iterations", sep=""))
  }    

  beta[ abs(beta) < coef.cutoff ] <- 0
  if(scalex){
	beta <- transform_coefs(beta,mu_x,sigma_x,intercept)
  }
  
  coefficients <- c(beta, zbeta)
  if(intercept){
    coefficients <- c(zbeta[1], beta, zbeta[-1])
  }

  return( coefficients )
}







QICD.group <- function(y, x, groups, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                 initial_beta=NULL, maxin=100, maxout=20, eps = 1e-05, coef.cutoff=1e-08,  
                 a=3.7, scalex=TRUE, norm=2, ...)
#y: response variable, length n vector
#x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
#groups: numeric vector of length ncol(x) with the group number of the coefficient (can be unique)
#tau is the quantile value
#lambda is the tuning parameter (numeric value > 0)
#intercept is a logical value,should intercept be fitted (default=TRUE) or set to zero (FALSE)
#maxin: maximum number of iterations for inside coordinate descent,default value is 100
#maxout: maximum number of iterations for outside MM step, default value is 20
#initial_beta: initial value for x-covariates, the default value is NULL (lasso estimates will be used)
#eps is the convergence threshold for coordinate descent and majorization minimization step
#penalty is the name of the penalty function ("SCAD", "MCP", "LASSO")
#coef.cutoff is threshold for determining nonzero coefficients
#initial_beta is vector containing initial values for intercept (if included) and coefficients.
#a is scale parameter, the default value is 3.7 (>2 for SCAD, >1 for MCP, not used in LASSO)

### Should be in the form (intercept, coefficients) intercept should be left out if intercept=FALSE
{
  cleanInputs(y, x, lambda, initial_beta, intercept, penalty, a)
  
  if(scalex){
	x <- scale(x)
	mu_x <- attributes(x)$`scaled:center`
	sigma_x <- attributes(x)$`scaled:scale`
  }

  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
      pentype <- as.integer(1)
  } else{
      pentype <- as.integer(2)
  }

  if( is.null(initial_beta) ){
    initial_beta <- LASSO.fit(y, x, tau, lambda, intercept, coef.cutoff) 
  }

  if( intercept ){
    beta <- initial_beta[-1]
    intval <- initial_beta[1]
  } else{
    beta <- initial_beta
    intval <- 0    
  }


  n         <- length(y)
  y         <- as.double(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps    <- as.double(eps)
  maxin     <- as.integer(maxin)

  if( length(lambda) == 1 ){
  	lambda <- rep( lambda, p)
  } else if ( length(lambda) != p ){
  	stop("lambda must be of length 1 or p")
  }
  lambda    <- as.double(lambda)

  i=0
  distance <- eps+1
  groupNorm <- rep(0, p)
  beta0 <- beta
  intval0 <- intval
  residuals <- as.double(y - x%*%beta - intval)

  badValues <- FALSE

  while( (i < maxout) & (distance >= eps) ){

    for(grps in unique(groups)){
	  if(norm==2){
		groupNorm[groups==grps] <- sqrt( sum( beta[groups==grps]^2) )
	  } else{
		groupNorm[groups==grps] <- sum( abs(beta[groups==grps]) )
	  }
    }
	#lambda <- n*lambda
    out <- .C("penderiv", as.double(groupNorm), p, a, lambda, pentype)
	#penweight <- as.double(out[[1]])
    penweight <- as.double( n*out[[1]] )

    out <- .C("QCD", xdoub, as.double(beta), as.double(intval), penweight, residuals,
                         nint, p, int, tau, eps, maxin)
    beta <- out[[2]]
    intval <- out[[3]]
    residuals <- as.double( out[[5]] )
    i <- i+1
    distance <- sqrt( sum((beta - beta0)^2) + (intval0 - intval)^2 )
    beta0 <- beta
    intval0 <- intval

    if( max(abs(beta)) > 1e10 ){
      badValues <- TRUE
      break
    }
  }

  if( badValues ){
      warning("Some coefficients diverged to infinity (bad results)")    
  }

  if(i == maxout & distance > eps){
      warning(paste("did not converge after ", maxout, " iterations", sep=""))
  }    

  beta[ abs(beta) < coef.cutoff ] <- 0
  coefficients <- beta
  if(intercept){
    coefficients <- c(intval, beta)
  }
  
  if(scalex){
	coefficients <- transform_coefs(coefficients,mu_x,sigma_x,intercept)
  }
  
  return( coefficients )
}


### Can fit models for many different lambdas using warm start.  
## If no intial_beta is provided, QICD.master will find an appropriate starting value and fit the model for each lambda
QICD.master <- function(y, x, z=NULL, groups=NULL, tau=.5, lambda, intercept=TRUE, penalty="SCAD", 
                 initial_beta, maxin=100, maxout=20, 
                 eps = 1e-05, coef.cutoff=1e-08, a=3.7, ...){
  if( !is.null(z) & !is.null(groups) )
    stop("Currently not supporting group penalty and some nonpenalized variables \n  z or groups (or both) must be NULL")

  lambda <- sort( unique(lambda) )
  if( lambda[1] <= 0)
    stop("lambda must be positive")
  nlam <- length(lambda)
  p <- ncol(x)
  q <- ifelse( is.null(z), 0, ncol(z) )
  penVars <- intercept + 1:p ### Keep track of penalized coefficients

  ### Get names of output correct
  out <- matrix(0, nrow=intercept+p+q, ncol=nlam) ### Each column is coefficients for 1 lambda
  rnms <- paste("x",1:p, sep="")
  if( intercept )
    rnms <- c("(Intercept)", rnms)
  if( !is.null(z) )
    rnms <- c(rnms, paste("z",1:q, sep=""))

  rownames(out) <- rnms
  colnames(out) <- 1:nlam

  ### Set correct QICD function
  QF <- match.call()
  if( is.null(z) & is.null(groups) ){ ### QICD
    QF[[1]] <- as.name("QICD")
  } else if (is.null(z) ){            ### QICD.group
    QF[[1]] <- as.name("QICD.group")
  } else {                            ### QICD.nonpen
    QF[[1]] <- as.name("QICD.nonpen")
  }
  
  #####################################
  ### Don't do warm start (for now) ###
  #####################################
  ### Do first iteration
  QF$lambda <- lambda[1]
  # if( penalty == "LASSO" & is.null(groups) & !is.null(initial_beta) ){     ### Use initial_beta as coefficients for smallest lambda when LASSO with no group penalty
  #     out[,1] <- initial_beta
  # } else if (penalty == "LASSO" & is.null(groups) & is.null(initial_beta)) {                                        ### SCAD, MCP, group LASSO penalty use initial_beta as starting values
  #     QF$maxout <- -1 
  #     out[,1] <- eval.parent(QF)
  # } else {
  #     out[,1] <- eval.parent(QF)
  # }
  
  # QF$maxout <- maxout
  #####################################
  out[,1] <- eval.parent(QF)

  ibeta <- out[,1]
  if( nlam > 1 ){
    for( i in 2:nlam ){

      if( max( abs(ibeta[penVars]) ) == 0 ){ ### If all penalized betas are 0, then we are done 
        out[,i:nlam] <- matrix(ibeta, nrow=length(ibeta), ncol=nlam-i+1) ### Assign coefficients for remaining lambdas to current value 
        break ### End the for loop
      }

      QF$lambda <- lambda[i]
      # QF$initial_beta <- ibeta # For now, use initial_beta as starting values for all lambda
      out[,i] <- eval.parent(QF)
      ibeta <- out[,i]
    }
  }

  return(list( coefficients=out, lambda=lambda ))
}



### LASSO estimates for initial values in QICD and QICD.group functions when initial_beta = NULL
LASSO.fit <- function(y, x, tau, lambda, intercept, coef.cutoff, weights=NULL)
{
  p <- ncol(x)
  n <- nrow(x)
  ynew <- c( y, rep(0, 2*p) )
  xnew <- rbind( x, diag(n*lambda, p), -diag(n*lambda, p) )
  if( intercept )
    xnew <- cbind( c( rep(1,n), rep(rep(0, 2*p)) ), xnew )

  if( !is.null(weights) ){
  	xnew[1:n,] <- xnew[1:n,] * weights
  	ynew[1:n]  <-  ynew[1:n] * weights
  }

  #Ben: had some problems with fnb, setting all to use br for now
  #if( n + 2*p < 500 ){ ### If problem is small, use "br" method
    out <- shortrq.fit.br(xnew, ynew, tau)
  #} else {             ### Else use "fn" method
  #  out <- shortrq.fit.fnb(xnew, ynew, tau)
  #}

  out[ abs(out) < coef.cutoff ] <- 0
  return(out)
}


### LASSO estimates for initial values in QICD.nonpen function when initial_beta = NULL
LASSO.fit.nonpen <- function(y, x, z, tau, lambda, intercept, coef.cutoff, weights=NULL)
{
  p <- ncol(x)
  pxz <- ncol(z) + p
  n <- nrow(x)
  ynew <- c( y, rep(0, 2*p) )
  xz <- cbind( x, z )
  xz <- rbind( xz, diag(n*lambda, p, pxz), -diag(n*lambda, p, pxz) )
  if( intercept )
    xz <- cbind( c( rep(1,n), rep(rep(0, 2*p)) ), xz )

  if( !is.null(weights) ){
  	xz[1:n,]   <-   xz[1:n,] * weights
  	ynew[1:n]  <-  ynew[1:n] * weights
  }

  #Ben: had some problems with fnb, setting all to use br for now
  #if( n + 2*p < 500 ){ ### If problem is small, use "br" method
    out <- shortrq.fit.br(xz, ynew, tau)
  #} else {             ### Else use "fn" method
  #  out <- shortrq.fit.fnb(xz, ynew, tau)
  #}

  out <- out

  out[ abs(out) < coef.cutoff ] <- 0
  return(out)
}







#################
### Cleaning up the functions
#################
### This function just checks to make sure that the inputs to the QICD functions are appropriate
# Easier to create this function to include in each QICD function
cleanInputs <- function(y, x, lambda, initial_beta=NULL, intercept=TRUE,
                        penalty, a, ...){
  if( is(x,"matrix") == FALSE){                                                                                    
    stop('x needs to be a matrix')
  }

  if( any(lambda <= 0)){
    stop("lambda must be positive")
  }

  # Make sure we use br or fn method ("?rq.fit.br" or "?rq.fit.fnc")
  # if( method != "br" & method != "fn"){
  #   stop("Incorrect method.  Choose br or fn")
  # }

  if( nrow(x) != length(y) ){
    stop('length of y and rows of x do not match')
  }

  if( !is.null(initial_beta) & (length(initial_beta) < (ncol(x) + intercept)) ){
    stop("initial_beta must contain initial value for intercept (if TRUE) and each coefficient")
  }

  if( penalty == "SCAD" ){
    if( a <= 2 )
      stop("a must be > 2 for SCAD penalty")
  } else if ( penalty == "MCP" ){
    if( a <= 1 )
      stop("a must be > 1 for MCP penalty")
  } else {
    if( penalty != "LASSO" )
      stop("wrong penalty function")
  }

  return(NULL)
}









