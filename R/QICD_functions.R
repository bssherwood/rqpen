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
    if (z$flag != 0)
        warning(switch(z$flag, "Solution may be nonunique", "Premature end - possible conditioning problem in x"))
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




QICD <- function(y, x, tau=.5, lambda, intercept=TRUE, maxin=100, maxout=20,
                eps = 1e-05, penalty="SCAD", a=3.7, coef.cutoff=1e-08, 
                initial_beta=NULL, ...)
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

  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
      pentype <- as.integer(1)
  } else{
      pentype <- as.integer(2)
  }

  if( intercept ){
    beta <- initial_beta[-1]
    intval <- initial_beta[1]
  } else{
    beta <- initial_beta
    intval <- 0    
  }


  y         <- as.double(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps    <- as.double(eps)
  maxin     <- as.integer(maxin)
  lambda    <- as.double(lambda)

  i=0
  distance <- eps+1
  groupl1 <- rep(0, p)
  beta0 <- beta
  intval0 <- intval
  residuals <- as.double(y - x%*%beta - intval)

  while( (i < maxout) & (distance >= eps) ){

    out <- .C("penderiv", as.double(beta), p, a, lambda, pentype)
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
  if(intercept)
    coefficients <- c(intval, beta)

  return( coefficients )
}



QICD.nonpen <- function(y, x, z, tau=.5, lambda, intercept=TRUE, maxin=100, maxout=20,
                eps = 1e-05, penalty="SCAD", a=3.7, coef.cutoff=1e-08, 
                initial_beta=NULL, method="br", ...)
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
  cleanInputs(y, x, lambda, initial_beta, intercept, penalty, a)

  if( class(z)!="matrix" ){                                                                                    
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

  if( intercept ){
    beta <- initial_beta[-1]
    zmat <- cbind(1,z)
  } else{
    beta <- initial_beta
    zmat <- z
  }

  # Check for Singularity of X since br fortran isn't very reliable about this
  if (qr(zmat)$rank < ncol(zmat))
      stop("z is a singular matrix (make sure intercept column is not included)")

  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps       <- as.double(eps)
  maxin     <- as.integer(maxin)
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

    while( (ii < maxin) & distance.inner >= thresh ){


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
  coefficients <- c(beta, zbeta)
  if(intercept)
    coefficients <- c(zbeta[1], beta, zbeta[-1])

  return( coefficients )
}







QICD.group <- function(y, x, groups, tau = 0.5, lambda, intercept = TRUE, 
    maxin = 100, maxout=20, eps = 1e-05, penalty = "SCAD", a = 3.7, coef.cutoff = 1e-08, 
    initial_beta = NULL, ...) 
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
#a is scale parameter, the default value is 3.7 (>2 for SCAD, >1 for MCP, not used in LASSO)
#coef.cutoff is threshold for determining nonzero coefficients
#initial_beta is vector containing initial values for intercept (if included) and coefficients.
### Should be in the form (intercept, coefficients) intercept should be left out if intercept=FALSE
{
  cleanInputs(y, x, lambda, initial_beta, intercept, penalty, a)

  # Set penalty function
  if(penalty == "SCAD"){
    pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
      pentype <- as.integer(1)
  } else{
      pentype <- as.integer(2)
  }

  if( intercept ){
    beta <- initial_beta[-1]
    intval <- initial_beta[1]
  } else{
    beta <- initial_beta
    intval <- 0    
  }


  y         <- as.double(y)
  xdoub     <- as.double(x)
  p         <- as.integer( ncol(x) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  nint      <- as.integer( length(y) )
  a         <- as.double(a)
  eps    <- as.double(eps)
  maxin     <- as.integer(maxin)
  lambda    <- as.double(lambda)

  i=0
  distance <- eps+1
  groupl1 <- rep(0, p)
  beta0 <- beta
  intval0 <- intval
  residuals <- as.double(y - x%*%beta - intval)

  while( (i < maxout) & (distance >= eps) ){

    for(grps in unique(groups)){
      groupl1[groups==grps] <- sum( abs(beta[groups==grps]) )
    }
    out <- .C("penderiv", as.double(groupl1), p, a, lambda, pentype)
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
  if(intercept)
    coefficients <- c(intval, beta)

  return( coefficients )
}













#################
### Cleaning up the functions
#################
### This function just checks to make sure that the inputs to the QICD functions are appropriate
# Easier to create this function to include in each QICD function
cleanInputs <- function(y, x, lambda, initial_beta=NULL, intercept=TRUE,
                        penalty, a, ...){
  if( class(x)!="matrix" ){                                                                                    
    stop('x needs to be a matrix')
  }

  if( lambda <= 0){
    stop("lambda must be positive")
  }

  # Make sure we use br or fn method ("?rq.fit.br" or "?rq.fit.fnc")
  # if( method != "br" & method != "fn"){
  #   stop("Incorrect method.  Choose br or fn")
  # }

  if( nrow(x) != length(y) ){
    stop('length of y and rows of x do not match')
  }

  if( !is.null(initial_beta) & (length(initial_beta) != (ncol(x) + intercept)) ){
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









