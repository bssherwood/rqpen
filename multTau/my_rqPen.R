myrq.lasso.fit <- function(x,y,tau=.5,lambda=NULL,weights=NULL,intercept=TRUE,
                         coef.cutoff=1e-08, method="br",penVars=NULL,scalex=TRUE, ...){
    # x is a n x p matrix without the intercept term
    # y is a n x 1 vector
    # lambda takes values of 1 or p
    # coef.cutoff is a threshold to set to zero. 
    # Choose the method used to estimate the coefficients ("br" or "fn")
    ### According to quantreg manual and my experience, "fn" is much faster for big n
    ### The "n" can grow rapidly using lin. prog. approach  
    # penVars - variables to be penalized, doesn't work if lambda has multiple entries (Ben: I think it does though it is a little bit strange to do)
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
    
    # if penVars is not null, then only variables in penVars will get penalty of lambda
    # at the end of this loop, lambda becomes length p vector
    if(is.null(penVars) !=TRUE){
        if(length(lambda)==1){
            mult_lambda <- rep(0,p)
            mult_lambda[penVars] <- lambda
            lambda <- mult_lambda
        } else{
            lambda[-penVars] <- 0
        }
    }
    
    lambda <- lambda*n # need this to account for the fact that rq does not normalize the objective function
    
    # this section of code creates augmented matrices.
    # pen_x is 2p x p: diagonal matrix of lambda, stacked on top of diag matrix of -lambda 
    # (unless some rows are all zero, in which case those rows are dropped)
    # aug_n is no. of rows in pen_x
    # aug_x is x stacked on top pen_x. if intercept is true, we add a column of 
    # 1s (and 0s) to it.
    # aug_y is y stacked on top of zeros
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
    
    # if weights are provided, give each row of the augmented part of the matrix 
    # weight 1
    # pass to rq
    if(is.null(weights)){
        model <- quantreg::rq(aug_y ~ aug_x+0, tau=tau, method=method)
    } else{
        if(length(weights) != n){
            stop("Length of weights does not match length of y")
        }
        orig_weights <- weights
        weights <- c(weights, rep(1,aug_n))
        model <- quantreg::rq(aug_y ~ aug_x+0, tau=tau, weights=weights, method=method)
    }
    
    # get coefficients from the returned output of rq
    # coefficients with abs value < coef.cutoff are thresholded to 0.
    p_star <- p+intercept
    if (length(tau) == 1) {
        coefs <- matrix(coefficients(model)[1:p_star], nrow = 1)
    } else {
        coefs <- t(coefficients(model)[1:p_star, ])
    }
    return_val <- NULL
    return_val$coefficients <- coefs
    if(is.null(colnames(x))){
        x_names <- paste("x",1:p,sep="")
    } else{
        x_names <- colnames(x)
    }
    if(intercept) {
        x_names <- c("intercept",x_names)
    }
    colnames(return_val$coefficients) <- x_names
    return_val$coefficients[abs(return_val$coefficients) < coef.cutoff] <- 0
    
    # if we standardized the variables, we need to put the coefficients back on
    # the original scale
    if(scalex){
        #need to update for penVars
        return_val$coefficients <- mytransform_coefs(return_val$coefficients,mu_x,sigma_x,intercept)
        if(intercept){
            fits <- cbind(1,original_x) %*% t(return_val$coefficients)
        } else{
            fits <- original_x %*% t(return_val$coefficients)
        }
        return_val$residuals <- y - fits
        return_val$PenRho <- sapply(1:length(tau), function(i) 
            sum(rqPen::check(return_val$residuals[, i], tau[i]))) +
            myget_coef_pen(return_val$coefficients,lambda,intercept,penVars)	 
    } else{
        return_val$PenRho <- model$rho
        return_val$residuals <- model$residuals[1:n, ]   
    }
    
    # computation of loss (i.e. w/o penalty)
    if(is.null(weights)){   
        return_val$rho <- sapply(1:length(tau), function(i) 
            sum(rqPen::check(return_val$residuals[, i], tau[i])))
    } else{
        return_val$rho <- sapply(1:length(tau), function(i) 
            sum(orig_weights * rqPen::check(return_val$residuals[, i], tau[i])))
    }
    
    return_val$tau <- tau
    return_val$n <- n                  
    return_val$intercept <- intercept
    class(return_val) <- c("rq.pen", "rqLASSO")
    return_val
}

mytransform_coefs <- function (coefs, mu_x, sigma_x, intercept = TRUE) {
    new_coefs <- coefs
    if (intercept) {
        intercept <- coefs[, 1]
        for (j in 2:ncol(coefs)) {
            new_coefs[, j] <- coefs[, j]/sigma_x[j - 1]
            intercept <- intercept - coefs[, j] * mu_x[j - 1]/sigma_x[j - 1]
        }
        new_coefs[, 1] <- intercept
    }
    else {
        for (j in 1:ncol(coefs)) {
            new_coefs[, j] <- coefs[, j]/sigma_x[j]
        }
    }
    new_coefs
}

myget_coef_pen <- function (coefs, lambda, intercept, penVars, penalty = "LASSO") {
    if (intercept) {
        coefs <- coefs[, -1,  drop = F]
    }
    if (is.null(penVars) == FALSE) {
        coefs <- coefs[, penVars,  drop = F]
    }
    if (penalty == "LASSO") {
        rowSums(abs(coefs)) * lambda
    }
}