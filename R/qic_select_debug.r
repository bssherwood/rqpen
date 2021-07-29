rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)
library(hrqglas)


library(hqreg)
library(glmnet)

set.seed(1)

x <- matrix(rnorm(100*8,sd=10),ncol=8)

y <- 1 + x[,1] + 3*x[,3] - x[,8] + rt(100,3)
g <- c(1,1,1,1,2,2,3,3)
tvals <- c(.25,.75)

obj <- rq.pen(x,y,tau=tvals,penalty="SCAD",a=c(3,4,5))

septau <- FALSE
method <- "BIC"
weights <- NULL


	if(class(obj) == "cv.rq.pen.seq"){
		obj <- obj$fit
	} else if(class(obj) != "rq.pen.seq"){
		stop("obj must be of class rq.pen.seq or cv.rq.pen.seq")
	}
	if(is.null(weights)==FALSE & septau){
		warning("Weights are only used when septau is set to true.")
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
	nl <- length(obj$models[[1]]$lambda)
	
	qic_vals <- sapply(obj$models,qic,n,method)
	if(septau){
		minQIC <- apply(qic_vals,2,min)
		lambdaIndex <- apply(qic_vals,2,which.min)
		modelsInfo <- cbind(obj$modelsInfo,minQIC,lambdaIndex)
		modelsInfo <- data.table(modelsInfo)
		modelsInfo <- modelsInfo[, .SD[which.min(minQIC)],by=tau]
		
		coefs <- vector(mode="list", length=nt)
		for(i in 1:nt){
			coefs[[i]] <- coef(obj$models[[modelsInfo$modelIndex[i]]])[,modelsInfo$lambdaIndex[i]]
		}
	} else{
		gic <- matrix(rep(0,na*nl),ncol=nl)
		tqic_vals <- t(qic_vals)
		for(i in 1:na){
			subIC <- subset(tqic_vals, obj$modelsInfo$a==obj$a[i])
			gic[i,] <- weights %*% subIC
		}
		minIndex <- which(gic==min(gic),arr.ind=TRUE)
		returnA <- obj$a[minIndex[1]]
		modelsInfo <- subset(obj$modelsInfo, a==returnA)
		modelsInfo <- cbind(modelsInfo,minQIC=tqic_vals[modelsInfo$modelIndex,minIndex[2]],lambdaIndex=minIndex[2])
		coefs <- coef(obj, lambdaIndex=minIndex[2], modelsIndex=modelsInfo$modelIndex)		
	}
	coefIndex <- 1:nt
	modelsInfo <- cbind(modelsInfo, coefIndex)
	
	return_val <- list(coefficients = coefs, ic=qic_vals,modelsInfo=modelsInfo)
	class(return_val) <- "qic.select"
	return_val
	
	
groupTauResults <- function(cvErr, tauvals,a,avals,models,tauWeights){
# improve note: maybe code in for loop could be improved upon by checking at each iteration if the better cv value has been found or not
	nl <- length(models[[1]]$lambda)
	na <- length(a)
	gcve <- matrix(rep(0,na*nl),ncol=nl)
	for(i in 1:na){
		subErr <- subset(cvErr, avals==a[i])
		gcve[i,] <- tauWeights %*% subErr
	}
	minIndex <- which(gcve==min(gcve),arr.ind=TRUE)
	returnA <- a[minIndex[1]]
	modelIndex <- which(avals==returnA)
	targetModels <- models[modelIndex]
	tauvals <- sapply(targetModels,modelTau)
	lambdavals <- sapply(targetModels,modelLambda,minIndex[1,2])
	nz <- sapply(targetModels, modelNz, minIndex[1,2])
	minCv <- cvErr[modelIndex,minIndex[1,2]]
	data.table(tau=tauvals,lambda=lambdavals,a=returnA,minCv=minCv,lambdaIndex=minIndex[1,2],modelsIndex=modelIndex, nonzero=nz)
}
