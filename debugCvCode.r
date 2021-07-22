rm(list=ls(all=TRUE))
library(devtools)
devtools::unload("rqPen")
install_github("bssherwood/rqpen")
library(rqPen)


x <- matrix(rnorm(800),ncol=8)
y <- 1 + 3*x[,1] + 2*x[,3] - 5*x[,8] + rt(100,5)
tau <- .5
lambda=NULL
weights=NULL
penalty=c("LASSO","ridge","enet","aLASSO","SCAD","MCP")[1]
a=NULL
cvFunc=NULL
nfolds=10
foldid=NULL
nlambda=100
groupError=TRUE
cvSummary=mean
tauWeights=rep(1,length(tau))


#cv.rq.pen <- function(x,y,tau=.5,lambda=NULL,weights=NULL,penalty=c("LASSO","ridge","enet","aLASSO","SCAD","MCP"),a=NULL,cvFunc=NULL,nfolds=10,foldid=NULL,nlambda=100,groupError=TRUE,cvSummary=mean,tauWeights=rep(1,length(tau)),...){
#need to think about how to handle this for multi vs one tau. Also multi-a vs single a. Do the four types or something like that and then run the code
	n <- length(y)
	if(is.null(weights)==FALSE){
		stop("weights not currently implemented. Can use cv.rq.pen.old, but it supports fewer penalties and is slower.")
	}
	if(is.null(foldid)){
      foldid <- randomly_assign(n,nfolds)
    }
	fit <- rq.pen(x,y,tau,lambda=lambda,penalty=penalty,a=a)#,...)
	if(!groupError){
		indErrors <- vector(mode="list",length=nfolds)
	}
	nt <- length(tau)
	na <- length(fit$a)
	nl <- length(fit$models[[1]]$lambda)
	foldErrors <- fe2ndMoment <- matrix(rep(0,nt*na*nl),ncol=nl)
    for(i in 1:nfolds){
		train_x <- x[foldid!=i,]
		train_y <- y[foldid!=i]
		test_x <- x[foldid==i,,drop=FALSE]
		test_y <- y[foldid==i]
		trainModel <- rq.pen(train_x,train_y,tau,lambda=fit$lambda,penalty=penalty,a=fit$a)#,...)
		if(is.null(cvFunc)){
			testErrors <- check.errors(trainModel,train_x,train_y)
		} else{
			testErrors <- lapply(predict.errors(trainModel,test_x,test_y),cvFunc)
		}
		if(!groupError){
			indErrors[[i]] <- testErrors # will this be a problem? A list of a list? 
		}
		foldMeans <- do.call(rbind, lapply(testErrors,apply,2,cvSummary))
		foldErrors <- foldErrors + foldMeans
		fe2ndMoment <- fe2ndMoment + foldMeans^2
    }
	fe2ndMoment <- fe2ndMoment/nfolds
	foldErrors <- foldErrors/nfolds
	stdErr <- sqrt( (nfolds/(nfolds-1))*(fe2ndMoment - foldErrors^2))
	tauvals <- sapply(fit$models,modelTau)
	avals <- sapply(fit$models,modelA)
	#Then get the min by this index or something....
#}




