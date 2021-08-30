##
## Testing if weight methods work for different approaches
##
## Conclusion: Weights not implemented for QICD, this is a low priority. 
## 			   Weight are implemented for the group penalty with the LP algorithm

n <- 100
p <- 8
x <- matrix( rnorm(n*p), ncol=p)
wts <- 1/runif(n)
wts[which(wts>5)] <- 5

y <- 1 + x[,1] + x[,3] - x[,8] + rnorm(n)

#these get the same answer
lp_lasso <- rq.lasso.fit(x,y,weights=wts,lambda=.1)
quick_lasso <- LASSO.fit(y,x,.5,.1,TRUE,1e-08,wts)

#and the answers are different then 
lp_lasso_no_wts <- rq.lasso.fit(x,y,lambda=.1)
quick_lasso_no_wts <- LASSO.fit(y,x,.5,.1,TRUE,1e-08)

#weights not implemented for QICD
qicd_lasso_wts <- QICD(y,x,lambda=.1,weights=wts,penalty="LASSO")
qicd_lasso_no_wts <-  QICD(y,x,lambda=.1,penalty="LASSO")

#test the nc problems
lp_scad <- rq.nc.fit(x,y,weights=wts, lambda=.1, penalty="SCAD", alg="LP")
qicd_scad <- rq.nc.fit(x,y,weights=wts, lambda=.1, penalty="SCAD", alg="QICD")

#no weights, below answers are slightly different due to the algorithm. 
#can see that the QICD implementation gets the same results as when we use the weights because weights not implemented for QICD
lp_scad_no_wts   <- rq.nc.fit(x,y,lambda=.1, penalty="SCAD", alg="LP")
qicd_scad_no_wts <- rq.nc.fit(x,y, lambda=.1, penalty="SCAD", alg="QICD")


## Testing weights for group penalties
##
g <- rep(seq(1,4),each=2)

lp_glasso_wts <- rq.group.fit(x,y,lambda=.1,groups=g,weights=wts, alg="LP")
lp_glasso_no_wts <- rq.group.fit(x,y,lambda=.1,group=g, alg="LP")

lp_gscad_wts <- rq.group.fit(x,y,lambda=.1,groups=g,weights=wts, alg="LP", penalty="SCAD")
lp_gscad_no_wts <- rq.group.fit(x,y,lambda=.1,group=g, alg="LP",penalty="SCAD")

qicd_glasso_wts <- rq.group.fit(x,y,lambda=.1,groups=g,weights=wts, alg="QICD")
qicd_glasso_no_wts <- rq.group.fit(x,y,lambda=.1,group=g, alg="QICD")

qicd_gscad_wts <- rq.group.fit(x,y,lambda=.1,groups=g,weights=wts, alg="QICD", penalty="SCAD")
qicd_gscad_no_wts <- rq.group.fit(x,y,lambda=.1,group=g, alg="QICD",penalty="SCAD")
