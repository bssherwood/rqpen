###写一个函数输入原始的X矩阵，输出tilde_X
##X:n*p;tilde_X:n*np;
trans_data = function(X){
  n = nrow(X)
  p = ncol(X)
  #tilde_X部分
  tilde_X = matrix(data = rep(0,n*n*p),nrow = n,ncol = n*p)
  for (i in 1:n){
    tilde_X[i,1:(i*p)] = rep(X[i,],i)
  }
  return(tilde_X)
}



###用lasso和group lasso求解参数theta部分
library(glmnet)
library(gglasso)
library(scales)
library(grplasso)
library(grpreg)
library(rqPen)
library(sparsegl)
library(robustlm)



N = 100
x0 = rep(1,N)
x1 = rnorm(N,0,1)
epsi = rnorm(N,0,1)


X= cbind(x0,x1)
beta = c(0,1,2.4,-6,-1.1,2,0.5,0)
Y = c(beta[1]+beta[2]*x1[1:19] + epsi[1:19],
      beta[3]+beta[4]*x1[20:49] +epsi[20:49],
      beta[5]+beta[6]*x1[50:69] +epsi[50:69],
      beta[7]+beta[8]*x1[70:100] +epsi[70:100])
X_trans1 = trans_data(X)


#############LAD+组LASSO################
p = 2
gr1 = rep(1:N,rep(p,N))
cv_model <- cv.grpreg(X_trans1, Y, group = gr1, penalty = "grLasso", tau = 0.5, nfolds = 5)
optimal_lambda <- cv_model$lambda.min
model=rq.group.pen(X_trans1, Y, groups = gr1, penalty = "gLASSO", lambda = optimal_lambda,tau = 0.5)


##############lambda#################
# 定义 lambda 的范围和数量  nlambda <- 100  # 选择合适的数量
lamMax <- 0.1  # 选择合适的最大值
eps <- 0.001  # 选择合适的比例

lambda_seq <- seq(from = lamMax, to = eps * lamMax, length.out = nlambda)
lambda_seq <- lambda_seq[!is.na(lambda_seq)]  # 排除缺失值

# 然后继续你的代码
cv_errors <- rep(0, length(lambda_seq))
for (i in 1:length(lambda_seq)) {
  fit <- rq.group.pen(X_trans1, Y, tau = 0.5, penalty = 'gLASSO', groups = gr1, lambda = lambda_seq[i])
  cv_errors[i] <- fit[['cvm']]
}


# 寻找最小交叉验证误差对应的lambda值
best_lambda <- lambda_seq[which.min(cv_errors)]
print(best_lambda)

# 使用最佳lambda值拟合最终模型
Y=matrix(Y,ncol=1)
fit_final <- rq.group.pen(X_trans1, Y, tau = 0.5, penalty = 'gLASSO', groups = gr1, lambda =0.1)