//#include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericVector stlSort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
NumericVector findIndices(NumericVector x, int k) {
  NumericVector indices;
  
  for (int i = 0; i < x.size(); i++) {
    if (x[i] == k) {
      indices.push_back(i + 1);
    }
  }
  
  return indices;
}

/* Quantile loss function for augmented data (tau needs to be vectorized) */
// [[Rcpp::export]]
NumericVector rqLossAug(NumericVector r, NumericVector tau) {
  int n = r.size();
  NumericVector val(n);
  
  for (int i = 0; i < r.size(); i++) {
    val[i] = (abs(r[i]) + (2 * tau[i] - 1) * r[i]) / 2;
  }
  
  return val;
}


// /* tanh loss (not used for current work) */
// // [[Rcpp::export]]
// NumericVector tanh_loss(NumericVector r, double gamma) {
//   int n = r.size();
//   NumericVector val(n);
//   
//   for (int i = 0; i < r.size(); i++) {
//     val[i] = gamma*log(cosh(r[i]/gamma));
//   }
//   
//   return val;
// }

/* First order derivative w.r.t. residual*/
// [[Rcpp::export]]
NumericVector rqHuberDerivAug(NumericVector r, NumericVector tau, double gamma) {
  int n = r.size();
  NumericVector result(n);
  
  for(int i = 0; i < n; i++){
    if(std::abs(r[i]) <= gamma){
      result[i] = ( r[i]/gamma + (2*tau[i]-1))/2;
    } else{
      result[i] = ((r[i] > 0 ? 1: -1) + (2*tau[i]-1))/2;
      //result = (sign(r[i]) + (2*tau[i]-1))/2;
    }
  }
 
  return result;
}

/* Negative gradient of huberized quantile loss (w.r.t. beta) */
// [[Rcpp::export]]
NumericVector negGradientAug(NumericVector r, arma::vec weights, NumericVector tau, double gamma, arma::sp_mat x, int ntau) {
  int n = r.size();
  //int p = x.ncol();
  int p = x.n_cols;
  NumericVector grad(p);
  //NumericVector deriv_b(n); 
  NumericVector deriv = rqHuberDerivAug(r, tau, gamma);
  arma::vec deriv_arma(deriv.begin(), n, false); // Convert deriv to Armadillo vector
  
  for (int j = 0; j < p; j++) {
    //arma::vec col_j = weights_arma % x.col(j) % deriv_arma;
    //std::copy(col_j.begin(), col_j.end(), deriv_b.begin()); // Copy the result to the NumericVector
    grad[j] = accu(weights % x.col(j) % deriv_arma);	
  }
  
  return grad/n*ntau;
  
}

// /* l2norm */
// // [[Rcpp::export]]
// double l2_norm(NumericVector x) {
//   double result = 0.0;
//   int n = x.size();
//   
//   for (int i = 0; i < n; i++) {
//     result += x[i] * x[i];
//   }
//   
//   return sqrt(result);
// }

/* l2norm weighted */
// [[Rcpp::export]]
double weightedNorm(Rcpp::NumericVector x, Rcpp::NumericVector normweights) {
  double result = 0.0;
  int n = x.size();
  
  for (int i = 0; i < n; i++) {
    result += normweights[i] * x[i] * x[i];
  }
  
  return sqrt(result);
}


/* coordinate descent for solving beta */
// [[Rcpp::export]]
List solvebetaCpp(arma::sp_mat x, arma::vec y, int n, NumericVector tau, double gamma, arma::vec weights, 
                  NumericVector groupIndex, double lambdaj, NumericVector wlambda, NumericVector wtau, 
                  NumericVector eigenval, NumericVector betaini, int maxIter, double epsilon, int ntau){
  
  //int nAug = n*ntau;
  int nbeta = betaini.size();
  NumericVector int0 = betaini[Rcpp::Range(0,(ntau-1))];
  NumericVector beta0 = betaini[Rcpp::Range(ntau,(nbeta-1))];
  NumericVector beta1 = beta0;
  NumericVector int0Aug = rep_each(int0, n);
  NumericVector unique_groups = unique(groupIndex);
  NumericVector unique_groups1 = stlSort(unique_groups);
  //Rcout << "unique_groups1 " << unique_groups1 << std::endl;
  int nGroup = groupIndex.length()/ntau;
  
  
  arma::vec int0Aug_arma = as<arma::vec>(wrap(int0Aug));
  arma::vec beta0_arma = as<arma::vec>(wrap(beta0));
  arma::vec r0 = y-int0Aug_arma-x*beta0_arma;
  Rcpp::NumericVector r01  = as<NumericVector>(wrap(r0));
  
  double delta = 2; 
  double iter = 0;
  NumericVector beta0_all(nbeta);
  NumericVector beta1_all(nbeta);
  NumericVector u(ntau);
  
  while (delta>epsilon && iter<maxIter){
    iter++;
    arma::mat Xint = arma::kron(arma::eye<arma::mat>(ntau, ntau), ones(n));
    NumericVector int1 = int0 + negGradientAug(r01, weights, tau, gamma, sp_mat(Xint), ntau)*gamma;
    r01 = r01-rep_each(int1, n)+rep_each(int0, n);
    
    for (int i = 0; i < nGroup; i++) {
      int k = unique_groups1[i];
      //NumericVector ind = find_indices(groupIndex, k);
      IntegerVector ind = seq(i * ntau + 1, (i+1)*ntau);
      arma::sp_mat subX = x.cols(i*ntau, (i+1)*ntau-1);
      u = negGradientAug(r01, weights, tau, gamma, subX, ntau);
      NumericVector subBeta = beta0[Rcpp::Range(i * ntau, (i+1)*ntau-1)];
      //Rcout << "sumU " << sum(u) << std::endl;
      NumericVector temp1 = (u+subBeta*eigenval[k-1]);
      //double temp2 = l2_norm(temp1);
      double temp2 = weightedNorm(temp1, wtau);
      NumericVector beta1_k(ntau);
      if(temp2>lambdaj*wlambda[k-1]){
        beta1_k = temp1*(1-lambdaj*wlambda[k-1]/temp2)/eigenval[k-1];
      }
      beta1[ind-1] = beta1_k;
      NumericVector betadiff = beta1_k-subBeta;
      
      arma::vec betadiff_arma = as<arma::vec>(wrap(betadiff));
      arma::vec r0_arma = as<arma::vec>(wrap(r01));
      //arma::vec r0 = r0_arma-subX*betadiff_arma;
      //Rcpp::NumericVector r01  = as<NumericVector>(wrap(r0));
      r0 = r0_arma-subX*betadiff_arma;
      r01  = as<NumericVector>(wrap(r0));
      //Rcout << "sumR0 " << sum(r01) << std::endl;
      
    }
    
    beta0_all[Rcpp::Range(0,(ntau-1))] = int0;
    beta0_all[Rcpp::Range(ntau,(nbeta-1))] = beta0;
    beta1_all[Rcpp::Range(0,(ntau-1))] = int1;
    beta1_all[Rcpp::Range(ntau,(nbeta-1))] = beta1;
    delta = max(abs(beta0_all-beta1_all));
    
    //Rcout << "delta " << delta << std::endl;
    
    beta0 = beta1;
    int0 = int1;
    
  }
  
  int converge;
  if(iter< maxIter){
    converge= 1;
  }else{
    converge= 0;
  }
  
  return Rcpp::List::create(Named("beta_update") = beta1_all,
                            Named("converge") = converge,
                            Named("resid") = r0,
                            Named("n_iter") = iter,
                            Named("gradient") = -u);
  
}