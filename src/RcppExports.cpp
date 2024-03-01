// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// neg_gradient_aug
NumericVector neg_gradient_aug(NumericVector r, arma::vec weights, NumericVector tau, double gamma, arma::sp_mat x, int ntau);
RcppExport SEXP _rqPen_neg_gradient_aug(SEXP rSEXP, SEXP weightsSEXP, SEXP tauSEXP, SEXP gammaSEXP, SEXP xSEXP, SEXP ntauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ntau(ntauSEXP);
    rcpp_result_gen = Rcpp::wrap(neg_gradient_aug(r, weights, tau, gamma, x, ntau));
    return rcpp_result_gen;
END_RCPP
}
// weighted_norm
double weighted_norm(Rcpp::NumericVector x, Rcpp::NumericVector normweights);
RcppExport SEXP _rqPen_weighted_norm(SEXP xSEXP, SEXP normweightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type normweights(normweightsSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_norm(x, normweights));
    return rcpp_result_gen;
END_RCPP
}
// solvebetaCpp
List solvebetaCpp(arma::sp_mat x, arma::vec y, int n, NumericVector tau, double gamma, arma::vec weights, NumericVector groupIndex, double lambdaj, NumericVector wlambda, NumericVector wtau, NumericVector eigenval, NumericVector betaini, int maxIter, double epsilon, int ntau);
RcppExport SEXP _rqPen_solvebetaCpp(SEXP xSEXP, SEXP ySEXP, SEXP nSEXP, SEXP tauSEXP, SEXP gammaSEXP, SEXP weightsSEXP, SEXP groupIndexSEXP, SEXP lambdajSEXP, SEXP wlambdaSEXP, SEXP wtauSEXP, SEXP eigenvalSEXP, SEXP betainiSEXP, SEXP maxIterSEXP, SEXP epsilonSEXP, SEXP ntauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groupIndex(groupIndexSEXP);
    Rcpp::traits::input_parameter< double >::type lambdaj(lambdajSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wlambda(wlambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wtau(wtauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eigenval(eigenvalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betaini(betainiSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type ntau(ntauSEXP);
    rcpp_result_gen = Rcpp::wrap(solvebetaCpp(x, y, n, tau, gamma, weights, groupIndex, lambdaj, wlambda, wtau, eigenval, betaini, maxIter, epsilon, ntau));
    return rcpp_result_gen;
END_RCPP
}