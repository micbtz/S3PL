// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// IRTlikRcppA
arma::mat IRTlikRcppA(arma::vec par, arma::mat data, bool pen, arma::vec nodes, arma::vec weights, double lambda);
RcppExport SEXP _S3PL_IRTlikRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP penSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< bool >::type pen(penSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(IRTlikRcppA(par, data, pen, nodes, weights, lambda));
    return rcpp_result_gen;
END_RCPP
}
// gradIRTlikRcppA
arma::mat gradIRTlikRcppA(arma::vec par, arma::mat data, bool pen, arma::vec nodes, arma::vec weights, double lambda);
RcppExport SEXP _S3PL_gradIRTlikRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP penSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< bool >::type pen(penSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gradIRTlikRcppA(par, data, pen, nodes, weights, lambda));
    return rcpp_result_gen;
END_RCPP
}
// grad_i_IRTlikRcppA
arma::mat grad_i_IRTlikRcppA(arma::vec par, arma::mat data, bool pen, arma::vec nodes, arma::vec weights, double lambda);
RcppExport SEXP _S3PL_grad_i_IRTlikRcppA(SEXP parSEXP, SEXP dataSEXP, SEXP penSEXP, SEXP nodesSEXP, SEXP weightsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< bool >::type pen(penSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_i_IRTlikRcppA(par, data, pen, nodes, weights, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_S3PL_IRTlikRcppA", (DL_FUNC) &_S3PL_IRTlikRcppA, 6},
    {"_S3PL_gradIRTlikRcppA", (DL_FUNC) &_S3PL_gradIRTlikRcppA, 6},
    {"_S3PL_grad_i_IRTlikRcppA", (DL_FUNC) &_S3PL_grad_i_IRTlikRcppA, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_S3PL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
