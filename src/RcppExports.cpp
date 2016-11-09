// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cumsum_bounded
NumericVector cumsum_bounded(NumericVector x, double low, double high);
RcppExport SEXP genieR_cumsum_bounded(SEXP xSEXP, SEXP lowSEXP, SEXP highSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type low(lowSEXP);
    Rcpp::traits::input_parameter< double >::type high(highSEXP);
    rcpp_result_gen = Rcpp::wrap(cumsum_bounded(x, low, high));
    return rcpp_result_gen;
END_RCPP
}
// negllc_const
double negllc_const(double N0, NumericVector t, NumericVector A);
RcppExport SEXP genieR_negllc_const(SEXP N0SEXP, SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type N0(N0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negllc_const(N0, t, A));
    return rcpp_result_gen;
END_RCPP
}