// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// negllc_const
double negllc_const(double N0, NumericVector t, NumericVector A);
RcppExport SEXP _genieR_negllc_const(SEXP N0SEXP, SEXP tSEXP, SEXP ASEXP) {
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
// negllc_expan
double negllc_expan(NumericVector parr, NumericVector t, NumericVector A);
RcppExport SEXP _genieR_negllc_expan(SEXP parrSEXP, SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parr(parrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negllc_expan(parr, t, A));
    return rcpp_result_gen;
END_RCPP
}
// negllc_expo
double negllc_expo(NumericVector parr, NumericVector t, NumericVector A);
RcppExport SEXP _genieR_negllc_expo(SEXP parrSEXP, SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parr(parrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negllc_expo(parr, t, A));
    return rcpp_result_gen;
END_RCPP
}
// negllc_log
double negllc_log(NumericVector parr, NumericVector t, NumericVector A);
RcppExport SEXP _genieR_negllc_log(SEXP parrSEXP, SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parr(parrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negllc_log(parr, t, A));
    return rcpp_result_gen;
END_RCPP
}
// negllc_pexpan
double negllc_pexpan(NumericVector parr, NumericVector t, NumericVector A);
RcppExport SEXP _genieR_negllc_pexpan(SEXP parrSEXP, SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parr(parrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negllc_pexpan(parr, t, A));
    return rcpp_result_gen;
END_RCPP
}
// negllc_plog
double negllc_plog(NumericVector parr, NumericVector t, NumericVector A);
RcppExport SEXP _genieR_negllc_plog(SEXP parrSEXP, SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parr(parrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negllc_plog(parr, t, A));
    return rcpp_result_gen;
END_RCPP
}
// negllc_step
double negllc_step(NumericVector parr, NumericVector t, NumericVector A);
RcppExport SEXP _genieR_negllc_step(SEXP parrSEXP, SEXP tSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parr(parrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(negllc_step(parr, t, A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_genieR_negllc_const", (DL_FUNC) &_genieR_negllc_const, 3},
    {"_genieR_negllc_expan", (DL_FUNC) &_genieR_negllc_expan, 3},
    {"_genieR_negllc_expo", (DL_FUNC) &_genieR_negllc_expo, 3},
    {"_genieR_negllc_log", (DL_FUNC) &_genieR_negllc_log, 3},
    {"_genieR_negllc_pexpan", (DL_FUNC) &_genieR_negllc_pexpan, 3},
    {"_genieR_negllc_plog", (DL_FUNC) &_genieR_negllc_plog, 3},
    {"_genieR_negllc_step", (DL_FUNC) &_genieR_negllc_step, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_genieR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
