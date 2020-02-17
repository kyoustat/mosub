// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fast_loss_prj
arma::vec fast_loss_prj(int nS, int dS, int mS, arma::mat PS, arma::mat xS, arma::vec muS);
RcppExport SEXP _mosub_fast_loss_prj(SEXP nSSEXP, SEXP dSSEXP, SEXP mSSEXP, SEXP PSSEXP, SEXP xSSEXP, SEXP muSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nS(nSSEXP);
    Rcpp::traits::input_parameter< int >::type dS(dSSEXP);
    Rcpp::traits::input_parameter< int >::type mS(mSSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type PS(PSSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xS(xSSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type muS(muSSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_loss_prj(nS, dS, mS, PS, xS, muS));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _mosub_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _mosub_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _mosub_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _mosub_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mosub_fast_loss_prj", (DL_FUNC) &_mosub_fast_loss_prj, 6},
    {"_mosub_rcpparma_hello_world", (DL_FUNC) &_mosub_rcpparma_hello_world, 0},
    {"_mosub_rcpparma_outerproduct", (DL_FUNC) &_mosub_rcpparma_outerproduct, 1},
    {"_mosub_rcpparma_innerproduct", (DL_FUNC) &_mosub_rcpparma_innerproduct, 1},
    {"_mosub_rcpparma_bothproducts", (DL_FUNC) &_mosub_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mosub(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
