// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_lacunarity
NumericVector rcpp_lacunarity(const NumericMatrix& mat, const IntegerVector& r_vec, const int fun, const int ncores, const bool display_progress);
RcppExport SEXP _CGEI_rcpp_lacunarity(SEXP matSEXP, SEXP r_vecSEXP, SEXP funSEXP, SEXP ncoresSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type r_vec(r_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type fun(funSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lacunarity(mat, r_vec, fun, ncores, display_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CGEI_rcpp_lacunarity", (DL_FUNC) &_CGEI_rcpp_lacunarity, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_CGEI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}