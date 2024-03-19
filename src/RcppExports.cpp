// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// IDW_cpp
NumericVector IDW_cpp(S4& rast, const NumericVector& x, const NumericVector& sf_x, const NumericVector& sf_y, const NumericVector& sf_z, const size_t n, const double b, const double radius, const bool na_only, const int ncores, const bool display_progress);
RcppExport SEXP _CGEI_IDW_cpp(SEXP rastSEXP, SEXP xSEXP, SEXP sf_xSEXP, SEXP sf_ySEXP, SEXP sf_zSEXP, SEXP nSEXP, SEXP bSEXP, SEXP radiusSEXP, SEXP na_onlySEXP, SEXP ncoresSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4& >::type rast(rastSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sf_x(sf_xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sf_y(sf_ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sf_z(sf_zSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< const bool >::type na_only(na_onlySEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(IDW_cpp(rast, x, sf_x, sf_y, sf_z, n, b, radius, na_only, ncores, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// LoS_reference
Rcpp::IntegerVector LoS_reference(const int x0_ref, const int y0_ref, const int r, const int nc_ref);
RcppExport SEXP _CGEI_LoS_reference(SEXP x0_refSEXP, SEXP y0_refSEXP, SEXP rSEXP, SEXP nc_refSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type x0_ref(x0_refSEXP);
    Rcpp::traits::input_parameter< const int >::type y0_ref(y0_refSEXP);
    Rcpp::traits::input_parameter< const int >::type r(rSEXP);
    Rcpp::traits::input_parameter< const int >::type nc_ref(nc_refSEXP);
    rcpp_result_gen = Rcpp::wrap(LoS_reference(x0_ref, y0_ref, r, nc_ref));
    return rcpp_result_gen;
END_RCPP
}
// focal_sum
NumericMatrix focal_sum(S4& x, const NumericMatrix& x_mat, const NumericMatrix& lac, const bool na_rm, const int ncores, const bool display_progress);
RcppExport SEXP _CGEI_focal_sum(SEXP xSEXP, SEXP x_matSEXP, SEXP lacSEXP, SEXP na_rmSEXP, SEXP ncoresSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type lac(lacSEXP);
    Rcpp::traits::input_parameter< const bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(focal_sum(x, x_mat, lac, na_rm, ncores, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lacunarity
NumericVector rcpp_lacunarity(Rcpp::S4& x, const Rcpp::NumericVector& x_values, const IntegerVector& r_vec, const int fun, const int ncores, const bool display_progress);
RcppExport SEXP _CGEI_rcpp_lacunarity(SEXP xSEXP, SEXP x_valuesSEXP, SEXP r_vecSEXP, SEXP funSEXP, SEXP ncoresSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x_values(x_valuesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type r_vec(r_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type fun(funSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lacunarity(x, x_values, r_vec, fun, ncores, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// VGVI_cpp
std::vector<double> VGVI_cpp(Rcpp::S4& dsm, const Rcpp::NumericVector& dsm_values, Rcpp::S4& greenspace, const Rcpp::NumericVector& greenspace_values, const Rcpp::IntegerVector& x0, const Rcpp::IntegerVector& y0, const Rcpp::NumericVector& h0, const int radius, const int fun, const double m, const double b, const int ncores, const bool display_progress);
RcppExport SEXP _CGEI_VGVI_cpp(SEXP dsmSEXP, SEXP dsm_valuesSEXP, SEXP greenspaceSEXP, SEXP greenspace_valuesSEXP, SEXP x0SEXP, SEXP y0SEXP, SEXP h0SEXP, SEXP radiusSEXP, SEXP funSEXP, SEXP mSEXP, SEXP bSEXP, SEXP ncoresSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type dsm(dsmSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type dsm_values(dsm_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4& >::type greenspace(greenspaceSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type greenspace_values(greenspace_valuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< const int >::type fun(funSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(VGVI_cpp(dsm, dsm_values, greenspace, greenspace_values, x0, y0, h0, radius, fun, m, b, ncores, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// VVI_cpp
Rcpp::List VVI_cpp(Rcpp::S4& dsm, const Rcpp::NumericVector& dsm_values, const Rcpp::IntegerVector& x0, const Rcpp::IntegerVector& y0, const Rcpp::NumericVector& h0, const int radius, const int ncores, const bool display_progress);
RcppExport SEXP _CGEI_VVI_cpp(SEXP dsmSEXP, SEXP dsm_valuesSEXP, SEXP x0SEXP, SEXP y0SEXP, SEXP h0SEXP, SEXP radiusSEXP, SEXP ncoresSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type dsm(dsmSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type dsm_values(dsm_valuesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type y0(y0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< const int >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(VVI_cpp(dsm, dsm_values, x0, y0, h0, radius, ncores, display_progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CGEI_IDW_cpp", (DL_FUNC) &_CGEI_IDW_cpp, 11},
    {"_CGEI_LoS_reference", (DL_FUNC) &_CGEI_LoS_reference, 4},
    {"_CGEI_focal_sum", (DL_FUNC) &_CGEI_focal_sum, 6},
    {"_CGEI_rcpp_lacunarity", (DL_FUNC) &_CGEI_rcpp_lacunarity, 6},
    {"_CGEI_VGVI_cpp", (DL_FUNC) &_CGEI_VGVI_cpp, 13},
    {"_CGEI_VVI_cpp", (DL_FUNC) &_CGEI_VVI_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_CGEI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
