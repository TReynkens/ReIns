#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ReIns_spliceEM_shape_adj_Rexport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ReIns_spliceEM_shape_red(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ReIns_spliceEM_splicefit_raw_Rexport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ReIns_stdf_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP ReIns_stdf2_cpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"ReIns_spliceEM_shape_adj_Rexport",     (DL_FUNC) &ReIns_spliceEM_shape_adj_Rexport,     19},
  {"ReIns_spliceEM_shape_red",             (DL_FUNC) &ReIns_spliceEM_shape_red,             22},
  {"ReIns_spliceEM_splicefit_raw_Rexport", (DL_FUNC) &ReIns_spliceEM_splicefit_raw_Rexport, 19},
  {"ReIns_stdf_cpp",                       (DL_FUNC) &ReIns_stdf_cpp,                        4},
  {"ReIns_stdf2_cpp",                      (DL_FUNC) &ReIns_stdf2_cpp,                       3},
  {NULL, NULL, 0}
};

void R_init_ReIns(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}