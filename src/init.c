#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ReIns_spliceEM_shape_adj_Rexport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ReIns_spliceEM_shape_red(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ReIns_spliceEM_splicefit_raw_Rexport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ReIns_stdf_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ReIns_stdf2_cpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_ReIns_spliceEM_shape_adj_Rexport",     (DL_FUNC) &_ReIns_spliceEM_shape_adj_Rexport,     19},
  {"_ReIns_spliceEM_shape_red",             (DL_FUNC) &_ReIns_spliceEM_shape_red,             22},
  {"_ReIns_spliceEM_splicefit_raw_Rexport", (DL_FUNC) &_ReIns_spliceEM_splicefit_raw_Rexport, 19},
  {"_ReIns_stdf_cpp",                       (DL_FUNC) &_ReIns_stdf_cpp,                        4},
  {"_ReIns_stdf2_cpp",                      (DL_FUNC) &_ReIns_stdf2_cpp,                       3},
  {NULL, NULL, 0}
};

void R_init_ReIns(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}