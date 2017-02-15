#include <R_ext/Rdynload.h>
#include "RcppExports.h"

// Short command
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

// Register Call functions with R
static R_CallMethodDef callMethods[] = {
  CALLDEF(ReIns_spliceEM_splicefit_raw_Rexport, 19),
  CALLDEF(ReIns_spliceEM_shape_adj_Rexport, 19),
  CALLDEF(ReIns_spliceEM_shape_red, 22),
  CALLDEF(ReIns_stdf_cpp, 4),
  CALLDEF(ReIns_stdf2_cpp, 3),
  {NULL, NULL, 0}
};


// Register functions with R and disable symbol search
void R_init_ReIns(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}