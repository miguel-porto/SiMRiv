#include "SiMRiv.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[]  = {
  {"_simulate_individuals", (DL_FUNC) &_simulate_individuals, 6},
  {"stepRasterAccumulator", (DL_FUNC) &stepRasterAccumulator, 3},
  {NULL, NULL, 0}
};

void R_init_SiMRiv(DllInfo *info)
{

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}
