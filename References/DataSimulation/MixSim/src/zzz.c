#include <R.h>
#include <R_ext/Rdynload.h>

#include "zzz.h"

static const R_CMethodDef cMethods[] = {
	{"runExactOverlap", (DL_FUNC) &runExactOverlap, 12},
	{"runOmegaClust", (DL_FUNC) &runOmegaClust, 22},
	{"runOmegaBarOmegaMax", (DL_FUNC) &runOmegaBarOmegaMax, 18},
	{"runAdjRand", (DL_FUNC) &runAdjRand, 8},
	{"runProAgree", (DL_FUNC) &runProAgree, 6},
	{"runVarInf", (DL_FUNC) &runVarInf, 6},
	{"runPerms", (DL_FUNC) &runPerms, 3},

	/* Finish R_CMethodDef. */
	{NULL, NULL, 0}
};
/* End of the callMethods[]. */


void R_init_MixSim(DllInfo *info){
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
} /* End of R_init_MixSim(). */
