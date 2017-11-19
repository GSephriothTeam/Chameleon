#include <R.h>
#include <Rinternals.h>

void runExactOverlap(int (*p1), int (*K1), double *Pi, double *Mu1, double *S1,
	double *pars, int (*lim1), double *OmegaMap1, double (*BarOmega),
	double (*MaxOmega), double (*EigOmega), int *rcMax);

void runOmegaClust(double (*Omega1), int (*method1), int (*p1), int (*K1), 
	double (*PiLow1), double (*Lbound1), double (*Ubound1), double (*emax1),
	double *pars, int (*lim1), int (*resN1), int (*sph1), int (*hom1),
	double *Pi, double *Mu1, double *S1, double *OmegaMap1,
	double (*BarOmega),
	double (*MaxOmega), double (*EigOmega), int *rcMax, int (*fail));

void runOmegaBarOmegaMax(int (*p1), int (*K1), double (*PiLow1),
	double (*Lbound1),
	double (*Ubound1), double (*emax1), double *pars, int (*lim1),
	int (*resN1),
	int (*sph1), double *Pi, double *Mu1, double *S1, double *OmegaMap1,
	double (*BarOmega), double (*MaxOmega),	int *rcMax, int (*fail));

void runAdjRand(int (*n), int (*K1), int (*K2), int *id1, int *id2,
	double (*Rand),
	double (*aRand), double (*F));

void runProAgree(int (*n), int (*K1), int (*K2), int *id1, int *id2,
	double (*maxPro));

void runVarInf(int (*n), int (*K1), int (*K2), int *id1, int *id2,
	double (*VI));

void runPerms(int (*n1), int (*permN1), int *perms);

