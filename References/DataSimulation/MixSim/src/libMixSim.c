
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#include <R.h>

#include "overlap.h"


void runExactOverlap(int (*p1), int (*K1), double *Pi, double *Mu1, double *S1,
	double *pars, int (*lim1), double *OmegaMap1, double (*BarOmega), double (*MaxOmega),
	double (*EigOmega), int *rcMax){

	double **Mu, **OmegaMap;
	double ***S;

	int p, K, lim;
	double BarOmega1, MaxOmega1, EigOmega1;


	p = (*p1);
	K = (*K1);
	lim = (*lim1);


	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);
	MAKE_MATRIX(OmegaMap, K, K);

	array1to2(K, p, Mu1, Mu);
	array1to3(K, p, p, S1, S);


	ExactOverlap(p, K, Pi, Mu, S, pars, lim, OmegaMap, &BarOmega1, &MaxOmega1, &EigOmega1, rcMax);


	(*BarOmega) = BarOmega1;
	(*MaxOmega) = MaxOmega1;
	(*EigOmega) = EigOmega1;
	
	array2to1(K, K, OmegaMap1, OmegaMap);


	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);
	FREE_MATRIX(OmegaMap);

}



void runOmegaClust(double (*Omega1), int (*method1), int (*p1), int (*K1), 
	double (*PiLow1), double (*Lbound1), double (*Ubound1), double (*emax1),
	double *pars, int (*lim1), int (*resN1), int (*sph1), int (*hom1),
	double *Pi, double *Mu1, double *S1, double *OmegaMap1,	double (*BarOmega),
	double (*MaxOmega), double (*EigOmega), int *rcMax, int (*fail)){


	double **Mu, **OmegaMap;
	double ***S;

	int fail1, p, K, lim, method, resN, sph, hom;
	double BarOmega1, MaxOmega1, EigOmega1, Omega, PiLow, Lbound, Ubound, emax;

	GetRNGstate();

	p = (*p1);
	K = (*K1);		


	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);
	MAKE_MATRIX(OmegaMap, K, K);

	fail1 = (*fail);
	lim = (*lim1);
	method = (*method1);
	resN = (*resN1);
	sph = (*sph1);
	hom = (*hom1);
	
	Omega = (*Omega1);
	PiLow = (*PiLow1);
	Lbound = (*Lbound1);
	Ubound = (*Ubound1);
	emax = (*emax1);

	OmegaClust(Omega, method, p, K, PiLow, Lbound, Ubound, emax, pars, lim, resN,
		sph, hom, Pi, Mu, S, OmegaMap, &BarOmega1, &MaxOmega1, &EigOmega1, rcMax, &fail1);

	(*BarOmega) = BarOmega1;
	(*MaxOmega) = MaxOmega1;
	(*EigOmega) = EigOmega1;
	(*fail) = fail1;

	array2to1(K, p, Mu1, Mu);
	array3to1(K, p, p, S1, S);
	array2to1(K, K, OmegaMap1, OmegaMap);


	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);
	FREE_MATRIX(OmegaMap);

	PutRNGstate();

}




void runOmegaBarOmegaMax(int (*p1), int (*K1), double (*PiLow1), double (*Lbound1),
	double (*Ubound1), double (*emax1), double *pars, int (*lim1), int (*resN1),
	int (*sph1), double *Pi, double *Mu1, double *S1, double *OmegaMap1,
	double (*BarOmega), double (*MaxOmega),	int *rcMax, int (*fail)){


	double **Mu, **OmegaMap;
	double ***S;

	int fail1, p, K, lim, resN, sph;
	double BarOmega1, MaxOmega1, PiLow, Lbound, Ubound, emax;


	GetRNGstate();


	p = (*p1);
	K = (*K1);		


	MAKE_MATRIX(Mu, K, p);
	MAKE_3ARRAY(S, K, p, p);
	MAKE_MATRIX(OmegaMap, K, K);

	fail1 = (*fail);
	lim = (*lim1);
	resN = (*resN1);
	sph = (*sph1);
	
	PiLow = (*PiLow1);
	Lbound = (*Lbound1);	
	Ubound = (*Ubound1);
	emax = (*emax1);
	
	BarOmega1 = (*BarOmega);
	MaxOmega1 = (*MaxOmega);

	OmegaBarOmegaMax(p, K, PiLow, Lbound, Ubound, emax, pars, lim, resN, sph, Pi, Mu, S, OmegaMap, &BarOmega1, &MaxOmega1, rcMax, &fail1);

	(*BarOmega) = BarOmega1;
	(*MaxOmega) = MaxOmega1;
	(*fail) = fail1;

	array2to1(K, p, Mu1, Mu);
	array3to1(K, p, p, S1, S);
	array2to1(K, K, OmegaMap1, OmegaMap);


	FREE_MATRIX(Mu);
	FREE_3ARRAY(S);
	FREE_MATRIX(OmegaMap);

	PutRNGstate();

}

		


void runAdjRand(int (*n), int (*K1), int (*K2), int *id1, int *id2, double (*Rand), double (*aRand), double (*F)){

	double Rand1, aRand1, F1;

	int n1, K11, K21;

	n1 = (*n);
	K11 = (*K1);
	K21 = (*K2);		

	Rand1 = (*Rand);
	aRand1 = (*aRand);
	F1 = (*F);
	
	RRand(n1, K11, K21, id1, id2, &Rand1, &aRand1, &F1);

	(*Rand) = Rand1;
	(*aRand) = aRand1;
	(*F) = F1;

}




void runProAgree(int (*n), int (*K1), int (*K2), int *id1, int *id2, double (*maxPro)){

	double maxPro1;

	int n1, K11, K21;

	n1 = (*n);
	K11 = (*K1);
	K21 = (*K2);		

	maxPro1 = (*maxPro);
	
	proAgree(n1, K11, K21, id1, id2, &maxPro1);

	(*maxPro) = maxPro1;

}




void runVarInf(int (*n), int (*K1), int (*K2), int *id1, int *id2, double (*VI)){

	double VI1;

	int n1, K11, K21;

	n1 = (*n);
	K11 = (*K1);
	K21 = (*K2);		

	VI1 = (*VI);
	
	VIindex(n1, K11, K21, id1, id2, &VI1);

	(*VI) = VI1;

}




void runPerms(int (*n1), int (*permN1), int *perms){

	int permN, n;
	int **permMat;

	permN = (*permN1);
	n = (*n1);

	MAKE_MATRIX(permMat, permN, n);
	
	array1to2i(permN, n, perms, permMat);

	AllPerms(n, permMat);
	
	array2to1i(permN, n, perms, permMat);
	
	FREE_MATRIX(permMat);
	
}
