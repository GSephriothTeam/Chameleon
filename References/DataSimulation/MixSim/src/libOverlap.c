 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#define Inf 1e+140

#include "overlap.h"

/* WCC */
#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


/* computes parameters needed for computing overlap
 * p  - dimensionality
 * K  - number of components
 * Pi - mixing proportions
 * Mu - mean vectors
 * S  - covariance matrices
 * li, di, const1 - parameters needed for computing overlap (see theory of method)
 */

void ComputePars(int p, int K, double *Pi, double **Mu, double ***S, double ***li, double ***di, double **const1){

	int i, j, k;

	double dtmt;

	double *m1, *m2, *Eig, *detS;
	double **Ga, **L, **Ga2, **L2, **Si;
	double ***Sinv, ***Sh;


	MAKE_VECTOR(detS, K);
	MAKE_3ARRAY(Sinv, K, p, p);
	MAKE_3ARRAY(Sh, K, p, p);
	
	MAKE_VECTOR(m1, p);
	MAKE_VECTOR(m2, p);
	MAKE_VECTOR(Eig, p);
	MAKE_MATRIX(Ga, p, p);	
	MAKE_MATRIX(L, p, p);
	MAKE_MATRIX(Ga2, p, p);	
	MAKE_MATRIX(L2, p, p);
	MAKE_MATRIX(Si, p, p);


	for (k=0; k<K; k++){

		cpy1(S, k, p, p, Ga);
		#ifdef __HAVE_R_
			EigValDec(p, Eig, Ga, &dtmt);
		#else
			cephes_symmeigens_down(p, Eig, Ga, &dtmt);
		#endif
		detS[k] = dtmt;

		Anull(L, p, p);
		for (i=0; i<p; i++){
			L[i][i] = pow(Eig[i], 0.5);
		}
		XAXt2(Ga, p, L, Sh, k);

		for (i=0; i<p; i++){
			L[i][i] = 1 / Eig[i];
		}
		XAXt2(Ga, p, L, Sinv, k);

	}



	for (i=0; i<(K-1); i++){
		for (j=i+1; j<K; j++){

			cpy1(Sh, i, p, p, Ga);
			cpy1(Sinv, j, p, p, L);
			cpy1(Sinv, i, p, p, L2);
			for (k=0; k<p; k++){
				m1[k] = Mu[i][k] - Mu[j][k];
				m2[k] = -m1[k];
			}


			XAXt(Ga, p, L, Si);

			#ifdef __HAVE_R_
				EigValDec(p, Eig, Si, &dtmt);
			#else
				cephes_symmeigens_down(p, Eig, Si, &dtmt);
			#endif
			for (k=0; k<p; k++){
				li[i][j][k] = Eig[k];
			}
			

			multiply(L2, p, p, Ga, p, p, Ga2);
			matxvec(Ga2, p, p, m1, p, Eig);
			tA(Si, p, p, Ga2);
			matxvec(Ga2, p, p, Eig, p, m1);
			for (k=0; k<p; k++){
				di[i][j][k] = m1[k];
			}



			cpy1(Sh, j, p, p, Ga2);


			XAXt(Ga2, p, L2, Si);

			#ifdef __HAVE_R_
				EigValDec(p, Eig, Si, &dtmt);
			#else
				cephes_symmeigens_down(p, Eig, Si, &dtmt);
			#endif
			for (k=0; k<p; k++){
				li[j][i][k] = Eig[k];
			}

			multiply(L, p, p, Ga2, p, p, Ga);
			matxvec(Ga, p, p, m2, p, Eig);
			tA(Si, p, p, Ga);
			matxvec(Ga, p, p, Eig, p, m2);
			for (k=0; k<p; k++){
				di[j][i][k] = m2[k];
			}


			const1[i][j] = log((Pi[j]*Pi[j]) / (Pi[i]*Pi[i]) * detS[i]/detS[j]);
			const1[j][i] = -const1[i][j];

		}

	}

	


	FREE_VECTOR(detS);
	FREE_3ARRAY(Sinv);
	FREE_3ARRAY(Sh);
	
	FREE_VECTOR(m1);
	FREE_VECTOR(m2);
	FREE_VECTOR(Eig);
	FREE_MATRIX(Ga);
	FREE_MATRIX(L);
	FREE_MATRIX(Ga2);
	FREE_MATRIX(L2);
	FREE_MATRIX(Si);

}


/* calculates the map of misclassificatons
 * c  - inflation parameter
 * p  - dimensionality
 * K  - number of components
 * li, di, const1 - parameters needed for computing overlap (see theory of method)
 * fix - fixed clusters that do not participate in inflation/deflation
 * pars, lim - parameters for qfc function
 * asympt - flag for regular or asymptotic overlap
 * OmegaMap - map of misclassification probabilities
 * BarOmega - average overlap
 * MaxOmega - maximum overlap
 * rcMax - contains the pair of components producing the highest overlap
 */


void GetOmegaMap(double c, int p, int K, double ***li, double ***di, double **const1, int *fix, double *pars, int lim, double asympt, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax){

	int i, j, k, hom;
	double Cnst1, t, s, TotalOmega, OmegaOverlap;
	double acc, sigma;
	
	double *Li, *Di, *ncp, *coef, *ldprod, *const2;
	double *trace;
	int ifault;
	int *df;

	MAKE_VECTOR(Li, p);
	MAKE_VECTOR(Di, p);
	MAKE_VECTOR(coef, p);
	MAKE_VECTOR(ldprod, p);
	MAKE_VECTOR(ncp, p);
	MAKE_VECTOR(const2, p);
	MAKE_VECTOR(df, p);

	MAKE_VECTOR(trace, 7);

	acc = pars[1];

	TotalOmega = 0.0;
	(*MaxOmega) = 0.0;

   	for (k=0; k<p; k++){
		df[k] = 1;
	}

	i = 0;
	j = 1;


	/* check if clusters are homogeneous */
	hom = 1;
	for (k=0; k<p; k++){
		if (li[0][1][k] != li[1][0][k]) hom = 0;
	}


	if (hom == 1){ /* homogeneous clusters */

		if (asympt == 0){
	
			while (i < (K-1)){

				Cnst1 = 0.0;
				for (k=0; k<p; k++){
					Di[k] = di[i][j][k] / pow(c, 0.5);
					Cnst1 = Cnst1 + Di[k] * Di[k];
					coef[k] = 0.0;
					ncp[k] = 0.0;
				}

				t = const1[i][j] - Cnst1;
				sigma = 2 * pow(Cnst1, 0.5);

				OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);



				Cnst1 = 0.0;
				for (k=0; k<p; k++){
					Di[k] = di[j][i][k] / pow(c, 0.5);
					Cnst1 = Cnst1 + Di[k] * Di[k];
					coef[k] = 0.0;
					ncp[k] = 0.0;					
				}

				t = const1[j][i] - Cnst1;
				sigma = 2 * pow(Cnst1, 0.5);

				OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
	


				OmegaOverlap = OmegaMap[i][j] + OmegaMap[j][i];
				TotalOmega = TotalOmega + OmegaOverlap;

				if (OmegaOverlap > (*MaxOmega)){
					(*MaxOmega) = OmegaOverlap;
					rcMax[0] = i;
					rcMax[1] = j;
				}

			
				if (j < (K - 1)){
					j = j + 1;
				} else {
					i = i + 1;
					j = i + 1;
				}

			}

		}




		if (asympt == 1){
			
			while (i < (K-1)){

				if (const1[i][j] > 0){
					OmegaMap[i][j] = 1;
					OmegaMap[j][i] = 0;
				}
					
				if (const1[i][j] < 0){
					OmegaMap[i][j] = 0;
					OmegaMap[j][i] = 1;
				}
				
				if (const1[i][j] == 0){
					OmegaMap[i][j] = 0.5;
					OmegaMap[j][i] = 0.5;
				}


				OmegaOverlap = OmegaMap[i][j] + OmegaMap[j][i];
				TotalOmega = TotalOmega + OmegaOverlap;

				if (OmegaOverlap > (*MaxOmega)){
					(*MaxOmega) = OmegaOverlap;
					rcMax[0] = i;
					rcMax[1] = j;
				}

			
				if (j < (K - 1)){
					j = j + 1;
				} else {
					i = i + 1;
					j = i + 1;
				}

			}

		}
		
	}	




	if (hom == 0){ /* heterogeneous clusters */

		sigma = 0.0;

		if (asympt == 0){
	
			while (i < (K-1)){

				if (fix[i] == 1){
				
					for (k=0; k<p; k++){
						Di[k] = di[i][j][k];
					}

					if (fix[j] == 1){
						for (k=0; k<p; k++){
							Li[k] = li[i][j][k];
						}
						Cnst1 = const1[i][j];
					} else {
						for (k=0; k<p; k++){
							Li[k] = li[i][j][k] / c;
						}	
						Cnst1 = const1[i][j] - p * log(c);
					}

				} else {

					for (k=0; k<p; k++){
						Di[k] = di[i][j][k] / pow(c, 0.5);
					}

					if (fix[j] == 1){
						for (k=0; k<p; k++){
							Li[k] = c * li[i][j][k];
						}
						Cnst1 = const1[i][j] + p * log(c);
					} else {
						for (k=0; k<p; k++){
							Li[k] = li[i][j][k];
						}	
						Cnst1 = const1[i][j];
					}

				}


				s = 0;
				for (k=0; k<p; k++){
					coef[k] = Li[k] - 1.0;
					ldprod[k] = Li[k] * Di[k];
					const2[k] = ldprod[k] * Di[k] / coef[k];
					s = s + const2[k];
					ncp[k] = pow(ldprod[k] / coef[k], 2);
				}
				t = s + Cnst1;

				OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);



				if (fix[j] == 1){
				
					for (k=0; k<p; k++){
						Di[k] = di[j][i][k];
					}

					if (fix[i] == 1){
						for (k=0; k<p; k++){
							Li[k] = li[j][i][k];
						}
						Cnst1 = const1[j][i];
					} else {
			      		for (k=0; k<p; k++){
							Li[k] = li[j][i][k] / c;
						}	
						Cnst1 = const1[j][i] - p * log(c);
					}

				} else {

					for (k=0; k<p; k++){
						Di[k] = di[j][i][k] / pow(c, 0.5);
					}

					if (fix[i] == 1){
						for (k=0; k<p; k++){
							Li[k] = c * li[j][i][k];
						}
						Cnst1 = const1[j][i] + p * log(c);
					} else {
			      		for (k=0; k<p; k++){
							Li[k] = li[j][i][k];
						}	
						Cnst1 = const1[j][i];
					}

				}

				s = 0;
				for (k=0; k<p; k++){
					coef[k] = Li[k] - 1.0;
					ldprod[k] = Li[k] * Di[k];
					const2[k] = ldprod[k] * Di[k] / coef[k];
					s = s + const2[k];
					ncp[k] = pow(ldprod[k] / coef[k], 2);
				}
				t = s + Cnst1;

				OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
	


				OmegaOverlap = OmegaMap[i][j] + OmegaMap[j][i];
				TotalOmega = TotalOmega + OmegaOverlap;

				if (OmegaOverlap > (*MaxOmega)){
					(*MaxOmega) = OmegaOverlap;
					rcMax[0] = i;
					rcMax[1] = j;
				}

			
				if (j < (K - 1)){
					j = j + 1;
				} else {
					i = i + 1;
					j = i + 1;
				}

			}

		}




		if (asympt == 1){

			while (i < (K-1)){

				if (fix[i] == 1){
				
					if (fix[j] == 1){
					
						s = 0;
						for (k=0; k<p; k++){
							Di[k] = di[i][j][k];
							Li[k] = li[i][j][k];
							coef[k] = Li[k] - 1.0;
							ldprod[k] = Li[k] * Di[k];
							const2[k] = ldprod[k] * Di[k] / coef[k];
							s = s + const2[k];
							ncp[k] = pow(ldprod[k] / coef[k], 2);
						}
						t = s + const1[i][j];

						OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);

					} else {

						OmegaMap[i][j] = 0.0;

					}

				} else {

					if (fix[j] == 1){

						OmegaMap[i][j] = 0.0;

					} else {
			      		for (k=0; k<p; k++){
							coef[k] = li[i][j][k] - 1.0;
							ncp[k] = 0.0;
						}
						t = const1[i][j];

						OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);					

					}
				}




				if (fix[j] == 1){
				
					if (fix[i] == 1){
					
						s = 0;
						for (k=0; k<p; k++){
							Di[k] = di[j][i][k];
							Li[k] = li[j][i][k];
							coef[k] = Li[k] - 1.0;
							ldprod[k] = Li[k] * Di[k];
							const2[k] = ldprod[k] * Di[k] / coef[k];
							s = s + const2[k];
							ncp[k] = pow(ldprod[k] / coef[k], 2);
						}
						t = s + const1[j][i];
				
						OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);

					} else {

						OmegaMap[j][i] = 0.0;

					}

				} else {

					if (fix[i] == 1){

						OmegaMap[j][i] = 0.0;

					} else {
						for (k=0; k<p; k++){
							coef[k] = li[j][i][k] - 1.0;
							ncp[k] = 0.0;
						}
						t = const1[j][i];
					
						OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);

					}
				}


				OmegaOverlap = OmegaMap[i][j] + OmegaMap[j][i];
				TotalOmega = TotalOmega + OmegaOverlap;

				if (OmegaOverlap > (*MaxOmega)){
					(*MaxOmega) = OmegaOverlap;
					rcMax[0] = i;
					rcMax[1] = j;
				}

			
				if (j < (K - 1)){
					j = j + 1;
				} else {
					i = i + 1;
					j = i + 1;
				}

			}

		}
		
	}	

	(*BarOmega) = TotalOmega / (K * (K - 1) / 2.0);

	for (k=0; k<K; k++){
		OmegaMap[k][k] = 1.0;
	}


	FREE_VECTOR(Li);
	FREE_VECTOR(Di);
	FREE_VECTOR(coef);
	FREE_VECTOR(ldprod);
	FREE_VECTOR(ncp);
	FREE_VECTOR(const2);
	FREE_VECTOR(df);

	FREE_VECTOR(trace);
	
}



/* computes the measure of overlap based on the largest eigenvalue
 * K  - number of components
 * OmegaMap - map of misclassification probabilities 
 */

double GetEigOmega(int K, double **OmegaMap){
	
	int i, j;
	double dtmt, eigOm;
	double *Eig;
	double **W2;

	MAKE_VECTOR(Eig, K);
	MAKE_MATRIX(W2, K, K);
	
	for (i=1; i<K; i++){
		for (j=0; j<i; j++){
			W2[i][j] = (OmegaMap[i][j] + OmegaMap[j][i]);
			W2[j][i] = W2[i][j];
		}
	}

	for (i=0; i<K; i++){
		W2[i][i] = 1.0;
	}

	#ifdef __HAVE_R_
		EigValDec(K, Eig, W2, &dtmt);
	#else
		cephes_symmeigens_down(K, Eig, W2, &dtmt);
	#endif


	eigOm = (Eig[K-1] - 1) / (K - 1);
	
	FREE_MATRIX(W2);
	FREE_VECTOR(Eig);

	return eigOm;
	
}


/* computes the exact overlap
 * p  - dimensionality
 * K  - number of components
 * Pi - mixing proportions
 * Mu - mean vectors
 * S  - covariance matrices
 * pars, lim - parameters for qfc function
 * OmegaMap - map of misclassification probabilities
 * BarOmega - average overlap
 * MaxOmega - maximum overlap
 * rcMax - contains the pair of components producing the highest overlap
 */

void ExactOverlap(int p, int K, double *Pi, double **Mu, double ***S, double *pars, int lim,
	double **OmegaMap, double (*BarOmega), double (*MaxOmega), double (*EigOmega), int *rcMax){

	double c, Balpha, Malpha, Ealpha;
	int asympt;

	int  *fix;
	double **const1;
	double ***li, ***di;

	MAKE_VECTOR(fix, K);

	MAKE_3ARRAY(li, K, K, p);
	MAKE_3ARRAY(di, K, K, p);
	MAKE_MATRIX(const1, K, K);

	Balpha = (*BarOmega);
	Malpha = (*MaxOmega);
	Ealpha = (*EigOmega);

	ComputePars(p, K, Pi, Mu, S, li, di, const1);

	anulli(fix, K);

	c = 1.0;
	asympt = 0;
	GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);
	Ealpha = GetEigOmega(K, OmegaMap);
	
	(*BarOmega) = Balpha;
	(*MaxOmega) = Malpha;
	(*EigOmega) = Ealpha;

	FREE_VECTOR(fix);

	FREE_3ARRAY(li);
	FREE_3ARRAY(di);
	FREE_MATRIX(const1);

}



/* FIND MULTIPLIER C ON THE INTERVAL (lower, upper)
 * lower - lower bound of the interval
 * upper - upper bound of the interval
 * Omega - overlap value
 * method - average or maximum overlap
 * p  - dimensionality
 * K  - number of components
 * li, di, const1 - parameters needed for computing overlap (see theory of method)
 * fix - fixed clusters that do not participate in inflation/deflation
 * pars, lim - parameters for qfc function
 * c  - inflation parameter
 * OmegaMap - map of misclassification probabilities
 * BarOmega - average overlap
 * MaxOmega - maximum overlap
 * EigOmega - eigenvalue overlap
 * rcMax - contains the pair of components producing the highest overlap
 */


void FindC(double lower, double upper, double Omega, int method, int p, int K, double ***li, double ***di, double **const1, int *fix, double *pars, int lim, double (*c), double **OmegaMap, double (*BarOmega), double (*MaxOmega), double (*EigOmega), int *rcMax){

	double diff, eps;
	int sch, asympt, stopIter;

	eps = pars[0];

	diff = Inf;
	stopIter = 1000;

	sch = 0;

	while (fabs(diff) > eps){

		(*c) = (lower + upper) / 2.0;

		asympt = 0;
		GetOmegaMap((*c), p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &(*BarOmega), &(*MaxOmega), rcMax);
	
		if (method == 0){

			if ((*BarOmega) < Omega){ /* clusters are too far */
				lower = (*c);
			} else {
				upper = (*c);
			}
			
			diff = (*BarOmega) - Omega;

		} else {
				
			if (method == 1){

				if ((*MaxOmega) < Omega){ /* clusters are too far */
					lower = (*c);
				} else {
					upper = (*c);
				}
			
				diff = (*MaxOmega) - Omega;			

			} else {

				(*EigOmega) = GetEigOmega(K, OmegaMap);

				if ((*EigOmega) < Omega){ /* clusters are too far */
					lower = (*c);
				} else {
					upper = (*c);
				}
			
				diff = (*EigOmega) - Omega;				
			
			}
			
		}

		sch = sch + 1;

		if (sch == stopIter){
			(*c) = 0.0;
/*			printf("Error: required overlap was not reached in %i iterations...\n", stopIter); */
			break;
		}

	}
	
	if ((method == 0) | (method == 1)) (*EigOmega) = GetEigOmega(K, OmegaMap);

}


/* FIND MULTIPLIER C ON THE INTERVAL (lower, upper)
 * lower - lower bound of the interval
 * upper - upper bound of the interval
 * Omega - overlap value
 * method - average or maximum overlap
 * p  - dimensionality
 * K  - number of components
 * li, di, const1 - parameters needed for computing overlap (see theory of method)
 * fix - fixed clusters that do not participate in inflation/deflation
 * pars, lim - parameters for qfc function
 * c  - inflation parameter
 * OmegaMap - map of misclassification probabilities
 * EigOmega - average overlap
 * rcMax - contains the pair of components producing the highest overlap
 */

/*
void FindCeig(double lower, double upper, double Omega, int method, int p, int K, double ***li, double ***di, double **const1, int *fix, double *pars, int lim, double (*c), double **OmegaMap, double (*BarOmega), double (*MaxOmega), double (*EigOmega), int *rcMax){

	double diff, eps;
	int sch, asympt, stopIter;

	eps = pars[0];

	diff = Inf;
	stopIter = 1000;

	sch = 0;

	while (fabs(diff) > eps){

		(*c) = (lower + upper) / 2.0;

		asympt = 0;
		GetOmegaMap((*c), p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &(*BarOmega), &(*MaxOmega), rcMax);
	
		(*EigOmega) = GetEigOmega(K, OmegaMap);

		if ((*EigOmega) < Omega){ // clusters are too far
			lower = (*c);
		} else {
			upper = (*c);
		}
			
		diff = (*EigOmega) - Omega;				
			
		sch = sch + 1;

		if (sch == stopIter){
			(*c) = 0.0;
			break;
		}

	}
	
}
*/

/* run the procedure when average or maximum overlap is specified
 * 
 * Omega - overlap value
 * method - average or maximum overlap
 * p  - dimensionality
 * K  - number of components
 * PiLow - smallest mixing proportion allowed
 * Lbound - lower bound for uniform hypercube at which mean vectors at simulated
 * Ubound - upper bound for uniform hypercube at which mean vectors at simulated
 * emax - maximum eccentricity
 * pars, lim - parameters for qfc function
 * resN - number of resamplings allowed
 * sph - sperical covariance matrices
 * hom - homogeneous covariance matrices
 * Pi - mixing proportions
 * Mu - mean vectors
 * S  - covariance matrices
 * OmegaMap - map of misclassification probabilities
 * BarOmega - average overlap
 * MaxOmega - maximum overlap
 * rcMax - contains the pair of components producing the highest overlap
 * fail - flag indicating if the process failed
 */
void OmegaClust(double Omega, int method, int p, int K, double PiLow, double Lbound,
	double Ubound,	double emax, double *pars, int lim, int resN, int sph, int hom,
	double *Pi, double **Mu, double ***S, double **OmegaMap, double (*BarOmega),
	double (*MaxOmega), double (*EigOmega), int *rcMax, int (*fail)){

	int asympt, sch;
	double c, diff, lower, upper, eps, Balpha, Malpha, Ealpha;

	int *fix;
	double **const1;
	double ***li, ***di;


	MAKE_3ARRAY(li, K, K, p);
	MAKE_3ARRAY(di, K, K, p);
	MAKE_MATRIX(const1, K, K);

	MAKE_VECTOR(fix, K);

	anulli(fix, K);

	eps = pars[0];

	Balpha = (*BarOmega);
	Malpha = (*MaxOmega);
	Ealpha = (*EigOmega);

	sch = 0;

	do{

		(*fail) = 0;

		/* generate parameters */
/*		printf("Simulating dataset...\n"); */

		genPi(K, PiLow, Pi);

		genMu(p, K, Mu, Lbound, Ubound);
		if (sph == 0){
			genSigmaEcc(p, K, emax, S, hom);
		} else {
			genSphSigma(p, K, S, hom);
		}


		/* prepare parameters */
		
		ComputePars(p, K, Pi, Mu, S, li, di, const1);

		/* check if desired overlap is reachable */

		asympt = 1;
		c = 0.0;
		GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);
		Ealpha = GetEigOmega(K, OmegaMap);

		if (method == 0){
			diff = Balpha - Omega;
		} else {
			if (method == 1){
				diff = Malpha - Omega;
			} else {
				diff = Ealpha - Omega;
			}
		}

		if (diff < -eps){ /* not reachable */
/*			printf("Warning: the desired overlap cannot be reached...\n"); */
			(*fail) = 1;
		} else {

			lower = 0.0;
			upper = pow(2.0, 2);

			do{
				FindC(lower, upper, Omega, method, p, K, li, di, const1, fix, pars, lim, &c, OmegaMap, &Balpha, &Malpha, &Ealpha, rcMax);

				lower = upper;
				upper = upper * upper;
				if (upper > 1000000){
/*					printf("Warning: the desired overlap cannot be reached...\n");
					printf("Simulating another dataset...\n"); */
					(*fail) = 1;			
					break;
				} /* (upper > 1000000) prevents nonstopping loops */

			} while(c == 0);

		}


		if ((*fail) == 0){

			/* correct covariances by multiplier C */

			cxS(p, K, S, c);

			break;

		}

		sch = sch + 1;

		if (sch == resN){
			// WCC: R doesn't like printf(...).
#ifndef __HAVE_R_
			printf("Error: the desired overlap has not been reached in %i simulations...\n", resN);
			printf("Increase the number of simulations allowed (option resN) or change the value of overlap...\n");
#else
			Rprintf("Error: the desired overlap has not been reached in %i simulations...\n", resN);
			Rprintf("Increase the number of simulations allowed (option resN) or change the value of overlap...\n");
#endif
			(*fail) = 1;
			break;
		}

	} while ((*fail) == 1);

	(*BarOmega) = Balpha;
	(*MaxOmega) = Malpha;
	(*EigOmega) = Ealpha;

	FREE_3ARRAY(li);
	FREE_3ARRAY(di);
	FREE_MATRIX(const1);

	FREE_VECTOR(fix);


}




/* run the procedure when average and maximum overlaps are both specified
 * 
 * p  - dimensionality
 * K  - number of components
 * PiLow - smallest mixing proportion allowed
 * Lbound - lower bound for uniform hypercube at which mean vectors at simulated
 * Ubound - upper bound for uniform hypercube at which mean vectors at simulated
 * emax - maximum eccentricity
 * pars, lim - parameters for qfc function
 * resN - number of resamplings allowed
 * sph - sperical covariance matrices
 * Pi - mixing proportions
 * Mu - mean vectors
 * S  - covariance matrices
 * OmegaMap - map of misclassification probabilities
 * BarOmega - average overlap
 * MaxOmega - maximum overlap
 * rcMax - contains the pair of components producing the highest overlap
 * fail - flag indicating if the process failed
 */
void OmegaBarOmegaMax(int p, int K, double PiLow, double Lbound, double Ubound, double emax,
	double *pars, int lim, int resN, int sph, double *Pi, double **Mu, double ***S,
	double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax, int (*fail)){


	int i, j, k, asympt, sch, method;
	double c, diff, lower = 0, upper = 0, eps, Balpha, Malpha, Ealpha;

	int *fix, *fix2;
	double **const1, **const12, **OmegaMap2;
	double ***li, ***di, ***li2, ***di2;

	MAKE_3ARRAY(li, K, K, p);
	MAKE_3ARRAY(di, K, K, p);
	MAKE_MATRIX(const1, K, K);
	MAKE_3ARRAY(li2, 2, 2, p);
	MAKE_3ARRAY(di2, 2, 2, p);
	MAKE_MATRIX(const12, 2, 2);
	MAKE_MATRIX(OmegaMap2, 2, 2);

	MAKE_VECTOR(fix, K);
	MAKE_VECTOR(fix2, 2);

	anulli(fix, K);
	anulli(fix2, 2);

	eps = pars[0];

	Balpha = (*BarOmega);
	Malpha = (*MaxOmega);
	
	Ealpha = GetEigOmega(K, OmegaMap);

	(*fail) = 1;

	if ((Malpha < Balpha) | (Malpha > Balpha * K * (K - 1) / 2.0)){ /* wrong parameters*/

#ifndef __HAVE_R_
		printf("Error: incorrect values of average and maximum overlaps...\n");
		printf("Both conditions should hold:\n1. MaxOverlap > AverOverlap\n2. MaxOverlap < AverOverlap * K (K - 1) / 2\n");
#else
		Rprintf("Error: incorrect values of average and maximum overlaps...\n");
		Rprintf("Both conditions should hold:\n1. MaxOverlap > AverOverlap\n2. MaxOverlap < AverOverlap * K (K - 1) / 2\n");
#endif

	} else {


		sch = 0;
	
		do{

			/* generate parameters */
/*			printf("Simulating dataset...\n");	*/

			genPi(K, PiLow, Pi);
			genMu(p, K, Mu, Lbound, Ubound);
			if (sph == 0){
				genSigmaEcc(p, K, emax, S, 0);
			} else {
				genSphSigma(p, K, S, 0);
			}

			/* prepare parameters */
		
			ComputePars(p, K, Pi, Mu, S, li, di, const1);

			/* check if maximum overlap is reachable */

			asympt = 1;
			c = 0.0;
			GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);

			diff = Malpha - (*MaxOmega);

			if (diff >= -eps){ /* reachable */

				lower = 0.0;
				upper = pow(2.0, 10);

				do{

					/* find C for two currently largest clusters */
					
					for (i=0; i<2; i++){
						for (j=0; j<2; j++){
							for (k=0; k<p; k++){
								li2[i][j][k] = li[rcMax[i]][rcMax[j]][k];
								di2[i][j][k] = di[rcMax[i]][rcMax[j]][k];
							}
							const12[i][j] = const1[rcMax[i]][rcMax[j]];
						}
					}				

					Malpha = (*MaxOmega);
					method = 1;
					FindC(lower, upper, Malpha, method, p, 2, li2, di2, const12, fix2, pars, lim, &c, OmegaMap2, &Balpha, &Malpha, &Ealpha, rcMax);

					if (c == 0){ /* abnormal termination */
/*						printf("Warning: the desired overlap cannot be reached...\n"); */
						(*fail) = 1;
						break;
					}

					asympt = 0;
					GetOmegaMap(c, p, K, li, di, const1, fix, pars, lim, asympt, OmegaMap, &Balpha, &Malpha, rcMax);
					upper = c;

			       		diff = Balpha - (*BarOmega);
					if (diff < -eps){ /* BarOmega is not reachable */
/*						printf("Warning: the desired overlap cannot be reached...\n"); */
						(*fail) = 1;
						break;
					}

					diff = Malpha - (*MaxOmega);
					if (diff < eps){ /* MaxOmega has been reached */
						(*fail) = 0;
						break;
					}

				} while ((*fail) != 0);

			}	



			if ((*fail) == 0){ /* OmegaMax is reached and OmegaBar is reachable */

				/* correct covariances by multiplier C */

				cxS(p, K, S, c);

				ComputePars(p, K, Pi, Mu, S, li, di, const1);

				fix[rcMax[0]] = 1;
				fix[rcMax[1]] = 1;
				upper = 1;

				Balpha = (*BarOmega);
				method = 0;
				FindC(lower, upper, Balpha, method, p, K, li, di, const1, fix, pars, lim, &c, OmegaMap, &Balpha, &Malpha, &Ealpha, rcMax);

				/* correct covariances by multiplier C */

				for (k=0; k<K; k++){
					for (i=0; i<p; i++){
						for (j=0; j<p; j++){
							if (fix[k] == 0){
								S[k][i][j] = c * S[k][i][j];
							}
						}
					}
				}

				break;

			}

			
			sch = sch + 1;

			if (sch == resN){
#ifndef __HAVE_R_
				printf("Error: the desired overlap has not been reached in %i simulations...\n", resN);
				printf("Increase the number of simulations allowed (option resN) or change the value of overlap...\n");
#else
				Rprintf("Error: the desired overlap has not been reached in %i simulations...\n", resN);
				Rprintf("Increase the number of simulations allowed (option resN) or change the value of overlap...\n");
#endif
				(*fail) = 1;
				break;
			}

		} while ((*fail) == 1);

	}

	(*BarOmega) = Balpha;
	(*MaxOmega) = Malpha;

	FREE_3ARRAY(li);
	FREE_3ARRAY(di);
	FREE_MATRIX(const1);
	FREE_3ARRAY(li2);
	FREE_3ARRAY(di2);
	FREE_MATRIX(const12);
	FREE_MATRIX(OmegaMap2);

	FREE_VECTOR(fix);
	FREE_VECTOR(fix2);	


}

