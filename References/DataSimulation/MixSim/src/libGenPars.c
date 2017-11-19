
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#include "overlap.h"

/* WCC */
#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


/* genSigma : generates covariance matrix based on (p + 1) observations (unstable covariance matrix)
 * Parameters:
 * 		p - number of dimensions
 * 		VC - variance-covariance matrix
 */

void genSigma(int p, double **VC){

	int i,j,k,n;
	double **x, *mu;

	n = p + 1;

	MAKE_MATRIX(x, n, p);
	MAKE_VECTOR(mu, p);

	anull(mu, p);

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			/* WCC */
			#ifdef __HAVE_R_
				x[i][j] = rnorm(0.0, 1.0);
			#else
				x[i][j] = rnor(0.0, 1.0);
			#endif
			mu[j] = mu[j] + x[i][j];
		}
        }

	for (j=0;j<p;j++){
		mu[j] = mu[j] / n;
	}	
	
	Anull(VC, p, p);

	for (i=0; i<n; i++){
		for (j=0; j<p; j++){
			for (k=0; k<p; k++){
				VC[j][k] = VC[j][k] + (x[i][j] - mu[j]) * (x[i][k] - mu[k]);
			}
		}
        }

	for (j=0; j<p; j++){
		for (k=0; k<p; k++){
			VC[j][k] = VC[j][k] / (n - 1);
		}
	}

	FREE_MATRIX(x);
	FREE_VECTOR(mu);

}


/* genSigmaEcc : generates covariance matrix with prespecified eccentricity
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		emax - maximum eccentricity
 * 		S - set of variance-covariance matrices
 */

void genSigmaEcc(int p, int K, double emax, double ***S, int hom){

	int i, k;

	double dtmt, minL, maxL, e;
	double *Eig;
	double **VC, **L, **R;

	MAKE_VECTOR(Eig, p);
	MAKE_MATRIX(VC, p, p);
	MAKE_MATRIX(L, p, p);
	MAKE_MATRIX(R, p, p);

	if (hom == 0){ /* heterogeneous clusters */

		for (k=0; k<K; k++){
		
			genSigma(p, VC);
			cpy2(VC, p, p, S, k);
							
			#ifdef __HAVE_R_
				EigValDec(p, Eig, VC, &dtmt);
			#else
				cephes_symmeigens_down(p, Eig, VC, &dtmt);
			#endif

			i = vecMin(Eig, p, &minL);
			i = vecMax(Eig, p, &maxL);

			e = pow(1 - minL / maxL, 0.5);

			if (e > emax){

				Anull(L, p, p);

				for (i=0; i<p; i++){
					Eig[i] = maxL * (1 - emax * emax * (maxL - Eig[i]) / (maxL - minL));
					L[i][i] = Eig[i];
				}

				XAXt(VC, p, L, R);
				cpy2(R, p, p, S, k);

			}

		}

	} else { /* homogeneous clusters */
	
		genSigma(p, VC);
		for (k=0; k<K; k++){
			cpy2(VC, p, p, S, k);
		}
	
		#ifdef __HAVE_R_
			EigValDec(p, Eig, VC, &dtmt);
		#else
			cephes_symmeigens_down(p, Eig, VC, &dtmt);
		#endif
	
		i = vecMin(Eig, p, &minL);
		i = vecMax(Eig, p, &maxL);
		
		e = pow(1 - minL / maxL, 0.5);
		
		if (e > emax){

			Anull(L, p, p);

			for (i=0; i<p; i++){
				Eig[i] = maxL * (1 - emax * emax * (maxL - Eig[i]) / (maxL - minL));
				L[i][i] = Eig[i];
			}

			XAXt(VC, p, L, R);

			for (k=0; k<K; k++){
				cpy2(R, p, p, S, k);
			}
		
		}
	
	}


	FREE_MATRIX(VC);
	FREE_MATRIX(L);
	FREE_MATRIX(R);

}



/* genSphSigma : generates spherical covariance matrix
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		S - set of variance-covariance matrices
 */

void genSphSigma(int p, int K, double ***S, int hom){

	int i, k;
	double r;
	double **L;

	MAKE_MATRIX(L, p, p);		

	Anull(L, p, p);
	
	/* WCC */
	#ifdef __HAVE_R_
		r = runif(0.0, 1.0);
	#else
		r = runir(0.0, 1.0);
	#endif
	
	for (k=0; k<K; k++){

		/* WCC */
		#ifdef __HAVE_R_
			if (hom == 0) r = runif(0.0, 1.0);
		#else
			if (hom == 0) r = runir(0.0, 1.0);
		#endif

		for (i=0; i<p; i++){
			L[i][i] = r;
		}

		cpy2(L, p, p, S, k);
	
	}

	FREE_MATRIX(L);

}



/* genSphSigma : generates matrix of means
 * Parameters:
 * 		p - number of dimensions
 * 		K - number of components
 * 		Mu - set of mean vectors
 * 		Lbound - lower bound for the hypercube
 * 		Ubound - upper bound for the hypercube
 */

void genMu(int p, int K, double **Mu, double Lbound, double Ubound){
		
	int i, k;
	
	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
		
			/* WCC */
			#ifdef __HAVE_R_
				Mu[k][i] = runif(Lbound, Ubound);
			#else
				Mu[k][i] = runir(Lbound, Ubound);
			#endif

		}
	}

}



/* genPi : generates mixing proportions
 * Parameters:
 * 		K - number of components
 * 		PiLow - smallest possible mixing proportion
 * 		Pi - vector of mixing proportions
 */

void genPi(int K, double PiLow, double *Pi){

	int flag, k;
	double s;

	flag = 0;


	if ((PiLow >= 1) | (PiLow <= 0)){
/*		printf("Warning: PiLow is out of range... generated equal mixing proportions...\n"); */
		for (k=0; k<K; k++){
			Pi[k] = 1.0 / K;
		}
	} else {
		s = 0.0;
		for (k=0; k<K; k++){
			/* WCC */
			#ifdef __HAVE_R_
				Pi[k] = rgamma(1.0, 1.0);
			#else
				Pi[k] = rgamma(1.0);
			#endif
			s += Pi[k];
		}
		for (k=0; k<K; k++){
			Pi[k] = PiLow + Pi[k] / s * (1 - K * PiLow);
			if (Pi[k] < PiLow){
				flag = 1;
				break;
			}
		}
		if (flag == 1){
/*			printf("Warning: PiLow is too high... generated equal mixing proportions...\n"); */
			for (k=0; k<K; k++){
				Pi[k] = 1.0 / K;
			}
		}
		
		
	}

}
