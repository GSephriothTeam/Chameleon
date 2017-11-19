

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "array.h"
#define inf 1e+40;

#include "overlap.h"

/* This procedure computes Adjusted Rand index
 * 
 *  Parameters:
 *  N - number of observations
 * 	TRUK - true number of clusters
 * 	PREDK - estimated number of clusters
 * 	trcl - true classification vector
 * 	prcl - estimated classification vector
 * 	Rand - Rand index
 * 	adjRand - Adjusted Rand index
 * 	Eindex - E index
 */

void RRand(int N, int TRUK, int PREDK, int *trcl, int *prcl, double *Rand, double *adjRand, double *F){
  
	// int i, j, n[TRUK][PREDK];
	int i, j;
	double sumtr[TRUK], sumpr[PREDK], sumprsq, sumtrsq, sumsq, discordant, sumtrprsq;
	double term1, term2, term3, W1, W2;
	double nij2sum, nidot2sum, ndotj2sum;

	/* Fails because n is a 2D array when *TRUK and *PREDK are too large.
	for (i = 0; i < TRUK; i++){
		for (j = 0; j < PREDK; j++){
			n[i][j] = 0;
		}
	}
	*/

  /* Fixed by WCC. */
  int **n;
  n = (int **) malloc(TRUK * sizeof(int *));
  if(n == NULL) {
    error("Memory allocation fails at n!\n");
  }
  for (i=0;i<TRUK;i++) {
    n[i] = (int *) malloc(PREDK * sizeof(int));
    if(n[i] == NULL) {
      error("Memory allocation fails at n[i]!\n");
    }
    for (j=0;j<PREDK;j++) {
      n[i][j]=0;
    }
  }
  
	for (i = 0; i < N; i++){
		n[trcl[i]][prcl[i]] += 1;
	}

	sumtrsq = 0.;
	for (i = 0; i < TRUK; i++){
		sumtr[i] = 0.;
		for (j = 0; j < PREDK; j++){
			sumtr[i] += n[i][j];
		}
		sumtrsq += sumtr[i] * sumtr[i];
	}
  
	sumprsq = 0.;
	for (j = 0; j < PREDK; j++){
		sumpr[j] = 0.;
		for (i = 0; i < TRUK; i++){
			sumpr[j] += (double)n[i][j];
		}
		sumprsq += sumpr[j] * sumpr[j];
	}

	sumtrprsq = 0.;
	for (i = 0; i < TRUK; i++){
		for (j = 0; j < PREDK; j++){
			sumtrprsq += sumtr[i] * sumtr[i] * sumpr[j] * sumpr[j];
		}
	}

//	(*Eindex) = sumtrprsq / (N * ((double)N - 1) + N * (double)N / (N - 1)) - (sumprsq + sumtrsq) / (N - 1);
//	(*Eindex) *= 2.;
//	(*Eindex) /= N * ((double)N - 1);
  
	sumsq = 0.;
	for (i = 0; i < TRUK; i++){
		for (j = 0; j < PREDK; j++){
			sumsq += (double)n[i][j] * n[i][j];
		}
	}

	nij2sum = 0.;
	for (i = 0; i < TRUK; i++){
		for (j = 0; j < PREDK; j++){
			nij2sum += (double)n[i][j] * (n[i][j] - 1) / 2.0;
		}
	}

	nidot2sum = 0.;
	for (i = 0; i < TRUK; i++){
		nidot2sum += (double)sumtr[i] * (sumtr[i] - 1) / 2.0;
	}

	ndotj2sum = 0.;
	for (i = 0; i < PREDK; i++){
		ndotj2sum += (double)sumpr[i] * (sumpr[i] - 1) / 2.0;
	}

	W1 = nij2sum / nidot2sum;
	W2 = nij2sum / ndotj2sum;
	
	(*F) = pow(W1 * W2, 0.5);
	
	discordant = 0.5 * (sumtrsq + sumprsq) - sumsq;

	(*Rand) = 1.0 - discordant / ((double)N * ((double)N - 1.) / 2.);

	term3 = nidot2sum * ndotj2sum / ((double)N * ((double)N - 1.) / 2.);

	term1 = nij2sum - term3;

	term2 = (nidot2sum + ndotj2sum) / 2 - term3;

	(*adjRand) = term1 / term2;

  /* Free 2D array pointers. */
  for (i=0;i<TRUK;i++){
    free(n[i]);
  }
  free(n);

}




void proAgree(int n, int K1, int K2, int *id1, int *id2, double *maxPro){

	int sch, i, j, v, w, finish, flag, ind;
	int size, h, idsmall, trclass;
	double curr;

    int *cn;
	double **pat;

	if (K1 < K2){
		size = K1;
		idsmall = 1;
	} else {
		size = K2;
		idsmall = 2;
	}
		
	curr = 0.0;
	
	sch = 0;
	i = 0;
	j = -1;
	flag = 0;
	finish = 0;
	ind = 0;

	MAKE_MATRIX(pat, size, size);
	for (v=0; v<size; v++){
		for (w=0; w<size; w++){
			pat[v][w] = 0;
		}
	}

	MAKE_VECTOR(cn,size);
	for (v=0; v<size; v++){
		cn[v] = 0;
	}
  
	while (finish == 0){
    
		if (j != (size-1)){
			j = j+1;
		} else {
			if (flag == 1){
				j = 0;
				i = i + 1;
				flag = 0;
			}
		}
    
		if (pat[i][j] == 0){
			for (v=0; v<size; v++){
				pat[i][v] = 1;
				pat[v][j] = 1;
			}
      
			sch = sch + 1;
			cn[sch-1] = j;
			flag = 1;
		}

		if ((sch == size) & (flag == 1)){
      
//			for (v=0; v<size; v++){
//				printf(" %i", cn[v]);
//			}
//			printf(" \n");

//			####################################

			trclass = 0;

			if (idsmall == 1){
				for (h=0; h<n; h++){				
					if (cn[id1[h]] == id2[h]) trclass++;
				}
			} else {
				for (h=0; h<n; h++){				
					if (cn[id2[h]] == id1[h]) trclass++;
				}
			}

			curr = (double)trclass / n;

			if ((*maxPro) < curr) (*maxPro) = curr;
			
//			printf("    %lf\n", curr);

//			####################################

			ind++;
			flag = 0;
			sch = sch - 1;
			i = i - 1;
			j = cn[sch-1];
			sch = sch - 1;
      
			for (v=0; v<size; v++){
				for (w=0; w<size; w++){
					pat[v][w] = 0;
				}
			}

			for (v=0; v<sch; v++){
				for (w=0; w<size; w++){
					pat[v][w] = 1;
					pat[w][cn[v]] = 1;
				}
			}    
      
		}



		if ((j == (size-1)) & (flag == 0)){
			i = i - 1;
			
			sch = sch - 1;

			if (sch >= 0){

				j = cn[sch];

				for (v=0; v<size; v++){
					for (w=0; w<size; w++){
						pat[v][w] = 0;
					}
				}

				if (sch > 0){
					for (v=0; v<sch; v++){
						for (w=0; w<size; w++){
							pat[v][w] = 1;
							pat[w][cn[v]] = 1;
						}
					}    
				}

			}

			if (i >= 0){
				pat[i][j] = 1;
			}

		}

		if (sch == -1){
			finish = 1;
		}

	}

	FREE_MATRIX(pat);
	FREE_VECTOR(cn);

}


void VIindex(int n, int K1, int K2, int *id1, int *id2, double *VI){
	
	int k, i;
	double H1, H2, H12;
	double *P1, *P2, **P12;
	
	
	MAKE_VECTOR(P1, K1);
	MAKE_VECTOR(P2, K2);
	MAKE_MATRIX(P12, K1, K2);	


	anull(P1, K1);
	anull(P2, K2);
	Anull(P12, K1, K2);
	for (i=0; i<n; i++){
		P1[id1[i]]++;
		P2[id2[i]]++;
		P12[id1[i]][id2[i]]++;
	}
		
	H1 = 0.0;
	for (k=0; k<K1; k++){
		P1[k] = P1[k] / n;
		H1 = H1 - P1[k] * log(P1[k]);
	}

	H2 = 0.0;
	for (k=0; k<K2; k++){
		P2[k] = P2[k] / n;
		H2 = H2 - P2[k] * log(P2[k]);
	}

	H12 = 0.0;
	for (k=0; k<K1; k++){
		for (i=0; i<K2; i++){
			P12[k][i] = P12[k][i] / n;
			if (P12[k][i] != 0.0) H12 = H12 + P12[k][i] * log(P12[k][i] / P1[k] / P2[i]);
		}
	}
	
	(*VI) = H1 + H2 - 2 * H12;
		
	
	FREE_VECTOR(P1);
	FREE_VECTOR(P2);
	FREE_MATRIX(P12);	
		
}






