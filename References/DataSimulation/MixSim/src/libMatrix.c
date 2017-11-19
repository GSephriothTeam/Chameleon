


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "array.h"

#include "overlap.h"


/* Multiplies matrices a and b and puts the result in c which should be
 pre-allocated */

void multiply(double **a, int arows, int acols,
	      double **b, int brows, int bcols, double **c)
{
  int i, j, k;
  
  for (i=0; i<arows; i++)
    for (j=0; j<bcols; j++) {
      c[i][j] = 0;
      for (k=0; k<acols; k++)
	c[i][j] += a[i][k] * b[k][j];
    }

}

/* Multiplies matrices a and b and puts the result in c[m] which should be
 pre-allocated */

void multiply2(double **a, int arows, int acols,
	       double **b, int brows, int bcols, double ***c, int m)
{
  int i, j, k;
  
  for (i=0; i<arows; i++)
    for (j=0; j<bcols; j++) {
      c[m][i][j] = 0;
      for (k=0; k<acols; k++)
	c[m][i][j] += a[i][k] * b[k][j];
    }

}

/* Multiplies matrix a and vector x and puts the result in y which should be
 pre-allocated */

void matxvec(double **a, int arows, int acols,
		double *x, int xrows, double *y)
{
  int i, k;
  
  for (i=0; i<arows; i++){
    y[i] = 0;
    for (k=0; k<acols; k++){
      y[i] += a[i][k] * x[k];
    }
  }

}


/*copies matrix A to matrix B*/

void cpy(double **a, int nrows, int ncols, double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[i][j];
    }
  }

}


/*copies matrix A[k] to matrix B */

void cpy1(double ***a, int k, int nrows, int ncols, double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[k][i][j];
    }
  }

}

/*copies matrix A to matrix B[k] */

void cpy2(double **a, int nrows, int ncols, double ***b, int k)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[k][i][j]=a[i][j];
    }
  }

}


/* finds the smallest element in x */

int vecMin(double *x, int p, double (*min)){

	int i, minN;

	(*min) = x[0];
	minN = 0;

	for (i=0;i<p;i++){
		if (x[i] < (*min)){
			(*min) = x[i];
			minN = i;
		}
	}

	return minN;

}


/* finds the largest element in x */

int vecMax(double *x, int p, double (*max)){

	int i, maxN;

	(*max) = x[0];
	maxN = 0;

	for (i=0;i<p;i++){
		if (x[i] > (*max)){
			(*max) = x[i];
			maxN = i;
		}
	}

	return maxN;

}


/* multiplies x by x' */

int vec11vecSQ(double *y, int p, double **Res){
	
	int i,j;

	for (i=0;i<p;i++){
		for (j=0;j<p;j++){
			Res[i][j] = y[i]*y[j];
		}
	}
	
	return 0;
}


/* multiplies x' by y */

double vecNNvec(int p, double *y, double *x){
	
	int i;
	double Res;

	Res = 0;
	for (i=0;i<p;i++){
		Res = Res + y[i]*x[i];
	}

	return Res;
}


/* subtracts matrices */

int mat_(int a, int b,double **Res, double **Y){
	
	int i,j;

	for (i=0;i<a;i++){
		for (j=0;j<b;j++){
			Res[i][j] = Res[i][j] - Y[i][j];
		}
	}
	
	return 0;
}


/* finds sums of rows in matrix */ 

int vecsum(int a, int b,double **OO, double *Res){
	
	int i,j;

	for (i=0;i<a;i++){
		Res[i] = 0;
		for (j=0;j<b;j++){
			Res[i] = Res[i] + OO[i][j];
		}
	}
	
	return 0;
}


/* multiplies OO by OO' */

int MatrixProd(double **OO, int p, int m, double **Res){
     
     int i,j,k;

     for (i=0; i<p; i++){
         for (j=0; j<p; j++){
             Res[i][j]=0;
             for (k=0; k<m; k++){
                 Res[i][j]=Res[i][j]+OO[i][k]*OO[j][k];
             }
         }     
     }
     
     return 0;
}


/* computes Kronecker product */

int Kronecker(double **A, int a1, int a2, double **B, int b1, int b2, double **Res){

  int inda1, inda2, indb1, indb2, indRes1, indRes2;
  int i;
  int n;


  n = a1 * b1 * a2 * b2;

  indRes1 = 0;
  indRes2 = -1;

  inda1 = 0;
  inda2 = 0;
  indb1 = 0;
  indb2 = -1;

  for (i=0; i<n; i++){

    indb2++;
    indRes2++;

    if (indb2 == b2){

      indb2 = 0;
      inda2++;

      if (inda2 == a2){
	
	inda2 = 0;
	indb1++;
	indRes1++;
	indRes2 = 0;

	if (indb1 == b1){

	  indb1 = 0;
	  inda1++;

	}

      }
    }

    Res[indRes1][indRes2] = A[inda1][inda2] * B[indb1][indb2];

  }

  return 0;

}


/* computes G matrix */

int Gmat(int p, int m, double **Res){
     
     int a,b,i,i1,i2,n,ind;

	 n = 0;

     for (a=0; a<p; a++){
         for (b=0; b<p; b++){
         	
         	if (a < b){
         		i1 = b;
         		i2 = a;
         	} else {
         		i1 = a;
         		i2 = b;
         	}
         	
         	ind = m - (p - i2) * (p - i2 + 1) / 2 + i1 - i2;
         	
         	for (i=0; i<m; i++){
         	
         		if (i != ind ){
         			Res[n][i] = 0;	
         		} else {
         			Res[n][i] = 1;
         		}
         		
         	}
         	
         	n++;

         }     
     }
     
     return 0;
}


/* provides transpose */

void tA(double **A, int a, int b, double **Res){

	int i,j;

   	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			Res[i][j] = A[j][i];
		}
	}
	
}


/* computes product of three matrices */

int ZXY(double **Z, int az, int bz, double **X, int ax, int bx, double **Y, int ay, int by, double **Res){

	double **Res1;

	MAKE_MATRIX(Res1, az, bx);	

	multiply(Z, az, bz, X, ax, bx, Res1);
	multiply(Res1, az, bx, Y, ay, by, Res);

	FREE_MATRIX(Res1);
 
        return 0;
    
}


/* Computes X %*% A %*% t(X) */

void XAXt(double **X, int p, double **A, double **Res){

	double **Res1, **Res2;

	MAKE_MATRIX(Res1, p, p);	
	MAKE_MATRIX(Res2, p, p);

	tA(X, p, p, Res2);

	multiply(X, p, p, A, p, p, Res1);
	multiply(Res1, p, p, Res2, p, p, Res);

	FREE_MATRIX(Res1);
 	FREE_MATRIX(Res2);

}


/* Computes X %*% A %*% t(X) and writes results into Res[k,,] */

void XAXt2(double **X, int p, double **A, double ***Res, int k){

	double **Res1, **Res2;

	MAKE_MATRIX(Res1, p, p);	
	MAKE_MATRIX(Res2, p, p);

	tA(X, p, p, Res2);

	multiply(X, p, p, A, p, p, Res1);
	multiply2(Res1, p, p, Res2, p, p, Res, k);

	FREE_MATRIX(Res1);
 	FREE_MATRIX(Res2);

}

/* provides 0-matrix */

void Anull(double **X, int ax, int bx){
     
     int i, j;

     for (i=0; i<ax; i++){
         for (j=0; j<bx; j++){
		X[i][j] = 0.0;
	 }
     }
}


/* provides 0-vector */

void anull(double *x, int p){
     
     int i;

     for (i=0; i<p; i++){
	     x[i] = 0.0;
     }
}


/* provides 0-matrix of integers */

void Anulli(int **X, int ax, int bx){
     
     int i, j;

     for (i=0; i<ax; i++){
         for (j=0; j<bx; j++){
		X[i][j] = 0;
	 }
     }
}


/* provides 0-vector of integers */

void anulli(int *x, int p){
     
     int i;

     for (i=0; i<p; i++){
	     x[i] = 0;
     }
}


/* transforms matrix into vector */

int asvector(double **X, int ax, int bx, double *ResVec){
     
    int i,j,k;

    k = 0;

    for (i=0; i<ax; i++){
         for (j=0; j<bx; j++){
		 	ResVec[k] = X[i][j];
		 	k++;
         }
     }

    return 0;
    
}


/* multiplies 3d array by constant  */

void cxS(int p, int K, double ***S, double c){

	int i, j, k;

	for (k=0; k<K; k++){
		for (i=0; i<p; i++){
			for (j=0; j<p; j++){
				S[k][i][j] = c * S[k][i][j];
			}
		}
	}

}



void array1to2(int a, int b, double *y, double **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			x[i][j] = y[k];
			k++;

		}
	}
	
	
}

void array1to2i(int a, int b, int *y, int **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			x[i][j] = y[k];
			k++;

		}
	}
	
	
}

void array1to3(int a, int b, int c, double *y, double ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				x[i][j][m] = y[k];
				k++;
			
			}
		}
	}
	
	
}


void array2to1(int a, int b, double *y, double **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			y[k] = x[i][j];
			k++;

		}
	}
	
	
}

void array2to1i(int a, int b, int *y, int **x){

	int i, j, k;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){

			y[k] = x[i][j];
			k++;

		}
	}
	
	
}

void array3to1(int a, int b, int c, double *y, double ***x){

	int i, j, k, m;

	k = 0;
	for (i=0; i<a; i++){
		for (j=0; j<b; j++){
			for (m=0; m<c; m++){

				y[k] = x[i][j][m];
				k++;
			
			}
		}
	}
	
	
}


void AllPerms(int size,int **perms){

	int sch, i, j, v, w, finish, flag, ind;
	double **pat;
	int *cn;

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

	MAKE_VECTOR(cn, size);
	for (v=0; v<size; v++){
		cn[v] = 0;
	}
  

	while (finish == 0){
    
		if (j != (size-1)){
			j = j+1;
		} else {
			if (flag == 1){
				j = 0;
				i = i+1;
				flag = 0;
			}
		}
    
		if (pat[i][j] == 0){
			for (v=0; v<size; v++){
				pat[i][v]=1;
				pat[v][j]=1;
			}
      
			sch = sch + 1;
			cn[sch-1] = j;
			flag = 1;
		}

		if ((sch == size) & (flag == 1)){
      
			for (v=0; v<size; v++){
				perms[ind][v] = cn[v];
			}

			ind++;
			flag = 0;
			sch = sch - 1;
			i = i - 1;
			j = cn[sch-1];
			sch = sch-1;
      
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


		if ((j == size - 1) & (flag == 0)){
			i = i - 1;
			
			sch = sch-1;

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


int Factorial(int a){
    int i;
    int res;
    
    res=1;
    for (i=1; i<(a+1); i++){
        res=res*i;
    }
    
    return res;
}

