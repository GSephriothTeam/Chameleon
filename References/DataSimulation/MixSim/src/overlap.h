
#ifndef OVERLAP_H
#define OVERLAP_H

int CMixSim(int p, int K, double BarOmega, double MaxOmega, int sph, int hom, double emax, double PiLow, double Lbound, double Ubound, int resN, double *pars, int lim, int g, int nf, int v, int w,	int o, char *Path, char *PIfname, char *MUfname, char *Sfname, char *overmap, char *overbarmax,	char *dataX, char *Nksizes, char *IDtrue, char *Lambda);
void ExactOverlap(int p, int K, double *Pi, double **Mu, double ***S, double *pars, int lim, double **OmegaMap, double (*BarOmega), double (*MaxOmega), double (*EigOmega), int *rcMax);
void OmegaClust(double Omega, int method, int p, int K, double PiLow, double Lbound, double Ubound, double emax, double *pars, int lim, int resN, int sph, int hom, double *Pi, double **Mu, double ***S, double **OmegaMap, double (*BarOmega), double (*MaxOmega), double (*EigOmega), int *rcMax, int (*fail));
void OmegaBarOmegaMax(int p, int K, double PiLow, double Lbound, double Ubound, double emax, double *pars, int lim, int resN, int sph, double *Pi, double **Mu, double ***S, double **OmegaMap, double (*BarOmega), double (*MaxOmega), int *rcMax, int (*fail));
void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);
int vecMin(double *x, int p, double (*min));
int vecMax(double *x, int p, double (*max));
void Anull(double **X, int ax, int bx);
void anulli(int *x, int p);
void anull(double *X, int p);
void XAXt(double **X, int p, double **A, double **Res);
void XAXt2(double **X, int p, double **A, double ***Res, int k);
void cpy(double **a, int nrows, int ncols, double **b);
void cpy1(double ***a, int k, int nrows, int ncols, double **b);
void cpy2(double **a, int nrows, int ncols, double ***b, int k);
void tA(double **A, int a, int b, double **Res);
double vecNNvec(int p, double *y, double *x);
void matxvec(double **a, int arows, int acols, double *x, int xrows, double *y);
void multiply(double **a, int arows, int acols, double **b, int brows, int bcols, double **c);

/* WCC: "libEVD.c" and "libEVD_LAPACK.c" */
#ifndef __HAVE_R_
	void cephes_symmeigens_down(int p, double *eval, double **A, double (*determinant));
#else
	void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);
	void EigValDec(int size, double *W, double **A, double (*determinant));
#endif


void cxS(int p, int K, double ***S, double c);
double qfc(double* lb1, double* nc1, int* n1, int *r1in, double *sigmain, double *c1in, int *lim1in, double *accin, double* trace, int* ifault);
void genPi(int K, double PiLow, double *Pi);
void genMu(int p, int K, double **Mu, double Lbound, double Ubound);
void genSigmaEcc(int p, int K, double emax, double ***S, int hom);
void genSphSigma(int p, int K, double ***S, int hom);

/* WCC: Delete "libPrint.c". */
#ifndef __HAVE_R_
	void fprintOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax, int fileN, char *overmap, char *overbarmax);
	void printOverlap(int K, double **OmegaMap, double BarOmega, double MaxOmega, int *rcMax);
	void fprintParameters(int p, int K, double *Pi, double **Mu, double ***S, char *PIfname, char *MUfname, char *Sfname, int fileN);
	void printParameters(int p, int K, double *Pi, double **Mu, double ***S);
	void freadParameters(int p, int K, double *Pi, double **Mu, double ***S, char *PIfname, char *MUfname, char *Sfname, int fileN);
	void fprintData(int n, int p, int K, double **x, int *Nk, int fileN, char *dataX, char *Nksizes, char *IDtrue);
#endif

void array1to2(int a, int b, double *y, double **x);
void array1to2i(int a, int b, int *y, int **x);
void array1to3(int a, int b, int c, double *y, double ***x);
void array2to1(int a, int b, double *y, double **x);
void array2to1i(int a, int b, int *y, int **x);
void array3to1(int a, int b, int c, double *y, double ***x);
void array3to1(int a, int b, int c, double *y, double ***x);
void AllPerms(int size, int **perms);
int Factorial(int a);
	
/* WCC: Delete "libRandom.c". */
#ifndef __HAVE_R_
	void setseed(unsigned int s);
	long genseed(void);
	double runir(double a,double b);
	double rnor(double mu,double sd);
	double rgamma(double alpha);
	void rmulti(int n, int K, double *pi, int *Nk);
	void rMVN(int n, int p, double *MU, double **S, double **X);
	float gammln(float xx);
	void gser(float *gamser, float a, float x, float *gln);
	float gammq(float a, float x);
#endif

void RRand(int N, int TRUK, int PREDK, int *trcl, int *prcl, double *Rand, double *adjRand, double *F);
void proAgree(int n, int K1, int K2, int *id1, int *id2, double *proMax);
void VIindex(int n, int K1, int K2, int *id1, int *id2, double *VI);

#endif /* OVERLAP_H */
