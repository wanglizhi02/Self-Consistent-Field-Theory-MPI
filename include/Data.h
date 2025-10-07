#ifndef __Data_h
#define __Data_h

#include "Head.h"

typedef struct{
	double a, b;
	int n;
} Interval;

typedef struct myVec
{
	double* Data;
	int* Index;
}mySortVec;

extern int nprocs, myrank;
extern double chi;
extern int phase;
extern int ItMax;
extern int DIM, DimPhy, DimCpt;
extern double dsMax;
extern double Dx, Dy;
extern int Ndeg;
extern ptrdiff_t *NCpt;
extern int cplxDofs, realDofs;
extern int Nblock, Nspecies, Nblend;
extern ptrdiff_t alloc_local, local_n0, local_0_start;
extern int *n, *displs;
extern int a;
extern double *singQ;
extern double **dirBox, **rcpBox, **gradB;
extern double **ProjMatrix;
extern int **indKspace, **indKspacelocal;
extern double *Gsquare;
extern double *Gsquarelocal;
extern fftw_complex **rho;
extern double **rhoglobal, **rhoreal;
extern fftw_complex **fieldW, *fieldWplus;
extern fftw_complex **gradW;
extern fftw_complex *fftw_Ctmp;
extern double *fftw_Rtmp;
extern fftw_plan PlanC2R, PlanR2C;
extern double fA, fB;

extern int *Nsplit;
extern int NA, NB, Ns;
extern Interval *Ranges;
extern fftw_complex **frdQ, **bakQ;
extern fftw_complex **hatq_qplus;
extern double  Hamilton;
extern double  FreeEnergy;
extern double internalEnergy, entropicEnergy1, entropicEnergy2; 
#endif
