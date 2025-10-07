#include "Head.h"
#include "Data.h"

int nprocs, myrank;
double chi;
int phase;
int DIM, DimPhy, DimCpt;
double dsMax;
double Dx, Dy;
int Ndeg;
int ItMax;
ptrdiff_t *NCpt;
int cplxDofs, realDofs;
int Nblock, Nspecies, Nblend;

ptrdiff_t alloc_local, local_n0, local_0_start;
int *n, *displs;
int a;
double *singQ;
double **dirBox, **rcpBox, **gradB;
double **ProjMatrix;

int **indKspace, **indKspacelocal;


double *Gsquare;
double *Gsquarelocal;

fftw_complex **rho;
double **rhoglobal, **rhoreal;

fftw_complex **fieldW, *fieldWplus;
fftw_complex **gradW;

fftw_complex *fftw_Ctmp;
double *fftw_Rtmp;


fftw_plan PlanC2R, PlanR2C;

double fA, fB;
int *Nsplit;
int NA, NB, Ns;
Interval *Ranges;

fftw_complex **frdQ, **bakQ;
fftw_complex **hatq_qplus;

double  Hamilton;
double  FreeEnergy;
double internalEnergy, entropicEnergy1, entropicEnergy2;
