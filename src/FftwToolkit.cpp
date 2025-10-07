#include "Head.h"
#include "Data.h"
#include "BasicFunc.h"
#include "FftwToolkit.h"
#include "mpi.h"                                                                
#include "fftw3-mpi.h" 

void getIndextheta(int *kspace, int n,  int ndeg)
{
	 int k1 = 0;	

	for(int i = 0; i < n; i++)
	{
			if(k1 > ndeg/2)
				kspace[i] = k1 - ndeg;
			else
				kspace[i] = k1;
		k1 ++;
	}
}


void getIndex(int **kspace, int n, int dim, ptrdiff_t *ndeg)
{
	int *k = (int *)malloc(sizeof(int)*dim);
	for(int i = 0; i < dim; i++) k[i] = 0;	

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < dim-1; j++)
		{
			if(k[j] > ndeg[j]/2)
				kspace[i][j] = k[j] - ndeg[j];
			else
				kspace[i][j] = k[j];
		}
		kspace[i][dim-1] = k[dim-1];
		k[dim-1] ++;
		if(k[dim-1] > ndeg[dim-1]/2)
		{
			k[dim-1] = 0;
			k[dim-2] ++;
		}
		for(int jj = dim-2; jj > 0; jj--)
		{
			if(k[jj] > ndeg[jj]-1)
			{
				k[jj] = 0;
				k[jj-1] ++;
			}
		}
	}
	free(k);
}
void FftwC2R(double *Rrslt, fftw_complex *Corig)
{
	memcpy(fftw_Ctmp, Corig, sizeof(fftw_complex)*alloc_local);
	fftw_execute(PlanC2R);
	memcpy(Rrslt, fftw_Rtmp, sizeof(double)*(2*alloc_local));
}

void FftwR2C(fftw_complex *Crslt, double *Rorig)
{
	memcpy(fftw_Rtmp, Rorig, sizeof(double)*(2*alloc_local));
	fftw_execute(PlanR2C);
	FuncsLinear1Cplx(fftw_Ctmp, alloc_local, 1.0/realDofs, fftw_Ctmp);
	memcpy(Crslt, fftw_Ctmp, sizeof(fftw_complex)*alloc_local);
}

void hatConv(fftw_complex *rslt, fftw_complex *src1, fftw_complex *src2)
{
	double *Rsrc1 = (double *)malloc(sizeof(double) *(2*alloc_local));
	double *Rsrc2 = (double *)malloc(sizeof(double) *(2*alloc_local));

	FftwC2R(Rsrc1, src1);
	FftwC2R(Rsrc2, src2);

	for(int i = 0; i < 2*alloc_local; i++)
        Rsrc1[i] *= Rsrc2[i];
	FftwR2C(rslt, Rsrc1);

	free(Rsrc1);
	free(Rsrc2);
}

void Conv(fftw_complex *rslt, double *src1, double *src2)
{
	for(int i = 0; i < 2*alloc_local; i++)
        src1[i] *= src2[i];
	FftwR2C(rslt, src1);
}


double Intergral_space(fftw_complex *src1, fftw_complex *src2)
{
	double val = 0.0;
	double *Rsrc1 = (double *)malloc(sizeof(double) *(2*alloc_local));
	double *Rsrc2 = (double *)malloc(sizeof(double) *(2*alloc_local));
	fftw_complex *rslt = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *alloc_local);

	FftwC2R(Rsrc1, src1);
	FftwC2R(Rsrc2, src2);

	for(int i = 0; i < 2*alloc_local; i++)
        	Rsrc1[i] *= Rsrc2[i];
	FftwR2C(rslt, Rsrc1);
	val = rslt[0][0];
	fftw_free(rslt);
	free(Rsrc1);
	free(Rsrc2);
	return val;
}

double Intergral_space_real(double *src1, double *src2)
{
	double val = 0.0;
	fftw_complex *rslt = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *alloc_local);
	double *Rsrc = (double *)malloc(sizeof(double) *2*alloc_local);

	for(int i = 0; i < 2*alloc_local; i++)
        	Rsrc[i] =src1[i] * src2[i];
	
	FftwR2C(rslt, Rsrc);
	
	val = rslt[0][0];
	fftw_free(rslt);
	free(Rsrc);
	return val;
}


double innerProduct_space(fftw_complex *lft, fftw_complex *rht)
{
	double *Rsrc1 = (double *)malloc(sizeof(double) *(2*alloc_local));
	double *Rsrc2 = (double *)malloc(sizeof(double) *(2*alloc_local));
	
	FftwC2R(Rsrc1, lft);
	FftwC2R(Rsrc2, rht);


	double rslt = 0.0;
	for(int i = 0; i < (2*alloc_local); i++)
		rslt += Rsrc1[i]*Rsrc2[i];
	return rslt;
}

/*
double Intergral_u(double *src1, double *src2)
{
	double val = 0.0;
	double *Rsrc1 = (double *)malloc(sizeof(double) *Ntheta);
	double *Rsrc2 = (double *)malloc(sizeof(double) *Ntheta);

	fftw_complex *rslt = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) *cplxNtheta);
	memcpy(Rsrc1, src1, sizeof(double) *Ntheta);
	memcpy(Rsrc2, src2, sizeof(double) *Ntheta);
	for(int i = 0; i < Ntheta; i++)
        	Rsrc1[i] *= Rsrc2[i];
	FftwR1C(rslt, Rsrc1);

	val = rslt[0][0]*2*PI;
	free(Rsrc1);
	free(Rsrc2);
	fftw_free(rslt);
	return val;
}
*/
double hatVecsMultiply(fftw_complex *lhs, fftw_complex *rhs)
{
	double *Rlhs = (double *)malloc(sizeof(double) *realDofs);
	double *Rrhs = (double *)malloc(sizeof(double) *realDofs);
	
	FftwC2R(Rlhs, lhs);
	FftwC2R(Rrhs, rhs);
	double rslt = 0.0;
	for(int i = 0; i < realDofs; i++)
		rslt += Rlhs[i]*Rrhs[i];		

	free(Rlhs);
	free(Rrhs);
	return rslt;
}

double InnerProd(fftw_complex *lhs, fftw_complex *rhs, int cplxDofs, int n)
{
	double rslt = 0;

    	int tmp = n/2+1;
	for(int i = 0; i < cplxDofs; i++)
        	rslt += ((i%tmp == 0) || ((i+1)%tmp == 0)) ? lhs[i][0]*rhs[i][0]+lhs[i][1]*rhs[i][1] : 2*(lhs[i][0]*rhs[i][0]+lhs[i][1]*rhs[i][1]);

	return rslt;
}

