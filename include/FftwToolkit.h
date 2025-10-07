#ifndef __FftwToolkit_h
#define __FftwToolkit_h

#include "Head.h"

//extern void getIndextheta(int *kspace, int n, int ndeg);
extern void getIndex(int **kspace, int n, int dim, ptrdiff_t *ndeg);
extern void FftwC2R(double *Rrslt, fftw_complex *Corig);
extern void FftwR2C(fftw_complex *Crslt, double *Rorig);
//extern void FftwC1R(double *Rrslt, fftw_complex *Corig);
//extern void FftwR1C(fftw_complex *Crslt, double *Rorig);
extern double hatVecsMultiply(fftw_complex *lhs, fftw_complex *rhs);
extern double InnerProd(fftw_complex *lhs, fftw_complex *rhs, int cplxNtheta, int n);
extern void hatConv(fftw_complex *rslt, fftw_complex *src1, fftw_complex *src2);
extern void Conv(fftw_complex *rslt, double *src1, double *src2);
extern double Intergral_space(fftw_complex *src1, fftw_complex *src2);
extern double Intergral_space_real(double *src1, double *src2);
//extern double Intergral_u(double *src1, double *src2);

extern double innerProduct_space(fftw_complex *lft, fftw_complex *rht);
#endif
