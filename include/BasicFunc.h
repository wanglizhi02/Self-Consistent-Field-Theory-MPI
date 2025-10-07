#ifndef __BasicFunc_h
#define __BasicFunc_h

#include "Head.h"

extern void printVecCplx(fftw_complex *src);
extern double normCplxInfty(fftw_complex *src, int n);
extern double normRealInfty(double *src, int n);
extern void eliminate(fftw_complex *src, int n);

extern void FuncsLinear1Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1);

extern void FuncsLinear2Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2);

extern void FuncsLinear3Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3);

extern void FuncsLinear4Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4);

void FuncsLinear5Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4,
					  const double a5, const fftw_complex *F5);

extern void FuncAddToCplx(fftw_complex *rslt, int n,
			   const double a1, const fftw_complex *F1);

extern void FuncCplxAddAConst(fftw_complex *rslt, int n, const double a);

extern void setCplxZero(fftw_complex rslt);

extern void FuncsLinear1Real(double *rslt, int n,
					  const double a1, const double *F1);

extern void FuncsLinear2Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2);

extern void FuncsLinear3Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3);

extern void FuncsLinear4Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4);

extern void FuncsLinear5Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4,
					  const double a5, const double *F5);


extern void FuncAddToReal(double *rslt, int n,
				   const double a1, const double *F1);

extern void MatCopy(double **dst, double **src, int n, int m);
extern void MatsAdd(double **dst, double d, double **src, double s, int n, int m);
extern double MatsDot(double **lhs, double **rhs, int n, int m);
extern void MatPrint(double **matrix, int n, int m);
extern double innerProduct(double *lft, double *rht, int n);

#endif
