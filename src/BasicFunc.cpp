#include "Head.h"
#include "Data.h"
#include "BasicFunc.h"

void printVecCplx(fftw_complex *src)
{
	for(int k = 0; k < alloc_local; k++)
	{
		printf("src[%d][0] = %f \t src[%d][1] = %f\n", k, src[k][0], k, src[k][1]);
	}
}

double normCplxInfty(fftw_complex *src, int n)
{
	double tmp;
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = src[i][0]*src[i][0]+src[i][1]*src[i][1];
		tmp = sqrt(tmp);
		rslt = (rslt > tmp ? rslt : tmp);//（条件表达式）？（条件为真时的表达式）：（条件为假时的表达式）
	}
	return rslt;
}

double normRealInfty(double *src, int n)
{
	double tmp;
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
	{
		tmp = src[i];
		rslt = (rslt > tmp ? rslt : tmp);
	}
	return rslt;
}

void eliminate(fftw_complex *src, int n)
{
	for(int i = 0; i < n; i++)
	{
		if(fabs(src[i][0]) < 1.0e-15)
			src[i][0] = 0.0;
		if(fabs(src[i][1]) < 1.0e-15)
			src[i][1] = 0.0;
	}
}

void FuncsLinear1Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0];
		rslt[i][1] = a1*F1[i][1];
	}
}

void FuncsLinear2Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1];
	}
}

void FuncsLinear3Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0] + a3*F3[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1] + a3*F3[i][1];
	}
}

void FuncsLinear4Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0] + a3*F3[i][0] + a4*F4[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1] + a3*F3[i][1] + a4*F4[i][1];
	}
}

void FuncsLinear5Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4,
					  const double a5, const fftw_complex *F5)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] = a1*F1[i][0] + a2*F2[i][0] + a3*F3[i][0] + a4*F4[i][0] + a5*F5[i][0];
		rslt[i][1] = a1*F1[i][1] + a2*F2[i][1] + a3*F3[i][1] + a4*F4[i][1] + a5*F5[i][1];
	}
}

void FuncAddToCplx(fftw_complex *rslt, int n,
			   const double a1, const fftw_complex *F1)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] += a1*F1[i][0];
		rslt[i][1] += a1*F1[i][1];
	}
}

void FuncCplxAddAConst(fftw_complex *rslt, int n,
					   const double a)
{
	for(int i = 0; i < n; i++)
	{
		rslt[i][0] += a;
		rslt[i][1] += a;
	}
}

void setCplxZero(fftw_complex rslt)
{
	rslt[0] = 0.0; rslt[1] = 0.0;
}

void FuncsLinear1Real(double *rslt, int n,
					  const double a1, const double *F1)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i];
}

void FuncsLinear2Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i] + a2*F2[i];
}

void FuncsLinear3Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i] + a2*F2[i] + a3*F3[i];
}

void FuncsLinear4Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i] + a2*F2[i] + a3*F3[i] + a4*F4[i];
}

void FuncsLinear5Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4,
					  const double a5, const double *F5)
{
	for(int i = 0; i < n; i++)
		rslt[i] = a1*F1[i] + a2*F2[i] + a3*F3[i] + a4*F4[i] +a5*F5[i];
}


void FuncAddToReal(double *rslt, int n,
				   const double a1, const double *F1)
{
	for(int i = 0; i < n; i++)
		rslt[i] += a1*F1[i];
}

void MatCopy(double **dst, double **src, int n, int m)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			dst[i][j] = src[i][j];
}

void MatsAdd(double **dst, double d, double **src, double s, int n, int m)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			dst[i][j] = d*dst[i][j] + s*src[i][j];
}

double MatsDot(double **lhs, double **rhs, int n, int m)
{
	double rslt = 0;
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			rslt += lhs[i][j]*rhs[i][j];
	return rslt;
}

void MatPrint(double **matrix, int n, int m)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			printf("matrix[%d][%d] = %.15f  ", i, j, matrix[i][j]);
		}
		printf("\n");
	}
}

//void MatPrintLdouble(long double **matrix, int n, int m)
//{
//    for(int i = 0; i < n; i++)
//    {
//        for(int j = 0; j < m; j++)
//        {
//            printf("matrix[%d][%d] = %Lf   ", i, j, matrix[i][j]);
//        }
//        printf("\n");
//    }
//}

double innerProduct(double *lft, double *rht, int n)
{
	double rslt = 0.0;
	for(int i = 0; i < n; i++)
		rslt += lft[i]*rht[i];
	return rslt;
}

