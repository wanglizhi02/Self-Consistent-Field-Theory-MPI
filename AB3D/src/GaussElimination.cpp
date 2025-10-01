#include "Head.h"
#include "GaussElimination.h"
#include "fftw3-mpi.h" 

#include "Data.h"
#include "Mytimer.h"
#include "Initialization.h"
#include "BasicFunc.h"
#include "MDESolver.h"
#include "FftwToolkit.h"
#include "IntegrationToolkit.h"
#include "AndersonMix.h"
#include "SCFTBaseAB.h"

#include "config.h"


/*
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_block.h>
*/                          
void GaussElimination(double *u, double **A_matrix, double *rhs, int n)
{
	for(int s=0; s<n-1; s++)
	{
		CheckIsHaveSolution(A_matrix, n, s);
		Pivot(A_matrix, rhs, n, s);
		for(int i=s+1; i<n; i++)
		{
			double c=A_matrix[i][s]/A_matrix[s][s];
			rhs[i]-=rhs[s]*c;
			for(int j=s+1; j<n; j++)
			{
				A_matrix[i][j]-=A_matrix[s][j]*c;
			}
		}
	}
//    for(int i = 0; i < n; i++)
//    {
//        for(int j = 0; j < n; j ++)
//        {
//            printf("A_MATRIX[%d][%d] = %e\t ", i, j , A_matrix[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
	////////////////////////////////back replace;
	int Vsize=n;
	u[Vsize-1]=rhs[Vsize-1]/A_matrix[Vsize-1][Vsize-1];
	for(int i=Vsize-2; i>=0; i--)
	{
		double S=0.0;
		for(int j=i+1; j<Vsize; j++)
		{
			S+=A_matrix[i][j]*u[j];
		}
		u[i]=(rhs[i]-S)/A_matrix[i][i];
	}
}

void NumberSwap(double *a, double *b)
{
	double tmp = *a;
	*a = *b;
	*b = tmp;
}

void CheckIsHaveSolution(double **A_matrix, int n, int k)
{
	int i=0;
	for(i=k; i<n; i++)
	{
		if(A_matrix[i][k] != 0.0)
			break;
	}
	if(i==n)
	{
		std::cout<<"No solution or no unique solution!"<<std::endl;
		std::cout<<"Please input CONTRAL+C to stop!"<<std::endl;
		getchar();
	} 
}

void Pivot(double **A_matrix, double *Rhs, int n, int k)
{
	int t=k;
	for(int j=k; j<n; j++)
	{
		if(fabs(A_matrix[j][k])>fabs(A_matrix[t][k]))
		{
			t=j;
		}
	}
	if(t!=k)
	{
		NumberSwap(Rhs+k, Rhs+t);
		for(int j=k; j<n; j++)
		{
			NumberSwap(*(A_matrix+k)+j, *(A_matrix+t)+j);
		}
	}
}
/*
void Solve_LU(double *u, double **A_matrix, double *rhs, int n)
{
    double *U = (double *)malloc(sizeof(double)*(n*n));
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            U[i*n+j] = A_matrix[i][j];
        }
    }
    gsl_matrix_view A;
    gsl_vector_view b, x;
    gsl_permutation *p;
    int s;

    A = gsl_matrix_view_array(U, n, n);
    b = gsl_vector_view_array(rhs, n);
    x = gsl_vector_view_array(u, n);

    p = gsl_permutation_alloc(n);

    gsl_linalg_LU_decomp(&A.matrix, p, &s);
    gsl_linalg_LU_solve(&A.matrix, p, &b.vector, &x.vector);

    gsl_permutation_free(p);
}
*/
void Solve_LU(double *u, double **A_matrix, double *rhs, int n)
// 实现对矩阵A的LU分解，L为下三角矩阵
{
    double **L = (double **)malloc(sizeof(double*)*n);
    for(int i=0; i<n; i++)
		L[i] = (double *)malloc(sizeof(double)*n);

    double **U = (double **)malloc(sizeof(double*)*n);
    for(int i=0; i<n; i++)
		U[i] = (double *)malloc(sizeof(double)*n);

	double *V = (double *)malloc(sizeof(double)*n);
    
	for(int i=0; i<n; i++)
	{
		 L[i][i] = 1.0;
	}
    
	for(int k=0; k<n; k++)
	{
		for(int j=k; j<n; j++)
		{
			U[k][j] = A_matrix[k][j]; 
			for(int m=0; m<k-1; m++)
				U[k][j] -= L[k][m]*U[m][j];
		}
		for(int i=k+1; i<n; i++)
		{
			L[i][k] = A_matrix[i][k];
			for(int m=0; m<k-1; m++)
				L[i][k] -= L[i][m]* U[m][k];
			
			L[i][k] = L[i][k]/U[k][k];
		}
	}

	V[0] = rhs[0];
    for(int i=1; i<n; i++)
	{
		for(int j=0; i<i-1; i++)
		{
			rhs[i] -=- L[i][j]*V[j];
		}
	 V[i] = rhs[i];
	}
	u[n-1]=V[n-1]/U[n-1][n-1];

    for(int i=n-2; i>=0; i--)
	{
		for(int j=n-1; j>=i+1; j--)
		{
			V[i] = V[i] - U[i][j] * u[j];
		}
     u[i] = V[i]/U[i][i];
	}

    for(int i=0; i<n; i++)
	{
		free(L[i]);
		free(U[i]);
	}
	
		free(L);
		free(U);
		free(V);

}
	
