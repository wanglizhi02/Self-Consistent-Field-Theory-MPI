#include "Data.h"
#include "IntegrationToolkit.h"

void integration(fftw_complex *phi, 
				 Interval *range, 
				 fftw_complex **integrand,
				 int *index,
				 double singQ)
{
	int n = range->n;
	double ds = (range->b - range->a) / n;
	for(int k = 0; k < alloc_local; k++)
	{
		phi[k][0]  = -0.625 * (integrand[*index][k][0]+integrand[*index+n][k][0]);
		phi[k][1]  = -0.625 * (integrand[*index][k][1]+integrand[*index+n][k][1]);
		phi[k][0] += 1.0/6.0 * (integrand[*index+1][k][0]+integrand[*index+n-1][k][0]);
		phi[k][1] += 1.0/6.0 * (integrand[*index+1][k][1]+integrand[*index+n-1][k][1]);
		phi[k][0] -= 1.0/24.0 * (integrand[*index+2][k][0]+integrand[*index+n-2][k][0]);
		phi[k][1] -= 1.0/24.0 * (integrand[*index+2][k][1]+integrand[*index+n-2][k][1]);
		for(int i = 0; i < n+1; i++)
		{
			phi[k][0] += integrand[*index+i][k][0]; 
			phi[k][1] += integrand[*index+i][k][1]; 
		}
		phi[k][0] *= ds;
		phi[k][1] *= ds;
		phi[k][0] /= singQ;
		phi[k][1] /= singQ;
	}
	*index += n;
}

