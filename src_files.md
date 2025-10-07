\subsubsection{AndersonMix.cpp}
\begin{lstlisting}
#include "mpi.h"
#include "fftw3-mpi.h"


#include "Head.h"
#include "Data.h"
#include "Mytimer.h"
#include "Initialization.h"
#include "BasicFunc.h"
#include "MDESolver.h"
#include "FftwToolkit.h"
#include "IntegrationToolkit.h"
#include "SCFTBaseAB.h"
#include "AndersonMix.h"
#include "GaussElimination.h"
#include "config.h"
#include <string>
void updateMatrix(double **srcMatrix, fftw_complex ***dw, fftw_complex ***dM, double ***dL,  int n)
{
	for(int i = n-2; i >= 0; i--)
	{
		for(int j = n-2; j >= 0; j--)
		{
			srcMatrix[i+1][j+1] = srcMatrix[i][j];
		}
	}
	for(int i = 0; i < n; i++)
	{
		srcMatrix[0][i] = 0.0;
		for(int j = 0; j < Nspecies; j++)
		{
			srcMatrix[0][i] += Intergral_space(dw[j][0], dw[j][i]);
		}
		for(int j = 0; j < 4; j++)
		{
			srcMatrix[0][i] += Intergral_space(dM[j][0], dM[j][i]);
		}
		for(int j = 0; j < DIM; j++)
		{
			srcMatrix[0][i] += dL[0][j][j]*dL[i][j][j];
		}
		srcMatrix[i][0] = srcMatrix[0][i];
	}
	if(myrank==0)
	{
		printf("\t\t ----- assemble Matrix order = %d ------\n", n);
		MatPrint(srcMatrix, n, n);
		printf("\t\t ----- assemble Matrix order = %d ------\n", n);
	}

}

void assembAnderMatrix(double **srcMatrix, double **rsltMatrix, double *rigCoeff, int n)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = i; j < n; j++)
		{
			rsltMatrix[i][j] = srcMatrix[0][0] - srcMatrix[i+1][0] - srcMatrix[0][j+1] + srcMatrix[i+1][j+1];	
		}
		rigCoeff[i] = srcMatrix[0][0] - srcMatrix[0][i+1];
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < i; j++)
		{
			rsltMatrix[i][j] = rsltMatrix[j][i];
		}
	}
	if (myrank==0)
	{
		printf("\t\t ----- U ------\n");
		MatPrint(rsltMatrix, n, n);
		printf("\n");
	}

}

void getAnderCoeff(double *u, double **A_matrix, double *rhs, int n)
{
	Solve_LU(u, A_matrix, rhs, n);
	for(int i = 0; i < n; i++)
		if (myrank==0)
        printf("rank=%d, C[%d]=%lf\n", myrank, i, u[i]); 
}

double andersonMixing(double tol)
{
    mytimer_t timer;
    mytimer_t timerGauss;
	double hamDiff = 10.0;
	double hamiltonNew = 1.0;
	double hamiltonOld = 1.0;
	double *resGradmu = (double *)malloc(sizeof(double)*Nspecies);
	double *resGradM = (double *)malloc(sizeof(double)*4);
	for(int i = 0; i < Nspecies; i++) 
	{
		resGradmu[i] = 0.0;
	}
	for(int i = 0; i < 4; i++) 
		resGradM[i] = 0.0;
	double resInfty = 1.0;
	double resmuinfty = 1.0;
	double resMinfty = 1.0;
	int nr;

    std::vector<double> res(Nspecies);
	double error = 1.0;
	int iterator = 0;

//~~~~~~~~~~~~~~~~~~~~~~~Memeory allocation~~~~~~~~~~~~~~~~~~~~~~~~~
	fftw_complex ***fieldmuold;  // 场mu的历史数据, 复空间的数据
	fieldmuold = (fftw_complex ***) fftw_malloc(sizeof(fftw_complex**) *Nspecies);
	for(int i = 0; i < Nspecies; i++)
	{
		fieldmuold[i] = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *(anderM+1));
		for(int j = 0; j < anderM+1; j++)
		{
			fieldmuold[i][j] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(alloc_local));
		}
	}
	
	fftw_complex ***gradmuold;  // 梯度W的历史数据, 复空间的数据
	gradmuold = (fftw_complex ***) fftw_malloc(sizeof(fftw_complex**) *Nspecies);
	for(int i = 0; i < Nspecies; i++)
	{
		gradmuold[i] = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *(anderM+1));
		for(int j = 0; j < anderM+1; j++)
		{
			gradmuold[i][j] = (fftw_complex *)malloc(sizeof(fftw_complex)*(alloc_local));
		}
	}
	

	fftw_complex **dmu, **mutmp, **gradmutmp;
	dmu = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *Nspecies);
	mutmp = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *Nspecies);
	gradmutmp = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *Nspecies);
	for(int i = 0; i < Nspecies; i++)
	{
		dmu[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex*) *alloc_local);
		mutmp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex*) *alloc_local);
		gradmutmp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex*) *alloc_local);
	}
	

	fftw_complex ***hatMold;  // 场mu的历史数据, 复空间的数据
	hatMold = (fftw_complex ***) fftw_malloc(sizeof(fftw_complex**) *4);
	for(int i = 0; i < 4; i++)
	{
		hatMold[i] = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *(anderM+1));
		for(int j = 0; j < anderM+1; j++)
		{
			hatMold[i][j] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*(alloc_local));
		}
	}
	
	fftw_complex ***gradMold;  // M梯度的历史数据, 复空间的数据
	gradMold = (fftw_complex ***) fftw_malloc(sizeof(fftw_complex**) *4);
	for(int i = 0; i < 4; i++)
	{
		gradMold[i] = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *(anderM+1));
		for(int j = 0; j < anderM+1; j++)
		{
			gradMold[i][j] = (fftw_complex *)malloc(sizeof(fftw_complex)*(alloc_local));
		}
	}
	
	fftw_complex **dM, **hatMtmp, **gradMtmp;  // M梯度的历史数据, 复空间的数据
	dM = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *4);
	hatMtmp = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *4);
	gradMtmp = (fftw_complex **) fftw_malloc(sizeof(fftw_complex*) *4);
	for(int i = 0; i < 4; i++)
	{
		dM[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *alloc_local);
		hatMtmp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *alloc_local);
		gradMtmp[i] = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *alloc_local);
	}
	
	double ***Bold;  // 计算区域的历史信息
	Bold = (double ***) malloc(sizeof(double**) *(anderM+1));
	for(int i = 0; i < anderM+1; i++)
	{
		Bold[i] = (double **) malloc(sizeof(double*) *DIM);
		for(int j = 0; j < DIM; j++)
		{
			Bold[i][j] = (double *) malloc(sizeof(double) *DIM);
		}
	}

	double ***gradBold;  // 计算区域的历史信息
	gradBold = (double ***) malloc(sizeof(double**) *anderM+1);
	for(int i = 0; i < anderM+1; i++)
	{
		gradBold[i] = (double **) malloc(sizeof(double*) *DIM);
		for(int j = 0; j < DIM; j++)
		{
			gradBold[i][j] = (double *) malloc(sizeof(double) *DIM);
		}
	}
	
	double **dB, **Btmp, **gradBtmp;  // 计算区域的历史信息
	dB = (double **) malloc(sizeof(double*) *DIM);
	Btmp = (double **) malloc(sizeof(double*) *DIM);
	gradBtmp = (double **) malloc(sizeof(double*) *DIM);
	for(int i = 0; i < DIM; i++)
	{
			dB[i] = (double *) malloc(sizeof(double) *DIM);
			Btmp[i] = (double *) malloc(sizeof(double) *DIM);
			gradBtmp[i] = (double *) malloc(sizeof(double) *DIM);
	}


	double **assembMatrix = (double **)malloc(sizeof(double* ) *(anderM+1));
	for(int i = 0; i < anderM+1; i++)
	{
		assembMatrix[i] = (double *)malloc(sizeof(double)*(anderM+1));
	}
	for(int i = 0; i < anderM+1; i++)
		for(int j = 0; j < anderM+1; j++)
			assembMatrix[i][j] = 0.0;

	double **LestMatrix = (double **)malloc(sizeof(double*) *anderM);
	for(int i = 0; i < anderM; i++)
		LestMatrix[i] = (double *)malloc(sizeof(double)*anderM);
	for(int i = 0; i < anderM; i++)
		for(int j = 0; j < anderM; j++)
			LestMatrix[i][j] = 0.0;


	double *rigTerm = (double *)malloc(sizeof(double)*anderM);
	double *coeff = (double *)malloc(sizeof(double)*anderM);
	for(int i = 0; i < anderM; i++) rigTerm[i] = 0.0;
	for(int i = 0; i < anderM; i++) coeff[i] = 0.0;
    
	double **rhoReallocal, **RhoReal, **MReal;
	rhoReallocal = (double **)malloc(sizeof(double*) *Nspecies);
	RhoReal = (double **)malloc(sizeof(double*) *Nspecies);                          
	MReal = (double **)malloc(sizeof(double*) *4);                          
    for(int i = 0; i < Nspecies; i++)                                           
    {                                                                           
        rhoReallocal[i] = (double *)malloc(sizeof(double)*(2*alloc_local));               
        RhoReal[i] = (double *)malloc(sizeof(double)*(2*alloc_local*nprocs));               
    }
	for(int i = 0; i < 4; i++)                                           
    {                                                                           
		MReal[i] = (double *)malloc(sizeof(double)*(2*alloc_local*nprocs));               
    }
	////一些准备量
    getUvalues();	
	getGsquare();
	getUsquare();
	getGK();
	getuumI2();

	do{
		iterator++;
		nr = std::min(iterator, anderM);
		
		if (myrank==0)
			printf("iter=%d\n", iterator); 
		
		MutoWTransform(fieldW, fieldmu);
		
		timer.reset();
		timer.start();
		updatePropagator();
		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();

		if (myrank==0)
			printf("\t\t time cost of updatePropagator : %f seconds\n", timer.get_current_time());

		timer.reset();
		timer.start();
		updateQ();
		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();
		if (myrank==0)
			printf("\t\t time cost of updateQ          : %f seconds\n", timer.get_current_time());
		
		timer.reset();
		timer.start();
		updateOrderParameter();
		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();
		if (myrank==0)
			printf("\t\t time cost of updateOrderParameter    : %f seconds\n", timer.get_current_time());
	
		timer.reset();
		timer.start();
		hamiltonNew=updateHamilton();
		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();
		if (myrank==0)
			printf("\t\t time cost of updateHamilton   : %f seconds\n", timer.get_current_time());

		hamDiff = fabs(hamiltonNew - hamiltonOld);
		hamiltonOld = hamiltonNew;

		if (myrank==0)
		{
			printf("Step %d: Energy = %.15e, entropicEnergy = %.15e, internalEnergy = %.15e\n", iterator, hamiltonNew, entropicEnergy, internalEnergy);
			printf("\tQ = %.15e,  error = %.15e, diffham = %.15e\n ", singQ[0], error, hamDiff);
		}
		
	
		///存储所有历史的场和区域
		for(int i = 0; i < Nspecies; i++)
		{
			for(int j = anderM; j > 0; j--)
			{
				memcpy(fieldmuold[i][j], fieldmuold[i][j-1], sizeof(fftw_complex)*(alloc_local));
			}
		memcpy(fieldmuold[i][0], fieldmu[i], sizeof(fftw_complex)*(alloc_local));
		}
		for(int i = 0; i < 4; i++)
		{
			for(int j = anderM; j > 0; j--)
			{
				memcpy(hatMold[i][j], hatMold[i][j-1], sizeof(fftw_complex)*(alloc_local));
			}
		memcpy(hatMold[i][0], hatM[i], sizeof(fftw_complex)*(alloc_local));
		}
		for(int i = anderM; i > 0; i--)
		{
			MatCopy(Bold[i], Bold[i-1], DimPhy, DimCpt);
		}
		MatCopy(Bold[0], rcpBox, DimPhy, DimCpt);
	
       	timer.reset();
		timer.start();
		if(iterator < anderM+5)
		{
			if (myrank==0)
			printf("    Simple Mixing ");
			explicitEuler(resGradmu, resGradM);
		}
		else
		{
			if (myrank==0)
			printf("    AndersonMixing is start ");
			
			assembAnderMatrix(assembMatrix, LestMatrix, rigTerm, nr);

			timer.reset();
			timer.start();
			for(int i = 0; i < nr; i++)
			{
				if(LestMatrix[i][i] <= 0)
				{
					if (myrank==0)
					{
						printf("Anderson iteration: matrix[%d][%d] = %f\n", i, i, LestMatrix[i][i]);
						printf("Cannot launch Anderson iteration\n");
					}
					break;
				}
			}
			getAnderCoeff(coeff, LestMatrix, rigTerm, nr);

		   for(int i = 0; i < nr; i++)
				MPI_Bcast(&coeff[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			timer.pause();
			
			if (myrank==0)
				printf("  LU solver use time   : %f seconds\n", timerGauss.get_current_time());
	
			double gamma = 1.0-std::pow(0.9, iterator);
			//mu
			for(int i = 0; i < Nspecies; i++)
			{
				for(int j = 0; j < alloc_local; j++)
				{
					mutmp[i][j][0] = fieldmu[i][j][0];
					mutmp[i][j][1] = fieldmu[i][j][1];
					
					gradmutmp[i][j][0] = gradmu[i][j][0];
					gradmutmp[i][j][1] = gradmu[i][j][1];
					
					for(int k = 0; k < nr; k++)
					{
						mutmp[i][j][0] += coeff[k]*(fieldmuold[i][k][j][0]-fieldmu[i][j][0]);
						mutmp[i][j][1] += coeff[k]*(fieldmuold[i][k][j][1]-fieldmu[i][j][1]);
						
						gradmutmp[i][j][0] += coeff[k]*(gradmuold[i][k][j][0]-dmu[i][j][0]);
						gradmutmp[i][j][1] += coeff[k]*(gradmuold[i][k][j][1]-dmu[i][j][1]);
					}
				}
			}
			//M
			for(int i = 0; i < 4; i++)
			{
				for(int j = 0; j < (alloc_local); j++)
				{
					hatMtmp[i][j][0] = hatM[i][j][0];
					hatMtmp[i][j][1] = hatM[i][j][1];
					
					gradMtmp[i][j][0] = dM[i][j][0];
					gradMtmp[i][j][1] = dM[i][j][1];
					for(int k = 0; k < nr; k++)
					{
						hatMtmp[i][j][0] += coeff[k]*(hatMold[i][k][j][0]-hatM[i][j][0]);
						hatMtmp[i][j][1] += coeff[k]*(hatMold[i][k][j][1]-hatM[i][j][1]);

						gradMtmp[i][j][0] += coeff[k]*(gradMold[i][k][j][0]-dM[i][j][0]);
						gradMtmp[i][j][1] += coeff[k]*(gradMold[i][k][j][1]-dM[i][j][1]);
					}
				}
			}
			//B
			for(int i = 0; i < DIM; i++)
			{
				Btmp[i][i]= rcpBox[i][i];
				gradBtmp[i][i] = dB[i][i];
				for(int k = 0; k < nr; k++)
					{
						Btmp[i][i] += coeff[k]*(Bold[k][i][i]-rcpBox[i][i]);
						gradBtmp[i][i] += coeff[k]*(gradBold[k][i][i]-dB[i][i]);
					}
			}
			
			gamma = 1.0- std::pow(0.9, iterator);

			for(int i = 0; i < Nspecies; i++)
				FuncsLinear2Cplx(fieldmu[i], alloc_local, 1.0,  mutmp[i], gamma, gradmutmp[i]); 
				
			for(int i = 0; i < 4; i++)
				FuncsLinear2Cplx(hatM[i], alloc_local, 1.0,  hatMtmp[i], gamma, gradMtmp[i]); 

			for(int i = 0; i < DIM; i++)
				rcpBox[i][i] = Btmp[i][i] + gamma*gradBtmp[i][i];
			
			for(int i=0; i<DIM; i++)
				dirBox[i][i] = 2*PI/rcpBox[i][i];
   
			if(myrank==0)
			{
				printf("\t\t\t === rcpBox === \n");
				MatPrint(rcpBox, DimCpt, DimCpt);
				printf("\n");
				printf("\t\t\t === dirBox === \n");
				MatPrint(dirBox, DimCpt, DimCpt);
				printf("\n");
			}
		}//anderson end 
		// updated d 
		for(int i = 0; i < Nspecies; i++)
			FuncsLinear2Cplx(dmu[i], alloc_local, 1.0, fieldmu[i], -1.0, fieldmuold[i][0]);
		
		for(int i = 0; i < 4; i++)
			FuncsLinear2Cplx(dM[i], alloc_local, 1.0, hatM[i], -1.0, hatMold[i][0]);
		
		for(int i = 0; i < DIM; i++)
			dB[i][i] = rcpBox[i][i]-Bold[0][i][i];

		///存储所有历史偏差
		for(int i = 0; i < Nspecies; i++)
		{
			for(int j = anderM; j > 0; j--)
			{
				memcpy(gradmuold[i][j], gradmuold[i][j-1], sizeof(fftw_complex)*(alloc_local));
			}
			memcpy(gradmuold[i][0], dmu[i], sizeof(fftw_complex)*(alloc_local));
		}
	
		for(int i = 0; i < 4; i++)
		{
			for(int j = anderM; j > 0; j--)
			{
				memcpy(gradMold[i][j], gradMold[i][j-1], sizeof(fftw_complex)*(alloc_local));
			}
		memcpy(gradMold[i][0], dM[i], sizeof(fftw_complex)*(alloc_local));
		}
		
		for(int j = anderM; j > 0; j--)
		{
			MatCopy(gradBold[j], gradBold[j-1], DIM, DIM);
		}
		MatCopy(gradBold[0], dB, DIM, DIM);

		double tmpd=0;
		double tmpw=0;
		for(int j = 0; j < Nspecies; j++)
		{
			tmpd += Intergral_space(dmu[j], dmu[j]);
			tmpw += Intergral_space(fieldmu[j], fieldmu[j]);
		}
		for(int j = 0; j < 4; j++)
		{
			tmpd += Intergral_space(dM[j], dM[j]);
			tmpw += Intergral_space(hatM[j], hatM[j]);
		}
		
		for(int j = 0; j < DIM; j++)
		{
			tmpd+= dB[j][j]*dB[j][j];
			tmpw+= rcpBox[j][j]*rcpBox[j][j];
		}
		
		timer.pause();
		if (myrank==0)
		printf("\t updateField use time  : %f seconds\n", timer.get_current_time());

		error = std::sqrt(tmpd/tmpw);

		timer.reset();
		timer.start();
		updateMatrix(assembMatrix, gradmuold, gradMold, gradBold, nr+1);
		timer.pause();

		if (myrank==0)
		printf("\t\t time cost of assembling matrix: %f seconds\n", timer.get_current_time());

	// save data
		for(int i=0;i< Nspecies; i++)
		{
			FftwC2R(rhoReallocal[i], rho[i]);
		}
		
		for(int i=0;i< 4; i++)
			FftwC2R(M[i], hatM[i]);
		
/*
		for(int i=0;i< Nspecies; i++)
		{
			MPI_Gatherv(rhoReallocal[i], a, MPI_DOUBLE, RhoReal[i], n, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		for(int i=0;i< 4; i++)
			MPI_Gatherv(M[i], a, MPI_DOUBLE, MReal[i], n, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	     if (myrank == 0)
		  {
		   if(iterator%10== 0||iterator ==1)
			{
				write_rho(RhoReal, MReal, iterator);
				char fname[255];
				sprintf(fname, "./results-%d/rst-%f-%d", phase, fB, iterator);
				printf("savedat=%s",fname);
				FILE *fp=fopen(fname,"a");
				fprintf(fp, "%.15e, %.15e\n", FreeEnergy, singQ[0]);
				fclose(fp);
		   }
		   }	*/     
		       }while(hamDiff > tol && iterator < 5000);

	free(resGradM);
	free(resGradmu);
	
	for(int i = 0; i < Nspecies; i++)
	{
		for(int j = 0; j < anderM+1; j++)
		{
			fftw_free(fieldmuold[i][j]);
			fftw_free(gradmuold[i][j]);
		}
		fftw_free(fieldmuold[i]);
		fftw_free(gradmuold[i]);
	}
	fftw_free(fieldmuold);
	fftw_free(gradmuold);
	for(int i = 0; i < Nspecies; i++)
	{
		fftw_free(mutmp[i]);
		fftw_free(gradmutmp[i]);
		fftw_free(dmu[i]);
	}
	fftw_free(mutmp);
	fftw_free(gradmutmp);
	fftw_free(dmu);

	
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < anderM+1; j++)
		{
			fftw_free(hatMold[i][j]);
			fftw_free(gradMold[i][j]);
		}
		fftw_free(hatMold[i]);
		fftw_free(gradMold[i]);
	}
	fftw_free(hatMold);
	fftw_free(gradMold);
	
	for(int i = 0; i < 4; i++)
	{
		fftw_free(hatMtmp[i]);
		fftw_free(gradMtmp[i]);
		fftw_free(dM[i]);
	}
	fftw_free(hatMtmp);
	fftw_free(gradMtmp);
	fftw_free(dM);

	for(int i = 0; i < DIM; i++)
	{
		for(int j = 0; j < anderM+1; j++)
		{
			free(Bold[i][j]);
			free(gradBold[i][j]);
		}
		free(Bold[i]);
		free(gradBold[i]);
	}
	free(Bold);
	free(gradBold);
	
	for(int i = 0; i < 4; i++)
	{
		free(Btmp[i]);
		free(gradBtmp[i]);
		free(dB[i]);
	}
	free(Btmp);
	free(gradBtmp);
	free(dB);


	for(int i = 0; i < anderM; i++)
	{
		free(LestMatrix[i]);
	}
	free(LestMatrix);
	free(rigTerm);
	free(coeff);
	
	for(int i = 0; i < anderM+1; i++)
	{
		free(assembMatrix[i]);
	}
	free(assembMatrix);
    for(int i = 0; i < Nspecies; i++)
	{
		free(rhoReallocal[i]);
		free(RhoReal[i]);
	}
    for(int i = 0; i < 4; i++)
		free(MReal[i]);
	free(MReal);

	return hamiltonNew;
}
\end{lstlisting}

\subsubsection{BasicFunc.cpp}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{Data.cpp}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{DisplayResults.cpp}
\begin{lstlisting}
#include <algorithm>    // std::sort
#include "Data.h"
#include "Mytimer.h"
#include "BasicFunc.h"
#include "FftwToolkit.h"
#include "DisplayResults.h"

void dispDensity(fftw_complex *Csrc, std::string &densityName, int step)
{
   // printf("\n=== Output plotted data, please waiting ...");
	mytimer_t timer;
    std::string name;

	timer.reset();
	timer.start();

	//double chiN = Ndeg * chiAB;
	FILE *densityfile;
        name = densityName + std::string(".dat");
	densityfile = fopen(name.c_str(), "w");
	
	int *Nexpand = (int *)malloc(sizeof(int)*DimCpt);
        for(int i = 0; i < DimCpt; i++) {Nexpand[0] = 41;Nexpand[1] = 41;}
	int cplxDofsExpand = 1;
	int realDofsExpand = 1;
	for(int i = 0; i < DimCpt-1; i++)
	{
		cplxDofsExpand *= Nexpand[i];
		realDofsExpand *= Nexpand[i];
	}
	cplxDofsExpand *= (Nexpand[DimCpt-1]/2+1);
	realDofsExpand *= Nexpand[DimCpt-1];
//printf("cplxDofsExpand=%d\t realDofsExpand=%d\n",cplxDofsExpand,realDofsExpand);
	fftw_complex *phi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofsExpand);
	for(int i = 0; i < cplxDofsExpand; i++)
		setCplxZero(phi[i]);

	double *realDensity = (double *) malloc(sizeof(double) * realDofsExpand);
	for(int i = 0; i < realDofsExpand; i++) realDensity[i] = 0.0;
	fftw_plan p_c2r;
	p_c2r = fftw_plan_dft_c2r(DimCpt, Nexpand, phi, realDensity, FFTW_MEASURE);
////   2D images
	if(DimCpt == 2 && DimPhy == 2)
	{
//        FuncsLinear1Cplx(phi, cplxDofsExpand, 1.0, Csrc);
		expand(phi, Nexpand, Csrc, NCpt, DimCpt);
		int n0 = Nexpand[0]; int n1 = Nexpand[1];
		fftw_execute(p_c2r);
		for(int i = 0; i <n0; i++)
		{
			for(int j = 0; j <n1; j++)
			{
				fprintf(densityfile, "%.20f\t", realDensity[(i%n0)*n1+(j%n1)]);
	//			printf( "realDensity=%.20f\n", realDensity[(i%n0)*n1+(j%n1)]);
			}
			fprintf(densityfile,"\n");
		}
	}

	if(DimCpt == 3 && DimPhy == 3)
	{
		expand(phi, Nexpand, Csrc, NCpt, DimCpt);

		int n0 = Nexpand[0]; int n1 = Nexpand[1]; int n2 = Nexpand[2];

		fftw_execute(p_c2r);

		for(int i=0;i<=n0;i++)	fprintf(densityfile,"%f\t",i*dirBox[0][0]/(double)n0);
		fprintf(densityfile,"\n");
		
		for(int i=0;i<=n1;i++)	fprintf(densityfile,"%f\t",i*dirBox[1][1]/(double)n1);
		fprintf(densityfile,"\n");

		for(int i=0;i<=n2;i++)	fprintf(densityfile,"%f\t",i*dirBox[2][2]/(double)n2);
		fprintf(densityfile,"\n");
		
		for(int i=0;i<=n0;i++)
		{
			for(int j=0;j<=n1;j++)
			{
				for(int k=0;k<=n2;k++)
				{
					fprintf(densityfile,"%f ",realDensity[k%n2+(j%n1)*n2+(i%n0)*n2*n1]);
				}
				fprintf(densityfile,"\n");
			}
			fprintf(densityfile,"\n");
		}

		fprintf(densityfile,"\n");
	}

	timer.pause();
	printf("\t\t time cost of output plotted data, step %d : %f seconds\n", step, timer.get_current_time());
//    printf("\n=== Output plotted data, Step %d: %ld seconds, \t %f minutes, \t %f hours\n\n", 
//            step, timeSum, double(timeSum/60.0), double(timeSum/3600.0));

	fftw_destroy_plan(p_c2r);
	free(realDensity);
	free(Nexpand);
	fftw_free(phi);
	fclose(densityfile);
}

void dispDensityPara(fftw_complex **Csrc, std::string &densityName, int step)
{
	mytimer_t timer;
    std::string name;

	timer.reset();

	FILE *densityfile;
        name = densityName + std::string(".dat");
	densityfile = fopen(name.c_str(), "w");
 double *MU  = (double *)malloc(sizeof(double)*realDofs);
 double *MU1 = (double *)malloc(sizeof(double)*realDofs);
 double *MU2 = (double *)malloc(sizeof(double)*realDofs);
 double *MU3 = (double *)malloc(sizeof(double)*realDofs);
	
	int *Nexpand = (int *)malloc(sizeof(int)*DimCpt);
        for(int i = 0; i < DimCpt; i++) {Nexpand[0] = 41;Nexpand[1] = 41;}
	int cplxDofsExpand = 1;
	int realDofsExpand = 1;
	for(int i = 0; i < DimCpt-1; i++)
	{
		cplxDofsExpand *= Nexpand[i];
		realDofsExpand *= Nexpand[i];
	}
	cplxDofsExpand *= (Nexpand[DimCpt-1]/2+1);
	realDofsExpand *= Nexpand[DimCpt-1];
//printf("cplxDofsExpand=%d\t realDofsExpand=%d\n",cplxDofsExpand,realDofsExpand);

 FftwC2R(MU,  hatM[0]);
 FftwC2R(MU1, hatM[1]);
 FftwC2R(MU2, hatM[2]);
 FftwC2R(MU3, hatM[3]);
                        for(int ir=0; ir<realDofs;ir++)
                        {
                                fprintf(densityfile,"%.17e\t %.17e\t %.17e\t %.17e\n",MU[ir],MU1[ir],MU2[ir],MU3[ir]);
                              //  if((ir+1)%41 == 0)
                               // fprintf(fp1,"\n");
                        }
	timer.pause();
	printf("\t\t time cost of output plotted data, step %d : %f seconds\n", step, timer.get_current_time());
	fclose(densityfile);
	free(MU);
	free(MU1);
	free(MU2);
	free(MU3);
}


void dispDensitycoil(Interval *range, fftw_complex **coil, std::string &densityName, int step)
{
	mytimer_t timer;
    std::string name;

	timer.reset();
//printf("range=%d\n",range->n);
	FILE *densityfile;
        name = densityName + std::string(".dat");
	densityfile = fopen(name.c_str(), "w");
 double **realqcoil  = (double **)malloc(sizeof(double*)*(range->n+1));
 for(int is=0; is<=range->n; is++)
 { 
     realqcoil[is]  = (double *)malloc(sizeof(double)*realDofs);
 }
// double *MU1 = (double *)malloc(sizeof(double)*realDofs);
	
	int *Nexpand = (int *)malloc(sizeof(int)*DimCpt);
        for(int i = 0; i < DimCpt; i++) {Nexpand[0] = 41;Nexpand[1] = 41;}
	int cplxDofsExpand = 1;
	int realDofsExpand = 1;
	for(int i = 0; i < DimCpt-1; i++)
	{
		cplxDofsExpand *= Nexpand[i];
		realDofsExpand *= Nexpand[i];
	}
	cplxDofsExpand *= (Nexpand[DimCpt-1]/2+1);
	realDofsExpand *= Nexpand[DimCpt-1];
//printf("cplxDofsExpand=%d\t realDofsExpand=%d\n",cplxDofsExpand,realDofsExpand);
for(int is=0; is<=range->n; is++)
{
    FftwC2R(realqcoil[is],coil[is] );
}

             for(int is=0; is<=range->n;is++)
                {
                        for(int ir=0; ir<realDofs;ir++)
                        {
                                fprintf(densityfile,"%.17e\t",realqcoil[is][ir]);
                               if((ir+1)%41 == 0)
                                fprintf(densityfile,"\n");
                        }
                                fprintf(densityfile,"\n");
				}
	timer.pause();
	printf("\t\t time cost of output plotted data, step %d : %f seconds\n", step, timer.get_current_time());
	fclose(densityfile);
	for(int is=0; is<=range->n; is++)
	free(realqcoil[is]);
	free(realqcoil);
}



void dispDensityrod(Interval *range, fftw_complex ***rod, std::string &densityName, int step)
{
	mytimer_t timer;
    std::string name;

	timer.reset();
//printf("range=%d\n",range->n);
	FILE *densityfile;
        name = densityName + std::string(".dat");
	densityfile = fopen(name.c_str(), "w");
 double ***realrod  = (double ***)malloc(sizeof(double**)*Ntheta);
 for(int iu=0; iu<Ntheta; iu++)
 {
            realrod[iu]  = (double **)malloc(sizeof(double*)*(range->n+1));
      for(int is=0; is<=range->n; is++)
         { 
            realrod[iu][is]  = (double *)malloc(sizeof(double)*realDofs);
          }
 }
	
	int *Nexpand = (int *)malloc(sizeof(int)*DimCpt);
        for(int i = 0; i < DimCpt; i++) {Nexpand[0] = 41;Nexpand[1] = 41;}
	int cplxDofsExpand = 1;
	int realDofsExpand = 1;
	for(int i = 0; i < DimCpt-1; i++)
	{
		cplxDofsExpand *= Nexpand[i];
		realDofsExpand *= Nexpand[i];
	}
	cplxDofsExpand *= (Nexpand[DimCpt-1]/2+1);
	realDofsExpand *= Nexpand[DimCpt-1];
//printf("cplxDofsExpand=%d\t realDofsExpand=%d\n",cplxDofsExpand,realDofsExpand);
  for(int iu=0; iu<Ntheta; iu++)
  {
	for(int is=0; is<=range->n; is++)
       {
         FftwC2R(realrod[iu][is],rod[iu][is] );
       }
  }
for(int iu=0;iu<Ntheta;iu++)
{
             for(int is=0; is<=range->n;is++)
                {
                        for(int ir=0; ir<realDofs;ir++)
                        {
                                fprintf(densityfile,"%.17e\t",realrod[iu][is][ir]);
                               if((ir+1)%41 == 0)
                                fprintf(densityfile,"\n");
                        }
                                fprintf(densityfile,"\n");
				}
                                fprintf(densityfile,"\n");
}
	timer.pause();
	printf("\t\t time cost of output plotted data, step %d : %f seconds\n", step, timer.get_current_time());
	fclose(densityfile);
	for(int iu=0;iu<Ntheta;iu++)
{
    	for(int is=0; is<=range->n; is++)
		{
        	free(realrod[iu][is]);
		}
        	free(realrod[iu]);
}
	free(realrod);
}
void dispDensity_Lx(fftw_complex *Csrc, std::string &densityName, int step)
{
//    printf("\n=== Output plotted data, please waiting ...");
        mytimer_t timer;
    std::string name;

        timer.reset();
        timer.start();

        double chiN = Ndeg * chiAB;
        FILE *densityfile;
        name = densityName + std::string(".dat");
        densityfile = fopen(name.c_str(), "w");

        int *Nexpand = (int *)malloc(sizeof(int)*DimCpt);
        for(int i = 0; i < DimCpt; i++) {Nexpand[0] = 64;Nexpand[1] = 64;}
        int cplxDofsExpand = 1;
        int realDofsExpand = 1;
        for(int i = 0; i < DimCpt-1; i++)
        {
                cplxDofsExpand *= Nexpand[i];
                realDofsExpand *= Nexpand[i];
        }
        cplxDofsExpand *= (Nexpand[DimCpt-1]/2+1);
        realDofsExpand *= Nexpand[DimCpt-1];

        fftw_complex *phi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofsExpand);
        for(int i = 0; i < cplxDofsExpand; i++)
                setCplxZero(phi[i]);

        double *realDensity = (double *) malloc(sizeof(double) * realDofsExpand);
        for(int i = 0; i < realDofsExpand; i++) realDensity[i] = 0.0;
        fftw_plan p_c2r;
        p_c2r = fftw_plan_dft_c2r(DimCpt, Nexpand, phi, realDensity, FFTW_MEASURE);
////   2D images
        if(DimCpt == 2 && DimPhy == 2)
        {
                expand(phi, Nexpand, Csrc, NCpt, DimCpt);
                int n0 = Nexpand[0]; int n1 = Nexpand[1];
                fftw_execute(p_c2r);

                for(int i=0;i<=n0;i++)	fprintf(densityfile,"%f\t",i*dirBox[0][0]/(double)n0);
                fprintf(densityfile,"\n");
        }

        if(DimCpt == 3 && DimPhy == 3)
        {
                expand(phi, Nexpand, Csrc, NCpt, DimCpt);

                int n0 = Nexpand[0]; int n1 = Nexpand[1]; int n2 = Nexpand[2];

                fftw_execute(p_c2r);

                for(int i=0;i<=n0;i++)	fprintf(densityfile,"%f\t",i*dirBox[0][0]/(double)n0);
                fprintf(densityfile,"\n");

                for(int i=0;i<=n1;i++)	fprintf(densityfile,"%f\t",i*dirBox[1][1]/(double)n1);
                fprintf(densityfile,"\n");

                for(int i=0;i<=n2;i++)	fprintf(densityfile,"%f\t",i*dirBox[2][2]/(double)n2);
                fprintf(densityfile,"\n");

                for(int i=0;i<=n0;i++)
                {
                        for(int j=0;j<=n1;j++)
                        {
                                for(int k=0;k<=n2;k++)
                                {
                                        fprintf(densityfile,"%f ",realDensity[k%n2+(j%n1)*n2+(i%n0)*n2*n1]);
                                }
                                fprintf(densityfile,"\n");
                        }
                        fprintf(densityfile,"\n");
                }

                fprintf(densityfile,"\n");
        }

        timer.pause();
        printf("\t\t time cost of output plotted data, step %d : %f seconds\n", step, timer.get_current_time());
//    printf("\n=== Output plotted data, Step %d: %ld seconds, \t %f minutes, \t %f hours\n\n",
//            step, timeSum, double(timeSum/60.0), double(timeSum/3600.0));

        fftw_destroy_plan(p_c2r);
        free(realDensity);
        free(Nexpand);
        fftw_free(phi);
        fclose(densityfile);
}
void dispDensity_Ly(fftw_complex *Csrc, std::string &densityName, int step)
{
//    printf("\n=== Output plotted data, please waiting ...");
        mytimer_t timer;
    std::string name;

        timer.reset();
        timer.start();

        double chiN = Ndeg * chiAB;
        FILE *densityfile;
        name = densityName + std::string(".dat");
        densityfile = fopen(name.c_str(), "w");

        int *Nexpand = (int *)malloc(sizeof(int)*DimCpt);
        for(int i = 0; i < DimCpt; i++) {Nexpand[0] = 64;Nexpand[1] = 64;}
        int cplxDofsExpand = 1;
        int realDofsExpand = 1;
        for(int i = 0; i < DimCpt-1; i++)
        {
                cplxDofsExpand *= Nexpand[i];
                realDofsExpand *= Nexpand[i];
        }
        cplxDofsExpand *= (Nexpand[DimCpt-1]/2+1);
        realDofsExpand *= Nexpand[DimCpt-1];

        fftw_complex *phi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*cplxDofsExpand);
        for(int i = 0; i < cplxDofsExpand; i++)
                setCplxZero(phi[i]);

        double *realDensity = (double *) malloc(sizeof(double) * realDofsExpand);
        for(int i = 0; i < realDofsExpand; i++) realDensity[i] = 0.0;
        fftw_plan p_c2r;
        p_c2r = fftw_plan_dft_c2r(DimCpt, Nexpand, phi, realDensity, FFTW_MEASURE);
////   2D images
        if(DimCpt == 2 && DimPhy == 2)
        {
//        FuncsLinear1Cplx(phi, cplxDofsExpand, 1.0, Csrc);
                expand(phi, Nexpand, Csrc, NCpt, DimCpt);
                int n0 = Nexpand[0]; int n1 = Nexpand[1];
                fftw_execute(p_c2r);
                for(int i=0;i<=n1;i++)	fprintf(densityfile,"%f\t",i*dirBox[1][1]/(double)n1);
                fprintf(densityfile,"\n");
        }
////   Actually for 3D images
//    if(DimCpt == 2 && DimPhy == 2)
//    {
//        expand(phi, Nexpand, Csrc, NCpt, DimCpt);
//        int n0 = Nexpand[0]; int n1 = Nexpand[1]; int n2 = Nexpand[0];
//        int n2 = 1;
//        fftw_execute(p_c2r);
//
//        for(int i=0;i<=n0;i++)	fprintf(densityfile,"%f\t",i*dirBox[0][0]/(double)n0);
//        fprintf(densityfile,"\n");
//
//        for(int i=0;i<=n1;i++)	fprintf(densityfile,"%f\t",i*dirBox[1][1]/(double)n1);
//        fprintf(densityfile,"\n");
//
//        for(int i=0;i<=n2;i++)	fprintf(densityfile,"%f\t",i*dirBox[0][0]/(double)n2);
//        fprintf(densityfile,"\n");
//
//        for(int k = 0; k <= n2; k++)
//        {
//            for(int i = 0; i <= n0; i++)
//            {
//                for(int j = 0; j <= n1; j++)
//                {
//                    fprintf(densityfile, "%f\t", realDensity[(i%n0)*n1+(j%n1)]);
//                }
//                fprintf(densityfile,"\n");
//            }
//            fprintf(densityfile,"\n");
//        }
//        fprintf(densityfile,"\n");
//    }

        if(DimCpt == 3 && DimPhy == 3)
        {
                expand(phi, Nexpand, Csrc, NCpt, DimCpt);

                int n0 = Nexpand[0]; int n1 = Nexpand[1]; int n2 = Nexpand[2];

                fftw_execute(p_c2r);

                for(int i=0;i<=n0;i++)	fprintf(densityfile,"%f\t",i*dirBox[0][0]/(double)n0);
                fprintf(densityfile,"\n");

                for(int i=0;i<=n1;i++)	fprintf(densityfile,"%f\t",i*dirBox[1][1]/(double)n1);
                fprintf(densityfile,"\n");

                for(int i=0;i<=n2;i++)	fprintf(densityfile,"%f\t",i*dirBox[2][2]/(double)n2);
                fprintf(densityfile,"\n");

                for(int i=0;i<=n0;i++)
                {
                        for(int j=0;j<=n1;j++)
                        {
                                for(int k=0;k<=n2;k++)
                                {
                                        fprintf(densityfile,"%f ",realDensity[k%n2+(j%n1)*n2+(i%n0)*n2*n1]);
                                }
                                fprintf(densityfile,"\n");
                        }
                        fprintf(densityfile,"\n");
                }

                fprintf(densityfile,"\n");
        }

        timer.pause();
        printf("\t\t time cost of output plotted data, step %d : %f seconds\n", step, timer.get_current_time());
//    printf("\n=== Output plotted data, Step %d: %ld seconds, \t %f minutes, \t %f hours\n\n",
//            step, timeSum, double(timeSum/60.0), double(timeSum/3600.0));

        fftw_destroy_plan(p_c2r);
        free(realDensity);
        free(Nexpand);
        fftw_free(phi);
        fclose(densityfile);
}

void dispPlaneWave(fftw_complex *Csrc, int n, std::string & vecName)
{
    std::string name;
	FILE *fprojPlane;
    name = vecName + std::string(".txt");
	fprojPlane = fopen(name.c_str(), "w");
	
	double *planeTmp = (double *)malloc(sizeof(double)*DimPhy);
	int index;
	for(int i = 0; i < n; i ++)
	{
		double tmp = std::sqrt(Csrc[i][0]*Csrc[i][0]+Csrc[i][1]*Csrc[i][1]);
		if(tmp > 1.0e-04) index = 1;
		else index = 0;
		
		for(int kk = 0; kk < DimPhy; kk ++)
		{
			planeTmp[kk] = 0.0;
		}

		for(int j = 0; j < DimPhy; j++)
		{
			double mnt = 0.0;
			for(int ii = 0; ii < DimPhy; ii++)
			{
				mnt += indKspace[i][ii]*rcpBox[ii][j]; 
			}
			for(int jj = 0; jj < DimPhy; jj++)
			{
				planeTmp[jj] += ProjMatrix[jj][j]*mnt;
			}
		}
		for(int j = 0; j < DimPhy; j++)
		{
			fprintf(fprojPlane, "%f\t ", planeTmp[j]);
		}
//            fprintf(fprojPlane, "%d\t ", index);
		fprintf(fprojPlane, "%e\t ", tmp);
		fprintf(fprojPlane, "\n");

		if(indKspace[i][DimPhy-1]>0)
		{
			for(int kk = 0; kk < DimPhy; kk ++)
			{
				planeTmp[kk] = 0.0;
			}
			for(int j = 0; j < DimPhy; j++)
			{
				double mnt = 0.0;
				for(int ii = 0; ii < DimPhy; ii++)
				{
					mnt += -1.0*indKspace[i][ii]*rcpBox[ii][j]; 
				}
				for(int jj = 0; jj < DimPhy; jj++)
				{
					planeTmp[jj] += ProjMatrix[jj][j]*mnt;
				}
			}
			for(int j = 0; j < DimPhy; j++)
			{
				fprintf(fprojPlane, "%f\t ", planeTmp[j]);
			}
//                fprintf(fprojPlane, "%d\t ", index);
			fprintf(fprojPlane, "%e\t ", tmp);
			fprintf(fprojPlane, "\n");
		}
	}
	free(planeTmp);
	fclose(fprojPlane);
}

bool myComp(mySortVec a, mySortVec b)
{
	return (a.Data[0]*a.Data[0]+a.Data[1]*a.Data[1] > b.Data[0]*b.Data[0]+b.Data[1]*b.Data[1]);
}

void dispCplxDensity(fftw_complex *Csrc, int **kspace, int n, int dim, std::string & name)
{

    name += std::string(".txt");
	FILE* cplxDensity = fopen(name.c_str(), "w");
	fprintf(cplxDensity, "%d\n", n);

	mySortVec *myVector = (mySortVec*)malloc(sizeof(mySortVec)*n);
	for(int k = 0; k < n; k++)
	{
		myVector[k].Data  = (double *)malloc(sizeof(double)*2);
		myVector[k].Index = (int *)malloc(sizeof(int)*DimCpt);
		for(int i = 0; i < dim; i++)
			myVector[k].Index[i] = kspace[k][i];
		for(int j = 0; j < 2; j++)
			myVector[k].Data[j] = Csrc[k][j];
	}

////      sort 
	std::sort(myVector, myVector+n, myComp);
	for(int k = 0; k < n; k++)
	{
		for(int i = 0; i < dim; i++)
		{
			fprintf(cplxDensity, "%d\t", myVector[k].Index[i]);
		}
		fprintf(cplxDensity, "%e\t%e\n", myVector[k].Data[0], myVector[k].Data[1]);
	}
////      unsort 
//    double TOL = 1.0e-14;
//    for(int k = 0; k < n; k++)
//    {
//        double tmp = std::sqrt(Csrc[k][0]*Csrc[k][0]+Csrc[k][1]*Csrc[k][1]);
//        if(tmp>TOL)
//        {
//            for(int i = 0; i < dim; i++)
//            {
//                fprintf(cplxDensity, "%d\t", kspace[k][i]);
//            }
//            fprintf(cplxDensity, "%.15f\t%.15f\n", Csrc[k][0], Csrc[k][1]);
//        }
//    }

	fclose(cplxDensity);
	for(int i = 0; i < n; i++)
	{
		free(myVector[i].Index);
		free(myVector[i].Data);
	}
	free(myVector);
}

void dispRealDensity(double *Rsrc, int **kspace, int n, int dim, std::string & name)
{
    name += std::string(".txt");
	FILE* fdensity = fopen(name.c_str(), "w");

	for(int k = 0; k < n; k++)
	{
		for(int i = 0; i < dim; i++)
		{
			fprintf(fdensity, "%d\t", kspace[k][i]);
		}
		fprintf(fdensity, "%e\n", Rsrc[k]);
	}
	fclose(fdensity);
}
\end{lstlisting}

\subsubsection{FftwToolkit.cpp}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{GaussElimination.cpp}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{Initialization.cpp}
\begin{lstlisting}
#include "mpi.h"
#include "fftw3-mpi.h" 

#include "Head.h"
#include "Data.h"
#include "Mytimer.h"
#include "BasicFunc.h"
#include "FftwToolkit.h"
#include "SCFTBaseAB.h"
#include "Initialization.h"

void initialize()
{
	mytimer_t timer;

	// allocate memory 
	timer.reset();
	timer.start();
	memAllocation_no_time();
	timer.pause();
	if (myrank==0)
		printf("\t\t time cost of memory allocation : %f seconds\n", timer.get_current_time());
	timer.reset();
	timer.start();
	getIndex(indKspace, cplxDofs, DimCpt, NCpt);
	timer.pause();
	if (myrank==0)
	printf("\t\t time cost of getIndex : %f seconds\n\n", timer.get_current_time());
}

void memAllocation_no_time()
{
	singQ = (double *)malloc(sizeof(double)*Nblend);
	for(int i = 0; i < Nblend; i++) 
    {
		singQ[i] = 0.0;
	}

	// Initialize computational box
	dirBox = (double **)malloc(sizeof(double*)*DimCpt);
	rcpBox = (double **)malloc(sizeof(double*)*DimCpt);
	for(int i = 0; i < DimCpt; i++)
	{
		dirBox[i] = (double *)malloc(sizeof(double)*DimCpt);
		rcpBox[i] = (double *)malloc(sizeof(double)*DimCpt);
	}
	for(int i = 0; i < DimCpt; i++)
	{
		for(int j = 0; j < DimCpt; j++)
		{
			dirBox[i][j] = 0.0; 	rcpBox[i][j] = 0.0;
		}
	}

    gradB = (double **)malloc(sizeof(double*)*DimCpt);
    for(int i = 0; i < DimCpt; i++)
    {
        gradB[i] = (double *)malloc(sizeof(double)*DimCpt);
    }
    for(int i = 0; i < DimCpt; i++)
    {
        for(int j = 0; j < DimCpt; j++)
        {
            gradB[i][j] = 0.0;
        }
    }


	// projective matrix
	ProjMatrix = (double **)malloc(sizeof(double*)*DimPhy);
	for(int i = 0; i < DimPhy; i++)
		ProjMatrix[i] = (double *)malloc(sizeof(double)*DimCpt);

	for(int i = 0; i < DimPhy; i++)
		for(int j = 0; j < DimCpt; j++)
			ProjMatrix[i][j] = 0.0;

	for(int i = 0; i < DimPhy; i++)
			ProjMatrix[i][i] = 1.0;
	
    indKspace = (int **)malloc(sizeof(int*)*cplxDofs);
	for(int i = 0; i < cplxDofs; i++)
		indKspace[i] = (int *)malloc(sizeof(int)*DimCpt);
	for(int i = 0; i < cplxDofs; i++)
		for(int j = 0; j < DimCpt; j++)
			indKspace[i][j] = 0;
	
    indKspacelocal = (int **)malloc(sizeof(int*)*alloc_local);
	for(int i = 0; i < alloc_local; i++)
		indKspacelocal[i] = (int *)malloc(sizeof(int)*DimCpt);
	
	for(int i = 0; i < alloc_local; i++)
		for(int j = 0; j < DimCpt; j++)
			indKspacelocal[i][j] = 0;
	
	// |K|^2, |K|
	Gsquare = (double *)malloc(sizeof(double) *cplxDofs);
	Gsquarelocal = (double *)malloc(sizeof(double) *alloc_local);

	for(int i = 0; i < alloc_local; i++)
		Gsquarelocal[i]=0;

	for(int i = 0; i < cplxDofs; i++)
		Gsquare[i]=0;
	// density and field functions
    rho = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*Nspecies);
    rhoglobal = (double **)malloc(sizeof(double*)*Nspecies);
    rhoreal = (double **)malloc(sizeof(double*)*Nspecies);

    gradW = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*Nspecies);
    fieldW = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*Nspecies);
	
	fieldWplus  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
	
	for(int i = 0; i < Nspecies; i++)
	{
	    fieldW[i]  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
	    gradW[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
        rho[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
        rhoglobal[i] = (double *)malloc(sizeof(double)*(2*cplxDofs));
        rhoreal[i] = (double *)malloc(sizeof(double)*(2*alloc_local));
    }
	for(int i = 0; i < Nspecies; i++)
	{
		for(int j = 0; j < alloc_local; j++)
		{
			setCplxZero(fieldW[i][j]);
			setCplxZero(gradW[i][j]);
         	setCplxZero(rho[i][j]);
        }
    }

	for(int j = 0; j < alloc_local; j++)
	{
		fieldWplus[j][0]=0.0;
		fieldWplus[j][1]=0.0;
	}

// temp arrays for FFTW
	fftw_Ctmp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
	fftw_Rtmp = (double *)malloc(sizeof(double) *(2*alloc_local));
// fftw_plan
    PlanR2C = fftw_mpi_plan_dft_r2c(DIM, NCpt, fftw_Rtmp, fftw_Ctmp, MPI_COMM_WORLD, FFTW_ESTIMATE);
    PlanC2R = fftw_mpi_plan_dft_c2r(DIM, NCpt, fftw_Ctmp, fftw_Rtmp, MPI_COMM_WORLD, FFTW_ESTIMATE);	
}

void memAllocation()
{
	Nsplit = (int *)malloc(sizeof(int)*Nblock);
	for(int i = 0; i < Nblock; i++) Nsplit[i] = 0;
    Nsplit[0] = round((fA+1.0e-8) / dsMax);
    Nsplit[1] = round((fB+1.0e-8) / dsMax);

	Ranges = (Interval *)malloc(sizeof(Interval)*Nblock);
	Ranges[0].a = 0.0;
	Ranges[0].b = fA;

	Ranges[1].a = Ranges[0].b;
	Ranges[1].b = Ranges[0].b+fB;

	for(int i = 0; i < Nblock; i++)
		Ranges[i].n = Nsplit[i];

	Ns=0;
	for(int i = 0; i < Nblock; i++) 
		Ns += Nsplit[i];
	Ns = Ns+1;

	// forward and backward propagators
	frdQ = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*Ns);
	bakQ = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*Ns);
	hatq_qplus = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*Ns);
	for(int i = 0; i < Ns; i++)
	{
		frdQ[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
		bakQ[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
		hatq_qplus[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
	}
	for(int i=0; i < Ns; i++)
	{
		for(int j=0; j < alloc_local; j++)
		{
			setCplxZero(frdQ[i][j]);
			setCplxZero(bakQ[i][j]);
			setCplxZero(hatq_qplus[i][j]);
		}
	}
}


void getGsquare()
{
	for(int i = 0; i < cplxDofs; i++)
		Gsquare[i]=0;
	double *tmp = (double *)malloc(sizeof(double)*DimPhy);
	for(int i = 0; i < DimPhy; i++)
		tmp[i] = 0.0;
	for(int i = 0; i < cplxDofs; i++)
	{
		tmp[0] = 0.0; tmp[1]=0.0; tmp[2]=0.0;
		for(int k = 0; k < DimPhy; k++) tmp[k] = 0.0;
		for(int ii = 0; ii < DimCpt; ii++)
		{
			double mnt = 0.0;
			for(int jj = 0; jj < DimCpt; jj++)
			{
				mnt += indKspace[i][jj]*rcpBox[jj][ii]; 
			}
			for(int kk = 0; kk < DimPhy; kk++) 
			{
				tmp[kk] += ProjMatrix[kk][ii]*mnt;
			}
		}
		for(int kk = 0; kk < DimPhy; kk++)
			Gsquare[i] += pow(tmp[kk], 2);
	}

   for (int i = 0; i < local_n0; i++)
   {
	   for (int j = 0; j < NCpt[1]; j++)
	   {
			for (int k = 0; k < NCpt[2]/2+1; k++)
			{
				int index = i*NCpt[1]*(NCpt[2]/2+1)+j*(NCpt[2]/2+1)+k;
				int index1 = (i+local_0_start)*NCpt[1]*(NCpt[2]/2+1)+j*(NCpt[2]/2+1)+k;
				Gsquarelocal[index] = Gsquare[index1];
			}
		}
    }         
	free(tmp);
}


void getRecipLattice(double **dBox, double **rBox)
{
	if(DimCpt==2)
	{
		double det = dBox[0][0]*dBox[1][1]-dBox[0][1]*dBox[1][0];
		rBox[0][0] =  2*PI*dBox[1][1] / det;
		rBox[1][0] = -2*PI*dBox[1][0] / det;
		rBox[0][1] = -2*PI*dBox[0][1] / det;
		rBox[1][1] =  2*PI*dBox[0][0] / det;
	}
	if(DimCpt == 3)
	{
		double aone[3], atwo[3], athree[3];
		double bone[3], btwo[3], bthree[3];
		double cone[3], ctwo[3], cthree[3];
		double volume[3];

		for(int i = 0; i < 3; i++)
		{
			aone[i]   = dBox[0][i];
			atwo[i]   = dBox[1][i];
			athree[i] = dBox[2][i];
		}

		cone[0]=atwo[1]*athree[2]-atwo[2]*athree[1];
		cone[1]=atwo[2]*athree[0]-atwo[0]*athree[2];
		cone[2]=atwo[0]*athree[1]-atwo[1]*athree[0];
		volume[0]=aone[0]*cone[0]+aone[1]*cone[1]+aone[2]*cone[2];
		bone[0]=2*PI/volume[0]*cone[0];
		bone[1]=2*PI/volume[0]*cone[1];
		bone[2]=2*PI/volume[0]*cone[2];    
	    
		ctwo[0]=athree[1]*aone[2]-athree[2]*aone[1];
		ctwo[1]=athree[2]*aone[0]-athree[0]*aone[2];
		ctwo[2]=athree[0]*aone[1]-athree[1]*aone[0];
		volume[1]=atwo[0]*ctwo[0]+atwo[1]*ctwo[1]+atwo[2]*ctwo[2];
		btwo[0]=2*PI/volume[1]*ctwo[0];
		btwo[1]=2*PI/volume[1]*ctwo[1];
		btwo[2]=2*PI/volume[1]*ctwo[2];  
	        
		cthree[0]=aone[1]*atwo[2]-aone[2]*atwo[1];
		cthree[1]=aone[2]*atwo[0]-aone[0]*atwo[2];
		cthree[2]=aone[0]*atwo[1]-aone[1]*atwo[0];
		volume[2]=athree[0]*cthree[0]+athree[1]*cthree[1]+athree[2]*cthree[2];
		bthree[0]=2*PI/volume[2]*cthree[0];
		bthree[1]=2*PI/volume[2]*cthree[1];
		bthree[2]=2*PI/volume[2]*cthree[2];   

		for(int i = 0; i < 3; i++)
		{
			rBox[0][i] = bone[i];
			rBox[1][i] = btwo[i]; 
			rBox[2][i] = bthree[i];
		}
	}
}

void writeRealData(double *field, const char *fname)
{
    int ftmp;
    FILE *fp; // = fopen(fname, "r");
    if((fp = fopen(fname, "r")) == NULL)
    {   
         printf("Cannot open file.\n");
         exit(1);
    }

    int initDofx=NCpt[0];
    int initDofy=NCpt[1];
    int initDofz=NCpt[2];

    for(int i = 0; i < initDofx; i++)
    {   
        for(int j = 0; j < initDofy; j++)
        {
			for(int k = 0; k < initDofz; k++)
			{
				ftmp = fscanf(fp,"%lf", &(field[i*(initDofy)*(initDofz)+j*(initDofz)+k]));
			}
        }
    }	
    fclose(fp);
}
            
void initFieldFourier(fftw_complex *field, const char *fname)
{
    int ftmp;
    int initDof;
    FILE *fp; // = fopen(fname, "r");

    if((fp = fopen(fname, "r")) == NULL)
    {   
         printf("Cannot open file.\n");
         exit(1);
    }   

    ftmp = fscanf(fp, "%d", &initDof);

    int **fin = (int **)malloc(sizeof(int*)*initDof);
    for(int i = 0; i < initDof; i++)
        fin[i] = (int *)malloc(sizeof(int)*DimCpt);

    fftw_complex *fieldInput = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*initDof);
    for(int i = 0; i < initDof; i++)
    {   
        for(int j = 0; j < DimCpt; j++)
        {   
            ftmp = fscanf(fp, "%d", &(fin[i][j]));
        }   
        ftmp = fscanf(fp, "%lf", &(fieldInput[i][0]));
        ftmp = fscanf(fp, "%lf", &(fieldInput[i][1]));
    }   
    for(int i = 0; i < cplxDofs; i ++)
    {
        for(int j = 0; j < initDof; j++)
        {
            if(DimCpt == 2)
            {
                if(fin[j][0]==indKspace[i][0] && fin[j][1]==indKspace[i][1])
                {
                    field[i][0] = fieldInput[j][0];
                    field[i][1] = fieldInput[j][1];
                }
            }
            
            if(DimCpt == 3)
            {
                if(fin[j][0]==indKspace[i][0] && fin[j][1]==indKspace[i][1] && fin[j][2]==indKspace[i][2])
                {
                    field[i][0] = fieldInput[j][0];
                    field[i][1] = fieldInput[j][1];
                }
            }
        }
    }

    fclose(fp);

    for(int i = 0; i < initDof; i++) free(fin[i]);
    free(fin);
    fftw_free(fieldInput);

}
\end{lstlisting}

\subsubsection{IntegrationToolkit.cpp}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{MDESolver.cpp}
\begin{lstlisting}
#include "mpi.h"                                                                
#include "fftw3-mpi.h" 
#include "Head.h"
#include "Data.h"
#include "BasicFunc.h"
#include "MDESolver.h"
#include "Initialization.h"
#include "BasicFunc.h"
#include "MDESolver.h"
#include "FftwToolkit.h"
#include "IntegrationToolkit.h"

#include "config.h"


void MDESolver2Order(Interval *range,
                     double  *WReal,
                     double **Q)
//  second-order operator-splitting scheme
{
    double ds = ((*range).b - range->a) / range->n;
    fftw_complex *Q_Ctmp;
    Q_Ctmp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
    double *Qreal  = (double *)malloc(sizeof(double) *(2*alloc_local));

	for(int i = 0; i < range->n; i++)
    {
        memcpy(Qreal, Q[i], sizeof(double)*(2*alloc_local));

        for(int j = 0; j < (2*alloc_local); j++)
            Qreal[j] *= std::exp(-WReal[j]*ds/2.0);
        FftwR2C(Q_Ctmp, Qreal);

        for(int k = 0; k < alloc_local; k ++)
        {
            Q_Ctmp[k][0] *= std::exp(-Gsquarelocal[k]*ds);
            Q_Ctmp[k][1] *= std::exp(-Gsquarelocal[k]*ds);
		}
        FftwC2R(Qreal, Q_Ctmp);
        for(int j = 0; j < 2*alloc_local; j ++)
            Qreal[j] *= std::exp(-WReal[j]*ds/2.0);
		
		memcpy(Q[i+1], Qreal, sizeof(double)*(2*alloc_local));
    }
    fftw_free(Q_Ctmp);
    free(Qreal);
}
void MDESolver4Adams(Interval *range,
                                     fftw_complex *hatW,
                                     fftw_complex **hatQ,
									 int *index)
{
        double ds = ((*range).b - range->a) / range->n;
        fftw_complex *Q_Ctmp;
        Q_Ctmp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
        double *Qreal  = (double *)malloc(sizeof(double) *(2*alloc_local));
        double *WReal= (double*)malloc(sizeof(double) *(2*alloc_local));

        fftw_complex **temp;
        temp = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*3);
        for(int i = 0; i < 3; i++)
                temp[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
        fftw_complex **tmp;
        tmp = (fftw_complex **)fftw_malloc(sizeof(fftw_complex*)*2);
        for(int i = 0; i < 2; i++)
                tmp[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
        fftw_complex *qConv = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
        fftw_complex *qRhs = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);
        fftw_complex *wqConv = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);

		FftwC2R(WReal,hatW);

       
        for(int i = 0; i < range->n; i++)
        {
                if(i < 3)
                {
                        for(int t = 0; t < 2; t++)
                        {
				 		memcpy(tmp[t], hatQ[*index+i], sizeof(fftw_complex)*alloc_local);
                                double P = std::pow(2.0, t);
                                for(int index1 = 0; index1 < P; index1++)
                                {
                                        for(int k = 0; k < alloc_local; k ++)
                                        {
                                                tmp[t][k][0] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                                tmp[t][k][1] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                        }
                                        FftwC2R(Qreal, tmp[t]);

                                        for(int j = 0; j < 2*alloc_local; j ++)
                                                Qreal[j] *= std::exp(-WReal[j]*ds/P);
                                        FftwR2C(tmp[t], Qreal);
                                        for(int k = 0; k < alloc_local; k ++)
                                        {
                                                tmp[t][k][0] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                                tmp[t][k][1] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                        }
                                }

                        }
                        FuncsLinear2Cplx(hatQ[*index+i+1],  alloc_local, -1.0/3.0, tmp[0], 4.0/3.0, tmp[1]);
		}
		else
                {
                        FuncsLinear4Cplx(qConv, alloc_local, 4.0, hatQ[*index+i], -6.0, hatQ[*index+i-1], 4.0, hatQ[*index+i-2], -1.0, hatQ[*index+i-3]);

                        FuncsLinear4Cplx(qRhs,  alloc_local, 4.0, hatQ[*index+i], -3.0, hatQ[*index+i-1], 4.0/3.0, hatQ[*index+i-2], -1.0/4.0, hatQ[*index+i-3]);
                        hatConv(wqConv, hatW, qConv);
                        for(int k = 0; k < alloc_local; k++)
                        {
                                hatQ[*index+i+1][k][0] = (qRhs[k][0] - ds*wqConv[k][0]) / (25.0/12 + Gsquarelocal[k]*ds);
                                hatQ[*index+i+1][k][1] = (qRhs[k][1] - ds*wqConv[k][1]) / (25.0/12 + Gsquarelocal[k]*ds);
                        }
                }
        }
		*index += range->n;
       for(int i = 0; i < 3; i++) fftw_free(temp[i]);
        fftw_free(temp);
	        for(int i = 0; i < 2; i++) fftw_free(tmp[i]);
        fftw_free(tmp);
        fftw_free(Q_Ctmp);
        free(WReal);
        free(Qreal);
        fftw_free(qConv);
        fftw_free(wqConv);
        fftw_free(qRhs);
}
\end{lstlisting}

\subsubsection{MemFree.cpp}
\begin{lstlisting}
#include "Head.h"
#include "Data.h"
#include "MemFree.h"

void memReleaser()
{
	free(NCpt);
	//some pra

	free(n);
	free(displs);
	
	free(singQ);
	
	for(int i = 0; i < DimCpt; i++)
	{
		free(dirBox[i]); free(rcpBox[i]); free(gradB[i]);
	}
	free(dirBox); free(rcpBox);free(gradB);
	
	for(int i = 0; i < DimPhy; i++) 
	      free(ProjMatrix[i]);
	free(ProjMatrix);

	for(int i = 0; i < cplxDofs; i++)
		free(indKspace[i]);
	free(indKspace);

	for(int i = 0; i < alloc_local; i++)
		free(indKspacelocal[i]);
	free(indKspacelocal);
	

	free(Gsquare);
	free(Gsquarelocal);

	free(Nsplit);
	for(int is=0;is<Ns;is++)
	{
		fftw_free(frdQ[is]);
		fftw_free(bakQ[is]);
		fftw_free(hatq_qplus[is]);
	}
	fftw_free(frdQ);
	fftw_free(bakQ);
	fftw_free(hatq_qplus);
///order para
	for(int i = 0; i < Nspecies; i++)
	{

		fftw_free(fieldW[i]);
		fftw_free(gradW[i]);
		fftw_free(rho[i]);
		free(rhoglobal[i]);
		free(rhoreal[i]);
	}
	fftw_free(fieldW);
	fftw_free(fieldWplus);
	fftw_free(gradW);
	fftw_free(rho);
	free(rhoglobal);
	free(rhoreal);


	fftw_free(fftw_Ctmp);
	free(fftw_Rtmp);

	//  destroy fftw plan
	fftw_destroy_plan(PlanC2R);
	fftw_destroy_plan(PlanR2C);
	}
\end{lstlisting}

\subsubsection{SCFTBaseAB.cpp}
\begin{lstlisting}
#include "mpi.h"                                                                
#include "fftw3-mpi.h" 

#include "Head.h"
#include "time.h"
#include "Data.h"
#include "Mytimer.h"
#include "Initialization.h"
#include "BasicFunc.h"
#include "MDESolver.h"
#include "FftwToolkit.h"
#include "IntegrationToolkit.h"
#include "SCFTBaseAB.h"

#include "config.h"

void write_rho(double **RhoReal,  int iter)
{                                                                               
    char fname[255];                                                            
     sprintf(fname, "./results-%d/rho/rho.[%.4f.%.4f].[%.2f].[%d].txt", phase, fA, fB, chi, iter);
    FILE *fp=fopen(fname,"w");     
    for(int i = 0; i < realDofs; i++)
	{
				fprintf(fp,"%lf %lf \n", RhoReal[0][i], RhoReal[1][i]);
	}
    fclose(fp);                                                                 
} 

void cut(double *rholocalold, double *rholocalnew)
{
	for (int i = 0; i < local_n0; i++)
      {
         for (int j = 0; j < NCpt[1]; j++)
         {
         	for (int k = 0; k < NCpt[2]; k++)
			{
             	int index = i*NCpt[1]*(2*(NCpt[2]/2+1))+j*(2*(NCpt[2]/2+1))+k;
             	int index1 = i*NCpt[1]*NCpt[2]+j*NCpt[2]+k;
              	rholocalnew[index1] = rholocalold[index];
			}
         }
      }

}

void updatePropagator()
{
	if (myrank==0)
	{
		frdQ[0][0][0]=1.0;frdQ[0][0][1]=0.0;
	}
	int it_q =0;

	MDESolver4Adams(&Ranges[0], fieldW[0], frdQ, &it_q);// A
	MDESolver4Adams(&Ranges[1], fieldW[1], frdQ, &it_q);// A

	if (myrank==0)
	{
		bakQ[0][0][0]=1.0;	bakQ[0][0][1]=0.0;
	}

	int it_qplus =0;
	MDESolver4Adams(&Ranges[1], fieldW[1], bakQ, &it_qplus);// fB
	MDESolver4Adams(&Ranges[0], fieldW[0], bakQ, &it_qplus);// BA
}

	
double updateQ()
{
	double Q;
	if (myrank == 0)
	{
		Q = frdQ[Ns-1][0][0];
		printf("Q1 = %.15e\n", Q);                
	}
	MPI_Bcast(&Q, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	return Q;
}

void updateOrderParameter()
{
	for(int i = 0; i < Ns; i++)
	{
		int it_inv1 = Ns-1-i;
		hatConv(hatq_qplus[i], frdQ[i], bakQ[it_inv1]);
	}
	
	int index = 0;
	for(int i = 0; i < Nspecies; i++)
		integration(rho[i], &Ranges[i], hatq_qplus, &index, singQ[0]);
}

void get_gradB(double dh)
{    
     double **oldB; 
     oldB = (double **)malloc(sizeof(double*)*DimCpt);
     for(int i = 0; i < DimCpt; i++)
     {   
         oldB[i] = (double *)malloc(sizeof(double)*DimCpt);
     }
     for(int i = 0; i < DimCpt; i++)
     {
         for(int j = 0; j < DimCpt; j++)
         {   
             oldB[i][j] = 0.0;
         }
     }

    MatCopy(oldB, rcpBox, DimCpt, DimCpt);

    for(int i = 0; i < DimCpt; i++)
    {   
	rcpBox[i][i] += dh;
    }
    double FR = calPartialF();
    MatCopy(rcpBox, oldB, DimCpt, DimCpt);
    
    for(int i = 0; i < DimCpt; i++)
    {   
	rcpBox[i][i] -= dh;
    }
    double FL = calPartialF();
    MatCopy(rcpBox, oldB, DimCpt, DimCpt);
    
    if(myrank==0)
    printf("[FR, FL] = [%.15e, \t %.15e]\n", FR, FL);
    
    for(int i = 0; i < DimCpt; i++)
	gradB[i][i] = (FR-FL)/(2.0*dh);
    
    
     for(int i = 0; i < DimCpt; i++)
	     free(oldB[i]);
    free(oldB);

}


void SimpleMixing(double *resGradW, double *resGradB)
{
	 double dt1 = 0.1;
	 double dt2 = 0.1;
	 double dt3 = 1.0e-06;
	
	 //update w
	 FuncsLinear3Cplx(gradW[0], alloc_local, chi*Ndeg, rho[1], 1.0, fieldWplus, -1.0, fieldW[0]);
	 FuncsLinear3Cplx(gradW[1], alloc_local, chi*Ndeg, rho[0], 1.0, fieldWplus, -1.0, fieldW[1]);
	 
	 FuncAddToCplx(fieldW[0], alloc_local, dt1, gradW[0]);
	 FuncAddToCplx(fieldW[1], alloc_local, dt2, gradW[1]);
	
	 //update w_plus
	 FuncsLinear2Cplx(fieldWplus, alloc_local, 0.5, fieldW[0], 0.5, fieldW[1]);
	 if (myrank ==0)
	 {
		fieldWplus[0][0] -= 0.5*chi*Ndeg;
	 }
 
	 double res1 = normCplxInfty(gradW[0], alloc_local);
	double res2 = normCplxInfty(gradW[1], alloc_local);
	MPI_Reduce(&res1, &resGradW[0], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);      
	MPI_Reduce(&res2, &resGradW[1], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);      
	//update domain
	double dh = 1e-05;
	get_gradB(dh);
	MatsAdd(rcpBox, 1.0, gradB, dt3, DimCpt, DimCpt);

	if(myrank==0)
	{
		printf("\n\t\t\t === Direct Box === \n");
		getRecipLattice(rcpBox, dirBox);
		MatPrint(dirBox, DimCpt, DimCpt);
		printf("\n");
	}
	  
	double resB1 = normRealInfty(gradB[0], DimCpt);
	double resB2 = normRealInfty(gradB[1], DimCpt);
	double resB3 = normRealInfty(gradB[2], DimCpt);
	MPI_Reduce(&resB1, &resGradB[0], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);      
	MPI_Reduce(&resB2, &resGradB[1], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);  
	MPI_Reduce(&resB3, &resGradB[2], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);  
}

double updateHamilton()
{
	double tmpAB;
	
	tmpAB = Intergral_space(rho[0], rho[1]);

	double tmpA, tmpB, tmpWplus;
	tmpA = Intergral_space(fieldW[0], rho[0]);
	tmpB = Intergral_space(fieldW[1], rho[1]);

	fftw_complex *rhotmp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*alloc_local);	
	FuncsLinear2Cplx(rhotmp, alloc_local, 1.0, rho[0], 1.0, rho[1]);
	tmpWplus = Intergral_space(fieldWplus, rhotmp);

	if(myrank==0)
	{
		internalEnergy = chi*Ndeg*tmpAB - tmpA -tmpB;
	//	internalEnergy = chi*Ndeg*tmpAB - tmpA -tmpB + tmpWplus - fieldWplus[0][0];
	    entropicEnergy1= - std::log(singQ[0]);
		Hamilton = internalEnergy + entropicEnergy1;
	printf("H = %.15e\t %15e\t \n", internalEnergy, chi*Ndeg*tmpAB);                
	}
	 
	MPI_Bcast(&Hamilton, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	fftw_free(rhotmp);
    return Hamilton;
}
double evlSaddle(double tol)
{
	mytimer_t timer;
	int iterator = 0;
	double error;
	double Energy, oldEnergy, diffEnergy;
	oldEnergy = 100.0;
	Energy = 100.0;
	double *resGradW  = (double *)malloc(sizeof(double)*Nspecies);
	double *resGradB = (double *)malloc(sizeof(double)*DimCpt);
	for(int i = 0; i < Nspecies; i++) 
	{
		resGradW[i] = 0.0;
	}
	for(int i = 0; i < DimCpt; i++) 
		resGradB[i] = 0.0;
	
	double resinftyW = 1.0;
	double resinftyB = 1.0;
	double **RhoReal,  **rhorealnew;
    RhoReal = (double **)malloc(sizeof(double*) *Nspecies);
    rhorealnew = (double **)malloc(sizeof(double*) *Nspecies);
     for(int i = 0; i < Nspecies; i++)
    {
        RhoReal[i] = (double *)malloc(sizeof(double)*2*cplxDofs);
        rhorealnew[i] = (double *)malloc(sizeof(double)*(local_n0*NCpt[1]*NCpt[2]));
    }
   
	do{
		getGsquare();
		iterator++; 
		timer.reset();
		timer.start();
		updatePropagator();
		timer.pause();

		if (myrank==0)
		printf("\t\t time cost of updatePropagator : %f seconds\n", timer.get_current_time());
		
		MPI_Barrier(MPI_COMM_WORLD);

		timer.reset();
		timer.start();
		singQ[0]=updateQ();
		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();
		if (myrank==0)
        printf("\t\t time cost of updateQ          : %f seconds\n", timer.get_current_time());
		
		timer.reset();
		timer.start();
		updateOrderParameter();

		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();
		if (myrank==0)
        printf("\t\t time cost of updateOrderParameter    : %f seconds\n", timer.get_current_time());
     
		timer.reset();
		timer.start();
		updateField(resGradW,resGradB);

		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();
		if (myrank==0)
        printf("\t\t time cost of updateField      : %f seconds\n", timer.get_current_time());
		
		timer.reset();
		timer.start();
		Energy = updateHamilton();


		MPI_Barrier(MPI_COMM_WORLD);
		timer.pause();
		if (myrank==0)
        printf("\t\t time cost of updateHamilton   : %f seconds\n", timer.get_current_time());

		if (myrank==0)
		{
			resinftyW = normRealInfty(resGradW, Nspecies);
			resinftyB  = normRealInfty(resGradB, DimCpt);
		}
		// save data
		for(int i=0;i< Nspecies; i++)
			FftwC2R(rhoreal[i], rho[i]);

		for(int i=0;i< Nspecies; i++)
			cut(rhoreal[i], rhorealnew[i]);
 
		n = (int *)malloc(sizeof(int)*nprocs);

		a = local_n0*NCpt[1]*NCpt[2];
		MPI_Allgather(&a, 1, MPI_INT, n, 1, MPI_INT, MPI_COMM_WORLD);

		displs= (int *)malloc(sizeof(int)*nprocs);//每个进程的偏移量
		displs[0] = 0;
		for(int i=1;i< nprocs; i++)
		{
			displs[i] = displs[i-1]+ n[i-1];
		}

		for(int i=0;i< Nspecies; i++)
			MPI_Gatherv(rhorealnew[i], a, MPI_DOUBLE, RhoReal[i], n, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		diffEnergy = fabs(Energy-oldEnergy);
		if(iterator==1 ||iterator%500== 0)
		{
			if (myrank == 0)
				write_rho(RhoReal, iterator);//	
		}
		if (myrank == 0)
		{
			   char fname[255];
			   sprintf(fname, "./results-%d/rst1.[%.4f.%.4f]-[%.2f].dat", phase, fA, fB, chi);
			   FILE *fp=fopen(fname,"a");
	           fprintf(fp, "%.4f\t %.4f\t  %.6e\t  %.6e\t  %.6e\t   %.7e\t  %.7e\t  %.7e\t  %.7e\t\n", dirBox[0][0], dirBox[1][1], diffEnergy, resinftyW, resinftyB, singQ[0], Energy, internalEnergy, entropicEnergy1);
			   fclose(fp);
		}
	  	if (myrank==0)
        {
			printf("ITERATOR %d:  singQ = %.15e\t, Energy=%.15e\t,  diffhm =%.15e, resinftyW = %.15e,  resinftyB=%.15e\n", iterator, singQ[0],  Energy, diffEnergy, resinftyW,  resinftyB);	
            printf( "%.5f\t %.5f\t %.15e\t %.15e\n", dirBox[0][0], dirBox[1][1], internalEnergy, entropicEnergy1);
        }	
        oldEnergy = Energy;
		
	    if (iterator > ItMax)
			break;
           
		if(diffEnergy > 10000)
			break;
	}while(diffEnergy> tol);
	if (myrank == 0)
	{
       		write_rho(RhoReal,  iterator);
	}
	if (myrank ==0)
	{
		char fname2[255];                                                           
		sprintf(fname2, "./results-%d/Data/Data.[%.4f.%.4f]-[%.2f].dat", phase, fA, fB, chi);
    	FILE *fp2=fopen(fname2,"w");     
		fprintf(fp2," %ld\t %ld\t  %.5f\t %d\t  %d\n", NCpt[0], NCpt[1], dsMax, phase, iterator);
		fprintf(fp2,"%.4f\t %.4f\t  \n", fA, fB);
	    fprintf(fp2,"%.2f\t \n", chi);
		fprintf(fp2,"%.4f\t %.4f\n", dirBox[0][0], dirBox[1][1]);
	    fprintf(fp2, "%.15e\t %.15e\n",diffEnergy, singQ[0]);
	    fprintf(fp2, "%.15e\t %.15e\t\n",resinftyW,  resinftyB);
        fprintf(fp2, "%.15e\t %.15e\t %.15e\n",Energy, internalEnergy, entropicEnergy1);
        fclose(fp2);  
	}

	free(resGradW);
	free(resGradB);
	for(int i = 0; i < Nspecies; i++)
	{
		free(RhoReal[i]);
		free(rhorealnew[i]);
	}
	free(RhoReal);
	free(rhorealnew);
   	return Energy;
}


void updateField(double *resGradW, double *resGradB)
{
	SimpleMixing(resGradW, resGradB);
}

double calPartialF()
{
	double partialF = 0;
	getGsquare();
	updatePropagator();
	double Q = updateQ();
	partialF -= std::log(Q);
	return partialF;
}
void SCFTiteration(double tol)
{
	evlSaddle(tol);

}
\end{lstlisting}
