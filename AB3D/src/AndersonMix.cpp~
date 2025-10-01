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
