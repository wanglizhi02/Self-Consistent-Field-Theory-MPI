#include "mpi.h"                                                                
#include "fftw3-mpi.h" 

#include "Head.h"
#include "Data.h"
#include "Initialization.h"
#include "BasicFunc.h"
#include "MemFree.h"
#include "SCFTBaseAB.h"
#include "FftwToolkit.h"
#include "Mytimer.h"

#include "config.h"

#include <iostream>
#include <string>

#include "time.h"

void fixedParameters()
{
	// Computational dimensions
	DIM = 3;
	DimPhy = 3; //real space;
	DimCpt = 3; //cplx space;

	NCpt = (ptrdiff_t *)malloc(sizeof(ptrdiff_t)*DimCpt);
	// for(int i = 0; i < DimCpt-1; i++)
	// 	NCpt[i] = 414;

	// NCpt[DimCpt-1] = 414;
	int pp=64;
	NCpt[0] = pp;
	NCpt[1] = pp;
	NCpt[2] = pp;
	// NCpt[0] = 510;
	// NCpt[1] = 76;
	// NCpt[2] = 76;
	
	// 总共的自由度
	realDofs = cplxDofs = 1;
	for(int i = 0; i < DimCpt-1; i++)
	{
		realDofs *= NCpt[i];
		cplxDofs *= NCpt[i];
	}
	realDofs *= NCpt[DimCpt-1];
	cplxDofs *= (NCpt[DimCpt-1]/2+1);
	
	if (myrank == 0)
	{
	printf("\t ------- Discrete Modes ------- \n");
	for(int i = 0; i < DimCpt; i++)
		printf("\t NCpt[%d] = %ld", i, NCpt[i]);
	printf("\n");
	printf("\t cplxDofs = %d,\t realDofs = %d\n",  cplxDofs, realDofs);
	}

	Nspecies = 2;
	Nblend = 1;
	Nblock = 2;



}
void Distribute(double *rho, double *rholocal)
{
    for (int i = 0; i < local_n0; i++)
      {
         for (int j = 0; j < NCpt[1]; j++)
         {
         	for (int k = 0; k < NCpt[2]; k++)
			{
             	int index = i*NCpt[1]*(2*(NCpt[2]/2+1))+j*(2*(NCpt[2]/2+1))+k;
             	int index1 = (i+local_0_start)*NCpt[1]*NCpt[2]+j*NCpt[2]+k;
              	rholocal[index] = rho[index1];
			}
         }
      }
}

void initRho()
{

	mytimer_t timer;

	// allocate memory 
	timer.reset();
	timer.start();
	
	//init field
    char fname[150];
	sprintf(fname, "/beegfs/home/phasetransition/wlz/AB3D/initData/DG/DGA_64.txt");
	// sprintf(fname, "/beegfs/home/phasetransition/wlz/AB3D/initData/rho_DG_128_128_128.txt");
	writeRealData(rhoglobal[0], fname);
	
	Distribute(rhoglobal[0], rhoreal[0]);
	
	FftwR2C(fieldW[0], rhoreal[0]);
	FuncsLinear1Cplx(fieldW[1], alloc_local, -1.0, fieldW[0]);
	
	FftwC2R(rhoreal[1],fieldW[0]);
	
	timer.pause();
	if (myrank==0)
		printf("\t\t time cost of initRho : %f seconds\n", timer.get_current_time());
	

} 
void initdirBox()
{
	dirBox[0][0] = 8.5;
	dirBox[1][1] = 8.5;
	dirBox[2][2] = 8.5;
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv); ///初始化                                           
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);//获取总进程数                       
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);//获取当前进程的进程号              

	if (myrank==0)
		printf("\n\n ==============================   START PROGRAM ======================================\n\n");

	phase = atoi(argv[1]);
	fA =atof(argv[2]);
	fB = 1.0-fA;
	chi = atof(argv[3]);
	double dN = atof(argv[4]); 
    dsMax = 1.0/dN;

	ItMax = atof(argv[5]);//Max Iteration 
	Ndeg = atoi(argv[6]);
	double TOL = 1.0e-08;
    
	mytimer_t timer;
	timer.reset();
	timer.start();
	long int timeSum;
    // time_t start, finish;
    // start = time(NULL);

	fixedParameters();
	
	fftw_mpi_init();	
	
	//  Allocate memory  --------------------------------------------
	alloc_local = fftw_mpi_local_size_3d(NCpt[0], NCpt[1], NCpt[2]/2+1, MPI_COMM_WORLD, &local_n0, &local_0_start);
	//alloc_local = fftw_mpi_local_size(DIM, NCpt, MPI_COMM_WORLD, &local_n0, &local_0_start);

	printf("n0[%d]=%ld, %ld\n", myrank, local_n0, alloc_local);
	MPI_Barrier(MPI_COMM_WORLD);

	
	initialize();
	//  Projective matrix --------------------------------------------
	for(int i = 0; i < DimPhy; i++)
		ProjMatrix[i][i] = 1.0;

   	//   Initial fields --------------------------------------
  	initRho();
  
    //  Computational domain ----------------------------------------
    initdirBox();
    getRecipLattice(dirBox, rcpBox);
	// MatPrint(rcpBox, DimCpt, DimCpt);
    memAllocation();
    if (myrank==0)
    {
		// Model parameters ---------------------------------------------
		printf("\n\n\t********************************* PARAMETERS ************************************* \n");
		printf("\t Discrete points in s:[Nsplit:A-B]=[%d-%d]\n", Nsplit[0],Nsplit[1]);
		printf("\t*\t\t   Nspecies = %d, \t Nblend = %d, \t Nblock = %d              *\n", Nspecies, Nblend, Nblock);
		printf("\t*\t\t[chiAB-fA-fB] = [%f-%f-%f]      *\n", chi, fA, fB);
		printf(" \t cplxDofs = %d,\t realDofs = %d, \t phase = %d\n", cplxDofs, realDofs, phase);
    }
	//   SCFT Iteration ----------------------------------------------
    SCFTiteration(TOL);
	
	if (myrank ==0)
	{
	char fname2[255];                                                            
	sprintf(fname2, "./results-%d/energy/energy.[phase=%d].[fB=%.2f]-chiAB=%.2f.dat", phase, phase, fB, chi);
	FILE *fp2=fopen(fname2,"a");     
    fprintf(fp2, "%.3f \t %.15e\t %.15e\t %.15e\n",fB,internalEnergy, entropicEnergy1, Hamilton);
	fclose(fp2);  
	}
//  Releasr memory  ----------------------------------------------
	memReleaser();
//  --------------------------------------------------------------
	// finish = time(NULL);
    // timeSum = finish - start;
	// if (myrank==0)
	// {
	// 	printf("\n=== Time cost: %ld seconds, \t %f minutes, \t %f hours\n\n", timeSum, double(timeSum/60.0), double(timeSum/3600.0));
	// 	printf("\n\n ======================================   END PROGRAM  ======================================\n\n");
	// }
	timer.pause();
	 if (myrank ==0){
		double mytime=timer.get_current_time();
		printf("\n===Mytime Time cost: %f seconds, \t %f minutes, \t %f hours\n\n", mytime, double(mytime/60.0), double(mytime/3600.0));
		printf("\n\n ======================================   END PROGRAM  ======================================\n\n");
	}


	MPI_Finalize();
	
	return 1;
}

