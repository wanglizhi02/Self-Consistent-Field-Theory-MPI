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

