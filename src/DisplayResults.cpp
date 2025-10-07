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

