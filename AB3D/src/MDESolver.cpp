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
                                     fftw_complex *hatW,
                                     fftw_complex **hatQ,
									 int *index)
//  second-order operator-splitting scheme
{
        double ds = ((*range).b - range->a) / range->n;
      
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
        FftwC2R(WReal,hatW);
   
        //  for(int i = 0; i < range->n; i++)
         for(int i = 0; i <range->n; i++)
        {
                // for(int t = 0; t < 2; t++)
                for(int t = 0; t < 1; t++)
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
                FuncsLinear1Cplx(hatQ[*index+i+1],  alloc_local, 1.0, tmp[0]);
            
        }
        *index += range->n;
       for(int i = 0; i < 3; i++) fftw_free(temp[i]);
        fftw_free(temp);
	        for(int i = 0; i < 2; i++) fftw_free(tmp[i]);
        fftw_free(tmp);
        free(WReal);
        free(Qreal);


}


void MDESolver2OrderExtrapolation(Interval *range,
                                     fftw_complex *hatW,
                                     fftw_complex **hatQ,
									 int *index)
//  second-order operator-splitting scheme
{
        double ds = ((*range).b - range->a) / range->n;
      
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
        FftwC2R(WReal,hatW);
        	// WReal
	// for(int i = 0; i < 20; i++)
	// {
	// 	printf("WReal %.20f\n", WReal[i]);
	// }
        // printf("=======================\n");
        //  for(int i = 0; i < range->n; i++)
         for(int i = 0; i <range->n; i++)
        {
                // for(int t = 0; t < 2; t++)
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
        *index += range->n;
       for(int i = 0; i < 3; i++) fftw_free(temp[i]);
        fftw_free(temp);
	        for(int i = 0; i < 2; i++) fftw_free(tmp[i]);
        fftw_free(tmp);
        free(WReal);
        free(Qreal);


}


void MDESolver2OrderExtrapolationprintf(Interval *range,
                                     fftw_complex *hatW,
                                     fftw_complex **hatQ,
									 int *index)
//  second-order operator-splitting scheme
{
        printf("MDESolver2Orderprintf\n");
        bool ifprintf=0;

        double ds = ((*range).b - range->a) / range->n;
      
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
         if(ifprintf)
         {
                for(int j = 0; j < 20; j ++)
                {
                        printf("hatW %.20f %.20f\n",hatW[j][0], hatW[j][1]);
                }
                printf("========================================\n");
         }      
        

        FftwC2R(WReal,hatW);
         if(ifprintf)
         {
                for(int j = 0; j < 50; j ++)
                {
                        printf("Wreal %.20f \n",WReal[j]);
                }
                printf("========================================\n");
                for(int k = 0; k < 20; k ++)
                {
                        printf("init data %.20f %.20f\n",hatQ[*index][0],hatQ[*index][1]);
                }
                printf("========================================\n");
         }

        
        //  for(int i = 0; i < range->n; i++)
         for(int i = 0; i <range->n; i++)
        {
                // printf("i = %d\n",i);
                // for(int t = 0; t < 2; t++)
                for(int t = 0; t < 2; t++)
                {
                        memcpy(tmp[t], hatQ[*index+i], sizeof(fftw_complex)*alloc_local);
                        double P = std::pow(2.0, t);
                        // printf("P = %f\n",P);
                        for(int index1 = 0; index1 < P; index1++)
                        {
                                if(ifprintf)
                                {
                                        for(int k = 0; k < 20; k ++)
                                        {
                                                printf("Gsquare %.20f\n",Gsquare[k]);
                                        }
                                        printf("========================================\n");

                                }
                                 
                                for(int k = 0; k < alloc_local; k ++)
                                {
                                        tmp[t][k][0] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                        tmp[t][k][1] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                }
                                if(ifprintf)
                                {
                                            for(int k = 0; k < 20; k ++)
                                        {
                                                printf("tmp1 %.20f %.20f\n",tmp[t][k][0],tmp[t][k][1]);
                                        }
                                        printf("========================================\n");

                                }
                            

                                FftwC2R(Qreal, tmp[t]);

                                if(ifprintf)
                                {
                                        for(int j = 0; j < 20; j ++)
                                        {
                                                printf("Qreal1 %.20f \n",Qreal[j]);
                                        }
                                        printf("========================================\n");

                                }

                                
                                
                                for(int j = 0; j < 2*alloc_local; j ++)
                                        Qreal[j] *= std::exp(-WReal[j]*ds/P);

                                if(ifprintf)
                                {
                                        // int local_n0n=NCpt[0];
                                        // int M=NCpt[1];
                                        // int N=NCpt[2];
                                        // printf("local_n0n %d M %d N %d \n",local_n0n,M,N);      
                                        // int mm=local_n0n*M*(N/2+1)*2;
                                        // mm=mm-2;
                                        // for(int ii=0;ii<20;ii++)
                                        // {
                                        // printf("Qreal2 = %.20g \n",Qreal[mm-1-ii]);
                                        // }
                                          for(int j = 0; j < 20; j ++)
                                        {
                                                printf("Qreal2 %.20f \n",Qreal[j]);
                                        }
                                        printf("========================================\n");
                                }


                                FftwR2C(tmp[t], Qreal); 
                                 
                                if(ifprintf)
                                {
                                        for(int k = 0; k < 100; k ++)
                                        {
                                                printf("tmp2 %.20f %.20f\n",tmp[t][k][0],tmp[t][k][1]);
                                        }
                                        printf("========================================\n");

                                }


                                for(int k = 0; k < alloc_local; k ++)
                                {
                                        tmp[t][k][0] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                        tmp[t][k][1] *= std::exp(-Gsquarelocal[k]*ds/(2.0*P));
                                }
                                if(ifprintf)
                                {
                                         for(int k = 0; k < 20; k ++)
                                        {
                                                printf("tmp3 %.20f %.20f\n",tmp[t][k][0],tmp[t][k][1]);
                                        }
                                        printf("========================================\n");

                                }
                        }

                }

                FuncsLinear2Cplx(hatQ[*index+i+1],  alloc_local, -1.0/3.0, tmp[0], 4.0/3.0, tmp[1]);
                ifprintf=1;
                if(ifprintf)
                {
                        printf("%d\n",i);
                        for(int k = 0; k < 20; k ++)
                        {
                                printf("tmp0 %.20f %.20f\n",tmp[0][k][0],tmp[0][k][1]);

                        }
                        printf("========================================\n");
                        for(int k = 0; k < 20; k ++)
                        {
                                printf("tmp1 %.20f %.20f\n",tmp[1][k][0],tmp[1][k][1]);
                        }
                        printf("========================================\n");
                         for(int k = 0; k < 20; k ++)
                        {
                              
                                printf("final %.20f %.20f\n",hatQ[*index+i+1][k][0],hatQ[*index+i+1][k][1]);
                        }
                        printf("========================================\n");
                }
                ifprintf=0;

        }

        *index += range->n;
       for(int i = 0; i < 3; i++) fftw_free(temp[i]);
        fftw_free(temp);
	        for(int i = 0; i < 2; i++) fftw_free(tmp[i]);
        fftw_free(tmp);
        free(WReal);
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
                // for(int i=0;i<500;i++)
                // {
                //         printf("fieldW[0]=%.20f %.20f\n",hatW[i][0],hatW[i][1]);
                // }
                // printf("--------------------------------");
		// FftwC2R(WReal,hatW);

                // for(int i=0;i<500;i++)
                // {
                //         printf("wReal=%.20f\n",WReal[i]);
                // }
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


