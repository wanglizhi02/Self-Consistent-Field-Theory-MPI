#ifndef __MDESolver_h
#define __MDESolver_h

#include "Data.h"
//extern double getUdotK(int indexu, int indexr);

extern void MDESolver2Order(Interval *range,  
					 double *WReal, 
					 double **Q);

extern void MDESolver4Adams(Interval *range,
                                     fftw_complex *hatW,
                                     fftw_complex **hatQ,
									 int *index);
#endif
