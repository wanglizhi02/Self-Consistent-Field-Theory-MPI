#ifndef __IntegrationToolkit_h
#define __IntegrationToolkit_h

#include <string.h>
#include "Data.h"

extern void integration(fftw_complex *phi, 
					 Interval *range, 
					 fftw_complex **integrand,
					 int *index,
					 double singQ);
#endif
