#ifndef ab_h
#define ab_h
#include "parameter.h"

void fixedParameters();

void write_rho(double **RhoReal, int iter, double chiAB, int phase, int realDofs, double mu);

void write_hatrho(fftw_complex *rho, double mu);

void initRho(int phase);


#endif // end of ab_h
