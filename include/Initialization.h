#ifndef __Initialization_h
#define __Initialization_h

extern void initialize();
extern void memAllocation();
extern void memAllocation_no_time();
extern void getGK();
extern void getGsquare();
extern void getRecipLattice(double **dBox, double **rBox);
extern void writeRealData(double *field, const char *fname);
extern void initFieldFourier(fftw_complex *field, const char *fname);
#endif
