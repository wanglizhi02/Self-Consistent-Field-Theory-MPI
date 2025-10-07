#ifndef __DisplayResults_h
#define __DisplayResults_h

#include "Head.h"
extern void dispDensityPara(fftw_complex **Csrc, std::string &densityName, int step);
extern void dispDensitycoil(Interval *range, fftw_complex **coil, std::string &densityName, int step);
extern void dispDensityrod(Interval *range, fftw_complex ***rod, std::string &densityName, int step);
extern void dispDensity(fftw_complex *Csrc, std::string & densityName, const int step);
extern void dispDensity_Lx(fftw_complex *Csrc, std::string & densityName, const int step);
extern void dispDensity_Ly(fftw_complex *Csrc, std::string & densityName, const int step);
extern void dispCplxDensity(fftw_complex *Csrc, int **kspace, int n, int dim, std::string & name);
extern void dispRealDensity(double *Rsrc, int **kspace, int n, int dim, std::string & name);
extern void dispPlaneWave(fftw_complex *Csrc, int n, std::string & vecName);
extern bool myComp(mySortVec a, mySortVec b);


#endif
