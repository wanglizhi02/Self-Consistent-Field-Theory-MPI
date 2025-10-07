\subsubsection{ab.h}
\begin{lstlisting}
#ifndef ab_h
#define ab_h
#include "parameter.h"

void fixedParameters();

void write_rho(double **RhoReal, int iter, double chiAB, int phase, int realDofs, double mu);

void write_hatrho(fftw_complex *rho, double mu);

void initRho(int phase);


#endif // end of ab_h
\end{lstlisting}

\subsubsection{AndersonMix.h}
\begin{lstlisting}
#ifndef __AndersonMix_h
#define __AndersonMix_h

// andersion mixing operators
extern double andersonMixing(double tol);
extern void updateMatrix(double **srcMatrix, fftw_complex ***dw, fftw_complex ***dM, double ***dL, int n);
extern void assembAnderMatrix(double **srcMatrix, double **rsltMatrix, double *rigCoeff, int n);
extern void anderInnProd(double **paraD, double *paraV, double ***gradWold, int n);
//extern void assembAnderMatrix(double **matrix, double *rht, double ***f, double ***df, int n);
extern void getAnderCoeff(double *u, double **A_matrix, double *rhs, int n);

#endif
\end{lstlisting}

\subsubsection{BasicFunc.h}
\begin{lstlisting}
#ifndef __BasicFunc_h
#define __BasicFunc_h

#include "Head.h"

extern void printVecCplx(fftw_complex *src);
extern double normCplxInfty(fftw_complex *src, int n);
extern double normRealInfty(double *src, int n);
extern void eliminate(fftw_complex *src, int n);

extern void FuncsLinear1Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1);

extern void FuncsLinear2Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2);

extern void FuncsLinear3Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3);

extern void FuncsLinear4Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4);

void FuncsLinear5Cplx(fftw_complex *rslt, int n,
					  const double a1, const fftw_complex *F1,
					  const double a2, const fftw_complex *F2,
					  const double a3, const fftw_complex *F3,
					  const double a4, const fftw_complex *F4,
					  const double a5, const fftw_complex *F5);

extern void FuncAddToCplx(fftw_complex *rslt, int n,
			   const double a1, const fftw_complex *F1);

extern void FuncCplxAddAConst(fftw_complex *rslt, int n, const double a);

extern void setCplxZero(fftw_complex rslt);

extern void FuncsLinear1Real(double *rslt, int n,
					  const double a1, const double *F1);

extern void FuncsLinear2Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2);

extern void FuncsLinear3Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3);

extern void FuncsLinear4Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4);

extern void FuncsLinear5Real(double *rslt, int n,
					  const double a1, const double *F1,
					  const double a2, const double *F2,
					  const double a3, const double *F3,
					  const double a4, const double *F4,
					  const double a5, const double *F5);


extern void FuncAddToReal(double *rslt, int n,
				   const double a1, const double *F1);

extern void MatCopy(double **dst, double **src, int n, int m);
extern void MatsAdd(double **dst, double d, double **src, double s, int n, int m);
extern double MatsDot(double **lhs, double **rhs, int n, int m);
extern void MatPrint(double **matrix, int n, int m);
extern double innerProduct(double *lft, double *rht, int n);

#endif
\end{lstlisting}

\subsubsection{config.h}
\begin{lstlisting}
#ifndef config_h
#define config_h

#define DATA_DIR  "@PROJECT_SOURCE_DIR@/initData/"
#define RESULT_DIR "@PROJECT_SOURCE_DIR@/result/"
#define INPUT_DIR "@PROJECT_SOURCE_DIR@/input/"

#endif // end of config_h
\end{lstlisting}

\subsubsection{Data.h}
\begin{lstlisting}
#ifndef __Data_h
#define __Data_h

#include "Head.h"

typedef struct{
	double a, b;
	int n;
} Interval;

typedef struct myVec
{
	double* Data;
	int* Index;
}mySortVec;

extern int nprocs, myrank;
extern double chi;
extern int phase;
extern int ItMax;
extern int DIM, DimPhy, DimCpt;
extern double dsMax;
extern double Dx, Dy;
extern int Ndeg;
extern ptrdiff_t *NCpt;
extern int cplxDofs, realDofs;
extern int Nblock, Nspecies, Nblend;
extern ptrdiff_t alloc_local, local_n0, local_0_start;
extern int *n, *displs;
extern int a;
extern double *singQ;
extern double **dirBox, **rcpBox, **gradB;
extern double **ProjMatrix;
extern int **indKspace, **indKspacelocal;
extern double *Gsquare;
extern double *Gsquarelocal;
extern fftw_complex **rho;
extern double **rhoglobal, **rhoreal;
extern fftw_complex **fieldW, *fieldWplus;
extern fftw_complex **gradW;
extern fftw_complex *fftw_Ctmp;
extern double *fftw_Rtmp;
extern fftw_plan PlanC2R, PlanR2C;
extern double fA, fB;

extern int *Nsplit;
extern int NA, NB, Ns;
extern Interval *Ranges;
extern fftw_complex **frdQ, **bakQ;
extern fftw_complex **hatq_qplus;
extern double  Hamilton;
extern double  FreeEnergy;
extern double internalEnergy, entropicEnergy1, entropicEnergy2; 
#endif
\end{lstlisting}

\subsubsection{DisplayResults.h}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{FftwToolkit.h}
\begin{lstlisting}
#ifndef __FftwToolkit_h
#define __FftwToolkit_h

#include "Head.h"

//extern void getIndextheta(int *kspace, int n, int ndeg);
extern void getIndex(int **kspace, int n, int dim, ptrdiff_t *ndeg);
extern void FftwC2R(double *Rrslt, fftw_complex *Corig);
extern void FftwR2C(fftw_complex *Crslt, double *Rorig);
//extern void FftwC1R(double *Rrslt, fftw_complex *Corig);
//extern void FftwR1C(fftw_complex *Crslt, double *Rorig);
extern double hatVecsMultiply(fftw_complex *lhs, fftw_complex *rhs);
extern double InnerProd(fftw_complex *lhs, fftw_complex *rhs, int cplxNtheta, int n);
extern void hatConv(fftw_complex *rslt, fftw_complex *src1, fftw_complex *src2);
extern void Conv(fftw_complex *rslt, double *src1, double *src2);
extern double Intergral_space(fftw_complex *src1, fftw_complex *src2);
extern double Intergral_space_real(double *src1, double *src2);
//extern double Intergral_u(double *src1, double *src2);

extern double innerProduct_space(fftw_complex *lft, fftw_complex *rht);
#endif
\end{lstlisting}

\subsubsection{GaussElimination.h}
\begin{lstlisting}
#ifndef __GaussElimination_h
#define __GaussElimination_h

extern void GaussElimination(double *u, double **A_matrix, double *rhs, int n);
extern void NumberSwap(double *a, double *b);
extern void CheckIsHaveSolution(double **A_matrix, int n, int k);
extern void Pivot(double **A_matrix, double *Rhs, int n, int k);

extern void Solve_LU(double *u, double **A_matrix, double *rhs, int n);
#endif
\end{lstlisting}

\subsubsection{Head.h}
\begin{lstlisting}
// Macro definition

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <string.h>
#include <assert.h>

#include "fftw3.h"

#define PI 3.14159265358979323846
#define HELLO printf("Hello!\n");
\end{lstlisting}

\subsubsection{Initialization.h}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{IntegrationToolkit.h}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{MDESolver.h}
\begin{lstlisting}
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
\end{lstlisting}

\subsubsection{MemFree.h}
\begin{lstlisting}
#ifndef __MemFree_h
#define __MemFree_h

extern void memReleaser();

#endif
\end{lstlisting}

\subsubsection{Mytimer.h}
\begin{lstlisting}
#ifndef __mytimer_h_
#define __mytimer_h_

#include <sys/time.h>

class mytimer_t {
    private:
	long time_total;
	long time_previous;

	struct timeval time_start;

	struct timeval time_end;

    public:
	mytimer_t() : time_total(0), time_previous(0) {}

	void start() { gettimeofday(&time_start, NULL); }

	void pause() {
	    gettimeofday(&time_end, NULL);
	    time_previous = time_total;
	    time_total += (time_end.tv_sec - time_start.tv_sec) * 1000000 + (time_end.tv_usec - time_start.tv_usec);
	}

	void reset() { time_total = 0; time_previous = 0;}

	double get_current_time() { return time_total * 1e-6; }
	double get_previous_time() { return time_previous * 1e-6; }
};

#endif
\end{lstlisting}

\subsubsection{parameter.h}
\begin{lstlisting}
#ifndef parameter_h
#define parameter_h
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <stdio.h>
#include <sstream>


struct Parameter
{
    std::string pattern;
    double chiAB;
    std::array<int, 2> degree;
    double dsmax;
    double u;
    std::array<int, 3> dim;
    std::array<double, 3> fA1;
    std::array<double, 3> fA2;
    std::array<double, 9> domain;
    std::string phiA;
    std::string phiB;
    std::string filename;

    /* constructor
     * read params from file fname
     */
    Parameter(const char * fname)
    {
        filename = fname;
    }

    void read()
    {
        std::ifstream file(filename.c_str());
        std::string line;
        while(std::getline(file, line))
        {
            auto found = line.find(":");
            if(found != std::string::npos)
            {
                auto name = line.substr(0, found);
                remove_space(name);
                std::cout << name << std::endl;
                auto data = line.substr(found+1);
                if(name == "pattern")
                {
                    pattern = data;
                    remove_space(pattern);
                }
                else if(name == "dim")
                {
                    read_array(data, dim, name);
                }
                else if(name == "dsmax" )
                {
                    read_data(data, dsmax);
                }
 /*               else if(name == "u" )
                {
                    read_data(data, u);
                }
*/
                else if(name == "degree")
                {
                    read_array(data, degree, name);
                }
                else if(name == "chiAB")
                {
                    read_data(data, chiAB);
                }
                else if(name == "fA1")
                {
                    read_array(data, fA1, name);
                }
                else if(name == "fA2")
                {
                    read_array(data, fA2, name);
                }
                else if(name == "domain")
                {
                    read_array(data, domain, name);
                }
                else if(name == "phiA")
                {
                    phiA = data;
                    remove_space(phiA);
                }
                else if(name == "phiB")
                {
                    phiB = data;
                    remove_space(phiB);
                }
            }
        }
    }

    void remove_space(std::string & str)
    {
        str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
        return;
    }

    /* read a array*/
    template<typename Array>
    void read_array(std::string & line, Array & a, std::string & name)
    { 
        std::istringstream iss(line);
        int i = 0;
        while(iss >> a[i])
        {
            i++;
            if(i == a.size())
                break;
        }
        if(i < a.size())
        {
            std::cout << "Warning! " 
                << name 
                << " parameter need " 
                << a.size() 
                << " values, but you just input " 
                << i 
                << " values "
                << std::endl;
            exit(-1);
        }
    }

    /* read one values*/
    template<typename Type>
    void read_data(std::string & line, Type & val)
    {
        std::istringstream iss(line);
        iss >> val;
    }
/*
    "pattern": "HHC",
    "chiBC": 0.30,
    "chiAC": 0.30,
    "chiAB": 0.30,
    "fA": [0.2, 0.4, 100],
    "fB": [0.2, 0.4, 100],
    "domain": [4.8, 4.8, 4.8],
    "phiA": "phiA.[30.00.30.00.30.00].[0.19.0.62.0.19]-[32]-[7.57397744]-[-0.57672435]-[4.88641].8.txt",
    "phiB": "phiB.[30.00.30.00.30.00].[0.19.0.62.0.19]-[32]-[7.57397744]-[-0.57672435]-[4.88641].8.txt",
    "phiC": "phiC.[30.00.30.00.30.00].[0.19.0.62.0.19]-[32]-[7.57397744]-[-0.57672435]-[4.88641].8.txt"
*/

    /* show the struct */
    void print()
    {
        std::cout << " degree : " << degree[0] << " " <<degree[1] << std::endl;
        std::cout << "  dsmax : " << dsmax << std::endl;
//        std::cout << "      u : " << u << std::endl;
        std::cout << "    dim : " << dim[0] << " " << dim[1] << " " << dim[2] << std::endl;
        std::cout << "pattern : " << pattern << std::endl;
        std::cout << "  chiAB : " << chiAB << std::endl;
        std::cout << "     fA1 : " << fA1[0] << " " << fA1[1] << " " << fA1[2] << std::endl;
        std::cout << "     fA2 : " << fA2[0] << " " << fA2[1] << " " << fA2[2] << std::endl;
        std::cout << " domain : " 
            << domain[0] << " "
            << domain[1] << " "
            << domain[2] << " "
            << domain[3] << " "
            << domain[4] << " "
            << domain[5] << " "
            << domain[6] << " "
            << domain[7] << " "
            << domain[8] << std::endl;
        std::cout << "   phiA : " << phiA << std::endl;
        std::cout << "   phiB : " << phiB << std::endl;
        printf("test\n");
    }
};
#endif // end of parameter_h
\end{lstlisting}

\subsubsection{SCFTBaseAB.h}
\begin{lstlisting}
#ifndef __SCFTBaseAB_h
#define __SCFTBaseAB_h

extern void updatePropagator();
extern double updateQ();
extern void updateOrderParameter();
extern void updateField(double *resGradmu,  double *resGradB);
extern double updateHamilton();
//extern double tmp_H();
extern void get_gradB(double dh);
//extern void savefinals(int step, double hamilton, double size);
extern double evlSaddle(double tol);
extern double calPartialF();
extern void SCFTiteration(double tol);

//extern void semiImplicit(double *resGrad);
extern void simpleMixing(double *resGradW, double *resGradB);

extern void write_rho(double **RhoReal,  int iter);
extern void cut(double *rholocalold, double *rholocalnew);
#endif
\end{lstlisting}
