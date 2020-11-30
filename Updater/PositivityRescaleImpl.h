#include <algorithm> 
#include <cmath>
#include <stdio.h>
#define EPSILON std::numeric_limits<double>::epsilon()
#define DBL_MIN std::numeric_limits<double>::min()
extern "C" { 
double findMinNodalValue(const double *fIn, int ndim, int polyOrder); 
double findMinNodalRatio(const double *fNum, const double *fDenom, double fac, int ndim, int polyOrder);
double rescale(const double *fIn, double *fOut, int ndim, int polyOrder, int numBasis, int *idx, double tCurr);
double calcVolTermRescale(const double tCurr, const double dt, const double *fIn, const double weight, const double *fRhsSurf, double *fRhsVol, int ndim, int polyOrder, int numBasis, int *idx, int printWarnings);
double rescaleVolTerm(const double tCurr, const double dt, const double *fIn, const double weight, const double *fRhsSurf, double *fRhsVol, int ndim, int polyOrder, int numBasis, int *idx, int printWarnings);
bool check(const double *fIn, int ndim, int polyOrder, int *idx, double tCurr, int rkIdx);
}; 
