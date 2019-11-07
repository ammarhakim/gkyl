#include <algorithm> 
#include <cmath>
#include <stdio.h>
#define EPSILON std::numeric_limits<double>::epsilon()
extern "C" { 
double findMinNodalValue(const double *fIn, int ndim); 
double findMinNodalRatio(const double *fNum, const double *fDenom, double fac, int ndim);
double rescale(const double *fIn, double *fOut, int ndim, int numBasis, int *idx, double tCurr);
double rescaleVolTerm(const double *fOutSurf, const double dt, double *fVol, int ndim, int numBasis, int *idx);
bool check(const double *fIn, int ndim, int numBasis, int *idx, double tCurr, int rkIdx);
}; 
