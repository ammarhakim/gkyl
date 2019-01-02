#include <algorithm> 
#include <cmath>
#define EPSILON std::numeric_limits<double>::epsilon()
extern "C" { 
double findMinNodalValue(const double *fIn, int ndim); 
double rescale(const double *fIn, double *fOut, int ndim, int numBasis, int *idx, double tCurr);
}; 
