#include <algorithm> 
#include <cmath>
#define EPSILON std::numeric_limits<double>::epsilon()
extern "C" { 
double limiter(const double *fIn, const double *fHat, double *fOut, int ndim, int numBasis);
}; 
