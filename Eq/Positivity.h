/* define C functions needed for positivity implementation */
#include <cmath>
#include <limits>
#define extraType PATCHFIT
#define EPSILON std::numeric_limits<double>::epsilon()

double patchFit(double r, double x, double CFL);
double limTheta(double r, double x, double CFL);
