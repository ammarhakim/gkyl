/* define C functions needed for positivity implementation */
#include <cmath>
#include <limits>
#define NONE 0
#define LINEAR 1
#define EXP 2
#define EXP0 3
#define PATCHFIT 4
#define EPSILON std::numeric_limits<double>::epsilon()

double patchFit(double r, double x);
double limTheta(double r, double x);
