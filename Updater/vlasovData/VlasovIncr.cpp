#include <VlasovIncrDecl.h>

void vlasovIncr(unsigned n, const double *aIn, double a, double *aOut) {
  for (unsigned i=0; i<n; ++i)
    aOut[i] = a*aIn[i];
}
