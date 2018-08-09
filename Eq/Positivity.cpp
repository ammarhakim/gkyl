/* define C functions needed for positivity implementation */
#include <Positivity.h>

double patchFit(double r, double x, double CFL) {
   double val = 0.0;
   if (x > 0) {
      if (r<2.2) {
         val = fmax(0.0, fmin(1.0/CFL, std::exp(2*r*x/3)*(1+r*x/3)));
      } else {
         val = fmin(1.0/CFL, 6/(3-fmin(2.999, std::abs(r))));
      }
   } else {
      if (r>-2.2) {
         val = fmax(0.0, fmin(1.0/CFL, std::exp(2*r*x/3)*(1+r*x/3)));
      } else {
         val = fmin(1.0/CFL, 6/(3+fmin(2.999, std::abs(r))));
      }
   }
   return val;
}

double limTheta(double r, double x, double CFL) {
#if extraType == NONE
   return 1 + r*x;
#elif extraType == LINEAR
   return fmax(0.0, 1+r*x);
#elif extraType == EXP
   return fmin(1.0/CFL, std::exp(r*x));
#elif extraType == EXP0
   return fmax(0.0, fmin(1.0/CFL, std::exp(2*r*x/3)*(1+r*x/3)));
#elif extraType == PATCHFIT
   return patchFit(r, x, CFL);
#else
   return 0.0;
#endif
}
