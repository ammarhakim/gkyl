/* define C functions needed for positivity implementation */
#include <Positivity.h>
#include <iostream>
#define extraType PATCHFIT

double patchFit(double r, double x) {
   double val = 0.0;
   if (x*r<2.2) {
      val = std::exp(2.*r*x/3.)*(1.+r*x/3.);
   } else {
      val = 6.0/(3.0-std::min(2.999, std::abs(r)));
   }
   return val;
}

double limTheta(double r, double x) {
#if extraType == NONE
   return 1 + r*x;
#elif extraType == EXP
   return std::exp(r*x);
#elif extraType == EXP0
   return std::exp(2*r*x/3)*(1+r*x/3);
#elif extraType == COHEN
   long double beta = r/3*(3-r*r/9)/(1-std::min(0.999,r*r/9));
   return (beta*r/3+1)*std::exp(beta*x)/std::cosh(beta);
#elif extraType == KROGER
   double r3 = r/3.;
   long double beta = (3.*r3-r3*r3*r3*(6. + r3*r3 - 2.*r3*r3*r3*r3)/5.)/(1.-std::min(0.999,r3*r3));
   return (beta*r/3+1)*std::exp(beta*x)/std::cosh(beta);
#elif extraType == PATCHFIT
   return patchFit(r, x);
#else
   return 0.0;
#endif
}
