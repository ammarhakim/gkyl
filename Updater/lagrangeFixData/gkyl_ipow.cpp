#include <VlasovLagrangeFixDecl.h> 

double gkyl_ipow(double base, int exp) {
  double result = 1;
  int sign = 1;
  if (exp < 0) {
    sign = -1;
    exp = -1*exp;
  }
  while (1) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    if (!exp)
      break;
    base *= base;
  }
  if (sign == -1)
    result = 1/result;
  return result;
}
