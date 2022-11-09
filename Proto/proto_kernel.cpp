#include <math.h>

// This demos how to include a C/C++ kernel in Proto and have it compiled
// by the build system.

double proto_kernel(double a, double b) {
  return a+b;
}
