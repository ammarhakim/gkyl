#include <test_cfuncs.h>


double calcSum(int n, double *v)
{
  double sum = 0.0;
  for (unsigned i=0; i<n; ++i)
    sum += v[i];
  return sum;
}
