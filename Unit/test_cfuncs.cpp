#include <test_cfuncs.h>
#include <stdlib.h>

double calcSum(int n, double *v)
{
  double sum = 0.0;
  for (unsigned i=0; i<n; ++i)
    sum += v[i];
  return sum;
}

double addValues(loc_t *v)
{
  return v->x + v->y + v->z;
}

void setValues(int n, int ix, double *v)
{
  v[0]  = ix+1;
  v[1]  = ix+2;
  v[2]  = ix+3;
}

void allocValues(int sz, void *v)
{
  v = malloc(sz);
}
