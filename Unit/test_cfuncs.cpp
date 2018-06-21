#include <test_cfuncs.h>
#include <classwrap.h>
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

void* new_Adder(int n)
{
  Adder *a = new Adder(n);
  return reinterpret_cast<void*>(a);
}

int add_Adder(Adder *a, int x)
{
  return a->add(x);
}

int sub_Adder(Adder *a, int x)
{
  return a->sub(x);
}

void setFuncPointer_Adder(Adder *a, int (*v)(int))
{ 
  a->setFuncPointer(v);
}

int incr_Adder(Adder *a, int x)
{
  return a->incr(x);
}

void
setMetricFuncPointer_Adder(Adder *a, void (*gfunc)(double *xc, double *g))
{
  a->setMetricFuncPointer(gfunc);
}

void
printg_Adder(Adder *a, double r, double t)
{
  a->printg(r, t);
}
