#ifndef GKYL_TEST_CFUNCS_H
#define GKYL_TEST_CFUNCS_H

#include <iostream>

class Adder;

extern "C"
{
    typedef struct { double x, y, z; } loc_t;
    double calcSum(int n, double *v);
    double addValues(loc_t *v);
    void setValues(int n, int ix, double *v);
    void allocValues(int sz, void *v);

    void *new_Adder(int n);
    int add_Adder(Adder *a, int x);
    int sub_Adder(Adder *a, int x);
    void setFuncPointer_Adder(Adder *a, int (*v)(int));
    int incr_Adder(Adder *a, int x);

    void setMetricFuncPointer_Adder(Adder *a, void (*gfunc)(double *xc, double *g));
    void printg_Adder(Adder *a, double r, double t);
}

class Adder {
  public:
    Adder(int n)
      : n(n) {
    }

    int add(int x) {
      return n+x;
    }

    int sub(int x) {
      return n-x;
    }

    void setFuncPointer(int (*v)(int)) {
      f = v;
    }

    void setMetricFuncPointer(void (*gfunc)(double *xc, double *g)) {
      this->gfunc = gfunc;
    }    

    int incr(int x) {
      return n + f(x);
    }

    void printg(double r, double t) {
      double xc[3] = {0.0, r, t};
      double g[4];
      gfunc(xc, g);
      std::cout << g[1] << " " << g[2] << " " << g[3] << std::endl;
    }

  private:
    int n;
    int (*f)(int);
    void (*gfunc)(double *, double *);
};

#endif // GKYL_TEST_CFUNCS_H
