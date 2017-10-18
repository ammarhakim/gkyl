#ifndef GKYL_TEST_CFUNCS_H
#define GKYL_TEST_CFUNCS_H

extern "C"
{
    typedef struct { double x, y, z; } loc_t;
    double calcSum(int n, double *v);
    double addValues(loc_t *v);
    void setValues(int n, int ix, double *v);
    void allocValues(int sz, void *v);
}

#endif // GKYL_TEST_CFUNCS_H
