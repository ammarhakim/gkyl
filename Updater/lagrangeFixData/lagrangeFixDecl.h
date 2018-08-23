#ifndef LAGRANGE_FIX_H 
#define LAGRANGE_FIX_H 

extern "C" { 

double gkyl_ipow(double base, int exp);

void lagrangeFixSer1x1v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

void lagrangeFixMax1x1v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f);

} 
#endif