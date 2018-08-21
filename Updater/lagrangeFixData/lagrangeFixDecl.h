#ifndef LAGRANGE_FIX_H 
#define LAGRANGE_FIX_H 

extern "C" { 
void lagrangeFixSer1x1v1p(const double *dm0, const double *dm1, const double *dm2, double *f, const double vc, const double L, const double Nv); 

void lagrangeFixMax1x1v1p(const double *dm0, const double *dm1, const double *dm2, double *f, const double vc, const double L, const double Nv); 

} 
#endif 
