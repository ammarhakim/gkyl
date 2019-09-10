#include <math.h>
#include <fpoKernelsDecl.h>

void fpoDiagKernelP1(const double *dv, const double *vc, const double *f, const double *h, double *out) {
  out[0] += 0.5*(1.732050807568877*dv[1]*f[2]*h[3]+1.732050807568877*f[0]*dv[1]*h[1]);
  out[1] += 0.5*(1.732050807568877*dv[0]*f[1]*h[3]+1.732050807568877*dv[0]*f[0]*h[2]);
  out[2] += 0.125*((5.0*dv[0]*dv[1]*f[3]+6.928203230275509*vc[0]*dv[1]*f[2]+6.928203230275509*dv[0]*f[1]*vc[1])*h[3]+(3.0*dv[0]*dv[1]*f[2]+6.928203230275509*dv[0]*f[0]*vc[1])*h[2]+(3.0*dv[0]*dv[1]*f[1]+6.928203230275509*f[0]*vc[0]*dv[1])*h[1]+dv[0]*f[0]*h[0]*dv[1]);
}