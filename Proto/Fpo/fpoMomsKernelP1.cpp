#include <math.h>
#include <fpoKernelsDecl.h>

void fpoMomsKernelP1(const double *dv, const double *vc, const double *f, double *out) {
  out[0] += 0.5*dv[0]*f[0]*dv[1];
  out[1] += 0.08333333333333333*(1.732050807568877*dv[1]*f[1]*gkyl_ipow(dv[0],2)+6.0*dv[0]*f[0]*vc[0]*dv[1]);
  out[2] += 0.08333333333333333*(1.732050807568877*dv[0]*f[2]*gkyl_ipow(dv[1],2)+6.0*dv[0]*f[0]*dv[1]*vc[1]);
  out[3] += 0.04166666666666666*(12.0*dv[0]*f[0]*dv[1]*gkyl_ipow(vc[1],2)+dv[0]*f[0]*gkyl_ipow(dv[1],3)+6.928203230275509*dv[0]*vc[1]*f[2]*gkyl_ipow(dv[1],2)+dv[1]*(12.0*dv[0]*f[0]*gkyl_ipow(vc[0],2)+f[0]*gkyl_ipow(dv[0],3))+6.928203230275509*vc[0]*dv[1]*f[1]*gkyl_ipow(dv[0],2));
}