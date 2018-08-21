#include <math.h> 
#include <lagrangeFixDecl.h> 

void lagrangeFixSer1x1v1p(const double *dm0, const double *dm1, const double*dm2, double *f, const double vc, const double L, const double Nv) {
  const double *dm0s = dm0 - 1;
  const double *dm1s = dm1 - 1;
  const double *dm2s = dm2 - 1;
  f[0] = f[0] +  -(1.0*((84.85281374238573*dm0s[1]*L*L-1018.233764908629*dm2s[1])*Nv*Nv*Nv*Nv*(vc*vc)+(67.8822509939086*dm1s[1]*L*L-67.8822509939086*dm1s[1]*L*L*Nv*Nv*Nv*Nv)*vc+(84.85281374238573*dm2s[1]*L*L-12.72792206135786*dm0s[1]*L*L*L*L)*Nv*Nv*Nv*Nv+(7.071067811865476*dm0s[1]*L*L*L*L-84.85281374238573*dm2s[1]*L*L)*Nv*Nv+5.656854249492382*dm0s[1]*L*L*L*L))/(4.0*(L*L*L*L*L)*Nv*Nv*Nv*Nv-4.0*(L*L*L*L*L));
  f[1] = f[1] +  -(1.0*((84.85281374238573*dm0s[2]*L*L-1018.233764908629*dm2s[2])*Nv*Nv*Nv*Nv*(vc*vc)+(67.8822509939086*dm1s[2]*L*L-67.8822509939086*dm1s[2]*L*L*Nv*Nv*Nv*Nv)*vc+(84.85281374238573*dm2s[2]*L*L-12.72792206135786*dm0s[2]*L*L*L*L)*Nv*Nv*Nv*Nv+(7.071067811865476*dm0s[2]*L*L*L*L-84.85281374238573*dm2s[2]*L*L)*Nv*Nv+5.656854249492382*dm0s[2]*L*L*L*L))/(4.0*(L*L*L*L*L)*Nv*Nv*Nv*Nv-4.0*(L*L*L*L*L));
  f[2] = f[2] +  -(1.0*((30.0*dm0s[1]*L*L-360.0*dm2s[1])*Nv*Nv*Nv*Nv*vc-12.0*dm1s[1]*L*L*Nv*Nv*Nv*Nv+12.0*dm1s[1]*L*L))/(2.449489742783178*L*L*L*L*Nv*Nv*Nv*Nv*Nv-2.449489742783178*L*L*L*L*Nv);
  f[3] = f[3] +  -(1.0*((30.0*dm0s[2]*L*L-360.0*dm2s[2])*Nv*Nv*Nv*Nv*vc-12.0*dm1s[2]*L*L*Nv*Nv*Nv*Nv+12.0*dm1s[2]*L*L))/(2.449489742783178*L*L*L*L*Nv*Nv*Nv*Nv*Nv-2.449489742783178*L*L*L*L*Nv);
}
