#include <math.h> 
#include <lagrangeFixDecl.h> 
#include <gkyl_ipow.h> 

void lagrangeFixSer1x1v1p(double *dm0, double *dm1, double *dm2, double *L, double *Nv, double *vc, double *f) {
  f[0] = f[0] +  -1.0*((84.85281374238573*dm0[0]*gkyl_ipow(L[0],2)-1018.233764908629*dm2[0])*gkyl_ipow(Nv[0],4)*gkyl_ipow(vc[0],2)+vc[0]*(67.8822509939086*dm1[0]*gkyl_ipow(L[0],2)-67.8822509939086*dm1[0]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],4))+(84.85281374238573*dm2[0]*gkyl_ipow(L[0],2)-12.72792206135786*dm0[0]*gkyl_ipow(L[0],4))*gkyl_ipow(Nv[0],4)+(7.071067811865476*dm0[0]*gkyl_ipow(L[0],4)-84.85281374238573*dm2[0]*gkyl_ipow(L[0],2))*gkyl_ipow(Nv[0],2)+5.656854249492382*dm0[0]*gkyl_ipow(L[0],4))*gkyl_ipow(4.0*gkyl_ipow(L[0],5)*gkyl_ipow(Nv[0],4)-4.0*gkyl_ipow(L[0],5),-1);
  f[1] = f[1] +  -1.0*((84.85281374238573*dm0[1]*gkyl_ipow(L[0],2)-1018.233764908629*dm2[1])*gkyl_ipow(Nv[0],4)*gkyl_ipow(vc[0],2)+vc[0]*(67.8822509939086*dm1[1]*gkyl_ipow(L[0],2)-67.8822509939086*dm1[1]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],4))+(84.85281374238573*dm2[1]*gkyl_ipow(L[0],2)-12.72792206135786*dm0[1]*gkyl_ipow(L[0],4))*gkyl_ipow(Nv[0],4)+(7.071067811865476*dm0[1]*gkyl_ipow(L[0],4)-84.85281374238573*dm2[1]*gkyl_ipow(L[0],2))*gkyl_ipow(Nv[0],2)+5.656854249492382*dm0[1]*gkyl_ipow(L[0],4))*gkyl_ipow(4.0*gkyl_ipow(L[0],5)*gkyl_ipow(Nv[0],4)-4.0*gkyl_ipow(L[0],5),-1);
  f[2] = f[2] +  -1.0*(vc[0]*(30.0*dm0[0]*gkyl_ipow(L[0],2)-360.0*dm2[0])*gkyl_ipow(Nv[0],4)-12.0*dm1[0]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],4)+12.0*dm1[0]*gkyl_ipow(L[0],2))*gkyl_ipow(2.449489742783178*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],5)-2.449489742783178*Nv[0]*gkyl_ipow(L[0],4),-1);
  f[3] = f[3] +  -1.0*(vc[0]*(30.0*dm0[1]*gkyl_ipow(L[0],2)-360.0*dm2[1])*gkyl_ipow(Nv[0],4)-12.0*dm1[1]*gkyl_ipow(L[0],2)*gkyl_ipow(Nv[0],4)+12.0*dm1[1]*gkyl_ipow(L[0],2))*gkyl_ipow(2.449489742783178*gkyl_ipow(L[0],4)*gkyl_ipow(Nv[0],5)-2.449489742783178*Nv[0]*gkyl_ipow(L[0],4),-1);
}
