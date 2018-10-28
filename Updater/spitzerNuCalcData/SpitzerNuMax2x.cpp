#include <SpitzerNuModDecl.h> 
#include <math.h> 
#include <../../Lib/gkyl_ipow.h> 

void SpitzerNuCellAv2xMax_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[2]:   number density. 
  // vtSq[2]: squared thermal speed, sqrt(T/m). 
  // nu[2]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.125*gkyl_ipow(vtSq[0],3)); 
 
} 
 
void SpitzerNuCellAv2xMax_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[2]:   number density. 
  // vtSq[2]: squared thermal speed, sqrt(T/m). 
  // nu[2]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.125*gkyl_ipow(vtSq[0],3)); 
 
} 
 
void SpitzerNuCellAv2xMax_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[2]:   number density. 
  // vtSq[2]: squared thermal speed, sqrt(T/m). 
  // nu[2]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.125*gkyl_ipow(vtSq[0],3)); 
 
} 
 
