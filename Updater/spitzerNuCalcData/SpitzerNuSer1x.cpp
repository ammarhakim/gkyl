#include <SpitzerNuModDecl.h> 
#include <math.h> 
#include <../../Lib/gkyl_ipow.h> 

void SpitzerNuCellAv1xSer_P1(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[1]:   number density. 
  // vtSq[1]: squared thermal speed, sqrt(T/m). 
  // nu[1]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.3535533905932737*gkyl_ipow(vtSq[0],3)); 
 
} 
 
void SpitzerNuCellAv1xSer_P2(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[1]:   number density. 
  // vtSq[1]: squared thermal speed, sqrt(T/m). 
  // nu[1]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.3535533905932737*gkyl_ipow(vtSq[0],3)); 
 
} 
 
void SpitzerNuCellAv1xSer_P3(const double normNu, const double rmR3d2, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[1]:   number density. 
  // vtSq[1]: squared thermal speed, sqrt(T/m). 
  // nu[1]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.3535533905932737*gkyl_ipow(vtSq[0],3)); 
 
} 
 
