#include <SpitzerNuModDecl.h> 
#include <math.h> 
#include <../../Lib/gkyl_ipow.h> 

void SpitzerNuCellAvScale3xMax_P1(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[4]:   number density. 
  // vtSq[4]: squared thermal speed, sqrt(T/m). 
  // nu[4]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.0441941738241592*gkyl_ipow(vtSq[0],3)); 
 
} 
 
void SpitzerNuCellAvScale3xMax_P2(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[10]:   number density. 
  // vtSq[10]: squared thermal speed, sqrt(T/m). 
  // nu[10]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.0441941738241592*gkyl_ipow(vtSq[0],3)); 
 
} 
 
void SpitzerNuCellAvScale3xMax_P3(const double normNu, const double rmR3d2, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu) 
{ 
  // nuNorm:  collisionality normalized by (T_0^(3/2)/n_0). 
  // rmR3d2:  reciprocal of mass raised to the (3/2) power. 
  // m0[20]:   number density. 
  // vtSq[20]: squared thermal speed, sqrt(T/m). 
  // nu[20]:   collisionality. 
 
  nu[0] = (m0[0]*normNu*rmR3d2)/sqrt(0.0441941738241592*gkyl_ipow(vtSq[0],3)); 
 
} 
 
void SpitzerNuCellAvBuild3xMax_P1(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu) 
{ 
  // elemCharge: elementary charge. 
  // mass:       mass of this species. 
  // massSq:     mass of this species squared. 
  // chargeR4:   charge of this species raised to the power of 4. 
  // eps0Sq:     vacuum permittivity squared. 
  // m0[4]:      number density. 
  // vtSq[4]:    squared thermal speed, sqrt(T/m). 
  // nu[4]:      collisionality. 
 
  double n0 = 0.3535533905932738*m0[0]; 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  double T0 = (0.3535533905932738*vtSq[0]*mass)/elemCharge; 
  double logLambda;
  if (T0 < 50.0) {
    logLambda = (-1.15*log10(9.999999999999999e-7*n0))+3.45*log10(T0)+23.4; 
  } else {
    logLambda = (-1.15*log10(9.999999999999999e-7*n0))+2.3*log10(T0)+25.3; 
  }
  nu[0] = (0.05986237404172218*chargeR4*logLambda*n0)/(eps0Sq*massSq*sqrt(gkyl_ipow(vtSq0,3))); 
 
} 
 
void SpitzerNuCellAvBuild3xMax_P2(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu) 
{ 
  // elemCharge: elementary charge. 
  // mass:       mass of this species. 
  // massSq:     mass of this species squared. 
  // chargeR4:   charge of this species raised to the power of 4. 
  // eps0Sq:     vacuum permittivity squared. 
  // m0[10]:      number density. 
  // vtSq[10]:    squared thermal speed, sqrt(T/m). 
  // nu[10]:      collisionality. 
 
  double n0 = 0.3535533905932738*m0[0]; 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  double T0 = (0.3535533905932738*vtSq[0]*mass)/elemCharge; 
  double logLambda;
  if (T0 < 50.0) {
    logLambda = (-1.15*log10(9.999999999999999e-7*n0))+3.45*log10(T0)+23.4; 
  } else {
    logLambda = (-1.15*log10(9.999999999999999e-7*n0))+2.3*log10(T0)+25.3; 
  }
  nu[0] = (0.05986237404172218*chargeR4*logLambda*n0)/(eps0Sq*massSq*sqrt(gkyl_ipow(vtSq0,3))); 
 
} 
 
void SpitzerNuCellAvBuild3xMax_P3(const double elemCharge, const double mass, const double massSq, const double chargeR4, const double eps0Sq, const double *m0, const double *vtSq, double *nu) 
{ 
  // elemCharge: elementary charge. 
  // mass:       mass of this species. 
  // massSq:     mass of this species squared. 
  // chargeR4:   charge of this species raised to the power of 4. 
  // eps0Sq:     vacuum permittivity squared. 
  // m0[20]:      number density. 
  // vtSq[20]:    squared thermal speed, sqrt(T/m). 
  // nu[20]:      collisionality. 
 
  double n0 = 0.3535533905932738*m0[0]; 
  double vtSq0 = 0.3535533905932738*vtSq[0]; 
  double T0 = (0.3535533905932738*vtSq[0]*mass)/elemCharge; 
  double logLambda;
  if (T0 < 50.0) {
    logLambda = (-1.15*log10(9.999999999999999e-7*n0))+3.45*log10(T0)+23.4; 
  } else {
    logLambda = (-1.15*log10(9.999999999999999e-7*n0))+2.3*log10(T0)+25.3; 
  }
  nu[0] = (0.05986237404172218*chargeR4*logLambda*n0)/(eps0Sq*massSq*sqrt(gkyl_ipow(vtSq0,3))); 
 
} 
 
