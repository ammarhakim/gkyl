#include <SpitzerNuModDecl.h> 
#include <math.h> 
#include <../../Lib/gkyl_ipow.h> 

void SpitzerNuCellAvScale2xMax_P1(const double elemCharge, const double eps0, const double hBar, const double nuFrac, const double qA, const double massA, const double *m0A, const double *vtSqA, const double qB, const double massB, const double *m0B, const double *vtSqB, const double normNu, const double *Bmag, double *nu) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // eps0:       vacuum permittivity. 
  // hBar:       Planck's constant h, divided by 2pi. 
  // nuFrac:     scaling factor. 
  // qA:         charge of species A. 
  // massA:      mass of species A. 
  // m0A[3]:     number density of species A. 
  // vtSqA[3]:   squared thermal speed, sqrt(T/m), of species A. 
  // qB:         charge of species B. 
  // massB:      mass of species B. 
  // m0B[3]:     number density of species B. 
  // vtSqB[3]:   squared thermal speed, sqrt(T/m), of species B. 
  // Bmag[3]:    magnetic field magnitude. 
  // nu[3]:      collisionality. 
 
  if ((m0B[0]>0.0) && (vtSqA[0]>0.0) && (vtSqB[0]>0.0)) {
    nu[0] = (m0B[0]*normNu*nuFrac)/sqrt(0.125*gkyl_ipow(vtSqB[0],3)+0.375*vtSqA[0]*gkyl_ipow(vtSqB[0],2)+0.375*gkyl_ipow(vtSqA[0],2)*vtSqB[0]+0.125*gkyl_ipow(vtSqA[0],3)); 
  } else {
    nu[0] = 0.0;
  }
 
} 
 
void SpitzerNuCellAvScale2xMax_P2(const double elemCharge, const double eps0, const double hBar, const double nuFrac, const double qA, const double massA, const double *m0A, const double *vtSqA, const double qB, const double massB, const double *m0B, const double *vtSqB, const double normNu, const double *Bmag, double *nu) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // eps0:       vacuum permittivity. 
  // hBar:       Planck's constant h, divided by 2pi. 
  // nuFrac:     scaling factor. 
  // qA:         charge of species A. 
  // massA:      mass of species A. 
  // m0A[6]:     number density of species A. 
  // vtSqA[6]:   squared thermal speed, sqrt(T/m), of species A. 
  // qB:         charge of species B. 
  // massB:      mass of species B. 
  // m0B[6]:     number density of species B. 
  // vtSqB[6]:   squared thermal speed, sqrt(T/m), of species B. 
  // Bmag[6]:    magnetic field magnitude. 
  // nu[6]:      collisionality. 
 
  if ((m0B[0]>0.0) && (vtSqA[0]>0.0) && (vtSqB[0]>0.0)) {
    nu[0] = (m0B[0]*normNu*nuFrac)/sqrt(0.125*gkyl_ipow(vtSqB[0],3)+0.375*vtSqA[0]*gkyl_ipow(vtSqB[0],2)+0.375*gkyl_ipow(vtSqA[0],2)*vtSqB[0]+0.125*gkyl_ipow(vtSqA[0],3)); 
  } else {
    nu[0] = 0.0;
  }
 
} 
 
void SpitzerNuCellAvScale2xMax_P3(const double elemCharge, const double eps0, const double hBar, const double nuFrac, const double qA, const double massA, const double *m0A, const double *vtSqA, const double qB, const double massB, const double *m0B, const double *vtSqB, const double normNu, const double *Bmag, double *nu) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // eps0:       vacuum permittivity. 
  // hBar:       Planck's constant h, divided by 2pi. 
  // nuFrac:     scaling factor. 
  // qA:         charge of species A. 
  // massA:      mass of species A. 
  // m0A[10]:     number density of species A. 
  // vtSqA[10]:   squared thermal speed, sqrt(T/m), of species A. 
  // qB:         charge of species B. 
  // massB:      mass of species B. 
  // m0B[10]:     number density of species B. 
  // vtSqB[10]:   squared thermal speed, sqrt(T/m), of species B. 
  // Bmag[10]:    magnetic field magnitude. 
  // nu[10]:      collisionality. 
 
  if ((m0B[0]>0.0) && (vtSqA[0]>0.0) && (vtSqB[0]>0.0)) {
    nu[0] = (m0B[0]*normNu*nuFrac)/sqrt(0.125*gkyl_ipow(vtSqB[0],3)+0.375*vtSqA[0]*gkyl_ipow(vtSqB[0],2)+0.375*gkyl_ipow(vtSqA[0],2)*vtSqB[0]+0.125*gkyl_ipow(vtSqA[0],3)); 
  } else {
    nu[0] = 0.0;
  }
 
} 
 
void SpitzerNuCellAvBuild2xMax_P1(const double elemCharge, const double eps0, const double hBar, const double nuFrac, const double qA, const double massA, const double *m0A, const double *vtSqA, const double qB, const double massB, const double *m0B, const double *vtSqB, const double normNu, const double *Bmag, double *nu) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // eps0:       vacuum permittivity. 
  // hBar:       Planck's constant h, divided by 2pi. 
  // nuFrac:     scaling factor. 
  // qA:         charge of species A. 
  // massA:      mass of species A. 
  // m0A[3]:     number density of species A. 
  // vtSqA[3]:   squared thermal speed, sqrt(T/m), of species A. 
  // qB:         charge of species B. 
  // massB:      mass of species B. 
  // m0B[3]:     number density of species B. 
  // vtSqB[3]:   squared thermal speed, sqrt(T/m), of species B. 
  // Bmag[3]:    magnetic field magnitude. 
  // nu[3]:      collisionality. 
 
  double nA0 = 0.5*m0A[0]; 
  double vtSqA0 = 0.5*vtSqA[0]; 
  double nB0 = 0.5*m0B[0]; 
  double vtSqB0 = 0.5*vtSqB[0]; 
  double Bmag0 = 0.5*Bmag[0]; 
  double logLambda;
  double massAB = (massA*massB)/(massB+massA); 
  double uRel = 3.0*vtSqB0+3.0*vtSqA0; 
  double rMaxA = 1.0/sqrt((gkyl_ipow(Bmag0,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massB,2)*vtSqB0+3.0*gkyl_ipow(massB,2)*vtSqA0)+(nB0*gkyl_ipow(qB,2))/(eps0*massB*vtSqB0+3.0*eps0*massB*vtSqA0)+(0.25*nA0*gkyl_ipow(qA,2))/(eps0*massA*vtSqA0)+(0.25*gkyl_ipow(Bmag0,2)*gkyl_ipow(qA,2))/(gkyl_ipow(massA,2)*vtSqA0)); 
  double rMaxB = 1.0/sqrt((gkyl_ipow(Bmag0,2)*gkyl_ipow(qA,2))/(3.0*gkyl_ipow(massA,2)*vtSqB0+gkyl_ipow(massA,2)*vtSqA0)+(nA0*gkyl_ipow(qA,2))/(3.0*eps0*massA*vtSqB0+eps0*massA*vtSqA0)+(0.25*nB0*gkyl_ipow(qB,2))/(eps0*massB*vtSqB0)+(0.25*gkyl_ipow(Bmag0,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massB,2)*vtSqB0)); 
  double rMin = std::max((0.07957747154594767*std::abs(qA*qB))/(eps0*massAB*gkyl_ipow(uRel,2)),(0.3032653298563167*hBar)/(massAB*uRel)); 
  logLambda = 0.5*log(rMaxB/rMin)+0.5*log(rMaxA/rMin); 
  if ((m0B[0]>0.0) && (vtSqA[0]>0.0) && (vtSqB[0]>0.0)) {
    nu[0] = ((0.04232909062282731*logLambda*nB0*nuFrac*gkyl_ipow(qA,2)*gkyl_ipow(qB,2))/(massA*massB*sqrt(gkyl_ipow(vtSqB0,3)+3.0*vtSqA0*gkyl_ipow(vtSqB0,2)+3.0*gkyl_ipow(vtSqA0,2)*vtSqB0+gkyl_ipow(vtSqA0,3)))+(0.04232909062282731*logLambda*nB0*nuFrac*gkyl_ipow(qA,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massA,2)*sqrt(gkyl_ipow(vtSqB0,3)+3.0*vtSqA0*gkyl_ipow(vtSqB0,2)+3.0*gkyl_ipow(vtSqA0,2)*vtSqB0+gkyl_ipow(vtSqA0,3))))/gkyl_ipow(eps0,2); 
  } else {
    nu[0] = 0.0;
  }
 
} 
 
void SpitzerNuCellAvBuild2xMax_P2(const double elemCharge, const double eps0, const double hBar, const double nuFrac, const double qA, const double massA, const double *m0A, const double *vtSqA, const double qB, const double massB, const double *m0B, const double *vtSqB, const double normNu, const double *Bmag, double *nu) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // eps0:       vacuum permittivity. 
  // hBar:       Planck's constant h, divided by 2pi. 
  // nuFrac:     scaling factor. 
  // qA:         charge of species A. 
  // massA:      mass of species A. 
  // m0A[6]:     number density of species A. 
  // vtSqA[6]:   squared thermal speed, sqrt(T/m), of species A. 
  // qB:         charge of species B. 
  // massB:      mass of species B. 
  // m0B[6]:     number density of species B. 
  // vtSqB[6]:   squared thermal speed, sqrt(T/m), of species B. 
  // Bmag[6]:    magnetic field magnitude. 
  // nu[6]:      collisionality. 
 
  double nA0 = 0.5*m0A[0]; 
  double vtSqA0 = 0.5*vtSqA[0]; 
  double nB0 = 0.5*m0B[0]; 
  double vtSqB0 = 0.5*vtSqB[0]; 
  double Bmag0 = 0.5*Bmag[0]; 
  double logLambda;
  double massAB = (massA*massB)/(massB+massA); 
  double uRel = 3.0*vtSqB0+3.0*vtSqA0; 
  double rMaxA = 1.0/sqrt((gkyl_ipow(Bmag0,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massB,2)*vtSqB0+3.0*gkyl_ipow(massB,2)*vtSqA0)+(nB0*gkyl_ipow(qB,2))/(eps0*massB*vtSqB0+3.0*eps0*massB*vtSqA0)+(0.25*nA0*gkyl_ipow(qA,2))/(eps0*massA*vtSqA0)+(0.25*gkyl_ipow(Bmag0,2)*gkyl_ipow(qA,2))/(gkyl_ipow(massA,2)*vtSqA0)); 
  double rMaxB = 1.0/sqrt((gkyl_ipow(Bmag0,2)*gkyl_ipow(qA,2))/(3.0*gkyl_ipow(massA,2)*vtSqB0+gkyl_ipow(massA,2)*vtSqA0)+(nA0*gkyl_ipow(qA,2))/(3.0*eps0*massA*vtSqB0+eps0*massA*vtSqA0)+(0.25*nB0*gkyl_ipow(qB,2))/(eps0*massB*vtSqB0)+(0.25*gkyl_ipow(Bmag0,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massB,2)*vtSqB0)); 
  double rMin = std::max((0.07957747154594767*std::abs(qA*qB))/(eps0*massAB*gkyl_ipow(uRel,2)),(0.3032653298563167*hBar)/(massAB*uRel)); 
  logLambda = 0.5*log(rMaxB/rMin)+0.5*log(rMaxA/rMin); 
  if ((m0B[0]>0.0) && (vtSqA[0]>0.0) && (vtSqB[0]>0.0)) {
    nu[0] = ((0.04232909062282731*logLambda*nB0*nuFrac*gkyl_ipow(qA,2)*gkyl_ipow(qB,2))/(massA*massB*sqrt(gkyl_ipow(vtSqB0,3)+3.0*vtSqA0*gkyl_ipow(vtSqB0,2)+3.0*gkyl_ipow(vtSqA0,2)*vtSqB0+gkyl_ipow(vtSqA0,3)))+(0.04232909062282731*logLambda*nB0*nuFrac*gkyl_ipow(qA,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massA,2)*sqrt(gkyl_ipow(vtSqB0,3)+3.0*vtSqA0*gkyl_ipow(vtSqB0,2)+3.0*gkyl_ipow(vtSqA0,2)*vtSqB0+gkyl_ipow(vtSqA0,3))))/gkyl_ipow(eps0,2); 
  } else {
    nu[0] = 0.0;
  }
 
} 
 
void SpitzerNuCellAvBuild2xMax_P3(const double elemCharge, const double eps0, const double hBar, const double nuFrac, const double qA, const double massA, const double *m0A, const double *vtSqA, const double qB, const double massB, const double *m0B, const double *vtSqB, const double normNu, const double *Bmag, double *nu) 
{ 
  // elemCharge: elementary charge (J - eV conversion factor). 
  // eps0:       vacuum permittivity. 
  // hBar:       Planck's constant h, divided by 2pi. 
  // nuFrac:     scaling factor. 
  // qA:         charge of species A. 
  // massA:      mass of species A. 
  // m0A[10]:     number density of species A. 
  // vtSqA[10]:   squared thermal speed, sqrt(T/m), of species A. 
  // qB:         charge of species B. 
  // massB:      mass of species B. 
  // m0B[10]:     number density of species B. 
  // vtSqB[10]:   squared thermal speed, sqrt(T/m), of species B. 
  // Bmag[10]:    magnetic field magnitude. 
  // nu[10]:      collisionality. 
 
  double nA0 = 0.5*m0A[0]; 
  double vtSqA0 = 0.5*vtSqA[0]; 
  double nB0 = 0.5*m0B[0]; 
  double vtSqB0 = 0.5*vtSqB[0]; 
  double Bmag0 = 0.5*Bmag[0]; 
  double logLambda;
  double massAB = (massA*massB)/(massB+massA); 
  double uRel = 3.0*vtSqB0+3.0*vtSqA0; 
  double rMaxA = 1.0/sqrt((gkyl_ipow(Bmag0,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massB,2)*vtSqB0+3.0*gkyl_ipow(massB,2)*vtSqA0)+(nB0*gkyl_ipow(qB,2))/(eps0*massB*vtSqB0+3.0*eps0*massB*vtSqA0)+(0.25*nA0*gkyl_ipow(qA,2))/(eps0*massA*vtSqA0)+(0.25*gkyl_ipow(Bmag0,2)*gkyl_ipow(qA,2))/(gkyl_ipow(massA,2)*vtSqA0)); 
  double rMaxB = 1.0/sqrt((gkyl_ipow(Bmag0,2)*gkyl_ipow(qA,2))/(3.0*gkyl_ipow(massA,2)*vtSqB0+gkyl_ipow(massA,2)*vtSqA0)+(nA0*gkyl_ipow(qA,2))/(3.0*eps0*massA*vtSqB0+eps0*massA*vtSqA0)+(0.25*nB0*gkyl_ipow(qB,2))/(eps0*massB*vtSqB0)+(0.25*gkyl_ipow(Bmag0,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massB,2)*vtSqB0)); 
  double rMin = std::max((0.07957747154594767*std::abs(qA*qB))/(eps0*massAB*gkyl_ipow(uRel,2)),(0.3032653298563167*hBar)/(massAB*uRel)); 
  logLambda = 0.5*log(rMaxB/rMin)+0.5*log(rMaxA/rMin); 
  if ((m0B[0]>0.0) && (vtSqA[0]>0.0) && (vtSqB[0]>0.0)) {
    nu[0] = ((0.04232909062282731*logLambda*nB0*nuFrac*gkyl_ipow(qA,2)*gkyl_ipow(qB,2))/(massA*massB*sqrt(gkyl_ipow(vtSqB0,3)+3.0*vtSqA0*gkyl_ipow(vtSqB0,2)+3.0*gkyl_ipow(vtSqA0,2)*vtSqB0+gkyl_ipow(vtSqA0,3)))+(0.04232909062282731*logLambda*nB0*nuFrac*gkyl_ipow(qA,2)*gkyl_ipow(qB,2))/(gkyl_ipow(massA,2)*sqrt(gkyl_ipow(vtSqB0,3)+3.0*vtSqA0*gkyl_ipow(vtSqB0,2)+3.0*gkyl_ipow(vtSqA0,2)*vtSqB0+gkyl_ipow(vtSqA0,3))))/gkyl_ipow(eps0,2); 
  } else {
    nu[0] = 0.0;
  }
 
} 
 
