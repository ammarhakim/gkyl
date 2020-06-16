#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void VmProdCXcellAvMax3x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[6]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double xSqAv = (0.125*pow(u[8],2))/vtSqAv-(0.7071067811865475*w[5]*u[8])/vtSqAv+pow(w[5],2)/vtSqAv+pow(w[4],2)/vtSqAv-(0.7071067811865475*u[4]*w[4])/vtSqAv+(0.125*pow(u[4],2))/vtSqAv+pow(w[3],2)/vtSqAv-(0.7071067811865475*u[0]*w[3])/vtSqAv+(0.125*pow(u[0],2))/vtSqAv; 
  double vrelCX = sqrt(vtSqAv)*sqrt(xSqAv+1.273239544735163); 
 
  prodCX[0] = 0.3535533905932737*fOther[3]*m0[3]*vrelCX+0.3535533905932737*fOther[2]*m0[2]*vrelCX+0.3535533905932737*fOther[1]*m0[1]*vrelCX+0.3535533905932737*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.3535533905932737*fOther[0]*m0[1]*vrelCX+0.3535533905932737*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.3535533905932737*fOther[0]*m0[2]*vrelCX+0.3535533905932737*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.3535533905932737*fOther[0]*m0[3]*vrelCX+0.3535533905932737*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.3535533905932737*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.3535533905932737*m0[0]*fOther[5]*vrelCX; 
  prodCX[6] = 0.3535533905932737*m0[0]*fOther[6]*vrelCX; 
 
} 
