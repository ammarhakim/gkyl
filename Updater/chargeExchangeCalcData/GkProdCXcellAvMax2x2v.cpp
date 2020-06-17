#include <RelativeVelocityModDecl.h> 
#include <math.h> 
void GkProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[2]:     density. 
  // uPar[2]:    velocity. 
  // vtSq[2]:   squared thermal speed, sqrt(T/m). 
  // fOther:     distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = pow(w[2],2)/vtSqAv-(1.0*uPar[0]*w[2])/vtSqAv+(0.25*pow(uPar[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.5*fOther[2]*m0[2]*vrelCX+0.5*fOther[1]*m0[1]*vrelCX+0.5*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.5*fOther[0]*m0[1]*vrelCX+0.5*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.5*fOther[0]*m0[2]*vrelCX+0.5*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.5*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.5*m0[0]*fOther[4]*vrelCX; 
 
} 
void GkProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:      cell-center coordinates. 
  // m0[2]:     density. 
  // uPar[2]:    velocity. 
  // vtSq[2]:   squared thermal speed, sqrt(T/m). 
  // fOther:     distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = pow(w[2],2)/vtSqAv-(1.0*uPar[0]*w[2])/vtSqAv+(0.25*pow(uPar[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.5*m0[5]*fOther[12]*vrelCX+0.5*m0[4]*fOther[11]*vrelCX+0.5*m0[3]*fOther[5]*vrelCX+0.5*fOther[2]*m0[2]*vrelCX+0.5*fOther[1]*m0[1]*vrelCX+0.5*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.4472135954999579*m0[1]*fOther[11]*vrelCX+0.5*m0[2]*fOther[5]*vrelCX+0.4472135954999579*fOther[1]*m0[4]*vrelCX+0.5*fOther[2]*m0[3]*vrelCX+0.5*fOther[0]*m0[1]*vrelCX+0.5*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.4472135954999579*m0[2]*fOther[12]*vrelCX+0.4472135954999579*fOther[2]*m0[5]*vrelCX+0.5*m0[1]*fOther[5]*vrelCX+0.5*fOther[1]*m0[3]*vrelCX+0.5*fOther[0]*m0[2]*vrelCX+0.5*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.5*m0[2]*fOther[7]*vrelCX+0.5*m0[1]*fOther[6]*vrelCX+0.5*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.5*m0[2]*fOther[9]*vrelCX+0.5*m0[1]*fOther[8]*vrelCX+0.5*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.4472135954999579*m0[3]*fOther[12]*vrelCX+0.4472135954999579*m0[3]*fOther[11]*vrelCX+0.4472135954999579*fOther[5]*m0[5]*vrelCX+0.4472135954999579*m0[4]*fOther[5]*vrelCX+0.5*m0[0]*fOther[5]*vrelCX+0.5*fOther[0]*m0[3]*vrelCX+0.5*fOther[1]*m0[2]*vrelCX+0.5*m0[1]*fOther[2]*vrelCX; 
  prodCX[6] = 0.5*m0[3]*fOther[7]*vrelCX+0.4472135954999579*m0[4]*fOther[6]*vrelCX+0.5*m0[0]*fOther[6]*vrelCX+0.5*m0[1]*fOther[3]*vrelCX; 
  prodCX[7] = 0.4472135954999579*m0[5]*fOther[7]*vrelCX+0.5*m0[0]*fOther[7]*vrelCX+0.5*m0[3]*fOther[6]*vrelCX+0.5*m0[2]*fOther[3]*vrelCX; 
  prodCX[8] = 0.5*m0[3]*fOther[9]*vrelCX+0.4472135954999579*m0[4]*fOther[8]*vrelCX+0.5*m0[0]*fOther[8]*vrelCX+0.5*m0[1]*fOther[4]*vrelCX; 
  prodCX[9] = 0.4472135954999579*m0[5]*fOther[9]*vrelCX+0.5*m0[0]*fOther[9]*vrelCX+0.5*m0[3]*fOther[8]*vrelCX+0.5*m0[2]*fOther[4]*vrelCX; 
  prodCX[10] = 0.5*m0[0]*fOther[10]*vrelCX; 
  prodCX[11] = 0.31943828249997*m0[4]*fOther[11]*vrelCX+0.5*m0[0]*fOther[11]*vrelCX+0.4472135954999579*m0[3]*fOther[5]*vrelCX+0.5*fOther[0]*m0[4]*vrelCX+0.4472135954999579*fOther[1]*m0[1]*vrelCX; 
  prodCX[12] = 0.31943828249997*m0[5]*fOther[12]*vrelCX+0.5*m0[0]*fOther[12]*vrelCX+0.5*fOther[0]*m0[5]*vrelCX+0.4472135954999579*m0[3]*fOther[5]*vrelCX+0.4472135954999579*fOther[2]*m0[2]*vrelCX; 
  prodCX[13] = 0.5*m0[0]*fOther[13]*vrelCX; 
  prodCX[14] = 0.5*m0[0]*fOther[14]*vrelCX; 
 
} 
