#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void GkProdCXcellAvMax3x2v_P1(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[3]:     density. 
  // uPar[3]:    velocity. 
  // vtSq[3]:   squared thermal speed, sqrt(T/m). 
  // fOther:     distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double xSqAv = pow(w[3],2)/vtSqAv-(0.7071067811865475*uPar[0]*w[3])/vtSqAv+(0.125*pow(uPar[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.3535533905932737*fOther[3]*m0[3]*vrelCX+0.3535533905932737*fOther[2]*m0[2]*vrelCX+0.3535533905932737*fOther[1]*m0[1]*vrelCX+0.3535533905932737*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.3535533905932737*fOther[0]*m0[1]*vrelCX+0.3535533905932737*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.3535533905932737*fOther[0]*m0[2]*vrelCX+0.3535533905932737*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.3535533905932737*fOther[0]*m0[3]*vrelCX+0.3535533905932737*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.3535533905932737*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.3535533905932737*m0[0]*fOther[5]*vrelCX; 
 
} 
void GkProdCXcellAvMax3x2v_P2(const double *w, const double *m0, const double *uPar, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[5]:      cell-center coordinates. 
  // m0[3]:     density. 
  // uPar[3]:    velocity. 
  // vtSq[3]:   squared thermal speed, sqrt(T/m). 
  // fOther:     distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.3535533905932738*vtSq[0]; 
  double xSqAv = pow(w[3],2)/vtSqAv-(0.7071067811865475*uPar[0]*w[3])/vtSqAv+(0.125*pow(uPar[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.3535533905932737*m0[9]*fOther[18]*vrelCX+0.3535533905932737*m0[8]*fOther[17]*vrelCX+0.3535533905932737*m0[7]*fOther[16]*vrelCX+0.3535533905932737*m0[6]*fOther[8]*vrelCX+0.3535533905932737*m0[5]*fOther[7]*vrelCX+0.3535533905932737*m0[4]*fOther[6]*vrelCX+0.3535533905932737*fOther[3]*m0[3]*vrelCX+0.3535533905932737*fOther[2]*m0[2]*vrelCX+0.3535533905932737*fOther[1]*m0[1]*vrelCX+0.3535533905932737*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.3162277660168379*m0[1]*fOther[16]*vrelCX+0.3162277660168379*fOther[1]*m0[7]*vrelCX+0.3535533905932737*m0[3]*fOther[7]*vrelCX+0.3535533905932737*m0[2]*fOther[6]*vrelCX+0.3535533905932737*fOther[3]*m0[5]*vrelCX+0.3535533905932737*fOther[2]*m0[4]*vrelCX+0.3535533905932737*fOther[0]*m0[1]*vrelCX+0.3535533905932737*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.3162277660168379*m0[2]*fOther[17]*vrelCX+0.3162277660168379*fOther[2]*m0[8]*vrelCX+0.3535533905932737*m0[3]*fOther[8]*vrelCX+0.3535533905932737*fOther[3]*m0[6]*vrelCX+0.3535533905932737*m0[1]*fOther[6]*vrelCX+0.3535533905932737*fOther[1]*m0[4]*vrelCX+0.3535533905932737*fOther[0]*m0[2]*vrelCX+0.3535533905932737*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.3162277660168379*m0[3]*fOther[18]*vrelCX+0.3162277660168379*fOther[3]*m0[9]*vrelCX+0.3535533905932737*m0[2]*fOther[8]*vrelCX+0.3535533905932737*m0[1]*fOther[7]*vrelCX+0.3535533905932737*fOther[2]*m0[6]*vrelCX+0.3535533905932737*fOther[1]*m0[5]*vrelCX+0.3535533905932737*fOther[0]*m0[3]*vrelCX+0.3535533905932737*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.3535533905932737*m0[3]*fOther[11]*vrelCX+0.3535533905932737*m0[2]*fOther[10]*vrelCX+0.3535533905932737*m0[1]*fOther[9]*vrelCX+0.3535533905932737*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.3535533905932737*m0[3]*fOther[14]*vrelCX+0.3535533905932737*m0[2]*fOther[13]*vrelCX+0.3535533905932737*m0[1]*fOther[12]*vrelCX+0.3535533905932737*m0[0]*fOther[5]*vrelCX; 
  prodCX[6] = 0.3162277660168379*m0[4]*fOther[17]*vrelCX+0.3162277660168379*m0[4]*fOther[16]*vrelCX+0.3162277660168379*fOther[6]*m0[8]*vrelCX+0.3535533905932737*m0[5]*fOther[8]*vrelCX+0.3162277660168379*fOther[6]*m0[7]*vrelCX+0.3535533905932737*m0[6]*fOther[7]*vrelCX+0.3535533905932737*m0[0]*fOther[6]*vrelCX+0.3535533905932737*fOther[0]*m0[4]*vrelCX+0.3535533905932737*fOther[1]*m0[2]*vrelCX+0.3535533905932737*m0[1]*fOther[2]*vrelCX; 
  prodCX[7] = 0.3162277660168379*m0[5]*fOther[18]*vrelCX+0.3162277660168379*m0[5]*fOther[16]*vrelCX+0.3162277660168379*fOther[7]*m0[9]*vrelCX+0.3535533905932737*m0[4]*fOther[8]*vrelCX+0.3162277660168379*fOther[7]*m0[7]*vrelCX+0.3535533905932737*m0[0]*fOther[7]*vrelCX+0.3535533905932737*fOther[6]*m0[6]*vrelCX+0.3535533905932737*fOther[0]*m0[5]*vrelCX+0.3535533905932737*fOther[1]*m0[3]*vrelCX+0.3535533905932737*m0[1]*fOther[3]*vrelCX; 
  prodCX[8] = 0.3162277660168379*m0[6]*fOther[18]*vrelCX+0.3162277660168379*m0[6]*fOther[17]*vrelCX+0.3162277660168379*fOther[8]*m0[9]*vrelCX+0.3162277660168379*fOther[8]*m0[8]*vrelCX+0.3535533905932737*m0[0]*fOther[8]*vrelCX+0.3535533905932737*m0[4]*fOther[7]*vrelCX+0.3535533905932737*fOther[0]*m0[6]*vrelCX+0.3535533905932737*m0[5]*fOther[6]*vrelCX+0.3535533905932737*fOther[2]*m0[3]*vrelCX+0.3535533905932737*m0[2]*fOther[3]*vrelCX; 
  prodCX[9] = 0.3535533905932737*m0[5]*fOther[11]*vrelCX+0.3535533905932737*m0[4]*fOther[10]*vrelCX+0.3162277660168379*m0[7]*fOther[9]*vrelCX+0.3535533905932737*m0[0]*fOther[9]*vrelCX+0.3535533905932737*m0[1]*fOther[4]*vrelCX; 
  prodCX[10] = 0.3535533905932737*m0[6]*fOther[11]*vrelCX+0.3162277660168379*m0[8]*fOther[10]*vrelCX+0.3535533905932737*m0[0]*fOther[10]*vrelCX+0.3535533905932737*m0[4]*fOther[9]*vrelCX+0.3535533905932737*m0[2]*fOther[4]*vrelCX; 
  prodCX[11] = 0.3162277660168379*m0[9]*fOther[11]*vrelCX+0.3535533905932737*m0[0]*fOther[11]*vrelCX+0.3535533905932737*m0[6]*fOther[10]*vrelCX+0.3535533905932737*m0[5]*fOther[9]*vrelCX+0.3535533905932737*m0[3]*fOther[4]*vrelCX; 
  prodCX[12] = 0.3535533905932737*m0[5]*fOther[14]*vrelCX+0.3535533905932737*m0[4]*fOther[13]*vrelCX+0.3162277660168379*m0[7]*fOther[12]*vrelCX+0.3535533905932737*m0[0]*fOther[12]*vrelCX+0.3535533905932737*m0[1]*fOther[5]*vrelCX; 
  prodCX[13] = 0.3535533905932737*m0[6]*fOther[14]*vrelCX+0.3162277660168379*m0[8]*fOther[13]*vrelCX+0.3535533905932737*m0[0]*fOther[13]*vrelCX+0.3535533905932737*m0[4]*fOther[12]*vrelCX+0.3535533905932737*m0[2]*fOther[5]*vrelCX; 
  prodCX[14] = 0.3162277660168379*m0[9]*fOther[14]*vrelCX+0.3535533905932737*m0[0]*fOther[14]*vrelCX+0.3535533905932737*m0[6]*fOther[13]*vrelCX+0.3535533905932737*m0[5]*fOther[12]*vrelCX+0.3535533905932737*m0[3]*fOther[5]*vrelCX; 
  prodCX[15] = 0.3535533905932737*m0[0]*fOther[15]*vrelCX; 
  prodCX[16] = 0.2258769757263128*m0[7]*fOther[16]*vrelCX+0.3535533905932737*m0[0]*fOther[16]*vrelCX+0.3535533905932737*fOther[0]*m0[7]*vrelCX+0.3162277660168379*m0[5]*fOther[7]*vrelCX+0.3162277660168379*m0[4]*fOther[6]*vrelCX+0.3162277660168379*fOther[1]*m0[1]*vrelCX; 
  prodCX[17] = 0.2258769757263128*m0[8]*fOther[17]*vrelCX+0.3535533905932737*m0[0]*fOther[17]*vrelCX+0.3535533905932737*fOther[0]*m0[8]*vrelCX+0.3162277660168379*m0[6]*fOther[8]*vrelCX+0.3162277660168379*m0[4]*fOther[6]*vrelCX+0.3162277660168379*fOther[2]*m0[2]*vrelCX; 
  prodCX[18] = 0.2258769757263128*m0[9]*fOther[18]*vrelCX+0.3535533905932737*m0[0]*fOther[18]*vrelCX+0.3535533905932737*fOther[0]*m0[9]*vrelCX+0.3162277660168379*m0[6]*fOther[8]*vrelCX+0.3162277660168379*m0[5]*fOther[7]*vrelCX+0.3162277660168379*fOther[3]*m0[3]*vrelCX; 
  prodCX[19] = 0.3535533905932737*m0[0]*fOther[19]*vrelCX; 
  prodCX[20] = 0.3535533905932737*m0[0]*fOther[20]*vrelCX; 
 
} 
