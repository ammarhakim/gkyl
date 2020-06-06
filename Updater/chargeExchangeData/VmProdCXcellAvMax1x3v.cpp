#include <ChargeExchangeModDecl.h> 
#include <math.h> 
void VmProdCXcellAvMax1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = (0.5*pow(u[4],2))/vtSqAv-(1.414213562373095*w[3]*u[4])/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv-(1.414213562373095*u[2]*w[2])/vtSqAv+(0.5*pow(u[2],2))/vtSqAv+pow(w[1],2)/vtSqAv-(1.414213562373095*u[0]*w[1])/vtSqAv+(0.5*pow(u[0],2))/vtSqAv; 
  double vrelCX = sqrt(vtSqAv)*sqrt(xSqAv+1.273239544735163); 
 
  prodCX[0] = 0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.7071067811865475*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.7071067811865475*m0[0]*fOther[4]*vrelCX; 
 
} 
void VmProdCXcellAvMax1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = (0.5*pow(u[6],2))/vtSqAv-(1.414213562373095*w[3]*u[6])/vtSqAv+pow(w[3],2)/vtSqAv+(0.5*pow(u[3],2))/vtSqAv-(1.414213562373095*w[2]*u[3])/vtSqAv+pow(w[2],2)/vtSqAv+pow(w[1],2)/vtSqAv-(1.414213562373095*u[0]*w[1])/vtSqAv+(0.5*pow(u[0],2))/vtSqAv; 
  double vrelCX = sqrt(vtSqAv)*sqrt(xSqAv+1.273239544735163); 
 
  prodCX[0] = 0.7071067811865475*m0[2]*fOther[11]*vrelCX+0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.6324555320336759*m0[1]*fOther[11]*vrelCX+0.6324555320336759*fOther[1]*m0[2]*vrelCX+0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865475*m0[1]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.7071067811865475*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.7071067811865475*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.6324555320336759*m0[2]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[5]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
  prodCX[6] = 0.6324555320336759*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[6]*vrelCX+0.7071067811865475*m0[1]*fOther[3]*vrelCX; 
  prodCX[7] = 0.7071067811865475*m0[0]*fOther[7]*vrelCX; 
  prodCX[8] = 0.6324555320336759*m0[2]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[8]*vrelCX+0.7071067811865475*m0[1]*fOther[4]*vrelCX; 
  prodCX[9] = 0.7071067811865475*m0[0]*fOther[9]*vrelCX; 
  prodCX[10] = 0.7071067811865475*m0[0]*fOther[10]*vrelCX; 
  prodCX[11] = 0.4517539514526256*m0[2]*fOther[11]*vrelCX+0.7071067811865475*m0[0]*fOther[11]*vrelCX+0.7071067811865475*fOther[0]*m0[2]*vrelCX+0.6324555320336759*fOther[1]*m0[1]*vrelCX; 
  prodCX[12] = 0.7071067811865475*m0[0]*fOther[12]*vrelCX; 
  prodCX[13] = 0.7071067811865475*m0[0]*fOther[13]*vrelCX; 
  prodCX[14] = 0.7071067811865475*m0[0]*fOther[14]*vrelCX; 
 
} 
void VmProdCXcellAvMax1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = (0.5*pow(u[8],2))/vtSqAv-(1.414213562373095*w[3]*u[8])/vtSqAv+(0.5*pow(u[4],2))/vtSqAv-(1.414213562373095*w[2]*u[4])/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv+pow(w[1],2)/vtSqAv-(1.414213562373095*u[0]*w[1])/vtSqAv+(0.5*pow(u[0],2))/vtSqAv; 
  double vrelCX = sqrt(vtSqAv)*sqrt(xSqAv+1.273239544735163); 
 
  prodCX[0] = 0.7071067811865475*m0[3]*fOther[31]*vrelCX+0.7071067811865475*m0[2]*fOther[11]*vrelCX+0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.6210590034081186*m0[2]*fOther[31]*vrelCX+0.6210590034081186*m0[3]*fOther[11]*vrelCX+0.6324555320336759*m0[1]*fOther[11]*vrelCX+0.6324555320336759*fOther[1]*m0[2]*vrelCX+0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865475*m0[2]*fOther[19]*vrelCX+0.7071067811865475*m0[1]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.7071067811865475*m0[2]*fOther[21]*vrelCX+0.7071067811865475*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.7071067811865475*m0[2]*fOther[25]*vrelCX+0.7071067811865475*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.6210590034081187*m0[3]*fOther[19]*vrelCX+0.632455532033676*m0[1]*fOther[19]*vrelCX+0.6324555320336759*m0[2]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[5]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
  prodCX[6] = 0.6210590034081187*m0[3]*fOther[21]*vrelCX+0.632455532033676*m0[1]*fOther[21]*vrelCX+0.6324555320336759*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[6]*vrelCX+0.7071067811865475*m0[1]*fOther[3]*vrelCX; 
  prodCX[7] = 0.7071067811865475*m0[1]*fOther[15]*vrelCX+0.7071067811865475*m0[0]*fOther[7]*vrelCX; 
  prodCX[8] = 0.6210590034081187*m0[3]*fOther[25]*vrelCX+0.632455532033676*m0[1]*fOther[25]*vrelCX+0.6324555320336759*m0[2]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[8]*vrelCX+0.7071067811865475*m0[1]*fOther[4]*vrelCX; 
  prodCX[9] = 0.7071067811865475*m0[1]*fOther[16]*vrelCX+0.7071067811865475*m0[0]*fOther[9]*vrelCX; 
  prodCX[10] = 0.7071067811865475*m0[1]*fOther[17]*vrelCX+0.7071067811865475*m0[0]*fOther[10]*vrelCX; 
  prodCX[11] = 0.421637021355784*m0[3]*fOther[31]*vrelCX+0.6210590034081186*m0[1]*fOther[31]*vrelCX+0.4517539514526256*m0[2]*fOther[11]*vrelCX+0.7071067811865475*m0[0]*fOther[11]*vrelCX+0.6210590034081186*fOther[1]*m0[3]*vrelCX+0.7071067811865475*fOther[0]*m0[2]*vrelCX+0.6324555320336759*fOther[1]*m0[1]*vrelCX; 
  prodCX[12] = 0.7071067811865475*m0[1]*fOther[20]*vrelCX+0.7071067811865475*m0[0]*fOther[12]*vrelCX; 
  prodCX[13] = 0.7071067811865475*m0[1]*fOther[23]*vrelCX+0.7071067811865475*m0[0]*fOther[13]*vrelCX; 
  prodCX[14] = 0.7071067811865475*m0[1]*fOther[28]*vrelCX+0.7071067811865475*m0[0]*fOther[14]*vrelCX; 
  prodCX[15] = 0.6324555320336759*m0[2]*fOther[15]*vrelCX+0.7071067811865475*m0[0]*fOther[15]*vrelCX+0.7071067811865475*m0[1]*fOther[7]*vrelCX; 
  prodCX[16] = 0.6324555320336759*m0[2]*fOther[16]*vrelCX+0.7071067811865475*m0[0]*fOther[16]*vrelCX+0.7071067811865475*m0[1]*fOther[9]*vrelCX; 
  prodCX[17] = 0.6324555320336759*m0[2]*fOther[17]*vrelCX+0.7071067811865475*m0[0]*fOther[17]*vrelCX+0.7071067811865475*m0[1]*fOther[10]*vrelCX; 
  prodCX[18] = 0.7071067811865475*m0[0]*fOther[18]*vrelCX; 
  prodCX[19] = 0.4517539514526256*m0[2]*fOther[19]*vrelCX+0.7071067811865475*m0[0]*fOther[19]*vrelCX+0.6210590034081187*m0[3]*fOther[5]*vrelCX+0.632455532033676*m0[1]*fOther[5]*vrelCX+0.7071067811865475*fOther[2]*m0[2]*vrelCX; 
  prodCX[20] = 0.6324555320336759*m0[2]*fOther[20]*vrelCX+0.7071067811865475*m0[0]*fOther[20]*vrelCX+0.7071067811865475*m0[1]*fOther[12]*vrelCX; 
  prodCX[21] = 0.4517539514526256*m0[2]*fOther[21]*vrelCX+0.7071067811865475*m0[0]*fOther[21]*vrelCX+0.6210590034081187*m0[3]*fOther[6]*vrelCX+0.632455532033676*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[2]*fOther[3]*vrelCX; 
  prodCX[22] = 0.7071067811865475*m0[0]*fOther[22]*vrelCX; 
  prodCX[23] = 0.6324555320336759*m0[2]*fOther[23]*vrelCX+0.7071067811865475*m0[0]*fOther[23]*vrelCX+0.7071067811865475*m0[1]*fOther[13]*vrelCX; 
  prodCX[24] = 0.7071067811865475*m0[0]*fOther[24]*vrelCX; 
  prodCX[25] = 0.4517539514526256*m0[2]*fOther[25]*vrelCX+0.7071067811865475*m0[0]*fOther[25]*vrelCX+0.6210590034081187*m0[3]*fOther[8]*vrelCX+0.632455532033676*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[2]*fOther[4]*vrelCX; 
  prodCX[26] = 0.7071067811865475*m0[0]*fOther[26]*vrelCX; 
  prodCX[27] = 0.7071067811865475*m0[0]*fOther[27]*vrelCX; 
  prodCX[28] = 0.6324555320336759*m0[2]*fOther[28]*vrelCX+0.7071067811865475*m0[0]*fOther[28]*vrelCX+0.7071067811865475*m0[1]*fOther[14]*vrelCX; 
  prodCX[29] = 0.7071067811865475*m0[0]*fOther[29]*vrelCX; 
  prodCX[30] = 0.7071067811865475*m0[0]*fOther[30]*vrelCX; 
  prodCX[31] = 0.421637021355784*m0[2]*fOther[31]*vrelCX+0.7071067811865475*m0[0]*fOther[31]*vrelCX+0.421637021355784*m0[3]*fOther[11]*vrelCX+0.6210590034081186*m0[1]*fOther[11]*vrelCX+0.7071067811865475*fOther[0]*m0[3]*vrelCX+0.6210590034081186*fOther[1]*m0[2]*vrelCX; 
  prodCX[32] = 0.7071067811865475*m0[0]*fOther[32]*vrelCX; 
  prodCX[33] = 0.7071067811865475*m0[0]*fOther[33]*vrelCX; 
  prodCX[34] = 0.7071067811865475*m0[0]*fOther[34]*vrelCX; 
 
} 
