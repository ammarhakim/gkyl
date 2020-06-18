#include <RelativeVelocityModDecl.h> 
#include <math.h> 
void VmProdCXcellAvSer1x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = (0.5*pow(u[4],2))/vtSqAv-(1.414213562373095*w[3]*u[4])/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv-(1.414213562373095*u[2]*w[2])/vtSqAv+(0.5*pow(u[2],2))/vtSqAv+pow(w[1],2)/vtSqAv-(1.414213562373095*u[0]*w[1])/vtSqAv+(0.5*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865475*m0[1]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.7071067811865475*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.7071067811865475*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.7071067811865475*m0[0]*fOther[5]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
  prodCX[6] = 0.7071067811865475*m0[0]*fOther[6]*vrelCX+0.7071067811865475*m0[1]*fOther[3]*vrelCX; 
  prodCX[7] = 0.7071067811865475*m0[1]*fOther[11]*vrelCX+0.7071067811865475*m0[0]*fOther[7]*vrelCX; 
  prodCX[8] = 0.7071067811865475*m0[0]*fOther[8]*vrelCX+0.7071067811865475*m0[1]*fOther[4]*vrelCX; 
  prodCX[9] = 0.7071067811865475*m0[1]*fOther[12]*vrelCX+0.7071067811865475*m0[0]*fOther[9]*vrelCX; 
  prodCX[10] = 0.7071067811865475*m0[1]*fOther[13]*vrelCX+0.7071067811865475*m0[0]*fOther[10]*vrelCX; 
  prodCX[11] = 0.7071067811865475*m0[0]*fOther[11]*vrelCX+0.7071067811865475*m0[1]*fOther[7]*vrelCX; 
  prodCX[12] = 0.7071067811865475*m0[0]*fOther[12]*vrelCX+0.7071067811865475*m0[1]*fOther[9]*vrelCX; 
  prodCX[13] = 0.7071067811865475*m0[0]*fOther[13]*vrelCX+0.7071067811865475*m0[1]*fOther[10]*vrelCX; 
  prodCX[14] = 0.7071067811865475*m0[1]*fOther[15]*vrelCX+0.7071067811865475*m0[0]*fOther[14]*vrelCX; 
  prodCX[15] = 0.7071067811865475*m0[0]*fOther[15]*vrelCX+0.7071067811865475*m0[1]*fOther[14]*vrelCX; 
 
} 
void VmProdCXcellAvSer1x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = (0.5*pow(u[6],2))/vtSqAv-(1.414213562373095*w[3]*u[6])/vtSqAv+pow(w[3],2)/vtSqAv+(0.5*pow(u[3],2))/vtSqAv-(1.414213562373095*w[2]*u[3])/vtSqAv+pow(w[2],2)/vtSqAv+pow(w[1],2)/vtSqAv-(1.414213562373095*u[0]*w[1])/vtSqAv+(0.5*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.7071067811865475*m0[2]*fOther[11]*vrelCX+0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.6324555320336759*m0[1]*fOther[11]*vrelCX+0.6324555320336759*fOther[1]*m0[2]*vrelCX+0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865475*m0[2]*fOther[19]*vrelCX+0.7071067811865475*m0[1]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.7071067811865475*m0[2]*fOther[21]*vrelCX+0.7071067811865475*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.7071067811865475*m0[2]*fOther[25]*vrelCX+0.7071067811865475*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.632455532033676*m0[1]*fOther[19]*vrelCX+0.6324555320336759*m0[2]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[5]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
  prodCX[6] = 0.632455532033676*m0[1]*fOther[21]*vrelCX+0.6324555320336759*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[6]*vrelCX+0.7071067811865475*m0[1]*fOther[3]*vrelCX; 
  prodCX[7] = 0.7071067811865475*m0[2]*fOther[32]*vrelCX+0.7071067811865475*m0[1]*fOther[15]*vrelCX+0.7071067811865475*m0[0]*fOther[7]*vrelCX; 
  prodCX[8] = 0.632455532033676*m0[1]*fOther[25]*vrelCX+0.6324555320336759*m0[2]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[8]*vrelCX+0.7071067811865475*m0[1]*fOther[4]*vrelCX; 
  prodCX[9] = 0.7071067811865475*m0[2]*fOther[35]*vrelCX+0.7071067811865475*m0[1]*fOther[16]*vrelCX+0.7071067811865475*m0[0]*fOther[9]*vrelCX; 
  prodCX[10] = 0.7071067811865475*m0[2]*fOther[37]*vrelCX+0.7071067811865475*m0[1]*fOther[17]*vrelCX+0.7071067811865475*m0[0]*fOther[10]*vrelCX; 
  prodCX[11] = 0.4517539514526256*m0[2]*fOther[11]*vrelCX+0.7071067811865475*m0[0]*fOther[11]*vrelCX+0.7071067811865475*fOther[0]*m0[2]*vrelCX+0.6324555320336759*fOther[1]*m0[1]*vrelCX; 
  prodCX[12] = 0.7071067811865475*m0[1]*fOther[20]*vrelCX+0.7071067811865475*m0[0]*fOther[12]*vrelCX; 
  prodCX[13] = 0.7071067811865475*m0[1]*fOther[23]*vrelCX+0.7071067811865475*m0[0]*fOther[13]*vrelCX; 
  prodCX[14] = 0.7071067811865475*m0[1]*fOther[28]*vrelCX+0.7071067811865475*m0[0]*fOther[14]*vrelCX; 
  prodCX[15] = 0.6324555320336759*m0[1]*fOther[32]*vrelCX+0.6324555320336759*m0[2]*fOther[15]*vrelCX+0.7071067811865475*m0[0]*fOther[15]*vrelCX+0.7071067811865475*m0[1]*fOther[7]*vrelCX; 
  prodCX[16] = 0.6324555320336759*m0[1]*fOther[35]*vrelCX+0.6324555320336759*m0[2]*fOther[16]*vrelCX+0.7071067811865475*m0[0]*fOther[16]*vrelCX+0.7071067811865475*m0[1]*fOther[9]*vrelCX; 
  prodCX[17] = 0.6324555320336759*m0[1]*fOther[37]*vrelCX+0.6324555320336759*m0[2]*fOther[17]*vrelCX+0.7071067811865475*m0[0]*fOther[17]*vrelCX+0.7071067811865475*m0[1]*fOther[10]*vrelCX; 
  prodCX[18] = 0.7071067811865475*m0[2]*fOther[44]*vrelCX+0.7071067811865475*m0[1]*fOther[31]*vrelCX+0.7071067811865475*m0[0]*fOther[18]*vrelCX; 
  prodCX[19] = 0.4517539514526256*m0[2]*fOther[19]*vrelCX+0.7071067811865475*m0[0]*fOther[19]*vrelCX+0.632455532033676*m0[1]*fOther[5]*vrelCX+0.7071067811865475*fOther[2]*m0[2]*vrelCX; 
  prodCX[20] = 0.6324555320336759*m0[2]*fOther[20]*vrelCX+0.7071067811865475*m0[0]*fOther[20]*vrelCX+0.7071067811865475*m0[1]*fOther[12]*vrelCX; 
  prodCX[21] = 0.4517539514526256*m0[2]*fOther[21]*vrelCX+0.7071067811865475*m0[0]*fOther[21]*vrelCX+0.632455532033676*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[2]*fOther[3]*vrelCX; 
  prodCX[22] = 0.7071067811865475*m0[1]*fOther[33]*vrelCX+0.7071067811865475*m0[0]*fOther[22]*vrelCX; 
  prodCX[23] = 0.6324555320336759*m0[2]*fOther[23]*vrelCX+0.7071067811865475*m0[0]*fOther[23]*vrelCX+0.7071067811865475*m0[1]*fOther[13]*vrelCX; 
  prodCX[24] = 0.7071067811865475*m0[1]*fOther[34]*vrelCX+0.7071067811865475*m0[0]*fOther[24]*vrelCX; 
  prodCX[25] = 0.4517539514526256*m0[2]*fOther[25]*vrelCX+0.7071067811865475*m0[0]*fOther[25]*vrelCX+0.632455532033676*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[2]*fOther[4]*vrelCX; 
  prodCX[26] = 0.7071067811865475*m0[1]*fOther[36]*vrelCX+0.7071067811865475*m0[0]*fOther[26]*vrelCX; 
  prodCX[27] = 0.7071067811865475*m0[1]*fOther[39]*vrelCX+0.7071067811865475*m0[0]*fOther[27]*vrelCX; 
  prodCX[28] = 0.6324555320336759*m0[2]*fOther[28]*vrelCX+0.7071067811865475*m0[0]*fOther[28]*vrelCX+0.7071067811865475*m0[1]*fOther[14]*vrelCX; 
  prodCX[29] = 0.7071067811865475*m0[1]*fOther[41]*vrelCX+0.7071067811865475*m0[0]*fOther[29]*vrelCX; 
  prodCX[30] = 0.7071067811865475*m0[1]*fOther[42]*vrelCX+0.7071067811865475*m0[0]*fOther[30]*vrelCX; 
  prodCX[31] = 0.632455532033676*m0[1]*fOther[44]*vrelCX+0.6324555320336759*m0[2]*fOther[31]*vrelCX+0.7071067811865475*m0[0]*fOther[31]*vrelCX+0.7071067811865475*m0[1]*fOther[18]*vrelCX; 
  prodCX[32] = 0.4517539514526256*m0[2]*fOther[32]*vrelCX+0.7071067811865475*m0[0]*fOther[32]*vrelCX+0.6324555320336759*m0[1]*fOther[15]*vrelCX+0.7071067811865475*m0[2]*fOther[7]*vrelCX; 
  prodCX[33] = 0.6324555320336759*m0[2]*fOther[33]*vrelCX+0.7071067811865475*m0[0]*fOther[33]*vrelCX+0.7071067811865475*m0[1]*fOther[22]*vrelCX; 
  prodCX[34] = 0.6324555320336759*m0[2]*fOther[34]*vrelCX+0.7071067811865475*m0[0]*fOther[34]*vrelCX+0.7071067811865475*m0[1]*fOther[24]*vrelCX; 
  prodCX[35] = 0.4517539514526256*m0[2]*fOther[35]*vrelCX+0.7071067811865475*m0[0]*fOther[35]*vrelCX+0.6324555320336759*m0[1]*fOther[16]*vrelCX+0.7071067811865475*m0[2]*fOther[9]*vrelCX; 
  prodCX[36] = 0.6324555320336759*m0[2]*fOther[36]*vrelCX+0.7071067811865475*m0[0]*fOther[36]*vrelCX+0.7071067811865475*m0[1]*fOther[26]*vrelCX; 
  prodCX[37] = 0.4517539514526256*m0[2]*fOther[37]*vrelCX+0.7071067811865475*m0[0]*fOther[37]*vrelCX+0.6324555320336759*m0[1]*fOther[17]*vrelCX+0.7071067811865475*m0[2]*fOther[10]*vrelCX; 
  prodCX[38] = 0.7071067811865475*m0[1]*fOther[45]*vrelCX+0.7071067811865475*m0[0]*fOther[38]*vrelCX; 
  prodCX[39] = 0.6324555320336759*m0[2]*fOther[39]*vrelCX+0.7071067811865475*m0[0]*fOther[39]*vrelCX+0.7071067811865475*m0[1]*fOther[27]*vrelCX; 
  prodCX[40] = 0.7071067811865475*m0[1]*fOther[46]*vrelCX+0.7071067811865475*m0[0]*fOther[40]*vrelCX; 
  prodCX[41] = 0.6324555320336759*m0[2]*fOther[41]*vrelCX+0.7071067811865475*m0[0]*fOther[41]*vrelCX+0.7071067811865475*m0[1]*fOther[29]*vrelCX; 
  prodCX[42] = 0.6324555320336759*m0[2]*fOther[42]*vrelCX+0.7071067811865475*m0[0]*fOther[42]*vrelCX+0.7071067811865475*m0[1]*fOther[30]*vrelCX; 
  prodCX[43] = 0.7071067811865475*m0[1]*fOther[47]*vrelCX+0.7071067811865475*m0[0]*fOther[43]*vrelCX; 
  prodCX[44] = 0.4517539514526256*m0[2]*fOther[44]*vrelCX+0.7071067811865475*m0[0]*fOther[44]*vrelCX+0.632455532033676*m0[1]*fOther[31]*vrelCX+0.7071067811865475*m0[2]*fOther[18]*vrelCX; 
  prodCX[45] = 0.6324555320336759*m0[2]*fOther[45]*vrelCX+0.7071067811865475*m0[0]*fOther[45]*vrelCX+0.7071067811865475*m0[1]*fOther[38]*vrelCX; 
  prodCX[46] = 0.6324555320336759*m0[2]*fOther[46]*vrelCX+0.7071067811865475*m0[0]*fOther[46]*vrelCX+0.7071067811865475*m0[1]*fOther[40]*vrelCX; 
  prodCX[47] = 0.6324555320336759*m0[2]*fOther[47]*vrelCX+0.7071067811865475*m0[0]*fOther[47]*vrelCX+0.7071067811865475*m0[1]*fOther[43]*vrelCX; 
 
} 
void VmProdCXcellAvSer1x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.7071067811865476*vtSq[0]; 
  double xSqAv = (0.5*pow(u[8],2))/vtSqAv-(1.414213562373095*w[3]*u[8])/vtSqAv+(0.5*pow(u[4],2))/vtSqAv-(1.414213562373095*w[2]*u[4])/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv+pow(w[1],2)/vtSqAv-(1.414213562373095*u[0]*w[1])/vtSqAv+(0.5*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.7071067811865475*m0[3]*fOther[31]*vrelCX+0.7071067811865475*m0[2]*fOther[11]*vrelCX+0.7071067811865475*fOther[1]*m0[1]*vrelCX+0.7071067811865475*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.6210590034081186*m0[2]*fOther[31]*vrelCX+0.6210590034081186*m0[3]*fOther[11]*vrelCX+0.6324555320336759*m0[1]*fOther[11]*vrelCX+0.6324555320336759*fOther[1]*m0[2]*vrelCX+0.7071067811865475*fOther[0]*m0[1]*vrelCX+0.7071067811865475*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.7071067811865474*m0[3]*fOther[48]*vrelCX+0.7071067811865475*m0[2]*fOther[19]*vrelCX+0.7071067811865475*m0[1]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.7071067811865474*m0[3]*fOther[50]*vrelCX+0.7071067811865475*m0[2]*fOther[21]*vrelCX+0.7071067811865475*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.7071067811865474*m0[3]*fOther[54]*vrelCX+0.7071067811865475*m0[2]*fOther[25]*vrelCX+0.7071067811865475*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.6210590034081187*m0[2]*fOther[48]*vrelCX+0.6210590034081187*m0[3]*fOther[19]*vrelCX+0.632455532033676*m0[1]*fOther[19]*vrelCX+0.6324555320336759*m0[2]*fOther[5]*vrelCX+0.7071067811865475*m0[0]*fOther[5]*vrelCX+0.7071067811865475*m0[1]*fOther[2]*vrelCX; 
  prodCX[6] = 0.6210590034081187*m0[2]*fOther[50]*vrelCX+0.6210590034081187*m0[3]*fOther[21]*vrelCX+0.632455532033676*m0[1]*fOther[21]*vrelCX+0.6324555320336759*m0[2]*fOther[6]*vrelCX+0.7071067811865475*m0[0]*fOther[6]*vrelCX+0.7071067811865475*m0[1]*fOther[3]*vrelCX; 
  prodCX[7] = 0.7071067811865475*m0[3]*fOther[64]*vrelCX+0.7071067811865475*m0[2]*fOther[36]*vrelCX+0.7071067811865475*m0[1]*fOther[15]*vrelCX+0.7071067811865475*m0[0]*fOther[7]*vrelCX; 
  prodCX[8] = 0.6210590034081187*m0[2]*fOther[54]*vrelCX+0.6210590034081187*m0[3]*fOther[25]*vrelCX+0.632455532033676*m0[1]*fOther[25]*vrelCX+0.6324555320336759*m0[2]*fOther[8]*vrelCX+0.7071067811865475*m0[0]*fOther[8]*vrelCX+0.7071067811865475*m0[1]*fOther[4]*vrelCX; 
  prodCX[9] = 0.7071067811865475*m0[3]*fOther[67]*vrelCX+0.7071067811865475*m0[2]*fOther[39]*vrelCX+0.7071067811865475*m0[1]*fOther[16]*vrelCX+0.7071067811865475*m0[0]*fOther[9]*vrelCX; 
  prodCX[10] = 0.7071067811865475*m0[3]*fOther[69]*vrelCX+0.7071067811865475*m0[2]*fOther[41]*vrelCX+0.7071067811865475*m0[1]*fOther[17]*vrelCX+0.7071067811865475*m0[0]*fOther[10]*vrelCX; 
  prodCX[11] = 0.421637021355784*m0[3]*fOther[31]*vrelCX+0.6210590034081186*m0[1]*fOther[31]*vrelCX+0.4517539514526256*m0[2]*fOther[11]*vrelCX+0.7071067811865475*m0[0]*fOther[11]*vrelCX+0.6210590034081186*fOther[1]*m0[3]*vrelCX+0.7071067811865475*fOther[0]*m0[2]*vrelCX+0.6324555320336759*fOther[1]*m0[1]*vrelCX; 
  prodCX[12] = 0.7071067811865475*m0[1]*fOther[20]*vrelCX+0.7071067811865475*m0[0]*fOther[12]*vrelCX; 
  prodCX[13] = 0.7071067811865475*m0[1]*fOther[23]*vrelCX+0.7071067811865475*m0[0]*fOther[13]*vrelCX; 
  prodCX[14] = 0.7071067811865475*m0[1]*fOther[28]*vrelCX+0.7071067811865475*m0[0]*fOther[14]*vrelCX; 
  prodCX[15] = 0.6210590034081186*m0[2]*fOther[64]*vrelCX+0.6210590034081186*m0[3]*fOther[36]*vrelCX+0.6324555320336759*m0[1]*fOther[36]*vrelCX+0.6324555320336759*m0[2]*fOther[15]*vrelCX+0.7071067811865475*m0[0]*fOther[15]*vrelCX+0.7071067811865475*m0[1]*fOther[7]*vrelCX; 
  prodCX[16] = 0.6210590034081186*m0[2]*fOther[67]*vrelCX+0.6210590034081186*m0[3]*fOther[39]*vrelCX+0.6324555320336759*m0[1]*fOther[39]*vrelCX+0.6324555320336759*m0[2]*fOther[16]*vrelCX+0.7071067811865475*m0[0]*fOther[16]*vrelCX+0.7071067811865475*m0[1]*fOther[9]*vrelCX; 
  prodCX[17] = 0.6210590034081186*m0[2]*fOther[69]*vrelCX+0.6210590034081186*m0[3]*fOther[41]*vrelCX+0.6324555320336759*m0[1]*fOther[41]*vrelCX+0.6324555320336759*m0[2]*fOther[17]*vrelCX+0.7071067811865475*m0[0]*fOther[17]*vrelCX+0.7071067811865475*m0[1]*fOther[10]*vrelCX; 
  prodCX[18] = 0.7071067811865474*m0[3]*fOther[76]*vrelCX+0.7071067811865475*m0[2]*fOther[60]*vrelCX+0.7071067811865475*m0[1]*fOther[35]*vrelCX+0.7071067811865475*m0[0]*fOther[18]*vrelCX; 
  prodCX[19] = 0.4216370213557839*m0[3]*fOther[48]*vrelCX+0.6210590034081187*m0[1]*fOther[48]*vrelCX+0.4517539514526256*m0[2]*fOther[19]*vrelCX+0.7071067811865475*m0[0]*fOther[19]*vrelCX+0.6210590034081187*m0[3]*fOther[5]*vrelCX+0.632455532033676*m0[1]*fOther[5]*vrelCX+0.7071067811865475*fOther[2]*m0[2]*vrelCX; 
  prodCX[20] = 0.6324555320336759*m0[2]*fOther[20]*vrelCX+0.7071067811865475*m0[0]*fOther[20]*vrelCX+0.7071067811865475*m0[1]*fOther[12]*vrelCX; 
  prodCX[21] = 0.4216370213557839*m0[3]*fOther[50]*vrelCX+0.6210590034081187*m0[1]*fOther[50]*vrelCX+0.4517539514526256*m0[2]*fOther[21]*vrelCX+0.7071067811865475*m0[0]*fOther[21]*vrelCX+0.6210590034081187*m0[3]*fOther[6]*vrelCX+0.632455532033676*m0[1]*fOther[6]*vrelCX+0.7071067811865475*m0[2]*fOther[3]*vrelCX; 
  prodCX[22] = 0.7071067811865475*m0[1]*fOther[37]*vrelCX+0.7071067811865475*m0[0]*fOther[22]*vrelCX; 
  prodCX[23] = 0.6324555320336759*m0[2]*fOther[23]*vrelCX+0.7071067811865475*m0[0]*fOther[23]*vrelCX+0.7071067811865475*m0[1]*fOther[13]*vrelCX; 
  prodCX[24] = 0.7071067811865475*m0[1]*fOther[38]*vrelCX+0.7071067811865475*m0[0]*fOther[24]*vrelCX; 
  prodCX[25] = 0.4216370213557839*m0[3]*fOther[54]*vrelCX+0.6210590034081187*m0[1]*fOther[54]*vrelCX+0.4517539514526256*m0[2]*fOther[25]*vrelCX+0.7071067811865475*m0[0]*fOther[25]*vrelCX+0.6210590034081187*m0[3]*fOther[8]*vrelCX+0.632455532033676*m0[1]*fOther[8]*vrelCX+0.7071067811865475*m0[2]*fOther[4]*vrelCX; 
  prodCX[26] = 0.7071067811865475*m0[1]*fOther[40]*vrelCX+0.7071067811865475*m0[0]*fOther[26]*vrelCX; 
  prodCX[27] = 0.7071067811865475*m0[1]*fOther[43]*vrelCX+0.7071067811865475*m0[0]*fOther[27]*vrelCX; 
  prodCX[28] = 0.6324555320336759*m0[2]*fOther[28]*vrelCX+0.7071067811865475*m0[0]*fOther[28]*vrelCX+0.7071067811865475*m0[1]*fOther[14]*vrelCX; 
  prodCX[29] = 0.7071067811865475*m0[1]*fOther[45]*vrelCX+0.7071067811865475*m0[0]*fOther[29]*vrelCX; 
  prodCX[30] = 0.7071067811865475*m0[1]*fOther[46]*vrelCX+0.7071067811865475*m0[0]*fOther[30]*vrelCX; 
  prodCX[31] = 0.421637021355784*m0[2]*fOther[31]*vrelCX+0.7071067811865475*m0[0]*fOther[31]*vrelCX+0.421637021355784*m0[3]*fOther[11]*vrelCX+0.6210590034081186*m0[1]*fOther[11]*vrelCX+0.7071067811865475*fOther[0]*m0[3]*vrelCX+0.6210590034081186*fOther[1]*m0[2]*vrelCX; 
  prodCX[32] = 0.7071067811865474*m0[1]*fOther[49]*vrelCX+0.7071067811865475*m0[0]*fOther[32]*vrelCX; 
  prodCX[33] = 0.7071067811865474*m0[1]*fOther[52]*vrelCX+0.7071067811865475*m0[0]*fOther[33]*vrelCX; 
  prodCX[34] = 0.7071067811865474*m0[1]*fOther[57]*vrelCX+0.7071067811865475*m0[0]*fOther[34]*vrelCX; 
  prodCX[35] = 0.6210590034081187*m0[2]*fOther[76]*vrelCX+0.6210590034081187*m0[3]*fOther[60]*vrelCX+0.632455532033676*m0[1]*fOther[60]*vrelCX+0.6324555320336759*m0[2]*fOther[35]*vrelCX+0.7071067811865475*m0[0]*fOther[35]*vrelCX+0.7071067811865475*m0[1]*fOther[18]*vrelCX; 
  prodCX[36] = 0.421637021355784*m0[3]*fOther[64]*vrelCX+0.6210590034081186*m0[1]*fOther[64]*vrelCX+0.4517539514526256*m0[2]*fOther[36]*vrelCX+0.7071067811865475*m0[0]*fOther[36]*vrelCX+0.6210590034081186*m0[3]*fOther[15]*vrelCX+0.6324555320336759*m0[1]*fOther[15]*vrelCX+0.7071067811865475*m0[2]*fOther[7]*vrelCX; 
  prodCX[37] = 0.6324555320336759*m0[2]*fOther[37]*vrelCX+0.7071067811865475*m0[0]*fOther[37]*vrelCX+0.7071067811865475*m0[1]*fOther[22]*vrelCX; 
  prodCX[38] = 0.6324555320336759*m0[2]*fOther[38]*vrelCX+0.7071067811865475*m0[0]*fOther[38]*vrelCX+0.7071067811865475*m0[1]*fOther[24]*vrelCX; 
  prodCX[39] = 0.421637021355784*m0[3]*fOther[67]*vrelCX+0.6210590034081186*m0[1]*fOther[67]*vrelCX+0.4517539514526256*m0[2]*fOther[39]*vrelCX+0.7071067811865475*m0[0]*fOther[39]*vrelCX+0.6210590034081186*m0[3]*fOther[16]*vrelCX+0.6324555320336759*m0[1]*fOther[16]*vrelCX+0.7071067811865475*m0[2]*fOther[9]*vrelCX; 
  prodCX[40] = 0.6324555320336759*m0[2]*fOther[40]*vrelCX+0.7071067811865475*m0[0]*fOther[40]*vrelCX+0.7071067811865475*m0[1]*fOther[26]*vrelCX; 
  prodCX[41] = 0.421637021355784*m0[3]*fOther[69]*vrelCX+0.6210590034081186*m0[1]*fOther[69]*vrelCX+0.4517539514526256*m0[2]*fOther[41]*vrelCX+0.7071067811865475*m0[0]*fOther[41]*vrelCX+0.6210590034081186*m0[3]*fOther[17]*vrelCX+0.6324555320336759*m0[1]*fOther[17]*vrelCX+0.7071067811865475*m0[2]*fOther[10]*vrelCX; 
  prodCX[42] = 0.7071067811865475*m0[1]*fOther[61]*vrelCX+0.7071067811865475*m0[0]*fOther[42]*vrelCX; 
  prodCX[43] = 0.6324555320336759*m0[2]*fOther[43]*vrelCX+0.7071067811865475*m0[0]*fOther[43]*vrelCX+0.7071067811865475*m0[1]*fOther[27]*vrelCX; 
  prodCX[44] = 0.7071067811865475*m0[1]*fOther[62]*vrelCX+0.7071067811865475*m0[0]*fOther[44]*vrelCX; 
  prodCX[45] = 0.6324555320336759*m0[2]*fOther[45]*vrelCX+0.7071067811865475*m0[0]*fOther[45]*vrelCX+0.7071067811865475*m0[1]*fOther[29]*vrelCX; 
  prodCX[46] = 0.6324555320336759*m0[2]*fOther[46]*vrelCX+0.7071067811865475*m0[0]*fOther[46]*vrelCX+0.7071067811865475*m0[1]*fOther[30]*vrelCX; 
  prodCX[47] = 0.7071067811865475*m0[1]*fOther[63]*vrelCX+0.7071067811865475*m0[0]*fOther[47]*vrelCX; 
  prodCX[48] = 0.421637021355784*m0[2]*fOther[48]*vrelCX+0.7071067811865475*m0[0]*fOther[48]*vrelCX+0.4216370213557839*m0[3]*fOther[19]*vrelCX+0.6210590034081187*m0[1]*fOther[19]*vrelCX+0.6210590034081187*m0[2]*fOther[5]*vrelCX+0.7071067811865474*fOther[2]*m0[3]*vrelCX; 
  prodCX[49] = 0.6324555320336759*m0[2]*fOther[49]*vrelCX+0.7071067811865475*m0[0]*fOther[49]*vrelCX+0.7071067811865474*m0[1]*fOther[32]*vrelCX; 
  prodCX[50] = 0.421637021355784*m0[2]*fOther[50]*vrelCX+0.7071067811865475*m0[0]*fOther[50]*vrelCX+0.4216370213557839*m0[3]*fOther[21]*vrelCX+0.6210590034081187*m0[1]*fOther[21]*vrelCX+0.6210590034081187*m0[2]*fOther[6]*vrelCX+0.7071067811865474*fOther[3]*m0[3]*vrelCX; 
  prodCX[51] = 0.7071067811865474*m0[1]*fOther[65]*vrelCX+0.7071067811865475*m0[0]*fOther[51]*vrelCX; 
  prodCX[52] = 0.6324555320336759*m0[2]*fOther[52]*vrelCX+0.7071067811865475*m0[0]*fOther[52]*vrelCX+0.7071067811865474*m0[1]*fOther[33]*vrelCX; 
  prodCX[53] = 0.7071067811865474*m0[1]*fOther[66]*vrelCX+0.7071067811865475*m0[0]*fOther[53]*vrelCX; 
  prodCX[54] = 0.421637021355784*m0[2]*fOther[54]*vrelCX+0.7071067811865475*m0[0]*fOther[54]*vrelCX+0.4216370213557839*m0[3]*fOther[25]*vrelCX+0.6210590034081187*m0[1]*fOther[25]*vrelCX+0.6210590034081187*m0[2]*fOther[8]*vrelCX+0.7071067811865474*m0[3]*fOther[4]*vrelCX; 
  prodCX[55] = 0.7071067811865474*m0[1]*fOther[68]*vrelCX+0.7071067811865475*m0[0]*fOther[55]*vrelCX; 
  prodCX[56] = 0.7071067811865474*m0[1]*fOther[71]*vrelCX+0.7071067811865475*m0[0]*fOther[56]*vrelCX; 
  prodCX[57] = 0.6324555320336759*m0[2]*fOther[57]*vrelCX+0.7071067811865475*m0[0]*fOther[57]*vrelCX+0.7071067811865474*m0[1]*fOther[34]*vrelCX; 
  prodCX[58] = 0.7071067811865474*m0[1]*fOther[73]*vrelCX+0.7071067811865475*m0[0]*fOther[58]*vrelCX; 
  prodCX[59] = 0.7071067811865474*m0[1]*fOther[74]*vrelCX+0.7071067811865475*m0[0]*fOther[59]*vrelCX; 
  prodCX[60] = 0.4216370213557839*m0[3]*fOther[76]*vrelCX+0.6210590034081187*m0[1]*fOther[76]*vrelCX+0.4517539514526256*m0[2]*fOther[60]*vrelCX+0.7071067811865475*m0[0]*fOther[60]*vrelCX+0.6210590034081187*m0[3]*fOther[35]*vrelCX+0.632455532033676*m0[1]*fOther[35]*vrelCX+0.7071067811865475*m0[2]*fOther[18]*vrelCX; 
  prodCX[61] = 0.6324555320336759*m0[2]*fOther[61]*vrelCX+0.7071067811865475*m0[0]*fOther[61]*vrelCX+0.7071067811865475*m0[1]*fOther[42]*vrelCX; 
  prodCX[62] = 0.6324555320336759*m0[2]*fOther[62]*vrelCX+0.7071067811865475*m0[0]*fOther[62]*vrelCX+0.7071067811865475*m0[1]*fOther[44]*vrelCX; 
  prodCX[63] = 0.6324555320336759*m0[2]*fOther[63]*vrelCX+0.7071067811865475*m0[0]*fOther[63]*vrelCX+0.7071067811865475*m0[1]*fOther[47]*vrelCX; 
  prodCX[64] = 0.421637021355784*m0[2]*fOther[64]*vrelCX+0.7071067811865475*m0[0]*fOther[64]*vrelCX+0.421637021355784*m0[3]*fOther[36]*vrelCX+0.6210590034081186*m0[1]*fOther[36]*vrelCX+0.6210590034081186*m0[2]*fOther[15]*vrelCX+0.7071067811865475*m0[3]*fOther[7]*vrelCX; 
  prodCX[65] = 0.6324555320336759*m0[2]*fOther[65]*vrelCX+0.7071067811865475*m0[0]*fOther[65]*vrelCX+0.7071067811865474*m0[1]*fOther[51]*vrelCX; 
  prodCX[66] = 0.6324555320336759*m0[2]*fOther[66]*vrelCX+0.7071067811865475*m0[0]*fOther[66]*vrelCX+0.7071067811865474*m0[1]*fOther[53]*vrelCX; 
  prodCX[67] = 0.421637021355784*m0[2]*fOther[67]*vrelCX+0.7071067811865475*m0[0]*fOther[67]*vrelCX+0.421637021355784*m0[3]*fOther[39]*vrelCX+0.6210590034081186*m0[1]*fOther[39]*vrelCX+0.6210590034081186*m0[2]*fOther[16]*vrelCX+0.7071067811865475*m0[3]*fOther[9]*vrelCX; 
  prodCX[68] = 0.6324555320336759*m0[2]*fOther[68]*vrelCX+0.7071067811865475*m0[0]*fOther[68]*vrelCX+0.7071067811865474*m0[1]*fOther[55]*vrelCX; 
  prodCX[69] = 0.421637021355784*m0[2]*fOther[69]*vrelCX+0.7071067811865475*m0[0]*fOther[69]*vrelCX+0.421637021355784*m0[3]*fOther[41]*vrelCX+0.6210590034081186*m0[1]*fOther[41]*vrelCX+0.6210590034081186*m0[2]*fOther[17]*vrelCX+0.7071067811865475*m0[3]*fOther[10]*vrelCX; 
  prodCX[70] = 0.7071067811865474*m0[1]*fOther[77]*vrelCX+0.7071067811865475*m0[0]*fOther[70]*vrelCX; 
  prodCX[71] = 0.6324555320336759*m0[2]*fOther[71]*vrelCX+0.7071067811865475*m0[0]*fOther[71]*vrelCX+0.7071067811865474*m0[1]*fOther[56]*vrelCX; 
  prodCX[72] = 0.7071067811865474*m0[1]*fOther[78]*vrelCX+0.7071067811865475*m0[0]*fOther[72]*vrelCX; 
  prodCX[73] = 0.6324555320336759*m0[2]*fOther[73]*vrelCX+0.7071067811865475*m0[0]*fOther[73]*vrelCX+0.7071067811865474*m0[1]*fOther[58]*vrelCX; 
  prodCX[74] = 0.6324555320336759*m0[2]*fOther[74]*vrelCX+0.7071067811865475*m0[0]*fOther[74]*vrelCX+0.7071067811865474*m0[1]*fOther[59]*vrelCX; 
  prodCX[75] = 0.7071067811865474*m0[1]*fOther[79]*vrelCX+0.7071067811865475*m0[0]*fOther[75]*vrelCX; 
  prodCX[76] = 0.421637021355784*m0[2]*fOther[76]*vrelCX+0.7071067811865475*m0[0]*fOther[76]*vrelCX+0.4216370213557839*m0[3]*fOther[60]*vrelCX+0.6210590034081187*m0[1]*fOther[60]*vrelCX+0.6210590034081187*m0[2]*fOther[35]*vrelCX+0.7071067811865474*m0[3]*fOther[18]*vrelCX; 
  prodCX[77] = 0.6324555320336759*m0[2]*fOther[77]*vrelCX+0.7071067811865475*m0[0]*fOther[77]*vrelCX+0.7071067811865474*m0[1]*fOther[70]*vrelCX; 
  prodCX[78] = 0.6324555320336759*m0[2]*fOther[78]*vrelCX+0.7071067811865475*m0[0]*fOther[78]*vrelCX+0.7071067811865474*m0[1]*fOther[72]*vrelCX; 
  prodCX[79] = 0.6324555320336759*m0[2]*fOther[79]*vrelCX+0.7071067811865475*m0[0]*fOther[79]*vrelCX+0.7071067811865474*m0[1]*fOther[75]*vrelCX; 
 
} 
