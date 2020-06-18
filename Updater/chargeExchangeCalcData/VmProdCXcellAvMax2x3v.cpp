#include <RelativeVelocityModDecl.h> 
#include <math.h> 
void VmProdCXcellAvMax2x3v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[5]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = (0.25*pow(u[6],2))/vtSqAv-(1.0*w[4]*u[6])/vtSqAv+pow(w[4],2)/vtSqAv+pow(w[3],2)/vtSqAv-(1.0*u[3]*w[3])/vtSqAv+(0.25*pow(u[3],2))/vtSqAv+pow(w[2],2)/vtSqAv-(1.0*u[0]*w[2])/vtSqAv+(0.25*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.5*fOther[2]*m0[2]*vrelCX+0.5*fOther[1]*m0[1]*vrelCX+0.5*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.5*fOther[0]*m0[1]*vrelCX+0.5*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.5*fOther[0]*m0[2]*vrelCX+0.5*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.5*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.5*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.5*m0[0]*fOther[5]*vrelCX; 
 
} 
void VmProdCXcellAvMax2x3v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[5]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = (0.25*pow(u[12],2))/vtSqAv-(1.0*w[4]*u[12])/vtSqAv+(0.25*pow(u[6],2))/vtSqAv-(1.0*w[3]*u[6])/vtSqAv+pow(w[4],2)/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv-(1.0*u[0]*w[2])/vtSqAv+(0.25*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.5*m0[5]*fOther[17]*vrelCX+0.5*m0[4]*fOther[16]*vrelCX+0.5*m0[3]*fOther[6]*vrelCX+0.5*fOther[2]*m0[2]*vrelCX+0.5*fOther[1]*m0[1]*vrelCX+0.5*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.4472135954999579*m0[1]*fOther[16]*vrelCX+0.5*m0[2]*fOther[6]*vrelCX+0.4472135954999579*fOther[1]*m0[4]*vrelCX+0.5*fOther[2]*m0[3]*vrelCX+0.5*fOther[0]*m0[1]*vrelCX+0.5*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.4472135954999579*m0[2]*fOther[17]*vrelCX+0.5*m0[1]*fOther[6]*vrelCX+0.4472135954999579*fOther[2]*m0[5]*vrelCX+0.5*fOther[1]*m0[3]*vrelCX+0.5*fOther[0]*m0[2]*vrelCX+0.5*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.5*m0[2]*fOther[8]*vrelCX+0.5*m0[1]*fOther[7]*vrelCX+0.5*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.5*m0[2]*fOther[10]*vrelCX+0.5*m0[1]*fOther[9]*vrelCX+0.5*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.5*m0[2]*fOther[13]*vrelCX+0.5*m0[1]*fOther[12]*vrelCX+0.5*m0[0]*fOther[5]*vrelCX; 
  prodCX[6] = 0.4472135954999579*m0[3]*fOther[17]*vrelCX+0.4472135954999579*m0[3]*fOther[16]*vrelCX+0.4472135954999579*m0[5]*fOther[6]*vrelCX+0.4472135954999579*m0[4]*fOther[6]*vrelCX+0.5*m0[0]*fOther[6]*vrelCX+0.5*fOther[0]*m0[3]*vrelCX+0.5*fOther[1]*m0[2]*vrelCX+0.5*m0[1]*fOther[2]*vrelCX; 
  prodCX[7] = 0.5*m0[3]*fOther[8]*vrelCX+0.4472135954999579*m0[4]*fOther[7]*vrelCX+0.5*m0[0]*fOther[7]*vrelCX+0.5*m0[1]*fOther[3]*vrelCX; 
  prodCX[8] = 0.4472135954999579*m0[5]*fOther[8]*vrelCX+0.5*m0[0]*fOther[8]*vrelCX+0.5*m0[3]*fOther[7]*vrelCX+0.5*m0[2]*fOther[3]*vrelCX; 
  prodCX[9] = 0.5*m0[3]*fOther[10]*vrelCX+0.4472135954999579*m0[4]*fOther[9]*vrelCX+0.5*m0[0]*fOther[9]*vrelCX+0.5*m0[1]*fOther[4]*vrelCX; 
  prodCX[10] = 0.4472135954999579*m0[5]*fOther[10]*vrelCX+0.5*m0[0]*fOther[10]*vrelCX+0.5*m0[3]*fOther[9]*vrelCX+0.5*m0[2]*fOther[4]*vrelCX; 
  prodCX[11] = 0.5*m0[0]*fOther[11]*vrelCX; 
  prodCX[12] = 0.5*m0[3]*fOther[13]*vrelCX+0.4472135954999579*m0[4]*fOther[12]*vrelCX+0.5*m0[0]*fOther[12]*vrelCX+0.5*m0[1]*fOther[5]*vrelCX; 
  prodCX[13] = 0.4472135954999579*m0[5]*fOther[13]*vrelCX+0.5*m0[0]*fOther[13]*vrelCX+0.5*m0[3]*fOther[12]*vrelCX+0.5*m0[2]*fOther[5]*vrelCX; 
  prodCX[14] = 0.5*m0[0]*fOther[14]*vrelCX; 
  prodCX[15] = 0.5*m0[0]*fOther[15]*vrelCX; 
  prodCX[16] = 0.31943828249997*m0[4]*fOther[16]*vrelCX+0.5*m0[0]*fOther[16]*vrelCX+0.4472135954999579*m0[3]*fOther[6]*vrelCX+0.5*fOther[0]*m0[4]*vrelCX+0.4472135954999579*fOther[1]*m0[1]*vrelCX; 
  prodCX[17] = 0.31943828249997*m0[5]*fOther[17]*vrelCX+0.5*m0[0]*fOther[17]*vrelCX+0.4472135954999579*m0[3]*fOther[6]*vrelCX+0.5*fOther[0]*m0[5]*vrelCX+0.4472135954999579*fOther[2]*m0[2]*vrelCX; 
  prodCX[18] = 0.5*m0[0]*fOther[18]*vrelCX; 
  prodCX[19] = 0.5*m0[0]*fOther[19]*vrelCX; 
  prodCX[20] = 0.5*m0[0]*fOther[20]*vrelCX; 
 
} 
void VmProdCXcellAvMax2x3v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[5]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = (0.25*pow(u[20],2))/vtSqAv-(1.0*w[4]*u[20])/vtSqAv+(0.25*pow(u[10],2))/vtSqAv-(1.0*w[3]*u[10])/vtSqAv+pow(w[4],2)/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv-(1.0*u[0]*w[2])/vtSqAv+(0.25*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.5*m0[9]*fOther[52]*vrelCX+0.5*m0[8]*fOther[51]*vrelCX+0.5*m0[7]*fOther[32]*vrelCX+0.5*m0[6]*fOther[31]*vrelCX+0.5*m0[5]*fOther[17]*vrelCX+0.5*m0[4]*fOther[16]*vrelCX+0.5*m0[3]*fOther[6]*vrelCX+0.5*fOther[2]*m0[2]*vrelCX+0.5*fOther[1]*m0[1]*vrelCX+0.5*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.4391550328268398*m0[4]*fOther[51]*vrelCX+0.5000000000000001*m0[5]*fOther[32]*vrelCX+0.447213595499958*m0[3]*fOther[31]*vrelCX+0.5000000000000001*m0[7]*fOther[17]*vrelCX+0.4391550328268398*m0[8]*fOther[16]*vrelCX+0.4472135954999579*m0[1]*fOther[16]*vrelCX+0.447213595499958*fOther[6]*m0[6]*vrelCX+0.5*m0[2]*fOther[6]*vrelCX+0.4472135954999579*fOther[1]*m0[4]*vrelCX+0.5*fOther[2]*m0[3]*vrelCX+0.5*fOther[0]*m0[1]*vrelCX+0.5*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.4391550328268398*m0[5]*fOther[52]*vrelCX+0.447213595499958*m0[3]*fOther[32]*vrelCX+0.5000000000000001*m0[4]*fOther[31]*vrelCX+0.4391550328268398*m0[9]*fOther[17]*vrelCX+0.4472135954999579*m0[2]*fOther[17]*vrelCX+0.5000000000000001*m0[6]*fOther[16]*vrelCX+0.447213595499958*fOther[6]*m0[7]*vrelCX+0.5*m0[1]*fOther[6]*vrelCX+0.4472135954999579*fOther[2]*m0[5]*vrelCX+0.5*fOther[1]*m0[3]*vrelCX+0.5*fOther[0]*m0[2]*vrelCX+0.5*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.5000000000000001*m0[5]*fOther[34]*vrelCX+0.5000000000000001*m0[4]*fOther[33]*vrelCX+0.5*m0[3]*fOther[21]*vrelCX+0.5*m0[2]*fOther[8]*vrelCX+0.5*m0[1]*fOther[7]*vrelCX+0.5*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.5000000000000001*m0[5]*fOther[38]*vrelCX+0.5000000000000001*m0[4]*fOther[37]*vrelCX+0.5*m0[3]*fOther[22]*vrelCX+0.5*m0[2]*fOther[10]*vrelCX+0.5*m0[1]*fOther[9]*vrelCX+0.5*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.5000000000000001*m0[5]*fOther[44]*vrelCX+0.5000000000000001*m0[4]*fOther[43]*vrelCX+0.5*m0[3]*fOther[25]*vrelCX+0.5*m0[2]*fOther[13]*vrelCX+0.5*m0[1]*fOther[12]*vrelCX+0.5*m0[0]*fOther[5]*vrelCX; 
  prodCX[6] = 0.4391550328268399*m0[7]*fOther[52]*vrelCX+0.4391550328268399*m0[6]*fOther[51]*vrelCX+0.4391550328268399*m0[9]*fOther[32]*vrelCX+0.4*m0[6]*fOther[32]*vrelCX+0.447213595499958*m0[2]*fOther[32]*vrelCX+0.4391550328268399*m0[8]*fOther[31]*vrelCX+0.4*m0[7]*fOther[31]*vrelCX+0.447213595499958*m0[1]*fOther[31]*vrelCX+0.4472135954999579*m0[3]*fOther[17]*vrelCX+0.4472135954999579*m0[3]*fOther[16]*vrelCX+0.447213595499958*fOther[2]*m0[7]*vrelCX+0.447213595499958*fOther[1]*m0[6]*vrelCX+0.4472135954999579*m0[5]*fOther[6]*vrelCX+0.4472135954999579*m0[4]*fOther[6]*vrelCX+0.5*m0[0]*fOther[6]*vrelCX+0.5*fOther[0]*m0[3]*vrelCX+0.5*fOther[1]*m0[2]*vrelCX+0.5*m0[1]*fOther[2]*vrelCX; 
  prodCX[7] = 0.5*m0[7]*fOther[34]*vrelCX+0.4391550328268399*m0[8]*fOther[33]*vrelCX+0.447213595499958*m0[1]*fOther[33]*vrelCX+0.447213595499958*m0[6]*fOther[21]*vrelCX+0.5*m0[2]*fOther[21]*vrelCX+0.5*m0[3]*fOther[8]*vrelCX+0.4472135954999579*m0[4]*fOther[7]*vrelCX+0.5*m0[0]*fOther[7]*vrelCX+0.5*m0[1]*fOther[3]*vrelCX; 
  prodCX[8] = 0.4391550328268399*m0[9]*fOther[34]*vrelCX+0.447213595499958*m0[2]*fOther[34]*vrelCX+0.5*m0[6]*fOther[33]*vrelCX+0.447213595499958*m0[7]*fOther[21]*vrelCX+0.5*m0[1]*fOther[21]*vrelCX+0.4472135954999579*m0[5]*fOther[8]*vrelCX+0.5*m0[0]*fOther[8]*vrelCX+0.5*m0[3]*fOther[7]*vrelCX+0.5*m0[2]*fOther[3]*vrelCX; 
  prodCX[9] = 0.5*m0[7]*fOther[38]*vrelCX+0.4391550328268399*m0[8]*fOther[37]*vrelCX+0.447213595499958*m0[1]*fOther[37]*vrelCX+0.447213595499958*m0[6]*fOther[22]*vrelCX+0.5*m0[2]*fOther[22]*vrelCX+0.5*m0[3]*fOther[10]*vrelCX+0.4472135954999579*m0[4]*fOther[9]*vrelCX+0.5*m0[0]*fOther[9]*vrelCX+0.5*m0[1]*fOther[4]*vrelCX; 
  prodCX[10] = 0.4391550328268399*m0[9]*fOther[38]*vrelCX+0.447213595499958*m0[2]*fOther[38]*vrelCX+0.5*m0[6]*fOther[37]*vrelCX+0.447213595499958*m0[7]*fOther[22]*vrelCX+0.5*m0[1]*fOther[22]*vrelCX+0.4472135954999579*m0[5]*fOther[10]*vrelCX+0.5*m0[0]*fOther[10]*vrelCX+0.5*m0[3]*fOther[9]*vrelCX+0.5*m0[2]*fOther[4]*vrelCX; 
  prodCX[11] = 0.5*m0[2]*fOther[24]*vrelCX+0.5*m0[1]*fOther[23]*vrelCX+0.5*m0[0]*fOther[11]*vrelCX; 
  prodCX[12] = 0.5*m0[7]*fOther[44]*vrelCX+0.4391550328268399*m0[8]*fOther[43]*vrelCX+0.447213595499958*m0[1]*fOther[43]*vrelCX+0.447213595499958*m0[6]*fOther[25]*vrelCX+0.5*m0[2]*fOther[25]*vrelCX+0.5*m0[3]*fOther[13]*vrelCX+0.4472135954999579*m0[4]*fOther[12]*vrelCX+0.5*m0[0]*fOther[12]*vrelCX+0.5*m0[1]*fOther[5]*vrelCX; 
  prodCX[13] = 0.4391550328268399*m0[9]*fOther[44]*vrelCX+0.447213595499958*m0[2]*fOther[44]*vrelCX+0.5*m0[6]*fOther[43]*vrelCX+0.447213595499958*m0[7]*fOther[25]*vrelCX+0.5*m0[1]*fOther[25]*vrelCX+0.4472135954999579*m0[5]*fOther[13]*vrelCX+0.5*m0[0]*fOther[13]*vrelCX+0.5*m0[3]*fOther[12]*vrelCX+0.5*m0[2]*fOther[5]*vrelCX; 
  prodCX[14] = 0.5*m0[2]*fOther[27]*vrelCX+0.5*m0[1]*fOther[26]*vrelCX+0.5*m0[0]*fOther[14]*vrelCX; 
  prodCX[15] = 0.5*m0[2]*fOther[29]*vrelCX+0.5*m0[1]*fOther[28]*vrelCX+0.5*m0[0]*fOther[15]*vrelCX; 
  prodCX[16] = 0.2981423969999719*m0[8]*fOther[51]*vrelCX+0.4391550328268398*m0[1]*fOther[51]*vrelCX+0.4472135954999579*m0[7]*fOther[32]*vrelCX+0.31943828249997*m0[6]*fOther[31]*vrelCX+0.5000000000000001*m0[2]*fOther[31]*vrelCX+0.31943828249997*m0[4]*fOther[16]*vrelCX+0.5*m0[0]*fOther[16]*vrelCX+0.4391550328268398*fOther[1]*m0[8]*vrelCX+0.5000000000000001*fOther[2]*m0[6]*vrelCX+0.4472135954999579*m0[3]*fOther[6]*vrelCX+0.5*fOther[0]*m0[4]*vrelCX+0.4472135954999579*fOther[1]*m0[1]*vrelCX; 
  prodCX[17] = 0.2981423969999719*m0[9]*fOther[52]*vrelCX+0.4391550328268398*m0[2]*fOther[52]*vrelCX+0.31943828249997*m0[7]*fOther[32]*vrelCX+0.5000000000000001*m0[1]*fOther[32]*vrelCX+0.4472135954999579*m0[6]*fOther[31]*vrelCX+0.31943828249997*m0[5]*fOther[17]*vrelCX+0.5*m0[0]*fOther[17]*vrelCX+0.4391550328268398*fOther[2]*m0[9]*vrelCX+0.5000000000000001*fOther[1]*m0[7]*vrelCX+0.4472135954999579*m0[3]*fOther[6]*vrelCX+0.5*fOther[0]*m0[5]*vrelCX+0.4472135954999579*fOther[2]*m0[2]*vrelCX; 
  prodCX[18] = 0.5000000000000001*m0[2]*fOther[36]*vrelCX+0.5000000000000001*m0[1]*fOther[35]*vrelCX+0.5*m0[0]*fOther[18]*vrelCX; 
  prodCX[19] = 0.5000000000000001*m0[2]*fOther[41]*vrelCX+0.5000000000000001*m0[1]*fOther[40]*vrelCX+0.5*m0[0]*fOther[19]*vrelCX; 
  prodCX[20] = 0.5000000000000001*m0[2]*fOther[48]*vrelCX+0.5000000000000001*m0[1]*fOther[47]*vrelCX+0.5*m0[0]*fOther[20]*vrelCX; 
  prodCX[21] = 0.447213595499958*m0[3]*fOther[34]*vrelCX+0.447213595499958*m0[3]*fOther[33]*vrelCX+0.4472135954999579*m0[5]*fOther[21]*vrelCX+0.4472135954999579*m0[4]*fOther[21]*vrelCX+0.5*m0[0]*fOther[21]*vrelCX+0.447213595499958*m0[7]*fOther[8]*vrelCX+0.5*m0[1]*fOther[8]*vrelCX+0.447213595499958*m0[6]*fOther[7]*vrelCX+0.5*m0[2]*fOther[7]*vrelCX+0.5*fOther[3]*m0[3]*vrelCX; 
  prodCX[22] = 0.447213595499958*m0[3]*fOther[38]*vrelCX+0.447213595499958*m0[3]*fOther[37]*vrelCX+0.4472135954999579*m0[5]*fOther[22]*vrelCX+0.4472135954999579*m0[4]*fOther[22]*vrelCX+0.5*m0[0]*fOther[22]*vrelCX+0.447213595499958*m0[7]*fOther[10]*vrelCX+0.5*m0[1]*fOther[10]*vrelCX+0.447213595499958*m0[6]*fOther[9]*vrelCX+0.5*m0[2]*fOther[9]*vrelCX+0.5*m0[3]*fOther[4]*vrelCX; 
  prodCX[23] = 0.5*m0[3]*fOther[24]*vrelCX+0.4472135954999579*m0[4]*fOther[23]*vrelCX+0.5*m0[0]*fOther[23]*vrelCX+0.5*m0[1]*fOther[11]*vrelCX; 
  prodCX[24] = 0.4472135954999579*m0[5]*fOther[24]*vrelCX+0.5*m0[0]*fOther[24]*vrelCX+0.5*m0[3]*fOther[23]*vrelCX+0.5*m0[2]*fOther[11]*vrelCX; 
  prodCX[25] = 0.447213595499958*m0[3]*fOther[44]*vrelCX+0.447213595499958*m0[3]*fOther[43]*vrelCX+0.4472135954999579*m0[5]*fOther[25]*vrelCX+0.4472135954999579*m0[4]*fOther[25]*vrelCX+0.5*m0[0]*fOther[25]*vrelCX+0.447213595499958*m0[7]*fOther[13]*vrelCX+0.5*m0[1]*fOther[13]*vrelCX+0.447213595499958*m0[6]*fOther[12]*vrelCX+0.5*m0[2]*fOther[12]*vrelCX+0.5*m0[3]*fOther[5]*vrelCX; 
  prodCX[26] = 0.5*m0[3]*fOther[27]*vrelCX+0.4472135954999579*m0[4]*fOther[26]*vrelCX+0.5*m0[0]*fOther[26]*vrelCX+0.5*m0[1]*fOther[14]*vrelCX; 
  prodCX[27] = 0.4472135954999579*m0[5]*fOther[27]*vrelCX+0.5*m0[0]*fOther[27]*vrelCX+0.5*m0[3]*fOther[26]*vrelCX+0.5*m0[2]*fOther[14]*vrelCX; 
  prodCX[28] = 0.5*m0[3]*fOther[29]*vrelCX+0.4472135954999579*m0[4]*fOther[28]*vrelCX+0.5*m0[0]*fOther[28]*vrelCX+0.5*m0[1]*fOther[15]*vrelCX; 
  prodCX[29] = 0.4472135954999579*m0[5]*fOther[29]*vrelCX+0.5*m0[0]*fOther[29]*vrelCX+0.5*m0[3]*fOther[28]*vrelCX+0.5*m0[2]*fOther[15]*vrelCX; 
  prodCX[30] = 0.5*m0[0]*fOther[30]*vrelCX; 
  prodCX[31] = 0.4391550328268399*m0[3]*fOther[51]*vrelCX+0.4*m0[3]*fOther[32]*vrelCX+0.4472135954999579*m0[5]*fOther[31]*vrelCX+0.31943828249997*m0[4]*fOther[31]*vrelCX+0.5*m0[0]*fOther[31]*vrelCX+0.4472135954999579*m0[6]*fOther[17]*vrelCX+0.31943828249997*m0[6]*fOther[16]*vrelCX+0.5000000000000001*m0[2]*fOther[16]*vrelCX+0.4391550328268399*fOther[6]*m0[8]*vrelCX+0.4*fOther[6]*m0[7]*vrelCX+0.5*fOther[0]*m0[6]*vrelCX+0.447213595499958*m0[1]*fOther[6]*vrelCX+0.5000000000000001*fOther[2]*m0[4]*vrelCX+0.447213595499958*fOther[1]*m0[3]*vrelCX; 
  prodCX[32] = 0.4391550328268399*m0[3]*fOther[52]*vrelCX+0.31943828249997*m0[5]*fOther[32]*vrelCX+0.4472135954999579*m0[4]*fOther[32]*vrelCX+0.5*m0[0]*fOther[32]*vrelCX+0.4*m0[3]*fOther[31]*vrelCX+0.31943828249997*m0[7]*fOther[17]*vrelCX+0.5000000000000001*m0[1]*fOther[17]*vrelCX+0.4472135954999579*m0[7]*fOther[16]*vrelCX+0.4391550328268399*fOther[6]*m0[9]*vrelCX+0.5*fOther[0]*m0[7]*vrelCX+0.4*fOther[6]*m0[6]*vrelCX+0.447213595499958*m0[2]*fOther[6]*vrelCX+0.5000000000000001*fOther[1]*m0[5]*vrelCX+0.447213595499958*fOther[2]*m0[3]*vrelCX; 
  prodCX[33] = 0.31943828249997*m0[4]*fOther[33]*vrelCX+0.5*m0[0]*fOther[33]*vrelCX+0.447213595499958*m0[3]*fOther[21]*vrelCX+0.4391550328268399*fOther[7]*m0[8]*vrelCX+0.5*m0[6]*fOther[8]*vrelCX+0.447213595499958*m0[1]*fOther[7]*vrelCX+0.5000000000000001*fOther[3]*m0[4]*vrelCX; 
  prodCX[34] = 0.31943828249997*m0[5]*fOther[34]*vrelCX+0.5*m0[0]*fOther[34]*vrelCX+0.447213595499958*m0[3]*fOther[21]*vrelCX+0.4391550328268399*fOther[8]*m0[9]*vrelCX+0.447213595499958*m0[2]*fOther[8]*vrelCX+0.5*fOther[7]*m0[7]*vrelCX+0.5000000000000001*fOther[3]*m0[5]*vrelCX; 
  prodCX[35] = 0.5*m0[3]*fOther[36]*vrelCX+0.4472135954999579*m0[4]*fOther[35]*vrelCX+0.5*m0[0]*fOther[35]*vrelCX+0.5000000000000001*m0[1]*fOther[18]*vrelCX; 
  prodCX[36] = 0.4472135954999579*m0[5]*fOther[36]*vrelCX+0.5*m0[0]*fOther[36]*vrelCX+0.5*m0[3]*fOther[35]*vrelCX+0.5000000000000001*m0[2]*fOther[18]*vrelCX; 
  prodCX[37] = 0.31943828249997*m0[4]*fOther[37]*vrelCX+0.5*m0[0]*fOther[37]*vrelCX+0.447213595499958*m0[3]*fOther[22]*vrelCX+0.5*m0[6]*fOther[10]*vrelCX+0.4391550328268399*m0[8]*fOther[9]*vrelCX+0.447213595499958*m0[1]*fOther[9]*vrelCX+0.5000000000000001*fOther[4]*m0[4]*vrelCX; 
  prodCX[38] = 0.31943828249997*m0[5]*fOther[38]*vrelCX+0.5*m0[0]*fOther[38]*vrelCX+0.447213595499958*m0[3]*fOther[22]*vrelCX+0.4391550328268399*m0[9]*fOther[10]*vrelCX+0.447213595499958*m0[2]*fOther[10]*vrelCX+0.5*m0[7]*fOther[9]*vrelCX+0.5000000000000001*fOther[4]*m0[5]*vrelCX; 
  prodCX[39] = 0.5*m0[0]*fOther[39]*vrelCX; 
  prodCX[40] = 0.5*m0[3]*fOther[41]*vrelCX+0.4472135954999579*m0[4]*fOther[40]*vrelCX+0.5*m0[0]*fOther[40]*vrelCX+0.5000000000000001*m0[1]*fOther[19]*vrelCX; 
  prodCX[41] = 0.4472135954999579*m0[5]*fOther[41]*vrelCX+0.5*m0[0]*fOther[41]*vrelCX+0.5*m0[3]*fOther[40]*vrelCX+0.5000000000000001*m0[2]*fOther[19]*vrelCX; 
  prodCX[42] = 0.5*m0[0]*fOther[42]*vrelCX; 
  prodCX[43] = 0.31943828249997*m0[4]*fOther[43]*vrelCX+0.5*m0[0]*fOther[43]*vrelCX+0.447213595499958*m0[3]*fOther[25]*vrelCX+0.5*m0[6]*fOther[13]*vrelCX+0.4391550328268399*m0[8]*fOther[12]*vrelCX+0.447213595499958*m0[1]*fOther[12]*vrelCX+0.5000000000000001*m0[4]*fOther[5]*vrelCX; 
  prodCX[44] = 0.31943828249997*m0[5]*fOther[44]*vrelCX+0.5*m0[0]*fOther[44]*vrelCX+0.447213595499958*m0[3]*fOther[25]*vrelCX+0.4391550328268399*m0[9]*fOther[13]*vrelCX+0.447213595499958*m0[2]*fOther[13]*vrelCX+0.5*m0[7]*fOther[12]*vrelCX+0.5000000000000001*fOther[5]*m0[5]*vrelCX; 
  prodCX[45] = 0.5*m0[0]*fOther[45]*vrelCX; 
  prodCX[46] = 0.5*m0[0]*fOther[46]*vrelCX; 
  prodCX[47] = 0.5*m0[3]*fOther[48]*vrelCX+0.4472135954999579*m0[4]*fOther[47]*vrelCX+0.5*m0[0]*fOther[47]*vrelCX+0.5000000000000001*m0[1]*fOther[20]*vrelCX; 
  prodCX[48] = 0.4472135954999579*m0[5]*fOther[48]*vrelCX+0.5*m0[0]*fOther[48]*vrelCX+0.5*m0[3]*fOther[47]*vrelCX+0.5000000000000001*m0[2]*fOther[20]*vrelCX; 
  prodCX[49] = 0.5*m0[0]*fOther[49]*vrelCX; 
  prodCX[50] = 0.5*m0[0]*fOther[50]*vrelCX; 
  prodCX[51] = 0.2981423969999719*m0[4]*fOther[51]*vrelCX+0.5*m0[0]*fOther[51]*vrelCX+0.4391550328268399*m0[3]*fOther[31]*vrelCX+0.2981423969999719*m0[8]*fOther[16]*vrelCX+0.4391550328268398*m0[1]*fOther[16]*vrelCX+0.5*fOther[0]*m0[8]*vrelCX+0.4391550328268399*fOther[6]*m0[6]*vrelCX+0.4391550328268398*fOther[1]*m0[4]*vrelCX; 
  prodCX[52] = 0.2981423969999719*m0[5]*fOther[52]*vrelCX+0.5*m0[0]*fOther[52]*vrelCX+0.4391550328268399*m0[3]*fOther[32]*vrelCX+0.2981423969999719*m0[9]*fOther[17]*vrelCX+0.4391550328268398*m0[2]*fOther[17]*vrelCX+0.5*fOther[0]*m0[9]*vrelCX+0.4391550328268399*fOther[6]*m0[7]*vrelCX+0.4391550328268398*fOther[2]*m0[5]*vrelCX; 
  prodCX[53] = 0.5*m0[0]*fOther[53]*vrelCX; 
  prodCX[54] = 0.5*m0[0]*fOther[54]*vrelCX; 
  prodCX[55] = 0.5*m0[0]*fOther[55]*vrelCX; 
 
} 
