#include <RelativeVelocityModDecl.h> 
#include <math.h> 
void VmProdCXcellAvMax2x2v_P1(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = pow(w[3],2)/vtSqAv-(1.0*u[3]*w[3])/vtSqAv+(0.25*pow(u[3],2))/vtSqAv+pow(w[2],2)/vtSqAv-(1.0*u[0]*w[2])/vtSqAv+(0.25*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.5*fOther[2]*m0[2]*vrelCX+0.5*fOther[1]*m0[1]*vrelCX+0.5*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.5*fOther[0]*m0[1]*vrelCX+0.5*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.5*fOther[0]*m0[2]*vrelCX+0.5*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.5*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.5*m0[0]*fOther[4]*vrelCX; 
 
} 
void VmProdCXcellAvMax2x2v_P2(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = (0.25*pow(u[6],2))/vtSqAv-(1.0*w[3]*u[6])/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv-(1.0*u[0]*w[2])/vtSqAv+(0.25*pow(u[0],2))/vtSqAv; 
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
void VmProdCXcellAvMax2x2v_P3(const double *w, const double *m0, const double *u, const double *vtSq, const double *fOther, double *prodCX) 
{ 
  // w[4]:   cell-center coordinates. 
  // m0:      density. 
  // u:       velocity. 
  // vtSq:    squared thermal speed, sqrt(T/m). 
  // fOther:    distribution function of other CX species. 
  // prodCX:  produce of v^*, m0, and f in Pauls CX model. 
 
  double vtSqAv = 0.5*vtSq[0]; 
  double xSqAv = (0.25*pow(u[10],2))/vtSqAv-(1.0*w[3]*u[10])/vtSqAv+pow(w[3],2)/vtSqAv+pow(w[2],2)/vtSqAv-(1.0*u[0]*w[2])/vtSqAv+(0.25*pow(u[0],2))/vtSqAv; 
  double vrelCX = 0.5641895835477563*sqrt(vtSqAv)*sqrt(3.141592653589793*xSqAv+4.0); 
 
  prodCX[0] = 0.5*m0[9]*fOther[32]*vrelCX+0.5*m0[8]*fOther[31]*vrelCX+0.5*m0[7]*fOther[20]*vrelCX+0.5*m0[6]*fOther[19]*vrelCX+0.5*m0[5]*fOther[12]*vrelCX+0.5*m0[4]*fOther[11]*vrelCX+0.5*m0[3]*fOther[5]*vrelCX+0.5*fOther[2]*m0[2]*vrelCX+0.5*fOther[1]*m0[1]*vrelCX+0.5*fOther[0]*m0[0]*vrelCX; 
  prodCX[1] = 0.4391550328268398*m0[4]*fOther[31]*vrelCX+0.5000000000000001*m0[5]*fOther[20]*vrelCX+0.447213595499958*m0[3]*fOther[19]*vrelCX+0.5000000000000001*m0[7]*fOther[12]*vrelCX+0.4391550328268398*m0[8]*fOther[11]*vrelCX+0.4472135954999579*m0[1]*fOther[11]*vrelCX+0.447213595499958*fOther[5]*m0[6]*vrelCX+0.5*m0[2]*fOther[5]*vrelCX+0.4472135954999579*fOther[1]*m0[4]*vrelCX+0.5*fOther[2]*m0[3]*vrelCX+0.5*fOther[0]*m0[1]*vrelCX+0.5*m0[0]*fOther[1]*vrelCX; 
  prodCX[2] = 0.4391550328268398*m0[5]*fOther[32]*vrelCX+0.447213595499958*m0[3]*fOther[20]*vrelCX+0.5000000000000001*m0[4]*fOther[19]*vrelCX+0.4391550328268398*m0[9]*fOther[12]*vrelCX+0.4472135954999579*m0[2]*fOther[12]*vrelCX+0.5000000000000001*m0[6]*fOther[11]*vrelCX+0.447213595499958*fOther[5]*m0[7]*vrelCX+0.4472135954999579*fOther[2]*m0[5]*vrelCX+0.5*m0[1]*fOther[5]*vrelCX+0.5*fOther[1]*m0[3]*vrelCX+0.5*fOther[0]*m0[2]*vrelCX+0.5*m0[0]*fOther[2]*vrelCX; 
  prodCX[3] = 0.5000000000000001*m0[5]*fOther[22]*vrelCX+0.5000000000000001*m0[4]*fOther[21]*vrelCX+0.5*m0[3]*fOther[15]*vrelCX+0.5*m0[2]*fOther[7]*vrelCX+0.5*m0[1]*fOther[6]*vrelCX+0.5*m0[0]*fOther[3]*vrelCX; 
  prodCX[4] = 0.5000000000000001*m0[5]*fOther[26]*vrelCX+0.5000000000000001*m0[4]*fOther[25]*vrelCX+0.5*m0[3]*fOther[16]*vrelCX+0.5*m0[2]*fOther[9]*vrelCX+0.5*m0[1]*fOther[8]*vrelCX+0.5*m0[0]*fOther[4]*vrelCX; 
  prodCX[5] = 0.4391550328268399*m0[7]*fOther[32]*vrelCX+0.4391550328268399*m0[6]*fOther[31]*vrelCX+0.4391550328268399*m0[9]*fOther[20]*vrelCX+0.4*m0[6]*fOther[20]*vrelCX+0.447213595499958*m0[2]*fOther[20]*vrelCX+0.4391550328268399*m0[8]*fOther[19]*vrelCX+0.4*m0[7]*fOther[19]*vrelCX+0.447213595499958*m0[1]*fOther[19]*vrelCX+0.4472135954999579*m0[3]*fOther[12]*vrelCX+0.4472135954999579*m0[3]*fOther[11]*vrelCX+0.447213595499958*fOther[2]*m0[7]*vrelCX+0.447213595499958*fOther[1]*m0[6]*vrelCX+0.4472135954999579*fOther[5]*m0[5]*vrelCX+0.4472135954999579*m0[4]*fOther[5]*vrelCX+0.5*m0[0]*fOther[5]*vrelCX+0.5*fOther[0]*m0[3]*vrelCX+0.5*fOther[1]*m0[2]*vrelCX+0.5*m0[1]*fOther[2]*vrelCX; 
  prodCX[6] = 0.5*m0[7]*fOther[22]*vrelCX+0.4391550328268399*m0[8]*fOther[21]*vrelCX+0.447213595499958*m0[1]*fOther[21]*vrelCX+0.447213595499958*m0[6]*fOther[15]*vrelCX+0.5*m0[2]*fOther[15]*vrelCX+0.5*m0[3]*fOther[7]*vrelCX+0.4472135954999579*m0[4]*fOther[6]*vrelCX+0.5*m0[0]*fOther[6]*vrelCX+0.5*m0[1]*fOther[3]*vrelCX; 
  prodCX[7] = 0.4391550328268399*m0[9]*fOther[22]*vrelCX+0.447213595499958*m0[2]*fOther[22]*vrelCX+0.5*m0[6]*fOther[21]*vrelCX+0.447213595499958*m0[7]*fOther[15]*vrelCX+0.5*m0[1]*fOther[15]*vrelCX+0.4472135954999579*m0[5]*fOther[7]*vrelCX+0.5*m0[0]*fOther[7]*vrelCX+0.5*m0[3]*fOther[6]*vrelCX+0.5*m0[2]*fOther[3]*vrelCX; 
  prodCX[8] = 0.5*m0[7]*fOther[26]*vrelCX+0.4391550328268399*m0[8]*fOther[25]*vrelCX+0.447213595499958*m0[1]*fOther[25]*vrelCX+0.447213595499958*m0[6]*fOther[16]*vrelCX+0.5*m0[2]*fOther[16]*vrelCX+0.5*m0[3]*fOther[9]*vrelCX+0.4472135954999579*m0[4]*fOther[8]*vrelCX+0.5*m0[0]*fOther[8]*vrelCX+0.5*m0[1]*fOther[4]*vrelCX; 
  prodCX[9] = 0.4391550328268399*m0[9]*fOther[26]*vrelCX+0.447213595499958*m0[2]*fOther[26]*vrelCX+0.5*m0[6]*fOther[25]*vrelCX+0.447213595499958*m0[7]*fOther[16]*vrelCX+0.5*m0[1]*fOther[16]*vrelCX+0.4472135954999579*m0[5]*fOther[9]*vrelCX+0.5*m0[0]*fOther[9]*vrelCX+0.5*m0[3]*fOther[8]*vrelCX+0.5*m0[2]*fOther[4]*vrelCX; 
  prodCX[10] = 0.5*m0[2]*fOther[18]*vrelCX+0.5*m0[1]*fOther[17]*vrelCX+0.5*m0[0]*fOther[10]*vrelCX; 
  prodCX[11] = 0.2981423969999719*m0[8]*fOther[31]*vrelCX+0.4391550328268398*m0[1]*fOther[31]*vrelCX+0.4472135954999579*m0[7]*fOther[20]*vrelCX+0.31943828249997*m0[6]*fOther[19]*vrelCX+0.5000000000000001*m0[2]*fOther[19]*vrelCX+0.31943828249997*m0[4]*fOther[11]*vrelCX+0.5*m0[0]*fOther[11]*vrelCX+0.4391550328268398*fOther[1]*m0[8]*vrelCX+0.5000000000000001*fOther[2]*m0[6]*vrelCX+0.4472135954999579*m0[3]*fOther[5]*vrelCX+0.5*fOther[0]*m0[4]*vrelCX+0.4472135954999579*fOther[1]*m0[1]*vrelCX; 
  prodCX[12] = 0.2981423969999719*m0[9]*fOther[32]*vrelCX+0.4391550328268398*m0[2]*fOther[32]*vrelCX+0.31943828249997*m0[7]*fOther[20]*vrelCX+0.5000000000000001*m0[1]*fOther[20]*vrelCX+0.4472135954999579*m0[6]*fOther[19]*vrelCX+0.31943828249997*m0[5]*fOther[12]*vrelCX+0.5*m0[0]*fOther[12]*vrelCX+0.4391550328268398*fOther[2]*m0[9]*vrelCX+0.5000000000000001*fOther[1]*m0[7]*vrelCX+0.5*fOther[0]*m0[5]*vrelCX+0.4472135954999579*m0[3]*fOther[5]*vrelCX+0.4472135954999579*fOther[2]*m0[2]*vrelCX; 
  prodCX[13] = 0.5000000000000001*m0[2]*fOther[24]*vrelCX+0.5000000000000001*m0[1]*fOther[23]*vrelCX+0.5*m0[0]*fOther[13]*vrelCX; 
  prodCX[14] = 0.5000000000000001*m0[2]*fOther[29]*vrelCX+0.5000000000000001*m0[1]*fOther[28]*vrelCX+0.5*m0[0]*fOther[14]*vrelCX; 
  prodCX[15] = 0.447213595499958*m0[3]*fOther[22]*vrelCX+0.447213595499958*m0[3]*fOther[21]*vrelCX+0.4472135954999579*m0[5]*fOther[15]*vrelCX+0.4472135954999579*m0[4]*fOther[15]*vrelCX+0.5*m0[0]*fOther[15]*vrelCX+0.447213595499958*fOther[7]*m0[7]*vrelCX+0.5*m0[1]*fOther[7]*vrelCX+0.447213595499958*fOther[6]*m0[6]*vrelCX+0.5*m0[2]*fOther[6]*vrelCX+0.5*fOther[3]*m0[3]*vrelCX; 
  prodCX[16] = 0.447213595499958*m0[3]*fOther[26]*vrelCX+0.447213595499958*m0[3]*fOther[25]*vrelCX+0.4472135954999579*m0[5]*fOther[16]*vrelCX+0.4472135954999579*m0[4]*fOther[16]*vrelCX+0.5*m0[0]*fOther[16]*vrelCX+0.447213595499958*m0[7]*fOther[9]*vrelCX+0.5*m0[1]*fOther[9]*vrelCX+0.447213595499958*m0[6]*fOther[8]*vrelCX+0.5*m0[2]*fOther[8]*vrelCX+0.5*m0[3]*fOther[4]*vrelCX; 
  prodCX[17] = 0.5*m0[3]*fOther[18]*vrelCX+0.4472135954999579*m0[4]*fOther[17]*vrelCX+0.5*m0[0]*fOther[17]*vrelCX+0.5*m0[1]*fOther[10]*vrelCX; 
  prodCX[18] = 0.4472135954999579*m0[5]*fOther[18]*vrelCX+0.5*m0[0]*fOther[18]*vrelCX+0.5*m0[3]*fOther[17]*vrelCX+0.5*m0[2]*fOther[10]*vrelCX; 
  prodCX[19] = 0.4391550328268399*m0[3]*fOther[31]*vrelCX+0.4*m0[3]*fOther[20]*vrelCX+0.4472135954999579*m0[5]*fOther[19]*vrelCX+0.31943828249997*m0[4]*fOther[19]*vrelCX+0.5*m0[0]*fOther[19]*vrelCX+0.4472135954999579*m0[6]*fOther[12]*vrelCX+0.31943828249997*m0[6]*fOther[11]*vrelCX+0.5000000000000001*m0[2]*fOther[11]*vrelCX+0.4391550328268399*fOther[5]*m0[8]*vrelCX+0.4*fOther[5]*m0[7]*vrelCX+0.5*fOther[0]*m0[6]*vrelCX+0.447213595499958*m0[1]*fOther[5]*vrelCX+0.5000000000000001*fOther[2]*m0[4]*vrelCX+0.447213595499958*fOther[1]*m0[3]*vrelCX; 
  prodCX[20] = 0.4391550328268399*m0[3]*fOther[32]*vrelCX+0.31943828249997*m0[5]*fOther[20]*vrelCX+0.4472135954999579*m0[4]*fOther[20]*vrelCX+0.5*m0[0]*fOther[20]*vrelCX+0.4*m0[3]*fOther[19]*vrelCX+0.31943828249997*m0[7]*fOther[12]*vrelCX+0.5000000000000001*m0[1]*fOther[12]*vrelCX+0.4472135954999579*m0[7]*fOther[11]*vrelCX+0.4391550328268399*fOther[5]*m0[9]*vrelCX+0.5*fOther[0]*m0[7]*vrelCX+0.4*fOther[5]*m0[6]*vrelCX+0.5000000000000001*fOther[1]*m0[5]*vrelCX+0.447213595499958*m0[2]*fOther[5]*vrelCX+0.447213595499958*fOther[2]*m0[3]*vrelCX; 
  prodCX[21] = 0.31943828249997*m0[4]*fOther[21]*vrelCX+0.5*m0[0]*fOther[21]*vrelCX+0.447213595499958*m0[3]*fOther[15]*vrelCX+0.4391550328268399*fOther[6]*m0[8]*vrelCX+0.5*m0[6]*fOther[7]*vrelCX+0.447213595499958*m0[1]*fOther[6]*vrelCX+0.5000000000000001*fOther[3]*m0[4]*vrelCX; 
  prodCX[22] = 0.31943828249997*m0[5]*fOther[22]*vrelCX+0.5*m0[0]*fOther[22]*vrelCX+0.447213595499958*m0[3]*fOther[15]*vrelCX+0.4391550328268399*fOther[7]*m0[9]*vrelCX+0.5*fOther[6]*m0[7]*vrelCX+0.447213595499958*m0[2]*fOther[7]*vrelCX+0.5000000000000001*fOther[3]*m0[5]*vrelCX; 
  prodCX[23] = 0.5*m0[3]*fOther[24]*vrelCX+0.4472135954999579*m0[4]*fOther[23]*vrelCX+0.5*m0[0]*fOther[23]*vrelCX+0.5000000000000001*m0[1]*fOther[13]*vrelCX; 
  prodCX[24] = 0.4472135954999579*m0[5]*fOther[24]*vrelCX+0.5*m0[0]*fOther[24]*vrelCX+0.5*m0[3]*fOther[23]*vrelCX+0.5000000000000001*m0[2]*fOther[13]*vrelCX; 
  prodCX[25] = 0.31943828249997*m0[4]*fOther[25]*vrelCX+0.5*m0[0]*fOther[25]*vrelCX+0.447213595499958*m0[3]*fOther[16]*vrelCX+0.5*m0[6]*fOther[9]*vrelCX+0.4391550328268399*fOther[8]*m0[8]*vrelCX+0.447213595499958*m0[1]*fOther[8]*vrelCX+0.5000000000000001*fOther[4]*m0[4]*vrelCX; 
  prodCX[26] = 0.31943828249997*m0[5]*fOther[26]*vrelCX+0.5*m0[0]*fOther[26]*vrelCX+0.447213595499958*m0[3]*fOther[16]*vrelCX+0.4391550328268399*fOther[9]*m0[9]*vrelCX+0.447213595499958*m0[2]*fOther[9]*vrelCX+0.5*m0[7]*fOther[8]*vrelCX+0.5000000000000001*fOther[4]*m0[5]*vrelCX; 
  prodCX[27] = 0.5*m0[0]*fOther[27]*vrelCX; 
  prodCX[28] = 0.5*m0[3]*fOther[29]*vrelCX+0.4472135954999579*m0[4]*fOther[28]*vrelCX+0.5*m0[0]*fOther[28]*vrelCX+0.5000000000000001*m0[1]*fOther[14]*vrelCX; 
  prodCX[29] = 0.4472135954999579*m0[5]*fOther[29]*vrelCX+0.5*m0[0]*fOther[29]*vrelCX+0.5*m0[3]*fOther[28]*vrelCX+0.5000000000000001*m0[2]*fOther[14]*vrelCX; 
  prodCX[30] = 0.5*m0[0]*fOther[30]*vrelCX; 
  prodCX[31] = 0.2981423969999719*m0[4]*fOther[31]*vrelCX+0.5*m0[0]*fOther[31]*vrelCX+0.4391550328268399*m0[3]*fOther[19]*vrelCX+0.2981423969999719*m0[8]*fOther[11]*vrelCX+0.4391550328268398*m0[1]*fOther[11]*vrelCX+0.5*fOther[0]*m0[8]*vrelCX+0.4391550328268399*fOther[5]*m0[6]*vrelCX+0.4391550328268398*fOther[1]*m0[4]*vrelCX; 
  prodCX[32] = 0.2981423969999719*m0[5]*fOther[32]*vrelCX+0.5*m0[0]*fOther[32]*vrelCX+0.4391550328268399*m0[3]*fOther[20]*vrelCX+0.2981423969999719*m0[9]*fOther[12]*vrelCX+0.4391550328268398*m0[2]*fOther[12]*vrelCX+0.5*fOther[0]*m0[9]*vrelCX+0.4391550328268399*fOther[5]*m0[7]*vrelCX+0.4391550328268398*fOther[2]*m0[5]*vrelCX; 
  prodCX[33] = 0.5*m0[0]*fOther[33]*vrelCX; 
  prodCX[34] = 0.5*m0[0]*fOther[34]*vrelCX; 
 
} 
