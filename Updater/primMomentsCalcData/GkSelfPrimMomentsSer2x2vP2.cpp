#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkSelfPrimMoments2x2vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(3.872983346207417*(m0[7]+m0[6])-2.23606797749979*(m0[5]+m0[4])-3.0*m0[3]+1.732050807568877*(m0[2]+m0[1])-1.0*m0[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*m0[7]-3.872983346207417*m0[6]+2.23606797749979*(m0[5]+m0[4])-3.0*m0[3]-1.732050807568877*m0[2]+1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (-0.5*(3.872983346207417*m0[7]-3.872983346207417*m0[6]-2.23606797749979*(m0[5]+m0[4])+3.0*m0[3]-1.732050807568877*m0[2]+1.732050807568877*m0[1]-1.0*m0[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*(m0[7]+m0[6])+2.23606797749979*(m0[5]+m0[4])+3.0*m0[3]+1.732050807568877*(m0[2]+m0[1])+m0[0]) < 0) cellAvg = true; 
  if (-0.5*(3.872983346207417*(m2[7]+m2[6])-2.23606797749979*(m2[5]+m2[4])-3.0*m2[3]+1.732050807568877*(m2[2]+m2[1])-1.0*m2[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*m2[7]-3.872983346207417*m2[6]+2.23606797749979*(m2[5]+m2[4])-3.0*m2[3]-1.732050807568877*m2[2]+1.732050807568877*m2[1]+m2[0]) < 0) cellAvg = true; 
  if (-0.5*(3.872983346207417*m2[7]-3.872983346207417*m2[6]-2.23606797749979*(m2[5]+m2[4])+3.0*m2[3]-1.732050807568877*m2[2]+1.732050807568877*m2[1]-1.0*m2[0]) < 0) cellAvg = true; 
  if (0.5*(3.872983346207417*(m2[7]+m2[6])+2.23606797749979*(m2[5]+m2[4])+3.0*m2[3]+1.732050807568877*(m2[2]+m2[1])+m2[0]) < 0) cellAvg = true; 
 
  double m0r[8]; 
  double m1r[8]; 
  double cMr[8]; 
  double cEr[8]; 
  double m2r[8]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m0r[6] = 0.0; 
    m0r[7] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    cEr[4] = 0.0; 
    cEr[5] = 0.0; 
    cEr[6] = 0.0; 
    cEr[7] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
    m2r[6] = 0.0; 
    m2r[7] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
    m0r[6] = m0[6]; 
    m0r[7] = m0[7]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m1r[6] = m1[6]; 
    m1r[7] = m1[7]; 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cMr[3] = cM[3]; 
    cMr[4] = cM[4]; 
    cMr[5] = cM[5]; 
    cMr[6] = cM[6]; 
    cMr[7] = cM[7]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    cEr[3] = cE[3]; 
    cEr[4] = cE[4]; 
    cEr[5] = cE[5]; 
    cEr[6] = cE[6]; 
    cEr[7] = cE[7]; 
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
    m2r[6] = m2[6]; 
    m2r[7] = m2[7]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(16,16); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(0,3) = 0.5*m0r[3]; 
  data->AEM_S(0,4) = 0.5*m0r[4]; 
  data->AEM_S(0,5) = 0.5*m0r[5]; 
  data->AEM_S(0,6) = 0.5*m0r[6]; 
  data->AEM_S(0,7) = 0.5*m0r[7]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(1,2) = 0.5*m0r[3]; 
  data->AEM_S(1,3) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(1,4) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(1,5) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(1,6) = 0.447213595499958*m0r[3]; 
  data->AEM_S(1,7) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,1) = 0.5*m0r[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(2,3) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(2,4) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(2,5) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(2,6) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(2,7) = 0.447213595499958*m0r[3]; 
  data->AEM_S(3,0) = 0.5*m0r[3]; 
  data->AEM_S(3,1) = 0.447213595499958*m0r[6]+0.5*m0r[2]; 
  data->AEM_S(3,2) = 0.447213595499958*m0r[7]+0.5*m0r[1]; 
  data->AEM_S(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(3,4) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(3,5) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(3,6) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(3,7) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(4,0) = 0.5*m0r[4]; 
  data->AEM_S(4,1) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(4,2) = 0.5000000000000001*m0r[6]; 
  data->AEM_S(4,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(4,6) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(4,7) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(5,0) = 0.5*m0r[5]; 
  data->AEM_S(5,1) = 0.5000000000000001*m0r[7]; 
  data->AEM_S(5,2) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(5,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(5,6) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(5,7) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(6,0) = 0.5*m0r[6]; 
  data->AEM_S(6,1) = 0.447213595499958*m0r[3]; 
  data->AEM_S(6,2) = 0.5000000000000001*m0r[4]; 
  data->AEM_S(6,3) = 0.4*m0r[7]+0.447213595499958*m0r[1]; 
  data->AEM_S(6,4) = 0.31943828249997*m0r[6]+0.5000000000000001*m0r[2]; 
  data->AEM_S(6,5) = 0.4472135954999579*m0r[6]; 
  data->AEM_S(6,6) = 0.4472135954999579*m0r[5]+0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(6,7) = 0.4*m0r[3]; 
  data->AEM_S(7,0) = 0.5*m0r[7]; 
  data->AEM_S(7,1) = 0.5000000000000001*m0r[5]; 
  data->AEM_S(7,2) = 0.447213595499958*m0r[3]; 
  data->AEM_S(7,3) = 0.4*m0r[6]+0.447213595499958*m0r[2]; 
  data->AEM_S(7,4) = 0.4472135954999579*m0r[7]; 
  data->AEM_S(7,5) = 0.31943828249997*m0r[7]+0.5000000000000001*m0r[1]; 
  data->AEM_S(7,6) = 0.4*m0r[3]; 
  data->AEM_S(7,7) = 0.31943828249997*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,8) = -0.5*cMr[0]; 
  data->AEM_S(0,9) = -0.5*cMr[1]; 
  data->AEM_S(0,10) = -0.5*cMr[2]; 
  data->AEM_S(0,11) = -0.5*cMr[3]; 
  data->AEM_S(0,12) = -0.5*cMr[4]; 
  data->AEM_S(0,13) = -0.5*cMr[5]; 
  data->AEM_S(0,14) = -0.5*cMr[6]; 
  data->AEM_S(0,15) = -0.5*cMr[7]; 
  data->AEM_S(1,8) = -0.5*cMr[1]; 
  data->AEM_S(1,9) = (-0.4472135954999579*cMr[4])-0.5*cMr[0]; 
  data->AEM_S(1,10) = -0.5*cMr[3]; 
  data->AEM_S(1,11) = (-0.447213595499958*cMr[6])-0.5*cMr[2]; 
  data->AEM_S(1,12) = -0.4472135954999579*cMr[1]; 
  data->AEM_S(1,13) = -0.5000000000000001*cMr[7]; 
  data->AEM_S(1,14) = -0.447213595499958*cMr[3]; 
  data->AEM_S(1,15) = -0.5000000000000001*cMr[5]; 
  data->AEM_S(2,8) = -0.5*cMr[2]; 
  data->AEM_S(2,9) = -0.5*cMr[3]; 
  data->AEM_S(2,10) = (-0.4472135954999579*cMr[5])-0.5*cMr[0]; 
  data->AEM_S(2,11) = (-0.447213595499958*cMr[7])-0.5*cMr[1]; 
  data->AEM_S(2,12) = -0.5000000000000001*cMr[6]; 
  data->AEM_S(2,13) = -0.4472135954999579*cMr[2]; 
  data->AEM_S(2,14) = -0.5000000000000001*cMr[4]; 
  data->AEM_S(2,15) = -0.447213595499958*cMr[3]; 
  data->AEM_S(3,8) = -0.5*cMr[3]; 
  data->AEM_S(3,9) = (-0.447213595499958*cMr[6])-0.5*cMr[2]; 
  data->AEM_S(3,10) = (-0.447213595499958*cMr[7])-0.5*cMr[1]; 
  data->AEM_S(3,11) = (-0.4472135954999579*cMr[5])-0.4472135954999579*cMr[4]-0.5*cMr[0]; 
  data->AEM_S(3,12) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(3,13) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(3,14) = (-0.4*cMr[7])-0.447213595499958*cMr[1]; 
  data->AEM_S(3,15) = (-0.4*cMr[6])-0.447213595499958*cMr[2]; 
  data->AEM_S(4,8) = -0.5*cMr[4]; 
  data->AEM_S(4,9) = -0.4472135954999579*cMr[1]; 
  data->AEM_S(4,10) = -0.5000000000000001*cMr[6]; 
  data->AEM_S(4,11) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(4,12) = (-0.31943828249997*cMr[4])-0.5*cMr[0]; 
  data->AEM_S(4,14) = (-0.31943828249997*cMr[6])-0.5000000000000001*cMr[2]; 
  data->AEM_S(4,15) = -0.4472135954999579*cMr[7]; 
  data->AEM_S(5,8) = -0.5*cMr[5]; 
  data->AEM_S(5,9) = -0.5000000000000001*cMr[7]; 
  data->AEM_S(5,10) = -0.4472135954999579*cMr[2]; 
  data->AEM_S(5,11) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(5,13) = (-0.31943828249997*cMr[5])-0.5*cMr[0]; 
  data->AEM_S(5,14) = -0.4472135954999579*cMr[6]; 
  data->AEM_S(5,15) = (-0.31943828249997*cMr[7])-0.5000000000000001*cMr[1]; 
  data->AEM_S(6,8) = -0.5*cMr[6]; 
  data->AEM_S(6,9) = -0.447213595499958*cMr[3]; 
  data->AEM_S(6,10) = -0.5000000000000001*cMr[4]; 
  data->AEM_S(6,11) = (-0.4*cMr[7])-0.447213595499958*cMr[1]; 
  data->AEM_S(6,12) = (-0.31943828249997*cMr[6])-0.5000000000000001*cMr[2]; 
  data->AEM_S(6,13) = -0.4472135954999579*cMr[6]; 
  data->AEM_S(6,14) = (-0.4472135954999579*cMr[5])-0.31943828249997*cMr[4]-0.5*cMr[0]; 
  data->AEM_S(6,15) = -0.4*cMr[3]; 
  data->AEM_S(7,8) = -0.5*cMr[7]; 
  data->AEM_S(7,9) = -0.5000000000000001*cMr[5]; 
  data->AEM_S(7,10) = -0.447213595499958*cMr[3]; 
  data->AEM_S(7,11) = (-0.4*cMr[6])-0.447213595499958*cMr[2]; 
  data->AEM_S(7,12) = -0.4472135954999579*cMr[7]; 
  data->AEM_S(7,13) = (-0.31943828249997*cMr[7])-0.5000000000000001*cMr[1]; 
  data->AEM_S(7,14) = -0.4*cMr[3]; 
  data->AEM_S(7,15) = (-0.31943828249997*cMr[5])-0.4472135954999579*cMr[4]-0.5*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(8,0) = 0.5*m1r[0]; 
  data->AEM_S(8,1) = 0.5*m1r[1]; 
  data->AEM_S(8,2) = 0.5*m1r[2]; 
  data->AEM_S(8,3) = 0.5*m1r[3]; 
  data->AEM_S(8,4) = 0.5*m1r[4]; 
  data->AEM_S(8,5) = 0.5*m1r[5]; 
  data->AEM_S(8,6) = 0.5*m1r[6]; 
  data->AEM_S(8,7) = 0.5*m1r[7]; 
  data->AEM_S(9,0) = 0.5*m1r[1]; 
  data->AEM_S(9,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(9,2) = 0.5*m1r[3]; 
  data->AEM_S(9,3) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  data->AEM_S(9,4) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(9,5) = 0.5000000000000001*m1r[7]; 
  data->AEM_S(9,6) = 0.447213595499958*m1r[3]; 
  data->AEM_S(9,7) = 0.5000000000000001*m1r[5]; 
  data->AEM_S(10,0) = 0.5*m1r[2]; 
  data->AEM_S(10,1) = 0.5*m1r[3]; 
  data->AEM_S(10,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  data->AEM_S(10,3) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  data->AEM_S(10,4) = 0.5000000000000001*m1r[6]; 
  data->AEM_S(10,5) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(10,6) = 0.5000000000000001*m1r[4]; 
  data->AEM_S(10,7) = 0.447213595499958*m1r[3]; 
  data->AEM_S(11,0) = 0.5*m1r[3]; 
  data->AEM_S(11,1) = 0.447213595499958*m1r[6]+0.5*m1r[2]; 
  data->AEM_S(11,2) = 0.447213595499958*m1r[7]+0.5*m1r[1]; 
  data->AEM_S(11,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(11,4) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(11,5) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(11,6) = 0.4*m1r[7]+0.447213595499958*m1r[1]; 
  data->AEM_S(11,7) = 0.4*m1r[6]+0.447213595499958*m1r[2]; 
  data->AEM_S(12,0) = 0.5*m1r[4]; 
  data->AEM_S(12,1) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(12,2) = 0.5000000000000001*m1r[6]; 
  data->AEM_S(12,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(12,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(12,6) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  data->AEM_S(12,7) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(13,0) = 0.5*m1r[5]; 
  data->AEM_S(13,1) = 0.5000000000000001*m1r[7]; 
  data->AEM_S(13,2) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(13,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(13,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
  data->AEM_S(13,6) = 0.4472135954999579*m1r[6]; 
  data->AEM_S(13,7) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  data->AEM_S(14,0) = 0.5*m1r[6]; 
  data->AEM_S(14,1) = 0.447213595499958*m1r[3]; 
  data->AEM_S(14,2) = 0.5000000000000001*m1r[4]; 
  data->AEM_S(14,3) = 0.4*m1r[7]+0.447213595499958*m1r[1]; 
  data->AEM_S(14,4) = 0.31943828249997*m1r[6]+0.5000000000000001*m1r[2]; 
  data->AEM_S(14,5) = 0.4472135954999579*m1r[6]; 
  data->AEM_S(14,6) = 0.4472135954999579*m1r[5]+0.31943828249997*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(14,7) = 0.4*m1r[3]; 
  data->AEM_S(15,0) = 0.5*m1r[7]; 
  data->AEM_S(15,1) = 0.5000000000000001*m1r[5]; 
  data->AEM_S(15,2) = 0.447213595499958*m1r[3]; 
  data->AEM_S(15,3) = 0.4*m1r[6]+0.447213595499958*m1r[2]; 
  data->AEM_S(15,4) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(15,5) = 0.31943828249997*m1r[7]+0.5000000000000001*m1r[1]; 
  data->AEM_S(15,6) = 0.4*m1r[3]; 
  data->AEM_S(15,7) = 0.31943828249997*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(8,8) = 1.5*m0r[0]-0.5*cEr[0]; 
  data->AEM_S(8,9) = 1.5*m0r[1]-0.5*cEr[1]; 
  data->AEM_S(8,10) = 1.5*m0r[2]-0.5*cEr[2]; 
  data->AEM_S(8,11) = 1.5*m0r[3]-0.5*cEr[3]; 
  data->AEM_S(8,12) = 1.5*m0r[4]-0.5*cEr[4]; 
  data->AEM_S(8,13) = 1.5*m0r[5]-0.5*cEr[5]; 
  data->AEM_S(8,14) = 1.5*m0r[6]-0.5*cEr[6]; 
  data->AEM_S(8,15) = 1.5*m0r[7]-0.5*cEr[7]; 
  data->AEM_S(9,8) = 1.5*m0r[1]-0.5*cEr[1]; 
  data->AEM_S(9,9) = 1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]; 
  data->AEM_S(9,10) = 1.5*m0r[3]-0.5*cEr[3]; 
  data->AEM_S(9,11) = 1.341640786499874*m0r[6]-0.447213595499958*cEr[6]+1.5*m0r[2]-0.5*cEr[2]; 
  data->AEM_S(9,12) = 1.341640786499874*m0r[1]-0.4472135954999579*cEr[1]; 
  data->AEM_S(9,13) = 1.5*m0r[7]-0.5000000000000001*cEr[7]; 
  data->AEM_S(9,14) = 1.341640786499874*m0r[3]-0.447213595499958*cEr[3]; 
  data->AEM_S(9,15) = 1.5*m0r[5]-0.5000000000000001*cEr[5]; 
  data->AEM_S(10,8) = 1.5*m0r[2]-0.5*cEr[2]; 
  data->AEM_S(10,9) = 1.5*m0r[3]-0.5*cEr[3]; 
  data->AEM_S(10,10) = 1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+1.5*m0r[0]-0.5*cEr[0]; 
  data->AEM_S(10,11) = 1.341640786499874*m0r[7]-0.447213595499958*cEr[7]+1.5*m0r[1]-0.5*cEr[1]; 
  data->AEM_S(10,12) = 1.5*m0r[6]-0.5000000000000001*cEr[6]; 
  data->AEM_S(10,13) = 1.341640786499874*m0r[2]-0.4472135954999579*cEr[2]; 
  data->AEM_S(10,14) = 1.5*m0r[4]-0.5000000000000001*cEr[4]; 
  data->AEM_S(10,15) = 1.341640786499874*m0r[3]-0.447213595499958*cEr[3]; 
  data->AEM_S(11,8) = 1.5*m0r[3]-0.5*cEr[3]; 
  data->AEM_S(11,9) = 1.341640786499874*m0r[6]-0.447213595499958*cEr[6]+1.5*m0r[2]-0.5*cEr[2]; 
  data->AEM_S(11,10) = 1.341640786499874*m0r[7]-0.447213595499958*cEr[7]+1.5*m0r[1]-0.5*cEr[1]; 
  data->AEM_S(11,11) = 1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]; 
  data->AEM_S(11,12) = 1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(11,13) = 1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(11,14) = 1.2*m0r[7]-0.4*cEr[7]+1.341640786499874*m0r[1]-0.447213595499958*cEr[1]; 
  data->AEM_S(11,15) = 1.2*m0r[6]-0.4*cEr[6]+1.341640786499874*m0r[2]-0.447213595499958*cEr[2]; 
  data->AEM_S(12,8) = 1.5*m0r[4]-0.5*cEr[4]; 
  data->AEM_S(12,9) = 1.341640786499874*m0r[1]-0.4472135954999579*cEr[1]; 
  data->AEM_S(12,10) = 1.5*m0r[6]-0.5000000000000001*cEr[6]; 
  data->AEM_S(12,11) = 1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(12,12) = 0.9583148474999099*m0r[4]-0.31943828249997*cEr[4]+1.5*m0r[0]-0.5*cEr[0]; 
  data->AEM_S(12,14) = 0.9583148474999099*m0r[6]-0.31943828249997*cEr[6]+1.5*m0r[2]-0.5000000000000001*cEr[2]; 
  data->AEM_S(12,15) = 1.341640786499874*m0r[7]-0.4472135954999579*cEr[7]; 
  data->AEM_S(13,8) = 1.5*m0r[5]-0.5*cEr[5]; 
  data->AEM_S(13,9) = 1.5*m0r[7]-0.5000000000000001*cEr[7]; 
  data->AEM_S(13,10) = 1.341640786499874*m0r[2]-0.4472135954999579*cEr[2]; 
  data->AEM_S(13,11) = 1.341640786499874*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(13,13) = 0.9583148474999099*m0r[5]-0.31943828249997*cEr[5]+1.5*m0r[0]-0.5*cEr[0]; 
  data->AEM_S(13,14) = 1.341640786499874*m0r[6]-0.4472135954999579*cEr[6]; 
  data->AEM_S(13,15) = 0.9583148474999099*m0r[7]-0.31943828249997*cEr[7]+1.5*m0r[1]-0.5000000000000001*cEr[1]; 
  data->AEM_S(14,8) = 1.5*m0r[6]-0.5*cEr[6]; 
  data->AEM_S(14,9) = 1.341640786499874*m0r[3]-0.447213595499958*cEr[3]; 
  data->AEM_S(14,10) = 1.5*m0r[4]-0.5000000000000001*cEr[4]; 
  data->AEM_S(14,11) = 1.2*m0r[7]-0.4*cEr[7]+1.341640786499874*m0r[1]-0.447213595499958*cEr[1]; 
  data->AEM_S(14,12) = 0.9583148474999099*m0r[6]-0.31943828249997*cEr[6]+1.5*m0r[2]-0.5000000000000001*cEr[2]; 
  data->AEM_S(14,13) = 1.341640786499874*m0r[6]-0.4472135954999579*cEr[6]; 
  data->AEM_S(14,14) = 1.341640786499874*m0r[5]-0.4472135954999579*cEr[5]+0.9583148474999099*m0r[4]-0.31943828249997*cEr[4]+1.5*m0r[0]-0.5*cEr[0]; 
  data->AEM_S(14,15) = 1.2*m0r[3]-0.4*cEr[3]; 
  data->AEM_S(15,8) = 1.5*m0r[7]-0.5*cEr[7]; 
  data->AEM_S(15,9) = 1.5*m0r[5]-0.5000000000000001*cEr[5]; 
  data->AEM_S(15,10) = 1.341640786499874*m0r[3]-0.447213595499958*cEr[3]; 
  data->AEM_S(15,11) = 1.2*m0r[6]-0.4*cEr[6]+1.341640786499874*m0r[2]-0.447213595499958*cEr[2]; 
  data->AEM_S(15,12) = 1.341640786499874*m0r[7]-0.4472135954999579*cEr[7]; 
  data->AEM_S(15,13) = 0.9583148474999099*m0r[7]-0.31943828249997*cEr[7]+1.5*m0r[1]-0.5000000000000001*cEr[1]; 
  data->AEM_S(15,14) = 1.2*m0r[3]-0.4*cEr[3]; 
  data->AEM_S(15,15) = 0.9583148474999099*m0r[5]-0.31943828249997*cEr[5]+1.341640786499874*m0r[4]-0.4472135954999579*cEr[4]+1.5*m0r[0]-0.5*cEr[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,8,1) = data->u_S.segment<8>(0); 
 
  Eigen::Map<VectorXd>(vtSq,8,1) = data->u_S.segment<8>(8); 
 
} 
 