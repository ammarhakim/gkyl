#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmSelfPrimMoments2x2vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.5*(2.23606797749979*(m0[5]+m0[4])+3.0*m0[3]-1.732050807568877*(m0[2]+m0[1])+m0[0]) < 0) cellAvg = true; 
  if (0.5*(2.23606797749979*(m0[5]+m0[4])-3.0*m0[3]-1.732050807568877*m0[2]+1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.5*(2.23606797749979*(m0[5]+m0[4])-3.0*m0[3]+1.732050807568877*m0[2]-1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.5*(2.23606797749979*(m0[5]+m0[4])+3.0*m0[3]+1.732050807568877*(m0[2]+m0[1])+m0[0]) < 0) cellAvg = true; 
 
  double m0r[6]; 
  double m1r[12]; 
  double cMr[12]; 
  double cEr[6]; 
  double m2r[6]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    m1r[6] = m1[6]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    cMr[6] = cM[6]; 
    cMr[7] = 0.0; 
    cMr[8] = 0.0; 
    cMr[9] = 0.0; 
    cMr[10] = 0.0; 
    cMr[11] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    cEr[4] = 0.0; 
    cEr[5] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m1r[6] = m1[6]; 
    m1r[7] = m1[7]; 
    m1r[8] = m1[8]; 
    m1r[9] = m1[9]; 
    m1r[10] = m1[10]; 
    m1r[11] = m1[11]; 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cMr[3] = cM[3]; 
    cMr[4] = cM[4]; 
    cMr[5] = cM[5]; 
    cMr[6] = cM[6]; 
    cMr[7] = cM[7]; 
    cMr[8] = cM[8]; 
    cMr[9] = cM[9]; 
    cMr[10] = cM[10]; 
    cMr[11] = cM[11]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    cEr[3] = cE[3]; 
    cEr[4] = cE[4]; 
    cEr[5] = cE[5]; 
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(18,18); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(0,3) = 0.5*m0r[3]; 
  data->AEM_S(0,4) = 0.5*m0r[4]; 
  data->AEM_S(0,5) = 0.5*m0r[5]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(1,2) = 0.5*m0r[3]; 
  data->AEM_S(1,3) = 0.5*m0r[2]; 
  data->AEM_S(1,4) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,1) = 0.5*m0r[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(2,3) = 0.5*m0r[1]; 
  data->AEM_S(2,5) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(3,0) = 0.5*m0r[3]; 
  data->AEM_S(3,1) = 0.5*m0r[2]; 
  data->AEM_S(3,2) = 0.5*m0r[1]; 
  data->AEM_S(3,3) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(3,4) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(3,5) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(4,0) = 0.5*m0r[4]; 
  data->AEM_S(4,1) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(4,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(4,4) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(5,0) = 0.5*m0r[5]; 
  data->AEM_S(5,2) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(5,3) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(5,5) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,12) = -0.5*cMr[0]; 
  data->AEM_S(0,13) = -0.5*cMr[1]; 
  data->AEM_S(0,14) = -0.5*cMr[2]; 
  data->AEM_S(0,15) = -0.5*cMr[3]; 
  data->AEM_S(0,16) = -0.5*cMr[4]; 
  data->AEM_S(0,17) = -0.5*cMr[5]; 
  data->AEM_S(1,12) = -0.5*cMr[1]; 
  data->AEM_S(1,13) = (-0.4472135954999579*cMr[4])-0.5*cMr[0]; 
  data->AEM_S(1,14) = -0.5*cMr[3]; 
  data->AEM_S(1,15) = -0.5*cMr[2]; 
  data->AEM_S(1,16) = -0.4472135954999579*cMr[1]; 
  data->AEM_S(2,12) = -0.5*cMr[2]; 
  data->AEM_S(2,13) = -0.5*cMr[3]; 
  data->AEM_S(2,14) = (-0.4472135954999579*cMr[5])-0.5*cMr[0]; 
  data->AEM_S(2,15) = -0.5*cMr[1]; 
  data->AEM_S(2,17) = -0.4472135954999579*cMr[2]; 
  data->AEM_S(3,12) = -0.5*cMr[3]; 
  data->AEM_S(3,13) = -0.5*cMr[2]; 
  data->AEM_S(3,14) = -0.5*cMr[1]; 
  data->AEM_S(3,15) = (-0.4472135954999579*cMr[5])-0.4472135954999579*cMr[4]-0.5*cMr[0]; 
  data->AEM_S(3,16) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(3,17) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(4,12) = -0.5*cMr[4]; 
  data->AEM_S(4,13) = -0.4472135954999579*cMr[1]; 
  data->AEM_S(4,15) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(4,16) = (-0.31943828249997*cMr[4])-0.5*cMr[0]; 
  data->AEM_S(5,12) = -0.5*cMr[5]; 
  data->AEM_S(5,14) = -0.4472135954999579*cMr[2]; 
  data->AEM_S(5,15) = -0.4472135954999579*cMr[3]; 
  data->AEM_S(5,17) = (-0.31943828249997*cMr[5])-0.5*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(12,0) = 0.5*m1r[0]; 
  data->AEM_S(12,1) = 0.5*m1r[1]; 
  data->AEM_S(12,2) = 0.5*m1r[2]; 
  data->AEM_S(12,3) = 0.5*m1r[3]; 
  data->AEM_S(12,4) = 0.5*m1r[4]; 
  data->AEM_S(12,5) = 0.5*m1r[5]; 
  data->AEM_S(13,0) = 0.5*m1r[1]; 
  data->AEM_S(13,1) = 0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(13,2) = 0.5*m1r[3]; 
  data->AEM_S(13,3) = 0.5*m1r[2]; 
  data->AEM_S(13,4) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(14,0) = 0.5*m1r[2]; 
  data->AEM_S(14,1) = 0.5*m1r[3]; 
  data->AEM_S(14,2) = 0.4472135954999579*m1r[5]+0.5*m1r[0]; 
  data->AEM_S(14,3) = 0.5*m1r[1]; 
  data->AEM_S(14,5) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(15,0) = 0.5*m1r[3]; 
  data->AEM_S(15,1) = 0.5*m1r[2]; 
  data->AEM_S(15,2) = 0.5*m1r[1]; 
  data->AEM_S(15,3) = 0.4472135954999579*m1r[5]+0.4472135954999579*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(15,4) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(15,5) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(16,0) = 0.5*m1r[4]; 
  data->AEM_S(16,1) = 0.4472135954999579*m1r[1]; 
  data->AEM_S(16,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(16,4) = 0.31943828249997*m1r[4]+0.5*m1r[0]; 
  data->AEM_S(17,0) = 0.5*m1r[5]; 
  data->AEM_S(17,2) = 0.4472135954999579*m1r[2]; 
  data->AEM_S(17,3) = 0.4472135954999579*m1r[3]; 
  data->AEM_S(17,5) = 0.31943828249997*m1r[5]+0.5*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(6,6) = 0.5*m0r[0]; 
  data->AEM_S(6,7) = 0.5*m0r[1]; 
  data->AEM_S(6,8) = 0.5*m0r[2]; 
  data->AEM_S(6,9) = 0.5*m0r[3]; 
  data->AEM_S(6,10) = 0.5*m0r[4]; 
  data->AEM_S(6,11) = 0.5*m0r[5]; 
  data->AEM_S(7,6) = 0.5*m0r[1]; 
  data->AEM_S(7,7) = 0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(7,8) = 0.5*m0r[3]; 
  data->AEM_S(7,9) = 0.5*m0r[2]; 
  data->AEM_S(7,10) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(8,6) = 0.5*m0r[2]; 
  data->AEM_S(8,7) = 0.5*m0r[3]; 
  data->AEM_S(8,8) = 0.4472135954999579*m0r[5]+0.5*m0r[0]; 
  data->AEM_S(8,9) = 0.5*m0r[1]; 
  data->AEM_S(8,11) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(9,6) = 0.5*m0r[3]; 
  data->AEM_S(9,7) = 0.5*m0r[2]; 
  data->AEM_S(9,8) = 0.5*m0r[1]; 
  data->AEM_S(9,9) = 0.4472135954999579*m0r[5]+0.4472135954999579*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(9,10) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(9,11) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(10,6) = 0.5*m0r[4]; 
  data->AEM_S(10,7) = 0.4472135954999579*m0r[1]; 
  data->AEM_S(10,9) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(10,10) = 0.31943828249997*m0r[4]+0.5*m0r[0]; 
  data->AEM_S(11,6) = 0.5*m0r[5]; 
  data->AEM_S(11,8) = 0.4472135954999579*m0r[2]; 
  data->AEM_S(11,9) = 0.4472135954999579*m0r[3]; 
  data->AEM_S(11,11) = 0.31943828249997*m0r[5]+0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(6,12) = -0.5*cMr[6]; 
  data->AEM_S(6,13) = -0.5*cMr[7]; 
  data->AEM_S(6,14) = -0.5*cMr[8]; 
  data->AEM_S(6,15) = -0.5*cMr[9]; 
  data->AEM_S(6,16) = -0.5*cMr[10]; 
  data->AEM_S(6,17) = -0.5*cMr[11]; 
  data->AEM_S(7,12) = -0.5*cMr[7]; 
  data->AEM_S(7,13) = (-0.4472135954999579*cMr[10])-0.5*cMr[6]; 
  data->AEM_S(7,14) = -0.5*cMr[9]; 
  data->AEM_S(7,15) = -0.5*cMr[8]; 
  data->AEM_S(7,16) = -0.4472135954999579*cMr[7]; 
  data->AEM_S(8,12) = -0.5*cMr[8]; 
  data->AEM_S(8,13) = -0.5*cMr[9]; 
  data->AEM_S(8,14) = (-0.4472135954999579*cMr[11])-0.5*cMr[6]; 
  data->AEM_S(8,15) = -0.5*cMr[7]; 
  data->AEM_S(8,17) = -0.4472135954999579*cMr[8]; 
  data->AEM_S(9,12) = -0.5*cMr[9]; 
  data->AEM_S(9,13) = -0.5*cMr[8]; 
  data->AEM_S(9,14) = -0.5*cMr[7]; 
  data->AEM_S(9,15) = (-0.4472135954999579*cMr[11])-0.4472135954999579*cMr[10]-0.5*cMr[6]; 
  data->AEM_S(9,16) = -0.4472135954999579*cMr[9]; 
  data->AEM_S(9,17) = -0.4472135954999579*cMr[9]; 
  data->AEM_S(10,12) = -0.5*cMr[10]; 
  data->AEM_S(10,13) = -0.4472135954999579*cMr[7]; 
  data->AEM_S(10,15) = -0.4472135954999579*cMr[9]; 
  data->AEM_S(10,16) = (-0.31943828249997*cMr[10])-0.5*cMr[6]; 
  data->AEM_S(11,12) = -0.5*cMr[11]; 
  data->AEM_S(11,14) = -0.4472135954999579*cMr[8]; 
  data->AEM_S(11,15) = -0.4472135954999579*cMr[9]; 
  data->AEM_S(11,17) = (-0.31943828249997*cMr[11])-0.5*cMr[6]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(12,6) = 0.5*m1r[6]; 
  data->AEM_S(12,7) = 0.5*m1r[7]; 
  data->AEM_S(12,8) = 0.5*m1r[8]; 
  data->AEM_S(12,9) = 0.5*m1r[9]; 
  data->AEM_S(12,10) = 0.5*m1r[10]; 
  data->AEM_S(12,11) = 0.5*m1r[11]; 
  data->AEM_S(13,6) = 0.5*m1r[7]; 
  data->AEM_S(13,7) = 0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  data->AEM_S(13,8) = 0.5*m1r[9]; 
  data->AEM_S(13,9) = 0.5*m1r[8]; 
  data->AEM_S(13,10) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(14,6) = 0.5*m1r[8]; 
  data->AEM_S(14,7) = 0.5*m1r[9]; 
  data->AEM_S(14,8) = 0.4472135954999579*m1r[11]+0.5*m1r[6]; 
  data->AEM_S(14,9) = 0.5*m1r[7]; 
  data->AEM_S(14,11) = 0.4472135954999579*m1r[8]; 
  data->AEM_S(15,6) = 0.5*m1r[9]; 
  data->AEM_S(15,7) = 0.5*m1r[8]; 
  data->AEM_S(15,8) = 0.5*m1r[7]; 
  data->AEM_S(15,9) = 0.4472135954999579*m1r[11]+0.4472135954999579*m1r[10]+0.5*m1r[6]; 
  data->AEM_S(15,10) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(15,11) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(16,6) = 0.5*m1r[10]; 
  data->AEM_S(16,7) = 0.4472135954999579*m1r[7]; 
  data->AEM_S(16,9) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(16,10) = 0.31943828249997*m1r[10]+0.5*m1r[6]; 
  data->AEM_S(17,6) = 0.5*m1r[11]; 
  data->AEM_S(17,8) = 0.4472135954999579*m1r[8]; 
  data->AEM_S(17,9) = 0.4472135954999579*m1r[9]; 
  data->AEM_S(17,11) = 0.31943828249997*m1r[11]+0.5*m1r[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(12,12) = m0r[0]-0.5*cEr[0]; 
  data->AEM_S(12,13) = m0r[1]-0.5*cEr[1]; 
  data->AEM_S(12,14) = m0r[2]-0.5*cEr[2]; 
  data->AEM_S(12,15) = m0r[3]-0.5*cEr[3]; 
  data->AEM_S(12,16) = m0r[4]-0.5*cEr[4]; 
  data->AEM_S(12,17) = m0r[5]-0.5*cEr[5]; 
  data->AEM_S(13,12) = m0r[1]-0.5*cEr[1]; 
  data->AEM_S(13,13) = 0.8944271909999159*m0r[4]-0.4472135954999579*cEr[4]+m0r[0]-0.5*cEr[0]; 
  data->AEM_S(13,14) = m0r[3]-0.5*cEr[3]; 
  data->AEM_S(13,15) = m0r[2]-0.5*cEr[2]; 
  data->AEM_S(13,16) = 0.8944271909999159*m0r[1]-0.4472135954999579*cEr[1]; 
  data->AEM_S(14,12) = m0r[2]-0.5*cEr[2]; 
  data->AEM_S(14,13) = m0r[3]-0.5*cEr[3]; 
  data->AEM_S(14,14) = 0.8944271909999159*m0r[5]-0.4472135954999579*cEr[5]+m0r[0]-0.5*cEr[0]; 
  data->AEM_S(14,15) = m0r[1]-0.5*cEr[1]; 
  data->AEM_S(14,17) = 0.8944271909999159*m0r[2]-0.4472135954999579*cEr[2]; 
  data->AEM_S(15,12) = m0r[3]-0.5*cEr[3]; 
  data->AEM_S(15,13) = m0r[2]-0.5*cEr[2]; 
  data->AEM_S(15,14) = m0r[1]-0.5*cEr[1]; 
  data->AEM_S(15,15) = 0.8944271909999159*m0r[5]-0.4472135954999579*cEr[5]+0.8944271909999159*m0r[4]-0.4472135954999579*cEr[4]+m0r[0]-0.5*cEr[0]; 
  data->AEM_S(15,16) = 0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(15,17) = 0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(16,12) = m0r[4]-0.5*cEr[4]; 
  data->AEM_S(16,13) = 0.8944271909999159*m0r[1]-0.4472135954999579*cEr[1]; 
  data->AEM_S(16,15) = 0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(16,16) = 0.6388765649999399*m0r[4]-0.31943828249997*cEr[4]+m0r[0]-0.5*cEr[0]; 
  data->AEM_S(17,12) = m0r[5]-0.5*cEr[5]; 
  data->AEM_S(17,14) = 0.8944271909999159*m0r[2]-0.4472135954999579*cEr[2]; 
  data->AEM_S(17,15) = 0.8944271909999159*m0r[3]-0.4472135954999579*cEr[3]; 
  data->AEM_S(17,17) = 0.6388765649999399*m0r[5]-0.31943828249997*cEr[5]+m0r[0]-0.5*cEr[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<6,6>(0,6).setZero(); 
  data->AEM_S.block<6,6>(6,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,12,1) = data->u_S.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSq,6,1) = data->u_S.segment<6>(12); 
 
} 
 
