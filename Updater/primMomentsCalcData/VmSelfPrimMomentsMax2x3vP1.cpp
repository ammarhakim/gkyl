#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmSelfPrimMoments2x3vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(1.732050807568877*(m0[2]+m0[1])-1.0*m0[0]) < 0) cellAvg = true; 
  if (-0.5*(1.732050807568877*m0[2]-1.732050807568877*m0[1]-1.0*m0[0]) < 0) cellAvg = true; 
  if (0.5*(1.732050807568877*m0[2]-1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.5*(1.732050807568877*(m0[2]+m0[1])+m0[0]) < 0) cellAvg = true; 
 
  double m0r[3]; 
  double m1r[9]; 
  double cMr[9]; 
  double cEr[3]; 
  double m0Sr[3]; 
  double m1Sr[9]; 
  double m2Sr[3]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    m1r[3] = m1[3]; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    cMr[3] = cM[3]; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    m1r[6] = m1[6]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    cMr[6] = cM[6]; 
    cMr[7] = 0.0; 
    cMr[8] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = 0.0; 
    m1Sr[5] = 0.0; 
    m1Sr[6] = m1S[6]; 
    m1Sr[7] = 0.0; 
    m1Sr[8] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
    m2Sr[2] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m1r[6] = m1[6]; 
    m1r[7] = m1[7]; 
    m1r[8] = m1[8]; 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cMr[3] = cM[3]; 
    cMr[4] = cM[4]; 
    cMr[5] = cM[5]; 
    cMr[6] = cM[6]; 
    cMr[7] = cM[7]; 
    cMr[8] = cM[8]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = m1S[5]; 
    m1Sr[6] = m1S[6]; 
    m1Sr[7] = m1S[7]; 
    m1Sr[8] = m1S[8]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
    m2Sr[2] = m2S[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(12,12); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.5*m0r[0]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,2) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,9) = -0.5*cMr[0]; 
  data->AEM_S(0,10) = -0.5*cMr[1]; 
  data->AEM_S(0,11) = -0.5*cMr[2]; 
  data->AEM_S(1,9) = -0.5*cMr[1]; 
  data->AEM_S(1,10) = -0.5*cMr[0]; 
  data->AEM_S(2,9) = -0.5*cMr[2]; 
  data->AEM_S(2,11) = -0.5*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(9,0) = 0.5*m1Sr[0]; 
  data->AEM_S(9,1) = 0.5*m1Sr[1]; 
  data->AEM_S(9,2) = 0.5*m1Sr[2]; 
  data->AEM_S(10,0) = 0.5*m1Sr[1]; 
  data->AEM_S(10,1) = 0.5*m1Sr[0]; 
  data->AEM_S(11,0) = 0.5*m1Sr[2]; 
  data->AEM_S(11,2) = 0.5*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(3,3) = 0.5*m0r[0]; 
  data->AEM_S(3,4) = 0.5*m0r[1]; 
  data->AEM_S(3,5) = 0.5*m0r[2]; 
  data->AEM_S(4,3) = 0.5*m0r[1]; 
  data->AEM_S(4,4) = 0.5*m0r[0]; 
  data->AEM_S(5,3) = 0.5*m0r[2]; 
  data->AEM_S(5,5) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(3,9) = -0.5*cMr[3]; 
  data->AEM_S(3,10) = -0.5*cMr[4]; 
  data->AEM_S(3,11) = -0.5*cMr[5]; 
  data->AEM_S(4,9) = -0.5*cMr[4]; 
  data->AEM_S(4,10) = -0.5*cMr[3]; 
  data->AEM_S(5,9) = -0.5*cMr[5]; 
  data->AEM_S(5,11) = -0.5*cMr[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(9,3) = 0.5*m1Sr[3]; 
  data->AEM_S(9,4) = 0.5*m1Sr[4]; 
  data->AEM_S(9,5) = 0.5*m1Sr[5]; 
  data->AEM_S(10,3) = 0.5*m1Sr[4]; 
  data->AEM_S(10,4) = 0.5*m1Sr[3]; 
  data->AEM_S(11,3) = 0.5*m1Sr[5]; 
  data->AEM_S(11,5) = 0.5*m1Sr[3]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(6,6) = 0.5*m0r[0]; 
  data->AEM_S(6,7) = 0.5*m0r[1]; 
  data->AEM_S(6,8) = 0.5*m0r[2]; 
  data->AEM_S(7,6) = 0.5*m0r[1]; 
  data->AEM_S(7,7) = 0.5*m0r[0]; 
  data->AEM_S(8,6) = 0.5*m0r[2]; 
  data->AEM_S(8,8) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(6,9) = -0.5*cMr[6]; 
  data->AEM_S(6,10) = -0.5*cMr[7]; 
  data->AEM_S(6,11) = -0.5*cMr[8]; 
  data->AEM_S(7,9) = -0.5*cMr[7]; 
  data->AEM_S(7,10) = -0.5*cMr[6]; 
  data->AEM_S(8,9) = -0.5*cMr[8]; 
  data->AEM_S(8,11) = -0.5*cMr[6]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(9,6) = 0.5*m1Sr[6]; 
  data->AEM_S(9,7) = 0.5*m1Sr[7]; 
  data->AEM_S(9,8) = 0.5*m1Sr[8]; 
  data->AEM_S(10,6) = 0.5*m1Sr[7]; 
  data->AEM_S(10,7) = 0.5*m1Sr[6]; 
  data->AEM_S(11,6) = 0.5*m1Sr[8]; 
  data->AEM_S(11,8) = 0.5*m1Sr[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(9,9) = 0.5*m0Sr[0]-0.5*cEr[0]; 
  data->AEM_S(9,10) = 0.5*m0Sr[1]-0.5*cEr[1]; 
  data->AEM_S(9,11) = 0.5*m0Sr[2]-0.5*cEr[2]; 
  data->AEM_S(10,9) = 0.5*m0Sr[1]-0.5*cEr[1]; 
  data->AEM_S(10,10) = 0.5*m0Sr[0]-0.5*cEr[0]; 
  data->AEM_S(11,9) = 0.5*m0Sr[2]-0.5*cEr[2]; 
  data->AEM_S(11,11) = 0.5*m0Sr[0]-0.5*cEr[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,6>(0,3).setZero(); 
  data->AEM_S.block<6,3>(3,0).setZero(); 
  data->AEM_S.block<3,3>(3,6).setZero(); 
  data->AEM_S.block<3,3>(6,3).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m2Sr[0],m2Sr[1],m2Sr[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,9,1) = data->u_S.segment<9>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(9); 
 
} 
 
