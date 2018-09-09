#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments2x3vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool m0Avg = false;
  if ((-0.8660254037844386*m0[2])-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if ((-0.8660254037844386*m0[2])-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if ((-0.8660254037844386*m0[2])+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if ((-0.8660254037844386*m0[2])+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
 
  double m0s[3]; 
  if (m0Avg) { 
    m0s[0] = m0[0]; 
    m0s[1] = 0.0; 
    m0s[2] = 0.0; 
  } else { 
    m0s[0] = m0[0]; 
    m0s[1] = m0[1]; 
    m0s[2] = m0[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(12,12); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(12);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(12);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0s[0]; 
  BigAEM(0,1) = 0.5*m0s[1]; 
  BigAEM(0,2) = 0.5*m0s[2]; 
  BigAEM(1,0) = 0.5*m0s[1]; 
  BigAEM(1,1) = 0.5*m0s[0]; 
  BigAEM(2,0) = 0.5*m0s[2]; 
  BigAEM(2,2) = 0.5*m0s[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,9) = -0.5*cM[0]; 
  BigAEM(0,10) = -0.5*cM[1]; 
  BigAEM(0,11) = -0.5*cM[2]; 
  BigAEM(1,9) = -0.5*cM[1]; 
  BigAEM(1,10) = -0.5*cM[0]; 
  BigAEM(2,9) = -0.5*cM[2]; 
  BigAEM(2,11) = -0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(9,0) = 0.5*m1[0]; 
  BigAEM(9,1) = 0.5*m1[1]; 
  BigAEM(9,2) = 0.5*m1[2]; 
  BigAEM(10,0) = 0.5*m1[1]; 
  BigAEM(10,1) = 0.5*m1[0]; 
  BigAEM(11,0) = 0.5*m1[2]; 
  BigAEM(11,2) = 0.5*m1[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(3,3) = 0.5*m0s[0]; 
  BigAEM(3,4) = 0.5*m0s[1]; 
  BigAEM(3,5) = 0.5*m0s[2]; 
  BigAEM(4,3) = 0.5*m0s[1]; 
  BigAEM(4,4) = 0.5*m0s[0]; 
  BigAEM(5,3) = 0.5*m0s[2]; 
  BigAEM(5,5) = 0.5*m0s[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(3,9) = -0.5*cM[3]; 
  BigAEM(3,10) = -0.5*cM[4]; 
  BigAEM(3,11) = -0.5*cM[5]; 
  BigAEM(4,9) = -0.5*cM[4]; 
  BigAEM(4,10) = -0.5*cM[3]; 
  BigAEM(5,9) = -0.5*cM[5]; 
  BigAEM(5,11) = -0.5*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(9,3) = 0.5*m1[3]; 
  BigAEM(9,4) = 0.5*m1[4]; 
  BigAEM(9,5) = 0.5*m1[5]; 
  BigAEM(10,3) = 0.5*m1[4]; 
  BigAEM(10,4) = 0.5*m1[3]; 
  BigAEM(11,3) = 0.5*m1[5]; 
  BigAEM(11,5) = 0.5*m1[3]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  BigAEM(6,6) = 0.5*m0s[0]; 
  BigAEM(6,7) = 0.5*m0s[1]; 
  BigAEM(6,8) = 0.5*m0s[2]; 
  BigAEM(7,6) = 0.5*m0s[1]; 
  BigAEM(7,7) = 0.5*m0s[0]; 
  BigAEM(8,6) = 0.5*m0s[2]; 
  BigAEM(8,8) = 0.5*m0s[0]; 
 
  // ....... Block from correction to uZ .......... // 
  BigAEM(6,9) = -0.5*cM[6]; 
  BigAEM(6,10) = -0.5*cM[7]; 
  BigAEM(6,11) = -0.5*cM[8]; 
  BigAEM(7,9) = -0.5*cM[7]; 
  BigAEM(7,10) = -0.5*cM[6]; 
  BigAEM(8,9) = -0.5*cM[8]; 
  BigAEM(8,11) = -0.5*cM[6]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  BigAEM(9,6) = 0.5*m1[6]; 
  BigAEM(9,7) = 0.5*m1[7]; 
  BigAEM(9,8) = 0.5*m1[8]; 
  BigAEM(10,6) = 0.5*m1[7]; 
  BigAEM(10,7) = 0.5*m1[6]; 
  BigAEM(11,6) = 0.5*m1[8]; 
  BigAEM(11,8) = 0.5*m1[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(9,9) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
  BigAEM(9,10) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(9,11) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(10,9) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(10,10) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
  BigAEM(11,9) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(11,11) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<3,6>(0,3).setZero(); 
  BigAEM.block<6,3>(3,0).setZero(); 
  BigAEM.block<3,3>(3,6).setZero(); 
  BigAEM.block<3,3>(6,3).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m1[2],m1[3],m1[4],m1[5],m1[6],m1[7],m1[8],m2[0],m2[1],m2[2]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,9,1) = xEV.segment<9>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV.segment<3>(9); 
 
} 
 
void SelfPrimMoments2x3vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool m0Avg = false;
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if (1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
 
  double m0s[6]; 
  if (m0Avg) { 
    m0s[0] = m0[0]; 
    m0s[1] = 0.0; 
    m0s[2] = 0.0; 
    m0s[3] = 0.0; 
    m0s[4] = 0.0; 
    m0s[5] = 0.0; 
  } else { 
    m0s[0] = m0[0]; 
    m0s[1] = m0[1]; 
    m0s[2] = m0[2]; 
    m0s[3] = m0[3]; 
    m0s[4] = m0[4]; 
    m0s[5] = m0[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(24,24); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(24);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(24);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0s[0]; 
  BigAEM(0,1) = 0.5*m0s[1]; 
  BigAEM(0,2) = 0.5*m0s[2]; 
  BigAEM(0,3) = 0.5*m0s[3]; 
  BigAEM(0,4) = 0.5*m0s[4]; 
  BigAEM(0,5) = 0.5*m0s[5]; 
  BigAEM(1,0) = 0.5*m0s[1]; 
  BigAEM(1,1) = 0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(1,2) = 0.5*m0s[3]; 
  BigAEM(1,3) = 0.5*m0s[2]; 
  BigAEM(1,4) = 0.4472135954999579*m0s[1]; 
  BigAEM(2,0) = 0.5*m0s[2]; 
  BigAEM(2,1) = 0.5*m0s[3]; 
  BigAEM(2,2) = 0.4472135954999579*m0s[5]+0.5*m0s[0]; 
  BigAEM(2,3) = 0.5*m0s[1]; 
  BigAEM(2,5) = 0.4472135954999579*m0s[2]; 
  BigAEM(3,0) = 0.5*m0s[3]; 
  BigAEM(3,1) = 0.5*m0s[2]; 
  BigAEM(3,2) = 0.5*m0s[1]; 
  BigAEM(3,3) = 0.4472135954999579*m0s[5]+0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(3,4) = 0.4472135954999579*m0s[3]; 
  BigAEM(3,5) = 0.4472135954999579*m0s[3]; 
  BigAEM(4,0) = 0.5*m0s[4]; 
  BigAEM(4,1) = 0.4472135954999579*m0s[1]; 
  BigAEM(4,3) = 0.4472135954999579*m0s[3]; 
  BigAEM(4,4) = 0.31943828249997*m0s[4]+0.5*m0s[0]; 
  BigAEM(5,0) = 0.5*m0s[5]; 
  BigAEM(5,2) = 0.4472135954999579*m0s[2]; 
  BigAEM(5,3) = 0.4472135954999579*m0s[3]; 
  BigAEM(5,5) = 0.31943828249997*m0s[5]+0.5*m0s[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,18) = -0.5*cM[0]; 
  BigAEM(0,19) = -0.5*cM[1]; 
  BigAEM(0,20) = -0.5*cM[2]; 
  BigAEM(0,21) = -0.5*cM[3]; 
  BigAEM(0,22) = -0.5*cM[4]; 
  BigAEM(0,23) = -0.5*cM[5]; 
  BigAEM(1,18) = -0.5*cM[1]; 
  BigAEM(1,19) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,20) = -0.5*cM[3]; 
  BigAEM(1,21) = -0.5*cM[2]; 
  BigAEM(1,22) = -0.4472135954999579*cM[1]; 
  BigAEM(2,18) = -0.5*cM[2]; 
  BigAEM(2,19) = -0.5*cM[3]; 
  BigAEM(2,20) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,21) = -0.5*cM[1]; 
  BigAEM(2,23) = -0.4472135954999579*cM[2]; 
  BigAEM(3,18) = -0.5*cM[3]; 
  BigAEM(3,19) = -0.5*cM[2]; 
  BigAEM(3,20) = -0.5*cM[1]; 
  BigAEM(3,21) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,22) = -0.4472135954999579*cM[3]; 
  BigAEM(3,23) = -0.4472135954999579*cM[3]; 
  BigAEM(4,18) = -0.5*cM[4]; 
  BigAEM(4,19) = -0.4472135954999579*cM[1]; 
  BigAEM(4,21) = -0.4472135954999579*cM[3]; 
  BigAEM(4,22) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(5,18) = -0.5*cM[5]; 
  BigAEM(5,20) = -0.4472135954999579*cM[2]; 
  BigAEM(5,21) = -0.4472135954999579*cM[3]; 
  BigAEM(5,23) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(18,0) = 0.5*m1[0]; 
  BigAEM(18,1) = 0.5*m1[1]; 
  BigAEM(18,2) = 0.5*m1[2]; 
  BigAEM(18,3) = 0.5*m1[3]; 
  BigAEM(18,4) = 0.5*m1[4]; 
  BigAEM(18,5) = 0.5*m1[5]; 
  BigAEM(19,0) = 0.5*m1[1]; 
  BigAEM(19,1) = 0.4472135954999579*m1[4]+0.5*m1[0]; 
  BigAEM(19,2) = 0.5*m1[3]; 
  BigAEM(19,3) = 0.5*m1[2]; 
  BigAEM(19,4) = 0.4472135954999579*m1[1]; 
  BigAEM(20,0) = 0.5*m1[2]; 
  BigAEM(20,1) = 0.5*m1[3]; 
  BigAEM(20,2) = 0.4472135954999579*m1[5]+0.5*m1[0]; 
  BigAEM(20,3) = 0.5*m1[1]; 
  BigAEM(20,5) = 0.4472135954999579*m1[2]; 
  BigAEM(21,0) = 0.5*m1[3]; 
  BigAEM(21,1) = 0.5*m1[2]; 
  BigAEM(21,2) = 0.5*m1[1]; 
  BigAEM(21,3) = 0.4472135954999579*m1[5]+0.4472135954999579*m1[4]+0.5*m1[0]; 
  BigAEM(21,4) = 0.4472135954999579*m1[3]; 
  BigAEM(21,5) = 0.4472135954999579*m1[3]; 
  BigAEM(22,0) = 0.5*m1[4]; 
  BigAEM(22,1) = 0.4472135954999579*m1[1]; 
  BigAEM(22,3) = 0.4472135954999579*m1[3]; 
  BigAEM(22,4) = 0.31943828249997*m1[4]+0.5*m1[0]; 
  BigAEM(23,0) = 0.5*m1[5]; 
  BigAEM(23,2) = 0.4472135954999579*m1[2]; 
  BigAEM(23,3) = 0.4472135954999579*m1[3]; 
  BigAEM(23,5) = 0.31943828249997*m1[5]+0.5*m1[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(6,6) = 0.5*m0s[0]; 
  BigAEM(6,7) = 0.5*m0s[1]; 
  BigAEM(6,8) = 0.5*m0s[2]; 
  BigAEM(6,9) = 0.5*m0s[3]; 
  BigAEM(6,10) = 0.5*m0s[4]; 
  BigAEM(6,11) = 0.5*m0s[5]; 
  BigAEM(7,6) = 0.5*m0s[1]; 
  BigAEM(7,7) = 0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(7,8) = 0.5*m0s[3]; 
  BigAEM(7,9) = 0.5*m0s[2]; 
  BigAEM(7,10) = 0.4472135954999579*m0s[1]; 
  BigAEM(8,6) = 0.5*m0s[2]; 
  BigAEM(8,7) = 0.5*m0s[3]; 
  BigAEM(8,8) = 0.4472135954999579*m0s[5]+0.5*m0s[0]; 
  BigAEM(8,9) = 0.5*m0s[1]; 
  BigAEM(8,11) = 0.4472135954999579*m0s[2]; 
  BigAEM(9,6) = 0.5*m0s[3]; 
  BigAEM(9,7) = 0.5*m0s[2]; 
  BigAEM(9,8) = 0.5*m0s[1]; 
  BigAEM(9,9) = 0.4472135954999579*m0s[5]+0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(9,10) = 0.4472135954999579*m0s[3]; 
  BigAEM(9,11) = 0.4472135954999579*m0s[3]; 
  BigAEM(10,6) = 0.5*m0s[4]; 
  BigAEM(10,7) = 0.4472135954999579*m0s[1]; 
  BigAEM(10,9) = 0.4472135954999579*m0s[3]; 
  BigAEM(10,10) = 0.31943828249997*m0s[4]+0.5*m0s[0]; 
  BigAEM(11,6) = 0.5*m0s[5]; 
  BigAEM(11,8) = 0.4472135954999579*m0s[2]; 
  BigAEM(11,9) = 0.4472135954999579*m0s[3]; 
  BigAEM(11,11) = 0.31943828249997*m0s[5]+0.5*m0s[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(6,18) = -0.5*cM[6]; 
  BigAEM(6,19) = -0.5*cM[7]; 
  BigAEM(6,20) = -0.5*cM[8]; 
  BigAEM(6,21) = -0.5*cM[9]; 
  BigAEM(6,22) = -0.5*cM[10]; 
  BigAEM(6,23) = -0.5*cM[11]; 
  BigAEM(7,18) = -0.5*cM[7]; 
  BigAEM(7,19) = (-0.4472135954999579*cM[10])-0.5*cM[6]; 
  BigAEM(7,20) = -0.5*cM[9]; 
  BigAEM(7,21) = -0.5*cM[8]; 
  BigAEM(7,22) = -0.4472135954999579*cM[7]; 
  BigAEM(8,18) = -0.5*cM[8]; 
  BigAEM(8,19) = -0.5*cM[9]; 
  BigAEM(8,20) = (-0.4472135954999579*cM[11])-0.5*cM[6]; 
  BigAEM(8,21) = -0.5*cM[7]; 
  BigAEM(8,23) = -0.4472135954999579*cM[8]; 
  BigAEM(9,18) = -0.5*cM[9]; 
  BigAEM(9,19) = -0.5*cM[8]; 
  BigAEM(9,20) = -0.5*cM[7]; 
  BigAEM(9,21) = (-0.4472135954999579*cM[11])-0.4472135954999579*cM[10]-0.5*cM[6]; 
  BigAEM(9,22) = -0.4472135954999579*cM[9]; 
  BigAEM(9,23) = -0.4472135954999579*cM[9]; 
  BigAEM(10,18) = -0.5*cM[10]; 
  BigAEM(10,19) = -0.4472135954999579*cM[7]; 
  BigAEM(10,21) = -0.4472135954999579*cM[9]; 
  BigAEM(10,22) = (-0.31943828249997*cM[10])-0.5*cM[6]; 
  BigAEM(11,18) = -0.5*cM[11]; 
  BigAEM(11,20) = -0.4472135954999579*cM[8]; 
  BigAEM(11,21) = -0.4472135954999579*cM[9]; 
  BigAEM(11,23) = (-0.31943828249997*cM[11])-0.5*cM[6]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(18,6) = 0.5*m1[6]; 
  BigAEM(18,7) = 0.5*m1[7]; 
  BigAEM(18,8) = 0.5*m1[8]; 
  BigAEM(18,9) = 0.5*m1[9]; 
  BigAEM(18,10) = 0.5*m1[10]; 
  BigAEM(18,11) = 0.5*m1[11]; 
  BigAEM(19,6) = 0.5*m1[7]; 
  BigAEM(19,7) = 0.4472135954999579*m1[10]+0.5*m1[6]; 
  BigAEM(19,8) = 0.5*m1[9]; 
  BigAEM(19,9) = 0.5*m1[8]; 
  BigAEM(19,10) = 0.4472135954999579*m1[7]; 
  BigAEM(20,6) = 0.5*m1[8]; 
  BigAEM(20,7) = 0.5*m1[9]; 
  BigAEM(20,8) = 0.4472135954999579*m1[11]+0.5*m1[6]; 
  BigAEM(20,9) = 0.5*m1[7]; 
  BigAEM(20,11) = 0.4472135954999579*m1[8]; 
  BigAEM(21,6) = 0.5*m1[9]; 
  BigAEM(21,7) = 0.5*m1[8]; 
  BigAEM(21,8) = 0.5*m1[7]; 
  BigAEM(21,9) = 0.4472135954999579*m1[11]+0.4472135954999579*m1[10]+0.5*m1[6]; 
  BigAEM(21,10) = 0.4472135954999579*m1[9]; 
  BigAEM(21,11) = 0.4472135954999579*m1[9]; 
  BigAEM(22,6) = 0.5*m1[10]; 
  BigAEM(22,7) = 0.4472135954999579*m1[7]; 
  BigAEM(22,9) = 0.4472135954999579*m1[9]; 
  BigAEM(22,10) = 0.31943828249997*m1[10]+0.5*m1[6]; 
  BigAEM(23,6) = 0.5*m1[11]; 
  BigAEM(23,8) = 0.4472135954999579*m1[8]; 
  BigAEM(23,9) = 0.4472135954999579*m1[9]; 
  BigAEM(23,11) = 0.31943828249997*m1[11]+0.5*m1[6]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  BigAEM(12,12) = 0.5*m0s[0]; 
  BigAEM(12,13) = 0.5*m0s[1]; 
  BigAEM(12,14) = 0.5*m0s[2]; 
  BigAEM(12,15) = 0.5*m0s[3]; 
  BigAEM(12,16) = 0.5*m0s[4]; 
  BigAEM(12,17) = 0.5*m0s[5]; 
  BigAEM(13,12) = 0.5*m0s[1]; 
  BigAEM(13,13) = 0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(13,14) = 0.5*m0s[3]; 
  BigAEM(13,15) = 0.5*m0s[2]; 
  BigAEM(13,16) = 0.4472135954999579*m0s[1]; 
  BigAEM(14,12) = 0.5*m0s[2]; 
  BigAEM(14,13) = 0.5*m0s[3]; 
  BigAEM(14,14) = 0.4472135954999579*m0s[5]+0.5*m0s[0]; 
  BigAEM(14,15) = 0.5*m0s[1]; 
  BigAEM(14,17) = 0.4472135954999579*m0s[2]; 
  BigAEM(15,12) = 0.5*m0s[3]; 
  BigAEM(15,13) = 0.5*m0s[2]; 
  BigAEM(15,14) = 0.5*m0s[1]; 
  BigAEM(15,15) = 0.4472135954999579*m0s[5]+0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(15,16) = 0.4472135954999579*m0s[3]; 
  BigAEM(15,17) = 0.4472135954999579*m0s[3]; 
  BigAEM(16,12) = 0.5*m0s[4]; 
  BigAEM(16,13) = 0.4472135954999579*m0s[1]; 
  BigAEM(16,15) = 0.4472135954999579*m0s[3]; 
  BigAEM(16,16) = 0.31943828249997*m0s[4]+0.5*m0s[0]; 
  BigAEM(17,12) = 0.5*m0s[5]; 
  BigAEM(17,14) = 0.4472135954999579*m0s[2]; 
  BigAEM(17,15) = 0.4472135954999579*m0s[3]; 
  BigAEM(17,17) = 0.31943828249997*m0s[5]+0.5*m0s[0]; 
 
  // ....... Block from correction to uZ .......... // 
  BigAEM(12,18) = -0.5*cM[12]; 
  BigAEM(12,19) = -0.5*cM[13]; 
  BigAEM(12,20) = -0.5*cM[14]; 
  BigAEM(12,21) = -0.5*cM[15]; 
  BigAEM(12,22) = -0.5*cM[16]; 
  BigAEM(12,23) = -0.5*cM[17]; 
  BigAEM(13,18) = -0.5*cM[13]; 
  BigAEM(13,19) = (-0.4472135954999579*cM[16])-0.5*cM[12]; 
  BigAEM(13,20) = -0.5*cM[15]; 
  BigAEM(13,21) = -0.5*cM[14]; 
  BigAEM(13,22) = -0.4472135954999579*cM[13]; 
  BigAEM(14,18) = -0.5*cM[14]; 
  BigAEM(14,19) = -0.5*cM[15]; 
  BigAEM(14,20) = (-0.4472135954999579*cM[17])-0.5*cM[12]; 
  BigAEM(14,21) = -0.5*cM[13]; 
  BigAEM(14,23) = -0.4472135954999579*cM[14]; 
  BigAEM(15,18) = -0.5*cM[15]; 
  BigAEM(15,19) = -0.5*cM[14]; 
  BigAEM(15,20) = -0.5*cM[13]; 
  BigAEM(15,21) = (-0.4472135954999579*cM[17])-0.4472135954999579*cM[16]-0.5*cM[12]; 
  BigAEM(15,22) = -0.4472135954999579*cM[15]; 
  BigAEM(15,23) = -0.4472135954999579*cM[15]; 
  BigAEM(16,18) = -0.5*cM[16]; 
  BigAEM(16,19) = -0.4472135954999579*cM[13]; 
  BigAEM(16,21) = -0.4472135954999579*cM[15]; 
  BigAEM(16,22) = (-0.31943828249997*cM[16])-0.5*cM[12]; 
  BigAEM(17,18) = -0.5*cM[17]; 
  BigAEM(17,20) = -0.4472135954999579*cM[14]; 
  BigAEM(17,21) = -0.4472135954999579*cM[15]; 
  BigAEM(17,23) = (-0.31943828249997*cM[17])-0.5*cM[12]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  BigAEM(18,12) = 0.5*m1[12]; 
  BigAEM(18,13) = 0.5*m1[13]; 
  BigAEM(18,14) = 0.5*m1[14]; 
  BigAEM(18,15) = 0.5*m1[15]; 
  BigAEM(18,16) = 0.5*m1[16]; 
  BigAEM(18,17) = 0.5*m1[17]; 
  BigAEM(19,12) = 0.5*m1[13]; 
  BigAEM(19,13) = 0.4472135954999579*m1[16]+0.5*m1[12]; 
  BigAEM(19,14) = 0.5*m1[15]; 
  BigAEM(19,15) = 0.5*m1[14]; 
  BigAEM(19,16) = 0.4472135954999579*m1[13]; 
  BigAEM(20,12) = 0.5*m1[14]; 
  BigAEM(20,13) = 0.5*m1[15]; 
  BigAEM(20,14) = 0.4472135954999579*m1[17]+0.5*m1[12]; 
  BigAEM(20,15) = 0.5*m1[13]; 
  BigAEM(20,17) = 0.4472135954999579*m1[14]; 
  BigAEM(21,12) = 0.5*m1[15]; 
  BigAEM(21,13) = 0.5*m1[14]; 
  BigAEM(21,14) = 0.5*m1[13]; 
  BigAEM(21,15) = 0.4472135954999579*m1[17]+0.4472135954999579*m1[16]+0.5*m1[12]; 
  BigAEM(21,16) = 0.4472135954999579*m1[15]; 
  BigAEM(21,17) = 0.4472135954999579*m1[15]; 
  BigAEM(22,12) = 0.5*m1[16]; 
  BigAEM(22,13) = 0.4472135954999579*m1[13]; 
  BigAEM(22,15) = 0.4472135954999579*m1[15]; 
  BigAEM(22,16) = 0.31943828249997*m1[16]+0.5*m1[12]; 
  BigAEM(23,12) = 0.5*m1[17]; 
  BigAEM(23,14) = 0.4472135954999579*m1[14]; 
  BigAEM(23,15) = 0.4472135954999579*m1[15]; 
  BigAEM(23,17) = 0.31943828249997*m1[17]+0.5*m1[12]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(18,18) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
  BigAEM(18,19) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(18,20) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(18,21) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(18,22) = 0.5*m0s[4]*pVdim-0.5*cE[4]; 
  BigAEM(18,23) = 0.5*m0s[5]*pVdim-0.5*cE[5]; 
  BigAEM(19,18) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(19,19) = 0.4472135954999579*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[4]-0.5*cE[0]; 
  BigAEM(19,20) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(19,21) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(19,22) = 0.4472135954999579*m0s[1]*pVdim-0.4472135954999579*cE[1]; 
  BigAEM(20,18) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(20,19) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(20,20) = 0.4472135954999579*m0s[5]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[5]-0.5*cE[0]; 
  BigAEM(20,21) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(20,23) = 0.4472135954999579*m0s[2]*pVdim-0.4472135954999579*cE[2]; 
  BigAEM(21,18) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(21,19) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(21,20) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(21,21) = 0.4472135954999579*m0s[5]*pVdim+0.4472135954999579*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[5]-0.4472135954999579*cE[4]-0.5*cE[0]; 
  BigAEM(21,22) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(21,23) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(22,18) = 0.5*m0s[4]*pVdim-0.5*cE[4]; 
  BigAEM(22,19) = 0.4472135954999579*m0s[1]*pVdim-0.4472135954999579*cE[1]; 
  BigAEM(22,21) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(22,22) = 0.31943828249997*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.31943828249997*cE[4]-0.5*cE[0]; 
  BigAEM(23,18) = 0.5*m0s[5]*pVdim-0.5*cE[5]; 
  BigAEM(23,20) = 0.4472135954999579*m0s[2]*pVdim-0.4472135954999579*cE[2]; 
  BigAEM(23,21) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(23,23) = 0.31943828249997*m0s[5]*pVdim+0.5*m0s[0]*pVdim-0.31943828249997*cE[5]-0.5*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<6,12>(0,6).setZero(); 
  BigAEM.block<12,6>(6,0).setZero(); 
  BigAEM.block<6,6>(6,12).setZero(); 
  BigAEM.block<6,6>(12,6).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m1[2],m1[3],m1[4],m1[5],m1[6],m1[7],m1[8],m1[9],m1[10],m1[11],m1[12],m1[13],m1[14],m1[15],m1[16],m1[17],m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,18,1) = xEV.segment<18>(0); 
 
  Eigen::Map<VectorXd>(vtSq,6,1) = xEV.segment<6>(18); 
 
} 
 
void BoundaryIntegral2x3vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  out[0] += (2.449489742783178*fvmin[3]*dS+2.449489742783178*fvmax[3]*dS-1.414213562373095*fvmin[0]*dS+1.414213562373095*fvmax[0]*dS)*intFac; 
  out[1] += (1.414213562373095*fvmax[1]*dS-1.414213562373095*fvmin[1]*dS)*intFac; 
  out[2] += (1.414213562373095*fvmax[2]*dS-1.414213562373095*fvmin[2]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x3vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  out[0] += ((-3.16227766016838*fvmin[18]*dS)+3.16227766016838*fvmax[18]*dS+2.449489742783178*fvmin[3]*dS+2.449489742783178*fvmax[3]*dS-1.414213562373095*fvmin[0]*dS+1.414213562373095*fvmax[0]*dS)*intFac; 
  out[1] += (2.449489742783178*fvmin[7]*dS+2.449489742783178*fvmax[7]*dS-1.414213562373095*fvmin[1]*dS+1.414213562373095*fvmax[1]*dS)*intFac; 
  out[2] += (2.449489742783178*fvmin[8]*dS+2.449489742783178*fvmax[8]*dS-1.414213562373095*fvmin[2]*dS+1.414213562373095*fvmax[2]*dS)*intFac; 
  out[3] += (1.414213562373095*fvmax[6]*dS-1.414213562373095*fvmin[6]*dS)*intFac; 
  out[4] += (1.414213562373095*fvmax[16]*dS-1.414213562373095*fvmin[16]*dS)*intFac; 
  out[5] += (1.414213562373095*fvmax[17]*dS-1.414213562373095*fvmin[17]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x3vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  out[0] += intFac*(2.449489742783178*fvmin[3]*dS*vmin-1.414213562373095*fvmin[0]*dS*vmin+2.449489742783178*fvmax[3]*dS*vmax+1.414213562373095*fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.414213562373095*fvmax[1]*dS*vmax-1.414213562373095*fvmin[1]*dS*vmin); 
  out[2] += intFac*(1.414213562373095*fvmax[2]*dS*vmax-1.414213562373095*fvmin[2]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x3vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[3]*dxv[4]; 
 
  out[0] += intFac*((-3.16227766016838*fvmin[18]*dS*vmin)+2.449489742783178*fvmin[3]*dS*vmin-1.414213562373095*fvmin[0]*dS*vmin+3.16227766016838*fvmax[18]*dS*vmax+2.449489742783178*fvmax[3]*dS*vmax+1.414213562373095*fvmax[0]*dS*vmax); 
  out[1] += intFac*(2.449489742783178*fvmin[7]*dS*vmin-1.414213562373095*fvmin[1]*dS*vmin+2.449489742783178*fvmax[7]*dS*vmax+1.414213562373095*fvmax[1]*dS*vmax); 
  out[2] += intFac*(2.449489742783178*fvmin[8]*dS*vmin-1.414213562373095*fvmin[2]*dS*vmin+2.449489742783178*fvmax[8]*dS*vmax+1.414213562373095*fvmax[2]*dS*vmax); 
  out[3] += intFac*(1.414213562373095*fvmax[6]*dS*vmax-1.414213562373095*fvmin[6]*dS*vmin); 
  out[4] += intFac*(1.414213562373095*fvmax[16]*dS*vmax-1.414213562373095*fvmin[16]*dS*vmin); 
  out[5] += intFac*(1.414213562373095*fvmax[17]*dS*vmax-1.414213562373095*fvmin[17]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x3vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  out[3] += (2.449489742783178*fvmin[4]*dS+2.449489742783178*fvmax[4]*dS-1.414213562373095*fvmin[0]*dS+1.414213562373095*fvmax[0]*dS)*intFac; 
  out[4] += (1.414213562373095*fvmax[1]*dS-1.414213562373095*fvmin[1]*dS)*intFac; 
  out[5] += (1.414213562373095*fvmax[2]*dS-1.414213562373095*fvmin[2]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x3vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  out[6] += ((-3.16227766016838*fvmin[19]*dS)+3.16227766016838*fvmax[19]*dS+2.449489742783178*fvmin[4]*dS+2.449489742783178*fvmax[4]*dS-1.414213562373095*fvmin[0]*dS+1.414213562373095*fvmax[0]*dS)*intFac; 
  out[7] += (2.449489742783178*fvmin[9]*dS+2.449489742783178*fvmax[9]*dS-1.414213562373095*fvmin[1]*dS+1.414213562373095*fvmax[1]*dS)*intFac; 
  out[8] += (2.449489742783178*fvmin[10]*dS+2.449489742783178*fvmax[10]*dS-1.414213562373095*fvmin[2]*dS+1.414213562373095*fvmax[2]*dS)*intFac; 
  out[9] += (1.414213562373095*fvmax[6]*dS-1.414213562373095*fvmin[6]*dS)*intFac; 
  out[10] += (1.414213562373095*fvmax[16]*dS-1.414213562373095*fvmin[16]*dS)*intFac; 
  out[11] += (1.414213562373095*fvmax[17]*dS-1.414213562373095*fvmin[17]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x3vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  out[0] += intFac*(2.449489742783178*fvmin[4]*dS*vmin-1.414213562373095*fvmin[0]*dS*vmin+2.449489742783178*fvmax[4]*dS*vmax+1.414213562373095*fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.414213562373095*fvmax[1]*dS*vmax-1.414213562373095*fvmin[1]*dS*vmin); 
  out[2] += intFac*(1.414213562373095*fvmax[2]*dS*vmax-1.414213562373095*fvmin[2]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x3vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[4]; 
 
  out[0] += intFac*((-3.16227766016838*fvmin[19]*dS*vmin)+2.449489742783178*fvmin[4]*dS*vmin-1.414213562373095*fvmin[0]*dS*vmin+3.16227766016838*fvmax[19]*dS*vmax+2.449489742783178*fvmax[4]*dS*vmax+1.414213562373095*fvmax[0]*dS*vmax); 
  out[1] += intFac*(2.449489742783178*fvmin[9]*dS*vmin-1.414213562373095*fvmin[1]*dS*vmin+2.449489742783178*fvmax[9]*dS*vmax+1.414213562373095*fvmax[1]*dS*vmax); 
  out[2] += intFac*(2.449489742783178*fvmin[10]*dS*vmin-1.414213562373095*fvmin[2]*dS*vmin+2.449489742783178*fvmax[10]*dS*vmax+1.414213562373095*fvmax[2]*dS*vmax); 
  out[3] += intFac*(1.414213562373095*fvmax[6]*dS*vmax-1.414213562373095*fvmin[6]*dS*vmin); 
  out[4] += intFac*(1.414213562373095*fvmax[16]*dS*vmax-1.414213562373095*fvmin[16]*dS*vmin); 
  out[5] += intFac*(1.414213562373095*fvmax[17]*dS*vmax-1.414213562373095*fvmin[17]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x3vMax_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  out[6] += (2.449489742783178*fvmin[5]*dS+2.449489742783178*fvmax[5]*dS-1.414213562373095*fvmin[0]*dS+1.414213562373095*fvmax[0]*dS)*intFac; 
  out[7] += (1.414213562373095*fvmax[1]*dS-1.414213562373095*fvmin[1]*dS)*intFac; 
  out[8] += (1.414213562373095*fvmax[2]*dS-1.414213562373095*fvmin[2]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x3vMax_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  out[12] += ((-3.16227766016838*fvmin[20]*dS)+3.16227766016838*fvmax[20]*dS+2.449489742783178*fvmin[5]*dS+2.449489742783178*fvmax[5]*dS-1.414213562373095*fvmin[0]*dS+1.414213562373095*fvmax[0]*dS)*intFac; 
  out[13] += (2.449489742783178*fvmin[12]*dS+2.449489742783178*fvmax[12]*dS-1.414213562373095*fvmin[1]*dS+1.414213562373095*fvmax[1]*dS)*intFac; 
  out[14] += (2.449489742783178*fvmin[13]*dS+2.449489742783178*fvmax[13]*dS-1.414213562373095*fvmin[2]*dS+1.414213562373095*fvmax[2]*dS)*intFac; 
  out[15] += (1.414213562373095*fvmax[6]*dS-1.414213562373095*fvmin[6]*dS)*intFac; 
  out[16] += (1.414213562373095*fvmax[16]*dS-1.414213562373095*fvmin[16]*dS)*intFac; 
  out[17] += (1.414213562373095*fvmax[17]*dS-1.414213562373095*fvmin[17]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x3vMax_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  out[0] += intFac*(2.449489742783178*fvmin[5]*dS*vmin-1.414213562373095*fvmin[0]*dS*vmin+2.449489742783178*fvmax[5]*dS*vmax+1.414213562373095*fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.414213562373095*fvmax[1]*dS*vmax-1.414213562373095*fvmin[1]*dS*vmin); 
  out[2] += intFac*(1.414213562373095*fvmax[2]*dS*vmax-1.414213562373095*fvmin[2]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x3vMax_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]; 
 
  out[0] += intFac*((-3.16227766016838*fvmin[20]*dS*vmin)+2.449489742783178*fvmin[5]*dS*vmin-1.414213562373095*fvmin[0]*dS*vmin+3.16227766016838*fvmax[20]*dS*vmax+2.449489742783178*fvmax[5]*dS*vmax+1.414213562373095*fvmax[0]*dS*vmax); 
  out[1] += intFac*(2.449489742783178*fvmin[12]*dS*vmin-1.414213562373095*fvmin[1]*dS*vmin+2.449489742783178*fvmax[12]*dS*vmax+1.414213562373095*fvmax[1]*dS*vmax); 
  out[2] += intFac*(2.449489742783178*fvmin[13]*dS*vmin-1.414213562373095*fvmin[2]*dS*vmin+2.449489742783178*fvmax[13]*dS*vmax+1.414213562373095*fvmax[2]*dS*vmax); 
  out[3] += intFac*(1.414213562373095*fvmax[6]*dS*vmax-1.414213562373095*fvmin[6]*dS*vmin); 
  out[4] += intFac*(1.414213562373095*fvmax[16]*dS*vmax-1.414213562373095*fvmin[16]*dS*vmin); 
  out[5] += intFac*(1.414213562373095*fvmax[17]*dS*vmax-1.414213562373095*fvmin[17]*dS*vmin); 
 
} 
 
