#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments1x2vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool m0Avg = false;
  if (0.7071067811865475*m0[0]-1.224744871391589*m0[1] < 0) { 
    m0Avg = true;
  }
  if (1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    m0Avg = true;
  }
 
  double m0s[2]; 
  double m1s[4]; 
  double m2s[2]; 
  if (m0Avg) { 
    m0s[0] = m0[0]; 
    m0s[1] = 0.0; 
    m1s[0] = m1[0]; 
    m1s[1] = 0.0; 
    m1s[2] = m1[2]; 
    m1s[3] = 0.0; 
    m2s[0] = m2[0]; 
    m2s[1] = 0.0; 
  } else { 
    m0s[0] = m0[0]; 
    m0s[1] = m0[1]; 
    m1s[0] = m1[0]; 
    m1s[1] = m1[1]; 
    m1s[2] = m1[2]; 
    m1s[3] = m1[3]; 
    m2s[0] = m2[0]; 
    m2s[1] = m2[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(6,6); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(6);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(6);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0s[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0s[1]; 
  BigAEM(1,0) = 0.7071067811865475*m0s[1]; 
  BigAEM(1,1) = 0.7071067811865475*m0s[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,4) = -0.7071067811865475*cM[0]; 
  BigAEM(0,5) = -0.7071067811865475*cM[1]; 
  BigAEM(1,4) = -0.7071067811865475*cM[1]; 
  BigAEM(1,5) = -0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(4,0) = 0.7071067811865475*m1s[0]; 
  BigAEM(4,1) = 0.7071067811865475*m1s[1]; 
  BigAEM(5,0) = 0.7071067811865475*m1s[1]; 
  BigAEM(5,1) = 0.7071067811865475*m1s[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(2,2) = 0.7071067811865475*m0s[0]; 
  BigAEM(2,3) = 0.7071067811865475*m0s[1]; 
  BigAEM(3,2) = 0.7071067811865475*m0s[1]; 
  BigAEM(3,3) = 0.7071067811865475*m0s[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(2,4) = -0.7071067811865475*cM[2]; 
  BigAEM(2,5) = -0.7071067811865475*cM[3]; 
  BigAEM(3,4) = -0.7071067811865475*cM[3]; 
  BigAEM(3,5) = -0.7071067811865475*cM[2]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(4,2) = 0.7071067811865475*m1s[2]; 
  BigAEM(4,3) = 0.7071067811865475*m1s[3]; 
  BigAEM(5,2) = 0.7071067811865475*m1s[3]; 
  BigAEM(5,3) = 0.7071067811865475*m1s[2]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(4,4) = 0.7071067811865475*m0s[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(4,5) = 0.7071067811865475*m0s[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(5,4) = 0.7071067811865475*m0s[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(5,5) = 0.7071067811865475*m0s[0]*pVdim-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<2,2>(0,2).setZero(); 
  BigAEM.block<2,2>(2,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1s[0],m1s[1],m1s[2],m1s[3],m2s[0],m2s[1]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,4,1) = xEV.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = xEV.segment<2>(4); 
 
} 
 
void SelfPrimMoments1x2vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool m0Avg = false;
  if (1.58113883008419*m0[2]-1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    m0Avg = true;
  }
  if (1.58113883008419*m0[2]+1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    m0Avg = true;
  }
 
  double m0s[3]; 
  double m1s[6]; 
  double m2s[3]; 
  if (m0Avg) { 
    m0s[0] = m0[0]; 
    m0s[1] = 0.0; 
    m0s[2] = 0.0; 
    m1s[0] = m1[0]; 
    m1s[1] = 0.0; 
    m1s[2] = 0.0; 
    m1s[3] = m1[3]; 
    m1s[4] = 0.0; 
    m1s[5] = 0.0; 
    m2s[0] = m2[0]; 
    m2s[1] = 0.0; 
    m2s[2] = 0.0; 
  } else { 
    m0s[0] = m0[0]; 
    m0s[1] = m0[1]; 
    m0s[2] = m0[2]; 
    m1s[0] = m1[0]; 
    m1s[1] = m1[1]; 
    m1s[2] = m1[2]; 
    m1s[3] = m1[3]; 
    m1s[4] = m1[4]; 
    m1s[5] = m1[5]; 
    m2s[0] = m2[0]; 
    m2s[1] = m2[1]; 
    m2s[2] = m2[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(9,9); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(9);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(9);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0s[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0s[1]; 
  BigAEM(0,2) = 0.7071067811865475*m0s[2]; 
  BigAEM(1,0) = 0.7071067811865475*m0s[1]; 
  BigAEM(1,1) = 0.6324555320336759*m0s[2]+0.7071067811865475*m0s[0]; 
  BigAEM(1,2) = 0.6324555320336759*m0s[1]; 
  BigAEM(2,0) = 0.7071067811865475*m0s[2]; 
  BigAEM(2,1) = 0.6324555320336759*m0s[1]; 
  BigAEM(2,2) = 0.4517539514526256*m0s[2]+0.7071067811865475*m0s[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,6) = -0.7071067811865475*cM[0]; 
  BigAEM(0,7) = -0.7071067811865475*cM[1]; 
  BigAEM(0,8) = -0.7071067811865475*cM[2]; 
  BigAEM(1,6) = -0.7071067811865475*cM[1]; 
  BigAEM(1,7) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(1,8) = -0.6324555320336759*cM[1]; 
  BigAEM(2,6) = -0.7071067811865475*cM[2]; 
  BigAEM(2,7) = -0.6324555320336759*cM[1]; 
  BigAEM(2,8) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(6,0) = 0.7071067811865475*m1s[0]; 
  BigAEM(6,1) = 0.7071067811865475*m1s[1]; 
  BigAEM(6,2) = 0.7071067811865475*m1s[2]; 
  BigAEM(7,0) = 0.7071067811865475*m1s[1]; 
  BigAEM(7,1) = 0.6324555320336759*m1s[2]+0.7071067811865475*m1s[0]; 
  BigAEM(7,2) = 0.6324555320336759*m1s[1]; 
  BigAEM(8,0) = 0.7071067811865475*m1s[2]; 
  BigAEM(8,1) = 0.6324555320336759*m1s[1]; 
  BigAEM(8,2) = 0.4517539514526256*m1s[2]+0.7071067811865475*m1s[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(3,3) = 0.7071067811865475*m0s[0]; 
  BigAEM(3,4) = 0.7071067811865475*m0s[1]; 
  BigAEM(3,5) = 0.7071067811865475*m0s[2]; 
  BigAEM(4,3) = 0.7071067811865475*m0s[1]; 
  BigAEM(4,4) = 0.6324555320336759*m0s[2]+0.7071067811865475*m0s[0]; 
  BigAEM(4,5) = 0.6324555320336759*m0s[1]; 
  BigAEM(5,3) = 0.7071067811865475*m0s[2]; 
  BigAEM(5,4) = 0.6324555320336759*m0s[1]; 
  BigAEM(5,5) = 0.4517539514526256*m0s[2]+0.7071067811865475*m0s[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(3,6) = -0.7071067811865475*cM[3]; 
  BigAEM(3,7) = -0.7071067811865475*cM[4]; 
  BigAEM(3,8) = -0.7071067811865475*cM[5]; 
  BigAEM(4,6) = -0.7071067811865475*cM[4]; 
  BigAEM(4,7) = (-0.6324555320336759*cM[5])-0.7071067811865475*cM[3]; 
  BigAEM(4,8) = -0.6324555320336759*cM[4]; 
  BigAEM(5,6) = -0.7071067811865475*cM[5]; 
  BigAEM(5,7) = -0.6324555320336759*cM[4]; 
  BigAEM(5,8) = (-0.4517539514526256*cM[5])-0.7071067811865475*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(6,3) = 0.7071067811865475*m1s[3]; 
  BigAEM(6,4) = 0.7071067811865475*m1s[4]; 
  BigAEM(6,5) = 0.7071067811865475*m1s[5]; 
  BigAEM(7,3) = 0.7071067811865475*m1s[4]; 
  BigAEM(7,4) = 0.6324555320336759*m1s[5]+0.7071067811865475*m1s[3]; 
  BigAEM(7,5) = 0.6324555320336759*m1s[4]; 
  BigAEM(8,3) = 0.7071067811865475*m1s[5]; 
  BigAEM(8,4) = 0.6324555320336759*m1s[4]; 
  BigAEM(8,5) = 0.4517539514526256*m1s[5]+0.7071067811865475*m1s[3]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(6,6) = 0.7071067811865475*m0s[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(6,7) = 0.7071067811865475*m0s[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(6,8) = 0.7071067811865475*m0s[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(7,6) = 0.7071067811865475*m0s[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(7,7) = 0.6324555320336759*m0s[2]*pVdim+0.7071067811865475*m0s[0]*pVdim-0.6324555320336759*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(7,8) = 0.6324555320336759*m0s[1]*pVdim-0.6324555320336759*cE[1]; 
  BigAEM(8,6) = 0.7071067811865475*m0s[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(8,7) = 0.6324555320336759*m0s[1]*pVdim-0.6324555320336759*cE[1]; 
  BigAEM(8,8) = 0.4517539514526256*m0s[2]*pVdim+0.7071067811865475*m0s[0]*pVdim-0.4517539514526256*cE[2]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<3,3>(0,3).setZero(); 
  BigAEM.block<3,3>(3,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1s[0],m1s[1],m1s[2],m1s[3],m1s[4],m1s[5],m2s[0],m2s[1],m2s[2]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,6,1) = xEV.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV.segment<3>(6); 
 
} 
 
void SelfPrimMoments1x2vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool m0Avg = false;
  if ((-1.870828693386971*m0[3])+1.58113883008419*m0[2]-1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    m0Avg = true;
  }
  if (1.870828693386971*m0[3]+1.58113883008419*m0[2]+1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    m0Avg = true;
  }
 
  double m0s[4]; 
  double m1s[8]; 
  double m2s[4]; 
  if (m0Avg) { 
    m0s[0] = m0[0]; 
    m0s[1] = 0.0; 
    m0s[2] = 0.0; 
    m0s[3] = 0.0; 
    m1s[0] = m1[0]; 
    m1s[1] = 0.0; 
    m1s[2] = 0.0; 
    m1s[3] = 0.0; 
    m1s[4] = m1[4]; 
    m1s[5] = 0.0; 
    m1s[6] = 0.0; 
    m1s[7] = 0.0; 
    m2s[0] = m2[0]; 
    m2s[1] = 0.0; 
    m2s[2] = 0.0; 
    m2s[3] = 0.0; 
  } else { 
    m0s[0] = m0[0]; 
    m0s[1] = m0[1]; 
    m0s[2] = m0[2]; 
    m0s[3] = m0[3]; 
    m1s[0] = m1[0]; 
    m1s[1] = m1[1]; 
    m1s[2] = m1[2]; 
    m1s[3] = m1[3]; 
    m1s[4] = m1[4]; 
    m1s[5] = m1[5]; 
    m1s[6] = m1[6]; 
    m1s[7] = m1[7]; 
    m2s[0] = m2[0]; 
    m2s[1] = m2[1]; 
    m2s[2] = m2[2]; 
    m2s[3] = m2[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(12,12); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(12);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(12);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0s[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0s[1]; 
  BigAEM(0,2) = 0.7071067811865475*m0s[2]; 
  BigAEM(0,3) = 0.7071067811865475*m0s[3]; 
  BigAEM(1,0) = 0.7071067811865475*m0s[1]; 
  BigAEM(1,1) = 0.6324555320336759*m0s[2]+0.7071067811865475*m0s[0]; 
  BigAEM(1,2) = 0.6210590034081186*m0s[3]+0.6324555320336759*m0s[1]; 
  BigAEM(1,3) = 0.6210590034081186*m0s[2]; 
  BigAEM(2,0) = 0.7071067811865475*m0s[2]; 
  BigAEM(2,1) = 0.6210590034081186*m0s[3]+0.6324555320336759*m0s[1]; 
  BigAEM(2,2) = 0.4517539514526256*m0s[2]+0.7071067811865475*m0s[0]; 
  BigAEM(2,3) = 0.421637021355784*m0s[3]+0.6210590034081186*m0s[1]; 
  BigAEM(3,0) = 0.7071067811865475*m0s[3]; 
  BigAEM(3,1) = 0.6210590034081186*m0s[2]; 
  BigAEM(3,2) = 0.421637021355784*m0s[3]+0.6210590034081186*m0s[1]; 
  BigAEM(3,3) = 0.421637021355784*m0s[2]+0.7071067811865475*m0s[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,8) = -0.7071067811865475*cM[0]; 
  BigAEM(0,9) = -0.7071067811865475*cM[1]; 
  BigAEM(0,10) = -0.7071067811865475*cM[2]; 
  BigAEM(0,11) = -0.7071067811865475*cM[3]; 
  BigAEM(1,8) = -0.7071067811865475*cM[1]; 
  BigAEM(1,9) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(1,10) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  BigAEM(1,11) = -0.6210590034081186*cM[2]; 
  BigAEM(2,8) = -0.7071067811865475*cM[2]; 
  BigAEM(2,9) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  BigAEM(2,10) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(2,11) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  BigAEM(3,8) = -0.7071067811865475*cM[3]; 
  BigAEM(3,9) = -0.6210590034081186*cM[2]; 
  BigAEM(3,10) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  BigAEM(3,11) = (-0.421637021355784*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(8,0) = 0.7071067811865475*m1s[0]; 
  BigAEM(8,1) = 0.7071067811865475*m1s[1]; 
  BigAEM(8,2) = 0.7071067811865475*m1s[2]; 
  BigAEM(8,3) = 0.7071067811865475*m1s[3]; 
  BigAEM(9,0) = 0.7071067811865475*m1s[1]; 
  BigAEM(9,1) = 0.6324555320336759*m1s[2]+0.7071067811865475*m1s[0]; 
  BigAEM(9,2) = 0.6210590034081186*m1s[3]+0.6324555320336759*m1s[1]; 
  BigAEM(9,3) = 0.6210590034081186*m1s[2]; 
  BigAEM(10,0) = 0.7071067811865475*m1s[2]; 
  BigAEM(10,1) = 0.6210590034081186*m1s[3]+0.6324555320336759*m1s[1]; 
  BigAEM(10,2) = 0.4517539514526256*m1s[2]+0.7071067811865475*m1s[0]; 
  BigAEM(10,3) = 0.421637021355784*m1s[3]+0.6210590034081186*m1s[1]; 
  BigAEM(11,0) = 0.7071067811865475*m1s[3]; 
  BigAEM(11,1) = 0.6210590034081186*m1s[2]; 
  BigAEM(11,2) = 0.421637021355784*m1s[3]+0.6210590034081186*m1s[1]; 
  BigAEM(11,3) = 0.421637021355784*m1s[2]+0.7071067811865475*m1s[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(4,4) = 0.7071067811865475*m0s[0]; 
  BigAEM(4,5) = 0.7071067811865475*m0s[1]; 
  BigAEM(4,6) = 0.7071067811865475*m0s[2]; 
  BigAEM(4,7) = 0.7071067811865475*m0s[3]; 
  BigAEM(5,4) = 0.7071067811865475*m0s[1]; 
  BigAEM(5,5) = 0.6324555320336759*m0s[2]+0.7071067811865475*m0s[0]; 
  BigAEM(5,6) = 0.6210590034081186*m0s[3]+0.6324555320336759*m0s[1]; 
  BigAEM(5,7) = 0.6210590034081186*m0s[2]; 
  BigAEM(6,4) = 0.7071067811865475*m0s[2]; 
  BigAEM(6,5) = 0.6210590034081186*m0s[3]+0.6324555320336759*m0s[1]; 
  BigAEM(6,6) = 0.4517539514526256*m0s[2]+0.7071067811865475*m0s[0]; 
  BigAEM(6,7) = 0.421637021355784*m0s[3]+0.6210590034081186*m0s[1]; 
  BigAEM(7,4) = 0.7071067811865475*m0s[3]; 
  BigAEM(7,5) = 0.6210590034081186*m0s[2]; 
  BigAEM(7,6) = 0.421637021355784*m0s[3]+0.6210590034081186*m0s[1]; 
  BigAEM(7,7) = 0.421637021355784*m0s[2]+0.7071067811865475*m0s[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(4,8) = -0.7071067811865475*cM[4]; 
  BigAEM(4,9) = -0.7071067811865475*cM[5]; 
  BigAEM(4,10) = -0.7071067811865475*cM[6]; 
  BigAEM(4,11) = -0.7071067811865475*cM[7]; 
  BigAEM(5,8) = -0.7071067811865475*cM[5]; 
  BigAEM(5,9) = (-0.6324555320336759*cM[6])-0.7071067811865475*cM[4]; 
  BigAEM(5,10) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  BigAEM(5,11) = -0.6210590034081186*cM[6]; 
  BigAEM(6,8) = -0.7071067811865475*cM[6]; 
  BigAEM(6,9) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  BigAEM(6,10) = (-0.4517539514526256*cM[6])-0.7071067811865475*cM[4]; 
  BigAEM(6,11) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  BigAEM(7,8) = -0.7071067811865475*cM[7]; 
  BigAEM(7,9) = -0.6210590034081186*cM[6]; 
  BigAEM(7,10) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  BigAEM(7,11) = (-0.421637021355784*cM[6])-0.7071067811865475*cM[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(8,4) = 0.7071067811865475*m1s[4]; 
  BigAEM(8,5) = 0.7071067811865475*m1s[5]; 
  BigAEM(8,6) = 0.7071067811865475*m1s[6]; 
  BigAEM(8,7) = 0.7071067811865475*m1s[7]; 
  BigAEM(9,4) = 0.7071067811865475*m1s[5]; 
  BigAEM(9,5) = 0.6324555320336759*m1s[6]+0.7071067811865475*m1s[4]; 
  BigAEM(9,6) = 0.6210590034081186*m1s[7]+0.6324555320336759*m1s[5]; 
  BigAEM(9,7) = 0.6210590034081186*m1s[6]; 
  BigAEM(10,4) = 0.7071067811865475*m1s[6]; 
  BigAEM(10,5) = 0.6210590034081186*m1s[7]+0.6324555320336759*m1s[5]; 
  BigAEM(10,6) = 0.4517539514526256*m1s[6]+0.7071067811865475*m1s[4]; 
  BigAEM(10,7) = 0.421637021355784*m1s[7]+0.6210590034081186*m1s[5]; 
  BigAEM(11,4) = 0.7071067811865475*m1s[7]; 
  BigAEM(11,5) = 0.6210590034081186*m1s[6]; 
  BigAEM(11,6) = 0.421637021355784*m1s[7]+0.6210590034081186*m1s[5]; 
  BigAEM(11,7) = 0.421637021355784*m1s[6]+0.7071067811865475*m1s[4]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(8,8) = 0.7071067811865475*m0s[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(8,9) = 0.7071067811865475*m0s[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(8,10) = 0.7071067811865475*m0s[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(8,11) = 0.7071067811865475*m0s[3]*pVdim-0.7071067811865475*cE[3]; 
  BigAEM(9,8) = 0.7071067811865475*m0s[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(9,9) = 0.6324555320336759*m0s[2]*pVdim+0.7071067811865475*m0s[0]*pVdim-0.6324555320336759*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(9,10) = 0.6210590034081186*m0s[3]*pVdim+0.6324555320336759*m0s[1]*pVdim-0.6210590034081186*cE[3]-0.6324555320336759*cE[1]; 
  BigAEM(9,11) = 0.6210590034081186*m0s[2]*pVdim-0.6210590034081186*cE[2]; 
  BigAEM(10,8) = 0.7071067811865475*m0s[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(10,9) = 0.6210590034081186*m0s[3]*pVdim+0.6324555320336759*m0s[1]*pVdim-0.6210590034081186*cE[3]-0.6324555320336759*cE[1]; 
  BigAEM(10,10) = 0.4517539514526256*m0s[2]*pVdim+0.7071067811865475*m0s[0]*pVdim-0.4517539514526256*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(10,11) = 0.421637021355784*m0s[3]*pVdim+0.6210590034081186*m0s[1]*pVdim-0.421637021355784*cE[3]-0.6210590034081186*cE[1]; 
  BigAEM(11,8) = 0.7071067811865475*m0s[3]*pVdim-0.7071067811865475*cE[3]; 
  BigAEM(11,9) = 0.6210590034081186*m0s[2]*pVdim-0.6210590034081186*cE[2]; 
  BigAEM(11,10) = 0.421637021355784*m0s[3]*pVdim+0.6210590034081186*m0s[1]*pVdim-0.421637021355784*cE[3]-0.6210590034081186*cE[1]; 
  BigAEM(11,11) = 0.421637021355784*m0s[2]*pVdim+0.7071067811865475*m0s[0]*pVdim-0.421637021355784*cE[2]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<4,4>(0,4).setZero(); 
  BigAEM.block<4,4>(4,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1s[0],m1s[1],m1s[2],m1s[3],m1s[4],m1s[5],m1s[6],m1s[7],m2s[0],m2s[1],m2s[2],m2s[3]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,8,1) = xEV.segment<8>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = xEV.segment<4>(8); 
 
} 
 
void BoundaryIntegral1x2vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[4], fvmin[4]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += (1.732050807568877*fvmin[2]*dS+1.732050807568877*fvmax[2]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[1] += (fvmax[1]*dS-1.0*fvmin[1]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x2vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[10], fvmin[10]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += ((-2.23606797749979*fvmin[8]*dS)+2.23606797749979*fvmax[8]*dS+1.732050807568877*fvmin[2]*dS+1.732050807568877*fvmax[2]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[1] += (1.732050807568877*fvmin[4]*dS+1.732050807568877*fvmax[4]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[2] += (fvmax[7]*dS-1.0*fvmin[7]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x2vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[20], fvmin[20]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += (2.645751311064591*fvmin[18]*dS+2.645751311064591*fvmax[18]*dS-2.23606797749979*fvmin[8]*dS+2.23606797749979*fvmax[8]*dS+1.732050807568877*fvmin[2]*dS+1.732050807568877*fvmax[2]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[1] += ((-2.23606797749979*fvmin[12]*dS)+2.23606797749979*fvmax[12]*dS+1.732050807568877*fvmin[4]*dS+1.732050807568877*fvmax[4]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[2] += (1.732050807568877*fvmin[11]*dS+1.732050807568877*fvmax[11]*dS-1.0*fvmin[7]*dS+fvmax[7]*dS)*intFac; 
  out[3] += (fvmax[17]*dS-1.0*fvmin[17]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x2vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[4], fvmin[4]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += intFac*(1.732050807568877*fvmin[2]*dS*vmin-1.0*fvmin[0]*dS*vmin+1.732050807568877*fvmax[2]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(fvmax[1]*dS*vmax-1.0*fvmin[1]*dS*vmin); 
 
} 
 
void BoundaryIntegral1x2vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[10], fvmin[10]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += intFac*((-2.23606797749979*fvmin[8]*dS*vmin)+1.732050807568877*fvmin[2]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.23606797749979*fvmax[8]*dS*vmax+1.732050807568877*fvmax[2]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.732050807568877*fvmin[4]*dS*vmin-1.0*fvmin[1]*dS*vmin+1.732050807568877*fvmax[4]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*(fvmax[7]*dS*vmax-1.0*fvmin[7]*dS*vmin); 
 
} 
 
void BoundaryIntegral1x2vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[20], fvmin[20]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += intFac*(2.645751311064591*fvmin[18]*dS*vmin-2.23606797749979*fvmin[8]*dS*vmin+1.732050807568877*fvmin[2]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.645751311064591*fvmax[18]*dS*vmax+2.23606797749979*fvmax[8]*dS*vmax+1.732050807568877*fvmax[2]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*((-2.23606797749979*fvmin[12]*dS*vmin)+1.732050807568877*fvmin[4]*dS*vmin-1.0*fvmin[1]*dS*vmin+2.23606797749979*fvmax[12]*dS*vmax+1.732050807568877*fvmax[4]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*(1.732050807568877*fvmin[11]*dS*vmin-1.0*fvmin[7]*dS*vmin+1.732050807568877*fvmax[11]*dS*vmax+fvmax[7]*dS*vmax); 
  out[3] += intFac*(fvmax[17]*dS*vmax-1.0*fvmin[17]*dS*vmin); 
 
} 
 
void BoundaryIntegral1x2vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[4], fvmin[4]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  out[2] += (1.732050807568877*fvmin[3]*dS+1.732050807568877*fvmax[3]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[3] += (fvmax[1]*dS-1.0*fvmin[1]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x2vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[10], fvmin[10]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  out[3] += ((-2.23606797749979*fvmin[9]*dS)+2.23606797749979*fvmax[9]*dS+1.732050807568877*fvmin[3]*dS+1.732050807568877*fvmax[3]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[4] += (1.732050807568877*fvmin[5]*dS+1.732050807568877*fvmax[5]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[5] += (fvmax[7]*dS-1.0*fvmin[7]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x2vMax_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[20], fvmin[20]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  out[4] += (2.645751311064591*fvmin[19]*dS+2.645751311064591*fvmax[19]*dS-2.23606797749979*fvmin[9]*dS+2.23606797749979*fvmax[9]*dS+1.732050807568877*fvmin[3]*dS+1.732050807568877*fvmax[3]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[5] += ((-2.23606797749979*fvmin[15]*dS)+2.23606797749979*fvmax[15]*dS+1.732050807568877*fvmin[5]*dS+1.732050807568877*fvmax[5]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[6] += (1.732050807568877*fvmin[13]*dS+1.732050807568877*fvmax[13]*dS-1.0*fvmin[7]*dS+fvmax[7]*dS)*intFac; 
  out[7] += (fvmax[17]*dS-1.0*fvmin[17]*dS)*intFac; 
 
} 
 
void BoundaryIntegral1x2vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[4], fvmin[4]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  out[0] += intFac*(1.732050807568877*fvmin[3]*dS*vmin-1.0*fvmin[0]*dS*vmin+1.732050807568877*fvmax[3]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(fvmax[1]*dS*vmax-1.0*fvmin[1]*dS*vmin); 
 
} 
 
void BoundaryIntegral1x2vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[10], fvmin[10]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  out[0] += intFac*((-2.23606797749979*fvmin[9]*dS*vmin)+1.732050807568877*fvmin[3]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.23606797749979*fvmax[9]*dS*vmax+1.732050807568877*fvmax[3]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.732050807568877*fvmin[5]*dS*vmin-1.0*fvmin[1]*dS*vmin+1.732050807568877*fvmax[5]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*(fvmax[7]*dS*vmax-1.0*fvmin[7]*dS*vmin); 
 
} 
 
void BoundaryIntegral1x2vMax_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[3]:             cell length in each direciton. 
  // fvmax[20], fvmin[20]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[1]; 
 
  out[0] += intFac*(2.645751311064591*fvmin[19]*dS*vmin-2.23606797749979*fvmin[9]*dS*vmin+1.732050807568877*fvmin[3]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.645751311064591*fvmax[19]*dS*vmax+2.23606797749979*fvmax[9]*dS*vmax+1.732050807568877*fvmax[3]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*((-2.23606797749979*fvmin[15]*dS*vmin)+1.732050807568877*fvmin[5]*dS*vmin-1.0*fvmin[1]*dS*vmin+2.23606797749979*fvmax[15]*dS*vmax+1.732050807568877*fvmax[5]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*(1.732050807568877*fvmin[13]*dS*vmin-1.0*fvmin[7]*dS*vmin+1.732050807568877*fvmax[13]*dS*vmax+fvmax[7]*dS*vmax); 
  out[3] += intFac*(fvmax[17]*dS*vmax-1.0*fvmin[17]*dS*vmin); 
 
} 
 
