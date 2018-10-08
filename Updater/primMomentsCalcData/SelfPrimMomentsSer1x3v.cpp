#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments1x3vSer_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0S, const double *m1S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // m0S,m1S:  star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*m0[0]-1.224744871391589*m0[1] < 0) { 
    cellAvg = true;
  }
  if (1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[2]; 
  double m1r[6]; 
  double m2r[2]; 
  double m0Sr[2]; 
  double m1Sr[6]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = m1[2]; 
    m1r[3] = 0.0; 
    m1r[4] = m1[4]; 
    m1r[5] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = 0.0; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m1r[4] = m1[4]; 
    m1r[5] = m1[5]; 
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = m1S[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(8,8); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(8);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(8);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0r[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0r[1]; 
  BigAEM(1,0) = 0.7071067811865475*m0r[1]; 
  BigAEM(1,1) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,6) = -0.7071067811865475*cM[0]; 
  BigAEM(0,7) = -0.7071067811865475*cM[1]; 
  BigAEM(1,6) = -0.7071067811865475*cM[1]; 
  BigAEM(1,7) = -0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(6,0) = 0.7071067811865475*m1Sr[0]; 
  BigAEM(6,1) = 0.7071067811865475*m1Sr[1]; 
  BigAEM(7,0) = 0.7071067811865475*m1Sr[1]; 
  BigAEM(7,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(2,2) = 0.7071067811865475*m0r[0]; 
  BigAEM(2,3) = 0.7071067811865475*m0r[1]; 
  BigAEM(3,2) = 0.7071067811865475*m0r[1]; 
  BigAEM(3,3) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(2,6) = -0.7071067811865475*cM[2]; 
  BigAEM(2,7) = -0.7071067811865475*cM[3]; 
  BigAEM(3,6) = -0.7071067811865475*cM[3]; 
  BigAEM(3,7) = -0.7071067811865475*cM[2]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(6,2) = 0.7071067811865475*m1Sr[2]; 
  BigAEM(6,3) = 0.7071067811865475*m1Sr[3]; 
  BigAEM(7,2) = 0.7071067811865475*m1Sr[3]; 
  BigAEM(7,3) = 0.7071067811865475*m1Sr[2]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  BigAEM(4,4) = 0.7071067811865475*m0r[0]; 
  BigAEM(4,5) = 0.7071067811865475*m0r[1]; 
  BigAEM(5,4) = 0.7071067811865475*m0r[1]; 
  BigAEM(5,5) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  BigAEM(4,6) = -0.7071067811865475*cM[4]; 
  BigAEM(4,7) = -0.7071067811865475*cM[5]; 
  BigAEM(5,6) = -0.7071067811865475*cM[5]; 
  BigAEM(5,7) = -0.7071067811865475*cM[4]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  BigAEM(6,4) = 0.7071067811865475*m1Sr[4]; 
  BigAEM(6,5) = 0.7071067811865475*m1Sr[5]; 
  BigAEM(7,4) = 0.7071067811865475*m1Sr[5]; 
  BigAEM(7,5) = 0.7071067811865475*m1Sr[4]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(6,6) = 0.7071067811865475*m0Sr[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(6,7) = 0.7071067811865475*m0Sr[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(7,6) = 0.7071067811865475*m0Sr[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(7,7) = 0.7071067811865475*m0Sr[0]*pVdim-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<2,4>(0,2).setZero(); 
  BigAEM.block<4,2>(2,0).setZero(); 
  BigAEM.block<2,2>(2,4).setZero(); 
  BigAEM.block<2,2>(4,2).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m2r[0],m2r[1]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,6,1) = xEV.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = xEV.segment<2>(6); 
 
} 
 
void SelfPrimMoments1x3vSer_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0S, const double *m1S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // m0S,m1S:  star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.58113883008419*m0[2]-1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.58113883008419*m0[2]+1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[3]; 
  double m1r[9]; 
  double m2r[3]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = m1[3]; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = m1[6]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(12,12); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(12);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(12);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0r[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0r[1]; 
  BigAEM(0,2) = 0.7071067811865475*m0r[2]; 
  BigAEM(1,0) = 0.7071067811865475*m0r[1]; 
  BigAEM(1,1) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(1,2) = 0.6324555320336759*m0r[1]; 
  BigAEM(2,0) = 0.7071067811865475*m0r[2]; 
  BigAEM(2,1) = 0.6324555320336759*m0r[1]; 
  BigAEM(2,2) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,9) = -0.7071067811865475*cM[0]; 
  BigAEM(0,10) = -0.7071067811865475*cM[1]; 
  BigAEM(0,11) = -0.7071067811865475*cM[2]; 
  BigAEM(1,9) = -0.7071067811865475*cM[1]; 
  BigAEM(1,10) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(1,11) = -0.6324555320336759*cM[1]; 
  BigAEM(2,9) = -0.7071067811865475*cM[2]; 
  BigAEM(2,10) = -0.6324555320336759*cM[1]; 
  BigAEM(2,11) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(9,0) = 0.7071067811865475*m1r[0]; 
  BigAEM(9,1) = 0.7071067811865475*m1r[1]; 
  BigAEM(9,2) = 0.7071067811865475*m1r[2]; 
  BigAEM(10,0) = 0.7071067811865475*m1r[1]; 
  BigAEM(10,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  BigAEM(10,2) = 0.6324555320336759*m1r[1]; 
  BigAEM(11,0) = 0.7071067811865475*m1r[2]; 
  BigAEM(11,1) = 0.6324555320336759*m1r[1]; 
  BigAEM(11,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(3,3) = 0.7071067811865475*m0r[0]; 
  BigAEM(3,4) = 0.7071067811865475*m0r[1]; 
  BigAEM(3,5) = 0.7071067811865475*m0r[2]; 
  BigAEM(4,3) = 0.7071067811865475*m0r[1]; 
  BigAEM(4,4) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(4,5) = 0.6324555320336759*m0r[1]; 
  BigAEM(5,3) = 0.7071067811865475*m0r[2]; 
  BigAEM(5,4) = 0.6324555320336759*m0r[1]; 
  BigAEM(5,5) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(3,9) = -0.7071067811865475*cM[3]; 
  BigAEM(3,10) = -0.7071067811865475*cM[4]; 
  BigAEM(3,11) = -0.7071067811865475*cM[5]; 
  BigAEM(4,9) = -0.7071067811865475*cM[4]; 
  BigAEM(4,10) = (-0.6324555320336759*cM[5])-0.7071067811865475*cM[3]; 
  BigAEM(4,11) = -0.6324555320336759*cM[4]; 
  BigAEM(5,9) = -0.7071067811865475*cM[5]; 
  BigAEM(5,10) = -0.6324555320336759*cM[4]; 
  BigAEM(5,11) = (-0.4517539514526256*cM[5])-0.7071067811865475*cM[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(9,3) = 0.7071067811865475*m1r[3]; 
  BigAEM(9,4) = 0.7071067811865475*m1r[4]; 
  BigAEM(9,5) = 0.7071067811865475*m1r[5]; 
  BigAEM(10,3) = 0.7071067811865475*m1r[4]; 
  BigAEM(10,4) = 0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]; 
  BigAEM(10,5) = 0.6324555320336759*m1r[4]; 
  BigAEM(11,3) = 0.7071067811865475*m1r[5]; 
  BigAEM(11,4) = 0.6324555320336759*m1r[4]; 
  BigAEM(11,5) = 0.4517539514526256*m1r[5]+0.7071067811865475*m1r[3]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  BigAEM(6,6) = 0.7071067811865475*m0r[0]; 
  BigAEM(6,7) = 0.7071067811865475*m0r[1]; 
  BigAEM(6,8) = 0.7071067811865475*m0r[2]; 
  BigAEM(7,6) = 0.7071067811865475*m0r[1]; 
  BigAEM(7,7) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(7,8) = 0.6324555320336759*m0r[1]; 
  BigAEM(8,6) = 0.7071067811865475*m0r[2]; 
  BigAEM(8,7) = 0.6324555320336759*m0r[1]; 
  BigAEM(8,8) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  BigAEM(6,9) = -0.7071067811865475*cM[6]; 
  BigAEM(6,10) = -0.7071067811865475*cM[7]; 
  BigAEM(6,11) = -0.7071067811865475*cM[8]; 
  BigAEM(7,9) = -0.7071067811865475*cM[7]; 
  BigAEM(7,10) = (-0.6324555320336759*cM[8])-0.7071067811865475*cM[6]; 
  BigAEM(7,11) = -0.6324555320336759*cM[7]; 
  BigAEM(8,9) = -0.7071067811865475*cM[8]; 
  BigAEM(8,10) = -0.6324555320336759*cM[7]; 
  BigAEM(8,11) = (-0.4517539514526256*cM[8])-0.7071067811865475*cM[6]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  BigAEM(9,6) = 0.7071067811865475*m1r[6]; 
  BigAEM(9,7) = 0.7071067811865475*m1r[7]; 
  BigAEM(9,8) = 0.7071067811865475*m1r[8]; 
  BigAEM(10,6) = 0.7071067811865475*m1r[7]; 
  BigAEM(10,7) = 0.6324555320336759*m1r[8]+0.7071067811865475*m1r[6]; 
  BigAEM(10,8) = 0.6324555320336759*m1r[7]; 
  BigAEM(11,6) = 0.7071067811865475*m1r[8]; 
  BigAEM(11,7) = 0.6324555320336759*m1r[7]; 
  BigAEM(11,8) = 0.4517539514526256*m1r[8]+0.7071067811865475*m1r[6]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(9,9) = 0.7071067811865475*m0r[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(9,10) = 0.7071067811865475*m0r[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(9,11) = 0.7071067811865475*m0r[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(10,9) = 0.7071067811865475*m0r[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(10,10) = 0.6324555320336759*m0r[2]*pVdim+0.7071067811865475*m0r[0]*pVdim-0.6324555320336759*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(10,11) = 0.6324555320336759*m0r[1]*pVdim-0.6324555320336759*cE[1]; 
  BigAEM(11,9) = 0.7071067811865475*m0r[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(11,10) = 0.6324555320336759*m0r[1]*pVdim-0.6324555320336759*cE[1]; 
  BigAEM(11,11) = 0.4517539514526256*m0r[2]*pVdim+0.7071067811865475*m0r[0]*pVdim-0.4517539514526256*cE[2]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<3,6>(0,3).setZero(); 
  BigAEM.block<6,3>(3,0).setZero(); 
  BigAEM.block<3,3>(3,6).setZero(); 
  BigAEM.block<3,3>(6,3).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m2r[0],m2r[1],m2r[2]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,9,1) = xEV.segment<9>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV.segment<3>(9); 
 
} 
 
void SelfPrimMoments1x3vSer_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0S, const double *m1S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // m0S,m1S:  star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-1.870828693386971*m0[3])+1.58113883008419*m0[2]-1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
  if (1.870828693386971*m0[3]+1.58113883008419*m0[2]+1.224744871391589*m0[1]+0.7071067811865475*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[4]; 
  double m1r[12]; 
  double m2r[4]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = m1[4]; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    m1r[8] = m1[8]; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(16,16); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(16);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(16);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.7071067811865475*m0r[0]; 
  BigAEM(0,1) = 0.7071067811865475*m0r[1]; 
  BigAEM(0,2) = 0.7071067811865475*m0r[2]; 
  BigAEM(0,3) = 0.7071067811865475*m0r[3]; 
  BigAEM(1,0) = 0.7071067811865475*m0r[1]; 
  BigAEM(1,1) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(1,2) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  BigAEM(1,3) = 0.6210590034081186*m0r[2]; 
  BigAEM(2,0) = 0.7071067811865475*m0r[2]; 
  BigAEM(2,1) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  BigAEM(2,2) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(2,3) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  BigAEM(3,0) = 0.7071067811865475*m0r[3]; 
  BigAEM(3,1) = 0.6210590034081186*m0r[2]; 
  BigAEM(3,2) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  BigAEM(3,3) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,12) = -0.7071067811865475*cM[0]; 
  BigAEM(0,13) = -0.7071067811865475*cM[1]; 
  BigAEM(0,14) = -0.7071067811865475*cM[2]; 
  BigAEM(0,15) = -0.7071067811865475*cM[3]; 
  BigAEM(1,12) = -0.7071067811865475*cM[1]; 
  BigAEM(1,13) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(1,14) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  BigAEM(1,15) = -0.6210590034081186*cM[2]; 
  BigAEM(2,12) = -0.7071067811865475*cM[2]; 
  BigAEM(2,13) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  BigAEM(2,14) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
  BigAEM(2,15) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  BigAEM(3,12) = -0.7071067811865475*cM[3]; 
  BigAEM(3,13) = -0.6210590034081186*cM[2]; 
  BigAEM(3,14) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  BigAEM(3,15) = (-0.421637021355784*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(12,0) = 0.7071067811865475*m1r[0]; 
  BigAEM(12,1) = 0.7071067811865475*m1r[1]; 
  BigAEM(12,2) = 0.7071067811865475*m1r[2]; 
  BigAEM(12,3) = 0.7071067811865475*m1r[3]; 
  BigAEM(13,0) = 0.7071067811865475*m1r[1]; 
  BigAEM(13,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  BigAEM(13,2) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  BigAEM(13,3) = 0.6210590034081186*m1r[2]; 
  BigAEM(14,0) = 0.7071067811865475*m1r[2]; 
  BigAEM(14,1) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  BigAEM(14,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
  BigAEM(14,3) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  BigAEM(15,0) = 0.7071067811865475*m1r[3]; 
  BigAEM(15,1) = 0.6210590034081186*m1r[2]; 
  BigAEM(15,2) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  BigAEM(15,3) = 0.421637021355784*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  BigAEM(4,4) = 0.7071067811865475*m0r[0]; 
  BigAEM(4,5) = 0.7071067811865475*m0r[1]; 
  BigAEM(4,6) = 0.7071067811865475*m0r[2]; 
  BigAEM(4,7) = 0.7071067811865475*m0r[3]; 
  BigAEM(5,4) = 0.7071067811865475*m0r[1]; 
  BigAEM(5,5) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(5,6) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  BigAEM(5,7) = 0.6210590034081186*m0r[2]; 
  BigAEM(6,4) = 0.7071067811865475*m0r[2]; 
  BigAEM(6,5) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  BigAEM(6,6) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(6,7) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  BigAEM(7,4) = 0.7071067811865475*m0r[3]; 
  BigAEM(7,5) = 0.6210590034081186*m0r[2]; 
  BigAEM(7,6) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  BigAEM(7,7) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  BigAEM(4,12) = -0.7071067811865475*cM[4]; 
  BigAEM(4,13) = -0.7071067811865475*cM[5]; 
  BigAEM(4,14) = -0.7071067811865475*cM[6]; 
  BigAEM(4,15) = -0.7071067811865475*cM[7]; 
  BigAEM(5,12) = -0.7071067811865475*cM[5]; 
  BigAEM(5,13) = (-0.6324555320336759*cM[6])-0.7071067811865475*cM[4]; 
  BigAEM(5,14) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  BigAEM(5,15) = -0.6210590034081186*cM[6]; 
  BigAEM(6,12) = -0.7071067811865475*cM[6]; 
  BigAEM(6,13) = (-0.6210590034081186*cM[7])-0.6324555320336759*cM[5]; 
  BigAEM(6,14) = (-0.4517539514526256*cM[6])-0.7071067811865475*cM[4]; 
  BigAEM(6,15) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  BigAEM(7,12) = -0.7071067811865475*cM[7]; 
  BigAEM(7,13) = -0.6210590034081186*cM[6]; 
  BigAEM(7,14) = (-0.421637021355784*cM[7])-0.6210590034081186*cM[5]; 
  BigAEM(7,15) = (-0.421637021355784*cM[6])-0.7071067811865475*cM[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  BigAEM(12,4) = 0.7071067811865475*m1r[4]; 
  BigAEM(12,5) = 0.7071067811865475*m1r[5]; 
  BigAEM(12,6) = 0.7071067811865475*m1r[6]; 
  BigAEM(12,7) = 0.7071067811865475*m1r[7]; 
  BigAEM(13,4) = 0.7071067811865475*m1r[5]; 
  BigAEM(13,5) = 0.6324555320336759*m1r[6]+0.7071067811865475*m1r[4]; 
  BigAEM(13,6) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  BigAEM(13,7) = 0.6210590034081186*m1r[6]; 
  BigAEM(14,4) = 0.7071067811865475*m1r[6]; 
  BigAEM(14,5) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  BigAEM(14,6) = 0.4517539514526256*m1r[6]+0.7071067811865475*m1r[4]; 
  BigAEM(14,7) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  BigAEM(15,4) = 0.7071067811865475*m1r[7]; 
  BigAEM(15,5) = 0.6210590034081186*m1r[6]; 
  BigAEM(15,6) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  BigAEM(15,7) = 0.421637021355784*m1r[6]+0.7071067811865475*m1r[4]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  BigAEM(8,8) = 0.7071067811865475*m0r[0]; 
  BigAEM(8,9) = 0.7071067811865475*m0r[1]; 
  BigAEM(8,10) = 0.7071067811865475*m0r[2]; 
  BigAEM(8,11) = 0.7071067811865475*m0r[3]; 
  BigAEM(9,8) = 0.7071067811865475*m0r[1]; 
  BigAEM(9,9) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(9,10) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  BigAEM(9,11) = 0.6210590034081186*m0r[2]; 
  BigAEM(10,8) = 0.7071067811865475*m0r[2]; 
  BigAEM(10,9) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  BigAEM(10,10) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  BigAEM(10,11) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  BigAEM(11,8) = 0.7071067811865475*m0r[3]; 
  BigAEM(11,9) = 0.6210590034081186*m0r[2]; 
  BigAEM(11,10) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  BigAEM(11,11) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  BigAEM(8,12) = -0.7071067811865475*cM[8]; 
  BigAEM(8,13) = -0.7071067811865475*cM[9]; 
  BigAEM(8,14) = -0.7071067811865475*cM[10]; 
  BigAEM(8,15) = -0.7071067811865475*cM[11]; 
  BigAEM(9,12) = -0.7071067811865475*cM[9]; 
  BigAEM(9,13) = (-0.6324555320336759*cM[10])-0.7071067811865475*cM[8]; 
  BigAEM(9,14) = (-0.6210590034081186*cM[11])-0.6324555320336759*cM[9]; 
  BigAEM(9,15) = -0.6210590034081186*cM[10]; 
  BigAEM(10,12) = -0.7071067811865475*cM[10]; 
  BigAEM(10,13) = (-0.6210590034081186*cM[11])-0.6324555320336759*cM[9]; 
  BigAEM(10,14) = (-0.4517539514526256*cM[10])-0.7071067811865475*cM[8]; 
  BigAEM(10,15) = (-0.421637021355784*cM[11])-0.6210590034081186*cM[9]; 
  BigAEM(11,12) = -0.7071067811865475*cM[11]; 
  BigAEM(11,13) = -0.6210590034081186*cM[10]; 
  BigAEM(11,14) = (-0.421637021355784*cM[11])-0.6210590034081186*cM[9]; 
  BigAEM(11,15) = (-0.421637021355784*cM[10])-0.7071067811865475*cM[8]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  BigAEM(12,8) = 0.7071067811865475*m1r[8]; 
  BigAEM(12,9) = 0.7071067811865475*m1r[9]; 
  BigAEM(12,10) = 0.7071067811865475*m1r[10]; 
  BigAEM(12,11) = 0.7071067811865475*m1r[11]; 
  BigAEM(13,8) = 0.7071067811865475*m1r[9]; 
  BigAEM(13,9) = 0.6324555320336759*m1r[10]+0.7071067811865475*m1r[8]; 
  BigAEM(13,10) = 0.6210590034081186*m1r[11]+0.6324555320336759*m1r[9]; 
  BigAEM(13,11) = 0.6210590034081186*m1r[10]; 
  BigAEM(14,8) = 0.7071067811865475*m1r[10]; 
  BigAEM(14,9) = 0.6210590034081186*m1r[11]+0.6324555320336759*m1r[9]; 
  BigAEM(14,10) = 0.4517539514526256*m1r[10]+0.7071067811865475*m1r[8]; 
  BigAEM(14,11) = 0.421637021355784*m1r[11]+0.6210590034081186*m1r[9]; 
  BigAEM(15,8) = 0.7071067811865475*m1r[11]; 
  BigAEM(15,9) = 0.6210590034081186*m1r[10]; 
  BigAEM(15,10) = 0.421637021355784*m1r[11]+0.6210590034081186*m1r[9]; 
  BigAEM(15,11) = 0.421637021355784*m1r[10]+0.7071067811865475*m1r[8]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(12,12) = 0.7071067811865475*m0r[0]*pVdim-0.7071067811865475*cE[0]; 
  BigAEM(12,13) = 0.7071067811865475*m0r[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(12,14) = 0.7071067811865475*m0r[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(12,15) = 0.7071067811865475*m0r[3]*pVdim-0.7071067811865475*cE[3]; 
  BigAEM(13,12) = 0.7071067811865475*m0r[1]*pVdim-0.7071067811865475*cE[1]; 
  BigAEM(13,13) = 0.6324555320336759*m0r[2]*pVdim+0.7071067811865475*m0r[0]*pVdim-0.6324555320336759*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(13,14) = 0.6210590034081186*m0r[3]*pVdim+0.6324555320336759*m0r[1]*pVdim-0.6210590034081186*cE[3]-0.6324555320336759*cE[1]; 
  BigAEM(13,15) = 0.6210590034081186*m0r[2]*pVdim-0.6210590034081186*cE[2]; 
  BigAEM(14,12) = 0.7071067811865475*m0r[2]*pVdim-0.7071067811865475*cE[2]; 
  BigAEM(14,13) = 0.6210590034081186*m0r[3]*pVdim+0.6324555320336759*m0r[1]*pVdim-0.6210590034081186*cE[3]-0.6324555320336759*cE[1]; 
  BigAEM(14,14) = 0.4517539514526256*m0r[2]*pVdim+0.7071067811865475*m0r[0]*pVdim-0.4517539514526256*cE[2]-0.7071067811865475*cE[0]; 
  BigAEM(14,15) = 0.421637021355784*m0r[3]*pVdim+0.6210590034081186*m0r[1]*pVdim-0.421637021355784*cE[3]-0.6210590034081186*cE[1]; 
  BigAEM(15,12) = 0.7071067811865475*m0r[3]*pVdim-0.7071067811865475*cE[3]; 
  BigAEM(15,13) = 0.6210590034081186*m0r[2]*pVdim-0.6210590034081186*cE[2]; 
  BigAEM(15,14) = 0.421637021355784*m0r[3]*pVdim+0.6210590034081186*m0r[1]*pVdim-0.421637021355784*cE[3]-0.6210590034081186*cE[1]; 
  BigAEM(15,15) = 0.421637021355784*m0r[2]*pVdim+0.7071067811865475*m0r[0]*pVdim-0.421637021355784*cE[2]-0.7071067811865475*cE[0]; 
 
  // Set other entries to 0. // 
  BigAEM.block<4,8>(0,4).setZero(); 
  BigAEM.block<8,4>(4,0).setZero(); 
  BigAEM.block<4,4>(4,8).setZero(); 
  BigAEM.block<4,4>(8,4).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m1r[10],m1r[11],m2r[0],m2r[1],m2r[2],m2r[3]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,12,1) = xEV.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = xEV.segment<4>(12); 
 
} 
 
void StarMoments1x3vSer_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =1 for VmLBO, =2pi/m for GkLBO. 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[2]*dxvl[3]*intFac*(wr[1]-wl[1]); 
 
  out[0] += ((-0.8164965809277261*fr[2])+0.8164965809277261*fl[2]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[5])+0.8164965809277261*fl[5]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void StarMoments1x3vSer_VY(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =1 for VmLBO, =2pi/m for GkLBO. 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[3]*intFac*(wr[2]-wl[2]); 
 
  out[0] += ((-0.8164965809277261*fr[3])+0.8164965809277261*fl[3]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[6])+0.8164965809277261*fl[6]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void StarMoments1x3vSer_VZ(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =1 for VmLBO, =2pi/m for GkLBO. 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[4]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.25*dxvl[1]*dxvl[2]*intFac*(wr[3]-wl[3]); 
 
  out[0] += ((-0.8164965809277261*fr[4])+0.8164965809277261*fl[4]+0.7071067811865475*fr[0]+0.7071067811865475*fl[0])*dS; 
  out[1] += ((-0.8164965809277261*fr[8])+0.8164965809277261*fl[8]+0.7071067811865475*fr[1]+0.7071067811865475*fl[1])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[16], fvmin[16]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]*intFac; 
 
  out[0] += (2.449489742783178*fvmin[2]+2.449489742783178*fvmax[2]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[1] += (2.449489742783178*fvmin[5]+2.449489742783178*fvmax[5]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[48], fvmin[48]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]*intFac; 
 
  out[0] += ((-3.16227766016838*fvmin[12])+3.16227766016838*fvmax[12]+2.449489742783178*fvmin[2]+2.449489742783178*fvmax[2]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[1] += ((-3.16227766016838*fvmin[20])+3.16227766016838*fvmax[20]+2.449489742783178*fvmin[5]+2.449489742783178*fvmax[5]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
  out[2] += (2.449489742783178*fvmin[19]+2.449489742783178*fvmax[19]-1.414213562373095*fvmin[11]+1.414213562373095*fvmax[11])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[80], fvmin[80]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]*intFac; 
 
  out[0] += (3.741657386773942*fvmin[32]+3.741657386773942*fvmax[32]-3.16227766016838*fvmin[12]+3.16227766016838*fvmax[12]+2.449489742783178*fvmin[2]+2.449489742783178*fvmax[2]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[1] += (3.741657386773942*fvmin[49]+3.741657386773942*fvmax[49]-3.16227766016838*fvmin[20]+3.16227766016838*fvmax[20]+2.449489742783178*fvmin[5]+2.449489742783178*fvmax[5]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
  out[2] += (2.449489742783178*fvmin[19]+2.449489742783178*fvmax[19]-1.414213562373095*fvmin[11]+1.414213562373095*fvmax[11])*dS; 
  out[3] += (2.449489742783178*fvmin[48]+2.449489742783178*fvmax[48]-1.414213562373095*fvmin[31]+1.414213562373095*fvmax[31])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[16], fvmin[16]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]*intFac; 
 
  out[0] += dS*(2.449489742783178*fvmin[2]*vmin-1.414213562373095*fvmin[0]*vmin+2.449489742783178*fvmax[2]*vmax+1.414213562373095*fvmax[0]*vmax+1.224744871391589*dxv[1]*fvmin[2]-1.224744871391589*dxv[1]*fvmax[2]-0.7071067811865475*fvmin[0]*dxv[1]-0.7071067811865475*fvmax[0]*dxv[1]); 
  out[1] += dS*(2.449489742783178*fvmin[5]*vmin-1.414213562373095*fvmin[1]*vmin+2.449489742783178*fvmax[5]*vmax+1.414213562373095*fvmax[1]*vmax+1.224744871391589*dxv[1]*fvmin[5]-1.224744871391589*dxv[1]*fvmax[5]-0.7071067811865475*dxv[1]*fvmin[1]-0.7071067811865475*dxv[1]*fvmax[1]); 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[48], fvmin[48]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]*intFac; 
 
  out[0] += dS*((-3.16227766016838*fvmin[12]*vmin)+2.449489742783178*fvmin[2]*vmin-1.414213562373095*fvmin[0]*vmin+3.16227766016838*fvmax[12]*vmax+2.449489742783178*fvmax[2]*vmax+1.414213562373095*fvmax[0]*vmax); 
  out[1] += dS*((-3.16227766016838*fvmin[20]*vmin)+2.449489742783178*fvmin[5]*vmin-1.414213562373095*fvmin[1]*vmin+3.16227766016838*fvmax[20]*vmax+2.449489742783178*fvmax[5]*vmax+1.414213562373095*fvmax[1]*vmax); 
  out[2] += dS*(2.449489742783178*fvmin[19]*vmin-1.414213562373095*fvmin[11]*vmin+2.449489742783178*fvmax[19]*vmax+1.414213562373095*fvmax[11]*vmax); 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[80], fvmin[80]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[2]*dxv[3]*intFac; 
 
  out[0] += dS*(3.741657386773942*fvmin[32]*vmin-3.16227766016838*fvmin[12]*vmin+2.449489742783178*fvmin[2]*vmin-1.414213562373095*fvmin[0]*vmin+3.741657386773942*fvmax[32]*vmax+3.16227766016838*fvmax[12]*vmax+2.449489742783178*fvmax[2]*vmax+1.414213562373095*fvmax[0]*vmax); 
  out[1] += dS*(3.741657386773942*fvmin[49]*vmin-3.16227766016838*fvmin[20]*vmin+2.449489742783178*fvmin[5]*vmin-1.414213562373095*fvmin[1]*vmin+3.741657386773942*fvmax[49]*vmax+3.16227766016838*fvmax[20]*vmax+2.449489742783178*fvmax[5]*vmax+1.414213562373095*fvmax[1]*vmax); 
  out[2] += dS*(2.449489742783178*fvmin[19]*vmin-1.414213562373095*fvmin[11]*vmin+2.449489742783178*fvmax[19]*vmax+1.414213562373095*fvmax[11]*vmax); 
  out[3] += dS*(2.449489742783178*fvmin[48]*vmin-1.414213562373095*fvmin[31]*vmin+2.449489742783178*fvmax[48]*vmax+1.414213562373095*fvmax[31]*vmax); 
 
} 
 
void BoundaryIntegral1x3vSer_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[16], fvmin[16]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]*intFac; 
 
  out[2] += (2.449489742783178*fvmin[3]+2.449489742783178*fvmax[3]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[3] += (2.449489742783178*fvmin[6]+2.449489742783178*fvmax[6]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[48], fvmin[48]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]*intFac; 
 
  out[3] += ((-3.16227766016838*fvmin[13])+3.16227766016838*fvmax[13]+2.449489742783178*fvmin[3]+2.449489742783178*fvmax[3]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[4] += ((-3.16227766016838*fvmin[23])+3.16227766016838*fvmax[23]+2.449489742783178*fvmin[6]+2.449489742783178*fvmax[6]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
  out[5] += (2.449489742783178*fvmin[21]+2.449489742783178*fvmax[21]-1.414213562373095*fvmin[11]+1.414213562373095*fvmax[11])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[80], fvmin[80]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]*intFac; 
 
  out[4] += (3.741657386773942*fvmin[33]+3.741657386773942*fvmax[33]-3.16227766016838*fvmin[13]+3.16227766016838*fvmax[13]+2.449489742783178*fvmin[3]+2.449489742783178*fvmax[3]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[5] += (3.741657386773942*fvmin[52]+3.741657386773942*fvmax[52]-3.16227766016838*fvmin[23]+3.16227766016838*fvmax[23]+2.449489742783178*fvmin[6]+2.449489742783178*fvmax[6]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
  out[6] += (2.449489742783178*fvmin[21]+2.449489742783178*fvmax[21]-1.414213562373095*fvmin[11]+1.414213562373095*fvmax[11])*dS; 
  out[7] += (2.449489742783178*fvmin[50]+2.449489742783178*fvmax[50]-1.414213562373095*fvmin[31]+1.414213562373095*fvmax[31])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[16], fvmin[16]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]*intFac; 
 
  out[0] += dS*(2.449489742783178*fvmin[3]*vmin-1.414213562373095*fvmin[0]*vmin+2.449489742783178*fvmax[3]*vmax+1.414213562373095*fvmax[0]*vmax+1.224744871391589*dxv[2]*fvmin[3]-1.224744871391589*dxv[2]*fvmax[3]-0.7071067811865475*fvmin[0]*dxv[2]-0.7071067811865475*fvmax[0]*dxv[2]); 
  out[1] += dS*(2.449489742783178*fvmin[6]*vmin-1.414213562373095*fvmin[1]*vmin+2.449489742783178*fvmax[6]*vmax+1.414213562373095*fvmax[1]*vmax+1.224744871391589*dxv[2]*fvmin[6]-1.224744871391589*dxv[2]*fvmax[6]-0.7071067811865475*fvmin[1]*dxv[2]-0.7071067811865475*fvmax[1]*dxv[2]); 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[48], fvmin[48]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]*intFac; 
 
  out[0] += dS*((-3.16227766016838*fvmin[13]*vmin)+2.449489742783178*fvmin[3]*vmin-1.414213562373095*fvmin[0]*vmin+3.16227766016838*fvmax[13]*vmax+2.449489742783178*fvmax[3]*vmax+1.414213562373095*fvmax[0]*vmax); 
  out[1] += dS*((-3.16227766016838*fvmin[23]*vmin)+2.449489742783178*fvmin[6]*vmin-1.414213562373095*fvmin[1]*vmin+3.16227766016838*fvmax[23]*vmax+2.449489742783178*fvmax[6]*vmax+1.414213562373095*fvmax[1]*vmax); 
  out[2] += dS*(2.449489742783178*fvmin[21]*vmin-1.414213562373095*fvmin[11]*vmin+2.449489742783178*fvmax[21]*vmax+1.414213562373095*fvmax[11]*vmax); 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[80], fvmin[80]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[3]*intFac; 
 
  out[0] += dS*(3.741657386773942*fvmin[33]*vmin-3.16227766016838*fvmin[13]*vmin+2.449489742783178*fvmin[3]*vmin-1.414213562373095*fvmin[0]*vmin+3.741657386773942*fvmax[33]*vmax+3.16227766016838*fvmax[13]*vmax+2.449489742783178*fvmax[3]*vmax+1.414213562373095*fvmax[0]*vmax); 
  out[1] += dS*(3.741657386773942*fvmin[52]*vmin-3.16227766016838*fvmin[23]*vmin+2.449489742783178*fvmin[6]*vmin-1.414213562373095*fvmin[1]*vmin+3.741657386773942*fvmax[52]*vmax+3.16227766016838*fvmax[23]*vmax+2.449489742783178*fvmax[6]*vmax+1.414213562373095*fvmax[1]*vmax); 
  out[2] += dS*(2.449489742783178*fvmin[21]*vmin-1.414213562373095*fvmin[11]*vmin+2.449489742783178*fvmax[21]*vmax+1.414213562373095*fvmax[11]*vmax); 
  out[3] += dS*(2.449489742783178*fvmin[50]*vmin-1.414213562373095*fvmin[31]*vmin+2.449489742783178*fvmax[50]*vmax+1.414213562373095*fvmax[31]*vmax); 
 
} 
 
void BoundaryIntegral1x3vSer_F_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[16], fvmin[16]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]*intFac; 
 
  out[4] += (2.449489742783178*fvmin[4]+2.449489742783178*fvmax[4]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[5] += (2.449489742783178*fvmin[8]+2.449489742783178*fvmax[8]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_F_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[48], fvmin[48]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]*intFac; 
 
  out[6] += ((-3.16227766016838*fvmin[14])+3.16227766016838*fvmax[14]+2.449489742783178*fvmin[4]+2.449489742783178*fvmax[4]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[7] += ((-3.16227766016838*fvmin[28])+3.16227766016838*fvmax[28]+2.449489742783178*fvmin[8]+2.449489742783178*fvmax[8]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
  out[8] += (2.449489742783178*fvmin[25]+2.449489742783178*fvmax[25]-1.414213562373095*fvmin[11]+1.414213562373095*fvmax[11])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_F_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[80], fvmin[80]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]*intFac; 
 
  out[8] += (3.741657386773942*fvmin[34]+3.741657386773942*fvmax[34]-3.16227766016838*fvmin[14]+3.16227766016838*fvmax[14]+2.449489742783178*fvmin[4]+2.449489742783178*fvmax[4]-1.414213562373095*fvmin[0]+1.414213562373095*fvmax[0])*dS; 
  out[9] += (3.741657386773942*fvmin[57]+3.741657386773942*fvmax[57]-3.16227766016838*fvmin[28]+3.16227766016838*fvmax[28]+2.449489742783178*fvmin[8]+2.449489742783178*fvmax[8]-1.414213562373095*fvmin[1]+1.414213562373095*fvmax[1])*dS; 
  out[10] += (2.449489742783178*fvmin[25]+2.449489742783178*fvmax[25]-1.414213562373095*fvmin[11]+1.414213562373095*fvmax[11])*dS; 
  out[11] += (2.449489742783178*fvmin[54]+2.449489742783178*fvmax[54]-1.414213562373095*fvmin[31]+1.414213562373095*fvmax[31])*dS; 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VZ_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[16], fvmin[16]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]*intFac; 
 
  out[0] += dS*(2.449489742783178*fvmin[4]*vmin-1.414213562373095*fvmin[0]*vmin+2.449489742783178*fvmax[4]*vmax+1.414213562373095*fvmax[0]*vmax+1.224744871391589*dxv[3]*fvmin[4]-1.224744871391589*dxv[3]*fvmax[4]-0.7071067811865475*fvmin[0]*dxv[3]-0.7071067811865475*fvmax[0]*dxv[3]); 
  out[1] += dS*(2.449489742783178*fvmin[8]*vmin-1.414213562373095*fvmin[1]*vmin+2.449489742783178*fvmax[8]*vmax+1.414213562373095*fvmax[1]*vmax+1.224744871391589*dxv[3]*fvmin[8]-1.224744871391589*dxv[3]*fvmax[8]-0.7071067811865475*fvmin[1]*dxv[3]-0.7071067811865475*fvmax[1]*dxv[3]); 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VZ_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[48], fvmin[48]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]*intFac; 
 
  out[0] += dS*((-3.16227766016838*fvmin[14]*vmin)+2.449489742783178*fvmin[4]*vmin-1.414213562373095*fvmin[0]*vmin+3.16227766016838*fvmax[14]*vmax+2.449489742783178*fvmax[4]*vmax+1.414213562373095*fvmax[0]*vmax); 
  out[1] += dS*((-3.16227766016838*fvmin[28]*vmin)+2.449489742783178*fvmin[8]*vmin-1.414213562373095*fvmin[1]*vmin+3.16227766016838*fvmax[28]*vmax+2.449489742783178*fvmax[8]*vmax+1.414213562373095*fvmax[1]*vmax); 
  out[2] += dS*(2.449489742783178*fvmin[25]*vmin-1.414213562373095*fvmin[11]*vmin+2.449489742783178*fvmax[25]*vmax+1.414213562373095*fvmax[11]*vmax); 
 
} 
 
void BoundaryIntegral1x3vSer_vF_VZ_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[80], fvmin[80]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.25*dxv[1]*dxv[2]*intFac; 
 
  out[0] += dS*(3.741657386773942*fvmin[34]*vmin-3.16227766016838*fvmin[14]*vmin+2.449489742783178*fvmin[4]*vmin-1.414213562373095*fvmin[0]*vmin+3.741657386773942*fvmax[34]*vmax+3.16227766016838*fvmax[14]*vmax+2.449489742783178*fvmax[4]*vmax+1.414213562373095*fvmax[0]*vmax); 
  out[1] += dS*(3.741657386773942*fvmin[57]*vmin-3.16227766016838*fvmin[28]*vmin+2.449489742783178*fvmin[8]*vmin-1.414213562373095*fvmin[1]*vmin+3.741657386773942*fvmax[57]*vmax+3.16227766016838*fvmax[28]*vmax+2.449489742783178*fvmax[8]*vmax+1.414213562373095*fvmax[1]*vmax); 
  out[2] += dS*(2.449489742783178*fvmin[25]*vmin-1.414213562373095*fvmin[11]*vmin+2.449489742783178*fvmax[25]*vmax+1.414213562373095*fvmax[11]*vmax); 
  out[3] += dS*(2.449489742783178*fvmin[54]*vmin-1.414213562373095*fvmin[31]*vmin+2.449489742783178*fvmax[54]*vmax+1.414213562373095*fvmax[31]*vmax); 
 
} 
 
