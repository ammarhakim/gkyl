#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments3x1vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0S, const double *m1S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // m0S,m1S:  star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0[3])-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[4]; 
  double m1r[4]; 
  double m2r[4]; 
  double m0Sr[4]; 
  double m1Sr[4]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m0Sr[3] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
    m1Sr[3] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m0Sr[3] = m0S[3]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(8,8); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(8);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(8);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.3535533905932737*m0r[0]; 
  BigAEM(0,1) = 0.3535533905932737*m0r[1]; 
  BigAEM(0,2) = 0.3535533905932737*m0r[2]; 
  BigAEM(0,3) = 0.3535533905932737*m0r[3]; 
  BigAEM(1,0) = 0.3535533905932737*m0r[1]; 
  BigAEM(1,1) = 0.3535533905932737*m0r[0]; 
  BigAEM(2,0) = 0.3535533905932737*m0r[2]; 
  BigAEM(2,2) = 0.3535533905932737*m0r[0]; 
  BigAEM(3,0) = 0.3535533905932737*m0r[3]; 
  BigAEM(3,3) = 0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,4) = -0.3535533905932737*cM[0]; 
  BigAEM(0,5) = -0.3535533905932737*cM[1]; 
  BigAEM(0,6) = -0.3535533905932737*cM[2]; 
  BigAEM(0,7) = -0.3535533905932737*cM[3]; 
  BigAEM(1,4) = -0.3535533905932737*cM[1]; 
  BigAEM(1,5) = -0.3535533905932737*cM[0]; 
  BigAEM(2,4) = -0.3535533905932737*cM[2]; 
  BigAEM(2,6) = -0.3535533905932737*cM[0]; 
  BigAEM(3,4) = -0.3535533905932737*cM[3]; 
  BigAEM(3,7) = -0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(4,0) = 0.3535533905932737*m1Sr[0]; 
  BigAEM(4,1) = 0.3535533905932737*m1Sr[1]; 
  BigAEM(4,2) = 0.3535533905932737*m1Sr[2]; 
  BigAEM(4,3) = 0.3535533905932737*m1Sr[3]; 
  BigAEM(5,0) = 0.3535533905932737*m1Sr[1]; 
  BigAEM(5,1) = 0.3535533905932737*m1Sr[0]; 
  BigAEM(6,0) = 0.3535533905932737*m1Sr[2]; 
  BigAEM(6,2) = 0.3535533905932737*m1Sr[0]; 
  BigAEM(7,0) = 0.3535533905932737*m1Sr[3]; 
  BigAEM(7,3) = 0.3535533905932737*m1Sr[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(4,4) = 0.3535533905932737*m0Sr[0]*pVdim-0.3535533905932737*cE[0]; 
  BigAEM(4,5) = 0.3535533905932737*m0Sr[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(4,6) = 0.3535533905932737*m0Sr[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(4,7) = 0.3535533905932737*m0Sr[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(5,4) = 0.3535533905932737*m0Sr[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(5,5) = 0.3535533905932737*m0Sr[0]*pVdim-0.3535533905932737*cE[0]; 
  BigAEM(6,4) = 0.3535533905932737*m0Sr[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(6,6) = 0.3535533905932737*m0Sr[0]*pVdim-0.3535533905932737*cE[0]; 
  BigAEM(7,4) = 0.3535533905932737*m0Sr[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(7,7) = 0.3535533905932737*m0Sr[0]*pVdim-0.3535533905932737*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m2r[0],m2r[1],m2r[2],m2r[3]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,4,1) = xEV.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = xEV.segment<4>(4); 
 
} 
 
void SelfPrimMoments3x1vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *m0S, const double *m1S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // m0S,m1S:  star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]+1.060660171779821*m0[5]+1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]-0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0[9]+0.7905694150420947*m0[8]+0.7905694150420947*m0[7]+1.060660171779821*m0[6]-1.060660171779821*m0[5]-1.060660171779821*m0[4]-0.6123724356957944*m0[3]-0.6123724356957944*m0[2]+0.6123724356957944*m0[1]+0.3535533905932737*m0[0] < 0) { 
    cellAvg = true;
  }
 
  double m0r[10]; 
  double m1r[10]; 
  double m2r[10]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m0r[4] = 0.0; 
    m0r[5] = 0.0; 
    m0r[6] = 0.0; 
    m0r[7] = 0.0; 
    m0r[8] = 0.0; 
    m0r[9] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    m1r[9] = 0.0; 
    m2r[0] = m2[0]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    m2r[4] = 0.0; 
    m2r[5] = 0.0; 
    m2r[6] = 0.0; 
    m2r[7] = 0.0; 
    m2r[8] = 0.0; 
    m2r[9] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m0r[2] = m0[2]; 
    m0r[3] = m0[3]; 
    m0r[4] = m0[4]; 
    m0r[5] = m0[5]; 
    m0r[6] = m0[6]; 
    m0r[7] = m0[7]; 
    m0r[8] = m0[8]; 
    m0r[9] = m0[9]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
    m2r[4] = m2[4]; 
    m2r[5] = m2[5]; 
    m2r[6] = m2[6]; 
    m2r[7] = m2[7]; 
    m2r[8] = m2[8]; 
    m2r[9] = m2[9]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(20,20); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(20);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(20);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.3535533905932737*m0r[0]; 
  BigAEM(0,1) = 0.3535533905932737*m0r[1]; 
  BigAEM(0,2) = 0.3535533905932737*m0r[2]; 
  BigAEM(0,3) = 0.3535533905932737*m0r[3]; 
  BigAEM(0,4) = 0.3535533905932737*m0r[4]; 
  BigAEM(0,5) = 0.3535533905932737*m0r[5]; 
  BigAEM(0,6) = 0.3535533905932737*m0r[6]; 
  BigAEM(0,7) = 0.3535533905932737*m0r[7]; 
  BigAEM(0,8) = 0.3535533905932737*m0r[8]; 
  BigAEM(0,9) = 0.3535533905932737*m0r[9]; 
  BigAEM(1,0) = 0.3535533905932737*m0r[1]; 
  BigAEM(1,1) = 0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(1,2) = 0.3535533905932737*m0r[4]; 
  BigAEM(1,3) = 0.3535533905932737*m0r[5]; 
  BigAEM(1,4) = 0.3535533905932737*m0r[2]; 
  BigAEM(1,5) = 0.3535533905932737*m0r[3]; 
  BigAEM(1,7) = 0.3162277660168379*m0r[1]; 
  BigAEM(2,0) = 0.3535533905932737*m0r[2]; 
  BigAEM(2,1) = 0.3535533905932737*m0r[4]; 
  BigAEM(2,2) = 0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  BigAEM(2,3) = 0.3535533905932737*m0r[6]; 
  BigAEM(2,4) = 0.3535533905932737*m0r[1]; 
  BigAEM(2,6) = 0.3535533905932737*m0r[3]; 
  BigAEM(2,8) = 0.3162277660168379*m0r[2]; 
  BigAEM(3,0) = 0.3535533905932737*m0r[3]; 
  BigAEM(3,1) = 0.3535533905932737*m0r[5]; 
  BigAEM(3,2) = 0.3535533905932737*m0r[6]; 
  BigAEM(3,3) = 0.3162277660168379*m0r[9]+0.3535533905932737*m0r[0]; 
  BigAEM(3,5) = 0.3535533905932737*m0r[1]; 
  BigAEM(3,6) = 0.3535533905932737*m0r[2]; 
  BigAEM(3,9) = 0.3162277660168379*m0r[3]; 
  BigAEM(4,0) = 0.3535533905932737*m0r[4]; 
  BigAEM(4,1) = 0.3535533905932737*m0r[2]; 
  BigAEM(4,2) = 0.3535533905932737*m0r[1]; 
  BigAEM(4,4) = 0.3162277660168379*m0r[8]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(4,5) = 0.3535533905932737*m0r[6]; 
  BigAEM(4,6) = 0.3535533905932737*m0r[5]; 
  BigAEM(4,7) = 0.3162277660168379*m0r[4]; 
  BigAEM(4,8) = 0.3162277660168379*m0r[4]; 
  BigAEM(5,0) = 0.3535533905932737*m0r[5]; 
  BigAEM(5,1) = 0.3535533905932737*m0r[3]; 
  BigAEM(5,3) = 0.3535533905932737*m0r[1]; 
  BigAEM(5,4) = 0.3535533905932737*m0r[6]; 
  BigAEM(5,5) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(5,6) = 0.3535533905932737*m0r[4]; 
  BigAEM(5,7) = 0.3162277660168379*m0r[5]; 
  BigAEM(5,9) = 0.3162277660168379*m0r[5]; 
  BigAEM(6,0) = 0.3535533905932737*m0r[6]; 
  BigAEM(6,2) = 0.3535533905932737*m0r[3]; 
  BigAEM(6,3) = 0.3535533905932737*m0r[2]; 
  BigAEM(6,4) = 0.3535533905932737*m0r[5]; 
  BigAEM(6,5) = 0.3535533905932737*m0r[4]; 
  BigAEM(6,6) = 0.3162277660168379*m0r[9]+0.3162277660168379*m0r[8]+0.3535533905932737*m0r[0]; 
  BigAEM(6,8) = 0.3162277660168379*m0r[6]; 
  BigAEM(6,9) = 0.3162277660168379*m0r[6]; 
  BigAEM(7,0) = 0.3535533905932737*m0r[7]; 
  BigAEM(7,1) = 0.3162277660168379*m0r[1]; 
  BigAEM(7,4) = 0.3162277660168379*m0r[4]; 
  BigAEM(7,5) = 0.3162277660168379*m0r[5]; 
  BigAEM(7,7) = 0.2258769757263128*m0r[7]+0.3535533905932737*m0r[0]; 
  BigAEM(8,0) = 0.3535533905932737*m0r[8]; 
  BigAEM(8,2) = 0.3162277660168379*m0r[2]; 
  BigAEM(8,4) = 0.3162277660168379*m0r[4]; 
  BigAEM(8,6) = 0.3162277660168379*m0r[6]; 
  BigAEM(8,8) = 0.2258769757263128*m0r[8]+0.3535533905932737*m0r[0]; 
  BigAEM(9,0) = 0.3535533905932737*m0r[9]; 
  BigAEM(9,3) = 0.3162277660168379*m0r[3]; 
  BigAEM(9,5) = 0.3162277660168379*m0r[5]; 
  BigAEM(9,6) = 0.3162277660168379*m0r[6]; 
  BigAEM(9,9) = 0.2258769757263128*m0r[9]+0.3535533905932737*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,10) = -0.3535533905932737*cM[0]; 
  BigAEM(0,11) = -0.3535533905932737*cM[1]; 
  BigAEM(0,12) = -0.3535533905932737*cM[2]; 
  BigAEM(0,13) = -0.3535533905932737*cM[3]; 
  BigAEM(0,14) = -0.3535533905932737*cM[4]; 
  BigAEM(0,15) = -0.3535533905932737*cM[5]; 
  BigAEM(0,16) = -0.3535533905932737*cM[6]; 
  BigAEM(0,17) = -0.3535533905932737*cM[7]; 
  BigAEM(0,18) = -0.3535533905932737*cM[8]; 
  BigAEM(0,19) = -0.3535533905932737*cM[9]; 
  BigAEM(1,10) = -0.3535533905932737*cM[1]; 
  BigAEM(1,11) = (-0.3162277660168379*cM[7])-0.3535533905932737*cM[0]; 
  BigAEM(1,12) = -0.3535533905932737*cM[4]; 
  BigAEM(1,13) = -0.3535533905932737*cM[5]; 
  BigAEM(1,14) = -0.3535533905932737*cM[2]; 
  BigAEM(1,15) = -0.3535533905932737*cM[3]; 
  BigAEM(1,17) = -0.3162277660168379*cM[1]; 
  BigAEM(2,10) = -0.3535533905932737*cM[2]; 
  BigAEM(2,11) = -0.3535533905932737*cM[4]; 
  BigAEM(2,12) = (-0.3162277660168379*cM[8])-0.3535533905932737*cM[0]; 
  BigAEM(2,13) = -0.3535533905932737*cM[6]; 
  BigAEM(2,14) = -0.3535533905932737*cM[1]; 
  BigAEM(2,16) = -0.3535533905932737*cM[3]; 
  BigAEM(2,18) = -0.3162277660168379*cM[2]; 
  BigAEM(3,10) = -0.3535533905932737*cM[3]; 
  BigAEM(3,11) = -0.3535533905932737*cM[5]; 
  BigAEM(3,12) = -0.3535533905932737*cM[6]; 
  BigAEM(3,13) = (-0.3162277660168379*cM[9])-0.3535533905932737*cM[0]; 
  BigAEM(3,15) = -0.3535533905932737*cM[1]; 
  BigAEM(3,16) = -0.3535533905932737*cM[2]; 
  BigAEM(3,19) = -0.3162277660168379*cM[3]; 
  BigAEM(4,10) = -0.3535533905932737*cM[4]; 
  BigAEM(4,11) = -0.3535533905932737*cM[2]; 
  BigAEM(4,12) = -0.3535533905932737*cM[1]; 
  BigAEM(4,14) = (-0.3162277660168379*cM[8])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  BigAEM(4,15) = -0.3535533905932737*cM[6]; 
  BigAEM(4,16) = -0.3535533905932737*cM[5]; 
  BigAEM(4,17) = -0.3162277660168379*cM[4]; 
  BigAEM(4,18) = -0.3162277660168379*cM[4]; 
  BigAEM(5,10) = -0.3535533905932737*cM[5]; 
  BigAEM(5,11) = -0.3535533905932737*cM[3]; 
  BigAEM(5,13) = -0.3535533905932737*cM[1]; 
  BigAEM(5,14) = -0.3535533905932737*cM[6]; 
  BigAEM(5,15) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[7]-0.3535533905932737*cM[0]; 
  BigAEM(5,16) = -0.3535533905932737*cM[4]; 
  BigAEM(5,17) = -0.3162277660168379*cM[5]; 
  BigAEM(5,19) = -0.3162277660168379*cM[5]; 
  BigAEM(6,10) = -0.3535533905932737*cM[6]; 
  BigAEM(6,12) = -0.3535533905932737*cM[3]; 
  BigAEM(6,13) = -0.3535533905932737*cM[2]; 
  BigAEM(6,14) = -0.3535533905932737*cM[5]; 
  BigAEM(6,15) = -0.3535533905932737*cM[4]; 
  BigAEM(6,16) = (-0.3162277660168379*cM[9])-0.3162277660168379*cM[8]-0.3535533905932737*cM[0]; 
  BigAEM(6,18) = -0.3162277660168379*cM[6]; 
  BigAEM(6,19) = -0.3162277660168379*cM[6]; 
  BigAEM(7,10) = -0.3535533905932737*cM[7]; 
  BigAEM(7,11) = -0.3162277660168379*cM[1]; 
  BigAEM(7,14) = -0.3162277660168379*cM[4]; 
  BigAEM(7,15) = -0.3162277660168379*cM[5]; 
  BigAEM(7,17) = (-0.2258769757263128*cM[7])-0.3535533905932737*cM[0]; 
  BigAEM(8,10) = -0.3535533905932737*cM[8]; 
  BigAEM(8,12) = -0.3162277660168379*cM[2]; 
  BigAEM(8,14) = -0.3162277660168379*cM[4]; 
  BigAEM(8,16) = -0.3162277660168379*cM[6]; 
  BigAEM(8,18) = (-0.2258769757263128*cM[8])-0.3535533905932737*cM[0]; 
  BigAEM(9,10) = -0.3535533905932737*cM[9]; 
  BigAEM(9,13) = -0.3162277660168379*cM[3]; 
  BigAEM(9,15) = -0.3162277660168379*cM[5]; 
  BigAEM(9,16) = -0.3162277660168379*cM[6]; 
  BigAEM(9,19) = (-0.2258769757263128*cM[9])-0.3535533905932737*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(10,0) = 0.3535533905932737*m1r[0]; 
  BigAEM(10,1) = 0.3535533905932737*m1r[1]; 
  BigAEM(10,2) = 0.3535533905932737*m1r[2]; 
  BigAEM(10,3) = 0.3535533905932737*m1r[3]; 
  BigAEM(10,4) = 0.3535533905932737*m1r[4]; 
  BigAEM(10,5) = 0.3535533905932737*m1r[5]; 
  BigAEM(10,6) = 0.3535533905932737*m1r[6]; 
  BigAEM(10,7) = 0.3535533905932737*m1r[7]; 
  BigAEM(10,8) = 0.3535533905932737*m1r[8]; 
  BigAEM(10,9) = 0.3535533905932737*m1r[9]; 
  BigAEM(11,0) = 0.3535533905932737*m1r[1]; 
  BigAEM(11,1) = 0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(11,2) = 0.3535533905932737*m1r[4]; 
  BigAEM(11,3) = 0.3535533905932737*m1r[5]; 
  BigAEM(11,4) = 0.3535533905932737*m1r[2]; 
  BigAEM(11,5) = 0.3535533905932737*m1r[3]; 
  BigAEM(11,7) = 0.3162277660168379*m1r[1]; 
  BigAEM(12,0) = 0.3535533905932737*m1r[2]; 
  BigAEM(12,1) = 0.3535533905932737*m1r[4]; 
  BigAEM(12,2) = 0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  BigAEM(12,3) = 0.3535533905932737*m1r[6]; 
  BigAEM(12,4) = 0.3535533905932737*m1r[1]; 
  BigAEM(12,6) = 0.3535533905932737*m1r[3]; 
  BigAEM(12,8) = 0.3162277660168379*m1r[2]; 
  BigAEM(13,0) = 0.3535533905932737*m1r[3]; 
  BigAEM(13,1) = 0.3535533905932737*m1r[5]; 
  BigAEM(13,2) = 0.3535533905932737*m1r[6]; 
  BigAEM(13,3) = 0.3162277660168379*m1r[9]+0.3535533905932737*m1r[0]; 
  BigAEM(13,5) = 0.3535533905932737*m1r[1]; 
  BigAEM(13,6) = 0.3535533905932737*m1r[2]; 
  BigAEM(13,9) = 0.3162277660168379*m1r[3]; 
  BigAEM(14,0) = 0.3535533905932737*m1r[4]; 
  BigAEM(14,1) = 0.3535533905932737*m1r[2]; 
  BigAEM(14,2) = 0.3535533905932737*m1r[1]; 
  BigAEM(14,4) = 0.3162277660168379*m1r[8]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(14,5) = 0.3535533905932737*m1r[6]; 
  BigAEM(14,6) = 0.3535533905932737*m1r[5]; 
  BigAEM(14,7) = 0.3162277660168379*m1r[4]; 
  BigAEM(14,8) = 0.3162277660168379*m1r[4]; 
  BigAEM(15,0) = 0.3535533905932737*m1r[5]; 
  BigAEM(15,1) = 0.3535533905932737*m1r[3]; 
  BigAEM(15,3) = 0.3535533905932737*m1r[1]; 
  BigAEM(15,4) = 0.3535533905932737*m1r[6]; 
  BigAEM(15,5) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(15,6) = 0.3535533905932737*m1r[4]; 
  BigAEM(15,7) = 0.3162277660168379*m1r[5]; 
  BigAEM(15,9) = 0.3162277660168379*m1r[5]; 
  BigAEM(16,0) = 0.3535533905932737*m1r[6]; 
  BigAEM(16,2) = 0.3535533905932737*m1r[3]; 
  BigAEM(16,3) = 0.3535533905932737*m1r[2]; 
  BigAEM(16,4) = 0.3535533905932737*m1r[5]; 
  BigAEM(16,5) = 0.3535533905932737*m1r[4]; 
  BigAEM(16,6) = 0.3162277660168379*m1r[9]+0.3162277660168379*m1r[8]+0.3535533905932737*m1r[0]; 
  BigAEM(16,8) = 0.3162277660168379*m1r[6]; 
  BigAEM(16,9) = 0.3162277660168379*m1r[6]; 
  BigAEM(17,0) = 0.3535533905932737*m1r[7]; 
  BigAEM(17,1) = 0.3162277660168379*m1r[1]; 
  BigAEM(17,4) = 0.3162277660168379*m1r[4]; 
  BigAEM(17,5) = 0.3162277660168379*m1r[5]; 
  BigAEM(17,7) = 0.2258769757263128*m1r[7]+0.3535533905932737*m1r[0]; 
  BigAEM(18,0) = 0.3535533905932737*m1r[8]; 
  BigAEM(18,2) = 0.3162277660168379*m1r[2]; 
  BigAEM(18,4) = 0.3162277660168379*m1r[4]; 
  BigAEM(18,6) = 0.3162277660168379*m1r[6]; 
  BigAEM(18,8) = 0.2258769757263128*m1r[8]+0.3535533905932737*m1r[0]; 
  BigAEM(19,0) = 0.3535533905932737*m1r[9]; 
  BigAEM(19,3) = 0.3162277660168379*m1r[3]; 
  BigAEM(19,5) = 0.3162277660168379*m1r[5]; 
  BigAEM(19,6) = 0.3162277660168379*m1r[6]; 
  BigAEM(19,9) = 0.2258769757263128*m1r[9]+0.3535533905932737*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(10,10) = 0.3535533905932737*m0r[0]*pVdim-0.3535533905932737*cE[0]; 
  BigAEM(10,11) = 0.3535533905932737*m0r[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(10,12) = 0.3535533905932737*m0r[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(10,13) = 0.3535533905932737*m0r[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(10,14) = 0.3535533905932737*m0r[4]*pVdim-0.3535533905932737*cE[4]; 
  BigAEM(10,15) = 0.3535533905932737*m0r[5]*pVdim-0.3535533905932737*cE[5]; 
  BigAEM(10,16) = 0.3535533905932737*m0r[6]*pVdim-0.3535533905932737*cE[6]; 
  BigAEM(10,17) = 0.3535533905932737*m0r[7]*pVdim-0.3535533905932737*cE[7]; 
  BigAEM(10,18) = 0.3535533905932737*m0r[8]*pVdim-0.3535533905932737*cE[8]; 
  BigAEM(10,19) = 0.3535533905932737*m0r[9]*pVdim-0.3535533905932737*cE[9]; 
  BigAEM(11,10) = 0.3535533905932737*m0r[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(11,11) = 0.3162277660168379*m0r[7]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.3162277660168379*cE[7]-0.3535533905932737*cE[0]; 
  BigAEM(11,12) = 0.3535533905932737*m0r[4]*pVdim-0.3535533905932737*cE[4]; 
  BigAEM(11,13) = 0.3535533905932737*m0r[5]*pVdim-0.3535533905932737*cE[5]; 
  BigAEM(11,14) = 0.3535533905932737*m0r[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(11,15) = 0.3535533905932737*m0r[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(11,17) = 0.3162277660168379*m0r[1]*pVdim-0.3162277660168379*cE[1]; 
  BigAEM(12,10) = 0.3535533905932737*m0r[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(12,11) = 0.3535533905932737*m0r[4]*pVdim-0.3535533905932737*cE[4]; 
  BigAEM(12,12) = 0.3162277660168379*m0r[8]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.3162277660168379*cE[8]-0.3535533905932737*cE[0]; 
  BigAEM(12,13) = 0.3535533905932737*m0r[6]*pVdim-0.3535533905932737*cE[6]; 
  BigAEM(12,14) = 0.3535533905932737*m0r[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(12,16) = 0.3535533905932737*m0r[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(12,18) = 0.3162277660168379*m0r[2]*pVdim-0.3162277660168379*cE[2]; 
  BigAEM(13,10) = 0.3535533905932737*m0r[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(13,11) = 0.3535533905932737*m0r[5]*pVdim-0.3535533905932737*cE[5]; 
  BigAEM(13,12) = 0.3535533905932737*m0r[6]*pVdim-0.3535533905932737*cE[6]; 
  BigAEM(13,13) = 0.3162277660168379*m0r[9]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.3162277660168379*cE[9]-0.3535533905932737*cE[0]; 
  BigAEM(13,15) = 0.3535533905932737*m0r[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(13,16) = 0.3535533905932737*m0r[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(13,19) = 0.3162277660168379*m0r[3]*pVdim-0.3162277660168379*cE[3]; 
  BigAEM(14,10) = 0.3535533905932737*m0r[4]*pVdim-0.3535533905932737*cE[4]; 
  BigAEM(14,11) = 0.3535533905932737*m0r[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(14,12) = 0.3535533905932737*m0r[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(14,14) = 0.3162277660168379*m0r[8]*pVdim+0.3162277660168379*m0r[7]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.3162277660168379*cE[8]-0.3162277660168379*cE[7]-0.3535533905932737*cE[0]; 
  BigAEM(14,15) = 0.3535533905932737*m0r[6]*pVdim-0.3535533905932737*cE[6]; 
  BigAEM(14,16) = 0.3535533905932737*m0r[5]*pVdim-0.3535533905932737*cE[5]; 
  BigAEM(14,17) = 0.3162277660168379*m0r[4]*pVdim-0.3162277660168379*cE[4]; 
  BigAEM(14,18) = 0.3162277660168379*m0r[4]*pVdim-0.3162277660168379*cE[4]; 
  BigAEM(15,10) = 0.3535533905932737*m0r[5]*pVdim-0.3535533905932737*cE[5]; 
  BigAEM(15,11) = 0.3535533905932737*m0r[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(15,13) = 0.3535533905932737*m0r[1]*pVdim-0.3535533905932737*cE[1]; 
  BigAEM(15,14) = 0.3535533905932737*m0r[6]*pVdim-0.3535533905932737*cE[6]; 
  BigAEM(15,15) = 0.3162277660168379*m0r[9]*pVdim+0.3162277660168379*m0r[7]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.3162277660168379*cE[9]-0.3162277660168379*cE[7]-0.3535533905932737*cE[0]; 
  BigAEM(15,16) = 0.3535533905932737*m0r[4]*pVdim-0.3535533905932737*cE[4]; 
  BigAEM(15,17) = 0.3162277660168379*m0r[5]*pVdim-0.3162277660168379*cE[5]; 
  BigAEM(15,19) = 0.3162277660168379*m0r[5]*pVdim-0.3162277660168379*cE[5]; 
  BigAEM(16,10) = 0.3535533905932737*m0r[6]*pVdim-0.3535533905932737*cE[6]; 
  BigAEM(16,12) = 0.3535533905932737*m0r[3]*pVdim-0.3535533905932737*cE[3]; 
  BigAEM(16,13) = 0.3535533905932737*m0r[2]*pVdim-0.3535533905932737*cE[2]; 
  BigAEM(16,14) = 0.3535533905932737*m0r[5]*pVdim-0.3535533905932737*cE[5]; 
  BigAEM(16,15) = 0.3535533905932737*m0r[4]*pVdim-0.3535533905932737*cE[4]; 
  BigAEM(16,16) = 0.3162277660168379*m0r[9]*pVdim+0.3162277660168379*m0r[8]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.3162277660168379*cE[9]-0.3162277660168379*cE[8]-0.3535533905932737*cE[0]; 
  BigAEM(16,18) = 0.3162277660168379*m0r[6]*pVdim-0.3162277660168379*cE[6]; 
  BigAEM(16,19) = 0.3162277660168379*m0r[6]*pVdim-0.3162277660168379*cE[6]; 
  BigAEM(17,10) = 0.3535533905932737*m0r[7]*pVdim-0.3535533905932737*cE[7]; 
  BigAEM(17,11) = 0.3162277660168379*m0r[1]*pVdim-0.3162277660168379*cE[1]; 
  BigAEM(17,14) = 0.3162277660168379*m0r[4]*pVdim-0.3162277660168379*cE[4]; 
  BigAEM(17,15) = 0.3162277660168379*m0r[5]*pVdim-0.3162277660168379*cE[5]; 
  BigAEM(17,17) = 0.2258769757263128*m0r[7]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.2258769757263128*cE[7]-0.3535533905932737*cE[0]; 
  BigAEM(18,10) = 0.3535533905932737*m0r[8]*pVdim-0.3535533905932737*cE[8]; 
  BigAEM(18,12) = 0.3162277660168379*m0r[2]*pVdim-0.3162277660168379*cE[2]; 
  BigAEM(18,14) = 0.3162277660168379*m0r[4]*pVdim-0.3162277660168379*cE[4]; 
  BigAEM(18,16) = 0.3162277660168379*m0r[6]*pVdim-0.3162277660168379*cE[6]; 
  BigAEM(18,18) = 0.2258769757263128*m0r[8]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.2258769757263128*cE[8]-0.3535533905932737*cE[0]; 
  BigAEM(19,10) = 0.3535533905932737*m0r[9]*pVdim-0.3535533905932737*cE[9]; 
  BigAEM(19,13) = 0.3162277660168379*m0r[3]*pVdim-0.3162277660168379*cE[3]; 
  BigAEM(19,15) = 0.3162277660168379*m0r[5]*pVdim-0.3162277660168379*cE[5]; 
  BigAEM(19,16) = 0.3162277660168379*m0r[6]*pVdim-0.3162277660168379*cE[6]; 
  BigAEM(19,19) = 0.2258769757263128*m0r[9]*pVdim+0.3535533905932737*m0r[0]*pVdim-0.2258769757263128*cE[9]-0.3535533905932737*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m1r[8],m1r[9],m2r[0],m2r[1],m2r[2],m2r[3],m2r[4],m2r[5],m2r[6],m2r[7],m2r[8],m2r[9]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,10,1) = xEV.segment<10>(0); 
 
  Eigen::Map<VectorXd>(vtSq,10,1) = xEV.segment<10>(10); 
 
} 
 
void StarMoments3x1vMax_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =1 for VmLBO, =2pi/m for GkLBO. 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[4]*intFac*(wr[3]-wl[3]); 
 
  out[0] += ((-0.5773502691896258*fr[4])+0.5773502691896258*fl[4]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += (0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += (0.5*fr[2]+0.5*fl[2])*dS; 
  out[3] += (0.5*fr[3]+0.5*fl[3])*dS; 
 
} 
 
void StarMoments3x1vMax_VY(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =1 for VmLBO, =2pi/m for GkLBO. 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[5]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 0.5*dxvl[3]*intFac*(wr[4]-wl[4]); 
 
  out[0] += ((-0.5773502691896258*fr[5])+0.5773502691896258*fl[5]+0.5*fr[0]+0.5*fl[0])*dS; 
  out[1] += (0.5*fr[1]+0.5*fl[1])*dS; 
  out[2] += (0.5*fr[2]+0.5*fl[2])*dS; 
  out[3] += (0.5*fr[3]+0.5*fl[3])*dS; 
 
} 
 
void BoundaryIntegral3x2vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  out[0] += (1.732050807568877*fvmin[4]+1.732050807568877*fvmax[4]-1.0*fvmin[0]+fvmax[0])*dS; 
  out[1] += (fvmax[1]-1.0*fvmin[1])*dS; 
  out[2] += (fvmax[2]-1.0*fvmin[2])*dS; 
  out[3] += (fvmax[3]-1.0*fvmin[3])*dS; 
 
} 
 
void BoundaryIntegral3x2vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  out[0] += ((-2.23606797749979*fvmin[19])+2.23606797749979*fvmax[19]+1.732050807568877*fvmin[4]+1.732050807568877*fvmax[4]-1.0*fvmin[0]+fvmax[0])*dS; 
  out[1] += (1.732050807568877*fvmin[9]+1.732050807568877*fvmax[9]-1.0*fvmin[1]+fvmax[1])*dS; 
  out[2] += (1.732050807568877*fvmin[10]+1.732050807568877*fvmax[10]-1.0*fvmin[2]+fvmax[2])*dS; 
  out[3] += (1.732050807568877*fvmin[11]+1.732050807568877*fvmax[11]-1.0*fvmin[3]+fvmax[3])*dS; 
  out[4] += (fvmax[6]-1.0*fvmin[6])*dS; 
  out[5] += (fvmax[7]-1.0*fvmin[7])*dS; 
  out[6] += (fvmax[8]-1.0*fvmin[8])*dS; 
  out[7] += (fvmax[16]-1.0*fvmin[16])*dS; 
  out[8] += (fvmax[17]-1.0*fvmin[17])*dS; 
  out[9] += (fvmax[18]-1.0*fvmin[18])*dS; 
 
} 
 
void BoundaryIntegral3x2vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  out[0] += dS*(1.732050807568877*fvmin[4]*vmin-1.0*fvmin[0]*vmin+1.732050807568877*fvmax[4]*vmax+fvmax[0]*vmax+0.8660254037844386*dxv[3]*fvmin[4]-0.8660254037844386*dxv[3]*fvmax[4]-0.5*fvmin[0]*dxv[3]-0.5*fvmax[0]*dxv[3]); 
  out[1] += dS*((-1.0*fvmin[1]*vmin)+fvmax[1]*vmax-0.5*fvmin[1]*dxv[3]-0.5*fvmax[1]*dxv[3]); 
  out[2] += dS*((-1.0*fvmin[2]*vmin)+fvmax[2]*vmax-0.5*fvmin[2]*dxv[3]-0.5*fvmax[2]*dxv[3]); 
  out[3] += dS*((-1.0*fvmin[3]*vmin)+fvmax[3]*vmax-0.5*dxv[3]*fvmin[3]-0.5*dxv[3]*fvmax[3]); 
 
} 
 
void BoundaryIntegral3x2vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[4]*intFac; 
 
  out[0] += dS*((-2.23606797749979*fvmin[19]*vmin)+1.732050807568877*fvmin[4]*vmin-1.0*fvmin[0]*vmin+2.23606797749979*fvmax[19]*vmax+1.732050807568877*fvmax[4]*vmax+fvmax[0]*vmax); 
  out[1] += dS*(1.732050807568877*fvmin[9]*vmin-1.0*fvmin[1]*vmin+1.732050807568877*fvmax[9]*vmax+fvmax[1]*vmax); 
  out[2] += dS*(1.732050807568877*fvmin[10]*vmin-1.0*fvmin[2]*vmin+1.732050807568877*fvmax[10]*vmax+fvmax[2]*vmax); 
  out[3] += dS*(1.732050807568877*fvmin[11]*vmin-1.0*fvmin[3]*vmin+1.732050807568877*fvmax[11]*vmax+fvmax[3]*vmax); 
  out[4] += dS*(fvmax[6]*vmax-1.0*fvmin[6]*vmin); 
  out[5] += dS*(fvmax[7]*vmax-1.0*fvmin[7]*vmin); 
  out[6] += dS*(fvmax[8]*vmax-1.0*fvmin[8]*vmin); 
  out[7] += dS*(fvmax[16]*vmax-1.0*fvmin[16]*vmin); 
  out[8] += dS*(fvmax[17]*vmax-1.0*fvmin[17]*vmin); 
  out[9] += dS*(fvmax[18]*vmax-1.0*fvmin[18]*vmin); 
 
} 
 
void BoundaryIntegral3x2vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  out[4] += (1.732050807568877*fvmin[5]+1.732050807568877*fvmax[5]-1.0*fvmin[0]+fvmax[0])*dS; 
  out[5] += (fvmax[1]-1.0*fvmin[1])*dS; 
  out[6] += (fvmax[2]-1.0*fvmin[2])*dS; 
  out[7] += (fvmax[3]-1.0*fvmin[3])*dS; 
 
} 
 
void BoundaryIntegral3x2vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  out[10] += ((-2.23606797749979*fvmin[20])+2.23606797749979*fvmax[20]+1.732050807568877*fvmin[5]+1.732050807568877*fvmax[5]-1.0*fvmin[0]+fvmax[0])*dS; 
  out[11] += (1.732050807568877*fvmin[12]+1.732050807568877*fvmax[12]-1.0*fvmin[1]+fvmax[1])*dS; 
  out[12] += (1.732050807568877*fvmin[13]+1.732050807568877*fvmax[13]-1.0*fvmin[2]+fvmax[2])*dS; 
  out[13] += (1.732050807568877*fvmin[14]+1.732050807568877*fvmax[14]-1.0*fvmin[3]+fvmax[3])*dS; 
  out[14] += (fvmax[6]-1.0*fvmin[6])*dS; 
  out[15] += (fvmax[7]-1.0*fvmin[7])*dS; 
  out[16] += (fvmax[8]-1.0*fvmin[8])*dS; 
  out[17] += (fvmax[16]-1.0*fvmin[16])*dS; 
  out[18] += (fvmax[17]-1.0*fvmin[17])*dS; 
  out[19] += (fvmax[18]-1.0*fvmin[18])*dS; 
 
} 
 
void BoundaryIntegral3x2vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[6], fvmin[6]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  out[0] += dS*(1.732050807568877*fvmin[5]*vmin-1.0*fvmin[0]*vmin+1.732050807568877*fvmax[5]*vmax+fvmax[0]*vmax+0.8660254037844386*dxv[4]*fvmin[5]-0.8660254037844386*dxv[4]*fvmax[5]-0.5*fvmin[0]*dxv[4]-0.5*fvmax[0]*dxv[4]); 
  out[1] += dS*((-1.0*fvmin[1]*vmin)+fvmax[1]*vmax-0.5*fvmin[1]*dxv[4]-0.5*fvmax[1]*dxv[4]); 
  out[2] += dS*((-1.0*fvmin[2]*vmin)+fvmax[2]*vmax-0.5*fvmin[2]*dxv[4]-0.5*fvmax[2]*dxv[4]); 
  out[3] += dS*((-1.0*fvmin[3]*vmin)+fvmax[3]*vmax-0.5*fvmin[3]*dxv[4]-0.5*fvmax[3]*dxv[4]); 
 
} 
 
void BoundaryIntegral3x2vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[5]:             cell length in each direciton. 
  // fvmax[21], fvmin[21]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]*intFac; 
 
  out[0] += dS*((-2.23606797749979*fvmin[20]*vmin)+1.732050807568877*fvmin[5]*vmin-1.0*fvmin[0]*vmin+2.23606797749979*fvmax[20]*vmax+1.732050807568877*fvmax[5]*vmax+fvmax[0]*vmax); 
  out[1] += dS*(1.732050807568877*fvmin[12]*vmin-1.0*fvmin[1]*vmin+1.732050807568877*fvmax[12]*vmax+fvmax[1]*vmax); 
  out[2] += dS*(1.732050807568877*fvmin[13]*vmin-1.0*fvmin[2]*vmin+1.732050807568877*fvmax[13]*vmax+fvmax[2]*vmax); 
  out[3] += dS*(1.732050807568877*fvmin[14]*vmin-1.0*fvmin[3]*vmin+1.732050807568877*fvmax[14]*vmax+fvmax[3]*vmax); 
  out[4] += dS*(fvmax[6]*vmax-1.0*fvmin[6]*vmin); 
  out[5] += dS*(fvmax[7]*vmax-1.0*fvmin[7]*vmin); 
  out[6] += dS*(fvmax[8]*vmax-1.0*fvmin[8]*vmin); 
  out[7] += dS*(fvmax[16]*vmax-1.0*fvmin[16]*vmin); 
  out[8] += dS*(fvmax[17]*vmax-1.0*fvmin[17]*vmin); 
  out[9] += dS*(fvmax[18]*vmax-1.0*fvmin[18]*vmin); 
 
} 
 
