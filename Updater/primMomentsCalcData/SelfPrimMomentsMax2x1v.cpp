#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments2x1vMax_P1(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(6,6); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(6);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(6);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0s[0]; 
  BigAEM(0,1) = 0.5*m0s[1]; 
  BigAEM(0,2) = 0.5*m0s[2]; 
  BigAEM(1,0) = 0.5*m0s[1]; 
  BigAEM(1,1) = 0.5*m0s[0]; 
  BigAEM(2,0) = 0.5*m0s[2]; 
  BigAEM(2,2) = 0.5*m0s[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,3) = -0.5*cM[0]; 
  BigAEM(0,4) = -0.5*cM[1]; 
  BigAEM(0,5) = -0.5*cM[2]; 
  BigAEM(1,3) = -0.5*cM[1]; 
  BigAEM(1,4) = -0.5*cM[0]; 
  BigAEM(2,3) = -0.5*cM[2]; 
  BigAEM(2,5) = -0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(3,0) = 0.5*m1[0]; 
  BigAEM(3,1) = 0.5*m1[1]; 
  BigAEM(3,2) = 0.5*m1[2]; 
  BigAEM(4,0) = 0.5*m1[1]; 
  BigAEM(4,1) = 0.5*m1[0]; 
  BigAEM(5,0) = 0.5*m1[2]; 
  BigAEM(5,2) = 0.5*m1[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(3,3) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
  BigAEM(3,4) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(3,5) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(4,3) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(4,4) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
  BigAEM(5,3) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(5,5) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m1[2],m2[0],m2[1],m2[2]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,3,1) = xEV.segment<3>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV.segment<3>(3); 
 
} 
 
void SelfPrimMoments2x1vMax_P2(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
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
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(12,12); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(12);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(12);  
 
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
  BigAEM(0,6) = -0.5*cM[0]; 
  BigAEM(0,7) = -0.5*cM[1]; 
  BigAEM(0,8) = -0.5*cM[2]; 
  BigAEM(0,9) = -0.5*cM[3]; 
  BigAEM(0,10) = -0.5*cM[4]; 
  BigAEM(0,11) = -0.5*cM[5]; 
  BigAEM(1,6) = -0.5*cM[1]; 
  BigAEM(1,7) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,8) = -0.5*cM[3]; 
  BigAEM(1,9) = -0.5*cM[2]; 
  BigAEM(1,10) = -0.4472135954999579*cM[1]; 
  BigAEM(2,6) = -0.5*cM[2]; 
  BigAEM(2,7) = -0.5*cM[3]; 
  BigAEM(2,8) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,9) = -0.5*cM[1]; 
  BigAEM(2,11) = -0.4472135954999579*cM[2]; 
  BigAEM(3,6) = -0.5*cM[3]; 
  BigAEM(3,7) = -0.5*cM[2]; 
  BigAEM(3,8) = -0.5*cM[1]; 
  BigAEM(3,9) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,10) = -0.4472135954999579*cM[3]; 
  BigAEM(3,11) = -0.4472135954999579*cM[3]; 
  BigAEM(4,6) = -0.5*cM[4]; 
  BigAEM(4,7) = -0.4472135954999579*cM[1]; 
  BigAEM(4,9) = -0.4472135954999579*cM[3]; 
  BigAEM(4,10) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(5,6) = -0.5*cM[5]; 
  BigAEM(5,8) = -0.4472135954999579*cM[2]; 
  BigAEM(5,9) = -0.4472135954999579*cM[3]; 
  BigAEM(5,11) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(6,0) = 0.5*m1[0]; 
  BigAEM(6,1) = 0.5*m1[1]; 
  BigAEM(6,2) = 0.5*m1[2]; 
  BigAEM(6,3) = 0.5*m1[3]; 
  BigAEM(6,4) = 0.5*m1[4]; 
  BigAEM(6,5) = 0.5*m1[5]; 
  BigAEM(7,0) = 0.5*m1[1]; 
  BigAEM(7,1) = 0.4472135954999579*m1[4]+0.5*m1[0]; 
  BigAEM(7,2) = 0.5*m1[3]; 
  BigAEM(7,3) = 0.5*m1[2]; 
  BigAEM(7,4) = 0.4472135954999579*m1[1]; 
  BigAEM(8,0) = 0.5*m1[2]; 
  BigAEM(8,1) = 0.5*m1[3]; 
  BigAEM(8,2) = 0.4472135954999579*m1[5]+0.5*m1[0]; 
  BigAEM(8,3) = 0.5*m1[1]; 
  BigAEM(8,5) = 0.4472135954999579*m1[2]; 
  BigAEM(9,0) = 0.5*m1[3]; 
  BigAEM(9,1) = 0.5*m1[2]; 
  BigAEM(9,2) = 0.5*m1[1]; 
  BigAEM(9,3) = 0.4472135954999579*m1[5]+0.4472135954999579*m1[4]+0.5*m1[0]; 
  BigAEM(9,4) = 0.4472135954999579*m1[3]; 
  BigAEM(9,5) = 0.4472135954999579*m1[3]; 
  BigAEM(10,0) = 0.5*m1[4]; 
  BigAEM(10,1) = 0.4472135954999579*m1[1]; 
  BigAEM(10,3) = 0.4472135954999579*m1[3]; 
  BigAEM(10,4) = 0.31943828249997*m1[4]+0.5*m1[0]; 
  BigAEM(11,0) = 0.5*m1[5]; 
  BigAEM(11,2) = 0.4472135954999579*m1[2]; 
  BigAEM(11,3) = 0.4472135954999579*m1[3]; 
  BigAEM(11,5) = 0.31943828249997*m1[5]+0.5*m1[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(6,6) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
  BigAEM(6,7) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(6,8) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(6,9) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(6,10) = 0.5*m0s[4]*pVdim-0.5*cE[4]; 
  BigAEM(6,11) = 0.5*m0s[5]*pVdim-0.5*cE[5]; 
  BigAEM(7,6) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(7,7) = 0.4472135954999579*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[4]-0.5*cE[0]; 
  BigAEM(7,8) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(7,9) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(7,10) = 0.4472135954999579*m0s[1]*pVdim-0.4472135954999579*cE[1]; 
  BigAEM(8,6) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(8,7) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(8,8) = 0.4472135954999579*m0s[5]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[5]-0.5*cE[0]; 
  BigAEM(8,9) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(8,11) = 0.4472135954999579*m0s[2]*pVdim-0.4472135954999579*cE[2]; 
  BigAEM(9,6) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(9,7) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(9,8) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(9,9) = 0.4472135954999579*m0s[5]*pVdim+0.4472135954999579*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[5]-0.4472135954999579*cE[4]-0.5*cE[0]; 
  BigAEM(9,10) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(9,11) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(10,6) = 0.5*m0s[4]*pVdim-0.5*cE[4]; 
  BigAEM(10,7) = 0.4472135954999579*m0s[1]*pVdim-0.4472135954999579*cE[1]; 
  BigAEM(10,9) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(10,10) = 0.31943828249997*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.31943828249997*cE[4]-0.5*cE[0]; 
  BigAEM(11,6) = 0.5*m0s[5]*pVdim-0.5*cE[5]; 
  BigAEM(11,8) = 0.4472135954999579*m0s[2]*pVdim-0.4472135954999579*cE[2]; 
  BigAEM(11,9) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(11,11) = 0.31943828249997*m0s[5]*pVdim+0.5*m0s[0]*pVdim-0.31943828249997*cE[5]-0.5*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m1[2],m1[3],m1[4],m1[5],m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,6,1) = xEV.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,6,1) = xEV.segment<6>(6); 
 
} 
 
void SelfPrimMoments2x1vMax_P3(const int pVdim, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE: vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool m0Avg = false;
  if ((-1.322875655532295*m0[9])-1.322875655532295*m0[8]-1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if ((-1.322875655532295*m0[9])-1.322875655532295*m0[8]-1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]+1.5*m0[3]-0.8660254037844386*m0[2]-0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if ((-1.322875655532295*m0[9])+1.322875655532295*m0[8]+1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
  if ((-1.322875655532295*m0[9])+1.322875655532295*m0[8]+1.936491673103709*m0[7]-1.936491673103709*m0[6]+1.118033988749895*m0[5]+1.118033988749895*m0[4]-1.5*m0[3]-0.8660254037844386*m0[2]+0.8660254037844386*m0[1]+0.5*m0[0] < 0) { 
    m0Avg = true;
  }
 
  double m0s[10]; 
  if (m0Avg) { 
    m0s[0] = m0[0]; 
    m0s[1] = 0.0; 
    m0s[2] = 0.0; 
    m0s[3] = 0.0; 
    m0s[4] = 0.0; 
    m0s[5] = 0.0; 
    m0s[6] = 0.0; 
    m0s[7] = 0.0; 
    m0s[8] = 0.0; 
    m0s[9] = 0.0; 
  } else { 
    m0s[0] = m0[0]; 
    m0s[1] = m0[1]; 
    m0s[2] = m0[2]; 
    m0s[3] = m0[3]; 
    m0s[4] = m0[4]; 
    m0s[5] = m0[5]; 
    m0s[6] = m0[6]; 
    m0s[7] = m0[7]; 
    m0s[8] = m0[8]; 
    m0s[9] = m0[9]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd BigAEM = Eigen::MatrixXd::Zero(20,20); 
  Eigen::VectorXd bEV = Eigen::VectorXd::Zero(20);  
  Eigen::VectorXd xEV = Eigen::VectorXd::Zero(20);  
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  BigAEM(0,0) = 0.5*m0s[0]; 
  BigAEM(0,1) = 0.5*m0s[1]; 
  BigAEM(0,2) = 0.5*m0s[2]; 
  BigAEM(0,3) = 0.5*m0s[3]; 
  BigAEM(0,4) = 0.5*m0s[4]; 
  BigAEM(0,5) = 0.5*m0s[5]; 
  BigAEM(0,6) = 0.5*m0s[6]; 
  BigAEM(0,7) = 0.5*m0s[7]; 
  BigAEM(0,8) = 0.5*m0s[8]; 
  BigAEM(0,9) = 0.5*m0s[9]; 
  BigAEM(1,0) = 0.5*m0s[1]; 
  BigAEM(1,1) = 0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(1,2) = 0.5*m0s[3]; 
  BigAEM(1,3) = 0.447213595499958*m0s[6]+0.5*m0s[2]; 
  BigAEM(1,4) = 0.4391550328268398*m0s[8]+0.4472135954999579*m0s[1]; 
  BigAEM(1,5) = 0.5000000000000001*m0s[7]; 
  BigAEM(1,6) = 0.447213595499958*m0s[3]; 
  BigAEM(1,7) = 0.5000000000000001*m0s[5]; 
  BigAEM(1,8) = 0.4391550328268398*m0s[4]; 
  BigAEM(2,0) = 0.5*m0s[2]; 
  BigAEM(2,1) = 0.5*m0s[3]; 
  BigAEM(2,2) = 0.4472135954999579*m0s[5]+0.5*m0s[0]; 
  BigAEM(2,3) = 0.447213595499958*m0s[7]+0.5*m0s[1]; 
  BigAEM(2,4) = 0.5000000000000001*m0s[6]; 
  BigAEM(2,5) = 0.4391550328268398*m0s[9]+0.4472135954999579*m0s[2]; 
  BigAEM(2,6) = 0.5000000000000001*m0s[4]; 
  BigAEM(2,7) = 0.447213595499958*m0s[3]; 
  BigAEM(2,9) = 0.4391550328268398*m0s[5]; 
  BigAEM(3,0) = 0.5*m0s[3]; 
  BigAEM(3,1) = 0.447213595499958*m0s[6]+0.5*m0s[2]; 
  BigAEM(3,2) = 0.447213595499958*m0s[7]+0.5*m0s[1]; 
  BigAEM(3,3) = 0.4472135954999579*m0s[5]+0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(3,4) = 0.4472135954999579*m0s[3]; 
  BigAEM(3,5) = 0.4472135954999579*m0s[3]; 
  BigAEM(3,6) = 0.4391550328268399*m0s[8]+0.4*m0s[7]+0.447213595499958*m0s[1]; 
  BigAEM(3,7) = 0.4391550328268399*m0s[9]+0.4*m0s[6]+0.447213595499958*m0s[2]; 
  BigAEM(3,8) = 0.4391550328268399*m0s[6]; 
  BigAEM(3,9) = 0.4391550328268399*m0s[7]; 
  BigAEM(4,0) = 0.5*m0s[4]; 
  BigAEM(4,1) = 0.4391550328268398*m0s[8]+0.4472135954999579*m0s[1]; 
  BigAEM(4,2) = 0.5000000000000001*m0s[6]; 
  BigAEM(4,3) = 0.4472135954999579*m0s[3]; 
  BigAEM(4,4) = 0.31943828249997*m0s[4]+0.5*m0s[0]; 
  BigAEM(4,6) = 0.31943828249997*m0s[6]+0.5000000000000001*m0s[2]; 
  BigAEM(4,7) = 0.4472135954999579*m0s[7]; 
  BigAEM(4,8) = 0.2981423969999719*m0s[8]+0.4391550328268398*m0s[1]; 
  BigAEM(5,0) = 0.5*m0s[5]; 
  BigAEM(5,1) = 0.5000000000000001*m0s[7]; 
  BigAEM(5,2) = 0.4391550328268398*m0s[9]+0.4472135954999579*m0s[2]; 
  BigAEM(5,3) = 0.4472135954999579*m0s[3]; 
  BigAEM(5,5) = 0.31943828249997*m0s[5]+0.5*m0s[0]; 
  BigAEM(5,6) = 0.4472135954999579*m0s[6]; 
  BigAEM(5,7) = 0.31943828249997*m0s[7]+0.5000000000000001*m0s[1]; 
  BigAEM(5,9) = 0.2981423969999719*m0s[9]+0.4391550328268398*m0s[2]; 
  BigAEM(6,0) = 0.5*m0s[6]; 
  BigAEM(6,1) = 0.447213595499958*m0s[3]; 
  BigAEM(6,2) = 0.5000000000000001*m0s[4]; 
  BigAEM(6,3) = 0.4391550328268399*m0s[8]+0.4*m0s[7]+0.447213595499958*m0s[1]; 
  BigAEM(6,4) = 0.31943828249997*m0s[6]+0.5000000000000001*m0s[2]; 
  BigAEM(6,5) = 0.4472135954999579*m0s[6]; 
  BigAEM(6,6) = 0.4472135954999579*m0s[5]+0.31943828249997*m0s[4]+0.5*m0s[0]; 
  BigAEM(6,7) = 0.4*m0s[3]; 
  BigAEM(6,8) = 0.4391550328268399*m0s[3]; 
  BigAEM(7,0) = 0.5*m0s[7]; 
  BigAEM(7,1) = 0.5000000000000001*m0s[5]; 
  BigAEM(7,2) = 0.447213595499958*m0s[3]; 
  BigAEM(7,3) = 0.4391550328268399*m0s[9]+0.4*m0s[6]+0.447213595499958*m0s[2]; 
  BigAEM(7,4) = 0.4472135954999579*m0s[7]; 
  BigAEM(7,5) = 0.31943828249997*m0s[7]+0.5000000000000001*m0s[1]; 
  BigAEM(7,6) = 0.4*m0s[3]; 
  BigAEM(7,7) = 0.31943828249997*m0s[5]+0.4472135954999579*m0s[4]+0.5*m0s[0]; 
  BigAEM(7,9) = 0.4391550328268399*m0s[3]; 
  BigAEM(8,0) = 0.5*m0s[8]; 
  BigAEM(8,1) = 0.4391550328268398*m0s[4]; 
  BigAEM(8,3) = 0.4391550328268399*m0s[6]; 
  BigAEM(8,4) = 0.2981423969999719*m0s[8]+0.4391550328268398*m0s[1]; 
  BigAEM(8,6) = 0.4391550328268399*m0s[3]; 
  BigAEM(8,8) = 0.2981423969999719*m0s[4]+0.5*m0s[0]; 
  BigAEM(9,0) = 0.5*m0s[9]; 
  BigAEM(9,2) = 0.4391550328268398*m0s[5]; 
  BigAEM(9,3) = 0.4391550328268399*m0s[7]; 
  BigAEM(9,5) = 0.2981423969999719*m0s[9]+0.4391550328268398*m0s[2]; 
  BigAEM(9,7) = 0.4391550328268399*m0s[3]; 
  BigAEM(9,9) = 0.2981423969999719*m0s[5]+0.5*m0s[0]; 
 
  // ....... Block from correction to uX .......... // 
  BigAEM(0,10) = -0.5*cM[0]; 
  BigAEM(0,11) = -0.5*cM[1]; 
  BigAEM(0,12) = -0.5*cM[2]; 
  BigAEM(0,13) = -0.5*cM[3]; 
  BigAEM(0,14) = -0.5*cM[4]; 
  BigAEM(0,15) = -0.5*cM[5]; 
  BigAEM(0,16) = -0.5*cM[6]; 
  BigAEM(0,17) = -0.5*cM[7]; 
  BigAEM(0,18) = -0.5*cM[8]; 
  BigAEM(0,19) = -0.5*cM[9]; 
  BigAEM(1,10) = -0.5*cM[1]; 
  BigAEM(1,11) = (-0.4472135954999579*cM[4])-0.5*cM[0]; 
  BigAEM(1,12) = -0.5*cM[3]; 
  BigAEM(1,13) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(1,14) = (-0.4391550328268398*cM[8])-0.4472135954999579*cM[1]; 
  BigAEM(1,15) = -0.5000000000000001*cM[7]; 
  BigAEM(1,16) = -0.447213595499958*cM[3]; 
  BigAEM(1,17) = -0.5000000000000001*cM[5]; 
  BigAEM(1,18) = -0.4391550328268398*cM[4]; 
  BigAEM(2,10) = -0.5*cM[2]; 
  BigAEM(2,11) = -0.5*cM[3]; 
  BigAEM(2,12) = (-0.4472135954999579*cM[5])-0.5*cM[0]; 
  BigAEM(2,13) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(2,14) = -0.5000000000000001*cM[6]; 
  BigAEM(2,15) = (-0.4391550328268398*cM[9])-0.4472135954999579*cM[2]; 
  BigAEM(2,16) = -0.5000000000000001*cM[4]; 
  BigAEM(2,17) = -0.447213595499958*cM[3]; 
  BigAEM(2,19) = -0.4391550328268398*cM[5]; 
  BigAEM(3,10) = -0.5*cM[3]; 
  BigAEM(3,11) = (-0.447213595499958*cM[6])-0.5*cM[2]; 
  BigAEM(3,12) = (-0.447213595499958*cM[7])-0.5*cM[1]; 
  BigAEM(3,13) = (-0.4472135954999579*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(3,14) = -0.4472135954999579*cM[3]; 
  BigAEM(3,15) = -0.4472135954999579*cM[3]; 
  BigAEM(3,16) = (-0.4391550328268399*cM[8])-0.4*cM[7]-0.447213595499958*cM[1]; 
  BigAEM(3,17) = (-0.4391550328268399*cM[9])-0.4*cM[6]-0.447213595499958*cM[2]; 
  BigAEM(3,18) = -0.4391550328268399*cM[6]; 
  BigAEM(3,19) = -0.4391550328268399*cM[7]; 
  BigAEM(4,10) = -0.5*cM[4]; 
  BigAEM(4,11) = (-0.4391550328268398*cM[8])-0.4472135954999579*cM[1]; 
  BigAEM(4,12) = -0.5000000000000001*cM[6]; 
  BigAEM(4,13) = -0.4472135954999579*cM[3]; 
  BigAEM(4,14) = (-0.31943828249997*cM[4])-0.5*cM[0]; 
  BigAEM(4,16) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(4,17) = -0.4472135954999579*cM[7]; 
  BigAEM(4,18) = (-0.2981423969999719*cM[8])-0.4391550328268398*cM[1]; 
  BigAEM(5,10) = -0.5*cM[5]; 
  BigAEM(5,11) = -0.5000000000000001*cM[7]; 
  BigAEM(5,12) = (-0.4391550328268398*cM[9])-0.4472135954999579*cM[2]; 
  BigAEM(5,13) = -0.4472135954999579*cM[3]; 
  BigAEM(5,15) = (-0.31943828249997*cM[5])-0.5*cM[0]; 
  BigAEM(5,16) = -0.4472135954999579*cM[6]; 
  BigAEM(5,17) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(5,19) = (-0.2981423969999719*cM[9])-0.4391550328268398*cM[2]; 
  BigAEM(6,10) = -0.5*cM[6]; 
  BigAEM(6,11) = -0.447213595499958*cM[3]; 
  BigAEM(6,12) = -0.5000000000000001*cM[4]; 
  BigAEM(6,13) = (-0.4391550328268399*cM[8])-0.4*cM[7]-0.447213595499958*cM[1]; 
  BigAEM(6,14) = (-0.31943828249997*cM[6])-0.5000000000000001*cM[2]; 
  BigAEM(6,15) = -0.4472135954999579*cM[6]; 
  BigAEM(6,16) = (-0.4472135954999579*cM[5])-0.31943828249997*cM[4]-0.5*cM[0]; 
  BigAEM(6,17) = -0.4*cM[3]; 
  BigAEM(6,18) = -0.4391550328268399*cM[3]; 
  BigAEM(7,10) = -0.5*cM[7]; 
  BigAEM(7,11) = -0.5000000000000001*cM[5]; 
  BigAEM(7,12) = -0.447213595499958*cM[3]; 
  BigAEM(7,13) = (-0.4391550328268399*cM[9])-0.4*cM[6]-0.447213595499958*cM[2]; 
  BigAEM(7,14) = -0.4472135954999579*cM[7]; 
  BigAEM(7,15) = (-0.31943828249997*cM[7])-0.5000000000000001*cM[1]; 
  BigAEM(7,16) = -0.4*cM[3]; 
  BigAEM(7,17) = (-0.31943828249997*cM[5])-0.4472135954999579*cM[4]-0.5*cM[0]; 
  BigAEM(7,19) = -0.4391550328268399*cM[3]; 
  BigAEM(8,10) = -0.5*cM[8]; 
  BigAEM(8,11) = -0.4391550328268398*cM[4]; 
  BigAEM(8,13) = -0.4391550328268399*cM[6]; 
  BigAEM(8,14) = (-0.2981423969999719*cM[8])-0.4391550328268398*cM[1]; 
  BigAEM(8,16) = -0.4391550328268399*cM[3]; 
  BigAEM(8,18) = (-0.2981423969999719*cM[4])-0.5*cM[0]; 
  BigAEM(9,10) = -0.5*cM[9]; 
  BigAEM(9,12) = -0.4391550328268398*cM[5]; 
  BigAEM(9,13) = -0.4391550328268399*cM[7]; 
  BigAEM(9,15) = (-0.2981423969999719*cM[9])-0.4391550328268398*cM[2]; 
  BigAEM(9,17) = -0.4391550328268399*cM[3]; 
  BigAEM(9,19) = (-0.2981423969999719*cM[5])-0.5*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  BigAEM(10,0) = 0.5*m1[0]; 
  BigAEM(10,1) = 0.5*m1[1]; 
  BigAEM(10,2) = 0.5*m1[2]; 
  BigAEM(10,3) = 0.5*m1[3]; 
  BigAEM(10,4) = 0.5*m1[4]; 
  BigAEM(10,5) = 0.5*m1[5]; 
  BigAEM(10,6) = 0.5*m1[6]; 
  BigAEM(10,7) = 0.5*m1[7]; 
  BigAEM(10,8) = 0.5*m1[8]; 
  BigAEM(10,9) = 0.5*m1[9]; 
  BigAEM(11,0) = 0.5*m1[1]; 
  BigAEM(11,1) = 0.4472135954999579*m1[4]+0.5*m1[0]; 
  BigAEM(11,2) = 0.5*m1[3]; 
  BigAEM(11,3) = 0.447213595499958*m1[6]+0.5*m1[2]; 
  BigAEM(11,4) = 0.4391550328268398*m1[8]+0.4472135954999579*m1[1]; 
  BigAEM(11,5) = 0.5000000000000001*m1[7]; 
  BigAEM(11,6) = 0.447213595499958*m1[3]; 
  BigAEM(11,7) = 0.5000000000000001*m1[5]; 
  BigAEM(11,8) = 0.4391550328268398*m1[4]; 
  BigAEM(12,0) = 0.5*m1[2]; 
  BigAEM(12,1) = 0.5*m1[3]; 
  BigAEM(12,2) = 0.4472135954999579*m1[5]+0.5*m1[0]; 
  BigAEM(12,3) = 0.447213595499958*m1[7]+0.5*m1[1]; 
  BigAEM(12,4) = 0.5000000000000001*m1[6]; 
  BigAEM(12,5) = 0.4391550328268398*m1[9]+0.4472135954999579*m1[2]; 
  BigAEM(12,6) = 0.5000000000000001*m1[4]; 
  BigAEM(12,7) = 0.447213595499958*m1[3]; 
  BigAEM(12,9) = 0.4391550328268398*m1[5]; 
  BigAEM(13,0) = 0.5*m1[3]; 
  BigAEM(13,1) = 0.447213595499958*m1[6]+0.5*m1[2]; 
  BigAEM(13,2) = 0.447213595499958*m1[7]+0.5*m1[1]; 
  BigAEM(13,3) = 0.4472135954999579*m1[5]+0.4472135954999579*m1[4]+0.5*m1[0]; 
  BigAEM(13,4) = 0.4472135954999579*m1[3]; 
  BigAEM(13,5) = 0.4472135954999579*m1[3]; 
  BigAEM(13,6) = 0.4391550328268399*m1[8]+0.4*m1[7]+0.447213595499958*m1[1]; 
  BigAEM(13,7) = 0.4391550328268399*m1[9]+0.4*m1[6]+0.447213595499958*m1[2]; 
  BigAEM(13,8) = 0.4391550328268399*m1[6]; 
  BigAEM(13,9) = 0.4391550328268399*m1[7]; 
  BigAEM(14,0) = 0.5*m1[4]; 
  BigAEM(14,1) = 0.4391550328268398*m1[8]+0.4472135954999579*m1[1]; 
  BigAEM(14,2) = 0.5000000000000001*m1[6]; 
  BigAEM(14,3) = 0.4472135954999579*m1[3]; 
  BigAEM(14,4) = 0.31943828249997*m1[4]+0.5*m1[0]; 
  BigAEM(14,6) = 0.31943828249997*m1[6]+0.5000000000000001*m1[2]; 
  BigAEM(14,7) = 0.4472135954999579*m1[7]; 
  BigAEM(14,8) = 0.2981423969999719*m1[8]+0.4391550328268398*m1[1]; 
  BigAEM(15,0) = 0.5*m1[5]; 
  BigAEM(15,1) = 0.5000000000000001*m1[7]; 
  BigAEM(15,2) = 0.4391550328268398*m1[9]+0.4472135954999579*m1[2]; 
  BigAEM(15,3) = 0.4472135954999579*m1[3]; 
  BigAEM(15,5) = 0.31943828249997*m1[5]+0.5*m1[0]; 
  BigAEM(15,6) = 0.4472135954999579*m1[6]; 
  BigAEM(15,7) = 0.31943828249997*m1[7]+0.5000000000000001*m1[1]; 
  BigAEM(15,9) = 0.2981423969999719*m1[9]+0.4391550328268398*m1[2]; 
  BigAEM(16,0) = 0.5*m1[6]; 
  BigAEM(16,1) = 0.447213595499958*m1[3]; 
  BigAEM(16,2) = 0.5000000000000001*m1[4]; 
  BigAEM(16,3) = 0.4391550328268399*m1[8]+0.4*m1[7]+0.447213595499958*m1[1]; 
  BigAEM(16,4) = 0.31943828249997*m1[6]+0.5000000000000001*m1[2]; 
  BigAEM(16,5) = 0.4472135954999579*m1[6]; 
  BigAEM(16,6) = 0.4472135954999579*m1[5]+0.31943828249997*m1[4]+0.5*m1[0]; 
  BigAEM(16,7) = 0.4*m1[3]; 
  BigAEM(16,8) = 0.4391550328268399*m1[3]; 
  BigAEM(17,0) = 0.5*m1[7]; 
  BigAEM(17,1) = 0.5000000000000001*m1[5]; 
  BigAEM(17,2) = 0.447213595499958*m1[3]; 
  BigAEM(17,3) = 0.4391550328268399*m1[9]+0.4*m1[6]+0.447213595499958*m1[2]; 
  BigAEM(17,4) = 0.4472135954999579*m1[7]; 
  BigAEM(17,5) = 0.31943828249997*m1[7]+0.5000000000000001*m1[1]; 
  BigAEM(17,6) = 0.4*m1[3]; 
  BigAEM(17,7) = 0.31943828249997*m1[5]+0.4472135954999579*m1[4]+0.5*m1[0]; 
  BigAEM(17,9) = 0.4391550328268399*m1[3]; 
  BigAEM(18,0) = 0.5*m1[8]; 
  BigAEM(18,1) = 0.4391550328268398*m1[4]; 
  BigAEM(18,3) = 0.4391550328268399*m1[6]; 
  BigAEM(18,4) = 0.2981423969999719*m1[8]+0.4391550328268398*m1[1]; 
  BigAEM(18,6) = 0.4391550328268399*m1[3]; 
  BigAEM(18,8) = 0.2981423969999719*m1[4]+0.5*m1[0]; 
  BigAEM(19,0) = 0.5*m1[9]; 
  BigAEM(19,2) = 0.4391550328268398*m1[5]; 
  BigAEM(19,3) = 0.4391550328268399*m1[7]; 
  BigAEM(19,5) = 0.2981423969999719*m1[9]+0.4391550328268398*m1[2]; 
  BigAEM(19,7) = 0.4391550328268399*m1[3]; 
  BigAEM(19,9) = 0.2981423969999719*m1[5]+0.5*m1[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  BigAEM(10,10) = 0.5*m0s[0]*pVdim-0.5*cE[0]; 
  BigAEM(10,11) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(10,12) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(10,13) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(10,14) = 0.5*m0s[4]*pVdim-0.5*cE[4]; 
  BigAEM(10,15) = 0.5*m0s[5]*pVdim-0.5*cE[5]; 
  BigAEM(10,16) = 0.5*m0s[6]*pVdim-0.5*cE[6]; 
  BigAEM(10,17) = 0.5*m0s[7]*pVdim-0.5*cE[7]; 
  BigAEM(10,18) = 0.5*m0s[8]*pVdim-0.5*cE[8]; 
  BigAEM(10,19) = 0.5*m0s[9]*pVdim-0.5*cE[9]; 
  BigAEM(11,10) = 0.5*m0s[1]*pVdim-0.5*cE[1]; 
  BigAEM(11,11) = 0.4472135954999579*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[4]-0.5*cE[0]; 
  BigAEM(11,12) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(11,13) = 0.447213595499958*m0s[6]*pVdim+0.5*m0s[2]*pVdim-0.447213595499958*cE[6]-0.5*cE[2]; 
  BigAEM(11,14) = 0.4391550328268398*m0s[8]*pVdim+0.4472135954999579*m0s[1]*pVdim-0.4391550328268398*cE[8]-0.4472135954999579*cE[1]; 
  BigAEM(11,15) = 0.5000000000000001*m0s[7]*pVdim-0.5000000000000001*cE[7]; 
  BigAEM(11,16) = 0.447213595499958*m0s[3]*pVdim-0.447213595499958*cE[3]; 
  BigAEM(11,17) = 0.5000000000000001*m0s[5]*pVdim-0.5000000000000001*cE[5]; 
  BigAEM(11,18) = 0.4391550328268398*m0s[4]*pVdim-0.4391550328268398*cE[4]; 
  BigAEM(12,10) = 0.5*m0s[2]*pVdim-0.5*cE[2]; 
  BigAEM(12,11) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(12,12) = 0.4472135954999579*m0s[5]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[5]-0.5*cE[0]; 
  BigAEM(12,13) = 0.447213595499958*m0s[7]*pVdim+0.5*m0s[1]*pVdim-0.447213595499958*cE[7]-0.5*cE[1]; 
  BigAEM(12,14) = 0.5000000000000001*m0s[6]*pVdim-0.5000000000000001*cE[6]; 
  BigAEM(12,15) = 0.4391550328268398*m0s[9]*pVdim+0.4472135954999579*m0s[2]*pVdim-0.4391550328268398*cE[9]-0.4472135954999579*cE[2]; 
  BigAEM(12,16) = 0.5000000000000001*m0s[4]*pVdim-0.5000000000000001*cE[4]; 
  BigAEM(12,17) = 0.447213595499958*m0s[3]*pVdim-0.447213595499958*cE[3]; 
  BigAEM(12,19) = 0.4391550328268398*m0s[5]*pVdim-0.4391550328268398*cE[5]; 
  BigAEM(13,10) = 0.5*m0s[3]*pVdim-0.5*cE[3]; 
  BigAEM(13,11) = 0.447213595499958*m0s[6]*pVdim+0.5*m0s[2]*pVdim-0.447213595499958*cE[6]-0.5*cE[2]; 
  BigAEM(13,12) = 0.447213595499958*m0s[7]*pVdim+0.5*m0s[1]*pVdim-0.447213595499958*cE[7]-0.5*cE[1]; 
  BigAEM(13,13) = 0.4472135954999579*m0s[5]*pVdim+0.4472135954999579*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[5]-0.4472135954999579*cE[4]-0.5*cE[0]; 
  BigAEM(13,14) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(13,15) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(13,16) = 0.4391550328268399*m0s[8]*pVdim+0.4*m0s[7]*pVdim+0.447213595499958*m0s[1]*pVdim-0.4391550328268399*cE[8]-0.4*cE[7]-0.447213595499958*cE[1]; 
  BigAEM(13,17) = 0.4391550328268399*m0s[9]*pVdim+0.4*m0s[6]*pVdim+0.447213595499958*m0s[2]*pVdim-0.4391550328268399*cE[9]-0.4*cE[6]-0.447213595499958*cE[2]; 
  BigAEM(13,18) = 0.4391550328268399*m0s[6]*pVdim-0.4391550328268399*cE[6]; 
  BigAEM(13,19) = 0.4391550328268399*m0s[7]*pVdim-0.4391550328268399*cE[7]; 
  BigAEM(14,10) = 0.5*m0s[4]*pVdim-0.5*cE[4]; 
  BigAEM(14,11) = 0.4391550328268398*m0s[8]*pVdim+0.4472135954999579*m0s[1]*pVdim-0.4391550328268398*cE[8]-0.4472135954999579*cE[1]; 
  BigAEM(14,12) = 0.5000000000000001*m0s[6]*pVdim-0.5000000000000001*cE[6]; 
  BigAEM(14,13) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(14,14) = 0.31943828249997*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.31943828249997*cE[4]-0.5*cE[0]; 
  BigAEM(14,16) = 0.31943828249997*m0s[6]*pVdim+0.5000000000000001*m0s[2]*pVdim-0.31943828249997*cE[6]-0.5000000000000001*cE[2]; 
  BigAEM(14,17) = 0.4472135954999579*m0s[7]*pVdim-0.4472135954999579*cE[7]; 
  BigAEM(14,18) = 0.2981423969999719*m0s[8]*pVdim+0.4391550328268398*m0s[1]*pVdim-0.2981423969999719*cE[8]-0.4391550328268398*cE[1]; 
  BigAEM(15,10) = 0.5*m0s[5]*pVdim-0.5*cE[5]; 
  BigAEM(15,11) = 0.5000000000000001*m0s[7]*pVdim-0.5000000000000001*cE[7]; 
  BigAEM(15,12) = 0.4391550328268398*m0s[9]*pVdim+0.4472135954999579*m0s[2]*pVdim-0.4391550328268398*cE[9]-0.4472135954999579*cE[2]; 
  BigAEM(15,13) = 0.4472135954999579*m0s[3]*pVdim-0.4472135954999579*cE[3]; 
  BigAEM(15,15) = 0.31943828249997*m0s[5]*pVdim+0.5*m0s[0]*pVdim-0.31943828249997*cE[5]-0.5*cE[0]; 
  BigAEM(15,16) = 0.4472135954999579*m0s[6]*pVdim-0.4472135954999579*cE[6]; 
  BigAEM(15,17) = 0.31943828249997*m0s[7]*pVdim+0.5000000000000001*m0s[1]*pVdim-0.31943828249997*cE[7]-0.5000000000000001*cE[1]; 
  BigAEM(15,19) = 0.2981423969999719*m0s[9]*pVdim+0.4391550328268398*m0s[2]*pVdim-0.2981423969999719*cE[9]-0.4391550328268398*cE[2]; 
  BigAEM(16,10) = 0.5*m0s[6]*pVdim-0.5*cE[6]; 
  BigAEM(16,11) = 0.447213595499958*m0s[3]*pVdim-0.447213595499958*cE[3]; 
  BigAEM(16,12) = 0.5000000000000001*m0s[4]*pVdim-0.5000000000000001*cE[4]; 
  BigAEM(16,13) = 0.4391550328268399*m0s[8]*pVdim+0.4*m0s[7]*pVdim+0.447213595499958*m0s[1]*pVdim-0.4391550328268399*cE[8]-0.4*cE[7]-0.447213595499958*cE[1]; 
  BigAEM(16,14) = 0.31943828249997*m0s[6]*pVdim+0.5000000000000001*m0s[2]*pVdim-0.31943828249997*cE[6]-0.5000000000000001*cE[2]; 
  BigAEM(16,15) = 0.4472135954999579*m0s[6]*pVdim-0.4472135954999579*cE[6]; 
  BigAEM(16,16) = 0.4472135954999579*m0s[5]*pVdim+0.31943828249997*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.4472135954999579*cE[5]-0.31943828249997*cE[4]-0.5*cE[0]; 
  BigAEM(16,17) = 0.4*m0s[3]*pVdim-0.4*cE[3]; 
  BigAEM(16,18) = 0.4391550328268399*m0s[3]*pVdim-0.4391550328268399*cE[3]; 
  BigAEM(17,10) = 0.5*m0s[7]*pVdim-0.5*cE[7]; 
  BigAEM(17,11) = 0.5000000000000001*m0s[5]*pVdim-0.5000000000000001*cE[5]; 
  BigAEM(17,12) = 0.447213595499958*m0s[3]*pVdim-0.447213595499958*cE[3]; 
  BigAEM(17,13) = 0.4391550328268399*m0s[9]*pVdim+0.4*m0s[6]*pVdim+0.447213595499958*m0s[2]*pVdim-0.4391550328268399*cE[9]-0.4*cE[6]-0.447213595499958*cE[2]; 
  BigAEM(17,14) = 0.4472135954999579*m0s[7]*pVdim-0.4472135954999579*cE[7]; 
  BigAEM(17,15) = 0.31943828249997*m0s[7]*pVdim+0.5000000000000001*m0s[1]*pVdim-0.31943828249997*cE[7]-0.5000000000000001*cE[1]; 
  BigAEM(17,16) = 0.4*m0s[3]*pVdim-0.4*cE[3]; 
  BigAEM(17,17) = 0.31943828249997*m0s[5]*pVdim+0.4472135954999579*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.31943828249997*cE[5]-0.4472135954999579*cE[4]-0.5*cE[0]; 
  BigAEM(17,19) = 0.4391550328268399*m0s[3]*pVdim-0.4391550328268399*cE[3]; 
  BigAEM(18,10) = 0.5*m0s[8]*pVdim-0.5*cE[8]; 
  BigAEM(18,11) = 0.4391550328268398*m0s[4]*pVdim-0.4391550328268398*cE[4]; 
  BigAEM(18,13) = 0.4391550328268399*m0s[6]*pVdim-0.4391550328268399*cE[6]; 
  BigAEM(18,14) = 0.2981423969999719*m0s[8]*pVdim+0.4391550328268398*m0s[1]*pVdim-0.2981423969999719*cE[8]-0.4391550328268398*cE[1]; 
  BigAEM(18,16) = 0.4391550328268399*m0s[3]*pVdim-0.4391550328268399*cE[3]; 
  BigAEM(18,18) = 0.2981423969999719*m0s[4]*pVdim+0.5*m0s[0]*pVdim-0.2981423969999719*cE[4]-0.5*cE[0]; 
  BigAEM(19,10) = 0.5*m0s[9]*pVdim-0.5*cE[9]; 
  BigAEM(19,12) = 0.4391550328268398*m0s[5]*pVdim-0.4391550328268398*cE[5]; 
  BigAEM(19,13) = 0.4391550328268399*m0s[7]*pVdim-0.4391550328268399*cE[7]; 
  BigAEM(19,15) = 0.2981423969999719*m0s[9]*pVdim+0.4391550328268398*m0s[2]*pVdim-0.2981423969999719*cE[9]-0.4391550328268398*cE[2]; 
  BigAEM(19,17) = 0.4391550328268399*m0s[3]*pVdim-0.4391550328268399*cE[3]; 
  BigAEM(19,19) = 0.2981423969999719*m0s[5]*pVdim+0.5*m0s[0]*pVdim-0.2981423969999719*cE[5]-0.5*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  bEV << m1[0],m1[1],m1[2],m1[3],m1[4],m1[5],m1[6],m1[7],m1[8],m1[9],m2[0],m2[1],m2[2],m2[3],m2[4],m2[5],m2[6],m2[7],m2[8],m2[9]; 
 
  xEV = BigAEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(u,10,1) = xEV.segment<10>(0); 
 
  Eigen::Map<VectorXd>(vtSq,10,1) = xEV.segment<10>(10); 
 
} 
 
void BoundaryIntegral2x1vMax_F_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[5], fvmin[5]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  out[0] += (1.732050807568877*fvmin[3]*dS+1.732050807568877*fvmax[3]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[1] += (fvmax[1]*dS-1.0*fvmin[1]*dS)*intFac; 
  out[2] += (fvmax[2]*dS-1.0*fvmin[2]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x1vMax_F_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[15], fvmin[15]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  out[0] += ((-2.23606797749979*fvmin[13]*dS)+2.23606797749979*fvmax[13]*dS+1.732050807568877*fvmin[3]*dS+1.732050807568877*fvmax[3]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[1] += (1.732050807568877*fvmin[6]*dS+1.732050807568877*fvmax[6]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[2] += (1.732050807568877*fvmin[7]*dS+1.732050807568877*fvmax[7]*dS-1.0*fvmin[2]*dS+fvmax[2]*dS)*intFac; 
  out[3] += (fvmax[5]*dS-1.0*fvmin[5]*dS)*intFac; 
  out[4] += (fvmax[11]*dS-1.0*fvmin[11]*dS)*intFac; 
  out[5] += (fvmax[12]*dS-1.0*fvmin[12]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x1vMax_F_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[35], fvmin[35]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  out[0] += (2.645751311064591*fvmin[33]*dS+2.645751311064591*fvmax[33]*dS-2.23606797749979*fvmin[13]*dS+2.23606797749979*fvmax[13]*dS+1.732050807568877*fvmin[3]*dS+1.732050807568877*fvmax[3]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[1] += ((-2.23606797749979*fvmin[23]*dS)+2.23606797749979*fvmax[23]*dS+1.732050807568877*fvmin[6]*dS+1.732050807568877*fvmax[6]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[2] += ((-2.23606797749979*fvmin[24]*dS)+2.23606797749979*fvmax[24]*dS+1.732050807568877*fvmin[7]*dS+1.732050807568877*fvmax[7]*dS-1.0*fvmin[2]*dS+fvmax[2]*dS)*intFac; 
  out[3] += (1.732050807568877*fvmin[15]*dS+1.732050807568877*fvmax[15]*dS-1.0*fvmin[5]*dS+fvmax[5]*dS)*intFac; 
  out[4] += (1.732050807568877*fvmin[21]*dS+1.732050807568877*fvmax[21]*dS-1.0*fvmin[11]*dS+fvmax[11]*dS)*intFac; 
  out[5] += (1.732050807568877*fvmin[22]*dS+1.732050807568877*fvmax[22]*dS-1.0*fvmin[12]*dS+fvmax[12]*dS)*intFac; 
  out[6] += (fvmax[19]*dS-1.0*fvmin[19]*dS)*intFac; 
  out[7] += (fvmax[20]*dS-1.0*fvmin[20]*dS)*intFac; 
  out[8] += (fvmax[31]*dS-1.0*fvmin[31]*dS)*intFac; 
  out[9] += (fvmax[32]*dS-1.0*fvmin[32]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x1vMax_vF_VX_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[5], fvmin[5]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  out[0] += intFac*(1.732050807568877*fvmin[3]*dS*vmin-1.0*fvmin[0]*dS*vmin+1.732050807568877*fvmax[3]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(fvmax[1]*dS*vmax-1.0*fvmin[1]*dS*vmin); 
  out[2] += intFac*(fvmax[2]*dS*vmax-1.0*fvmin[2]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x1vMax_vF_VX_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[15], fvmin[15]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  out[0] += intFac*((-2.23606797749979*fvmin[13]*dS*vmin)+1.732050807568877*fvmin[3]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.23606797749979*fvmax[13]*dS*vmax+1.732050807568877*fvmax[3]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.732050807568877*fvmin[6]*dS*vmin-1.0*fvmin[1]*dS*vmin+1.732050807568877*fvmax[6]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*(1.732050807568877*fvmin[7]*dS*vmin-1.0*fvmin[2]*dS*vmin+1.732050807568877*fvmax[7]*dS*vmax+fvmax[2]*dS*vmax); 
  out[3] += intFac*(fvmax[5]*dS*vmax-1.0*fvmin[5]*dS*vmin); 
  out[4] += intFac*(fvmax[11]*dS*vmax-1.0*fvmin[11]*dS*vmin); 
  out[5] += intFac*(fvmax[12]*dS*vmax-1.0*fvmin[12]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x1vMax_vF_VX_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[35], fvmin[35]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[3]; 
 
  out[0] += intFac*(2.645751311064591*fvmin[33]*dS*vmin-2.23606797749979*fvmin[13]*dS*vmin+1.732050807568877*fvmin[3]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.645751311064591*fvmax[33]*dS*vmax+2.23606797749979*fvmax[13]*dS*vmax+1.732050807568877*fvmax[3]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*((-2.23606797749979*fvmin[23]*dS*vmin)+1.732050807568877*fvmin[6]*dS*vmin-1.0*fvmin[1]*dS*vmin+2.23606797749979*fvmax[23]*dS*vmax+1.732050807568877*fvmax[6]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*((-2.23606797749979*fvmin[24]*dS*vmin)+1.732050807568877*fvmin[7]*dS*vmin-1.0*fvmin[2]*dS*vmin+2.23606797749979*fvmax[24]*dS*vmax+1.732050807568877*fvmax[7]*dS*vmax+fvmax[2]*dS*vmax); 
  out[3] += intFac*(1.732050807568877*fvmin[15]*dS*vmin-1.0*fvmin[5]*dS*vmin+1.732050807568877*fvmax[15]*dS*vmax+fvmax[5]*dS*vmax); 
  out[4] += intFac*(1.732050807568877*fvmin[21]*dS*vmin-1.0*fvmin[11]*dS*vmin+1.732050807568877*fvmax[21]*dS*vmax+fvmax[11]*dS*vmax); 
  out[5] += intFac*(1.732050807568877*fvmin[22]*dS*vmin-1.0*fvmin[12]*dS*vmin+1.732050807568877*fvmax[22]*dS*vmax+fvmax[12]*dS*vmax); 
  out[6] += intFac*(fvmax[19]*dS*vmax-1.0*fvmin[19]*dS*vmin); 
  out[7] += intFac*(fvmax[20]*dS*vmax-1.0*fvmin[20]*dS*vmin); 
  out[8] += intFac*(fvmax[31]*dS*vmax-1.0*fvmin[31]*dS*vmin); 
  out[9] += intFac*(fvmax[32]*dS*vmax-1.0*fvmin[32]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x1vMax_F_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[5], fvmin[5]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[3] += (1.732050807568877*fvmin[4]*dS+1.732050807568877*fvmax[4]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[4] += (fvmax[1]*dS-1.0*fvmin[1]*dS)*intFac; 
  out[5] += (fvmax[2]*dS-1.0*fvmin[2]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x1vMax_F_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[15], fvmin[15]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[6] += ((-2.23606797749979*fvmin[14]*dS)+2.23606797749979*fvmax[14]*dS+1.732050807568877*fvmin[4]*dS+1.732050807568877*fvmax[4]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[7] += (1.732050807568877*fvmin[8]*dS+1.732050807568877*fvmax[8]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[8] += (1.732050807568877*fvmin[9]*dS+1.732050807568877*fvmax[9]*dS-1.0*fvmin[2]*dS+fvmax[2]*dS)*intFac; 
  out[9] += (fvmax[5]*dS-1.0*fvmin[5]*dS)*intFac; 
  out[10] += (fvmax[11]*dS-1.0*fvmin[11]*dS)*intFac; 
  out[11] += (fvmax[12]*dS-1.0*fvmin[12]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x1vMax_F_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[35], fvmin[35]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[10] += (2.645751311064591*fvmin[34]*dS+2.645751311064591*fvmax[34]*dS-2.23606797749979*fvmin[14]*dS+2.23606797749979*fvmax[14]*dS+1.732050807568877*fvmin[4]*dS+1.732050807568877*fvmax[4]*dS-1.0*fvmin[0]*dS+fvmax[0]*dS)*intFac; 
  out[11] += ((-2.23606797749979*fvmin[28]*dS)+2.23606797749979*fvmax[28]*dS+1.732050807568877*fvmin[8]*dS+1.732050807568877*fvmax[8]*dS-1.0*fvmin[1]*dS+fvmax[1]*dS)*intFac; 
  out[12] += ((-2.23606797749979*fvmin[29]*dS)+2.23606797749979*fvmax[29]*dS+1.732050807568877*fvmin[9]*dS+1.732050807568877*fvmax[9]*dS-1.0*fvmin[2]*dS+fvmax[2]*dS)*intFac; 
  out[13] += (1.732050807568877*fvmin[16]*dS+1.732050807568877*fvmax[16]*dS-1.0*fvmin[5]*dS+fvmax[5]*dS)*intFac; 
  out[14] += (1.732050807568877*fvmin[25]*dS+1.732050807568877*fvmax[25]*dS-1.0*fvmin[11]*dS+fvmax[11]*dS)*intFac; 
  out[15] += (1.732050807568877*fvmin[26]*dS+1.732050807568877*fvmax[26]*dS-1.0*fvmin[12]*dS+fvmax[12]*dS)*intFac; 
  out[16] += (fvmax[19]*dS-1.0*fvmin[19]*dS)*intFac; 
  out[17] += (fvmax[20]*dS-1.0*fvmin[20]*dS)*intFac; 
  out[18] += (fvmax[31]*dS-1.0*fvmin[31]*dS)*intFac; 
  out[19] += (fvmax[32]*dS-1.0*fvmin[32]*dS)*intFac; 
 
} 
 
void BoundaryIntegral2x1vMax_vF_VY_P1(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[5], fvmin[5]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += intFac*(1.732050807568877*fvmin[4]*dS*vmin-1.0*fvmin[0]*dS*vmin+1.732050807568877*fvmax[4]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(fvmax[1]*dS*vmax-1.0*fvmin[1]*dS*vmin); 
  out[2] += intFac*(fvmax[2]*dS*vmax-1.0*fvmin[2]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x1vMax_vF_VY_P2(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[15], fvmin[15]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += intFac*((-2.23606797749979*fvmin[14]*dS*vmin)+1.732050807568877*fvmin[4]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.23606797749979*fvmax[14]*dS*vmax+1.732050807568877*fvmax[4]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*(1.732050807568877*fvmin[8]*dS*vmin-1.0*fvmin[1]*dS*vmin+1.732050807568877*fvmax[8]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*(1.732050807568877*fvmin[9]*dS*vmin-1.0*fvmin[2]*dS*vmin+1.732050807568877*fvmax[9]*dS*vmax+fvmax[2]*dS*vmax); 
  out[3] += intFac*(fvmax[5]*dS*vmax-1.0*fvmin[5]*dS*vmin); 
  out[4] += intFac*(fvmax[11]*dS*vmax-1.0*fvmin[11]*dS*vmin); 
  out[5] += intFac*(fvmax[12]*dS*vmax-1.0*fvmin[12]*dS*vmin); 
 
} 
 
void BoundaryIntegral2x1vMax_vF_VY_P3(const double intFac, const double vmin, const double vmax, const double *dxv, const double *fvmin, const double *fvmax, double *out) 
{ 
  // intFac:             =1 for VmLBO, =2pi/m or 4pi/m for GkLBO. 
  // vmax, vmin:         maximum and minimum velocity of the velocity grid. 
  // dxv[4]:             cell length in each direciton. 
  // fvmax[35], fvmin[35]: distribution function at the velocity boundaries. 
  // out:                int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 0.5*dxv[2]; 
 
  out[0] += intFac*(2.645751311064591*fvmin[34]*dS*vmin-2.23606797749979*fvmin[14]*dS*vmin+1.732050807568877*fvmin[4]*dS*vmin-1.0*fvmin[0]*dS*vmin+2.645751311064591*fvmax[34]*dS*vmax+2.23606797749979*fvmax[14]*dS*vmax+1.732050807568877*fvmax[4]*dS*vmax+fvmax[0]*dS*vmax); 
  out[1] += intFac*((-2.23606797749979*fvmin[28]*dS*vmin)+1.732050807568877*fvmin[8]*dS*vmin-1.0*fvmin[1]*dS*vmin+2.23606797749979*fvmax[28]*dS*vmax+1.732050807568877*fvmax[8]*dS*vmax+fvmax[1]*dS*vmax); 
  out[2] += intFac*((-2.23606797749979*fvmin[29]*dS*vmin)+1.732050807568877*fvmin[9]*dS*vmin-1.0*fvmin[2]*dS*vmin+2.23606797749979*fvmax[29]*dS*vmax+1.732050807568877*fvmax[9]*dS*vmax+fvmax[2]*dS*vmax); 
  out[3] += intFac*(1.732050807568877*fvmin[16]*dS*vmin-1.0*fvmin[5]*dS*vmin+1.732050807568877*fvmax[16]*dS*vmax+fvmax[5]*dS*vmax); 
  out[4] += intFac*(1.732050807568877*fvmin[25]*dS*vmin-1.0*fvmin[11]*dS*vmin+1.732050807568877*fvmax[25]*dS*vmax+fvmax[11]*dS*vmax); 
  out[5] += intFac*(1.732050807568877*fvmin[26]*dS*vmin-1.0*fvmin[12]*dS*vmin+1.732050807568877*fvmax[26]*dS*vmax+fvmax[12]*dS*vmax); 
  out[6] += intFac*(fvmax[19]*dS*vmax-1.0*fvmin[19]*dS*vmin); 
  out[7] += intFac*(fvmax[20]*dS*vmax-1.0*fvmin[20]*dS*vmin); 
  out[8] += intFac*(fvmax[31]*dS*vmax-1.0*fvmin[31]*dS*vmin); 
  out[9] += intFac*(fvmax[32]*dS*vmax-1.0*fvmin[32]*dS*vmin); 
 
} 
 
