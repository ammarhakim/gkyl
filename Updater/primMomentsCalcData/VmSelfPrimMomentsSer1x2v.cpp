#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmSelfPrimMoments1x2vSer_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(2.449489742783178*m0[1]-1.414213562373095*m0[0]) < 0) cellAvg = true; 
  if (0.5*(2.449489742783178*m0[1]+1.414213562373095*m0[0]) < 0) cellAvg = true; 
 
  double m0r[2]; 
  double m1r[4]; 
  double cMr[4]; 
  double cEr[2]; 
  double m0Sr[2]; 
  double m1Sr[4]; 
  double m2Sr[2]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    m1r[2] = m1[2]; 
    m1r[3] = 0.0; 
    cMr[2] = cM[2]; 
    cMr[3] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    m1r[2] = m1[2]; 
    m1r[3] = m1[3]; 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cMr[3] = cM[3]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(6,6); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,4) = -0.7071067811865475*cMr[0]; 
  data->AEM_S(0,5) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,4) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,5) = -0.7071067811865475*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(4,0) = 0.7071067811865475*m1Sr[0]; 
  data->AEM_S(4,1) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(5,0) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(5,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(2,3) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,3) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(2,4) = -0.7071067811865475*cMr[2]; 
  data->AEM_S(2,5) = -0.7071067811865475*cMr[3]; 
  data->AEM_S(3,4) = -0.7071067811865475*cMr[3]; 
  data->AEM_S(3,5) = -0.7071067811865475*cMr[2]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(4,2) = 0.7071067811865475*m1Sr[2]; 
  data->AEM_S(4,3) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(5,2) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(5,3) = 0.7071067811865475*m1Sr[2]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(4,4) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(4,5) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(5,4) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(5,5) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cEr[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<2,2>(0,2).setZero(); 
  data->AEM_S.block<2,2>(2,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m2Sr[0],m2Sr[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,4,1) = data->u_S.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = data->u_S.segment<2>(4); 
 
} 
 
void VmSelfPrimMoments1x2vSer_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*(2.23606797749979*m0[2]-1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*m0[2]+1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
 
  double m0r[3]; 
  double m1r[6]; 
  double cMr[6]; 
  double cEr[3]; 
  double m2r[3]; 
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
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
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
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cMr[3] = cM[3]; 
    cMr[4] = cM[4]; 
    cMr[5] = cM[5]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(9,9); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(0,2) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(1,2) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(2,0) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(2,1) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(2,2) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,6) = -0.7071067811865475*cMr[0]; 
  data->AEM_S(0,7) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(0,8) = -0.7071067811865475*cMr[2]; 
  data->AEM_S(1,6) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,7) = (-0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]; 
  data->AEM_S(1,8) = -0.6324555320336759*cMr[1]; 
  data->AEM_S(2,6) = -0.7071067811865475*cMr[2]; 
  data->AEM_S(2,7) = -0.6324555320336759*cMr[1]; 
  data->AEM_S(2,8) = (-0.4517539514526256*cMr[2])-0.7071067811865475*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(6,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(6,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(6,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(7,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(7,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(7,2) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(8,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(8,1) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(8,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(3,3) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(3,4) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,5) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(4,3) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(4,4) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(4,5) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(5,3) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(5,4) = 0.6324555320336759*m0r[1]; 
  data->AEM_S(5,5) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(3,6) = -0.7071067811865475*cMr[3]; 
  data->AEM_S(3,7) = -0.7071067811865475*cMr[4]; 
  data->AEM_S(3,8) = -0.7071067811865475*cMr[5]; 
  data->AEM_S(4,6) = -0.7071067811865475*cMr[4]; 
  data->AEM_S(4,7) = (-0.6324555320336759*cMr[5])-0.7071067811865475*cMr[3]; 
  data->AEM_S(4,8) = -0.6324555320336759*cMr[4]; 
  data->AEM_S(5,6) = -0.7071067811865475*cMr[5]; 
  data->AEM_S(5,7) = -0.6324555320336759*cMr[4]; 
  data->AEM_S(5,8) = (-0.4517539514526256*cMr[5])-0.7071067811865475*cMr[3]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(6,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(6,4) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(6,5) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(7,3) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(7,4) = 0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]; 
  data->AEM_S(7,5) = 0.6324555320336759*m1r[4]; 
  data->AEM_S(8,3) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(8,4) = 0.6324555320336759*m1r[4]; 
  data->AEM_S(8,5) = 0.4517539514526256*m1r[5]+0.7071067811865475*m1r[3]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(6,6) = 1.414213562373095*m0r[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(6,7) = 1.414213562373095*m0r[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(6,8) = 1.414213562373095*m0r[2]-0.7071067811865475*cEr[2]; 
  data->AEM_S(7,6) = 1.414213562373095*m0r[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(7,7) = 1.264911064067352*m0r[2]-0.6324555320336759*cEr[2]+1.414213562373095*m0r[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(7,8) = 1.264911064067352*m0r[1]-0.6324555320336759*cEr[1]; 
  data->AEM_S(8,6) = 1.414213562373095*m0r[2]-0.7071067811865475*cEr[2]; 
  data->AEM_S(8,7) = 1.264911064067352*m0r[1]-0.6324555320336759*cEr[1]; 
  data->AEM_S(8,8) = 0.9035079029052515*m0r[2]-0.4517539514526256*cEr[2]+1.414213562373095*m0r[0]-0.7071067811865475*cEr[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,3>(0,3).setZero(); 
  data->AEM_S.block<3,3>(3,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m2r[0],m2r[1],m2r[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,6,1) = data->u_S.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(6); 
 
} 
 
void VmSelfPrimMoments1x2vSer_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (-0.5*(3.741657386773942*m0[3]-3.16227766016838*m0[2]+2.449489742783178*m0[1]-1.414213562373095*m0[0]) < 0) cellAvg = true; 
  if (0.5*(3.741657386773942*m0[3]+3.16227766016838*m0[2]+2.449489742783178*m0[1]+1.414213562373095*m0[0]) < 0) cellAvg = true; 
 
  double m0r[4]; 
  double m1r[8]; 
  double cMr[8]; 
  double cEr[4]; 
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
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    m1r[4] = m1[4]; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    cMr[4] = cM[4]; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(12,12); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(0,2) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(0,3) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(1,2) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(1,3) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(2,0) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(2,1) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(2,2) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(2,3) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(3,0) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(3,1) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(3,2) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(3,3) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,8) = -0.7071067811865475*cMr[0]; 
  data->AEM_S(0,9) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(0,10) = -0.7071067811865475*cMr[2]; 
  data->AEM_S(0,11) = -0.7071067811865475*cMr[3]; 
  data->AEM_S(1,8) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,9) = (-0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]; 
  data->AEM_S(1,10) = (-0.6210590034081186*cMr[3])-0.6324555320336759*cMr[1]; 
  data->AEM_S(1,11) = -0.6210590034081186*cMr[2]; 
  data->AEM_S(2,8) = -0.7071067811865475*cMr[2]; 
  data->AEM_S(2,9) = (-0.6210590034081186*cMr[3])-0.6324555320336759*cMr[1]; 
  data->AEM_S(2,10) = (-0.4517539514526256*cMr[2])-0.7071067811865475*cMr[0]; 
  data->AEM_S(2,11) = (-0.421637021355784*cMr[3])-0.6210590034081186*cMr[1]; 
  data->AEM_S(3,8) = -0.7071067811865475*cMr[3]; 
  data->AEM_S(3,9) = -0.6210590034081186*cMr[2]; 
  data->AEM_S(3,10) = (-0.421637021355784*cMr[3])-0.6210590034081186*cMr[1]; 
  data->AEM_S(3,11) = (-0.421637021355784*cMr[2])-0.7071067811865475*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(8,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(8,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(8,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(8,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(9,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(9,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(9,2) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(9,3) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(10,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(10,1) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(10,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(10,3) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(11,0) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(11,1) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(11,2) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(11,3) = 0.421637021355784*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(4,4) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(4,5) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(4,6) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(4,7) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(5,4) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(5,5) = 0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(5,6) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(5,7) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(6,4) = 0.7071067811865475*m0r[2]; 
  data->AEM_S(6,5) = 0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]; 
  data->AEM_S(6,6) = 0.4517539514526256*m0r[2]+0.7071067811865475*m0r[0]; 
  data->AEM_S(6,7) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(7,4) = 0.7071067811865475*m0r[3]; 
  data->AEM_S(7,5) = 0.6210590034081186*m0r[2]; 
  data->AEM_S(7,6) = 0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]; 
  data->AEM_S(7,7) = 0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(4,8) = -0.7071067811865475*cMr[4]; 
  data->AEM_S(4,9) = -0.7071067811865475*cMr[5]; 
  data->AEM_S(4,10) = -0.7071067811865475*cMr[6]; 
  data->AEM_S(4,11) = -0.7071067811865475*cMr[7]; 
  data->AEM_S(5,8) = -0.7071067811865475*cMr[5]; 
  data->AEM_S(5,9) = (-0.6324555320336759*cMr[6])-0.7071067811865475*cMr[4]; 
  data->AEM_S(5,10) = (-0.6210590034081186*cMr[7])-0.6324555320336759*cMr[5]; 
  data->AEM_S(5,11) = -0.6210590034081186*cMr[6]; 
  data->AEM_S(6,8) = -0.7071067811865475*cMr[6]; 
  data->AEM_S(6,9) = (-0.6210590034081186*cMr[7])-0.6324555320336759*cMr[5]; 
  data->AEM_S(6,10) = (-0.4517539514526256*cMr[6])-0.7071067811865475*cMr[4]; 
  data->AEM_S(6,11) = (-0.421637021355784*cMr[7])-0.6210590034081186*cMr[5]; 
  data->AEM_S(7,8) = -0.7071067811865475*cMr[7]; 
  data->AEM_S(7,9) = -0.6210590034081186*cMr[6]; 
  data->AEM_S(7,10) = (-0.421637021355784*cMr[7])-0.6210590034081186*cMr[5]; 
  data->AEM_S(7,11) = (-0.421637021355784*cMr[6])-0.7071067811865475*cMr[4]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(8,4) = 0.7071067811865475*m1r[4]; 
  data->AEM_S(8,5) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(8,6) = 0.7071067811865475*m1r[6]; 
  data->AEM_S(8,7) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(9,4) = 0.7071067811865475*m1r[5]; 
  data->AEM_S(9,5) = 0.6324555320336759*m1r[6]+0.7071067811865475*m1r[4]; 
  data->AEM_S(9,6) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  data->AEM_S(9,7) = 0.6210590034081186*m1r[6]; 
  data->AEM_S(10,4) = 0.7071067811865475*m1r[6]; 
  data->AEM_S(10,5) = 0.6210590034081186*m1r[7]+0.6324555320336759*m1r[5]; 
  data->AEM_S(10,6) = 0.4517539514526256*m1r[6]+0.7071067811865475*m1r[4]; 
  data->AEM_S(10,7) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  data->AEM_S(11,4) = 0.7071067811865475*m1r[7]; 
  data->AEM_S(11,5) = 0.6210590034081186*m1r[6]; 
  data->AEM_S(11,6) = 0.421637021355784*m1r[7]+0.6210590034081186*m1r[5]; 
  data->AEM_S(11,7) = 0.421637021355784*m1r[6]+0.7071067811865475*m1r[4]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(8,8) = 1.414213562373095*m0r[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(8,9) = 1.414213562373095*m0r[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(8,10) = 1.414213562373095*m0r[2]-0.7071067811865475*cEr[2]; 
  data->AEM_S(8,11) = 1.414213562373095*m0r[3]-0.7071067811865475*cEr[3]; 
  data->AEM_S(9,8) = 1.414213562373095*m0r[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(9,9) = 1.264911064067352*m0r[2]-0.6324555320336759*cEr[2]+1.414213562373095*m0r[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(9,10) = 1.242118006816237*m0r[3]-0.6210590034081186*cEr[3]+1.264911064067352*m0r[1]-0.6324555320336759*cEr[1]; 
  data->AEM_S(9,11) = 1.242118006816237*m0r[2]-0.6210590034081186*cEr[2]; 
  data->AEM_S(10,8) = 1.414213562373095*m0r[2]-0.7071067811865475*cEr[2]; 
  data->AEM_S(10,9) = 1.242118006816237*m0r[3]-0.6210590034081186*cEr[3]+1.264911064067352*m0r[1]-0.6324555320336759*cEr[1]; 
  data->AEM_S(10,10) = 0.9035079029052515*m0r[2]-0.4517539514526256*cEr[2]+1.414213562373095*m0r[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(10,11) = 0.8432740427115681*m0r[3]-0.421637021355784*cEr[3]+1.242118006816237*m0r[1]-0.6210590034081186*cEr[1]; 
  data->AEM_S(11,8) = 1.414213562373095*m0r[3]-0.7071067811865475*cEr[3]; 
  data->AEM_S(11,9) = 1.242118006816237*m0r[2]-0.6210590034081186*cEr[2]; 
  data->AEM_S(11,10) = 0.8432740427115681*m0r[3]-0.421637021355784*cEr[3]+1.242118006816237*m0r[1]-0.6210590034081186*cEr[1]; 
  data->AEM_S(11,11) = 0.8432740427115681*m0r[2]-0.421637021355784*cEr[2]+1.414213562373095*m0r[0]-0.7071067811865475*cEr[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<4,4>(0,4).setZero(); 
  data->AEM_S.block<4,4>(4,0).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m1r[6],m1r[7],m2r[0],m2r[1],m2r[2],m2r[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,8,1) = data->u_S.segment<8>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = data->u_S.segment<4>(8); 
 
} 
 
