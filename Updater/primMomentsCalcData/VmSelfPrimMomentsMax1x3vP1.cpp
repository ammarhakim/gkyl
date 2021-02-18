#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmSelfPrimMoments1x3vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
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
  if (-0.5*(2.449489742783178*m2S[1]-1.414213562373095*m2S[0]) < 0) cellAvg = true; 
  if (0.5*(2.449489742783178*m2S[1]+1.414213562373095*m2S[0]) < 0) cellAvg = true; 
 
  double m0r[2]; 
  double m1r[6]; 
  double cMr[6]; 
  double cEr[2]; 
  double m0Sr[2]; 
  double m1Sr[6]; 
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
    m1r[4] = m1[4]; 
    m1r[5] = 0.0; 
    cMr[4] = cM[4]; 
    cMr[5] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = 0.0; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
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
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m1Sr[3] = m1S[3]; 
    m1Sr[4] = m1S[4]; 
    m1Sr[5] = m1S[5]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(8,8); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,6) = -0.7071067811865475*cMr[0]; 
  data->AEM_S(0,7) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,6) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,7) = -0.7071067811865475*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(6,0) = 0.7071067811865475*m1Sr[0]; 
  data->AEM_S(6,1) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(7,0) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(7,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from weak multiply of uY and m0  .......... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(2,3) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(3,3) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uY .......... // 
  data->AEM_S(2,6) = -0.7071067811865475*cMr[2]; 
  data->AEM_S(2,7) = -0.7071067811865475*cMr[3]; 
  data->AEM_S(3,6) = -0.7071067811865475*cMr[3]; 
  data->AEM_S(3,7) = -0.7071067811865475*cMr[2]; 
 
  // ....... Block from weak multiply of uY and m1Y  .......... // 
  data->AEM_S(6,2) = 0.7071067811865475*m1Sr[2]; 
  data->AEM_S(6,3) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(7,2) = 0.7071067811865475*m1Sr[3]; 
  data->AEM_S(7,3) = 0.7071067811865475*m1Sr[2]; 
 
  // ....... Block from weak multiply of uZ and m0  .......... // 
  data->AEM_S(4,4) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(4,5) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(5,4) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(5,5) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uZ .......... // 
  data->AEM_S(4,6) = -0.7071067811865475*cMr[4]; 
  data->AEM_S(4,7) = -0.7071067811865475*cMr[5]; 
  data->AEM_S(5,6) = -0.7071067811865475*cMr[5]; 
  data->AEM_S(5,7) = -0.7071067811865475*cMr[4]; 
 
  // ....... Block from weak multiply of uZ and m1Z  .......... // 
  data->AEM_S(6,4) = 0.7071067811865475*m1Sr[4]; 
  data->AEM_S(6,5) = 0.7071067811865475*m1Sr[5]; 
  data->AEM_S(7,4) = 0.7071067811865475*m1Sr[5]; 
  data->AEM_S(7,5) = 0.7071067811865475*m1Sr[4]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(6,6) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(6,7) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(7,6) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(7,7) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cEr[0]; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<2,4>(0,2).setZero(); 
  data->AEM_S.block<4,2>(2,0).setZero(); 
  data->AEM_S.block<2,2>(2,4).setZero(); 
  data->AEM_S.block<2,2>(4,2).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m1r[4],m1r[5],m2Sr[0],m2Sr[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,6,1) = data->u_S.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = data->u_S.segment<2>(6); 
 
} 
 
