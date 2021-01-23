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
  if (-0.5*(2.449489742783178*m2S[1]-1.414213562373095*m2S[0]) < 0) cellAvg = true; 
  if (0.5*(2.449489742783178*m2S[1]+1.414213562373095*m2S[0]) < 0) cellAvg = true; 
 
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
 
