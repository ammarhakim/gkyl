#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmSelfPrimMoments1x1vTensor_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
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
  double m1r[2]; 
  double cMr[2]; 
  double cEr[2]; 
  double m0Sr[2]; 
  double m1Sr[2]; 
  double m2Sr[2]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    cMr[0] = cM[0]; 
    cMr[1] = 0.0; 
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = 0.0; 
  } else { 
    m0r[0] = m0[0]; 
    m0r[1] = m0[1]; 
    m1r[0] = m1[0]; 
    m1r[1] = m1[1]; 
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(4,4); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0r[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,0) = 0.7071067811865475*m0r[1]; 
  data->AEM_S(1,1) = 0.7071067811865475*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,2) = -0.7071067811865475*cMr[0]; 
  data->AEM_S(0,3) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,2) = -0.7071067811865475*cMr[1]; 
  data->AEM_S(1,3) = -0.7071067811865475*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(2,0) = 0.7071067811865475*m1Sr[0]; 
  data->AEM_S(2,1) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(3,0) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(3,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cEr[0]; 
  data->AEM_S(2,3) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m0Sr[1]-0.7071067811865475*cEr[1]; 
  data->AEM_S(3,3) = 0.7071067811865475*m0Sr[0]-0.7071067811865475*cEr[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m2Sr[0],m2Sr[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,2,1) = data->u_S.segment<2>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = data->u_S.segment<2>(2); 
 
} 
 
