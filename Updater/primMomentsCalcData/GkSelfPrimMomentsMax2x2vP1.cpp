#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkSelfPrimMoments2x2vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
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
  if (-0.5*(1.732050807568877*(m2S[2]+m2S[1])-1.0*m2S[0]) < 0) cellAvg = true; 
  if (-0.5*(1.732050807568877*m2S[2]-1.732050807568877*m2S[1]-1.0*m2S[0]) < 0) cellAvg = true; 
  if (0.5*(1.732050807568877*m2S[2]-1.732050807568877*m2S[1]+m2S[0]) < 0) cellAvg = true; 
  if (0.5*(1.732050807568877*(m2S[2]+m2S[1])+m2S[0]) < 0) cellAvg = true; 
 
  double m0r[3]; 
  double m1r[3]; 
  double cMr[3]; 
  double cEr[3]; 
  double m0Sr[3]; 
  double m1Sr[3]; 
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
    cEr[0] = cE[0]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = 0.0; 
    m0Sr[2] = 0.0; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = 0.0; 
    m1Sr[2] = 0.0; 
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
    cMr[0] = cM[0]; 
    cMr[1] = cM[1]; 
    cMr[2] = cM[2]; 
    cEr[0] = cE[0]; 
    cEr[1] = cE[1]; 
    cEr[2] = cE[2]; 
    m0Sr[0] = m0S[0]; 
    m0Sr[1] = m0S[1]; 
    m0Sr[2] = m0S[2]; 
    m1Sr[0] = m1S[0]; 
    m1Sr[1] = m1S[1]; 
    m1Sr[2] = m1S[2]; 
    m2Sr[0] = m2S[0]; 
    m2Sr[1] = m2S[1]; 
    m2Sr[2] = m2S[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(6,6); 
 
  // ....... Block from weak multiply of uX and m0  .......... // 
  data->AEM_S(0,0) = 0.5*m0r[0]; 
  data->AEM_S(0,1) = 0.5*m0r[1]; 
  data->AEM_S(0,2) = 0.5*m0r[2]; 
  data->AEM_S(1,0) = 0.5*m0r[1]; 
  data->AEM_S(1,1) = 0.5*m0r[0]; 
  data->AEM_S(2,0) = 0.5*m0r[2]; 
  data->AEM_S(2,2) = 0.5*m0r[0]; 
 
  // ....... Block from correction to uX .......... // 
  data->AEM_S(0,3) = -0.5*cMr[0]; 
  data->AEM_S(0,4) = -0.5*cMr[1]; 
  data->AEM_S(0,5) = -0.5*cMr[2]; 
  data->AEM_S(1,3) = -0.5*cMr[1]; 
  data->AEM_S(1,4) = -0.5*cMr[0]; 
  data->AEM_S(2,3) = -0.5*cMr[2]; 
  data->AEM_S(2,5) = -0.5*cMr[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(3,0) = 0.5*m1Sr[0]; 
  data->AEM_S(3,1) = 0.5*m1Sr[1]; 
  data->AEM_S(3,2) = 0.5*m1Sr[2]; 
  data->AEM_S(4,0) = 0.5*m1Sr[1]; 
  data->AEM_S(4,1) = 0.5*m1Sr[0]; 
  data->AEM_S(5,0) = 0.5*m1Sr[2]; 
  data->AEM_S(5,2) = 0.5*m1Sr[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(3,3) = m0r[0]+0.5*m0Sr[0]-0.5*cEr[0]; 
  data->AEM_S(3,4) = m0r[1]+0.5*m0Sr[1]-0.5*cEr[1]; 
  data->AEM_S(3,5) = m0r[2]+0.5*m0Sr[2]-0.5*cEr[2]; 
  data->AEM_S(4,3) = m0r[1]+0.5*m0Sr[1]-0.5*cEr[1]; 
  data->AEM_S(4,4) = m0r[0]+0.5*m0Sr[0]-0.5*cEr[0]; 
  data->AEM_S(5,3) = m0r[2]+0.5*m0Sr[2]-0.5*cEr[2]; 
  data->AEM_S(5,5) = m0r[0]+0.5*m0Sr[0]-0.5*cEr[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m2Sr[0],m2Sr[1],m2Sr[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,3,1) = data->u_S.segment<3>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(3); 
 
} 
 
