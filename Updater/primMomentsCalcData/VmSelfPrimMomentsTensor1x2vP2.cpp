#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmSelfPrimMoments1x2vTensor_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
  // cM, cE:   vtSq*cM and vtSq*cE are corrections to u and vtSq, respectively. 
  // u:        velocity. 
  // vtSq:     squared thermal speed, sqrt(T/m). 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*(2.23606797749979*m0[2]-1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*m0[2]+1.732050807568877*m0[1]+m0[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*m2[2]-1.732050807568877*m2[1]+m2[0]) < 0) cellAvg = true; 
  if (0.7071067811865475*(2.23606797749979*m2[2]+1.732050807568877*m2[1]+m2[0]) < 0) cellAvg = true; 
 
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
 
