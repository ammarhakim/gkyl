#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkSelfPrimMoments1x2vMax_P1(binOpData_t *data, const double *m0, const double *m1, const double *m0S, const double *m1S, const double *m2S, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1:       moments of the distribution function. 
  // m0S,m1S,m1S: star moments (only used for piecewise linear). 
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
  double m1r[2]; 
  double m0Sr[2]; 
  double m1Sr[2]; 
  double m2Sr[2]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
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
  data->AEM_S(0,2) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,3) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,2) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,3) = -0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(2,0) = 0.7071067811865475*m1Sr[0]; 
  data->AEM_S(2,1) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(3,0) = 0.7071067811865475*m1Sr[1]; 
  data->AEM_S(3,1) = 0.7071067811865475*m1Sr[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(2,2) = 1.414213562373095*m0r[0]+0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(2,3) = 1.414213562373095*m0r[1]+0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(3,2) = 1.414213562373095*m0r[1]+0.7071067811865475*m0Sr[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(3,3) = 1.414213562373095*m0r[0]+0.7071067811865475*m0Sr[0]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m2Sr[0],m2Sr[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,2,1) = data->u_S.segment<2>(0); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = data->u_S.segment<2>(2); 
 
} 
 
void GkSelfPrimMoments1x2vMax_P2(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
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
  double m1r[3]; 
  double m2r[3]; 
  if (cellAvg) { 
    m0r[0] = m0[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = m1[0]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(6,6); 
 
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
  data->AEM_S(0,3) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,4) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,5) = -0.7071067811865475*cM[2]; 
  data->AEM_S(1,3) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,4) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,5) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,3) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,4) = -0.6324555320336759*cM[1]; 
  data->AEM_S(2,5) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(3,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(3,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(3,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(4,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(4,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(4,2) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(5,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(5,1) = 0.6324555320336759*m1r[1]; 
  data->AEM_S(5,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(3,3) = 2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(3,4) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(3,5) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(4,3) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(4,4) = 1.897366596101028*m0r[2]-0.6324555320336759*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(4,5) = 1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(5,3) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(5,4) = 1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(5,5) = 1.355261854357877*m0r[2]-0.4517539514526256*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m2r[0],m2r[1],m2r[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,3,1) = data->u_S.segment<3>(0); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = data->u_S.segment<3>(3); 
 
} 
 
void GkSelfPrimMoments1x2vMax_P3(binOpData_t *data, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE, double *u, double *vtSq) 
{ 
  // m0,m1,m2: moments of the distribution function. 
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
  double m1r[4]; 
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
    m2r[0] = m2[0]; 
    m2r[1] = m2[1]; 
    m2r[2] = m2[2]; 
    m2r[3] = m2[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(8,8); 
 
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
  data->AEM_S(0,4) = -0.7071067811865475*cM[0]; 
  data->AEM_S(0,5) = -0.7071067811865475*cM[1]; 
  data->AEM_S(0,6) = -0.7071067811865475*cM[2]; 
  data->AEM_S(0,7) = -0.7071067811865475*cM[3]; 
  data->AEM_S(1,4) = -0.7071067811865475*cM[1]; 
  data->AEM_S(1,5) = (-0.6324555320336759*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(1,6) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(1,7) = -0.6210590034081186*cM[2]; 
  data->AEM_S(2,4) = -0.7071067811865475*cM[2]; 
  data->AEM_S(2,5) = (-0.6210590034081186*cM[3])-0.6324555320336759*cM[1]; 
  data->AEM_S(2,6) = (-0.4517539514526256*cM[2])-0.7071067811865475*cM[0]; 
  data->AEM_S(2,7) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,4) = -0.7071067811865475*cM[3]; 
  data->AEM_S(3,5) = -0.6210590034081186*cM[2]; 
  data->AEM_S(3,6) = (-0.421637021355784*cM[3])-0.6210590034081186*cM[1]; 
  data->AEM_S(3,7) = (-0.421637021355784*cM[2])-0.7071067811865475*cM[0]; 
 
  // ....... Block from weak multiply of uX and m1X  .......... // 
  data->AEM_S(4,0) = 0.7071067811865475*m1r[0]; 
  data->AEM_S(4,1) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(4,2) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(4,3) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(5,0) = 0.7071067811865475*m1r[1]; 
  data->AEM_S(5,1) = 0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(5,2) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(5,3) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(6,0) = 0.7071067811865475*m1r[2]; 
  data->AEM_S(6,1) = 0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]; 
  data->AEM_S(6,2) = 0.4517539514526256*m1r[2]+0.7071067811865475*m1r[0]; 
  data->AEM_S(6,3) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(7,0) = 0.7071067811865475*m1r[3]; 
  data->AEM_S(7,1) = 0.6210590034081186*m1r[2]; 
  data->AEM_S(7,2) = 0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]; 
  data->AEM_S(7,3) = 0.421637021355784*m1r[2]+0.7071067811865475*m1r[0]; 
 
  // ....... Block from correction to vtSq .......... // 
  data->AEM_S(4,4) = 2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(4,5) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(4,6) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(4,7) = 2.121320343559642*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(5,4) = 2.121320343559642*m0r[1]-0.7071067811865475*cE[1]; 
  data->AEM_S(5,5) = 1.897366596101028*m0r[2]-0.6324555320336759*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(5,6) = 1.863177010224355*m0r[3]-0.6210590034081186*cE[3]+1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(5,7) = 1.863177010224355*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(6,4) = 2.121320343559642*m0r[2]-0.7071067811865475*cE[2]; 
  data->AEM_S(6,5) = 1.863177010224355*m0r[3]-0.6210590034081186*cE[3]+1.897366596101028*m0r[1]-0.6324555320336759*cE[1]; 
  data->AEM_S(6,6) = 1.355261854357877*m0r[2]-0.4517539514526256*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
  data->AEM_S(6,7) = 1.264911064067352*m0r[3]-0.421637021355784*cE[3]+1.863177010224355*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(7,4) = 2.121320343559642*m0r[3]-0.7071067811865475*cE[3]; 
  data->AEM_S(7,5) = 1.863177010224355*m0r[2]-0.6210590034081186*cE[2]; 
  data->AEM_S(7,6) = 1.264911064067352*m0r[3]-0.421637021355784*cE[3]+1.863177010224355*m0r[1]-0.6210590034081186*cE[1]; 
  data->AEM_S(7,7) = 1.264911064067352*m0r[2]-0.421637021355784*cE[2]+2.121320343559642*m0r[0]-0.7071067811865475*cE[0]; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << m1r[0],m1r[1],m1r[2],m1r[3],m2r[0],m2r[1],m2r[2],m2r[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(u,4,1) = data->u_S.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = data->u_S.segment<4>(4); 
 
} 
 
