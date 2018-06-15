#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments1x1vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq) 
{ 
  // m0,m1,m2:     moments of the distribution function. 
  // fvmax, fvmin: distribution function at the velocity boundaries. 
  // vmax, vmin:   maximum and minimum velocity of the velocity grid. 
  // u:            velocity. 
  // vtSq:         squared thermal speed, sqrt(T/m). 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(2,2); 
  Eigen::VectorXd bEV(2); 
  Eigen::VectorXd xEV(2); 
 
  // ....... Compute u through weak division m1/m0 .......... // 
  AEM(0,0) = 0.7071067811865475*m0[0]; 
  AEM(0,1) = 0.7071067811865475*m0[1]; 
  AEM(1,0) = 0.7071067811865475*m0[1]; 
  AEM(1,1) = 0.7071067811865475*m0[0]; 
 
  for(unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int b0 = 2*vd; 
    bEV << m1[b0],m1[b0+1]; 
 
    xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
    Eigen::Map<VectorXd>(u+vd*2,2,1) = xEV; 
  } 
 
  // ....... Get kinetic energy density via weak dot product u.m1 .......... // 
  double kinEnergyDens[2]; 
  for (unsigned short int k=0; k<2; k++) 
  { 
    kinEnergyDens[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    kinEnergyDens[0] += 0.7071067811865475*m1[a0+1]*u[a0+1]+0.7071067811865475*m1[a0]*u[a0]; 
    kinEnergyDens[1] += 0.7071067811865475*m1[a0]*u[a0+1]+0.7071067811865475*u[a0]*m1[a0+1]; 
  } 
 
  // ....... Thermal energy density: M2-u.M1 .......... // 
  double thEnergyDens[2]; 
  for (unsigned short int i=0; i<2; i++) 
  { 
    thEnergyDens[i] = m2[i] - kinEnergyDens[i]; 
  } 
 
  // ....... M0-sum_i int dS_i (v_i*f)|^(vmax_i)_(vmin_i) .......... // 
  double m0c[2]; 
  m0c[0] = vmin[0]*(1.224744871391589*fvmin[2]-0.7071067811865475*fvmin[0])+vmax[0]*((-1.224744871391589*fvmax[2])-0.7071067811865475*fvmax[0])+m0[0]; 
  m0c[1] = vmin[0]*(1.224744871391589*fvmin[3]-0.7071067811865475*fvmin[1])+vmax[0]*((-1.224744871391589*fvmax[3])-0.7071067811865475*fvmax[1])+m0[1]; 
 
  // ....... Compute vtSq through weak division thEnergyDens/m0c .......... // 
  AEM(0,0) = 0.7071067811865475*m0c[0]; 
  AEM(0,1) = 0.7071067811865475*m0c[1]; 
  AEM(1,0) = 0.7071067811865475*m0c[1]; 
  AEM(1,1) = 0.7071067811865475*m0c[0]; 
 
  bEV << thEnergyDens[0],thEnergyDens[1]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(vtSq,2,1) = xEV; 
 
} 
 
void SelfPrimMoments1x1vSer_P2(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq) 
{ 
  // m0,m1,m2:     moments of the distribution function. 
  // fvmax, fvmin: distribution function at the velocity boundaries. 
  // vmax, vmin:   maximum and minimum velocity of the velocity grid. 
  // u:            velocity. 
  // vtSq:         squared thermal speed, sqrt(T/m). 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(3,3); 
  Eigen::VectorXd bEV(3); 
  Eigen::VectorXd xEV(3); 
 
  // ....... Compute u through weak division m1/m0 .......... // 
  AEM(0,0) = 0.7071067811865475*m0[0]; 
  AEM(0,1) = 0.7071067811865475*m0[1]; 
  AEM(0,2) = 0.7071067811865475*m0[2]; 
  AEM(1,0) = 0.7071067811865475*m0[1]; 
  AEM(1,1) = 0.6324555320336759*m0[2]+0.7071067811865475*m0[0]; 
  AEM(1,2) = 0.6324555320336759*m0[1]; 
  AEM(2,0) = 0.7071067811865475*m0[2]; 
  AEM(2,1) = 0.6324555320336759*m0[1]; 
  AEM(2,2) = 0.4517539514526256*m0[2]+0.7071067811865475*m0[0]; 
 
  for(unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int b0 = 3*vd; 
    bEV << m1[b0],m1[b0+1],m1[b0+2]; 
 
    xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
    Eigen::Map<VectorXd>(u+vd*3,3,1) = xEV; 
  } 
 
  // ....... Get kinetic energy density via weak dot product u.m1 .......... // 
  double kinEnergyDens[3]; 
  for (unsigned short int k=0; k<3; k++) 
  { 
    kinEnergyDens[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    kinEnergyDens[0] += 0.7071067811865475*m1[a0+2]*u[a0+2]+0.7071067811865475*m1[a0+1]*u[a0+1]+0.7071067811865475*m1[a0]*u[a0]; 
    kinEnergyDens[1] += 0.6324555320336759*m1[a0+1]*u[a0+2]+0.6324555320336759*u[a0+1]*m1[a0+2]+0.7071067811865475*m1[a0]*u[a0+1]+0.7071067811865475*u[a0]*m1[a0+1]; 
    kinEnergyDens[2] += 0.4517539514526256*m1[a0+2]*u[a0+2]+0.7071067811865475*m1[a0]*u[a0+2]+0.7071067811865475*u[a0]*m1[a0+2]+0.6324555320336759*m1[a0+1]*u[a0+1]; 
  } 
 
  // ....... Thermal energy density: M2-u.M1 .......... // 
  double thEnergyDens[3]; 
  for (unsigned short int i=0; i<3; i++) 
  { 
    thEnergyDens[i] = m2[i] - kinEnergyDens[i]; 
  } 
 
  // ....... M0-sum_i int dS_i (v_i*f)|^(vmax_i)_(vmin_i) .......... // 
  double m0c[3]; 
  m0c[0] = vmin[0]*((-1.58113883008419*fvmin[5])+1.224744871391589*fvmin[2]-0.7071067811865475*fvmin[0])+vmax[0]*((-1.58113883008419*fvmax[5])-1.224744871391589*fvmax[2]-0.7071067811865475*fvmax[0])+m0[0]; 
  m0c[1] = vmin[0]*((-1.58113883008419*fvmin[7])+1.224744871391589*fvmin[3]-0.7071067811865475*fvmin[1])+vmax[0]*((-1.58113883008419*fvmax[7])-1.224744871391589*fvmax[3]-0.7071067811865475*fvmax[1])+m0[1]; 
  m0c[2] = vmin[0]*(1.224744871391589*fvmin[6]-0.7071067811865475*fvmin[4])+vmax[0]*((-1.224744871391589*fvmax[6])-0.7071067811865475*fvmax[4])+m0[2]; 
 
  // ....... Compute vtSq through weak division thEnergyDens/m0c .......... // 
  AEM(0,0) = 0.7071067811865475*m0c[0]; 
  AEM(0,1) = 0.7071067811865475*m0c[1]; 
  AEM(0,2) = 0.7071067811865475*m0c[2]; 
  AEM(1,0) = 0.7071067811865475*m0c[1]; 
  AEM(1,1) = 0.6324555320336759*m0c[2]+0.7071067811865475*m0c[0]; 
  AEM(1,2) = 0.6324555320336759*m0c[1]; 
  AEM(2,0) = 0.7071067811865475*m0c[2]; 
  AEM(2,1) = 0.6324555320336759*m0c[1]; 
  AEM(2,2) = 0.4517539514526256*m0c[2]+0.7071067811865475*m0c[0]; 
 
  bEV << thEnergyDens[0],thEnergyDens[1],thEnergyDens[2]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(vtSq,3,1) = xEV; 
 
} 
 
void SelfPrimMoments1x1vSer_P3(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq) 
{ 
  // m0,m1,m2:     moments of the distribution function. 
  // fvmax, fvmin: distribution function at the velocity boundaries. 
  // vmax, vmin:   maximum and minimum velocity of the velocity grid. 
  // u:            velocity. 
  // vtSq:         squared thermal speed, sqrt(T/m). 
 
  // Declare Eigen matrix and vectors for weak division. 
  Eigen::MatrixXd AEM(4,4); 
  Eigen::VectorXd bEV(4); 
  Eigen::VectorXd xEV(4); 
 
  // ....... Compute u through weak division m1/m0 .......... // 
  AEM(0,0) = 0.7071067811865475*m0[0]; 
  AEM(0,1) = 0.7071067811865475*m0[1]; 
  AEM(0,2) = 0.7071067811865475*m0[2]; 
  AEM(0,3) = 0.7071067811865475*m0[3]; 
  AEM(1,0) = 0.7071067811865475*m0[1]; 
  AEM(1,1) = 0.6324555320336759*m0[2]+0.7071067811865475*m0[0]; 
  AEM(1,2) = 0.6210590034081186*m0[3]+0.6324555320336759*m0[1]; 
  AEM(1,3) = 0.6210590034081186*m0[2]; 
  AEM(2,0) = 0.7071067811865475*m0[2]; 
  AEM(2,1) = 0.6210590034081186*m0[3]+0.6324555320336759*m0[1]; 
  AEM(2,2) = 0.4517539514526256*m0[2]+0.7071067811865475*m0[0]; 
  AEM(2,3) = 0.421637021355784*m0[3]+0.6210590034081186*m0[1]; 
  AEM(3,0) = 0.7071067811865475*m0[3]; 
  AEM(3,1) = 0.6210590034081186*m0[2]; 
  AEM(3,2) = 0.421637021355784*m0[3]+0.6210590034081186*m0[1]; 
  AEM(3,3) = 0.421637021355784*m0[2]+0.7071067811865475*m0[0]; 
 
  for(unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int b0 = 4*vd; 
    bEV << m1[b0],m1[b0+1],m1[b0+2],m1[b0+3]; 
 
    xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
    Eigen::Map<VectorXd>(u+vd*4,4,1) = xEV; 
  } 
 
  // ....... Get kinetic energy density via weak dot product u.m1 .......... // 
  double kinEnergyDens[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    kinEnergyDens[k] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    kinEnergyDens[0] += 0.7071067811865475*m1[a0+3]*u[a0+3]+0.7071067811865475*m1[a0+2]*u[a0+2]+0.7071067811865475*m1[a0+1]*u[a0+1]+0.7071067811865475*m1[a0]*u[a0]; 
    kinEnergyDens[1] += 0.6210590034081186*m1[a0+2]*u[a0+3]+0.6210590034081186*u[a0+2]*m1[a0+3]+0.6324555320336759*m1[a0+1]*u[a0+2]+0.6324555320336759*u[a0+1]*m1[a0+2]+0.7071067811865475*m1[a0]*u[a0+1]+0.7071067811865475*u[a0]*m1[a0+1]; 
    kinEnergyDens[2] += 0.421637021355784*m1[a0+3]*u[a0+3]+0.6210590034081186*m1[a0+1]*u[a0+3]+0.6210590034081186*u[a0+1]*m1[a0+3]+0.4517539514526256*m1[a0+2]*u[a0+2]+0.7071067811865475*m1[a0]*u[a0+2]+0.7071067811865475*u[a0]*m1[a0+2]+0.6324555320336759*m1[a0+1]*u[a0+1]; 
    kinEnergyDens[3] += 0.421637021355784*m1[a0+2]*u[a0+3]+0.7071067811865475*m1[a0]*u[a0+3]+0.421637021355784*u[a0+2]*m1[a0+3]+0.7071067811865475*u[a0]*m1[a0+3]+0.6210590034081186*m1[a0+1]*u[a0+2]+0.6210590034081186*u[a0+1]*m1[a0+2]; 
  } 
 
  // ....... Thermal energy density: M2-u.M1 .......... // 
  double thEnergyDens[4]; 
  for (unsigned short int i=0; i<4; i++) 
  { 
    thEnergyDens[i] = m2[i] - kinEnergyDens[i]; 
  } 
 
  // ....... M0-sum_i int dS_i (v_i*f)|^(vmax_i)_(vmin_i) .......... // 
  double m0c[4]; 
  m0c[0] = vmin[0]*(1.870828693386971*fvmin[9]-1.58113883008419*fvmin[5]+1.224744871391589*fvmin[2]-0.7071067811865475*fvmin[0])+vmax[0]*((-1.870828693386971*fvmax[9])-1.58113883008419*fvmax[5]-1.224744871391589*fvmax[2]-0.7071067811865475*fvmax[0])+m0[0]; 
  m0c[1] = vmin[0]*(1.870828693386971*fvmin[11]-1.58113883008419*fvmin[7]+1.224744871391589*fvmin[3]-0.7071067811865475*fvmin[1])+vmax[0]*((-1.870828693386971*fvmax[11])-1.58113883008419*fvmax[7]-1.224744871391589*fvmax[3]-0.7071067811865475*fvmax[1])+m0[1]; 
  m0c[2] = vmin[0]*(1.224744871391589*fvmin[6]-0.7071067811865475*fvmin[4])+vmax[0]*((-1.224744871391589*fvmax[6])-0.7071067811865475*fvmax[4])+m0[2]; 
  m0c[3] = vmin[0]*(1.224744871391589*fvmin[10]-0.7071067811865475*fvmin[8])+vmax[0]*((-1.224744871391589*fvmax[10])-0.7071067811865475*fvmax[8])+m0[3]; 
 
  // ....... Compute vtSq through weak division thEnergyDens/m0c .......... // 
  AEM(0,0) = 0.7071067811865475*m0c[0]; 
  AEM(0,1) = 0.7071067811865475*m0c[1]; 
  AEM(0,2) = 0.7071067811865475*m0c[2]; 
  AEM(0,3) = 0.7071067811865475*m0c[3]; 
  AEM(1,0) = 0.7071067811865475*m0c[1]; 
  AEM(1,1) = 0.6324555320336759*m0c[2]+0.7071067811865475*m0c[0]; 
  AEM(1,2) = 0.6210590034081186*m0c[3]+0.6324555320336759*m0c[1]; 
  AEM(1,3) = 0.6210590034081186*m0c[2]; 
  AEM(2,0) = 0.7071067811865475*m0c[2]; 
  AEM(2,1) = 0.6210590034081186*m0c[3]+0.6324555320336759*m0c[1]; 
  AEM(2,2) = 0.4517539514526256*m0c[2]+0.7071067811865475*m0c[0]; 
  AEM(2,3) = 0.421637021355784*m0c[3]+0.6210590034081186*m0c[1]; 
  AEM(3,0) = 0.7071067811865475*m0c[3]; 
  AEM(3,1) = 0.6210590034081186*m0c[2]; 
  AEM(3,2) = 0.421637021355784*m0c[3]+0.6210590034081186*m0c[1]; 
  AEM(3,3) = 0.421637021355784*m0c[2]+0.7071067811865475*m0c[0]; 
 
  bEV << thEnergyDens[0],thEnergyDens[1],thEnergyDens[2],thEnergyDens[3]; 
 
  xEV = AEM.colPivHouseholderQr().solve(bEV); 
 
  Eigen::Map<VectorXd>(vtSq,4,1) = xEV; 
 
} 
 
