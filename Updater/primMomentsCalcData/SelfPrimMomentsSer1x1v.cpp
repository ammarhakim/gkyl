#include <math.h> 
#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void SelfPrimMoments1x1vSer_P1(const double *m0, const double *m1, const double *m2, const double *fvmin, const double *fvmax, const double *vmin, const double *vmax, double *u, double *vtSq) 
{ 
  // m0,m1,m2:     moments of the distribution function. 
  // fvmax, fvmin: distribution function at the velocity boundaries. 
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
 
  // ....... M0-(v*f)|^(+vmax)_(-vmax) .......... // 
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
 
