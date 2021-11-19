#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkLBOGCrossPrimMoments1x1vSer_P1(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // m0S,m1S,m1S:        star moments (only used for piecewise linear). 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg;
  cellAvg = false;
  if (0.7071067811865475*m0Self[0]-1.224744871391589*m0Self[1] < 0) cellAvg = true; 
  if (0.7071067811865475*m2SSelf[0]-1.224744871391589*m2SSelf[1] < 0) cellAvg = true; 
  if (1.224744871391589*m0Self[1]+0.7071067811865475*m0Self[0] < 0) cellAvg = true; 
  if (1.224744871391589*m2SSelf[1]+0.7071067811865475*m2SSelf[0] < 0) cellAvg = true; 
 
  double m0rSelf[2]; 
  double m1rSelf[2]; 
  double m2rSelf[2]; 
  double cMrSelf[2]; 
  double cErSelf[2]; 
  double m0SrSelf[2]; 
  double m1SrSelf[2]; 
  double m2SrSelf[2]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = 0.0; 
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = 0.0; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = 0.0; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = m1Self[1]; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = m2Self[1]; 
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = cMSelf[1]; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = cESelf[1]; 
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = m0SSelf[1]; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = m1SSelf[1]; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = m2SSelf[1]; 
  } 
 
  cellAvg = false;
  if (0.7071067811865475*m0Other[0]-1.224744871391589*m0Other[1] < 0) cellAvg = true; 
  if (0.7071067811865475*m2SOther[0]-1.224744871391589*m2SOther[1] < 0) cellAvg = true; 
  if (1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) cellAvg = true; 
  if (1.224744871391589*m2SOther[1]+0.7071067811865475*m2SOther[0] < 0) cellAvg = true; 
 
  double m0rOther[2]; 
  double m1rOther[2]; 
  double m2rOther[2]; 
  double cMrOther[2]; 
  double cErOther[2]; 
  double m0SrOther[2]; 
  double m1SrOther[2]; 
  double m2SrOther[2]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = 0.0; 
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = 0.0; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = 0.0; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = m1Other[1]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = cMOther[1]; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = cEOther[1]; 
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = m0SOther[1]; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = m1SOther[1]; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak system. 
  data->AEM_S.setZero(); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double mnuM1sum[2]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    mnuM1sum[vd] = 0.0; 
  } 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,2) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
  data->AEM_S(0,3) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(1,2) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(1,3) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,4) = 0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(0,5) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,4) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,5) = 0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,6) = -0.7071067811865475*cMrOther[0]*mnuOther; 
  data->AEM_S(0,7) = -0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(1,6) = -0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(1,7) = -0.7071067811865475*cMrOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(2,0) = 0.7071067811865475*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(2,1) = 0.7071067811865475*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(3,0) = 0.7071067811865475*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(3,1) = 0.7071067811865475*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(2,4) = 0.7071067811865475*m1SrOther[0]*mnuOther; 
  data->AEM_S(2,5) = 0.7071067811865475*m1SrOther[1]*mnuOther; 
  data->AEM_S(3,4) = 0.7071067811865475*m1SrOther[1]*mnuOther; 
  data->AEM_S(3,5) = 0.7071067811865475*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(3,2) = 0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(3,3) = 0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(2,6) = 0.7071067811865475*m0SrOther[0]*mnuOther-0.7071067811865475*cErOther[0]*mnuOther; 
  data->AEM_S(2,7) = 0.7071067811865475*m0SrOther[1]*mnuOther-0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(3,6) = 0.7071067811865475*m0SrOther[1]*mnuOther-0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(3,7) = 0.7071067811865475*m0SrOther[0]*mnuOther-0.7071067811865475*cErOther[0]*mnuOther; 
 
  double mnuM2sum[2]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2SrSelf[0]*mnuSelf+m2SrOther[0]*mnuOther; 
  mnuM2sum[1] = m2SrSelf[1]*mnuSelf+m2SrOther[1]*mnuOther; 
 
  double m1Relax[2]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  double m1EffD[2]; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(4,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,1) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(4,2) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
  data->AEM_S(4,3) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(5,2) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(5,3) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(4,4) = -0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(4,5) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(5,4) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(5,5) = -0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(4,6) = 0.7071067811865475*cMrOther[0]*mnuOther; 
  data->AEM_S(4,7) = 0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(5,6) = 0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(5,7) = 0.7071067811865475*cMrOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(6,0) = (-0.5*m0rSelf[1]*uSelf[1]*mnuSelf)-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(6,1) = (-0.5*m0rSelf[0]*uSelf[1]*mnuSelf)+0.7071067811865475*m1SrSelf[1]*mnuSelf-0.5*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,0) = (-0.5*m0rSelf[0]*uSelf[1]*mnuSelf)+0.7071067811865475*m1SrSelf[1]*mnuSelf-0.5*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,1) = (-0.9*m0rSelf[1]*uSelf[1]*mnuSelf)-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(6,4) = 0.5*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1SrOther[0]*mnuOther; 
  data->AEM_S(6,5) = 0.5*m0rOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1SrOther[1]*mnuOther+0.5*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(7,4) = 0.5*m0rOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1SrOther[1]*mnuOther+0.5*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(7,5) = 0.9*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1SrOther[0]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfX-m0Self*m1OtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[0] = 0.7071067811865475*m0rOther[1]*m1rSelf[1]-0.7071067811865475*m0rSelf[1]*m1rOther[1]+0.7071067811865475*m0rOther[0]*m1rSelf[0]-0.7071067811865475*m0rSelf[0]*m1rOther[0]; 
  m1EffD[1] = 0.7071067811865475*m0rOther[0]*m1rSelf[1]-0.7071067811865475*m0rSelf[0]*m1rOther[1]-0.7071067811865475*m1rOther[0]*m0rSelf[1]+0.7071067811865475*m1rSelf[0]*m0rOther[1]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(2,2); 
  dataDiv->AEM_S(0,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf+0.7071067811865475*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf+0.7071067811865475*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf+0.7071067811865475*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.7071067811865475*m0rSelf[0]*mnuSelf+0.7071067811865475*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[0],m1EffD[1]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+0,2,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (-2.0*m1EffD[0]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (-2.0*m1EffD[1]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
 
  double ucMSelf[2]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMrSelf[a0+1]*uSelf[a0+1]+0.7071067811865475*cMrSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.7071067811865475*cMrSelf[a0]*uSelf[a0+1]+0.7071067811865475*uSelf[a0]*cMrSelf[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(6,2) = 0.7071067811865475*ucMSelf[0]*mnuSelf+0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
  data->AEM_S(6,3) = 0.7071067811865475*ucMSelf[1]*mnuSelf+0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(7,2) = 0.7071067811865475*ucMSelf[1]*mnuSelf+0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(7,3) = 0.7071067811865475*ucMSelf[0]*mnuSelf+0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
 
  double ucMOther[2]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.7071067811865475*cMrOther[a0+1]*uOther[a0+1]+0.7071067811865475*cMrOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.7071067811865475*cMrOther[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*cMrOther[a0+1]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(6,6) = (-0.7071067811865475*ucMOther[0]*mnuOther)-0.7071067811865475*m0SrOther[0]*mnuOther+0.7071067811865475*cErOther[0]*mnuOther; 
  data->AEM_S(6,7) = (-0.7071067811865475*ucMOther[1]*mnuOther)-0.7071067811865475*m0SrOther[1]*mnuOther+0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(7,6) = (-0.7071067811865475*ucMOther[1]*mnuOther)-0.7071067811865475*m0SrOther[1]*mnuOther+0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(7,7) = (-0.7071067811865475*ucMOther[0]*mnuOther)-0.7071067811865475*m0SrOther[0]*mnuOther+0.7071067811865475*cErOther[0]*mnuOther; 
 
  double kinESelf[2]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinESelf[0] += 0.7071067811865475*m1rSelf[a0+1]*uSelf[a0+1]+0.7071067811865475*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.7071067811865475*m1rSelf[a0]*uSelf[a0+1]+0.7071067811865475*uSelf[a0]*m1rSelf[a0+1]; 
  } 
 
  double kinEOther[2]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    kinEOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.7071067811865475*m1rOther[a0+1]*uOther[a0+1]+0.7071067811865475*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.7071067811865475*m1rOther[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*m1rOther[a0+1]; 
  } 
 
  double relKinE[2]; 
  // zero out array with dot product of uSelf-uOther and m1EffD. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.7071067811865475*m1EffD[a0+1]*uSelf[a0+1]-0.7071067811865475*m1EffD[a0+1]*uOther[a0+1]+0.7071067811865475*m1EffD[a0]*uSelf[a0]-0.7071067811865475*m1EffD[a0]*uOther[a0]; 
    relKinE[1] += 0.7071067811865475*m1EffD[a0]*uSelf[a0+1]-0.7071067811865475*m1EffD[a0]*uOther[a0+1]+0.7071067811865475*uSelf[a0]*m1EffD[a0+1]-0.7071067811865475*uOther[a0]*m1EffD[a0+1]; 
  } 
 
  // Divide m0Other*(m2Self-kinESelf) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Other and m2Self-uSelf.m1Self. 
  double m0OtherThESelf[2]; 
  m0OtherThESelf[0] = 0.7071067811865475*m0rOther[1]*m2rSelf[1]-0.7071067811865475*kinESelf[1]*m0rOther[1]+0.7071067811865475*m0rOther[0]*m2rSelf[0]-0.7071067811865475*kinESelf[0]*m0rOther[0]; 
  m0OtherThESelf[1] = 0.7071067811865475*m0rOther[0]*m2rSelf[1]+0.7071067811865475*m2rSelf[0]*m0rOther[1]-0.7071067811865475*kinESelf[0]*m0rOther[1]-0.7071067811865475*m0rOther[0]*kinESelf[1]; 
  dataDiv->BEV_S << m0OtherThESelf[0],m0OtherThESelf[1]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthSelf[2]; 
  Eigen::Map<VectorXd>(effEthSelf,2,1) = dataDiv->u_S; 
 
  // Divide m0Self*(m2Other-kinEOther) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Self and m2Other-uOther.m1Other. 
  double m0SelfThEOther[2]; 
  m0SelfThEOther[0] = 0.7071067811865475*m0rSelf[1]*m2rOther[1]-0.7071067811865475*kinEOther[1]*m0rSelf[1]+0.7071067811865475*m0rSelf[0]*m2rOther[0]-0.7071067811865475*kinEOther[0]*m0rSelf[0]; 
  m0SelfThEOther[1] = 0.7071067811865475*m0rSelf[0]*m2rOther[1]+0.7071067811865475*m2rOther[0]*m0rSelf[1]-0.7071067811865475*kinEOther[0]*m0rSelf[1]-0.7071067811865475*m0rSelf[0]*kinEOther[1]; 
  dataDiv->BEV_S << m0SelfThEOther[0],m0SelfThEOther[1]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthOther[2]; 
  Eigen::Map<VectorXd>(effEthOther,2,1) = dataDiv->u_S; 
 
  double m2Relax[2]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(1.0*relKinE[0]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[0]*mSelf)/(mSelf+mOther)+(1.0*relKinE[0]*mOther)/(mSelf+mOther)+(2.0*effEthOther[0]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2SrSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2SrOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(1.0*relKinE[1]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[1]*mSelf)/(mSelf+mOther)+(1.0*relKinE[1]*mOther)/(mSelf+mOther)+(2.0*effEthOther[1]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2SrSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2SrOther[1])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM2sum[0],mnuM2sum[1],m1Relax[0],m1Relax[1],m2Relax[0],m2Relax[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,2,1) = data->u_S.segment<2>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,2,1) = data->u_S.segment<2>(2); 
 
  Eigen::Map<VectorXd>(uCrossOther,2,1) = data->u_S.segment<2>(4); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,2,1) = data->u_S.segment<2>(6); 
 
} 
 
void GkLBOECrossPrimMoments1x1vSer_P1(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // m0S,m1S,m1S:        star moments (only used for piecewise linear). 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg;
  cellAvg = false;
  if (0.7071067811865475*m0Self[0]-1.224744871391589*m0Self[1] < 0) cellAvg = true; 
  if (0.7071067811865475*m2SSelf[0]-1.224744871391589*m2SSelf[1] < 0) cellAvg = true; 
  if (1.224744871391589*m0Self[1]+0.7071067811865475*m0Self[0] < 0) cellAvg = true; 
  if (1.224744871391589*m2SSelf[1]+0.7071067811865475*m2SSelf[0] < 0) cellAvg = true; 
 
  double m0rSelf[2]; 
  double m1rSelf[2]; 
  double m2rSelf[2]; 
  double cMrSelf[2]; 
  double cErSelf[2]; 
  double m0SrSelf[2]; 
  double m1SrSelf[2]; 
  double m2SrSelf[2]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = 0.0; 
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = 0.0; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = 0.0; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = m1Self[1]; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = m2Self[1]; 
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = cMSelf[1]; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = cESelf[1]; 
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = m0SSelf[1]; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = m1SSelf[1]; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = m2SSelf[1]; 
  } 
 
  cellAvg = false;
  if (0.7071067811865475*m0Other[0]-1.224744871391589*m0Other[1] < 0) cellAvg = true; 
  if (0.7071067811865475*m2SOther[0]-1.224744871391589*m2SOther[1] < 0) cellAvg = true; 
  if (1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) cellAvg = true; 
  if (1.224744871391589*m2SOther[1]+0.7071067811865475*m2SOther[0] < 0) cellAvg = true; 
 
  double m0rOther[2]; 
  double m1rOther[2]; 
  double m2rOther[2]; 
  double cMrOther[2]; 
  double cErOther[2]; 
  double m0SrOther[2]; 
  double m1SrOther[2]; 
  double m2SrOther[2]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = 0.0; 
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = 0.0; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = 0.0; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = m1Other[1]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = cMOther[1]; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = cEOther[1]; 
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = m0SOther[1]; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = m1SOther[1]; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak system. 
  data->AEM_S.setZero(); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double mnuM1sum[2]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    mnuM1sum[vd] = 0.0; 
  } 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,2) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
  data->AEM_S(0,3) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(1,2) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(1,3) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,0) += 0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(0,1) += 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,0) += 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,1) += 0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,2) += -(0.7071067811865475*cMrOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(0,3) += -(0.7071067811865475*cMrOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(1,2) += -(0.7071067811865475*cMrOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(1,3) += -(0.7071067811865475*cMrOther[0]*mSelf*mnuOther)/mOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(2,0) = 0.7071067811865475*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(2,1) = 0.7071067811865475*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(3,0) = 0.7071067811865475*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(3,1) = 0.7071067811865475*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(2,0) += 0.7071067811865475*m1SrOther[0]*mnuOther; 
  data->AEM_S(2,1) += 0.7071067811865475*m1SrOther[1]*mnuOther; 
  data->AEM_S(3,0) += 0.7071067811865475*m1SrOther[1]*mnuOther; 
  data->AEM_S(3,1) += 0.7071067811865475*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(2,2) = 0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(3,2) = 0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(3,3) = 0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(2,2) += (0.7071067811865475*m0SrOther[0]*mSelf*mnuOther)/mOther-(0.7071067811865475*cErOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(2,3) += (0.7071067811865475*m0SrOther[1]*mSelf*mnuOther)/mOther-(0.7071067811865475*cErOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(3,2) += (0.7071067811865475*m0SrOther[1]*mSelf*mnuOther)/mOther-(0.7071067811865475*cErOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(3,3) += (0.7071067811865475*m0SrOther[0]*mSelf*mnuOther)/mOther-(0.7071067811865475*cErOther[0]*mSelf*mnuOther)/mOther; 
 
  double mnuM2sum[2]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2SrSelf[0]*mnuSelf+m2SrOther[0]*mnuOther; 
  mnuM2sum[1] = m2SrSelf[1]*mnuSelf+m2SrOther[1]*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM2sum[0],mnuM2sum[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,2,1) = data->u_S.segment<2>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,2,1) = data->u_S.segment<2>(2); 
 
  uCrossOther[0] = uCrossSelf[0]; 
  uCrossOther[1] = uCrossSelf[1]; 
  vtSqCrossOther[0] = (vtSqCrossSelf[0]*mSelf)/mOther; 
  vtSqCrossOther[1] = (vtSqCrossSelf[1]*mSelf)/mOther; 
} 
 
