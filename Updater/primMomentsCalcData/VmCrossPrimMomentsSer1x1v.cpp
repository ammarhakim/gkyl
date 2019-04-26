#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmCrossPrimMoments1x1vSer_P1(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m0SOther, const double *m1SOther, const double *m2SOther, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // m0S,m1S,m1S:        star moments (only used for piecewise linear). 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*m0Self[0]-1.224744871391589*m0Self[1] < 0) { 
    cellAvg = true;
  }
  if (1.224744871391589*m0Self[1]+0.7071067811865475*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[2]; 
  double m1rSelf[2]; 
  double m0SrSelf[2]; 
  double m1SrSelf[2]; 
  double m2SrSelf[2]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
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
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = m0SSelf[1]; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = m1SSelf[1]; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = m2SSelf[1]; 
  } 
 
  if (0.7071067811865475*m0Other[0]-1.224744871391589*m0Other[1] < 0) { 
    cellAvg = true;
  }
  if (1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[2]; 
  double m1rOther[2]; 
  double m0SrOther[2]; 
  double m1SrOther[2]; 
  double m2SrOther[2]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
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
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = m0SOther[1]; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = m1SOther[1]; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(8,8); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
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
  data->AEM_S(0,2) = -0.7071067811865475*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,3) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,2) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,3) = -0.7071067811865475*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,4) = 0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(0,5) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,4) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,5) = 0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,6) = -0.7071067811865475*cMOther[0]*mnuOther; 
  data->AEM_S(0,7) = -0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(1,6) = -0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(1,7) = -0.7071067811865475*cMOther[0]*mnuOther; 
 
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
  data->AEM_S(2,2) = 0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(3,2) = 0.7071067811865475*m0SrSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(3,3) = 0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(2,6) = 0.7071067811865475*m0SrOther[0]*mnuOther-0.7071067811865475*cEOther[0]*mnuOther; 
  data->AEM_S(2,7) = 0.7071067811865475*m0SrOther[1]*mnuOther-0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(3,6) = 0.7071067811865475*m0SrOther[1]*mnuOther-0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(3,7) = 0.7071067811865475*m0SrOther[0]*mnuOther-0.7071067811865475*cEOther[0]*mnuOther; 
 
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
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(4,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,1) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(4,2) = -0.7071067811865475*cMSelf[0]*mnuSelf; 
  data->AEM_S(4,3) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(5,2) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(5,3) = -0.7071067811865475*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(4,4) = -0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(4,5) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(5,4) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(5,5) = -0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(4,6) = 0.7071067811865475*cMOther[0]*mnuOther; 
  data->AEM_S(4,7) = 0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(5,6) = 0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(5,7) = 0.7071067811865475*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(6,0) = (-0.5*m0SrSelf[1]*uSelf[1]*mnuSelf)-0.5*m0SrSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(6,1) = (-0.5*m0SrSelf[0]*uSelf[1]*mnuSelf)+0.7071067811865475*m1SrSelf[1]*mnuSelf-0.5*uSelf[0]*m0SrSelf[1]*mnuSelf; 
  data->AEM_S(7,0) = (-0.5*m0SrSelf[0]*uSelf[1]*mnuSelf)+0.7071067811865475*m1SrSelf[1]*mnuSelf-0.5*uSelf[0]*m0SrSelf[1]*mnuSelf; 
  data->AEM_S(7,1) = (-0.9*m0SrSelf[1]*uSelf[1]*mnuSelf)-0.5*m0SrSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(6,4) = 0.5*m0SrOther[1]*uOther[1]*mnuOther+0.5*m0SrOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1SrOther[0]*mnuOther; 
  data->AEM_S(6,5) = 0.5*m0SrOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1SrOther[1]*mnuOther+0.5*uOther[0]*m0SrOther[1]*mnuOther; 
  data->AEM_S(7,4) = 0.5*m0SrOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1SrOther[1]*mnuOther+0.5*uOther[0]*m0SrOther[1]*mnuOther; 
  data->AEM_S(7,5) = 0.9*m0SrOther[1]*uOther[1]*mnuOther+0.5*m0SrOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(6,2) = 0.5*cMSelf[1]*uSelf[1]*mnuSelf+0.5*cMSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
  data->AEM_S(6,3) = 0.5*cMSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m0SrSelf[1]*mnuSelf+0.5*uSelf[0]*cMSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(7,2) = 0.5*cMSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m0SrSelf[1]*mnuSelf+0.5*uSelf[0]*cMSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(7,3) = 0.5*cMSelf[1]*uSelf[1]*mnuSelf+0.5*cMSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m0SrSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(6,6) = (-0.5*cMOther[1]*uOther[1]*mnuOther)-0.5*cMOther[0]*uOther[0]*mnuOther-0.7071067811865475*m0SrOther[0]*mnuOther+0.7071067811865475*cEOther[0]*mnuOther; 
  data->AEM_S(6,7) = (-0.5*cMOther[0]*uOther[1]*mnuOther)-0.7071067811865475*m0SrOther[1]*mnuOther-0.5*uOther[0]*cMOther[1]*mnuOther+0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(7,6) = (-0.5*cMOther[0]*uOther[1]*mnuOther)-0.7071067811865475*m0SrOther[1]*mnuOther-0.5*uOther[0]*cMOther[1]*mnuOther+0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(7,7) = (-0.5*cMOther[1]*uOther[1]*mnuOther)-0.5*cMOther[0]*uOther[0]*mnuOther-0.7071067811865475*m0SrOther[0]*mnuOther+0.7071067811865475*cEOther[0]*mnuOther; 
 
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
    kinESelf[0] += 0.7071067811865475*m1SrSelf[a0+1]*uSelf[a0+1]+0.7071067811865475*m1SrSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.7071067811865475*m1SrSelf[a0]*uSelf[a0+1]+0.7071067811865475*uSelf[a0]*m1SrSelf[a0+1]; 
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
    kinEOther[0] += 0.7071067811865475*m1SrOther[a0+1]*uOther[a0+1]+0.7071067811865475*m1SrOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.7071067811865475*m1SrOther[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*m1SrOther[a0+1]; 
  } 
 
  double relKinE[2]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.7071067811865475*m1SrSelf[a0+1]*uSelf[a0+1]-0.7071067811865475*m1SrOther[a0+1]*uSelf[a0+1]-0.7071067811865475*m1SrSelf[a0+1]*uOther[a0+1]+0.7071067811865475*m1SrOther[a0+1]*uOther[a0+1]+0.7071067811865475*m1SrSelf[a0]*uSelf[a0]-0.7071067811865475*m1SrOther[a0]*uSelf[a0]-0.7071067811865475*m1SrSelf[a0]*uOther[a0]+0.7071067811865475*m1SrOther[a0]*uOther[a0]; 
    relKinE[1] += 0.7071067811865475*m1SrSelf[a0]*uSelf[a0+1]-0.7071067811865475*m1SrOther[a0]*uSelf[a0+1]-0.7071067811865475*m1SrSelf[a0]*uOther[a0+1]+0.7071067811865475*m1SrOther[a0]*uOther[a0+1]+0.7071067811865475*uSelf[a0]*m1SrSelf[a0+1]-0.7071067811865475*uOther[a0]*m1SrSelf[a0+1]-0.7071067811865475*uSelf[a0]*m1SrOther[a0+1]+0.7071067811865475*uOther[a0]*m1SrOther[a0+1]; 
  } 
 
  double m2Relax[2]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(0.5*relKinE[0]*mSelf)/(mSelf+mOther))-(1.0*m2SrSelf[0]*mSelf)/(mSelf+mOther)+(kinESelf[0]*mSelf)/(mSelf+mOther)+(0.5*relKinE[0]*mOther)/(mSelf+mOther)+(m2SrOther[0]*mOther)/(mSelf+mOther)-(1.0*kinEOther[0]*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2SrOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(0.5*relKinE[1]*mSelf)/(mSelf+mOther))-(1.0*m2SrSelf[1]*mSelf)/(mSelf+mOther)+(kinESelf[1]*mSelf)/(mSelf+mOther)+(0.5*relKinE[1]*mOther)/(mSelf+mOther)+(m2SrOther[1]*mOther)/(mSelf+mOther)-(1.0*kinEOther[1]*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2SrOther[1])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM2sum[0],mnuM2sum[1],m1Relax[0],m1Relax[1],m2Relax[0],m2Relax[1]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,2,1) = data->u_S.segment<2>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,2,1) = data->u_S.segment<2>(2); 
 
  Eigen::Map<VectorXd>(uCrossOther,2,1) = data->u_S.segment<2>(4); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,2,1) = data->u_S.segment<2>(6); 
 
} 
 
void VmCrossPrimMoments1x1vSer_P2(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.58113883008419*m0Self[2]-1.224744871391589*m0Self[1]+0.7071067811865475*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.58113883008419*m0Self[2]+1.224744871391589*m0Self[1]+0.7071067811865475*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[3]; 
  double m1rSelf[3]; 
  double m2rSelf[3]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    m1rSelf[2] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = m1Self[1]; 
    m1rSelf[2] = m1Self[2]; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = m2Self[1]; 
    m2rSelf[2] = m2Self[2]; 
  } 
 
  if (1.58113883008419*m0Other[2]-1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.58113883008419*m0Other[2]+1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[3]; 
  double m1rOther[3]; 
  double m2rOther[3]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    m1rOther[2] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = m1Other[1]; 
    m1rOther[2] = m1Other[2]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(12,12); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
  double mnuM1sum[3]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    mnuM1sum[vd] = 0.0; 
  } 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(0,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(0,2) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(1,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(1,2) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,0) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,1) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,2) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,3) = -0.7071067811865475*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,4) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,5) = -0.7071067811865475*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,3) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,4) = (-0.6324555320336759*cMSelf[2]*mnuSelf)-0.7071067811865475*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,5) = -0.6324555320336759*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,3) = -0.7071067811865475*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,4) = -0.6324555320336759*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,5) = (-0.4517539514526256*cMSelf[2]*mnuSelf)-0.7071067811865475*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,6) = 0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(0,7) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(0,8) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(1,6) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,7) = 0.6324555320336759*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(1,8) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(2,6) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(2,7) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(2,8) = 0.4517539514526256*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,9) = -0.7071067811865475*cMOther[0]*mnuOther; 
  data->AEM_S(0,10) = -0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(0,11) = -0.7071067811865475*cMOther[2]*mnuOther; 
  data->AEM_S(1,9) = -0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(1,10) = (-0.6324555320336759*cMOther[2]*mnuOther)-0.7071067811865475*cMOther[0]*mnuOther; 
  data->AEM_S(1,11) = -0.6324555320336759*cMOther[1]*mnuOther; 
  data->AEM_S(2,9) = -0.7071067811865475*cMOther[2]*mnuOther; 
  data->AEM_S(2,10) = -0.6324555320336759*cMOther[1]*mnuOther; 
  data->AEM_S(2,11) = (-0.4517539514526256*cMOther[2]*mnuOther)-0.7071067811865475*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(3,0) = 0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(3,1) = 0.7071067811865475*m1rSelf[1]*mnuSelf; 
  data->AEM_S(3,2) = 0.7071067811865475*m1rSelf[2]*mnuSelf; 
  data->AEM_S(4,0) = 0.7071067811865475*m1rSelf[1]*mnuSelf; 
  data->AEM_S(4,1) = 0.6324555320336759*m1rSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(4,2) = 0.6324555320336759*m1rSelf[1]*mnuSelf; 
  data->AEM_S(5,0) = 0.7071067811865475*m1rSelf[2]*mnuSelf; 
  data->AEM_S(5,1) = 0.6324555320336759*m1rSelf[1]*mnuSelf; 
  data->AEM_S(5,2) = 0.4517539514526256*m1rSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(3,6) = 0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(3,7) = 0.7071067811865475*m1rOther[1]*mnuOther; 
  data->AEM_S(3,8) = 0.7071067811865475*m1rOther[2]*mnuOther; 
  data->AEM_S(4,6) = 0.7071067811865475*m1rOther[1]*mnuOther; 
  data->AEM_S(4,7) = 0.6324555320336759*m1rOther[2]*mnuOther+0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(4,8) = 0.6324555320336759*m1rOther[1]*mnuOther; 
  data->AEM_S(5,6) = 0.7071067811865475*m1rOther[2]*mnuOther; 
  data->AEM_S(5,7) = 0.6324555320336759*m1rOther[1]*mnuOther; 
  data->AEM_S(5,8) = 0.4517539514526256*m1rOther[2]*mnuOther+0.7071067811865475*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(3,3) = 0.7071067811865475*m0rSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
  data->AEM_S(3,4) = 0.7071067811865475*m0rSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(3,5) = 0.7071067811865475*m0rSelf[2]*mnuSelf-0.7071067811865475*cESelf[2]*mnuSelf; 
  data->AEM_S(4,3) = 0.7071067811865475*m0rSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(4,4) = 0.6324555320336759*m0rSelf[2]*mnuSelf-0.6324555320336759*cESelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
  data->AEM_S(4,5) = 0.6324555320336759*m0rSelf[1]*mnuSelf-0.6324555320336759*cESelf[1]*mnuSelf; 
  data->AEM_S(5,3) = 0.7071067811865475*m0rSelf[2]*mnuSelf-0.7071067811865475*cESelf[2]*mnuSelf; 
  data->AEM_S(5,4) = 0.6324555320336759*m0rSelf[1]*mnuSelf-0.6324555320336759*cESelf[1]*mnuSelf; 
  data->AEM_S(5,5) = 0.4517539514526256*m0rSelf[2]*mnuSelf-0.4517539514526256*cESelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(3,9) = 0.7071067811865475*m0rOther[0]*mnuOther-0.7071067811865475*cEOther[0]*mnuOther; 
  data->AEM_S(3,10) = 0.7071067811865475*m0rOther[1]*mnuOther-0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(3,11) = 0.7071067811865475*m0rOther[2]*mnuOther-0.7071067811865475*cEOther[2]*mnuOther; 
  data->AEM_S(4,9) = 0.7071067811865475*m0rOther[1]*mnuOther-0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(4,10) = 0.6324555320336759*m0rOther[2]*mnuOther-0.6324555320336759*cEOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther-0.7071067811865475*cEOther[0]*mnuOther; 
  data->AEM_S(4,11) = 0.6324555320336759*m0rOther[1]*mnuOther-0.6324555320336759*cEOther[1]*mnuOther; 
  data->AEM_S(5,9) = 0.7071067811865475*m0rOther[2]*mnuOther-0.7071067811865475*cEOther[2]*mnuOther; 
  data->AEM_S(5,10) = 0.6324555320336759*m0rOther[1]*mnuOther-0.6324555320336759*cEOther[1]*mnuOther; 
  data->AEM_S(5,11) = 0.4517539514526256*m0rOther[2]*mnuOther-0.4517539514526256*cEOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther-0.7071067811865475*cEOther[0]*mnuOther; 
 
  double mnuM2sum[3]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2rSelf[0]*mnuSelf+m2rOther[0]*mnuOther; 
  mnuM2sum[1] = m2rSelf[1]*mnuSelf+m2rOther[1]*mnuOther; 
  mnuM2sum[2] = m2rSelf[2]*mnuSelf+m2rOther[2]*mnuOther; 
 
  double m1Relax[3]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(6,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,2) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,1) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(7,2) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,0) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,1) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,2) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(6,3) = -0.7071067811865475*cMSelf[0]*mnuSelf; 
  data->AEM_S(6,4) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,5) = -0.7071067811865475*cMSelf[2]*mnuSelf; 
  data->AEM_S(7,3) = -0.7071067811865475*cMSelf[1]*mnuSelf; 
  data->AEM_S(7,4) = (-0.6324555320336759*cMSelf[2]*mnuSelf)-0.7071067811865475*cMSelf[0]*mnuSelf; 
  data->AEM_S(7,5) = -0.6324555320336759*cMSelf[1]*mnuSelf; 
  data->AEM_S(8,3) = -0.7071067811865475*cMSelf[2]*mnuSelf; 
  data->AEM_S(8,4) = -0.6324555320336759*cMSelf[1]*mnuSelf; 
  data->AEM_S(8,5) = (-0.4517539514526256*cMSelf[2]*mnuSelf)-0.7071067811865475*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(6,6) = -0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(6,7) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(6,8) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(7,6) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(7,7) = (-0.6324555320336759*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(7,8) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(8,6) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(8,7) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(8,8) = (-0.4517539514526256*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(6,9) = 0.7071067811865475*cMOther[0]*mnuOther; 
  data->AEM_S(6,10) = 0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(6,11) = 0.7071067811865475*cMOther[2]*mnuOther; 
  data->AEM_S(7,9) = 0.7071067811865475*cMOther[1]*mnuOther; 
  data->AEM_S(7,10) = 0.6324555320336759*cMOther[2]*mnuOther+0.7071067811865475*cMOther[0]*mnuOther; 
  data->AEM_S(7,11) = 0.6324555320336759*cMOther[1]*mnuOther; 
  data->AEM_S(8,9) = 0.7071067811865475*cMOther[2]*mnuOther; 
  data->AEM_S(8,10) = 0.6324555320336759*cMOther[1]*mnuOther; 
  data->AEM_S(8,11) = 0.4517539514526256*cMOther[2]*mnuOther+0.7071067811865475*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(9,0) = (-0.5*m0rSelf[2]*uSelf[2]*mnuSelf)-0.5*m0rSelf[1]*uSelf[1]*mnuSelf-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(9,1) = (-0.4472135954999579*m0rSelf[1]*uSelf[2]*mnuSelf)-0.4472135954999579*uSelf[1]*m0rSelf[2]*mnuSelf-0.5*m0rSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m1rSelf[1]*mnuSelf-0.5*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,2) = (-0.31943828249997*m0rSelf[2]*uSelf[2]*mnuSelf)-0.5*m0rSelf[0]*uSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[2]*mnuSelf-0.5*uSelf[0]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(10,0) = (-0.4472135954999579*m0rSelf[1]*uSelf[2]*mnuSelf)-0.4472135954999579*uSelf[1]*m0rSelf[2]*mnuSelf-0.5*m0rSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m1rSelf[1]*mnuSelf-0.5*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(10,1) = (-0.7857142857142857*m0rSelf[2]*uSelf[2]*mnuSelf)-0.4472135954999579*m0rSelf[0]*uSelf[2]*mnuSelf+0.6324555320336759*m1rSelf[2]*mnuSelf-0.4472135954999579*uSelf[0]*m0rSelf[2]*mnuSelf-0.9*m0rSelf[1]*uSelf[1]*mnuSelf-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(10,2) = (-0.7857142857142857*m0rSelf[1]*uSelf[2]*mnuSelf)-0.7857142857142857*uSelf[1]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[1]*mnuSelf+0.6324555320336759*m1rSelf[1]*mnuSelf-0.4472135954999579*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,0) = (-0.31943828249997*m0rSelf[2]*uSelf[2]*mnuSelf)-0.5*m0rSelf[0]*uSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[2]*mnuSelf-0.5*uSelf[0]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(11,1) = (-0.7857142857142857*m0rSelf[1]*uSelf[2]*mnuSelf)-0.7857142857142857*uSelf[1]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[1]*mnuSelf+0.6324555320336759*m1rSelf[1]*mnuSelf-0.4472135954999579*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,2) = (-1.071428571428571*m0rSelf[2]*uSelf[2]*mnuSelf)-0.31943828249997*m0rSelf[0]*uSelf[2]*mnuSelf+0.4517539514526256*m1rSelf[2]*mnuSelf-0.31943828249997*uSelf[0]*m0rSelf[2]*mnuSelf-0.7857142857142857*m0rSelf[1]*uSelf[1]*mnuSelf-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(9,6) = 0.5*m0rOther[2]*uOther[2]*mnuOther+0.5*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(9,7) = 0.4472135954999579*m0rOther[1]*uOther[2]*mnuOther+0.4472135954999579*uOther[1]*m0rOther[2]*mnuOther+0.5*m0rOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1rOther[1]*mnuOther+0.5*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(9,8) = 0.31943828249997*m0rOther[2]*uOther[2]*mnuOther+0.5*m0rOther[0]*uOther[2]*mnuOther-0.7071067811865475*m1rOther[2]*mnuOther+0.5*uOther[0]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(10,6) = 0.4472135954999579*m0rOther[1]*uOther[2]*mnuOther+0.4472135954999579*uOther[1]*m0rOther[2]*mnuOther+0.5*m0rOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1rOther[1]*mnuOther+0.5*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(10,7) = 0.7857142857142857*m0rOther[2]*uOther[2]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[2]*mnuOther-0.6324555320336759*m1rOther[2]*mnuOther+0.4472135954999579*uOther[0]*m0rOther[2]*mnuOther+0.9*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(10,8) = 0.7857142857142857*m0rOther[1]*uOther[2]*mnuOther+0.7857142857142857*uOther[1]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[1]*mnuOther-0.6324555320336759*m1rOther[1]*mnuOther+0.4472135954999579*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(11,6) = 0.31943828249997*m0rOther[2]*uOther[2]*mnuOther+0.5*m0rOther[0]*uOther[2]*mnuOther-0.7071067811865475*m1rOther[2]*mnuOther+0.5*uOther[0]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(11,7) = 0.7857142857142857*m0rOther[1]*uOther[2]*mnuOther+0.7857142857142857*uOther[1]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[1]*mnuOther-0.6324555320336759*m1rOther[1]*mnuOther+0.4472135954999579*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(11,8) = 1.071428571428571*m0rOther[2]*uOther[2]*mnuOther+0.31943828249997*m0rOther[0]*uOther[2]*mnuOther-0.4517539514526256*m1rOther[2]*mnuOther+0.31943828249997*uOther[0]*m0rOther[2]*mnuOther+0.7857142857142857*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(9,3) = 0.5*cMSelf[2]*uSelf[2]*mnuSelf+0.5*cMSelf[1]*uSelf[1]*mnuSelf+0.5*cMSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
  data->AEM_S(9,4) = 0.4472135954999579*cMSelf[1]*uSelf[2]*mnuSelf+0.4472135954999579*uSelf[1]*cMSelf[2]*mnuSelf+0.5*cMSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m0rSelf[1]*mnuSelf+0.5*uSelf[0]*cMSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(9,5) = 0.31943828249997*cMSelf[2]*uSelf[2]*mnuSelf+0.5*cMSelf[0]*uSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[2]*mnuSelf+0.5*uSelf[0]*cMSelf[2]*mnuSelf-0.7071067811865475*cESelf[2]*mnuSelf+0.4472135954999579*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(10,3) = 0.4472135954999579*cMSelf[1]*uSelf[2]*mnuSelf+0.4472135954999579*uSelf[1]*cMSelf[2]*mnuSelf+0.5*cMSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m0rSelf[1]*mnuSelf+0.5*uSelf[0]*cMSelf[1]*mnuSelf-0.7071067811865475*cESelf[1]*mnuSelf; 
  data->AEM_S(10,4) = 0.7857142857142857*cMSelf[2]*uSelf[2]*mnuSelf+0.4472135954999579*cMSelf[0]*uSelf[2]*mnuSelf+0.6324555320336759*m0rSelf[2]*mnuSelf+0.4472135954999579*uSelf[0]*cMSelf[2]*mnuSelf-0.6324555320336759*cESelf[2]*mnuSelf+0.9*cMSelf[1]*uSelf[1]*mnuSelf+0.5*cMSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
  data->AEM_S(10,5) = 0.4*cMSelf[1]*uSelf[2]*mnuSelf+0.4*uSelf[1]*cMSelf[2]*mnuSelf+0.4472135954999579*cMSelf[0]*uSelf[1]*mnuSelf+0.6324555320336759*m0rSelf[1]*mnuSelf+0.4472135954999579*uSelf[0]*cMSelf[1]*mnuSelf-0.6324555320336759*cESelf[1]*mnuSelf; 
  data->AEM_S(11,3) = 0.31943828249997*cMSelf[2]*uSelf[2]*mnuSelf+0.5*cMSelf[0]*uSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[2]*mnuSelf+0.5*uSelf[0]*cMSelf[2]*mnuSelf-0.7071067811865475*cESelf[2]*mnuSelf+0.4472135954999579*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(11,4) = 0.4*cMSelf[1]*uSelf[2]*mnuSelf+0.4*uSelf[1]*cMSelf[2]*mnuSelf+0.4472135954999579*cMSelf[0]*uSelf[1]*mnuSelf+0.6324555320336759*m0rSelf[1]*mnuSelf+0.4472135954999579*uSelf[0]*cMSelf[1]*mnuSelf-0.6324555320336759*cESelf[1]*mnuSelf; 
  data->AEM_S(11,5) = 0.7040816326530612*cMSelf[2]*uSelf[2]*mnuSelf+0.31943828249997*cMSelf[0]*uSelf[2]*mnuSelf+0.4517539514526256*m0rSelf[2]*mnuSelf+0.31943828249997*uSelf[0]*cMSelf[2]*mnuSelf-0.4517539514526256*cESelf[2]*mnuSelf+0.7857142857142857*cMSelf[1]*uSelf[1]*mnuSelf+0.5*cMSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf-0.7071067811865475*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(9,9) = (-0.5*cMOther[2]*uOther[2]*mnuOther)-0.5*cMOther[1]*uOther[1]*mnuOther-0.5*cMOther[0]*uOther[0]*mnuOther-0.7071067811865475*m0rOther[0]*mnuOther+0.7071067811865475*cEOther[0]*mnuOther; 
  data->AEM_S(9,10) = (-0.4472135954999579*cMOther[1]*uOther[2]*mnuOther)-0.4472135954999579*uOther[1]*cMOther[2]*mnuOther-0.5*cMOther[0]*uOther[1]*mnuOther-0.7071067811865475*m0rOther[1]*mnuOther-0.5*uOther[0]*cMOther[1]*mnuOther+0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(9,11) = (-0.31943828249997*cMOther[2]*uOther[2]*mnuOther)-0.5*cMOther[0]*uOther[2]*mnuOther-0.7071067811865475*m0rOther[2]*mnuOther-0.5*uOther[0]*cMOther[2]*mnuOther+0.7071067811865475*cEOther[2]*mnuOther-0.4472135954999579*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(10,9) = (-0.4472135954999579*cMOther[1]*uOther[2]*mnuOther)-0.4472135954999579*uOther[1]*cMOther[2]*mnuOther-0.5*cMOther[0]*uOther[1]*mnuOther-0.7071067811865475*m0rOther[1]*mnuOther-0.5*uOther[0]*cMOther[1]*mnuOther+0.7071067811865475*cEOther[1]*mnuOther; 
  data->AEM_S(10,10) = (-0.7857142857142857*cMOther[2]*uOther[2]*mnuOther)-0.4472135954999579*cMOther[0]*uOther[2]*mnuOther-0.6324555320336759*m0rOther[2]*mnuOther-0.4472135954999579*uOther[0]*cMOther[2]*mnuOther+0.6324555320336759*cEOther[2]*mnuOther-0.9*cMOther[1]*uOther[1]*mnuOther-0.5*cMOther[0]*uOther[0]*mnuOther-0.7071067811865475*m0rOther[0]*mnuOther+0.7071067811865475*cEOther[0]*mnuOther; 
  data->AEM_S(10,11) = (-0.4*cMOther[1]*uOther[2]*mnuOther)-0.4*uOther[1]*cMOther[2]*mnuOther-0.4472135954999579*cMOther[0]*uOther[1]*mnuOther-0.6324555320336759*m0rOther[1]*mnuOther-0.4472135954999579*uOther[0]*cMOther[1]*mnuOther+0.6324555320336759*cEOther[1]*mnuOther; 
  data->AEM_S(11,9) = (-0.31943828249997*cMOther[2]*uOther[2]*mnuOther)-0.5*cMOther[0]*uOther[2]*mnuOther-0.7071067811865475*m0rOther[2]*mnuOther-0.5*uOther[0]*cMOther[2]*mnuOther+0.7071067811865475*cEOther[2]*mnuOther-0.4472135954999579*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(11,10) = (-0.4*cMOther[1]*uOther[2]*mnuOther)-0.4*uOther[1]*cMOther[2]*mnuOther-0.4472135954999579*cMOther[0]*uOther[1]*mnuOther-0.6324555320336759*m0rOther[1]*mnuOther-0.4472135954999579*uOther[0]*cMOther[1]*mnuOther+0.6324555320336759*cEOther[1]*mnuOther; 
  data->AEM_S(11,11) = (-0.7040816326530612*cMOther[2]*uOther[2]*mnuOther)-0.31943828249997*cMOther[0]*uOther[2]*mnuOther-0.4517539514526256*m0rOther[2]*mnuOther-0.31943828249997*uOther[0]*cMOther[2]*mnuOther+0.4517539514526256*cEOther[2]*mnuOther-0.7857142857142857*cMOther[1]*uOther[1]*mnuOther-0.5*cMOther[0]*uOther[0]*mnuOther-0.7071067811865475*m0rOther[0]*mnuOther+0.7071067811865475*cEOther[0]*mnuOther; 
 
  double kinESelf[3]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinESelf[0] += 0.7071067811865475*m1rSelf[a0+2]*uSelf[a0+2]+0.7071067811865475*m1rSelf[a0+1]*uSelf[a0+1]+0.7071067811865475*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.6324555320336759*m1rSelf[a0+1]*uSelf[a0+2]+0.6324555320336759*uSelf[a0+1]*m1rSelf[a0+2]+0.7071067811865475*m1rSelf[a0]*uSelf[a0+1]+0.7071067811865475*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.4517539514526256*m1rSelf[a0+2]*uSelf[a0+2]+0.7071067811865475*m1rSelf[a0]*uSelf[a0+2]+0.7071067811865475*uSelf[a0]*m1rSelf[a0+2]+0.6324555320336759*m1rSelf[a0+1]*uSelf[a0+1]; 
  } 
 
  double kinEOther[3]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    kinEOther[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.7071067811865475*m1rOther[a0+2]*uOther[a0+2]+0.7071067811865475*m1rOther[a0+1]*uOther[a0+1]+0.7071067811865475*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.6324555320336759*m1rOther[a0+1]*uOther[a0+2]+0.6324555320336759*uOther[a0+1]*m1rOther[a0+2]+0.7071067811865475*m1rOther[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.4517539514526256*m1rOther[a0+2]*uOther[a0+2]+0.7071067811865475*m1rOther[a0]*uOther[a0+2]+0.7071067811865475*uOther[a0]*m1rOther[a0+2]+0.6324555320336759*m1rOther[a0+1]*uOther[a0+1]; 
  } 
 
  double relKinE[3]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.7071067811865475*m1rSelf[a0+2]*uSelf[a0+2]-0.7071067811865475*m1rOther[a0+2]*uSelf[a0+2]-0.7071067811865475*m1rSelf[a0+2]*uOther[a0+2]+0.7071067811865475*m1rOther[a0+2]*uOther[a0+2]+0.7071067811865475*m1rSelf[a0+1]*uSelf[a0+1]-0.7071067811865475*m1rOther[a0+1]*uSelf[a0+1]-0.7071067811865475*m1rSelf[a0+1]*uOther[a0+1]+0.7071067811865475*m1rOther[a0+1]*uOther[a0+1]+0.7071067811865475*m1rSelf[a0]*uSelf[a0]-0.7071067811865475*m1rOther[a0]*uSelf[a0]-0.7071067811865475*m1rSelf[a0]*uOther[a0]+0.7071067811865475*m1rOther[a0]*uOther[a0]; 
    relKinE[1] += 0.6324555320336759*m1rSelf[a0+1]*uSelf[a0+2]-0.6324555320336759*m1rOther[a0+1]*uSelf[a0+2]-0.6324555320336759*m1rSelf[a0+1]*uOther[a0+2]+0.6324555320336759*m1rOther[a0+1]*uOther[a0+2]+0.6324555320336759*uSelf[a0+1]*m1rSelf[a0+2]-0.6324555320336759*uOther[a0+1]*m1rSelf[a0+2]-0.6324555320336759*uSelf[a0+1]*m1rOther[a0+2]+0.6324555320336759*uOther[a0+1]*m1rOther[a0+2]+0.7071067811865475*m1rSelf[a0]*uSelf[a0+1]-0.7071067811865475*m1rOther[a0]*uSelf[a0+1]-0.7071067811865475*m1rSelf[a0]*uOther[a0+1]+0.7071067811865475*m1rOther[a0]*uOther[a0+1]+0.7071067811865475*uSelf[a0]*m1rSelf[a0+1]-0.7071067811865475*uOther[a0]*m1rSelf[a0+1]-0.7071067811865475*uSelf[a0]*m1rOther[a0+1]+0.7071067811865475*uOther[a0]*m1rOther[a0+1]; 
    relKinE[2] += 0.4517539514526256*m1rSelf[a0+2]*uSelf[a0+2]-0.4517539514526256*m1rOther[a0+2]*uSelf[a0+2]+0.7071067811865475*m1rSelf[a0]*uSelf[a0+2]-0.7071067811865475*m1rOther[a0]*uSelf[a0+2]-0.4517539514526256*m1rSelf[a0+2]*uOther[a0+2]+0.4517539514526256*m1rOther[a0+2]*uOther[a0+2]-0.7071067811865475*m1rSelf[a0]*uOther[a0+2]+0.7071067811865475*m1rOther[a0]*uOther[a0+2]+0.7071067811865475*uSelf[a0]*m1rSelf[a0+2]-0.7071067811865475*uOther[a0]*m1rSelf[a0+2]-0.7071067811865475*uSelf[a0]*m1rOther[a0+2]+0.7071067811865475*uOther[a0]*m1rOther[a0+2]+0.6324555320336759*m1rSelf[a0+1]*uSelf[a0+1]-0.6324555320336759*m1rOther[a0+1]*uSelf[a0+1]-0.6324555320336759*m1rSelf[a0+1]*uOther[a0+1]+0.6324555320336759*m1rOther[a0+1]*uOther[a0+1]; 
  } 
 
  double m2Relax[3]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(0.5*relKinE[0]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[0]*mSelf)/(mSelf+mOther)+(kinESelf[0]*mSelf)/(mSelf+mOther)+(0.5*relKinE[0]*mOther)/(mSelf+mOther)+(m2rOther[0]*mOther)/(mSelf+mOther)-(1.0*kinEOther[0]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2rOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(0.5*relKinE[1]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[1]*mSelf)/(mSelf+mOther)+(kinESelf[1]*mSelf)/(mSelf+mOther)+(0.5*relKinE[1]*mOther)/(mSelf+mOther)+(m2rOther[1]*mOther)/(mSelf+mOther)-(1.0*kinEOther[1]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2rOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(0.5*relKinE[2]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[2]*mSelf)/(mSelf+mOther)+(kinESelf[2]*mSelf)/(mSelf+mOther)+(0.5*relKinE[2]*mOther)/(mSelf+mOther)+(m2rOther[2]*mOther)/(mSelf+mOther)-(1.0*kinEOther[2]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2rOther[2])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],m1Relax[0],m1Relax[1],m1Relax[2],m2Relax[0],m2Relax[1],m2Relax[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,3,1) = data->u_S.segment<3>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,3,1) = data->u_S.segment<3>(3); 
 
  Eigen::Map<VectorXd>(uCrossOther,3,1) = data->u_S.segment<3>(6); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,3,1) = data->u_S.segment<3>(9); 
 
} 
 
