#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkCrossPrimMoments2x2vMax_P1(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  if ((-0.8660254037844386*m0Self[2])-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Self[2])-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Self[2])+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Self[2])+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[3]; 
  double m1rSelf[3]; 
  double m2rSelf[3]; 
  double m0SrSelf[3]; 
  double m1SrSelf[3]; 
  double m2SrSelf[3]; 
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
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = 0.0; 
    m0SrSelf[2] = 0.0; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = 0.0; 
    m1SrSelf[2] = 0.0; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = 0.0; 
    m2SrSelf[2] = 0.0; 
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
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = m0SSelf[1]; 
    m0SrSelf[2] = m0SSelf[2]; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = m1SSelf[1]; 
    m1SrSelf[2] = m1SSelf[2]; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = m2SSelf[1]; 
    m2SrSelf[2] = m2SSelf[2]; 
  } 
 
  if ((-0.8660254037844386*m0Other[2])-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Other[2])-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Other[2])+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Other[2])+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[3]; 
  double m1rOther[3]; 
  double m2rOther[3]; 
  double m0SrOther[3]; 
  double m1SrOther[3]; 
  double m2SrOther[3]; 
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
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = 0.0; 
    m0SrOther[2] = 0.0; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = 0.0; 
    m1SrOther[2] = 0.0; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = 0.0; 
    m2SrOther[2] = 0.0; 
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
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = m0SOther[1]; 
    m0SrOther[2] = m0SOther[2]; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = m1SOther[1]; 
    m1SrOther[2] = m1SOther[2]; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
    m2SrOther[2] = m2SOther[2]; 
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
  data->AEM_S(0,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(0,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,2) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,3) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,4) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,5) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,3) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,4) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,5) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,6) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,7) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,8) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,6) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,7) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,6) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,8) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,9) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,10) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,11) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,9) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,10) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,9) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,11) = -0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(3,0) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(3,1) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(3,2) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(4,0) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(4,1) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(5,0) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(5,2) = 0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(3,6) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(3,7) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(3,8) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(4,6) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(4,7) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(5,6) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(5,8) = 0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(3,3) = m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(3,4) = m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(3,5) = m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(4,3) = m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(4,4) = m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(5,3) = m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(5,5) = m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(3,9) = m0rOther[0]*mnuOther+0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(3,10) = m0rOther[1]*mnuOther+0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(3,11) = m0rOther[2]*mnuOther+0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(4,9) = m0rOther[1]*mnuOther+0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(4,10) = m0rOther[0]*mnuOther+0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(5,9) = m0rOther[2]*mnuOther+0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(5,11) = m0rOther[0]*mnuOther+0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[3]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2SrSelf[0]*mnuSelf+m2SrOther[0]*mnuOther; 
  mnuM2sum[1] = m2SrSelf[1]*mnuSelf+m2SrOther[1]*mnuOther; 
  mnuM2sum[2] = m2SrSelf[2]*mnuSelf+m2SrOther[2]*mnuOther; 
 
  double m1Relax[3]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(6,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,2) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(6,3) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(6,4) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,5) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(7,3) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(7,4) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(8,3) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(8,5) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(6,6) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,7) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(6,8) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(7,6) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,7) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,6) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,8) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(6,9) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(6,10) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(6,11) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(7,9) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(7,10) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(8,9) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(8,11) = 0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(9,0) = (-0.25*m0rSelf[2]*uSelf[2]*mnuSelf)-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(9,1) = (-0.25*m0rSelf[0]*uSelf[1]*mnuSelf)+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,2) = (-0.25*m0rSelf[0]*uSelf[2]*mnuSelf)+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,0) = (-0.25*m0rSelf[0]*uSelf[1]*mnuSelf)+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(10,1) = (-0.25*m0rSelf[2]*uSelf[2]*mnuSelf)-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(10,2) = (-0.25*m0rSelf[1]*uSelf[2]*mnuSelf)-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,0) = (-0.25*m0rSelf[0]*uSelf[2]*mnuSelf)+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,1) = (-0.25*m0rSelf[1]*uSelf[2]*mnuSelf)-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,2) = (-0.45*m0rSelf[2]*uSelf[2]*mnuSelf)-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(9,6) = 0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(9,7) = 0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(9,8) = 0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(10,6) = 0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(10,7) = 0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(10,8) = 0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(11,6) = 0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(11,7) = 0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(11,8) = 0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
 
  double ucMSelf[3]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.5*cMSelf[a0+2]*uSelf[a0+2]+0.5*cMSelf[a0+1]*uSelf[a0+1]+0.5*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.5*cMSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.5*cMSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*cMSelf[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(9,3) = 0.5*ucMSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(9,4) = 0.5*ucMSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(9,5) = 0.5*ucMSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(10,3) = 0.5*ucMSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(10,4) = 0.5*ucMSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(11,3) = 0.5*ucMSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(11,5) = 0.5*ucMSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  double ucMOther[3]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.5*cMOther[a0+2]*uOther[a0+2]+0.5*cMOther[a0+1]*uOther[a0+1]+0.5*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.5*cMOther[a0]*uOther[a0+1]+0.5*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.5*cMOther[a0]*uOther[a0+2]+0.5*uOther[a0]*cMOther[a0+2]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(9,9) = (-0.5*ucMOther[0]*mnuOther)-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(9,10) = (-0.5*ucMOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(9,11) = (-0.5*ucMOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(10,9) = (-0.5*ucMOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(10,10) = (-0.5*ucMOther[0]*mnuOther)-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(11,9) = (-0.5*ucMOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(11,11) = (-0.5*ucMOther[0]*mnuOther)-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
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
    kinESelf[0] += 0.5*m1rSelf[a0+2]*uSelf[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.5*m1rSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.5*m1rSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]; 
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
    kinEOther[0] += 0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.5*m1rOther[a0]*uOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.5*m1rOther[a0]*uOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
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
    relKinE[0] += 0.5*m1rSelf[a0+2]*uSelf[a0+2]-0.5*m1rOther[a0+2]*uSelf[a0+2]-0.5*m1rSelf[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]-0.5*m1rOther[a0+1]*uSelf[a0+1]-0.5*m1rSelf[a0+1]*uOther[a0+1]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]-0.5*m1rOther[a0]*uSelf[a0]-0.5*m1rSelf[a0]*uOther[a0]+0.5*m1rOther[a0]*uOther[a0]; 
    relKinE[1] += 0.5*m1rSelf[a0]*uSelf[a0+1]-0.5*m1rOther[a0]*uSelf[a0+1]-0.5*m1rSelf[a0]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]-0.5*uOther[a0]*m1rSelf[a0+1]-0.5*uSelf[a0]*m1rOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    relKinE[2] += 0.5*m1rSelf[a0]*uSelf[a0+2]-0.5*m1rOther[a0]*uSelf[a0+2]-0.5*m1rSelf[a0]*uOther[a0+2]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]-0.5*uOther[a0]*m1rSelf[a0+2]-0.5*uSelf[a0]*m1rOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
  } 
 
  double m2Relax[3]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(0.5*relKinE[0]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[0]*mSelf)/(mSelf+mOther)+(kinESelf[0]*mSelf)/(mSelf+mOther)+(0.5*relKinE[0]*mOther)/(mSelf+mOther)+(m2rOther[0]*mOther)/(mSelf+mOther)-(1.0*kinEOther[0]*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2SrOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(0.5*relKinE[1]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[1]*mSelf)/(mSelf+mOther)+(kinESelf[1]*mSelf)/(mSelf+mOther)+(0.5*relKinE[1]*mOther)/(mSelf+mOther)+(m2rOther[1]*mOther)/(mSelf+mOther)-(1.0*kinEOther[1]*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2SrOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(0.5*relKinE[2]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[2]*mSelf)/(mSelf+mOther)+(kinESelf[2]*mSelf)/(mSelf+mOther)+(0.5*relKinE[2]*mOther)/(mSelf+mOther)+(m2rOther[2]*mOther)/(mSelf+mOther)-(1.0*kinEOther[2]*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2SrOther[2])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],m1Relax[0],m1Relax[1],m1Relax[2],m2Relax[0],m2Relax[1],m2Relax[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,3,1) = data->u_S.segment<3>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,3,1) = data->u_S.segment<3>(3); 
 
  Eigen::Map<VectorXd>(uCrossOther,3,1) = data->u_S.segment<3>(6); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,3,1) = data->u_S.segment<3>(9); 
 
} 
 
void GkCrossPrimMoments2x2vMax_P2(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[6]; 
  double m1rSelf[6]; 
  double m2rSelf[6]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m0rSelf[3] = 0.0; 
    m0rSelf[4] = 0.0; 
    m0rSelf[5] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    m1rSelf[2] = 0.0; 
    m1rSelf[3] = 0.0; 
    m1rSelf[4] = 0.0; 
    m1rSelf[5] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    m2rSelf[3] = 0.0; 
    m2rSelf[4] = 0.0; 
    m2rSelf[5] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
    m0rSelf[3] = m0Self[3]; 
    m0rSelf[4] = m0Self[4]; 
    m0rSelf[5] = m0Self[5]; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = m1Self[1]; 
    m1rSelf[2] = m1Self[2]; 
    m1rSelf[3] = m1Self[3]; 
    m1rSelf[4] = m1Self[4]; 
    m1rSelf[5] = m1Self[5]; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = m2Self[1]; 
    m2rSelf[2] = m2Self[2]; 
    m2rSelf[3] = m2Self[3]; 
    m2rSelf[4] = m2Self[4]; 
    m2rSelf[5] = m2Self[5]; 
  } 
 
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[6]; 
  double m1rOther[6]; 
  double m2rOther[6]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m0rOther[3] = 0.0; 
    m0rOther[4] = 0.0; 
    m0rOther[5] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    m1rOther[2] = 0.0; 
    m1rOther[3] = 0.0; 
    m1rOther[4] = 0.0; 
    m1rOther[5] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    m2rOther[3] = 0.0; 
    m2rOther[4] = 0.0; 
    m2rOther[5] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
    m0rOther[3] = m0Other[3]; 
    m0rOther[4] = m0Other[4]; 
    m0rOther[5] = m0Other[5]; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = m1Other[1]; 
    m1rOther[2] = m1Other[2]; 
    m1rOther[3] = m1Other[3]; 
    m1rOther[4] = m1Other[4]; 
    m1rOther[5] = m1Other[5]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
    m2rOther[3] = m2Other[3]; 
    m2rOther[4] = m2Other[4]; 
    m2rOther[5] = m2Other[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(24,24); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
  double mnuM1sum[6]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    mnuM1sum[vd] = 0.0; 
  } 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(0,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(0,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(0,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(0,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(0,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(1,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(1,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(2,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,1) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,2) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(3,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(4,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(4,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(4,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(5,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(5,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(5,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(5,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,6) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,7) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,8) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,9) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,10) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,11) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(1,6) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,7) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,8) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,9) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,10) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,6) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,7) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,8) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,9) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,11) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,6) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,7) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,8) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,9) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,10) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,11) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,6) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,7) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,9) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,10) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,6) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,8) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(5,9) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,11) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,12) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,13) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,14) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,15) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(0,16) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(0,17) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(1,12) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,13) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,14) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,15) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,16) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(2,12) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,13) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,14) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,15) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(2,17) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(3,12) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,13) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,14) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,15) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,16) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,17) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,12) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(4,13) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(4,15) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,16) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,12) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(5,14) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(5,15) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(5,17) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,18) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,19) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,20) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,21) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(0,22) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(0,23) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(1,18) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,19) = (-0.4472135954999579*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,20) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,21) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,22) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(2,18) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,19) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,20) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,21) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(2,23) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(3,18) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,19) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,20) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,21) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(3,22) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,23) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,18) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,19) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(4,21) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,22) = (-0.31943828249997*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(5,18) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,20) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(5,21) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(5,23) = (-0.31943828249997*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(6,0) = 0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(6,1) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(6,2) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(6,3) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(6,4) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(6,5) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(7,0) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(7,1) = 0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(7,2) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(7,3) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(7,4) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(8,0) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(8,1) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(8,2) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(8,3) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(8,5) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(9,0) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(9,1) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(9,2) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(9,3) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(9,4) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(9,5) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(10,0) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(10,1) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(10,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(10,4) = 0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(11,0) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(11,2) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(11,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(11,5) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(6,12) = 0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(6,13) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(6,14) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(6,15) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(6,16) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(6,17) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(7,12) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(7,13) = 0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(7,14) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(7,15) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(7,16) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(8,12) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(8,13) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(8,14) = 0.4472135954999579*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(8,15) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(8,17) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(9,12) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(9,13) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(9,14) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(9,15) = 0.4472135954999579*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(9,16) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(9,17) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(10,12) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(10,13) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(10,15) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(10,16) = 0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(11,12) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(11,14) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(11,15) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(11,17) = 0.31943828249997*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(6,6) = 1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(6,8) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(6,9) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(6,10) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(6,11) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(7,6) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(7,7) = 1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(7,8) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(7,9) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(7,10) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(8,6) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(8,7) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(8,8) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(8,9) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(8,11) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(9,6) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(9,7) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(9,8) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(9,9) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(9,10) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(9,11) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(10,6) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(10,7) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(10,9) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(10,10) = 0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(11,6) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(11,8) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(11,9) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(11,11) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(6,18) = 1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(6,19) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(6,20) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(6,21) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(6,22) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(6,23) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(7,18) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(7,19) = 1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(7,20) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(7,21) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(7,22) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(8,18) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(8,19) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(8,20) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(8,21) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(8,23) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(9,18) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(9,19) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(9,20) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(9,21) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(9,22) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(9,23) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(10,18) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(10,19) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(10,21) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(10,22) = 0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(11,18) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(11,20) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(11,21) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(11,23) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[6]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2rSelf[0]*mnuSelf+m2rOther[0]*mnuOther; 
  mnuM2sum[1] = m2rSelf[1]*mnuSelf+m2rOther[1]*mnuOther; 
  mnuM2sum[2] = m2rSelf[2]*mnuSelf+m2rOther[2]*mnuOther; 
  mnuM2sum[3] = m2rSelf[3]*mnuSelf+m2rOther[3]*mnuOther; 
  mnuM2sum[4] = m2rSelf[4]*mnuSelf+m2rOther[4]*mnuOther; 
  mnuM2sum[5] = m2rSelf[5]*mnuSelf+m2rOther[5]*mnuOther; 
 
  double m1Relax[6]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(12,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(12,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(12,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(13,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(13,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(13,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(14,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(14,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,1) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,2) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(15,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(16,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(17,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(17,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(12,6) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(12,7) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(12,8) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(12,9) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(12,10) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(12,11) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(13,6) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(13,7) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(13,8) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(13,9) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(13,10) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(14,6) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(14,7) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(14,8) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(14,9) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(14,11) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(15,6) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(15,7) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(15,8) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(15,9) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(15,10) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(15,11) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(16,6) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(16,7) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(16,9) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(16,10) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(17,6) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(17,8) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(17,9) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(17,11) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(12,12) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(12,13) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(12,14) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(12,15) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(12,16) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(12,17) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(13,12) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(13,13) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(13,14) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(13,15) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(13,16) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(14,12) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(14,13) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(14,14) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(14,15) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(14,17) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(15,12) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(15,13) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(15,14) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(15,15) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(15,16) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(15,17) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(16,12) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(16,13) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(16,15) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(16,16) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(17,12) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(17,14) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(17,15) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(17,17) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(12,18) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(12,19) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(12,20) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(12,21) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(12,22) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(12,23) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(13,18) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(13,19) = 0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(13,20) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(13,21) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(13,22) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(14,18) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(14,19) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(14,20) = 0.4472135954999579*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(14,21) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(14,23) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(15,18) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(15,19) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(15,20) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(15,21) = 0.4472135954999579*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(15,22) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(15,23) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(16,18) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(16,19) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(16,21) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(16,22) = 0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(17,18) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(17,20) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(17,21) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(17,23) = 0.31943828249997*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(18,0) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(18,1) = (-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(18,2) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,3) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,4) = (-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf)-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(18,5) = (-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(19,0) = (-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,1) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(19,2) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,3) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,4) = (-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf)-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,5) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(20,0) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,1) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,2) = (-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf)-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(20,3) = (-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,4) = (-0.25*m0rSelf[2]*uSelf[4]*mnuSelf)-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(20,5) = (-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,0) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,1) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,2) = (-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,3) = (-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf)-0.2*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.2*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.81*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(21,4) = (-0.2*m0rSelf[3]*uSelf[5]*mnuSelf)-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,5) = (-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,0) = (-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf)-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(22,1) = (-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf)-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,2) = (-0.25*m0rSelf[2]*uSelf[4]*mnuSelf)-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,3) = (-0.2*m0rSelf[3]*uSelf[5]*mnuSelf)-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,4) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(22,5) = (-0.25*m0rSelf[4]*uSelf[5]*mnuSelf)-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(23,0) = (-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(23,1) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,2) = (-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,3) = (-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,4) = (-0.25*m0rSelf[4]*uSelf[5]*mnuSelf)-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(23,5) = (-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf)-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(18,12) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(18,13) = 0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(18,14) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(18,15) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(18,16) = 0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(18,17) = 0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(19,12) = 0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(19,13) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(19,14) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(19,15) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(19,16) = 0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(19,17) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(20,12) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(20,13) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(20,14) = 0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(20,15) = 0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(20,16) = 0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(20,17) = 0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(21,12) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(21,13) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(21,14) = 0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(21,15) = 0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.2*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.2*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(21,16) = 0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(21,17) = 0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(22,12) = 0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(22,13) = 0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(22,14) = 0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(22,15) = 0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(22,16) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(22,17) = 0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(23,12) = 0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(23,13) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(23,14) = 0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(23,15) = 0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(23,16) = 0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(23,17) = 0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += (m1rOther[3]-1.0*m1rSelf[3])*betaGreenep1*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (m1rOther[4]-1.0*m1rSelf[4])*betaGreenep1*mnuSelf+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (m1rOther[5]-1.0*m1rSelf[5])*betaGreenep1*mnuSelf+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
 
  double ucMSelf[6]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.5*cMSelf[a0+5]*uSelf[a0+5]+0.5*cMSelf[a0+4]*uSelf[a0+4]+0.5*cMSelf[a0+3]*uSelf[a0+3]+0.5*cMSelf[a0+2]*uSelf[a0+2]+0.5*cMSelf[a0+1]*uSelf[a0+1]+0.5*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.4472135954999579*cMSelf[a0+1]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+1]*cMSelf[a0+4]+0.5*cMSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.4472135954999579*cMSelf[a0+2]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+2]*cMSelf[a0+5]+0.5*cMSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*cMSelf[a0+2]; 
    ucMSelf[3] += 0.4472135954999579*cMSelf[a0+3]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+3]*cMSelf[a0+5]+0.4472135954999579*cMSelf[a0+3]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+3]*cMSelf[a0+4]+0.5*cMSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*cMSelf[a0+3]+0.5*cMSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*cMSelf[a0+2]; 
    ucMSelf[4] += 0.31943828249997*cMSelf[a0+4]*uSelf[a0+4]+0.5*cMSelf[a0]*uSelf[a0+4]+0.5*uSelf[a0]*cMSelf[a0+4]+0.4472135954999579*cMSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*cMSelf[a0+1]*uSelf[a0+1]; 
    ucMSelf[5] += 0.31943828249997*cMSelf[a0+5]*uSelf[a0+5]+0.5*cMSelf[a0]*uSelf[a0+5]+0.5*uSelf[a0]*cMSelf[a0+5]+0.4472135954999579*cMSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*cMSelf[a0+2]*uSelf[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(18,6) = 0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(18,7) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(18,8) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(18,9) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(18,10) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(18,11) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(19,6) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(19,7) = 0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(19,8) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(19,9) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(19,10) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(20,6) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(20,7) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(20,8) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(20,9) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(20,11) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(21,6) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(21,7) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(21,8) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(21,9) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(21,10) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(21,11) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(22,6) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(22,7) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(22,9) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(22,10) = 0.31943828249997*ucMSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(23,6) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(23,8) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(23,9) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(23,11) = 0.31943828249997*ucMSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  double ucMOther[6]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.5*cMOther[a0+5]*uOther[a0+5]+0.5*cMOther[a0+4]*uOther[a0+4]+0.5*cMOther[a0+3]*uOther[a0+3]+0.5*cMOther[a0+2]*uOther[a0+2]+0.5*cMOther[a0+1]*uOther[a0+1]+0.5*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.4472135954999579*cMOther[a0+1]*uOther[a0+4]+0.4472135954999579*uOther[a0+1]*cMOther[a0+4]+0.5*cMOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+1]+0.5*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.4472135954999579*cMOther[a0+2]*uOther[a0+5]+0.4472135954999579*uOther[a0+2]*cMOther[a0+5]+0.5*cMOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+2]+0.5*uOther[a0]*cMOther[a0+2]; 
    ucMOther[3] += 0.4472135954999579*cMOther[a0+3]*uOther[a0+5]+0.4472135954999579*uOther[a0+3]*cMOther[a0+5]+0.4472135954999579*cMOther[a0+3]*uOther[a0+4]+0.4472135954999579*uOther[a0+3]*cMOther[a0+4]+0.5*cMOther[a0]*uOther[a0+3]+0.5*uOther[a0]*cMOther[a0+3]+0.5*cMOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*cMOther[a0+2]; 
    ucMOther[4] += 0.31943828249997*cMOther[a0+4]*uOther[a0+4]+0.5*cMOther[a0]*uOther[a0+4]+0.5*uOther[a0]*cMOther[a0+4]+0.4472135954999579*cMOther[a0+3]*uOther[a0+3]+0.4472135954999579*cMOther[a0+1]*uOther[a0+1]; 
    ucMOther[5] += 0.31943828249997*cMOther[a0+5]*uOther[a0+5]+0.5*cMOther[a0]*uOther[a0+5]+0.5*uOther[a0]*cMOther[a0+5]+0.4472135954999579*cMOther[a0+3]*uOther[a0+3]+0.4472135954999579*cMOther[a0+2]*uOther[a0+2]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(18,18) = (-0.5*ucMOther[0]*mnuOther)-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(18,19) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(18,20) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(18,21) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(18,22) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(18,23) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(19,18) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(19,19) = (-0.4472135954999579*ucMOther[4]*mnuOther)-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(19,20) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(19,21) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(19,22) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(20,18) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(20,19) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(20,20) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(20,21) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(20,23) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(21,18) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(21,19) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(21,20) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(21,21) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.4472135954999579*ucMOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(21,22) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(21,23) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(22,18) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(22,19) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(22,21) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(22,22) = (-0.31943828249997*ucMOther[4]*mnuOther)-0.9583148474999099*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(23,18) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(23,20) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(23,21) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(23,23) = (-0.31943828249997*ucMOther[5]*mnuOther)-0.9583148474999099*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[6]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinESelf[0] += 0.5*m1rSelf[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0+3]*uSelf[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+1]*m1rSelf[a0+4]+0.5*m1rSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+2]*m1rSelf[a0+5]+0.5*m1rSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]; 
    kinESelf[3] += 0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]; 
    kinESelf[4] += 0.31943828249997*m1rSelf[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+4]+0.5*uSelf[a0]*m1rSelf[a0+4]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+1]; 
    kinESelf[5] += 0.31943828249997*m1rSelf[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0]*uSelf[a0+5]+0.5*uSelf[a0]*m1rSelf[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+2]; 
  } 
 
  double kinEOther[6]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    kinEOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.5*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.4472135954999579*m1rOther[a0+1]*uOther[a0+4]+0.4472135954999579*uOther[a0+1]*m1rOther[a0+4]+0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.4472135954999579*m1rOther[a0+2]*uOther[a0+5]+0.4472135954999579*uOther[a0+2]*m1rOther[a0+5]+0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    kinEOther[3] += 0.4472135954999579*m1rOther[a0+3]*uOther[a0+5]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+4]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
    kinEOther[4] += 0.31943828249997*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+4]+0.5*uOther[a0]*m1rOther[a0+4]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+1]; 
    kinEOther[5] += 0.31943828249997*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rOther[a0]*uOther[a0+5]+0.5*uOther[a0]*m1rOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+2]; 
  } 
 
  double relKinE[6]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.5*m1rSelf[a0+5]*uSelf[a0+5]-0.5*m1rOther[a0+5]*uSelf[a0+5]-0.5*m1rSelf[a0+5]*uOther[a0+5]+0.5*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rSelf[a0+4]*uSelf[a0+4]-0.5*m1rOther[a0+4]*uSelf[a0+4]-0.5*m1rSelf[a0+4]*uOther[a0+4]+0.5*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rSelf[a0+3]*uSelf[a0+3]-0.5*m1rOther[a0+3]*uSelf[a0+3]-0.5*m1rSelf[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]-0.5*m1rOther[a0+2]*uSelf[a0+2]-0.5*m1rSelf[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]-0.5*m1rOther[a0+1]*uSelf[a0+1]-0.5*m1rSelf[a0+1]*uOther[a0+1]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]-0.5*m1rOther[a0]*uSelf[a0]-0.5*m1rSelf[a0]*uOther[a0]+0.5*m1rOther[a0]*uOther[a0]; 
    relKinE[1] += 0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+4]-0.4472135954999579*m1rOther[a0+1]*uSelf[a0+4]-0.4472135954999579*m1rSelf[a0+1]*uOther[a0+4]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+4]+0.4472135954999579*uSelf[a0+1]*m1rSelf[a0+4]-0.4472135954999579*uOther[a0+1]*m1rSelf[a0+4]-0.4472135954999579*uSelf[a0+1]*m1rOther[a0+4]+0.4472135954999579*uOther[a0+1]*m1rOther[a0+4]+0.5*m1rSelf[a0+2]*uSelf[a0+3]-0.5*m1rOther[a0+2]*uSelf[a0+3]-0.5*m1rSelf[a0+2]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]-0.5*uOther[a0+2]*m1rSelf[a0+3]-0.5*uSelf[a0+2]*m1rOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]-0.5*m1rOther[a0]*uSelf[a0+1]-0.5*m1rSelf[a0]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]-0.5*uOther[a0]*m1rSelf[a0+1]-0.5*uSelf[a0]*m1rOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    relKinE[2] += 0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+5]-0.4472135954999579*m1rOther[a0+2]*uSelf[a0+5]-0.4472135954999579*m1rSelf[a0+2]*uOther[a0+5]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+5]+0.4472135954999579*uSelf[a0+2]*m1rSelf[a0+5]-0.4472135954999579*uOther[a0+2]*m1rSelf[a0+5]-0.4472135954999579*uSelf[a0+2]*m1rOther[a0+5]+0.4472135954999579*uOther[a0+2]*m1rOther[a0+5]+0.5*m1rSelf[a0+1]*uSelf[a0+3]-0.5*m1rOther[a0+1]*uSelf[a0+3]-0.5*m1rSelf[a0+1]*uOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]-0.5*uOther[a0+1]*m1rSelf[a0+3]-0.5*uSelf[a0+1]*m1rOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]-0.5*m1rOther[a0]*uSelf[a0+2]-0.5*m1rSelf[a0]*uOther[a0+2]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]-0.5*uOther[a0]*m1rSelf[a0+2]-0.5*uSelf[a0]*m1rOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    relKinE[3] += 0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+5]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+5]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+5]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+5]-0.4472135954999579*uOther[a0+3]*m1rSelf[a0+5]-0.4472135954999579*uSelf[a0+3]*m1rOther[a0+5]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+4]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+4]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+4]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+4]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+4]-0.4472135954999579*uOther[a0+3]*m1rSelf[a0+4]-0.4472135954999579*uSelf[a0+3]*m1rOther[a0+4]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+3]-0.5*m1rOther[a0]*uSelf[a0+3]-0.5*m1rSelf[a0]*uOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]-0.5*uOther[a0]*m1rSelf[a0+3]-0.5*uSelf[a0]*m1rOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]-0.5*m1rOther[a0+1]*uSelf[a0+2]-0.5*m1rSelf[a0+1]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]-0.5*uOther[a0+1]*m1rSelf[a0+2]-0.5*uSelf[a0+1]*m1rOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
    relKinE[4] += 0.31943828249997*m1rSelf[a0+4]*uSelf[a0+4]-0.31943828249997*m1rOther[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+4]-0.5*m1rOther[a0]*uSelf[a0+4]-0.31943828249997*m1rSelf[a0+4]*uOther[a0+4]+0.31943828249997*m1rOther[a0+4]*uOther[a0+4]-0.5*m1rSelf[a0]*uOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+4]+0.5*uSelf[a0]*m1rSelf[a0+4]-0.5*uOther[a0]*m1rSelf[a0+4]-0.5*uSelf[a0]*m1rOther[a0+4]+0.5*uOther[a0]*m1rOther[a0+4]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+1]-0.4472135954999579*m1rOther[a0+1]*uSelf[a0+1]-0.4472135954999579*m1rSelf[a0+1]*uOther[a0+1]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+1]; 
    relKinE[5] += 0.31943828249997*m1rSelf[a0+5]*uSelf[a0+5]-0.31943828249997*m1rOther[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0]*uSelf[a0+5]-0.5*m1rOther[a0]*uSelf[a0+5]-0.31943828249997*m1rSelf[a0+5]*uOther[a0+5]+0.31943828249997*m1rOther[a0+5]*uOther[a0+5]-0.5*m1rSelf[a0]*uOther[a0+5]+0.5*m1rOther[a0]*uOther[a0+5]+0.5*uSelf[a0]*m1rSelf[a0+5]-0.5*uOther[a0]*m1rSelf[a0+5]-0.5*uSelf[a0]*m1rOther[a0+5]+0.5*uOther[a0]*m1rOther[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+2]-0.4472135954999579*m1rOther[a0+2]*uSelf[a0+2]-0.4472135954999579*m1rSelf[a0+2]*uOther[a0+2]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+2]; 
  } 
 
  double m2Relax[6]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(0.5*relKinE[0]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[0]*mSelf)/(mSelf+mOther)+(kinESelf[0]*mSelf)/(mSelf+mOther)+(0.5*relKinE[0]*mOther)/(mSelf+mOther)+(m2rOther[0]*mOther)/(mSelf+mOther)-(1.0*kinEOther[0]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2rOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(0.5*relKinE[1]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[1]*mSelf)/(mSelf+mOther)+(kinESelf[1]*mSelf)/(mSelf+mOther)+(0.5*relKinE[1]*mOther)/(mSelf+mOther)+(m2rOther[1]*mOther)/(mSelf+mOther)-(1.0*kinEOther[1]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2rOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(0.5*relKinE[2]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[2]*mSelf)/(mSelf+mOther)+(kinESelf[2]*mSelf)/(mSelf+mOther)+(0.5*relKinE[2]*mOther)/(mSelf+mOther)+(m2rOther[2]*mOther)/(mSelf+mOther)-(1.0*kinEOther[2]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2rOther[2])*mnuOther; 
  m2Relax[3] = betaGreenep1*((-(0.5*relKinE[3]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[3]*mSelf)/(mSelf+mOther)+(kinESelf[3]*mSelf)/(mSelf+mOther)+(0.5*relKinE[3]*mOther)/(mSelf+mOther)+(m2rOther[3]*mOther)/(mSelf+mOther)-(1.0*kinEOther[3]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[3]-1.0*kinESelf[3])*mnuSelf+(kinEOther[3]-1.0*m2rOther[3])*mnuOther; 
  m2Relax[4] = betaGreenep1*((-(0.5*relKinE[4]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[4]*mSelf)/(mSelf+mOther)+(kinESelf[4]*mSelf)/(mSelf+mOther)+(0.5*relKinE[4]*mOther)/(mSelf+mOther)+(m2rOther[4]*mOther)/(mSelf+mOther)-(1.0*kinEOther[4]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[4]-1.0*kinESelf[4])*mnuSelf+(kinEOther[4]-1.0*m2rOther[4])*mnuOther; 
  m2Relax[5] = betaGreenep1*((-(0.5*relKinE[5]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[5]*mSelf)/(mSelf+mOther)+(kinESelf[5]*mSelf)/(mSelf+mOther)+(0.5*relKinE[5]*mOther)/(mSelf+mOther)+(m2rOther[5]*mOther)/(mSelf+mOther)-(1.0*kinEOther[5]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[5]-1.0*kinESelf[5])*mnuSelf+(kinEOther[5]-1.0*m2rOther[5])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,6,1) = data->u_S.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,6,1) = data->u_S.segment<6>(6); 
 
  Eigen::Map<VectorXd>(uCrossOther,6,1) = data->u_S.segment<6>(12); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,6,1) = data->u_S.segment<6>(18); 
 
} 
 
