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
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(9,3) = 0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(9,4) = 0.25*cMSelf[0]*uSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(9,5) = 0.25*cMSelf[0]*uSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(10,3) = 0.25*cMSelf[0]*uSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(10,4) = 0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(11,3) = 0.25*cMSelf[0]*uSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(11,5) = 0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(9,9) = (-0.25*cMOther[2]*uOther[2]*mnuOther)-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(9,10) = (-0.25*cMOther[0]*uOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(9,11) = (-0.25*cMOther[0]*uOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(10,9) = (-0.25*cMOther[0]*uOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(10,10) = (-0.25*cMOther[2]*uOther[2]*mnuOther)-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(11,9) = (-0.25*cMOther[0]*uOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(11,11) = (-0.25*cMOther[2]*uOther[2]*mnuOther)-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
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
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(18,6) = 0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[3]*uSelf[3]*mnuSelf+0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(18,7) = 0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.25*cMSelf[2]*uSelf[3]*mnuSelf+0.25*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(18,8) = 0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.25*cMSelf[1]*uSelf[3]*mnuSelf+0.25*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(18,9) = 0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(18,10) = 0.159719141249985*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf+0.25*uSelf[0]*cMSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(18,11) = 0.159719141249985*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[0]*uSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf+0.25*uSelf[0]*cMSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(19,6) = 0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.25*cMSelf[2]*uSelf[3]*mnuSelf+0.25*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(19,7) = 0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.3928571428571428*cMSelf[4]*uSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.45*cMSelf[3]*uSelf[3]*mnuSelf+0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.45*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(19,8) = 0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(19,9) = 0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.25*cMSelf[1]*uSelf[3]*mnuSelf+0.25*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(19,10) = 0.2*cMSelf[1]*uSelf[4]*mnuSelf+0.2*uSelf[1]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(20,6) = 0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.25*cMSelf[1]*uSelf[3]*mnuSelf+0.25*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(20,7) = 0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(20,8) = 0.3928571428571428*cMSelf[5]*uSelf[5]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.45*cMSelf[3]*uSelf[3]*mnuSelf+0.45*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(20,9) = 0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.25*cMSelf[2]*uSelf[3]*mnuSelf+0.25*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(20,11) = 0.2*cMSelf[2]*uSelf[5]*mnuSelf+0.2*uSelf[2]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(21,6) = 0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(21,7) = 0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.25*cMSelf[1]*uSelf[3]*mnuSelf+0.25*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(21,8) = 0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.25*cMSelf[2]*uSelf[3]*mnuSelf+0.25*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(21,9) = 0.3928571428571428*cMSelf[5]*uSelf[5]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.3928571428571428*cMSelf[4]*uSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.65*cMSelf[3]*uSelf[3]*mnuSelf+0.45*cMSelf[2]*uSelf[2]*mnuSelf+0.45*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(21,10) = 0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(21,11) = 0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(22,6) = 0.159719141249985*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf+0.25*uSelf[0]*cMSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(22,7) = 0.2*cMSelf[1]*uSelf[4]*mnuSelf+0.2*uSelf[1]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(22,9) = 0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(22,10) = 0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.3520408163265306*cMSelf[4]*uSelf[4]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.3928571428571428*cMSelf[3]*uSelf[3]*mnuSelf+0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.3928571428571428*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(23,6) = 0.159719141249985*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[0]*uSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf+0.25*uSelf[0]*cMSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(23,8) = 0.2*cMSelf[2]*uSelf[5]*mnuSelf+0.2*uSelf[2]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(23,9) = 0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(23,11) = 0.3520408163265306*cMSelf[5]*uSelf[5]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.3928571428571428*cMSelf[3]*uSelf[3]*mnuSelf+0.3928571428571428*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(18,18) = (-0.25*cMOther[5]*uOther[5]*mnuOther)-0.25*cMOther[4]*uOther[4]*mnuOther-0.25*cMOther[3]*uOther[3]*mnuOther-0.25*cMOther[2]*uOther[2]*mnuOther-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(18,19) = (-0.223606797749979*cMOther[1]*uOther[4]*mnuOther)-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.25*cMOther[2]*uOther[3]*mnuOther-0.25*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(18,20) = (-0.223606797749979*cMOther[2]*uOther[5]*mnuOther)-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.25*cMOther[1]*uOther[3]*mnuOther-0.25*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(18,21) = (-0.223606797749979*cMOther[3]*uOther[5]*mnuOther)-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(18,22) = (-0.159719141249985*cMOther[4]*uOther[4]*mnuOther)-0.25*cMOther[0]*uOther[4]*mnuOther-1.5*m0rOther[4]*mnuOther-0.25*uOther[0]*cMOther[4]*mnuOther+0.5*cEOther[4]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(18,23) = (-0.159719141249985*cMOther[5]*uOther[5]*mnuOther)-0.25*cMOther[0]*uOther[5]*mnuOther-1.5*m0rOther[5]*mnuOther-0.25*uOther[0]*cMOther[5]*mnuOther+0.5*cEOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(19,18) = (-0.223606797749979*cMOther[1]*uOther[4]*mnuOther)-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.25*cMOther[2]*uOther[3]*mnuOther-0.25*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(19,19) = (-0.25*cMOther[5]*uOther[5]*mnuOther)-0.3928571428571428*cMOther[4]*uOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther-0.223606797749979*uOther[0]*cMOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.45*cMOther[3]*uOther[3]*mnuOther-0.25*cMOther[2]*uOther[2]*mnuOther-0.45*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(19,20) = (-0.223606797749979*cMOther[3]*uOther[5]*mnuOther)-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(19,21) = (-0.223606797749979*cMOther[2]*uOther[5]*mnuOther)-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.25*cMOther[1]*uOther[3]*mnuOther-0.25*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(19,22) = (-0.2*cMOther[1]*uOther[4]*mnuOther)-0.2*uOther[1]*cMOther[4]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther-0.223606797749979*uOther[0]*cMOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(20,18) = (-0.223606797749979*cMOther[2]*uOther[5]*mnuOther)-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.25*cMOther[1]*uOther[3]*mnuOther-0.25*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(20,19) = (-0.223606797749979*cMOther[3]*uOther[5]*mnuOther)-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(20,20) = (-0.3928571428571428*cMOther[5]*uOther[5]*mnuOther)-0.223606797749979*cMOther[0]*uOther[5]*mnuOther-1.341640786499874*m0rOther[5]*mnuOther-0.223606797749979*uOther[0]*cMOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.25*cMOther[4]*uOther[4]*mnuOther-0.45*cMOther[3]*uOther[3]*mnuOther-0.45*cMOther[2]*uOther[2]*mnuOther-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(20,21) = (-0.223606797749979*cMOther[1]*uOther[4]*mnuOther)-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.25*cMOther[2]*uOther[3]*mnuOther-0.25*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(20,23) = (-0.2*cMOther[2]*uOther[5]*mnuOther)-0.2*uOther[2]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther-0.223606797749979*uOther[0]*cMOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(21,18) = (-0.223606797749979*cMOther[3]*uOther[5]*mnuOther)-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(21,19) = (-0.223606797749979*cMOther[2]*uOther[5]*mnuOther)-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.25*cMOther[1]*uOther[3]*mnuOther-0.25*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(21,20) = (-0.223606797749979*cMOther[1]*uOther[4]*mnuOther)-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.25*cMOther[2]*uOther[3]*mnuOther-0.25*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(21,21) = (-0.3928571428571428*cMOther[5]*uOther[5]*mnuOther)-0.223606797749979*cMOther[0]*uOther[5]*mnuOther-1.341640786499874*m0rOther[5]*mnuOther-0.223606797749979*uOther[0]*cMOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.3928571428571428*cMOther[4]*uOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther-0.223606797749979*uOther[0]*cMOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.65*cMOther[3]*uOther[3]*mnuOther-0.45*cMOther[2]*uOther[2]*mnuOther-0.45*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(21,22) = (-0.2*cMOther[3]*uOther[5]*mnuOther)-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(21,23) = (-0.2*cMOther[3]*uOther[5]*mnuOther)-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(22,18) = (-0.159719141249985*cMOther[4]*uOther[4]*mnuOther)-0.25*cMOther[0]*uOther[4]*mnuOther-1.5*m0rOther[4]*mnuOther-0.25*uOther[0]*cMOther[4]*mnuOther+0.5*cEOther[4]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(22,19) = (-0.2*cMOther[1]*uOther[4]*mnuOther)-0.2*uOther[1]*cMOther[4]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther-0.223606797749979*uOther[0]*cMOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(22,21) = (-0.2*cMOther[3]*uOther[5]*mnuOther)-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(22,22) = (-0.25*cMOther[5]*uOther[5]*mnuOther)-0.3520408163265306*cMOther[4]*uOther[4]*mnuOther-0.159719141249985*cMOther[0]*uOther[4]*mnuOther-0.9583148474999099*m0rOther[4]*mnuOther-0.159719141249985*uOther[0]*cMOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.3928571428571428*cMOther[3]*uOther[3]*mnuOther-0.25*cMOther[2]*uOther[2]*mnuOther-0.3928571428571428*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(23,18) = (-0.159719141249985*cMOther[5]*uOther[5]*mnuOther)-0.25*cMOther[0]*uOther[5]*mnuOther-1.5*m0rOther[5]*mnuOther-0.25*uOther[0]*cMOther[5]*mnuOther+0.5*cEOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(23,20) = (-0.2*cMOther[2]*uOther[5]*mnuOther)-0.2*uOther[2]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther-0.223606797749979*uOther[0]*cMOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(23,21) = (-0.2*cMOther[3]*uOther[5]*mnuOther)-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(23,23) = (-0.3520408163265306*cMOther[5]*uOther[5]*mnuOther)-0.159719141249985*cMOther[0]*uOther[5]*mnuOther-0.9583148474999099*m0rOther[5]*mnuOther-0.159719141249985*uOther[0]*cMOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.25*cMOther[4]*uOther[4]*mnuOther-0.3928571428571428*cMOther[3]*uOther[3]*mnuOther-0.3928571428571428*cMOther[2]*uOther[2]*mnuOther-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
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
 
void GkCrossPrimMoments2x2vMax_P3(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-1.322875655532295*m0Self[9])-1.322875655532295*m0Self[8]-1.936491673103709*m0Self[7]-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0Self[9])-1.322875655532295*m0Self[8]-1.936491673103709*m0Self[7]-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0Self[9])+1.322875655532295*m0Self[8]+1.936491673103709*m0Self[7]-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0Self[9])+1.322875655532295*m0Self[8]+1.936491673103709*m0Self[7]-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[10]; 
  double m1rSelf[10]; 
  double m2rSelf[10]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m0rSelf[3] = 0.0; 
    m0rSelf[4] = 0.0; 
    m0rSelf[5] = 0.0; 
    m0rSelf[6] = 0.0; 
    m0rSelf[7] = 0.0; 
    m0rSelf[8] = 0.0; 
    m0rSelf[9] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    m1rSelf[2] = 0.0; 
    m1rSelf[3] = 0.0; 
    m1rSelf[4] = 0.0; 
    m1rSelf[5] = 0.0; 
    m1rSelf[6] = 0.0; 
    m1rSelf[7] = 0.0; 
    m1rSelf[8] = 0.0; 
    m1rSelf[9] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    m2rSelf[3] = 0.0; 
    m2rSelf[4] = 0.0; 
    m2rSelf[5] = 0.0; 
    m2rSelf[6] = 0.0; 
    m2rSelf[7] = 0.0; 
    m2rSelf[8] = 0.0; 
    m2rSelf[9] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
    m0rSelf[3] = m0Self[3]; 
    m0rSelf[4] = m0Self[4]; 
    m0rSelf[5] = m0Self[5]; 
    m0rSelf[6] = m0Self[6]; 
    m0rSelf[7] = m0Self[7]; 
    m0rSelf[8] = m0Self[8]; 
    m0rSelf[9] = m0Self[9]; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = m1Self[1]; 
    m1rSelf[2] = m1Self[2]; 
    m1rSelf[3] = m1Self[3]; 
    m1rSelf[4] = m1Self[4]; 
    m1rSelf[5] = m1Self[5]; 
    m1rSelf[6] = m1Self[6]; 
    m1rSelf[7] = m1Self[7]; 
    m1rSelf[8] = m1Self[8]; 
    m1rSelf[9] = m1Self[9]; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = m2Self[1]; 
    m2rSelf[2] = m2Self[2]; 
    m2rSelf[3] = m2Self[3]; 
    m2rSelf[4] = m2Self[4]; 
    m2rSelf[5] = m2Self[5]; 
    m2rSelf[6] = m2Self[6]; 
    m2rSelf[7] = m2Self[7]; 
    m2rSelf[8] = m2Self[8]; 
    m2rSelf[9] = m2Self[9]; 
  } 
 
  if ((-1.322875655532295*m0Other[9])-1.322875655532295*m0Other[8]-1.936491673103709*m0Other[7]-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0Other[9])-1.322875655532295*m0Other[8]-1.936491673103709*m0Other[7]-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0Other[9])+1.322875655532295*m0Other[8]+1.936491673103709*m0Other[7]-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.322875655532295*m0Other[9])+1.322875655532295*m0Other[8]+1.936491673103709*m0Other[7]-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[10]; 
  double m1rOther[10]; 
  double m2rOther[10]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m0rOther[3] = 0.0; 
    m0rOther[4] = 0.0; 
    m0rOther[5] = 0.0; 
    m0rOther[6] = 0.0; 
    m0rOther[7] = 0.0; 
    m0rOther[8] = 0.0; 
    m0rOther[9] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    m1rOther[2] = 0.0; 
    m1rOther[3] = 0.0; 
    m1rOther[4] = 0.0; 
    m1rOther[5] = 0.0; 
    m1rOther[6] = 0.0; 
    m1rOther[7] = 0.0; 
    m1rOther[8] = 0.0; 
    m1rOther[9] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    m2rOther[3] = 0.0; 
    m2rOther[4] = 0.0; 
    m2rOther[5] = 0.0; 
    m2rOther[6] = 0.0; 
    m2rOther[7] = 0.0; 
    m2rOther[8] = 0.0; 
    m2rOther[9] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
    m0rOther[3] = m0Other[3]; 
    m0rOther[4] = m0Other[4]; 
    m0rOther[5] = m0Other[5]; 
    m0rOther[6] = m0Other[6]; 
    m0rOther[7] = m0Other[7]; 
    m0rOther[8] = m0Other[8]; 
    m0rOther[9] = m0Other[9]; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = m1Other[1]; 
    m1rOther[2] = m1Other[2]; 
    m1rOther[3] = m1Other[3]; 
    m1rOther[4] = m1Other[4]; 
    m1rOther[5] = m1Other[5]; 
    m1rOther[6] = m1Other[6]; 
    m1rOther[7] = m1Other[7]; 
    m1rOther[8] = m1Other[8]; 
    m1rOther[9] = m1Other[9]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
    m2rOther[3] = m2Other[3]; 
    m2rOther[4] = m2Other[4]; 
    m2rOther[5] = m2Other[5]; 
    m2rOther[6] = m2Other[6]; 
    m2rOther[7] = m2Other[7]; 
    m2rOther[8] = m2Other[8]; 
    m2rOther[9] = m2Other[9]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(40,40); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
  double mnuM1sum[10]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<10; vd++) 
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
  data->AEM_S(0,6) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(0,7) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(0,8) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(0,9) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(1,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,3) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(1,4) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,5) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(1,6) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,7) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(1,8) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(2,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,4) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(2,5) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,6) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(2,7) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(2,9) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(3,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,1) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,2) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(3,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,6) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,7) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,8) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(3,9) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(4,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(4,1) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,2) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(4,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(4,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,6) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,7) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(4,8) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(5,1) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(5,2) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(5,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(5,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(5,6) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(5,7) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,9) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,0) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(6,1) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(6,2) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(6,3) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,4) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,5) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(6,6) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(6,8) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,0) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(7,1) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(7,2) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,3) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,4) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(7,5) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,6) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,7) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(7,9) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(8,0) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(8,1) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(8,3) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(8,4) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,6) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(8,8) = 0.2981423969999719*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(9,0) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(9,2) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(9,3) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(9,5) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(9,7) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,9) = 0.2981423969999719*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,10) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,11) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,12) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,13) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,14) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,15) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(0,16) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(0,17) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(0,18) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(0,19) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(1,10) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,11) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,12) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,13) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,14) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,15) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(1,16) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,17) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(1,18) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(2,10) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,11) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,12) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,13) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,14) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(2,15) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,16) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(2,17) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,19) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(3,10) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,11) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,12) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,13) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,14) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,15) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,16) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,17) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,18) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(3,19) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(4,10) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,11) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,12) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(4,13) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,14) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(4,16) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(4,17) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(4,18) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(5,10) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,11) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,12) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(5,13) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,15) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,16) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(5,17) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(5,19) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(6,10) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,11) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(6,12) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(6,13) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,14) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(6,15) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,16) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(6,17) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(6,18) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,10) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,11) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(7,12) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,13) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(7,14) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,15) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(7,16) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,17) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(7,19) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(8,10) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,11) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(8,13) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(8,14) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(8,16) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(8,18) = (-0.2981423969999719*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(9,10) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,12) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(9,13) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(9,15) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(9,17) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(9,19) = (-0.2981423969999719*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,20) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,21) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,22) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,23) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(0,24) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(0,25) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(0,26) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(0,27) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(0,28) = 0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(0,29) = 0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(1,20) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,21) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,22) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,23) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,24) = 0.4391550328268398*m0rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(1,25) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(1,26) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(1,27) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(1,28) = 0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(2,20) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,21) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,22) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,23) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(2,24) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(2,25) = 0.4391550328268398*m0rOther[9]*mnuOther+0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(2,26) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(2,27) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(2,29) = 0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(3,20) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,21) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,22) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,23) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,24) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,25) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,26) = 0.4391550328268399*m0rOther[8]*mnuOther+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(3,27) = 0.4391550328268399*m0rOther[9]*mnuOther+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(3,28) = 0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(3,29) = 0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(4,20) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(4,21) = 0.4391550328268398*m0rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(4,22) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(4,23) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,24) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(4,26) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(4,27) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(4,28) = 0.2981423969999719*m0rOther[8]*mnuOther+0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(5,20) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(5,21) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(5,22) = 0.4391550328268398*m0rOther[9]*mnuOther+0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(5,23) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(5,25) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,26) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(5,27) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(5,29) = 0.2981423969999719*m0rOther[9]*mnuOther+0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(6,20) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(6,21) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(6,22) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(6,23) = 0.4391550328268399*m0rOther[8]*mnuOther+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(6,24) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(6,25) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(6,26) = 0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,27) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(6,28) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(7,20) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(7,21) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(7,22) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(7,23) = 0.4391550328268399*m0rOther[9]*mnuOther+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(7,24) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(7,25) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(7,26) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(7,27) = 0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(7,29) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(8,20) = 0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(8,21) = 0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(8,23) = 0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(8,24) = 0.2981423969999719*m0rOther[8]*mnuOther+0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(8,26) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(8,28) = 0.2981423969999719*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,20) = 0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(9,22) = 0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(9,23) = 0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(9,25) = 0.2981423969999719*m0rOther[9]*mnuOther+0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(9,27) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(9,29) = 0.2981423969999719*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,30) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,31) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,32) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,33) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(0,34) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(0,35) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(0,36) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(0,37) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(0,38) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(0,39) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(1,30) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,31) = (-0.4472135954999579*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,32) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,33) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,34) = (-0.4391550328268398*cMOther[8]*mnuOther)-0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(1,35) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(1,36) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(1,37) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(1,38) = -0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(2,30) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,31) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,32) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,33) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(2,34) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(2,35) = (-0.4391550328268398*cMOther[9]*mnuOther)-0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(2,36) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(2,37) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(2,39) = -0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(3,30) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,31) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,32) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,33) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(3,34) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,35) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,36) = (-0.4391550328268399*cMOther[8]*mnuOther)-0.4*cMOther[7]*mnuOther-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(3,37) = (-0.4391550328268399*cMOther[9]*mnuOther)-0.4*cMOther[6]*mnuOther-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(3,38) = -0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(3,39) = -0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(4,30) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,31) = (-0.4391550328268398*cMOther[8]*mnuOther)-0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(4,32) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(4,33) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,34) = (-0.31943828249997*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(4,36) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(4,37) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(4,38) = (-0.2981423969999719*cMOther[8]*mnuOther)-0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(5,30) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,31) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(5,32) = (-0.4391550328268398*cMOther[9]*mnuOther)-0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(5,33) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(5,35) = (-0.31943828249997*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(5,36) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(5,37) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(5,39) = (-0.2981423969999719*cMOther[9]*mnuOther)-0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(6,30) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,31) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(6,32) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(6,33) = (-0.4391550328268399*cMOther[8]*mnuOther)-0.4*cMOther[7]*mnuOther-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(6,34) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(6,35) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(6,36) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.31943828249997*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(6,37) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(6,38) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(7,30) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,31) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(7,32) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(7,33) = (-0.4391550328268399*cMOther[9]*mnuOther)-0.4*cMOther[6]*mnuOther-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(7,34) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(7,35) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(7,36) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(7,37) = (-0.31943828249997*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(7,39) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(8,30) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(8,31) = -0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(8,33) = -0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(8,34) = (-0.2981423969999719*cMOther[8]*mnuOther)-0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(8,36) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(8,38) = (-0.2981423969999719*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(9,30) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(9,32) = -0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(9,33) = -0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(9,35) = (-0.2981423969999719*cMOther[9]*mnuOther)-0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(9,37) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(9,39) = (-0.2981423969999719*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(10,0) = 0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(10,1) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(10,2) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(10,3) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(10,4) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(10,5) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(10,6) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(10,7) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(10,8) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(10,9) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(11,0) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(11,1) = 0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(11,2) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(11,3) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(11,4) = 0.4391550328268398*m1rSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(11,5) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(11,6) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(11,7) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(11,8) = 0.4391550328268398*m1rSelf[4]*mnuSelf; 
  data->AEM_S(12,0) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(12,1) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(12,2) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(12,3) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(12,4) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(12,5) = 0.4391550328268398*m1rSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(12,6) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(12,7) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(12,9) = 0.4391550328268398*m1rSelf[5]*mnuSelf; 
  data->AEM_S(13,0) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(13,1) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(13,2) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(13,3) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(13,4) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(13,5) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(13,6) = 0.4391550328268399*m1rSelf[8]*mnuSelf+0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(13,7) = 0.4391550328268399*m1rSelf[9]*mnuSelf+0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(13,8) = 0.4391550328268399*m1rSelf[6]*mnuSelf; 
  data->AEM_S(13,9) = 0.4391550328268399*m1rSelf[7]*mnuSelf; 
  data->AEM_S(14,0) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(14,1) = 0.4391550328268398*m1rSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(14,2) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(14,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(14,4) = 0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(14,6) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(14,7) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(14,8) = 0.2981423969999719*m1rSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf; 
  data->AEM_S(15,0) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(15,1) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(15,2) = 0.4391550328268398*m1rSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(15,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(15,5) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(15,6) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(15,7) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(15,9) = 0.2981423969999719*m1rSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf; 
  data->AEM_S(16,0) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(16,1) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(16,2) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(16,3) = 0.4391550328268399*m1rSelf[8]*mnuSelf+0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(16,4) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(16,5) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(16,6) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(16,7) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(16,8) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(17,0) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(17,1) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(17,2) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(17,3) = 0.4391550328268399*m1rSelf[9]*mnuSelf+0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(17,4) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(17,5) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(17,6) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(17,7) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(17,9) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(18,0) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(18,1) = 0.4391550328268398*m1rSelf[4]*mnuSelf; 
  data->AEM_S(18,3) = 0.4391550328268399*m1rSelf[6]*mnuSelf; 
  data->AEM_S(18,4) = 0.2981423969999719*m1rSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf; 
  data->AEM_S(18,6) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(18,8) = 0.2981423969999719*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(19,0) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(19,2) = 0.4391550328268398*m1rSelf[5]*mnuSelf; 
  data->AEM_S(19,3) = 0.4391550328268399*m1rSelf[7]*mnuSelf; 
  data->AEM_S(19,5) = 0.2981423969999719*m1rSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf; 
  data->AEM_S(19,7) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(19,9) = 0.2981423969999719*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(10,20) = 0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(10,21) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(10,22) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(10,23) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(10,24) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(10,25) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(10,26) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(10,27) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(10,28) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(10,29) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(11,20) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(11,21) = 0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(11,22) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(11,23) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(11,24) = 0.4391550328268398*m1rOther[8]*mnuOther+0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(11,25) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(11,26) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(11,27) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(11,28) = 0.4391550328268398*m1rOther[4]*mnuOther; 
  data->AEM_S(12,20) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(12,21) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(12,22) = 0.4472135954999579*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(12,23) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(12,24) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(12,25) = 0.4391550328268398*m1rOther[9]*mnuOther+0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(12,26) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(12,27) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(12,29) = 0.4391550328268398*m1rOther[5]*mnuOther; 
  data->AEM_S(13,20) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(13,21) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(13,22) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(13,23) = 0.4472135954999579*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(13,24) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(13,25) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(13,26) = 0.4391550328268399*m1rOther[8]*mnuOther+0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(13,27) = 0.4391550328268399*m1rOther[9]*mnuOther+0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(13,28) = 0.4391550328268399*m1rOther[6]*mnuOther; 
  data->AEM_S(13,29) = 0.4391550328268399*m1rOther[7]*mnuOther; 
  data->AEM_S(14,20) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(14,21) = 0.4391550328268398*m1rOther[8]*mnuOther+0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(14,22) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(14,23) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(14,24) = 0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(14,26) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(14,27) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(14,28) = 0.2981423969999719*m1rOther[8]*mnuOther+0.4391550328268398*m1rOther[1]*mnuOther; 
  data->AEM_S(15,20) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(15,21) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(15,22) = 0.4391550328268398*m1rOther[9]*mnuOther+0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(15,23) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(15,25) = 0.31943828249997*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(15,26) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(15,27) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(15,29) = 0.2981423969999719*m1rOther[9]*mnuOther+0.4391550328268398*m1rOther[2]*mnuOther; 
  data->AEM_S(16,20) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(16,21) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(16,22) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(16,23) = 0.4391550328268399*m1rOther[8]*mnuOther+0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(16,24) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(16,25) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(16,26) = 0.4472135954999579*m1rOther[5]*mnuOther+0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(16,27) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(16,28) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(17,20) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(17,21) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(17,22) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(17,23) = 0.4391550328268399*m1rOther[9]*mnuOther+0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(17,24) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(17,25) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(17,26) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(17,27) = 0.31943828249997*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(17,29) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(18,20) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(18,21) = 0.4391550328268398*m1rOther[4]*mnuOther; 
  data->AEM_S(18,23) = 0.4391550328268399*m1rOther[6]*mnuOther; 
  data->AEM_S(18,24) = 0.2981423969999719*m1rOther[8]*mnuOther+0.4391550328268398*m1rOther[1]*mnuOther; 
  data->AEM_S(18,26) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(18,28) = 0.2981423969999719*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(19,20) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(19,22) = 0.4391550328268398*m1rOther[5]*mnuOther; 
  data->AEM_S(19,23) = 0.4391550328268399*m1rOther[7]*mnuOther; 
  data->AEM_S(19,25) = 0.2981423969999719*m1rOther[9]*mnuOther+0.4391550328268398*m1rOther[2]*mnuOther; 
  data->AEM_S(19,27) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(19,29) = 0.2981423969999719*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
  mnuM1sum[9] += m1rSelf[9]*mnuSelf+m1rOther[9]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(10,10) = 1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(10,11) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(10,12) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(10,13) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(10,14) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(10,15) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(10,16) = 1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(10,17) = 1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(10,18) = 1.5*m0rSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf; 
  data->AEM_S(10,19) = 1.5*m0rSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf; 
  data->AEM_S(11,10) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(11,11) = 1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(11,12) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(11,13) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(11,14) = 1.317465098480519*m0rSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(11,15) = 1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(11,16) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(11,17) = 1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(11,18) = 1.317465098480519*m0rSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf; 
  data->AEM_S(12,10) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(12,11) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(12,12) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(12,14) = 1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(12,15) = 1.317465098480519*m0rSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(12,16) = 1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(12,17) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(12,19) = 1.317465098480519*m0rSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf; 
  data->AEM_S(13,10) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(13,11) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(13,12) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(13,15) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(13,16) = 1.31746509848052*m0rSelf[8]*mnuSelf-0.4391550328268399*cESelf[8]*mnuSelf+1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(13,17) = 1.31746509848052*m0rSelf[9]*mnuSelf-0.4391550328268399*cESelf[9]*mnuSelf+1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(13,18) = 1.31746509848052*m0rSelf[6]*mnuSelf-0.4391550328268399*cESelf[6]*mnuSelf; 
  data->AEM_S(13,19) = 1.31746509848052*m0rSelf[7]*mnuSelf-0.4391550328268399*cESelf[7]*mnuSelf; 
  data->AEM_S(14,10) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(14,11) = 1.317465098480519*m0rSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(14,12) = 1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(14,13) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(14,14) = 0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(14,16) = 0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(14,17) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(14,18) = 0.8944271909999159*m0rSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+1.317465098480519*m0rSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(15,10) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(15,11) = 1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(15,12) = 1.317465098480519*m0rSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(15,13) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(15,15) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(15,16) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(15,17) = 0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(15,19) = 0.8944271909999159*m0rSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+1.317465098480519*m0rSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(16,10) = 1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(16,11) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(16,12) = 1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(16,13) = 1.31746509848052*m0rSelf[8]*mnuSelf-0.4391550328268399*cESelf[8]*mnuSelf+1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(16,14) = 0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(16,15) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(16,16) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(16,17) = 1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(16,18) = 1.31746509848052*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(17,10) = 1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(17,11) = 1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(17,12) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(17,13) = 1.31746509848052*m0rSelf[9]*mnuSelf-0.4391550328268399*cESelf[9]*mnuSelf+1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(17,14) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(17,15) = 0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(17,16) = 1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(17,17) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(17,19) = 1.31746509848052*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(18,10) = 1.5*m0rSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf; 
  data->AEM_S(18,11) = 1.317465098480519*m0rSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf; 
  data->AEM_S(18,13) = 1.31746509848052*m0rSelf[6]*mnuSelf-0.4391550328268399*cESelf[6]*mnuSelf; 
  data->AEM_S(18,14) = 0.8944271909999159*m0rSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+1.317465098480519*m0rSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(18,16) = 1.31746509848052*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(18,18) = 0.8944271909999159*m0rSelf[4]*mnuSelf-0.2981423969999719*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(19,10) = 1.5*m0rSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf; 
  data->AEM_S(19,12) = 1.317465098480519*m0rSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf; 
  data->AEM_S(19,13) = 1.31746509848052*m0rSelf[7]*mnuSelf-0.4391550328268399*cESelf[7]*mnuSelf; 
  data->AEM_S(19,15) = 0.8944271909999159*m0rSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+1.317465098480519*m0rSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(19,17) = 1.31746509848052*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(19,19) = 0.8944271909999159*m0rSelf[5]*mnuSelf-0.2981423969999719*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(10,30) = 1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(10,31) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(10,32) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(10,33) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(10,34) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(10,35) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(10,36) = 1.5*m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(10,37) = 1.5*m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(10,38) = 1.5*m0rOther[8]*mnuOther-0.5*cEOther[8]*mnuOther; 
  data->AEM_S(10,39) = 1.5*m0rOther[9]*mnuOther-0.5*cEOther[9]*mnuOther; 
  data->AEM_S(11,30) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(11,31) = 1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(11,32) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(11,33) = 1.341640786499874*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(11,34) = 1.317465098480519*m0rOther[8]*mnuOther-0.4391550328268398*cEOther[8]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(11,35) = 1.5*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(11,36) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(11,37) = 1.5*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(11,38) = 1.317465098480519*m0rOther[4]*mnuOther-0.4391550328268398*cEOther[4]*mnuOther; 
  data->AEM_S(12,30) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(12,31) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(12,32) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(12,33) = 1.341640786499874*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(12,34) = 1.5*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(12,35) = 1.317465098480519*m0rOther[9]*mnuOther-0.4391550328268398*cEOther[9]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(12,36) = 1.5*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(12,37) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(12,39) = 1.317465098480519*m0rOther[5]*mnuOther-0.4391550328268398*cEOther[5]*mnuOther; 
  data->AEM_S(13,30) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(13,31) = 1.341640786499874*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(13,32) = 1.341640786499874*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(13,33) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(13,34) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(13,35) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(13,36) = 1.31746509848052*m0rOther[8]*mnuOther-0.4391550328268399*cEOther[8]*mnuOther+1.2*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(13,37) = 1.31746509848052*m0rOther[9]*mnuOther-0.4391550328268399*cEOther[9]*mnuOther+1.2*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(13,38) = 1.31746509848052*m0rOther[6]*mnuOther-0.4391550328268399*cEOther[6]*mnuOther; 
  data->AEM_S(13,39) = 1.31746509848052*m0rOther[7]*mnuOther-0.4391550328268399*cEOther[7]*mnuOther; 
  data->AEM_S(14,30) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(14,31) = 1.317465098480519*m0rOther[8]*mnuOther-0.4391550328268398*cEOther[8]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(14,32) = 1.5*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(14,33) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(14,34) = 0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(14,36) = 0.9583148474999099*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(14,37) = 1.341640786499874*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(14,38) = 0.8944271909999159*m0rOther[8]*mnuOther-0.2981423969999719*cEOther[8]*mnuOther+1.317465098480519*m0rOther[1]*mnuOther-0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(15,30) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(15,31) = 1.5*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(15,32) = 1.317465098480519*m0rOther[9]*mnuOther-0.4391550328268398*cEOther[9]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(15,33) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(15,35) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(15,36) = 1.341640786499874*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(15,37) = 0.9583148474999099*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(15,39) = 0.8944271909999159*m0rOther[9]*mnuOther-0.2981423969999719*cEOther[9]*mnuOther+1.317465098480519*m0rOther[2]*mnuOther-0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(16,30) = 1.5*m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(16,31) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(16,32) = 1.5*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(16,33) = 1.31746509848052*m0rOther[8]*mnuOther-0.4391550328268399*cEOther[8]*mnuOther+1.2*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(16,34) = 0.9583148474999099*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(16,35) = 1.341640786499874*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(16,36) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(16,37) = 1.2*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(16,38) = 1.31746509848052*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(17,30) = 1.5*m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(17,31) = 1.5*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(17,32) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(17,33) = 1.31746509848052*m0rOther[9]*mnuOther-0.4391550328268399*cEOther[9]*mnuOther+1.2*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(17,34) = 1.341640786499874*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(17,35) = 0.9583148474999099*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(17,36) = 1.2*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(17,37) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(17,39) = 1.31746509848052*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(18,30) = 1.5*m0rOther[8]*mnuOther-0.5*cEOther[8]*mnuOther; 
  data->AEM_S(18,31) = 1.317465098480519*m0rOther[4]*mnuOther-0.4391550328268398*cEOther[4]*mnuOther; 
  data->AEM_S(18,33) = 1.31746509848052*m0rOther[6]*mnuOther-0.4391550328268399*cEOther[6]*mnuOther; 
  data->AEM_S(18,34) = 0.8944271909999159*m0rOther[8]*mnuOther-0.2981423969999719*cEOther[8]*mnuOther+1.317465098480519*m0rOther[1]*mnuOther-0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(18,36) = 1.31746509848052*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(18,38) = 0.8944271909999159*m0rOther[4]*mnuOther-0.2981423969999719*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(19,30) = 1.5*m0rOther[9]*mnuOther-0.5*cEOther[9]*mnuOther; 
  data->AEM_S(19,32) = 1.317465098480519*m0rOther[5]*mnuOther-0.4391550328268398*cEOther[5]*mnuOther; 
  data->AEM_S(19,33) = 1.31746509848052*m0rOther[7]*mnuOther-0.4391550328268399*cEOther[7]*mnuOther; 
  data->AEM_S(19,35) = 0.8944271909999159*m0rOther[9]*mnuOther-0.2981423969999719*cEOther[9]*mnuOther+1.317465098480519*m0rOther[2]*mnuOther-0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(19,37) = 1.31746509848052*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(19,39) = 0.8944271909999159*m0rOther[5]*mnuOther-0.2981423969999719*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[10]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2rSelf[0]*mnuSelf+m2rOther[0]*mnuOther; 
  mnuM2sum[1] = m2rSelf[1]*mnuSelf+m2rOther[1]*mnuOther; 
  mnuM2sum[2] = m2rSelf[2]*mnuSelf+m2rOther[2]*mnuOther; 
  mnuM2sum[3] = m2rSelf[3]*mnuSelf+m2rOther[3]*mnuOther; 
  mnuM2sum[4] = m2rSelf[4]*mnuSelf+m2rOther[4]*mnuOther; 
  mnuM2sum[5] = m2rSelf[5]*mnuSelf+m2rOther[5]*mnuOther; 
  mnuM2sum[6] = m2rSelf[6]*mnuSelf+m2rOther[6]*mnuOther; 
  mnuM2sum[7] = m2rSelf[7]*mnuSelf+m2rOther[7]*mnuOther; 
  mnuM2sum[8] = m2rSelf[8]*mnuSelf+m2rOther[8]*mnuOther; 
  mnuM2sum[9] = m2rSelf[9]*mnuSelf+m2rOther[9]*mnuOther; 
 
  double m1Relax[10]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(20,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(20,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(20,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(20,6) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(20,7) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(20,8) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(20,9) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(21,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(21,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,3) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,4) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,5) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(21,6) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,7) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(21,8) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(22,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(22,3) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,4) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(22,5) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,6) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(22,7) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,9) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(23,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,1) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,2) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(23,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,6) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,7) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,8) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(23,9) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(24,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(24,1) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(24,2) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(24,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(24,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(24,6) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,7) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(24,8) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(25,1) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(25,2) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(25,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(25,6) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(25,7) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,9) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,0) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(26,1) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(26,2) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(26,3) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(26,4) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,5) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(26,6) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(26,7) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(26,8) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,0) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(27,1) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(27,2) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,3) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,4) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(27,5) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,6) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,7) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(27,9) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(28,0) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(28,1) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(28,3) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(28,4) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(28,6) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(28,8) = 0.2981423969999719*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(29,0) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(29,2) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(29,3) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(29,5) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,7) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,9) = 0.2981423969999719*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(20,10) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(20,11) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(20,12) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(20,13) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(20,14) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(20,15) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(20,16) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(20,17) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(20,18) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(20,19) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(21,10) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(21,11) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(21,12) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(21,13) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(21,14) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(21,15) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(21,16) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(21,17) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(21,18) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(22,10) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(22,11) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(22,12) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(22,13) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(22,14) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(22,15) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(22,16) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(22,17) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(22,19) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(23,10) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(23,11) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(23,12) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(23,13) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(23,14) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(23,15) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(23,16) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(23,17) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(23,18) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(23,19) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(24,10) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(24,11) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(24,12) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(24,13) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(24,14) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(24,16) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(24,17) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(24,18) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(25,10) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(25,11) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(25,12) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(25,13) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(25,15) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(25,16) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(25,17) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(25,19) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(26,10) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(26,11) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(26,12) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(26,13) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(26,14) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(26,15) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(26,16) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(26,17) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(26,18) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(27,10) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(27,11) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(27,12) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(27,13) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(27,14) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(27,15) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(27,16) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(27,17) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(27,19) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(28,10) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(28,11) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(28,13) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(28,14) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(28,16) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(28,18) = (-0.2981423969999719*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(29,10) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(29,12) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(29,13) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(29,15) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(29,17) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(29,19) = (-0.2981423969999719*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(20,20) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(20,21) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(20,22) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(20,23) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(20,24) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(20,25) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(20,26) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(20,27) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(20,28) = -0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(20,29) = -0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(21,20) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(21,21) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(21,22) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(21,23) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(21,24) = (-0.4391550328268398*m0rOther[8]*mnuOther)-0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(21,25) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(21,26) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(21,27) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(21,28) = -0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(22,20) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(22,21) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(22,22) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(22,23) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(22,24) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(22,25) = (-0.4391550328268398*m0rOther[9]*mnuOther)-0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(22,26) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(22,27) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(22,29) = -0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(23,20) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(23,21) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(23,22) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(23,23) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(23,24) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(23,25) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(23,26) = (-0.4391550328268399*m0rOther[8]*mnuOther)-0.4*m0rOther[7]*mnuOther-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(23,27) = (-0.4391550328268399*m0rOther[9]*mnuOther)-0.4*m0rOther[6]*mnuOther-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(23,28) = -0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(23,29) = -0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(24,20) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(24,21) = (-0.4391550328268398*m0rOther[8]*mnuOther)-0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(24,22) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(24,23) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(24,24) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(24,26) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(24,27) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(24,28) = (-0.2981423969999719*m0rOther[8]*mnuOther)-0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(25,20) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(25,21) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(25,22) = (-0.4391550328268398*m0rOther[9]*mnuOther)-0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(25,23) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(25,25) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(25,26) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(25,27) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(25,29) = (-0.2981423969999719*m0rOther[9]*mnuOther)-0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(26,20) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(26,21) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(26,22) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(26,23) = (-0.4391550328268399*m0rOther[8]*mnuOther)-0.4*m0rOther[7]*mnuOther-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(26,24) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(26,25) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(26,26) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.31943828249997*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(26,27) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(26,28) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(27,20) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(27,21) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(27,22) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(27,23) = (-0.4391550328268399*m0rOther[9]*mnuOther)-0.4*m0rOther[6]*mnuOther-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(27,24) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(27,25) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(27,26) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(27,27) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(27,29) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(28,20) = -0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(28,21) = -0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(28,23) = -0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(28,24) = (-0.2981423969999719*m0rOther[8]*mnuOther)-0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(28,26) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(28,28) = (-0.2981423969999719*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(29,20) = -0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(29,22) = -0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(29,23) = -0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(29,25) = (-0.2981423969999719*m0rOther[9]*mnuOther)-0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(29,27) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(29,29) = (-0.2981423969999719*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(20,30) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(20,31) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(20,32) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(20,33) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(20,34) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(20,35) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(20,36) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(20,37) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(20,38) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(20,39) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(21,30) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(21,31) = 0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(21,32) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(21,33) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(21,34) = 0.4391550328268398*cMOther[8]*mnuOther+0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(21,35) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(21,36) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(21,37) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(21,38) = 0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(22,30) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(22,31) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(22,32) = 0.4472135954999579*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(22,33) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(22,34) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(22,35) = 0.4391550328268398*cMOther[9]*mnuOther+0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(22,36) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(22,37) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(22,39) = 0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(23,30) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(23,31) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(23,32) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(23,33) = 0.4472135954999579*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(23,34) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(23,35) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(23,36) = 0.4391550328268399*cMOther[8]*mnuOther+0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(23,37) = 0.4391550328268399*cMOther[9]*mnuOther+0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(23,38) = 0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(23,39) = 0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(24,30) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(24,31) = 0.4391550328268398*cMOther[8]*mnuOther+0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(24,32) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(24,33) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(24,34) = 0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(24,36) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(24,37) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(24,38) = 0.2981423969999719*cMOther[8]*mnuOther+0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(25,30) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(25,31) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(25,32) = 0.4391550328268398*cMOther[9]*mnuOther+0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(25,33) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(25,35) = 0.31943828249997*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(25,36) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(25,37) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(25,39) = 0.2981423969999719*cMOther[9]*mnuOther+0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(26,30) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(26,31) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(26,32) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(26,33) = 0.4391550328268399*cMOther[8]*mnuOther+0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(26,34) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(26,35) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(26,36) = 0.4472135954999579*cMOther[5]*mnuOther+0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(26,37) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(26,38) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(27,30) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(27,31) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(27,32) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(27,33) = 0.4391550328268399*cMOther[9]*mnuOther+0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(27,34) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(27,35) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(27,36) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(27,37) = 0.31943828249997*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(27,39) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(28,30) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(28,31) = 0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(28,33) = 0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(28,34) = 0.2981423969999719*cMOther[8]*mnuOther+0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(28,36) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(28,38) = 0.2981423969999719*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(29,30) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(29,32) = 0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(29,33) = 0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(29,35) = 0.2981423969999719*cMOther[9]*mnuOther+0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(29,37) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(29,39) = 0.2981423969999719*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(30,0) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.25*m0rSelf[7]*uSelf[7]*mnuSelf-0.25*m0rSelf[6]*uSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(30,1) = (-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(30,2) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,3) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,4) = (-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(30,5) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[2]*uSelf[9]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(30,6) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(30,7) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(30,8) = (-0.149071198499986*m0rSelf[4]*uSelf[8]*mnuSelf)-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.149071198499986*uSelf[4]*m0rSelf[8]*mnuSelf-0.25*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(30,9) = (-0.149071198499986*m0rSelf[5]*uSelf[9]*mnuSelf)-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.149071198499986*uSelf[5]*m0rSelf[9]*mnuSelf-0.25*uSelf[0]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(31,0) = (-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,1) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3833333333333334*m0rSelf[8]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[8]*mnuSelf-0.45*m0rSelf[7]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(31,2) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,3) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,4) = (-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,5) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,6) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,7) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[9]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(31,8) = (-0.2499586742703185*m0rSelf[8]*uSelf[8]*mnuSelf)-0.3833333333333334*m0rSelf[1]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[1]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[4]*mnuSelf+0.4391550328268398*m1rSelf[4]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(31,9) = (-0.149071198499986*m0rSelf[7]*uSelf[9]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.149071198499986*uSelf[7]*m0rSelf[9]*mnuSelf-0.25*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(32,0) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,1) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,2) = (-0.3833333333333334*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[2]*uSelf[9]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[9]*mnuSelf-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.45*m0rSelf[6]*uSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(32,3) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(32,4) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(32,5) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[9]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,6) = (-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(32,7) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.21957751641342*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,8) = (-0.149071198499986*m0rSelf[6]*uSelf[8]*mnuSelf)-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.149071198499986*uSelf[6]*m0rSelf[8]*mnuSelf-0.25*uSelf[2]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(32,9) = (-0.2499586742703185*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3833333333333334*m0rSelf[2]*uSelf[9]*mnuSelf-0.3833333333333334*uSelf[2]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[5]*mnuSelf+0.4391550328268398*m1rSelf[5]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(33,0) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,1) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,2) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,3) = (-0.3833333333333334*m0rSelf[9]*uSelf[9]*mnuSelf)-0.175662013130736*m0rSelf[6]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[9]*mnuSelf-0.175662013130736*uSelf[6]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[8]*mnuSelf-0.175662013130736*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[8]*mnuSelf-0.175662013130736*uSelf[7]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[8]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[6]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.2*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.2*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.81*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(33,4) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,5) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,6) = (-0.175662013130736*m0rSelf[3]*uSelf[9]*mnuSelf)-0.175662013130736*uSelf[3]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[8]*mnuSelf+0.43915503282684*m1rSelf[8]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,7) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[4]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[9]*mnuSelf+0.43915503282684*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[9]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[8]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[8]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,8) = (-0.3833333333333334*m0rSelf[3]*uSelf[8]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[8]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[7]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[6]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[6]*mnuSelf+0.43915503282684*m1rSelf[6]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[6]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[6]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,9) = (-0.3833333333333334*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[7]*mnuSelf+0.43915503282684*m1rSelf[7]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[7]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[7]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[6]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,0) = (-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(34,1) = (-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(34,2) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,3) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(34,4) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.4621212121212121*m0rSelf[8]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[6]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(34,5) = (-0.2195775164134199*m0rSelf[6]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[6]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[7]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(34,6) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(34,7) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[8]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,8) = (-0.4621212121212121*m0rSelf[4]*uSelf[8]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[8]*mnuSelf+0.2981423969999719*m1rSelf[8]*mnuSelf-0.4621212121212121*uSelf[4]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[5]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[1]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(34,9) = (-0.25*m0rSelf[4]*uSelf[9]*mnuSelf)-0.25*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[7]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[6]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(35,0) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[2]*uSelf[9]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(35,1) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,2) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[9]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(35,3) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(35,4) = (-0.2195775164134199*m0rSelf[6]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[6]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[7]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(35,5) = (-0.4621212121212121*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3273268353539885*m0rSelf[2]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[9]*mnuSelf-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(35,6) = (-0.21957751641342*m0rSelf[4]*uSelf[9]*mnuSelf)-0.21957751641342*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,7) = (-0.3273268353539885*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3273268353539885*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(35,8) = (-0.25*m0rSelf[5]*uSelf[8]*mnuSelf)-0.25*uSelf[5]*m0rSelf[8]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(35,9) = (-0.4621212121212121*m0rSelf[5]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[9]*mnuSelf+0.2981423969999719*m1rSelf[9]*mnuSelf-0.4621212121212121*uSelf[5]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[7]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[2]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,0) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(36,1) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,2) = (-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(36,3) = (-0.175662013130736*m0rSelf[3]*uSelf[9]*mnuSelf)-0.175662013130736*uSelf[3]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[8]*mnuSelf+0.43915503282684*m1rSelf[8]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(36,4) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,5) = (-0.21957751641342*m0rSelf[4]*uSelf[9]*mnuSelf)-0.21957751641342*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(36,6) = (-0.3833333333333334*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1254728665219542*m0rSelf[6]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[9]*mnuSelf-0.1254728665219542*uSelf[6]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[9]*mnuSelf-0.4621212121212121*m0rSelf[8]*uSelf[8]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[8]*mnuSelf-0.29277002188456*uSelf[7]*m0rSelf[8]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[8]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[7]*mnuSelf-0.9642857142857143*m0rSelf[6]*uSelf[6]*mnuSelf-0.2874944542499729*m0rSelf[2]*uSelf[6]*mnuSelf-0.2874944542499729*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(36,7) = (-0.1928571428571429*m0rSelf[8]*uSelf[9]*mnuSelf)-0.29277002188456*m0rSelf[7]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[9]*mnuSelf-0.1928571428571429*uSelf[8]*m0rSelf[9]*mnuSelf-0.29277002188456*uSelf[7]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[9]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[8]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[8]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,8) = (-0.1928571428571429*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[7]*m0rSelf[9]*mnuSelf-0.4621212121212121*m0rSelf[6]*uSelf[8]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[8]*mnuSelf-0.4621212121212121*uSelf[6]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[2]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[5]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[5]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,9) = (-0.3833333333333334*m0rSelf[6]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[6]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[7]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[7]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.1254728665219543*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(37,0) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(37,1) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[9]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(37,2) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.21957751641342*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(37,3) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[4]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[9]*mnuSelf+0.43915503282684*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[9]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[8]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[8]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(37,4) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[8]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(37,5) = (-0.3273268353539885*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3273268353539885*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(37,6) = (-0.1928571428571429*m0rSelf[8]*uSelf[9]*mnuSelf)-0.29277002188456*m0rSelf[7]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[9]*mnuSelf-0.1928571428571429*uSelf[8]*m0rSelf[9]*mnuSelf-0.29277002188456*uSelf[7]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[9]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[8]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[8]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(37,7) = (-0.4621212121212121*m0rSelf[9]*uSelf[9]*mnuSelf)-0.29277002188456*m0rSelf[6]*uSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[9]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[9]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[8]*mnuSelf-0.1254728665219542*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[8]*mnuSelf-0.1254728665219542*uSelf[7]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[8]*mnuSelf-0.9642857142857143*m0rSelf[7]*uSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[1]*uSelf[7]*mnuSelf-0.2874944542499729*uSelf[1]*m0rSelf[7]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[6]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(37,8) = (-0.1928571428571429*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[6]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[7]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[7]*m0rSelf[8]*mnuSelf-0.1254728665219543*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(37,9) = (-0.4621212121212121*m0rSelf[7]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[1]*uSelf[9]*mnuSelf-0.4621212121212121*uSelf[7]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[1]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[4]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(38,0) = (-0.149071198499986*m0rSelf[4]*uSelf[8]*mnuSelf)-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.149071198499986*uSelf[4]*m0rSelf[8]*mnuSelf-0.25*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(38,1) = (-0.2499586742703185*m0rSelf[8]*uSelf[8]*mnuSelf)-0.3833333333333334*m0rSelf[1]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[1]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[4]*mnuSelf+0.4391550328268398*m1rSelf[4]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(38,2) = (-0.149071198499986*m0rSelf[6]*uSelf[8]*mnuSelf)-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.149071198499986*uSelf[6]*m0rSelf[8]*mnuSelf-0.25*uSelf[2]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(38,3) = (-0.3833333333333334*m0rSelf[3]*uSelf[8]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[8]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[7]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[6]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[6]*mnuSelf+0.43915503282684*m1rSelf[6]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[6]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[6]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(38,4) = (-0.4621212121212121*m0rSelf[4]*uSelf[8]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[8]*mnuSelf+0.2981423969999719*m1rSelf[8]*mnuSelf-0.4621212121212121*uSelf[4]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[5]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[1]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(38,5) = (-0.25*m0rSelf[5]*uSelf[8]*mnuSelf)-0.25*uSelf[5]*m0rSelf[8]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(38,6) = (-0.1928571428571429*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[7]*m0rSelf[9]*mnuSelf-0.4621212121212121*m0rSelf[6]*uSelf[8]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[8]*mnuSelf-0.4621212121212121*uSelf[6]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[2]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[5]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[5]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(38,7) = (-0.1928571428571429*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[6]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[7]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[7]*m0rSelf[8]*mnuSelf-0.1254728665219543*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(38,8) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.5898601398601399*m0rSelf[8]*uSelf[8]*mnuSelf-0.2499586742703185*m0rSelf[1]*uSelf[8]*mnuSelf-0.2499586742703185*uSelf[1]*m0rSelf[8]*mnuSelf-0.3833333333333334*m0rSelf[7]*uSelf[7]*mnuSelf-0.4621212121212121*m0rSelf[6]*uSelf[6]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[6]*mnuSelf-0.149071198499986*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.4621212121212121*m0rSelf[4]*uSelf[4]*mnuSelf-0.149071198499986*m0rSelf[0]*uSelf[4]*mnuSelf+0.2981423969999719*m1rSelf[4]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[4]*mnuSelf-0.3833333333333334*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3833333333333334*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(38,9) = (-0.25*m0rSelf[8]*uSelf[9]*mnuSelf)-0.25*uSelf[8]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[7]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[7]*mnuSelf; 
  data->AEM_S(39,0) = (-0.149071198499986*m0rSelf[5]*uSelf[9]*mnuSelf)-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.149071198499986*uSelf[5]*m0rSelf[9]*mnuSelf-0.25*uSelf[0]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(39,1) = (-0.149071198499986*m0rSelf[7]*uSelf[9]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.149071198499986*uSelf[7]*m0rSelf[9]*mnuSelf-0.25*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(39,2) = (-0.2499586742703185*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3833333333333334*m0rSelf[2]*uSelf[9]*mnuSelf-0.3833333333333334*uSelf[2]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[5]*mnuSelf+0.4391550328268398*m1rSelf[5]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(39,3) = (-0.3833333333333334*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[7]*mnuSelf+0.43915503282684*m1rSelf[7]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[7]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[7]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[6]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(39,4) = (-0.25*m0rSelf[4]*uSelf[9]*mnuSelf)-0.25*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[7]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[6]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(39,5) = (-0.4621212121212121*m0rSelf[5]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[9]*mnuSelf+0.2981423969999719*m1rSelf[9]*mnuSelf-0.4621212121212121*uSelf[5]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[7]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[2]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(39,6) = (-0.3833333333333334*m0rSelf[6]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[6]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[7]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[7]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.1254728665219543*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(39,7) = (-0.4621212121212121*m0rSelf[7]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[1]*uSelf[9]*mnuSelf-0.4621212121212121*uSelf[7]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[1]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[4]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(39,8) = (-0.25*m0rSelf[8]*uSelf[9]*mnuSelf)-0.25*uSelf[8]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[7]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[7]*mnuSelf; 
  data->AEM_S(39,9) = (-0.5898601398601399*m0rSelf[9]*uSelf[9]*mnuSelf)-0.2499586742703185*m0rSelf[2]*uSelf[9]*mnuSelf-0.2499586742703185*uSelf[2]*m0rSelf[9]*mnuSelf-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.4621212121212121*m0rSelf[7]*uSelf[7]*mnuSelf-0.149071198499986*m0rSelf[1]*uSelf[7]*mnuSelf-0.149071198499986*uSelf[1]*m0rSelf[7]*mnuSelf-0.3833333333333334*m0rSelf[6]*uSelf[6]*mnuSelf-0.4621212121212121*m0rSelf[5]*uSelf[5]*mnuSelf-0.149071198499986*m0rSelf[0]*uSelf[5]*mnuSelf+0.2981423969999719*m1rSelf[5]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3833333333333334*m0rSelf[3]*uSelf[3]*mnuSelf-0.3833333333333334*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(30,20) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.25*m0rOther[7]*uOther[7]*mnuOther+0.25*m0rOther[6]*uOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(30,21) = 0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(30,22) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,23) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,24) = 0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[8]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(30,25) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[9]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(30,26) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(30,27) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(30,28) = 0.149071198499986*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.149071198499986*uOther[4]*m0rOther[8]*mnuOther+0.25*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[6]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[4]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[4]*mnuOther; 
  data->AEM_S(30,29) = 0.149071198499986*m0rOther[5]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.149071198499986*uOther[5]*m0rOther[9]*mnuOther+0.25*uOther[0]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[7]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[5]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[5]*mnuOther; 
  data->AEM_S(31,20) = 0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(31,21) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[8]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[8]*mnuOther+0.45*m0rOther[7]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(31,22) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,23) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,24) = 0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[8]*mnuOther-0.4391550328268398*m1rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(31,25) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(31,26) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.21957751641342*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.21957751641342*uOther[2]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,27) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.21957751641342*m0rOther[2]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.21957751641342*uOther[2]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(31,28) = 0.2499586742703185*m0rOther[8]*uOther[8]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[8]*mnuOther+0.3833333333333334*uOther[1]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[6]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[4]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[4]*mnuOther-0.4391550328268398*m1rOther[4]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(31,29) = 0.149071198499986*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.149071198499986*uOther[7]*m0rOther[9]*mnuOther+0.25*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[5]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[5]*mnuOther; 
  data->AEM_S(32,20) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(32,21) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(32,22) = 0.3833333333333334*m0rOther[9]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[9]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.45*m0rOther[6]*uOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(32,23) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(32,24) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(32,25) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[9]*mnuOther-0.4391550328268398*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[9]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(32,26) = 0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.21957751641342*m0rOther[1]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.21957751641342*uOther[1]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(32,27) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.21957751641342*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.21957751641342*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(32,28) = 0.149071198499986*m0rOther[6]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.149071198499986*uOther[6]*m0rOther[8]*mnuOther+0.25*uOther[2]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[6]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[4]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(32,29) = 0.2499586742703185*m0rOther[9]*uOther[9]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[9]*mnuOther+0.3833333333333334*uOther[2]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[7]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[5]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[5]*mnuOther-0.4391550328268398*m1rOther[5]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(33,20) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,21) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,22) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(33,23) = 0.3833333333333334*m0rOther[9]*uOther[9]*mnuOther+0.175662013130736*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[9]*mnuOther+0.175662013130736*uOther[6]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[8]*mnuOther+0.175662013130736*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[8]*mnuOther+0.175662013130736*uOther[7]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[8]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[7]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[7]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[6]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[6]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.2*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.2*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(33,24) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,25) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,26) = 0.175662013130736*m0rOther[3]*uOther[9]*mnuOther+0.175662013130736*uOther[3]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.21957751641342*m0rOther[0]*uOther[8]*mnuOther-0.43915503282684*m1rOther[8]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.21957751641342*uOther[0]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(33,27) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*m0rOther[0]*uOther[9]*mnuOther-0.43915503282684*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[9]*mnuOther+0.21957751641342*uOther[0]*m0rOther[9]*mnuOther+0.175662013130736*m0rOther[3]*uOther[8]*mnuOther+0.175662013130736*uOther[3]*m0rOther[8]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,28) = 0.3833333333333334*m0rOther[3]*uOther[8]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[8]*mnuOther+0.175662013130736*m0rOther[3]*uOther[7]*mnuOther+0.175662013130736*uOther[3]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[6]*mnuOther+0.21957751641342*m0rOther[0]*uOther[6]*mnuOther-0.43915503282684*m1rOther[6]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[6]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[6]*mnuOther+0.21957751641342*uOther[0]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[4]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[3]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(33,29) = 0.3833333333333334*m0rOther[3]*uOther[9]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*m0rOther[0]*uOther[7]*mnuOther-0.43915503282684*m1rOther[7]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[7]*mnuOther+0.21957751641342*uOther[0]*m0rOther[7]*mnuOther+0.175662013130736*m0rOther[3]*uOther[6]*mnuOther+0.175662013130736*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[5]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[3]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(34,20) = 0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[8]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(34,21) = 0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[8]*mnuOther-0.4391550328268398*m1rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(34,22) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(34,23) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(34,24) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[8]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[8]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[2]*uOther[6]*mnuOther+0.159719141249985*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(34,25) = 0.2195775164134199*m0rOther[6]*uOther[9]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[8]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(34,26) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[8]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(34,27) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.21957751641342*m0rOther[5]*uOther[8]*mnuOther+0.21957751641342*uOther[5]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(34,28) = 0.4621212121212121*m0rOther[4]*uOther[8]*mnuOther+0.149071198499986*m0rOther[0]*uOther[8]*mnuOther-0.2981423969999719*m1rOther[8]*mnuOther+0.4621212121212121*uOther[4]*m0rOther[8]*mnuOther+0.149071198499986*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[7]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[6]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[4]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[4]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[3]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[1]*mnuOther-0.4391550328268398*m1rOther[1]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(34,29) = 0.25*m0rOther[4]*uOther[9]*mnuOther+0.25*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[7]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[7]*mnuOther+0.21957751641342*m0rOther[5]*uOther[6]*mnuOther+0.21957751641342*uOther[5]*m0rOther[6]*mnuOther; 
  data->AEM_S(35,20) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[9]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(35,21) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(35,22) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[9]*mnuOther-0.4391550328268398*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[9]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(35,23) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(35,24) = 0.2195775164134199*m0rOther[6]*uOther[9]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[8]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(35,25) = 0.4621212121212121*m0rOther[9]*uOther[9]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[9]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*uOther[1]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(35,26) = 0.21957751641342*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(35,27) = 0.3273268353539885*m0rOther[3]*uOther[9]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(35,28) = 0.25*m0rOther[5]*uOther[8]*mnuOther+0.25*uOther[5]*m0rOther[8]*mnuOther+0.21957751641342*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*uOther[4]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[6]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[6]*mnuOther; 
  data->AEM_S(35,29) = 0.4621212121212121*m0rOther[5]*uOther[9]*mnuOther+0.149071198499986*m0rOther[0]*uOther[9]*mnuOther-0.2981423969999719*m1rOther[9]*mnuOther+0.4621212121212121*uOther[5]*m0rOther[9]*mnuOther+0.149071198499986*uOther[0]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[7]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[6]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[5]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[5]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[3]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[2]*mnuOther-0.4391550328268398*m1rOther[2]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(36,20) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(36,21) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.21957751641342*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.21957751641342*uOther[2]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(36,22) = 0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.21957751641342*m0rOther[1]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.21957751641342*uOther[1]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(36,23) = 0.175662013130736*m0rOther[3]*uOther[9]*mnuOther+0.175662013130736*uOther[3]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.21957751641342*m0rOther[0]*uOther[8]*mnuOther-0.43915503282684*m1rOther[8]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.21957751641342*uOther[0]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(36,24) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[8]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(36,25) = 0.21957751641342*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(36,26) = 0.3833333333333334*m0rOther[9]*uOther[9]*mnuOther+0.1254728665219542*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[9]*mnuOther+0.1254728665219542*uOther[6]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[9]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[8]*mnuOther+0.29277002188456*m0rOther[7]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[8]*mnuOther+0.29277002188456*uOther[7]*m0rOther[8]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[8]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[7]*mnuOther+0.351382110749967*uOther[1]*m0rOther[7]*mnuOther+0.9642857142857143*m0rOther[6]*uOther[6]*mnuOther+0.2874944542499729*m0rOther[2]*uOther[6]*mnuOther+0.2874944542499729*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(36,27) = 0.1928571428571429*m0rOther[8]*uOther[9]*mnuOther+0.29277002188456*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[9]*mnuOther+0.1928571428571429*uOther[8]*m0rOther[9]*mnuOther+0.29277002188456*uOther[7]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[9]*mnuOther+0.29277002188456*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[8]*mnuOther+0.29277002188456*uOther[6]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[8]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(36,28) = 0.1928571428571429*m0rOther[7]*uOther[9]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[9]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[8]*mnuOther+0.149071198499986*m0rOther[2]*uOther[8]*mnuOther+0.4621212121212121*uOther[6]*m0rOther[8]*mnuOther+0.149071198499986*uOther[2]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[6]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[6]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[5]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[5]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[4]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(36,29) = 0.3833333333333334*m0rOther[6]*uOther[9]*mnuOther+0.3833333333333334*uOther[6]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[8]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.1254728665219543*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(37,20) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(37,21) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.21957751641342*m0rOther[2]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.21957751641342*uOther[2]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(37,22) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.21957751641342*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.21957751641342*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(37,23) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*m0rOther[0]*uOther[9]*mnuOther-0.43915503282684*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[9]*mnuOther+0.21957751641342*uOther[0]*m0rOther[9]*mnuOther+0.175662013130736*m0rOther[3]*uOther[8]*mnuOther+0.175662013130736*uOther[3]*m0rOther[8]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(37,24) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.21957751641342*m0rOther[5]*uOther[8]*mnuOther+0.21957751641342*uOther[5]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(37,25) = 0.3273268353539885*m0rOther[3]*uOther[9]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(37,26) = 0.1928571428571429*m0rOther[8]*uOther[9]*mnuOther+0.29277002188456*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[9]*mnuOther+0.1928571428571429*uOther[8]*m0rOther[9]*mnuOther+0.29277002188456*uOther[7]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[9]*mnuOther+0.29277002188456*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[8]*mnuOther+0.29277002188456*uOther[6]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[8]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(37,27) = 0.4621212121212121*m0rOther[9]*uOther[9]*mnuOther+0.29277002188456*m0rOther[6]*uOther[9]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[9]*mnuOther+0.29277002188456*uOther[6]*m0rOther[9]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[8]*mnuOther+0.1254728665219542*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[8]*mnuOther+0.1254728665219542*uOther[7]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[8]*mnuOther+0.9642857142857143*m0rOther[7]*uOther[7]*mnuOther+0.2874944542499729*m0rOther[1]*uOther[7]*mnuOther+0.2874944542499729*uOther[1]*m0rOther[7]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[6]*mnuOther+0.351382110749967*m0rOther[2]*uOther[6]*mnuOther+0.351382110749967*uOther[2]*m0rOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(37,28) = 0.1928571428571429*m0rOther[6]*uOther[9]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[8]*mnuOther+0.3833333333333334*uOther[7]*m0rOther[8]*mnuOther+0.1254728665219543*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.29277002188456*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(37,29) = 0.4621212121212121*m0rOther[7]*uOther[9]*mnuOther+0.149071198499986*m0rOther[1]*uOther[9]*mnuOther+0.4621212121212121*uOther[7]*m0rOther[9]*mnuOther+0.149071198499986*uOther[1]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[8]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[6]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[5]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[4]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(38,20) = 0.149071198499986*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.149071198499986*uOther[4]*m0rOther[8]*mnuOther+0.25*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[6]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[4]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[4]*mnuOther; 
  data->AEM_S(38,21) = 0.2499586742703185*m0rOther[8]*uOther[8]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[8]*mnuOther+0.3833333333333334*uOther[1]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[6]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[4]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[4]*mnuOther-0.4391550328268398*m1rOther[4]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(38,22) = 0.149071198499986*m0rOther[6]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.149071198499986*uOther[6]*m0rOther[8]*mnuOther+0.25*uOther[2]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[6]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[4]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(38,23) = 0.3833333333333334*m0rOther[3]*uOther[8]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[8]*mnuOther+0.175662013130736*m0rOther[3]*uOther[7]*mnuOther+0.175662013130736*uOther[3]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[6]*mnuOther+0.21957751641342*m0rOther[0]*uOther[6]*mnuOther-0.43915503282684*m1rOther[6]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[6]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[6]*mnuOther+0.21957751641342*uOther[0]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[4]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[3]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(38,24) = 0.4621212121212121*m0rOther[4]*uOther[8]*mnuOther+0.149071198499986*m0rOther[0]*uOther[8]*mnuOther-0.2981423969999719*m1rOther[8]*mnuOther+0.4621212121212121*uOther[4]*m0rOther[8]*mnuOther+0.149071198499986*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[7]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[6]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[4]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[4]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[3]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[1]*mnuOther-0.4391550328268398*m1rOther[1]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(38,25) = 0.25*m0rOther[5]*uOther[8]*mnuOther+0.25*uOther[5]*m0rOther[8]*mnuOther+0.21957751641342*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*uOther[4]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[6]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[6]*mnuOther; 
  data->AEM_S(38,26) = 0.1928571428571429*m0rOther[7]*uOther[9]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[9]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[8]*mnuOther+0.149071198499986*m0rOther[2]*uOther[8]*mnuOther+0.4621212121212121*uOther[6]*m0rOther[8]*mnuOther+0.149071198499986*uOther[2]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[6]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[6]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[5]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[5]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[4]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(38,27) = 0.1928571428571429*m0rOther[6]*uOther[9]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[8]*mnuOther+0.3833333333333334*uOther[7]*m0rOther[8]*mnuOther+0.1254728665219543*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.29277002188456*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(38,28) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.5898601398601399*m0rOther[8]*uOther[8]*mnuOther+0.2499586742703185*m0rOther[1]*uOther[8]*mnuOther+0.2499586742703185*uOther[1]*m0rOther[8]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[7]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[6]*mnuOther+0.149071198499986*m0rOther[2]*uOther[6]*mnuOther+0.149071198499986*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.4621212121212121*m0rOther[4]*uOther[4]*mnuOther+0.149071198499986*m0rOther[0]*uOther[4]*mnuOther-0.2981423969999719*m1rOther[4]*mnuOther+0.149071198499986*uOther[0]*m0rOther[4]*mnuOther+0.3833333333333334*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(38,29) = 0.25*m0rOther[8]*uOther[9]*mnuOther+0.25*uOther[8]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[7]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[7]*mnuOther; 
  data->AEM_S(39,20) = 0.149071198499986*m0rOther[5]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.149071198499986*uOther[5]*m0rOther[9]*mnuOther+0.25*uOther[0]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[7]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[5]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[5]*mnuOther; 
  data->AEM_S(39,21) = 0.149071198499986*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.149071198499986*uOther[7]*m0rOther[9]*mnuOther+0.25*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[5]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[5]*mnuOther; 
  data->AEM_S(39,22) = 0.2499586742703185*m0rOther[9]*uOther[9]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[9]*mnuOther+0.3833333333333334*uOther[2]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[7]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[5]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[5]*mnuOther-0.4391550328268398*m1rOther[5]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(39,23) = 0.3833333333333334*m0rOther[3]*uOther[9]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*m0rOther[0]*uOther[7]*mnuOther-0.43915503282684*m1rOther[7]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[7]*mnuOther+0.21957751641342*uOther[0]*m0rOther[7]*mnuOther+0.175662013130736*m0rOther[3]*uOther[6]*mnuOther+0.175662013130736*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[5]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[3]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(39,24) = 0.25*m0rOther[4]*uOther[9]*mnuOther+0.25*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[7]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[7]*mnuOther+0.21957751641342*m0rOther[5]*uOther[6]*mnuOther+0.21957751641342*uOther[5]*m0rOther[6]*mnuOther; 
  data->AEM_S(39,25) = 0.4621212121212121*m0rOther[5]*uOther[9]*mnuOther+0.149071198499986*m0rOther[0]*uOther[9]*mnuOther-0.2981423969999719*m1rOther[9]*mnuOther+0.4621212121212121*uOther[5]*m0rOther[9]*mnuOther+0.149071198499986*uOther[0]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[7]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[6]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[5]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[5]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[3]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[2]*mnuOther-0.4391550328268398*m1rOther[2]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(39,26) = 0.3833333333333334*m0rOther[6]*uOther[9]*mnuOther+0.3833333333333334*uOther[6]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[8]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.1254728665219543*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(39,27) = 0.4621212121212121*m0rOther[7]*uOther[9]*mnuOther+0.149071198499986*m0rOther[1]*uOther[9]*mnuOther+0.4621212121212121*uOther[7]*m0rOther[9]*mnuOther+0.149071198499986*uOther[1]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[8]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[6]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[5]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[4]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(39,28) = 0.25*m0rOther[8]*uOther[9]*mnuOther+0.25*uOther[8]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[7]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[7]*mnuOther; 
  data->AEM_S(39,29) = 0.5898601398601399*m0rOther[9]*uOther[9]*mnuOther+0.2499586742703185*m0rOther[2]*uOther[9]*mnuOther+0.2499586742703185*uOther[2]*m0rOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.4621212121212121*m0rOther[7]*uOther[7]*mnuOther+0.149071198499986*m0rOther[1]*uOther[7]*mnuOther+0.149071198499986*uOther[1]*m0rOther[7]*mnuOther+0.3833333333333334*m0rOther[6]*uOther[6]*mnuOther+0.4621212121212121*m0rOther[5]*uOther[5]*mnuOther+0.149071198499986*m0rOther[0]*uOther[5]*mnuOther-0.2981423969999719*m1rOther[5]*mnuOther+0.149071198499986*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3833333333333334*m0rOther[3]*uOther[3]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += (m1rOther[3]-1.0*m1rSelf[3])*betaGreenep1*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (m1rOther[4]-1.0*m1rSelf[4])*betaGreenep1*mnuSelf+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (m1rOther[5]-1.0*m1rSelf[5])*betaGreenep1*mnuSelf+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
  m1Relax[6] += (m1rOther[6]-1.0*m1rSelf[6])*betaGreenep1*mnuSelf+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (m1rOther[7]-1.0*m1rSelf[7])*betaGreenep1*mnuSelf+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
  m1Relax[8] += (m1rOther[8]-1.0*m1rSelf[8])*betaGreenep1*mnuSelf+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
  m1Relax[9] += (m1rOther[9]-1.0*m1rSelf[9])*betaGreenep1*mnuSelf+m1rSelf[9]*mnuSelf-1.0*m1rOther[9]*mnuOther; 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(30,10) = 0.25*cMSelf[9]*uSelf[9]*mnuSelf+0.25*cMSelf[8]*uSelf[8]*mnuSelf+0.25*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[6]*uSelf[6]*mnuSelf+0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[3]*uSelf[3]*mnuSelf+0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(30,11) = 0.2195775164134199*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[8]*mnuSelf+0.2500000000000001*cMSelf[5]*uSelf[7]*mnuSelf+0.2500000000000001*uSelf[5]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.25*cMSelf[2]*uSelf[3]*mnuSelf+0.25*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(30,12) = 0.2195775164134199*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[7]*mnuSelf+0.2500000000000001*cMSelf[4]*uSelf[6]*mnuSelf+0.2500000000000001*uSelf[4]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.25*cMSelf[1]*uSelf[3]*mnuSelf+0.25*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(30,13) = 0.2195775164134199*cMSelf[7]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[7]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[6]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[6]*cMSelf[8]*mnuSelf+0.2*cMSelf[6]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[7]*mnuSelf+0.2*uSelf[6]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(30,14) = 0.149071198499986*cMSelf[8]*uSelf[8]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[8]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[6]*uSelf[6]*mnuSelf+0.2500000000000001*cMSelf[2]*uSelf[6]*mnuSelf+0.2500000000000001*uSelf[2]*cMSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf+0.25*uSelf[0]*cMSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(30,15) = 0.149071198499986*cMSelf[9]*uSelf[9]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[7]*uSelf[7]*mnuSelf+0.2500000000000001*cMSelf[1]*uSelf[7]*mnuSelf+0.2500000000000001*uSelf[1]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[0]*uSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf+0.25*uSelf[0]*cMSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(30,16) = 0.2195775164134199*cMSelf[3]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[8]*mnuSelf+0.2*cMSelf[3]*uSelf[7]*mnuSelf+0.2*uSelf[3]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[5]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[6]*mnuSelf+0.25*cMSelf[0]*uSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf+0.223606797749979*uSelf[5]*cMSelf[6]*mnuSelf+0.159719141249985*uSelf[4]*cMSelf[6]*mnuSelf+0.25*uSelf[0]*cMSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf+0.2500000000000001*cMSelf[2]*uSelf[4]*mnuSelf+0.2500000000000001*uSelf[2]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(30,17) = 0.2195775164134199*cMSelf[3]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[4]*uSelf[7]*mnuSelf+0.25*cMSelf[0]*uSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf+0.159719141249985*uSelf[5]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[4]*cMSelf[7]*mnuSelf+0.25*uSelf[0]*cMSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf+0.2*cMSelf[3]*uSelf[6]*mnuSelf+0.2*uSelf[3]*cMSelf[6]*mnuSelf+0.2500000000000001*cMSelf[1]*uSelf[5]*mnuSelf+0.2500000000000001*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(30,18) = 0.149071198499986*cMSelf[4]*uSelf[8]*mnuSelf+0.25*cMSelf[0]*uSelf[8]*mnuSelf+1.5*m0rSelf[8]*mnuSelf+0.149071198499986*uSelf[4]*cMSelf[8]*mnuSelf+0.25*uSelf[0]*cMSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf+0.2195775164134199*cMSelf[3]*uSelf[6]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[6]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[4]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[4]*mnuSelf; 
  data->AEM_S(30,19) = 0.149071198499986*cMSelf[5]*uSelf[9]*mnuSelf+0.25*cMSelf[0]*uSelf[9]*mnuSelf+1.5*m0rSelf[9]*mnuSelf+0.149071198499986*uSelf[5]*cMSelf[9]*mnuSelf+0.25*uSelf[0]*cMSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf+0.2195775164134199*cMSelf[3]*uSelf[7]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[7]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[5]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[5]*mnuSelf; 
  data->AEM_S(31,10) = 0.2195775164134199*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[8]*mnuSelf+0.2500000000000001*cMSelf[5]*uSelf[7]*mnuSelf+0.2500000000000001*uSelf[5]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.25*cMSelf[2]*uSelf[3]*mnuSelf+0.25*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(31,11) = 0.25*cMSelf[9]*uSelf[9]*mnuSelf+0.3833333333333334*cMSelf[8]*uSelf[8]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[8]*mnuSelf+0.45*cMSelf[7]*uSelf[7]*mnuSelf+0.3928571428571428*cMSelf[6]*uSelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[6]*mnuSelf+0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.3928571428571428*cMSelf[4]*uSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.45*cMSelf[3]*uSelf[3]*mnuSelf+0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.45*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(31,12) = 0.2195775164134199*cMSelf[7]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[7]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[6]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[6]*cMSelf[8]*mnuSelf+0.2*cMSelf[6]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[7]*mnuSelf+0.2*uSelf[6]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(31,13) = 0.2195775164134199*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[8]*mnuSelf+0.4024922359499621*cMSelf[3]*uSelf[7]*mnuSelf+0.4024922359499621*uSelf[3]*cMSelf[7]*mnuSelf+0.2*cMSelf[5]*uSelf[6]*mnuSelf+0.3928571428571429*cMSelf[4]*uSelf[6]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf+0.2*uSelf[5]*cMSelf[6]*mnuSelf+0.3928571428571429*uSelf[4]*cMSelf[6]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[4]*mnuSelf+0.45*cMSelf[1]*uSelf[3]*mnuSelf+0.45*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(31,14) = 0.3273268353539885*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[8]*mnuSelf+1.317465098480519*m0rSelf[8]*mnuSelf+0.3273268353539885*uSelf[4]*cMSelf[8]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+0.223606797749979*cMSelf[5]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[5]*cMSelf[7]*mnuSelf+0.3928571428571429*cMSelf[3]*uSelf[6]*mnuSelf+0.3928571428571429*uSelf[3]*cMSelf[6]*mnuSelf+0.3928571428571428*cMSelf[1]*uSelf[4]*mnuSelf+0.3928571428571428*uSelf[1]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(31,15) = 0.2195775164134199*cMSelf[3]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[4]*uSelf[7]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf+0.159719141249985*uSelf[5]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[4]*cMSelf[7]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf+0.2*cMSelf[3]*uSelf[6]*mnuSelf+0.2*uSelf[3]*cMSelf[6]*mnuSelf+0.25*cMSelf[1]*uSelf[5]*mnuSelf+0.25*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(31,16) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999832*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999832*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(31,17) = 0.149071198499986*cMSelf[9]*uSelf[9]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[1]*uSelf[7]*mnuSelf+0.25*uSelf[1]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[5]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(31,18) = 0.1309307341415954*cMSelf[8]*uSelf[8]*mnuSelf+0.1928571428571429*cMSelf[1]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[1]*cMSelf[8]*mnuSelf+0.1963961012123931*cMSelf[7]*uSelf[7]*mnuSelf+0.1402829294374236*cMSelf[6]*uSelf[6]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[6]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[6]*mnuSelf+0.1402829294374236*cMSelf[4]*uSelf[4]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[4]*mnuSelf+1.317465098480519*m0rSelf[4]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[3]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(32,10) = 0.2195775164134199*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[7]*mnuSelf+0.2500000000000001*cMSelf[4]*uSelf[6]*mnuSelf+0.2500000000000001*uSelf[4]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.25*cMSelf[1]*uSelf[3]*mnuSelf+0.25*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(32,11) = 0.2195775164134199*cMSelf[7]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[7]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[6]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[6]*cMSelf[8]*mnuSelf+0.2*cMSelf[6]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[7]*mnuSelf+0.2*uSelf[6]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(32,12) = 0.3833333333333334*cMSelf[9]*uSelf[9]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[9]*mnuSelf+0.25*cMSelf[8]*uSelf[8]*mnuSelf+0.3928571428571428*cMSelf[7]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[7]*mnuSelf+0.45*cMSelf[6]*uSelf[6]*mnuSelf+0.3928571428571428*cMSelf[5]*uSelf[5]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.45*cMSelf[3]*uSelf[3]*mnuSelf+0.45*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(32,13) = 0.1963961012123931*cMSelf[3]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[8]*mnuSelf+0.3928571428571429*cMSelf[5]*uSelf[7]*mnuSelf+0.2*cMSelf[4]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf+0.3928571428571429*uSelf[5]*cMSelf[7]*mnuSelf+0.2*uSelf[4]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+0.4024922359499621*cMSelf[3]*uSelf[6]*mnuSelf+0.4024922359499621*uSelf[3]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.45*cMSelf[2]*uSelf[3]*mnuSelf+0.45*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(32,14) = 0.2195775164134199*cMSelf[3]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[8]*mnuSelf+0.2*cMSelf[3]*uSelf[7]*mnuSelf+0.2*uSelf[3]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[5]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[6]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf+0.223606797749979*uSelf[5]*cMSelf[6]*mnuSelf+0.159719141249985*uSelf[4]*cMSelf[6]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf+0.25*cMSelf[2]*uSelf[4]*mnuSelf+0.25*uSelf[2]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(32,15) = 0.3273268353539885*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[9]*mnuSelf+1.317465098480519*m0rSelf[9]*mnuSelf+0.3273268353539885*uSelf[5]*cMSelf[9]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+0.3928571428571429*cMSelf[3]*uSelf[7]*mnuSelf+0.3928571428571429*uSelf[3]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[4]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[4]*cMSelf[6]*mnuSelf+0.3928571428571428*cMSelf[2]*uSelf[5]*mnuSelf+0.3928571428571428*uSelf[2]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(32,16) = 0.149071198499986*cMSelf[8]*uSelf[8]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[8]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[6]*uSelf[6]*mnuSelf+0.25*cMSelf[2]*uSelf[6]*mnuSelf+0.25*uSelf[2]*cMSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[4]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(32,17) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999832*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999832*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(32,19) = 0.1309307341415954*cMSelf[9]*uSelf[9]*mnuSelf+0.1928571428571429*cMSelf[2]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[2]*cMSelf[9]*mnuSelf+0.1402829294374236*cMSelf[7]*uSelf[7]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[7]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[6]*mnuSelf+0.1402829294374236*cMSelf[5]*uSelf[5]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[5]*mnuSelf+1.317465098480519*m0rSelf[5]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[3]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(33,10) = 0.2195775164134199*cMSelf[7]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[7]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[6]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[6]*cMSelf[8]*mnuSelf+0.2*cMSelf[6]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[7]*mnuSelf+0.2*uSelf[6]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[3]*cMSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf+0.25*uSelf[0]*cMSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf+0.25*cMSelf[1]*uSelf[2]*mnuSelf+0.25*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(33,11) = 0.2195775164134199*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[8]*mnuSelf+0.4024922359499621*cMSelf[3]*uSelf[7]*mnuSelf+0.4024922359499621*uSelf[3]*cMSelf[7]*mnuSelf+0.2*cMSelf[5]*uSelf[6]*mnuSelf+0.3928571428571429*cMSelf[4]*uSelf[6]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf+0.2*uSelf[5]*cMSelf[6]*mnuSelf+0.3928571428571429*uSelf[4]*cMSelf[6]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[4]*mnuSelf+0.45*cMSelf[1]*uSelf[3]*mnuSelf+0.45*uSelf[1]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.25*uSelf[0]*cMSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(33,12) = 0.1963961012123931*cMSelf[3]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[8]*mnuSelf+0.3928571428571429*cMSelf[5]*uSelf[7]*mnuSelf+0.2*cMSelf[4]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf+0.3928571428571429*uSelf[5]*cMSelf[7]*mnuSelf+0.2*uSelf[4]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+0.4024922359499621*cMSelf[3]*uSelf[6]*mnuSelf+0.4024922359499621*uSelf[3]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.45*cMSelf[2]*uSelf[3]*mnuSelf+0.45*uSelf[2]*cMSelf[3]*mnuSelf+0.25*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.25*uSelf[0]*cMSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(33,13) = 0.3833333333333334*cMSelf[9]*uSelf[9]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[9]*mnuSelf+0.3833333333333334*cMSelf[8]*uSelf[8]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[8]*mnuSelf+0.5928571428571429*cMSelf[7]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[7]*mnuSelf+0.5928571428571429*cMSelf[6]*uSelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[6]*mnuSelf+0.3928571428571428*cMSelf[5]*uSelf[5]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.3928571428571428*cMSelf[4]*uSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.65*cMSelf[3]*uSelf[3]*mnuSelf+0.45*cMSelf[2]*uSelf[2]*mnuSelf+0.45*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(33,14) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999831*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999831*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(33,15) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999831*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999831*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(33,16) = 0.175662013130736*cMSelf[3]*uSelf[9]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[9]*mnuSelf+0.3273268353539885*cMSelf[4]*uSelf[8]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[8]*mnuSelf+1.31746509848052*m0rSelf[8]*mnuSelf+0.3273268353539885*uSelf[4]*cMSelf[8]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[8]*mnuSelf-0.43915503282684*cESelf[8]*mnuSelf+0.3513821107499669*cMSelf[5]*uSelf[7]*mnuSelf+0.1788854381999831*cMSelf[4]*uSelf[7]*mnuSelf+0.2*cMSelf[0]*uSelf[7]*mnuSelf+1.2*m0rSelf[7]*mnuSelf+0.3513821107499669*uSelf[5]*cMSelf[7]*mnuSelf+0.1788854381999831*uSelf[4]*cMSelf[7]*mnuSelf+0.2*uSelf[0]*cMSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.5528571428571428*cMSelf[3]*uSelf[6]*mnuSelf+0.5528571428571428*uSelf[3]*cMSelf[6]*mnuSelf+0.2*cMSelf[1]*uSelf[5]*mnuSelf+0.2*uSelf[1]*cMSelf[5]*mnuSelf+0.3928571428571429*cMSelf[1]*uSelf[4]*mnuSelf+0.3928571428571429*uSelf[1]*cMSelf[4]*mnuSelf+0.4024922359499621*cMSelf[2]*uSelf[3]*mnuSelf+0.4024922359499621*uSelf[2]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(33,17) = 0.3273268353539885*cMSelf[5]*uSelf[9]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[9]*mnuSelf+1.31746509848052*m0rSelf[9]*mnuSelf+0.3273268353539885*uSelf[5]*cMSelf[9]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[9]*mnuSelf-0.43915503282684*cESelf[9]*mnuSelf+0.175662013130736*cMSelf[3]*uSelf[8]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[8]*mnuSelf+0.5528571428571428*cMSelf[3]*uSelf[7]*mnuSelf+0.5528571428571428*uSelf[3]*cMSelf[7]*mnuSelf+0.1788854381999831*cMSelf[5]*uSelf[6]*mnuSelf+0.3513821107499669*cMSelf[4]*uSelf[6]*mnuSelf+0.2*cMSelf[0]*uSelf[6]*mnuSelf+1.2*m0rSelf[6]*mnuSelf+0.1788854381999831*uSelf[5]*cMSelf[6]*mnuSelf+0.3513821107499669*uSelf[4]*cMSelf[6]*mnuSelf+0.2*uSelf[0]*cMSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.3928571428571429*cMSelf[2]*uSelf[5]*mnuSelf+0.3928571428571429*uSelf[2]*cMSelf[5]*mnuSelf+0.2*cMSelf[2]*uSelf[4]*mnuSelf+0.2*uSelf[2]*cMSelf[4]*mnuSelf+0.4024922359499621*cMSelf[1]*uSelf[3]*mnuSelf+0.4024922359499621*uSelf[1]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(33,18) = 0.1928571428571429*cMSelf[3]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[3]*cMSelf[8]*mnuSelf+0.175662013130736*cMSelf[3]*uSelf[7]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[5]*uSelf[6]*mnuSelf+0.1402829294374237*cMSelf[4]*uSelf[6]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[6]*mnuSelf+1.31746509848052*m0rSelf[6]*mnuSelf+0.1963961012123931*uSelf[5]*cMSelf[6]*mnuSelf+0.1402829294374237*uSelf[4]*cMSelf[6]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[6]*mnuSelf-0.43915503282684*cESelf[6]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[4]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[4]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[3]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(33,19) = 0.1928571428571429*cMSelf[3]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[3]*cMSelf[9]*mnuSelf+0.1402829294374237*cMSelf[5]*uSelf[7]*mnuSelf+0.1963961012123931*cMSelf[4]*uSelf[7]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[7]*mnuSelf+1.31746509848052*m0rSelf[7]*mnuSelf+0.1402829294374237*uSelf[5]*cMSelf[7]*mnuSelf+0.1963961012123931*uSelf[4]*cMSelf[7]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[7]*mnuSelf-0.43915503282684*cESelf[7]*mnuSelf+0.175662013130736*cMSelf[3]*uSelf[6]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[6]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[5]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[5]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[3]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(34,10) = 0.149071198499986*cMSelf[8]*uSelf[8]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[8]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[6]*uSelf[6]*mnuSelf+0.2500000000000001*cMSelf[2]*uSelf[6]*mnuSelf+0.2500000000000001*uSelf[2]*cMSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[0]*uSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf+0.25*uSelf[0]*cMSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(34,11) = 0.3273268353539885*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[8]*mnuSelf+1.317465098480519*m0rSelf[8]*mnuSelf+0.3273268353539885*uSelf[4]*cMSelf[8]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+0.223606797749979*cMSelf[5]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[5]*cMSelf[7]*mnuSelf+0.3928571428571429*cMSelf[3]*uSelf[6]*mnuSelf+0.3928571428571429*uSelf[3]*cMSelf[6]*mnuSelf+0.3928571428571428*cMSelf[1]*uSelf[4]*mnuSelf+0.3928571428571428*uSelf[1]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(34,12) = 0.2195775164134199*cMSelf[3]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[8]*mnuSelf+0.2*cMSelf[3]*uSelf[7]*mnuSelf+0.2*uSelf[3]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[5]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[6]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf+0.223606797749979*uSelf[5]*cMSelf[6]*mnuSelf+0.159719141249985*uSelf[4]*cMSelf[6]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf+0.25*cMSelf[2]*uSelf[4]*mnuSelf+0.25*uSelf[2]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(34,13) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999831*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999831*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(34,14) = 0.25*cMSelf[9]*uSelf[9]*mnuSelf+0.3452380952380952*cMSelf[8]*uSelf[8]*mnuSelf+0.1402829294374236*cMSelf[1]*uSelf[8]*mnuSelf+0.1402829294374236*uSelf[1]*cMSelf[8]*mnuSelf+0.3928571428571428*cMSelf[7]*uSelf[7]*mnuSelf+0.3520408163265306*cMSelf[6]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[2]*uSelf[6]*mnuSelf+0.159719141249985*uSelf[2]*cMSelf[6]*mnuSelf+0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.3520408163265306*cMSelf[4]*uSelf[4]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.3928571428571428*cMSelf[3]*uSelf[3]*mnuSelf+0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.3928571428571428*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(34,16) = 0.2195775164134199*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[9]*mnuSelf+0.1402829294374237*cMSelf[3]*uSelf[8]*mnuSelf+0.1402829294374237*uSelf[3]*cMSelf[8]*mnuSelf+0.3513821107499669*cMSelf[3]*uSelf[7]*mnuSelf+0.3513821107499669*uSelf[3]*cMSelf[7]*mnuSelf+0.1428571428571428*cMSelf[5]*uSelf[6]*mnuSelf+0.3520408163265306*cMSelf[4]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[6]*mnuSelf+0.9583148474999099*m0rSelf[6]*mnuSelf+0.1428571428571428*uSelf[5]*cMSelf[6]*mnuSelf+0.3520408163265306*uSelf[4]*cMSelf[6]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.159719141249985*cMSelf[2]*uSelf[4]*mnuSelf+0.159719141249985*uSelf[2]*cMSelf[4]*mnuSelf+0.3928571428571429*cMSelf[1]*uSelf[3]*mnuSelf+0.3928571428571429*uSelf[1]*cMSelf[3]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(34,17) = 0.1963961012123931*cMSelf[3]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[9]*mnuSelf+0.1428571428571428*cMSelf[5]*uSelf[7]*mnuSelf+0.2*cMSelf[4]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf+0.1428571428571428*uSelf[5]*cMSelf[7]*mnuSelf+0.2*uSelf[4]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[6]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[5]*mnuSelf+0.2*cMSelf[2]*uSelf[3]*mnuSelf+0.2*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(34,18) = 0.2817460317460317*cMSelf[4]*uSelf[8]*mnuSelf+0.149071198499986*cMSelf[0]*uSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[8]*mnuSelf+0.2817460317460317*uSelf[4]*cMSelf[8]*mnuSelf+0.149071198499986*uSelf[0]*cMSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+0.2195775164134199*cMSelf[5]*uSelf[7]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[7]*mnuSelf+0.3273268353539885*cMSelf[3]*uSelf[6]*mnuSelf+0.3273268353539885*uSelf[3]*cMSelf[6]*mnuSelf+0.3273268353539885*cMSelf[1]*uSelf[4]*mnuSelf+0.3273268353539885*uSelf[1]*cMSelf[4]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[3]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[3]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[1]*mnuSelf+1.317465098480519*m0rSelf[1]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(35,10) = 0.149071198499986*cMSelf[9]*uSelf[9]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[7]*uSelf[7]*mnuSelf+0.2500000000000001*cMSelf[1]*uSelf[7]*mnuSelf+0.2500000000000001*uSelf[1]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[0]*uSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf+0.25*uSelf[0]*cMSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(35,11) = 0.2195775164134199*cMSelf[3]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[4]*uSelf[7]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf+0.159719141249985*uSelf[5]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[4]*cMSelf[7]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf+0.2*cMSelf[3]*uSelf[6]*mnuSelf+0.2*uSelf[3]*cMSelf[6]*mnuSelf+0.25*cMSelf[1]*uSelf[5]*mnuSelf+0.25*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(35,12) = 0.3273268353539885*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[9]*mnuSelf+1.317465098480519*m0rSelf[9]*mnuSelf+0.3273268353539885*uSelf[5]*cMSelf[9]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+0.3928571428571429*cMSelf[3]*uSelf[7]*mnuSelf+0.3928571428571429*uSelf[3]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[4]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[4]*cMSelf[6]*mnuSelf+0.3928571428571428*cMSelf[2]*uSelf[5]*mnuSelf+0.3928571428571428*uSelf[2]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(35,13) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999831*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999831*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(35,15) = 0.3452380952380952*cMSelf[9]*uSelf[9]*mnuSelf+0.1402829294374236*cMSelf[2]*uSelf[9]*mnuSelf+0.1402829294374236*uSelf[2]*cMSelf[9]*mnuSelf+0.25*cMSelf[8]*uSelf[8]*mnuSelf+0.3520408163265306*cMSelf[7]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[1]*uSelf[7]*mnuSelf+0.159719141249985*uSelf[1]*cMSelf[7]*mnuSelf+0.3928571428571428*cMSelf[6]*uSelf[6]*mnuSelf+0.3520408163265306*cMSelf[5]*uSelf[5]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.3928571428571428*cMSelf[3]*uSelf[3]*mnuSelf+0.3928571428571428*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(35,16) = 0.1963961012123931*cMSelf[3]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[8]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[7]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[7]*mnuSelf+0.2*cMSelf[5]*uSelf[6]*mnuSelf+0.1428571428571428*cMSelf[4]*uSelf[6]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf+0.2*uSelf[5]*cMSelf[6]*mnuSelf+0.1428571428571428*uSelf[4]*cMSelf[6]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[4]*mnuSelf+0.2*cMSelf[1]*uSelf[3]*mnuSelf+0.2*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(35,17) = 0.1402829294374237*cMSelf[3]*uSelf[9]*mnuSelf+0.1402829294374237*uSelf[3]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[8]*mnuSelf+0.3520408163265306*cMSelf[5]*uSelf[7]*mnuSelf+0.1428571428571428*cMSelf[4]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[7]*mnuSelf+0.9583148474999099*m0rSelf[7]*mnuSelf+0.3520408163265306*uSelf[5]*cMSelf[7]*mnuSelf+0.1428571428571428*uSelf[4]*cMSelf[7]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+0.3513821107499669*cMSelf[3]*uSelf[6]*mnuSelf+0.3513821107499669*uSelf[3]*cMSelf[6]*mnuSelf+0.159719141249985*cMSelf[1]*uSelf[5]*mnuSelf+0.159719141249985*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.3928571428571429*cMSelf[2]*uSelf[3]*mnuSelf+0.3928571428571429*uSelf[2]*cMSelf[3]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(35,19) = 0.2817460317460317*cMSelf[5]*uSelf[9]*mnuSelf+0.149071198499986*cMSelf[0]*uSelf[9]*mnuSelf+0.8944271909999159*m0rSelf[9]*mnuSelf+0.2817460317460317*uSelf[5]*cMSelf[9]*mnuSelf+0.149071198499986*uSelf[0]*cMSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+0.3273268353539885*cMSelf[3]*uSelf[7]*mnuSelf+0.3273268353539885*uSelf[3]*cMSelf[7]*mnuSelf+0.2195775164134199*cMSelf[4]*uSelf[6]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[6]*mnuSelf+0.3273268353539885*cMSelf[2]*uSelf[5]*mnuSelf+0.3273268353539885*uSelf[2]*cMSelf[5]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[3]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[3]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[2]*mnuSelf+1.317465098480519*m0rSelf[2]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(36,10) = 0.2195775164134199*cMSelf[3]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[8]*mnuSelf+0.2*cMSelf[3]*uSelf[7]*mnuSelf+0.2*uSelf[3]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[5]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[6]*mnuSelf+0.25*cMSelf[0]*uSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf+0.223606797749979*uSelf[5]*cMSelf[6]*mnuSelf+0.159719141249985*uSelf[4]*cMSelf[6]*mnuSelf+0.25*uSelf[0]*cMSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf+0.2500000000000001*cMSelf[2]*uSelf[4]*mnuSelf+0.2500000000000001*uSelf[2]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(36,11) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999832*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999832*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(36,12) = 0.149071198499986*cMSelf[8]*uSelf[8]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[8]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[6]*uSelf[6]*mnuSelf+0.25*cMSelf[2]*uSelf[6]*mnuSelf+0.25*uSelf[2]*cMSelf[6]*mnuSelf+0.159719141249985*cMSelf[4]*uSelf[4]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(36,13) = 0.175662013130736*cMSelf[3]*uSelf[9]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[9]*mnuSelf+0.3273268353539885*cMSelf[4]*uSelf[8]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[8]*mnuSelf+1.31746509848052*m0rSelf[8]*mnuSelf+0.3273268353539885*uSelf[4]*cMSelf[8]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[8]*mnuSelf-0.43915503282684*cESelf[8]*mnuSelf+0.3513821107499669*cMSelf[5]*uSelf[7]*mnuSelf+0.1788854381999831*cMSelf[4]*uSelf[7]*mnuSelf+0.2*cMSelf[0]*uSelf[7]*mnuSelf+1.2*m0rSelf[7]*mnuSelf+0.3513821107499669*uSelf[5]*cMSelf[7]*mnuSelf+0.1788854381999831*uSelf[4]*cMSelf[7]*mnuSelf+0.2*uSelf[0]*cMSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.5528571428571428*cMSelf[3]*uSelf[6]*mnuSelf+0.5528571428571428*uSelf[3]*cMSelf[6]*mnuSelf+0.2*cMSelf[1]*uSelf[5]*mnuSelf+0.2*uSelf[1]*cMSelf[5]*mnuSelf+0.3928571428571429*cMSelf[1]*uSelf[4]*mnuSelf+0.3928571428571429*uSelf[1]*cMSelf[4]*mnuSelf+0.4024922359499621*cMSelf[2]*uSelf[3]*mnuSelf+0.4024922359499621*uSelf[2]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(36,14) = 0.2195775164134199*cMSelf[5]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[9]*mnuSelf+0.1402829294374237*cMSelf[3]*uSelf[8]*mnuSelf+0.1402829294374237*uSelf[3]*cMSelf[8]*mnuSelf+0.3513821107499669*cMSelf[3]*uSelf[7]*mnuSelf+0.3513821107499669*uSelf[3]*cMSelf[7]*mnuSelf+0.1428571428571428*cMSelf[5]*uSelf[6]*mnuSelf+0.3520408163265306*cMSelf[4]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[6]*mnuSelf+0.9583148474999099*m0rSelf[6]*mnuSelf+0.1428571428571428*uSelf[5]*cMSelf[6]*mnuSelf+0.3520408163265306*uSelf[4]*cMSelf[6]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[5]*mnuSelf+0.159719141249985*cMSelf[2]*uSelf[4]*mnuSelf+0.159719141249985*uSelf[2]*cMSelf[4]*mnuSelf+0.3928571428571429*cMSelf[1]*uSelf[3]*mnuSelf+0.3928571428571429*uSelf[1]*cMSelf[3]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(36,15) = 0.1963961012123931*cMSelf[3]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[8]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[7]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[7]*mnuSelf+0.2*cMSelf[5]*uSelf[6]*mnuSelf+0.1428571428571428*cMSelf[4]*uSelf[6]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf+0.2*uSelf[5]*cMSelf[6]*mnuSelf+0.1428571428571428*uSelf[4]*cMSelf[6]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[4]*mnuSelf+0.2*cMSelf[1]*uSelf[3]*mnuSelf+0.2*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(36,16) = 0.3833333333333334*cMSelf[9]*uSelf[9]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[9]*mnuSelf+0.3452380952380952*cMSelf[8]*uSelf[8]*mnuSelf+0.1402829294374236*cMSelf[1]*uSelf[8]*mnuSelf+0.1402829294374236*uSelf[1]*cMSelf[8]*mnuSelf+0.5357142857142857*cMSelf[7]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[7]*mnuSelf+0.5520408163265306*cMSelf[6]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[2]*uSelf[6]*mnuSelf+0.159719141249985*uSelf[2]*cMSelf[6]*mnuSelf+0.3928571428571428*cMSelf[5]*uSelf[5]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.3520408163265306*cMSelf[4]*uSelf[4]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.5928571428571429*cMSelf[3]*uSelf[3]*mnuSelf+0.45*cMSelf[2]*uSelf[2]*mnuSelf+0.3928571428571428*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(36,17) = 0.175662013130736*cMSelf[7]*uSelf[9]*mnuSelf+0.175662013130736*uSelf[7]*cMSelf[9]*mnuSelf+0.175662013130736*cMSelf[6]*uSelf[8]*mnuSelf+0.175662013130736*uSelf[6]*cMSelf[8]*mnuSelf+0.16*cMSelf[6]*uSelf[7]*mnuSelf+0.1788854381999832*cMSelf[2]*uSelf[7]*mnuSelf+0.16*uSelf[6]*cMSelf[7]*mnuSelf+0.1788854381999832*uSelf[2]*cMSelf[7]*mnuSelf+0.1788854381999832*cMSelf[1]*uSelf[6]*mnuSelf+0.1788854381999832*uSelf[1]*cMSelf[6]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[5]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[5]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[4]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[4]*mnuSelf+0.2*cMSelf[0]*uSelf[3]*mnuSelf+1.2*m0rSelf[3]*mnuSelf+0.2*uSelf[0]*cMSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf+0.2*cMSelf[1]*uSelf[2]*mnuSelf+0.2*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(36,18) = 0.1928571428571429*cMSelf[7]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[7]*cMSelf[9]*mnuSelf+0.1928571428571429*cMSelf[6]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[6]*cMSelf[8]*mnuSelf+0.175662013130736*cMSelf[6]*uSelf[7]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[7]*mnuSelf+0.175662013130736*uSelf[6]*cMSelf[7]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[6]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[6]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[5]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[5]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[4]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[4]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[3]*mnuSelf+1.31746509848052*m0rSelf[3]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf+0.21957751641342*cMSelf[1]*uSelf[2]*mnuSelf+0.21957751641342*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(37,10) = 0.2195775164134199*cMSelf[3]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[4]*uSelf[7]*mnuSelf+0.25*cMSelf[0]*uSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf+0.159719141249985*uSelf[5]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[4]*cMSelf[7]*mnuSelf+0.25*uSelf[0]*cMSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf+0.2*cMSelf[3]*uSelf[6]*mnuSelf+0.2*uSelf[3]*cMSelf[6]*mnuSelf+0.2500000000000001*cMSelf[1]*uSelf[5]*mnuSelf+0.2500000000000001*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[3]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(37,11) = 0.149071198499986*cMSelf[9]*uSelf[9]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[9]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[9]*mnuSelf+0.159719141249985*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[1]*uSelf[7]*mnuSelf+0.25*uSelf[1]*cMSelf[7]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[6]*mnuSelf+0.159719141249985*cMSelf[5]*uSelf[5]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf+0.223606797749979*cMSelf[3]*uSelf[3]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(37,12) = 0.1963961012123931*cMSelf[7]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[7]*cMSelf[9]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[6]*cMSelf[8]*mnuSelf+0.1788854381999832*cMSelf[6]*uSelf[7]*mnuSelf+0.2*cMSelf[2]*uSelf[7]*mnuSelf+0.1788854381999832*uSelf[6]*cMSelf[7]*mnuSelf+0.2*uSelf[2]*cMSelf[7]*mnuSelf+0.2*cMSelf[1]*uSelf[6]*mnuSelf+0.2*uSelf[1]*cMSelf[6]*mnuSelf+0.2*cMSelf[3]*uSelf[5]*mnuSelf+0.2*uSelf[3]*cMSelf[5]*mnuSelf+0.2*cMSelf[3]*uSelf[4]*mnuSelf+0.2*uSelf[3]*cMSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[2]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(37,13) = 0.3273268353539885*cMSelf[5]*uSelf[9]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[9]*mnuSelf+1.31746509848052*m0rSelf[9]*mnuSelf+0.3273268353539885*uSelf[5]*cMSelf[9]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[9]*mnuSelf-0.43915503282684*cESelf[9]*mnuSelf+0.175662013130736*cMSelf[3]*uSelf[8]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[8]*mnuSelf+0.5528571428571428*cMSelf[3]*uSelf[7]*mnuSelf+0.5528571428571428*uSelf[3]*cMSelf[7]*mnuSelf+0.1788854381999831*cMSelf[5]*uSelf[6]*mnuSelf+0.3513821107499669*cMSelf[4]*uSelf[6]*mnuSelf+0.2*cMSelf[0]*uSelf[6]*mnuSelf+1.2*m0rSelf[6]*mnuSelf+0.1788854381999831*uSelf[5]*cMSelf[6]*mnuSelf+0.3513821107499669*uSelf[4]*cMSelf[6]*mnuSelf+0.2*uSelf[0]*cMSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.3928571428571429*cMSelf[2]*uSelf[5]*mnuSelf+0.3928571428571429*uSelf[2]*cMSelf[5]*mnuSelf+0.2*cMSelf[2]*uSelf[4]*mnuSelf+0.2*uSelf[2]*cMSelf[4]*mnuSelf+0.4024922359499621*cMSelf[1]*uSelf[3]*mnuSelf+0.4024922359499621*uSelf[1]*cMSelf[3]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(37,14) = 0.1963961012123931*cMSelf[3]*uSelf[9]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[9]*mnuSelf+0.1428571428571428*cMSelf[5]*uSelf[7]*mnuSelf+0.2*cMSelf[4]*uSelf[7]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf+0.1428571428571428*uSelf[5]*cMSelf[7]*mnuSelf+0.2*uSelf[4]*cMSelf[7]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[6]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[6]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[5]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[5]*mnuSelf+0.2*cMSelf[2]*uSelf[3]*mnuSelf+0.2*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(37,15) = 0.1402829294374237*cMSelf[3]*uSelf[9]*mnuSelf+0.1402829294374237*uSelf[3]*cMSelf[9]*mnuSelf+0.2195775164134199*cMSelf[4]*uSelf[8]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[8]*mnuSelf+0.3520408163265306*cMSelf[5]*uSelf[7]*mnuSelf+0.1428571428571428*cMSelf[4]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[7]*mnuSelf+0.9583148474999099*m0rSelf[7]*mnuSelf+0.3520408163265306*uSelf[5]*cMSelf[7]*mnuSelf+0.1428571428571428*uSelf[4]*cMSelf[7]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+0.3513821107499669*cMSelf[3]*uSelf[6]*mnuSelf+0.3513821107499669*uSelf[3]*cMSelf[6]*mnuSelf+0.159719141249985*cMSelf[1]*uSelf[5]*mnuSelf+0.159719141249985*uSelf[1]*cMSelf[5]*mnuSelf+0.223606797749979*cMSelf[1]*uSelf[4]*mnuSelf+0.223606797749979*uSelf[1]*cMSelf[4]*mnuSelf+0.3928571428571429*cMSelf[2]*uSelf[3]*mnuSelf+0.3928571428571429*uSelf[2]*cMSelf[3]*mnuSelf+0.2500000000000001*cMSelf[0]*uSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf+0.2500000000000001*uSelf[0]*cMSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(37,16) = 0.175662013130736*cMSelf[7]*uSelf[9]*mnuSelf+0.175662013130736*uSelf[7]*cMSelf[9]*mnuSelf+0.175662013130736*cMSelf[6]*uSelf[8]*mnuSelf+0.175662013130736*uSelf[6]*cMSelf[8]*mnuSelf+0.16*cMSelf[6]*uSelf[7]*mnuSelf+0.1788854381999832*cMSelf[2]*uSelf[7]*mnuSelf+0.16*uSelf[6]*cMSelf[7]*mnuSelf+0.1788854381999832*uSelf[2]*cMSelf[7]*mnuSelf+0.1788854381999832*cMSelf[1]*uSelf[6]*mnuSelf+0.1788854381999832*uSelf[1]*cMSelf[6]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[5]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[5]*mnuSelf+0.1788854381999831*cMSelf[3]*uSelf[4]*mnuSelf+0.1788854381999831*uSelf[3]*cMSelf[4]*mnuSelf+0.2*cMSelf[0]*uSelf[3]*mnuSelf+1.2*m0rSelf[3]*mnuSelf+0.2*uSelf[0]*cMSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf+0.2*cMSelf[1]*uSelf[2]*mnuSelf+0.2*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(37,17) = 0.3452380952380952*cMSelf[9]*uSelf[9]*mnuSelf+0.1402829294374236*cMSelf[2]*uSelf[9]*mnuSelf+0.1402829294374236*uSelf[2]*cMSelf[9]*mnuSelf+0.3833333333333334*cMSelf[8]*uSelf[8]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[8]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[8]*mnuSelf+0.5520408163265306*cMSelf[7]*uSelf[7]*mnuSelf+0.159719141249985*cMSelf[1]*uSelf[7]*mnuSelf+0.159719141249985*uSelf[1]*cMSelf[7]*mnuSelf+0.5357142857142857*cMSelf[6]*uSelf[6]*mnuSelf+0.223606797749979*cMSelf[2]*uSelf[6]*mnuSelf+0.223606797749979*uSelf[2]*cMSelf[6]*mnuSelf+0.3520408163265306*cMSelf[5]*uSelf[5]*mnuSelf+0.159719141249985*cMSelf[0]*uSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf+0.159719141249985*uSelf[0]*cMSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.3928571428571428*cMSelf[4]*uSelf[4]*mnuSelf+0.223606797749979*cMSelf[0]*uSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf+0.223606797749979*uSelf[0]*cMSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5928571428571429*cMSelf[3]*uSelf[3]*mnuSelf+0.3928571428571428*cMSelf[2]*uSelf[2]*mnuSelf+0.45*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(37,19) = 0.1928571428571429*cMSelf[7]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[7]*cMSelf[9]*mnuSelf+0.1928571428571429*cMSelf[6]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[6]*cMSelf[8]*mnuSelf+0.175662013130736*cMSelf[6]*uSelf[7]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[7]*mnuSelf+0.175662013130736*uSelf[6]*cMSelf[7]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[6]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[6]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[5]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[5]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[4]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[4]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[3]*mnuSelf+1.31746509848052*m0rSelf[3]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf+0.21957751641342*cMSelf[1]*uSelf[2]*mnuSelf+0.21957751641342*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(38,10) = 0.149071198499986*cMSelf[4]*uSelf[8]*mnuSelf+0.25*cMSelf[0]*uSelf[8]*mnuSelf+1.5*m0rSelf[8]*mnuSelf+0.149071198499986*uSelf[4]*cMSelf[8]*mnuSelf+0.25*uSelf[0]*cMSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf+0.2195775164134199*cMSelf[3]*uSelf[6]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[6]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[4]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[4]*mnuSelf; 
  data->AEM_S(38,11) = 0.1309307341415954*cMSelf[8]*uSelf[8]*mnuSelf+0.1928571428571429*cMSelf[1]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[1]*cMSelf[8]*mnuSelf+0.1963961012123931*cMSelf[7]*uSelf[7]*mnuSelf+0.1402829294374236*cMSelf[6]*uSelf[6]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[6]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[6]*mnuSelf+0.1402829294374236*cMSelf[4]*uSelf[4]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[4]*mnuSelf+1.317465098480519*m0rSelf[4]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[3]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(38,13) = 0.1928571428571429*cMSelf[3]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[3]*cMSelf[8]*mnuSelf+0.175662013130736*cMSelf[3]*uSelf[7]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[5]*uSelf[6]*mnuSelf+0.1402829294374237*cMSelf[4]*uSelf[6]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[6]*mnuSelf+1.31746509848052*m0rSelf[6]*mnuSelf+0.1963961012123931*uSelf[5]*cMSelf[6]*mnuSelf+0.1402829294374237*uSelf[4]*cMSelf[6]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[6]*mnuSelf-0.43915503282684*cESelf[6]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[4]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[4]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[3]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[3]*mnuSelf; 
  data->AEM_S(38,14) = 0.2817460317460317*cMSelf[4]*uSelf[8]*mnuSelf+0.149071198499986*cMSelf[0]*uSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[8]*mnuSelf+0.2817460317460317*uSelf[4]*cMSelf[8]*mnuSelf+0.149071198499986*uSelf[0]*cMSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+0.2195775164134199*cMSelf[5]*uSelf[7]*mnuSelf+0.2195775164134199*uSelf[5]*cMSelf[7]*mnuSelf+0.3273268353539885*cMSelf[3]*uSelf[6]*mnuSelf+0.3273268353539885*uSelf[3]*cMSelf[6]*mnuSelf+0.3273268353539885*cMSelf[1]*uSelf[4]*mnuSelf+0.3273268353539885*uSelf[1]*cMSelf[4]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[3]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[3]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[1]*mnuSelf+1.317465098480519*m0rSelf[1]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(38,16) = 0.1928571428571429*cMSelf[7]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[7]*cMSelf[9]*mnuSelf+0.1928571428571429*cMSelf[6]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[6]*cMSelf[8]*mnuSelf+0.175662013130736*cMSelf[6]*uSelf[7]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[7]*mnuSelf+0.175662013130736*uSelf[6]*cMSelf[7]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[6]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[6]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[5]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[5]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[4]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[4]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[3]*mnuSelf+1.31746509848052*m0rSelf[3]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf+0.21957751641342*cMSelf[1]*uSelf[2]*mnuSelf+0.21957751641342*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(38,18) = 0.25*cMSelf[9]*uSelf[9]*mnuSelf+0.3388888888888889*cMSelf[8]*uSelf[8]*mnuSelf+0.1309307341415954*cMSelf[1]*uSelf[8]*mnuSelf+0.1309307341415954*uSelf[1]*cMSelf[8]*mnuSelf+0.3833333333333334*cMSelf[7]*uSelf[7]*mnuSelf+0.3452380952380952*cMSelf[6]*uSelf[6]*mnuSelf+0.149071198499986*cMSelf[2]*uSelf[6]*mnuSelf+0.149071198499986*uSelf[2]*cMSelf[6]*mnuSelf+0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.3452380952380952*cMSelf[4]*uSelf[4]*mnuSelf+0.149071198499986*cMSelf[0]*uSelf[4]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf+0.149071198499986*uSelf[0]*cMSelf[4]*mnuSelf-0.2981423969999719*cESelf[4]*mnuSelf+0.3833333333333334*cMSelf[3]*uSelf[3]*mnuSelf+0.25*cMSelf[2]*uSelf[2]*mnuSelf+0.3833333333333334*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(39,10) = 0.149071198499986*cMSelf[5]*uSelf[9]*mnuSelf+0.25*cMSelf[0]*uSelf[9]*mnuSelf+1.5*m0rSelf[9]*mnuSelf+0.149071198499986*uSelf[5]*cMSelf[9]*mnuSelf+0.25*uSelf[0]*cMSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf+0.2195775164134199*cMSelf[3]*uSelf[7]*mnuSelf+0.2195775164134199*uSelf[3]*cMSelf[7]*mnuSelf+0.2195775164134199*cMSelf[2]*uSelf[5]*mnuSelf+0.2195775164134199*uSelf[2]*cMSelf[5]*mnuSelf; 
  data->AEM_S(39,12) = 0.1309307341415954*cMSelf[9]*uSelf[9]*mnuSelf+0.1928571428571429*cMSelf[2]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[2]*cMSelf[9]*mnuSelf+0.1402829294374236*cMSelf[7]*uSelf[7]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[7]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[6]*uSelf[6]*mnuSelf+0.1402829294374236*cMSelf[5]*uSelf[5]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[5]*mnuSelf+1.317465098480519*m0rSelf[5]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[3]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(39,13) = 0.1928571428571429*cMSelf[3]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[3]*cMSelf[9]*mnuSelf+0.1402829294374237*cMSelf[5]*uSelf[7]*mnuSelf+0.1963961012123931*cMSelf[4]*uSelf[7]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[7]*mnuSelf+1.31746509848052*m0rSelf[7]*mnuSelf+0.1402829294374237*uSelf[5]*cMSelf[7]*mnuSelf+0.1963961012123931*uSelf[4]*cMSelf[7]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[7]*mnuSelf-0.43915503282684*cESelf[7]*mnuSelf+0.175662013130736*cMSelf[3]*uSelf[6]*mnuSelf+0.175662013130736*uSelf[3]*cMSelf[6]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[5]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[5]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[3]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[3]*mnuSelf; 
  data->AEM_S(39,15) = 0.2817460317460317*cMSelf[5]*uSelf[9]*mnuSelf+0.149071198499986*cMSelf[0]*uSelf[9]*mnuSelf+0.8944271909999159*m0rSelf[9]*mnuSelf+0.2817460317460317*uSelf[5]*cMSelf[9]*mnuSelf+0.149071198499986*uSelf[0]*cMSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+0.3273268353539885*cMSelf[3]*uSelf[7]*mnuSelf+0.3273268353539885*uSelf[3]*cMSelf[7]*mnuSelf+0.2195775164134199*cMSelf[4]*uSelf[6]*mnuSelf+0.2195775164134199*uSelf[4]*cMSelf[6]*mnuSelf+0.3273268353539885*cMSelf[2]*uSelf[5]*mnuSelf+0.3273268353539885*uSelf[2]*cMSelf[5]*mnuSelf+0.2195775164134199*cMSelf[1]*uSelf[3]*mnuSelf+0.2195775164134199*uSelf[1]*cMSelf[3]*mnuSelf+0.2195775164134199*cMSelf[0]*uSelf[2]*mnuSelf+1.317465098480519*m0rSelf[2]*mnuSelf+0.2195775164134199*uSelf[0]*cMSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(39,17) = 0.1928571428571429*cMSelf[7]*uSelf[9]*mnuSelf+0.1928571428571429*uSelf[7]*cMSelf[9]*mnuSelf+0.1928571428571429*cMSelf[6]*uSelf[8]*mnuSelf+0.1928571428571429*uSelf[6]*cMSelf[8]*mnuSelf+0.175662013130736*cMSelf[6]*uSelf[7]*mnuSelf+0.1963961012123931*cMSelf[2]*uSelf[7]*mnuSelf+0.175662013130736*uSelf[6]*cMSelf[7]*mnuSelf+0.1963961012123931*uSelf[2]*cMSelf[7]*mnuSelf+0.1963961012123931*cMSelf[1]*uSelf[6]*mnuSelf+0.1963961012123931*uSelf[1]*cMSelf[6]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[5]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[5]*mnuSelf+0.1963961012123931*cMSelf[3]*uSelf[4]*mnuSelf+0.1963961012123931*uSelf[3]*cMSelf[4]*mnuSelf+0.21957751641342*cMSelf[0]*uSelf[3]*mnuSelf+1.31746509848052*m0rSelf[3]*mnuSelf+0.21957751641342*uSelf[0]*cMSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf+0.21957751641342*cMSelf[1]*uSelf[2]*mnuSelf+0.21957751641342*uSelf[1]*cMSelf[2]*mnuSelf; 
  data->AEM_S(39,19) = 0.3388888888888889*cMSelf[9]*uSelf[9]*mnuSelf+0.1309307341415954*cMSelf[2]*uSelf[9]*mnuSelf+0.1309307341415954*uSelf[2]*cMSelf[9]*mnuSelf+0.25*cMSelf[8]*uSelf[8]*mnuSelf+0.3452380952380952*cMSelf[7]*uSelf[7]*mnuSelf+0.149071198499986*cMSelf[1]*uSelf[7]*mnuSelf+0.149071198499986*uSelf[1]*cMSelf[7]*mnuSelf+0.3833333333333334*cMSelf[6]*uSelf[6]*mnuSelf+0.3452380952380952*cMSelf[5]*uSelf[5]*mnuSelf+0.149071198499986*cMSelf[0]*uSelf[5]*mnuSelf+0.8944271909999159*m0rSelf[5]*mnuSelf+0.149071198499986*uSelf[0]*cMSelf[5]*mnuSelf-0.2981423969999719*cESelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.3833333333333334*cMSelf[3]*uSelf[3]*mnuSelf+0.3833333333333334*cMSelf[2]*uSelf[2]*mnuSelf+0.25*cMSelf[1]*uSelf[1]*mnuSelf+0.25*cMSelf[0]*uSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(30,30) = (-0.25*cMOther[9]*uOther[9]*mnuOther)-0.25*cMOther[8]*uOther[8]*mnuOther-0.25*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[6]*uOther[6]*mnuOther-0.25*cMOther[5]*uOther[5]*mnuOther-0.25*cMOther[4]*uOther[4]*mnuOther-0.25*cMOther[3]*uOther[3]*mnuOther-0.25*cMOther[2]*uOther[2]*mnuOther-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(30,31) = (-0.2195775164134199*cMOther[4]*uOther[8]*mnuOther)-0.2195775164134199*uOther[4]*cMOther[8]*mnuOther-0.2500000000000001*cMOther[5]*uOther[7]*mnuOther-0.2500000000000001*uOther[5]*cMOther[7]*mnuOther-0.223606797749979*cMOther[3]*uOther[6]*mnuOther-0.223606797749979*uOther[3]*cMOther[6]*mnuOther-0.223606797749979*cMOther[1]*uOther[4]*mnuOther-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.25*cMOther[2]*uOther[3]*mnuOther-0.25*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(30,32) = (-0.2195775164134199*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*uOther[5]*cMOther[9]*mnuOther-0.223606797749979*cMOther[3]*uOther[7]*mnuOther-0.223606797749979*uOther[3]*cMOther[7]*mnuOther-0.2500000000000001*cMOther[4]*uOther[6]*mnuOther-0.2500000000000001*uOther[4]*cMOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[5]*mnuOther-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.25*cMOther[1]*uOther[3]*mnuOther-0.25*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(30,33) = (-0.2195775164134199*cMOther[7]*uOther[9]*mnuOther)-0.2195775164134199*uOther[7]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[6]*uOther[8]*mnuOther-0.2195775164134199*uOther[6]*cMOther[8]*mnuOther-0.2*cMOther[6]*uOther[7]*mnuOther-0.223606797749979*cMOther[2]*uOther[7]*mnuOther-0.2*uOther[6]*cMOther[7]*mnuOther-0.223606797749979*uOther[2]*cMOther[7]*mnuOther-0.223606797749979*cMOther[1]*uOther[6]*mnuOther-0.223606797749979*uOther[1]*cMOther[6]*mnuOther-0.223606797749979*cMOther[3]*uOther[5]*mnuOther-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(30,34) = (-0.149071198499986*cMOther[8]*uOther[8]*mnuOther)-0.2195775164134199*cMOther[1]*uOther[8]*mnuOther-0.2195775164134199*uOther[1]*cMOther[8]*mnuOther-0.223606797749979*cMOther[7]*uOther[7]*mnuOther-0.159719141249985*cMOther[6]*uOther[6]*mnuOther-0.2500000000000001*cMOther[2]*uOther[6]*mnuOther-0.2500000000000001*uOther[2]*cMOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[4]*mnuOther-0.25*cMOther[0]*uOther[4]*mnuOther-1.5*m0rOther[4]*mnuOther-0.25*uOther[0]*cMOther[4]*mnuOther+0.5*cEOther[4]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(30,35) = (-0.149071198499986*cMOther[9]*uOther[9]*mnuOther)-0.2195775164134199*cMOther[2]*uOther[9]*mnuOther-0.2195775164134199*uOther[2]*cMOther[9]*mnuOther-0.159719141249985*cMOther[7]*uOther[7]*mnuOther-0.2500000000000001*cMOther[1]*uOther[7]*mnuOther-0.2500000000000001*uOther[1]*cMOther[7]*mnuOther-0.223606797749979*cMOther[6]*uOther[6]*mnuOther-0.159719141249985*cMOther[5]*uOther[5]*mnuOther-0.25*cMOther[0]*uOther[5]*mnuOther-1.5*m0rOther[5]*mnuOther-0.25*uOther[0]*cMOther[5]*mnuOther+0.5*cEOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(30,36) = (-0.2195775164134199*cMOther[3]*uOther[8]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[8]*mnuOther-0.2*cMOther[3]*uOther[7]*mnuOther-0.2*uOther[3]*cMOther[7]*mnuOther-0.223606797749979*cMOther[5]*uOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[6]*mnuOther-0.25*cMOther[0]*uOther[6]*mnuOther-1.5*m0rOther[6]*mnuOther-0.223606797749979*uOther[5]*cMOther[6]*mnuOther-0.159719141249985*uOther[4]*cMOther[6]*mnuOther-0.25*uOther[0]*cMOther[6]*mnuOther+0.5*cEOther[6]*mnuOther-0.2500000000000001*cMOther[2]*uOther[4]*mnuOther-0.2500000000000001*uOther[2]*cMOther[4]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(30,37) = (-0.2195775164134199*cMOther[3]*uOther[9]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[9]*mnuOther-0.159719141249985*cMOther[5]*uOther[7]*mnuOther-0.223606797749979*cMOther[4]*uOther[7]*mnuOther-0.25*cMOther[0]*uOther[7]*mnuOther-1.5*m0rOther[7]*mnuOther-0.159719141249985*uOther[5]*cMOther[7]*mnuOther-0.223606797749979*uOther[4]*cMOther[7]*mnuOther-0.25*uOther[0]*cMOther[7]*mnuOther+0.5*cEOther[7]*mnuOther-0.2*cMOther[3]*uOther[6]*mnuOther-0.2*uOther[3]*cMOther[6]*mnuOther-0.2500000000000001*cMOther[1]*uOther[5]*mnuOther-0.2500000000000001*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(30,38) = (-0.149071198499986*cMOther[4]*uOther[8]*mnuOther)-0.25*cMOther[0]*uOther[8]*mnuOther-1.5*m0rOther[8]*mnuOther-0.149071198499986*uOther[4]*cMOther[8]*mnuOther-0.25*uOther[0]*cMOther[8]*mnuOther+0.5*cEOther[8]*mnuOther-0.2195775164134199*cMOther[3]*uOther[6]*mnuOther-0.2195775164134199*uOther[3]*cMOther[6]*mnuOther-0.2195775164134199*cMOther[1]*uOther[4]*mnuOther-0.2195775164134199*uOther[1]*cMOther[4]*mnuOther; 
  data->AEM_S(30,39) = (-0.149071198499986*cMOther[5]*uOther[9]*mnuOther)-0.25*cMOther[0]*uOther[9]*mnuOther-1.5*m0rOther[9]*mnuOther-0.149071198499986*uOther[5]*cMOther[9]*mnuOther-0.25*uOther[0]*cMOther[9]*mnuOther+0.5*cEOther[9]*mnuOther-0.2195775164134199*cMOther[3]*uOther[7]*mnuOther-0.2195775164134199*uOther[3]*cMOther[7]*mnuOther-0.2195775164134199*cMOther[2]*uOther[5]*mnuOther-0.2195775164134199*uOther[2]*cMOther[5]*mnuOther; 
  data->AEM_S(31,30) = (-0.2195775164134199*cMOther[4]*uOther[8]*mnuOther)-0.2195775164134199*uOther[4]*cMOther[8]*mnuOther-0.2500000000000001*cMOther[5]*uOther[7]*mnuOther-0.2500000000000001*uOther[5]*cMOther[7]*mnuOther-0.223606797749979*cMOther[3]*uOther[6]*mnuOther-0.223606797749979*uOther[3]*cMOther[6]*mnuOther-0.223606797749979*cMOther[1]*uOther[4]*mnuOther-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.25*cMOther[2]*uOther[3]*mnuOther-0.25*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(31,31) = (-0.25*cMOther[9]*uOther[9]*mnuOther)-0.3833333333333334*cMOther[8]*uOther[8]*mnuOther-0.1963961012123931*cMOther[1]*uOther[8]*mnuOther-0.1963961012123931*uOther[1]*cMOther[8]*mnuOther-0.45*cMOther[7]*uOther[7]*mnuOther-0.3928571428571428*cMOther[6]*uOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[6]*mnuOther-0.223606797749979*uOther[2]*cMOther[6]*mnuOther-0.25*cMOther[5]*uOther[5]*mnuOther-0.3928571428571428*cMOther[4]*uOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther-0.223606797749979*uOther[0]*cMOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.45*cMOther[3]*uOther[3]*mnuOther-0.25*cMOther[2]*uOther[2]*mnuOther-0.45*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(31,32) = (-0.2195775164134199*cMOther[7]*uOther[9]*mnuOther)-0.2195775164134199*uOther[7]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[6]*uOther[8]*mnuOther-0.2195775164134199*uOther[6]*cMOther[8]*mnuOther-0.2*cMOther[6]*uOther[7]*mnuOther-0.223606797749979*cMOther[2]*uOther[7]*mnuOther-0.2*uOther[6]*cMOther[7]*mnuOther-0.223606797749979*uOther[2]*cMOther[7]*mnuOther-0.223606797749979*cMOther[1]*uOther[6]*mnuOther-0.223606797749979*uOther[1]*cMOther[6]*mnuOther-0.223606797749979*cMOther[3]*uOther[5]*mnuOther-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(31,33) = (-0.2195775164134199*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*uOther[5]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[3]*uOther[8]*mnuOther-0.1963961012123931*uOther[3]*cMOther[8]*mnuOther-0.4024922359499621*cMOther[3]*uOther[7]*mnuOther-0.4024922359499621*uOther[3]*cMOther[7]*mnuOther-0.2*cMOther[5]*uOther[6]*mnuOther-0.3928571428571429*cMOther[4]*uOther[6]*mnuOther-0.223606797749979*cMOther[0]*uOther[6]*mnuOther-1.341640786499874*m0rOther[6]*mnuOther-0.2*uOther[5]*cMOther[6]*mnuOther-0.3928571428571429*uOther[4]*cMOther[6]*mnuOther-0.223606797749979*uOther[0]*cMOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[5]*mnuOther-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.223606797749979*cMOther[2]*uOther[4]*mnuOther-0.223606797749979*uOther[2]*cMOther[4]*mnuOther-0.45*cMOther[1]*uOther[3]*mnuOther-0.45*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(31,34) = (-0.3273268353539885*cMOther[4]*uOther[8]*mnuOther)-0.2195775164134199*cMOther[0]*uOther[8]*mnuOther-1.317465098480519*m0rOther[8]*mnuOther-0.3273268353539885*uOther[4]*cMOther[8]*mnuOther-0.2195775164134199*uOther[0]*cMOther[8]*mnuOther+0.4391550328268398*cEOther[8]*mnuOther-0.223606797749979*cMOther[5]*uOther[7]*mnuOther-0.223606797749979*uOther[5]*cMOther[7]*mnuOther-0.3928571428571429*cMOther[3]*uOther[6]*mnuOther-0.3928571428571429*uOther[3]*cMOther[6]*mnuOther-0.3928571428571428*cMOther[1]*uOther[4]*mnuOther-0.3928571428571428*uOther[1]*cMOther[4]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther-0.223606797749979*uOther[0]*cMOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(31,35) = (-0.2195775164134199*cMOther[3]*uOther[9]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[9]*mnuOther-0.159719141249985*cMOther[5]*uOther[7]*mnuOther-0.223606797749979*cMOther[4]*uOther[7]*mnuOther-0.2500000000000001*cMOther[0]*uOther[7]*mnuOther-1.5*m0rOther[7]*mnuOther-0.159719141249985*uOther[5]*cMOther[7]*mnuOther-0.223606797749979*uOther[4]*cMOther[7]*mnuOther-0.2500000000000001*uOther[0]*cMOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther-0.2*cMOther[3]*uOther[6]*mnuOther-0.2*uOther[3]*cMOther[6]*mnuOther-0.25*cMOther[1]*uOther[5]*mnuOther-0.25*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(31,36) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999832*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999832*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(31,37) = (-0.149071198499986*cMOther[9]*uOther[9]*mnuOther)-0.2195775164134199*cMOther[2]*uOther[9]*mnuOther-0.2195775164134199*uOther[2]*cMOther[9]*mnuOther-0.159719141249985*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[1]*uOther[7]*mnuOther-0.25*uOther[1]*cMOther[7]*mnuOther-0.223606797749979*cMOther[6]*uOther[6]*mnuOther-0.159719141249985*cMOther[5]*uOther[5]*mnuOther-0.2500000000000001*cMOther[0]*uOther[5]*mnuOther-1.5*m0rOther[5]*mnuOther-0.2500000000000001*uOther[0]*cMOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(31,38) = (-0.1309307341415954*cMOther[8]*uOther[8]*mnuOther)-0.1928571428571429*cMOther[1]*uOther[8]*mnuOther-0.1928571428571429*uOther[1]*cMOther[8]*mnuOther-0.1963961012123931*cMOther[7]*uOther[7]*mnuOther-0.1402829294374236*cMOther[6]*uOther[6]*mnuOther-0.2195775164134199*cMOther[2]*uOther[6]*mnuOther-0.2195775164134199*uOther[2]*cMOther[6]*mnuOther-0.1402829294374236*cMOther[4]*uOther[4]*mnuOther-0.2195775164134199*cMOther[0]*uOther[4]*mnuOther-1.317465098480519*m0rOther[4]*mnuOther-0.2195775164134199*uOther[0]*cMOther[4]*mnuOther+0.4391550328268398*cEOther[4]*mnuOther-0.1963961012123931*cMOther[3]*uOther[3]*mnuOther-0.1963961012123931*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(32,30) = (-0.2195775164134199*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*uOther[5]*cMOther[9]*mnuOther-0.223606797749979*cMOther[3]*uOther[7]*mnuOther-0.223606797749979*uOther[3]*cMOther[7]*mnuOther-0.2500000000000001*cMOther[4]*uOther[6]*mnuOther-0.2500000000000001*uOther[4]*cMOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[5]*mnuOther-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.25*cMOther[1]*uOther[3]*mnuOther-0.25*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(32,31) = (-0.2195775164134199*cMOther[7]*uOther[9]*mnuOther)-0.2195775164134199*uOther[7]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[6]*uOther[8]*mnuOther-0.2195775164134199*uOther[6]*cMOther[8]*mnuOther-0.2*cMOther[6]*uOther[7]*mnuOther-0.223606797749979*cMOther[2]*uOther[7]*mnuOther-0.2*uOther[6]*cMOther[7]*mnuOther-0.223606797749979*uOther[2]*cMOther[7]*mnuOther-0.223606797749979*cMOther[1]*uOther[6]*mnuOther-0.223606797749979*uOther[1]*cMOther[6]*mnuOther-0.223606797749979*cMOther[3]*uOther[5]*mnuOther-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(32,32) = (-0.3833333333333334*cMOther[9]*uOther[9]*mnuOther)-0.1963961012123931*cMOther[2]*uOther[9]*mnuOther-0.1963961012123931*uOther[2]*cMOther[9]*mnuOther-0.25*cMOther[8]*uOther[8]*mnuOther-0.3928571428571428*cMOther[7]*uOther[7]*mnuOther-0.223606797749979*cMOther[1]*uOther[7]*mnuOther-0.223606797749979*uOther[1]*cMOther[7]*mnuOther-0.45*cMOther[6]*uOther[6]*mnuOther-0.3928571428571428*cMOther[5]*uOther[5]*mnuOther-0.223606797749979*cMOther[0]*uOther[5]*mnuOther-1.341640786499874*m0rOther[5]*mnuOther-0.223606797749979*uOther[0]*cMOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.25*cMOther[4]*uOther[4]*mnuOther-0.45*cMOther[3]*uOther[3]*mnuOther-0.45*cMOther[2]*uOther[2]*mnuOther-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(32,33) = (-0.1963961012123931*cMOther[3]*uOther[9]*mnuOther)-0.1963961012123931*uOther[3]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[4]*uOther[8]*mnuOther-0.2195775164134199*uOther[4]*cMOther[8]*mnuOther-0.3928571428571429*cMOther[5]*uOther[7]*mnuOther-0.2*cMOther[4]*uOther[7]*mnuOther-0.223606797749979*cMOther[0]*uOther[7]*mnuOther-1.341640786499874*m0rOther[7]*mnuOther-0.3928571428571429*uOther[5]*cMOther[7]*mnuOther-0.2*uOther[4]*cMOther[7]*mnuOther-0.223606797749979*uOther[0]*cMOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-0.4024922359499621*cMOther[3]*uOther[6]*mnuOther-0.4024922359499621*uOther[3]*cMOther[6]*mnuOther-0.223606797749979*cMOther[1]*uOther[5]*mnuOther-0.223606797749979*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[4]*mnuOther-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.45*cMOther[2]*uOther[3]*mnuOther-0.45*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(32,34) = (-0.2195775164134199*cMOther[3]*uOther[8]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[8]*mnuOther-0.2*cMOther[3]*uOther[7]*mnuOther-0.2*uOther[3]*cMOther[7]*mnuOther-0.223606797749979*cMOther[5]*uOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[6]*mnuOther-0.2500000000000001*cMOther[0]*uOther[6]*mnuOther-1.5*m0rOther[6]*mnuOther-0.223606797749979*uOther[5]*cMOther[6]*mnuOther-0.159719141249985*uOther[4]*cMOther[6]*mnuOther-0.2500000000000001*uOther[0]*cMOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther-0.25*cMOther[2]*uOther[4]*mnuOther-0.25*uOther[2]*cMOther[4]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(32,35) = (-0.3273268353539885*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*cMOther[0]*uOther[9]*mnuOther-1.317465098480519*m0rOther[9]*mnuOther-0.3273268353539885*uOther[5]*cMOther[9]*mnuOther-0.2195775164134199*uOther[0]*cMOther[9]*mnuOther+0.4391550328268398*cEOther[9]*mnuOther-0.3928571428571429*cMOther[3]*uOther[7]*mnuOther-0.3928571428571429*uOther[3]*cMOther[7]*mnuOther-0.223606797749979*cMOther[4]*uOther[6]*mnuOther-0.223606797749979*uOther[4]*cMOther[6]*mnuOther-0.3928571428571428*cMOther[2]*uOther[5]*mnuOther-0.3928571428571428*uOther[2]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther-0.223606797749979*uOther[0]*cMOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(32,36) = (-0.149071198499986*cMOther[8]*uOther[8]*mnuOther)-0.2195775164134199*cMOther[1]*uOther[8]*mnuOther-0.2195775164134199*uOther[1]*cMOther[8]*mnuOther-0.223606797749979*cMOther[7]*uOther[7]*mnuOther-0.159719141249985*cMOther[6]*uOther[6]*mnuOther-0.25*cMOther[2]*uOther[6]*mnuOther-0.25*uOther[2]*cMOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[4]*mnuOther-0.2500000000000001*cMOther[0]*uOther[4]*mnuOther-1.5*m0rOther[4]*mnuOther-0.2500000000000001*uOther[0]*cMOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(32,37) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999832*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999832*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(32,39) = (-0.1309307341415954*cMOther[9]*uOther[9]*mnuOther)-0.1928571428571429*cMOther[2]*uOther[9]*mnuOther-0.1928571428571429*uOther[2]*cMOther[9]*mnuOther-0.1402829294374236*cMOther[7]*uOther[7]*mnuOther-0.2195775164134199*cMOther[1]*uOther[7]*mnuOther-0.2195775164134199*uOther[1]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[6]*uOther[6]*mnuOther-0.1402829294374236*cMOther[5]*uOther[5]*mnuOther-0.2195775164134199*cMOther[0]*uOther[5]*mnuOther-1.317465098480519*m0rOther[5]*mnuOther-0.2195775164134199*uOther[0]*cMOther[5]*mnuOther+0.4391550328268398*cEOther[5]*mnuOther-0.1963961012123931*cMOther[3]*uOther[3]*mnuOther-0.1963961012123931*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(33,30) = (-0.2195775164134199*cMOther[7]*uOther[9]*mnuOther)-0.2195775164134199*uOther[7]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[6]*uOther[8]*mnuOther-0.2195775164134199*uOther[6]*cMOther[8]*mnuOther-0.2*cMOther[6]*uOther[7]*mnuOther-0.223606797749979*cMOther[2]*uOther[7]*mnuOther-0.2*uOther[6]*cMOther[7]*mnuOther-0.223606797749979*uOther[2]*cMOther[7]*mnuOther-0.223606797749979*cMOther[1]*uOther[6]*mnuOther-0.223606797749979*uOther[1]*cMOther[6]*mnuOther-0.223606797749979*cMOther[3]*uOther[5]*mnuOther-0.223606797749979*uOther[3]*cMOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[4]*mnuOther-0.223606797749979*uOther[3]*cMOther[4]*mnuOther-0.25*cMOther[0]*uOther[3]*mnuOther-1.5*m0rOther[3]*mnuOther-0.25*uOther[0]*cMOther[3]*mnuOther+0.5*cEOther[3]*mnuOther-0.25*cMOther[1]*uOther[2]*mnuOther-0.25*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(33,31) = (-0.2195775164134199*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*uOther[5]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[3]*uOther[8]*mnuOther-0.1963961012123931*uOther[3]*cMOther[8]*mnuOther-0.4024922359499621*cMOther[3]*uOther[7]*mnuOther-0.4024922359499621*uOther[3]*cMOther[7]*mnuOther-0.2*cMOther[5]*uOther[6]*mnuOther-0.3928571428571429*cMOther[4]*uOther[6]*mnuOther-0.223606797749979*cMOther[0]*uOther[6]*mnuOther-1.341640786499874*m0rOther[6]*mnuOther-0.2*uOther[5]*cMOther[6]*mnuOther-0.3928571428571429*uOther[4]*cMOther[6]*mnuOther-0.223606797749979*uOther[0]*cMOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[5]*mnuOther-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.223606797749979*cMOther[2]*uOther[4]*mnuOther-0.223606797749979*uOther[2]*cMOther[4]*mnuOther-0.45*cMOther[1]*uOther[3]*mnuOther-0.45*uOther[1]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.25*uOther[0]*cMOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(33,32) = (-0.1963961012123931*cMOther[3]*uOther[9]*mnuOther)-0.1963961012123931*uOther[3]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[4]*uOther[8]*mnuOther-0.2195775164134199*uOther[4]*cMOther[8]*mnuOther-0.3928571428571429*cMOther[5]*uOther[7]*mnuOther-0.2*cMOther[4]*uOther[7]*mnuOther-0.223606797749979*cMOther[0]*uOther[7]*mnuOther-1.341640786499874*m0rOther[7]*mnuOther-0.3928571428571429*uOther[5]*cMOther[7]*mnuOther-0.2*uOther[4]*cMOther[7]*mnuOther-0.223606797749979*uOther[0]*cMOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-0.4024922359499621*cMOther[3]*uOther[6]*mnuOther-0.4024922359499621*uOther[3]*cMOther[6]*mnuOther-0.223606797749979*cMOther[1]*uOther[5]*mnuOther-0.223606797749979*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[4]*mnuOther-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.45*cMOther[2]*uOther[3]*mnuOther-0.45*uOther[2]*cMOther[3]*mnuOther-0.25*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.25*uOther[0]*cMOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(33,33) = (-0.3833333333333334*cMOther[9]*uOther[9]*mnuOther)-0.1963961012123931*cMOther[2]*uOther[9]*mnuOther-0.1963961012123931*uOther[2]*cMOther[9]*mnuOther-0.3833333333333334*cMOther[8]*uOther[8]*mnuOther-0.1963961012123931*cMOther[1]*uOther[8]*mnuOther-0.1963961012123931*uOther[1]*cMOther[8]*mnuOther-0.5928571428571429*cMOther[7]*uOther[7]*mnuOther-0.223606797749979*cMOther[1]*uOther[7]*mnuOther-0.223606797749979*uOther[1]*cMOther[7]*mnuOther-0.5928571428571429*cMOther[6]*uOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[6]*mnuOther-0.223606797749979*uOther[2]*cMOther[6]*mnuOther-0.3928571428571428*cMOther[5]*uOther[5]*mnuOther-0.223606797749979*cMOther[0]*uOther[5]*mnuOther-1.341640786499874*m0rOther[5]*mnuOther-0.223606797749979*uOther[0]*cMOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.3928571428571428*cMOther[4]*uOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther-0.223606797749979*uOther[0]*cMOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.65*cMOther[3]*uOther[3]*mnuOther-0.45*cMOther[2]*uOther[2]*mnuOther-0.45*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(33,34) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999831*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999831*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(33,35) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999831*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999831*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(33,36) = (-0.175662013130736*cMOther[3]*uOther[9]*mnuOther)-0.175662013130736*uOther[3]*cMOther[9]*mnuOther-0.3273268353539885*cMOther[4]*uOther[8]*mnuOther-0.21957751641342*cMOther[0]*uOther[8]*mnuOther-1.31746509848052*m0rOther[8]*mnuOther-0.3273268353539885*uOther[4]*cMOther[8]*mnuOther-0.21957751641342*uOther[0]*cMOther[8]*mnuOther+0.43915503282684*cEOther[8]*mnuOther-0.3513821107499669*cMOther[5]*uOther[7]*mnuOther-0.1788854381999831*cMOther[4]*uOther[7]*mnuOther-0.2*cMOther[0]*uOther[7]*mnuOther-1.2*m0rOther[7]*mnuOther-0.3513821107499669*uOther[5]*cMOther[7]*mnuOther-0.1788854381999831*uOther[4]*cMOther[7]*mnuOther-0.2*uOther[0]*cMOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.5528571428571428*cMOther[3]*uOther[6]*mnuOther-0.5528571428571428*uOther[3]*cMOther[6]*mnuOther-0.2*cMOther[1]*uOther[5]*mnuOther-0.2*uOther[1]*cMOther[5]*mnuOther-0.3928571428571429*cMOther[1]*uOther[4]*mnuOther-0.3928571428571429*uOther[1]*cMOther[4]*mnuOther-0.4024922359499621*cMOther[2]*uOther[3]*mnuOther-0.4024922359499621*uOther[2]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther-0.223606797749979*uOther[0]*cMOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(33,37) = (-0.3273268353539885*cMOther[5]*uOther[9]*mnuOther)-0.21957751641342*cMOther[0]*uOther[9]*mnuOther-1.31746509848052*m0rOther[9]*mnuOther-0.3273268353539885*uOther[5]*cMOther[9]*mnuOther-0.21957751641342*uOther[0]*cMOther[9]*mnuOther+0.43915503282684*cEOther[9]*mnuOther-0.175662013130736*cMOther[3]*uOther[8]*mnuOther-0.175662013130736*uOther[3]*cMOther[8]*mnuOther-0.5528571428571428*cMOther[3]*uOther[7]*mnuOther-0.5528571428571428*uOther[3]*cMOther[7]*mnuOther-0.1788854381999831*cMOther[5]*uOther[6]*mnuOther-0.3513821107499669*cMOther[4]*uOther[6]*mnuOther-0.2*cMOther[0]*uOther[6]*mnuOther-1.2*m0rOther[6]*mnuOther-0.1788854381999831*uOther[5]*cMOther[6]*mnuOther-0.3513821107499669*uOther[4]*cMOther[6]*mnuOther-0.2*uOther[0]*cMOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.3928571428571429*cMOther[2]*uOther[5]*mnuOther-0.3928571428571429*uOther[2]*cMOther[5]*mnuOther-0.2*cMOther[2]*uOther[4]*mnuOther-0.2*uOther[2]*cMOther[4]*mnuOther-0.4024922359499621*cMOther[1]*uOther[3]*mnuOther-0.4024922359499621*uOther[1]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther-0.223606797749979*uOther[0]*cMOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(33,38) = (-0.1928571428571429*cMOther[3]*uOther[8]*mnuOther)-0.1928571428571429*uOther[3]*cMOther[8]*mnuOther-0.175662013130736*cMOther[3]*uOther[7]*mnuOther-0.175662013130736*uOther[3]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[5]*uOther[6]*mnuOther-0.1402829294374237*cMOther[4]*uOther[6]*mnuOther-0.21957751641342*cMOther[0]*uOther[6]*mnuOther-1.31746509848052*m0rOther[6]*mnuOther-0.1963961012123931*uOther[5]*cMOther[6]*mnuOther-0.1402829294374237*uOther[4]*cMOther[6]*mnuOther-0.21957751641342*uOther[0]*cMOther[6]*mnuOther+0.43915503282684*cEOther[6]*mnuOther-0.2195775164134199*cMOther[2]*uOther[4]*mnuOther-0.2195775164134199*uOther[2]*cMOther[4]*mnuOther-0.1963961012123931*cMOther[1]*uOther[3]*mnuOther-0.1963961012123931*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(33,39) = (-0.1928571428571429*cMOther[3]*uOther[9]*mnuOther)-0.1928571428571429*uOther[3]*cMOther[9]*mnuOther-0.1402829294374237*cMOther[5]*uOther[7]*mnuOther-0.1963961012123931*cMOther[4]*uOther[7]*mnuOther-0.21957751641342*cMOther[0]*uOther[7]*mnuOther-1.31746509848052*m0rOther[7]*mnuOther-0.1402829294374237*uOther[5]*cMOther[7]*mnuOther-0.1963961012123931*uOther[4]*cMOther[7]*mnuOther-0.21957751641342*uOther[0]*cMOther[7]*mnuOther+0.43915503282684*cEOther[7]*mnuOther-0.175662013130736*cMOther[3]*uOther[6]*mnuOther-0.175662013130736*uOther[3]*cMOther[6]*mnuOther-0.2195775164134199*cMOther[1]*uOther[5]*mnuOther-0.2195775164134199*uOther[1]*cMOther[5]*mnuOther-0.1963961012123931*cMOther[2]*uOther[3]*mnuOther-0.1963961012123931*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(34,30) = (-0.149071198499986*cMOther[8]*uOther[8]*mnuOther)-0.2195775164134199*cMOther[1]*uOther[8]*mnuOther-0.2195775164134199*uOther[1]*cMOther[8]*mnuOther-0.223606797749979*cMOther[7]*uOther[7]*mnuOther-0.159719141249985*cMOther[6]*uOther[6]*mnuOther-0.2500000000000001*cMOther[2]*uOther[6]*mnuOther-0.2500000000000001*uOther[2]*cMOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[4]*mnuOther-0.25*cMOther[0]*uOther[4]*mnuOther-1.5*m0rOther[4]*mnuOther-0.25*uOther[0]*cMOther[4]*mnuOther+0.5*cEOther[4]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(34,31) = (-0.3273268353539885*cMOther[4]*uOther[8]*mnuOther)-0.2195775164134199*cMOther[0]*uOther[8]*mnuOther-1.317465098480519*m0rOther[8]*mnuOther-0.3273268353539885*uOther[4]*cMOther[8]*mnuOther-0.2195775164134199*uOther[0]*cMOther[8]*mnuOther+0.4391550328268398*cEOther[8]*mnuOther-0.223606797749979*cMOther[5]*uOther[7]*mnuOther-0.223606797749979*uOther[5]*cMOther[7]*mnuOther-0.3928571428571429*cMOther[3]*uOther[6]*mnuOther-0.3928571428571429*uOther[3]*cMOther[6]*mnuOther-0.3928571428571428*cMOther[1]*uOther[4]*mnuOther-0.3928571428571428*uOther[1]*cMOther[4]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther-0.223606797749979*uOther[0]*cMOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(34,32) = (-0.2195775164134199*cMOther[3]*uOther[8]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[8]*mnuOther-0.2*cMOther[3]*uOther[7]*mnuOther-0.2*uOther[3]*cMOther[7]*mnuOther-0.223606797749979*cMOther[5]*uOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[6]*mnuOther-0.2500000000000001*cMOther[0]*uOther[6]*mnuOther-1.5*m0rOther[6]*mnuOther-0.223606797749979*uOther[5]*cMOther[6]*mnuOther-0.159719141249985*uOther[4]*cMOther[6]*mnuOther-0.2500000000000001*uOther[0]*cMOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther-0.25*cMOther[2]*uOther[4]*mnuOther-0.25*uOther[2]*cMOther[4]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(34,33) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999831*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999831*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(34,34) = (-0.25*cMOther[9]*uOther[9]*mnuOther)-0.3452380952380952*cMOther[8]*uOther[8]*mnuOther-0.1402829294374236*cMOther[1]*uOther[8]*mnuOther-0.1402829294374236*uOther[1]*cMOther[8]*mnuOther-0.3928571428571428*cMOther[7]*uOther[7]*mnuOther-0.3520408163265306*cMOther[6]*uOther[6]*mnuOther-0.159719141249985*cMOther[2]*uOther[6]*mnuOther-0.159719141249985*uOther[2]*cMOther[6]*mnuOther-0.25*cMOther[5]*uOther[5]*mnuOther-0.3520408163265306*cMOther[4]*uOther[4]*mnuOther-0.159719141249985*cMOther[0]*uOther[4]*mnuOther-0.9583148474999099*m0rOther[4]*mnuOther-0.159719141249985*uOther[0]*cMOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.3928571428571428*cMOther[3]*uOther[3]*mnuOther-0.25*cMOther[2]*uOther[2]*mnuOther-0.3928571428571428*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(34,36) = (-0.2195775164134199*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*uOther[5]*cMOther[9]*mnuOther-0.1402829294374237*cMOther[3]*uOther[8]*mnuOther-0.1402829294374237*uOther[3]*cMOther[8]*mnuOther-0.3513821107499669*cMOther[3]*uOther[7]*mnuOther-0.3513821107499669*uOther[3]*cMOther[7]*mnuOther-0.1428571428571428*cMOther[5]*uOther[6]*mnuOther-0.3520408163265306*cMOther[4]*uOther[6]*mnuOther-0.159719141249985*cMOther[0]*uOther[6]*mnuOther-0.9583148474999099*m0rOther[6]*mnuOther-0.1428571428571428*uOther[5]*cMOther[6]*mnuOther-0.3520408163265306*uOther[4]*cMOther[6]*mnuOther-0.159719141249985*uOther[0]*cMOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[5]*mnuOther-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.159719141249985*cMOther[2]*uOther[4]*mnuOther-0.159719141249985*uOther[2]*cMOther[4]*mnuOther-0.3928571428571429*cMOther[1]*uOther[3]*mnuOther-0.3928571428571429*uOther[1]*cMOther[3]*mnuOther-0.2500000000000001*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.2500000000000001*uOther[0]*cMOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(34,37) = (-0.1963961012123931*cMOther[3]*uOther[9]*mnuOther)-0.1963961012123931*uOther[3]*cMOther[9]*mnuOther-0.1428571428571428*cMOther[5]*uOther[7]*mnuOther-0.2*cMOther[4]*uOther[7]*mnuOther-0.223606797749979*cMOther[0]*uOther[7]*mnuOther-1.341640786499874*m0rOther[7]*mnuOther-0.1428571428571428*uOther[5]*cMOther[7]*mnuOther-0.2*uOther[4]*cMOther[7]*mnuOther-0.223606797749979*uOther[0]*cMOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther-0.1788854381999831*cMOther[3]*uOther[6]*mnuOther-0.1788854381999831*uOther[3]*cMOther[6]*mnuOther-0.223606797749979*cMOther[1]*uOther[5]*mnuOther-0.223606797749979*uOther[1]*cMOther[5]*mnuOther-0.2*cMOther[2]*uOther[3]*mnuOther-0.2*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(34,38) = (-0.2817460317460317*cMOther[4]*uOther[8]*mnuOther)-0.149071198499986*cMOther[0]*uOther[8]*mnuOther-0.8944271909999159*m0rOther[8]*mnuOther-0.2817460317460317*uOther[4]*cMOther[8]*mnuOther-0.149071198499986*uOther[0]*cMOther[8]*mnuOther+0.2981423969999719*cEOther[8]*mnuOther-0.2195775164134199*cMOther[5]*uOther[7]*mnuOther-0.2195775164134199*uOther[5]*cMOther[7]*mnuOther-0.3273268353539885*cMOther[3]*uOther[6]*mnuOther-0.3273268353539885*uOther[3]*cMOther[6]*mnuOther-0.3273268353539885*cMOther[1]*uOther[4]*mnuOther-0.3273268353539885*uOther[1]*cMOther[4]*mnuOther-0.2195775164134199*cMOther[2]*uOther[3]*mnuOther-0.2195775164134199*uOther[2]*cMOther[3]*mnuOther-0.2195775164134199*cMOther[0]*uOther[1]*mnuOther-1.317465098480519*m0rOther[1]*mnuOther-0.2195775164134199*uOther[0]*cMOther[1]*mnuOther+0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(35,30) = (-0.149071198499986*cMOther[9]*uOther[9]*mnuOther)-0.2195775164134199*cMOther[2]*uOther[9]*mnuOther-0.2195775164134199*uOther[2]*cMOther[9]*mnuOther-0.159719141249985*cMOther[7]*uOther[7]*mnuOther-0.2500000000000001*cMOther[1]*uOther[7]*mnuOther-0.2500000000000001*uOther[1]*cMOther[7]*mnuOther-0.223606797749979*cMOther[6]*uOther[6]*mnuOther-0.159719141249985*cMOther[5]*uOther[5]*mnuOther-0.25*cMOther[0]*uOther[5]*mnuOther-1.5*m0rOther[5]*mnuOther-0.25*uOther[0]*cMOther[5]*mnuOther+0.5*cEOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(35,31) = (-0.2195775164134199*cMOther[3]*uOther[9]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[9]*mnuOther-0.159719141249985*cMOther[5]*uOther[7]*mnuOther-0.223606797749979*cMOther[4]*uOther[7]*mnuOther-0.2500000000000001*cMOther[0]*uOther[7]*mnuOther-1.5*m0rOther[7]*mnuOther-0.159719141249985*uOther[5]*cMOther[7]*mnuOther-0.223606797749979*uOther[4]*cMOther[7]*mnuOther-0.2500000000000001*uOther[0]*cMOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther-0.2*cMOther[3]*uOther[6]*mnuOther-0.2*uOther[3]*cMOther[6]*mnuOther-0.25*cMOther[1]*uOther[5]*mnuOther-0.25*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(35,32) = (-0.3273268353539885*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*cMOther[0]*uOther[9]*mnuOther-1.317465098480519*m0rOther[9]*mnuOther-0.3273268353539885*uOther[5]*cMOther[9]*mnuOther-0.2195775164134199*uOther[0]*cMOther[9]*mnuOther+0.4391550328268398*cEOther[9]*mnuOther-0.3928571428571429*cMOther[3]*uOther[7]*mnuOther-0.3928571428571429*uOther[3]*cMOther[7]*mnuOther-0.223606797749979*cMOther[4]*uOther[6]*mnuOther-0.223606797749979*uOther[4]*cMOther[6]*mnuOther-0.3928571428571428*cMOther[2]*uOther[5]*mnuOther-0.3928571428571428*uOther[2]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther-0.223606797749979*uOther[0]*cMOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(35,33) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999831*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999831*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(35,35) = (-0.3452380952380952*cMOther[9]*uOther[9]*mnuOther)-0.1402829294374236*cMOther[2]*uOther[9]*mnuOther-0.1402829294374236*uOther[2]*cMOther[9]*mnuOther-0.25*cMOther[8]*uOther[8]*mnuOther-0.3520408163265306*cMOther[7]*uOther[7]*mnuOther-0.159719141249985*cMOther[1]*uOther[7]*mnuOther-0.159719141249985*uOther[1]*cMOther[7]*mnuOther-0.3928571428571428*cMOther[6]*uOther[6]*mnuOther-0.3520408163265306*cMOther[5]*uOther[5]*mnuOther-0.159719141249985*cMOther[0]*uOther[5]*mnuOther-0.9583148474999099*m0rOther[5]*mnuOther-0.159719141249985*uOther[0]*cMOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.25*cMOther[4]*uOther[4]*mnuOther-0.3928571428571428*cMOther[3]*uOther[3]*mnuOther-0.3928571428571428*cMOther[2]*uOther[2]*mnuOther-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(35,36) = (-0.1963961012123931*cMOther[3]*uOther[8]*mnuOther)-0.1963961012123931*uOther[3]*cMOther[8]*mnuOther-0.1788854381999831*cMOther[3]*uOther[7]*mnuOther-0.1788854381999831*uOther[3]*cMOther[7]*mnuOther-0.2*cMOther[5]*uOther[6]*mnuOther-0.1428571428571428*cMOther[4]*uOther[6]*mnuOther-0.223606797749979*cMOther[0]*uOther[6]*mnuOther-1.341640786499874*m0rOther[6]*mnuOther-0.2*uOther[5]*cMOther[6]*mnuOther-0.1428571428571428*uOther[4]*cMOther[6]*mnuOther-0.223606797749979*uOther[0]*cMOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[4]*mnuOther-0.223606797749979*uOther[2]*cMOther[4]*mnuOther-0.2*cMOther[1]*uOther[3]*mnuOther-0.2*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(35,37) = (-0.1402829294374237*cMOther[3]*uOther[9]*mnuOther)-0.1402829294374237*uOther[3]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[4]*uOther[8]*mnuOther-0.2195775164134199*uOther[4]*cMOther[8]*mnuOther-0.3520408163265306*cMOther[5]*uOther[7]*mnuOther-0.1428571428571428*cMOther[4]*uOther[7]*mnuOther-0.159719141249985*cMOther[0]*uOther[7]*mnuOther-0.9583148474999099*m0rOther[7]*mnuOther-0.3520408163265306*uOther[5]*cMOther[7]*mnuOther-0.1428571428571428*uOther[4]*cMOther[7]*mnuOther-0.159719141249985*uOther[0]*cMOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-0.3513821107499669*cMOther[3]*uOther[6]*mnuOther-0.3513821107499669*uOther[3]*cMOther[6]*mnuOther-0.159719141249985*cMOther[1]*uOther[5]*mnuOther-0.159719141249985*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[4]*mnuOther-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.3928571428571429*cMOther[2]*uOther[3]*mnuOther-0.3928571428571429*uOther[2]*cMOther[3]*mnuOther-0.2500000000000001*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.2500000000000001*uOther[0]*cMOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(35,39) = (-0.2817460317460317*cMOther[5]*uOther[9]*mnuOther)-0.149071198499986*cMOther[0]*uOther[9]*mnuOther-0.8944271909999159*m0rOther[9]*mnuOther-0.2817460317460317*uOther[5]*cMOther[9]*mnuOther-0.149071198499986*uOther[0]*cMOther[9]*mnuOther+0.2981423969999719*cEOther[9]*mnuOther-0.3273268353539885*cMOther[3]*uOther[7]*mnuOther-0.3273268353539885*uOther[3]*cMOther[7]*mnuOther-0.2195775164134199*cMOther[4]*uOther[6]*mnuOther-0.2195775164134199*uOther[4]*cMOther[6]*mnuOther-0.3273268353539885*cMOther[2]*uOther[5]*mnuOther-0.3273268353539885*uOther[2]*cMOther[5]*mnuOther-0.2195775164134199*cMOther[1]*uOther[3]*mnuOther-0.2195775164134199*uOther[1]*cMOther[3]*mnuOther-0.2195775164134199*cMOther[0]*uOther[2]*mnuOther-1.317465098480519*m0rOther[2]*mnuOther-0.2195775164134199*uOther[0]*cMOther[2]*mnuOther+0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(36,30) = (-0.2195775164134199*cMOther[3]*uOther[8]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[8]*mnuOther-0.2*cMOther[3]*uOther[7]*mnuOther-0.2*uOther[3]*cMOther[7]*mnuOther-0.223606797749979*cMOther[5]*uOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[6]*mnuOther-0.25*cMOther[0]*uOther[6]*mnuOther-1.5*m0rOther[6]*mnuOther-0.223606797749979*uOther[5]*cMOther[6]*mnuOther-0.159719141249985*uOther[4]*cMOther[6]*mnuOther-0.25*uOther[0]*cMOther[6]*mnuOther+0.5*cEOther[6]*mnuOther-0.2500000000000001*cMOther[2]*uOther[4]*mnuOther-0.2500000000000001*uOther[2]*cMOther[4]*mnuOther-0.223606797749979*cMOther[1]*uOther[3]*mnuOther-0.223606797749979*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(36,31) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999832*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999832*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(36,32) = (-0.149071198499986*cMOther[8]*uOther[8]*mnuOther)-0.2195775164134199*cMOther[1]*uOther[8]*mnuOther-0.2195775164134199*uOther[1]*cMOther[8]*mnuOther-0.223606797749979*cMOther[7]*uOther[7]*mnuOther-0.159719141249985*cMOther[6]*uOther[6]*mnuOther-0.25*cMOther[2]*uOther[6]*mnuOther-0.25*uOther[2]*cMOther[6]*mnuOther-0.159719141249985*cMOther[4]*uOther[4]*mnuOther-0.2500000000000001*cMOther[0]*uOther[4]*mnuOther-1.5*m0rOther[4]*mnuOther-0.2500000000000001*uOther[0]*cMOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(36,33) = (-0.175662013130736*cMOther[3]*uOther[9]*mnuOther)-0.175662013130736*uOther[3]*cMOther[9]*mnuOther-0.3273268353539885*cMOther[4]*uOther[8]*mnuOther-0.21957751641342*cMOther[0]*uOther[8]*mnuOther-1.31746509848052*m0rOther[8]*mnuOther-0.3273268353539885*uOther[4]*cMOther[8]*mnuOther-0.21957751641342*uOther[0]*cMOther[8]*mnuOther+0.43915503282684*cEOther[8]*mnuOther-0.3513821107499669*cMOther[5]*uOther[7]*mnuOther-0.1788854381999831*cMOther[4]*uOther[7]*mnuOther-0.2*cMOther[0]*uOther[7]*mnuOther-1.2*m0rOther[7]*mnuOther-0.3513821107499669*uOther[5]*cMOther[7]*mnuOther-0.1788854381999831*uOther[4]*cMOther[7]*mnuOther-0.2*uOther[0]*cMOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.5528571428571428*cMOther[3]*uOther[6]*mnuOther-0.5528571428571428*uOther[3]*cMOther[6]*mnuOther-0.2*cMOther[1]*uOther[5]*mnuOther-0.2*uOther[1]*cMOther[5]*mnuOther-0.3928571428571429*cMOther[1]*uOther[4]*mnuOther-0.3928571428571429*uOther[1]*cMOther[4]*mnuOther-0.4024922359499621*cMOther[2]*uOther[3]*mnuOther-0.4024922359499621*uOther[2]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther-0.223606797749979*uOther[0]*cMOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(36,34) = (-0.2195775164134199*cMOther[5]*uOther[9]*mnuOther)-0.2195775164134199*uOther[5]*cMOther[9]*mnuOther-0.1402829294374237*cMOther[3]*uOther[8]*mnuOther-0.1402829294374237*uOther[3]*cMOther[8]*mnuOther-0.3513821107499669*cMOther[3]*uOther[7]*mnuOther-0.3513821107499669*uOther[3]*cMOther[7]*mnuOther-0.1428571428571428*cMOther[5]*uOther[6]*mnuOther-0.3520408163265306*cMOther[4]*uOther[6]*mnuOther-0.159719141249985*cMOther[0]*uOther[6]*mnuOther-0.9583148474999099*m0rOther[6]*mnuOther-0.1428571428571428*uOther[5]*cMOther[6]*mnuOther-0.3520408163265306*uOther[4]*cMOther[6]*mnuOther-0.159719141249985*uOther[0]*cMOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[5]*mnuOther-0.223606797749979*uOther[2]*cMOther[5]*mnuOther-0.159719141249985*cMOther[2]*uOther[4]*mnuOther-0.159719141249985*uOther[2]*cMOther[4]*mnuOther-0.3928571428571429*cMOther[1]*uOther[3]*mnuOther-0.3928571428571429*uOther[1]*cMOther[3]*mnuOther-0.2500000000000001*cMOther[0]*uOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther-0.2500000000000001*uOther[0]*cMOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(36,35) = (-0.1963961012123931*cMOther[3]*uOther[8]*mnuOther)-0.1963961012123931*uOther[3]*cMOther[8]*mnuOther-0.1788854381999831*cMOther[3]*uOther[7]*mnuOther-0.1788854381999831*uOther[3]*cMOther[7]*mnuOther-0.2*cMOther[5]*uOther[6]*mnuOther-0.1428571428571428*cMOther[4]*uOther[6]*mnuOther-0.223606797749979*cMOther[0]*uOther[6]*mnuOther-1.341640786499874*m0rOther[6]*mnuOther-0.2*uOther[5]*cMOther[6]*mnuOther-0.1428571428571428*uOther[4]*cMOther[6]*mnuOther-0.223606797749979*uOther[0]*cMOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[4]*mnuOther-0.223606797749979*uOther[2]*cMOther[4]*mnuOther-0.2*cMOther[1]*uOther[3]*mnuOther-0.2*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(36,36) = (-0.3833333333333334*cMOther[9]*uOther[9]*mnuOther)-0.1963961012123931*cMOther[2]*uOther[9]*mnuOther-0.1963961012123931*uOther[2]*cMOther[9]*mnuOther-0.3452380952380952*cMOther[8]*uOther[8]*mnuOther-0.1402829294374236*cMOther[1]*uOther[8]*mnuOther-0.1402829294374236*uOther[1]*cMOther[8]*mnuOther-0.5357142857142857*cMOther[7]*uOther[7]*mnuOther-0.223606797749979*cMOther[1]*uOther[7]*mnuOther-0.223606797749979*uOther[1]*cMOther[7]*mnuOther-0.5520408163265306*cMOther[6]*uOther[6]*mnuOther-0.159719141249985*cMOther[2]*uOther[6]*mnuOther-0.159719141249985*uOther[2]*cMOther[6]*mnuOther-0.3928571428571428*cMOther[5]*uOther[5]*mnuOther-0.223606797749979*cMOther[0]*uOther[5]*mnuOther-1.341640786499874*m0rOther[5]*mnuOther-0.223606797749979*uOther[0]*cMOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.3520408163265306*cMOther[4]*uOther[4]*mnuOther-0.159719141249985*cMOther[0]*uOther[4]*mnuOther-0.9583148474999099*m0rOther[4]*mnuOther-0.159719141249985*uOther[0]*cMOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.5928571428571429*cMOther[3]*uOther[3]*mnuOther-0.45*cMOther[2]*uOther[2]*mnuOther-0.3928571428571428*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(36,37) = (-0.175662013130736*cMOther[7]*uOther[9]*mnuOther)-0.175662013130736*uOther[7]*cMOther[9]*mnuOther-0.175662013130736*cMOther[6]*uOther[8]*mnuOther-0.175662013130736*uOther[6]*cMOther[8]*mnuOther-0.16*cMOther[6]*uOther[7]*mnuOther-0.1788854381999832*cMOther[2]*uOther[7]*mnuOther-0.16*uOther[6]*cMOther[7]*mnuOther-0.1788854381999832*uOther[2]*cMOther[7]*mnuOther-0.1788854381999832*cMOther[1]*uOther[6]*mnuOther-0.1788854381999832*uOther[1]*cMOther[6]*mnuOther-0.1788854381999831*cMOther[3]*uOther[5]*mnuOther-0.1788854381999831*uOther[3]*cMOther[5]*mnuOther-0.1788854381999831*cMOther[3]*uOther[4]*mnuOther-0.1788854381999831*uOther[3]*cMOther[4]*mnuOther-0.2*cMOther[0]*uOther[3]*mnuOther-1.2*m0rOther[3]*mnuOther-0.2*uOther[0]*cMOther[3]*mnuOther+0.4*cEOther[3]*mnuOther-0.2*cMOther[1]*uOther[2]*mnuOther-0.2*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(36,38) = (-0.1928571428571429*cMOther[7]*uOther[9]*mnuOther)-0.1928571428571429*uOther[7]*cMOther[9]*mnuOther-0.1928571428571429*cMOther[6]*uOther[8]*mnuOther-0.1928571428571429*uOther[6]*cMOther[8]*mnuOther-0.175662013130736*cMOther[6]*uOther[7]*mnuOther-0.1963961012123931*cMOther[2]*uOther[7]*mnuOther-0.175662013130736*uOther[6]*cMOther[7]*mnuOther-0.1963961012123931*uOther[2]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[1]*uOther[6]*mnuOther-0.1963961012123931*uOther[1]*cMOther[6]*mnuOther-0.1963961012123931*cMOther[3]*uOther[5]*mnuOther-0.1963961012123931*uOther[3]*cMOther[5]*mnuOther-0.1963961012123931*cMOther[3]*uOther[4]*mnuOther-0.1963961012123931*uOther[3]*cMOther[4]*mnuOther-0.21957751641342*cMOther[0]*uOther[3]*mnuOther-1.31746509848052*m0rOther[3]*mnuOther-0.21957751641342*uOther[0]*cMOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther-0.21957751641342*cMOther[1]*uOther[2]*mnuOther-0.21957751641342*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(37,30) = (-0.2195775164134199*cMOther[3]*uOther[9]*mnuOther)-0.2195775164134199*uOther[3]*cMOther[9]*mnuOther-0.159719141249985*cMOther[5]*uOther[7]*mnuOther-0.223606797749979*cMOther[4]*uOther[7]*mnuOther-0.25*cMOther[0]*uOther[7]*mnuOther-1.5*m0rOther[7]*mnuOther-0.159719141249985*uOther[5]*cMOther[7]*mnuOther-0.223606797749979*uOther[4]*cMOther[7]*mnuOther-0.25*uOther[0]*cMOther[7]*mnuOther+0.5*cEOther[7]*mnuOther-0.2*cMOther[3]*uOther[6]*mnuOther-0.2*uOther[3]*cMOther[6]*mnuOther-0.2500000000000001*cMOther[1]*uOther[5]*mnuOther-0.2500000000000001*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[2]*uOther[3]*mnuOther-0.223606797749979*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(37,31) = (-0.149071198499986*cMOther[9]*uOther[9]*mnuOther)-0.2195775164134199*cMOther[2]*uOther[9]*mnuOther-0.2195775164134199*uOther[2]*cMOther[9]*mnuOther-0.159719141249985*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[1]*uOther[7]*mnuOther-0.25*uOther[1]*cMOther[7]*mnuOther-0.223606797749979*cMOther[6]*uOther[6]*mnuOther-0.159719141249985*cMOther[5]*uOther[5]*mnuOther-0.2500000000000001*cMOther[0]*uOther[5]*mnuOther-1.5*m0rOther[5]*mnuOther-0.2500000000000001*uOther[0]*cMOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther-0.223606797749979*cMOther[3]*uOther[3]*mnuOther-0.223606797749979*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(37,32) = (-0.1963961012123931*cMOther[7]*uOther[9]*mnuOther)-0.1963961012123931*uOther[7]*cMOther[9]*mnuOther-0.1963961012123931*cMOther[6]*uOther[8]*mnuOther-0.1963961012123931*uOther[6]*cMOther[8]*mnuOther-0.1788854381999832*cMOther[6]*uOther[7]*mnuOther-0.2*cMOther[2]*uOther[7]*mnuOther-0.1788854381999832*uOther[6]*cMOther[7]*mnuOther-0.2*uOther[2]*cMOther[7]*mnuOther-0.2*cMOther[1]*uOther[6]*mnuOther-0.2*uOther[1]*cMOther[6]*mnuOther-0.2*cMOther[3]*uOther[5]*mnuOther-0.2*uOther[3]*cMOther[5]*mnuOther-0.2*cMOther[3]*uOther[4]*mnuOther-0.2*uOther[3]*cMOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[3]*mnuOther-1.341640786499874*m0rOther[3]*mnuOther-0.223606797749979*uOther[0]*cMOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther-0.223606797749979*cMOther[1]*uOther[2]*mnuOther-0.223606797749979*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(37,33) = (-0.3273268353539885*cMOther[5]*uOther[9]*mnuOther)-0.21957751641342*cMOther[0]*uOther[9]*mnuOther-1.31746509848052*m0rOther[9]*mnuOther-0.3273268353539885*uOther[5]*cMOther[9]*mnuOther-0.21957751641342*uOther[0]*cMOther[9]*mnuOther+0.43915503282684*cEOther[9]*mnuOther-0.175662013130736*cMOther[3]*uOther[8]*mnuOther-0.175662013130736*uOther[3]*cMOther[8]*mnuOther-0.5528571428571428*cMOther[3]*uOther[7]*mnuOther-0.5528571428571428*uOther[3]*cMOther[7]*mnuOther-0.1788854381999831*cMOther[5]*uOther[6]*mnuOther-0.3513821107499669*cMOther[4]*uOther[6]*mnuOther-0.2*cMOther[0]*uOther[6]*mnuOther-1.2*m0rOther[6]*mnuOther-0.1788854381999831*uOther[5]*cMOther[6]*mnuOther-0.3513821107499669*uOther[4]*cMOther[6]*mnuOther-0.2*uOther[0]*cMOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.3928571428571429*cMOther[2]*uOther[5]*mnuOther-0.3928571428571429*uOther[2]*cMOther[5]*mnuOther-0.2*cMOther[2]*uOther[4]*mnuOther-0.2*uOther[2]*cMOther[4]*mnuOther-0.4024922359499621*cMOther[1]*uOther[3]*mnuOther-0.4024922359499621*uOther[1]*cMOther[3]*mnuOther-0.223606797749979*cMOther[0]*uOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther-0.223606797749979*uOther[0]*cMOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(37,34) = (-0.1963961012123931*cMOther[3]*uOther[9]*mnuOther)-0.1963961012123931*uOther[3]*cMOther[9]*mnuOther-0.1428571428571428*cMOther[5]*uOther[7]*mnuOther-0.2*cMOther[4]*uOther[7]*mnuOther-0.223606797749979*cMOther[0]*uOther[7]*mnuOther-1.341640786499874*m0rOther[7]*mnuOther-0.1428571428571428*uOther[5]*cMOther[7]*mnuOther-0.2*uOther[4]*cMOther[7]*mnuOther-0.223606797749979*uOther[0]*cMOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther-0.1788854381999831*cMOther[3]*uOther[6]*mnuOther-0.1788854381999831*uOther[3]*cMOther[6]*mnuOther-0.223606797749979*cMOther[1]*uOther[5]*mnuOther-0.223606797749979*uOther[1]*cMOther[5]*mnuOther-0.2*cMOther[2]*uOther[3]*mnuOther-0.2*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(37,35) = (-0.1402829294374237*cMOther[3]*uOther[9]*mnuOther)-0.1402829294374237*uOther[3]*cMOther[9]*mnuOther-0.2195775164134199*cMOther[4]*uOther[8]*mnuOther-0.2195775164134199*uOther[4]*cMOther[8]*mnuOther-0.3520408163265306*cMOther[5]*uOther[7]*mnuOther-0.1428571428571428*cMOther[4]*uOther[7]*mnuOther-0.159719141249985*cMOther[0]*uOther[7]*mnuOther-0.9583148474999099*m0rOther[7]*mnuOther-0.3520408163265306*uOther[5]*cMOther[7]*mnuOther-0.1428571428571428*uOther[4]*cMOther[7]*mnuOther-0.159719141249985*uOther[0]*cMOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-0.3513821107499669*cMOther[3]*uOther[6]*mnuOther-0.3513821107499669*uOther[3]*cMOther[6]*mnuOther-0.159719141249985*cMOther[1]*uOther[5]*mnuOther-0.159719141249985*uOther[1]*cMOther[5]*mnuOther-0.223606797749979*cMOther[1]*uOther[4]*mnuOther-0.223606797749979*uOther[1]*cMOther[4]*mnuOther-0.3928571428571429*cMOther[2]*uOther[3]*mnuOther-0.3928571428571429*uOther[2]*cMOther[3]*mnuOther-0.2500000000000001*cMOther[0]*uOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther-0.2500000000000001*uOther[0]*cMOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(37,36) = (-0.175662013130736*cMOther[7]*uOther[9]*mnuOther)-0.175662013130736*uOther[7]*cMOther[9]*mnuOther-0.175662013130736*cMOther[6]*uOther[8]*mnuOther-0.175662013130736*uOther[6]*cMOther[8]*mnuOther-0.16*cMOther[6]*uOther[7]*mnuOther-0.1788854381999832*cMOther[2]*uOther[7]*mnuOther-0.16*uOther[6]*cMOther[7]*mnuOther-0.1788854381999832*uOther[2]*cMOther[7]*mnuOther-0.1788854381999832*cMOther[1]*uOther[6]*mnuOther-0.1788854381999832*uOther[1]*cMOther[6]*mnuOther-0.1788854381999831*cMOther[3]*uOther[5]*mnuOther-0.1788854381999831*uOther[3]*cMOther[5]*mnuOther-0.1788854381999831*cMOther[3]*uOther[4]*mnuOther-0.1788854381999831*uOther[3]*cMOther[4]*mnuOther-0.2*cMOther[0]*uOther[3]*mnuOther-1.2*m0rOther[3]*mnuOther-0.2*uOther[0]*cMOther[3]*mnuOther+0.4*cEOther[3]*mnuOther-0.2*cMOther[1]*uOther[2]*mnuOther-0.2*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(37,37) = (-0.3452380952380952*cMOther[9]*uOther[9]*mnuOther)-0.1402829294374236*cMOther[2]*uOther[9]*mnuOther-0.1402829294374236*uOther[2]*cMOther[9]*mnuOther-0.3833333333333334*cMOther[8]*uOther[8]*mnuOther-0.1963961012123931*cMOther[1]*uOther[8]*mnuOther-0.1963961012123931*uOther[1]*cMOther[8]*mnuOther-0.5520408163265306*cMOther[7]*uOther[7]*mnuOther-0.159719141249985*cMOther[1]*uOther[7]*mnuOther-0.159719141249985*uOther[1]*cMOther[7]*mnuOther-0.5357142857142857*cMOther[6]*uOther[6]*mnuOther-0.223606797749979*cMOther[2]*uOther[6]*mnuOther-0.223606797749979*uOther[2]*cMOther[6]*mnuOther-0.3520408163265306*cMOther[5]*uOther[5]*mnuOther-0.159719141249985*cMOther[0]*uOther[5]*mnuOther-0.9583148474999099*m0rOther[5]*mnuOther-0.159719141249985*uOther[0]*cMOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.3928571428571428*cMOther[4]*uOther[4]*mnuOther-0.223606797749979*cMOther[0]*uOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther-0.223606797749979*uOther[0]*cMOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5928571428571429*cMOther[3]*uOther[3]*mnuOther-0.3928571428571428*cMOther[2]*uOther[2]*mnuOther-0.45*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(37,39) = (-0.1928571428571429*cMOther[7]*uOther[9]*mnuOther)-0.1928571428571429*uOther[7]*cMOther[9]*mnuOther-0.1928571428571429*cMOther[6]*uOther[8]*mnuOther-0.1928571428571429*uOther[6]*cMOther[8]*mnuOther-0.175662013130736*cMOther[6]*uOther[7]*mnuOther-0.1963961012123931*cMOther[2]*uOther[7]*mnuOther-0.175662013130736*uOther[6]*cMOther[7]*mnuOther-0.1963961012123931*uOther[2]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[1]*uOther[6]*mnuOther-0.1963961012123931*uOther[1]*cMOther[6]*mnuOther-0.1963961012123931*cMOther[3]*uOther[5]*mnuOther-0.1963961012123931*uOther[3]*cMOther[5]*mnuOther-0.1963961012123931*cMOther[3]*uOther[4]*mnuOther-0.1963961012123931*uOther[3]*cMOther[4]*mnuOther-0.21957751641342*cMOther[0]*uOther[3]*mnuOther-1.31746509848052*m0rOther[3]*mnuOther-0.21957751641342*uOther[0]*cMOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther-0.21957751641342*cMOther[1]*uOther[2]*mnuOther-0.21957751641342*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(38,30) = (-0.149071198499986*cMOther[4]*uOther[8]*mnuOther)-0.25*cMOther[0]*uOther[8]*mnuOther-1.5*m0rOther[8]*mnuOther-0.149071198499986*uOther[4]*cMOther[8]*mnuOther-0.25*uOther[0]*cMOther[8]*mnuOther+0.5*cEOther[8]*mnuOther-0.2195775164134199*cMOther[3]*uOther[6]*mnuOther-0.2195775164134199*uOther[3]*cMOther[6]*mnuOther-0.2195775164134199*cMOther[1]*uOther[4]*mnuOther-0.2195775164134199*uOther[1]*cMOther[4]*mnuOther; 
  data->AEM_S(38,31) = (-0.1309307341415954*cMOther[8]*uOther[8]*mnuOther)-0.1928571428571429*cMOther[1]*uOther[8]*mnuOther-0.1928571428571429*uOther[1]*cMOther[8]*mnuOther-0.1963961012123931*cMOther[7]*uOther[7]*mnuOther-0.1402829294374236*cMOther[6]*uOther[6]*mnuOther-0.2195775164134199*cMOther[2]*uOther[6]*mnuOther-0.2195775164134199*uOther[2]*cMOther[6]*mnuOther-0.1402829294374236*cMOther[4]*uOther[4]*mnuOther-0.2195775164134199*cMOther[0]*uOther[4]*mnuOther-1.317465098480519*m0rOther[4]*mnuOther-0.2195775164134199*uOther[0]*cMOther[4]*mnuOther+0.4391550328268398*cEOther[4]*mnuOther-0.1963961012123931*cMOther[3]*uOther[3]*mnuOther-0.1963961012123931*cMOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(38,33) = (-0.1928571428571429*cMOther[3]*uOther[8]*mnuOther)-0.1928571428571429*uOther[3]*cMOther[8]*mnuOther-0.175662013130736*cMOther[3]*uOther[7]*mnuOther-0.175662013130736*uOther[3]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[5]*uOther[6]*mnuOther-0.1402829294374237*cMOther[4]*uOther[6]*mnuOther-0.21957751641342*cMOther[0]*uOther[6]*mnuOther-1.31746509848052*m0rOther[6]*mnuOther-0.1963961012123931*uOther[5]*cMOther[6]*mnuOther-0.1402829294374237*uOther[4]*cMOther[6]*mnuOther-0.21957751641342*uOther[0]*cMOther[6]*mnuOther+0.43915503282684*cEOther[6]*mnuOther-0.2195775164134199*cMOther[2]*uOther[4]*mnuOther-0.2195775164134199*uOther[2]*cMOther[4]*mnuOther-0.1963961012123931*cMOther[1]*uOther[3]*mnuOther-0.1963961012123931*uOther[1]*cMOther[3]*mnuOther; 
  data->AEM_S(38,34) = (-0.2817460317460317*cMOther[4]*uOther[8]*mnuOther)-0.149071198499986*cMOther[0]*uOther[8]*mnuOther-0.8944271909999159*m0rOther[8]*mnuOther-0.2817460317460317*uOther[4]*cMOther[8]*mnuOther-0.149071198499986*uOther[0]*cMOther[8]*mnuOther+0.2981423969999719*cEOther[8]*mnuOther-0.2195775164134199*cMOther[5]*uOther[7]*mnuOther-0.2195775164134199*uOther[5]*cMOther[7]*mnuOther-0.3273268353539885*cMOther[3]*uOther[6]*mnuOther-0.3273268353539885*uOther[3]*cMOther[6]*mnuOther-0.3273268353539885*cMOther[1]*uOther[4]*mnuOther-0.3273268353539885*uOther[1]*cMOther[4]*mnuOther-0.2195775164134199*cMOther[2]*uOther[3]*mnuOther-0.2195775164134199*uOther[2]*cMOther[3]*mnuOther-0.2195775164134199*cMOther[0]*uOther[1]*mnuOther-1.317465098480519*m0rOther[1]*mnuOther-0.2195775164134199*uOther[0]*cMOther[1]*mnuOther+0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(38,36) = (-0.1928571428571429*cMOther[7]*uOther[9]*mnuOther)-0.1928571428571429*uOther[7]*cMOther[9]*mnuOther-0.1928571428571429*cMOther[6]*uOther[8]*mnuOther-0.1928571428571429*uOther[6]*cMOther[8]*mnuOther-0.175662013130736*cMOther[6]*uOther[7]*mnuOther-0.1963961012123931*cMOther[2]*uOther[7]*mnuOther-0.175662013130736*uOther[6]*cMOther[7]*mnuOther-0.1963961012123931*uOther[2]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[1]*uOther[6]*mnuOther-0.1963961012123931*uOther[1]*cMOther[6]*mnuOther-0.1963961012123931*cMOther[3]*uOther[5]*mnuOther-0.1963961012123931*uOther[3]*cMOther[5]*mnuOther-0.1963961012123931*cMOther[3]*uOther[4]*mnuOther-0.1963961012123931*uOther[3]*cMOther[4]*mnuOther-0.21957751641342*cMOther[0]*uOther[3]*mnuOther-1.31746509848052*m0rOther[3]*mnuOther-0.21957751641342*uOther[0]*cMOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther-0.21957751641342*cMOther[1]*uOther[2]*mnuOther-0.21957751641342*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(38,38) = (-0.25*cMOther[9]*uOther[9]*mnuOther)-0.3388888888888889*cMOther[8]*uOther[8]*mnuOther-0.1309307341415954*cMOther[1]*uOther[8]*mnuOther-0.1309307341415954*uOther[1]*cMOther[8]*mnuOther-0.3833333333333334*cMOther[7]*uOther[7]*mnuOther-0.3452380952380952*cMOther[6]*uOther[6]*mnuOther-0.149071198499986*cMOther[2]*uOther[6]*mnuOther-0.149071198499986*uOther[2]*cMOther[6]*mnuOther-0.25*cMOther[5]*uOther[5]*mnuOther-0.3452380952380952*cMOther[4]*uOther[4]*mnuOther-0.149071198499986*cMOther[0]*uOther[4]*mnuOther-0.8944271909999159*m0rOther[4]*mnuOther-0.149071198499986*uOther[0]*cMOther[4]*mnuOther+0.2981423969999719*cEOther[4]*mnuOther-0.3833333333333334*cMOther[3]*uOther[3]*mnuOther-0.25*cMOther[2]*uOther[2]*mnuOther-0.3833333333333334*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(39,30) = (-0.149071198499986*cMOther[5]*uOther[9]*mnuOther)-0.25*cMOther[0]*uOther[9]*mnuOther-1.5*m0rOther[9]*mnuOther-0.149071198499986*uOther[5]*cMOther[9]*mnuOther-0.25*uOther[0]*cMOther[9]*mnuOther+0.5*cEOther[9]*mnuOther-0.2195775164134199*cMOther[3]*uOther[7]*mnuOther-0.2195775164134199*uOther[3]*cMOther[7]*mnuOther-0.2195775164134199*cMOther[2]*uOther[5]*mnuOther-0.2195775164134199*uOther[2]*cMOther[5]*mnuOther; 
  data->AEM_S(39,32) = (-0.1309307341415954*cMOther[9]*uOther[9]*mnuOther)-0.1928571428571429*cMOther[2]*uOther[9]*mnuOther-0.1928571428571429*uOther[2]*cMOther[9]*mnuOther-0.1402829294374236*cMOther[7]*uOther[7]*mnuOther-0.2195775164134199*cMOther[1]*uOther[7]*mnuOther-0.2195775164134199*uOther[1]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[6]*uOther[6]*mnuOther-0.1402829294374236*cMOther[5]*uOther[5]*mnuOther-0.2195775164134199*cMOther[0]*uOther[5]*mnuOther-1.317465098480519*m0rOther[5]*mnuOther-0.2195775164134199*uOther[0]*cMOther[5]*mnuOther+0.4391550328268398*cEOther[5]*mnuOther-0.1963961012123931*cMOther[3]*uOther[3]*mnuOther-0.1963961012123931*cMOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(39,33) = (-0.1928571428571429*cMOther[3]*uOther[9]*mnuOther)-0.1928571428571429*uOther[3]*cMOther[9]*mnuOther-0.1402829294374237*cMOther[5]*uOther[7]*mnuOther-0.1963961012123931*cMOther[4]*uOther[7]*mnuOther-0.21957751641342*cMOther[0]*uOther[7]*mnuOther-1.31746509848052*m0rOther[7]*mnuOther-0.1402829294374237*uOther[5]*cMOther[7]*mnuOther-0.1963961012123931*uOther[4]*cMOther[7]*mnuOther-0.21957751641342*uOther[0]*cMOther[7]*mnuOther+0.43915503282684*cEOther[7]*mnuOther-0.175662013130736*cMOther[3]*uOther[6]*mnuOther-0.175662013130736*uOther[3]*cMOther[6]*mnuOther-0.2195775164134199*cMOther[1]*uOther[5]*mnuOther-0.2195775164134199*uOther[1]*cMOther[5]*mnuOther-0.1963961012123931*cMOther[2]*uOther[3]*mnuOther-0.1963961012123931*uOther[2]*cMOther[3]*mnuOther; 
  data->AEM_S(39,35) = (-0.2817460317460317*cMOther[5]*uOther[9]*mnuOther)-0.149071198499986*cMOther[0]*uOther[9]*mnuOther-0.8944271909999159*m0rOther[9]*mnuOther-0.2817460317460317*uOther[5]*cMOther[9]*mnuOther-0.149071198499986*uOther[0]*cMOther[9]*mnuOther+0.2981423969999719*cEOther[9]*mnuOther-0.3273268353539885*cMOther[3]*uOther[7]*mnuOther-0.3273268353539885*uOther[3]*cMOther[7]*mnuOther-0.2195775164134199*cMOther[4]*uOther[6]*mnuOther-0.2195775164134199*uOther[4]*cMOther[6]*mnuOther-0.3273268353539885*cMOther[2]*uOther[5]*mnuOther-0.3273268353539885*uOther[2]*cMOther[5]*mnuOther-0.2195775164134199*cMOther[1]*uOther[3]*mnuOther-0.2195775164134199*uOther[1]*cMOther[3]*mnuOther-0.2195775164134199*cMOther[0]*uOther[2]*mnuOther-1.317465098480519*m0rOther[2]*mnuOther-0.2195775164134199*uOther[0]*cMOther[2]*mnuOther+0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(39,37) = (-0.1928571428571429*cMOther[7]*uOther[9]*mnuOther)-0.1928571428571429*uOther[7]*cMOther[9]*mnuOther-0.1928571428571429*cMOther[6]*uOther[8]*mnuOther-0.1928571428571429*uOther[6]*cMOther[8]*mnuOther-0.175662013130736*cMOther[6]*uOther[7]*mnuOther-0.1963961012123931*cMOther[2]*uOther[7]*mnuOther-0.175662013130736*uOther[6]*cMOther[7]*mnuOther-0.1963961012123931*uOther[2]*cMOther[7]*mnuOther-0.1963961012123931*cMOther[1]*uOther[6]*mnuOther-0.1963961012123931*uOther[1]*cMOther[6]*mnuOther-0.1963961012123931*cMOther[3]*uOther[5]*mnuOther-0.1963961012123931*uOther[3]*cMOther[5]*mnuOther-0.1963961012123931*cMOther[3]*uOther[4]*mnuOther-0.1963961012123931*uOther[3]*cMOther[4]*mnuOther-0.21957751641342*cMOther[0]*uOther[3]*mnuOther-1.31746509848052*m0rOther[3]*mnuOther-0.21957751641342*uOther[0]*cMOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther-0.21957751641342*cMOther[1]*uOther[2]*mnuOther-0.21957751641342*uOther[1]*cMOther[2]*mnuOther; 
  data->AEM_S(39,39) = (-0.3388888888888889*cMOther[9]*uOther[9]*mnuOther)-0.1309307341415954*cMOther[2]*uOther[9]*mnuOther-0.1309307341415954*uOther[2]*cMOther[9]*mnuOther-0.25*cMOther[8]*uOther[8]*mnuOther-0.3452380952380952*cMOther[7]*uOther[7]*mnuOther-0.149071198499986*cMOther[1]*uOther[7]*mnuOther-0.149071198499986*uOther[1]*cMOther[7]*mnuOther-0.3833333333333334*cMOther[6]*uOther[6]*mnuOther-0.3452380952380952*cMOther[5]*uOther[5]*mnuOther-0.149071198499986*cMOther[0]*uOther[5]*mnuOther-0.8944271909999159*m0rOther[5]*mnuOther-0.149071198499986*uOther[0]*cMOther[5]*mnuOther+0.2981423969999719*cEOther[5]*mnuOther-0.25*cMOther[4]*uOther[4]*mnuOther-0.3833333333333334*cMOther[3]*uOther[3]*mnuOther-0.3833333333333334*cMOther[2]*uOther[2]*mnuOther-0.25*cMOther[1]*uOther[1]*mnuOther-0.25*cMOther[0]*uOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[10]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinESelf[0] += 0.5*m1rSelf[a0+9]*uSelf[a0+9]+0.5*m1rSelf[a0+8]*uSelf[a0+8]+0.5*m1rSelf[a0+7]*uSelf[a0+7]+0.5*m1rSelf[a0+6]*uSelf[a0+6]+0.5*m1rSelf[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0+3]*uSelf[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.4391550328268398*m1rSelf[a0+4]*uSelf[a0+8]+0.4391550328268398*uSelf[a0+4]*m1rSelf[a0+8]+0.5000000000000001*m1rSelf[a0+5]*uSelf[a0+7]+0.5000000000000001*uSelf[a0+5]*m1rSelf[a0+7]+0.447213595499958*m1rSelf[a0+3]*uSelf[a0+6]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+6]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+1]*m1rSelf[a0+4]+0.5*m1rSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.4391550328268398*m1rSelf[a0+5]*uSelf[a0+9]+0.4391550328268398*uSelf[a0+5]*m1rSelf[a0+9]+0.447213595499958*m1rSelf[a0+3]*uSelf[a0+7]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+7]+0.5000000000000001*m1rSelf[a0+4]*uSelf[a0+6]+0.5000000000000001*uSelf[a0+4]*m1rSelf[a0+6]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+2]*m1rSelf[a0+5]+0.5*m1rSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]; 
    kinESelf[3] += 0.4391550328268399*m1rSelf[a0+7]*uSelf[a0+9]+0.4391550328268399*uSelf[a0+7]*m1rSelf[a0+9]+0.4391550328268399*m1rSelf[a0+6]*uSelf[a0+8]+0.4391550328268399*uSelf[a0+6]*m1rSelf[a0+8]+0.4*m1rSelf[a0+6]*uSelf[a0+7]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+7]+0.4*uSelf[a0+6]*m1rSelf[a0+7]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+7]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+6]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+6]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]; 
    kinESelf[4] += 0.2981423969999719*m1rSelf[a0+8]*uSelf[a0+8]+0.4391550328268398*m1rSelf[a0+1]*uSelf[a0+8]+0.4391550328268398*uSelf[a0+1]*m1rSelf[a0+8]+0.4472135954999579*m1rSelf[a0+7]*uSelf[a0+7]+0.31943828249997*m1rSelf[a0+6]*uSelf[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+6]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+4]+0.5*uSelf[a0]*m1rSelf[a0+4]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+1]; 
    kinESelf[5] += 0.2981423969999719*m1rSelf[a0+9]*uSelf[a0+9]+0.4391550328268398*m1rSelf[a0+2]*uSelf[a0+9]+0.4391550328268398*uSelf[a0+2]*m1rSelf[a0+9]+0.31943828249997*m1rSelf[a0+7]*uSelf[a0+7]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+7]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+7]+0.4472135954999579*m1rSelf[a0+6]*uSelf[a0+6]+0.31943828249997*m1rSelf[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0]*uSelf[a0+5]+0.5*uSelf[a0]*m1rSelf[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+2]; 
    kinESelf[6] += 0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+8]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+8]+0.4*m1rSelf[a0+3]*uSelf[a0+7]+0.4*uSelf[a0+3]*m1rSelf[a0+7]+0.4472135954999579*m1rSelf[a0+5]*uSelf[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+6]+0.5*m1rSelf[a0]*uSelf[a0+6]+0.4472135954999579*uSelf[a0+5]*m1rSelf[a0+6]+0.31943828249997*uSelf[a0+4]*m1rSelf[a0+6]+0.5*uSelf[a0]*m1rSelf[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+4]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+4]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+3]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+3]; 
    kinESelf[7] += 0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+9]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+9]+0.31943828249997*m1rSelf[a0+5]*uSelf[a0+7]+0.4472135954999579*m1rSelf[a0+4]*uSelf[a0+7]+0.5*m1rSelf[a0]*uSelf[a0+7]+0.31943828249997*uSelf[a0+5]*m1rSelf[a0+7]+0.4472135954999579*uSelf[a0+4]*m1rSelf[a0+7]+0.5*uSelf[a0]*m1rSelf[a0+7]+0.4*m1rSelf[a0+3]*uSelf[a0+6]+0.4*uSelf[a0+3]*m1rSelf[a0+6]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+5]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+5]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+3]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+3]; 
    kinESelf[8] += 0.2981423969999719*m1rSelf[a0+4]*uSelf[a0+8]+0.5*m1rSelf[a0]*uSelf[a0+8]+0.2981423969999719*uSelf[a0+4]*m1rSelf[a0+8]+0.5*uSelf[a0]*m1rSelf[a0+8]+0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+6]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+6]+0.4391550328268398*m1rSelf[a0+1]*uSelf[a0+4]+0.4391550328268398*uSelf[a0+1]*m1rSelf[a0+4]; 
    kinESelf[9] += 0.2981423969999719*m1rSelf[a0+5]*uSelf[a0+9]+0.5*m1rSelf[a0]*uSelf[a0+9]+0.2981423969999719*uSelf[a0+5]*m1rSelf[a0+9]+0.5*uSelf[a0]*m1rSelf[a0+9]+0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+7]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+7]+0.4391550328268398*m1rSelf[a0+2]*uSelf[a0+5]+0.4391550328268398*uSelf[a0+2]*m1rSelf[a0+5]; 
  } 
 
  double kinEOther[10]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    kinEOther[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.5*m1rOther[a0+9]*uOther[a0+9]+0.5*m1rOther[a0+8]*uOther[a0+8]+0.5*m1rOther[a0+7]*uOther[a0+7]+0.5*m1rOther[a0+6]*uOther[a0+6]+0.5*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.4391550328268398*m1rOther[a0+4]*uOther[a0+8]+0.4391550328268398*uOther[a0+4]*m1rOther[a0+8]+0.5000000000000001*m1rOther[a0+5]*uOther[a0+7]+0.5000000000000001*uOther[a0+5]*m1rOther[a0+7]+0.447213595499958*m1rOther[a0+3]*uOther[a0+6]+0.447213595499958*uOther[a0+3]*m1rOther[a0+6]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+4]+0.4472135954999579*uOther[a0+1]*m1rOther[a0+4]+0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.4391550328268398*m1rOther[a0+5]*uOther[a0+9]+0.4391550328268398*uOther[a0+5]*m1rOther[a0+9]+0.447213595499958*m1rOther[a0+3]*uOther[a0+7]+0.447213595499958*uOther[a0+3]*m1rOther[a0+7]+0.5000000000000001*m1rOther[a0+4]*uOther[a0+6]+0.5000000000000001*uOther[a0+4]*m1rOther[a0+6]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+5]+0.4472135954999579*uOther[a0+2]*m1rOther[a0+5]+0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    kinEOther[3] += 0.4391550328268399*m1rOther[a0+7]*uOther[a0+9]+0.4391550328268399*uOther[a0+7]*m1rOther[a0+9]+0.4391550328268399*m1rOther[a0+6]*uOther[a0+8]+0.4391550328268399*uOther[a0+6]*m1rOther[a0+8]+0.4*m1rOther[a0+6]*uOther[a0+7]+0.447213595499958*m1rOther[a0+2]*uOther[a0+7]+0.4*uOther[a0+6]*m1rOther[a0+7]+0.447213595499958*uOther[a0+2]*m1rOther[a0+7]+0.447213595499958*m1rOther[a0+1]*uOther[a0+6]+0.447213595499958*uOther[a0+1]*m1rOther[a0+6]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+5]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+4]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
    kinEOther[4] += 0.2981423969999719*m1rOther[a0+8]*uOther[a0+8]+0.4391550328268398*m1rOther[a0+1]*uOther[a0+8]+0.4391550328268398*uOther[a0+1]*m1rOther[a0+8]+0.4472135954999579*m1rOther[a0+7]*uOther[a0+7]+0.31943828249997*m1rOther[a0+6]*uOther[a0+6]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+6]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+6]+0.31943828249997*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+4]+0.5*uOther[a0]*m1rOther[a0+4]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+1]; 
    kinEOther[5] += 0.2981423969999719*m1rOther[a0+9]*uOther[a0+9]+0.4391550328268398*m1rOther[a0+2]*uOther[a0+9]+0.4391550328268398*uOther[a0+2]*m1rOther[a0+9]+0.31943828249997*m1rOther[a0+7]*uOther[a0+7]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+7]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+7]+0.4472135954999579*m1rOther[a0+6]*uOther[a0+6]+0.31943828249997*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rOther[a0]*uOther[a0+5]+0.5*uOther[a0]*m1rOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+2]; 
    kinEOther[6] += 0.4391550328268399*m1rOther[a0+3]*uOther[a0+8]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+8]+0.4*m1rOther[a0+3]*uOther[a0+7]+0.4*uOther[a0+3]*m1rOther[a0+7]+0.4472135954999579*m1rOther[a0+5]*uOther[a0+6]+0.31943828249997*m1rOther[a0+4]*uOther[a0+6]+0.5*m1rOther[a0]*uOther[a0+6]+0.4472135954999579*uOther[a0+5]*m1rOther[a0+6]+0.31943828249997*uOther[a0+4]*m1rOther[a0+6]+0.5*uOther[a0]*m1rOther[a0+6]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+4]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+4]+0.447213595499958*m1rOther[a0+1]*uOther[a0+3]+0.447213595499958*uOther[a0+1]*m1rOther[a0+3]; 
    kinEOther[7] += 0.4391550328268399*m1rOther[a0+3]*uOther[a0+9]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+9]+0.31943828249997*m1rOther[a0+5]*uOther[a0+7]+0.4472135954999579*m1rOther[a0+4]*uOther[a0+7]+0.5*m1rOther[a0]*uOther[a0+7]+0.31943828249997*uOther[a0+5]*m1rOther[a0+7]+0.4472135954999579*uOther[a0+4]*m1rOther[a0+7]+0.5*uOther[a0]*m1rOther[a0+7]+0.4*m1rOther[a0+3]*uOther[a0+6]+0.4*uOther[a0+3]*m1rOther[a0+6]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+5]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+5]+0.447213595499958*m1rOther[a0+2]*uOther[a0+3]+0.447213595499958*uOther[a0+2]*m1rOther[a0+3]; 
    kinEOther[8] += 0.2981423969999719*m1rOther[a0+4]*uOther[a0+8]+0.5*m1rOther[a0]*uOther[a0+8]+0.2981423969999719*uOther[a0+4]*m1rOther[a0+8]+0.5*uOther[a0]*m1rOther[a0+8]+0.4391550328268399*m1rOther[a0+3]*uOther[a0+6]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+6]+0.4391550328268398*m1rOther[a0+1]*uOther[a0+4]+0.4391550328268398*uOther[a0+1]*m1rOther[a0+4]; 
    kinEOther[9] += 0.2981423969999719*m1rOther[a0+5]*uOther[a0+9]+0.5*m1rOther[a0]*uOther[a0+9]+0.2981423969999719*uOther[a0+5]*m1rOther[a0+9]+0.5*uOther[a0]*m1rOther[a0+9]+0.4391550328268399*m1rOther[a0+3]*uOther[a0+7]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+7]+0.4391550328268398*m1rOther[a0+2]*uOther[a0+5]+0.4391550328268398*uOther[a0+2]*m1rOther[a0+5]; 
  } 
 
  double relKinE[10]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.5*m1rSelf[a0+9]*uSelf[a0+9]-0.5*m1rOther[a0+9]*uSelf[a0+9]-0.5*m1rSelf[a0+9]*uOther[a0+9]+0.5*m1rOther[a0+9]*uOther[a0+9]+0.5*m1rSelf[a0+8]*uSelf[a0+8]-0.5*m1rOther[a0+8]*uSelf[a0+8]-0.5*m1rSelf[a0+8]*uOther[a0+8]+0.5*m1rOther[a0+8]*uOther[a0+8]+0.5*m1rSelf[a0+7]*uSelf[a0+7]-0.5*m1rOther[a0+7]*uSelf[a0+7]-0.5*m1rSelf[a0+7]*uOther[a0+7]+0.5*m1rOther[a0+7]*uOther[a0+7]+0.5*m1rSelf[a0+6]*uSelf[a0+6]-0.5*m1rOther[a0+6]*uSelf[a0+6]-0.5*m1rSelf[a0+6]*uOther[a0+6]+0.5*m1rOther[a0+6]*uOther[a0+6]+0.5*m1rSelf[a0+5]*uSelf[a0+5]-0.5*m1rOther[a0+5]*uSelf[a0+5]-0.5*m1rSelf[a0+5]*uOther[a0+5]+0.5*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rSelf[a0+4]*uSelf[a0+4]-0.5*m1rOther[a0+4]*uSelf[a0+4]-0.5*m1rSelf[a0+4]*uOther[a0+4]+0.5*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rSelf[a0+3]*uSelf[a0+3]-0.5*m1rOther[a0+3]*uSelf[a0+3]-0.5*m1rSelf[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]-0.5*m1rOther[a0+2]*uSelf[a0+2]-0.5*m1rSelf[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]-0.5*m1rOther[a0+1]*uSelf[a0+1]-0.5*m1rSelf[a0+1]*uOther[a0+1]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]-0.5*m1rOther[a0]*uSelf[a0]-0.5*m1rSelf[a0]*uOther[a0]+0.5*m1rOther[a0]*uOther[a0]; 
    relKinE[1] += 0.4391550328268398*m1rSelf[a0+4]*uSelf[a0+8]-0.4391550328268398*m1rOther[a0+4]*uSelf[a0+8]-0.4391550328268398*m1rSelf[a0+4]*uOther[a0+8]+0.4391550328268398*m1rOther[a0+4]*uOther[a0+8]+0.4391550328268398*uSelf[a0+4]*m1rSelf[a0+8]-0.4391550328268398*uOther[a0+4]*m1rSelf[a0+8]-0.4391550328268398*uSelf[a0+4]*m1rOther[a0+8]+0.4391550328268398*uOther[a0+4]*m1rOther[a0+8]+0.5000000000000001*m1rSelf[a0+5]*uSelf[a0+7]-0.5000000000000001*m1rOther[a0+5]*uSelf[a0+7]-0.5000000000000001*m1rSelf[a0+5]*uOther[a0+7]+0.5000000000000001*m1rOther[a0+5]*uOther[a0+7]+0.5000000000000001*uSelf[a0+5]*m1rSelf[a0+7]-0.5000000000000001*uOther[a0+5]*m1rSelf[a0+7]-0.5000000000000001*uSelf[a0+5]*m1rOther[a0+7]+0.5000000000000001*uOther[a0+5]*m1rOther[a0+7]+0.447213595499958*m1rSelf[a0+3]*uSelf[a0+6]-0.447213595499958*m1rOther[a0+3]*uSelf[a0+6]-0.447213595499958*m1rSelf[a0+3]*uOther[a0+6]+0.447213595499958*m1rOther[a0+3]*uOther[a0+6]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+6]-0.447213595499958*uOther[a0+3]*m1rSelf[a0+6]-0.447213595499958*uSelf[a0+3]*m1rOther[a0+6]+0.447213595499958*uOther[a0+3]*m1rOther[a0+6]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+4]-0.4472135954999579*m1rOther[a0+1]*uSelf[a0+4]-0.4472135954999579*m1rSelf[a0+1]*uOther[a0+4]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+4]+0.4472135954999579*uSelf[a0+1]*m1rSelf[a0+4]-0.4472135954999579*uOther[a0+1]*m1rSelf[a0+4]-0.4472135954999579*uSelf[a0+1]*m1rOther[a0+4]+0.4472135954999579*uOther[a0+1]*m1rOther[a0+4]+0.5*m1rSelf[a0+2]*uSelf[a0+3]-0.5*m1rOther[a0+2]*uSelf[a0+3]-0.5*m1rSelf[a0+2]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]-0.5*uOther[a0+2]*m1rSelf[a0+3]-0.5*uSelf[a0+2]*m1rOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]-0.5*m1rOther[a0]*uSelf[a0+1]-0.5*m1rSelf[a0]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]-0.5*uOther[a0]*m1rSelf[a0+1]-0.5*uSelf[a0]*m1rOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    relKinE[2] += 0.4391550328268398*m1rSelf[a0+5]*uSelf[a0+9]-0.4391550328268398*m1rOther[a0+5]*uSelf[a0+9]-0.4391550328268398*m1rSelf[a0+5]*uOther[a0+9]+0.4391550328268398*m1rOther[a0+5]*uOther[a0+9]+0.4391550328268398*uSelf[a0+5]*m1rSelf[a0+9]-0.4391550328268398*uOther[a0+5]*m1rSelf[a0+9]-0.4391550328268398*uSelf[a0+5]*m1rOther[a0+9]+0.4391550328268398*uOther[a0+5]*m1rOther[a0+9]+0.447213595499958*m1rSelf[a0+3]*uSelf[a0+7]-0.447213595499958*m1rOther[a0+3]*uSelf[a0+7]-0.447213595499958*m1rSelf[a0+3]*uOther[a0+7]+0.447213595499958*m1rOther[a0+3]*uOther[a0+7]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+7]-0.447213595499958*uOther[a0+3]*m1rSelf[a0+7]-0.447213595499958*uSelf[a0+3]*m1rOther[a0+7]+0.447213595499958*uOther[a0+3]*m1rOther[a0+7]+0.5000000000000001*m1rSelf[a0+4]*uSelf[a0+6]-0.5000000000000001*m1rOther[a0+4]*uSelf[a0+6]-0.5000000000000001*m1rSelf[a0+4]*uOther[a0+6]+0.5000000000000001*m1rOther[a0+4]*uOther[a0+6]+0.5000000000000001*uSelf[a0+4]*m1rSelf[a0+6]-0.5000000000000001*uOther[a0+4]*m1rSelf[a0+6]-0.5000000000000001*uSelf[a0+4]*m1rOther[a0+6]+0.5000000000000001*uOther[a0+4]*m1rOther[a0+6]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+5]-0.4472135954999579*m1rOther[a0+2]*uSelf[a0+5]-0.4472135954999579*m1rSelf[a0+2]*uOther[a0+5]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+5]+0.4472135954999579*uSelf[a0+2]*m1rSelf[a0+5]-0.4472135954999579*uOther[a0+2]*m1rSelf[a0+5]-0.4472135954999579*uSelf[a0+2]*m1rOther[a0+5]+0.4472135954999579*uOther[a0+2]*m1rOther[a0+5]+0.5*m1rSelf[a0+1]*uSelf[a0+3]-0.5*m1rOther[a0+1]*uSelf[a0+3]-0.5*m1rSelf[a0+1]*uOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]-0.5*uOther[a0+1]*m1rSelf[a0+3]-0.5*uSelf[a0+1]*m1rOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]-0.5*m1rOther[a0]*uSelf[a0+2]-0.5*m1rSelf[a0]*uOther[a0+2]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]-0.5*uOther[a0]*m1rSelf[a0+2]-0.5*uSelf[a0]*m1rOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    relKinE[3] += 0.4391550328268399*m1rSelf[a0+7]*uSelf[a0+9]-0.4391550328268399*m1rOther[a0+7]*uSelf[a0+9]-0.4391550328268399*m1rSelf[a0+7]*uOther[a0+9]+0.4391550328268399*m1rOther[a0+7]*uOther[a0+9]+0.4391550328268399*uSelf[a0+7]*m1rSelf[a0+9]-0.4391550328268399*uOther[a0+7]*m1rSelf[a0+9]-0.4391550328268399*uSelf[a0+7]*m1rOther[a0+9]+0.4391550328268399*uOther[a0+7]*m1rOther[a0+9]+0.4391550328268399*m1rSelf[a0+6]*uSelf[a0+8]-0.4391550328268399*m1rOther[a0+6]*uSelf[a0+8]-0.4391550328268399*m1rSelf[a0+6]*uOther[a0+8]+0.4391550328268399*m1rOther[a0+6]*uOther[a0+8]+0.4391550328268399*uSelf[a0+6]*m1rSelf[a0+8]-0.4391550328268399*uOther[a0+6]*m1rSelf[a0+8]-0.4391550328268399*uSelf[a0+6]*m1rOther[a0+8]+0.4391550328268399*uOther[a0+6]*m1rOther[a0+8]+0.4*m1rSelf[a0+6]*uSelf[a0+7]-0.4*m1rOther[a0+6]*uSelf[a0+7]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+7]-0.447213595499958*m1rOther[a0+2]*uSelf[a0+7]-0.4*m1rSelf[a0+6]*uOther[a0+7]+0.4*m1rOther[a0+6]*uOther[a0+7]-0.447213595499958*m1rSelf[a0+2]*uOther[a0+7]+0.447213595499958*m1rOther[a0+2]*uOther[a0+7]+0.4*uSelf[a0+6]*m1rSelf[a0+7]-0.4*uOther[a0+6]*m1rSelf[a0+7]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+7]-0.447213595499958*uOther[a0+2]*m1rSelf[a0+7]-0.4*uSelf[a0+6]*m1rOther[a0+7]+0.4*uOther[a0+6]*m1rOther[a0+7]-0.447213595499958*uSelf[a0+2]*m1rOther[a0+7]+0.447213595499958*uOther[a0+2]*m1rOther[a0+7]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+6]-0.447213595499958*m1rOther[a0+1]*uSelf[a0+6]-0.447213595499958*m1rSelf[a0+1]*uOther[a0+6]+0.447213595499958*m1rOther[a0+1]*uOther[a0+6]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+6]-0.447213595499958*uOther[a0+1]*m1rSelf[a0+6]-0.447213595499958*uSelf[a0+1]*m1rOther[a0+6]+0.447213595499958*uOther[a0+1]*m1rOther[a0+6]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+5]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+5]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+5]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+5]-0.4472135954999579*uOther[a0+3]*m1rSelf[a0+5]-0.4472135954999579*uSelf[a0+3]*m1rOther[a0+5]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+4]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+4]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+4]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+4]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+4]-0.4472135954999579*uOther[a0+3]*m1rSelf[a0+4]-0.4472135954999579*uSelf[a0+3]*m1rOther[a0+4]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+3]-0.5*m1rOther[a0]*uSelf[a0+3]-0.5*m1rSelf[a0]*uOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]-0.5*uOther[a0]*m1rSelf[a0+3]-0.5*uSelf[a0]*m1rOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]-0.5*m1rOther[a0+1]*uSelf[a0+2]-0.5*m1rSelf[a0+1]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]-0.5*uOther[a0+1]*m1rSelf[a0+2]-0.5*uSelf[a0+1]*m1rOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
    relKinE[4] += 0.2981423969999719*m1rSelf[a0+8]*uSelf[a0+8]-0.2981423969999719*m1rOther[a0+8]*uSelf[a0+8]+0.4391550328268398*m1rSelf[a0+1]*uSelf[a0+8]-0.4391550328268398*m1rOther[a0+1]*uSelf[a0+8]-0.2981423969999719*m1rSelf[a0+8]*uOther[a0+8]+0.2981423969999719*m1rOther[a0+8]*uOther[a0+8]-0.4391550328268398*m1rSelf[a0+1]*uOther[a0+8]+0.4391550328268398*m1rOther[a0+1]*uOther[a0+8]+0.4391550328268398*uSelf[a0+1]*m1rSelf[a0+8]-0.4391550328268398*uOther[a0+1]*m1rSelf[a0+8]-0.4391550328268398*uSelf[a0+1]*m1rOther[a0+8]+0.4391550328268398*uOther[a0+1]*m1rOther[a0+8]+0.4472135954999579*m1rSelf[a0+7]*uSelf[a0+7]-0.4472135954999579*m1rOther[a0+7]*uSelf[a0+7]-0.4472135954999579*m1rSelf[a0+7]*uOther[a0+7]+0.4472135954999579*m1rOther[a0+7]*uOther[a0+7]+0.31943828249997*m1rSelf[a0+6]*uSelf[a0+6]-0.31943828249997*m1rOther[a0+6]*uSelf[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+6]-0.5000000000000001*m1rOther[a0+2]*uSelf[a0+6]-0.31943828249997*m1rSelf[a0+6]*uOther[a0+6]+0.31943828249997*m1rOther[a0+6]*uOther[a0+6]-0.5000000000000001*m1rSelf[a0+2]*uOther[a0+6]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+6]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+6]-0.5000000000000001*uOther[a0+2]*m1rSelf[a0+6]-0.5000000000000001*uSelf[a0+2]*m1rOther[a0+6]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+4]-0.31943828249997*m1rOther[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+4]-0.5*m1rOther[a0]*uSelf[a0+4]-0.31943828249997*m1rSelf[a0+4]*uOther[a0+4]+0.31943828249997*m1rOther[a0+4]*uOther[a0+4]-0.5*m1rSelf[a0]*uOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+4]+0.5*uSelf[a0]*m1rSelf[a0+4]-0.5*uOther[a0]*m1rSelf[a0+4]-0.5*uSelf[a0]*m1rOther[a0+4]+0.5*uOther[a0]*m1rOther[a0+4]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+1]-0.4472135954999579*m1rOther[a0+1]*uSelf[a0+1]-0.4472135954999579*m1rSelf[a0+1]*uOther[a0+1]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+1]; 
    relKinE[5] += 0.2981423969999719*m1rSelf[a0+9]*uSelf[a0+9]-0.2981423969999719*m1rOther[a0+9]*uSelf[a0+9]+0.4391550328268398*m1rSelf[a0+2]*uSelf[a0+9]-0.4391550328268398*m1rOther[a0+2]*uSelf[a0+9]-0.2981423969999719*m1rSelf[a0+9]*uOther[a0+9]+0.2981423969999719*m1rOther[a0+9]*uOther[a0+9]-0.4391550328268398*m1rSelf[a0+2]*uOther[a0+9]+0.4391550328268398*m1rOther[a0+2]*uOther[a0+9]+0.4391550328268398*uSelf[a0+2]*m1rSelf[a0+9]-0.4391550328268398*uOther[a0+2]*m1rSelf[a0+9]-0.4391550328268398*uSelf[a0+2]*m1rOther[a0+9]+0.4391550328268398*uOther[a0+2]*m1rOther[a0+9]+0.31943828249997*m1rSelf[a0+7]*uSelf[a0+7]-0.31943828249997*m1rOther[a0+7]*uSelf[a0+7]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+7]-0.5000000000000001*m1rOther[a0+1]*uSelf[a0+7]-0.31943828249997*m1rSelf[a0+7]*uOther[a0+7]+0.31943828249997*m1rOther[a0+7]*uOther[a0+7]-0.5000000000000001*m1rSelf[a0+1]*uOther[a0+7]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+7]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+7]-0.5000000000000001*uOther[a0+1]*m1rSelf[a0+7]-0.5000000000000001*uSelf[a0+1]*m1rOther[a0+7]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+7]+0.4472135954999579*m1rSelf[a0+6]*uSelf[a0+6]-0.4472135954999579*m1rOther[a0+6]*uSelf[a0+6]-0.4472135954999579*m1rSelf[a0+6]*uOther[a0+6]+0.4472135954999579*m1rOther[a0+6]*uOther[a0+6]+0.31943828249997*m1rSelf[a0+5]*uSelf[a0+5]-0.31943828249997*m1rOther[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0]*uSelf[a0+5]-0.5*m1rOther[a0]*uSelf[a0+5]-0.31943828249997*m1rSelf[a0+5]*uOther[a0+5]+0.31943828249997*m1rOther[a0+5]*uOther[a0+5]-0.5*m1rSelf[a0]*uOther[a0+5]+0.5*m1rOther[a0]*uOther[a0+5]+0.5*uSelf[a0]*m1rSelf[a0+5]-0.5*uOther[a0]*m1rSelf[a0+5]-0.5*uSelf[a0]*m1rOther[a0+5]+0.5*uOther[a0]*m1rOther[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+2]-0.4472135954999579*m1rOther[a0+2]*uSelf[a0+2]-0.4472135954999579*m1rSelf[a0+2]*uOther[a0+2]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+2]; 
    relKinE[6] += 0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+8]-0.4391550328268399*m1rOther[a0+3]*uSelf[a0+8]-0.4391550328268399*m1rSelf[a0+3]*uOther[a0+8]+0.4391550328268399*m1rOther[a0+3]*uOther[a0+8]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+8]-0.4391550328268399*uOther[a0+3]*m1rSelf[a0+8]-0.4391550328268399*uSelf[a0+3]*m1rOther[a0+8]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+8]+0.4*m1rSelf[a0+3]*uSelf[a0+7]-0.4*m1rOther[a0+3]*uSelf[a0+7]-0.4*m1rSelf[a0+3]*uOther[a0+7]+0.4*m1rOther[a0+3]*uOther[a0+7]+0.4*uSelf[a0+3]*m1rSelf[a0+7]-0.4*uOther[a0+3]*m1rSelf[a0+7]-0.4*uSelf[a0+3]*m1rOther[a0+7]+0.4*uOther[a0+3]*m1rOther[a0+7]+0.4472135954999579*m1rSelf[a0+5]*uSelf[a0+6]-0.4472135954999579*m1rOther[a0+5]*uSelf[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+6]-0.31943828249997*m1rOther[a0+4]*uSelf[a0+6]+0.5*m1rSelf[a0]*uSelf[a0+6]-0.5*m1rOther[a0]*uSelf[a0+6]-0.4472135954999579*m1rSelf[a0+5]*uOther[a0+6]+0.4472135954999579*m1rOther[a0+5]*uOther[a0+6]-0.31943828249997*m1rSelf[a0+4]*uOther[a0+6]+0.31943828249997*m1rOther[a0+4]*uOther[a0+6]-0.5*m1rSelf[a0]*uOther[a0+6]+0.5*m1rOther[a0]*uOther[a0+6]+0.4472135954999579*uSelf[a0+5]*m1rSelf[a0+6]-0.4472135954999579*uOther[a0+5]*m1rSelf[a0+6]+0.31943828249997*uSelf[a0+4]*m1rSelf[a0+6]-0.31943828249997*uOther[a0+4]*m1rSelf[a0+6]+0.5*uSelf[a0]*m1rSelf[a0+6]-0.5*uOther[a0]*m1rSelf[a0+6]-0.4472135954999579*uSelf[a0+5]*m1rOther[a0+6]+0.4472135954999579*uOther[a0+5]*m1rOther[a0+6]-0.31943828249997*uSelf[a0+4]*m1rOther[a0+6]+0.31943828249997*uOther[a0+4]*m1rOther[a0+6]-0.5*uSelf[a0]*m1rOther[a0+6]+0.5*uOther[a0]*m1rOther[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+4]-0.5000000000000001*m1rOther[a0+2]*uSelf[a0+4]-0.5000000000000001*m1rSelf[a0+2]*uOther[a0+4]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+4]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+4]-0.5000000000000001*uOther[a0+2]*m1rSelf[a0+4]-0.5000000000000001*uSelf[a0+2]*m1rOther[a0+4]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+4]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+3]-0.447213595499958*m1rOther[a0+1]*uSelf[a0+3]-0.447213595499958*m1rSelf[a0+1]*uOther[a0+3]+0.447213595499958*m1rOther[a0+1]*uOther[a0+3]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+3]-0.447213595499958*uOther[a0+1]*m1rSelf[a0+3]-0.447213595499958*uSelf[a0+1]*m1rOther[a0+3]+0.447213595499958*uOther[a0+1]*m1rOther[a0+3]; 
    relKinE[7] += 0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+9]-0.4391550328268399*m1rOther[a0+3]*uSelf[a0+9]-0.4391550328268399*m1rSelf[a0+3]*uOther[a0+9]+0.4391550328268399*m1rOther[a0+3]*uOther[a0+9]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+9]-0.4391550328268399*uOther[a0+3]*m1rSelf[a0+9]-0.4391550328268399*uSelf[a0+3]*m1rOther[a0+9]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+9]+0.31943828249997*m1rSelf[a0+5]*uSelf[a0+7]-0.31943828249997*m1rOther[a0+5]*uSelf[a0+7]+0.4472135954999579*m1rSelf[a0+4]*uSelf[a0+7]-0.4472135954999579*m1rOther[a0+4]*uSelf[a0+7]+0.5*m1rSelf[a0]*uSelf[a0+7]-0.5*m1rOther[a0]*uSelf[a0+7]-0.31943828249997*m1rSelf[a0+5]*uOther[a0+7]+0.31943828249997*m1rOther[a0+5]*uOther[a0+7]-0.4472135954999579*m1rSelf[a0+4]*uOther[a0+7]+0.4472135954999579*m1rOther[a0+4]*uOther[a0+7]-0.5*m1rSelf[a0]*uOther[a0+7]+0.5*m1rOther[a0]*uOther[a0+7]+0.31943828249997*uSelf[a0+5]*m1rSelf[a0+7]-0.31943828249997*uOther[a0+5]*m1rSelf[a0+7]+0.4472135954999579*uSelf[a0+4]*m1rSelf[a0+7]-0.4472135954999579*uOther[a0+4]*m1rSelf[a0+7]+0.5*uSelf[a0]*m1rSelf[a0+7]-0.5*uOther[a0]*m1rSelf[a0+7]-0.31943828249997*uSelf[a0+5]*m1rOther[a0+7]+0.31943828249997*uOther[a0+5]*m1rOther[a0+7]-0.4472135954999579*uSelf[a0+4]*m1rOther[a0+7]+0.4472135954999579*uOther[a0+4]*m1rOther[a0+7]-0.5*uSelf[a0]*m1rOther[a0+7]+0.5*uOther[a0]*m1rOther[a0+7]+0.4*m1rSelf[a0+3]*uSelf[a0+6]-0.4*m1rOther[a0+3]*uSelf[a0+6]-0.4*m1rSelf[a0+3]*uOther[a0+6]+0.4*m1rOther[a0+3]*uOther[a0+6]+0.4*uSelf[a0+3]*m1rSelf[a0+6]-0.4*uOther[a0+3]*m1rSelf[a0+6]-0.4*uSelf[a0+3]*m1rOther[a0+6]+0.4*uOther[a0+3]*m1rOther[a0+6]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+5]-0.5000000000000001*m1rOther[a0+1]*uSelf[a0+5]-0.5000000000000001*m1rSelf[a0+1]*uOther[a0+5]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+5]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+5]-0.5000000000000001*uOther[a0+1]*m1rSelf[a0+5]-0.5000000000000001*uSelf[a0+1]*m1rOther[a0+5]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+5]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+3]-0.447213595499958*m1rOther[a0+2]*uSelf[a0+3]-0.447213595499958*m1rSelf[a0+2]*uOther[a0+3]+0.447213595499958*m1rOther[a0+2]*uOther[a0+3]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+3]-0.447213595499958*uOther[a0+2]*m1rSelf[a0+3]-0.447213595499958*uSelf[a0+2]*m1rOther[a0+3]+0.447213595499958*uOther[a0+2]*m1rOther[a0+3]; 
    relKinE[8] += 0.2981423969999719*m1rSelf[a0+4]*uSelf[a0+8]-0.2981423969999719*m1rOther[a0+4]*uSelf[a0+8]+0.5*m1rSelf[a0]*uSelf[a0+8]-0.5*m1rOther[a0]*uSelf[a0+8]-0.2981423969999719*m1rSelf[a0+4]*uOther[a0+8]+0.2981423969999719*m1rOther[a0+4]*uOther[a0+8]-0.5*m1rSelf[a0]*uOther[a0+8]+0.5*m1rOther[a0]*uOther[a0+8]+0.2981423969999719*uSelf[a0+4]*m1rSelf[a0+8]-0.2981423969999719*uOther[a0+4]*m1rSelf[a0+8]+0.5*uSelf[a0]*m1rSelf[a0+8]-0.5*uOther[a0]*m1rSelf[a0+8]-0.2981423969999719*uSelf[a0+4]*m1rOther[a0+8]+0.2981423969999719*uOther[a0+4]*m1rOther[a0+8]-0.5*uSelf[a0]*m1rOther[a0+8]+0.5*uOther[a0]*m1rOther[a0+8]+0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+6]-0.4391550328268399*m1rOther[a0+3]*uSelf[a0+6]-0.4391550328268399*m1rSelf[a0+3]*uOther[a0+6]+0.4391550328268399*m1rOther[a0+3]*uOther[a0+6]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+6]-0.4391550328268399*uOther[a0+3]*m1rSelf[a0+6]-0.4391550328268399*uSelf[a0+3]*m1rOther[a0+6]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+6]+0.4391550328268398*m1rSelf[a0+1]*uSelf[a0+4]-0.4391550328268398*m1rOther[a0+1]*uSelf[a0+4]-0.4391550328268398*m1rSelf[a0+1]*uOther[a0+4]+0.4391550328268398*m1rOther[a0+1]*uOther[a0+4]+0.4391550328268398*uSelf[a0+1]*m1rSelf[a0+4]-0.4391550328268398*uOther[a0+1]*m1rSelf[a0+4]-0.4391550328268398*uSelf[a0+1]*m1rOther[a0+4]+0.4391550328268398*uOther[a0+1]*m1rOther[a0+4]; 
    relKinE[9] += 0.2981423969999719*m1rSelf[a0+5]*uSelf[a0+9]-0.2981423969999719*m1rOther[a0+5]*uSelf[a0+9]+0.5*m1rSelf[a0]*uSelf[a0+9]-0.5*m1rOther[a0]*uSelf[a0+9]-0.2981423969999719*m1rSelf[a0+5]*uOther[a0+9]+0.2981423969999719*m1rOther[a0+5]*uOther[a0+9]-0.5*m1rSelf[a0]*uOther[a0+9]+0.5*m1rOther[a0]*uOther[a0+9]+0.2981423969999719*uSelf[a0+5]*m1rSelf[a0+9]-0.2981423969999719*uOther[a0+5]*m1rSelf[a0+9]+0.5*uSelf[a0]*m1rSelf[a0+9]-0.5*uOther[a0]*m1rSelf[a0+9]-0.2981423969999719*uSelf[a0+5]*m1rOther[a0+9]+0.2981423969999719*uOther[a0+5]*m1rOther[a0+9]-0.5*uSelf[a0]*m1rOther[a0+9]+0.5*uOther[a0]*m1rOther[a0+9]+0.4391550328268399*m1rSelf[a0+3]*uSelf[a0+7]-0.4391550328268399*m1rOther[a0+3]*uSelf[a0+7]-0.4391550328268399*m1rSelf[a0+3]*uOther[a0+7]+0.4391550328268399*m1rOther[a0+3]*uOther[a0+7]+0.4391550328268399*uSelf[a0+3]*m1rSelf[a0+7]-0.4391550328268399*uOther[a0+3]*m1rSelf[a0+7]-0.4391550328268399*uSelf[a0+3]*m1rOther[a0+7]+0.4391550328268399*uOther[a0+3]*m1rOther[a0+7]+0.4391550328268398*m1rSelf[a0+2]*uSelf[a0+5]-0.4391550328268398*m1rOther[a0+2]*uSelf[a0+5]-0.4391550328268398*m1rSelf[a0+2]*uOther[a0+5]+0.4391550328268398*m1rOther[a0+2]*uOther[a0+5]+0.4391550328268398*uSelf[a0+2]*m1rSelf[a0+5]-0.4391550328268398*uOther[a0+2]*m1rSelf[a0+5]-0.4391550328268398*uSelf[a0+2]*m1rOther[a0+5]+0.4391550328268398*uOther[a0+2]*m1rOther[a0+5]; 
  } 
 
  double m2Relax[10]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(0.5*relKinE[0]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[0]*mSelf)/(mSelf+mOther)+(kinESelf[0]*mSelf)/(mSelf+mOther)+(0.5*relKinE[0]*mOther)/(mSelf+mOther)+(m2rOther[0]*mOther)/(mSelf+mOther)-(1.0*kinEOther[0]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2rOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(0.5*relKinE[1]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[1]*mSelf)/(mSelf+mOther)+(kinESelf[1]*mSelf)/(mSelf+mOther)+(0.5*relKinE[1]*mOther)/(mSelf+mOther)+(m2rOther[1]*mOther)/(mSelf+mOther)-(1.0*kinEOther[1]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2rOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(0.5*relKinE[2]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[2]*mSelf)/(mSelf+mOther)+(kinESelf[2]*mSelf)/(mSelf+mOther)+(0.5*relKinE[2]*mOther)/(mSelf+mOther)+(m2rOther[2]*mOther)/(mSelf+mOther)-(1.0*kinEOther[2]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2rOther[2])*mnuOther; 
  m2Relax[3] = betaGreenep1*((-(0.5*relKinE[3]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[3]*mSelf)/(mSelf+mOther)+(kinESelf[3]*mSelf)/(mSelf+mOther)+(0.5*relKinE[3]*mOther)/(mSelf+mOther)+(m2rOther[3]*mOther)/(mSelf+mOther)-(1.0*kinEOther[3]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[3]-1.0*kinESelf[3])*mnuSelf+(kinEOther[3]-1.0*m2rOther[3])*mnuOther; 
  m2Relax[4] = betaGreenep1*((-(0.5*relKinE[4]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[4]*mSelf)/(mSelf+mOther)+(kinESelf[4]*mSelf)/(mSelf+mOther)+(0.5*relKinE[4]*mOther)/(mSelf+mOther)+(m2rOther[4]*mOther)/(mSelf+mOther)-(1.0*kinEOther[4]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[4]-1.0*kinESelf[4])*mnuSelf+(kinEOther[4]-1.0*m2rOther[4])*mnuOther; 
  m2Relax[5] = betaGreenep1*((-(0.5*relKinE[5]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[5]*mSelf)/(mSelf+mOther)+(kinESelf[5]*mSelf)/(mSelf+mOther)+(0.5*relKinE[5]*mOther)/(mSelf+mOther)+(m2rOther[5]*mOther)/(mSelf+mOther)-(1.0*kinEOther[5]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[5]-1.0*kinESelf[5])*mnuSelf+(kinEOther[5]-1.0*m2rOther[5])*mnuOther; 
  m2Relax[6] = betaGreenep1*((-(0.5*relKinE[6]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[6]*mSelf)/(mSelf+mOther)+(kinESelf[6]*mSelf)/(mSelf+mOther)+(0.5*relKinE[6]*mOther)/(mSelf+mOther)+(m2rOther[6]*mOther)/(mSelf+mOther)-(1.0*kinEOther[6]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[6]-1.0*kinESelf[6])*mnuSelf+(kinEOther[6]-1.0*m2rOther[6])*mnuOther; 
  m2Relax[7] = betaGreenep1*((-(0.5*relKinE[7]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[7]*mSelf)/(mSelf+mOther)+(kinESelf[7]*mSelf)/(mSelf+mOther)+(0.5*relKinE[7]*mOther)/(mSelf+mOther)+(m2rOther[7]*mOther)/(mSelf+mOther)-(1.0*kinEOther[7]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[7]-1.0*kinESelf[7])*mnuSelf+(kinEOther[7]-1.0*m2rOther[7])*mnuOther; 
  m2Relax[8] = betaGreenep1*((-(0.5*relKinE[8]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[8]*mSelf)/(mSelf+mOther)+(kinESelf[8]*mSelf)/(mSelf+mOther)+(0.5*relKinE[8]*mOther)/(mSelf+mOther)+(m2rOther[8]*mOther)/(mSelf+mOther)-(1.0*kinEOther[8]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[8]-1.0*kinESelf[8])*mnuSelf+(kinEOther[8]-1.0*m2rOther[8])*mnuOther; 
  m2Relax[9] = betaGreenep1*((-(0.5*relKinE[9]*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[9]*mSelf)/(mSelf+mOther)+(kinESelf[9]*mSelf)/(mSelf+mOther)+(0.5*relKinE[9]*mOther)/(mSelf+mOther)+(m2rOther[9]*mOther)/(mSelf+mOther)-(1.0*kinEOther[9]*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[9]-1.0*kinESelf[9])*mnuSelf+(kinEOther[9]-1.0*m2rOther[9])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],mnuM2sum[6],mnuM2sum[7],mnuM2sum[8],mnuM2sum[9],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m1Relax[9],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5],m2Relax[6],m2Relax[7],m2Relax[8],m2Relax[9]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,10,1) = data->u_S.segment<10>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,10,1) = data->u_S.segment<10>(10); 
 
  Eigen::Map<VectorXd>(uCrossOther,10,1) = data->u_S.segment<10>(20); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,10,1) = data->u_S.segment<10>(30); 
 
} 
 
