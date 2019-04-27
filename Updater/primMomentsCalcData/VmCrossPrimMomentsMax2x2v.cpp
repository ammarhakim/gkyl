#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmCrossPrimMoments2x2vMax_P1(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[6]; 
  double m2rSelf[3]; 
  double m0SrSelf[3]; 
  double m1SrSelf[6]; 
  double m2SrSelf[3]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    m1rSelf[2] = 0.0; 
    m1rSelf[3] = m1Self[3]; 
    m1rSelf[4] = 0.0; 
    m1rSelf[5] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = 0.0; 
    m0SrSelf[2] = 0.0; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = 0.0; 
    m1SrSelf[2] = 0.0; 
    m1SrSelf[3] = m1SSelf[3]; 
    m1SrSelf[4] = 0.0; 
    m1SrSelf[5] = 0.0; 
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
    m1rSelf[3] = m1Self[3]; 
    m1rSelf[4] = m1Self[4]; 
    m1rSelf[5] = m1Self[5]; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = m2Self[1]; 
    m2rSelf[2] = m2Self[2]; 
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = m0SSelf[1]; 
    m0SrSelf[2] = m0SSelf[2]; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = m1SSelf[1]; 
    m1SrSelf[2] = m1SSelf[2]; 
    m1SrSelf[3] = m1SSelf[3]; 
    m1SrSelf[4] = m1SSelf[4]; 
    m1SrSelf[5] = m1SSelf[5]; 
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
  double m1rOther[6]; 
  double m2rOther[3]; 
  double m0SrOther[3]; 
  double m1SrOther[6]; 
  double m2SrOther[3]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    m1rOther[2] = 0.0; 
    m1rOther[3] = m1Other[3]; 
    m1rOther[4] = 0.0; 
    m1rOther[5] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = 0.0; 
    m0SrOther[2] = 0.0; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = 0.0; 
    m1SrOther[2] = 0.0; 
    m1SrOther[3] = m1SOther[3]; 
    m1SrOther[4] = 0.0; 
    m1SrOther[5] = 0.0; 
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
    m1rOther[3] = m1Other[3]; 
    m1rOther[4] = m1Other[4]; 
    m1rOther[5] = m1Other[5]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = m0SOther[1]; 
    m0SrOther[2] = m0SOther[2]; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = m1SOther[1]; 
    m1SrOther[2] = m1SOther[2]; 
    m1SrOther[3] = m1SOther[3]; 
    m1SrOther[4] = m1SOther[4]; 
    m1SrOther[5] = m1SOther[5]; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
    m2SrOther[2] = m2SOther[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(18,18); 
 
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
  data->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,2) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,6) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,7) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,8) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,6) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,7) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,6) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,8) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,9) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,10) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,11) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,9) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,10) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,9) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,11) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,15) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,16) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,17) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,15) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,16) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,15) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,17) = -0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(6,0) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(6,1) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(6,2) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(7,0) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(7,1) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(8,0) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(8,2) = 0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(6,9) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(6,10) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(6,11) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(7,9) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(7,10) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(8,9) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(8,11) = 0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(3,3) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(3,4) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,5) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,4) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(5,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(5,5) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(3,6) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,7) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(3,8) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(4,6) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,7) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,6) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,8) = -0.5*cMSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(3,12) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,13) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,14) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(4,12) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(4,13) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,12) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(5,14) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(3,15) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,16) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(3,17) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(4,15) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,16) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(5,15) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,17) = -0.5*cMOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(6,3) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(6,4) = 0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(6,5) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(7,3) = 0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(7,4) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(8,3) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(8,5) = 0.5*m1SrSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(6,12) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(6,13) = 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(6,14) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(7,12) = 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(7,13) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(8,12) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(8,14) = 0.5*m1SrOther[3]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(6,6) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(6,8) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(7,6) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(7,7) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(8,6) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(8,8) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(6,15) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(6,16) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(6,17) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(7,15) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(7,16) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(8,15) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(8,17) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[3]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2SrSelf[0]*mnuSelf+m2SrOther[0]*mnuOther; 
  mnuM2sum[1] = m2SrSelf[1]*mnuSelf+m2SrOther[1]*mnuOther; 
  mnuM2sum[2] = m2SrSelf[2]*mnuSelf+m2SrOther[2]*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,3>(0,3).setZero(); 
  data->AEM_S.block<3,3>(3,0).setZero(); 
  data->AEM_S.block<3,3>(0,12).setZero(); 
  data->AEM_S.block<3,3>(3,9).setZero(); 
 
  double m1Relax[6]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(9,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(9,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(10,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(11,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,2) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(9,6) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(9,7) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(9,8) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(10,6) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(10,7) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(11,6) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(11,8) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(9,9) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,10) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(9,11) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,9) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(10,10) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(11,9) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(11,11) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(9,15) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(9,16) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(9,17) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(10,15) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(10,16) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(11,15) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(11,17) = 0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(15,0) = (-0.25*m0rSelf[2]*uSelf[2]*mnuSelf)-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(15,1) = (-0.25*m0rSelf[0]*uSelf[1]*mnuSelf)+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,2) = (-0.25*m0rSelf[0]*uSelf[2]*mnuSelf)+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,0) = (-0.25*m0rSelf[0]*uSelf[1]*mnuSelf)+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,1) = (-0.25*m0rSelf[2]*uSelf[2]*mnuSelf)-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(16,2) = (-0.25*m0rSelf[1]*uSelf[2]*mnuSelf)-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,0) = (-0.25*m0rSelf[0]*uSelf[2]*mnuSelf)+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,1) = (-0.25*m0rSelf[1]*uSelf[2]*mnuSelf)-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,2) = (-0.45*m0rSelf[2]*uSelf[2]*mnuSelf)-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(15,9) = 0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(15,10) = 0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(15,11) = 0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(16,9) = 0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(16,10) = 0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(16,11) = 0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(17,9) = 0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(17,10) = 0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(17,11) = 0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(12,3) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,4) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,5) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,4) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(14,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,5) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(12,6) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(12,7) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(12,8) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(13,6) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(13,7) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(14,6) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(14,8) = -0.5*cMSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(12,12) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(12,13) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(12,14) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(13,12) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(13,13) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(14,12) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(14,14) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(12,15) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(12,16) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(12,17) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(13,15) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(13,16) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(14,15) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(14,17) = 0.5*cMOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(15,3) = (-0.25*m0rSelf[2]*uSelf[5]*mnuSelf)-0.25*m0rSelf[1]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(15,4) = (-0.25*m0rSelf[0]*uSelf[4]*mnuSelf)+0.5*m1SrSelf[4]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(15,5) = (-0.25*m0rSelf[0]*uSelf[5]*mnuSelf)+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf; 
  data->AEM_S(16,3) = (-0.25*m0rSelf[0]*uSelf[4]*mnuSelf)+0.5*m1SrSelf[4]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(16,4) = (-0.25*m0rSelf[2]*uSelf[5]*mnuSelf)-0.45*m0rSelf[1]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(16,5) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(17,3) = (-0.25*m0rSelf[0]*uSelf[5]*mnuSelf)+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf; 
  data->AEM_S(17,4) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(17,5) = (-0.45*m0rSelf[2]*uSelf[5]*mnuSelf)-0.25*m0rSelf[1]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1SrSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(15,12) = 0.25*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(15,13) = 0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(15,14) = 0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther; 
  data->AEM_S(16,12) = 0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(16,13) = 0.25*m0rOther[2]*uOther[5]*mnuOther+0.45*m0rOther[1]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(16,14) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(17,12) = 0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther; 
  data->AEM_S(17,13) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(17,14) = 0.45*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[3] += (m1rOther[3]-1.0*m1rSelf[3])*betaGreenep1*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (m1rOther[4]-1.0*m1rSelf[4])*betaGreenep1*mnuSelf+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (m1rOther[5]-1.0*m1rSelf[5])*betaGreenep1*mnuSelf+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(15,6) = 0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[3]*uSelf[3]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(15,7) = 0.25*cMSelf[3]*uSelf[4]*mnuSelf+0.25*uSelf[3]*cMSelf[4]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(15,8) = 0.25*cMSelf[3]*uSelf[5]*mnuSelf+0.25*uSelf[3]*cMSelf[5]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(16,6) = 0.25*cMSelf[3]*uSelf[4]*mnuSelf+0.25*uSelf[3]*cMSelf[4]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(16,7) = 0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[3]*uSelf[3]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(17,6) = 0.25*cMSelf[3]*uSelf[5]*mnuSelf+0.25*uSelf[3]*cMSelf[5]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(17,8) = 0.25*cMSelf[5]*uSelf[5]*mnuSelf+0.25*cMSelf[4]*uSelf[4]*mnuSelf+0.25*cMSelf[3]*uSelf[3]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(15,15) = (-0.25*cMOther[5]*uOther[5]*mnuOther)-0.25*cMOther[4]*uOther[4]*mnuOther-0.25*cMOther[3]*uOther[3]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(15,16) = (-0.25*cMOther[3]*uOther[4]*mnuOther)-0.25*uOther[3]*cMOther[4]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(15,17) = (-0.25*cMOther[3]*uOther[5]*mnuOther)-0.25*uOther[3]*cMOther[5]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(16,15) = (-0.25*cMOther[3]*uOther[4]*mnuOther)-0.25*uOther[3]*cMOther[4]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(16,16) = (-0.25*cMOther[5]*uOther[5]*mnuOther)-0.25*cMOther[4]*uOther[4]*mnuOther-0.25*cMOther[3]*uOther[3]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(17,15) = (-0.25*cMOther[3]*uOther[5]*mnuOther)-0.25*uOther[3]*cMOther[5]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(17,17) = (-0.25*cMOther[5]*uOther[5]*mnuOther)-0.25*cMOther[4]*uOther[4]*mnuOther-0.25*cMOther[3]*uOther[3]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[3]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,3>(9,3).setZero(); 
  data->AEM_S.block<3,3>(12,0).setZero(); 
  data->AEM_S.block<3,3>(9,12).setZero(); 
  data->AEM_S.block<3,3>(12,9).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m2Relax[0],m2Relax[1],m2Relax[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,6,1) = data->u_S.segment<6>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,3,1) = data->u_S.segment<3>(6); 
 
  Eigen::Map<VectorXd>(uCrossOther,6,1) = data->u_S.segment<6>(9); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,3,1) = data->u_S.segment<3>(15); 
 
} 
 
void VmCrossPrimMoments2x2vMax_P2(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[12]; 
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
    m1rSelf[6] = m1Self[6]; 
    m1rSelf[7] = 0.0; 
    m1rSelf[8] = 0.0; 
    m1rSelf[9] = 0.0; 
    m1rSelf[10] = 0.0; 
    m1rSelf[11] = 0.0; 
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
    m1rSelf[6] = m1Self[6]; 
    m1rSelf[7] = m1Self[7]; 
    m1rSelf[8] = m1Self[8]; 
    m1rSelf[9] = m1Self[9]; 
    m1rSelf[10] = m1Self[10]; 
    m1rSelf[11] = m1Self[11]; 
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
  double m1rOther[12]; 
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
    m1rOther[6] = m1Other[6]; 
    m1rOther[7] = 0.0; 
    m1rOther[8] = 0.0; 
    m1rOther[9] = 0.0; 
    m1rOther[10] = 0.0; 
    m1rOther[11] = 0.0; 
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
    m1rOther[6] = m1Other[6]; 
    m1rOther[7] = m1Other[7]; 
    m1rOther[8] = m1Other[8]; 
    m1rOther[9] = m1Other[9]; 
    m1rOther[10] = m1Other[10]; 
    m1rOther[11] = m1Other[11]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
    m2rOther[3] = m2Other[3]; 
    m2rOther[4] = m2Other[4]; 
    m2rOther[5] = m2Other[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(36,36); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
  double mnuM1sum[12]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<12; vd++) 
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
  data->AEM_S(0,12) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,13) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,14) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,15) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,16) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,17) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(1,12) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,13) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,14) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,15) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,16) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,12) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,13) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,14) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,15) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,17) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,12) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,13) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,14) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,15) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,16) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,17) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,12) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,13) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,15) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,16) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,12) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,14) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(5,15) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,17) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,18) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,19) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,20) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,21) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(0,22) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(0,23) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(1,18) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,19) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,20) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,21) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,22) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(2,18) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,19) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,20) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,21) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(2,23) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(3,18) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,19) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,20) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,21) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,22) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,23) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,18) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(4,19) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(4,21) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,22) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,18) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(5,20) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(5,21) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(5,23) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,30) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,31) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,32) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,33) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(0,34) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(0,35) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(1,30) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,31) = (-0.4472135954999579*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,32) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,33) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,34) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(2,30) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,31) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,32) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,33) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(2,35) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(3,30) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,31) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,32) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,33) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(3,34) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,35) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,30) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,31) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(4,33) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,34) = (-0.31943828249997*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(5,30) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,32) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(5,33) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(5,35) = (-0.31943828249997*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(12,0) = 0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(12,1) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(12,3) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(12,4) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(12,5) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(13,0) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(13,1) = 0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(13,2) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(13,3) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(13,4) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(14,0) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(14,1) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(14,2) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(14,3) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(14,5) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(15,0) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(15,1) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(15,2) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(15,3) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(15,4) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(15,5) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(16,0) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(16,1) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(16,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(16,4) = 0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(17,0) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(17,2) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(17,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(17,5) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(12,18) = 0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(12,19) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(12,20) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(12,21) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(12,22) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(12,23) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(13,18) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(13,19) = 0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(13,20) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(13,21) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(13,22) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(14,18) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(14,19) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(14,20) = 0.4472135954999579*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(14,21) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(14,23) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(15,18) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(15,19) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(15,20) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(15,21) = 0.4472135954999579*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(15,22) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(15,23) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(16,18) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(16,19) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(16,21) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(16,22) = 0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(17,18) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(17,20) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(17,21) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(17,23) = 0.31943828249997*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(6,6) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,9) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(6,10) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(6,11) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(7,6) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,7) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(7,8) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,9) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,10) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,6) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,7) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(8,8) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,9) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,11) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(9,6) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,7) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(9,8) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,9) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(9,10) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,11) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(10,6) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(10,7) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(10,9) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(10,10) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(11,6) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(11,8) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,9) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,11) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(6,12) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,13) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(6,14) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(6,15) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(6,16) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(6,17) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(7,12) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,13) = (-0.4472135954999579*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(7,14) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(7,15) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(7,16) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(8,12) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,13) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(8,14) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(8,15) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(8,17) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(9,12) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,13) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(9,14) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(9,15) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.4472135954999579*cMSelf[10]*mnuSelf-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(9,16) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,17) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(10,12) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,13) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(10,15) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(10,16) = (-0.31943828249997*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(11,12) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,14) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(11,15) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(11,17) = (-0.31943828249997*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(6,24) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,25) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(6,26) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(6,27) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(6,28) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(6,29) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(7,24) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,25) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(7,26) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(7,27) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(7,28) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(8,24) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,25) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(8,26) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,27) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(8,29) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(9,24) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,25) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(9,26) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(9,27) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,28) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(9,29) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(10,24) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(10,25) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(10,27) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(10,28) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(11,24) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(11,26) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(11,27) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(11,29) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(6,30) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,31) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(6,32) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(6,33) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(6,34) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(6,35) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(7,30) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,31) = (-0.4472135954999579*cMOther[10]*mnuOther)-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(7,32) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(7,33) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(7,34) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(8,30) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(8,31) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(8,32) = (-0.4472135954999579*cMOther[11]*mnuOther)-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(8,33) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(8,35) = -0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(9,30) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(9,31) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(9,32) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(9,33) = (-0.4472135954999579*cMOther[11]*mnuOther)-0.4472135954999579*cMOther[10]*mnuOther-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(9,34) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(9,35) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(10,30) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(10,31) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(10,33) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(10,34) = (-0.31943828249997*cMOther[10]*mnuOther)-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(11,30) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(11,32) = -0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(11,33) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(11,35) = (-0.31943828249997*cMOther[11]*mnuOther)-0.5*cMOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(12,6) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(12,7) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(12,8) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(12,9) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(12,10) = 0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(12,11) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(13,6) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(13,7) = 0.4472135954999579*m1rSelf[10]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(13,8) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(13,9) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(13,10) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(14,6) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(14,7) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(14,8) = 0.4472135954999579*m1rSelf[11]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(14,9) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(14,11) = 0.4472135954999579*m1rSelf[8]*mnuSelf; 
  data->AEM_S(15,6) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(15,7) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(15,8) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(15,9) = 0.4472135954999579*m1rSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(15,10) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(15,11) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(16,6) = 0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(16,7) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(16,9) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(16,10) = 0.31943828249997*m1rSelf[10]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(17,6) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(17,8) = 0.4472135954999579*m1rSelf[8]*mnuSelf; 
  data->AEM_S(17,9) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(17,11) = 0.31943828249997*m1rSelf[11]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(12,24) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(12,25) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(12,26) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(12,27) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(12,28) = 0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(12,29) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(13,24) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(13,25) = 0.4472135954999579*m1rOther[10]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(13,26) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(13,27) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(13,28) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(14,24) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(14,25) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(14,26) = 0.4472135954999579*m1rOther[11]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(14,27) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(14,29) = 0.4472135954999579*m1rOther[8]*mnuOther; 
  data->AEM_S(15,24) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(15,25) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(15,26) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(15,27) = 0.4472135954999579*m1rOther[11]*mnuOther+0.4472135954999579*m1rOther[10]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(15,28) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(15,29) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(16,24) = 0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(16,25) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(16,27) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(16,28) = 0.31943828249997*m1rOther[10]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(17,24) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(17,26) = 0.4472135954999579*m1rOther[8]*mnuOther; 
  data->AEM_S(17,27) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(17,29) = 0.31943828249997*m1rOther[11]*mnuOther+0.5*m1rOther[6]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
  mnuM1sum[9] += m1rSelf[9]*mnuSelf+m1rOther[9]*mnuOther; 
  mnuM1sum[10] += m1rSelf[10]*mnuSelf+m1rOther[10]*mnuOther; 
  mnuM1sum[11] += m1rSelf[11]*mnuSelf+m1rOther[11]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(12,12) = m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(12,13) = m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(12,14) = m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(12,15) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(12,16) = m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(12,17) = m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(13,12) = m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(13,14) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(13,15) = m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(13,16) = 0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(14,12) = m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(14,13) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(14,14) = 0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(14,15) = m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(14,17) = 0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(15,12) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(15,13) = m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(15,14) = m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(15,15) = 0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(15,16) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(15,17) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(16,12) = m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(16,13) = 0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(16,15) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(16,16) = 0.6388765649999399*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(17,12) = m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(17,14) = 0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(17,15) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(17,17) = 0.6388765649999399*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(12,30) = m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(12,31) = m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(12,32) = m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(12,33) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(12,34) = m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(12,35) = m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(13,30) = m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(13,31) = 0.8944271909999159*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(13,32) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(13,33) = m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(13,34) = 0.8944271909999159*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(14,30) = m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(14,31) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(14,32) = 0.8944271909999159*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(14,33) = m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(14,35) = 0.8944271909999159*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(15,30) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(15,31) = m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(15,32) = m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(15,33) = 0.8944271909999159*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+0.8944271909999159*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(15,34) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(15,35) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(16,30) = m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(16,31) = 0.8944271909999159*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(16,33) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(16,34) = 0.6388765649999399*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(17,30) = m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(17,32) = 0.8944271909999159*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(17,33) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(17,35) = 0.6388765649999399*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[6]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2rSelf[0]*mnuSelf+m2rOther[0]*mnuOther; 
  mnuM2sum[1] = m2rSelf[1]*mnuSelf+m2rOther[1]*mnuOther; 
  mnuM2sum[2] = m2rSelf[2]*mnuSelf+m2rOther[2]*mnuOther; 
  mnuM2sum[3] = m2rSelf[3]*mnuSelf+m2rOther[3]*mnuOther; 
  mnuM2sum[4] = m2rSelf[4]*mnuSelf+m2rOther[4]*mnuOther; 
  mnuM2sum[5] = m2rSelf[5]*mnuSelf+m2rOther[5]*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<6,6>(0,6).setZero(); 
  data->AEM_S.block<6,6>(6,0).setZero(); 
  data->AEM_S.block<6,6>(0,24).setZero(); 
  data->AEM_S.block<6,6>(6,18).setZero(); 
 
  double m1Relax[12]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<12; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(18,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(18,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(18,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(18,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(19,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(19,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(20,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,1) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,2) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(21,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(22,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(23,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(23,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(18,12) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(18,13) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(18,14) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(18,15) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(18,16) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(18,17) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(19,12) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(19,13) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(19,14) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(19,15) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(19,16) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(20,12) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(20,13) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(20,14) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(20,15) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(20,17) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(21,12) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(21,13) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(21,14) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(21,15) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(21,16) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(21,17) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(22,12) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(22,13) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(22,15) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(22,16) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(23,12) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(23,14) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(23,15) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(23,17) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(18,18) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(18,19) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(18,20) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(18,21) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(18,22) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(18,23) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(19,18) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(19,19) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(19,20) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(19,21) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(19,22) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(20,18) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(20,19) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(20,20) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(20,21) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(20,23) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(21,18) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(21,19) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(21,20) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(21,21) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(21,22) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(21,23) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(22,18) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(22,19) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(22,21) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(22,22) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(23,18) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(23,20) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(23,21) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(23,23) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(18,30) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(18,31) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(18,32) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(18,33) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(18,34) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(18,35) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(19,30) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(19,31) = 0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(19,32) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(19,33) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(19,34) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(20,30) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(20,31) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(20,32) = 0.4472135954999579*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(20,33) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(20,35) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(21,30) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(21,31) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(21,32) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(21,33) = 0.4472135954999579*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(21,34) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(21,35) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(22,30) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(22,31) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(22,33) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(22,34) = 0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(23,30) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(23,32) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(23,33) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(23,35) = 0.31943828249997*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(30,0) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(30,1) = (-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(30,2) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,3) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,4) = (-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf)-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(30,5) = (-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(31,0) = (-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,1) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(31,2) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,3) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,4) = (-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf)-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,5) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(32,0) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,1) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,2) = (-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf)-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(32,3) = (-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(32,4) = (-0.25*m0rSelf[2]*uSelf[4]*mnuSelf)-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(32,5) = (-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,0) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,1) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,2) = (-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,3) = (-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf)-0.2*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.2*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.81*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(33,4) = (-0.2*m0rSelf[3]*uSelf[5]*mnuSelf)-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,5) = (-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(34,0) = (-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf)-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(34,1) = (-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf)-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(34,2) = (-0.25*m0rSelf[2]*uSelf[4]*mnuSelf)-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,3) = (-0.2*m0rSelf[3]*uSelf[5]*mnuSelf)-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(34,4) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(34,5) = (-0.25*m0rSelf[4]*uSelf[5]*mnuSelf)-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(35,0) = (-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(35,1) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,2) = (-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(35,3) = (-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(35,4) = (-0.25*m0rSelf[4]*uSelf[5]*mnuSelf)-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(35,5) = (-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf)-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(30,18) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(30,19) = 0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(30,20) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,21) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,22) = 0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(30,23) = 0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(31,18) = 0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(31,19) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(31,20) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,21) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,22) = 0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(31,23) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(32,18) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(32,19) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(32,20) = 0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(32,21) = 0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(32,22) = 0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(32,23) = 0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,18) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,19) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,20) = 0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(33,21) = 0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.2*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.2*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(33,22) = 0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(33,23) = 0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(34,18) = 0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(34,19) = 0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(34,20) = 0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(34,21) = 0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(34,22) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(34,23) = 0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(35,18) = 0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(35,19) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(35,20) = 0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(35,21) = 0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(35,22) = 0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(35,23) = 0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += (m1rOther[3]-1.0*m1rSelf[3])*betaGreenep1*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (m1rOther[4]-1.0*m1rSelf[4])*betaGreenep1*mnuSelf+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (m1rOther[5]-1.0*m1rSelf[5])*betaGreenep1*mnuSelf+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(24,6) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(24,7) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(24,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,9) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(24,10) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(24,11) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(25,6) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,7) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(25,8) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,9) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(25,10) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(26,6) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,7) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(26,8) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(26,9) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(26,11) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,6) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,7) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,8) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,9) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(27,10) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,11) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(28,6) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(28,7) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(28,9) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(28,10) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(29,6) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(29,8) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,9) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,11) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(24,12) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(24,13) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(24,14) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(24,15) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(24,16) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(24,17) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(25,12) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(25,13) = (-0.4472135954999579*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(25,14) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(25,15) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(25,16) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(26,12) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(26,13) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(26,14) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(26,15) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(26,17) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(27,12) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(27,13) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(27,14) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(27,15) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.4472135954999579*cMSelf[10]*mnuSelf-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(27,16) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(27,17) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(28,12) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(28,13) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(28,15) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(28,16) = (-0.31943828249997*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(29,12) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(29,14) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(29,15) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(29,17) = (-0.31943828249997*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(24,24) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(24,25) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(24,26) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(24,27) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(24,28) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(24,29) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(25,24) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(25,25) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(25,26) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(25,27) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(25,28) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(26,24) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(26,25) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(26,26) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(26,27) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(26,29) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(27,24) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(27,25) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(27,26) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(27,27) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(27,28) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(27,29) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(28,24) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(28,25) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(28,27) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(28,28) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(29,24) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(29,26) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(29,27) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(29,29) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(24,30) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(24,31) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(24,32) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(24,33) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(24,34) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(24,35) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(25,30) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(25,31) = 0.4472135954999579*cMOther[10]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(25,32) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(25,33) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(25,34) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(26,30) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(26,31) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(26,32) = 0.4472135954999579*cMOther[11]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(26,33) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(26,35) = 0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(27,30) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(27,31) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(27,32) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(27,33) = 0.4472135954999579*cMOther[11]*mnuOther+0.4472135954999579*cMOther[10]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(27,34) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(27,35) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(28,30) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(28,31) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(28,33) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(28,34) = 0.31943828249997*cMOther[10]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(29,30) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(29,32) = 0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(29,33) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(29,35) = 0.31943828249997*cMOther[11]*mnuOther+0.5*cMOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(30,6) = (-0.25*m0rSelf[5]*uSelf[11]*mnuSelf)-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(30,7) = (-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf)-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(30,8) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.25*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(30,9) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(30,10) = (-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[4]*uSelf[6]*mnuSelf; 
  data->AEM_S(30,11) = (-0.159719141249985*m0rSelf[5]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[6]*mnuSelf; 
  data->AEM_S(31,6) = (-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf)-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(31,7) = (-0.25*m0rSelf[5]*uSelf[11]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(31,8) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(31,9) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.45*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(31,10) = (-0.3928571428571428*m0rSelf[1]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(31,11) = (-0.25*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[7]*mnuSelf; 
  data->AEM_S(32,6) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.25*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(32,7) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(32,8) = (-0.3928571428571428*m0rSelf[5]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.45*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(32,9) = (-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.45*m0rSelf[2]*uSelf[9]*mnuSelf-0.45*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(32,10) = (-0.25*m0rSelf[2]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[4]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf; 
  data->AEM_S(32,11) = (-0.3928571428571428*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(33,6) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(33,7) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.45*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(33,8) = (-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.45*m0rSelf[2]*uSelf[9]*mnuSelf-0.45*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(33,9) = (-0.3928571428571428*m0rSelf[5]*uSelf[11]*mnuSelf)-0.2*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.2*m0rSelf[5]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf-0.81*m0rSelf[3]*uSelf[9]*mnuSelf-0.45*m0rSelf[2]*uSelf[8]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(33,10) = (-0.2*m0rSelf[3]*uSelf[11]*mnuSelf)-0.3928571428571428*m0rSelf[3]*uSelf[10]*mnuSelf-0.2*m0rSelf[5]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(33,11) = (-0.3928571428571428*m0rSelf[3]*uSelf[11]*mnuSelf)-0.2*m0rSelf[3]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[9]*mnuSelf-0.2*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(34,6) = (-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[4]*uSelf[6]*mnuSelf; 
  data->AEM_S(34,7) = (-0.3928571428571428*m0rSelf[1]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(34,8) = (-0.25*m0rSelf[2]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[4]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf; 
  data->AEM_S(34,9) = (-0.2*m0rSelf[3]*uSelf[11]*mnuSelf)-0.3928571428571428*m0rSelf[3]*uSelf[10]*mnuSelf-0.2*m0rSelf[5]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(34,10) = (-0.25*m0rSelf[5]*uSelf[11]*mnuSelf)-0.5357142857142857*m0rSelf[4]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[10]*mnuSelf+0.31943828249997*m1rSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(34,11) = (-0.25*m0rSelf[4]*uSelf[11]*mnuSelf)-0.25*m0rSelf[5]*uSelf[10]*mnuSelf-0.2*m0rSelf[3]*uSelf[9]*mnuSelf; 
  data->AEM_S(35,6) = (-0.159719141249985*m0rSelf[5]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[6]*mnuSelf; 
  data->AEM_S(35,7) = (-0.25*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[7]*mnuSelf; 
  data->AEM_S(35,8) = (-0.3928571428571428*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(35,9) = (-0.3928571428571428*m0rSelf[3]*uSelf[11]*mnuSelf)-0.2*m0rSelf[3]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[9]*mnuSelf-0.2*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(35,10) = (-0.25*m0rSelf[4]*uSelf[11]*mnuSelf)-0.25*m0rSelf[5]*uSelf[10]*mnuSelf-0.2*m0rSelf[3]*uSelf[9]*mnuSelf; 
  data->AEM_S(35,11) = (-0.5357142857142857*m0rSelf[5]*uSelf[11]*mnuSelf)-0.159719141249985*m0rSelf[0]*uSelf[11]*mnuSelf+0.31943828249997*m1rSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(30,24) = 0.25*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(30,25) = 0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(30,26) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.25*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(30,27) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(30,28) = 0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[4]*uOther[6]*mnuOther; 
  data->AEM_S(30,29) = 0.159719141249985*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[6]*mnuOther; 
  data->AEM_S(31,24) = 0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(31,25) = 0.25*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.4472135954999579*m1rOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(31,26) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(31,27) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.45*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(31,28) = 0.3928571428571428*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(31,29) = 0.25*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[7]*mnuOther; 
  data->AEM_S(32,24) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.25*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(32,25) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(32,26) = 0.3928571428571428*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.45*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(32,27) = 0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.45*m0rOther[2]*uOther[9]*mnuOther+0.45*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(32,28) = 0.25*m0rOther[2]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[4]*uOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther; 
  data->AEM_S(32,29) = 0.3928571428571428*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[0]*uOther[8]*mnuOther-0.4472135954999579*m1rOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(33,24) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(33,25) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.45*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(33,26) = 0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.45*m0rOther[2]*uOther[9]*mnuOther+0.45*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(33,27) = 0.3928571428571428*m0rOther[5]*uOther[11]*mnuOther+0.2*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.2*m0rOther[5]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.4472135954999579*m1rOther[10]*mnuOther+0.81*m0rOther[3]*uOther[9]*mnuOther+0.45*m0rOther[2]*uOther[8]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(33,28) = 0.2*m0rOther[3]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[10]*mnuOther+0.2*m0rOther[5]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(33,29) = 0.3928571428571428*m0rOther[3]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[9]*mnuOther+0.2*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(34,24) = 0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[4]*uOther[6]*mnuOther; 
  data->AEM_S(34,25) = 0.3928571428571428*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(34,26) = 0.25*m0rOther[2]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[4]*uOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther; 
  data->AEM_S(34,27) = 0.2*m0rOther[3]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[10]*mnuOther+0.2*m0rOther[5]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(34,28) = 0.25*m0rOther[5]*uOther[11]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[10]*mnuOther+0.159719141249985*m0rOther[0]*uOther[10]*mnuOther-0.31943828249997*m1rOther[10]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(34,29) = 0.25*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[5]*uOther[10]*mnuOther+0.2*m0rOther[3]*uOther[9]*mnuOther; 
  data->AEM_S(35,24) = 0.159719141249985*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[6]*mnuOther; 
  data->AEM_S(35,25) = 0.25*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[7]*mnuOther; 
  data->AEM_S(35,26) = 0.3928571428571428*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[0]*uOther[8]*mnuOther-0.4472135954999579*m1rOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(35,27) = 0.3928571428571428*m0rOther[3]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[9]*mnuOther+0.2*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(35,28) = 0.25*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[5]*uOther[10]*mnuOther+0.2*m0rOther[3]*uOther[9]*mnuOther; 
  data->AEM_S(35,29) = 0.5357142857142857*m0rOther[5]*uOther[11]*mnuOther+0.159719141249985*m0rOther[0]*uOther[11]*mnuOther-0.31943828249997*m1rOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*m0rOther[5]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[6] += (m1rOther[6]-1.0*m1rSelf[6])*betaGreenep1*mnuSelf+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (m1rOther[7]-1.0*m1rSelf[7])*betaGreenep1*mnuSelf+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
  m1Relax[8] += (m1rOther[8]-1.0*m1rSelf[8])*betaGreenep1*mnuSelf+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
  m1Relax[9] += (m1rOther[9]-1.0*m1rSelf[9])*betaGreenep1*mnuSelf+m1rSelf[9]*mnuSelf-1.0*m1rOther[9]*mnuOther; 
  m1Relax[10] += (m1rOther[10]-1.0*m1rSelf[10])*betaGreenep1*mnuSelf+m1rSelf[10]*mnuSelf-1.0*m1rOther[10]*mnuOther; 
  m1Relax[11] += (m1rOther[11]-1.0*m1rSelf[11])*betaGreenep1*mnuSelf+m1rSelf[11]*mnuSelf-1.0*m1rOther[11]*mnuOther; 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(30,12) = 0.25*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.25*cMSelf[9]*uSelf[9]*mnuSelf+0.25*cMSelf[8]*uSelf[8]*mnuSelf+0.25*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[6]*uSelf[6]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(30,13) = 0.223606797749979*cMSelf[7]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[10]*mnuSelf+0.25*cMSelf[8]*uSelf[9]*mnuSelf+0.25*uSelf[8]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[7]*mnuSelf+0.25*uSelf[6]*cMSelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(30,14) = 0.223606797749979*cMSelf[8]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[8]*cMSelf[11]*mnuSelf+0.25*cMSelf[7]*uSelf[9]*mnuSelf+0.25*uSelf[7]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[8]*mnuSelf+0.25*uSelf[6]*cMSelf[8]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(30,15) = 0.223606797749979*cMSelf[9]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[10]*mnuSelf+0.25*cMSelf[6]*uSelf[9]*mnuSelf+0.25*uSelf[6]*cMSelf[9]*mnuSelf+0.25*cMSelf[7]*uSelf[8]*mnuSelf+0.25*uSelf[7]*cMSelf[8]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(30,16) = 0.159719141249985*cMSelf[10]*uSelf[10]*mnuSelf+0.25*cMSelf[6]*uSelf[10]*mnuSelf+0.25*uSelf[6]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[9]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[7]*mnuSelf+m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(30,17) = 0.159719141249985*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[6]*uSelf[11]*mnuSelf+0.25*uSelf[6]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[9]*mnuSelf+0.223606797749979*cMSelf[8]*uSelf[8]*mnuSelf+m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(31,12) = 0.223606797749979*cMSelf[7]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[10]*mnuSelf+0.25*cMSelf[8]*uSelf[9]*mnuSelf+0.25*uSelf[8]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[7]*mnuSelf+0.25*uSelf[6]*cMSelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(31,13) = 0.25*cMSelf[11]*uSelf[11]*mnuSelf+0.3928571428571428*cMSelf[10]*uSelf[10]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[10]*mnuSelf+0.45*cMSelf[9]*uSelf[9]*mnuSelf+0.25*cMSelf[8]*uSelf[8]*mnuSelf+0.45*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[6]*uSelf[6]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(31,14) = 0.223606797749979*cMSelf[9]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[10]*mnuSelf+0.25*cMSelf[6]*uSelf[9]*mnuSelf+0.25*uSelf[6]*cMSelf[9]*mnuSelf+0.25*cMSelf[7]*uSelf[8]*mnuSelf+0.25*uSelf[7]*cMSelf[8]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(31,15) = 0.223606797749979*cMSelf[8]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[8]*cMSelf[11]*mnuSelf+0.25*cMSelf[7]*uSelf[9]*mnuSelf+0.25*uSelf[7]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[8]*mnuSelf+0.25*uSelf[6]*cMSelf[8]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(31,16) = 0.2*cMSelf[7]*uSelf[10]*mnuSelf+0.2*uSelf[7]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[8]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[8]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[7]*mnuSelf+0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(32,12) = 0.223606797749979*cMSelf[8]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[8]*cMSelf[11]*mnuSelf+0.25*cMSelf[7]*uSelf[9]*mnuSelf+0.25*uSelf[7]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[8]*mnuSelf+0.25*uSelf[6]*cMSelf[8]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(32,13) = 0.223606797749979*cMSelf[9]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[10]*mnuSelf+0.25*cMSelf[6]*uSelf[9]*mnuSelf+0.25*uSelf[6]*cMSelf[9]*mnuSelf+0.25*cMSelf[7]*uSelf[8]*mnuSelf+0.25*uSelf[7]*cMSelf[8]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(32,14) = 0.3928571428571428*cMSelf[11]*uSelf[11]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.45*cMSelf[9]*uSelf[9]*mnuSelf+0.45*cMSelf[8]*uSelf[8]*mnuSelf+0.25*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[6]*uSelf[6]*mnuSelf+0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(32,15) = 0.223606797749979*cMSelf[7]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[10]*mnuSelf+0.25*cMSelf[8]*uSelf[9]*mnuSelf+0.25*uSelf[8]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[7]*mnuSelf+0.25*uSelf[6]*cMSelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(32,17) = 0.2*cMSelf[8]*uSelf[11]*mnuSelf+0.2*uSelf[8]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[8]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(33,12) = 0.223606797749979*cMSelf[9]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[9]*cMSelf[10]*mnuSelf+0.25*cMSelf[6]*uSelf[9]*mnuSelf+0.25*uSelf[6]*cMSelf[9]*mnuSelf+0.25*cMSelf[7]*uSelf[8]*mnuSelf+0.25*uSelf[7]*cMSelf[8]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(33,13) = 0.223606797749979*cMSelf[8]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[8]*cMSelf[11]*mnuSelf+0.25*cMSelf[7]*uSelf[9]*mnuSelf+0.25*uSelf[7]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[8]*mnuSelf+0.25*uSelf[6]*cMSelf[8]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(33,14) = 0.223606797749979*cMSelf[7]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[10]*mnuSelf+0.25*cMSelf[8]*uSelf[9]*mnuSelf+0.25*uSelf[8]*cMSelf[9]*mnuSelf+0.25*cMSelf[6]*uSelf[7]*mnuSelf+0.25*uSelf[6]*cMSelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(33,15) = 0.3928571428571428*cMSelf[11]*uSelf[11]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[11]*mnuSelf+0.3928571428571428*cMSelf[10]*uSelf[10]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[10]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[10]*mnuSelf+0.65*cMSelf[9]*uSelf[9]*mnuSelf+0.45*cMSelf[8]*uSelf[8]*mnuSelf+0.45*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[6]*uSelf[6]*mnuSelf+0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(33,16) = 0.2*cMSelf[9]*uSelf[11]*mnuSelf+0.2*uSelf[9]*cMSelf[11]*mnuSelf+0.2*cMSelf[9]*uSelf[10]*mnuSelf+0.2*uSelf[9]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[8]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(33,17) = 0.2*cMSelf[9]*uSelf[11]*mnuSelf+0.2*uSelf[9]*cMSelf[11]*mnuSelf+0.2*cMSelf[9]*uSelf[10]*mnuSelf+0.2*uSelf[9]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[8]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(34,12) = 0.159719141249985*cMSelf[10]*uSelf[10]*mnuSelf+0.25*cMSelf[6]*uSelf[10]*mnuSelf+0.25*uSelf[6]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[9]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[7]*mnuSelf+m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(34,13) = 0.2*cMSelf[7]*uSelf[10]*mnuSelf+0.2*uSelf[7]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[8]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[8]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[7]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[7]*mnuSelf+0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(34,15) = 0.2*cMSelf[9]*uSelf[11]*mnuSelf+0.2*uSelf[9]*cMSelf[11]*mnuSelf+0.2*cMSelf[9]*uSelf[10]*mnuSelf+0.2*uSelf[9]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[8]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(34,16) = 0.25*cMSelf[11]*uSelf[11]*mnuSelf+0.3520408163265306*cMSelf[10]*uSelf[10]*mnuSelf+0.159719141249985*cMSelf[6]*uSelf[10]*mnuSelf+0.159719141249985*uSelf[6]*cMSelf[10]*mnuSelf+0.3928571428571428*cMSelf[9]*uSelf[9]*mnuSelf+0.25*cMSelf[8]*uSelf[8]*mnuSelf+0.3928571428571428*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[6]*uSelf[6]*mnuSelf+0.6388765649999399*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(35,12) = 0.159719141249985*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[6]*uSelf[11]*mnuSelf+0.25*uSelf[6]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[9]*uSelf[9]*mnuSelf+0.223606797749979*cMSelf[8]*uSelf[8]*mnuSelf+m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(35,14) = 0.2*cMSelf[8]*uSelf[11]*mnuSelf+0.2*uSelf[8]*cMSelf[11]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[8]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(35,15) = 0.2*cMSelf[9]*uSelf[11]*mnuSelf+0.2*uSelf[9]*cMSelf[11]*mnuSelf+0.2*cMSelf[9]*uSelf[10]*mnuSelf+0.2*uSelf[9]*cMSelf[10]*mnuSelf+0.223606797749979*cMSelf[6]*uSelf[9]*mnuSelf+0.223606797749979*uSelf[6]*cMSelf[9]*mnuSelf+0.223606797749979*cMSelf[7]*uSelf[8]*mnuSelf+0.223606797749979*uSelf[7]*cMSelf[8]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(35,17) = 0.3520408163265306*cMSelf[11]*uSelf[11]*mnuSelf+0.159719141249985*cMSelf[6]*uSelf[11]*mnuSelf+0.159719141249985*uSelf[6]*cMSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.3928571428571428*cMSelf[9]*uSelf[9]*mnuSelf+0.3928571428571428*cMSelf[8]*uSelf[8]*mnuSelf+0.25*cMSelf[7]*uSelf[7]*mnuSelf+0.25*cMSelf[6]*uSelf[6]*mnuSelf+0.6388765649999399*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(30,30) = (-0.25*cMOther[11]*uOther[11]*mnuOther)-0.25*cMOther[10]*uOther[10]*mnuOther-0.25*cMOther[9]*uOther[9]*mnuOther-0.25*cMOther[8]*uOther[8]*mnuOther-0.25*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[6]*uOther[6]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(30,31) = (-0.223606797749979*cMOther[7]*uOther[10]*mnuOther)-0.223606797749979*uOther[7]*cMOther[10]*mnuOther-0.25*cMOther[8]*uOther[9]*mnuOther-0.25*uOther[8]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[7]*mnuOther-0.25*uOther[6]*cMOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(30,32) = (-0.223606797749979*cMOther[8]*uOther[11]*mnuOther)-0.223606797749979*uOther[8]*cMOther[11]*mnuOther-0.25*cMOther[7]*uOther[9]*mnuOther-0.25*uOther[7]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[8]*mnuOther-0.25*uOther[6]*cMOther[8]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(30,33) = (-0.223606797749979*cMOther[9]*uOther[11]*mnuOther)-0.223606797749979*uOther[9]*cMOther[11]*mnuOther-0.223606797749979*cMOther[9]*uOther[10]*mnuOther-0.223606797749979*uOther[9]*cMOther[10]*mnuOther-0.25*cMOther[6]*uOther[9]*mnuOther-0.25*uOther[6]*cMOther[9]*mnuOther-0.25*cMOther[7]*uOther[8]*mnuOther-0.25*uOther[7]*cMOther[8]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(30,34) = (-0.159719141249985*cMOther[10]*uOther[10]*mnuOther)-0.25*cMOther[6]*uOther[10]*mnuOther-0.25*uOther[6]*cMOther[10]*mnuOther-0.223606797749979*cMOther[9]*uOther[9]*mnuOther-0.223606797749979*cMOther[7]*uOther[7]*mnuOther-1.0*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(30,35) = (-0.159719141249985*cMOther[11]*uOther[11]*mnuOther)-0.25*cMOther[6]*uOther[11]*mnuOther-0.25*uOther[6]*cMOther[11]*mnuOther-0.223606797749979*cMOther[9]*uOther[9]*mnuOther-0.223606797749979*cMOther[8]*uOther[8]*mnuOther-1.0*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(31,30) = (-0.223606797749979*cMOther[7]*uOther[10]*mnuOther)-0.223606797749979*uOther[7]*cMOther[10]*mnuOther-0.25*cMOther[8]*uOther[9]*mnuOther-0.25*uOther[8]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[7]*mnuOther-0.25*uOther[6]*cMOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(31,31) = (-0.25*cMOther[11]*uOther[11]*mnuOther)-0.3928571428571428*cMOther[10]*uOther[10]*mnuOther-0.223606797749979*cMOther[6]*uOther[10]*mnuOther-0.223606797749979*uOther[6]*cMOther[10]*mnuOther-0.45*cMOther[9]*uOther[9]*mnuOther-0.25*cMOther[8]*uOther[8]*mnuOther-0.45*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[6]*uOther[6]*mnuOther-0.8944271909999159*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(31,32) = (-0.223606797749979*cMOther[9]*uOther[11]*mnuOther)-0.223606797749979*uOther[9]*cMOther[11]*mnuOther-0.223606797749979*cMOther[9]*uOther[10]*mnuOther-0.223606797749979*uOther[9]*cMOther[10]*mnuOther-0.25*cMOther[6]*uOther[9]*mnuOther-0.25*uOther[6]*cMOther[9]*mnuOther-0.25*cMOther[7]*uOther[8]*mnuOther-0.25*uOther[7]*cMOther[8]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(31,33) = (-0.223606797749979*cMOther[8]*uOther[11]*mnuOther)-0.223606797749979*uOther[8]*cMOther[11]*mnuOther-0.25*cMOther[7]*uOther[9]*mnuOther-0.25*uOther[7]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[8]*mnuOther-0.25*uOther[6]*cMOther[8]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(31,34) = (-0.2*cMOther[7]*uOther[10]*mnuOther)-0.2*uOther[7]*cMOther[10]*mnuOther-0.223606797749979*cMOther[8]*uOther[9]*mnuOther-0.223606797749979*uOther[8]*cMOther[9]*mnuOther-0.223606797749979*cMOther[6]*uOther[7]*mnuOther-0.223606797749979*uOther[6]*cMOther[7]*mnuOther-0.8944271909999159*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(32,30) = (-0.223606797749979*cMOther[8]*uOther[11]*mnuOther)-0.223606797749979*uOther[8]*cMOther[11]*mnuOther-0.25*cMOther[7]*uOther[9]*mnuOther-0.25*uOther[7]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[8]*mnuOther-0.25*uOther[6]*cMOther[8]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(32,31) = (-0.223606797749979*cMOther[9]*uOther[11]*mnuOther)-0.223606797749979*uOther[9]*cMOther[11]*mnuOther-0.223606797749979*cMOther[9]*uOther[10]*mnuOther-0.223606797749979*uOther[9]*cMOther[10]*mnuOther-0.25*cMOther[6]*uOther[9]*mnuOther-0.25*uOther[6]*cMOther[9]*mnuOther-0.25*cMOther[7]*uOther[8]*mnuOther-0.25*uOther[7]*cMOther[8]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(32,32) = (-0.3928571428571428*cMOther[11]*uOther[11]*mnuOther)-0.223606797749979*cMOther[6]*uOther[11]*mnuOther-0.223606797749979*uOther[6]*cMOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.45*cMOther[9]*uOther[9]*mnuOther-0.45*cMOther[8]*uOther[8]*mnuOther-0.25*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[6]*uOther[6]*mnuOther-0.8944271909999159*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(32,33) = (-0.223606797749979*cMOther[7]*uOther[10]*mnuOther)-0.223606797749979*uOther[7]*cMOther[10]*mnuOther-0.25*cMOther[8]*uOther[9]*mnuOther-0.25*uOther[8]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[7]*mnuOther-0.25*uOther[6]*cMOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(32,35) = (-0.2*cMOther[8]*uOther[11]*mnuOther)-0.2*uOther[8]*cMOther[11]*mnuOther-0.223606797749979*cMOther[7]*uOther[9]*mnuOther-0.223606797749979*uOther[7]*cMOther[9]*mnuOther-0.223606797749979*cMOther[6]*uOther[8]*mnuOther-0.223606797749979*uOther[6]*cMOther[8]*mnuOther-0.8944271909999159*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(33,30) = (-0.223606797749979*cMOther[9]*uOther[11]*mnuOther)-0.223606797749979*uOther[9]*cMOther[11]*mnuOther-0.223606797749979*cMOther[9]*uOther[10]*mnuOther-0.223606797749979*uOther[9]*cMOther[10]*mnuOther-0.25*cMOther[6]*uOther[9]*mnuOther-0.25*uOther[6]*cMOther[9]*mnuOther-0.25*cMOther[7]*uOther[8]*mnuOther-0.25*uOther[7]*cMOther[8]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(33,31) = (-0.223606797749979*cMOther[8]*uOther[11]*mnuOther)-0.223606797749979*uOther[8]*cMOther[11]*mnuOther-0.25*cMOther[7]*uOther[9]*mnuOther-0.25*uOther[7]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[8]*mnuOther-0.25*uOther[6]*cMOther[8]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(33,32) = (-0.223606797749979*cMOther[7]*uOther[10]*mnuOther)-0.223606797749979*uOther[7]*cMOther[10]*mnuOther-0.25*cMOther[8]*uOther[9]*mnuOther-0.25*uOther[8]*cMOther[9]*mnuOther-0.25*cMOther[6]*uOther[7]*mnuOther-0.25*uOther[6]*cMOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(33,33) = (-0.3928571428571428*cMOther[11]*uOther[11]*mnuOther)-0.223606797749979*cMOther[6]*uOther[11]*mnuOther-0.223606797749979*uOther[6]*cMOther[11]*mnuOther-0.3928571428571428*cMOther[10]*uOther[10]*mnuOther-0.223606797749979*cMOther[6]*uOther[10]*mnuOther-0.223606797749979*uOther[6]*cMOther[10]*mnuOther-0.65*cMOther[9]*uOther[9]*mnuOther-0.45*cMOther[8]*uOther[8]*mnuOther-0.45*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[6]*uOther[6]*mnuOther-0.8944271909999159*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.8944271909999159*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(33,34) = (-0.2*cMOther[9]*uOther[11]*mnuOther)-0.2*uOther[9]*cMOther[11]*mnuOther-0.2*cMOther[9]*uOther[10]*mnuOther-0.2*uOther[9]*cMOther[10]*mnuOther-0.223606797749979*cMOther[6]*uOther[9]*mnuOther-0.223606797749979*uOther[6]*cMOther[9]*mnuOther-0.223606797749979*cMOther[7]*uOther[8]*mnuOther-0.223606797749979*uOther[7]*cMOther[8]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(33,35) = (-0.2*cMOther[9]*uOther[11]*mnuOther)-0.2*uOther[9]*cMOther[11]*mnuOther-0.2*cMOther[9]*uOther[10]*mnuOther-0.2*uOther[9]*cMOther[10]*mnuOther-0.223606797749979*cMOther[6]*uOther[9]*mnuOther-0.223606797749979*uOther[6]*cMOther[9]*mnuOther-0.223606797749979*cMOther[7]*uOther[8]*mnuOther-0.223606797749979*uOther[7]*cMOther[8]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(34,30) = (-0.159719141249985*cMOther[10]*uOther[10]*mnuOther)-0.25*cMOther[6]*uOther[10]*mnuOther-0.25*uOther[6]*cMOther[10]*mnuOther-0.223606797749979*cMOther[9]*uOther[9]*mnuOther-0.223606797749979*cMOther[7]*uOther[7]*mnuOther-1.0*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(34,31) = (-0.2*cMOther[7]*uOther[10]*mnuOther)-0.2*uOther[7]*cMOther[10]*mnuOther-0.223606797749979*cMOther[8]*uOther[9]*mnuOther-0.223606797749979*uOther[8]*cMOther[9]*mnuOther-0.223606797749979*cMOther[6]*uOther[7]*mnuOther-0.223606797749979*uOther[6]*cMOther[7]*mnuOther-0.8944271909999159*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(34,33) = (-0.2*cMOther[9]*uOther[11]*mnuOther)-0.2*uOther[9]*cMOther[11]*mnuOther-0.2*cMOther[9]*uOther[10]*mnuOther-0.2*uOther[9]*cMOther[10]*mnuOther-0.223606797749979*cMOther[6]*uOther[9]*mnuOther-0.223606797749979*uOther[6]*cMOther[9]*mnuOther-0.223606797749979*cMOther[7]*uOther[8]*mnuOther-0.223606797749979*uOther[7]*cMOther[8]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(34,34) = (-0.25*cMOther[11]*uOther[11]*mnuOther)-0.3520408163265306*cMOther[10]*uOther[10]*mnuOther-0.159719141249985*cMOther[6]*uOther[10]*mnuOther-0.159719141249985*uOther[6]*cMOther[10]*mnuOther-0.3928571428571428*cMOther[9]*uOther[9]*mnuOther-0.25*cMOther[8]*uOther[8]*mnuOther-0.3928571428571428*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[6]*uOther[6]*mnuOther-0.6388765649999399*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(35,30) = (-0.159719141249985*cMOther[11]*uOther[11]*mnuOther)-0.25*cMOther[6]*uOther[11]*mnuOther-0.25*uOther[6]*cMOther[11]*mnuOther-0.223606797749979*cMOther[9]*uOther[9]*mnuOther-0.223606797749979*cMOther[8]*uOther[8]*mnuOther-1.0*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(35,32) = (-0.2*cMOther[8]*uOther[11]*mnuOther)-0.2*uOther[8]*cMOther[11]*mnuOther-0.223606797749979*cMOther[7]*uOther[9]*mnuOther-0.223606797749979*uOther[7]*cMOther[9]*mnuOther-0.223606797749979*cMOther[6]*uOther[8]*mnuOther-0.223606797749979*uOther[6]*cMOther[8]*mnuOther-0.8944271909999159*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(35,33) = (-0.2*cMOther[9]*uOther[11]*mnuOther)-0.2*uOther[9]*cMOther[11]*mnuOther-0.2*cMOther[9]*uOther[10]*mnuOther-0.2*uOther[9]*cMOther[10]*mnuOther-0.223606797749979*cMOther[6]*uOther[9]*mnuOther-0.223606797749979*uOther[6]*cMOther[9]*mnuOther-0.223606797749979*cMOther[7]*uOther[8]*mnuOther-0.223606797749979*uOther[7]*cMOther[8]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(35,35) = (-0.3520408163265306*cMOther[11]*uOther[11]*mnuOther)-0.159719141249985*cMOther[6]*uOther[11]*mnuOther-0.159719141249985*uOther[6]*cMOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.3928571428571428*cMOther[9]*uOther[9]*mnuOther-0.3928571428571428*cMOther[8]*uOther[8]*mnuOther-0.25*cMOther[7]*uOther[7]*mnuOther-0.25*cMOther[6]*uOther[6]*mnuOther-0.6388765649999399*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[6]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  // Set other entries to 0. // 
  data->AEM_S.block<6,6>(18,6).setZero(); 
  data->AEM_S.block<6,6>(24,0).setZero(); 
  data->AEM_S.block<6,6>(18,24).setZero(); 
  data->AEM_S.block<6,6>(24,18).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM1sum[10],mnuM1sum[11],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m1Relax[9],m1Relax[10],m1Relax[11],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,12,1) = data->u_S.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,6,1) = data->u_S.segment<6>(12); 
 
  Eigen::Map<VectorXd>(uCrossOther,12,1) = data->u_S.segment<12>(18); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,6,1) = data->u_S.segment<6>(30); 
 
} 
 
void VmCrossPrimMoments2x2vMax_P3(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[20]; 
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
    m1rSelf[10] = m1Self[10]; 
    m1rSelf[11] = 0.0; 
    m1rSelf[12] = 0.0; 
    m1rSelf[13] = 0.0; 
    m1rSelf[14] = 0.0; 
    m1rSelf[15] = 0.0; 
    m1rSelf[16] = 0.0; 
    m1rSelf[17] = 0.0; 
    m1rSelf[18] = 0.0; 
    m1rSelf[19] = 0.0; 
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
    m1rSelf[10] = m1Self[10]; 
    m1rSelf[11] = m1Self[11]; 
    m1rSelf[12] = m1Self[12]; 
    m1rSelf[13] = m1Self[13]; 
    m1rSelf[14] = m1Self[14]; 
    m1rSelf[15] = m1Self[15]; 
    m1rSelf[16] = m1Self[16]; 
    m1rSelf[17] = m1Self[17]; 
    m1rSelf[18] = m1Self[18]; 
    m1rSelf[19] = m1Self[19]; 
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
  double m1rOther[20]; 
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
    m1rOther[10] = m1Other[10]; 
    m1rOther[11] = 0.0; 
    m1rOther[12] = 0.0; 
    m1rOther[13] = 0.0; 
    m1rOther[14] = 0.0; 
    m1rOther[15] = 0.0; 
    m1rOther[16] = 0.0; 
    m1rOther[17] = 0.0; 
    m1rOther[18] = 0.0; 
    m1rOther[19] = 0.0; 
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
    m1rOther[10] = m1Other[10]; 
    m1rOther[11] = m1Other[11]; 
    m1rOther[12] = m1Other[12]; 
    m1rOther[13] = m1Other[13]; 
    m1rOther[14] = m1Other[14]; 
    m1rOther[15] = m1Other[15]; 
    m1rOther[16] = m1Other[16]; 
    m1rOther[17] = m1Other[17]; 
    m1rOther[18] = m1Other[18]; 
    m1rOther[19] = m1Other[19]; 
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
  data->AEM_S = Eigen::MatrixXd::Zero(60,60); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
  double mnuM1sum[20]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<20; vd++) 
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
  data->AEM_S(0,20) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,21) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,22) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,23) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,24) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,25) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(0,26) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(0,27) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(0,28) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(0,29) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(1,20) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,21) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,22) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,23) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,24) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,25) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(1,26) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,27) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(1,28) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(2,20) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,21) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,22) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,23) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,24) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(2,25) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,26) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(2,27) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,29) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(3,20) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,21) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,22) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,23) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,24) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,25) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,26) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,27) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,28) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(3,29) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(4,20) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,21) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,22) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(4,23) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,24) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(4,26) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(4,27) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(4,28) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(5,20) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,21) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,22) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(5,23) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,25) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,26) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(5,27) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(5,29) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(6,20) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,21) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(6,22) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(6,23) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,24) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(6,25) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,26) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(6,27) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(6,28) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,20) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,21) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(7,22) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,23) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(7,24) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,25) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(7,26) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,27) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(7,29) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(8,20) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,21) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(8,23) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(8,24) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(8,26) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(8,28) = (-0.2981423969999719*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(9,20) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,22) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(9,23) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(9,25) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(9,27) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(9,29) = (-0.2981423969999719*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,30) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,31) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,32) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,33) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(0,34) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(0,35) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(0,36) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(0,37) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(0,38) = 0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(0,39) = 0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(1,30) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,31) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,32) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,33) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,34) = 0.4391550328268398*m0rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(1,35) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(1,36) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(1,37) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(1,38) = 0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(2,30) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,31) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,32) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,33) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(2,34) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(2,35) = 0.4391550328268398*m0rOther[9]*mnuOther+0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(2,36) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(2,37) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(2,39) = 0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(3,30) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,31) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,32) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,33) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,34) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,35) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,36) = 0.4391550328268399*m0rOther[8]*mnuOther+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(3,37) = 0.4391550328268399*m0rOther[9]*mnuOther+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(3,38) = 0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(3,39) = 0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(4,30) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(4,31) = 0.4391550328268398*m0rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(4,32) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(4,33) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,34) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(4,36) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(4,37) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(4,38) = 0.2981423969999719*m0rOther[8]*mnuOther+0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(5,30) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(5,31) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(5,32) = 0.4391550328268398*m0rOther[9]*mnuOther+0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(5,33) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(5,35) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,36) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(5,37) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(5,39) = 0.2981423969999719*m0rOther[9]*mnuOther+0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(6,30) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(6,31) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(6,32) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(6,33) = 0.4391550328268399*m0rOther[8]*mnuOther+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(6,34) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(6,35) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(6,36) = 0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,37) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(6,38) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(7,30) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(7,31) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(7,32) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(7,33) = 0.4391550328268399*m0rOther[9]*mnuOther+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(7,34) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(7,35) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(7,36) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(7,37) = 0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(7,39) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(8,30) = 0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(8,31) = 0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(8,33) = 0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(8,34) = 0.2981423969999719*m0rOther[8]*mnuOther+0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(8,36) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(8,38) = 0.2981423969999719*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,30) = 0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(9,32) = 0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(9,33) = 0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(9,35) = 0.2981423969999719*m0rOther[9]*mnuOther+0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(9,37) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(9,39) = 0.2981423969999719*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,50) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,51) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,52) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,53) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(0,54) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(0,55) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(0,56) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(0,57) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(0,58) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(0,59) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(1,50) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,51) = (-0.4472135954999579*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,52) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,53) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,54) = (-0.4391550328268398*cMOther[8]*mnuOther)-0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(1,55) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(1,56) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(1,57) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(1,58) = -0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(2,50) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,51) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,52) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,53) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(2,54) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(2,55) = (-0.4391550328268398*cMOther[9]*mnuOther)-0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(2,56) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(2,57) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(2,59) = -0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(3,50) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,51) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,52) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,53) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(3,54) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,55) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,56) = (-0.4391550328268399*cMOther[8]*mnuOther)-0.4*cMOther[7]*mnuOther-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(3,57) = (-0.4391550328268399*cMOther[9]*mnuOther)-0.4*cMOther[6]*mnuOther-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(3,58) = -0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(3,59) = -0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(4,50) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,51) = (-0.4391550328268398*cMOther[8]*mnuOther)-0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(4,52) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(4,53) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,54) = (-0.31943828249997*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(4,56) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(4,57) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(4,58) = (-0.2981423969999719*cMOther[8]*mnuOther)-0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(5,50) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,51) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(5,52) = (-0.4391550328268398*cMOther[9]*mnuOther)-0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(5,53) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(5,55) = (-0.31943828249997*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(5,56) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(5,57) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(5,59) = (-0.2981423969999719*cMOther[9]*mnuOther)-0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(6,50) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,51) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(6,52) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(6,53) = (-0.4391550328268399*cMOther[8]*mnuOther)-0.4*cMOther[7]*mnuOther-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(6,54) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(6,55) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(6,56) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.31943828249997*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(6,57) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(6,58) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(7,50) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,51) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(7,52) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(7,53) = (-0.4391550328268399*cMOther[9]*mnuOther)-0.4*cMOther[6]*mnuOther-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(7,54) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(7,55) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(7,56) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(7,57) = (-0.31943828249997*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(7,59) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(8,50) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(8,51) = -0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(8,53) = -0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(8,54) = (-0.2981423969999719*cMOther[8]*mnuOther)-0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(8,56) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(8,58) = (-0.2981423969999719*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(9,50) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(9,52) = -0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(9,53) = -0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(9,55) = (-0.2981423969999719*cMOther[9]*mnuOther)-0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(9,57) = -0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(9,59) = (-0.2981423969999719*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(20,0) = 0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(20,1) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(20,2) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(20,3) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(20,4) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(20,5) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(20,6) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(20,7) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(20,8) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(20,9) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(21,0) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(21,1) = 0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(21,2) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(21,3) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(21,4) = 0.4391550328268398*m1rSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(21,5) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(21,6) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(21,7) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(21,8) = 0.4391550328268398*m1rSelf[4]*mnuSelf; 
  data->AEM_S(22,0) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(22,1) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(22,2) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(22,3) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(22,4) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(22,5) = 0.4391550328268398*m1rSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(22,6) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(22,7) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(22,9) = 0.4391550328268398*m1rSelf[5]*mnuSelf; 
  data->AEM_S(23,0) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(23,1) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(23,2) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(23,3) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(23,4) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(23,5) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(23,6) = 0.4391550328268399*m1rSelf[8]*mnuSelf+0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(23,7) = 0.4391550328268399*m1rSelf[9]*mnuSelf+0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(23,8) = 0.4391550328268399*m1rSelf[6]*mnuSelf; 
  data->AEM_S(23,9) = 0.4391550328268399*m1rSelf[7]*mnuSelf; 
  data->AEM_S(24,0) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(24,1) = 0.4391550328268398*m1rSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(24,2) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(24,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(24,4) = 0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(24,6) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(24,7) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(24,8) = 0.2981423969999719*m1rSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf; 
  data->AEM_S(25,0) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(25,1) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(25,2) = 0.4391550328268398*m1rSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(25,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(25,5) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(25,6) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(25,7) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(25,9) = 0.2981423969999719*m1rSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf; 
  data->AEM_S(26,0) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(26,1) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(26,2) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(26,3) = 0.4391550328268399*m1rSelf[8]*mnuSelf+0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(26,4) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(26,5) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(26,6) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(26,7) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(26,8) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(27,0) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(27,1) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(27,2) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(27,3) = 0.4391550328268399*m1rSelf[9]*mnuSelf+0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(27,4) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(27,5) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(27,6) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(27,7) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(27,9) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(28,0) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(28,1) = 0.4391550328268398*m1rSelf[4]*mnuSelf; 
  data->AEM_S(28,3) = 0.4391550328268399*m1rSelf[6]*mnuSelf; 
  data->AEM_S(28,4) = 0.2981423969999719*m1rSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf; 
  data->AEM_S(28,6) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(28,8) = 0.2981423969999719*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(29,0) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(29,2) = 0.4391550328268398*m1rSelf[5]*mnuSelf; 
  data->AEM_S(29,3) = 0.4391550328268399*m1rSelf[7]*mnuSelf; 
  data->AEM_S(29,5) = 0.2981423969999719*m1rSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf; 
  data->AEM_S(29,7) = 0.4391550328268399*m1rSelf[3]*mnuSelf; 
  data->AEM_S(29,9) = 0.2981423969999719*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(20,30) = 0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(20,31) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(20,32) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(20,33) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(20,34) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(20,35) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(20,36) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(20,37) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(20,38) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(20,39) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(21,30) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(21,31) = 0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(21,32) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(21,33) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(21,34) = 0.4391550328268398*m1rOther[8]*mnuOther+0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(21,35) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(21,36) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(21,37) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(21,38) = 0.4391550328268398*m1rOther[4]*mnuOther; 
  data->AEM_S(22,30) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(22,31) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(22,32) = 0.4472135954999579*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(22,33) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(22,34) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(22,35) = 0.4391550328268398*m1rOther[9]*mnuOther+0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(22,36) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(22,37) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(22,39) = 0.4391550328268398*m1rOther[5]*mnuOther; 
  data->AEM_S(23,30) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(23,31) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(23,32) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(23,33) = 0.4472135954999579*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(23,34) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(23,35) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(23,36) = 0.4391550328268399*m1rOther[8]*mnuOther+0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(23,37) = 0.4391550328268399*m1rOther[9]*mnuOther+0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(23,38) = 0.4391550328268399*m1rOther[6]*mnuOther; 
  data->AEM_S(23,39) = 0.4391550328268399*m1rOther[7]*mnuOther; 
  data->AEM_S(24,30) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(24,31) = 0.4391550328268398*m1rOther[8]*mnuOther+0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(24,32) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(24,33) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(24,34) = 0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(24,36) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(24,37) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(24,38) = 0.2981423969999719*m1rOther[8]*mnuOther+0.4391550328268398*m1rOther[1]*mnuOther; 
  data->AEM_S(25,30) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(25,31) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(25,32) = 0.4391550328268398*m1rOther[9]*mnuOther+0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(25,33) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(25,35) = 0.31943828249997*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(25,36) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(25,37) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(25,39) = 0.2981423969999719*m1rOther[9]*mnuOther+0.4391550328268398*m1rOther[2]*mnuOther; 
  data->AEM_S(26,30) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(26,31) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(26,32) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(26,33) = 0.4391550328268399*m1rOther[8]*mnuOther+0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(26,34) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(26,35) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(26,36) = 0.4472135954999579*m1rOther[5]*mnuOther+0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(26,37) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(26,38) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(27,30) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(27,31) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(27,32) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(27,33) = 0.4391550328268399*m1rOther[9]*mnuOther+0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(27,34) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(27,35) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(27,36) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(27,37) = 0.31943828249997*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(27,39) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(28,30) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(28,31) = 0.4391550328268398*m1rOther[4]*mnuOther; 
  data->AEM_S(28,33) = 0.4391550328268399*m1rOther[6]*mnuOther; 
  data->AEM_S(28,34) = 0.2981423969999719*m1rOther[8]*mnuOther+0.4391550328268398*m1rOther[1]*mnuOther; 
  data->AEM_S(28,36) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(28,38) = 0.2981423969999719*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(29,30) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(29,32) = 0.4391550328268398*m1rOther[5]*mnuOther; 
  data->AEM_S(29,33) = 0.4391550328268399*m1rOther[7]*mnuOther; 
  data->AEM_S(29,35) = 0.2981423969999719*m1rOther[9]*mnuOther+0.4391550328268398*m1rOther[2]*mnuOther; 
  data->AEM_S(29,37) = 0.4391550328268399*m1rOther[3]*mnuOther; 
  data->AEM_S(29,39) = 0.2981423969999719*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
 
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
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(10,10) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(10,11) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(10,12) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,13) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(10,14) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(10,15) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(10,16) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(10,17) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(10,18) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(10,19) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(11,10) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,11) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(11,12) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,13) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,14) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,15) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(11,16) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,17) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(11,18) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(12,10) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,11) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(12,12) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,14) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(12,15) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,16) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(12,17) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(12,19) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(13,10) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(13,11) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,12) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(13,15) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(13,16) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,17) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,18) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(13,19) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(14,10) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(14,11) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,12) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(14,13) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(14,14) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(14,16) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,17) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(14,18) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,10) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(15,11) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(15,12) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,13) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,15) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(15,16) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(15,17) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,19) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,10) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(16,11) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,12) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(16,13) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,14) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,15) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(16,16) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(16,17) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,18) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,10) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(17,11) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(17,12) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,13) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,14) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(17,15) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,16) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,17) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(17,19) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(18,10) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(18,11) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(18,13) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(18,14) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(18,16) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(18,18) = 0.2981423969999719*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(19,10) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(19,12) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(19,13) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(19,15) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,17) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,19) = 0.2981423969999719*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(10,20) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,21) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(10,22) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(10,23) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(10,24) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(10,25) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(10,26) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(10,27) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(10,28) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(10,29) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(11,20) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,21) = (-0.4472135954999579*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(11,22) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(11,23) = (-0.447213595499958*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(11,24) = (-0.4391550328268398*cMSelf[18]*mnuSelf)-0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,25) = -0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(11,26) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(11,27) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(11,28) = -0.4391550328268398*cMSelf[14]*mnuSelf; 
  data->AEM_S(12,20) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(12,21) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(12,22) = (-0.4472135954999579*cMSelf[15]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(12,23) = (-0.447213595499958*cMSelf[17]*mnuSelf)-0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(12,24) = -0.5000000000000001*cMSelf[16]*mnuSelf; 
  data->AEM_S(12,25) = (-0.4391550328268398*cMSelf[19]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf; 
  data->AEM_S(12,26) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(12,27) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(12,29) = -0.4391550328268398*cMSelf[15]*mnuSelf; 
  data->AEM_S(13,20) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(13,21) = (-0.447213595499958*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(13,22) = (-0.447213595499958*cMSelf[17]*mnuSelf)-0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(13,23) = (-0.4472135954999579*cMSelf[15]*mnuSelf)-0.4472135954999579*cMSelf[14]*mnuSelf-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(13,24) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(13,25) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(13,26) = (-0.4391550328268399*cMSelf[18]*mnuSelf)-0.4*cMSelf[17]*mnuSelf-0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(13,27) = (-0.4391550328268399*cMSelf[19]*mnuSelf)-0.4*cMSelf[16]*mnuSelf-0.447213595499958*cMSelf[12]*mnuSelf; 
  data->AEM_S(13,28) = -0.4391550328268399*cMSelf[16]*mnuSelf; 
  data->AEM_S(13,29) = -0.4391550328268399*cMSelf[17]*mnuSelf; 
  data->AEM_S(14,20) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(14,21) = (-0.4391550328268398*cMSelf[18]*mnuSelf)-0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(14,22) = -0.5000000000000001*cMSelf[16]*mnuSelf; 
  data->AEM_S(14,23) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(14,24) = (-0.31943828249997*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(14,26) = (-0.31943828249997*cMSelf[16]*mnuSelf)-0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(14,27) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(14,28) = (-0.2981423969999719*cMSelf[18]*mnuSelf)-0.4391550328268398*cMSelf[11]*mnuSelf; 
  data->AEM_S(15,20) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(15,21) = -0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(15,22) = (-0.4391550328268398*cMSelf[19]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf; 
  data->AEM_S(15,23) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(15,25) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(15,26) = -0.4472135954999579*cMSelf[16]*mnuSelf; 
  data->AEM_S(15,27) = (-0.31943828249997*cMSelf[17]*mnuSelf)-0.5000000000000001*cMSelf[11]*mnuSelf; 
  data->AEM_S(15,29) = (-0.2981423969999719*cMSelf[19]*mnuSelf)-0.4391550328268398*cMSelf[12]*mnuSelf; 
  data->AEM_S(16,20) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(16,21) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(16,22) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(16,23) = (-0.4391550328268399*cMSelf[18]*mnuSelf)-0.4*cMSelf[17]*mnuSelf-0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(16,24) = (-0.31943828249997*cMSelf[16]*mnuSelf)-0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(16,25) = -0.4472135954999579*cMSelf[16]*mnuSelf; 
  data->AEM_S(16,26) = (-0.4472135954999579*cMSelf[15]*mnuSelf)-0.31943828249997*cMSelf[14]*mnuSelf-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(16,27) = -0.4*cMSelf[13]*mnuSelf; 
  data->AEM_S(16,28) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(17,20) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(17,21) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(17,22) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(17,23) = (-0.4391550328268399*cMSelf[19]*mnuSelf)-0.4*cMSelf[16]*mnuSelf-0.447213595499958*cMSelf[12]*mnuSelf; 
  data->AEM_S(17,24) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(17,25) = (-0.31943828249997*cMSelf[17]*mnuSelf)-0.5000000000000001*cMSelf[11]*mnuSelf; 
  data->AEM_S(17,26) = -0.4*cMSelf[13]*mnuSelf; 
  data->AEM_S(17,27) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.4472135954999579*cMSelf[14]*mnuSelf-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(17,29) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(18,20) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(18,21) = -0.4391550328268398*cMSelf[14]*mnuSelf; 
  data->AEM_S(18,23) = -0.4391550328268399*cMSelf[16]*mnuSelf; 
  data->AEM_S(18,24) = (-0.2981423969999719*cMSelf[18]*mnuSelf)-0.4391550328268398*cMSelf[11]*mnuSelf; 
  data->AEM_S(18,26) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(18,28) = (-0.2981423969999719*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(19,20) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(19,22) = -0.4391550328268398*cMSelf[15]*mnuSelf; 
  data->AEM_S(19,23) = -0.4391550328268399*cMSelf[17]*mnuSelf; 
  data->AEM_S(19,25) = (-0.2981423969999719*cMSelf[19]*mnuSelf)-0.4391550328268398*cMSelf[12]*mnuSelf; 
  data->AEM_S(19,27) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(19,29) = (-0.2981423969999719*cMSelf[15]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(10,40) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(10,41) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(10,42) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,43) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(10,44) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(10,45) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(10,46) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(10,47) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(10,48) = 0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(10,49) = 0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(11,40) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,41) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(11,42) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(11,43) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(11,44) = 0.4391550328268398*m0rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(11,45) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(11,46) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(11,47) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(11,48) = 0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(12,40) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(12,41) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(12,42) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(12,43) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(12,44) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(12,45) = 0.4391550328268398*m0rOther[9]*mnuOther+0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(12,46) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(12,47) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(12,49) = 0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(13,40) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(13,41) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(13,42) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(13,43) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(13,44) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(13,45) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(13,46) = 0.4391550328268399*m0rOther[8]*mnuOther+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(13,47) = 0.4391550328268399*m0rOther[9]*mnuOther+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(13,48) = 0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(13,49) = 0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(14,40) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(14,41) = 0.4391550328268398*m0rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(14,42) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(14,43) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(14,44) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(14,46) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(14,47) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(14,48) = 0.2981423969999719*m0rOther[8]*mnuOther+0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(15,40) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(15,41) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(15,42) = 0.4391550328268398*m0rOther[9]*mnuOther+0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(15,43) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(15,45) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(15,46) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(15,47) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(15,49) = 0.2981423969999719*m0rOther[9]*mnuOther+0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(16,40) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(16,41) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(16,42) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(16,43) = 0.4391550328268399*m0rOther[8]*mnuOther+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(16,44) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(16,45) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(16,46) = 0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(16,47) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(16,48) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(17,40) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(17,41) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(17,42) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(17,43) = 0.4391550328268399*m0rOther[9]*mnuOther+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(17,44) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(17,45) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(17,46) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(17,47) = 0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(17,49) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(18,40) = 0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(18,41) = 0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(18,43) = 0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(18,44) = 0.2981423969999719*m0rOther[8]*mnuOther+0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(18,46) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(18,48) = 0.2981423969999719*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(19,40) = 0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(19,42) = 0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(19,43) = 0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(19,45) = 0.2981423969999719*m0rOther[9]*mnuOther+0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(19,47) = 0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(19,49) = 0.2981423969999719*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(10,50) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(10,51) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(10,52) = -0.5*cMOther[12]*mnuOther; 
  data->AEM_S(10,53) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(10,54) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(10,55) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(10,56) = -0.5*cMOther[16]*mnuOther; 
  data->AEM_S(10,57) = -0.5*cMOther[17]*mnuOther; 
  data->AEM_S(10,58) = -0.5*cMOther[18]*mnuOther; 
  data->AEM_S(10,59) = -0.5*cMOther[19]*mnuOther; 
  data->AEM_S(11,50) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(11,51) = (-0.4472135954999579*cMOther[14]*mnuOther)-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(11,52) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(11,53) = (-0.447213595499958*cMOther[16]*mnuOther)-0.5*cMOther[12]*mnuOther; 
  data->AEM_S(11,54) = (-0.4391550328268398*cMOther[18]*mnuOther)-0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(11,55) = -0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(11,56) = -0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(11,57) = -0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(11,58) = -0.4391550328268398*cMOther[14]*mnuOther; 
  data->AEM_S(12,50) = -0.5*cMOther[12]*mnuOther; 
  data->AEM_S(12,51) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(12,52) = (-0.4472135954999579*cMOther[15]*mnuOther)-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(12,53) = (-0.447213595499958*cMOther[17]*mnuOther)-0.5*cMOther[11]*mnuOther; 
  data->AEM_S(12,54) = -0.5000000000000001*cMOther[16]*mnuOther; 
  data->AEM_S(12,55) = (-0.4391550328268398*cMOther[19]*mnuOther)-0.4472135954999579*cMOther[12]*mnuOther; 
  data->AEM_S(12,56) = -0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(12,57) = -0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(12,59) = -0.4391550328268398*cMOther[15]*mnuOther; 
  data->AEM_S(13,50) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(13,51) = (-0.447213595499958*cMOther[16]*mnuOther)-0.5*cMOther[12]*mnuOther; 
  data->AEM_S(13,52) = (-0.447213595499958*cMOther[17]*mnuOther)-0.5*cMOther[11]*mnuOther; 
  data->AEM_S(13,53) = (-0.4472135954999579*cMOther[15]*mnuOther)-0.4472135954999579*cMOther[14]*mnuOther-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(13,54) = -0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(13,55) = -0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(13,56) = (-0.4391550328268399*cMOther[18]*mnuOther)-0.4*cMOther[17]*mnuOther-0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(13,57) = (-0.4391550328268399*cMOther[19]*mnuOther)-0.4*cMOther[16]*mnuOther-0.447213595499958*cMOther[12]*mnuOther; 
  data->AEM_S(13,58) = -0.4391550328268399*cMOther[16]*mnuOther; 
  data->AEM_S(13,59) = -0.4391550328268399*cMOther[17]*mnuOther; 
  data->AEM_S(14,50) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(14,51) = (-0.4391550328268398*cMOther[18]*mnuOther)-0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(14,52) = -0.5000000000000001*cMOther[16]*mnuOther; 
  data->AEM_S(14,53) = -0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(14,54) = (-0.31943828249997*cMOther[14]*mnuOther)-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(14,56) = (-0.31943828249997*cMOther[16]*mnuOther)-0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(14,57) = -0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(14,58) = (-0.2981423969999719*cMOther[18]*mnuOther)-0.4391550328268398*cMOther[11]*mnuOther; 
  data->AEM_S(15,50) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(15,51) = -0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(15,52) = (-0.4391550328268398*cMOther[19]*mnuOther)-0.4472135954999579*cMOther[12]*mnuOther; 
  data->AEM_S(15,53) = -0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(15,55) = (-0.31943828249997*cMOther[15]*mnuOther)-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(15,56) = -0.4472135954999579*cMOther[16]*mnuOther; 
  data->AEM_S(15,57) = (-0.31943828249997*cMOther[17]*mnuOther)-0.5000000000000001*cMOther[11]*mnuOther; 
  data->AEM_S(15,59) = (-0.2981423969999719*cMOther[19]*mnuOther)-0.4391550328268398*cMOther[12]*mnuOther; 
  data->AEM_S(16,50) = -0.5*cMOther[16]*mnuOther; 
  data->AEM_S(16,51) = -0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(16,52) = -0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(16,53) = (-0.4391550328268399*cMOther[18]*mnuOther)-0.4*cMOther[17]*mnuOther-0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(16,54) = (-0.31943828249997*cMOther[16]*mnuOther)-0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(16,55) = -0.4472135954999579*cMOther[16]*mnuOther; 
  data->AEM_S(16,56) = (-0.4472135954999579*cMOther[15]*mnuOther)-0.31943828249997*cMOther[14]*mnuOther-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(16,57) = -0.4*cMOther[13]*mnuOther; 
  data->AEM_S(16,58) = -0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(17,50) = -0.5*cMOther[17]*mnuOther; 
  data->AEM_S(17,51) = -0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(17,52) = -0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(17,53) = (-0.4391550328268399*cMOther[19]*mnuOther)-0.4*cMOther[16]*mnuOther-0.447213595499958*cMOther[12]*mnuOther; 
  data->AEM_S(17,54) = -0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(17,55) = (-0.31943828249997*cMOther[17]*mnuOther)-0.5000000000000001*cMOther[11]*mnuOther; 
  data->AEM_S(17,56) = -0.4*cMOther[13]*mnuOther; 
  data->AEM_S(17,57) = (-0.31943828249997*cMOther[15]*mnuOther)-0.4472135954999579*cMOther[14]*mnuOther-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(17,59) = -0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(18,50) = -0.5*cMOther[18]*mnuOther; 
  data->AEM_S(18,51) = -0.4391550328268398*cMOther[14]*mnuOther; 
  data->AEM_S(18,53) = -0.4391550328268399*cMOther[16]*mnuOther; 
  data->AEM_S(18,54) = (-0.2981423969999719*cMOther[18]*mnuOther)-0.4391550328268398*cMOther[11]*mnuOther; 
  data->AEM_S(18,56) = -0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(18,58) = (-0.2981423969999719*cMOther[14]*mnuOther)-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(19,50) = -0.5*cMOther[19]*mnuOther; 
  data->AEM_S(19,52) = -0.4391550328268398*cMOther[15]*mnuOther; 
  data->AEM_S(19,53) = -0.4391550328268399*cMOther[17]*mnuOther; 
  data->AEM_S(19,55) = (-0.2981423969999719*cMOther[19]*mnuOther)-0.4391550328268398*cMOther[12]*mnuOther; 
  data->AEM_S(19,57) = -0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(19,59) = (-0.2981423969999719*cMOther[15]*mnuOther)-0.5*cMOther[10]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(20,10) = 0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(20,11) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(20,12) = 0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(20,13) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(20,14) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(20,15) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(20,16) = 0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(20,17) = 0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(20,18) = 0.5*m1rSelf[18]*mnuSelf; 
  data->AEM_S(20,19) = 0.5*m1rSelf[19]*mnuSelf; 
  data->AEM_S(21,10) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(21,11) = 0.4472135954999579*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(21,12) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(21,13) = 0.447213595499958*m1rSelf[16]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(21,14) = 0.4391550328268398*m1rSelf[18]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf; 
  data->AEM_S(21,15) = 0.5000000000000001*m1rSelf[17]*mnuSelf; 
  data->AEM_S(21,16) = 0.447213595499958*m1rSelf[13]*mnuSelf; 
  data->AEM_S(21,17) = 0.5000000000000001*m1rSelf[15]*mnuSelf; 
  data->AEM_S(21,18) = 0.4391550328268398*m1rSelf[14]*mnuSelf; 
  data->AEM_S(22,10) = 0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(22,11) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(22,12) = 0.4472135954999579*m1rSelf[15]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(22,13) = 0.447213595499958*m1rSelf[17]*mnuSelf+0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(22,14) = 0.5000000000000001*m1rSelf[16]*mnuSelf; 
  data->AEM_S(22,15) = 0.4391550328268398*m1rSelf[19]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf; 
  data->AEM_S(22,16) = 0.5000000000000001*m1rSelf[14]*mnuSelf; 
  data->AEM_S(22,17) = 0.447213595499958*m1rSelf[13]*mnuSelf; 
  data->AEM_S(22,19) = 0.4391550328268398*m1rSelf[15]*mnuSelf; 
  data->AEM_S(23,10) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(23,11) = 0.447213595499958*m1rSelf[16]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(23,12) = 0.447213595499958*m1rSelf[17]*mnuSelf+0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(23,13) = 0.4472135954999579*m1rSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(23,14) = 0.4472135954999579*m1rSelf[13]*mnuSelf; 
  data->AEM_S(23,15) = 0.4472135954999579*m1rSelf[13]*mnuSelf; 
  data->AEM_S(23,16) = 0.4391550328268399*m1rSelf[18]*mnuSelf+0.4*m1rSelf[17]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf; 
  data->AEM_S(23,17) = 0.4391550328268399*m1rSelf[19]*mnuSelf+0.4*m1rSelf[16]*mnuSelf+0.447213595499958*m1rSelf[12]*mnuSelf; 
  data->AEM_S(23,18) = 0.4391550328268399*m1rSelf[16]*mnuSelf; 
  data->AEM_S(23,19) = 0.4391550328268399*m1rSelf[17]*mnuSelf; 
  data->AEM_S(24,10) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(24,11) = 0.4391550328268398*m1rSelf[18]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf; 
  data->AEM_S(24,12) = 0.5000000000000001*m1rSelf[16]*mnuSelf; 
  data->AEM_S(24,13) = 0.4472135954999579*m1rSelf[13]*mnuSelf; 
  data->AEM_S(24,14) = 0.31943828249997*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(24,16) = 0.31943828249997*m1rSelf[16]*mnuSelf+0.5000000000000001*m1rSelf[12]*mnuSelf; 
  data->AEM_S(24,17) = 0.4472135954999579*m1rSelf[17]*mnuSelf; 
  data->AEM_S(24,18) = 0.2981423969999719*m1rSelf[18]*mnuSelf+0.4391550328268398*m1rSelf[11]*mnuSelf; 
  data->AEM_S(25,10) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(25,11) = 0.5000000000000001*m1rSelf[17]*mnuSelf; 
  data->AEM_S(25,12) = 0.4391550328268398*m1rSelf[19]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf; 
  data->AEM_S(25,13) = 0.4472135954999579*m1rSelf[13]*mnuSelf; 
  data->AEM_S(25,15) = 0.31943828249997*m1rSelf[15]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(25,16) = 0.4472135954999579*m1rSelf[16]*mnuSelf; 
  data->AEM_S(25,17) = 0.31943828249997*m1rSelf[17]*mnuSelf+0.5000000000000001*m1rSelf[11]*mnuSelf; 
  data->AEM_S(25,19) = 0.2981423969999719*m1rSelf[19]*mnuSelf+0.4391550328268398*m1rSelf[12]*mnuSelf; 
  data->AEM_S(26,10) = 0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(26,11) = 0.447213595499958*m1rSelf[13]*mnuSelf; 
  data->AEM_S(26,12) = 0.5000000000000001*m1rSelf[14]*mnuSelf; 
  data->AEM_S(26,13) = 0.4391550328268399*m1rSelf[18]*mnuSelf+0.4*m1rSelf[17]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf; 
  data->AEM_S(26,14) = 0.31943828249997*m1rSelf[16]*mnuSelf+0.5000000000000001*m1rSelf[12]*mnuSelf; 
  data->AEM_S(26,15) = 0.4472135954999579*m1rSelf[16]*mnuSelf; 
  data->AEM_S(26,16) = 0.4472135954999579*m1rSelf[15]*mnuSelf+0.31943828249997*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(26,17) = 0.4*m1rSelf[13]*mnuSelf; 
  data->AEM_S(26,18) = 0.4391550328268399*m1rSelf[13]*mnuSelf; 
  data->AEM_S(27,10) = 0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(27,11) = 0.5000000000000001*m1rSelf[15]*mnuSelf; 
  data->AEM_S(27,12) = 0.447213595499958*m1rSelf[13]*mnuSelf; 
  data->AEM_S(27,13) = 0.4391550328268399*m1rSelf[19]*mnuSelf+0.4*m1rSelf[16]*mnuSelf+0.447213595499958*m1rSelf[12]*mnuSelf; 
  data->AEM_S(27,14) = 0.4472135954999579*m1rSelf[17]*mnuSelf; 
  data->AEM_S(27,15) = 0.31943828249997*m1rSelf[17]*mnuSelf+0.5000000000000001*m1rSelf[11]*mnuSelf; 
  data->AEM_S(27,16) = 0.4*m1rSelf[13]*mnuSelf; 
  data->AEM_S(27,17) = 0.31943828249997*m1rSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(27,19) = 0.4391550328268399*m1rSelf[13]*mnuSelf; 
  data->AEM_S(28,10) = 0.5*m1rSelf[18]*mnuSelf; 
  data->AEM_S(28,11) = 0.4391550328268398*m1rSelf[14]*mnuSelf; 
  data->AEM_S(28,13) = 0.4391550328268399*m1rSelf[16]*mnuSelf; 
  data->AEM_S(28,14) = 0.2981423969999719*m1rSelf[18]*mnuSelf+0.4391550328268398*m1rSelf[11]*mnuSelf; 
  data->AEM_S(28,16) = 0.4391550328268399*m1rSelf[13]*mnuSelf; 
  data->AEM_S(28,18) = 0.2981423969999719*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(29,10) = 0.5*m1rSelf[19]*mnuSelf; 
  data->AEM_S(29,12) = 0.4391550328268398*m1rSelf[15]*mnuSelf; 
  data->AEM_S(29,13) = 0.4391550328268399*m1rSelf[17]*mnuSelf; 
  data->AEM_S(29,15) = 0.2981423969999719*m1rSelf[19]*mnuSelf+0.4391550328268398*m1rSelf[12]*mnuSelf; 
  data->AEM_S(29,17) = 0.4391550328268399*m1rSelf[13]*mnuSelf; 
  data->AEM_S(29,19) = 0.2981423969999719*m1rSelf[15]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(20,40) = 0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(20,41) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(20,42) = 0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(20,43) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(20,44) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(20,45) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(20,46) = 0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(20,47) = 0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(20,48) = 0.5*m1rOther[18]*mnuOther; 
  data->AEM_S(20,49) = 0.5*m1rOther[19]*mnuOther; 
  data->AEM_S(21,40) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(21,41) = 0.4472135954999579*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(21,42) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(21,43) = 0.447213595499958*m1rOther[16]*mnuOther+0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(21,44) = 0.4391550328268398*m1rOther[18]*mnuOther+0.4472135954999579*m1rOther[11]*mnuOther; 
  data->AEM_S(21,45) = 0.5000000000000001*m1rOther[17]*mnuOther; 
  data->AEM_S(21,46) = 0.447213595499958*m1rOther[13]*mnuOther; 
  data->AEM_S(21,47) = 0.5000000000000001*m1rOther[15]*mnuOther; 
  data->AEM_S(21,48) = 0.4391550328268398*m1rOther[14]*mnuOther; 
  data->AEM_S(22,40) = 0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(22,41) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(22,42) = 0.4472135954999579*m1rOther[15]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(22,43) = 0.447213595499958*m1rOther[17]*mnuOther+0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(22,44) = 0.5000000000000001*m1rOther[16]*mnuOther; 
  data->AEM_S(22,45) = 0.4391550328268398*m1rOther[19]*mnuOther+0.4472135954999579*m1rOther[12]*mnuOther; 
  data->AEM_S(22,46) = 0.5000000000000001*m1rOther[14]*mnuOther; 
  data->AEM_S(22,47) = 0.447213595499958*m1rOther[13]*mnuOther; 
  data->AEM_S(22,49) = 0.4391550328268398*m1rOther[15]*mnuOther; 
  data->AEM_S(23,40) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(23,41) = 0.447213595499958*m1rOther[16]*mnuOther+0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(23,42) = 0.447213595499958*m1rOther[17]*mnuOther+0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(23,43) = 0.4472135954999579*m1rOther[15]*mnuOther+0.4472135954999579*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(23,44) = 0.4472135954999579*m1rOther[13]*mnuOther; 
  data->AEM_S(23,45) = 0.4472135954999579*m1rOther[13]*mnuOther; 
  data->AEM_S(23,46) = 0.4391550328268399*m1rOther[18]*mnuOther+0.4*m1rOther[17]*mnuOther+0.447213595499958*m1rOther[11]*mnuOther; 
  data->AEM_S(23,47) = 0.4391550328268399*m1rOther[19]*mnuOther+0.4*m1rOther[16]*mnuOther+0.447213595499958*m1rOther[12]*mnuOther; 
  data->AEM_S(23,48) = 0.4391550328268399*m1rOther[16]*mnuOther; 
  data->AEM_S(23,49) = 0.4391550328268399*m1rOther[17]*mnuOther; 
  data->AEM_S(24,40) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(24,41) = 0.4391550328268398*m1rOther[18]*mnuOther+0.4472135954999579*m1rOther[11]*mnuOther; 
  data->AEM_S(24,42) = 0.5000000000000001*m1rOther[16]*mnuOther; 
  data->AEM_S(24,43) = 0.4472135954999579*m1rOther[13]*mnuOther; 
  data->AEM_S(24,44) = 0.31943828249997*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(24,46) = 0.31943828249997*m1rOther[16]*mnuOther+0.5000000000000001*m1rOther[12]*mnuOther; 
  data->AEM_S(24,47) = 0.4472135954999579*m1rOther[17]*mnuOther; 
  data->AEM_S(24,48) = 0.2981423969999719*m1rOther[18]*mnuOther+0.4391550328268398*m1rOther[11]*mnuOther; 
  data->AEM_S(25,40) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(25,41) = 0.5000000000000001*m1rOther[17]*mnuOther; 
  data->AEM_S(25,42) = 0.4391550328268398*m1rOther[19]*mnuOther+0.4472135954999579*m1rOther[12]*mnuOther; 
  data->AEM_S(25,43) = 0.4472135954999579*m1rOther[13]*mnuOther; 
  data->AEM_S(25,45) = 0.31943828249997*m1rOther[15]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(25,46) = 0.4472135954999579*m1rOther[16]*mnuOther; 
  data->AEM_S(25,47) = 0.31943828249997*m1rOther[17]*mnuOther+0.5000000000000001*m1rOther[11]*mnuOther; 
  data->AEM_S(25,49) = 0.2981423969999719*m1rOther[19]*mnuOther+0.4391550328268398*m1rOther[12]*mnuOther; 
  data->AEM_S(26,40) = 0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(26,41) = 0.447213595499958*m1rOther[13]*mnuOther; 
  data->AEM_S(26,42) = 0.5000000000000001*m1rOther[14]*mnuOther; 
  data->AEM_S(26,43) = 0.4391550328268399*m1rOther[18]*mnuOther+0.4*m1rOther[17]*mnuOther+0.447213595499958*m1rOther[11]*mnuOther; 
  data->AEM_S(26,44) = 0.31943828249997*m1rOther[16]*mnuOther+0.5000000000000001*m1rOther[12]*mnuOther; 
  data->AEM_S(26,45) = 0.4472135954999579*m1rOther[16]*mnuOther; 
  data->AEM_S(26,46) = 0.4472135954999579*m1rOther[15]*mnuOther+0.31943828249997*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(26,47) = 0.4*m1rOther[13]*mnuOther; 
  data->AEM_S(26,48) = 0.4391550328268399*m1rOther[13]*mnuOther; 
  data->AEM_S(27,40) = 0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(27,41) = 0.5000000000000001*m1rOther[15]*mnuOther; 
  data->AEM_S(27,42) = 0.447213595499958*m1rOther[13]*mnuOther; 
  data->AEM_S(27,43) = 0.4391550328268399*m1rOther[19]*mnuOther+0.4*m1rOther[16]*mnuOther+0.447213595499958*m1rOther[12]*mnuOther; 
  data->AEM_S(27,44) = 0.4472135954999579*m1rOther[17]*mnuOther; 
  data->AEM_S(27,45) = 0.31943828249997*m1rOther[17]*mnuOther+0.5000000000000001*m1rOther[11]*mnuOther; 
  data->AEM_S(27,46) = 0.4*m1rOther[13]*mnuOther; 
  data->AEM_S(27,47) = 0.31943828249997*m1rOther[15]*mnuOther+0.4472135954999579*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(27,49) = 0.4391550328268399*m1rOther[13]*mnuOther; 
  data->AEM_S(28,40) = 0.5*m1rOther[18]*mnuOther; 
  data->AEM_S(28,41) = 0.4391550328268398*m1rOther[14]*mnuOther; 
  data->AEM_S(28,43) = 0.4391550328268399*m1rOther[16]*mnuOther; 
  data->AEM_S(28,44) = 0.2981423969999719*m1rOther[18]*mnuOther+0.4391550328268398*m1rOther[11]*mnuOther; 
  data->AEM_S(28,46) = 0.4391550328268399*m1rOther[13]*mnuOther; 
  data->AEM_S(28,48) = 0.2981423969999719*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(29,40) = 0.5*m1rOther[19]*mnuOther; 
  data->AEM_S(29,42) = 0.4391550328268398*m1rOther[15]*mnuOther; 
  data->AEM_S(29,43) = 0.4391550328268399*m1rOther[17]*mnuOther; 
  data->AEM_S(29,45) = 0.2981423969999719*m1rOther[19]*mnuOther+0.4391550328268398*m1rOther[12]*mnuOther; 
  data->AEM_S(29,47) = 0.4391550328268399*m1rOther[13]*mnuOther; 
  data->AEM_S(29,49) = 0.2981423969999719*m1rOther[15]*mnuOther+0.5*m1rOther[10]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[10] += m1rSelf[10]*mnuSelf+m1rOther[10]*mnuOther; 
  mnuM1sum[11] += m1rSelf[11]*mnuSelf+m1rOther[11]*mnuOther; 
  mnuM1sum[12] += m1rSelf[12]*mnuSelf+m1rOther[12]*mnuOther; 
  mnuM1sum[13] += m1rSelf[13]*mnuSelf+m1rOther[13]*mnuOther; 
  mnuM1sum[14] += m1rSelf[14]*mnuSelf+m1rOther[14]*mnuOther; 
  mnuM1sum[15] += m1rSelf[15]*mnuSelf+m1rOther[15]*mnuOther; 
  mnuM1sum[16] += m1rSelf[16]*mnuSelf+m1rOther[16]*mnuOther; 
  mnuM1sum[17] += m1rSelf[17]*mnuSelf+m1rOther[17]*mnuOther; 
  mnuM1sum[18] += m1rSelf[18]*mnuSelf+m1rOther[18]*mnuOther; 
  mnuM1sum[19] += m1rSelf[19]*mnuSelf+m1rOther[19]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(20,20) = m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(20,21) = m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(20,22) = m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(20,23) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(20,24) = m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(20,25) = m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(20,26) = m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(20,27) = m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(20,28) = m0rSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf; 
  data->AEM_S(20,29) = m0rSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf; 
  data->AEM_S(21,20) = m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(21,21) = 0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(21,22) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(21,23) = 0.8944271909999161*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(21,24) = 0.8783100656536796*m0rSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(21,25) = 1.0*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(21,26) = 0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(21,27) = 1.0*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(21,28) = 0.8783100656536796*m0rSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf; 
  data->AEM_S(22,20) = m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(22,21) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(22,22) = 0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(22,23) = 0.8944271909999161*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(22,24) = 1.0*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(22,25) = 0.8783100656536796*m0rSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(22,26) = 1.0*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(22,27) = 0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(22,29) = 0.8783100656536796*m0rSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf; 
  data->AEM_S(23,20) = m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(23,21) = 0.8944271909999161*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(23,22) = 0.8944271909999161*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(23,23) = 0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(23,24) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(23,25) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(23,26) = 0.8783100656536798*m0rSelf[8]*mnuSelf-0.4391550328268399*cESelf[8]*mnuSelf+0.8*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.8944271909999161*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(23,27) = 0.8783100656536798*m0rSelf[9]*mnuSelf-0.4391550328268399*cESelf[9]*mnuSelf+0.8*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.8944271909999161*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(23,28) = 0.8783100656536798*m0rSelf[6]*mnuSelf-0.4391550328268399*cESelf[6]*mnuSelf; 
  data->AEM_S(23,29) = 0.8783100656536798*m0rSelf[7]*mnuSelf-0.4391550328268399*cESelf[7]*mnuSelf; 
  data->AEM_S(24,20) = m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(24,21) = 0.8783100656536796*m0rSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(24,22) = 1.0*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(24,23) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(24,24) = 0.6388765649999399*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(24,26) = 0.6388765649999399*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.0*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(24,27) = 0.8944271909999159*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(24,28) = 0.5962847939999438*m0rSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+0.8783100656536796*m0rSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(25,20) = m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(25,21) = 1.0*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(25,22) = 0.8783100656536796*m0rSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(25,23) = 0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(25,25) = 0.6388765649999399*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(25,26) = 0.8944271909999159*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(25,27) = 0.6388765649999399*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.0*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(25,29) = 0.5962847939999438*m0rSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+0.8783100656536796*m0rSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(26,20) = m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(26,21) = 0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(26,22) = 1.0*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(26,23) = 0.8783100656536798*m0rSelf[8]*mnuSelf-0.4391550328268399*cESelf[8]*mnuSelf+0.8*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.8944271909999161*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(26,24) = 0.6388765649999399*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.0*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(26,25) = 0.8944271909999159*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(26,26) = 0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.6388765649999399*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(26,27) = 0.8*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(26,28) = 0.8783100656536798*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(27,20) = m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(27,21) = 1.0*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(27,22) = 0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(27,23) = 0.8783100656536798*m0rSelf[9]*mnuSelf-0.4391550328268399*cESelf[9]*mnuSelf+0.8*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.8944271909999161*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(27,24) = 0.8944271909999159*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(27,25) = 0.6388765649999399*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.0*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(27,26) = 0.8*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(27,27) = 0.6388765649999399*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(27,29) = 0.8783100656536798*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(28,20) = m0rSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf; 
  data->AEM_S(28,21) = 0.8783100656536796*m0rSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf; 
  data->AEM_S(28,23) = 0.8783100656536798*m0rSelf[6]*mnuSelf-0.4391550328268399*cESelf[6]*mnuSelf; 
  data->AEM_S(28,24) = 0.5962847939999438*m0rSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+0.8783100656536796*m0rSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(28,26) = 0.8783100656536798*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(28,28) = 0.5962847939999438*m0rSelf[4]*mnuSelf-0.2981423969999719*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(29,20) = m0rSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf; 
  data->AEM_S(29,22) = 0.8783100656536796*m0rSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf; 
  data->AEM_S(29,23) = 0.8783100656536798*m0rSelf[7]*mnuSelf-0.4391550328268399*cESelf[7]*mnuSelf; 
  data->AEM_S(29,25) = 0.5962847939999438*m0rSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+0.8783100656536796*m0rSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(29,27) = 0.8783100656536798*m0rSelf[3]*mnuSelf-0.4391550328268399*cESelf[3]*mnuSelf; 
  data->AEM_S(29,29) = 0.5962847939999438*m0rSelf[5]*mnuSelf-0.2981423969999719*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(20,50) = m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(20,51) = m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(20,52) = m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(20,53) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(20,54) = m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(20,55) = m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(20,56) = m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(20,57) = m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(20,58) = m0rOther[8]*mnuOther-0.5*cEOther[8]*mnuOther; 
  data->AEM_S(20,59) = m0rOther[9]*mnuOther-0.5*cEOther[9]*mnuOther; 
  data->AEM_S(21,50) = m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(21,51) = 0.8944271909999159*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(21,52) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(21,53) = 0.8944271909999161*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(21,54) = 0.8783100656536796*m0rOther[8]*mnuOther-0.4391550328268398*cEOther[8]*mnuOther+0.8944271909999159*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(21,55) = 1.0*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(21,56) = 0.8944271909999161*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(21,57) = 1.0*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(21,58) = 0.8783100656536796*m0rOther[4]*mnuOther-0.4391550328268398*cEOther[4]*mnuOther; 
  data->AEM_S(22,50) = m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(22,51) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(22,52) = 0.8944271909999159*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(22,53) = 0.8944271909999161*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(22,54) = 1.0*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(22,55) = 0.8783100656536796*m0rOther[9]*mnuOther-0.4391550328268398*cEOther[9]*mnuOther+0.8944271909999159*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(22,56) = 1.0*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(22,57) = 0.8944271909999161*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(22,59) = 0.8783100656536796*m0rOther[5]*mnuOther-0.4391550328268398*cEOther[5]*mnuOther; 
  data->AEM_S(23,50) = m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(23,51) = 0.8944271909999161*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(23,52) = 0.8944271909999161*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(23,53) = 0.8944271909999159*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+0.8944271909999159*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(23,54) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(23,55) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(23,56) = 0.8783100656536798*m0rOther[8]*mnuOther-0.4391550328268399*cEOther[8]*mnuOther+0.8*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+0.8944271909999161*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(23,57) = 0.8783100656536798*m0rOther[9]*mnuOther-0.4391550328268399*cEOther[9]*mnuOther+0.8*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+0.8944271909999161*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(23,58) = 0.8783100656536798*m0rOther[6]*mnuOther-0.4391550328268399*cEOther[6]*mnuOther; 
  data->AEM_S(23,59) = 0.8783100656536798*m0rOther[7]*mnuOther-0.4391550328268399*cEOther[7]*mnuOther; 
  data->AEM_S(24,50) = m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(24,51) = 0.8783100656536796*m0rOther[8]*mnuOther-0.4391550328268398*cEOther[8]*mnuOther+0.8944271909999159*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(24,52) = 1.0*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(24,53) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(24,54) = 0.6388765649999399*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(24,56) = 0.6388765649999399*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.0*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(24,57) = 0.8944271909999159*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(24,58) = 0.5962847939999438*m0rOther[8]*mnuOther-0.2981423969999719*cEOther[8]*mnuOther+0.8783100656536796*m0rOther[1]*mnuOther-0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(25,50) = m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(25,51) = 1.0*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(25,52) = 0.8783100656536796*m0rOther[9]*mnuOther-0.4391550328268398*cEOther[9]*mnuOther+0.8944271909999159*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(25,53) = 0.8944271909999159*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(25,55) = 0.6388765649999399*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(25,56) = 0.8944271909999159*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(25,57) = 0.6388765649999399*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.0*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(25,59) = 0.5962847939999438*m0rOther[9]*mnuOther-0.2981423969999719*cEOther[9]*mnuOther+0.8783100656536796*m0rOther[2]*mnuOther-0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(26,50) = m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(26,51) = 0.8944271909999161*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(26,52) = 1.0*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(26,53) = 0.8783100656536798*m0rOther[8]*mnuOther-0.4391550328268399*cEOther[8]*mnuOther+0.8*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+0.8944271909999161*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(26,54) = 0.6388765649999399*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.0*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(26,55) = 0.8944271909999159*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(26,56) = 0.8944271909999159*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+0.6388765649999399*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(26,57) = 0.8*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(26,58) = 0.8783100656536798*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(27,50) = m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(27,51) = 1.0*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(27,52) = 0.8944271909999161*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(27,53) = 0.8783100656536798*m0rOther[9]*mnuOther-0.4391550328268399*cEOther[9]*mnuOther+0.8*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+0.8944271909999161*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(27,54) = 0.8944271909999159*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(27,55) = 0.6388765649999399*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.0*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(27,56) = 0.8*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(27,57) = 0.6388765649999399*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+0.8944271909999159*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(27,59) = 0.8783100656536798*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(28,50) = m0rOther[8]*mnuOther-0.5*cEOther[8]*mnuOther; 
  data->AEM_S(28,51) = 0.8783100656536796*m0rOther[4]*mnuOther-0.4391550328268398*cEOther[4]*mnuOther; 
  data->AEM_S(28,53) = 0.8783100656536798*m0rOther[6]*mnuOther-0.4391550328268399*cEOther[6]*mnuOther; 
  data->AEM_S(28,54) = 0.5962847939999438*m0rOther[8]*mnuOther-0.2981423969999719*cEOther[8]*mnuOther+0.8783100656536796*m0rOther[1]*mnuOther-0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(28,56) = 0.8783100656536798*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(28,58) = 0.5962847939999438*m0rOther[4]*mnuOther-0.2981423969999719*cEOther[4]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(29,50) = m0rOther[9]*mnuOther-0.5*cEOther[9]*mnuOther; 
  data->AEM_S(29,52) = 0.8783100656536796*m0rOther[5]*mnuOther-0.4391550328268398*cEOther[5]*mnuOther; 
  data->AEM_S(29,53) = 0.8783100656536798*m0rOther[7]*mnuOther-0.4391550328268399*cEOther[7]*mnuOther; 
  data->AEM_S(29,55) = 0.5962847939999438*m0rOther[9]*mnuOther-0.2981423969999719*cEOther[9]*mnuOther+0.8783100656536796*m0rOther[2]*mnuOther-0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(29,57) = 0.8783100656536798*m0rOther[3]*mnuOther-0.4391550328268399*cEOther[3]*mnuOther; 
  data->AEM_S(29,59) = 0.5962847939999438*m0rOther[5]*mnuOther-0.2981423969999719*cEOther[5]*mnuOther+m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
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
 
  // Set other entries to 0. // 
  data->AEM_S.block<10,10>(0,10).setZero(); 
  data->AEM_S.block<10,10>(10,0).setZero(); 
  data->AEM_S.block<10,10>(0,40).setZero(); 
  data->AEM_S.block<10,10>(10,30).setZero(); 
 
  double m1Relax[20]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<20; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(30,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(30,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(30,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(30,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(30,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(30,6) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(30,7) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(30,8) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(30,9) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(31,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(31,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,3) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,4) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,5) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(31,6) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,7) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(31,8) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(32,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(32,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(32,3) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(32,4) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(32,5) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,6) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(32,7) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(32,9) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(33,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,1) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,2) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(33,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,6) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,7) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,8) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(33,9) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(34,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(34,1) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(34,2) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(34,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(34,6) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(34,7) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(34,8) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(35,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(35,1) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(35,2) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(35,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(35,6) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(35,7) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(35,9) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,0) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(36,1) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(36,2) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(36,3) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(36,4) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,5) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(36,6) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(36,7) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(36,8) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(37,0) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(37,1) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(37,2) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(37,3) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(37,4) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(37,5) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(37,6) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(37,7) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(37,9) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(38,0) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(38,1) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(38,3) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(38,4) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(38,6) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(38,8) = 0.2981423969999719*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(39,0) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(39,2) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(39,3) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(39,5) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(39,7) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(39,9) = 0.2981423969999719*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(30,20) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(30,21) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(30,22) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(30,23) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(30,24) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(30,25) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(30,26) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(30,27) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(30,28) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(30,29) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(31,20) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(31,21) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(31,22) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(31,23) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(31,24) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(31,25) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(31,26) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(31,27) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(31,28) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(32,20) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(32,21) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(32,22) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(32,23) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(32,24) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(32,25) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(32,26) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(32,27) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(32,29) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(33,20) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(33,21) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(33,22) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(33,23) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(33,24) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(33,25) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(33,26) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(33,27) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(33,28) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(33,29) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(34,20) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(34,21) = (-0.4391550328268398*cMSelf[8]*mnuSelf)-0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(34,22) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(34,23) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(34,24) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(34,26) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(34,27) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(34,28) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(35,20) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(35,21) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(35,22) = (-0.4391550328268398*cMSelf[9]*mnuSelf)-0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(35,23) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(35,25) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(35,26) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(35,27) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(35,29) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(36,20) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(36,21) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(36,22) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(36,23) = (-0.4391550328268399*cMSelf[8]*mnuSelf)-0.4*cMSelf[7]*mnuSelf-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(36,24) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(36,25) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(36,26) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(36,27) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(36,28) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(37,20) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(37,21) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(37,22) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(37,23) = (-0.4391550328268399*cMSelf[9]*mnuSelf)-0.4*cMSelf[6]*mnuSelf-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(37,24) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(37,25) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(37,26) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(37,27) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(37,29) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(38,20) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(38,21) = -0.4391550328268398*cMSelf[4]*mnuSelf; 
  data->AEM_S(38,23) = -0.4391550328268399*cMSelf[6]*mnuSelf; 
  data->AEM_S(38,24) = (-0.2981423969999719*cMSelf[8]*mnuSelf)-0.4391550328268398*cMSelf[1]*mnuSelf; 
  data->AEM_S(38,26) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(38,28) = (-0.2981423969999719*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(39,20) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(39,22) = -0.4391550328268398*cMSelf[5]*mnuSelf; 
  data->AEM_S(39,23) = -0.4391550328268399*cMSelf[7]*mnuSelf; 
  data->AEM_S(39,25) = (-0.2981423969999719*cMSelf[9]*mnuSelf)-0.4391550328268398*cMSelf[2]*mnuSelf; 
  data->AEM_S(39,27) = -0.4391550328268399*cMSelf[3]*mnuSelf; 
  data->AEM_S(39,29) = (-0.2981423969999719*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(30,30) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(30,31) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(30,32) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(30,33) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(30,34) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(30,35) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(30,36) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(30,37) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(30,38) = -0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(30,39) = -0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(31,30) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(31,31) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(31,32) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(31,33) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(31,34) = (-0.4391550328268398*m0rOther[8]*mnuOther)-0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(31,35) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(31,36) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(31,37) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(31,38) = -0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(32,30) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(32,31) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(32,32) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(32,33) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(32,34) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(32,35) = (-0.4391550328268398*m0rOther[9]*mnuOther)-0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(32,36) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(32,37) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(32,39) = -0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(33,30) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(33,31) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(33,32) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(33,33) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(33,34) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(33,35) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(33,36) = (-0.4391550328268399*m0rOther[8]*mnuOther)-0.4*m0rOther[7]*mnuOther-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(33,37) = (-0.4391550328268399*m0rOther[9]*mnuOther)-0.4*m0rOther[6]*mnuOther-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(33,38) = -0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(33,39) = -0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(34,30) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(34,31) = (-0.4391550328268398*m0rOther[8]*mnuOther)-0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(34,32) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(34,33) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(34,34) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(34,36) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(34,37) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(34,38) = (-0.2981423969999719*m0rOther[8]*mnuOther)-0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(35,30) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(35,31) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(35,32) = (-0.4391550328268398*m0rOther[9]*mnuOther)-0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(35,33) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(35,35) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(35,36) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(35,37) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(35,39) = (-0.2981423969999719*m0rOther[9]*mnuOther)-0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(36,30) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(36,31) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(36,32) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(36,33) = (-0.4391550328268399*m0rOther[8]*mnuOther)-0.4*m0rOther[7]*mnuOther-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(36,34) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(36,35) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(36,36) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.31943828249997*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(36,37) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(36,38) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(37,30) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(37,31) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(37,32) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(37,33) = (-0.4391550328268399*m0rOther[9]*mnuOther)-0.4*m0rOther[6]*mnuOther-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(37,34) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(37,35) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(37,36) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(37,37) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(37,39) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(38,30) = -0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(38,31) = -0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(38,33) = -0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(38,34) = (-0.2981423969999719*m0rOther[8]*mnuOther)-0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(38,36) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(38,38) = (-0.2981423969999719*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(39,30) = -0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(39,32) = -0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(39,33) = -0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(39,35) = (-0.2981423969999719*m0rOther[9]*mnuOther)-0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(39,37) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(39,39) = (-0.2981423969999719*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(30,50) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(30,51) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(30,52) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(30,53) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(30,54) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(30,55) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(30,56) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(30,57) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(30,58) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(30,59) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(31,50) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(31,51) = 0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(31,52) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(31,53) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(31,54) = 0.4391550328268398*cMOther[8]*mnuOther+0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(31,55) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(31,56) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(31,57) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(31,58) = 0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(32,50) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(32,51) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(32,52) = 0.4472135954999579*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(32,53) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(32,54) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(32,55) = 0.4391550328268398*cMOther[9]*mnuOther+0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(32,56) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(32,57) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(32,59) = 0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(33,50) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(33,51) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(33,52) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(33,53) = 0.4472135954999579*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(33,54) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(33,55) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(33,56) = 0.4391550328268399*cMOther[8]*mnuOther+0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(33,57) = 0.4391550328268399*cMOther[9]*mnuOther+0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(33,58) = 0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(33,59) = 0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(34,50) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(34,51) = 0.4391550328268398*cMOther[8]*mnuOther+0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(34,52) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(34,53) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(34,54) = 0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(34,56) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(34,57) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(34,58) = 0.2981423969999719*cMOther[8]*mnuOther+0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(35,50) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(35,51) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(35,52) = 0.4391550328268398*cMOther[9]*mnuOther+0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(35,53) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(35,55) = 0.31943828249997*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(35,56) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(35,57) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(35,59) = 0.2981423969999719*cMOther[9]*mnuOther+0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(36,50) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(36,51) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(36,52) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(36,53) = 0.4391550328268399*cMOther[8]*mnuOther+0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(36,54) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(36,55) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(36,56) = 0.4472135954999579*cMOther[5]*mnuOther+0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(36,57) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(36,58) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(37,50) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(37,51) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(37,52) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(37,53) = 0.4391550328268399*cMOther[9]*mnuOther+0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(37,54) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(37,55) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(37,56) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(37,57) = 0.31943828249997*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(37,59) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(38,50) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(38,51) = 0.4391550328268398*cMOther[4]*mnuOther; 
  data->AEM_S(38,53) = 0.4391550328268399*cMOther[6]*mnuOther; 
  data->AEM_S(38,54) = 0.2981423969999719*cMOther[8]*mnuOther+0.4391550328268398*cMOther[1]*mnuOther; 
  data->AEM_S(38,56) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(38,58) = 0.2981423969999719*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(39,50) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(39,52) = 0.4391550328268398*cMOther[5]*mnuOther; 
  data->AEM_S(39,53) = 0.4391550328268399*cMOther[7]*mnuOther; 
  data->AEM_S(39,55) = 0.2981423969999719*cMOther[9]*mnuOther+0.4391550328268398*cMOther[2]*mnuOther; 
  data->AEM_S(39,57) = 0.4391550328268399*cMOther[3]*mnuOther; 
  data->AEM_S(39,59) = 0.2981423969999719*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(50,0) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.25*m0rSelf[7]*uSelf[7]*mnuSelf-0.25*m0rSelf[6]*uSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(50,1) = (-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(50,2) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(50,3) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(50,4) = (-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(50,5) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[2]*uSelf[9]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(50,6) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(50,7) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(50,8) = (-0.149071198499986*m0rSelf[4]*uSelf[8]*mnuSelf)-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.149071198499986*uSelf[4]*m0rSelf[8]*mnuSelf-0.25*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(50,9) = (-0.149071198499986*m0rSelf[5]*uSelf[9]*mnuSelf)-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.149071198499986*uSelf[5]*m0rSelf[9]*mnuSelf-0.25*uSelf[0]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(51,0) = (-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(51,1) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3833333333333334*m0rSelf[8]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[8]*mnuSelf-0.45*m0rSelf[7]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(51,2) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(51,3) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(51,4) = (-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(51,5) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(51,6) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(51,7) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[9]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(51,8) = (-0.2499586742703185*m0rSelf[8]*uSelf[8]*mnuSelf)-0.3833333333333334*m0rSelf[1]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[1]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[4]*mnuSelf+0.4391550328268398*m1rSelf[4]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(51,9) = (-0.149071198499986*m0rSelf[7]*uSelf[9]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.149071198499986*uSelf[7]*m0rSelf[9]*mnuSelf-0.25*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(52,0) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(52,1) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(52,2) = (-0.3833333333333334*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[2]*uSelf[9]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[9]*mnuSelf-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.45*m0rSelf[6]*uSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(52,3) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(52,4) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(52,5) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[9]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(52,6) = (-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(52,7) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.21957751641342*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(52,8) = (-0.149071198499986*m0rSelf[6]*uSelf[8]*mnuSelf)-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.149071198499986*uSelf[6]*m0rSelf[8]*mnuSelf-0.25*uSelf[2]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(52,9) = (-0.2499586742703185*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3833333333333334*m0rSelf[2]*uSelf[9]*mnuSelf-0.3833333333333334*uSelf[2]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[5]*mnuSelf+0.4391550328268398*m1rSelf[5]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(53,0) = (-0.2195775164134199*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[6]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[6]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,1) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,2) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(53,3) = (-0.3833333333333334*m0rSelf[9]*uSelf[9]*mnuSelf)-0.175662013130736*m0rSelf[6]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[9]*mnuSelf-0.175662013130736*uSelf[6]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[8]*mnuSelf-0.175662013130736*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[8]*mnuSelf-0.175662013130736*uSelf[7]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[8]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[6]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.2*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.2*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.81*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(53,4) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,5) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,6) = (-0.175662013130736*m0rSelf[3]*uSelf[9]*mnuSelf)-0.175662013130736*uSelf[3]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[8]*mnuSelf+0.43915503282684*m1rSelf[8]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(53,7) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[4]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[9]*mnuSelf+0.43915503282684*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[9]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[8]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[8]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,8) = (-0.3833333333333334*m0rSelf[3]*uSelf[8]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[8]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[7]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[6]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[6]*mnuSelf+0.43915503282684*m1rSelf[6]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[6]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[6]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(53,9) = (-0.3833333333333334*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[7]*mnuSelf+0.43915503282684*m1rSelf[7]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[7]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[7]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[6]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(54,0) = (-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(54,1) = (-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[8]*mnuSelf+0.4391550328268398*m1rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(54,2) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(54,3) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(54,4) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.4621212121212121*m0rSelf[8]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[6]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(54,5) = (-0.2195775164134199*m0rSelf[6]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[6]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[7]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(54,6) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(54,7) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[8]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(54,8) = (-0.4621212121212121*m0rSelf[4]*uSelf[8]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[8]*mnuSelf+0.2981423969999719*m1rSelf[8]*mnuSelf-0.4621212121212121*uSelf[4]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[5]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[1]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(54,9) = (-0.25*m0rSelf[4]*uSelf[9]*mnuSelf)-0.25*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[7]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[6]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(55,0) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[2]*uSelf[9]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(55,1) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(55,2) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[9]*mnuSelf+0.4391550328268398*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[9]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(55,3) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(55,4) = (-0.2195775164134199*m0rSelf[6]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[6]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[7]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(55,5) = (-0.4621212121212121*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3273268353539885*m0rSelf[2]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[9]*mnuSelf-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(55,6) = (-0.21957751641342*m0rSelf[4]*uSelf[9]*mnuSelf)-0.21957751641342*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(55,7) = (-0.3273268353539885*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3273268353539885*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(55,8) = (-0.25*m0rSelf[5]*uSelf[8]*mnuSelf)-0.25*uSelf[5]*m0rSelf[8]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(55,9) = (-0.4621212121212121*m0rSelf[5]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[9]*mnuSelf+0.2981423969999719*m1rSelf[9]*mnuSelf-0.4621212121212121*uSelf[5]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[7]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[2]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,0) = (-0.2195775164134199*m0rSelf[3]*uSelf[8]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[8]*mnuSelf-0.2*m0rSelf[3]*uSelf[7]*mnuSelf-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(56,1) = (-0.1963961012123931*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[7]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[6]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,2) = (-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(56,3) = (-0.175662013130736*m0rSelf[3]*uSelf[9]*mnuSelf)-0.175662013130736*uSelf[3]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[8]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[8]*mnuSelf+0.43915503282684*m1rSelf[8]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[8]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[8]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(56,4) = (-0.2195775164134199*m0rSelf[5]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[5]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[8]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,5) = (-0.21957751641342*m0rSelf[4]*uSelf[9]*mnuSelf)-0.21957751641342*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[8]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(56,6) = (-0.3833333333333334*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1254728665219542*m0rSelf[6]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[9]*mnuSelf-0.1254728665219542*uSelf[6]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[9]*mnuSelf-0.4621212121212121*m0rSelf[8]*uSelf[8]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[8]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[8]*mnuSelf-0.29277002188456*uSelf[7]*m0rSelf[8]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[8]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[7]*mnuSelf-0.9642857142857143*m0rSelf[6]*uSelf[6]*mnuSelf-0.2874944542499729*m0rSelf[2]*uSelf[6]*mnuSelf-0.2874944542499729*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(56,7) = (-0.1928571428571429*m0rSelf[8]*uSelf[9]*mnuSelf)-0.29277002188456*m0rSelf[7]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[9]*mnuSelf-0.1928571428571429*uSelf[8]*m0rSelf[9]*mnuSelf-0.29277002188456*uSelf[7]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[9]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[8]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[8]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,8) = (-0.1928571428571429*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[7]*m0rSelf[9]*mnuSelf-0.4621212121212121*m0rSelf[6]*uSelf[8]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[8]*mnuSelf-0.4621212121212121*uSelf[6]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[2]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[5]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[5]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,9) = (-0.3833333333333334*m0rSelf[6]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[6]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[7]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[7]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.1254728665219543*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(57,0) = (-0.2195775164134199*m0rSelf[3]*uSelf[9]*mnuSelf)-0.2195775164134199*uSelf[3]*m0rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(57,1) = (-0.149071198499986*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[6]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[9]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[2]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[7]*m0rSelf[8]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(57,2) = (-0.3273268353539885*m0rSelf[7]*uSelf[9]*mnuSelf)-0.21957751641342*m0rSelf[1]*uSelf[9]*mnuSelf-0.3273268353539885*uSelf[7]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[8]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(57,3) = (-0.3273268353539885*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1963961012123931*m0rSelf[4]*uSelf[9]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[9]*mnuSelf+0.43915503282684*m1rSelf[9]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[9]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[9]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[8]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[8]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(57,4) = (-0.1963961012123931*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1963961012123931*uSelf[3]*m0rSelf[9]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[8]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[8]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(57,5) = (-0.3273268353539885*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3273268353539885*uSelf[3]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[8]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[8]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(57,6) = (-0.1928571428571429*m0rSelf[8]*uSelf[9]*mnuSelf)-0.29277002188456*m0rSelf[7]*uSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[9]*mnuSelf-0.1928571428571429*uSelf[8]*m0rSelf[9]*mnuSelf-0.29277002188456*uSelf[7]*m0rSelf[9]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[9]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[8]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[8]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(57,7) = (-0.4621212121212121*m0rSelf[9]*uSelf[9]*mnuSelf)-0.29277002188456*m0rSelf[6]*uSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[9]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[9]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[8]*mnuSelf-0.1254728665219542*m0rSelf[7]*uSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[8]*mnuSelf-0.1254728665219542*uSelf[7]*m0rSelf[8]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[8]*mnuSelf-0.9642857142857143*m0rSelf[7]*uSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[1]*uSelf[7]*mnuSelf-0.2874944542499729*uSelf[1]*m0rSelf[7]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[6]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(57,8) = (-0.1928571428571429*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[6]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[7]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[7]*m0rSelf[8]*mnuSelf-0.1254728665219543*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(57,9) = (-0.4621212121212121*m0rSelf[7]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[1]*uSelf[9]*mnuSelf-0.4621212121212121*uSelf[7]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[1]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[4]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(58,0) = (-0.149071198499986*m0rSelf[4]*uSelf[8]*mnuSelf)-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.149071198499986*uSelf[4]*m0rSelf[8]*mnuSelf-0.25*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(58,1) = (-0.2499586742703185*m0rSelf[8]*uSelf[8]*mnuSelf)-0.3833333333333334*m0rSelf[1]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[1]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[4]*mnuSelf+0.4391550328268398*m1rSelf[4]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(58,2) = (-0.149071198499986*m0rSelf[6]*uSelf[8]*mnuSelf)-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.149071198499986*uSelf[6]*m0rSelf[8]*mnuSelf-0.25*uSelf[2]*m0rSelf[8]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(58,3) = (-0.3833333333333334*m0rSelf[3]*uSelf[8]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[8]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[7]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[6]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[6]*mnuSelf+0.43915503282684*m1rSelf[6]*mnuSelf-0.1963961012123931*uSelf[5]*m0rSelf[6]*mnuSelf-0.3273268353539885*uSelf[4]*m0rSelf[6]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[4]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[4]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(58,4) = (-0.4621212121212121*m0rSelf[4]*uSelf[8]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[8]*mnuSelf+0.2981423969999719*m1rSelf[8]*mnuSelf-0.4621212121212121*uSelf[4]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[8]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[5]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[4]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[1]*mnuSelf+0.4391550328268398*m1rSelf[1]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(58,5) = (-0.25*m0rSelf[5]*uSelf[8]*mnuSelf)-0.25*uSelf[5]*m0rSelf[8]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(58,6) = (-0.1928571428571429*m0rSelf[7]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[7]*m0rSelf[9]*mnuSelf-0.4621212121212121*m0rSelf[6]*uSelf[8]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[8]*mnuSelf-0.4621212121212121*uSelf[6]*m0rSelf[8]*mnuSelf-0.149071198499986*uSelf[2]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[6]*mnuSelf-0.3273268353539885*uSelf[1]*m0rSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[5]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[5]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[4]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(58,7) = (-0.1928571428571429*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1928571428571429*uSelf[6]*m0rSelf[9]*mnuSelf-0.3833333333333334*m0rSelf[7]*uSelf[8]*mnuSelf-0.3833333333333334*uSelf[7]*m0rSelf[8]*mnuSelf-0.1254728665219543*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(58,8) = (-0.25*m0rSelf[9]*uSelf[9]*mnuSelf)-0.5898601398601399*m0rSelf[8]*uSelf[8]*mnuSelf-0.2499586742703185*m0rSelf[1]*uSelf[8]*mnuSelf-0.2499586742703185*uSelf[1]*m0rSelf[8]*mnuSelf-0.3833333333333334*m0rSelf[7]*uSelf[7]*mnuSelf-0.4621212121212121*m0rSelf[6]*uSelf[6]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[6]*mnuSelf-0.149071198499986*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.4621212121212121*m0rSelf[4]*uSelf[4]*mnuSelf-0.149071198499986*m0rSelf[0]*uSelf[4]*mnuSelf+0.2981423969999719*m1rSelf[4]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[4]*mnuSelf-0.3833333333333334*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3833333333333334*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(58,9) = (-0.25*m0rSelf[8]*uSelf[9]*mnuSelf)-0.25*uSelf[8]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[7]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[7]*mnuSelf; 
  data->AEM_S(59,0) = (-0.149071198499986*m0rSelf[5]*uSelf[9]*mnuSelf)-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.149071198499986*uSelf[5]*m0rSelf[9]*mnuSelf-0.25*uSelf[0]*m0rSelf[9]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(59,1) = (-0.149071198499986*m0rSelf[7]*uSelf[9]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.149071198499986*uSelf[7]*m0rSelf[9]*mnuSelf-0.25*uSelf[1]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[6]*m0rSelf[7]*mnuSelf-0.2195775164134199*uSelf[2]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[3]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(59,2) = (-0.2499586742703185*m0rSelf[9]*uSelf[9]*mnuSelf)-0.3833333333333334*m0rSelf[2]*uSelf[9]*mnuSelf-0.3833333333333334*uSelf[2]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[7]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[5]*mnuSelf+0.4391550328268398*m1rSelf[5]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[3]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(59,3) = (-0.3833333333333334*m0rSelf[3]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[3]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[7]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[7]*mnuSelf+0.43915503282684*m1rSelf[7]*mnuSelf-0.3273268353539885*uSelf[5]*m0rSelf[7]*mnuSelf-0.1963961012123931*uSelf[4]*m0rSelf[7]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[7]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[6]*mnuSelf-0.175662013130736*uSelf[3]*m0rSelf[6]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[5]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[3]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(59,4) = (-0.25*m0rSelf[4]*uSelf[9]*mnuSelf)-0.25*uSelf[4]*m0rSelf[9]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[7]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[6]*mnuSelf-0.21957751641342*uSelf[5]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(59,5) = (-0.4621212121212121*m0rSelf[5]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[9]*mnuSelf+0.2981423969999719*m1rSelf[9]*mnuSelf-0.4621212121212121*uSelf[5]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[9]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[7]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[7]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[6]*mnuSelf-0.2195775164134199*uSelf[4]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[5]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[3]*mnuSelf-0.2195775164134199*uSelf[1]*m0rSelf[3]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[2]*mnuSelf+0.4391550328268398*m1rSelf[2]*mnuSelf-0.2195775164134199*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(59,6) = (-0.3833333333333334*m0rSelf[6]*uSelf[9]*mnuSelf)-0.3833333333333334*uSelf[6]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[7]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[7]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[7]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[7]*mnuSelf-0.1254728665219543*m0rSelf[6]*uSelf[6]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[2]*m0rSelf[6]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[5]*mnuSelf-0.21957751641342*uSelf[4]*m0rSelf[5]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(59,7) = (-0.4621212121212121*m0rSelf[7]*uSelf[9]*mnuSelf)-0.149071198499986*m0rSelf[1]*uSelf[9]*mnuSelf-0.4621212121212121*uSelf[7]*m0rSelf[9]*mnuSelf-0.149071198499986*uSelf[1]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[8]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[8]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[7]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[7]*mnuSelf-0.29277002188456*uSelf[6]*m0rSelf[7]*mnuSelf-0.3273268353539885*uSelf[2]*m0rSelf[7]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[6]*mnuSelf-0.1963961012123931*uSelf[1]*m0rSelf[6]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[5]*mnuSelf-0.3273268353539885*uSelf[3]*m0rSelf[5]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[4]*mnuSelf-0.1963961012123931*uSelf[3]*m0rSelf[4]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[3]*mnuSelf+0.43915503282684*m1rSelf[3]*mnuSelf-0.21957751641342*uSelf[0]*m0rSelf[3]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[2]*mnuSelf-0.21957751641342*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(59,8) = (-0.25*m0rSelf[8]*uSelf[9]*mnuSelf)-0.25*uSelf[8]*m0rSelf[9]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[7]*mnuSelf-0.1928571428571429*uSelf[6]*m0rSelf[7]*mnuSelf; 
  data->AEM_S(59,9) = (-0.5898601398601399*m0rSelf[9]*uSelf[9]*mnuSelf)-0.2499586742703185*m0rSelf[2]*uSelf[9]*mnuSelf-0.2499586742703185*uSelf[2]*m0rSelf[9]*mnuSelf-0.25*m0rSelf[8]*uSelf[8]*mnuSelf-0.4621212121212121*m0rSelf[7]*uSelf[7]*mnuSelf-0.149071198499986*m0rSelf[1]*uSelf[7]*mnuSelf-0.149071198499986*uSelf[1]*m0rSelf[7]*mnuSelf-0.3833333333333334*m0rSelf[6]*uSelf[6]*mnuSelf-0.4621212121212121*m0rSelf[5]*uSelf[5]*mnuSelf-0.149071198499986*m0rSelf[0]*uSelf[5]*mnuSelf+0.2981423969999719*m1rSelf[5]*mnuSelf-0.149071198499986*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3833333333333334*m0rSelf[3]*uSelf[3]*mnuSelf-0.3833333333333334*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(50,30) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.25*m0rOther[7]*uOther[7]*mnuOther+0.25*m0rOther[6]*uOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(50,31) = 0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(50,32) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(50,33) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(50,34) = 0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[8]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(50,35) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[9]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(50,36) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(50,37) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(50,38) = 0.149071198499986*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.149071198499986*uOther[4]*m0rOther[8]*mnuOther+0.25*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[6]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[4]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[4]*mnuOther; 
  data->AEM_S(50,39) = 0.149071198499986*m0rOther[5]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.149071198499986*uOther[5]*m0rOther[9]*mnuOther+0.25*uOther[0]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[7]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[5]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[5]*mnuOther; 
  data->AEM_S(51,30) = 0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(51,31) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[8]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[8]*mnuOther+0.45*m0rOther[7]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(51,32) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(51,33) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(51,34) = 0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[8]*mnuOther-0.4391550328268398*m1rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(51,35) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(51,36) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.21957751641342*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.21957751641342*uOther[2]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(51,37) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.21957751641342*m0rOther[2]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.21957751641342*uOther[2]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(51,38) = 0.2499586742703185*m0rOther[8]*uOther[8]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[8]*mnuOther+0.3833333333333334*uOther[1]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[6]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[4]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[4]*mnuOther-0.4391550328268398*m1rOther[4]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(51,39) = 0.149071198499986*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.149071198499986*uOther[7]*m0rOther[9]*mnuOther+0.25*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[5]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[5]*mnuOther; 
  data->AEM_S(52,30) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(52,31) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(52,32) = 0.3833333333333334*m0rOther[9]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[9]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.45*m0rOther[6]*uOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(52,33) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(52,34) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(52,35) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[9]*mnuOther-0.4391550328268398*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[9]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(52,36) = 0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.21957751641342*m0rOther[1]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.21957751641342*uOther[1]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(52,37) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.21957751641342*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.21957751641342*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(52,38) = 0.149071198499986*m0rOther[6]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.149071198499986*uOther[6]*m0rOther[8]*mnuOther+0.25*uOther[2]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[6]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[4]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(52,39) = 0.2499586742703185*m0rOther[9]*uOther[9]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[9]*mnuOther+0.3833333333333334*uOther[2]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[7]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[5]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[5]*mnuOther-0.4391550328268398*m1rOther[5]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(53,30) = 0.2195775164134199*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[8]*mnuOther+0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(53,31) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(53,32) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(53,33) = 0.3833333333333334*m0rOther[9]*uOther[9]*mnuOther+0.175662013130736*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[9]*mnuOther+0.175662013130736*uOther[6]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[8]*mnuOther+0.175662013130736*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[8]*mnuOther+0.175662013130736*uOther[7]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[8]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[7]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[7]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[6]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[6]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.2*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.2*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(53,34) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(53,35) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(53,36) = 0.175662013130736*m0rOther[3]*uOther[9]*mnuOther+0.175662013130736*uOther[3]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.21957751641342*m0rOther[0]*uOther[8]*mnuOther-0.43915503282684*m1rOther[8]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.21957751641342*uOther[0]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(53,37) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*m0rOther[0]*uOther[9]*mnuOther-0.43915503282684*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[9]*mnuOther+0.21957751641342*uOther[0]*m0rOther[9]*mnuOther+0.175662013130736*m0rOther[3]*uOther[8]*mnuOther+0.175662013130736*uOther[3]*m0rOther[8]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(53,38) = 0.3833333333333334*m0rOther[3]*uOther[8]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[8]*mnuOther+0.175662013130736*m0rOther[3]*uOther[7]*mnuOther+0.175662013130736*uOther[3]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[6]*mnuOther+0.21957751641342*m0rOther[0]*uOther[6]*mnuOther-0.43915503282684*m1rOther[6]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[6]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[6]*mnuOther+0.21957751641342*uOther[0]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[4]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[3]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(53,39) = 0.3833333333333334*m0rOther[3]*uOther[9]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*m0rOther[0]*uOther[7]*mnuOther-0.43915503282684*m1rOther[7]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[7]*mnuOther+0.21957751641342*uOther[0]*m0rOther[7]*mnuOther+0.175662013130736*m0rOther[3]*uOther[6]*mnuOther+0.175662013130736*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[5]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[3]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(54,30) = 0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[8]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(54,31) = 0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[8]*mnuOther-0.4391550328268398*m1rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(54,32) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(54,33) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(54,34) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[8]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[8]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[2]*uOther[6]*mnuOther+0.159719141249985*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(54,35) = 0.2195775164134199*m0rOther[6]*uOther[9]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[8]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(54,36) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[8]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(54,37) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.21957751641342*m0rOther[5]*uOther[8]*mnuOther+0.21957751641342*uOther[5]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(54,38) = 0.4621212121212121*m0rOther[4]*uOther[8]*mnuOther+0.149071198499986*m0rOther[0]*uOther[8]*mnuOther-0.2981423969999719*m1rOther[8]*mnuOther+0.4621212121212121*uOther[4]*m0rOther[8]*mnuOther+0.149071198499986*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[7]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[6]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[4]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[4]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[3]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[1]*mnuOther-0.4391550328268398*m1rOther[1]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(54,39) = 0.25*m0rOther[4]*uOther[9]*mnuOther+0.25*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[7]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[7]*mnuOther+0.21957751641342*m0rOther[5]*uOther[6]*mnuOther+0.21957751641342*uOther[5]*m0rOther[6]*mnuOther; 
  data->AEM_S(55,30) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[9]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(55,31) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(55,32) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[9]*mnuOther-0.4391550328268398*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[9]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(55,33) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(55,34) = 0.2195775164134199*m0rOther[6]*uOther[9]*mnuOther+0.2195775164134199*uOther[6]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[8]*mnuOther+0.2195775164134199*uOther[7]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(55,35) = 0.4621212121212121*m0rOther[9]*uOther[9]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[9]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*uOther[1]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(55,36) = 0.21957751641342*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(55,37) = 0.3273268353539885*m0rOther[3]*uOther[9]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(55,38) = 0.25*m0rOther[5]*uOther[8]*mnuOther+0.25*uOther[5]*m0rOther[8]*mnuOther+0.21957751641342*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*uOther[4]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[6]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[6]*mnuOther; 
  data->AEM_S(55,39) = 0.4621212121212121*m0rOther[5]*uOther[9]*mnuOther+0.149071198499986*m0rOther[0]*uOther[9]*mnuOther-0.2981423969999719*m1rOther[9]*mnuOther+0.4621212121212121*uOther[5]*m0rOther[9]*mnuOther+0.149071198499986*uOther[0]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[7]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[6]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[5]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[5]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[3]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[2]*mnuOther-0.4391550328268398*m1rOther[2]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(56,30) = 0.2195775164134199*m0rOther[3]*uOther[8]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[8]*mnuOther+0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(56,31) = 0.1963961012123931*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[8]*mnuOther+0.21957751641342*m0rOther[2]*uOther[8]*mnuOther+0.3273268353539885*uOther[6]*m0rOther[8]*mnuOther+0.21957751641342*uOther[2]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(56,32) = 0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.149071198499986*m0rOther[8]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.21957751641342*m0rOther[1]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.21957751641342*uOther[1]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(56,33) = 0.175662013130736*m0rOther[3]*uOther[9]*mnuOther+0.175662013130736*uOther[3]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[8]*mnuOther+0.21957751641342*m0rOther[0]*uOther[8]*mnuOther-0.43915503282684*m1rOther[8]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[8]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[8]*mnuOther+0.21957751641342*uOther[0]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(56,34) = 0.2195775164134199*m0rOther[5]*uOther[9]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[8]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(56,35) = 0.21957751641342*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[8]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[8]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(56,36) = 0.3833333333333334*m0rOther[9]*uOther[9]*mnuOther+0.1254728665219542*m0rOther[6]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[9]*mnuOther+0.1254728665219542*uOther[6]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[9]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[8]*mnuOther+0.29277002188456*m0rOther[7]*uOther[8]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[8]*mnuOther+0.29277002188456*uOther[7]*m0rOther[8]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[8]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[7]*mnuOther+0.351382110749967*uOther[1]*m0rOther[7]*mnuOther+0.9642857142857143*m0rOther[6]*uOther[6]*mnuOther+0.2874944542499729*m0rOther[2]*uOther[6]*mnuOther+0.2874944542499729*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(56,37) = 0.1928571428571429*m0rOther[8]*uOther[9]*mnuOther+0.29277002188456*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[9]*mnuOther+0.1928571428571429*uOther[8]*m0rOther[9]*mnuOther+0.29277002188456*uOther[7]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[9]*mnuOther+0.29277002188456*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[8]*mnuOther+0.29277002188456*uOther[6]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[8]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(56,38) = 0.1928571428571429*m0rOther[7]*uOther[9]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[9]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[8]*mnuOther+0.149071198499986*m0rOther[2]*uOther[8]*mnuOther+0.4621212121212121*uOther[6]*m0rOther[8]*mnuOther+0.149071198499986*uOther[2]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[6]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[6]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[5]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[5]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[4]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(56,39) = 0.3833333333333334*m0rOther[6]*uOther[9]*mnuOther+0.3833333333333334*uOther[6]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[8]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.1254728665219543*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(57,30) = 0.2195775164134199*m0rOther[3]*uOther[9]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(57,31) = 0.149071198499986*m0rOther[9]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[9]*mnuOther+0.21957751641342*m0rOther[2]*uOther[9]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[9]*mnuOther+0.21957751641342*uOther[2]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*uOther[7]*m0rOther[8]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(57,32) = 0.3273268353539885*m0rOther[7]*uOther[9]*mnuOther+0.21957751641342*m0rOther[1]*uOther[9]*mnuOther+0.3273268353539885*uOther[7]*m0rOther[9]*mnuOther+0.21957751641342*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[8]*mnuOther+0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(57,33) = 0.3273268353539885*m0rOther[5]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[9]*mnuOther+0.21957751641342*m0rOther[0]*uOther[9]*mnuOther-0.43915503282684*m1rOther[9]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[9]*mnuOther+0.21957751641342*uOther[0]*m0rOther[9]*mnuOther+0.175662013130736*m0rOther[3]*uOther[8]*mnuOther+0.175662013130736*uOther[3]*m0rOther[8]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(57,34) = 0.1963961012123931*m0rOther[3]*uOther[9]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[9]*mnuOther+0.21957751641342*m0rOther[5]*uOther[8]*mnuOther+0.21957751641342*uOther[5]*m0rOther[8]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(57,35) = 0.3273268353539885*m0rOther[3]*uOther[9]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[8]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[8]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(57,36) = 0.1928571428571429*m0rOther[8]*uOther[9]*mnuOther+0.29277002188456*m0rOther[7]*uOther[9]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[9]*mnuOther+0.1928571428571429*uOther[8]*m0rOther[9]*mnuOther+0.29277002188456*uOther[7]*m0rOther[9]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[9]*mnuOther+0.29277002188456*m0rOther[6]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[8]*mnuOther+0.29277002188456*uOther[6]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[8]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(57,37) = 0.4621212121212121*m0rOther[9]*uOther[9]*mnuOther+0.29277002188456*m0rOther[6]*uOther[9]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[9]*mnuOther+0.29277002188456*uOther[6]*m0rOther[9]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[8]*mnuOther+0.1254728665219542*m0rOther[7]*uOther[8]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[8]*mnuOther+0.1254728665219542*uOther[7]*m0rOther[8]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[8]*mnuOther+0.9642857142857143*m0rOther[7]*uOther[7]*mnuOther+0.2874944542499729*m0rOther[1]*uOther[7]*mnuOther+0.2874944542499729*uOther[1]*m0rOther[7]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[6]*mnuOther+0.351382110749967*m0rOther[2]*uOther[6]*mnuOther+0.351382110749967*uOther[2]*m0rOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(57,38) = 0.1928571428571429*m0rOther[6]*uOther[9]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[8]*mnuOther+0.3833333333333334*uOther[7]*m0rOther[8]*mnuOther+0.1254728665219543*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.29277002188456*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(57,39) = 0.4621212121212121*m0rOther[7]*uOther[9]*mnuOther+0.149071198499986*m0rOther[1]*uOther[9]*mnuOther+0.4621212121212121*uOther[7]*m0rOther[9]*mnuOther+0.149071198499986*uOther[1]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[8]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[6]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[5]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[4]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(58,30) = 0.149071198499986*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.149071198499986*uOther[4]*m0rOther[8]*mnuOther+0.25*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[6]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[4]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[4]*mnuOther; 
  data->AEM_S(58,31) = 0.2499586742703185*m0rOther[8]*uOther[8]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[8]*mnuOther+0.3833333333333334*uOther[1]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[6]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[4]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[4]*mnuOther-0.4391550328268398*m1rOther[4]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(58,32) = 0.149071198499986*m0rOther[6]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.149071198499986*uOther[6]*m0rOther[8]*mnuOther+0.25*uOther[2]*m0rOther[8]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[6]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[4]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(58,33) = 0.3833333333333334*m0rOther[3]*uOther[8]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[8]*mnuOther+0.175662013130736*m0rOther[3]*uOther[7]*mnuOther+0.175662013130736*uOther[3]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[6]*mnuOther+0.21957751641342*m0rOther[0]*uOther[6]*mnuOther-0.43915503282684*m1rOther[6]*mnuOther+0.1963961012123931*uOther[5]*m0rOther[6]*mnuOther+0.3273268353539885*uOther[4]*m0rOther[6]*mnuOther+0.21957751641342*uOther[0]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[4]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[4]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[3]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(58,34) = 0.4621212121212121*m0rOther[4]*uOther[8]*mnuOther+0.149071198499986*m0rOther[0]*uOther[8]*mnuOther-0.2981423969999719*m1rOther[8]*mnuOther+0.4621212121212121*uOther[4]*m0rOther[8]*mnuOther+0.149071198499986*uOther[0]*m0rOther[8]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[7]*mnuOther+0.2195775164134199*uOther[5]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[6]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[4]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[4]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[3]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[1]*mnuOther-0.4391550328268398*m1rOther[1]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(58,35) = 0.25*m0rOther[5]*uOther[8]*mnuOther+0.25*uOther[5]*m0rOther[8]*mnuOther+0.21957751641342*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*uOther[4]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[6]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[6]*mnuOther; 
  data->AEM_S(58,36) = 0.1928571428571429*m0rOther[7]*uOther[9]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[9]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[8]*mnuOther+0.149071198499986*m0rOther[2]*uOther[8]*mnuOther+0.4621212121212121*uOther[6]*m0rOther[8]*mnuOther+0.149071198499986*uOther[2]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[7]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[6]*mnuOther+0.3273268353539885*uOther[1]*m0rOther[6]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[5]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[5]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[4]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(58,37) = 0.1928571428571429*m0rOther[6]*uOther[9]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[9]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[8]*mnuOther+0.3833333333333334*uOther[7]*m0rOther[8]*mnuOther+0.1254728665219543*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.29277002188456*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(58,38) = 0.25*m0rOther[9]*uOther[9]*mnuOther+0.5898601398601399*m0rOther[8]*uOther[8]*mnuOther+0.2499586742703185*m0rOther[1]*uOther[8]*mnuOther+0.2499586742703185*uOther[1]*m0rOther[8]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[7]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[6]*mnuOther+0.149071198499986*m0rOther[2]*uOther[6]*mnuOther+0.149071198499986*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.4621212121212121*m0rOther[4]*uOther[4]*mnuOther+0.149071198499986*m0rOther[0]*uOther[4]*mnuOther-0.2981423969999719*m1rOther[4]*mnuOther+0.149071198499986*uOther[0]*m0rOther[4]*mnuOther+0.3833333333333334*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(58,39) = 0.25*m0rOther[8]*uOther[9]*mnuOther+0.25*uOther[8]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[7]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[7]*mnuOther; 
  data->AEM_S(59,30) = 0.149071198499986*m0rOther[5]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.149071198499986*uOther[5]*m0rOther[9]*mnuOther+0.25*uOther[0]*m0rOther[9]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[7]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[5]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[5]*mnuOther; 
  data->AEM_S(59,31) = 0.149071198499986*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.149071198499986*uOther[7]*m0rOther[9]*mnuOther+0.25*uOther[1]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[7]*mnuOther+0.1963961012123931*uOther[6]*m0rOther[7]*mnuOther+0.2195775164134199*uOther[2]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[5]*mnuOther+0.2195775164134199*uOther[3]*m0rOther[5]*mnuOther; 
  data->AEM_S(59,32) = 0.2499586742703185*m0rOther[9]*uOther[9]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[9]*mnuOther+0.3833333333333334*uOther[2]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[7]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[7]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[6]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[5]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[5]*mnuOther-0.4391550328268398*m1rOther[5]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[3]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(59,33) = 0.3833333333333334*m0rOther[3]*uOther[9]*mnuOther+0.3833333333333334*uOther[3]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[7]*mnuOther+0.21957751641342*m0rOther[0]*uOther[7]*mnuOther-0.43915503282684*m1rOther[7]*mnuOther+0.3273268353539885*uOther[5]*m0rOther[7]*mnuOther+0.1963961012123931*uOther[4]*m0rOther[7]*mnuOther+0.21957751641342*uOther[0]*m0rOther[7]*mnuOther+0.175662013130736*m0rOther[3]*uOther[6]*mnuOther+0.175662013130736*uOther[3]*m0rOther[6]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[5]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[3]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(59,34) = 0.25*m0rOther[4]*uOther[9]*mnuOther+0.25*uOther[4]*m0rOther[9]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[7]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[7]*mnuOther+0.21957751641342*m0rOther[5]*uOther[6]*mnuOther+0.21957751641342*uOther[5]*m0rOther[6]*mnuOther; 
  data->AEM_S(59,35) = 0.4621212121212121*m0rOther[5]*uOther[9]*mnuOther+0.149071198499986*m0rOther[0]*uOther[9]*mnuOther-0.2981423969999719*m1rOther[9]*mnuOther+0.4621212121212121*uOther[5]*m0rOther[9]*mnuOther+0.149071198499986*uOther[0]*m0rOther[9]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[7]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[7]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[6]*mnuOther+0.2195775164134199*uOther[4]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[5]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[5]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[3]*mnuOther+0.2195775164134199*uOther[1]*m0rOther[3]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[2]*mnuOther-0.4391550328268398*m1rOther[2]*mnuOther+0.2195775164134199*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(59,36) = 0.3833333333333334*m0rOther[6]*uOther[9]*mnuOther+0.3833333333333334*uOther[6]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[8]*mnuOther+0.1928571428571429*uOther[7]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[7]*uOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[7]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[7]*mnuOther+0.1254728665219543*m0rOther[6]*uOther[6]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[6]*mnuOther+0.1963961012123931*uOther[2]*m0rOther[6]*mnuOther+0.21957751641342*m0rOther[4]*uOther[5]*mnuOther+0.21957751641342*uOther[4]*m0rOther[5]*mnuOther+0.175662013130736*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(59,37) = 0.4621212121212121*m0rOther[7]*uOther[9]*mnuOther+0.149071198499986*m0rOther[1]*uOther[9]*mnuOther+0.4621212121212121*uOther[7]*m0rOther[9]*mnuOther+0.149071198499986*uOther[1]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[8]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[8]*mnuOther+0.29277002188456*m0rOther[6]*uOther[7]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[7]*mnuOther+0.29277002188456*uOther[6]*m0rOther[7]*mnuOther+0.3273268353539885*uOther[2]*m0rOther[7]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[6]*mnuOther+0.1963961012123931*uOther[1]*m0rOther[6]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[5]*mnuOther+0.3273268353539885*uOther[3]*m0rOther[5]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[4]*mnuOther+0.1963961012123931*uOther[3]*m0rOther[4]*mnuOther+0.21957751641342*m0rOther[0]*uOther[3]*mnuOther-0.43915503282684*m1rOther[3]*mnuOther+0.21957751641342*uOther[0]*m0rOther[3]*mnuOther+0.21957751641342*m0rOther[1]*uOther[2]*mnuOther+0.21957751641342*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(59,38) = 0.25*m0rOther[8]*uOther[9]*mnuOther+0.25*uOther[8]*m0rOther[9]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[7]*mnuOther+0.1928571428571429*uOther[6]*m0rOther[7]*mnuOther; 
  data->AEM_S(59,39) = 0.5898601398601399*m0rOther[9]*uOther[9]*mnuOther+0.2499586742703185*m0rOther[2]*uOther[9]*mnuOther+0.2499586742703185*uOther[2]*m0rOther[9]*mnuOther+0.25*m0rOther[8]*uOther[8]*mnuOther+0.4621212121212121*m0rOther[7]*uOther[7]*mnuOther+0.149071198499986*m0rOther[1]*uOther[7]*mnuOther+0.149071198499986*uOther[1]*m0rOther[7]*mnuOther+0.3833333333333334*m0rOther[6]*uOther[6]*mnuOther+0.4621212121212121*m0rOther[5]*uOther[5]*mnuOther+0.149071198499986*m0rOther[0]*uOther[5]*mnuOther-0.2981423969999719*m1rOther[5]*mnuOther+0.149071198499986*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3833333333333334*m0rOther[3]*uOther[3]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
 
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
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(40,10) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(40,11) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(40,12) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(40,13) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(40,14) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(40,15) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(40,16) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(40,17) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(40,18) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(40,19) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(41,10) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(41,11) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(41,12) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(41,13) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(41,14) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(41,15) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(41,16) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(41,17) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(41,18) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(42,10) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(42,11) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(42,12) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(42,13) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(42,14) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(42,15) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(42,16) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(42,17) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(42,19) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(43,10) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,11) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(43,12) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(43,13) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(43,14) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,15) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,16) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(43,17) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(43,18) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(43,19) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(44,10) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(44,11) = 0.4391550328268398*m0rSelf[8]*mnuSelf+0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(44,12) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(44,13) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(44,14) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(44,16) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(44,17) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(44,18) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(45,10) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(45,11) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(45,12) = 0.4391550328268398*m0rSelf[9]*mnuSelf+0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(45,13) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(45,15) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(45,16) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(45,17) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(45,19) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(46,10) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(46,11) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(46,12) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(46,13) = 0.4391550328268399*m0rSelf[8]*mnuSelf+0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(46,14) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(46,15) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(46,16) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(46,17) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(46,18) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(47,10) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(47,11) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(47,12) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(47,13) = 0.4391550328268399*m0rSelf[9]*mnuSelf+0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(47,14) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(47,15) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(47,16) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(47,17) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(47,19) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(48,10) = 0.5*m0rSelf[8]*mnuSelf; 
  data->AEM_S(48,11) = 0.4391550328268398*m0rSelf[4]*mnuSelf; 
  data->AEM_S(48,13) = 0.4391550328268399*m0rSelf[6]*mnuSelf; 
  data->AEM_S(48,14) = 0.2981423969999719*m0rSelf[8]*mnuSelf+0.4391550328268398*m0rSelf[1]*mnuSelf; 
  data->AEM_S(48,16) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(48,18) = 0.2981423969999719*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(49,10) = 0.5*m0rSelf[9]*mnuSelf; 
  data->AEM_S(49,12) = 0.4391550328268398*m0rSelf[5]*mnuSelf; 
  data->AEM_S(49,13) = 0.4391550328268399*m0rSelf[7]*mnuSelf; 
  data->AEM_S(49,15) = 0.2981423969999719*m0rSelf[9]*mnuSelf+0.4391550328268398*m0rSelf[2]*mnuSelf; 
  data->AEM_S(49,17) = 0.4391550328268399*m0rSelf[3]*mnuSelf; 
  data->AEM_S(49,19) = 0.2981423969999719*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(40,20) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(40,21) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(40,22) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(40,23) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(40,24) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(40,25) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(40,26) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(40,27) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(40,28) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(40,29) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(41,20) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(41,21) = (-0.4472135954999579*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(41,22) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(41,23) = (-0.447213595499958*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(41,24) = (-0.4391550328268398*cMSelf[18]*mnuSelf)-0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(41,25) = -0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(41,26) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(41,27) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(41,28) = -0.4391550328268398*cMSelf[14]*mnuSelf; 
  data->AEM_S(42,20) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(42,21) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(42,22) = (-0.4472135954999579*cMSelf[15]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(42,23) = (-0.447213595499958*cMSelf[17]*mnuSelf)-0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(42,24) = -0.5000000000000001*cMSelf[16]*mnuSelf; 
  data->AEM_S(42,25) = (-0.4391550328268398*cMSelf[19]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf; 
  data->AEM_S(42,26) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(42,27) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(42,29) = -0.4391550328268398*cMSelf[15]*mnuSelf; 
  data->AEM_S(43,20) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(43,21) = (-0.447213595499958*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(43,22) = (-0.447213595499958*cMSelf[17]*mnuSelf)-0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(43,23) = (-0.4472135954999579*cMSelf[15]*mnuSelf)-0.4472135954999579*cMSelf[14]*mnuSelf-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(43,24) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(43,25) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(43,26) = (-0.4391550328268399*cMSelf[18]*mnuSelf)-0.4*cMSelf[17]*mnuSelf-0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(43,27) = (-0.4391550328268399*cMSelf[19]*mnuSelf)-0.4*cMSelf[16]*mnuSelf-0.447213595499958*cMSelf[12]*mnuSelf; 
  data->AEM_S(43,28) = -0.4391550328268399*cMSelf[16]*mnuSelf; 
  data->AEM_S(43,29) = -0.4391550328268399*cMSelf[17]*mnuSelf; 
  data->AEM_S(44,20) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(44,21) = (-0.4391550328268398*cMSelf[18]*mnuSelf)-0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(44,22) = -0.5000000000000001*cMSelf[16]*mnuSelf; 
  data->AEM_S(44,23) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(44,24) = (-0.31943828249997*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(44,26) = (-0.31943828249997*cMSelf[16]*mnuSelf)-0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(44,27) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(44,28) = (-0.2981423969999719*cMSelf[18]*mnuSelf)-0.4391550328268398*cMSelf[11]*mnuSelf; 
  data->AEM_S(45,20) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(45,21) = -0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(45,22) = (-0.4391550328268398*cMSelf[19]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf; 
  data->AEM_S(45,23) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(45,25) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(45,26) = -0.4472135954999579*cMSelf[16]*mnuSelf; 
  data->AEM_S(45,27) = (-0.31943828249997*cMSelf[17]*mnuSelf)-0.5000000000000001*cMSelf[11]*mnuSelf; 
  data->AEM_S(45,29) = (-0.2981423969999719*cMSelf[19]*mnuSelf)-0.4391550328268398*cMSelf[12]*mnuSelf; 
  data->AEM_S(46,20) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(46,21) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(46,22) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(46,23) = (-0.4391550328268399*cMSelf[18]*mnuSelf)-0.4*cMSelf[17]*mnuSelf-0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(46,24) = (-0.31943828249997*cMSelf[16]*mnuSelf)-0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(46,25) = -0.4472135954999579*cMSelf[16]*mnuSelf; 
  data->AEM_S(46,26) = (-0.4472135954999579*cMSelf[15]*mnuSelf)-0.31943828249997*cMSelf[14]*mnuSelf-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(46,27) = -0.4*cMSelf[13]*mnuSelf; 
  data->AEM_S(46,28) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(47,20) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(47,21) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(47,22) = -0.447213595499958*cMSelf[13]*mnuSelf; 
  data->AEM_S(47,23) = (-0.4391550328268399*cMSelf[19]*mnuSelf)-0.4*cMSelf[16]*mnuSelf-0.447213595499958*cMSelf[12]*mnuSelf; 
  data->AEM_S(47,24) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(47,25) = (-0.31943828249997*cMSelf[17]*mnuSelf)-0.5000000000000001*cMSelf[11]*mnuSelf; 
  data->AEM_S(47,26) = -0.4*cMSelf[13]*mnuSelf; 
  data->AEM_S(47,27) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.4472135954999579*cMSelf[14]*mnuSelf-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(47,29) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(48,20) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(48,21) = -0.4391550328268398*cMSelf[14]*mnuSelf; 
  data->AEM_S(48,23) = -0.4391550328268399*cMSelf[16]*mnuSelf; 
  data->AEM_S(48,24) = (-0.2981423969999719*cMSelf[18]*mnuSelf)-0.4391550328268398*cMSelf[11]*mnuSelf; 
  data->AEM_S(48,26) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(48,28) = (-0.2981423969999719*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(49,20) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(49,22) = -0.4391550328268398*cMSelf[15]*mnuSelf; 
  data->AEM_S(49,23) = -0.4391550328268399*cMSelf[17]*mnuSelf; 
  data->AEM_S(49,25) = (-0.2981423969999719*cMSelf[19]*mnuSelf)-0.4391550328268398*cMSelf[12]*mnuSelf; 
  data->AEM_S(49,27) = -0.4391550328268399*cMSelf[13]*mnuSelf; 
  data->AEM_S(49,29) = (-0.2981423969999719*cMSelf[15]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(40,40) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(40,41) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(40,42) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(40,43) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(40,44) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(40,45) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(40,46) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(40,47) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(40,48) = -0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(40,49) = -0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(41,40) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(41,41) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(41,42) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(41,43) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(41,44) = (-0.4391550328268398*m0rOther[8]*mnuOther)-0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(41,45) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(41,46) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(41,47) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(41,48) = -0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(42,40) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(42,41) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(42,42) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(42,43) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(42,44) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(42,45) = (-0.4391550328268398*m0rOther[9]*mnuOther)-0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(42,46) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(42,47) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(42,49) = -0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(43,40) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(43,41) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(43,42) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(43,43) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(43,44) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(43,45) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(43,46) = (-0.4391550328268399*m0rOther[8]*mnuOther)-0.4*m0rOther[7]*mnuOther-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(43,47) = (-0.4391550328268399*m0rOther[9]*mnuOther)-0.4*m0rOther[6]*mnuOther-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(43,48) = -0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(43,49) = -0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(44,40) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(44,41) = (-0.4391550328268398*m0rOther[8]*mnuOther)-0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(44,42) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(44,43) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(44,44) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(44,46) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(44,47) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(44,48) = (-0.2981423969999719*m0rOther[8]*mnuOther)-0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(45,40) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(45,41) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(45,42) = (-0.4391550328268398*m0rOther[9]*mnuOther)-0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(45,43) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(45,45) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(45,46) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(45,47) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(45,49) = (-0.2981423969999719*m0rOther[9]*mnuOther)-0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(46,40) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(46,41) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(46,42) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(46,43) = (-0.4391550328268399*m0rOther[8]*mnuOther)-0.4*m0rOther[7]*mnuOther-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(46,44) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(46,45) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(46,46) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.31943828249997*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(46,47) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(46,48) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(47,40) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(47,41) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(47,42) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(47,43) = (-0.4391550328268399*m0rOther[9]*mnuOther)-0.4*m0rOther[6]*mnuOther-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(47,44) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(47,45) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(47,46) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(47,47) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(47,49) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(48,40) = -0.5*m0rOther[8]*mnuOther; 
  data->AEM_S(48,41) = -0.4391550328268398*m0rOther[4]*mnuOther; 
  data->AEM_S(48,43) = -0.4391550328268399*m0rOther[6]*mnuOther; 
  data->AEM_S(48,44) = (-0.2981423969999719*m0rOther[8]*mnuOther)-0.4391550328268398*m0rOther[1]*mnuOther; 
  data->AEM_S(48,46) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(48,48) = (-0.2981423969999719*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(49,40) = -0.5*m0rOther[9]*mnuOther; 
  data->AEM_S(49,42) = -0.4391550328268398*m0rOther[5]*mnuOther; 
  data->AEM_S(49,43) = -0.4391550328268399*m0rOther[7]*mnuOther; 
  data->AEM_S(49,45) = (-0.2981423969999719*m0rOther[9]*mnuOther)-0.4391550328268398*m0rOther[2]*mnuOther; 
  data->AEM_S(49,47) = -0.4391550328268399*m0rOther[3]*mnuOther; 
  data->AEM_S(49,49) = (-0.2981423969999719*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(40,50) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(40,51) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(40,52) = 0.5*cMOther[12]*mnuOther; 
  data->AEM_S(40,53) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(40,54) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(40,55) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(40,56) = 0.5*cMOther[16]*mnuOther; 
  data->AEM_S(40,57) = 0.5*cMOther[17]*mnuOther; 
  data->AEM_S(40,58) = 0.5*cMOther[18]*mnuOther; 
  data->AEM_S(40,59) = 0.5*cMOther[19]*mnuOther; 
  data->AEM_S(41,50) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(41,51) = 0.4472135954999579*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(41,52) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(41,53) = 0.447213595499958*cMOther[16]*mnuOther+0.5*cMOther[12]*mnuOther; 
  data->AEM_S(41,54) = 0.4391550328268398*cMOther[18]*mnuOther+0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(41,55) = 0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(41,56) = 0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(41,57) = 0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(41,58) = 0.4391550328268398*cMOther[14]*mnuOther; 
  data->AEM_S(42,50) = 0.5*cMOther[12]*mnuOther; 
  data->AEM_S(42,51) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(42,52) = 0.4472135954999579*cMOther[15]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(42,53) = 0.447213595499958*cMOther[17]*mnuOther+0.5*cMOther[11]*mnuOther; 
  data->AEM_S(42,54) = 0.5000000000000001*cMOther[16]*mnuOther; 
  data->AEM_S(42,55) = 0.4391550328268398*cMOther[19]*mnuOther+0.4472135954999579*cMOther[12]*mnuOther; 
  data->AEM_S(42,56) = 0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(42,57) = 0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(42,59) = 0.4391550328268398*cMOther[15]*mnuOther; 
  data->AEM_S(43,50) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(43,51) = 0.447213595499958*cMOther[16]*mnuOther+0.5*cMOther[12]*mnuOther; 
  data->AEM_S(43,52) = 0.447213595499958*cMOther[17]*mnuOther+0.5*cMOther[11]*mnuOther; 
  data->AEM_S(43,53) = 0.4472135954999579*cMOther[15]*mnuOther+0.4472135954999579*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(43,54) = 0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(43,55) = 0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(43,56) = 0.4391550328268399*cMOther[18]*mnuOther+0.4*cMOther[17]*mnuOther+0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(43,57) = 0.4391550328268399*cMOther[19]*mnuOther+0.4*cMOther[16]*mnuOther+0.447213595499958*cMOther[12]*mnuOther; 
  data->AEM_S(43,58) = 0.4391550328268399*cMOther[16]*mnuOther; 
  data->AEM_S(43,59) = 0.4391550328268399*cMOther[17]*mnuOther; 
  data->AEM_S(44,50) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(44,51) = 0.4391550328268398*cMOther[18]*mnuOther+0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(44,52) = 0.5000000000000001*cMOther[16]*mnuOther; 
  data->AEM_S(44,53) = 0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(44,54) = 0.31943828249997*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(44,56) = 0.31943828249997*cMOther[16]*mnuOther+0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(44,57) = 0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(44,58) = 0.2981423969999719*cMOther[18]*mnuOther+0.4391550328268398*cMOther[11]*mnuOther; 
  data->AEM_S(45,50) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(45,51) = 0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(45,52) = 0.4391550328268398*cMOther[19]*mnuOther+0.4472135954999579*cMOther[12]*mnuOther; 
  data->AEM_S(45,53) = 0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(45,55) = 0.31943828249997*cMOther[15]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(45,56) = 0.4472135954999579*cMOther[16]*mnuOther; 
  data->AEM_S(45,57) = 0.31943828249997*cMOther[17]*mnuOther+0.5000000000000001*cMOther[11]*mnuOther; 
  data->AEM_S(45,59) = 0.2981423969999719*cMOther[19]*mnuOther+0.4391550328268398*cMOther[12]*mnuOther; 
  data->AEM_S(46,50) = 0.5*cMOther[16]*mnuOther; 
  data->AEM_S(46,51) = 0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(46,52) = 0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(46,53) = 0.4391550328268399*cMOther[18]*mnuOther+0.4*cMOther[17]*mnuOther+0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(46,54) = 0.31943828249997*cMOther[16]*mnuOther+0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(46,55) = 0.4472135954999579*cMOther[16]*mnuOther; 
  data->AEM_S(46,56) = 0.4472135954999579*cMOther[15]*mnuOther+0.31943828249997*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(46,57) = 0.4*cMOther[13]*mnuOther; 
  data->AEM_S(46,58) = 0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(47,50) = 0.5*cMOther[17]*mnuOther; 
  data->AEM_S(47,51) = 0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(47,52) = 0.447213595499958*cMOther[13]*mnuOther; 
  data->AEM_S(47,53) = 0.4391550328268399*cMOther[19]*mnuOther+0.4*cMOther[16]*mnuOther+0.447213595499958*cMOther[12]*mnuOther; 
  data->AEM_S(47,54) = 0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(47,55) = 0.31943828249997*cMOther[17]*mnuOther+0.5000000000000001*cMOther[11]*mnuOther; 
  data->AEM_S(47,56) = 0.4*cMOther[13]*mnuOther; 
  data->AEM_S(47,57) = 0.31943828249997*cMOther[15]*mnuOther+0.4472135954999579*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(47,59) = 0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(48,50) = 0.5*cMOther[18]*mnuOther; 
  data->AEM_S(48,51) = 0.4391550328268398*cMOther[14]*mnuOther; 
  data->AEM_S(48,53) = 0.4391550328268399*cMOther[16]*mnuOther; 
  data->AEM_S(48,54) = 0.2981423969999719*cMOther[18]*mnuOther+0.4391550328268398*cMOther[11]*mnuOther; 
  data->AEM_S(48,56) = 0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(48,58) = 0.2981423969999719*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(49,50) = 0.5*cMOther[19]*mnuOther; 
  data->AEM_S(49,52) = 0.4391550328268398*cMOther[15]*mnuOther; 
  data->AEM_S(49,53) = 0.4391550328268399*cMOther[17]*mnuOther; 
  data->AEM_S(49,55) = 0.2981423969999719*cMOther[19]*mnuOther+0.4391550328268398*cMOther[12]*mnuOther; 
  data->AEM_S(49,57) = 0.4391550328268399*cMOther[13]*mnuOther; 
  data->AEM_S(49,59) = 0.2981423969999719*cMOther[15]*mnuOther+0.5*cMOther[10]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(50,10) = (-0.25*m0rSelf[9]*uSelf[19]*mnuSelf)-0.25*m0rSelf[8]*uSelf[18]*mnuSelf-0.25*m0rSelf[7]*uSelf[17]*mnuSelf-0.25*m0rSelf[6]*uSelf[16]*mnuSelf-0.25*m0rSelf[5]*uSelf[15]*mnuSelf-0.25*m0rSelf[4]*uSelf[14]*mnuSelf-0.25*m0rSelf[3]*uSelf[13]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf-0.25*m0rSelf[1]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(50,11) = (-0.2195775164134199*m0rSelf[4]*uSelf[18]*mnuSelf)-0.2500000000000001*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[13]*mnuSelf-0.25*m0rSelf[2]*uSelf[13]*mnuSelf-0.25*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,12) = (-0.2195775164134199*m0rSelf[5]*uSelf[19]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf-0.25*m0rSelf[3]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,13) = (-0.2195775164134199*m0rSelf[7]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[6]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[17]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[16]*mnuSelf-0.2*m0rSelf[7]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[11]*mnuSelf-0.25*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,14) = (-0.149071198499986*m0rSelf[8]*uSelf[18]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,15) = (-0.149071198499986*m0rSelf[9]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[2]*uSelf[19]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[11]*mnuSelf-0.25*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,16) = (-0.2195775164134199*m0rSelf[3]*uSelf[18]*mnuSelf)-0.2*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[13]*mnuSelf-0.2*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.25*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,17) = (-0.2195775164134199*m0rSelf[3]*uSelf[19]*mnuSelf)-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.2*m0rSelf[3]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[13]*mnuSelf-0.2*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[11]*mnuSelf-0.25*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,18) = (-0.149071198499986*m0rSelf[4]*uSelf[18]*mnuSelf)-0.25*m0rSelf[0]*uSelf[18]*mnuSelf+0.5*m1rSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[16]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[8]*uSelf[10]*mnuSelf; 
  data->AEM_S(50,19) = (-0.149071198499986*m0rSelf[5]*uSelf[19]*mnuSelf)-0.25*m0rSelf[0]*uSelf[19]*mnuSelf+0.5*m1rSelf[19]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[17]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[12]*mnuSelf-0.25*m0rSelf[9]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,10) = (-0.2195775164134199*m0rSelf[4]*uSelf[18]*mnuSelf)-0.2500000000000001*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[13]*mnuSelf-0.25*m0rSelf[2]*uSelf[13]*mnuSelf-0.25*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,11) = (-0.25*m0rSelf[9]*uSelf[19]*mnuSelf)-0.3833333333333334*m0rSelf[8]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[18]*mnuSelf-0.45*m0rSelf[7]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf-0.25*m0rSelf[5]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf-0.45*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[12]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[11]*mnuSelf-0.45*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(51,12) = (-0.2195775164134199*m0rSelf[7]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[6]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[17]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[16]*mnuSelf-0.2*m0rSelf[7]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[11]*mnuSelf-0.25*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,13) = (-0.2195775164134199*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[3]*uSelf[18]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[17]*mnuSelf-0.2*m0rSelf[5]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[16]*mnuSelf+0.447213595499958*m1rSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[15]*mnuSelf-0.2*m0rSelf[6]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[13]*mnuSelf-0.45*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf-0.45*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,14) = (-0.3273268353539885*m0rSelf[4]*uSelf[18]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[18]*mnuSelf+0.4391550328268398*m1rSelf[18]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[17]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,15) = (-0.2195775164134199*m0rSelf[3]*uSelf[19]*mnuSelf)-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[17]*mnuSelf+0.5000000000000001*m1rSelf[17]*mnuSelf-0.2*m0rSelf[3]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[15]*mnuSelf-0.25*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[13]*mnuSelf-0.2*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.25*m0rSelf[5]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,16) = (-0.1963961012123931*m0rSelf[7]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[6]*uSelf[18]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[17]*mnuSelf-0.2*m0rSelf[2]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[16]*mnuSelf-0.2*m0rSelf[3]*uSelf[15]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[14]*mnuSelf-0.2*m0rSelf[5]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.447213595499958*m1rSelf[13]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,17) = (-0.149071198499986*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[6]*uSelf[19]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[17]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[17]*mnuSelf-0.45*m0rSelf[1]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[16]*mnuSelf-0.2*m0rSelf[2]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[15]*mnuSelf+0.5000000000000001*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[13]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[12]*mnuSelf-0.2*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.45*m0rSelf[7]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,18) = (-0.2499586742703185*m0rSelf[8]*uSelf[18]*mnuSelf)-0.3833333333333334*m0rSelf[1]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[14]*mnuSelf+0.4391550328268398*m1rSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[12]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[11]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(51,19) = (-0.149071198499986*m0rSelf[7]*uSelf[19]*mnuSelf)-0.25*m0rSelf[1]*uSelf[19]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[12]*mnuSelf-0.25*m0rSelf[9]*uSelf[11]*mnuSelf; 
  data->AEM_S(52,10) = (-0.2195775164134199*m0rSelf[5]*uSelf[19]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf-0.25*m0rSelf[3]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(52,11) = (-0.2195775164134199*m0rSelf[7]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[6]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[17]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[16]*mnuSelf-0.2*m0rSelf[7]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[11]*mnuSelf-0.25*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(52,12) = (-0.3833333333333334*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[2]*uSelf[19]*mnuSelf-0.25*m0rSelf[8]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf-0.45*m0rSelf[6]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.25*m0rSelf[4]*uSelf[14]*mnuSelf-0.45*m0rSelf[3]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[12]*mnuSelf-0.45*m0rSelf[2]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(52,13) = (-0.1963961012123931*m0rSelf[3]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[4]*uSelf[18]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[17]*mnuSelf-0.2*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.447213595499958*m1rSelf[17]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[14]*mnuSelf-0.2*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[13]*mnuSelf-0.45*m0rSelf[2]*uSelf[13]*mnuSelf-0.45*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(52,14) = (-0.2195775164134199*m0rSelf[3]*uSelf[18]*mnuSelf)-0.2*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[16]*mnuSelf+0.5000000000000001*m1rSelf[16]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[14]*mnuSelf-0.25*m0rSelf[2]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[13]*mnuSelf-0.2*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.25*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(52,15) = (-0.3273268353539885*m0rSelf[5]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[19]*mnuSelf+0.4391550328268398*m1rSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[12]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(52,16) = (-0.1963961012123931*m0rSelf[6]*uSelf[19]*mnuSelf)-0.149071198499986*m0rSelf[8]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[18]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[17]*mnuSelf-0.2*m0rSelf[1]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[16]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[16]*mnuSelf-0.45*m0rSelf[2]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[14]*mnuSelf+0.5000000000000001*m1rSelf[14]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[13]*mnuSelf-0.45*m0rSelf[6]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[11]*mnuSelf-0.2*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(52,17) = (-0.3273268353539885*m0rSelf[7]*uSelf[19]*mnuSelf)-0.21957751641342*m0rSelf[1]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[16]*mnuSelf-0.2*m0rSelf[1]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[13]*mnuSelf-0.2*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.447213595499958*m1rSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[11]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(52,18) = (-0.149071198499986*m0rSelf[6]*uSelf[18]*mnuSelf)-0.25*m0rSelf[2]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[17]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[8]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(52,19) = (-0.2499586742703185*m0rSelf[9]*uSelf[19]*mnuSelf)-0.3833333333333334*m0rSelf[2]*uSelf[19]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[15]*mnuSelf+0.4391550328268398*m1rSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[13]*mnuSelf-0.3833333333333334*m0rSelf[9]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,10) = (-0.2195775164134199*m0rSelf[7]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[6]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[17]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[16]*mnuSelf-0.2*m0rSelf[7]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[11]*mnuSelf-0.25*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,11) = (-0.2195775164134199*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[3]*uSelf[18]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[17]*mnuSelf-0.2*m0rSelf[5]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[16]*mnuSelf+0.447213595499958*m1rSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[15]*mnuSelf-0.2*m0rSelf[6]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[13]*mnuSelf-0.45*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf-0.45*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,12) = (-0.1963961012123931*m0rSelf[3]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[4]*uSelf[18]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[17]*mnuSelf-0.2*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.447213595499958*m1rSelf[17]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[14]*mnuSelf-0.2*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[13]*mnuSelf-0.45*m0rSelf[2]*uSelf[13]*mnuSelf-0.45*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,13) = (-0.3833333333333334*m0rSelf[9]*uSelf[19]*mnuSelf)-0.175662013130736*m0rSelf[6]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[19]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[18]*mnuSelf-0.175662013130736*m0rSelf[7]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[18]*mnuSelf-0.175662013130736*m0rSelf[8]*uSelf[17]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[17]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[17]*mnuSelf-0.175662013130736*m0rSelf[9]*uSelf[16]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[16]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[15]*mnuSelf-0.2*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.2*m0rSelf[5]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf-0.81*m0rSelf[3]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[12]*mnuSelf-0.45*m0rSelf[2]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[11]*mnuSelf-0.45*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(53,14) = (-0.1963961012123931*m0rSelf[7]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[6]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[17]*mnuSelf-0.2*m0rSelf[2]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[16]*mnuSelf-0.2*m0rSelf[3]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[14]*mnuSelf-0.2*m0rSelf[5]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,15) = (-0.3273268353539885*m0rSelf[7]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[17]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[16]*mnuSelf-0.2*m0rSelf[1]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[13]*mnuSelf-0.2*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[11]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,16) = (-0.175662013130736*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[5]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[18]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[18]*mnuSelf+0.43915503282684*m1rSelf[18]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[17]*mnuSelf-0.2*m0rSelf[0]*uSelf[17]*mnuSelf+0.4*m1rSelf[17]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[15]*mnuSelf-0.2*m0rSelf[1]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[9]*uSelf[13]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[12]*mnuSelf-0.2*m0rSelf[5]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[10]*mnuSelf-0.2*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,17) = (-0.3273268353539885*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[4]*uSelf[19]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[19]*mnuSelf+0.43915503282684*m1rSelf[19]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[18]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[16]*mnuSelf-0.2*m0rSelf[0]*uSelf[16]*mnuSelf+0.4*m1rSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[15]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[14]*mnuSelf-0.2*m0rSelf[2]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[8]*uSelf[13]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[12]*mnuSelf-0.2*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[12]*mnuSelf+0.447213595499958*m1rSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[10]*mnuSelf-0.2*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,18) = (-0.3833333333333334*m0rSelf[3]*uSelf[18]*mnuSelf)-0.175662013130736*m0rSelf[3]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[16]*mnuSelf+0.43915503282684*m1rSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[14]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[13]*mnuSelf-0.175662013130736*m0rSelf[7]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(53,19) = (-0.3833333333333334*m0rSelf[3]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[5]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[17]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[17]*mnuSelf+0.43915503282684*m1rSelf[17]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[14]*mnuSelf-0.3833333333333334*m0rSelf[9]*uSelf[13]*mnuSelf-0.175662013130736*m0rSelf[6]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,10) = (-0.149071198499986*m0rSelf[8]*uSelf[18]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,11) = (-0.3273268353539885*m0rSelf[4]*uSelf[18]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[18]*mnuSelf+0.4391550328268398*m1rSelf[18]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[17]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,12) = (-0.2195775164134199*m0rSelf[3]*uSelf[18]*mnuSelf)-0.2*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[16]*mnuSelf+0.5000000000000001*m1rSelf[16]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[14]*mnuSelf-0.25*m0rSelf[2]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[13]*mnuSelf-0.2*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.25*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,13) = (-0.1963961012123931*m0rSelf[7]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[6]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[17]*mnuSelf-0.2*m0rSelf[2]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[16]*mnuSelf-0.2*m0rSelf[3]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[14]*mnuSelf-0.2*m0rSelf[5]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,14) = (-0.25*m0rSelf[9]*uSelf[19]*mnuSelf)-0.4621212121212121*m0rSelf[8]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[17]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[16]*mnuSelf-0.25*m0rSelf[5]*uSelf[15]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[14]*mnuSelf+0.31943828249997*m1rSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[12]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[11]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(54,15) = (-0.2195775164134199*m0rSelf[6]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[7]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[17]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[16]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf-0.25*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[5]*uSelf[14]*mnuSelf-0.2*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[11]*mnuSelf; 
  data->AEM_S(54,16) = (-0.2195775164134199*m0rSelf[5]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[3]*uSelf[18]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[17]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[16]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[16]*mnuSelf+0.31943828249997*m1rSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[12]*mnuSelf+0.5000000000000001*m1rSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[11]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,17) = (-0.1963961012123931*m0rSelf[3]*uSelf[19]*mnuSelf)-0.21957751641342*m0rSelf[5]*uSelf[18]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.4472135954999579*m1rSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[13]*mnuSelf-0.2*m0rSelf[2]*uSelf[13]*mnuSelf-0.2*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,18) = (-0.4621212121212121*m0rSelf[4]*uSelf[18]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[18]*mnuSelf+0.2981423969999719*m1rSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[15]*mnuSelf-0.4621212121212121*m0rSelf[8]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[12]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[11]*mnuSelf+0.4391550328268398*m1rSelf[11]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[10]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(54,19) = (-0.25*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[3]*uSelf[17]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[6]*uSelf[15]*mnuSelf-0.25*m0rSelf[9]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[13]*mnuSelf; 
  data->AEM_S(55,10) = (-0.149071198499986*m0rSelf[9]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[2]*uSelf[19]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[11]*mnuSelf-0.25*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(55,11) = (-0.2195775164134199*m0rSelf[3]*uSelf[19]*mnuSelf)-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[17]*mnuSelf+0.5000000000000001*m1rSelf[17]*mnuSelf-0.2*m0rSelf[3]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[15]*mnuSelf-0.25*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[13]*mnuSelf-0.2*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.25*m0rSelf[5]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(55,12) = (-0.3273268353539885*m0rSelf[5]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[0]*uSelf[19]*mnuSelf+0.4391550328268398*m1rSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[12]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(55,13) = (-0.3273268353539885*m0rSelf[7]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[1]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[17]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[16]*mnuSelf-0.2*m0rSelf[1]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[13]*mnuSelf-0.2*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[11]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(55,14) = (-0.2195775164134199*m0rSelf[6]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[7]*uSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[17]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[16]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf-0.25*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[5]*uSelf[14]*mnuSelf-0.2*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[11]*mnuSelf; 
  data->AEM_S(55,15) = (-0.4621212121212121*m0rSelf[9]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[2]*uSelf[19]*mnuSelf-0.25*m0rSelf[8]*uSelf[18]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[16]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[15]*mnuSelf+0.31943828249997*m1rSelf[15]*mnuSelf-0.25*m0rSelf[4]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[11]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(55,16) = (-0.21957751641342*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[3]*uSelf[18]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[16]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[16]*mnuSelf+0.4472135954999579*m1rSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[13]*mnuSelf-0.2*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[12]*mnuSelf-0.2*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(55,17) = (-0.3273268353539885*m0rSelf[3]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[4]*uSelf[18]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[17]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[17]*mnuSelf+0.31943828249997*m1rSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[16]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[11]*mnuSelf+0.5000000000000001*m1rSelf[11]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(55,18) = (-0.25*m0rSelf[5]*uSelf[18]*mnuSelf)-0.21957751641342*m0rSelf[4]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[16]*mnuSelf-0.25*m0rSelf[8]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[7]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[13]*mnuSelf; 
  data->AEM_S(55,19) = (-0.4621212121212121*m0rSelf[5]*uSelf[19]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[19]*mnuSelf+0.2981423969999719*m1rSelf[19]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[16]*mnuSelf-0.4621212121212121*m0rSelf[9]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[12]*mnuSelf+0.4391550328268398*m1rSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[11]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[10]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,10) = (-0.2195775164134199*m0rSelf[3]*uSelf[18]*mnuSelf)-0.2*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[13]*mnuSelf-0.2*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.25*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,11) = (-0.1963961012123931*m0rSelf[7]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[6]*uSelf[18]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[17]*mnuSelf-0.2*m0rSelf[2]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[16]*mnuSelf-0.2*m0rSelf[3]*uSelf[15]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[14]*mnuSelf-0.2*m0rSelf[5]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.447213595499958*m1rSelf[13]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,12) = (-0.1963961012123931*m0rSelf[6]*uSelf[19]*mnuSelf)-0.149071198499986*m0rSelf[8]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[18]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[17]*mnuSelf-0.2*m0rSelf[1]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[16]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[16]*mnuSelf-0.45*m0rSelf[2]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[14]*mnuSelf+0.5000000000000001*m1rSelf[14]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[13]*mnuSelf-0.45*m0rSelf[6]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[11]*mnuSelf-0.2*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,13) = (-0.175662013130736*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[5]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[18]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[18]*mnuSelf+0.43915503282684*m1rSelf[18]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[17]*mnuSelf-0.2*m0rSelf[0]*uSelf[17]*mnuSelf+0.4*m1rSelf[17]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[15]*mnuSelf-0.2*m0rSelf[1]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[9]*uSelf[13]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[12]*mnuSelf-0.2*m0rSelf[5]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[10]*mnuSelf-0.2*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,14) = (-0.2195775164134199*m0rSelf[5]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[3]*uSelf[18]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[17]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[16]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[16]*mnuSelf+0.31943828249997*m1rSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[12]*mnuSelf+0.5000000000000001*m1rSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[11]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,15) = (-0.21957751641342*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[3]*uSelf[18]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[16]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[16]*mnuSelf+0.4472135954999579*m1rSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[13]*mnuSelf-0.2*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[12]*mnuSelf-0.2*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,16) = (-0.3833333333333334*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1254728665219542*m0rSelf[6]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[19]*mnuSelf-0.4621212121212121*m0rSelf[8]*uSelf[18]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[18]*mnuSelf-0.29277002188456*m0rSelf[8]*uSelf[17]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[17]*mnuSelf-0.1254728665219542*m0rSelf[9]*uSelf[16]*mnuSelf-0.9642857142857143*m0rSelf[6]*uSelf[16]*mnuSelf-0.2874944542499729*m0rSelf[2]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[14]*mnuSelf+0.31943828249997*m1rSelf[14]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[12]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[12]*mnuSelf-0.45*m0rSelf[2]*uSelf[12]*mnuSelf-0.3273268353539885*m0rSelf[8]*uSelf[11]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(56,17) = (-0.1928571428571429*m0rSelf[8]*uSelf[19]*mnuSelf)-0.29277002188456*m0rSelf[7]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[19]*mnuSelf-0.1928571428571429*m0rSelf[9]*uSelf[18]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[18]*mnuSelf-0.29277002188456*m0rSelf[9]*uSelf[17]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[8]*uSelf[16]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[13]*mnuSelf-0.2*m0rSelf[0]*uSelf[13]*mnuSelf+0.4*m1rSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[12]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[12]*mnuSelf-0.2*m0rSelf[1]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[11]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[11]*mnuSelf-0.2*m0rSelf[2]*uSelf[11]*mnuSelf-0.2*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,18) = (-0.1928571428571429*m0rSelf[7]*uSelf[19]*mnuSelf)-0.4621212121212121*m0rSelf[6]*uSelf[18]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[18]*mnuSelf-0.1928571428571429*m0rSelf[9]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[17]*mnuSelf-0.4621212121212121*m0rSelf[8]*uSelf[16]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[13]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[13]*mnuSelf+0.43915503282684*m1rSelf[13]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[12]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(56,19) = (-0.3833333333333334*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1928571428571429*m0rSelf[7]*uSelf[18]*mnuSelf-0.1928571428571429*m0rSelf[8]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[17]*mnuSelf-0.3833333333333334*m0rSelf[9]*uSelf[16]*mnuSelf-0.1254728665219543*m0rSelf[6]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[11]*mnuSelf; 
  data->AEM_S(57,10) = (-0.2195775164134199*m0rSelf[3]*uSelf[19]*mnuSelf)-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.2*m0rSelf[3]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[9]*uSelf[13]*mnuSelf-0.2*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[11]*mnuSelf-0.25*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(57,11) = (-0.149071198499986*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[6]*uSelf[19]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[17]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[17]*mnuSelf-0.45*m0rSelf[1]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[16]*mnuSelf-0.2*m0rSelf[2]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[15]*mnuSelf+0.5000000000000001*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[13]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[12]*mnuSelf-0.2*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.45*m0rSelf[7]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(57,12) = (-0.3273268353539885*m0rSelf[7]*uSelf[19]*mnuSelf)-0.21957751641342*m0rSelf[1]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[18]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[16]*mnuSelf-0.2*m0rSelf[1]*uSelf[16]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[13]*mnuSelf-0.2*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.447213595499958*m1rSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[11]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(57,13) = (-0.3273268353539885*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[4]*uSelf[19]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[19]*mnuSelf+0.43915503282684*m1rSelf[19]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[18]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[16]*mnuSelf-0.2*m0rSelf[0]*uSelf[16]*mnuSelf+0.4*m1rSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[15]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[14]*mnuSelf-0.2*m0rSelf[2]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[8]*uSelf[13]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[13]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[12]*mnuSelf-0.2*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[12]*mnuSelf+0.447213595499958*m1rSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[9]*uSelf[10]*mnuSelf-0.2*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(57,14) = (-0.1963961012123931*m0rSelf[3]*uSelf[19]*mnuSelf)-0.21957751641342*m0rSelf[5]*uSelf[18]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.4472135954999579*m1rSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[8]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[13]*mnuSelf-0.2*m0rSelf[2]*uSelf[13]*mnuSelf-0.2*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(57,15) = (-0.3273268353539885*m0rSelf[3]*uSelf[19]*mnuSelf)-0.2195775164134199*m0rSelf[4]*uSelf[18]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[17]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[17]*mnuSelf+0.31943828249997*m1rSelf[17]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[16]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[8]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[11]*mnuSelf+0.5000000000000001*m1rSelf[11]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(57,16) = (-0.1928571428571429*m0rSelf[8]*uSelf[19]*mnuSelf)-0.29277002188456*m0rSelf[7]*uSelf[19]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[19]*mnuSelf-0.1928571428571429*m0rSelf[9]*uSelf[18]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[18]*mnuSelf-0.29277002188456*m0rSelf[9]*uSelf[17]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[17]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[8]*uSelf[16]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[16]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[13]*mnuSelf-0.2*m0rSelf[0]*uSelf[13]*mnuSelf+0.4*m1rSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[12]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[12]*mnuSelf-0.2*m0rSelf[1]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[9]*uSelf[11]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[11]*mnuSelf-0.2*m0rSelf[2]*uSelf[11]*mnuSelf-0.2*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(57,17) = (-0.4621212121212121*m0rSelf[9]*uSelf[19]*mnuSelf)-0.29277002188456*m0rSelf[6]*uSelf[19]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[19]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[18]*mnuSelf-0.1254728665219542*m0rSelf[7]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[18]*mnuSelf-0.1254728665219542*m0rSelf[8]*uSelf[17]*mnuSelf-0.9642857142857143*m0rSelf[7]*uSelf[17]*mnuSelf-0.2874944542499729*m0rSelf[1]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[9]*uSelf[16]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[16]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[16]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[15]*mnuSelf+0.31943828249997*m1rSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[9]*uSelf[12]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[8]*uSelf[11]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[11]*mnuSelf-0.45*m0rSelf[1]*uSelf[11]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(57,18) = (-0.1928571428571429*m0rSelf[6]*uSelf[19]*mnuSelf)-0.3833333333333334*m0rSelf[7]*uSelf[18]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[17]*mnuSelf-0.1254728665219543*m0rSelf[7]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[17]*mnuSelf-0.1928571428571429*m0rSelf[9]*uSelf[16]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[11]*mnuSelf; 
  data->AEM_S(57,19) = (-0.4621212121212121*m0rSelf[7]*uSelf[19]*mnuSelf)-0.149071198499986*m0rSelf[1]*uSelf[19]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[18]*mnuSelf-0.4621212121212121*m0rSelf[9]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[17]*mnuSelf-0.1928571428571429*m0rSelf[8]*uSelf[16]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[13]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[13]*mnuSelf+0.43915503282684*m1rSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[12]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[11]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(58,10) = (-0.149071198499986*m0rSelf[4]*uSelf[18]*mnuSelf)-0.25*m0rSelf[0]*uSelf[18]*mnuSelf+0.5*m1rSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[16]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[8]*uSelf[10]*mnuSelf; 
  data->AEM_S(58,11) = (-0.2499586742703185*m0rSelf[8]*uSelf[18]*mnuSelf)-0.3833333333333334*m0rSelf[1]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[14]*mnuSelf+0.4391550328268398*m1rSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[12]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[11]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(58,12) = (-0.149071198499986*m0rSelf[6]*uSelf[18]*mnuSelf)-0.25*m0rSelf[2]*uSelf[18]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[17]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[8]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(58,13) = (-0.3833333333333334*m0rSelf[3]*uSelf[18]*mnuSelf)-0.175662013130736*m0rSelf[3]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[16]*mnuSelf+0.43915503282684*m1rSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[14]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[14]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[13]*mnuSelf-0.175662013130736*m0rSelf[7]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(58,14) = (-0.4621212121212121*m0rSelf[4]*uSelf[18]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[18]*mnuSelf+0.2981423969999719*m1rSelf[18]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[15]*mnuSelf-0.4621212121212121*m0rSelf[8]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[12]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[11]*mnuSelf+0.4391550328268398*m1rSelf[11]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[10]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(58,15) = (-0.25*m0rSelf[5]*uSelf[18]*mnuSelf)-0.21957751641342*m0rSelf[4]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[16]*mnuSelf-0.25*m0rSelf[8]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[7]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[13]*mnuSelf; 
  data->AEM_S(58,16) = (-0.1928571428571429*m0rSelf[7]*uSelf[19]*mnuSelf)-0.4621212121212121*m0rSelf[6]*uSelf[18]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[18]*mnuSelf-0.1928571428571429*m0rSelf[9]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[17]*mnuSelf-0.4621212121212121*m0rSelf[8]*uSelf[16]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[1]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[5]*uSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[4]*uSelf[13]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[13]*mnuSelf+0.43915503282684*m1rSelf[13]*mnuSelf-0.149071198499986*m0rSelf[8]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[12]*mnuSelf-0.3273268353539885*m0rSelf[6]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(58,17) = (-0.1928571428571429*m0rSelf[6]*uSelf[19]*mnuSelf)-0.3833333333333334*m0rSelf[7]*uSelf[18]*mnuSelf-0.3833333333333334*m0rSelf[8]*uSelf[17]*mnuSelf-0.1254728665219543*m0rSelf[7]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[17]*mnuSelf-0.1928571428571429*m0rSelf[9]*uSelf[16]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[11]*mnuSelf; 
  data->AEM_S(58,18) = (-0.25*m0rSelf[9]*uSelf[19]*mnuSelf)-0.5898601398601399*m0rSelf[8]*uSelf[18]*mnuSelf-0.2499586742703185*m0rSelf[1]*uSelf[18]*mnuSelf-0.3833333333333334*m0rSelf[7]*uSelf[17]*mnuSelf-0.4621212121212121*m0rSelf[6]*uSelf[16]*mnuSelf-0.149071198499986*m0rSelf[2]*uSelf[16]*mnuSelf-0.25*m0rSelf[5]*uSelf[15]*mnuSelf-0.4621212121212121*m0rSelf[4]*uSelf[14]*mnuSelf-0.149071198499986*m0rSelf[0]*uSelf[14]*mnuSelf+0.2981423969999719*m1rSelf[14]*mnuSelf-0.3833333333333334*m0rSelf[3]*uSelf[13]*mnuSelf-0.149071198499986*m0rSelf[6]*uSelf[12]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf-0.2499586742703185*m0rSelf[8]*uSelf[11]*mnuSelf-0.3833333333333334*m0rSelf[1]*uSelf[11]*mnuSelf-0.149071198499986*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(58,19) = (-0.25*m0rSelf[8]*uSelf[19]*mnuSelf)-0.25*m0rSelf[9]*uSelf[18]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[17]*mnuSelf-0.1928571428571429*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,10) = (-0.149071198499986*m0rSelf[5]*uSelf[19]*mnuSelf)-0.25*m0rSelf[0]*uSelf[19]*mnuSelf+0.5*m1rSelf[19]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[17]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[12]*mnuSelf-0.25*m0rSelf[9]*uSelf[10]*mnuSelf; 
  data->AEM_S(59,11) = (-0.149071198499986*m0rSelf[7]*uSelf[19]*mnuSelf)-0.25*m0rSelf[1]*uSelf[19]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[16]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[12]*mnuSelf-0.25*m0rSelf[9]*uSelf[11]*mnuSelf; 
  data->AEM_S(59,12) = (-0.2499586742703185*m0rSelf[9]*uSelf[19]*mnuSelf)-0.3833333333333334*m0rSelf[2]*uSelf[19]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[15]*mnuSelf+0.4391550328268398*m1rSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[13]*mnuSelf-0.3833333333333334*m0rSelf[9]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[7]*uSelf[11]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(59,13) = (-0.3833333333333334*m0rSelf[3]*uSelf[19]*mnuSelf)-0.3273268353539885*m0rSelf[5]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[17]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[17]*mnuSelf+0.43915503282684*m1rSelf[17]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[14]*mnuSelf-0.3833333333333334*m0rSelf[9]*uSelf[13]*mnuSelf-0.175662013130736*m0rSelf[6]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[5]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(59,14) = (-0.25*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1963961012123931*m0rSelf[3]*uSelf[17]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[6]*uSelf[15]*mnuSelf-0.25*m0rSelf[9]*uSelf[14]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[13]*mnuSelf; 
  data->AEM_S(59,15) = (-0.4621212121212121*m0rSelf[5]*uSelf[19]*mnuSelf)-0.149071198499986*m0rSelf[0]*uSelf[19]*mnuSelf+0.2981423969999719*m1rSelf[19]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[17]*mnuSelf-0.2195775164134199*m0rSelf[4]*uSelf[16]*mnuSelf-0.4621212121212121*m0rSelf[9]*uSelf[15]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[15]*mnuSelf-0.2195775164134199*m0rSelf[6]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[13]*mnuSelf-0.2195775164134199*m0rSelf[1]*uSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[0]*uSelf[12]*mnuSelf+0.4391550328268398*m1rSelf[12]*mnuSelf-0.2195775164134199*m0rSelf[3]*uSelf[11]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[10]*mnuSelf-0.2195775164134199*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(59,16) = (-0.3833333333333334*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1928571428571429*m0rSelf[7]*uSelf[18]*mnuSelf-0.1928571428571429*m0rSelf[8]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[17]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[17]*mnuSelf-0.3833333333333334*m0rSelf[9]*uSelf[16]*mnuSelf-0.1254728665219543*m0rSelf[6]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[2]*uSelf[16]*mnuSelf-0.21957751641342*m0rSelf[4]*uSelf[15]*mnuSelf-0.21957751641342*m0rSelf[5]*uSelf[14]*mnuSelf-0.175662013130736*m0rSelf[3]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[12]*mnuSelf-0.1963961012123931*m0rSelf[7]*uSelf[11]*mnuSelf; 
  data->AEM_S(59,17) = (-0.4621212121212121*m0rSelf[7]*uSelf[19]*mnuSelf)-0.149071198499986*m0rSelf[1]*uSelf[19]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[18]*mnuSelf-0.4621212121212121*m0rSelf[9]*uSelf[17]*mnuSelf-0.29277002188456*m0rSelf[6]*uSelf[17]*mnuSelf-0.3273268353539885*m0rSelf[2]*uSelf[17]*mnuSelf-0.1928571428571429*m0rSelf[8]*uSelf[16]*mnuSelf-0.29277002188456*m0rSelf[7]*uSelf[16]*mnuSelf-0.1963961012123931*m0rSelf[1]*uSelf[16]*mnuSelf-0.3273268353539885*m0rSelf[3]*uSelf[15]*mnuSelf-0.1963961012123931*m0rSelf[3]*uSelf[14]*mnuSelf-0.3273268353539885*m0rSelf[5]*uSelf[13]*mnuSelf-0.1963961012123931*m0rSelf[4]*uSelf[13]*mnuSelf-0.21957751641342*m0rSelf[0]*uSelf[13]*mnuSelf+0.43915503282684*m1rSelf[13]*mnuSelf-0.3273268353539885*m0rSelf[7]*uSelf[12]*mnuSelf-0.21957751641342*m0rSelf[1]*uSelf[12]*mnuSelf-0.149071198499986*m0rSelf[9]*uSelf[11]*mnuSelf-0.1963961012123931*m0rSelf[6]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[2]*uSelf[11]*mnuSelf-0.21957751641342*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(59,18) = (-0.25*m0rSelf[8]*uSelf[19]*mnuSelf)-0.25*m0rSelf[9]*uSelf[18]*mnuSelf-0.1928571428571429*m0rSelf[6]*uSelf[17]*mnuSelf-0.1928571428571429*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,19) = (-0.5898601398601399*m0rSelf[9]*uSelf[19]*mnuSelf)-0.2499586742703185*m0rSelf[2]*uSelf[19]*mnuSelf-0.25*m0rSelf[8]*uSelf[18]*mnuSelf-0.4621212121212121*m0rSelf[7]*uSelf[17]*mnuSelf-0.149071198499986*m0rSelf[1]*uSelf[17]*mnuSelf-0.3833333333333334*m0rSelf[6]*uSelf[16]*mnuSelf-0.4621212121212121*m0rSelf[5]*uSelf[15]*mnuSelf-0.149071198499986*m0rSelf[0]*uSelf[15]*mnuSelf+0.2981423969999719*m1rSelf[15]*mnuSelf-0.25*m0rSelf[4]*uSelf[14]*mnuSelf-0.3833333333333334*m0rSelf[3]*uSelf[13]*mnuSelf-0.2499586742703185*m0rSelf[9]*uSelf[12]*mnuSelf-0.3833333333333334*m0rSelf[2]*uSelf[12]*mnuSelf-0.149071198499986*m0rSelf[7]*uSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[11]*mnuSelf-0.149071198499986*m0rSelf[5]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(50,40) = 0.25*m0rOther[9]*uOther[19]*mnuOther+0.25*m0rOther[8]*uOther[18]*mnuOther+0.25*m0rOther[7]*uOther[17]*mnuOther+0.25*m0rOther[6]*uOther[16]*mnuOther+0.25*m0rOther[5]*uOther[15]*mnuOther+0.25*m0rOther[4]*uOther[14]*mnuOther+0.25*m0rOther[3]*uOther[13]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther+0.25*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(50,41) = 0.2195775164134199*m0rOther[4]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[6]*uOther[13]*mnuOther+0.25*m0rOther[2]*uOther[13]*mnuOther+0.25*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(50,42) = 0.2195775164134199*m0rOther[5]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[7]*uOther[13]*mnuOther+0.25*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther+0.25*m0rOther[3]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(50,43) = 0.2195775164134199*m0rOther[7]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[17]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[16]*mnuOther+0.2*m0rOther[7]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther+0.223606797749979*m0rOther[6]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(50,44) = 0.149071198499986*m0rOther[8]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[18]*mnuOther+0.223606797749979*m0rOther[7]*uOther[17]*mnuOther+0.159719141249985*m0rOther[6]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[16]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(50,45) = 0.149071198499986*m0rOther[9]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[19]*mnuOther+0.159719141249985*m0rOther[7]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[6]*uOther[16]*mnuOther+0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[11]*mnuOther+0.25*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(50,46) = 0.2195775164134199*m0rOther[3]*uOther[18]*mnuOther+0.2*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[5]*uOther[16]*mnuOther+0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther+0.223606797749979*m0rOther[6]*uOther[15]*mnuOther+0.159719141249985*m0rOther[6]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[13]*mnuOther+0.2*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.25*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(50,47) = 0.2195775164134199*m0rOther[3]*uOther[19]*mnuOther+0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther+0.159719141249985*m0rOther[7]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[7]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[13]*mnuOther+0.2*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(50,48) = 0.149071198499986*m0rOther[4]*uOther[18]*mnuOther+0.25*m0rOther[0]*uOther[18]*mnuOther-0.5*m1rOther[18]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[16]*mnuOther+0.149071198499986*m0rOther[8]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[8]*uOther[10]*mnuOther; 
  data->AEM_S(50,49) = 0.149071198499986*m0rOther[5]*uOther[19]*mnuOther+0.25*m0rOther[0]*uOther[19]*mnuOther-0.5*m1rOther[19]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[17]*mnuOther+0.149071198499986*m0rOther[9]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[12]*mnuOther+0.25*m0rOther[9]*uOther[10]*mnuOther; 
  data->AEM_S(51,40) = 0.2195775164134199*m0rOther[4]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[6]*uOther[13]*mnuOther+0.25*m0rOther[2]*uOther[13]*mnuOther+0.25*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(51,41) = 0.25*m0rOther[9]*uOther[19]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[18]*mnuOther+0.45*m0rOther[7]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[16]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther+0.25*m0rOther[5]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.4472135954999579*m1rOther[14]*mnuOther+0.45*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[6]*uOther[12]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[11]*mnuOther+0.45*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(51,42) = 0.2195775164134199*m0rOther[7]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[17]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[16]*mnuOther+0.2*m0rOther[7]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther+0.223606797749979*m0rOther[6]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(51,43) = 0.2195775164134199*m0rOther[5]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[18]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[17]*mnuOther+0.2*m0rOther[5]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[16]*mnuOther+0.223606797749979*m0rOther[0]*uOther[16]*mnuOther-0.447213595499958*m1rOther[16]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[15]*mnuOther+0.2*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[13]*mnuOther+0.45*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.223606797749979*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther+0.45*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[6]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(51,44) = 0.3273268353539885*m0rOther[4]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[18]*mnuOther-0.4391550328268398*m1rOther[18]*mnuOther+0.223606797749979*m0rOther[5]*uOther[17]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[16]*mnuOther+0.223606797749979*m0rOther[7]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(51,45) = 0.2195775164134199*m0rOther[3]*uOther[19]*mnuOther+0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[17]*mnuOther-0.5000000000000001*m1rOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther+0.159719141249985*m0rOther[7]*uOther[15]*mnuOther+0.25*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[7]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[13]*mnuOther+0.2*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.25*m0rOther[5]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(51,46) = 0.1963961012123931*m0rOther[7]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[18]*mnuOther+0.21957751641342*m0rOther[2]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[17]*mnuOther+0.351382110749967*m0rOther[6]*uOther[17]*mnuOther+0.2*m0rOther[2]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[16]*mnuOther+0.351382110749967*m0rOther[7]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[16]*mnuOther+0.2*m0rOther[3]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[14]*mnuOther+0.2*m0rOther[5]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.447213595499958*m1rOther[13]*mnuOther+0.21957751641342*m0rOther[8]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(51,47) = 0.149071198499986*m0rOther[9]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[19]*mnuOther+0.21957751641342*m0rOther[2]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[17]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[17]*mnuOther+0.45*m0rOther[1]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[16]*mnuOther+0.351382110749967*m0rOther[6]*uOther[16]*mnuOther+0.2*m0rOther[2]*uOther[16]*mnuOther+0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[15]*mnuOther-0.5000000000000001*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[13]*mnuOther+0.21957751641342*m0rOther[9]*uOther[12]*mnuOther+0.2*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.45*m0rOther[7]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(51,48) = 0.2499586742703185*m0rOther[8]*uOther[18]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[14]*mnuOther-0.4391550328268398*m1rOther[14]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[12]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[11]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(51,49) = 0.149071198499986*m0rOther[7]*uOther[19]*mnuOther+0.25*m0rOther[1]*uOther[19]*mnuOther+0.149071198499986*m0rOther[9]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[12]*mnuOther+0.25*m0rOther[9]*uOther[11]*mnuOther; 
  data->AEM_S(52,40) = 0.2195775164134199*m0rOther[5]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[7]*uOther[13]*mnuOther+0.25*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther+0.25*m0rOther[3]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(52,41) = 0.2195775164134199*m0rOther[7]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[17]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[16]*mnuOther+0.2*m0rOther[7]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther+0.223606797749979*m0rOther[6]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(52,42) = 0.3833333333333334*m0rOther[9]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[19]*mnuOther+0.25*m0rOther[8]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.45*m0rOther[6]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.25*m0rOther[4]*uOther[14]*mnuOther+0.45*m0rOther[3]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[12]*mnuOther+0.45*m0rOther[2]*uOther[12]*mnuOther+0.223606797749979*m0rOther[7]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(52,43) = 0.1963961012123931*m0rOther[3]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[18]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[17]*mnuOther+0.2*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.447213595499958*m1rOther[17]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[14]*mnuOther+0.2*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[13]*mnuOther+0.45*m0rOther[2]*uOther[13]*mnuOther+0.45*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(52,44) = 0.2195775164134199*m0rOther[3]*uOther[18]*mnuOther+0.2*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[5]*uOther[16]*mnuOther+0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[16]*mnuOther-0.5000000000000001*m1rOther[16]*mnuOther+0.223606797749979*m0rOther[6]*uOther[15]*mnuOther+0.159719141249985*m0rOther[6]*uOther[14]*mnuOther+0.25*m0rOther[2]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[13]*mnuOther+0.2*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.25*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(52,45) = 0.3273268353539885*m0rOther[5]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[19]*mnuOther-0.4391550328268398*m1rOther[19]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[15]*mnuOther+0.223606797749979*m0rOther[6]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[12]*mnuOther+0.223606797749979*m0rOther[0]*uOther[12]*mnuOther-0.4472135954999579*m1rOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(52,46) = 0.1963961012123931*m0rOther[6]*uOther[19]*mnuOther+0.149071198499986*m0rOther[8]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[18]*mnuOther+0.21957751641342*m0rOther[1]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[17]*mnuOther+0.351382110749967*m0rOther[7]*uOther[17]*mnuOther+0.2*m0rOther[1]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[16]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[16]*mnuOther+0.45*m0rOther[2]*uOther[16]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[14]*mnuOther-0.5000000000000001*m1rOther[14]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[13]*mnuOther+0.45*m0rOther[6]*uOther[12]*mnuOther+0.21957751641342*m0rOther[8]*uOther[11]*mnuOther+0.2*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(52,47) = 0.3273268353539885*m0rOther[7]*uOther[19]*mnuOther+0.21957751641342*m0rOther[1]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[17]*mnuOther+0.351382110749967*m0rOther[6]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[16]*mnuOther+0.351382110749967*m0rOther[7]*uOther[16]*mnuOther+0.2*m0rOther[1]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[13]*mnuOther+0.2*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.447213595499958*m1rOther[13]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.21957751641342*m0rOther[9]*uOther[11]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(52,48) = 0.149071198499986*m0rOther[6]*uOther[18]*mnuOther+0.25*m0rOther[2]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[17]*mnuOther+0.149071198499986*m0rOther[8]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[8]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(52,49) = 0.2499586742703185*m0rOther[9]*uOther[19]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[15]*mnuOther-0.4391550328268398*m1rOther[15]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[13]*mnuOther+0.3833333333333334*m0rOther[9]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(53,40) = 0.2195775164134199*m0rOther[7]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[17]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[16]*mnuOther+0.2*m0rOther[7]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther+0.223606797749979*m0rOther[6]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(53,41) = 0.2195775164134199*m0rOther[5]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[18]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[17]*mnuOther+0.2*m0rOther[5]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[16]*mnuOther+0.223606797749979*m0rOther[0]*uOther[16]*mnuOther-0.447213595499958*m1rOther[16]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[15]*mnuOther+0.2*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[13]*mnuOther+0.45*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.223606797749979*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther+0.45*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[6]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(53,42) = 0.1963961012123931*m0rOther[3]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[18]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[17]*mnuOther+0.2*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.447213595499958*m1rOther[17]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[14]*mnuOther+0.2*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[13]*mnuOther+0.45*m0rOther[2]*uOther[13]*mnuOther+0.45*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(53,43) = 0.3833333333333334*m0rOther[9]*uOther[19]*mnuOther+0.175662013130736*m0rOther[6]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[19]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[18]*mnuOther+0.175662013130736*m0rOther[7]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[18]*mnuOther+0.175662013130736*m0rOther[8]*uOther[17]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[17]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[17]*mnuOther+0.175662013130736*m0rOther[9]*uOther[16]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[16]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[15]*mnuOther+0.2*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.2*m0rOther[5]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.4472135954999579*m1rOther[14]*mnuOther+0.81*m0rOther[3]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[12]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[12]*mnuOther+0.45*m0rOther[2]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[11]*mnuOther+0.45*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.223606797749979*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(53,44) = 0.1963961012123931*m0rOther[7]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[17]*mnuOther+0.2*m0rOther[2]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[16]*mnuOther+0.2*m0rOther[3]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[14]*mnuOther+0.2*m0rOther[5]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(53,45) = 0.3273268353539885*m0rOther[7]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[17]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[16]*mnuOther+0.2*m0rOther[1]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[13]*mnuOther+0.2*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[11]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(53,46) = 0.175662013130736*m0rOther[3]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[18]*mnuOther+0.21957751641342*m0rOther[0]*uOther[18]*mnuOther-0.43915503282684*m1rOther[18]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[17]*mnuOther+0.2*m0rOther[0]*uOther[17]*mnuOther-0.4*m1rOther[17]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[15]*mnuOther+0.2*m0rOther[1]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[14]*mnuOther+0.175662013130736*m0rOther[9]*uOther[13]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[12]*mnuOther+0.2*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.447213595499958*m1rOther[11]*mnuOther+0.21957751641342*m0rOther[8]*uOther[10]*mnuOther+0.2*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(53,47) = 0.3273268353539885*m0rOther[5]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[19]*mnuOther+0.21957751641342*m0rOther[0]*uOther[19]*mnuOther-0.43915503282684*m1rOther[19]*mnuOther+0.175662013130736*m0rOther[3]*uOther[18]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[16]*mnuOther+0.2*m0rOther[0]*uOther[16]*mnuOther-0.4*m1rOther[16]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[15]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[14]*mnuOther+0.2*m0rOther[2]*uOther[14]*mnuOther+0.175662013130736*m0rOther[8]*uOther[13]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[12]*mnuOther+0.2*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[0]*uOther[12]*mnuOther-0.447213595499958*m1rOther[12]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[11]*mnuOther+0.21957751641342*m0rOther[9]*uOther[10]*mnuOther+0.2*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(53,48) = 0.3833333333333334*m0rOther[3]*uOther[18]*mnuOther+0.175662013130736*m0rOther[3]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[16]*mnuOther+0.21957751641342*m0rOther[0]*uOther[16]*mnuOther-0.43915503282684*m1rOther[16]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[14]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[13]*mnuOther+0.175662013130736*m0rOther[7]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[11]*mnuOther+0.21957751641342*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(53,49) = 0.3833333333333334*m0rOther[3]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[17]*mnuOther+0.21957751641342*m0rOther[0]*uOther[17]*mnuOther-0.43915503282684*m1rOther[17]*mnuOther+0.175662013130736*m0rOther[3]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[15]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[14]*mnuOther+0.3833333333333334*m0rOther[9]*uOther[13]*mnuOther+0.175662013130736*m0rOther[6]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[11]*mnuOther+0.21957751641342*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(54,40) = 0.149071198499986*m0rOther[8]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[18]*mnuOther+0.223606797749979*m0rOther[7]*uOther[17]*mnuOther+0.159719141249985*m0rOther[6]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[16]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(54,41) = 0.3273268353539885*m0rOther[4]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[18]*mnuOther-0.4391550328268398*m1rOther[18]*mnuOther+0.223606797749979*m0rOther[5]*uOther[17]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[16]*mnuOther+0.223606797749979*m0rOther[7]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(54,42) = 0.2195775164134199*m0rOther[3]*uOther[18]*mnuOther+0.2*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[5]*uOther[16]*mnuOther+0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[16]*mnuOther-0.5000000000000001*m1rOther[16]*mnuOther+0.223606797749979*m0rOther[6]*uOther[15]*mnuOther+0.159719141249985*m0rOther[6]*uOther[14]*mnuOther+0.25*m0rOther[2]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[13]*mnuOther+0.2*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.25*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(54,43) = 0.1963961012123931*m0rOther[7]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[17]*mnuOther+0.2*m0rOther[2]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[16]*mnuOther+0.2*m0rOther[3]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[14]*mnuOther+0.2*m0rOther[5]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(54,44) = 0.25*m0rOther[9]*uOther[19]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[17]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[16]*mnuOther+0.159719141249985*m0rOther[2]*uOther[16]*mnuOther+0.25*m0rOther[5]*uOther[15]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[14]*mnuOther+0.159719141249985*m0rOther[0]*uOther[14]*mnuOther-0.31943828249997*m1rOther[14]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[13]*mnuOther+0.159719141249985*m0rOther[6]*uOther[12]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[11]*mnuOther+0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(54,45) = 0.2195775164134199*m0rOther[6]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[17]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[16]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[16]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther+0.25*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[5]*uOther[14]*mnuOther+0.2*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[7]*uOther[11]*mnuOther; 
  data->AEM_S(54,46) = 0.2195775164134199*m0rOther[5]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[18]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[17]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[16]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[16]*mnuOther+0.159719141249985*m0rOther[0]*uOther[16]*mnuOther-0.31943828249997*m1rOther[16]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[14]*mnuOther+0.159719141249985*m0rOther[2]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.159719141249985*m0rOther[4]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[12]*mnuOther-0.5000000000000001*m1rOther[12]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[11]*mnuOther+0.159719141249985*m0rOther[6]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(54,47) = 0.1963961012123931*m0rOther[3]*uOther[19]*mnuOther+0.21957751641342*m0rOther[5]*uOther[18]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.4472135954999579*m1rOther[17]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[16]*mnuOther+0.21957751641342*m0rOther[8]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[13]*mnuOther+0.2*m0rOther[2]*uOther[13]*mnuOther+0.2*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(54,48) = 0.4621212121212121*m0rOther[4]*uOther[18]*mnuOther+0.149071198499986*m0rOther[0]*uOther[18]*mnuOther-0.2981423969999719*m1rOther[18]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[15]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[12]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[11]*mnuOther-0.4391550328268398*m1rOther[11]*mnuOther+0.149071198499986*m0rOther[8]*uOther[10]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(54,49) = 0.25*m0rOther[4]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[17]*mnuOther+0.21957751641342*m0rOther[5]*uOther[16]*mnuOther+0.21957751641342*m0rOther[6]*uOther[15]*mnuOther+0.25*m0rOther[9]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[13]*mnuOther; 
  data->AEM_S(55,40) = 0.149071198499986*m0rOther[9]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[19]*mnuOther+0.159719141249985*m0rOther[7]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[6]*uOther[16]*mnuOther+0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[11]*mnuOther+0.25*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(55,41) = 0.2195775164134199*m0rOther[3]*uOther[19]*mnuOther+0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[17]*mnuOther-0.5000000000000001*m1rOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther+0.159719141249985*m0rOther[7]*uOther[15]*mnuOther+0.25*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[7]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[13]*mnuOther+0.2*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.25*m0rOther[5]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(55,42) = 0.3273268353539885*m0rOther[5]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[19]*mnuOther-0.4391550328268398*m1rOther[19]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[15]*mnuOther+0.223606797749979*m0rOther[6]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[12]*mnuOther+0.223606797749979*m0rOther[0]*uOther[12]*mnuOther-0.4472135954999579*m1rOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(55,43) = 0.3273268353539885*m0rOther[7]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[17]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[16]*mnuOther+0.2*m0rOther[1]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[13]*mnuOther+0.2*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[11]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(55,44) = 0.2195775164134199*m0rOther[6]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[18]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[17]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[16]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[16]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther+0.25*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[5]*uOther[14]*mnuOther+0.2*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[7]*uOther[11]*mnuOther; 
  data->AEM_S(55,45) = 0.4621212121212121*m0rOther[9]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[19]*mnuOther+0.25*m0rOther[8]*uOther[18]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[17]*mnuOther+0.159719141249985*m0rOther[1]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[16]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[15]*mnuOther+0.159719141249985*m0rOther[0]*uOther[15]*mnuOther-0.31943828249997*m1rOther[15]*mnuOther+0.25*m0rOther[4]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[13]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[12]*mnuOther+0.159719141249985*m0rOther[7]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[11]*mnuOther+0.159719141249985*m0rOther[5]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(55,46) = 0.21957751641342*m0rOther[4]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[18]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[16]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[16]*mnuOther+0.223606797749979*m0rOther[0]*uOther[16]*mnuOther-0.4472135954999579*m1rOther[16]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[15]*mnuOther+0.21957751641342*m0rOther[9]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[13]*mnuOther+0.2*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[12]*mnuOther+0.2*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(55,47) = 0.3273268353539885*m0rOther[3]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[18]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[17]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[17]*mnuOther+0.159719141249985*m0rOther[0]*uOther[17]*mnuOther-0.31943828249997*m1rOther[17]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[16]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[15]*mnuOther+0.159719141249985*m0rOther[1]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[12]*mnuOther+0.159719141249985*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[11]*mnuOther-0.5000000000000001*m1rOther[11]*mnuOther+0.159719141249985*m0rOther[7]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(55,48) = 0.25*m0rOther[5]*uOther[18]*mnuOther+0.21957751641342*m0rOther[4]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[16]*mnuOther+0.25*m0rOther[8]*uOther[15]*mnuOther+0.21957751641342*m0rOther[7]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[13]*mnuOther; 
  data->AEM_S(55,49) = 0.4621212121212121*m0rOther[5]*uOther[19]*mnuOther+0.149071198499986*m0rOther[0]*uOther[19]*mnuOther-0.2981423969999719*m1rOther[19]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[16]*mnuOther+0.4621212121212121*m0rOther[9]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[13]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[12]*mnuOther-0.4391550328268398*m1rOther[12]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[11]*mnuOther+0.149071198499986*m0rOther[9]*uOther[10]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(56,40) = 0.2195775164134199*m0rOther[3]*uOther[18]*mnuOther+0.2*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[5]*uOther[16]*mnuOther+0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther+0.223606797749979*m0rOther[6]*uOther[15]*mnuOther+0.159719141249985*m0rOther[6]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[13]*mnuOther+0.2*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.25*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(56,41) = 0.1963961012123931*m0rOther[7]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[18]*mnuOther+0.21957751641342*m0rOther[2]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[17]*mnuOther+0.351382110749967*m0rOther[6]*uOther[17]*mnuOther+0.2*m0rOther[2]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[16]*mnuOther+0.351382110749967*m0rOther[7]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[16]*mnuOther+0.2*m0rOther[3]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[14]*mnuOther+0.2*m0rOther[5]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.447213595499958*m1rOther[13]*mnuOther+0.21957751641342*m0rOther[8]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(56,42) = 0.1963961012123931*m0rOther[6]*uOther[19]*mnuOther+0.149071198499986*m0rOther[8]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[18]*mnuOther+0.21957751641342*m0rOther[1]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[17]*mnuOther+0.351382110749967*m0rOther[7]*uOther[17]*mnuOther+0.2*m0rOther[1]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[16]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[16]*mnuOther+0.45*m0rOther[2]*uOther[16]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[14]*mnuOther-0.5000000000000001*m1rOther[14]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[13]*mnuOther+0.45*m0rOther[6]*uOther[12]*mnuOther+0.21957751641342*m0rOther[8]*uOther[11]*mnuOther+0.2*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(56,43) = 0.175662013130736*m0rOther[3]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[18]*mnuOther+0.21957751641342*m0rOther[0]*uOther[18]*mnuOther-0.43915503282684*m1rOther[18]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[17]*mnuOther+0.2*m0rOther[0]*uOther[17]*mnuOther-0.4*m1rOther[17]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[15]*mnuOther+0.2*m0rOther[1]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[14]*mnuOther+0.175662013130736*m0rOther[9]*uOther[13]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[12]*mnuOther+0.2*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.447213595499958*m1rOther[11]*mnuOther+0.21957751641342*m0rOther[8]*uOther[10]*mnuOther+0.2*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(56,44) = 0.2195775164134199*m0rOther[5]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[18]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[17]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[16]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[16]*mnuOther+0.159719141249985*m0rOther[0]*uOther[16]*mnuOther-0.31943828249997*m1rOther[16]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[14]*mnuOther+0.159719141249985*m0rOther[2]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.159719141249985*m0rOther[4]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[12]*mnuOther-0.5000000000000001*m1rOther[12]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[11]*mnuOther+0.159719141249985*m0rOther[6]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(56,45) = 0.21957751641342*m0rOther[4]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[18]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[16]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[16]*mnuOther+0.223606797749979*m0rOther[0]*uOther[16]*mnuOther-0.4472135954999579*m1rOther[16]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[15]*mnuOther+0.21957751641342*m0rOther[9]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[13]*mnuOther+0.2*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[12]*mnuOther+0.2*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(56,46) = 0.3833333333333334*m0rOther[9]*uOther[19]*mnuOther+0.1254728665219542*m0rOther[6]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[19]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[18]*mnuOther+0.29277002188456*m0rOther[7]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[18]*mnuOther+0.29277002188456*m0rOther[8]*uOther[17]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[17]*mnuOther+0.351382110749967*m0rOther[1]*uOther[17]*mnuOther+0.1254728665219542*m0rOther[9]*uOther[16]*mnuOther+0.9642857142857143*m0rOther[6]*uOther[16]*mnuOther+0.2874944542499729*m0rOther[2]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[14]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[14]*mnuOther+0.159719141249985*m0rOther[0]*uOther[14]*mnuOther-0.31943828249997*m1rOther[14]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[12]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[12]*mnuOther+0.45*m0rOther[2]*uOther[12]*mnuOther+0.3273268353539885*m0rOther[8]*uOther[11]*mnuOther+0.351382110749967*m0rOther[7]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(56,47) = 0.1928571428571429*m0rOther[8]*uOther[19]*mnuOther+0.29277002188456*m0rOther[7]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[19]*mnuOther+0.1928571428571429*m0rOther[9]*uOther[18]*mnuOther+0.29277002188456*m0rOther[6]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[18]*mnuOther+0.29277002188456*m0rOther[9]*uOther[17]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[17]*mnuOther+0.351382110749967*m0rOther[2]*uOther[17]*mnuOther+0.29277002188456*m0rOther[8]*uOther[16]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[16]*mnuOther+0.351382110749967*m0rOther[1]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[13]*mnuOther+0.2*m0rOther[0]*uOther[13]*mnuOther-0.4*m1rOther[13]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[12]*mnuOther+0.351382110749967*m0rOther[7]*uOther[12]*mnuOther+0.2*m0rOther[1]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[11]*mnuOther+0.351382110749967*m0rOther[6]*uOther[11]*mnuOther+0.2*m0rOther[2]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(56,48) = 0.1928571428571429*m0rOther[7]*uOther[19]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[18]*mnuOther+0.149071198499986*m0rOther[2]*uOther[18]*mnuOther+0.1928571428571429*m0rOther[9]*uOther[17]*mnuOther+0.29277002188456*m0rOther[6]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[17]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[16]*mnuOther+0.29277002188456*m0rOther[7]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[13]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[13]*mnuOther+0.21957751641342*m0rOther[0]*uOther[13]*mnuOther-0.43915503282684*m1rOther[13]*mnuOther+0.149071198499986*m0rOther[8]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[12]*mnuOther+0.21957751641342*m0rOther[1]*uOther[12]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[11]*mnuOther+0.21957751641342*m0rOther[2]*uOther[11]*mnuOther+0.21957751641342*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(56,49) = 0.3833333333333334*m0rOther[6]*uOther[19]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[18]*mnuOther+0.1928571428571429*m0rOther[8]*uOther[17]*mnuOther+0.29277002188456*m0rOther[7]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[17]*mnuOther+0.3833333333333334*m0rOther[9]*uOther[16]*mnuOther+0.1254728665219543*m0rOther[6]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[16]*mnuOther+0.21957751641342*m0rOther[4]*uOther[15]*mnuOther+0.21957751641342*m0rOther[5]*uOther[14]*mnuOther+0.175662013130736*m0rOther[3]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[11]*mnuOther; 
  data->AEM_S(57,40) = 0.2195775164134199*m0rOther[3]*uOther[19]*mnuOther+0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther+0.159719141249985*m0rOther[7]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[7]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[9]*uOther[13]*mnuOther+0.2*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(57,41) = 0.149071198499986*m0rOther[9]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[19]*mnuOther+0.21957751641342*m0rOther[2]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[17]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[17]*mnuOther+0.45*m0rOther[1]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[16]*mnuOther+0.351382110749967*m0rOther[6]*uOther[16]*mnuOther+0.2*m0rOther[2]*uOther[16]*mnuOther+0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[15]*mnuOther-0.5000000000000001*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[13]*mnuOther+0.21957751641342*m0rOther[9]*uOther[12]*mnuOther+0.2*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.45*m0rOther[7]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(57,42) = 0.3273268353539885*m0rOther[7]*uOther[19]*mnuOther+0.21957751641342*m0rOther[1]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[18]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[17]*mnuOther+0.351382110749967*m0rOther[6]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[16]*mnuOther+0.351382110749967*m0rOther[7]*uOther[16]*mnuOther+0.2*m0rOther[1]*uOther[16]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[13]*mnuOther+0.2*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.447213595499958*m1rOther[13]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.21957751641342*m0rOther[9]*uOther[11]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(57,43) = 0.3273268353539885*m0rOther[5]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[19]*mnuOther+0.21957751641342*m0rOther[0]*uOther[19]*mnuOther-0.43915503282684*m1rOther[19]*mnuOther+0.175662013130736*m0rOther[3]*uOther[18]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[17]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[16]*mnuOther+0.2*m0rOther[0]*uOther[16]*mnuOther-0.4*m1rOther[16]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[15]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[14]*mnuOther+0.2*m0rOther[2]*uOther[14]*mnuOther+0.175662013130736*m0rOther[8]*uOther[13]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[13]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[12]*mnuOther+0.2*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[0]*uOther[12]*mnuOther-0.447213595499958*m1rOther[12]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[11]*mnuOther+0.21957751641342*m0rOther[9]*uOther[10]*mnuOther+0.2*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(57,44) = 0.1963961012123931*m0rOther[3]*uOther[19]*mnuOther+0.21957751641342*m0rOther[5]*uOther[18]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.4472135954999579*m1rOther[17]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[16]*mnuOther+0.21957751641342*m0rOther[8]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[13]*mnuOther+0.2*m0rOther[2]*uOther[13]*mnuOther+0.2*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(57,45) = 0.3273268353539885*m0rOther[3]*uOther[19]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[18]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[17]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[17]*mnuOther+0.159719141249985*m0rOther[0]*uOther[17]*mnuOther-0.31943828249997*m1rOther[17]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[16]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[15]*mnuOther+0.159719141249985*m0rOther[1]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[8]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[12]*mnuOther+0.159719141249985*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[11]*mnuOther-0.5000000000000001*m1rOther[11]*mnuOther+0.159719141249985*m0rOther[7]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(57,46) = 0.1928571428571429*m0rOther[8]*uOther[19]*mnuOther+0.29277002188456*m0rOther[7]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[19]*mnuOther+0.1928571428571429*m0rOther[9]*uOther[18]*mnuOther+0.29277002188456*m0rOther[6]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[18]*mnuOther+0.29277002188456*m0rOther[9]*uOther[17]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[17]*mnuOther+0.351382110749967*m0rOther[2]*uOther[17]*mnuOther+0.29277002188456*m0rOther[8]*uOther[16]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[16]*mnuOther+0.351382110749967*m0rOther[1]*uOther[16]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[13]*mnuOther+0.2*m0rOther[0]*uOther[13]*mnuOther-0.4*m1rOther[13]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[12]*mnuOther+0.351382110749967*m0rOther[7]*uOther[12]*mnuOther+0.2*m0rOther[1]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[9]*uOther[11]*mnuOther+0.351382110749967*m0rOther[6]*uOther[11]*mnuOther+0.2*m0rOther[2]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(57,47) = 0.4621212121212121*m0rOther[9]*uOther[19]*mnuOther+0.29277002188456*m0rOther[6]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[19]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[18]*mnuOther+0.1254728665219542*m0rOther[7]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[18]*mnuOther+0.1254728665219542*m0rOther[8]*uOther[17]*mnuOther+0.9642857142857143*m0rOther[7]*uOther[17]*mnuOther+0.2874944542499729*m0rOther[1]*uOther[17]*mnuOther+0.29277002188456*m0rOther[9]*uOther[16]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[16]*mnuOther+0.351382110749967*m0rOther[2]*uOther[16]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[15]*mnuOther+0.159719141249985*m0rOther[0]*uOther[15]*mnuOther-0.31943828249997*m1rOther[15]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.4472135954999579*m1rOther[14]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[13]*mnuOther+0.3273268353539885*m0rOther[9]*uOther[12]*mnuOther+0.351382110749967*m0rOther[6]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[8]*uOther[11]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[11]*mnuOther+0.45*m0rOther[1]*uOther[11]*mnuOther+0.159719141249985*m0rOther[5]*uOther[10]*mnuOther+0.223606797749979*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(57,48) = 0.1928571428571429*m0rOther[6]*uOther[19]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[18]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[17]*mnuOther+0.1254728665219543*m0rOther[7]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[17]*mnuOther+0.1928571428571429*m0rOther[9]*uOther[16]*mnuOther+0.29277002188456*m0rOther[6]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[16]*mnuOther+0.21957751641342*m0rOther[4]*uOther[15]*mnuOther+0.21957751641342*m0rOther[5]*uOther[14]*mnuOther+0.175662013130736*m0rOther[3]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[11]*mnuOther; 
  data->AEM_S(57,49) = 0.4621212121212121*m0rOther[7]*uOther[19]*mnuOther+0.149071198499986*m0rOther[1]*uOther[19]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[18]*mnuOther+0.4621212121212121*m0rOther[9]*uOther[17]*mnuOther+0.29277002188456*m0rOther[6]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[17]*mnuOther+0.1928571428571429*m0rOther[8]*uOther[16]*mnuOther+0.29277002188456*m0rOther[7]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[15]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[13]*mnuOther+0.21957751641342*m0rOther[0]*uOther[13]*mnuOther-0.43915503282684*m1rOther[13]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[12]*mnuOther+0.21957751641342*m0rOther[1]*uOther[12]*mnuOther+0.149071198499986*m0rOther[9]*uOther[11]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[11]*mnuOther+0.21957751641342*m0rOther[2]*uOther[11]*mnuOther+0.21957751641342*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(58,40) = 0.149071198499986*m0rOther[4]*uOther[18]*mnuOther+0.25*m0rOther[0]*uOther[18]*mnuOther-0.5*m1rOther[18]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[16]*mnuOther+0.149071198499986*m0rOther[8]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[8]*uOther[10]*mnuOther; 
  data->AEM_S(58,41) = 0.2499586742703185*m0rOther[8]*uOther[18]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[14]*mnuOther-0.4391550328268398*m1rOther[14]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[12]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[11]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(58,42) = 0.149071198499986*m0rOther[6]*uOther[18]*mnuOther+0.25*m0rOther[2]*uOther[18]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[17]*mnuOther+0.149071198499986*m0rOther[8]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[8]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(58,43) = 0.3833333333333334*m0rOther[3]*uOther[18]*mnuOther+0.175662013130736*m0rOther[3]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[16]*mnuOther+0.21957751641342*m0rOther[0]*uOther[16]*mnuOther-0.43915503282684*m1rOther[16]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[14]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[14]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[13]*mnuOther+0.175662013130736*m0rOther[7]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[11]*mnuOther+0.21957751641342*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(58,44) = 0.4621212121212121*m0rOther[4]*uOther[18]*mnuOther+0.149071198499986*m0rOther[0]*uOther[18]*mnuOther-0.2981423969999719*m1rOther[18]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[15]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[12]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[11]*mnuOther-0.4391550328268398*m1rOther[11]*mnuOther+0.149071198499986*m0rOther[8]*uOther[10]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(58,45) = 0.25*m0rOther[5]*uOther[18]*mnuOther+0.21957751641342*m0rOther[4]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[16]*mnuOther+0.25*m0rOther[8]*uOther[15]*mnuOther+0.21957751641342*m0rOther[7]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[13]*mnuOther; 
  data->AEM_S(58,46) = 0.1928571428571429*m0rOther[7]*uOther[19]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[18]*mnuOther+0.149071198499986*m0rOther[2]*uOther[18]*mnuOther+0.1928571428571429*m0rOther[9]*uOther[17]*mnuOther+0.29277002188456*m0rOther[6]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[17]*mnuOther+0.4621212121212121*m0rOther[8]*uOther[16]*mnuOther+0.29277002188456*m0rOther[7]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[1]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[5]*uOther[13]*mnuOther+0.3273268353539885*m0rOther[4]*uOther[13]*mnuOther+0.21957751641342*m0rOther[0]*uOther[13]*mnuOther-0.43915503282684*m1rOther[13]*mnuOther+0.149071198499986*m0rOther[8]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[12]*mnuOther+0.21957751641342*m0rOther[1]*uOther[12]*mnuOther+0.3273268353539885*m0rOther[6]*uOther[11]*mnuOther+0.21957751641342*m0rOther[2]*uOther[11]*mnuOther+0.21957751641342*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(58,47) = 0.1928571428571429*m0rOther[6]*uOther[19]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[18]*mnuOther+0.3833333333333334*m0rOther[8]*uOther[17]*mnuOther+0.1254728665219543*m0rOther[7]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[17]*mnuOther+0.1928571428571429*m0rOther[9]*uOther[16]*mnuOther+0.29277002188456*m0rOther[6]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[16]*mnuOther+0.21957751641342*m0rOther[4]*uOther[15]*mnuOther+0.21957751641342*m0rOther[5]*uOther[14]*mnuOther+0.175662013130736*m0rOther[3]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[11]*mnuOther; 
  data->AEM_S(58,48) = 0.25*m0rOther[9]*uOther[19]*mnuOther+0.5898601398601399*m0rOther[8]*uOther[18]*mnuOther+0.2499586742703185*m0rOther[1]*uOther[18]*mnuOther+0.3833333333333334*m0rOther[7]*uOther[17]*mnuOther+0.4621212121212121*m0rOther[6]*uOther[16]*mnuOther+0.149071198499986*m0rOther[2]*uOther[16]*mnuOther+0.25*m0rOther[5]*uOther[15]*mnuOther+0.4621212121212121*m0rOther[4]*uOther[14]*mnuOther+0.149071198499986*m0rOther[0]*uOther[14]*mnuOther-0.2981423969999719*m1rOther[14]*mnuOther+0.3833333333333334*m0rOther[3]*uOther[13]*mnuOther+0.149071198499986*m0rOther[6]*uOther[12]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther+0.2499586742703185*m0rOther[8]*uOther[11]*mnuOther+0.3833333333333334*m0rOther[1]*uOther[11]*mnuOther+0.149071198499986*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(58,49) = 0.25*m0rOther[8]*uOther[19]*mnuOther+0.25*m0rOther[9]*uOther[18]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[17]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(59,40) = 0.149071198499986*m0rOther[5]*uOther[19]*mnuOther+0.25*m0rOther[0]*uOther[19]*mnuOther-0.5*m1rOther[19]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[17]*mnuOther+0.149071198499986*m0rOther[9]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[12]*mnuOther+0.25*m0rOther[9]*uOther[10]*mnuOther; 
  data->AEM_S(59,41) = 0.149071198499986*m0rOther[7]*uOther[19]*mnuOther+0.25*m0rOther[1]*uOther[19]*mnuOther+0.149071198499986*m0rOther[9]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[16]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[12]*mnuOther+0.25*m0rOther[9]*uOther[11]*mnuOther; 
  data->AEM_S(59,42) = 0.2499586742703185*m0rOther[9]*uOther[19]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[15]*mnuOther-0.4391550328268398*m1rOther[15]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[13]*mnuOther+0.3833333333333334*m0rOther[9]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[7]*uOther[11]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(59,43) = 0.3833333333333334*m0rOther[3]*uOther[19]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[17]*mnuOther+0.21957751641342*m0rOther[0]*uOther[17]*mnuOther-0.43915503282684*m1rOther[17]*mnuOther+0.175662013130736*m0rOther[3]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[15]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[14]*mnuOther+0.3833333333333334*m0rOther[9]*uOther[13]*mnuOther+0.175662013130736*m0rOther[6]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[5]*uOther[11]*mnuOther+0.21957751641342*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(59,44) = 0.25*m0rOther[4]*uOther[19]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[17]*mnuOther+0.21957751641342*m0rOther[5]*uOther[16]*mnuOther+0.21957751641342*m0rOther[6]*uOther[15]*mnuOther+0.25*m0rOther[9]*uOther[14]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[13]*mnuOther; 
  data->AEM_S(59,45) = 0.4621212121212121*m0rOther[5]*uOther[19]*mnuOther+0.149071198499986*m0rOther[0]*uOther[19]*mnuOther-0.2981423969999719*m1rOther[19]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[17]*mnuOther+0.2195775164134199*m0rOther[4]*uOther[16]*mnuOther+0.4621212121212121*m0rOther[9]*uOther[15]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[15]*mnuOther+0.2195775164134199*m0rOther[6]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[13]*mnuOther+0.2195775164134199*m0rOther[1]*uOther[13]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[12]*mnuOther+0.2195775164134199*m0rOther[0]*uOther[12]*mnuOther-0.4391550328268398*m1rOther[12]*mnuOther+0.2195775164134199*m0rOther[3]*uOther[11]*mnuOther+0.149071198499986*m0rOther[9]*uOther[10]*mnuOther+0.2195775164134199*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(59,46) = 0.3833333333333334*m0rOther[6]*uOther[19]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[18]*mnuOther+0.1928571428571429*m0rOther[8]*uOther[17]*mnuOther+0.29277002188456*m0rOther[7]*uOther[17]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[17]*mnuOther+0.3833333333333334*m0rOther[9]*uOther[16]*mnuOther+0.1254728665219543*m0rOther[6]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[2]*uOther[16]*mnuOther+0.21957751641342*m0rOther[4]*uOther[15]*mnuOther+0.21957751641342*m0rOther[5]*uOther[14]*mnuOther+0.175662013130736*m0rOther[3]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[12]*mnuOther+0.1963961012123931*m0rOther[7]*uOther[11]*mnuOther; 
  data->AEM_S(59,47) = 0.4621212121212121*m0rOther[7]*uOther[19]*mnuOther+0.149071198499986*m0rOther[1]*uOther[19]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[18]*mnuOther+0.4621212121212121*m0rOther[9]*uOther[17]*mnuOther+0.29277002188456*m0rOther[6]*uOther[17]*mnuOther+0.3273268353539885*m0rOther[2]*uOther[17]*mnuOther+0.1928571428571429*m0rOther[8]*uOther[16]*mnuOther+0.29277002188456*m0rOther[7]*uOther[16]*mnuOther+0.1963961012123931*m0rOther[1]*uOther[16]*mnuOther+0.3273268353539885*m0rOther[3]*uOther[15]*mnuOther+0.1963961012123931*m0rOther[3]*uOther[14]*mnuOther+0.3273268353539885*m0rOther[5]*uOther[13]*mnuOther+0.1963961012123931*m0rOther[4]*uOther[13]*mnuOther+0.21957751641342*m0rOther[0]*uOther[13]*mnuOther-0.43915503282684*m1rOther[13]*mnuOther+0.3273268353539885*m0rOther[7]*uOther[12]*mnuOther+0.21957751641342*m0rOther[1]*uOther[12]*mnuOther+0.149071198499986*m0rOther[9]*uOther[11]*mnuOther+0.1963961012123931*m0rOther[6]*uOther[11]*mnuOther+0.21957751641342*m0rOther[2]*uOther[11]*mnuOther+0.21957751641342*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(59,48) = 0.25*m0rOther[8]*uOther[19]*mnuOther+0.25*m0rOther[9]*uOther[18]*mnuOther+0.1928571428571429*m0rOther[6]*uOther[17]*mnuOther+0.1928571428571429*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(59,49) = 0.5898601398601399*m0rOther[9]*uOther[19]*mnuOther+0.2499586742703185*m0rOther[2]*uOther[19]*mnuOther+0.25*m0rOther[8]*uOther[18]*mnuOther+0.4621212121212121*m0rOther[7]*uOther[17]*mnuOther+0.149071198499986*m0rOther[1]*uOther[17]*mnuOther+0.3833333333333334*m0rOther[6]*uOther[16]*mnuOther+0.4621212121212121*m0rOther[5]*uOther[15]*mnuOther+0.149071198499986*m0rOther[0]*uOther[15]*mnuOther-0.2981423969999719*m1rOther[15]*mnuOther+0.25*m0rOther[4]*uOther[14]*mnuOther+0.3833333333333334*m0rOther[3]*uOther[13]*mnuOther+0.2499586742703185*m0rOther[9]*uOther[12]*mnuOther+0.3833333333333334*m0rOther[2]*uOther[12]*mnuOther+0.149071198499986*m0rOther[7]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[11]*mnuOther+0.149071198499986*m0rOther[5]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[10] += (m1rOther[10]-1.0*m1rSelf[10])*betaGreenep1*mnuSelf+m1rSelf[10]*mnuSelf-1.0*m1rOther[10]*mnuOther; 
  m1Relax[11] += (m1rOther[11]-1.0*m1rSelf[11])*betaGreenep1*mnuSelf+m1rSelf[11]*mnuSelf-1.0*m1rOther[11]*mnuOther; 
  m1Relax[12] += (m1rOther[12]-1.0*m1rSelf[12])*betaGreenep1*mnuSelf+m1rSelf[12]*mnuSelf-1.0*m1rOther[12]*mnuOther; 
  m1Relax[13] += (m1rOther[13]-1.0*m1rSelf[13])*betaGreenep1*mnuSelf+m1rSelf[13]*mnuSelf-1.0*m1rOther[13]*mnuOther; 
  m1Relax[14] += (m1rOther[14]-1.0*m1rSelf[14])*betaGreenep1*mnuSelf+m1rSelf[14]*mnuSelf-1.0*m1rOther[14]*mnuOther; 
  m1Relax[15] += (m1rOther[15]-1.0*m1rSelf[15])*betaGreenep1*mnuSelf+m1rSelf[15]*mnuSelf-1.0*m1rOther[15]*mnuOther; 
  m1Relax[16] += (m1rOther[16]-1.0*m1rSelf[16])*betaGreenep1*mnuSelf+m1rSelf[16]*mnuSelf-1.0*m1rOther[16]*mnuOther; 
  m1Relax[17] += (m1rOther[17]-1.0*m1rSelf[17])*betaGreenep1*mnuSelf+m1rSelf[17]*mnuSelf-1.0*m1rOther[17]*mnuOther; 
  m1Relax[18] += (m1rOther[18]-1.0*m1rSelf[18])*betaGreenep1*mnuSelf+m1rSelf[18]*mnuSelf-1.0*m1rOther[18]*mnuOther; 
  m1Relax[19] += (m1rOther[19]-1.0*m1rSelf[19])*betaGreenep1*mnuSelf+m1rSelf[19]*mnuSelf-1.0*m1rOther[19]*mnuOther; 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(50,20) = 0.25*cMSelf[19]*uSelf[19]*mnuSelf+0.25*cMSelf[18]*uSelf[18]*mnuSelf+0.25*cMSelf[17]*uSelf[17]*mnuSelf+0.25*cMSelf[16]*uSelf[16]*mnuSelf+0.25*cMSelf[15]*uSelf[15]*mnuSelf+0.25*cMSelf[14]*uSelf[14]*mnuSelf+0.25*cMSelf[13]*uSelf[13]*mnuSelf+0.25*cMSelf[12]*uSelf[12]*mnuSelf+0.25*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(50,21) = 0.2195775164134199*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[18]*mnuSelf+0.2500000000000001*cMSelf[15]*uSelf[17]*mnuSelf+0.2500000000000001*uSelf[15]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[14]*mnuSelf+0.25*cMSelf[12]*uSelf[13]*mnuSelf+0.25*uSelf[12]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[11]*mnuSelf+0.25*uSelf[10]*cMSelf[11]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(50,22) = 0.2195775164134199*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[19]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[17]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[17]*mnuSelf+0.2500000000000001*cMSelf[14]*uSelf[16]*mnuSelf+0.2500000000000001*uSelf[14]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[15]*mnuSelf+0.25*cMSelf[11]*uSelf[13]*mnuSelf+0.25*uSelf[11]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[12]*mnuSelf+0.25*uSelf[10]*cMSelf[12]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(50,23) = 0.2195775164134199*cMSelf[17]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[17]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[16]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[16]*cMSelf[18]*mnuSelf+0.2*cMSelf[16]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[17]*mnuSelf+0.2*uSelf[16]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[14]*mnuSelf+0.25*cMSelf[10]*uSelf[13]*mnuSelf+0.25*uSelf[10]*cMSelf[13]*mnuSelf+0.25*cMSelf[11]*uSelf[12]*mnuSelf+0.25*uSelf[11]*cMSelf[12]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(50,24) = 0.149071198499986*cMSelf[18]*uSelf[18]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[18]*mnuSelf+0.223606797749979*cMSelf[17]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[16]*uSelf[16]*mnuSelf+0.2500000000000001*cMSelf[12]*uSelf[16]*mnuSelf+0.2500000000000001*uSelf[12]*cMSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[14]*mnuSelf+0.25*cMSelf[10]*uSelf[14]*mnuSelf+0.25*uSelf[10]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[11]*mnuSelf+m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(50,25) = 0.149071198499986*cMSelf[19]*uSelf[19]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[17]*uSelf[17]*mnuSelf+0.2500000000000001*cMSelf[11]*uSelf[17]*mnuSelf+0.2500000000000001*uSelf[11]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[16]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[15]*mnuSelf+0.25*cMSelf[10]*uSelf[15]*mnuSelf+0.25*uSelf[10]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[12]*mnuSelf+m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(50,26) = 0.2195775164134199*cMSelf[13]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[18]*mnuSelf+0.2*cMSelf[13]*uSelf[17]*mnuSelf+0.2*uSelf[13]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[15]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[16]*mnuSelf+0.25*cMSelf[10]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[15]*cMSelf[16]*mnuSelf+0.159719141249985*uSelf[14]*cMSelf[16]*mnuSelf+0.25*uSelf[10]*cMSelf[16]*mnuSelf+0.2500000000000001*cMSelf[12]*uSelf[14]*mnuSelf+0.2500000000000001*uSelf[12]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[13]*mnuSelf+m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(50,27) = 0.2195775164134199*cMSelf[13]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[14]*uSelf[17]*mnuSelf+0.25*cMSelf[10]*uSelf[17]*mnuSelf+0.159719141249985*uSelf[15]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[14]*cMSelf[17]*mnuSelf+0.25*uSelf[10]*cMSelf[17]*mnuSelf+0.2*cMSelf[13]*uSelf[16]*mnuSelf+0.2*uSelf[13]*cMSelf[16]*mnuSelf+0.2500000000000001*cMSelf[11]*uSelf[15]*mnuSelf+0.2500000000000001*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[13]*mnuSelf+m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(50,28) = 0.149071198499986*cMSelf[14]*uSelf[18]*mnuSelf+0.25*cMSelf[10]*uSelf[18]*mnuSelf+0.149071198499986*uSelf[14]*cMSelf[18]*mnuSelf+0.25*uSelf[10]*cMSelf[18]*mnuSelf+0.2195775164134199*cMSelf[13]*uSelf[16]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[16]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[14]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[14]*mnuSelf+m0rSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf; 
  data->AEM_S(50,29) = 0.149071198499986*cMSelf[15]*uSelf[19]*mnuSelf+0.25*cMSelf[10]*uSelf[19]*mnuSelf+0.149071198499986*uSelf[15]*cMSelf[19]*mnuSelf+0.25*uSelf[10]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[13]*uSelf[17]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[17]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[15]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[15]*mnuSelf+m0rSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf; 
  data->AEM_S(51,20) = 0.2195775164134199*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[18]*mnuSelf+0.2500000000000001*cMSelf[15]*uSelf[17]*mnuSelf+0.2500000000000001*uSelf[15]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[14]*mnuSelf+0.25*cMSelf[12]*uSelf[13]*mnuSelf+0.25*uSelf[12]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[11]*mnuSelf+0.25*uSelf[10]*cMSelf[11]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(51,21) = 0.25*cMSelf[19]*uSelf[19]*mnuSelf+0.3833333333333334*cMSelf[18]*uSelf[18]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[18]*mnuSelf+0.45*cMSelf[17]*uSelf[17]*mnuSelf+0.3928571428571428*cMSelf[16]*uSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[16]*mnuSelf+0.25*cMSelf[15]*uSelf[15]*mnuSelf+0.3928571428571428*cMSelf[14]*uSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[14]*mnuSelf+0.45*cMSelf[13]*uSelf[13]*mnuSelf+0.25*cMSelf[12]*uSelf[12]*mnuSelf+0.45*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(51,22) = 0.2195775164134199*cMSelf[17]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[17]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[16]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[16]*cMSelf[18]*mnuSelf+0.2*cMSelf[16]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[17]*mnuSelf+0.2*uSelf[16]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[14]*mnuSelf+0.25*cMSelf[10]*uSelf[13]*mnuSelf+0.25*uSelf[10]*cMSelf[13]*mnuSelf+0.25*cMSelf[11]*uSelf[12]*mnuSelf+0.25*uSelf[11]*cMSelf[12]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(51,23) = 0.2195775164134199*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[18]*mnuSelf+0.4024922359499621*cMSelf[13]*uSelf[17]*mnuSelf+0.4024922359499621*uSelf[13]*cMSelf[17]*mnuSelf+0.2*cMSelf[15]*uSelf[16]*mnuSelf+0.3928571428571429*cMSelf[14]*uSelf[16]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[16]*mnuSelf+0.2*uSelf[15]*cMSelf[16]*mnuSelf+0.3928571428571429*uSelf[14]*cMSelf[16]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[14]*mnuSelf+0.45*cMSelf[11]*uSelf[13]*mnuSelf+0.45*uSelf[11]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[12]*mnuSelf+0.25*uSelf[10]*cMSelf[12]*mnuSelf+0.8944271909999161*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(51,24) = 0.3273268353539885*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[18]*mnuSelf+0.3273268353539885*uSelf[14]*cMSelf[18]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[18]*mnuSelf+0.223606797749979*cMSelf[15]*uSelf[17]*mnuSelf+0.223606797749979*uSelf[15]*cMSelf[17]*mnuSelf+0.3928571428571429*cMSelf[13]*uSelf[16]*mnuSelf+0.3928571428571429*uSelf[13]*cMSelf[16]*mnuSelf+0.3928571428571428*cMSelf[11]*uSelf[14]*mnuSelf+0.3928571428571428*uSelf[11]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[11]*mnuSelf+0.8783100656536796*m0rSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(51,25) = 0.2195775164134199*cMSelf[13]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[14]*uSelf[17]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[17]*mnuSelf+0.159719141249985*uSelf[15]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[14]*cMSelf[17]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[17]*mnuSelf+0.2*cMSelf[13]*uSelf[16]*mnuSelf+0.2*uSelf[13]*cMSelf[16]*mnuSelf+0.25*cMSelf[11]*uSelf[15]*mnuSelf+0.25*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[13]*mnuSelf+1.0*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(51,26) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999832*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999832*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(51,27) = 0.149071198499986*cMSelf[19]*uSelf[19]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[17]*uSelf[17]*mnuSelf+0.25*cMSelf[11]*uSelf[17]*mnuSelf+0.25*uSelf[11]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[16]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[15]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[15]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[12]*mnuSelf+1.0*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(51,28) = 0.1309307341415954*cMSelf[18]*uSelf[18]*mnuSelf+0.1928571428571429*cMSelf[11]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[11]*cMSelf[18]*mnuSelf+0.1963961012123931*cMSelf[17]*uSelf[17]*mnuSelf+0.1402829294374236*cMSelf[16]*uSelf[16]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[16]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[16]*mnuSelf+0.1402829294374236*cMSelf[14]*uSelf[14]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[14]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[14]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[13]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[11]*mnuSelf+0.8783100656536796*m0rSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf; 
  data->AEM_S(52,20) = 0.2195775164134199*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[19]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[17]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[17]*mnuSelf+0.2500000000000001*cMSelf[14]*uSelf[16]*mnuSelf+0.2500000000000001*uSelf[14]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[15]*mnuSelf+0.25*cMSelf[11]*uSelf[13]*mnuSelf+0.25*uSelf[11]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[12]*mnuSelf+0.25*uSelf[10]*cMSelf[12]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(52,21) = 0.2195775164134199*cMSelf[17]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[17]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[16]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[16]*cMSelf[18]*mnuSelf+0.2*cMSelf[16]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[17]*mnuSelf+0.2*uSelf[16]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[14]*mnuSelf+0.25*cMSelf[10]*uSelf[13]*mnuSelf+0.25*uSelf[10]*cMSelf[13]*mnuSelf+0.25*cMSelf[11]*uSelf[12]*mnuSelf+0.25*uSelf[11]*cMSelf[12]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(52,22) = 0.3833333333333334*cMSelf[19]*uSelf[19]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[19]*mnuSelf+0.25*cMSelf[18]*uSelf[18]*mnuSelf+0.3928571428571428*cMSelf[17]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[17]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[17]*mnuSelf+0.45*cMSelf[16]*uSelf[16]*mnuSelf+0.3928571428571428*cMSelf[15]*uSelf[15]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[15]*mnuSelf+0.25*cMSelf[14]*uSelf[14]*mnuSelf+0.45*cMSelf[13]*uSelf[13]*mnuSelf+0.45*cMSelf[12]*uSelf[12]*mnuSelf+0.25*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(52,23) = 0.1963961012123931*cMSelf[13]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[18]*mnuSelf+0.3928571428571429*cMSelf[15]*uSelf[17]*mnuSelf+0.2*cMSelf[14]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[17]*mnuSelf+0.3928571428571429*uSelf[15]*cMSelf[17]*mnuSelf+0.2*uSelf[14]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[17]*mnuSelf+0.4024922359499621*cMSelf[13]*uSelf[16]*mnuSelf+0.4024922359499621*uSelf[13]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[14]*mnuSelf+0.45*cMSelf[12]*uSelf[13]*mnuSelf+0.45*uSelf[12]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[11]*mnuSelf+0.25*uSelf[10]*cMSelf[11]*mnuSelf+0.8944271909999161*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(52,24) = 0.2195775164134199*cMSelf[13]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[18]*mnuSelf+0.2*cMSelf[13]*uSelf[17]*mnuSelf+0.2*uSelf[13]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[15]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[16]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[15]*cMSelf[16]*mnuSelf+0.159719141249985*uSelf[14]*cMSelf[16]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[16]*mnuSelf+0.25*cMSelf[12]*uSelf[14]*mnuSelf+0.25*uSelf[12]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[13]*mnuSelf+1.0*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(52,25) = 0.3273268353539885*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[19]*mnuSelf+0.3273268353539885*uSelf[15]*cMSelf[19]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[19]*mnuSelf+0.3928571428571429*cMSelf[13]*uSelf[17]*mnuSelf+0.3928571428571429*uSelf[13]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[14]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[14]*cMSelf[16]*mnuSelf+0.3928571428571428*cMSelf[12]*uSelf[15]*mnuSelf+0.3928571428571428*uSelf[12]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[12]*mnuSelf+0.8783100656536796*m0rSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(52,26) = 0.149071198499986*cMSelf[18]*uSelf[18]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[18]*mnuSelf+0.223606797749979*cMSelf[17]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[16]*uSelf[16]*mnuSelf+0.25*cMSelf[12]*uSelf[16]*mnuSelf+0.25*uSelf[12]*cMSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[14]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[14]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[11]*mnuSelf+1.0*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(52,27) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999832*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999832*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(52,29) = 0.1309307341415954*cMSelf[19]*uSelf[19]*mnuSelf+0.1928571428571429*cMSelf[12]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[12]*cMSelf[19]*mnuSelf+0.1402829294374236*cMSelf[17]*uSelf[17]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[17]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[16]*mnuSelf+0.1402829294374236*cMSelf[15]*uSelf[15]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[15]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[13]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[12]*mnuSelf+0.8783100656536796*m0rSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf; 
  data->AEM_S(53,20) = 0.2195775164134199*cMSelf[17]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[17]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[16]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[16]*cMSelf[18]*mnuSelf+0.2*cMSelf[16]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[17]*mnuSelf+0.2*uSelf[16]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[13]*cMSelf[14]*mnuSelf+0.25*cMSelf[10]*uSelf[13]*mnuSelf+0.25*uSelf[10]*cMSelf[13]*mnuSelf+0.25*cMSelf[11]*uSelf[12]*mnuSelf+0.25*uSelf[11]*cMSelf[12]*mnuSelf+m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(53,21) = 0.2195775164134199*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[18]*mnuSelf+0.4024922359499621*cMSelf[13]*uSelf[17]*mnuSelf+0.4024922359499621*uSelf[13]*cMSelf[17]*mnuSelf+0.2*cMSelf[15]*uSelf[16]*mnuSelf+0.3928571428571429*cMSelf[14]*uSelf[16]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[16]*mnuSelf+0.2*uSelf[15]*cMSelf[16]*mnuSelf+0.3928571428571429*uSelf[14]*cMSelf[16]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[14]*mnuSelf+0.45*cMSelf[11]*uSelf[13]*mnuSelf+0.45*uSelf[11]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[12]*mnuSelf+0.25*uSelf[10]*cMSelf[12]*mnuSelf+0.8944271909999161*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(53,22) = 0.1963961012123931*cMSelf[13]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[18]*mnuSelf+0.3928571428571429*cMSelf[15]*uSelf[17]*mnuSelf+0.2*cMSelf[14]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[17]*mnuSelf+0.3928571428571429*uSelf[15]*cMSelf[17]*mnuSelf+0.2*uSelf[14]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[17]*mnuSelf+0.4024922359499621*cMSelf[13]*uSelf[16]*mnuSelf+0.4024922359499621*uSelf[13]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[14]*mnuSelf+0.45*cMSelf[12]*uSelf[13]*mnuSelf+0.45*uSelf[12]*cMSelf[13]*mnuSelf+0.25*cMSelf[10]*uSelf[11]*mnuSelf+0.25*uSelf[10]*cMSelf[11]*mnuSelf+0.8944271909999161*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(53,23) = 0.3833333333333334*cMSelf[19]*uSelf[19]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[19]*mnuSelf+0.3833333333333334*cMSelf[18]*uSelf[18]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[18]*mnuSelf+0.5928571428571429*cMSelf[17]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[17]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[17]*mnuSelf+0.5928571428571429*cMSelf[16]*uSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[16]*mnuSelf+0.3928571428571428*cMSelf[15]*uSelf[15]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[15]*mnuSelf+0.3928571428571428*cMSelf[14]*uSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[14]*mnuSelf+0.65*cMSelf[13]*uSelf[13]*mnuSelf+0.45*cMSelf[12]*uSelf[12]*mnuSelf+0.45*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(53,24) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999831*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999831*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(53,25) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999831*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999831*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(53,26) = 0.175662013130736*cMSelf[13]*uSelf[19]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[19]*mnuSelf+0.3273268353539885*cMSelf[14]*uSelf[18]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[18]*mnuSelf+0.3273268353539885*uSelf[14]*cMSelf[18]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[18]*mnuSelf+0.3513821107499669*cMSelf[15]*uSelf[17]*mnuSelf+0.1788854381999831*cMSelf[14]*uSelf[17]*mnuSelf+0.2*cMSelf[10]*uSelf[17]*mnuSelf+0.3513821107499669*uSelf[15]*cMSelf[17]*mnuSelf+0.1788854381999831*uSelf[14]*cMSelf[17]*mnuSelf+0.2*uSelf[10]*cMSelf[17]*mnuSelf+0.5528571428571428*cMSelf[13]*uSelf[16]*mnuSelf+0.5528571428571428*uSelf[13]*cMSelf[16]*mnuSelf+0.2*cMSelf[11]*uSelf[15]*mnuSelf+0.2*uSelf[11]*cMSelf[15]*mnuSelf+0.3928571428571429*cMSelf[11]*uSelf[14]*mnuSelf+0.3928571428571429*uSelf[11]*cMSelf[14]*mnuSelf+0.4024922359499621*cMSelf[12]*uSelf[13]*mnuSelf+0.4024922359499621*uSelf[12]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[11]*mnuSelf+0.87831006565368*m0rSelf[8]*mnuSelf-0.43915503282684*cESelf[8]*mnuSelf+0.8*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.8944271909999161*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(53,27) = 0.3273268353539885*cMSelf[15]*uSelf[19]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[19]*mnuSelf+0.3273268353539885*uSelf[15]*cMSelf[19]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[19]*mnuSelf+0.175662013130736*cMSelf[13]*uSelf[18]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[18]*mnuSelf+0.5528571428571428*cMSelf[13]*uSelf[17]*mnuSelf+0.5528571428571428*uSelf[13]*cMSelf[17]*mnuSelf+0.1788854381999831*cMSelf[15]*uSelf[16]*mnuSelf+0.3513821107499669*cMSelf[14]*uSelf[16]*mnuSelf+0.2*cMSelf[10]*uSelf[16]*mnuSelf+0.1788854381999831*uSelf[15]*cMSelf[16]*mnuSelf+0.3513821107499669*uSelf[14]*cMSelf[16]*mnuSelf+0.2*uSelf[10]*cMSelf[16]*mnuSelf+0.3928571428571429*cMSelf[12]*uSelf[15]*mnuSelf+0.3928571428571429*uSelf[12]*cMSelf[15]*mnuSelf+0.2*cMSelf[12]*uSelf[14]*mnuSelf+0.2*uSelf[12]*cMSelf[14]*mnuSelf+0.4024922359499621*cMSelf[11]*uSelf[13]*mnuSelf+0.4024922359499621*uSelf[11]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[12]*mnuSelf+0.87831006565368*m0rSelf[9]*mnuSelf-0.43915503282684*cESelf[9]*mnuSelf+0.8*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.8944271909999161*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(53,28) = 0.1928571428571429*cMSelf[13]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[13]*cMSelf[18]*mnuSelf+0.175662013130736*cMSelf[13]*uSelf[17]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[15]*uSelf[16]*mnuSelf+0.1402829294374237*cMSelf[14]*uSelf[16]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[16]*mnuSelf+0.1963961012123931*uSelf[15]*cMSelf[16]*mnuSelf+0.1402829294374237*uSelf[14]*cMSelf[16]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[16]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[14]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[14]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[13]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[13]*mnuSelf+0.87831006565368*m0rSelf[6]*mnuSelf-0.43915503282684*cESelf[6]*mnuSelf; 
  data->AEM_S(53,29) = 0.1928571428571429*cMSelf[13]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[13]*cMSelf[19]*mnuSelf+0.1402829294374237*cMSelf[15]*uSelf[17]*mnuSelf+0.1963961012123931*cMSelf[14]*uSelf[17]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[17]*mnuSelf+0.1402829294374237*uSelf[15]*cMSelf[17]*mnuSelf+0.1963961012123931*uSelf[14]*cMSelf[17]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[17]*mnuSelf+0.175662013130736*cMSelf[13]*uSelf[16]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[16]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[15]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[13]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[13]*mnuSelf+0.87831006565368*m0rSelf[7]*mnuSelf-0.43915503282684*cESelf[7]*mnuSelf; 
  data->AEM_S(54,20) = 0.149071198499986*cMSelf[18]*uSelf[18]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[18]*mnuSelf+0.223606797749979*cMSelf[17]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[16]*uSelf[16]*mnuSelf+0.2500000000000001*cMSelf[12]*uSelf[16]*mnuSelf+0.2500000000000001*uSelf[12]*cMSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[14]*mnuSelf+0.25*cMSelf[10]*uSelf[14]*mnuSelf+0.25*uSelf[10]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[11]*mnuSelf+m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(54,21) = 0.3273268353539885*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[18]*mnuSelf+0.3273268353539885*uSelf[14]*cMSelf[18]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[18]*mnuSelf+0.223606797749979*cMSelf[15]*uSelf[17]*mnuSelf+0.223606797749979*uSelf[15]*cMSelf[17]*mnuSelf+0.3928571428571429*cMSelf[13]*uSelf[16]*mnuSelf+0.3928571428571429*uSelf[13]*cMSelf[16]*mnuSelf+0.3928571428571428*cMSelf[11]*uSelf[14]*mnuSelf+0.3928571428571428*uSelf[11]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[11]*mnuSelf+0.8783100656536796*m0rSelf[8]*mnuSelf-0.4391550328268398*cESelf[8]*mnuSelf+0.8944271909999159*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(54,22) = 0.2195775164134199*cMSelf[13]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[18]*mnuSelf+0.2*cMSelf[13]*uSelf[17]*mnuSelf+0.2*uSelf[13]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[15]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[16]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[15]*cMSelf[16]*mnuSelf+0.159719141249985*uSelf[14]*cMSelf[16]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[16]*mnuSelf+0.25*cMSelf[12]*uSelf[14]*mnuSelf+0.25*uSelf[12]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[13]*mnuSelf+1.0*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(54,23) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999831*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999831*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(54,24) = 0.25*cMSelf[19]*uSelf[19]*mnuSelf+0.3452380952380952*cMSelf[18]*uSelf[18]*mnuSelf+0.1402829294374236*cMSelf[11]*uSelf[18]*mnuSelf+0.1402829294374236*uSelf[11]*cMSelf[18]*mnuSelf+0.3928571428571428*cMSelf[17]*uSelf[17]*mnuSelf+0.3520408163265306*cMSelf[16]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[12]*uSelf[16]*mnuSelf+0.159719141249985*uSelf[12]*cMSelf[16]*mnuSelf+0.25*cMSelf[15]*uSelf[15]*mnuSelf+0.3520408163265306*cMSelf[14]*uSelf[14]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[14]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[14]*mnuSelf+0.3928571428571428*cMSelf[13]*uSelf[13]*mnuSelf+0.25*cMSelf[12]*uSelf[12]*mnuSelf+0.3928571428571428*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.6388765649999399*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(54,26) = 0.2195775164134199*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[19]*mnuSelf+0.1402829294374237*cMSelf[13]*uSelf[18]*mnuSelf+0.1402829294374237*uSelf[13]*cMSelf[18]*mnuSelf+0.3513821107499669*cMSelf[13]*uSelf[17]*mnuSelf+0.3513821107499669*uSelf[13]*cMSelf[17]*mnuSelf+0.1428571428571428*cMSelf[15]*uSelf[16]*mnuSelf+0.3520408163265306*cMSelf[14]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[16]*mnuSelf+0.1428571428571428*uSelf[15]*cMSelf[16]*mnuSelf+0.3520408163265306*uSelf[14]*cMSelf[16]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[15]*mnuSelf+0.159719141249985*cMSelf[12]*uSelf[14]*mnuSelf+0.159719141249985*uSelf[12]*cMSelf[14]*mnuSelf+0.3928571428571429*cMSelf[11]*uSelf[13]*mnuSelf+0.3928571428571429*uSelf[11]*cMSelf[13]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[12]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[12]*mnuSelf+0.6388765649999399*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.0*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(54,27) = 0.1963961012123931*cMSelf[13]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[19]*mnuSelf+0.1428571428571428*cMSelf[15]*uSelf[17]*mnuSelf+0.2*cMSelf[14]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[17]*mnuSelf+0.1428571428571428*uSelf[15]*cMSelf[17]*mnuSelf+0.2*uSelf[14]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[17]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[16]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[15]*mnuSelf+0.2*cMSelf[12]*uSelf[13]*mnuSelf+0.2*uSelf[12]*cMSelf[13]*mnuSelf+0.8944271909999159*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(54,28) = 0.2817460317460317*cMSelf[14]*uSelf[18]*mnuSelf+0.149071198499986*cMSelf[10]*uSelf[18]*mnuSelf+0.2817460317460317*uSelf[14]*cMSelf[18]*mnuSelf+0.149071198499986*uSelf[10]*cMSelf[18]*mnuSelf+0.2195775164134199*cMSelf[15]*uSelf[17]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[17]*mnuSelf+0.3273268353539885*cMSelf[13]*uSelf[16]*mnuSelf+0.3273268353539885*uSelf[13]*cMSelf[16]*mnuSelf+0.3273268353539885*cMSelf[11]*uSelf[14]*mnuSelf+0.3273268353539885*uSelf[11]*cMSelf[14]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[13]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[13]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[11]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[11]*mnuSelf+0.5962847939999438*m0rSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+0.8783100656536796*m0rSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(55,20) = 0.149071198499986*cMSelf[19]*uSelf[19]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[17]*uSelf[17]*mnuSelf+0.2500000000000001*cMSelf[11]*uSelf[17]*mnuSelf+0.2500000000000001*uSelf[11]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[16]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[15]*mnuSelf+0.25*cMSelf[10]*uSelf[15]*mnuSelf+0.25*uSelf[10]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[12]*mnuSelf+m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(55,21) = 0.2195775164134199*cMSelf[13]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[14]*uSelf[17]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[17]*mnuSelf+0.159719141249985*uSelf[15]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[14]*cMSelf[17]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[17]*mnuSelf+0.2*cMSelf[13]*uSelf[16]*mnuSelf+0.2*uSelf[13]*cMSelf[16]*mnuSelf+0.25*cMSelf[11]*uSelf[15]*mnuSelf+0.25*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[13]*mnuSelf+1.0*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(55,22) = 0.3273268353539885*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[19]*mnuSelf+0.3273268353539885*uSelf[15]*cMSelf[19]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[19]*mnuSelf+0.3928571428571429*cMSelf[13]*uSelf[17]*mnuSelf+0.3928571428571429*uSelf[13]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[14]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[14]*cMSelf[16]*mnuSelf+0.3928571428571428*cMSelf[12]*uSelf[15]*mnuSelf+0.3928571428571428*uSelf[12]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[12]*mnuSelf+0.8783100656536796*m0rSelf[9]*mnuSelf-0.4391550328268398*cESelf[9]*mnuSelf+0.8944271909999159*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(55,23) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999831*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999831*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999159*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(55,25) = 0.3452380952380952*cMSelf[19]*uSelf[19]*mnuSelf+0.1402829294374236*cMSelf[12]*uSelf[19]*mnuSelf+0.1402829294374236*uSelf[12]*cMSelf[19]*mnuSelf+0.25*cMSelf[18]*uSelf[18]*mnuSelf+0.3520408163265306*cMSelf[17]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[11]*uSelf[17]*mnuSelf+0.159719141249985*uSelf[11]*cMSelf[17]*mnuSelf+0.3928571428571428*cMSelf[16]*uSelf[16]*mnuSelf+0.3520408163265306*cMSelf[15]*uSelf[15]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[15]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[15]*mnuSelf+0.25*cMSelf[14]*uSelf[14]*mnuSelf+0.3928571428571428*cMSelf[13]*uSelf[13]*mnuSelf+0.3928571428571428*cMSelf[12]*uSelf[12]*mnuSelf+0.25*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.6388765649999399*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(55,26) = 0.1963961012123931*cMSelf[13]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[18]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[17]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[17]*mnuSelf+0.2*cMSelf[15]*uSelf[16]*mnuSelf+0.1428571428571428*cMSelf[14]*uSelf[16]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[16]*mnuSelf+0.2*uSelf[15]*cMSelf[16]*mnuSelf+0.1428571428571428*uSelf[14]*cMSelf[16]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[14]*mnuSelf+0.2*cMSelf[11]*uSelf[13]*mnuSelf+0.2*uSelf[11]*cMSelf[13]*mnuSelf+0.8944271909999159*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(55,27) = 0.1402829294374237*cMSelf[13]*uSelf[19]*mnuSelf+0.1402829294374237*uSelf[13]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[18]*mnuSelf+0.3520408163265306*cMSelf[15]*uSelf[17]*mnuSelf+0.1428571428571428*cMSelf[14]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[17]*mnuSelf+0.3520408163265306*uSelf[15]*cMSelf[17]*mnuSelf+0.1428571428571428*uSelf[14]*cMSelf[17]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[17]*mnuSelf+0.3513821107499669*cMSelf[13]*uSelf[16]*mnuSelf+0.3513821107499669*uSelf[13]*cMSelf[16]*mnuSelf+0.159719141249985*cMSelf[11]*uSelf[15]*mnuSelf+0.159719141249985*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[14]*mnuSelf+0.3928571428571429*cMSelf[12]*uSelf[13]*mnuSelf+0.3928571428571429*uSelf[12]*cMSelf[13]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[11]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[11]*mnuSelf+0.6388765649999399*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.0*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(55,29) = 0.2817460317460317*cMSelf[15]*uSelf[19]*mnuSelf+0.149071198499986*cMSelf[10]*uSelf[19]*mnuSelf+0.2817460317460317*uSelf[15]*cMSelf[19]*mnuSelf+0.149071198499986*uSelf[10]*cMSelf[19]*mnuSelf+0.3273268353539885*cMSelf[13]*uSelf[17]*mnuSelf+0.3273268353539885*uSelf[13]*cMSelf[17]*mnuSelf+0.2195775164134199*cMSelf[14]*uSelf[16]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[16]*mnuSelf+0.3273268353539885*cMSelf[12]*uSelf[15]*mnuSelf+0.3273268353539885*uSelf[12]*cMSelf[15]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[13]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[13]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[12]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[12]*mnuSelf+0.5962847939999438*m0rSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+0.8783100656536796*m0rSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(56,20) = 0.2195775164134199*cMSelf[13]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[18]*mnuSelf+0.2*cMSelf[13]*uSelf[17]*mnuSelf+0.2*uSelf[13]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[15]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[16]*mnuSelf+0.25*cMSelf[10]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[15]*cMSelf[16]*mnuSelf+0.159719141249985*uSelf[14]*cMSelf[16]*mnuSelf+0.25*uSelf[10]*cMSelf[16]*mnuSelf+0.2500000000000001*cMSelf[12]*uSelf[14]*mnuSelf+0.2500000000000001*uSelf[12]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[13]*mnuSelf+m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(56,21) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999832*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999832*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(56,22) = 0.149071198499986*cMSelf[18]*uSelf[18]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[18]*mnuSelf+0.223606797749979*cMSelf[17]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[16]*uSelf[16]*mnuSelf+0.25*cMSelf[12]*uSelf[16]*mnuSelf+0.25*uSelf[12]*cMSelf[16]*mnuSelf+0.159719141249985*cMSelf[14]*uSelf[14]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[14]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[11]*mnuSelf+1.0*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(56,23) = 0.175662013130736*cMSelf[13]*uSelf[19]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[19]*mnuSelf+0.3273268353539885*cMSelf[14]*uSelf[18]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[18]*mnuSelf+0.3273268353539885*uSelf[14]*cMSelf[18]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[18]*mnuSelf+0.3513821107499669*cMSelf[15]*uSelf[17]*mnuSelf+0.1788854381999831*cMSelf[14]*uSelf[17]*mnuSelf+0.2*cMSelf[10]*uSelf[17]*mnuSelf+0.3513821107499669*uSelf[15]*cMSelf[17]*mnuSelf+0.1788854381999831*uSelf[14]*cMSelf[17]*mnuSelf+0.2*uSelf[10]*cMSelf[17]*mnuSelf+0.5528571428571428*cMSelf[13]*uSelf[16]*mnuSelf+0.5528571428571428*uSelf[13]*cMSelf[16]*mnuSelf+0.2*cMSelf[11]*uSelf[15]*mnuSelf+0.2*uSelf[11]*cMSelf[15]*mnuSelf+0.3928571428571429*cMSelf[11]*uSelf[14]*mnuSelf+0.3928571428571429*uSelf[11]*cMSelf[14]*mnuSelf+0.4024922359499621*cMSelf[12]*uSelf[13]*mnuSelf+0.4024922359499621*uSelf[12]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[11]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[11]*mnuSelf+0.87831006565368*m0rSelf[8]*mnuSelf-0.43915503282684*cESelf[8]*mnuSelf+0.8*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.8944271909999161*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(56,24) = 0.2195775164134199*cMSelf[15]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[19]*mnuSelf+0.1402829294374237*cMSelf[13]*uSelf[18]*mnuSelf+0.1402829294374237*uSelf[13]*cMSelf[18]*mnuSelf+0.3513821107499669*cMSelf[13]*uSelf[17]*mnuSelf+0.3513821107499669*uSelf[13]*cMSelf[17]*mnuSelf+0.1428571428571428*cMSelf[15]*uSelf[16]*mnuSelf+0.3520408163265306*cMSelf[14]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[16]*mnuSelf+0.1428571428571428*uSelf[15]*cMSelf[16]*mnuSelf+0.3520408163265306*uSelf[14]*cMSelf[16]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[15]*mnuSelf+0.159719141249985*cMSelf[12]*uSelf[14]*mnuSelf+0.159719141249985*uSelf[12]*cMSelf[14]*mnuSelf+0.3928571428571429*cMSelf[11]*uSelf[13]*mnuSelf+0.3928571428571429*uSelf[11]*cMSelf[13]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[12]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[12]*mnuSelf+0.6388765649999399*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.0*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(56,25) = 0.1963961012123931*cMSelf[13]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[18]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[17]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[17]*mnuSelf+0.2*cMSelf[15]*uSelf[16]*mnuSelf+0.1428571428571428*cMSelf[14]*uSelf[16]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[16]*mnuSelf+0.2*uSelf[15]*cMSelf[16]*mnuSelf+0.1428571428571428*uSelf[14]*cMSelf[16]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[14]*mnuSelf+0.2*cMSelf[11]*uSelf[13]*mnuSelf+0.2*uSelf[11]*cMSelf[13]*mnuSelf+0.8944271909999159*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(56,26) = 0.3833333333333334*cMSelf[19]*uSelf[19]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[19]*mnuSelf+0.3452380952380952*cMSelf[18]*uSelf[18]*mnuSelf+0.1402829294374236*cMSelf[11]*uSelf[18]*mnuSelf+0.1402829294374236*uSelf[11]*cMSelf[18]*mnuSelf+0.5357142857142857*cMSelf[17]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[17]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[17]*mnuSelf+0.5520408163265306*cMSelf[16]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[12]*uSelf[16]*mnuSelf+0.159719141249985*uSelf[12]*cMSelf[16]*mnuSelf+0.3928571428571428*cMSelf[15]*uSelf[15]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[15]*mnuSelf+0.3520408163265306*cMSelf[14]*uSelf[14]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[14]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[14]*mnuSelf+0.5928571428571429*cMSelf[13]*uSelf[13]*mnuSelf+0.45*cMSelf[12]*uSelf[12]*mnuSelf+0.3928571428571428*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.8944271909999159*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.6388765649999399*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(56,27) = 0.175662013130736*cMSelf[17]*uSelf[19]*mnuSelf+0.175662013130736*uSelf[17]*cMSelf[19]*mnuSelf+0.175662013130736*cMSelf[16]*uSelf[18]*mnuSelf+0.175662013130736*uSelf[16]*cMSelf[18]*mnuSelf+0.16*cMSelf[16]*uSelf[17]*mnuSelf+0.1788854381999832*cMSelf[12]*uSelf[17]*mnuSelf+0.16*uSelf[16]*cMSelf[17]*mnuSelf+0.1788854381999832*uSelf[12]*cMSelf[17]*mnuSelf+0.1788854381999832*cMSelf[11]*uSelf[16]*mnuSelf+0.1788854381999832*uSelf[11]*cMSelf[16]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[15]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[15]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[14]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[14]*mnuSelf+0.2*cMSelf[10]*uSelf[13]*mnuSelf+0.2*uSelf[10]*cMSelf[13]*mnuSelf+0.2*cMSelf[11]*uSelf[12]*mnuSelf+0.2*uSelf[11]*cMSelf[12]*mnuSelf+0.8*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(56,28) = 0.1928571428571429*cMSelf[17]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[17]*cMSelf[19]*mnuSelf+0.1928571428571429*cMSelf[16]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[16]*cMSelf[18]*mnuSelf+0.175662013130736*cMSelf[16]*uSelf[17]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[17]*mnuSelf+0.175662013130736*uSelf[16]*cMSelf[17]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[16]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[16]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[15]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[14]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[14]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[13]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[13]*mnuSelf+0.21957751641342*cMSelf[11]*uSelf[12]*mnuSelf+0.21957751641342*uSelf[11]*cMSelf[12]*mnuSelf+0.87831006565368*m0rSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf; 
  data->AEM_S(57,20) = 0.2195775164134199*cMSelf[13]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[14]*uSelf[17]*mnuSelf+0.25*cMSelf[10]*uSelf[17]*mnuSelf+0.159719141249985*uSelf[15]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[14]*cMSelf[17]*mnuSelf+0.25*uSelf[10]*cMSelf[17]*mnuSelf+0.2*cMSelf[13]*uSelf[16]*mnuSelf+0.2*uSelf[13]*cMSelf[16]*mnuSelf+0.2500000000000001*cMSelf[11]*uSelf[15]*mnuSelf+0.2500000000000001*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[13]*mnuSelf+m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(57,21) = 0.149071198499986*cMSelf[19]*uSelf[19]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[19]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[19]*mnuSelf+0.159719141249985*cMSelf[17]*uSelf[17]*mnuSelf+0.25*cMSelf[11]*uSelf[17]*mnuSelf+0.25*uSelf[11]*cMSelf[17]*mnuSelf+0.223606797749979*cMSelf[16]*uSelf[16]*mnuSelf+0.159719141249985*cMSelf[15]*uSelf[15]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[15]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[13]*uSelf[13]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[12]*mnuSelf+1.0*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(57,22) = 0.1963961012123931*cMSelf[17]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[17]*cMSelf[19]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[16]*cMSelf[18]*mnuSelf+0.1788854381999832*cMSelf[16]*uSelf[17]*mnuSelf+0.2*cMSelf[12]*uSelf[17]*mnuSelf+0.1788854381999832*uSelf[16]*cMSelf[17]*mnuSelf+0.2*uSelf[12]*cMSelf[17]*mnuSelf+0.2*cMSelf[11]*uSelf[16]*mnuSelf+0.2*uSelf[11]*cMSelf[16]*mnuSelf+0.2*cMSelf[13]*uSelf[15]*mnuSelf+0.2*uSelf[13]*cMSelf[15]*mnuSelf+0.2*cMSelf[13]*uSelf[14]*mnuSelf+0.2*uSelf[13]*cMSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[13]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[12]*mnuSelf+0.8944271909999161*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(57,23) = 0.3273268353539885*cMSelf[15]*uSelf[19]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[19]*mnuSelf+0.3273268353539885*uSelf[15]*cMSelf[19]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[19]*mnuSelf+0.175662013130736*cMSelf[13]*uSelf[18]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[18]*mnuSelf+0.5528571428571428*cMSelf[13]*uSelf[17]*mnuSelf+0.5528571428571428*uSelf[13]*cMSelf[17]*mnuSelf+0.1788854381999831*cMSelf[15]*uSelf[16]*mnuSelf+0.3513821107499669*cMSelf[14]*uSelf[16]*mnuSelf+0.2*cMSelf[10]*uSelf[16]*mnuSelf+0.1788854381999831*uSelf[15]*cMSelf[16]*mnuSelf+0.3513821107499669*uSelf[14]*cMSelf[16]*mnuSelf+0.2*uSelf[10]*cMSelf[16]*mnuSelf+0.3928571428571429*cMSelf[12]*uSelf[15]*mnuSelf+0.3928571428571429*uSelf[12]*cMSelf[15]*mnuSelf+0.2*cMSelf[12]*uSelf[14]*mnuSelf+0.2*uSelf[12]*cMSelf[14]*mnuSelf+0.4024922359499621*cMSelf[11]*uSelf[13]*mnuSelf+0.4024922359499621*uSelf[11]*cMSelf[13]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[12]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[12]*mnuSelf+0.87831006565368*m0rSelf[9]*mnuSelf-0.43915503282684*cESelf[9]*mnuSelf+0.8*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.8944271909999161*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(57,24) = 0.1963961012123931*cMSelf[13]*uSelf[19]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[19]*mnuSelf+0.1428571428571428*cMSelf[15]*uSelf[17]*mnuSelf+0.2*cMSelf[14]*uSelf[17]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[17]*mnuSelf+0.1428571428571428*uSelf[15]*cMSelf[17]*mnuSelf+0.2*uSelf[14]*cMSelf[17]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[17]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[16]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[16]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[15]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[15]*mnuSelf+0.2*cMSelf[12]*uSelf[13]*mnuSelf+0.2*uSelf[12]*cMSelf[13]*mnuSelf+0.8944271909999159*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(57,25) = 0.1402829294374237*cMSelf[13]*uSelf[19]*mnuSelf+0.1402829294374237*uSelf[13]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[14]*uSelf[18]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[18]*mnuSelf+0.3520408163265306*cMSelf[15]*uSelf[17]*mnuSelf+0.1428571428571428*cMSelf[14]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[17]*mnuSelf+0.3520408163265306*uSelf[15]*cMSelf[17]*mnuSelf+0.1428571428571428*uSelf[14]*cMSelf[17]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[17]*mnuSelf+0.3513821107499669*cMSelf[13]*uSelf[16]*mnuSelf+0.3513821107499669*uSelf[13]*cMSelf[16]*mnuSelf+0.159719141249985*cMSelf[11]*uSelf[15]*mnuSelf+0.159719141249985*uSelf[11]*cMSelf[15]*mnuSelf+0.223606797749979*cMSelf[11]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[11]*cMSelf[14]*mnuSelf+0.3928571428571429*cMSelf[12]*uSelf[13]*mnuSelf+0.3928571428571429*uSelf[12]*cMSelf[13]*mnuSelf+0.2500000000000001*cMSelf[10]*uSelf[11]*mnuSelf+0.2500000000000001*uSelf[10]*cMSelf[11]*mnuSelf+0.6388765649999399*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.0*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(57,26) = 0.175662013130736*cMSelf[17]*uSelf[19]*mnuSelf+0.175662013130736*uSelf[17]*cMSelf[19]*mnuSelf+0.175662013130736*cMSelf[16]*uSelf[18]*mnuSelf+0.175662013130736*uSelf[16]*cMSelf[18]*mnuSelf+0.16*cMSelf[16]*uSelf[17]*mnuSelf+0.1788854381999832*cMSelf[12]*uSelf[17]*mnuSelf+0.16*uSelf[16]*cMSelf[17]*mnuSelf+0.1788854381999832*uSelf[12]*cMSelf[17]*mnuSelf+0.1788854381999832*cMSelf[11]*uSelf[16]*mnuSelf+0.1788854381999832*uSelf[11]*cMSelf[16]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[15]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[15]*mnuSelf+0.1788854381999831*cMSelf[13]*uSelf[14]*mnuSelf+0.1788854381999831*uSelf[13]*cMSelf[14]*mnuSelf+0.2*cMSelf[10]*uSelf[13]*mnuSelf+0.2*uSelf[10]*cMSelf[13]*mnuSelf+0.2*cMSelf[11]*uSelf[12]*mnuSelf+0.2*uSelf[11]*cMSelf[12]*mnuSelf+0.8*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(57,27) = 0.3452380952380952*cMSelf[19]*uSelf[19]*mnuSelf+0.1402829294374236*cMSelf[12]*uSelf[19]*mnuSelf+0.1402829294374236*uSelf[12]*cMSelf[19]*mnuSelf+0.3833333333333334*cMSelf[18]*uSelf[18]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[18]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[18]*mnuSelf+0.5520408163265306*cMSelf[17]*uSelf[17]*mnuSelf+0.159719141249985*cMSelf[11]*uSelf[17]*mnuSelf+0.159719141249985*uSelf[11]*cMSelf[17]*mnuSelf+0.5357142857142857*cMSelf[16]*uSelf[16]*mnuSelf+0.223606797749979*cMSelf[12]*uSelf[16]*mnuSelf+0.223606797749979*uSelf[12]*cMSelf[16]*mnuSelf+0.3520408163265306*cMSelf[15]*uSelf[15]*mnuSelf+0.159719141249985*cMSelf[10]*uSelf[15]*mnuSelf+0.159719141249985*uSelf[10]*cMSelf[15]*mnuSelf+0.3928571428571428*cMSelf[14]*uSelf[14]*mnuSelf+0.223606797749979*cMSelf[10]*uSelf[14]*mnuSelf+0.223606797749979*uSelf[10]*cMSelf[14]*mnuSelf+0.5928571428571429*cMSelf[13]*uSelf[13]*mnuSelf+0.3928571428571428*cMSelf[12]*uSelf[12]*mnuSelf+0.45*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.6388765649999399*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.8944271909999159*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(57,29) = 0.1928571428571429*cMSelf[17]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[17]*cMSelf[19]*mnuSelf+0.1928571428571429*cMSelf[16]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[16]*cMSelf[18]*mnuSelf+0.175662013130736*cMSelf[16]*uSelf[17]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[17]*mnuSelf+0.175662013130736*uSelf[16]*cMSelf[17]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[16]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[16]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[15]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[14]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[14]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[13]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[13]*mnuSelf+0.21957751641342*cMSelf[11]*uSelf[12]*mnuSelf+0.21957751641342*uSelf[11]*cMSelf[12]*mnuSelf+0.87831006565368*m0rSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf; 
  data->AEM_S(58,20) = 0.149071198499986*cMSelf[14]*uSelf[18]*mnuSelf+0.25*cMSelf[10]*uSelf[18]*mnuSelf+0.149071198499986*uSelf[14]*cMSelf[18]*mnuSelf+0.25*uSelf[10]*cMSelf[18]*mnuSelf+0.2195775164134199*cMSelf[13]*uSelf[16]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[16]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[14]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[14]*mnuSelf+m0rSelf[8]*mnuSelf-0.5*cESelf[8]*mnuSelf; 
  data->AEM_S(58,21) = 0.1309307341415954*cMSelf[18]*uSelf[18]*mnuSelf+0.1928571428571429*cMSelf[11]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[11]*cMSelf[18]*mnuSelf+0.1963961012123931*cMSelf[17]*uSelf[17]*mnuSelf+0.1402829294374236*cMSelf[16]*uSelf[16]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[16]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[16]*mnuSelf+0.1402829294374236*cMSelf[14]*uSelf[14]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[14]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[14]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[13]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[11]*mnuSelf+0.8783100656536796*m0rSelf[4]*mnuSelf-0.4391550328268398*cESelf[4]*mnuSelf; 
  data->AEM_S(58,23) = 0.1928571428571429*cMSelf[13]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[13]*cMSelf[18]*mnuSelf+0.175662013130736*cMSelf[13]*uSelf[17]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[15]*uSelf[16]*mnuSelf+0.1402829294374237*cMSelf[14]*uSelf[16]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[16]*mnuSelf+0.1963961012123931*uSelf[15]*cMSelf[16]*mnuSelf+0.1402829294374237*uSelf[14]*cMSelf[16]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[16]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[14]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[14]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[13]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[13]*mnuSelf+0.87831006565368*m0rSelf[6]*mnuSelf-0.43915503282684*cESelf[6]*mnuSelf; 
  data->AEM_S(58,24) = 0.2817460317460317*cMSelf[14]*uSelf[18]*mnuSelf+0.149071198499986*cMSelf[10]*uSelf[18]*mnuSelf+0.2817460317460317*uSelf[14]*cMSelf[18]*mnuSelf+0.149071198499986*uSelf[10]*cMSelf[18]*mnuSelf+0.2195775164134199*cMSelf[15]*uSelf[17]*mnuSelf+0.2195775164134199*uSelf[15]*cMSelf[17]*mnuSelf+0.3273268353539885*cMSelf[13]*uSelf[16]*mnuSelf+0.3273268353539885*uSelf[13]*cMSelf[16]*mnuSelf+0.3273268353539885*cMSelf[11]*uSelf[14]*mnuSelf+0.3273268353539885*uSelf[11]*cMSelf[14]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[13]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[13]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[11]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[11]*mnuSelf+0.5962847939999438*m0rSelf[8]*mnuSelf-0.2981423969999719*cESelf[8]*mnuSelf+0.8783100656536796*m0rSelf[1]*mnuSelf-0.4391550328268398*cESelf[1]*mnuSelf; 
  data->AEM_S(58,26) = 0.1928571428571429*cMSelf[17]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[17]*cMSelf[19]*mnuSelf+0.1928571428571429*cMSelf[16]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[16]*cMSelf[18]*mnuSelf+0.175662013130736*cMSelf[16]*uSelf[17]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[17]*mnuSelf+0.175662013130736*uSelf[16]*cMSelf[17]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[16]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[16]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[15]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[14]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[14]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[13]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[13]*mnuSelf+0.21957751641342*cMSelf[11]*uSelf[12]*mnuSelf+0.21957751641342*uSelf[11]*cMSelf[12]*mnuSelf+0.87831006565368*m0rSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf; 
  data->AEM_S(58,28) = 0.25*cMSelf[19]*uSelf[19]*mnuSelf+0.3388888888888889*cMSelf[18]*uSelf[18]*mnuSelf+0.1309307341415954*cMSelf[11]*uSelf[18]*mnuSelf+0.1309307341415954*uSelf[11]*cMSelf[18]*mnuSelf+0.3833333333333334*cMSelf[17]*uSelf[17]*mnuSelf+0.3452380952380952*cMSelf[16]*uSelf[16]*mnuSelf+0.149071198499986*cMSelf[12]*uSelf[16]*mnuSelf+0.149071198499986*uSelf[12]*cMSelf[16]*mnuSelf+0.25*cMSelf[15]*uSelf[15]*mnuSelf+0.3452380952380952*cMSelf[14]*uSelf[14]*mnuSelf+0.149071198499986*cMSelf[10]*uSelf[14]*mnuSelf+0.149071198499986*uSelf[10]*cMSelf[14]*mnuSelf+0.3833333333333334*cMSelf[13]*uSelf[13]*mnuSelf+0.25*cMSelf[12]*uSelf[12]*mnuSelf+0.3833333333333334*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.5962847939999438*m0rSelf[4]*mnuSelf-0.2981423969999719*cESelf[4]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(59,20) = 0.149071198499986*cMSelf[15]*uSelf[19]*mnuSelf+0.25*cMSelf[10]*uSelf[19]*mnuSelf+0.149071198499986*uSelf[15]*cMSelf[19]*mnuSelf+0.25*uSelf[10]*cMSelf[19]*mnuSelf+0.2195775164134199*cMSelf[13]*uSelf[17]*mnuSelf+0.2195775164134199*uSelf[13]*cMSelf[17]*mnuSelf+0.2195775164134199*cMSelf[12]*uSelf[15]*mnuSelf+0.2195775164134199*uSelf[12]*cMSelf[15]*mnuSelf+m0rSelf[9]*mnuSelf-0.5*cESelf[9]*mnuSelf; 
  data->AEM_S(59,22) = 0.1309307341415954*cMSelf[19]*uSelf[19]*mnuSelf+0.1928571428571429*cMSelf[12]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[12]*cMSelf[19]*mnuSelf+0.1402829294374236*cMSelf[17]*uSelf[17]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[17]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[16]*uSelf[16]*mnuSelf+0.1402829294374236*cMSelf[15]*uSelf[15]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[15]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[13]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[12]*mnuSelf+0.8783100656536796*m0rSelf[5]*mnuSelf-0.4391550328268398*cESelf[5]*mnuSelf; 
  data->AEM_S(59,23) = 0.1928571428571429*cMSelf[13]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[13]*cMSelf[19]*mnuSelf+0.1402829294374237*cMSelf[15]*uSelf[17]*mnuSelf+0.1963961012123931*cMSelf[14]*uSelf[17]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[17]*mnuSelf+0.1402829294374237*uSelf[15]*cMSelf[17]*mnuSelf+0.1963961012123931*uSelf[14]*cMSelf[17]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[17]*mnuSelf+0.175662013130736*cMSelf[13]*uSelf[16]*mnuSelf+0.175662013130736*uSelf[13]*cMSelf[16]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[15]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[13]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[13]*mnuSelf+0.87831006565368*m0rSelf[7]*mnuSelf-0.43915503282684*cESelf[7]*mnuSelf; 
  data->AEM_S(59,25) = 0.2817460317460317*cMSelf[15]*uSelf[19]*mnuSelf+0.149071198499986*cMSelf[10]*uSelf[19]*mnuSelf+0.2817460317460317*uSelf[15]*cMSelf[19]*mnuSelf+0.149071198499986*uSelf[10]*cMSelf[19]*mnuSelf+0.3273268353539885*cMSelf[13]*uSelf[17]*mnuSelf+0.3273268353539885*uSelf[13]*cMSelf[17]*mnuSelf+0.2195775164134199*cMSelf[14]*uSelf[16]*mnuSelf+0.2195775164134199*uSelf[14]*cMSelf[16]*mnuSelf+0.3273268353539885*cMSelf[12]*uSelf[15]*mnuSelf+0.3273268353539885*uSelf[12]*cMSelf[15]*mnuSelf+0.2195775164134199*cMSelf[11]*uSelf[13]*mnuSelf+0.2195775164134199*uSelf[11]*cMSelf[13]*mnuSelf+0.2195775164134199*cMSelf[10]*uSelf[12]*mnuSelf+0.2195775164134199*uSelf[10]*cMSelf[12]*mnuSelf+0.5962847939999438*m0rSelf[9]*mnuSelf-0.2981423969999719*cESelf[9]*mnuSelf+0.8783100656536796*m0rSelf[2]*mnuSelf-0.4391550328268398*cESelf[2]*mnuSelf; 
  data->AEM_S(59,27) = 0.1928571428571429*cMSelf[17]*uSelf[19]*mnuSelf+0.1928571428571429*uSelf[17]*cMSelf[19]*mnuSelf+0.1928571428571429*cMSelf[16]*uSelf[18]*mnuSelf+0.1928571428571429*uSelf[16]*cMSelf[18]*mnuSelf+0.175662013130736*cMSelf[16]*uSelf[17]*mnuSelf+0.1963961012123931*cMSelf[12]*uSelf[17]*mnuSelf+0.175662013130736*uSelf[16]*cMSelf[17]*mnuSelf+0.1963961012123931*uSelf[12]*cMSelf[17]*mnuSelf+0.1963961012123931*cMSelf[11]*uSelf[16]*mnuSelf+0.1963961012123931*uSelf[11]*cMSelf[16]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[15]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[15]*mnuSelf+0.1963961012123931*cMSelf[13]*uSelf[14]*mnuSelf+0.1963961012123931*uSelf[13]*cMSelf[14]*mnuSelf+0.21957751641342*cMSelf[10]*uSelf[13]*mnuSelf+0.21957751641342*uSelf[10]*cMSelf[13]*mnuSelf+0.21957751641342*cMSelf[11]*uSelf[12]*mnuSelf+0.21957751641342*uSelf[11]*cMSelf[12]*mnuSelf+0.87831006565368*m0rSelf[3]*mnuSelf-0.43915503282684*cESelf[3]*mnuSelf; 
  data->AEM_S(59,29) = 0.3388888888888889*cMSelf[19]*uSelf[19]*mnuSelf+0.1309307341415954*cMSelf[12]*uSelf[19]*mnuSelf+0.1309307341415954*uSelf[12]*cMSelf[19]*mnuSelf+0.25*cMSelf[18]*uSelf[18]*mnuSelf+0.3452380952380952*cMSelf[17]*uSelf[17]*mnuSelf+0.149071198499986*cMSelf[11]*uSelf[17]*mnuSelf+0.149071198499986*uSelf[11]*cMSelf[17]*mnuSelf+0.3833333333333334*cMSelf[16]*uSelf[16]*mnuSelf+0.3452380952380952*cMSelf[15]*uSelf[15]*mnuSelf+0.149071198499986*cMSelf[10]*uSelf[15]*mnuSelf+0.149071198499986*uSelf[10]*cMSelf[15]*mnuSelf+0.25*cMSelf[14]*uSelf[14]*mnuSelf+0.3833333333333334*cMSelf[13]*uSelf[13]*mnuSelf+0.3833333333333334*cMSelf[12]*uSelf[12]*mnuSelf+0.25*cMSelf[11]*uSelf[11]*mnuSelf+0.25*cMSelf[10]*uSelf[10]*mnuSelf+0.5962847939999438*m0rSelf[5]*mnuSelf-0.2981423969999719*cESelf[5]*mnuSelf+m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(50,50) = (-0.25*cMOther[19]*uOther[19]*mnuOther)-0.25*cMOther[18]*uOther[18]*mnuOther-0.25*cMOther[17]*uOther[17]*mnuOther-0.25*cMOther[16]*uOther[16]*mnuOther-0.25*cMOther[15]*uOther[15]*mnuOther-0.25*cMOther[14]*uOther[14]*mnuOther-0.25*cMOther[13]*uOther[13]*mnuOther-0.25*cMOther[12]*uOther[12]*mnuOther-0.25*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(50,51) = (-0.2195775164134199*cMOther[14]*uOther[18]*mnuOther)-0.2195775164134199*uOther[14]*cMOther[18]*mnuOther-0.2500000000000001*cMOther[15]*uOther[17]*mnuOther-0.2500000000000001*uOther[15]*cMOther[17]*mnuOther-0.223606797749979*cMOther[13]*uOther[16]*mnuOther-0.223606797749979*uOther[13]*cMOther[16]*mnuOther-0.223606797749979*cMOther[11]*uOther[14]*mnuOther-0.223606797749979*uOther[11]*cMOther[14]*mnuOther-0.25*cMOther[12]*uOther[13]*mnuOther-0.25*uOther[12]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[11]*mnuOther-0.25*uOther[10]*cMOther[11]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(50,52) = (-0.2195775164134199*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*uOther[15]*cMOther[19]*mnuOther-0.223606797749979*cMOther[13]*uOther[17]*mnuOther-0.223606797749979*uOther[13]*cMOther[17]*mnuOther-0.2500000000000001*cMOther[14]*uOther[16]*mnuOther-0.2500000000000001*uOther[14]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[15]*mnuOther-0.223606797749979*uOther[12]*cMOther[15]*mnuOther-0.25*cMOther[11]*uOther[13]*mnuOther-0.25*uOther[11]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[12]*mnuOther-0.25*uOther[10]*cMOther[12]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(50,53) = (-0.2195775164134199*cMOther[17]*uOther[19]*mnuOther)-0.2195775164134199*uOther[17]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[16]*uOther[18]*mnuOther-0.2195775164134199*uOther[16]*cMOther[18]*mnuOther-0.2*cMOther[16]*uOther[17]*mnuOther-0.223606797749979*cMOther[12]*uOther[17]*mnuOther-0.2*uOther[16]*cMOther[17]*mnuOther-0.223606797749979*uOther[12]*cMOther[17]*mnuOther-0.223606797749979*cMOther[11]*uOther[16]*mnuOther-0.223606797749979*uOther[11]*cMOther[16]*mnuOther-0.223606797749979*cMOther[13]*uOther[15]*mnuOther-0.223606797749979*uOther[13]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[14]*mnuOther-0.223606797749979*uOther[13]*cMOther[14]*mnuOther-0.25*cMOther[10]*uOther[13]*mnuOther-0.25*uOther[10]*cMOther[13]*mnuOther-0.25*cMOther[11]*uOther[12]*mnuOther-0.25*uOther[11]*cMOther[12]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(50,54) = (-0.149071198499986*cMOther[18]*uOther[18]*mnuOther)-0.2195775164134199*cMOther[11]*uOther[18]*mnuOther-0.2195775164134199*uOther[11]*cMOther[18]*mnuOther-0.223606797749979*cMOther[17]*uOther[17]*mnuOther-0.159719141249985*cMOther[16]*uOther[16]*mnuOther-0.2500000000000001*cMOther[12]*uOther[16]*mnuOther-0.2500000000000001*uOther[12]*cMOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[14]*mnuOther-0.25*cMOther[10]*uOther[14]*mnuOther-0.25*uOther[10]*cMOther[14]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[11]*mnuOther-1.0*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(50,55) = (-0.149071198499986*cMOther[19]*uOther[19]*mnuOther)-0.2195775164134199*cMOther[12]*uOther[19]*mnuOther-0.2195775164134199*uOther[12]*cMOther[19]*mnuOther-0.159719141249985*cMOther[17]*uOther[17]*mnuOther-0.2500000000000001*cMOther[11]*uOther[17]*mnuOther-0.2500000000000001*uOther[11]*cMOther[17]*mnuOther-0.223606797749979*cMOther[16]*uOther[16]*mnuOther-0.159719141249985*cMOther[15]*uOther[15]*mnuOther-0.25*cMOther[10]*uOther[15]*mnuOther-0.25*uOther[10]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[12]*uOther[12]*mnuOther-1.0*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(50,56) = (-0.2195775164134199*cMOther[13]*uOther[18]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[18]*mnuOther-0.2*cMOther[13]*uOther[17]*mnuOther-0.2*uOther[13]*cMOther[17]*mnuOther-0.223606797749979*cMOther[15]*uOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[16]*mnuOther-0.25*cMOther[10]*uOther[16]*mnuOther-0.223606797749979*uOther[15]*cMOther[16]*mnuOther-0.159719141249985*uOther[14]*cMOther[16]*mnuOther-0.25*uOther[10]*cMOther[16]*mnuOther-0.2500000000000001*cMOther[12]*uOther[14]*mnuOther-0.2500000000000001*uOther[12]*cMOther[14]*mnuOther-0.223606797749979*cMOther[11]*uOther[13]*mnuOther-0.223606797749979*uOther[11]*cMOther[13]*mnuOther-1.0*m0rOther[6]*mnuOther+0.5*cEOther[6]*mnuOther; 
  data->AEM_S(50,57) = (-0.2195775164134199*cMOther[13]*uOther[19]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[19]*mnuOther-0.159719141249985*cMOther[15]*uOther[17]*mnuOther-0.223606797749979*cMOther[14]*uOther[17]*mnuOther-0.25*cMOther[10]*uOther[17]*mnuOther-0.159719141249985*uOther[15]*cMOther[17]*mnuOther-0.223606797749979*uOther[14]*cMOther[17]*mnuOther-0.25*uOther[10]*cMOther[17]*mnuOther-0.2*cMOther[13]*uOther[16]*mnuOther-0.2*uOther[13]*cMOther[16]*mnuOther-0.2500000000000001*cMOther[11]*uOther[15]*mnuOther-0.2500000000000001*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[12]*uOther[13]*mnuOther-0.223606797749979*uOther[12]*cMOther[13]*mnuOther-1.0*m0rOther[7]*mnuOther+0.5*cEOther[7]*mnuOther; 
  data->AEM_S(50,58) = (-0.149071198499986*cMOther[14]*uOther[18]*mnuOther)-0.25*cMOther[10]*uOther[18]*mnuOther-0.149071198499986*uOther[14]*cMOther[18]*mnuOther-0.25*uOther[10]*cMOther[18]*mnuOther-0.2195775164134199*cMOther[13]*uOther[16]*mnuOther-0.2195775164134199*uOther[13]*cMOther[16]*mnuOther-0.2195775164134199*cMOther[11]*uOther[14]*mnuOther-0.2195775164134199*uOther[11]*cMOther[14]*mnuOther-1.0*m0rOther[8]*mnuOther+0.5*cEOther[8]*mnuOther; 
  data->AEM_S(50,59) = (-0.149071198499986*cMOther[15]*uOther[19]*mnuOther)-0.25*cMOther[10]*uOther[19]*mnuOther-0.149071198499986*uOther[15]*cMOther[19]*mnuOther-0.25*uOther[10]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[13]*uOther[17]*mnuOther-0.2195775164134199*uOther[13]*cMOther[17]*mnuOther-0.2195775164134199*cMOther[12]*uOther[15]*mnuOther-0.2195775164134199*uOther[12]*cMOther[15]*mnuOther-1.0*m0rOther[9]*mnuOther+0.5*cEOther[9]*mnuOther; 
  data->AEM_S(51,50) = (-0.2195775164134199*cMOther[14]*uOther[18]*mnuOther)-0.2195775164134199*uOther[14]*cMOther[18]*mnuOther-0.2500000000000001*cMOther[15]*uOther[17]*mnuOther-0.2500000000000001*uOther[15]*cMOther[17]*mnuOther-0.223606797749979*cMOther[13]*uOther[16]*mnuOther-0.223606797749979*uOther[13]*cMOther[16]*mnuOther-0.223606797749979*cMOther[11]*uOther[14]*mnuOther-0.223606797749979*uOther[11]*cMOther[14]*mnuOther-0.25*cMOther[12]*uOther[13]*mnuOther-0.25*uOther[12]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[11]*mnuOther-0.25*uOther[10]*cMOther[11]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(51,51) = (-0.25*cMOther[19]*uOther[19]*mnuOther)-0.3833333333333334*cMOther[18]*uOther[18]*mnuOther-0.1963961012123931*cMOther[11]*uOther[18]*mnuOther-0.1963961012123931*uOther[11]*cMOther[18]*mnuOther-0.45*cMOther[17]*uOther[17]*mnuOther-0.3928571428571428*cMOther[16]*uOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[16]*mnuOther-0.223606797749979*uOther[12]*cMOther[16]*mnuOther-0.25*cMOther[15]*uOther[15]*mnuOther-0.3928571428571428*cMOther[14]*uOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[14]*mnuOther-0.223606797749979*uOther[10]*cMOther[14]*mnuOther-0.45*cMOther[13]*uOther[13]*mnuOther-0.25*cMOther[12]*uOther[12]*mnuOther-0.45*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.8944271909999159*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(51,52) = (-0.2195775164134199*cMOther[17]*uOther[19]*mnuOther)-0.2195775164134199*uOther[17]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[16]*uOther[18]*mnuOther-0.2195775164134199*uOther[16]*cMOther[18]*mnuOther-0.2*cMOther[16]*uOther[17]*mnuOther-0.223606797749979*cMOther[12]*uOther[17]*mnuOther-0.2*uOther[16]*cMOther[17]*mnuOther-0.223606797749979*uOther[12]*cMOther[17]*mnuOther-0.223606797749979*cMOther[11]*uOther[16]*mnuOther-0.223606797749979*uOther[11]*cMOther[16]*mnuOther-0.223606797749979*cMOther[13]*uOther[15]*mnuOther-0.223606797749979*uOther[13]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[14]*mnuOther-0.223606797749979*uOther[13]*cMOther[14]*mnuOther-0.25*cMOther[10]*uOther[13]*mnuOther-0.25*uOther[10]*cMOther[13]*mnuOther-0.25*cMOther[11]*uOther[12]*mnuOther-0.25*uOther[11]*cMOther[12]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(51,53) = (-0.2195775164134199*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*uOther[15]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[13]*uOther[18]*mnuOther-0.1963961012123931*uOther[13]*cMOther[18]*mnuOther-0.4024922359499621*cMOther[13]*uOther[17]*mnuOther-0.4024922359499621*uOther[13]*cMOther[17]*mnuOther-0.2*cMOther[15]*uOther[16]*mnuOther-0.3928571428571429*cMOther[14]*uOther[16]*mnuOther-0.223606797749979*cMOther[10]*uOther[16]*mnuOther-0.2*uOther[15]*cMOther[16]*mnuOther-0.3928571428571429*uOther[14]*cMOther[16]*mnuOther-0.223606797749979*uOther[10]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[15]*mnuOther-0.223606797749979*uOther[12]*cMOther[15]*mnuOther-0.223606797749979*cMOther[12]*uOther[14]*mnuOther-0.223606797749979*uOther[12]*cMOther[14]*mnuOther-0.45*cMOther[11]*uOther[13]*mnuOther-0.45*uOther[11]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[12]*mnuOther-0.25*uOther[10]*cMOther[12]*mnuOther-0.8944271909999161*m0rOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(51,54) = (-0.3273268353539885*cMOther[14]*uOther[18]*mnuOther)-0.2195775164134199*cMOther[10]*uOther[18]*mnuOther-0.3273268353539885*uOther[14]*cMOther[18]*mnuOther-0.2195775164134199*uOther[10]*cMOther[18]*mnuOther-0.223606797749979*cMOther[15]*uOther[17]*mnuOther-0.223606797749979*uOther[15]*cMOther[17]*mnuOther-0.3928571428571429*cMOther[13]*uOther[16]*mnuOther-0.3928571428571429*uOther[13]*cMOther[16]*mnuOther-0.3928571428571428*cMOther[11]*uOther[14]*mnuOther-0.3928571428571428*uOther[11]*cMOther[14]*mnuOther-0.223606797749979*cMOther[12]*uOther[13]*mnuOther-0.223606797749979*uOther[12]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[11]*mnuOther-0.223606797749979*uOther[10]*cMOther[11]*mnuOther-0.8783100656536796*m0rOther[8]*mnuOther+0.4391550328268398*cEOther[8]*mnuOther-0.8944271909999159*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(51,55) = (-0.2195775164134199*cMOther[13]*uOther[19]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[19]*mnuOther-0.159719141249985*cMOther[15]*uOther[17]*mnuOther-0.223606797749979*cMOther[14]*uOther[17]*mnuOther-0.2500000000000001*cMOther[10]*uOther[17]*mnuOther-0.159719141249985*uOther[15]*cMOther[17]*mnuOther-0.223606797749979*uOther[14]*cMOther[17]*mnuOther-0.2500000000000001*uOther[10]*cMOther[17]*mnuOther-0.2*cMOther[13]*uOther[16]*mnuOther-0.2*uOther[13]*cMOther[16]*mnuOther-0.25*cMOther[11]*uOther[15]*mnuOther-0.25*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[12]*uOther[13]*mnuOther-0.223606797749979*uOther[12]*cMOther[13]*mnuOther-1.0*m0rOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(51,56) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999832*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999832*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999161*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(51,57) = (-0.149071198499986*cMOther[19]*uOther[19]*mnuOther)-0.2195775164134199*cMOther[12]*uOther[19]*mnuOther-0.2195775164134199*uOther[12]*cMOther[19]*mnuOther-0.159719141249985*cMOther[17]*uOther[17]*mnuOther-0.25*cMOther[11]*uOther[17]*mnuOther-0.25*uOther[11]*cMOther[17]*mnuOther-0.223606797749979*cMOther[16]*uOther[16]*mnuOther-0.159719141249985*cMOther[15]*uOther[15]*mnuOther-0.2500000000000001*cMOther[10]*uOther[15]*mnuOther-0.2500000000000001*uOther[10]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[12]*uOther[12]*mnuOther-1.0*m0rOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(51,58) = (-0.1309307341415954*cMOther[18]*uOther[18]*mnuOther)-0.1928571428571429*cMOther[11]*uOther[18]*mnuOther-0.1928571428571429*uOther[11]*cMOther[18]*mnuOther-0.1963961012123931*cMOther[17]*uOther[17]*mnuOther-0.1402829294374236*cMOther[16]*uOther[16]*mnuOther-0.2195775164134199*cMOther[12]*uOther[16]*mnuOther-0.2195775164134199*uOther[12]*cMOther[16]*mnuOther-0.1402829294374236*cMOther[14]*uOther[14]*mnuOther-0.2195775164134199*cMOther[10]*uOther[14]*mnuOther-0.2195775164134199*uOther[10]*cMOther[14]*mnuOther-0.1963961012123931*cMOther[13]*uOther[13]*mnuOther-0.1963961012123931*cMOther[11]*uOther[11]*mnuOther-0.8783100656536796*m0rOther[4]*mnuOther+0.4391550328268398*cEOther[4]*mnuOther; 
  data->AEM_S(52,50) = (-0.2195775164134199*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*uOther[15]*cMOther[19]*mnuOther-0.223606797749979*cMOther[13]*uOther[17]*mnuOther-0.223606797749979*uOther[13]*cMOther[17]*mnuOther-0.2500000000000001*cMOther[14]*uOther[16]*mnuOther-0.2500000000000001*uOther[14]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[15]*mnuOther-0.223606797749979*uOther[12]*cMOther[15]*mnuOther-0.25*cMOther[11]*uOther[13]*mnuOther-0.25*uOther[11]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[12]*mnuOther-0.25*uOther[10]*cMOther[12]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(52,51) = (-0.2195775164134199*cMOther[17]*uOther[19]*mnuOther)-0.2195775164134199*uOther[17]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[16]*uOther[18]*mnuOther-0.2195775164134199*uOther[16]*cMOther[18]*mnuOther-0.2*cMOther[16]*uOther[17]*mnuOther-0.223606797749979*cMOther[12]*uOther[17]*mnuOther-0.2*uOther[16]*cMOther[17]*mnuOther-0.223606797749979*uOther[12]*cMOther[17]*mnuOther-0.223606797749979*cMOther[11]*uOther[16]*mnuOther-0.223606797749979*uOther[11]*cMOther[16]*mnuOther-0.223606797749979*cMOther[13]*uOther[15]*mnuOther-0.223606797749979*uOther[13]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[14]*mnuOther-0.223606797749979*uOther[13]*cMOther[14]*mnuOther-0.25*cMOther[10]*uOther[13]*mnuOther-0.25*uOther[10]*cMOther[13]*mnuOther-0.25*cMOther[11]*uOther[12]*mnuOther-0.25*uOther[11]*cMOther[12]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(52,52) = (-0.3833333333333334*cMOther[19]*uOther[19]*mnuOther)-0.1963961012123931*cMOther[12]*uOther[19]*mnuOther-0.1963961012123931*uOther[12]*cMOther[19]*mnuOther-0.25*cMOther[18]*uOther[18]*mnuOther-0.3928571428571428*cMOther[17]*uOther[17]*mnuOther-0.223606797749979*cMOther[11]*uOther[17]*mnuOther-0.223606797749979*uOther[11]*cMOther[17]*mnuOther-0.45*cMOther[16]*uOther[16]*mnuOther-0.3928571428571428*cMOther[15]*uOther[15]*mnuOther-0.223606797749979*cMOther[10]*uOther[15]*mnuOther-0.223606797749979*uOther[10]*cMOther[15]*mnuOther-0.25*cMOther[14]*uOther[14]*mnuOther-0.45*cMOther[13]*uOther[13]*mnuOther-0.45*cMOther[12]*uOther[12]*mnuOther-0.25*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.8944271909999159*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(52,53) = (-0.1963961012123931*cMOther[13]*uOther[19]*mnuOther)-0.1963961012123931*uOther[13]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[14]*uOther[18]*mnuOther-0.2195775164134199*uOther[14]*cMOther[18]*mnuOther-0.3928571428571429*cMOther[15]*uOther[17]*mnuOther-0.2*cMOther[14]*uOther[17]*mnuOther-0.223606797749979*cMOther[10]*uOther[17]*mnuOther-0.3928571428571429*uOther[15]*cMOther[17]*mnuOther-0.2*uOther[14]*cMOther[17]*mnuOther-0.223606797749979*uOther[10]*cMOther[17]*mnuOther-0.4024922359499621*cMOther[13]*uOther[16]*mnuOther-0.4024922359499621*uOther[13]*cMOther[16]*mnuOther-0.223606797749979*cMOther[11]*uOther[15]*mnuOther-0.223606797749979*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[11]*uOther[14]*mnuOther-0.223606797749979*uOther[11]*cMOther[14]*mnuOther-0.45*cMOther[12]*uOther[13]*mnuOther-0.45*uOther[12]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[11]*mnuOther-0.25*uOther[10]*cMOther[11]*mnuOther-0.8944271909999161*m0rOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(52,54) = (-0.2195775164134199*cMOther[13]*uOther[18]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[18]*mnuOther-0.2*cMOther[13]*uOther[17]*mnuOther-0.2*uOther[13]*cMOther[17]*mnuOther-0.223606797749979*cMOther[15]*uOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[16]*mnuOther-0.2500000000000001*cMOther[10]*uOther[16]*mnuOther-0.223606797749979*uOther[15]*cMOther[16]*mnuOther-0.159719141249985*uOther[14]*cMOther[16]*mnuOther-0.2500000000000001*uOther[10]*cMOther[16]*mnuOther-0.25*cMOther[12]*uOther[14]*mnuOther-0.25*uOther[12]*cMOther[14]*mnuOther-0.223606797749979*cMOther[11]*uOther[13]*mnuOther-0.223606797749979*uOther[11]*cMOther[13]*mnuOther-1.0*m0rOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(52,55) = (-0.3273268353539885*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*cMOther[10]*uOther[19]*mnuOther-0.3273268353539885*uOther[15]*cMOther[19]*mnuOther-0.2195775164134199*uOther[10]*cMOther[19]*mnuOther-0.3928571428571429*cMOther[13]*uOther[17]*mnuOther-0.3928571428571429*uOther[13]*cMOther[17]*mnuOther-0.223606797749979*cMOther[14]*uOther[16]*mnuOther-0.223606797749979*uOther[14]*cMOther[16]*mnuOther-0.3928571428571428*cMOther[12]*uOther[15]*mnuOther-0.3928571428571428*uOther[12]*cMOther[15]*mnuOther-0.223606797749979*cMOther[11]*uOther[13]*mnuOther-0.223606797749979*uOther[11]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[12]*mnuOther-0.223606797749979*uOther[10]*cMOther[12]*mnuOther-0.8783100656536796*m0rOther[9]*mnuOther+0.4391550328268398*cEOther[9]*mnuOther-0.8944271909999159*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(52,56) = (-0.149071198499986*cMOther[18]*uOther[18]*mnuOther)-0.2195775164134199*cMOther[11]*uOther[18]*mnuOther-0.2195775164134199*uOther[11]*cMOther[18]*mnuOther-0.223606797749979*cMOther[17]*uOther[17]*mnuOther-0.159719141249985*cMOther[16]*uOther[16]*mnuOther-0.25*cMOther[12]*uOther[16]*mnuOther-0.25*uOther[12]*cMOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[14]*mnuOther-0.2500000000000001*cMOther[10]*uOther[14]*mnuOther-0.2500000000000001*uOther[10]*cMOther[14]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[11]*mnuOther-1.0*m0rOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(52,57) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999832*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999832*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999161*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(52,59) = (-0.1309307341415954*cMOther[19]*uOther[19]*mnuOther)-0.1928571428571429*cMOther[12]*uOther[19]*mnuOther-0.1928571428571429*uOther[12]*cMOther[19]*mnuOther-0.1402829294374236*cMOther[17]*uOther[17]*mnuOther-0.2195775164134199*cMOther[11]*uOther[17]*mnuOther-0.2195775164134199*uOther[11]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[16]*uOther[16]*mnuOther-0.1402829294374236*cMOther[15]*uOther[15]*mnuOther-0.2195775164134199*cMOther[10]*uOther[15]*mnuOther-0.2195775164134199*uOther[10]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[13]*uOther[13]*mnuOther-0.1963961012123931*cMOther[12]*uOther[12]*mnuOther-0.8783100656536796*m0rOther[5]*mnuOther+0.4391550328268398*cEOther[5]*mnuOther; 
  data->AEM_S(53,50) = (-0.2195775164134199*cMOther[17]*uOther[19]*mnuOther)-0.2195775164134199*uOther[17]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[16]*uOther[18]*mnuOther-0.2195775164134199*uOther[16]*cMOther[18]*mnuOther-0.2*cMOther[16]*uOther[17]*mnuOther-0.223606797749979*cMOther[12]*uOther[17]*mnuOther-0.2*uOther[16]*cMOther[17]*mnuOther-0.223606797749979*uOther[12]*cMOther[17]*mnuOther-0.223606797749979*cMOther[11]*uOther[16]*mnuOther-0.223606797749979*uOther[11]*cMOther[16]*mnuOther-0.223606797749979*cMOther[13]*uOther[15]*mnuOther-0.223606797749979*uOther[13]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[14]*mnuOther-0.223606797749979*uOther[13]*cMOther[14]*mnuOther-0.25*cMOther[10]*uOther[13]*mnuOther-0.25*uOther[10]*cMOther[13]*mnuOther-0.25*cMOther[11]*uOther[12]*mnuOther-0.25*uOther[11]*cMOther[12]*mnuOther-1.0*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(53,51) = (-0.2195775164134199*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*uOther[15]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[13]*uOther[18]*mnuOther-0.1963961012123931*uOther[13]*cMOther[18]*mnuOther-0.4024922359499621*cMOther[13]*uOther[17]*mnuOther-0.4024922359499621*uOther[13]*cMOther[17]*mnuOther-0.2*cMOther[15]*uOther[16]*mnuOther-0.3928571428571429*cMOther[14]*uOther[16]*mnuOther-0.223606797749979*cMOther[10]*uOther[16]*mnuOther-0.2*uOther[15]*cMOther[16]*mnuOther-0.3928571428571429*uOther[14]*cMOther[16]*mnuOther-0.223606797749979*uOther[10]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[15]*mnuOther-0.223606797749979*uOther[12]*cMOther[15]*mnuOther-0.223606797749979*cMOther[12]*uOther[14]*mnuOther-0.223606797749979*uOther[12]*cMOther[14]*mnuOther-0.45*cMOther[11]*uOther[13]*mnuOther-0.45*uOther[11]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[12]*mnuOther-0.25*uOther[10]*cMOther[12]*mnuOther-0.8944271909999161*m0rOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(53,52) = (-0.1963961012123931*cMOther[13]*uOther[19]*mnuOther)-0.1963961012123931*uOther[13]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[14]*uOther[18]*mnuOther-0.2195775164134199*uOther[14]*cMOther[18]*mnuOther-0.3928571428571429*cMOther[15]*uOther[17]*mnuOther-0.2*cMOther[14]*uOther[17]*mnuOther-0.223606797749979*cMOther[10]*uOther[17]*mnuOther-0.3928571428571429*uOther[15]*cMOther[17]*mnuOther-0.2*uOther[14]*cMOther[17]*mnuOther-0.223606797749979*uOther[10]*cMOther[17]*mnuOther-0.4024922359499621*cMOther[13]*uOther[16]*mnuOther-0.4024922359499621*uOther[13]*cMOther[16]*mnuOther-0.223606797749979*cMOther[11]*uOther[15]*mnuOther-0.223606797749979*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[11]*uOther[14]*mnuOther-0.223606797749979*uOther[11]*cMOther[14]*mnuOther-0.45*cMOther[12]*uOther[13]*mnuOther-0.45*uOther[12]*cMOther[13]*mnuOther-0.25*cMOther[10]*uOther[11]*mnuOther-0.25*uOther[10]*cMOther[11]*mnuOther-0.8944271909999161*m0rOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(53,53) = (-0.3833333333333334*cMOther[19]*uOther[19]*mnuOther)-0.1963961012123931*cMOther[12]*uOther[19]*mnuOther-0.1963961012123931*uOther[12]*cMOther[19]*mnuOther-0.3833333333333334*cMOther[18]*uOther[18]*mnuOther-0.1963961012123931*cMOther[11]*uOther[18]*mnuOther-0.1963961012123931*uOther[11]*cMOther[18]*mnuOther-0.5928571428571429*cMOther[17]*uOther[17]*mnuOther-0.223606797749979*cMOther[11]*uOther[17]*mnuOther-0.223606797749979*uOther[11]*cMOther[17]*mnuOther-0.5928571428571429*cMOther[16]*uOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[16]*mnuOther-0.223606797749979*uOther[12]*cMOther[16]*mnuOther-0.3928571428571428*cMOther[15]*uOther[15]*mnuOther-0.223606797749979*cMOther[10]*uOther[15]*mnuOther-0.223606797749979*uOther[10]*cMOther[15]*mnuOther-0.3928571428571428*cMOther[14]*uOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[14]*mnuOther-0.223606797749979*uOther[10]*cMOther[14]*mnuOther-0.65*cMOther[13]*uOther[13]*mnuOther-0.45*cMOther[12]*uOther[12]*mnuOther-0.45*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.8944271909999159*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.8944271909999159*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(53,54) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999831*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999831*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(53,55) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999831*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999831*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(53,56) = (-0.175662013130736*cMOther[13]*uOther[19]*mnuOther)-0.175662013130736*uOther[13]*cMOther[19]*mnuOther-0.3273268353539885*cMOther[14]*uOther[18]*mnuOther-0.21957751641342*cMOther[10]*uOther[18]*mnuOther-0.3273268353539885*uOther[14]*cMOther[18]*mnuOther-0.21957751641342*uOther[10]*cMOther[18]*mnuOther-0.3513821107499669*cMOther[15]*uOther[17]*mnuOther-0.1788854381999831*cMOther[14]*uOther[17]*mnuOther-0.2*cMOther[10]*uOther[17]*mnuOther-0.3513821107499669*uOther[15]*cMOther[17]*mnuOther-0.1788854381999831*uOther[14]*cMOther[17]*mnuOther-0.2*uOther[10]*cMOther[17]*mnuOther-0.5528571428571428*cMOther[13]*uOther[16]*mnuOther-0.5528571428571428*uOther[13]*cMOther[16]*mnuOther-0.2*cMOther[11]*uOther[15]*mnuOther-0.2*uOther[11]*cMOther[15]*mnuOther-0.3928571428571429*cMOther[11]*uOther[14]*mnuOther-0.3928571428571429*uOther[11]*cMOther[14]*mnuOther-0.4024922359499621*cMOther[12]*uOther[13]*mnuOther-0.4024922359499621*uOther[12]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[11]*mnuOther-0.223606797749979*uOther[10]*cMOther[11]*mnuOther-0.87831006565368*m0rOther[8]*mnuOther+0.43915503282684*cEOther[8]*mnuOther-0.8*m0rOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.8944271909999161*m0rOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(53,57) = (-0.3273268353539885*cMOther[15]*uOther[19]*mnuOther)-0.21957751641342*cMOther[10]*uOther[19]*mnuOther-0.3273268353539885*uOther[15]*cMOther[19]*mnuOther-0.21957751641342*uOther[10]*cMOther[19]*mnuOther-0.175662013130736*cMOther[13]*uOther[18]*mnuOther-0.175662013130736*uOther[13]*cMOther[18]*mnuOther-0.5528571428571428*cMOther[13]*uOther[17]*mnuOther-0.5528571428571428*uOther[13]*cMOther[17]*mnuOther-0.1788854381999831*cMOther[15]*uOther[16]*mnuOther-0.3513821107499669*cMOther[14]*uOther[16]*mnuOther-0.2*cMOther[10]*uOther[16]*mnuOther-0.1788854381999831*uOther[15]*cMOther[16]*mnuOther-0.3513821107499669*uOther[14]*cMOther[16]*mnuOther-0.2*uOther[10]*cMOther[16]*mnuOther-0.3928571428571429*cMOther[12]*uOther[15]*mnuOther-0.3928571428571429*uOther[12]*cMOther[15]*mnuOther-0.2*cMOther[12]*uOther[14]*mnuOther-0.2*uOther[12]*cMOther[14]*mnuOther-0.4024922359499621*cMOther[11]*uOther[13]*mnuOther-0.4024922359499621*uOther[11]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[12]*mnuOther-0.223606797749979*uOther[10]*cMOther[12]*mnuOther-0.87831006565368*m0rOther[9]*mnuOther+0.43915503282684*cEOther[9]*mnuOther-0.8*m0rOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.8944271909999161*m0rOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(53,58) = (-0.1928571428571429*cMOther[13]*uOther[18]*mnuOther)-0.1928571428571429*uOther[13]*cMOther[18]*mnuOther-0.175662013130736*cMOther[13]*uOther[17]*mnuOther-0.175662013130736*uOther[13]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[15]*uOther[16]*mnuOther-0.1402829294374237*cMOther[14]*uOther[16]*mnuOther-0.21957751641342*cMOther[10]*uOther[16]*mnuOther-0.1963961012123931*uOther[15]*cMOther[16]*mnuOther-0.1402829294374237*uOther[14]*cMOther[16]*mnuOther-0.21957751641342*uOther[10]*cMOther[16]*mnuOther-0.2195775164134199*cMOther[12]*uOther[14]*mnuOther-0.2195775164134199*uOther[12]*cMOther[14]*mnuOther-0.1963961012123931*cMOther[11]*uOther[13]*mnuOther-0.1963961012123931*uOther[11]*cMOther[13]*mnuOther-0.87831006565368*m0rOther[6]*mnuOther+0.43915503282684*cEOther[6]*mnuOther; 
  data->AEM_S(53,59) = (-0.1928571428571429*cMOther[13]*uOther[19]*mnuOther)-0.1928571428571429*uOther[13]*cMOther[19]*mnuOther-0.1402829294374237*cMOther[15]*uOther[17]*mnuOther-0.1963961012123931*cMOther[14]*uOther[17]*mnuOther-0.21957751641342*cMOther[10]*uOther[17]*mnuOther-0.1402829294374237*uOther[15]*cMOther[17]*mnuOther-0.1963961012123931*uOther[14]*cMOther[17]*mnuOther-0.21957751641342*uOther[10]*cMOther[17]*mnuOther-0.175662013130736*cMOther[13]*uOther[16]*mnuOther-0.175662013130736*uOther[13]*cMOther[16]*mnuOther-0.2195775164134199*cMOther[11]*uOther[15]*mnuOther-0.2195775164134199*uOther[11]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[12]*uOther[13]*mnuOther-0.1963961012123931*uOther[12]*cMOther[13]*mnuOther-0.87831006565368*m0rOther[7]*mnuOther+0.43915503282684*cEOther[7]*mnuOther; 
  data->AEM_S(54,50) = (-0.149071198499986*cMOther[18]*uOther[18]*mnuOther)-0.2195775164134199*cMOther[11]*uOther[18]*mnuOther-0.2195775164134199*uOther[11]*cMOther[18]*mnuOther-0.223606797749979*cMOther[17]*uOther[17]*mnuOther-0.159719141249985*cMOther[16]*uOther[16]*mnuOther-0.2500000000000001*cMOther[12]*uOther[16]*mnuOther-0.2500000000000001*uOther[12]*cMOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[14]*mnuOther-0.25*cMOther[10]*uOther[14]*mnuOther-0.25*uOther[10]*cMOther[14]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[11]*mnuOther-1.0*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(54,51) = (-0.3273268353539885*cMOther[14]*uOther[18]*mnuOther)-0.2195775164134199*cMOther[10]*uOther[18]*mnuOther-0.3273268353539885*uOther[14]*cMOther[18]*mnuOther-0.2195775164134199*uOther[10]*cMOther[18]*mnuOther-0.223606797749979*cMOther[15]*uOther[17]*mnuOther-0.223606797749979*uOther[15]*cMOther[17]*mnuOther-0.3928571428571429*cMOther[13]*uOther[16]*mnuOther-0.3928571428571429*uOther[13]*cMOther[16]*mnuOther-0.3928571428571428*cMOther[11]*uOther[14]*mnuOther-0.3928571428571428*uOther[11]*cMOther[14]*mnuOther-0.223606797749979*cMOther[12]*uOther[13]*mnuOther-0.223606797749979*uOther[12]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[11]*mnuOther-0.223606797749979*uOther[10]*cMOther[11]*mnuOther-0.8783100656536796*m0rOther[8]*mnuOther+0.4391550328268398*cEOther[8]*mnuOther-0.8944271909999159*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(54,52) = (-0.2195775164134199*cMOther[13]*uOther[18]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[18]*mnuOther-0.2*cMOther[13]*uOther[17]*mnuOther-0.2*uOther[13]*cMOther[17]*mnuOther-0.223606797749979*cMOther[15]*uOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[16]*mnuOther-0.2500000000000001*cMOther[10]*uOther[16]*mnuOther-0.223606797749979*uOther[15]*cMOther[16]*mnuOther-0.159719141249985*uOther[14]*cMOther[16]*mnuOther-0.2500000000000001*uOther[10]*cMOther[16]*mnuOther-0.25*cMOther[12]*uOther[14]*mnuOther-0.25*uOther[12]*cMOther[14]*mnuOther-0.223606797749979*cMOther[11]*uOther[13]*mnuOther-0.223606797749979*uOther[11]*cMOther[13]*mnuOther-1.0*m0rOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(54,53) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999831*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999831*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(54,54) = (-0.25*cMOther[19]*uOther[19]*mnuOther)-0.3452380952380952*cMOther[18]*uOther[18]*mnuOther-0.1402829294374236*cMOther[11]*uOther[18]*mnuOther-0.1402829294374236*uOther[11]*cMOther[18]*mnuOther-0.3928571428571428*cMOther[17]*uOther[17]*mnuOther-0.3520408163265306*cMOther[16]*uOther[16]*mnuOther-0.159719141249985*cMOther[12]*uOther[16]*mnuOther-0.159719141249985*uOther[12]*cMOther[16]*mnuOther-0.25*cMOther[15]*uOther[15]*mnuOther-0.3520408163265306*cMOther[14]*uOther[14]*mnuOther-0.159719141249985*cMOther[10]*uOther[14]*mnuOther-0.159719141249985*uOther[10]*cMOther[14]*mnuOther-0.3928571428571428*cMOther[13]*uOther[13]*mnuOther-0.25*cMOther[12]*uOther[12]*mnuOther-0.3928571428571428*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.6388765649999399*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(54,56) = (-0.2195775164134199*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*uOther[15]*cMOther[19]*mnuOther-0.1402829294374237*cMOther[13]*uOther[18]*mnuOther-0.1402829294374237*uOther[13]*cMOther[18]*mnuOther-0.3513821107499669*cMOther[13]*uOther[17]*mnuOther-0.3513821107499669*uOther[13]*cMOther[17]*mnuOther-0.1428571428571428*cMOther[15]*uOther[16]*mnuOther-0.3520408163265306*cMOther[14]*uOther[16]*mnuOther-0.159719141249985*cMOther[10]*uOther[16]*mnuOther-0.1428571428571428*uOther[15]*cMOther[16]*mnuOther-0.3520408163265306*uOther[14]*cMOther[16]*mnuOther-0.159719141249985*uOther[10]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[15]*mnuOther-0.223606797749979*uOther[12]*cMOther[15]*mnuOther-0.159719141249985*cMOther[12]*uOther[14]*mnuOther-0.159719141249985*uOther[12]*cMOther[14]*mnuOther-0.3928571428571429*cMOther[11]*uOther[13]*mnuOther-0.3928571428571429*uOther[11]*cMOther[13]*mnuOther-0.2500000000000001*cMOther[10]*uOther[12]*mnuOther-0.2500000000000001*uOther[10]*cMOther[12]*mnuOther-0.6388765649999399*m0rOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(54,57) = (-0.1963961012123931*cMOther[13]*uOther[19]*mnuOther)-0.1963961012123931*uOther[13]*cMOther[19]*mnuOther-0.1428571428571428*cMOther[15]*uOther[17]*mnuOther-0.2*cMOther[14]*uOther[17]*mnuOther-0.223606797749979*cMOther[10]*uOther[17]*mnuOther-0.1428571428571428*uOther[15]*cMOther[17]*mnuOther-0.2*uOther[14]*cMOther[17]*mnuOther-0.223606797749979*uOther[10]*cMOther[17]*mnuOther-0.1788854381999831*cMOther[13]*uOther[16]*mnuOther-0.1788854381999831*uOther[13]*cMOther[16]*mnuOther-0.223606797749979*cMOther[11]*uOther[15]*mnuOther-0.223606797749979*uOther[11]*cMOther[15]*mnuOther-0.2*cMOther[12]*uOther[13]*mnuOther-0.2*uOther[12]*cMOther[13]*mnuOther-0.8944271909999159*m0rOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(54,58) = (-0.2817460317460317*cMOther[14]*uOther[18]*mnuOther)-0.149071198499986*cMOther[10]*uOther[18]*mnuOther-0.2817460317460317*uOther[14]*cMOther[18]*mnuOther-0.149071198499986*uOther[10]*cMOther[18]*mnuOther-0.2195775164134199*cMOther[15]*uOther[17]*mnuOther-0.2195775164134199*uOther[15]*cMOther[17]*mnuOther-0.3273268353539885*cMOther[13]*uOther[16]*mnuOther-0.3273268353539885*uOther[13]*cMOther[16]*mnuOther-0.3273268353539885*cMOther[11]*uOther[14]*mnuOther-0.3273268353539885*uOther[11]*cMOther[14]*mnuOther-0.2195775164134199*cMOther[12]*uOther[13]*mnuOther-0.2195775164134199*uOther[12]*cMOther[13]*mnuOther-0.2195775164134199*cMOther[10]*uOther[11]*mnuOther-0.2195775164134199*uOther[10]*cMOther[11]*mnuOther-0.5962847939999438*m0rOther[8]*mnuOther+0.2981423969999719*cEOther[8]*mnuOther-0.8783100656536796*m0rOther[1]*mnuOther+0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(55,50) = (-0.149071198499986*cMOther[19]*uOther[19]*mnuOther)-0.2195775164134199*cMOther[12]*uOther[19]*mnuOther-0.2195775164134199*uOther[12]*cMOther[19]*mnuOther-0.159719141249985*cMOther[17]*uOther[17]*mnuOther-0.2500000000000001*cMOther[11]*uOther[17]*mnuOther-0.2500000000000001*uOther[11]*cMOther[17]*mnuOther-0.223606797749979*cMOther[16]*uOther[16]*mnuOther-0.159719141249985*cMOther[15]*uOther[15]*mnuOther-0.25*cMOther[10]*uOther[15]*mnuOther-0.25*uOther[10]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[12]*uOther[12]*mnuOther-1.0*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(55,51) = (-0.2195775164134199*cMOther[13]*uOther[19]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[19]*mnuOther-0.159719141249985*cMOther[15]*uOther[17]*mnuOther-0.223606797749979*cMOther[14]*uOther[17]*mnuOther-0.2500000000000001*cMOther[10]*uOther[17]*mnuOther-0.159719141249985*uOther[15]*cMOther[17]*mnuOther-0.223606797749979*uOther[14]*cMOther[17]*mnuOther-0.2500000000000001*uOther[10]*cMOther[17]*mnuOther-0.2*cMOther[13]*uOther[16]*mnuOther-0.2*uOther[13]*cMOther[16]*mnuOther-0.25*cMOther[11]*uOther[15]*mnuOther-0.25*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[12]*uOther[13]*mnuOther-0.223606797749979*uOther[12]*cMOther[13]*mnuOther-1.0*m0rOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(55,52) = (-0.3273268353539885*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*cMOther[10]*uOther[19]*mnuOther-0.3273268353539885*uOther[15]*cMOther[19]*mnuOther-0.2195775164134199*uOther[10]*cMOther[19]*mnuOther-0.3928571428571429*cMOther[13]*uOther[17]*mnuOther-0.3928571428571429*uOther[13]*cMOther[17]*mnuOther-0.223606797749979*cMOther[14]*uOther[16]*mnuOther-0.223606797749979*uOther[14]*cMOther[16]*mnuOther-0.3928571428571428*cMOther[12]*uOther[15]*mnuOther-0.3928571428571428*uOther[12]*cMOther[15]*mnuOther-0.223606797749979*cMOther[11]*uOther[13]*mnuOther-0.223606797749979*uOther[11]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[12]*mnuOther-0.223606797749979*uOther[10]*cMOther[12]*mnuOther-0.8783100656536796*m0rOther[9]*mnuOther+0.4391550328268398*cEOther[9]*mnuOther-0.8944271909999159*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(55,53) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999831*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999831*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999159*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(55,55) = (-0.3452380952380952*cMOther[19]*uOther[19]*mnuOther)-0.1402829294374236*cMOther[12]*uOther[19]*mnuOther-0.1402829294374236*uOther[12]*cMOther[19]*mnuOther-0.25*cMOther[18]*uOther[18]*mnuOther-0.3520408163265306*cMOther[17]*uOther[17]*mnuOther-0.159719141249985*cMOther[11]*uOther[17]*mnuOther-0.159719141249985*uOther[11]*cMOther[17]*mnuOther-0.3928571428571428*cMOther[16]*uOther[16]*mnuOther-0.3520408163265306*cMOther[15]*uOther[15]*mnuOther-0.159719141249985*cMOther[10]*uOther[15]*mnuOther-0.159719141249985*uOther[10]*cMOther[15]*mnuOther-0.25*cMOther[14]*uOther[14]*mnuOther-0.3928571428571428*cMOther[13]*uOther[13]*mnuOther-0.3928571428571428*cMOther[12]*uOther[12]*mnuOther-0.25*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.6388765649999399*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(55,56) = (-0.1963961012123931*cMOther[13]*uOther[18]*mnuOther)-0.1963961012123931*uOther[13]*cMOther[18]*mnuOther-0.1788854381999831*cMOther[13]*uOther[17]*mnuOther-0.1788854381999831*uOther[13]*cMOther[17]*mnuOther-0.2*cMOther[15]*uOther[16]*mnuOther-0.1428571428571428*cMOther[14]*uOther[16]*mnuOther-0.223606797749979*cMOther[10]*uOther[16]*mnuOther-0.2*uOther[15]*cMOther[16]*mnuOther-0.1428571428571428*uOther[14]*cMOther[16]*mnuOther-0.223606797749979*uOther[10]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[14]*mnuOther-0.223606797749979*uOther[12]*cMOther[14]*mnuOther-0.2*cMOther[11]*uOther[13]*mnuOther-0.2*uOther[11]*cMOther[13]*mnuOther-0.8944271909999159*m0rOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(55,57) = (-0.1402829294374237*cMOther[13]*uOther[19]*mnuOther)-0.1402829294374237*uOther[13]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[14]*uOther[18]*mnuOther-0.2195775164134199*uOther[14]*cMOther[18]*mnuOther-0.3520408163265306*cMOther[15]*uOther[17]*mnuOther-0.1428571428571428*cMOther[14]*uOther[17]*mnuOther-0.159719141249985*cMOther[10]*uOther[17]*mnuOther-0.3520408163265306*uOther[15]*cMOther[17]*mnuOther-0.1428571428571428*uOther[14]*cMOther[17]*mnuOther-0.159719141249985*uOther[10]*cMOther[17]*mnuOther-0.3513821107499669*cMOther[13]*uOther[16]*mnuOther-0.3513821107499669*uOther[13]*cMOther[16]*mnuOther-0.159719141249985*cMOther[11]*uOther[15]*mnuOther-0.159719141249985*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[11]*uOther[14]*mnuOther-0.223606797749979*uOther[11]*cMOther[14]*mnuOther-0.3928571428571429*cMOther[12]*uOther[13]*mnuOther-0.3928571428571429*uOther[12]*cMOther[13]*mnuOther-0.2500000000000001*cMOther[10]*uOther[11]*mnuOther-0.2500000000000001*uOther[10]*cMOther[11]*mnuOther-0.6388765649999399*m0rOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(55,59) = (-0.2817460317460317*cMOther[15]*uOther[19]*mnuOther)-0.149071198499986*cMOther[10]*uOther[19]*mnuOther-0.2817460317460317*uOther[15]*cMOther[19]*mnuOther-0.149071198499986*uOther[10]*cMOther[19]*mnuOther-0.3273268353539885*cMOther[13]*uOther[17]*mnuOther-0.3273268353539885*uOther[13]*cMOther[17]*mnuOther-0.2195775164134199*cMOther[14]*uOther[16]*mnuOther-0.2195775164134199*uOther[14]*cMOther[16]*mnuOther-0.3273268353539885*cMOther[12]*uOther[15]*mnuOther-0.3273268353539885*uOther[12]*cMOther[15]*mnuOther-0.2195775164134199*cMOther[11]*uOther[13]*mnuOther-0.2195775164134199*uOther[11]*cMOther[13]*mnuOther-0.2195775164134199*cMOther[10]*uOther[12]*mnuOther-0.2195775164134199*uOther[10]*cMOther[12]*mnuOther-0.5962847939999438*m0rOther[9]*mnuOther+0.2981423969999719*cEOther[9]*mnuOther-0.8783100656536796*m0rOther[2]*mnuOther+0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(56,50) = (-0.2195775164134199*cMOther[13]*uOther[18]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[18]*mnuOther-0.2*cMOther[13]*uOther[17]*mnuOther-0.2*uOther[13]*cMOther[17]*mnuOther-0.223606797749979*cMOther[15]*uOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[16]*mnuOther-0.25*cMOther[10]*uOther[16]*mnuOther-0.223606797749979*uOther[15]*cMOther[16]*mnuOther-0.159719141249985*uOther[14]*cMOther[16]*mnuOther-0.25*uOther[10]*cMOther[16]*mnuOther-0.2500000000000001*cMOther[12]*uOther[14]*mnuOther-0.2500000000000001*uOther[12]*cMOther[14]*mnuOther-0.223606797749979*cMOther[11]*uOther[13]*mnuOther-0.223606797749979*uOther[11]*cMOther[13]*mnuOther-1.0*m0rOther[6]*mnuOther+0.5*cEOther[6]*mnuOther; 
  data->AEM_S(56,51) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999832*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999832*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999161*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(56,52) = (-0.149071198499986*cMOther[18]*uOther[18]*mnuOther)-0.2195775164134199*cMOther[11]*uOther[18]*mnuOther-0.2195775164134199*uOther[11]*cMOther[18]*mnuOther-0.223606797749979*cMOther[17]*uOther[17]*mnuOther-0.159719141249985*cMOther[16]*uOther[16]*mnuOther-0.25*cMOther[12]*uOther[16]*mnuOther-0.25*uOther[12]*cMOther[16]*mnuOther-0.159719141249985*cMOther[14]*uOther[14]*mnuOther-0.2500000000000001*cMOther[10]*uOther[14]*mnuOther-0.2500000000000001*uOther[10]*cMOther[14]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[11]*mnuOther-1.0*m0rOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(56,53) = (-0.175662013130736*cMOther[13]*uOther[19]*mnuOther)-0.175662013130736*uOther[13]*cMOther[19]*mnuOther-0.3273268353539885*cMOther[14]*uOther[18]*mnuOther-0.21957751641342*cMOther[10]*uOther[18]*mnuOther-0.3273268353539885*uOther[14]*cMOther[18]*mnuOther-0.21957751641342*uOther[10]*cMOther[18]*mnuOther-0.3513821107499669*cMOther[15]*uOther[17]*mnuOther-0.1788854381999831*cMOther[14]*uOther[17]*mnuOther-0.2*cMOther[10]*uOther[17]*mnuOther-0.3513821107499669*uOther[15]*cMOther[17]*mnuOther-0.1788854381999831*uOther[14]*cMOther[17]*mnuOther-0.2*uOther[10]*cMOther[17]*mnuOther-0.5528571428571428*cMOther[13]*uOther[16]*mnuOther-0.5528571428571428*uOther[13]*cMOther[16]*mnuOther-0.2*cMOther[11]*uOther[15]*mnuOther-0.2*uOther[11]*cMOther[15]*mnuOther-0.3928571428571429*cMOther[11]*uOther[14]*mnuOther-0.3928571428571429*uOther[11]*cMOther[14]*mnuOther-0.4024922359499621*cMOther[12]*uOther[13]*mnuOther-0.4024922359499621*uOther[12]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[11]*mnuOther-0.223606797749979*uOther[10]*cMOther[11]*mnuOther-0.87831006565368*m0rOther[8]*mnuOther+0.43915503282684*cEOther[8]*mnuOther-0.8*m0rOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.8944271909999161*m0rOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(56,54) = (-0.2195775164134199*cMOther[15]*uOther[19]*mnuOther)-0.2195775164134199*uOther[15]*cMOther[19]*mnuOther-0.1402829294374237*cMOther[13]*uOther[18]*mnuOther-0.1402829294374237*uOther[13]*cMOther[18]*mnuOther-0.3513821107499669*cMOther[13]*uOther[17]*mnuOther-0.3513821107499669*uOther[13]*cMOther[17]*mnuOther-0.1428571428571428*cMOther[15]*uOther[16]*mnuOther-0.3520408163265306*cMOther[14]*uOther[16]*mnuOther-0.159719141249985*cMOther[10]*uOther[16]*mnuOther-0.1428571428571428*uOther[15]*cMOther[16]*mnuOther-0.3520408163265306*uOther[14]*cMOther[16]*mnuOther-0.159719141249985*uOther[10]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[15]*mnuOther-0.223606797749979*uOther[12]*cMOther[15]*mnuOther-0.159719141249985*cMOther[12]*uOther[14]*mnuOther-0.159719141249985*uOther[12]*cMOther[14]*mnuOther-0.3928571428571429*cMOther[11]*uOther[13]*mnuOther-0.3928571428571429*uOther[11]*cMOther[13]*mnuOther-0.2500000000000001*cMOther[10]*uOther[12]*mnuOther-0.2500000000000001*uOther[10]*cMOther[12]*mnuOther-0.6388765649999399*m0rOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-1.0*m0rOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(56,55) = (-0.1963961012123931*cMOther[13]*uOther[18]*mnuOther)-0.1963961012123931*uOther[13]*cMOther[18]*mnuOther-0.1788854381999831*cMOther[13]*uOther[17]*mnuOther-0.1788854381999831*uOther[13]*cMOther[17]*mnuOther-0.2*cMOther[15]*uOther[16]*mnuOther-0.1428571428571428*cMOther[14]*uOther[16]*mnuOther-0.223606797749979*cMOther[10]*uOther[16]*mnuOther-0.2*uOther[15]*cMOther[16]*mnuOther-0.1428571428571428*uOther[14]*cMOther[16]*mnuOther-0.223606797749979*uOther[10]*cMOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[14]*mnuOther-0.223606797749979*uOther[12]*cMOther[14]*mnuOther-0.2*cMOther[11]*uOther[13]*mnuOther-0.2*uOther[11]*cMOther[13]*mnuOther-0.8944271909999159*m0rOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(56,56) = (-0.3833333333333334*cMOther[19]*uOther[19]*mnuOther)-0.1963961012123931*cMOther[12]*uOther[19]*mnuOther-0.1963961012123931*uOther[12]*cMOther[19]*mnuOther-0.3452380952380952*cMOther[18]*uOther[18]*mnuOther-0.1402829294374236*cMOther[11]*uOther[18]*mnuOther-0.1402829294374236*uOther[11]*cMOther[18]*mnuOther-0.5357142857142857*cMOther[17]*uOther[17]*mnuOther-0.223606797749979*cMOther[11]*uOther[17]*mnuOther-0.223606797749979*uOther[11]*cMOther[17]*mnuOther-0.5520408163265306*cMOther[16]*uOther[16]*mnuOther-0.159719141249985*cMOther[12]*uOther[16]*mnuOther-0.159719141249985*uOther[12]*cMOther[16]*mnuOther-0.3928571428571428*cMOther[15]*uOther[15]*mnuOther-0.223606797749979*cMOther[10]*uOther[15]*mnuOther-0.223606797749979*uOther[10]*cMOther[15]*mnuOther-0.3520408163265306*cMOther[14]*uOther[14]*mnuOther-0.159719141249985*cMOther[10]*uOther[14]*mnuOther-0.159719141249985*uOther[10]*cMOther[14]*mnuOther-0.5928571428571429*cMOther[13]*uOther[13]*mnuOther-0.45*cMOther[12]*uOther[12]*mnuOther-0.3928571428571428*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.8944271909999159*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.6388765649999399*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(56,57) = (-0.175662013130736*cMOther[17]*uOther[19]*mnuOther)-0.175662013130736*uOther[17]*cMOther[19]*mnuOther-0.175662013130736*cMOther[16]*uOther[18]*mnuOther-0.175662013130736*uOther[16]*cMOther[18]*mnuOther-0.16*cMOther[16]*uOther[17]*mnuOther-0.1788854381999832*cMOther[12]*uOther[17]*mnuOther-0.16*uOther[16]*cMOther[17]*mnuOther-0.1788854381999832*uOther[12]*cMOther[17]*mnuOther-0.1788854381999832*cMOther[11]*uOther[16]*mnuOther-0.1788854381999832*uOther[11]*cMOther[16]*mnuOther-0.1788854381999831*cMOther[13]*uOther[15]*mnuOther-0.1788854381999831*uOther[13]*cMOther[15]*mnuOther-0.1788854381999831*cMOther[13]*uOther[14]*mnuOther-0.1788854381999831*uOther[13]*cMOther[14]*mnuOther-0.2*cMOther[10]*uOther[13]*mnuOther-0.2*uOther[10]*cMOther[13]*mnuOther-0.2*cMOther[11]*uOther[12]*mnuOther-0.2*uOther[11]*cMOther[12]*mnuOther-0.8*m0rOther[3]*mnuOther+0.4*cEOther[3]*mnuOther; 
  data->AEM_S(56,58) = (-0.1928571428571429*cMOther[17]*uOther[19]*mnuOther)-0.1928571428571429*uOther[17]*cMOther[19]*mnuOther-0.1928571428571429*cMOther[16]*uOther[18]*mnuOther-0.1928571428571429*uOther[16]*cMOther[18]*mnuOther-0.175662013130736*cMOther[16]*uOther[17]*mnuOther-0.1963961012123931*cMOther[12]*uOther[17]*mnuOther-0.175662013130736*uOther[16]*cMOther[17]*mnuOther-0.1963961012123931*uOther[12]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[11]*uOther[16]*mnuOther-0.1963961012123931*uOther[11]*cMOther[16]*mnuOther-0.1963961012123931*cMOther[13]*uOther[15]*mnuOther-0.1963961012123931*uOther[13]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[13]*uOther[14]*mnuOther-0.1963961012123931*uOther[13]*cMOther[14]*mnuOther-0.21957751641342*cMOther[10]*uOther[13]*mnuOther-0.21957751641342*uOther[10]*cMOther[13]*mnuOther-0.21957751641342*cMOther[11]*uOther[12]*mnuOther-0.21957751641342*uOther[11]*cMOther[12]*mnuOther-0.87831006565368*m0rOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther; 
  data->AEM_S(57,50) = (-0.2195775164134199*cMOther[13]*uOther[19]*mnuOther)-0.2195775164134199*uOther[13]*cMOther[19]*mnuOther-0.159719141249985*cMOther[15]*uOther[17]*mnuOther-0.223606797749979*cMOther[14]*uOther[17]*mnuOther-0.25*cMOther[10]*uOther[17]*mnuOther-0.159719141249985*uOther[15]*cMOther[17]*mnuOther-0.223606797749979*uOther[14]*cMOther[17]*mnuOther-0.25*uOther[10]*cMOther[17]*mnuOther-0.2*cMOther[13]*uOther[16]*mnuOther-0.2*uOther[13]*cMOther[16]*mnuOther-0.2500000000000001*cMOther[11]*uOther[15]*mnuOther-0.2500000000000001*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[12]*uOther[13]*mnuOther-0.223606797749979*uOther[12]*cMOther[13]*mnuOther-1.0*m0rOther[7]*mnuOther+0.5*cEOther[7]*mnuOther; 
  data->AEM_S(57,51) = (-0.149071198499986*cMOther[19]*uOther[19]*mnuOther)-0.2195775164134199*cMOther[12]*uOther[19]*mnuOther-0.2195775164134199*uOther[12]*cMOther[19]*mnuOther-0.159719141249985*cMOther[17]*uOther[17]*mnuOther-0.25*cMOther[11]*uOther[17]*mnuOther-0.25*uOther[11]*cMOther[17]*mnuOther-0.223606797749979*cMOther[16]*uOther[16]*mnuOther-0.159719141249985*cMOther[15]*uOther[15]*mnuOther-0.2500000000000001*cMOther[10]*uOther[15]*mnuOther-0.2500000000000001*uOther[10]*cMOther[15]*mnuOther-0.223606797749979*cMOther[13]*uOther[13]*mnuOther-0.223606797749979*cMOther[12]*uOther[12]*mnuOther-1.0*m0rOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(57,52) = (-0.1963961012123931*cMOther[17]*uOther[19]*mnuOther)-0.1963961012123931*uOther[17]*cMOther[19]*mnuOther-0.1963961012123931*cMOther[16]*uOther[18]*mnuOther-0.1963961012123931*uOther[16]*cMOther[18]*mnuOther-0.1788854381999832*cMOther[16]*uOther[17]*mnuOther-0.2*cMOther[12]*uOther[17]*mnuOther-0.1788854381999832*uOther[16]*cMOther[17]*mnuOther-0.2*uOther[12]*cMOther[17]*mnuOther-0.2*cMOther[11]*uOther[16]*mnuOther-0.2*uOther[11]*cMOther[16]*mnuOther-0.2*cMOther[13]*uOther[15]*mnuOther-0.2*uOther[13]*cMOther[15]*mnuOther-0.2*cMOther[13]*uOther[14]*mnuOther-0.2*uOther[13]*cMOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[13]*mnuOther-0.223606797749979*uOther[10]*cMOther[13]*mnuOther-0.223606797749979*cMOther[11]*uOther[12]*mnuOther-0.223606797749979*uOther[11]*cMOther[12]*mnuOther-0.8944271909999161*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(57,53) = (-0.3273268353539885*cMOther[15]*uOther[19]*mnuOther)-0.21957751641342*cMOther[10]*uOther[19]*mnuOther-0.3273268353539885*uOther[15]*cMOther[19]*mnuOther-0.21957751641342*uOther[10]*cMOther[19]*mnuOther-0.175662013130736*cMOther[13]*uOther[18]*mnuOther-0.175662013130736*uOther[13]*cMOther[18]*mnuOther-0.5528571428571428*cMOther[13]*uOther[17]*mnuOther-0.5528571428571428*uOther[13]*cMOther[17]*mnuOther-0.1788854381999831*cMOther[15]*uOther[16]*mnuOther-0.3513821107499669*cMOther[14]*uOther[16]*mnuOther-0.2*cMOther[10]*uOther[16]*mnuOther-0.1788854381999831*uOther[15]*cMOther[16]*mnuOther-0.3513821107499669*uOther[14]*cMOther[16]*mnuOther-0.2*uOther[10]*cMOther[16]*mnuOther-0.3928571428571429*cMOther[12]*uOther[15]*mnuOther-0.3928571428571429*uOther[12]*cMOther[15]*mnuOther-0.2*cMOther[12]*uOther[14]*mnuOther-0.2*uOther[12]*cMOther[14]*mnuOther-0.4024922359499621*cMOther[11]*uOther[13]*mnuOther-0.4024922359499621*uOther[11]*cMOther[13]*mnuOther-0.223606797749979*cMOther[10]*uOther[12]*mnuOther-0.223606797749979*uOther[10]*cMOther[12]*mnuOther-0.87831006565368*m0rOther[9]*mnuOther+0.43915503282684*cEOther[9]*mnuOther-0.8*m0rOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.8944271909999161*m0rOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(57,54) = (-0.1963961012123931*cMOther[13]*uOther[19]*mnuOther)-0.1963961012123931*uOther[13]*cMOther[19]*mnuOther-0.1428571428571428*cMOther[15]*uOther[17]*mnuOther-0.2*cMOther[14]*uOther[17]*mnuOther-0.223606797749979*cMOther[10]*uOther[17]*mnuOther-0.1428571428571428*uOther[15]*cMOther[17]*mnuOther-0.2*uOther[14]*cMOther[17]*mnuOther-0.223606797749979*uOther[10]*cMOther[17]*mnuOther-0.1788854381999831*cMOther[13]*uOther[16]*mnuOther-0.1788854381999831*uOther[13]*cMOther[16]*mnuOther-0.223606797749979*cMOther[11]*uOther[15]*mnuOther-0.223606797749979*uOther[11]*cMOther[15]*mnuOther-0.2*cMOther[12]*uOther[13]*mnuOther-0.2*uOther[12]*cMOther[13]*mnuOther-0.8944271909999159*m0rOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(57,55) = (-0.1402829294374237*cMOther[13]*uOther[19]*mnuOther)-0.1402829294374237*uOther[13]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[14]*uOther[18]*mnuOther-0.2195775164134199*uOther[14]*cMOther[18]*mnuOther-0.3520408163265306*cMOther[15]*uOther[17]*mnuOther-0.1428571428571428*cMOther[14]*uOther[17]*mnuOther-0.159719141249985*cMOther[10]*uOther[17]*mnuOther-0.3520408163265306*uOther[15]*cMOther[17]*mnuOther-0.1428571428571428*uOther[14]*cMOther[17]*mnuOther-0.159719141249985*uOther[10]*cMOther[17]*mnuOther-0.3513821107499669*cMOther[13]*uOther[16]*mnuOther-0.3513821107499669*uOther[13]*cMOther[16]*mnuOther-0.159719141249985*cMOther[11]*uOther[15]*mnuOther-0.159719141249985*uOther[11]*cMOther[15]*mnuOther-0.223606797749979*cMOther[11]*uOther[14]*mnuOther-0.223606797749979*uOther[11]*cMOther[14]*mnuOther-0.3928571428571429*cMOther[12]*uOther[13]*mnuOther-0.3928571428571429*uOther[12]*cMOther[13]*mnuOther-0.2500000000000001*cMOther[10]*uOther[11]*mnuOther-0.2500000000000001*uOther[10]*cMOther[11]*mnuOther-0.6388765649999399*m0rOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-1.0*m0rOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(57,56) = (-0.175662013130736*cMOther[17]*uOther[19]*mnuOther)-0.175662013130736*uOther[17]*cMOther[19]*mnuOther-0.175662013130736*cMOther[16]*uOther[18]*mnuOther-0.175662013130736*uOther[16]*cMOther[18]*mnuOther-0.16*cMOther[16]*uOther[17]*mnuOther-0.1788854381999832*cMOther[12]*uOther[17]*mnuOther-0.16*uOther[16]*cMOther[17]*mnuOther-0.1788854381999832*uOther[12]*cMOther[17]*mnuOther-0.1788854381999832*cMOther[11]*uOther[16]*mnuOther-0.1788854381999832*uOther[11]*cMOther[16]*mnuOther-0.1788854381999831*cMOther[13]*uOther[15]*mnuOther-0.1788854381999831*uOther[13]*cMOther[15]*mnuOther-0.1788854381999831*cMOther[13]*uOther[14]*mnuOther-0.1788854381999831*uOther[13]*cMOther[14]*mnuOther-0.2*cMOther[10]*uOther[13]*mnuOther-0.2*uOther[10]*cMOther[13]*mnuOther-0.2*cMOther[11]*uOther[12]*mnuOther-0.2*uOther[11]*cMOther[12]*mnuOther-0.8*m0rOther[3]*mnuOther+0.4*cEOther[3]*mnuOther; 
  data->AEM_S(57,57) = (-0.3452380952380952*cMOther[19]*uOther[19]*mnuOther)-0.1402829294374236*cMOther[12]*uOther[19]*mnuOther-0.1402829294374236*uOther[12]*cMOther[19]*mnuOther-0.3833333333333334*cMOther[18]*uOther[18]*mnuOther-0.1963961012123931*cMOther[11]*uOther[18]*mnuOther-0.1963961012123931*uOther[11]*cMOther[18]*mnuOther-0.5520408163265306*cMOther[17]*uOther[17]*mnuOther-0.159719141249985*cMOther[11]*uOther[17]*mnuOther-0.159719141249985*uOther[11]*cMOther[17]*mnuOther-0.5357142857142857*cMOther[16]*uOther[16]*mnuOther-0.223606797749979*cMOther[12]*uOther[16]*mnuOther-0.223606797749979*uOther[12]*cMOther[16]*mnuOther-0.3520408163265306*cMOther[15]*uOther[15]*mnuOther-0.159719141249985*cMOther[10]*uOther[15]*mnuOther-0.159719141249985*uOther[10]*cMOther[15]*mnuOther-0.3928571428571428*cMOther[14]*uOther[14]*mnuOther-0.223606797749979*cMOther[10]*uOther[14]*mnuOther-0.223606797749979*uOther[10]*cMOther[14]*mnuOther-0.5928571428571429*cMOther[13]*uOther[13]*mnuOther-0.3928571428571428*cMOther[12]*uOther[12]*mnuOther-0.45*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.6388765649999399*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.8944271909999159*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(57,59) = (-0.1928571428571429*cMOther[17]*uOther[19]*mnuOther)-0.1928571428571429*uOther[17]*cMOther[19]*mnuOther-0.1928571428571429*cMOther[16]*uOther[18]*mnuOther-0.1928571428571429*uOther[16]*cMOther[18]*mnuOther-0.175662013130736*cMOther[16]*uOther[17]*mnuOther-0.1963961012123931*cMOther[12]*uOther[17]*mnuOther-0.175662013130736*uOther[16]*cMOther[17]*mnuOther-0.1963961012123931*uOther[12]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[11]*uOther[16]*mnuOther-0.1963961012123931*uOther[11]*cMOther[16]*mnuOther-0.1963961012123931*cMOther[13]*uOther[15]*mnuOther-0.1963961012123931*uOther[13]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[13]*uOther[14]*mnuOther-0.1963961012123931*uOther[13]*cMOther[14]*mnuOther-0.21957751641342*cMOther[10]*uOther[13]*mnuOther-0.21957751641342*uOther[10]*cMOther[13]*mnuOther-0.21957751641342*cMOther[11]*uOther[12]*mnuOther-0.21957751641342*uOther[11]*cMOther[12]*mnuOther-0.87831006565368*m0rOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther; 
  data->AEM_S(58,50) = (-0.149071198499986*cMOther[14]*uOther[18]*mnuOther)-0.25*cMOther[10]*uOther[18]*mnuOther-0.149071198499986*uOther[14]*cMOther[18]*mnuOther-0.25*uOther[10]*cMOther[18]*mnuOther-0.2195775164134199*cMOther[13]*uOther[16]*mnuOther-0.2195775164134199*uOther[13]*cMOther[16]*mnuOther-0.2195775164134199*cMOther[11]*uOther[14]*mnuOther-0.2195775164134199*uOther[11]*cMOther[14]*mnuOther-1.0*m0rOther[8]*mnuOther+0.5*cEOther[8]*mnuOther; 
  data->AEM_S(58,51) = (-0.1309307341415954*cMOther[18]*uOther[18]*mnuOther)-0.1928571428571429*cMOther[11]*uOther[18]*mnuOther-0.1928571428571429*uOther[11]*cMOther[18]*mnuOther-0.1963961012123931*cMOther[17]*uOther[17]*mnuOther-0.1402829294374236*cMOther[16]*uOther[16]*mnuOther-0.2195775164134199*cMOther[12]*uOther[16]*mnuOther-0.2195775164134199*uOther[12]*cMOther[16]*mnuOther-0.1402829294374236*cMOther[14]*uOther[14]*mnuOther-0.2195775164134199*cMOther[10]*uOther[14]*mnuOther-0.2195775164134199*uOther[10]*cMOther[14]*mnuOther-0.1963961012123931*cMOther[13]*uOther[13]*mnuOther-0.1963961012123931*cMOther[11]*uOther[11]*mnuOther-0.8783100656536796*m0rOther[4]*mnuOther+0.4391550328268398*cEOther[4]*mnuOther; 
  data->AEM_S(58,53) = (-0.1928571428571429*cMOther[13]*uOther[18]*mnuOther)-0.1928571428571429*uOther[13]*cMOther[18]*mnuOther-0.175662013130736*cMOther[13]*uOther[17]*mnuOther-0.175662013130736*uOther[13]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[15]*uOther[16]*mnuOther-0.1402829294374237*cMOther[14]*uOther[16]*mnuOther-0.21957751641342*cMOther[10]*uOther[16]*mnuOther-0.1963961012123931*uOther[15]*cMOther[16]*mnuOther-0.1402829294374237*uOther[14]*cMOther[16]*mnuOther-0.21957751641342*uOther[10]*cMOther[16]*mnuOther-0.2195775164134199*cMOther[12]*uOther[14]*mnuOther-0.2195775164134199*uOther[12]*cMOther[14]*mnuOther-0.1963961012123931*cMOther[11]*uOther[13]*mnuOther-0.1963961012123931*uOther[11]*cMOther[13]*mnuOther-0.87831006565368*m0rOther[6]*mnuOther+0.43915503282684*cEOther[6]*mnuOther; 
  data->AEM_S(58,54) = (-0.2817460317460317*cMOther[14]*uOther[18]*mnuOther)-0.149071198499986*cMOther[10]*uOther[18]*mnuOther-0.2817460317460317*uOther[14]*cMOther[18]*mnuOther-0.149071198499986*uOther[10]*cMOther[18]*mnuOther-0.2195775164134199*cMOther[15]*uOther[17]*mnuOther-0.2195775164134199*uOther[15]*cMOther[17]*mnuOther-0.3273268353539885*cMOther[13]*uOther[16]*mnuOther-0.3273268353539885*uOther[13]*cMOther[16]*mnuOther-0.3273268353539885*cMOther[11]*uOther[14]*mnuOther-0.3273268353539885*uOther[11]*cMOther[14]*mnuOther-0.2195775164134199*cMOther[12]*uOther[13]*mnuOther-0.2195775164134199*uOther[12]*cMOther[13]*mnuOther-0.2195775164134199*cMOther[10]*uOther[11]*mnuOther-0.2195775164134199*uOther[10]*cMOther[11]*mnuOther-0.5962847939999438*m0rOther[8]*mnuOther+0.2981423969999719*cEOther[8]*mnuOther-0.8783100656536796*m0rOther[1]*mnuOther+0.4391550328268398*cEOther[1]*mnuOther; 
  data->AEM_S(58,56) = (-0.1928571428571429*cMOther[17]*uOther[19]*mnuOther)-0.1928571428571429*uOther[17]*cMOther[19]*mnuOther-0.1928571428571429*cMOther[16]*uOther[18]*mnuOther-0.1928571428571429*uOther[16]*cMOther[18]*mnuOther-0.175662013130736*cMOther[16]*uOther[17]*mnuOther-0.1963961012123931*cMOther[12]*uOther[17]*mnuOther-0.175662013130736*uOther[16]*cMOther[17]*mnuOther-0.1963961012123931*uOther[12]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[11]*uOther[16]*mnuOther-0.1963961012123931*uOther[11]*cMOther[16]*mnuOther-0.1963961012123931*cMOther[13]*uOther[15]*mnuOther-0.1963961012123931*uOther[13]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[13]*uOther[14]*mnuOther-0.1963961012123931*uOther[13]*cMOther[14]*mnuOther-0.21957751641342*cMOther[10]*uOther[13]*mnuOther-0.21957751641342*uOther[10]*cMOther[13]*mnuOther-0.21957751641342*cMOther[11]*uOther[12]*mnuOther-0.21957751641342*uOther[11]*cMOther[12]*mnuOther-0.87831006565368*m0rOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther; 
  data->AEM_S(58,58) = (-0.25*cMOther[19]*uOther[19]*mnuOther)-0.3388888888888889*cMOther[18]*uOther[18]*mnuOther-0.1309307341415954*cMOther[11]*uOther[18]*mnuOther-0.1309307341415954*uOther[11]*cMOther[18]*mnuOther-0.3833333333333334*cMOther[17]*uOther[17]*mnuOther-0.3452380952380952*cMOther[16]*uOther[16]*mnuOther-0.149071198499986*cMOther[12]*uOther[16]*mnuOther-0.149071198499986*uOther[12]*cMOther[16]*mnuOther-0.25*cMOther[15]*uOther[15]*mnuOther-0.3452380952380952*cMOther[14]*uOther[14]*mnuOther-0.149071198499986*cMOther[10]*uOther[14]*mnuOther-0.149071198499986*uOther[10]*cMOther[14]*mnuOther-0.3833333333333334*cMOther[13]*uOther[13]*mnuOther-0.25*cMOther[12]*uOther[12]*mnuOther-0.3833333333333334*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.5962847939999438*m0rOther[4]*mnuOther+0.2981423969999719*cEOther[4]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(59,50) = (-0.149071198499986*cMOther[15]*uOther[19]*mnuOther)-0.25*cMOther[10]*uOther[19]*mnuOther-0.149071198499986*uOther[15]*cMOther[19]*mnuOther-0.25*uOther[10]*cMOther[19]*mnuOther-0.2195775164134199*cMOther[13]*uOther[17]*mnuOther-0.2195775164134199*uOther[13]*cMOther[17]*mnuOther-0.2195775164134199*cMOther[12]*uOther[15]*mnuOther-0.2195775164134199*uOther[12]*cMOther[15]*mnuOther-1.0*m0rOther[9]*mnuOther+0.5*cEOther[9]*mnuOther; 
  data->AEM_S(59,52) = (-0.1309307341415954*cMOther[19]*uOther[19]*mnuOther)-0.1928571428571429*cMOther[12]*uOther[19]*mnuOther-0.1928571428571429*uOther[12]*cMOther[19]*mnuOther-0.1402829294374236*cMOther[17]*uOther[17]*mnuOther-0.2195775164134199*cMOther[11]*uOther[17]*mnuOther-0.2195775164134199*uOther[11]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[16]*uOther[16]*mnuOther-0.1402829294374236*cMOther[15]*uOther[15]*mnuOther-0.2195775164134199*cMOther[10]*uOther[15]*mnuOther-0.2195775164134199*uOther[10]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[13]*uOther[13]*mnuOther-0.1963961012123931*cMOther[12]*uOther[12]*mnuOther-0.8783100656536796*m0rOther[5]*mnuOther+0.4391550328268398*cEOther[5]*mnuOther; 
  data->AEM_S(59,53) = (-0.1928571428571429*cMOther[13]*uOther[19]*mnuOther)-0.1928571428571429*uOther[13]*cMOther[19]*mnuOther-0.1402829294374237*cMOther[15]*uOther[17]*mnuOther-0.1963961012123931*cMOther[14]*uOther[17]*mnuOther-0.21957751641342*cMOther[10]*uOther[17]*mnuOther-0.1402829294374237*uOther[15]*cMOther[17]*mnuOther-0.1963961012123931*uOther[14]*cMOther[17]*mnuOther-0.21957751641342*uOther[10]*cMOther[17]*mnuOther-0.175662013130736*cMOther[13]*uOther[16]*mnuOther-0.175662013130736*uOther[13]*cMOther[16]*mnuOther-0.2195775164134199*cMOther[11]*uOther[15]*mnuOther-0.2195775164134199*uOther[11]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[12]*uOther[13]*mnuOther-0.1963961012123931*uOther[12]*cMOther[13]*mnuOther-0.87831006565368*m0rOther[7]*mnuOther+0.43915503282684*cEOther[7]*mnuOther; 
  data->AEM_S(59,55) = (-0.2817460317460317*cMOther[15]*uOther[19]*mnuOther)-0.149071198499986*cMOther[10]*uOther[19]*mnuOther-0.2817460317460317*uOther[15]*cMOther[19]*mnuOther-0.149071198499986*uOther[10]*cMOther[19]*mnuOther-0.3273268353539885*cMOther[13]*uOther[17]*mnuOther-0.3273268353539885*uOther[13]*cMOther[17]*mnuOther-0.2195775164134199*cMOther[14]*uOther[16]*mnuOther-0.2195775164134199*uOther[14]*cMOther[16]*mnuOther-0.3273268353539885*cMOther[12]*uOther[15]*mnuOther-0.3273268353539885*uOther[12]*cMOther[15]*mnuOther-0.2195775164134199*cMOther[11]*uOther[13]*mnuOther-0.2195775164134199*uOther[11]*cMOther[13]*mnuOther-0.2195775164134199*cMOther[10]*uOther[12]*mnuOther-0.2195775164134199*uOther[10]*cMOther[12]*mnuOther-0.5962847939999438*m0rOther[9]*mnuOther+0.2981423969999719*cEOther[9]*mnuOther-0.8783100656536796*m0rOther[2]*mnuOther+0.4391550328268398*cEOther[2]*mnuOther; 
  data->AEM_S(59,57) = (-0.1928571428571429*cMOther[17]*uOther[19]*mnuOther)-0.1928571428571429*uOther[17]*cMOther[19]*mnuOther-0.1928571428571429*cMOther[16]*uOther[18]*mnuOther-0.1928571428571429*uOther[16]*cMOther[18]*mnuOther-0.175662013130736*cMOther[16]*uOther[17]*mnuOther-0.1963961012123931*cMOther[12]*uOther[17]*mnuOther-0.175662013130736*uOther[16]*cMOther[17]*mnuOther-0.1963961012123931*uOther[12]*cMOther[17]*mnuOther-0.1963961012123931*cMOther[11]*uOther[16]*mnuOther-0.1963961012123931*uOther[11]*cMOther[16]*mnuOther-0.1963961012123931*cMOther[13]*uOther[15]*mnuOther-0.1963961012123931*uOther[13]*cMOther[15]*mnuOther-0.1963961012123931*cMOther[13]*uOther[14]*mnuOther-0.1963961012123931*uOther[13]*cMOther[14]*mnuOther-0.21957751641342*cMOther[10]*uOther[13]*mnuOther-0.21957751641342*uOther[10]*cMOther[13]*mnuOther-0.21957751641342*cMOther[11]*uOther[12]*mnuOther-0.21957751641342*uOther[11]*cMOther[12]*mnuOther-0.87831006565368*m0rOther[3]*mnuOther+0.43915503282684*cEOther[3]*mnuOther; 
  data->AEM_S(59,59) = (-0.3388888888888889*cMOther[19]*uOther[19]*mnuOther)-0.1309307341415954*cMOther[12]*uOther[19]*mnuOther-0.1309307341415954*uOther[12]*cMOther[19]*mnuOther-0.25*cMOther[18]*uOther[18]*mnuOther-0.3452380952380952*cMOther[17]*uOther[17]*mnuOther-0.149071198499986*cMOther[11]*uOther[17]*mnuOther-0.149071198499986*uOther[11]*cMOther[17]*mnuOther-0.3833333333333334*cMOther[16]*uOther[16]*mnuOther-0.3452380952380952*cMOther[15]*uOther[15]*mnuOther-0.149071198499986*cMOther[10]*uOther[15]*mnuOther-0.149071198499986*uOther[10]*cMOther[15]*mnuOther-0.25*cMOther[14]*uOther[14]*mnuOther-0.3833333333333334*cMOther[13]*uOther[13]*mnuOther-0.3833333333333334*cMOther[12]*uOther[12]*mnuOther-0.25*cMOther[11]*uOther[11]*mnuOther-0.25*cMOther[10]*uOther[10]*mnuOther-0.5962847939999438*m0rOther[5]*mnuOther+0.2981423969999719*cEOther[5]*mnuOther-1.0*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[10]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  for (unsigned short int vd=0; vd<2; vd++) 
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
 
  // Set other entries to 0. // 
  data->AEM_S.block<10,10>(30,10).setZero(); 
  data->AEM_S.block<10,10>(40,0).setZero(); 
  data->AEM_S.block<10,10>(30,40).setZero(); 
  data->AEM_S.block<10,10>(40,30).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM1sum[10],mnuM1sum[11],mnuM1sum[12],mnuM1sum[13],mnuM1sum[14],mnuM1sum[15],mnuM1sum[16],mnuM1sum[17],mnuM1sum[18],mnuM1sum[19],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],mnuM2sum[6],mnuM2sum[7],mnuM2sum[8],mnuM2sum[9],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m1Relax[9],m1Relax[10],m1Relax[11],m1Relax[12],m1Relax[13],m1Relax[14],m1Relax[15],m1Relax[16],m1Relax[17],m1Relax[18],m1Relax[19],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5],m2Relax[6],m2Relax[7],m2Relax[8],m2Relax[9]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,20,1) = data->u_S.segment<20>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,10,1) = data->u_S.segment<10>(20); 
 
  Eigen::Map<VectorXd>(uCrossOther,20,1) = data->u_S.segment<20>(30); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,10,1) = data->u_S.segment<10>(50); 
 
} 
 
