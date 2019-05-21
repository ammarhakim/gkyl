#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void LBOCrossPrimMoments3x3vMax_P1(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Self[3])-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[4]; 
  double m1rSelf[12]; 
  double m2rSelf[4]; 
  double m0SrSelf[4]; 
  double m1SrSelf[12]; 
  double m2SrSelf[4]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m0rSelf[3] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    m1rSelf[2] = 0.0; 
    m1rSelf[3] = 0.0; 
    m1rSelf[4] = m1Self[4]; 
    m1rSelf[5] = 0.0; 
    m1rSelf[6] = 0.0; 
    m1rSelf[7] = 0.0; 
    m1rSelf[8] = m1Self[8]; 
    m1rSelf[9] = 0.0; 
    m1rSelf[10] = 0.0; 
    m1rSelf[11] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    m2rSelf[3] = 0.0; 
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = 0.0; 
    m0SrSelf[2] = 0.0; 
    m0SrSelf[3] = 0.0; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = 0.0; 
    m1SrSelf[2] = 0.0; 
    m1SrSelf[3] = 0.0; 
    m1SrSelf[4] = m1SSelf[4]; 
    m1SrSelf[5] = 0.0; 
    m1SrSelf[6] = 0.0; 
    m1SrSelf[7] = 0.0; 
    m1SrSelf[8] = m1SSelf[8]; 
    m1SrSelf[9] = 0.0; 
    m1SrSelf[10] = 0.0; 
    m1SrSelf[11] = 0.0; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = 0.0; 
    m2SrSelf[2] = 0.0; 
    m2SrSelf[3] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
    m0rSelf[3] = m0Self[3]; 
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
    m0SrSelf[0] = m0SSelf[0]; 
    m0SrSelf[1] = m0SSelf[1]; 
    m0SrSelf[2] = m0SSelf[2]; 
    m0SrSelf[3] = m0SSelf[3]; 
    m1SrSelf[0] = m1SSelf[0]; 
    m1SrSelf[1] = m1SSelf[1]; 
    m1SrSelf[2] = m1SSelf[2]; 
    m1SrSelf[3] = m1SSelf[3]; 
    m1SrSelf[4] = m1SSelf[4]; 
    m1SrSelf[5] = m1SSelf[5]; 
    m1SrSelf[6] = m1SSelf[6]; 
    m1SrSelf[7] = m1SSelf[7]; 
    m1SrSelf[8] = m1SSelf[8]; 
    m1SrSelf[9] = m1SSelf[9]; 
    m1SrSelf[10] = m1SSelf[10]; 
    m1SrSelf[11] = m1SSelf[11]; 
    m2SrSelf[0] = m2SSelf[0]; 
    m2SrSelf[1] = m2SSelf[1]; 
    m2SrSelf[2] = m2SSelf[2]; 
    m2SrSelf[3] = m2SSelf[3]; 
  } 
 
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.6123724356957944*m0Other[3])-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[4]; 
  double m1rOther[12]; 
  double m2rOther[4]; 
  double m0SrOther[4]; 
  double m1SrOther[12]; 
  double m2SrOther[4]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m0rOther[3] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    m1rOther[2] = 0.0; 
    m1rOther[3] = 0.0; 
    m1rOther[4] = m1Other[4]; 
    m1rOther[5] = 0.0; 
    m1rOther[6] = 0.0; 
    m1rOther[7] = 0.0; 
    m1rOther[8] = m1Other[8]; 
    m1rOther[9] = 0.0; 
    m1rOther[10] = 0.0; 
    m1rOther[11] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    m2rOther[3] = 0.0; 
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = 0.0; 
    m0SrOther[2] = 0.0; 
    m0SrOther[3] = 0.0; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = 0.0; 
    m1SrOther[2] = 0.0; 
    m1SrOther[3] = 0.0; 
    m1SrOther[4] = m1SOther[4]; 
    m1SrOther[5] = 0.0; 
    m1SrOther[6] = 0.0; 
    m1SrOther[7] = 0.0; 
    m1SrOther[8] = m1SOther[8]; 
    m1SrOther[9] = 0.0; 
    m1SrOther[10] = 0.0; 
    m1SrOther[11] = 0.0; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = 0.0; 
    m2SrOther[2] = 0.0; 
    m2SrOther[3] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
    m0rOther[3] = m0Other[3]; 
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
    m0SrOther[0] = m0SOther[0]; 
    m0SrOther[1] = m0SOther[1]; 
    m0SrOther[2] = m0SOther[2]; 
    m0SrOther[3] = m0SOther[3]; 
    m1SrOther[0] = m1SOther[0]; 
    m1SrOther[1] = m1SOther[1]; 
    m1SrOther[2] = m1SOther[2]; 
    m1SrOther[3] = m1SOther[3]; 
    m1SrOther[4] = m1SOther[4]; 
    m1SrOther[5] = m1SOther[5]; 
    m1SrOther[6] = m1SOther[6]; 
    m1SrOther[7] = m1SOther[7]; 
    m1SrOther[8] = m1SOther[8]; 
    m1SrOther[9] = m1SOther[9]; 
    m1SrOther[10] = m1SOther[10]; 
    m1SrOther[11] = m1SOther[11]; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
    m2SrOther[2] = m2SOther[2]; 
    m2SrOther[3] = m2SOther[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak system. 
  data->AEM_S = Eigen::MatrixXd::Zero(32,32); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double mnuM1sum[12]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<12; vd++) 
  { 
    mnuM1sum[vd] = 0.0; 
  } 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(0,0) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(0,2) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(0,3) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,0) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,0) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,2) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(3,0) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,3) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,12) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,13) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,14) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,15) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,12) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,13) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,12) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,14) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,12) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,15) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,16) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(0,17) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(0,18) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(0,19) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(1,16) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(1,17) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(2,16) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(2,18) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(3,16) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(3,19) = 0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,28) = -0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(0,29) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(0,30) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(0,31) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(1,28) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(1,29) = -0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(2,28) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(2,30) = -0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(3,28) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(3,31) = -0.3535533905932737*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(12,0) = 0.3535533905932737*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(12,1) = 0.3535533905932737*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = 0.3535533905932737*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(12,3) = 0.3535533905932737*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(13,0) = 0.3535533905932737*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(13,1) = 0.3535533905932737*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(14,0) = 0.3535533905932737*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(14,2) = 0.3535533905932737*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(15,0) = 0.3535533905932737*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(15,3) = 0.3535533905932737*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(12,16) = 0.3535533905932737*m1SrOther[0]*mnuOther; 
  data->AEM_S(12,17) = 0.3535533905932737*m1SrOther[1]*mnuOther; 
  data->AEM_S(12,18) = 0.3535533905932737*m1SrOther[2]*mnuOther; 
  data->AEM_S(12,19) = 0.3535533905932737*m1SrOther[3]*mnuOther; 
  data->AEM_S(13,16) = 0.3535533905932737*m1SrOther[1]*mnuOther; 
  data->AEM_S(13,17) = 0.3535533905932737*m1SrOther[0]*mnuOther; 
  data->AEM_S(14,16) = 0.3535533905932737*m1SrOther[2]*mnuOther; 
  data->AEM_S(14,18) = 0.3535533905932737*m1SrOther[0]*mnuOther; 
  data->AEM_S(15,16) = 0.3535533905932737*m1SrOther[3]*mnuOther; 
  data->AEM_S(15,19) = 0.3535533905932737*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(4,4) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,5) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,6) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,7) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(5,4) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,5) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,4) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,6) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(7,4) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,7) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(4,12) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,13) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(4,14) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(4,15) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,12) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,13) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(6,12) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,14) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(7,12) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,15) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(4,20) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(4,21) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(4,22) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(4,23) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(5,20) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(5,21) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(6,20) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(6,22) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(7,20) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(7,23) = 0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(4,28) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(4,29) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(4,30) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(4,31) = -0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(5,28) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(5,29) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(6,28) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(6,30) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(7,28) = -0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(7,31) = -0.3535533905932737*cMOther[4]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(12,4) = 0.3535533905932737*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(12,5) = 0.3535533905932737*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(12,6) = 0.3535533905932737*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(12,7) = 0.3535533905932737*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(13,4) = 0.3535533905932737*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(13,5) = 0.3535533905932737*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(14,4) = 0.3535533905932737*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(14,6) = 0.3535533905932737*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(15,4) = 0.3535533905932737*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(15,7) = 0.3535533905932737*m1SrSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(12,20) = 0.3535533905932737*m1SrOther[4]*mnuOther; 
  data->AEM_S(12,21) = 0.3535533905932737*m1SrOther[5]*mnuOther; 
  data->AEM_S(12,22) = 0.3535533905932737*m1SrOther[6]*mnuOther; 
  data->AEM_S(12,23) = 0.3535533905932737*m1SrOther[7]*mnuOther; 
  data->AEM_S(13,20) = 0.3535533905932737*m1SrOther[5]*mnuOther; 
  data->AEM_S(13,21) = 0.3535533905932737*m1SrOther[4]*mnuOther; 
  data->AEM_S(14,20) = 0.3535533905932737*m1SrOther[6]*mnuOther; 
  data->AEM_S(14,22) = 0.3535533905932737*m1SrOther[4]*mnuOther; 
  data->AEM_S(15,20) = 0.3535533905932737*m1SrOther[7]*mnuOther; 
  data->AEM_S(15,23) = 0.3535533905932737*m1SrOther[4]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(8,8) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,9) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,10) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,11) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,8) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,9) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(10,8) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,10) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(11,8) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,11) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(8,12) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,13) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(8,14) = -0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(8,15) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(9,12) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,13) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(10,12) = -0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,14) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(11,12) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,15) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(8,24) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(8,25) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(8,26) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(8,27) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(9,24) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(9,25) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(10,24) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(10,26) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(11,24) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(11,27) = 0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(8,28) = -0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(8,29) = -0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(8,30) = -0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(8,31) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(9,28) = -0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(9,29) = -0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(10,28) = -0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(10,30) = -0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(11,28) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(11,31) = -0.3535533905932737*cMOther[8]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfZ and uCrossSelfZ ... // 
  data->AEM_S(12,8) = 0.3535533905932737*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(12,9) = 0.3535533905932737*m1SrSelf[9]*mnuSelf; 
  data->AEM_S(12,10) = 0.3535533905932737*m1SrSelf[10]*mnuSelf; 
  data->AEM_S(12,11) = 0.3535533905932737*m1SrSelf[11]*mnuSelf; 
  data->AEM_S(13,8) = 0.3535533905932737*m1SrSelf[9]*mnuSelf; 
  data->AEM_S(13,9) = 0.3535533905932737*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(14,8) = 0.3535533905932737*m1SrSelf[10]*mnuSelf; 
  data->AEM_S(14,10) = 0.3535533905932737*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(15,8) = 0.3535533905932737*m1SrSelf[11]*mnuSelf; 
  data->AEM_S(15,11) = 0.3535533905932737*m1SrSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherZ and uCrossOtherZ ... // 
  data->AEM_S(12,24) = 0.3535533905932737*m1SrOther[8]*mnuOther; 
  data->AEM_S(12,25) = 0.3535533905932737*m1SrOther[9]*mnuOther; 
  data->AEM_S(12,26) = 0.3535533905932737*m1SrOther[10]*mnuOther; 
  data->AEM_S(12,27) = 0.3535533905932737*m1SrOther[11]*mnuOther; 
  data->AEM_S(13,24) = 0.3535533905932737*m1SrOther[9]*mnuOther; 
  data->AEM_S(13,25) = 0.3535533905932737*m1SrOther[8]*mnuOther; 
  data->AEM_S(14,24) = 0.3535533905932737*m1SrOther[10]*mnuOther; 
  data->AEM_S(14,26) = 0.3535533905932737*m1SrOther[8]*mnuOther; 
  data->AEM_S(15,24) = 0.3535533905932737*m1SrOther[11]*mnuOther; 
  data->AEM_S(15,27) = 0.3535533905932737*m1SrOther[8]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
  mnuM1sum[9] += m1rSelf[9]*mnuSelf+m1rOther[9]*mnuOther; 
  mnuM1sum[10] += m1rSelf[10]*mnuSelf+m1rOther[10]*mnuOther; 
  mnuM1sum[11] += m1rSelf[11]*mnuSelf+m1rOther[11]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(12,12) = 0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 0.3535533905932737*m0SrSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(12,14) = 0.3535533905932737*m0SrSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(12,15) = 0.3535533905932737*m0SrSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(13,12) = 0.3535533905932737*m0SrSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(14,12) = 0.3535533905932737*m0SrSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(14,14) = 0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(15,12) = 0.3535533905932737*m0SrSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(15,15) = 0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(12,28) = 0.3535533905932737*m0SrOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(12,29) = 0.3535533905932737*m0SrOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(12,30) = 0.3535533905932737*m0SrOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(12,31) = 0.3535533905932737*m0SrOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(13,28) = 0.3535533905932737*m0SrOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(13,29) = 0.3535533905932737*m0SrOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(14,28) = 0.3535533905932737*m0SrOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(14,30) = 0.3535533905932737*m0SrOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(15,28) = 0.3535533905932737*m0SrOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(15,31) = 0.3535533905932737*m0SrOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
 
  double mnuM2sum[4]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2SrSelf[0]*mnuSelf+m2SrOther[0]*mnuOther; 
  mnuM2sum[1] = m2SrSelf[1]*mnuSelf+m2SrOther[1]*mnuOther; 
  mnuM2sum[2] = m2SrSelf[2]*mnuSelf+m2SrOther[2]*mnuOther; 
  mnuM2sum[3] = m2SrSelf[3]*mnuSelf+m2SrOther[3]*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<4,8>(0,4).setZero(); 
  data->AEM_S.block<8,4>(4,0).setZero(); 
  data->AEM_S.block<4,4>(4,8).setZero(); 
  data->AEM_S.block<4,4>(8,4).setZero(); 
  data->AEM_S.block<4,8>(0,20).setZero(); 
  data->AEM_S.block<8,4>(4,16).setZero(); 
  data->AEM_S.block<4,4>(4,24).setZero(); 
  data->AEM_S.block<4,4>(8,20).setZero(); 
 
  double m1Relax[12]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<12; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  double m1EffD[12]; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(16,0) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(16,1) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,2) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,3) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,0) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,1) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,0) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,2) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(19,0) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,3) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(16,12) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(16,13) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(16,14) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(16,15) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(17,12) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(17,13) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(18,12) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(18,14) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(19,12) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(19,15) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(16,16) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(16,17) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(16,18) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(16,19) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(17,16) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(17,17) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(18,16) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(18,18) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(19,16) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(19,19) = -0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(16,28) = 0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(16,29) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(16,30) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(16,31) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(17,28) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(17,29) = 0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(18,28) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(18,30) = 0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(19,28) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(19,31) = 0.3535533905932737*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(28,0) = (-0.125*m0rSelf[3]*uSelf[3]*mnuSelf)-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(28,1) = (-0.125*m0rSelf[0]*uSelf[1]*mnuSelf)+0.3535533905932737*m1SrSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(28,2) = (-0.125*m0rSelf[0]*uSelf[2]*mnuSelf)+0.3535533905932737*m1SrSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(28,3) = (-0.125*m0rSelf[0]*uSelf[3]*mnuSelf)+0.3535533905932737*m1SrSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,0) = (-0.125*m0rSelf[0]*uSelf[1]*mnuSelf)+0.3535533905932737*m1SrSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(29,1) = (-0.125*m0rSelf[3]*uSelf[3]*mnuSelf)-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.225*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(29,2) = (-0.125*m0rSelf[1]*uSelf[2]*mnuSelf)-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,3) = (-0.125*m0rSelf[1]*uSelf[3]*mnuSelf)-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(30,0) = (-0.125*m0rSelf[0]*uSelf[2]*mnuSelf)+0.3535533905932737*m1SrSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,1) = (-0.125*m0rSelf[1]*uSelf[2]*mnuSelf)-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,2) = (-0.125*m0rSelf[3]*uSelf[3]*mnuSelf)-0.225*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(30,3) = (-0.125*m0rSelf[2]*uSelf[3]*mnuSelf)-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,0) = (-0.125*m0rSelf[0]*uSelf[3]*mnuSelf)+0.3535533905932737*m1SrSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,1) = (-0.125*m0rSelf[1]*uSelf[3]*mnuSelf)-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,2) = (-0.125*m0rSelf[2]*uSelf[3]*mnuSelf)-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,3) = (-0.225*m0rSelf[3]*uSelf[3]*mnuSelf)-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(28,16) = 0.125*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1SrOther[0]*mnuOther; 
  data->AEM_S(28,17) = 0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1SrOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(28,18) = 0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1SrOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(28,19) = 0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1SrOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(29,16) = 0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1SrOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(29,17) = 0.125*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.225*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1SrOther[0]*mnuOther; 
  data->AEM_S(29,18) = 0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(29,19) = 0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(30,16) = 0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1SrOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,17) = 0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,18) = 0.125*m0rOther[3]*uOther[3]*mnuOther+0.225*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1SrOther[0]*mnuOther; 
  data->AEM_S(30,19) = 0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(31,16) = 0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1SrOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(31,17) = 0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(31,18) = 0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(31,19) = 0.225*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1SrOther[0]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfX-m0Self*m1OtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[0] = 0.3535533905932737*m0rOther[3]*m1rSelf[3]-0.3535533905932737*m0rSelf[3]*m1rOther[3]+0.3535533905932737*m0rOther[2]*m1rSelf[2]-0.3535533905932737*m0rSelf[2]*m1rOther[2]+0.3535533905932737*m0rOther[1]*m1rSelf[1]-0.3535533905932737*m0rSelf[1]*m1rOther[1]+0.3535533905932737*m0rOther[0]*m1rSelf[0]-0.3535533905932737*m0rSelf[0]*m1rOther[0]; 
  m1EffD[1] = 0.3535533905932737*m0rOther[0]*m1rSelf[1]-0.3535533905932737*m0rSelf[0]*m1rOther[1]-0.3535533905932737*m1rOther[0]*m0rSelf[1]+0.3535533905932737*m1rSelf[0]*m0rOther[1]; 
  m1EffD[2] = 0.3535533905932737*m0rOther[0]*m1rSelf[2]-0.3535533905932737*m0rSelf[0]*m1rOther[2]-0.3535533905932737*m1rOther[0]*m0rSelf[2]+0.3535533905932737*m1rSelf[0]*m0rOther[2]; 
  m1EffD[3] = 0.3535533905932737*m0rOther[0]*m1rSelf[3]-0.3535533905932737*m0rSelf[0]*m1rOther[3]-0.3535533905932737*m1rOther[0]*m0rSelf[3]+0.3535533905932737*m1rSelf[0]*m0rOther[3]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(4,4); 
  dataDiv->AEM_S(0,0) = 0.3535533905932737*m0rSelf[0]*mnuSelf+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(0,2) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(0,3) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.3535533905932737*m0rSelf[0]*mnuSelf+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(2,0) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,2) = 0.3535533905932737*m0rSelf[0]*mnuSelf+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(3,0) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,3) = 0.3535533905932737*m0rSelf[0]*mnuSelf+0.3535533905932737*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[0],m1EffD[1],m1EffD[2],m1EffD[3]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+0,4,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (-2.0*m1EffD[0]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (-2.0*m1EffD[1]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (-2.0*m1EffD[2]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += (-2.0*m1EffD[3]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(20,4) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,5) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,6) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,7) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,4) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,5) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(22,4) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,6) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(23,4) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,7) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(20,12) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(20,13) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(20,14) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(20,15) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(21,12) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(21,13) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(22,12) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(22,14) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(23,12) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(23,15) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(20,20) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(20,21) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(20,22) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(20,23) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(21,20) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(21,21) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(22,20) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(22,22) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(23,20) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(23,23) = -0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(20,28) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(20,29) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(20,30) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(20,31) = 0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(21,28) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(21,29) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(22,28) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(22,30) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(23,28) = 0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(23,31) = 0.3535533905932737*cMOther[4]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(28,4) = (-0.125*m0rSelf[3]*uSelf[7]*mnuSelf)-0.125*m0rSelf[2]*uSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(28,5) = (-0.125*m0rSelf[0]*uSelf[5]*mnuSelf)+0.3535533905932737*m1SrSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[4]*mnuSelf; 
  data->AEM_S(28,6) = (-0.125*m0rSelf[0]*uSelf[6]*mnuSelf)+0.3535533905932737*m1SrSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(28,7) = (-0.125*m0rSelf[0]*uSelf[7]*mnuSelf)+0.3535533905932737*m1SrSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf; 
  data->AEM_S(29,4) = (-0.125*m0rSelf[0]*uSelf[5]*mnuSelf)+0.3535533905932737*m1SrSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[4]*mnuSelf; 
  data->AEM_S(29,5) = (-0.125*m0rSelf[3]*uSelf[7]*mnuSelf)-0.125*m0rSelf[2]*uSelf[6]*mnuSelf-0.225*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(29,6) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*m0rSelf[2]*uSelf[5]*mnuSelf; 
  data->AEM_S(29,7) = (-0.125*m0rSelf[1]*uSelf[7]*mnuSelf)-0.125*m0rSelf[3]*uSelf[5]*mnuSelf; 
  data->AEM_S(30,4) = (-0.125*m0rSelf[0]*uSelf[6]*mnuSelf)+0.3535533905932737*m1SrSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(30,5) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*m0rSelf[2]*uSelf[5]*mnuSelf; 
  data->AEM_S(30,6) = (-0.125*m0rSelf[3]*uSelf[7]*mnuSelf)-0.225*m0rSelf[2]*uSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(30,7) = (-0.125*m0rSelf[2]*uSelf[7]*mnuSelf)-0.125*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(31,4) = (-0.125*m0rSelf[0]*uSelf[7]*mnuSelf)+0.3535533905932737*m1SrSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf; 
  data->AEM_S(31,5) = (-0.125*m0rSelf[1]*uSelf[7]*mnuSelf)-0.125*m0rSelf[3]*uSelf[5]*mnuSelf; 
  data->AEM_S(31,6) = (-0.125*m0rSelf[2]*uSelf[7]*mnuSelf)-0.125*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(31,7) = (-0.225*m0rSelf[3]*uSelf[7]*mnuSelf)-0.125*m0rSelf[2]*uSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1SrSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(28,20) = 0.125*m0rOther[3]*uOther[7]*mnuOther+0.125*m0rOther[2]*uOther[6]*mnuOther+0.125*m0rOther[1]*uOther[5]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1SrOther[4]*mnuOther; 
  data->AEM_S(28,21) = 0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1SrOther[5]*mnuOther+0.125*m0rOther[1]*uOther[4]*mnuOther; 
  data->AEM_S(28,22) = 0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1SrOther[6]*mnuOther+0.125*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(28,23) = 0.125*m0rOther[0]*uOther[7]*mnuOther-0.3535533905932737*m1SrOther[7]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther; 
  data->AEM_S(29,20) = 0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1SrOther[5]*mnuOther+0.125*m0rOther[1]*uOther[4]*mnuOther; 
  data->AEM_S(29,21) = 0.125*m0rOther[3]*uOther[7]*mnuOther+0.125*m0rOther[2]*uOther[6]*mnuOther+0.225*m0rOther[1]*uOther[5]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1SrOther[4]*mnuOther; 
  data->AEM_S(29,22) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther; 
  data->AEM_S(29,23) = 0.125*m0rOther[1]*uOther[7]*mnuOther+0.125*m0rOther[3]*uOther[5]*mnuOther; 
  data->AEM_S(30,20) = 0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1SrOther[6]*mnuOther+0.125*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(30,21) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther; 
  data->AEM_S(30,22) = 0.125*m0rOther[3]*uOther[7]*mnuOther+0.225*m0rOther[2]*uOther[6]*mnuOther+0.125*m0rOther[1]*uOther[5]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1SrOther[4]*mnuOther; 
  data->AEM_S(30,23) = 0.125*m0rOther[2]*uOther[7]*mnuOther+0.125*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(31,20) = 0.125*m0rOther[0]*uOther[7]*mnuOther-0.3535533905932737*m1SrOther[7]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther; 
  data->AEM_S(31,21) = 0.125*m0rOther[1]*uOther[7]*mnuOther+0.125*m0rOther[3]*uOther[5]*mnuOther; 
  data->AEM_S(31,22) = 0.125*m0rOther[2]*uOther[7]*mnuOther+0.125*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(31,23) = 0.225*m0rOther[3]*uOther[7]*mnuOther+0.125*m0rOther[2]*uOther[6]*mnuOther+0.125*m0rOther[1]*uOther[5]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1SrOther[4]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfY-m0Self*m1OtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[4] = 0.3535533905932737*m0rOther[3]*m1rSelf[7]-0.3535533905932737*m0rSelf[3]*m1rOther[7]+0.3535533905932737*m0rOther[2]*m1rSelf[6]-0.3535533905932737*m0rSelf[2]*m1rOther[6]+0.3535533905932737*m0rOther[1]*m1rSelf[5]-0.3535533905932737*m0rSelf[1]*m1rOther[5]+0.3535533905932737*m0rOther[0]*m1rSelf[4]-0.3535533905932737*m0rSelf[0]*m1rOther[4]; 
  m1EffD[5] = 0.3535533905932737*m0rOther[0]*m1rSelf[5]-0.3535533905932737*m0rSelf[0]*m1rOther[5]+0.3535533905932737*m0rOther[1]*m1rSelf[4]-0.3535533905932737*m0rSelf[1]*m1rOther[4]; 
  m1EffD[6] = 0.3535533905932737*m0rOther[0]*m1rSelf[6]-0.3535533905932737*m0rSelf[0]*m1rOther[6]+0.3535533905932737*m0rOther[2]*m1rSelf[4]-0.3535533905932737*m0rSelf[2]*m1rOther[4]; 
  m1EffD[7] = 0.3535533905932737*m0rOther[0]*m1rSelf[7]-0.3535533905932737*m0rSelf[0]*m1rOther[7]+0.3535533905932737*m0rOther[3]*m1rSelf[4]-0.3535533905932737*m0rSelf[3]*m1rOther[4]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[4],m1EffD[5],m1EffD[6],m1EffD[7]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+4,4,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[4] += (-2.0*m1EffD[4]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (-2.0*m1EffD[5]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
  m1Relax[6] += (-2.0*m1EffD[6]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (-2.0*m1EffD[7]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(24,8) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(24,9) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(24,10) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,11) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,8) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,9) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(26,8) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,10) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(27,8) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,11) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(24,12) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(24,13) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(24,14) = -0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(24,15) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(25,12) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(25,13) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(26,12) = -0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(26,14) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(27,12) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(27,15) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(24,24) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(24,25) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(24,26) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(24,27) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(25,24) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(25,25) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(26,24) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(26,26) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(27,24) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(27,27) = -0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(24,28) = 0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(24,29) = 0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(24,30) = 0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(24,31) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(25,28) = 0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(25,29) = 0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(26,28) = 0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(26,30) = 0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(27,28) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(27,31) = 0.3535533905932737*cMOther[8]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfZ-uSelfZ*m0Self) and uCrossSelfZ ... // 
  data->AEM_S(28,8) = (-0.125*m0rSelf[3]*uSelf[11]*mnuSelf)-0.125*m0rSelf[2]*uSelf[10]*mnuSelf-0.125*m0rSelf[1]*uSelf[9]*mnuSelf-0.125*m0rSelf[0]*uSelf[8]*mnuSelf+0.3535533905932737*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(28,9) = (-0.125*m0rSelf[0]*uSelf[9]*mnuSelf)+0.3535533905932737*m1SrSelf[9]*mnuSelf-0.125*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(28,10) = (-0.125*m0rSelf[0]*uSelf[10]*mnuSelf)+0.3535533905932737*m1SrSelf[10]*mnuSelf-0.125*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(28,11) = (-0.125*m0rSelf[0]*uSelf[11]*mnuSelf)+0.3535533905932737*m1SrSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(29,8) = (-0.125*m0rSelf[0]*uSelf[9]*mnuSelf)+0.3535533905932737*m1SrSelf[9]*mnuSelf-0.125*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(29,9) = (-0.125*m0rSelf[3]*uSelf[11]*mnuSelf)-0.125*m0rSelf[2]*uSelf[10]*mnuSelf-0.225*m0rSelf[1]*uSelf[9]*mnuSelf-0.125*m0rSelf[0]*uSelf[8]*mnuSelf+0.3535533905932737*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(29,10) = (-0.125*m0rSelf[1]*uSelf[10]*mnuSelf)-0.125*m0rSelf[2]*uSelf[9]*mnuSelf; 
  data->AEM_S(29,11) = (-0.125*m0rSelf[1]*uSelf[11]*mnuSelf)-0.125*m0rSelf[3]*uSelf[9]*mnuSelf; 
  data->AEM_S(30,8) = (-0.125*m0rSelf[0]*uSelf[10]*mnuSelf)+0.3535533905932737*m1SrSelf[10]*mnuSelf-0.125*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(30,9) = (-0.125*m0rSelf[1]*uSelf[10]*mnuSelf)-0.125*m0rSelf[2]*uSelf[9]*mnuSelf; 
  data->AEM_S(30,10) = (-0.125*m0rSelf[3]*uSelf[11]*mnuSelf)-0.225*m0rSelf[2]*uSelf[10]*mnuSelf-0.125*m0rSelf[1]*uSelf[9]*mnuSelf-0.125*m0rSelf[0]*uSelf[8]*mnuSelf+0.3535533905932737*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(30,11) = (-0.125*m0rSelf[2]*uSelf[11]*mnuSelf)-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(31,8) = (-0.125*m0rSelf[0]*uSelf[11]*mnuSelf)+0.3535533905932737*m1SrSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(31,9) = (-0.125*m0rSelf[1]*uSelf[11]*mnuSelf)-0.125*m0rSelf[3]*uSelf[9]*mnuSelf; 
  data->AEM_S(31,10) = (-0.125*m0rSelf[2]*uSelf[11]*mnuSelf)-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(31,11) = (-0.225*m0rSelf[3]*uSelf[11]*mnuSelf)-0.125*m0rSelf[2]*uSelf[10]*mnuSelf-0.125*m0rSelf[1]*uSelf[9]*mnuSelf-0.125*m0rSelf[0]*uSelf[8]*mnuSelf+0.3535533905932737*m1SrSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherZ-uOtherZ*m0Other) and uCrossOtherZ ... // 
  data->AEM_S(28,24) = 0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther+0.125*m0rOther[1]*uOther[9]*mnuOther+0.125*m0rOther[0]*uOther[8]*mnuOther-0.3535533905932737*m1SrOther[8]*mnuOther; 
  data->AEM_S(28,25) = 0.125*m0rOther[0]*uOther[9]*mnuOther-0.3535533905932737*m1SrOther[9]*mnuOther+0.125*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(28,26) = 0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1SrOther[10]*mnuOther+0.125*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(28,27) = 0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1SrOther[11]*mnuOther+0.125*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(29,24) = 0.125*m0rOther[0]*uOther[9]*mnuOther-0.3535533905932737*m1SrOther[9]*mnuOther+0.125*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(29,25) = 0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther+0.225*m0rOther[1]*uOther[9]*mnuOther+0.125*m0rOther[0]*uOther[8]*mnuOther-0.3535533905932737*m1SrOther[8]*mnuOther; 
  data->AEM_S(29,26) = 0.125*m0rOther[1]*uOther[10]*mnuOther+0.125*m0rOther[2]*uOther[9]*mnuOther; 
  data->AEM_S(29,27) = 0.125*m0rOther[1]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[9]*mnuOther; 
  data->AEM_S(30,24) = 0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1SrOther[10]*mnuOther+0.125*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(30,25) = 0.125*m0rOther[1]*uOther[10]*mnuOther+0.125*m0rOther[2]*uOther[9]*mnuOther; 
  data->AEM_S(30,26) = 0.125*m0rOther[3]*uOther[11]*mnuOther+0.225*m0rOther[2]*uOther[10]*mnuOther+0.125*m0rOther[1]*uOther[9]*mnuOther+0.125*m0rOther[0]*uOther[8]*mnuOther-0.3535533905932737*m1SrOther[8]*mnuOther; 
  data->AEM_S(30,27) = 0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(31,24) = 0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1SrOther[11]*mnuOther+0.125*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(31,25) = 0.125*m0rOther[1]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[9]*mnuOther; 
  data->AEM_S(31,26) = 0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(31,27) = 0.225*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther+0.125*m0rOther[1]*uOther[9]*mnuOther+0.125*m0rOther[0]*uOther[8]*mnuOther-0.3535533905932737*m1SrOther[8]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfZ-m0Self*m1OtherZ)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[8] = 0.3535533905932737*m0rOther[3]*m1rSelf[11]-0.3535533905932737*m0rSelf[3]*m1rOther[11]+0.3535533905932737*m0rOther[2]*m1rSelf[10]-0.3535533905932737*m0rSelf[2]*m1rOther[10]+0.3535533905932737*m0rOther[1]*m1rSelf[9]-0.3535533905932737*m0rSelf[1]*m1rOther[9]+0.3535533905932737*m0rOther[0]*m1rSelf[8]-0.3535533905932737*m0rSelf[0]*m1rOther[8]; 
  m1EffD[9] = 0.3535533905932737*m0rOther[0]*m1rSelf[9]-0.3535533905932737*m0rSelf[0]*m1rOther[9]+0.3535533905932737*m0rOther[1]*m1rSelf[8]-0.3535533905932737*m0rSelf[1]*m1rOther[8]; 
  m1EffD[10] = 0.3535533905932737*m0rOther[0]*m1rSelf[10]-0.3535533905932737*m0rSelf[0]*m1rOther[10]+0.3535533905932737*m0rOther[2]*m1rSelf[8]-0.3535533905932737*m0rSelf[2]*m1rOther[8]; 
  m1EffD[11] = 0.3535533905932737*m0rOther[0]*m1rSelf[11]-0.3535533905932737*m0rSelf[0]*m1rOther[11]+0.3535533905932737*m0rOther[3]*m1rSelf[8]-0.3535533905932737*m0rSelf[3]*m1rOther[8]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[8],m1EffD[9],m1EffD[10],m1EffD[11]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+8,4,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 3 of momentum relaxation. 
  m1Relax[8] += (-2.0*m1EffD[8]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
  m1Relax[9] += (-2.0*m1EffD[9]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[9]*mnuSelf-1.0*m1rOther[9]*mnuOther; 
  m1Relax[10] += (-2.0*m1EffD[10]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[10]*mnuSelf-1.0*m1rOther[10]*mnuOther; 
  m1Relax[11] += (-2.0*m1EffD[11]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[11]*mnuSelf-1.0*m1rOther[11]*mnuOther; 
 
  double ucMSelf[4]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.3535533905932737*cMSelf[a0+3]*uSelf[a0+3]+0.3535533905932737*cMSelf[a0+2]*uSelf[a0+2]+0.3535533905932737*cMSelf[a0+1]*uSelf[a0+1]+0.3535533905932737*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.3535533905932737*cMSelf[a0]*uSelf[a0+1]+0.3535533905932737*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.3535533905932737*cMSelf[a0]*uSelf[a0+2]+0.3535533905932737*uSelf[a0]*cMSelf[a0+2]; 
    ucMSelf[3] += 0.3535533905932737*cMSelf[a0]*uSelf[a0+3]+0.3535533905932737*uSelf[a0]*cMSelf[a0+3]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(28,12) = 0.3535533905932737*ucMSelf[0]*mnuSelf+0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(28,13) = 0.3535533905932737*ucMSelf[1]*mnuSelf+0.3535533905932737*m0SrSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(28,14) = 0.3535533905932737*ucMSelf[2]*mnuSelf+0.3535533905932737*m0SrSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(28,15) = 0.3535533905932737*ucMSelf[3]*mnuSelf+0.3535533905932737*m0SrSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(29,12) = 0.3535533905932737*ucMSelf[1]*mnuSelf+0.3535533905932737*m0SrSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(29,13) = 0.3535533905932737*ucMSelf[0]*mnuSelf+0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(30,12) = 0.3535533905932737*ucMSelf[2]*mnuSelf+0.3535533905932737*m0SrSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(30,14) = 0.3535533905932737*ucMSelf[0]*mnuSelf+0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(31,12) = 0.3535533905932737*ucMSelf[3]*mnuSelf+0.3535533905932737*m0SrSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(31,15) = 0.3535533905932737*ucMSelf[0]*mnuSelf+0.3535533905932737*m0SrSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
 
  double ucMOther[4]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.3535533905932737*cMOther[a0+3]*uOther[a0+3]+0.3535533905932737*cMOther[a0+2]*uOther[a0+2]+0.3535533905932737*cMOther[a0+1]*uOther[a0+1]+0.3535533905932737*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.3535533905932737*cMOther[a0]*uOther[a0+1]+0.3535533905932737*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.3535533905932737*cMOther[a0]*uOther[a0+2]+0.3535533905932737*uOther[a0]*cMOther[a0+2]; 
    ucMOther[3] += 0.3535533905932737*cMOther[a0]*uOther[a0+3]+0.3535533905932737*uOther[a0]*cMOther[a0+3]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(28,28) = (-0.3535533905932737*ucMOther[0]*mnuOther)-0.3535533905932737*m0SrOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(28,29) = (-0.3535533905932737*ucMOther[1]*mnuOther)-0.3535533905932737*m0SrOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(28,30) = (-0.3535533905932737*ucMOther[2]*mnuOther)-0.3535533905932737*m0SrOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(28,31) = (-0.3535533905932737*ucMOther[3]*mnuOther)-0.3535533905932737*m0SrOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(29,28) = (-0.3535533905932737*ucMOther[1]*mnuOther)-0.3535533905932737*m0SrOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(29,29) = (-0.3535533905932737*ucMOther[0]*mnuOther)-0.3535533905932737*m0SrOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(30,28) = (-0.3535533905932737*ucMOther[2]*mnuOther)-0.3535533905932737*m0SrOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(30,30) = (-0.3535533905932737*ucMOther[0]*mnuOther)-0.3535533905932737*m0SrOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(31,28) = (-0.3535533905932737*ucMOther[3]*mnuOther)-0.3535533905932737*m0SrOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(31,31) = (-0.3535533905932737*ucMOther[0]*mnuOther)-0.3535533905932737*m0SrOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
 
  double kinESelf[4]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinESelf[0] += 0.3535533905932737*m1rSelf[a0+3]*uSelf[a0+3]+0.3535533905932737*m1rSelf[a0+2]*uSelf[a0+2]+0.3535533905932737*m1rSelf[a0+1]*uSelf[a0+1]+0.3535533905932737*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.3535533905932737*m1rSelf[a0]*uSelf[a0+1]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.3535533905932737*m1rSelf[a0]*uSelf[a0+2]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+2]; 
    kinESelf[3] += 0.3535533905932737*m1rSelf[a0]*uSelf[a0+3]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+3]; 
  } 
 
  double kinEOther[4]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    kinEOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.3535533905932737*m1rOther[a0+3]*uOther[a0+3]+0.3535533905932737*m1rOther[a0+2]*uOther[a0+2]+0.3535533905932737*m1rOther[a0+1]*uOther[a0+1]+0.3535533905932737*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.3535533905932737*m1rOther[a0]*uOther[a0+1]+0.3535533905932737*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.3535533905932737*m1rOther[a0]*uOther[a0+2]+0.3535533905932737*uOther[a0]*m1rOther[a0+2]; 
    kinEOther[3] += 0.3535533905932737*m1rOther[a0]*uOther[a0+3]+0.3535533905932737*uOther[a0]*m1rOther[a0+3]; 
  } 
 
  double relKinE[4]; 
  // zero out array with dot product of uSelf-uOther and m1EffD. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.3535533905932737*m1EffD[a0+3]*uSelf[a0+3]-0.3535533905932737*m1EffD[a0+3]*uOther[a0+3]+0.3535533905932737*m1EffD[a0+2]*uSelf[a0+2]-0.3535533905932737*m1EffD[a0+2]*uOther[a0+2]+0.3535533905932737*m1EffD[a0+1]*uSelf[a0+1]-0.3535533905932737*m1EffD[a0+1]*uOther[a0+1]+0.3535533905932737*m1EffD[a0]*uSelf[a0]-0.3535533905932737*m1EffD[a0]*uOther[a0]; 
    relKinE[1] += 0.3535533905932737*m1EffD[a0]*uSelf[a0+1]-0.3535533905932737*m1EffD[a0]*uOther[a0+1]+0.3535533905932737*uSelf[a0]*m1EffD[a0+1]-0.3535533905932737*uOther[a0]*m1EffD[a0+1]; 
    relKinE[2] += 0.3535533905932737*m1EffD[a0]*uSelf[a0+2]-0.3535533905932737*m1EffD[a0]*uOther[a0+2]+0.3535533905932737*uSelf[a0]*m1EffD[a0+2]-0.3535533905932737*uOther[a0]*m1EffD[a0+2]; 
    relKinE[3] += 0.3535533905932737*m1EffD[a0]*uSelf[a0+3]-0.3535533905932737*m1EffD[a0]*uOther[a0+3]+0.3535533905932737*uSelf[a0]*m1EffD[a0+3]-0.3535533905932737*uOther[a0]*m1EffD[a0+3]; 
  } 
 
  // Divide m0Other*(m2Self-kinESelf) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Other and m2Self-uSelf.m1Self. 
  double m0OtherThESelf[4]; 
  m0OtherThESelf[0] = 0.3535533905932737*m0rOther[3]*m2rSelf[3]-0.3535533905932737*kinESelf[3]*m0rOther[3]+0.3535533905932737*m0rOther[2]*m2rSelf[2]-0.3535533905932737*kinESelf[2]*m0rOther[2]+0.3535533905932737*m0rOther[1]*m2rSelf[1]-0.3535533905932737*kinESelf[1]*m0rOther[1]+0.3535533905932737*m0rOther[0]*m2rSelf[0]-0.3535533905932737*kinESelf[0]*m0rOther[0]; 
  m0OtherThESelf[1] = 0.3535533905932737*m0rOther[0]*m2rSelf[1]+0.3535533905932737*m2rSelf[0]*m0rOther[1]-0.3535533905932737*kinESelf[0]*m0rOther[1]-0.3535533905932737*m0rOther[0]*kinESelf[1]; 
  m0OtherThESelf[2] = 0.3535533905932737*m0rOther[0]*m2rSelf[2]+0.3535533905932737*m2rSelf[0]*m0rOther[2]-0.3535533905932737*kinESelf[0]*m0rOther[2]-0.3535533905932737*m0rOther[0]*kinESelf[2]; 
  m0OtherThESelf[3] = 0.3535533905932737*m0rOther[0]*m2rSelf[3]+0.3535533905932737*m2rSelf[0]*m0rOther[3]-0.3535533905932737*kinESelf[0]*m0rOther[3]-0.3535533905932737*m0rOther[0]*kinESelf[3]; 
  dataDiv->BEV_S << m0OtherThESelf[0],m0OtherThESelf[1],m0OtherThESelf[2],m0OtherThESelf[3]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthSelf[4]; 
  Eigen::Map<VectorXd>(effEthSelf,4,1) = dataDiv->u_S; 
 
  // Divide m0Self*(m2Other-kinEOther) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Self and m2Other-uOther.m1Other. 
  double m0SelfThEOther[4]; 
  m0SelfThEOther[0] = 0.3535533905932737*m0rSelf[3]*m2rOther[3]-0.3535533905932737*kinEOther[3]*m0rSelf[3]+0.3535533905932737*m0rSelf[2]*m2rOther[2]-0.3535533905932737*kinEOther[2]*m0rSelf[2]+0.3535533905932737*m0rSelf[1]*m2rOther[1]-0.3535533905932737*kinEOther[1]*m0rSelf[1]+0.3535533905932737*m0rSelf[0]*m2rOther[0]-0.3535533905932737*kinEOther[0]*m0rSelf[0]; 
  m0SelfThEOther[1] = 0.3535533905932737*m0rSelf[0]*m2rOther[1]+0.3535533905932737*m2rOther[0]*m0rSelf[1]-0.3535533905932737*kinEOther[0]*m0rSelf[1]-0.3535533905932737*m0rSelf[0]*kinEOther[1]; 
  m0SelfThEOther[2] = 0.3535533905932737*m0rSelf[0]*m2rOther[2]+0.3535533905932737*m2rOther[0]*m0rSelf[2]-0.3535533905932737*kinEOther[0]*m0rSelf[2]-0.3535533905932737*m0rSelf[0]*kinEOther[2]; 
  m0SelfThEOther[3] = 0.3535533905932737*m0rSelf[0]*m2rOther[3]+0.3535533905932737*m2rOther[0]*m0rSelf[3]-0.3535533905932737*kinEOther[0]*m0rSelf[3]-0.3535533905932737*m0rSelf[0]*kinEOther[3]; 
  dataDiv->BEV_S << m0SelfThEOther[0],m0SelfThEOther[1],m0SelfThEOther[2],m0SelfThEOther[3]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthOther[4]; 
  Eigen::Map<VectorXd>(effEthOther,4,1) = dataDiv->u_S; 
 
  double m2Relax[4]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(1.0*relKinE[0]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[0]*mSelf)/(mSelf+mOther)+(1.0*relKinE[0]*mOther)/(mSelf+mOther)+(2.0*effEthOther[0]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2SrSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2SrOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(1.0*relKinE[1]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[1]*mSelf)/(mSelf+mOther)+(1.0*relKinE[1]*mOther)/(mSelf+mOther)+(2.0*effEthOther[1]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2SrSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2SrOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(1.0*relKinE[2]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[2]*mSelf)/(mSelf+mOther)+(1.0*relKinE[2]*mOther)/(mSelf+mOther)+(2.0*effEthOther[2]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2SrSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2SrOther[2])*mnuOther; 
  m2Relax[3] = betaGreenep1*((-(1.0*relKinE[3]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[3]*mSelf)/(mSelf+mOther)+(1.0*relKinE[3]*mOther)/(mSelf+mOther)+(2.0*effEthOther[3]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2SrSelf[3]-1.0*kinESelf[3])*mnuSelf+(kinEOther[3]-1.0*m2SrOther[3])*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<4,8>(16,4).setZero(); 
  data->AEM_S.block<8,4>(20,0).setZero(); 
  data->AEM_S.block<4,4>(20,8).setZero(); 
  data->AEM_S.block<4,4>(24,4).setZero(); 
  data->AEM_S.block<4,8>(16,20).setZero(); 
  data->AEM_S.block<8,4>(20,16).setZero(); 
  data->AEM_S.block<4,4>(20,24).setZero(); 
  data->AEM_S.block<4,4>(24,20).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM1sum[10],mnuM1sum[11],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m1Relax[9],m1Relax[10],m1Relax[11],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,12,1) = data->u_S.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,4,1) = data->u_S.segment<4>(12); 
 
  Eigen::Map<VectorXd>(uCrossOther,12,1) = data->u_S.segment<12>(16); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,4,1) = data->u_S.segment<4>(28); 
 
} 
 
void VmLBOCrossPrimMoments3x3vMax_P2(binOpData_t *data, binOpData_t *dataDiv,const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]+1.060660171779821*m0Self[5]+1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]+1.060660171779821*m0Self[5]+1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]+1.060660171779821*m0Self[5]+1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]+1.060660171779821*m0Self[5]+1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]-0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]-1.060660171779821*m0Self[5]-1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]-1.060660171779821*m0Self[5]-1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]-1.060660171779821*m0Self[5]-1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Self[9]+0.7905694150420947*m0Self[8]+0.7905694150420947*m0Self[7]+1.060660171779821*m0Self[6]-1.060660171779821*m0Self[5]-1.060660171779821*m0Self[4]-0.6123724356957944*m0Self[3]-0.6123724356957944*m0Self[2]+0.6123724356957944*m0Self[1]+0.3535533905932737*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[10]; 
  double m1rSelf[30]; 
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
    m1rSelf[20] = m1Self[20]; 
    m1rSelf[21] = 0.0; 
    m1rSelf[22] = 0.0; 
    m1rSelf[23] = 0.0; 
    m1rSelf[24] = 0.0; 
    m1rSelf[25] = 0.0; 
    m1rSelf[26] = 0.0; 
    m1rSelf[27] = 0.0; 
    m1rSelf[28] = 0.0; 
    m1rSelf[29] = 0.0; 
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
    m1rSelf[20] = m1Self[20]; 
    m1rSelf[21] = m1Self[21]; 
    m1rSelf[22] = m1Self[22]; 
    m1rSelf[23] = m1Self[23]; 
    m1rSelf[24] = m1Self[24]; 
    m1rSelf[25] = m1Self[25]; 
    m1rSelf[26] = m1Self[26]; 
    m1rSelf[27] = m1Self[27]; 
    m1rSelf[28] = m1Self[28]; 
    m1rSelf[29] = m1Self[29]; 
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
 
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]+1.060660171779821*m0Other[5]+1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]+1.060660171779821*m0Other[5]+1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]+1.060660171779821*m0Other[5]+1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]+1.060660171779821*m0Other[5]+1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]-0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]-1.060660171779821*m0Other[5]-1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]-1.060660171779821*m0Other[5]-1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]-1.060660171779821*m0Other[5]-1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.7905694150420947*m0Other[9]+0.7905694150420947*m0Other[8]+0.7905694150420947*m0Other[7]+1.060660171779821*m0Other[6]-1.060660171779821*m0Other[5]-1.060660171779821*m0Other[4]-0.6123724356957944*m0Other[3]-0.6123724356957944*m0Other[2]+0.6123724356957944*m0Other[1]+0.3535533905932737*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[10]; 
  double m1rOther[30]; 
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
    m1rOther[20] = m1Other[20]; 
    m1rOther[21] = 0.0; 
    m1rOther[22] = 0.0; 
    m1rOther[23] = 0.0; 
    m1rOther[24] = 0.0; 
    m1rOther[25] = 0.0; 
    m1rOther[26] = 0.0; 
    m1rOther[27] = 0.0; 
    m1rOther[28] = 0.0; 
    m1rOther[29] = 0.0; 
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
    m1rOther[20] = m1Other[20]; 
    m1rOther[21] = m1Other[21]; 
    m1rOther[22] = m1Other[22]; 
    m1rOther[23] = m1Other[23]; 
    m1rOther[24] = m1Other[24]; 
    m1rOther[25] = m1Other[25]; 
    m1rOther[26] = m1Other[26]; 
    m1rOther[27] = m1Other[27]; 
    m1rOther[28] = m1Other[28]; 
    m1rOther[29] = m1Other[29]; 
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
 
  // Declare Eigen matrix and vectors for weak system. 
  data->AEM_S = Eigen::MatrixXd::Zero(80,80); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double mnuM1sum[30]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<30; vd++) 
  { 
    mnuM1sum[vd] = 0.0; 
  } 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(0,0) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(0,2) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(0,3) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(0,4) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(0,5) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(0,6) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(0,7) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(0,8) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(0,9) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(1,0) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(1,2) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(1,3) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(1,4) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(1,5) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,7) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,0) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,1) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(2,2) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(2,4) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,6) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(2,8) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,0) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,1) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(3,2) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(3,3) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(3,5) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,6) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,9) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(4,0) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(4,1) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,2) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,4) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,5) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(4,6) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(4,7) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(4,8) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(5,0) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(5,1) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(5,3) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,4) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(5,5) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(5,6) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(5,7) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(5,9) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(6,0) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(6,2) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(6,3) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,4) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(6,5) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(6,6) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,8) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(6,9) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(7,0) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(7,1) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,4) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(7,5) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(7,7) = 0.2258769757263128*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,0) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(8,2) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,4) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(8,6) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(8,8) = 0.2258769757263128*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(9,0) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(9,3) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,5) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(9,6) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(9,9) = 0.2258769757263128*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,30) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,31) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,32) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,33) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,34) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,35) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(0,36) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(0,37) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(0,38) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(0,39) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(1,30) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,31) = (-0.3162277660168379*cMSelf[7]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,32) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(1,33) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(1,34) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,35) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,37) = -0.3162277660168379*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,30) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,31) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(2,32) = (-0.3162277660168379*cMSelf[8]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,33) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(2,34) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,36) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,38) = -0.3162277660168379*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,30) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,31) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(3,32) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(3,33) = (-0.3162277660168379*cMSelf[9]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,35) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,36) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,39) = -0.3162277660168379*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,30) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,31) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(4,32) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,34) = (-0.3162277660168379*cMSelf[8]*mnuSelf)-0.3162277660168379*cMSelf[7]*mnuSelf-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(4,35) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(4,36) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(4,37) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,38) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(5,30) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,31) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,33) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(5,34) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(5,35) = (-0.3162277660168379*cMSelf[9]*mnuSelf)-0.3162277660168379*cMSelf[7]*mnuSelf-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,36) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(5,37) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,39) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(6,30) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,32) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(6,33) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(6,34) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(6,35) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(6,36) = (-0.3162277660168379*cMSelf[9]*mnuSelf)-0.3162277660168379*cMSelf[8]*mnuSelf-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(6,38) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,39) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(7,30) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,31) = -0.3162277660168379*cMSelf[1]*mnuSelf; 
  data->AEM_S(7,34) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(7,35) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(7,37) = (-0.2258769757263128*cMSelf[7]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(8,30) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,32) = -0.3162277660168379*cMSelf[2]*mnuSelf; 
  data->AEM_S(8,34) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(8,36) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(8,38) = (-0.2258769757263128*cMSelf[8]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(9,30) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,33) = -0.3162277660168379*cMSelf[3]*mnuSelf; 
  data->AEM_S(9,35) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(9,36) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(9,39) = (-0.2258769757263128*cMSelf[9]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,40) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(0,41) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(0,42) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(0,43) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(0,44) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(0,45) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(0,46) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(0,47) = 0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(0,48) = 0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(0,49) = 0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(1,40) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(1,41) = 0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(1,42) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(1,43) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(1,44) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(1,45) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(1,47) = 0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(2,40) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(2,41) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(2,42) = 0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(2,43) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(2,44) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(2,46) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(2,48) = 0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(3,40) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(3,41) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(3,42) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(3,43) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(3,45) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(3,46) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(3,49) = 0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(4,40) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(4,41) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(4,42) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(4,44) = 0.3162277660168379*m0rOther[8]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(4,45) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(4,46) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(4,47) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(4,48) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(5,40) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(5,41) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(5,43) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(5,44) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(5,45) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(5,46) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(5,47) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(5,49) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(6,40) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(6,42) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(6,43) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(6,44) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(6,45) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(6,46) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(6,48) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(6,49) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(7,40) = 0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(7,41) = 0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(7,44) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(7,45) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(7,47) = 0.2258769757263128*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(8,40) = 0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(8,42) = 0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(8,44) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(8,46) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(8,48) = 0.2258769757263128*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(9,40) = 0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(9,43) = 0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(9,45) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(9,46) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(9,49) = 0.2258769757263128*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,70) = -0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(0,71) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(0,72) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(0,73) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(0,74) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(0,75) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(0,76) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(0,77) = -0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(0,78) = -0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(0,79) = -0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(1,70) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(1,71) = (-0.3162277660168379*cMOther[7]*mnuOther)-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(1,72) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(1,73) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(1,74) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(1,75) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(1,77) = -0.3162277660168379*cMOther[1]*mnuOther; 
  data->AEM_S(2,70) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(2,71) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(2,72) = (-0.3162277660168379*cMOther[8]*mnuOther)-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(2,73) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(2,74) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(2,76) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(2,78) = -0.3162277660168379*cMOther[2]*mnuOther; 
  data->AEM_S(3,70) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(3,71) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(3,72) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(3,73) = (-0.3162277660168379*cMOther[9]*mnuOther)-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(3,75) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(3,76) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(3,79) = -0.3162277660168379*cMOther[3]*mnuOther; 
  data->AEM_S(4,70) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(4,71) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(4,72) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(4,74) = (-0.3162277660168379*cMOther[8]*mnuOther)-0.3162277660168379*cMOther[7]*mnuOther-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(4,75) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(4,76) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(4,77) = -0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(4,78) = -0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(5,70) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(5,71) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(5,73) = -0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(5,74) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(5,75) = (-0.3162277660168379*cMOther[9]*mnuOther)-0.3162277660168379*cMOther[7]*mnuOther-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(5,76) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(5,77) = -0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(5,79) = -0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(6,70) = -0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(6,72) = -0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(6,73) = -0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(6,74) = -0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(6,75) = -0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(6,76) = (-0.3162277660168379*cMOther[9]*mnuOther)-0.3162277660168379*cMOther[8]*mnuOther-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(6,78) = -0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(6,79) = -0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(7,70) = -0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(7,71) = -0.3162277660168379*cMOther[1]*mnuOther; 
  data->AEM_S(7,74) = -0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(7,75) = -0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(7,77) = (-0.2258769757263128*cMOther[7]*mnuOther)-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(8,70) = -0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(8,72) = -0.3162277660168379*cMOther[2]*mnuOther; 
  data->AEM_S(8,74) = -0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(8,76) = -0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(8,78) = (-0.2258769757263128*cMOther[8]*mnuOther)-0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(9,70) = -0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(9,73) = -0.3162277660168379*cMOther[3]*mnuOther; 
  data->AEM_S(9,75) = -0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(9,76) = -0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(9,79) = (-0.2258769757263128*cMOther[9]*mnuOther)-0.3535533905932737*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(30,0) = 0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(30,1) = 0.3535533905932737*m1rSelf[1]*mnuSelf; 
  data->AEM_S(30,2) = 0.3535533905932737*m1rSelf[2]*mnuSelf; 
  data->AEM_S(30,3) = 0.3535533905932737*m1rSelf[3]*mnuSelf; 
  data->AEM_S(30,4) = 0.3535533905932737*m1rSelf[4]*mnuSelf; 
  data->AEM_S(30,5) = 0.3535533905932737*m1rSelf[5]*mnuSelf; 
  data->AEM_S(30,6) = 0.3535533905932737*m1rSelf[6]*mnuSelf; 
  data->AEM_S(30,7) = 0.3535533905932737*m1rSelf[7]*mnuSelf; 
  data->AEM_S(30,8) = 0.3535533905932737*m1rSelf[8]*mnuSelf; 
  data->AEM_S(30,9) = 0.3535533905932737*m1rSelf[9]*mnuSelf; 
  data->AEM_S(31,0) = 0.3535533905932737*m1rSelf[1]*mnuSelf; 
  data->AEM_S(31,1) = 0.3162277660168379*m1rSelf[7]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(31,2) = 0.3535533905932737*m1rSelf[4]*mnuSelf; 
  data->AEM_S(31,3) = 0.3535533905932737*m1rSelf[5]*mnuSelf; 
  data->AEM_S(31,4) = 0.3535533905932737*m1rSelf[2]*mnuSelf; 
  data->AEM_S(31,5) = 0.3535533905932737*m1rSelf[3]*mnuSelf; 
  data->AEM_S(31,7) = 0.3162277660168379*m1rSelf[1]*mnuSelf; 
  data->AEM_S(32,0) = 0.3535533905932737*m1rSelf[2]*mnuSelf; 
  data->AEM_S(32,1) = 0.3535533905932737*m1rSelf[4]*mnuSelf; 
  data->AEM_S(32,2) = 0.3162277660168379*m1rSelf[8]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(32,3) = 0.3535533905932737*m1rSelf[6]*mnuSelf; 
  data->AEM_S(32,4) = 0.3535533905932737*m1rSelf[1]*mnuSelf; 
  data->AEM_S(32,6) = 0.3535533905932737*m1rSelf[3]*mnuSelf; 
  data->AEM_S(32,8) = 0.3162277660168379*m1rSelf[2]*mnuSelf; 
  data->AEM_S(33,0) = 0.3535533905932737*m1rSelf[3]*mnuSelf; 
  data->AEM_S(33,1) = 0.3535533905932737*m1rSelf[5]*mnuSelf; 
  data->AEM_S(33,2) = 0.3535533905932737*m1rSelf[6]*mnuSelf; 
  data->AEM_S(33,3) = 0.3162277660168379*m1rSelf[9]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(33,5) = 0.3535533905932737*m1rSelf[1]*mnuSelf; 
  data->AEM_S(33,6) = 0.3535533905932737*m1rSelf[2]*mnuSelf; 
  data->AEM_S(33,9) = 0.3162277660168379*m1rSelf[3]*mnuSelf; 
  data->AEM_S(34,0) = 0.3535533905932737*m1rSelf[4]*mnuSelf; 
  data->AEM_S(34,1) = 0.3535533905932737*m1rSelf[2]*mnuSelf; 
  data->AEM_S(34,2) = 0.3535533905932737*m1rSelf[1]*mnuSelf; 
  data->AEM_S(34,4) = 0.3162277660168379*m1rSelf[8]*mnuSelf+0.3162277660168379*m1rSelf[7]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(34,5) = 0.3535533905932737*m1rSelf[6]*mnuSelf; 
  data->AEM_S(34,6) = 0.3535533905932737*m1rSelf[5]*mnuSelf; 
  data->AEM_S(34,7) = 0.3162277660168379*m1rSelf[4]*mnuSelf; 
  data->AEM_S(34,8) = 0.3162277660168379*m1rSelf[4]*mnuSelf; 
  data->AEM_S(35,0) = 0.3535533905932737*m1rSelf[5]*mnuSelf; 
  data->AEM_S(35,1) = 0.3535533905932737*m1rSelf[3]*mnuSelf; 
  data->AEM_S(35,3) = 0.3535533905932737*m1rSelf[1]*mnuSelf; 
  data->AEM_S(35,4) = 0.3535533905932737*m1rSelf[6]*mnuSelf; 
  data->AEM_S(35,5) = 0.3162277660168379*m1rSelf[9]*mnuSelf+0.3162277660168379*m1rSelf[7]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(35,6) = 0.3535533905932737*m1rSelf[4]*mnuSelf; 
  data->AEM_S(35,7) = 0.3162277660168379*m1rSelf[5]*mnuSelf; 
  data->AEM_S(35,9) = 0.3162277660168379*m1rSelf[5]*mnuSelf; 
  data->AEM_S(36,0) = 0.3535533905932737*m1rSelf[6]*mnuSelf; 
  data->AEM_S(36,2) = 0.3535533905932737*m1rSelf[3]*mnuSelf; 
  data->AEM_S(36,3) = 0.3535533905932737*m1rSelf[2]*mnuSelf; 
  data->AEM_S(36,4) = 0.3535533905932737*m1rSelf[5]*mnuSelf; 
  data->AEM_S(36,5) = 0.3535533905932737*m1rSelf[4]*mnuSelf; 
  data->AEM_S(36,6) = 0.3162277660168379*m1rSelf[9]*mnuSelf+0.3162277660168379*m1rSelf[8]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(36,8) = 0.3162277660168379*m1rSelf[6]*mnuSelf; 
  data->AEM_S(36,9) = 0.3162277660168379*m1rSelf[6]*mnuSelf; 
  data->AEM_S(37,0) = 0.3535533905932737*m1rSelf[7]*mnuSelf; 
  data->AEM_S(37,1) = 0.3162277660168379*m1rSelf[1]*mnuSelf; 
  data->AEM_S(37,4) = 0.3162277660168379*m1rSelf[4]*mnuSelf; 
  data->AEM_S(37,5) = 0.3162277660168379*m1rSelf[5]*mnuSelf; 
  data->AEM_S(37,7) = 0.2258769757263128*m1rSelf[7]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(38,0) = 0.3535533905932737*m1rSelf[8]*mnuSelf; 
  data->AEM_S(38,2) = 0.3162277660168379*m1rSelf[2]*mnuSelf; 
  data->AEM_S(38,4) = 0.3162277660168379*m1rSelf[4]*mnuSelf; 
  data->AEM_S(38,6) = 0.3162277660168379*m1rSelf[6]*mnuSelf; 
  data->AEM_S(38,8) = 0.2258769757263128*m1rSelf[8]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(39,0) = 0.3535533905932737*m1rSelf[9]*mnuSelf; 
  data->AEM_S(39,3) = 0.3162277660168379*m1rSelf[3]*mnuSelf; 
  data->AEM_S(39,5) = 0.3162277660168379*m1rSelf[5]*mnuSelf; 
  data->AEM_S(39,6) = 0.3162277660168379*m1rSelf[6]*mnuSelf; 
  data->AEM_S(39,9) = 0.2258769757263128*m1rSelf[9]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(30,40) = 0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(30,41) = 0.3535533905932737*m1rOther[1]*mnuOther; 
  data->AEM_S(30,42) = 0.3535533905932737*m1rOther[2]*mnuOther; 
  data->AEM_S(30,43) = 0.3535533905932737*m1rOther[3]*mnuOther; 
  data->AEM_S(30,44) = 0.3535533905932737*m1rOther[4]*mnuOther; 
  data->AEM_S(30,45) = 0.3535533905932737*m1rOther[5]*mnuOther; 
  data->AEM_S(30,46) = 0.3535533905932737*m1rOther[6]*mnuOther; 
  data->AEM_S(30,47) = 0.3535533905932737*m1rOther[7]*mnuOther; 
  data->AEM_S(30,48) = 0.3535533905932737*m1rOther[8]*mnuOther; 
  data->AEM_S(30,49) = 0.3535533905932737*m1rOther[9]*mnuOther; 
  data->AEM_S(31,40) = 0.3535533905932737*m1rOther[1]*mnuOther; 
  data->AEM_S(31,41) = 0.3162277660168379*m1rOther[7]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(31,42) = 0.3535533905932737*m1rOther[4]*mnuOther; 
  data->AEM_S(31,43) = 0.3535533905932737*m1rOther[5]*mnuOther; 
  data->AEM_S(31,44) = 0.3535533905932737*m1rOther[2]*mnuOther; 
  data->AEM_S(31,45) = 0.3535533905932737*m1rOther[3]*mnuOther; 
  data->AEM_S(31,47) = 0.3162277660168379*m1rOther[1]*mnuOther; 
  data->AEM_S(32,40) = 0.3535533905932737*m1rOther[2]*mnuOther; 
  data->AEM_S(32,41) = 0.3535533905932737*m1rOther[4]*mnuOther; 
  data->AEM_S(32,42) = 0.3162277660168379*m1rOther[8]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(32,43) = 0.3535533905932737*m1rOther[6]*mnuOther; 
  data->AEM_S(32,44) = 0.3535533905932737*m1rOther[1]*mnuOther; 
  data->AEM_S(32,46) = 0.3535533905932737*m1rOther[3]*mnuOther; 
  data->AEM_S(32,48) = 0.3162277660168379*m1rOther[2]*mnuOther; 
  data->AEM_S(33,40) = 0.3535533905932737*m1rOther[3]*mnuOther; 
  data->AEM_S(33,41) = 0.3535533905932737*m1rOther[5]*mnuOther; 
  data->AEM_S(33,42) = 0.3535533905932737*m1rOther[6]*mnuOther; 
  data->AEM_S(33,43) = 0.3162277660168379*m1rOther[9]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(33,45) = 0.3535533905932737*m1rOther[1]*mnuOther; 
  data->AEM_S(33,46) = 0.3535533905932737*m1rOther[2]*mnuOther; 
  data->AEM_S(33,49) = 0.3162277660168379*m1rOther[3]*mnuOther; 
  data->AEM_S(34,40) = 0.3535533905932737*m1rOther[4]*mnuOther; 
  data->AEM_S(34,41) = 0.3535533905932737*m1rOther[2]*mnuOther; 
  data->AEM_S(34,42) = 0.3535533905932737*m1rOther[1]*mnuOther; 
  data->AEM_S(34,44) = 0.3162277660168379*m1rOther[8]*mnuOther+0.3162277660168379*m1rOther[7]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(34,45) = 0.3535533905932737*m1rOther[6]*mnuOther; 
  data->AEM_S(34,46) = 0.3535533905932737*m1rOther[5]*mnuOther; 
  data->AEM_S(34,47) = 0.3162277660168379*m1rOther[4]*mnuOther; 
  data->AEM_S(34,48) = 0.3162277660168379*m1rOther[4]*mnuOther; 
  data->AEM_S(35,40) = 0.3535533905932737*m1rOther[5]*mnuOther; 
  data->AEM_S(35,41) = 0.3535533905932737*m1rOther[3]*mnuOther; 
  data->AEM_S(35,43) = 0.3535533905932737*m1rOther[1]*mnuOther; 
  data->AEM_S(35,44) = 0.3535533905932737*m1rOther[6]*mnuOther; 
  data->AEM_S(35,45) = 0.3162277660168379*m1rOther[9]*mnuOther+0.3162277660168379*m1rOther[7]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(35,46) = 0.3535533905932737*m1rOther[4]*mnuOther; 
  data->AEM_S(35,47) = 0.3162277660168379*m1rOther[5]*mnuOther; 
  data->AEM_S(35,49) = 0.3162277660168379*m1rOther[5]*mnuOther; 
  data->AEM_S(36,40) = 0.3535533905932737*m1rOther[6]*mnuOther; 
  data->AEM_S(36,42) = 0.3535533905932737*m1rOther[3]*mnuOther; 
  data->AEM_S(36,43) = 0.3535533905932737*m1rOther[2]*mnuOther; 
  data->AEM_S(36,44) = 0.3535533905932737*m1rOther[5]*mnuOther; 
  data->AEM_S(36,45) = 0.3535533905932737*m1rOther[4]*mnuOther; 
  data->AEM_S(36,46) = 0.3162277660168379*m1rOther[9]*mnuOther+0.3162277660168379*m1rOther[8]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(36,48) = 0.3162277660168379*m1rOther[6]*mnuOther; 
  data->AEM_S(36,49) = 0.3162277660168379*m1rOther[6]*mnuOther; 
  data->AEM_S(37,40) = 0.3535533905932737*m1rOther[7]*mnuOther; 
  data->AEM_S(37,41) = 0.3162277660168379*m1rOther[1]*mnuOther; 
  data->AEM_S(37,44) = 0.3162277660168379*m1rOther[4]*mnuOther; 
  data->AEM_S(37,45) = 0.3162277660168379*m1rOther[5]*mnuOther; 
  data->AEM_S(37,47) = 0.2258769757263128*m1rOther[7]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(38,40) = 0.3535533905932737*m1rOther[8]*mnuOther; 
  data->AEM_S(38,42) = 0.3162277660168379*m1rOther[2]*mnuOther; 
  data->AEM_S(38,44) = 0.3162277660168379*m1rOther[4]*mnuOther; 
  data->AEM_S(38,46) = 0.3162277660168379*m1rOther[6]*mnuOther; 
  data->AEM_S(38,48) = 0.2258769757263128*m1rOther[8]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(39,40) = 0.3535533905932737*m1rOther[9]*mnuOther; 
  data->AEM_S(39,43) = 0.3162277660168379*m1rOther[3]*mnuOther; 
  data->AEM_S(39,45) = 0.3162277660168379*m1rOther[5]*mnuOther; 
  data->AEM_S(39,46) = 0.3162277660168379*m1rOther[6]*mnuOther; 
  data->AEM_S(39,49) = 0.2258769757263128*m1rOther[9]*mnuOther+0.3535533905932737*m1rOther[0]*mnuOther; 
 
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
  data->AEM_S(10,10) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(10,11) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(10,12) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,13) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(10,14) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(10,15) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(10,16) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(10,17) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(10,18) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(10,19) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(11,10) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,11) = 0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(11,12) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(11,13) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(11,14) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,15) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,17) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,10) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,11) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(12,12) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(12,14) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,16) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(12,18) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,10) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(13,11) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(13,12) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(13,13) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(13,15) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,16) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,19) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(14,10) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(14,11) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,12) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,14) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(14,15) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(14,16) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(14,17) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(14,18) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(15,10) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(15,11) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,13) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,14) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(15,15) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(15,16) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(15,17) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(15,19) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(16,10) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(16,12) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,13) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,14) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(16,15) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(16,16) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(16,18) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(16,19) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(17,10) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(17,11) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,14) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(17,15) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(17,17) = 0.2258769757263128*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,10) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(18,12) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,14) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(18,16) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(18,18) = 0.2258769757263128*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(19,10) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(19,13) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,15) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(19,16) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(19,19) = 0.2258769757263128*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(10,30) = -0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,31) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(10,32) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(10,33) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(10,34) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(10,35) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(10,36) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(10,37) = -0.3535533905932737*cMSelf[17]*mnuSelf; 
  data->AEM_S(10,38) = -0.3535533905932737*cMSelf[18]*mnuSelf; 
  data->AEM_S(10,39) = -0.3535533905932737*cMSelf[19]*mnuSelf; 
  data->AEM_S(11,30) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,31) = (-0.3162277660168379*cMSelf[17]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(11,32) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(11,33) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(11,34) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(11,35) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(11,37) = -0.3162277660168379*cMSelf[11]*mnuSelf; 
  data->AEM_S(12,30) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(12,31) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(12,32) = (-0.3162277660168379*cMSelf[18]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(12,33) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(12,34) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(12,36) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(12,38) = -0.3162277660168379*cMSelf[12]*mnuSelf; 
  data->AEM_S(13,30) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(13,31) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(13,32) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(13,33) = (-0.3162277660168379*cMSelf[19]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(13,35) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(13,36) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(13,39) = -0.3162277660168379*cMSelf[13]*mnuSelf; 
  data->AEM_S(14,30) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(14,31) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(14,32) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(14,34) = (-0.3162277660168379*cMSelf[18]*mnuSelf)-0.3162277660168379*cMSelf[17]*mnuSelf-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(14,35) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(14,36) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(14,37) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(14,38) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(15,30) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(15,31) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(15,33) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(15,34) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(15,35) = (-0.3162277660168379*cMSelf[19]*mnuSelf)-0.3162277660168379*cMSelf[17]*mnuSelf-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(15,36) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(15,37) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(15,39) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(16,30) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(16,32) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(16,33) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(16,34) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(16,35) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(16,36) = (-0.3162277660168379*cMSelf[19]*mnuSelf)-0.3162277660168379*cMSelf[18]*mnuSelf-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(16,38) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(16,39) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(17,30) = -0.3535533905932737*cMSelf[17]*mnuSelf; 
  data->AEM_S(17,31) = -0.3162277660168379*cMSelf[11]*mnuSelf; 
  data->AEM_S(17,34) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(17,35) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(17,37) = (-0.2258769757263128*cMSelf[17]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(18,30) = -0.3535533905932737*cMSelf[18]*mnuSelf; 
  data->AEM_S(18,32) = -0.3162277660168379*cMSelf[12]*mnuSelf; 
  data->AEM_S(18,34) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(18,36) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(18,38) = (-0.2258769757263128*cMSelf[18]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(19,30) = -0.3535533905932737*cMSelf[19]*mnuSelf; 
  data->AEM_S(19,33) = -0.3162277660168379*cMSelf[13]*mnuSelf; 
  data->AEM_S(19,35) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(19,36) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(19,39) = (-0.2258769757263128*cMSelf[19]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(10,50) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(10,51) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(10,52) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(10,53) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(10,54) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(10,55) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(10,56) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(10,57) = 0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(10,58) = 0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(10,59) = 0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(11,50) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(11,51) = 0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(11,52) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(11,53) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(11,54) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(11,55) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(11,57) = 0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(12,50) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(12,51) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(12,52) = 0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(12,53) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(12,54) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(12,56) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(12,58) = 0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(13,50) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(13,51) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(13,52) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(13,53) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(13,55) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(13,56) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(13,59) = 0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(14,50) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(14,51) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(14,52) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(14,54) = 0.3162277660168379*m0rOther[8]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(14,55) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(14,56) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(14,57) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(14,58) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(15,50) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(15,51) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(15,53) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(15,54) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(15,55) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(15,56) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(15,57) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(15,59) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(16,50) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(16,52) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(16,53) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(16,54) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(16,55) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(16,56) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(16,58) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(16,59) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(17,50) = 0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(17,51) = 0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(17,54) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(17,55) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(17,57) = 0.2258769757263128*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(18,50) = 0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(18,52) = 0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(18,54) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(18,56) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(18,58) = 0.2258769757263128*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(19,50) = 0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(19,53) = 0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(19,55) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(19,56) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(19,59) = 0.2258769757263128*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(10,70) = -0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(10,71) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(10,72) = -0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(10,73) = -0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(10,74) = -0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(10,75) = -0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(10,76) = -0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(10,77) = -0.3535533905932737*cMOther[17]*mnuOther; 
  data->AEM_S(10,78) = -0.3535533905932737*cMOther[18]*mnuOther; 
  data->AEM_S(10,79) = -0.3535533905932737*cMOther[19]*mnuOther; 
  data->AEM_S(11,70) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(11,71) = (-0.3162277660168379*cMOther[17]*mnuOther)-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(11,72) = -0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(11,73) = -0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(11,74) = -0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(11,75) = -0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(11,77) = -0.3162277660168379*cMOther[11]*mnuOther; 
  data->AEM_S(12,70) = -0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(12,71) = -0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(12,72) = (-0.3162277660168379*cMOther[18]*mnuOther)-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(12,73) = -0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(12,74) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(12,76) = -0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(12,78) = -0.3162277660168379*cMOther[12]*mnuOther; 
  data->AEM_S(13,70) = -0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(13,71) = -0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(13,72) = -0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(13,73) = (-0.3162277660168379*cMOther[19]*mnuOther)-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(13,75) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(13,76) = -0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(13,79) = -0.3162277660168379*cMOther[13]*mnuOther; 
  data->AEM_S(14,70) = -0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(14,71) = -0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(14,72) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(14,74) = (-0.3162277660168379*cMOther[18]*mnuOther)-0.3162277660168379*cMOther[17]*mnuOther-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(14,75) = -0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(14,76) = -0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(14,77) = -0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(14,78) = -0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(15,70) = -0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(15,71) = -0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(15,73) = -0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(15,74) = -0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(15,75) = (-0.3162277660168379*cMOther[19]*mnuOther)-0.3162277660168379*cMOther[17]*mnuOther-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(15,76) = -0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(15,77) = -0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(15,79) = -0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(16,70) = -0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(16,72) = -0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(16,73) = -0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(16,74) = -0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(16,75) = -0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(16,76) = (-0.3162277660168379*cMOther[19]*mnuOther)-0.3162277660168379*cMOther[18]*mnuOther-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(16,78) = -0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(16,79) = -0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(17,70) = -0.3535533905932737*cMOther[17]*mnuOther; 
  data->AEM_S(17,71) = -0.3162277660168379*cMOther[11]*mnuOther; 
  data->AEM_S(17,74) = -0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(17,75) = -0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(17,77) = (-0.2258769757263128*cMOther[17]*mnuOther)-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(18,70) = -0.3535533905932737*cMOther[18]*mnuOther; 
  data->AEM_S(18,72) = -0.3162277660168379*cMOther[12]*mnuOther; 
  data->AEM_S(18,74) = -0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(18,76) = -0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(18,78) = (-0.2258769757263128*cMOther[18]*mnuOther)-0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(19,70) = -0.3535533905932737*cMOther[19]*mnuOther; 
  data->AEM_S(19,73) = -0.3162277660168379*cMOther[13]*mnuOther; 
  data->AEM_S(19,75) = -0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(19,76) = -0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(19,79) = (-0.2258769757263128*cMOther[19]*mnuOther)-0.3535533905932737*cMOther[10]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(30,10) = 0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(30,11) = 0.3535533905932737*m1rSelf[11]*mnuSelf; 
  data->AEM_S(30,12) = 0.3535533905932737*m1rSelf[12]*mnuSelf; 
  data->AEM_S(30,13) = 0.3535533905932737*m1rSelf[13]*mnuSelf; 
  data->AEM_S(30,14) = 0.3535533905932737*m1rSelf[14]*mnuSelf; 
  data->AEM_S(30,15) = 0.3535533905932737*m1rSelf[15]*mnuSelf; 
  data->AEM_S(30,16) = 0.3535533905932737*m1rSelf[16]*mnuSelf; 
  data->AEM_S(30,17) = 0.3535533905932737*m1rSelf[17]*mnuSelf; 
  data->AEM_S(30,18) = 0.3535533905932737*m1rSelf[18]*mnuSelf; 
  data->AEM_S(30,19) = 0.3535533905932737*m1rSelf[19]*mnuSelf; 
  data->AEM_S(31,10) = 0.3535533905932737*m1rSelf[11]*mnuSelf; 
  data->AEM_S(31,11) = 0.3162277660168379*m1rSelf[17]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(31,12) = 0.3535533905932737*m1rSelf[14]*mnuSelf; 
  data->AEM_S(31,13) = 0.3535533905932737*m1rSelf[15]*mnuSelf; 
  data->AEM_S(31,14) = 0.3535533905932737*m1rSelf[12]*mnuSelf; 
  data->AEM_S(31,15) = 0.3535533905932737*m1rSelf[13]*mnuSelf; 
  data->AEM_S(31,17) = 0.3162277660168379*m1rSelf[11]*mnuSelf; 
  data->AEM_S(32,10) = 0.3535533905932737*m1rSelf[12]*mnuSelf; 
  data->AEM_S(32,11) = 0.3535533905932737*m1rSelf[14]*mnuSelf; 
  data->AEM_S(32,12) = 0.3162277660168379*m1rSelf[18]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(32,13) = 0.3535533905932737*m1rSelf[16]*mnuSelf; 
  data->AEM_S(32,14) = 0.3535533905932737*m1rSelf[11]*mnuSelf; 
  data->AEM_S(32,16) = 0.3535533905932737*m1rSelf[13]*mnuSelf; 
  data->AEM_S(32,18) = 0.3162277660168379*m1rSelf[12]*mnuSelf; 
  data->AEM_S(33,10) = 0.3535533905932737*m1rSelf[13]*mnuSelf; 
  data->AEM_S(33,11) = 0.3535533905932737*m1rSelf[15]*mnuSelf; 
  data->AEM_S(33,12) = 0.3535533905932737*m1rSelf[16]*mnuSelf; 
  data->AEM_S(33,13) = 0.3162277660168379*m1rSelf[19]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(33,15) = 0.3535533905932737*m1rSelf[11]*mnuSelf; 
  data->AEM_S(33,16) = 0.3535533905932737*m1rSelf[12]*mnuSelf; 
  data->AEM_S(33,19) = 0.3162277660168379*m1rSelf[13]*mnuSelf; 
  data->AEM_S(34,10) = 0.3535533905932737*m1rSelf[14]*mnuSelf; 
  data->AEM_S(34,11) = 0.3535533905932737*m1rSelf[12]*mnuSelf; 
  data->AEM_S(34,12) = 0.3535533905932737*m1rSelf[11]*mnuSelf; 
  data->AEM_S(34,14) = 0.3162277660168379*m1rSelf[18]*mnuSelf+0.3162277660168379*m1rSelf[17]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(34,15) = 0.3535533905932737*m1rSelf[16]*mnuSelf; 
  data->AEM_S(34,16) = 0.3535533905932737*m1rSelf[15]*mnuSelf; 
  data->AEM_S(34,17) = 0.3162277660168379*m1rSelf[14]*mnuSelf; 
  data->AEM_S(34,18) = 0.3162277660168379*m1rSelf[14]*mnuSelf; 
  data->AEM_S(35,10) = 0.3535533905932737*m1rSelf[15]*mnuSelf; 
  data->AEM_S(35,11) = 0.3535533905932737*m1rSelf[13]*mnuSelf; 
  data->AEM_S(35,13) = 0.3535533905932737*m1rSelf[11]*mnuSelf; 
  data->AEM_S(35,14) = 0.3535533905932737*m1rSelf[16]*mnuSelf; 
  data->AEM_S(35,15) = 0.3162277660168379*m1rSelf[19]*mnuSelf+0.3162277660168379*m1rSelf[17]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(35,16) = 0.3535533905932737*m1rSelf[14]*mnuSelf; 
  data->AEM_S(35,17) = 0.3162277660168379*m1rSelf[15]*mnuSelf; 
  data->AEM_S(35,19) = 0.3162277660168379*m1rSelf[15]*mnuSelf; 
  data->AEM_S(36,10) = 0.3535533905932737*m1rSelf[16]*mnuSelf; 
  data->AEM_S(36,12) = 0.3535533905932737*m1rSelf[13]*mnuSelf; 
  data->AEM_S(36,13) = 0.3535533905932737*m1rSelf[12]*mnuSelf; 
  data->AEM_S(36,14) = 0.3535533905932737*m1rSelf[15]*mnuSelf; 
  data->AEM_S(36,15) = 0.3535533905932737*m1rSelf[14]*mnuSelf; 
  data->AEM_S(36,16) = 0.3162277660168379*m1rSelf[19]*mnuSelf+0.3162277660168379*m1rSelf[18]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(36,18) = 0.3162277660168379*m1rSelf[16]*mnuSelf; 
  data->AEM_S(36,19) = 0.3162277660168379*m1rSelf[16]*mnuSelf; 
  data->AEM_S(37,10) = 0.3535533905932737*m1rSelf[17]*mnuSelf; 
  data->AEM_S(37,11) = 0.3162277660168379*m1rSelf[11]*mnuSelf; 
  data->AEM_S(37,14) = 0.3162277660168379*m1rSelf[14]*mnuSelf; 
  data->AEM_S(37,15) = 0.3162277660168379*m1rSelf[15]*mnuSelf; 
  data->AEM_S(37,17) = 0.2258769757263128*m1rSelf[17]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(38,10) = 0.3535533905932737*m1rSelf[18]*mnuSelf; 
  data->AEM_S(38,12) = 0.3162277660168379*m1rSelf[12]*mnuSelf; 
  data->AEM_S(38,14) = 0.3162277660168379*m1rSelf[14]*mnuSelf; 
  data->AEM_S(38,16) = 0.3162277660168379*m1rSelf[16]*mnuSelf; 
  data->AEM_S(38,18) = 0.2258769757263128*m1rSelf[18]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(39,10) = 0.3535533905932737*m1rSelf[19]*mnuSelf; 
  data->AEM_S(39,13) = 0.3162277660168379*m1rSelf[13]*mnuSelf; 
  data->AEM_S(39,15) = 0.3162277660168379*m1rSelf[15]*mnuSelf; 
  data->AEM_S(39,16) = 0.3162277660168379*m1rSelf[16]*mnuSelf; 
  data->AEM_S(39,19) = 0.2258769757263128*m1rSelf[19]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(30,50) = 0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(30,51) = 0.3535533905932737*m1rOther[11]*mnuOther; 
  data->AEM_S(30,52) = 0.3535533905932737*m1rOther[12]*mnuOther; 
  data->AEM_S(30,53) = 0.3535533905932737*m1rOther[13]*mnuOther; 
  data->AEM_S(30,54) = 0.3535533905932737*m1rOther[14]*mnuOther; 
  data->AEM_S(30,55) = 0.3535533905932737*m1rOther[15]*mnuOther; 
  data->AEM_S(30,56) = 0.3535533905932737*m1rOther[16]*mnuOther; 
  data->AEM_S(30,57) = 0.3535533905932737*m1rOther[17]*mnuOther; 
  data->AEM_S(30,58) = 0.3535533905932737*m1rOther[18]*mnuOther; 
  data->AEM_S(30,59) = 0.3535533905932737*m1rOther[19]*mnuOther; 
  data->AEM_S(31,50) = 0.3535533905932737*m1rOther[11]*mnuOther; 
  data->AEM_S(31,51) = 0.3162277660168379*m1rOther[17]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(31,52) = 0.3535533905932737*m1rOther[14]*mnuOther; 
  data->AEM_S(31,53) = 0.3535533905932737*m1rOther[15]*mnuOther; 
  data->AEM_S(31,54) = 0.3535533905932737*m1rOther[12]*mnuOther; 
  data->AEM_S(31,55) = 0.3535533905932737*m1rOther[13]*mnuOther; 
  data->AEM_S(31,57) = 0.3162277660168379*m1rOther[11]*mnuOther; 
  data->AEM_S(32,50) = 0.3535533905932737*m1rOther[12]*mnuOther; 
  data->AEM_S(32,51) = 0.3535533905932737*m1rOther[14]*mnuOther; 
  data->AEM_S(32,52) = 0.3162277660168379*m1rOther[18]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(32,53) = 0.3535533905932737*m1rOther[16]*mnuOther; 
  data->AEM_S(32,54) = 0.3535533905932737*m1rOther[11]*mnuOther; 
  data->AEM_S(32,56) = 0.3535533905932737*m1rOther[13]*mnuOther; 
  data->AEM_S(32,58) = 0.3162277660168379*m1rOther[12]*mnuOther; 
  data->AEM_S(33,50) = 0.3535533905932737*m1rOther[13]*mnuOther; 
  data->AEM_S(33,51) = 0.3535533905932737*m1rOther[15]*mnuOther; 
  data->AEM_S(33,52) = 0.3535533905932737*m1rOther[16]*mnuOther; 
  data->AEM_S(33,53) = 0.3162277660168379*m1rOther[19]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(33,55) = 0.3535533905932737*m1rOther[11]*mnuOther; 
  data->AEM_S(33,56) = 0.3535533905932737*m1rOther[12]*mnuOther; 
  data->AEM_S(33,59) = 0.3162277660168379*m1rOther[13]*mnuOther; 
  data->AEM_S(34,50) = 0.3535533905932737*m1rOther[14]*mnuOther; 
  data->AEM_S(34,51) = 0.3535533905932737*m1rOther[12]*mnuOther; 
  data->AEM_S(34,52) = 0.3535533905932737*m1rOther[11]*mnuOther; 
  data->AEM_S(34,54) = 0.3162277660168379*m1rOther[18]*mnuOther+0.3162277660168379*m1rOther[17]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(34,55) = 0.3535533905932737*m1rOther[16]*mnuOther; 
  data->AEM_S(34,56) = 0.3535533905932737*m1rOther[15]*mnuOther; 
  data->AEM_S(34,57) = 0.3162277660168379*m1rOther[14]*mnuOther; 
  data->AEM_S(34,58) = 0.3162277660168379*m1rOther[14]*mnuOther; 
  data->AEM_S(35,50) = 0.3535533905932737*m1rOther[15]*mnuOther; 
  data->AEM_S(35,51) = 0.3535533905932737*m1rOther[13]*mnuOther; 
  data->AEM_S(35,53) = 0.3535533905932737*m1rOther[11]*mnuOther; 
  data->AEM_S(35,54) = 0.3535533905932737*m1rOther[16]*mnuOther; 
  data->AEM_S(35,55) = 0.3162277660168379*m1rOther[19]*mnuOther+0.3162277660168379*m1rOther[17]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(35,56) = 0.3535533905932737*m1rOther[14]*mnuOther; 
  data->AEM_S(35,57) = 0.3162277660168379*m1rOther[15]*mnuOther; 
  data->AEM_S(35,59) = 0.3162277660168379*m1rOther[15]*mnuOther; 
  data->AEM_S(36,50) = 0.3535533905932737*m1rOther[16]*mnuOther; 
  data->AEM_S(36,52) = 0.3535533905932737*m1rOther[13]*mnuOther; 
  data->AEM_S(36,53) = 0.3535533905932737*m1rOther[12]*mnuOther; 
  data->AEM_S(36,54) = 0.3535533905932737*m1rOther[15]*mnuOther; 
  data->AEM_S(36,55) = 0.3535533905932737*m1rOther[14]*mnuOther; 
  data->AEM_S(36,56) = 0.3162277660168379*m1rOther[19]*mnuOther+0.3162277660168379*m1rOther[18]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(36,58) = 0.3162277660168379*m1rOther[16]*mnuOther; 
  data->AEM_S(36,59) = 0.3162277660168379*m1rOther[16]*mnuOther; 
  data->AEM_S(37,50) = 0.3535533905932737*m1rOther[17]*mnuOther; 
  data->AEM_S(37,51) = 0.3162277660168379*m1rOther[11]*mnuOther; 
  data->AEM_S(37,54) = 0.3162277660168379*m1rOther[14]*mnuOther; 
  data->AEM_S(37,55) = 0.3162277660168379*m1rOther[15]*mnuOther; 
  data->AEM_S(37,57) = 0.2258769757263128*m1rOther[17]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(38,50) = 0.3535533905932737*m1rOther[18]*mnuOther; 
  data->AEM_S(38,52) = 0.3162277660168379*m1rOther[12]*mnuOther; 
  data->AEM_S(38,54) = 0.3162277660168379*m1rOther[14]*mnuOther; 
  data->AEM_S(38,56) = 0.3162277660168379*m1rOther[16]*mnuOther; 
  data->AEM_S(38,58) = 0.2258769757263128*m1rOther[18]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(39,50) = 0.3535533905932737*m1rOther[19]*mnuOther; 
  data->AEM_S(39,53) = 0.3162277660168379*m1rOther[13]*mnuOther; 
  data->AEM_S(39,55) = 0.3162277660168379*m1rOther[15]*mnuOther; 
  data->AEM_S(39,56) = 0.3162277660168379*m1rOther[16]*mnuOther; 
  data->AEM_S(39,59) = 0.2258769757263128*m1rOther[19]*mnuOther+0.3535533905932737*m1rOther[10]*mnuOther; 
 
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
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(20,20) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,21) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,22) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,23) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(20,24) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(20,25) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(20,26) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(20,27) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(20,28) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(20,29) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(21,20) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,21) = 0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(21,22) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(21,23) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(21,24) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,25) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,27) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,20) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,21) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(22,22) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(22,23) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(22,24) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,26) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,28) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,20) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,21) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(23,22) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(23,23) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(23,25) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,26) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,29) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(24,20) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(24,21) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,22) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(24,24) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(24,25) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(24,26) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(24,27) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(24,28) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(25,20) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(25,21) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,23) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,24) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(25,25) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(25,26) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(25,27) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(25,29) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(26,20) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(26,22) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(26,23) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,24) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(26,25) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(26,26) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(26,28) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(26,29) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(27,20) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(27,21) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,24) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(27,25) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(27,27) = 0.2258769757263128*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(28,20) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(28,22) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(28,24) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(28,26) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(28,28) = 0.2258769757263128*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(29,20) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(29,23) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,25) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(29,26) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(29,29) = 0.2258769757263128*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(20,30) = -0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(20,31) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(20,32) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(20,33) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(20,34) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(20,35) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(20,36) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(20,37) = -0.3535533905932737*cMSelf[27]*mnuSelf; 
  data->AEM_S(20,38) = -0.3535533905932737*cMSelf[28]*mnuSelf; 
  data->AEM_S(20,39) = -0.3535533905932737*cMSelf[29]*mnuSelf; 
  data->AEM_S(21,30) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(21,31) = (-0.3162277660168379*cMSelf[27]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(21,32) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(21,33) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(21,34) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(21,35) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(21,37) = -0.3162277660168379*cMSelf[21]*mnuSelf; 
  data->AEM_S(22,30) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(22,31) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(22,32) = (-0.3162277660168379*cMSelf[28]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(22,33) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(22,34) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(22,36) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(22,38) = -0.3162277660168379*cMSelf[22]*mnuSelf; 
  data->AEM_S(23,30) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(23,31) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(23,32) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(23,33) = (-0.3162277660168379*cMSelf[29]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(23,35) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(23,36) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(23,39) = -0.3162277660168379*cMSelf[23]*mnuSelf; 
  data->AEM_S(24,30) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(24,31) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(24,32) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(24,34) = (-0.3162277660168379*cMSelf[28]*mnuSelf)-0.3162277660168379*cMSelf[27]*mnuSelf-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(24,35) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(24,36) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(24,37) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(24,38) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(25,30) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(25,31) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(25,33) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(25,34) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(25,35) = (-0.3162277660168379*cMSelf[29]*mnuSelf)-0.3162277660168379*cMSelf[27]*mnuSelf-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(25,36) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(25,37) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(25,39) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(26,30) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(26,32) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(26,33) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(26,34) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(26,35) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(26,36) = (-0.3162277660168379*cMSelf[29]*mnuSelf)-0.3162277660168379*cMSelf[28]*mnuSelf-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(26,38) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(26,39) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(27,30) = -0.3535533905932737*cMSelf[27]*mnuSelf; 
  data->AEM_S(27,31) = -0.3162277660168379*cMSelf[21]*mnuSelf; 
  data->AEM_S(27,34) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(27,35) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(27,37) = (-0.2258769757263128*cMSelf[27]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(28,30) = -0.3535533905932737*cMSelf[28]*mnuSelf; 
  data->AEM_S(28,32) = -0.3162277660168379*cMSelf[22]*mnuSelf; 
  data->AEM_S(28,34) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(28,36) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(28,38) = (-0.2258769757263128*cMSelf[28]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(29,30) = -0.3535533905932737*cMSelf[29]*mnuSelf; 
  data->AEM_S(29,33) = -0.3162277660168379*cMSelf[23]*mnuSelf; 
  data->AEM_S(29,35) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(29,36) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(29,39) = (-0.2258769757263128*cMSelf[29]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(20,60) = 0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(20,61) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(20,62) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(20,63) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(20,64) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(20,65) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(20,66) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(20,67) = 0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(20,68) = 0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(20,69) = 0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(21,60) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(21,61) = 0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(21,62) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(21,63) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(21,64) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(21,65) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(21,67) = 0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(22,60) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(22,61) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(22,62) = 0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(22,63) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(22,64) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(22,66) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(22,68) = 0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(23,60) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(23,61) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(23,62) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(23,63) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(23,65) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(23,66) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(23,69) = 0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(24,60) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(24,61) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(24,62) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(24,64) = 0.3162277660168379*m0rOther[8]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(24,65) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(24,66) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(24,67) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(24,68) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(25,60) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(25,61) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(25,63) = 0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(25,64) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(25,65) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(25,66) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(25,67) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(25,69) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(26,60) = 0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(26,62) = 0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(26,63) = 0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(26,64) = 0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(26,65) = 0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(26,66) = 0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(26,68) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(26,69) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(27,60) = 0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(27,61) = 0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(27,64) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(27,65) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(27,67) = 0.2258769757263128*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(28,60) = 0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(28,62) = 0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(28,64) = 0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(28,66) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(28,68) = 0.2258769757263128*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(29,60) = 0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(29,63) = 0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(29,65) = 0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(29,66) = 0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(29,69) = 0.2258769757263128*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(20,70) = -0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(20,71) = -0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(20,72) = -0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(20,73) = -0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(20,74) = -0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(20,75) = -0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(20,76) = -0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(20,77) = -0.3535533905932737*cMOther[27]*mnuOther; 
  data->AEM_S(20,78) = -0.3535533905932737*cMOther[28]*mnuOther; 
  data->AEM_S(20,79) = -0.3535533905932737*cMOther[29]*mnuOther; 
  data->AEM_S(21,70) = -0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(21,71) = (-0.3162277660168379*cMOther[27]*mnuOther)-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(21,72) = -0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(21,73) = -0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(21,74) = -0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(21,75) = -0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(21,77) = -0.3162277660168379*cMOther[21]*mnuOther; 
  data->AEM_S(22,70) = -0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(22,71) = -0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(22,72) = (-0.3162277660168379*cMOther[28]*mnuOther)-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(22,73) = -0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(22,74) = -0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(22,76) = -0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(22,78) = -0.3162277660168379*cMOther[22]*mnuOther; 
  data->AEM_S(23,70) = -0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(23,71) = -0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(23,72) = -0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(23,73) = (-0.3162277660168379*cMOther[29]*mnuOther)-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(23,75) = -0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(23,76) = -0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(23,79) = -0.3162277660168379*cMOther[23]*mnuOther; 
  data->AEM_S(24,70) = -0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(24,71) = -0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(24,72) = -0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(24,74) = (-0.3162277660168379*cMOther[28]*mnuOther)-0.3162277660168379*cMOther[27]*mnuOther-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(24,75) = -0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(24,76) = -0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(24,77) = -0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(24,78) = -0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(25,70) = -0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(25,71) = -0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(25,73) = -0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(25,74) = -0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(25,75) = (-0.3162277660168379*cMOther[29]*mnuOther)-0.3162277660168379*cMOther[27]*mnuOther-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(25,76) = -0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(25,77) = -0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(25,79) = -0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(26,70) = -0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(26,72) = -0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(26,73) = -0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(26,74) = -0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(26,75) = -0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(26,76) = (-0.3162277660168379*cMOther[29]*mnuOther)-0.3162277660168379*cMOther[28]*mnuOther-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(26,78) = -0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(26,79) = -0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(27,70) = -0.3535533905932737*cMOther[27]*mnuOther; 
  data->AEM_S(27,71) = -0.3162277660168379*cMOther[21]*mnuOther; 
  data->AEM_S(27,74) = -0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(27,75) = -0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(27,77) = (-0.2258769757263128*cMOther[27]*mnuOther)-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(28,70) = -0.3535533905932737*cMOther[28]*mnuOther; 
  data->AEM_S(28,72) = -0.3162277660168379*cMOther[22]*mnuOther; 
  data->AEM_S(28,74) = -0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(28,76) = -0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(28,78) = (-0.2258769757263128*cMOther[28]*mnuOther)-0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(29,70) = -0.3535533905932737*cMOther[29]*mnuOther; 
  data->AEM_S(29,73) = -0.3162277660168379*cMOther[23]*mnuOther; 
  data->AEM_S(29,75) = -0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(29,76) = -0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(29,79) = (-0.2258769757263128*cMOther[29]*mnuOther)-0.3535533905932737*cMOther[20]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfZ and uCrossSelfZ ... // 
  data->AEM_S(30,20) = 0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(30,21) = 0.3535533905932737*m1rSelf[21]*mnuSelf; 
  data->AEM_S(30,22) = 0.3535533905932737*m1rSelf[22]*mnuSelf; 
  data->AEM_S(30,23) = 0.3535533905932737*m1rSelf[23]*mnuSelf; 
  data->AEM_S(30,24) = 0.3535533905932737*m1rSelf[24]*mnuSelf; 
  data->AEM_S(30,25) = 0.3535533905932737*m1rSelf[25]*mnuSelf; 
  data->AEM_S(30,26) = 0.3535533905932737*m1rSelf[26]*mnuSelf; 
  data->AEM_S(30,27) = 0.3535533905932737*m1rSelf[27]*mnuSelf; 
  data->AEM_S(30,28) = 0.3535533905932737*m1rSelf[28]*mnuSelf; 
  data->AEM_S(30,29) = 0.3535533905932737*m1rSelf[29]*mnuSelf; 
  data->AEM_S(31,20) = 0.3535533905932737*m1rSelf[21]*mnuSelf; 
  data->AEM_S(31,21) = 0.3162277660168379*m1rSelf[27]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(31,22) = 0.3535533905932737*m1rSelf[24]*mnuSelf; 
  data->AEM_S(31,23) = 0.3535533905932737*m1rSelf[25]*mnuSelf; 
  data->AEM_S(31,24) = 0.3535533905932737*m1rSelf[22]*mnuSelf; 
  data->AEM_S(31,25) = 0.3535533905932737*m1rSelf[23]*mnuSelf; 
  data->AEM_S(31,27) = 0.3162277660168379*m1rSelf[21]*mnuSelf; 
  data->AEM_S(32,20) = 0.3535533905932737*m1rSelf[22]*mnuSelf; 
  data->AEM_S(32,21) = 0.3535533905932737*m1rSelf[24]*mnuSelf; 
  data->AEM_S(32,22) = 0.3162277660168379*m1rSelf[28]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(32,23) = 0.3535533905932737*m1rSelf[26]*mnuSelf; 
  data->AEM_S(32,24) = 0.3535533905932737*m1rSelf[21]*mnuSelf; 
  data->AEM_S(32,26) = 0.3535533905932737*m1rSelf[23]*mnuSelf; 
  data->AEM_S(32,28) = 0.3162277660168379*m1rSelf[22]*mnuSelf; 
  data->AEM_S(33,20) = 0.3535533905932737*m1rSelf[23]*mnuSelf; 
  data->AEM_S(33,21) = 0.3535533905932737*m1rSelf[25]*mnuSelf; 
  data->AEM_S(33,22) = 0.3535533905932737*m1rSelf[26]*mnuSelf; 
  data->AEM_S(33,23) = 0.3162277660168379*m1rSelf[29]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(33,25) = 0.3535533905932737*m1rSelf[21]*mnuSelf; 
  data->AEM_S(33,26) = 0.3535533905932737*m1rSelf[22]*mnuSelf; 
  data->AEM_S(33,29) = 0.3162277660168379*m1rSelf[23]*mnuSelf; 
  data->AEM_S(34,20) = 0.3535533905932737*m1rSelf[24]*mnuSelf; 
  data->AEM_S(34,21) = 0.3535533905932737*m1rSelf[22]*mnuSelf; 
  data->AEM_S(34,22) = 0.3535533905932737*m1rSelf[21]*mnuSelf; 
  data->AEM_S(34,24) = 0.3162277660168379*m1rSelf[28]*mnuSelf+0.3162277660168379*m1rSelf[27]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(34,25) = 0.3535533905932737*m1rSelf[26]*mnuSelf; 
  data->AEM_S(34,26) = 0.3535533905932737*m1rSelf[25]*mnuSelf; 
  data->AEM_S(34,27) = 0.3162277660168379*m1rSelf[24]*mnuSelf; 
  data->AEM_S(34,28) = 0.3162277660168379*m1rSelf[24]*mnuSelf; 
  data->AEM_S(35,20) = 0.3535533905932737*m1rSelf[25]*mnuSelf; 
  data->AEM_S(35,21) = 0.3535533905932737*m1rSelf[23]*mnuSelf; 
  data->AEM_S(35,23) = 0.3535533905932737*m1rSelf[21]*mnuSelf; 
  data->AEM_S(35,24) = 0.3535533905932737*m1rSelf[26]*mnuSelf; 
  data->AEM_S(35,25) = 0.3162277660168379*m1rSelf[29]*mnuSelf+0.3162277660168379*m1rSelf[27]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(35,26) = 0.3535533905932737*m1rSelf[24]*mnuSelf; 
  data->AEM_S(35,27) = 0.3162277660168379*m1rSelf[25]*mnuSelf; 
  data->AEM_S(35,29) = 0.3162277660168379*m1rSelf[25]*mnuSelf; 
  data->AEM_S(36,20) = 0.3535533905932737*m1rSelf[26]*mnuSelf; 
  data->AEM_S(36,22) = 0.3535533905932737*m1rSelf[23]*mnuSelf; 
  data->AEM_S(36,23) = 0.3535533905932737*m1rSelf[22]*mnuSelf; 
  data->AEM_S(36,24) = 0.3535533905932737*m1rSelf[25]*mnuSelf; 
  data->AEM_S(36,25) = 0.3535533905932737*m1rSelf[24]*mnuSelf; 
  data->AEM_S(36,26) = 0.3162277660168379*m1rSelf[29]*mnuSelf+0.3162277660168379*m1rSelf[28]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(36,28) = 0.3162277660168379*m1rSelf[26]*mnuSelf; 
  data->AEM_S(36,29) = 0.3162277660168379*m1rSelf[26]*mnuSelf; 
  data->AEM_S(37,20) = 0.3535533905932737*m1rSelf[27]*mnuSelf; 
  data->AEM_S(37,21) = 0.3162277660168379*m1rSelf[21]*mnuSelf; 
  data->AEM_S(37,24) = 0.3162277660168379*m1rSelf[24]*mnuSelf; 
  data->AEM_S(37,25) = 0.3162277660168379*m1rSelf[25]*mnuSelf; 
  data->AEM_S(37,27) = 0.2258769757263128*m1rSelf[27]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(38,20) = 0.3535533905932737*m1rSelf[28]*mnuSelf; 
  data->AEM_S(38,22) = 0.3162277660168379*m1rSelf[22]*mnuSelf; 
  data->AEM_S(38,24) = 0.3162277660168379*m1rSelf[24]*mnuSelf; 
  data->AEM_S(38,26) = 0.3162277660168379*m1rSelf[26]*mnuSelf; 
  data->AEM_S(38,28) = 0.2258769757263128*m1rSelf[28]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(39,20) = 0.3535533905932737*m1rSelf[29]*mnuSelf; 
  data->AEM_S(39,23) = 0.3162277660168379*m1rSelf[23]*mnuSelf; 
  data->AEM_S(39,25) = 0.3162277660168379*m1rSelf[25]*mnuSelf; 
  data->AEM_S(39,26) = 0.3162277660168379*m1rSelf[26]*mnuSelf; 
  data->AEM_S(39,29) = 0.2258769757263128*m1rSelf[29]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherZ and uCrossOtherZ ... // 
  data->AEM_S(30,60) = 0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(30,61) = 0.3535533905932737*m1rOther[21]*mnuOther; 
  data->AEM_S(30,62) = 0.3535533905932737*m1rOther[22]*mnuOther; 
  data->AEM_S(30,63) = 0.3535533905932737*m1rOther[23]*mnuOther; 
  data->AEM_S(30,64) = 0.3535533905932737*m1rOther[24]*mnuOther; 
  data->AEM_S(30,65) = 0.3535533905932737*m1rOther[25]*mnuOther; 
  data->AEM_S(30,66) = 0.3535533905932737*m1rOther[26]*mnuOther; 
  data->AEM_S(30,67) = 0.3535533905932737*m1rOther[27]*mnuOther; 
  data->AEM_S(30,68) = 0.3535533905932737*m1rOther[28]*mnuOther; 
  data->AEM_S(30,69) = 0.3535533905932737*m1rOther[29]*mnuOther; 
  data->AEM_S(31,60) = 0.3535533905932737*m1rOther[21]*mnuOther; 
  data->AEM_S(31,61) = 0.3162277660168379*m1rOther[27]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(31,62) = 0.3535533905932737*m1rOther[24]*mnuOther; 
  data->AEM_S(31,63) = 0.3535533905932737*m1rOther[25]*mnuOther; 
  data->AEM_S(31,64) = 0.3535533905932737*m1rOther[22]*mnuOther; 
  data->AEM_S(31,65) = 0.3535533905932737*m1rOther[23]*mnuOther; 
  data->AEM_S(31,67) = 0.3162277660168379*m1rOther[21]*mnuOther; 
  data->AEM_S(32,60) = 0.3535533905932737*m1rOther[22]*mnuOther; 
  data->AEM_S(32,61) = 0.3535533905932737*m1rOther[24]*mnuOther; 
  data->AEM_S(32,62) = 0.3162277660168379*m1rOther[28]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(32,63) = 0.3535533905932737*m1rOther[26]*mnuOther; 
  data->AEM_S(32,64) = 0.3535533905932737*m1rOther[21]*mnuOther; 
  data->AEM_S(32,66) = 0.3535533905932737*m1rOther[23]*mnuOther; 
  data->AEM_S(32,68) = 0.3162277660168379*m1rOther[22]*mnuOther; 
  data->AEM_S(33,60) = 0.3535533905932737*m1rOther[23]*mnuOther; 
  data->AEM_S(33,61) = 0.3535533905932737*m1rOther[25]*mnuOther; 
  data->AEM_S(33,62) = 0.3535533905932737*m1rOther[26]*mnuOther; 
  data->AEM_S(33,63) = 0.3162277660168379*m1rOther[29]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(33,65) = 0.3535533905932737*m1rOther[21]*mnuOther; 
  data->AEM_S(33,66) = 0.3535533905932737*m1rOther[22]*mnuOther; 
  data->AEM_S(33,69) = 0.3162277660168379*m1rOther[23]*mnuOther; 
  data->AEM_S(34,60) = 0.3535533905932737*m1rOther[24]*mnuOther; 
  data->AEM_S(34,61) = 0.3535533905932737*m1rOther[22]*mnuOther; 
  data->AEM_S(34,62) = 0.3535533905932737*m1rOther[21]*mnuOther; 
  data->AEM_S(34,64) = 0.3162277660168379*m1rOther[28]*mnuOther+0.3162277660168379*m1rOther[27]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(34,65) = 0.3535533905932737*m1rOther[26]*mnuOther; 
  data->AEM_S(34,66) = 0.3535533905932737*m1rOther[25]*mnuOther; 
  data->AEM_S(34,67) = 0.3162277660168379*m1rOther[24]*mnuOther; 
  data->AEM_S(34,68) = 0.3162277660168379*m1rOther[24]*mnuOther; 
  data->AEM_S(35,60) = 0.3535533905932737*m1rOther[25]*mnuOther; 
  data->AEM_S(35,61) = 0.3535533905932737*m1rOther[23]*mnuOther; 
  data->AEM_S(35,63) = 0.3535533905932737*m1rOther[21]*mnuOther; 
  data->AEM_S(35,64) = 0.3535533905932737*m1rOther[26]*mnuOther; 
  data->AEM_S(35,65) = 0.3162277660168379*m1rOther[29]*mnuOther+0.3162277660168379*m1rOther[27]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(35,66) = 0.3535533905932737*m1rOther[24]*mnuOther; 
  data->AEM_S(35,67) = 0.3162277660168379*m1rOther[25]*mnuOther; 
  data->AEM_S(35,69) = 0.3162277660168379*m1rOther[25]*mnuOther; 
  data->AEM_S(36,60) = 0.3535533905932737*m1rOther[26]*mnuOther; 
  data->AEM_S(36,62) = 0.3535533905932737*m1rOther[23]*mnuOther; 
  data->AEM_S(36,63) = 0.3535533905932737*m1rOther[22]*mnuOther; 
  data->AEM_S(36,64) = 0.3535533905932737*m1rOther[25]*mnuOther; 
  data->AEM_S(36,65) = 0.3535533905932737*m1rOther[24]*mnuOther; 
  data->AEM_S(36,66) = 0.3162277660168379*m1rOther[29]*mnuOther+0.3162277660168379*m1rOther[28]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(36,68) = 0.3162277660168379*m1rOther[26]*mnuOther; 
  data->AEM_S(36,69) = 0.3162277660168379*m1rOther[26]*mnuOther; 
  data->AEM_S(37,60) = 0.3535533905932737*m1rOther[27]*mnuOther; 
  data->AEM_S(37,61) = 0.3162277660168379*m1rOther[21]*mnuOther; 
  data->AEM_S(37,64) = 0.3162277660168379*m1rOther[24]*mnuOther; 
  data->AEM_S(37,65) = 0.3162277660168379*m1rOther[25]*mnuOther; 
  data->AEM_S(37,67) = 0.2258769757263128*m1rOther[27]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(38,60) = 0.3535533905932737*m1rOther[28]*mnuOther; 
  data->AEM_S(38,62) = 0.3162277660168379*m1rOther[22]*mnuOther; 
  data->AEM_S(38,64) = 0.3162277660168379*m1rOther[24]*mnuOther; 
  data->AEM_S(38,66) = 0.3162277660168379*m1rOther[26]*mnuOther; 
  data->AEM_S(38,68) = 0.2258769757263128*m1rOther[28]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(39,60) = 0.3535533905932737*m1rOther[29]*mnuOther; 
  data->AEM_S(39,63) = 0.3162277660168379*m1rOther[23]*mnuOther; 
  data->AEM_S(39,65) = 0.3162277660168379*m1rOther[25]*mnuOther; 
  data->AEM_S(39,66) = 0.3162277660168379*m1rOther[26]*mnuOther; 
  data->AEM_S(39,69) = 0.2258769757263128*m1rOther[29]*mnuOther+0.3535533905932737*m1rOther[20]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[20] += m1rSelf[20]*mnuSelf+m1rOther[20]*mnuOther; 
  mnuM1sum[21] += m1rSelf[21]*mnuSelf+m1rOther[21]*mnuOther; 
  mnuM1sum[22] += m1rSelf[22]*mnuSelf+m1rOther[22]*mnuOther; 
  mnuM1sum[23] += m1rSelf[23]*mnuSelf+m1rOther[23]*mnuOther; 
  mnuM1sum[24] += m1rSelf[24]*mnuSelf+m1rOther[24]*mnuOther; 
  mnuM1sum[25] += m1rSelf[25]*mnuSelf+m1rOther[25]*mnuOther; 
  mnuM1sum[26] += m1rSelf[26]*mnuSelf+m1rOther[26]*mnuOther; 
  mnuM1sum[27] += m1rSelf[27]*mnuSelf+m1rOther[27]*mnuOther; 
  mnuM1sum[28] += m1rSelf[28]*mnuSelf+m1rOther[28]*mnuOther; 
  mnuM1sum[29] += m1rSelf[29]*mnuSelf+m1rOther[29]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(30,30) = 1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(30,31) = 1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(30,32) = 1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(30,33) = 1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(30,34) = 1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(30,35) = 1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(30,36) = 1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(30,37) = 1.060660171779821*m0rSelf[7]*mnuSelf-0.3535533905932737*cESelf[7]*mnuSelf; 
  data->AEM_S(30,38) = 1.060660171779821*m0rSelf[8]*mnuSelf-0.3535533905932737*cESelf[8]*mnuSelf; 
  data->AEM_S(30,39) = 1.060660171779821*m0rSelf[9]*mnuSelf-0.3535533905932737*cESelf[9]*mnuSelf; 
  data->AEM_S(31,30) = 1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(31,31) = 0.9486832980505137*m0rSelf[7]*mnuSelf-0.3162277660168379*cESelf[7]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(31,32) = 1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(31,33) = 1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(31,34) = 1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(31,35) = 1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(31,37) = 0.9486832980505137*m0rSelf[1]*mnuSelf-0.3162277660168379*cESelf[1]*mnuSelf; 
  data->AEM_S(32,30) = 1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(32,31) = 1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(32,32) = 0.9486832980505137*m0rSelf[8]*mnuSelf-0.3162277660168379*cESelf[8]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(32,33) = 1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(32,34) = 1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(32,36) = 1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(32,38) = 0.9486832980505137*m0rSelf[2]*mnuSelf-0.3162277660168379*cESelf[2]*mnuSelf; 
  data->AEM_S(33,30) = 1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(33,31) = 1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(33,32) = 1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(33,33) = 0.9486832980505137*m0rSelf[9]*mnuSelf-0.3162277660168379*cESelf[9]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(33,35) = 1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(33,36) = 1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(33,39) = 0.9486832980505137*m0rSelf[3]*mnuSelf-0.3162277660168379*cESelf[3]*mnuSelf; 
  data->AEM_S(34,30) = 1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(34,31) = 1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(34,32) = 1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(34,34) = 0.9486832980505137*m0rSelf[8]*mnuSelf-0.3162277660168379*cESelf[8]*mnuSelf+0.9486832980505137*m0rSelf[7]*mnuSelf-0.3162277660168379*cESelf[7]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(34,35) = 1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(34,36) = 1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(34,37) = 0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(34,38) = 0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(35,30) = 1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(35,31) = 1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(35,33) = 1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(35,34) = 1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(35,35) = 0.9486832980505137*m0rSelf[9]*mnuSelf-0.3162277660168379*cESelf[9]*mnuSelf+0.9486832980505137*m0rSelf[7]*mnuSelf-0.3162277660168379*cESelf[7]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(35,36) = 1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(35,37) = 0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(35,39) = 0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(36,30) = 1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(36,32) = 1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(36,33) = 1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(36,34) = 1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(36,35) = 1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(36,36) = 0.9486832980505137*m0rSelf[9]*mnuSelf-0.3162277660168379*cESelf[9]*mnuSelf+0.9486832980505137*m0rSelf[8]*mnuSelf-0.3162277660168379*cESelf[8]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(36,38) = 0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(36,39) = 0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(37,30) = 1.060660171779821*m0rSelf[7]*mnuSelf-0.3535533905932737*cESelf[7]*mnuSelf; 
  data->AEM_S(37,31) = 0.9486832980505137*m0rSelf[1]*mnuSelf-0.3162277660168379*cESelf[1]*mnuSelf; 
  data->AEM_S(37,34) = 0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(37,35) = 0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(37,37) = 0.6776309271789384*m0rSelf[7]*mnuSelf-0.2258769757263128*cESelf[7]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(38,30) = 1.060660171779821*m0rSelf[8]*mnuSelf-0.3535533905932737*cESelf[8]*mnuSelf; 
  data->AEM_S(38,32) = 0.9486832980505137*m0rSelf[2]*mnuSelf-0.3162277660168379*cESelf[2]*mnuSelf; 
  data->AEM_S(38,34) = 0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(38,36) = 0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(38,38) = 0.6776309271789384*m0rSelf[8]*mnuSelf-0.2258769757263128*cESelf[8]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(39,30) = 1.060660171779821*m0rSelf[9]*mnuSelf-0.3535533905932737*cESelf[9]*mnuSelf; 
  data->AEM_S(39,33) = 0.9486832980505137*m0rSelf[3]*mnuSelf-0.3162277660168379*cESelf[3]*mnuSelf; 
  data->AEM_S(39,35) = 0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(39,36) = 0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(39,39) = 0.6776309271789384*m0rSelf[9]*mnuSelf-0.2258769757263128*cESelf[9]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(30,70) = 1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(30,71) = 1.060660171779821*m0rOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(30,72) = 1.060660171779821*m0rOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(30,73) = 1.060660171779821*m0rOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(30,74) = 1.060660171779821*m0rOther[4]*mnuOther-0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(30,75) = 1.060660171779821*m0rOther[5]*mnuOther-0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(30,76) = 1.060660171779821*m0rOther[6]*mnuOther-0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(30,77) = 1.060660171779821*m0rOther[7]*mnuOther-0.3535533905932737*cEOther[7]*mnuOther; 
  data->AEM_S(30,78) = 1.060660171779821*m0rOther[8]*mnuOther-0.3535533905932737*cEOther[8]*mnuOther; 
  data->AEM_S(30,79) = 1.060660171779821*m0rOther[9]*mnuOther-0.3535533905932737*cEOther[9]*mnuOther; 
  data->AEM_S(31,70) = 1.060660171779821*m0rOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(31,71) = 0.9486832980505137*m0rOther[7]*mnuOther-0.3162277660168379*cEOther[7]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(31,72) = 1.060660171779821*m0rOther[4]*mnuOther-0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(31,73) = 1.060660171779821*m0rOther[5]*mnuOther-0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(31,74) = 1.060660171779821*m0rOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(31,75) = 1.060660171779821*m0rOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(31,77) = 0.9486832980505137*m0rOther[1]*mnuOther-0.3162277660168379*cEOther[1]*mnuOther; 
  data->AEM_S(32,70) = 1.060660171779821*m0rOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(32,71) = 1.060660171779821*m0rOther[4]*mnuOther-0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(32,72) = 0.9486832980505137*m0rOther[8]*mnuOther-0.3162277660168379*cEOther[8]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(32,73) = 1.060660171779821*m0rOther[6]*mnuOther-0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(32,74) = 1.060660171779821*m0rOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(32,76) = 1.060660171779821*m0rOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(32,78) = 0.9486832980505137*m0rOther[2]*mnuOther-0.3162277660168379*cEOther[2]*mnuOther; 
  data->AEM_S(33,70) = 1.060660171779821*m0rOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(33,71) = 1.060660171779821*m0rOther[5]*mnuOther-0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(33,72) = 1.060660171779821*m0rOther[6]*mnuOther-0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(33,73) = 0.9486832980505137*m0rOther[9]*mnuOther-0.3162277660168379*cEOther[9]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(33,75) = 1.060660171779821*m0rOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(33,76) = 1.060660171779821*m0rOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(33,79) = 0.9486832980505137*m0rOther[3]*mnuOther-0.3162277660168379*cEOther[3]*mnuOther; 
  data->AEM_S(34,70) = 1.060660171779821*m0rOther[4]*mnuOther-0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(34,71) = 1.060660171779821*m0rOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(34,72) = 1.060660171779821*m0rOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(34,74) = 0.9486832980505137*m0rOther[8]*mnuOther-0.3162277660168379*cEOther[8]*mnuOther+0.9486832980505137*m0rOther[7]*mnuOther-0.3162277660168379*cEOther[7]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(34,75) = 1.060660171779821*m0rOther[6]*mnuOther-0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(34,76) = 1.060660171779821*m0rOther[5]*mnuOther-0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(34,77) = 0.9486832980505137*m0rOther[4]*mnuOther-0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(34,78) = 0.9486832980505137*m0rOther[4]*mnuOther-0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(35,70) = 1.060660171779821*m0rOther[5]*mnuOther-0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(35,71) = 1.060660171779821*m0rOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(35,73) = 1.060660171779821*m0rOther[1]*mnuOther-0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(35,74) = 1.060660171779821*m0rOther[6]*mnuOther-0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(35,75) = 0.9486832980505137*m0rOther[9]*mnuOther-0.3162277660168379*cEOther[9]*mnuOther+0.9486832980505137*m0rOther[7]*mnuOther-0.3162277660168379*cEOther[7]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(35,76) = 1.060660171779821*m0rOther[4]*mnuOther-0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(35,77) = 0.9486832980505137*m0rOther[5]*mnuOther-0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(35,79) = 0.9486832980505137*m0rOther[5]*mnuOther-0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(36,70) = 1.060660171779821*m0rOther[6]*mnuOther-0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(36,72) = 1.060660171779821*m0rOther[3]*mnuOther-0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(36,73) = 1.060660171779821*m0rOther[2]*mnuOther-0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(36,74) = 1.060660171779821*m0rOther[5]*mnuOther-0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(36,75) = 1.060660171779821*m0rOther[4]*mnuOther-0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(36,76) = 0.9486832980505137*m0rOther[9]*mnuOther-0.3162277660168379*cEOther[9]*mnuOther+0.9486832980505137*m0rOther[8]*mnuOther-0.3162277660168379*cEOther[8]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(36,78) = 0.9486832980505137*m0rOther[6]*mnuOther-0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(36,79) = 0.9486832980505137*m0rOther[6]*mnuOther-0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(37,70) = 1.060660171779821*m0rOther[7]*mnuOther-0.3535533905932737*cEOther[7]*mnuOther; 
  data->AEM_S(37,71) = 0.9486832980505137*m0rOther[1]*mnuOther-0.3162277660168379*cEOther[1]*mnuOther; 
  data->AEM_S(37,74) = 0.9486832980505137*m0rOther[4]*mnuOther-0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(37,75) = 0.9486832980505137*m0rOther[5]*mnuOther-0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(37,77) = 0.6776309271789384*m0rOther[7]*mnuOther-0.2258769757263128*cEOther[7]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(38,70) = 1.060660171779821*m0rOther[8]*mnuOther-0.3535533905932737*cEOther[8]*mnuOther; 
  data->AEM_S(38,72) = 0.9486832980505137*m0rOther[2]*mnuOther-0.3162277660168379*cEOther[2]*mnuOther; 
  data->AEM_S(38,74) = 0.9486832980505137*m0rOther[4]*mnuOther-0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(38,76) = 0.9486832980505137*m0rOther[6]*mnuOther-0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(38,78) = 0.6776309271789384*m0rOther[8]*mnuOther-0.2258769757263128*cEOther[8]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(39,70) = 1.060660171779821*m0rOther[9]*mnuOther-0.3535533905932737*cEOther[9]*mnuOther; 
  data->AEM_S(39,73) = 0.9486832980505137*m0rOther[3]*mnuOther-0.3162277660168379*cEOther[3]*mnuOther; 
  data->AEM_S(39,75) = 0.9486832980505137*m0rOther[5]*mnuOther-0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(39,76) = 0.9486832980505137*m0rOther[6]*mnuOther-0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(39,79) = 0.6776309271789384*m0rOther[9]*mnuOther-0.2258769757263128*cEOther[9]*mnuOther+1.060660171779821*m0rOther[0]*mnuOther-0.3535533905932737*cEOther[0]*mnuOther; 
 
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
  data->AEM_S.block<10,20>(0,10).setZero(); 
  data->AEM_S.block<20,10>(10,0).setZero(); 
  data->AEM_S.block<10,10>(10,20).setZero(); 
  data->AEM_S.block<10,10>(20,10).setZero(); 
  data->AEM_S.block<10,20>(0,50).setZero(); 
  data->AEM_S.block<20,10>(10,40).setZero(); 
  data->AEM_S.block<10,10>(10,60).setZero(); 
  data->AEM_S.block<10,10>(20,50).setZero(); 
 
  double m1Relax[30]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<30; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  double m1EffD[30]; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(40,0) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(40,1) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(40,2) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(40,3) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(40,4) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(40,5) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(40,6) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(40,7) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(40,8) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(40,9) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(41,0) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(41,1) = 0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(41,2) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(41,3) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(41,4) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(41,5) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(41,7) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(42,0) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(42,1) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(42,2) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(42,3) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(42,4) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(42,6) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(42,8) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(43,0) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,1) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(43,2) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(43,3) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(43,5) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(43,6) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(43,9) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(44,0) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(44,1) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(44,2) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(44,4) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(44,5) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(44,6) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(44,7) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(44,8) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(45,0) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(45,1) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(45,3) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(45,4) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(45,5) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(45,6) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(45,7) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(45,9) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(46,0) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(46,2) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(46,3) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(46,4) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(46,5) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(46,6) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(46,8) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(46,9) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(47,0) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(47,1) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(47,4) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(47,5) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(47,7) = 0.2258769757263128*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(48,0) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(48,2) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(48,4) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(48,6) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(48,8) = 0.2258769757263128*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(49,0) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(49,3) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(49,5) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(49,6) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(49,9) = 0.2258769757263128*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(40,30) = -0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(40,31) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(40,32) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(40,33) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(40,34) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(40,35) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(40,36) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(40,37) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(40,38) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(40,39) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(41,30) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(41,31) = (-0.3162277660168379*cMSelf[7]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(41,32) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(41,33) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(41,34) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(41,35) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(41,37) = -0.3162277660168379*cMSelf[1]*mnuSelf; 
  data->AEM_S(42,30) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(42,31) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(42,32) = (-0.3162277660168379*cMSelf[8]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(42,33) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(42,34) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(42,36) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(42,38) = -0.3162277660168379*cMSelf[2]*mnuSelf; 
  data->AEM_S(43,30) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(43,31) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(43,32) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(43,33) = (-0.3162277660168379*cMSelf[9]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(43,35) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(43,36) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(43,39) = -0.3162277660168379*cMSelf[3]*mnuSelf; 
  data->AEM_S(44,30) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(44,31) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(44,32) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(44,34) = (-0.3162277660168379*cMSelf[8]*mnuSelf)-0.3162277660168379*cMSelf[7]*mnuSelf-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(44,35) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(44,36) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(44,37) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(44,38) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(45,30) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(45,31) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(45,33) = -0.3535533905932737*cMSelf[1]*mnuSelf; 
  data->AEM_S(45,34) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(45,35) = (-0.3162277660168379*cMSelf[9]*mnuSelf)-0.3162277660168379*cMSelf[7]*mnuSelf-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(45,36) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(45,37) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(45,39) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(46,30) = -0.3535533905932737*cMSelf[6]*mnuSelf; 
  data->AEM_S(46,32) = -0.3535533905932737*cMSelf[3]*mnuSelf; 
  data->AEM_S(46,33) = -0.3535533905932737*cMSelf[2]*mnuSelf; 
  data->AEM_S(46,34) = -0.3535533905932737*cMSelf[5]*mnuSelf; 
  data->AEM_S(46,35) = -0.3535533905932737*cMSelf[4]*mnuSelf; 
  data->AEM_S(46,36) = (-0.3162277660168379*cMSelf[9]*mnuSelf)-0.3162277660168379*cMSelf[8]*mnuSelf-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(46,38) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(46,39) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(47,30) = -0.3535533905932737*cMSelf[7]*mnuSelf; 
  data->AEM_S(47,31) = -0.3162277660168379*cMSelf[1]*mnuSelf; 
  data->AEM_S(47,34) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(47,35) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(47,37) = (-0.2258769757263128*cMSelf[7]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(48,30) = -0.3535533905932737*cMSelf[8]*mnuSelf; 
  data->AEM_S(48,32) = -0.3162277660168379*cMSelf[2]*mnuSelf; 
  data->AEM_S(48,34) = -0.3162277660168379*cMSelf[4]*mnuSelf; 
  data->AEM_S(48,36) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(48,38) = (-0.2258769757263128*cMSelf[8]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
  data->AEM_S(49,30) = -0.3535533905932737*cMSelf[9]*mnuSelf; 
  data->AEM_S(49,33) = -0.3162277660168379*cMSelf[3]*mnuSelf; 
  data->AEM_S(49,35) = -0.3162277660168379*cMSelf[5]*mnuSelf; 
  data->AEM_S(49,36) = -0.3162277660168379*cMSelf[6]*mnuSelf; 
  data->AEM_S(49,39) = (-0.2258769757263128*cMSelf[9]*mnuSelf)-0.3535533905932737*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(40,40) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(40,41) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(40,42) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(40,43) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(40,44) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(40,45) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(40,46) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(40,47) = -0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(40,48) = -0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(40,49) = -0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(41,40) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(41,41) = (-0.3162277660168379*m0rOther[7]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(41,42) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(41,43) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(41,44) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(41,45) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(41,47) = -0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(42,40) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(42,41) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(42,42) = (-0.3162277660168379*m0rOther[8]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(42,43) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(42,44) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(42,46) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(42,48) = -0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(43,40) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(43,41) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(43,42) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(43,43) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(43,45) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(43,46) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(43,49) = -0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(44,40) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(44,41) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(44,42) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(44,44) = (-0.3162277660168379*m0rOther[8]*mnuOther)-0.3162277660168379*m0rOther[7]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(44,45) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(44,46) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(44,47) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(44,48) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(45,40) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(45,41) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(45,43) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(45,44) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(45,45) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3162277660168379*m0rOther[7]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(45,46) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(45,47) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(45,49) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(46,40) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(46,42) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(46,43) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(46,44) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(46,45) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(46,46) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3162277660168379*m0rOther[8]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(46,48) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(46,49) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(47,40) = -0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(47,41) = -0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(47,44) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(47,45) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(47,47) = (-0.2258769757263128*m0rOther[7]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(48,40) = -0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(48,42) = -0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(48,44) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(48,46) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(48,48) = (-0.2258769757263128*m0rOther[8]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(49,40) = -0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(49,43) = -0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(49,45) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(49,46) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(49,49) = (-0.2258769757263128*m0rOther[9]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(40,70) = 0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(40,71) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(40,72) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(40,73) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(40,74) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(40,75) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(40,76) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(40,77) = 0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(40,78) = 0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(40,79) = 0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(41,70) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(41,71) = 0.3162277660168379*cMOther[7]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(41,72) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(41,73) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(41,74) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(41,75) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(41,77) = 0.3162277660168379*cMOther[1]*mnuOther; 
  data->AEM_S(42,70) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(42,71) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(42,72) = 0.3162277660168379*cMOther[8]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(42,73) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(42,74) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(42,76) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(42,78) = 0.3162277660168379*cMOther[2]*mnuOther; 
  data->AEM_S(43,70) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(43,71) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(43,72) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(43,73) = 0.3162277660168379*cMOther[9]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(43,75) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(43,76) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(43,79) = 0.3162277660168379*cMOther[3]*mnuOther; 
  data->AEM_S(44,70) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(44,71) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(44,72) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(44,74) = 0.3162277660168379*cMOther[8]*mnuOther+0.3162277660168379*cMOther[7]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(44,75) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(44,76) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(44,77) = 0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(44,78) = 0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(45,70) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(45,71) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(45,73) = 0.3535533905932737*cMOther[1]*mnuOther; 
  data->AEM_S(45,74) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(45,75) = 0.3162277660168379*cMOther[9]*mnuOther+0.3162277660168379*cMOther[7]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(45,76) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(45,77) = 0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(45,79) = 0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(46,70) = 0.3535533905932737*cMOther[6]*mnuOther; 
  data->AEM_S(46,72) = 0.3535533905932737*cMOther[3]*mnuOther; 
  data->AEM_S(46,73) = 0.3535533905932737*cMOther[2]*mnuOther; 
  data->AEM_S(46,74) = 0.3535533905932737*cMOther[5]*mnuOther; 
  data->AEM_S(46,75) = 0.3535533905932737*cMOther[4]*mnuOther; 
  data->AEM_S(46,76) = 0.3162277660168379*cMOther[9]*mnuOther+0.3162277660168379*cMOther[8]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(46,78) = 0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(46,79) = 0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(47,70) = 0.3535533905932737*cMOther[7]*mnuOther; 
  data->AEM_S(47,71) = 0.3162277660168379*cMOther[1]*mnuOther; 
  data->AEM_S(47,74) = 0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(47,75) = 0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(47,77) = 0.2258769757263128*cMOther[7]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(48,70) = 0.3535533905932737*cMOther[8]*mnuOther; 
  data->AEM_S(48,72) = 0.3162277660168379*cMOther[2]*mnuOther; 
  data->AEM_S(48,74) = 0.3162277660168379*cMOther[4]*mnuOther; 
  data->AEM_S(48,76) = 0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(48,78) = 0.2258769757263128*cMOther[8]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
  data->AEM_S(49,70) = 0.3535533905932737*cMOther[9]*mnuOther; 
  data->AEM_S(49,73) = 0.3162277660168379*cMOther[3]*mnuOther; 
  data->AEM_S(49,75) = 0.3162277660168379*cMOther[5]*mnuOther; 
  data->AEM_S(49,76) = 0.3162277660168379*cMOther[6]*mnuOther; 
  data->AEM_S(49,79) = 0.2258769757263128*cMOther[9]*mnuOther+0.3535533905932737*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(70,0) = (-0.125*m0rSelf[9]*uSelf[9]*mnuSelf)-0.125*m0rSelf[8]*uSelf[8]*mnuSelf-0.125*m0rSelf[7]*uSelf[7]*mnuSelf-0.125*m0rSelf[6]*uSelf[6]*mnuSelf-0.125*m0rSelf[5]*uSelf[5]*mnuSelf-0.125*m0rSelf[4]*uSelf[4]*mnuSelf-0.125*m0rSelf[3]*uSelf[3]*mnuSelf-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(70,1) = (-0.1118033988749895*m0rSelf[1]*uSelf[7]*mnuSelf)-0.1118033988749895*uSelf[1]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[5]*mnuSelf-0.125*uSelf[3]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[4]*mnuSelf-0.125*uSelf[2]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[1]*mnuSelf+0.3535533905932737*m1rSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(70,2) = (-0.1118033988749895*m0rSelf[2]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[2]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[3]*uSelf[6]*mnuSelf-0.125*uSelf[3]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[4]*mnuSelf-0.125*uSelf[1]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[2]*mnuSelf+0.3535533905932737*m1rSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(70,3) = (-0.1118033988749895*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[3]*m0rSelf[9]*mnuSelf-0.125*m0rSelf[2]*uSelf[6]*mnuSelf-0.125*uSelf[2]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*uSelf[1]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[3]*mnuSelf+0.3535533905932737*m1rSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(70,4) = (-0.1118033988749895*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[4]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[5]*uSelf[6]*mnuSelf-0.125*uSelf[5]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1rSelf[4]*mnuSelf-0.125*uSelf[0]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[1]*uSelf[2]*mnuSelf-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(70,5) = (-0.1118033988749895*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[5]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[4]*uSelf[6]*mnuSelf-0.125*uSelf[4]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[5]*mnuSelf+0.3535533905932737*m1rSelf[5]*mnuSelf-0.125*uSelf[0]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[3]*mnuSelf-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(70,6) = (-0.1118033988749895*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[6]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[0]*uSelf[6]*mnuSelf+0.3535533905932737*m1rSelf[6]*mnuSelf-0.125*uSelf[0]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[4]*uSelf[5]*mnuSelf-0.125*uSelf[4]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[3]*mnuSelf-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(70,7) = (-0.07985957062499249*m0rSelf[7]*uSelf[7]*mnuSelf)-0.125*m0rSelf[0]*uSelf[7]*mnuSelf+0.3535533905932737*m1rSelf[7]*mnuSelf-0.125*uSelf[0]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(70,8) = (-0.07985957062499249*m0rSelf[8]*uSelf[8]*mnuSelf)-0.125*m0rSelf[0]*uSelf[8]*mnuSelf+0.3535533905932737*m1rSelf[8]*mnuSelf-0.125*uSelf[0]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(70,9) = (-0.07985957062499249*m0rSelf[9]*uSelf[9]*mnuSelf)-0.125*m0rSelf[0]*uSelf[9]*mnuSelf+0.3535533905932737*m1rSelf[9]*mnuSelf-0.125*uSelf[0]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(71,0) = (-0.1118033988749895*m0rSelf[1]*uSelf[7]*mnuSelf)-0.1118033988749895*uSelf[1]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[5]*mnuSelf-0.125*uSelf[3]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[4]*mnuSelf-0.125*uSelf[2]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[1]*mnuSelf+0.3535533905932737*m1rSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(71,1) = (-0.125*m0rSelf[9]*uSelf[9]*mnuSelf)-0.125*m0rSelf[8]*uSelf[8]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[7]*mnuSelf+0.3162277660168379*m1rSelf[7]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[6]*uSelf[6]*mnuSelf-0.225*m0rSelf[5]*uSelf[5]*mnuSelf-0.225*m0rSelf[4]*uSelf[4]*mnuSelf-0.125*m0rSelf[3]*uSelf[3]*mnuSelf-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.225*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(71,2) = (-0.1118033988749895*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[4]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[5]*uSelf[6]*mnuSelf-0.125*uSelf[5]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1rSelf[4]*mnuSelf-0.125*uSelf[0]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[1]*uSelf[2]*mnuSelf-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(71,3) = (-0.1118033988749895*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[5]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[4]*uSelf[6]*mnuSelf-0.125*uSelf[4]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[5]*mnuSelf+0.3535533905932737*m1rSelf[5]*mnuSelf-0.125*uSelf[0]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[3]*mnuSelf-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(71,4) = (-0.1118033988749895*m0rSelf[2]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[2]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[6]*mnuSelf-0.125*uSelf[3]*m0rSelf[6]*mnuSelf-0.225*m0rSelf[1]*uSelf[4]*mnuSelf-0.225*uSelf[1]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[2]*mnuSelf+0.3535533905932737*m1rSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(71,5) = (-0.1118033988749895*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[3]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[2]*uSelf[6]*mnuSelf-0.125*uSelf[2]*m0rSelf[6]*mnuSelf-0.225*m0rSelf[1]*uSelf[5]*mnuSelf-0.225*uSelf[1]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[3]*mnuSelf+0.3535533905932737*m1rSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(71,6) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*uSelf[1]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[5]*mnuSelf-0.125*uSelf[2]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf-0.125*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(71,7) = (-0.1964285714285714*m0rSelf[1]*uSelf[7]*mnuSelf)-0.1964285714285714*uSelf[1]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[1]*mnuSelf+0.3162277660168379*m1rSelf[1]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(71,8) = (-0.125*m0rSelf[1]*uSelf[8]*mnuSelf)-0.125*uSelf[1]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(71,9) = (-0.125*m0rSelf[1]*uSelf[9]*mnuSelf)-0.125*uSelf[1]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(72,0) = (-0.1118033988749895*m0rSelf[2]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[2]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[3]*uSelf[6]*mnuSelf-0.125*uSelf[3]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[4]*mnuSelf-0.125*uSelf[1]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[2]*mnuSelf+0.3535533905932737*m1rSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(72,1) = (-0.1118033988749895*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[4]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[5]*uSelf[6]*mnuSelf-0.125*uSelf[5]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1rSelf[4]*mnuSelf-0.125*uSelf[0]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[1]*uSelf[2]*mnuSelf-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(72,2) = (-0.125*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1964285714285714*m0rSelf[8]*uSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[8]*mnuSelf+0.3162277660168379*m1rSelf[8]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[7]*uSelf[7]*mnuSelf-0.225*m0rSelf[6]*uSelf[6]*mnuSelf-0.125*m0rSelf[5]*uSelf[5]*mnuSelf-0.225*m0rSelf[4]*uSelf[4]*mnuSelf-0.125*m0rSelf[3]*uSelf[3]*mnuSelf-0.225*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(72,3) = (-0.1118033988749895*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[6]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[0]*uSelf[6]*mnuSelf+0.3535533905932737*m1rSelf[6]*mnuSelf-0.125*uSelf[0]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[4]*uSelf[5]*mnuSelf-0.125*uSelf[4]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[3]*mnuSelf-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(72,4) = (-0.1118033988749895*m0rSelf[1]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[1]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[5]*mnuSelf-0.125*uSelf[3]*m0rSelf[5]*mnuSelf-0.225*m0rSelf[2]*uSelf[4]*mnuSelf-0.225*uSelf[2]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[1]*mnuSelf+0.3535533905932737*m1rSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(72,5) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*uSelf[1]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[5]*mnuSelf-0.125*uSelf[2]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf-0.125*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(72,6) = (-0.1118033988749895*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[3]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[8]*mnuSelf-0.225*m0rSelf[2]*uSelf[6]*mnuSelf-0.225*uSelf[2]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*uSelf[1]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[3]*mnuSelf+0.3535533905932737*m1rSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(72,7) = (-0.125*m0rSelf[2]*uSelf[7]*mnuSelf)-0.125*uSelf[2]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(72,8) = (-0.1964285714285714*m0rSelf[2]*uSelf[8]*mnuSelf)-0.1964285714285714*uSelf[2]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[2]*mnuSelf+0.3162277660168379*m1rSelf[2]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(72,9) = (-0.125*m0rSelf[2]*uSelf[9]*mnuSelf)-0.125*uSelf[2]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(73,0) = (-0.1118033988749895*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[3]*m0rSelf[9]*mnuSelf-0.125*m0rSelf[2]*uSelf[6]*mnuSelf-0.125*uSelf[2]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*uSelf[1]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[3]*mnuSelf+0.3535533905932737*m1rSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(73,1) = (-0.1118033988749895*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[5]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[4]*uSelf[6]*mnuSelf-0.125*uSelf[4]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[5]*mnuSelf+0.3535533905932737*m1rSelf[5]*mnuSelf-0.125*uSelf[0]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[3]*mnuSelf-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(73,2) = (-0.1118033988749895*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[6]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[0]*uSelf[6]*mnuSelf+0.3535533905932737*m1rSelf[6]*mnuSelf-0.125*uSelf[0]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[4]*uSelf[5]*mnuSelf-0.125*uSelf[4]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[3]*mnuSelf-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(73,3) = (-0.1964285714285714*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1118033988749895*m0rSelf[0]*uSelf[9]*mnuSelf+0.3162277660168379*m1rSelf[9]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[9]*mnuSelf-0.125*m0rSelf[8]*uSelf[8]*mnuSelf-0.125*m0rSelf[7]*uSelf[7]*mnuSelf-0.225*m0rSelf[6]*uSelf[6]*mnuSelf-0.225*m0rSelf[5]*uSelf[5]*mnuSelf-0.125*m0rSelf[4]*uSelf[4]*mnuSelf-0.225*m0rSelf[3]*uSelf[3]*mnuSelf-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(73,4) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*uSelf[1]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[5]*mnuSelf-0.125*uSelf[2]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf-0.125*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(73,5) = (-0.1118033988749895*m0rSelf[1]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[1]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[3]*uSelf[5]*mnuSelf-0.225*uSelf[3]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[4]*mnuSelf-0.125*uSelf[2]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[1]*mnuSelf+0.3535533905932737*m1rSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(73,6) = (-0.1118033988749895*m0rSelf[2]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[2]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[8]*mnuSelf-0.225*m0rSelf[3]*uSelf[6]*mnuSelf-0.225*uSelf[3]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[4]*mnuSelf-0.125*uSelf[1]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[2]*mnuSelf+0.3535533905932737*m1rSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(73,7) = (-0.125*m0rSelf[3]*uSelf[7]*mnuSelf)-0.125*uSelf[3]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(73,8) = (-0.125*m0rSelf[3]*uSelf[8]*mnuSelf)-0.125*uSelf[3]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(73,9) = (-0.1964285714285714*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1964285714285714*uSelf[3]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[3]*mnuSelf+0.3162277660168379*m1rSelf[3]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(74,0) = (-0.1118033988749895*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[4]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[5]*uSelf[6]*mnuSelf-0.125*uSelf[5]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1rSelf[4]*mnuSelf-0.125*uSelf[0]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[1]*uSelf[2]*mnuSelf-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(74,1) = (-0.1118033988749895*m0rSelf[2]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[2]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[6]*mnuSelf-0.125*uSelf[3]*m0rSelf[6]*mnuSelf-0.225*m0rSelf[1]*uSelf[4]*mnuSelf-0.225*uSelf[1]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[2]*mnuSelf+0.3535533905932737*m1rSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(74,2) = (-0.1118033988749895*m0rSelf[1]*uSelf[8]*mnuSelf)-0.1118033988749895*uSelf[1]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[3]*uSelf[5]*mnuSelf-0.125*uSelf[3]*m0rSelf[5]*mnuSelf-0.225*m0rSelf[2]*uSelf[4]*mnuSelf-0.225*uSelf[2]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[1]*mnuSelf+0.3535533905932737*m1rSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(74,3) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*uSelf[1]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[5]*mnuSelf-0.125*uSelf[2]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf-0.125*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(74,4) = (-0.125*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1964285714285714*m0rSelf[8]*uSelf[8]*mnuSelf-0.1*m0rSelf[7]*uSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[8]*mnuSelf+0.3162277660168379*m1rSelf[8]*mnuSelf-0.1*uSelf[7]*m0rSelf[8]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[8]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[7]*mnuSelf+0.3162277660168379*m1rSelf[7]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[6]*uSelf[6]*mnuSelf-0.225*m0rSelf[5]*uSelf[5]*mnuSelf-0.405*m0rSelf[4]*uSelf[4]*mnuSelf-0.125*m0rSelf[3]*uSelf[3]*mnuSelf-0.225*m0rSelf[2]*uSelf[2]*mnuSelf-0.225*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(74,5) = (-0.1118033988749895*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[6]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[0]*uSelf[6]*mnuSelf+0.3535533905932737*m1rSelf[6]*mnuSelf-0.125*uSelf[0]*m0rSelf[6]*mnuSelf-0.225*m0rSelf[4]*uSelf[5]*mnuSelf-0.225*uSelf[4]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[3]*mnuSelf-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(74,6) = (-0.1118033988749895*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[5]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[4]*uSelf[6]*mnuSelf-0.225*uSelf[4]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[5]*mnuSelf+0.3535533905932737*m1rSelf[5]*mnuSelf-0.125*uSelf[0]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[3]*mnuSelf-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(74,7) = (-0.1*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1*uSelf[4]*m0rSelf[8]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[7]*mnuSelf-0.1964285714285714*uSelf[4]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[4]*mnuSelf+0.3162277660168379*m1rSelf[4]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[2]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(74,8) = (-0.1964285714285714*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1964285714285714*uSelf[4]*m0rSelf[8]*mnuSelf-0.1*m0rSelf[4]*uSelf[7]*mnuSelf-0.1*uSelf[4]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[4]*mnuSelf+0.3162277660168379*m1rSelf[4]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[2]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(74,9) = (-0.125*m0rSelf[4]*uSelf[9]*mnuSelf)-0.125*uSelf[4]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(75,0) = (-0.1118033988749895*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[5]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[4]*uSelf[6]*mnuSelf-0.125*uSelf[4]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[5]*mnuSelf+0.3535533905932737*m1rSelf[5]*mnuSelf-0.125*uSelf[0]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[3]*mnuSelf-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(75,1) = (-0.1118033988749895*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[3]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[2]*uSelf[6]*mnuSelf-0.125*uSelf[2]*m0rSelf[6]*mnuSelf-0.225*m0rSelf[1]*uSelf[5]*mnuSelf-0.225*uSelf[1]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[3]*mnuSelf+0.3535533905932737*m1rSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(75,2) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*uSelf[1]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[5]*mnuSelf-0.125*uSelf[2]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf-0.125*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(75,3) = (-0.1118033988749895*m0rSelf[1]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[1]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[3]*uSelf[5]*mnuSelf-0.225*uSelf[3]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[4]*mnuSelf-0.125*uSelf[2]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[1]*mnuSelf+0.3535533905932737*m1rSelf[1]*mnuSelf-0.125*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(75,4) = (-0.1118033988749895*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[6]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[0]*uSelf[6]*mnuSelf+0.3535533905932737*m1rSelf[6]*mnuSelf-0.125*uSelf[0]*m0rSelf[6]*mnuSelf-0.225*m0rSelf[4]*uSelf[5]*mnuSelf-0.225*uSelf[4]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[3]*mnuSelf-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(75,5) = (-0.1964285714285714*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1*m0rSelf[7]*uSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[9]*mnuSelf+0.3162277660168379*m1rSelf[9]*mnuSelf-0.1*uSelf[7]*m0rSelf[9]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[9]*mnuSelf-0.125*m0rSelf[8]*uSelf[8]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[7]*mnuSelf+0.3162277660168379*m1rSelf[7]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[6]*uSelf[6]*mnuSelf-0.405*m0rSelf[5]*uSelf[5]*mnuSelf-0.225*m0rSelf[4]*uSelf[4]*mnuSelf-0.225*m0rSelf[3]*uSelf[3]*mnuSelf-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.225*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(75,6) = (-0.1118033988749895*m0rSelf[4]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[4]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[5]*uSelf[6]*mnuSelf-0.225*uSelf[5]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1rSelf[4]*mnuSelf-0.125*uSelf[0]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[1]*uSelf[2]*mnuSelf-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(75,7) = (-0.1*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1*uSelf[5]*m0rSelf[9]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[7]*mnuSelf-0.1964285714285714*uSelf[5]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[5]*mnuSelf+0.3162277660168379*m1rSelf[5]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(75,8) = (-0.125*m0rSelf[5]*uSelf[8]*mnuSelf)-0.125*uSelf[5]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(75,9) = (-0.1964285714285714*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1964285714285714*uSelf[5]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[5]*uSelf[7]*mnuSelf-0.1*uSelf[5]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[5]*mnuSelf+0.3162277660168379*m1rSelf[5]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(76,0) = (-0.1118033988749895*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[6]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[6]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[0]*uSelf[6]*mnuSelf+0.3535533905932737*m1rSelf[6]*mnuSelf-0.125*uSelf[0]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[4]*uSelf[5]*mnuSelf-0.125*uSelf[4]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[2]*uSelf[3]*mnuSelf-0.125*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(76,1) = (-0.125*m0rSelf[1]*uSelf[6]*mnuSelf)-0.125*uSelf[1]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[2]*uSelf[5]*mnuSelf-0.125*uSelf[2]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[3]*uSelf[4]*mnuSelf-0.125*uSelf[3]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(76,2) = (-0.1118033988749895*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[3]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[8]*mnuSelf-0.225*m0rSelf[2]*uSelf[6]*mnuSelf-0.225*uSelf[2]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[5]*mnuSelf-0.125*uSelf[1]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[0]*uSelf[3]*mnuSelf+0.3535533905932737*m1rSelf[3]*mnuSelf-0.125*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(76,3) = (-0.1118033988749895*m0rSelf[2]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[2]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[8]*mnuSelf-0.225*m0rSelf[3]*uSelf[6]*mnuSelf-0.225*uSelf[3]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[1]*uSelf[4]*mnuSelf-0.125*uSelf[1]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[0]*uSelf[2]*mnuSelf+0.3535533905932737*m1rSelf[2]*mnuSelf-0.125*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(76,4) = (-0.1118033988749895*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[5]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[4]*uSelf[6]*mnuSelf-0.225*uSelf[4]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[5]*mnuSelf+0.3535533905932737*m1rSelf[5]*mnuSelf-0.125*uSelf[0]*m0rSelf[5]*mnuSelf-0.125*m0rSelf[1]*uSelf[3]*mnuSelf-0.125*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(76,5) = (-0.1118033988749895*m0rSelf[4]*uSelf[9]*mnuSelf)-0.1118033988749895*uSelf[4]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[8]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[7]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[7]*mnuSelf-0.225*m0rSelf[5]*uSelf[6]*mnuSelf-0.225*uSelf[5]*m0rSelf[6]*mnuSelf-0.125*m0rSelf[0]*uSelf[4]*mnuSelf+0.3535533905932737*m1rSelf[4]*mnuSelf-0.125*uSelf[0]*m0rSelf[4]*mnuSelf-0.125*m0rSelf[1]*uSelf[2]*mnuSelf-0.125*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(76,6) = (-0.1964285714285714*m0rSelf[9]*uSelf[9]*mnuSelf)-0.1*m0rSelf[8]*uSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[9]*mnuSelf+0.3162277660168379*m1rSelf[9]*mnuSelf-0.1*uSelf[8]*m0rSelf[9]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[9]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[8]*mnuSelf+0.3162277660168379*m1rSelf[8]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[7]*uSelf[7]*mnuSelf-0.405*m0rSelf[6]*uSelf[6]*mnuSelf-0.225*m0rSelf[5]*uSelf[5]*mnuSelf-0.225*m0rSelf[4]*uSelf[4]*mnuSelf-0.225*m0rSelf[3]*uSelf[3]*mnuSelf-0.225*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(76,7) = (-0.125*m0rSelf[6]*uSelf[7]*mnuSelf)-0.125*uSelf[6]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(76,8) = (-0.1*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1*uSelf[6]*m0rSelf[9]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[8]*mnuSelf-0.1964285714285714*uSelf[6]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[6]*mnuSelf+0.3162277660168379*m1rSelf[6]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(76,9) = (-0.1964285714285714*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1964285714285714*uSelf[6]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[6]*uSelf[8]*mnuSelf-0.1*uSelf[6]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[6]*mnuSelf+0.3162277660168379*m1rSelf[6]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(77,0) = (-0.07985957062499249*m0rSelf[7]*uSelf[7]*mnuSelf)-0.125*m0rSelf[0]*uSelf[7]*mnuSelf+0.3535533905932737*m1rSelf[7]*mnuSelf-0.125*uSelf[0]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(77,1) = (-0.1964285714285714*m0rSelf[1]*uSelf[7]*mnuSelf)-0.1964285714285714*uSelf[1]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[1]*mnuSelf+0.3162277660168379*m1rSelf[1]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(77,2) = (-0.125*m0rSelf[2]*uSelf[7]*mnuSelf)-0.125*uSelf[2]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(77,3) = (-0.125*m0rSelf[3]*uSelf[7]*mnuSelf)-0.125*uSelf[3]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(77,4) = (-0.1*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1*uSelf[4]*m0rSelf[8]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[7]*mnuSelf-0.1964285714285714*uSelf[4]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[4]*mnuSelf+0.3162277660168379*m1rSelf[4]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[2]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(77,5) = (-0.1*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1*uSelf[5]*m0rSelf[9]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[7]*mnuSelf-0.1964285714285714*uSelf[5]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[5]*mnuSelf+0.3162277660168379*m1rSelf[5]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(77,6) = (-0.125*m0rSelf[6]*uSelf[7]*mnuSelf)-0.125*uSelf[6]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(77,7) = (-0.125*m0rSelf[9]*uSelf[9]*mnuSelf)-0.125*m0rSelf[8]*uSelf[8]*mnuSelf-0.2678571428571428*m0rSelf[7]*uSelf[7]*mnuSelf-0.07985957062499249*m0rSelf[0]*uSelf[7]*mnuSelf+0.2258769757263128*m1rSelf[7]*mnuSelf-0.07985957062499249*uSelf[0]*m0rSelf[7]*mnuSelf-0.125*m0rSelf[6]*uSelf[6]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[5]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[4]*mnuSelf-0.125*m0rSelf[3]*uSelf[3]*mnuSelf-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.1964285714285714*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(77,8) = (-0.125*m0rSelf[7]*uSelf[8]*mnuSelf)-0.125*uSelf[7]*m0rSelf[8]*mnuSelf-0.1*m0rSelf[4]*uSelf[4]*mnuSelf; 
  data->AEM_S(77,9) = (-0.125*m0rSelf[7]*uSelf[9]*mnuSelf)-0.125*uSelf[7]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[5]*uSelf[5]*mnuSelf; 
  data->AEM_S(78,0) = (-0.07985957062499249*m0rSelf[8]*uSelf[8]*mnuSelf)-0.125*m0rSelf[0]*uSelf[8]*mnuSelf+0.3535533905932737*m1rSelf[8]*mnuSelf-0.125*uSelf[0]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(78,1) = (-0.125*m0rSelf[1]*uSelf[8]*mnuSelf)-0.125*uSelf[1]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[4]*mnuSelf; 
  data->AEM_S(78,2) = (-0.1964285714285714*m0rSelf[2]*uSelf[8]*mnuSelf)-0.1964285714285714*uSelf[2]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[4]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[2]*mnuSelf+0.3162277660168379*m1rSelf[2]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(78,3) = (-0.125*m0rSelf[3]*uSelf[8]*mnuSelf)-0.125*uSelf[3]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(78,4) = (-0.1964285714285714*m0rSelf[4]*uSelf[8]*mnuSelf)-0.1964285714285714*uSelf[4]*m0rSelf[8]*mnuSelf-0.1*m0rSelf[4]*uSelf[7]*mnuSelf-0.1*uSelf[4]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[4]*mnuSelf+0.3162277660168379*m1rSelf[4]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[4]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[2]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(78,5) = (-0.125*m0rSelf[5]*uSelf[8]*mnuSelf)-0.125*uSelf[5]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(78,6) = (-0.1*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1*uSelf[6]*m0rSelf[9]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[8]*mnuSelf-0.1964285714285714*uSelf[6]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[6]*mnuSelf+0.3162277660168379*m1rSelf[6]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(78,7) = (-0.125*m0rSelf[7]*uSelf[8]*mnuSelf)-0.125*uSelf[7]*m0rSelf[8]*mnuSelf-0.1*m0rSelf[4]*uSelf[4]*mnuSelf; 
  data->AEM_S(78,8) = (-0.125*m0rSelf[9]*uSelf[9]*mnuSelf)-0.2678571428571428*m0rSelf[8]*uSelf[8]*mnuSelf-0.07985957062499249*m0rSelf[0]*uSelf[8]*mnuSelf+0.2258769757263128*m1rSelf[8]*mnuSelf-0.07985957062499249*uSelf[0]*m0rSelf[8]*mnuSelf-0.125*m0rSelf[7]*uSelf[7]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[6]*mnuSelf-0.125*m0rSelf[5]*uSelf[5]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[4]*mnuSelf-0.125*m0rSelf[3]*uSelf[3]*mnuSelf-0.1964285714285714*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
  data->AEM_S(78,9) = (-0.125*m0rSelf[8]*uSelf[9]*mnuSelf)-0.125*uSelf[8]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[6]*uSelf[6]*mnuSelf; 
  data->AEM_S(79,0) = (-0.07985957062499249*m0rSelf[9]*uSelf[9]*mnuSelf)-0.125*m0rSelf[0]*uSelf[9]*mnuSelf+0.3535533905932737*m1rSelf[9]*mnuSelf-0.125*uSelf[0]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(79,1) = (-0.125*m0rSelf[1]*uSelf[9]*mnuSelf)-0.125*uSelf[1]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[5]*mnuSelf; 
  data->AEM_S(79,2) = (-0.125*m0rSelf[2]*uSelf[9]*mnuSelf)-0.125*uSelf[2]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[3]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(79,3) = (-0.1964285714285714*m0rSelf[3]*uSelf[9]*mnuSelf)-0.1964285714285714*uSelf[3]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[3]*mnuSelf+0.3162277660168379*m1rSelf[3]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(79,4) = (-0.125*m0rSelf[4]*uSelf[9]*mnuSelf)-0.125*uSelf[4]*m0rSelf[9]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[5]*m0rSelf[6]*mnuSelf; 
  data->AEM_S(79,5) = (-0.1964285714285714*m0rSelf[5]*uSelf[9]*mnuSelf)-0.1964285714285714*uSelf[5]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[5]*uSelf[7]*mnuSelf-0.1*uSelf[5]*m0rSelf[7]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[6]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[5]*mnuSelf+0.3162277660168379*m1rSelf[5]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(79,6) = (-0.1964285714285714*m0rSelf[6]*uSelf[9]*mnuSelf)-0.1964285714285714*uSelf[6]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[6]*uSelf[8]*mnuSelf-0.1*uSelf[6]*m0rSelf[8]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[6]*mnuSelf+0.3162277660168379*m1rSelf[6]*mnuSelf-0.1118033988749895*uSelf[0]*m0rSelf[6]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[5]*mnuSelf-0.1118033988749895*uSelf[4]*m0rSelf[5]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[3]*mnuSelf-0.1118033988749895*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(79,7) = (-0.125*m0rSelf[7]*uSelf[9]*mnuSelf)-0.125*uSelf[7]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[5]*uSelf[5]*mnuSelf; 
  data->AEM_S(79,8) = (-0.125*m0rSelf[8]*uSelf[9]*mnuSelf)-0.125*uSelf[8]*m0rSelf[9]*mnuSelf-0.1*m0rSelf[6]*uSelf[6]*mnuSelf; 
  data->AEM_S(79,9) = (-0.2678571428571428*m0rSelf[9]*uSelf[9]*mnuSelf)-0.07985957062499249*m0rSelf[0]*uSelf[9]*mnuSelf+0.2258769757263128*m1rSelf[9]*mnuSelf-0.07985957062499249*uSelf[0]*m0rSelf[9]*mnuSelf-0.125*m0rSelf[8]*uSelf[8]*mnuSelf-0.125*m0rSelf[7]*uSelf[7]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[6]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[5]*mnuSelf-0.125*m0rSelf[4]*uSelf[4]*mnuSelf-0.1964285714285714*m0rSelf[3]*uSelf[3]*mnuSelf-0.125*m0rSelf[2]*uSelf[2]*mnuSelf-0.125*m0rSelf[1]*uSelf[1]*mnuSelf-0.125*m0rSelf[0]*uSelf[0]*mnuSelf+0.3535533905932737*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(70,40) = 0.125*m0rOther[9]*uOther[9]*mnuOther+0.125*m0rOther[8]*uOther[8]*mnuOther+0.125*m0rOther[7]*uOther[7]*mnuOther+0.125*m0rOther[6]*uOther[6]*mnuOther+0.125*m0rOther[5]*uOther[5]*mnuOther+0.125*m0rOther[4]*uOther[4]*mnuOther+0.125*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(70,41) = 0.1118033988749895*m0rOther[1]*uOther[7]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[7]*mnuOther+0.125*m0rOther[3]*uOther[5]*mnuOther+0.125*uOther[3]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[4]*mnuOther+0.125*uOther[2]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1rOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(70,42) = 0.1118033988749895*m0rOther[2]*uOther[8]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[8]*mnuOther+0.125*m0rOther[3]*uOther[6]*mnuOther+0.125*uOther[3]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[4]*mnuOther+0.125*uOther[1]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1rOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(70,43) = 0.1118033988749895*m0rOther[3]*uOther[9]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[9]*mnuOther+0.125*m0rOther[2]*uOther[6]*mnuOther+0.125*uOther[2]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[5]*mnuOther+0.125*uOther[1]*m0rOther[5]*mnuOther+0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1rOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(70,44) = 0.1118033988749895*m0rOther[4]*uOther[8]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[7]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[7]*mnuOther+0.125*m0rOther[5]*uOther[6]*mnuOther+0.125*uOther[5]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1rOther[4]*mnuOther+0.125*uOther[0]*m0rOther[4]*mnuOther+0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(70,45) = 0.1118033988749895*m0rOther[5]*uOther[9]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[7]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[7]*mnuOther+0.125*m0rOther[4]*uOther[6]*mnuOther+0.125*uOther[4]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1rOther[5]*mnuOther+0.125*uOther[0]*m0rOther[5]*mnuOther+0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(70,46) = 0.1118033988749895*m0rOther[6]*uOther[9]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[8]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[8]*mnuOther+0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1rOther[6]*mnuOther+0.125*uOther[0]*m0rOther[6]*mnuOther+0.125*m0rOther[4]*uOther[5]*mnuOther+0.125*uOther[4]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(70,47) = 0.07985957062499249*m0rOther[7]*uOther[7]*mnuOther+0.125*m0rOther[0]*uOther[7]*mnuOther-0.3535533905932737*m1rOther[7]*mnuOther+0.125*uOther[0]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[5]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[4]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(70,48) = 0.07985957062499249*m0rOther[8]*uOther[8]*mnuOther+0.125*m0rOther[0]*uOther[8]*mnuOther-0.3535533905932737*m1rOther[8]*mnuOther+0.125*uOther[0]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[6]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[4]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(70,49) = 0.07985957062499249*m0rOther[9]*uOther[9]*mnuOther+0.125*m0rOther[0]*uOther[9]*mnuOther-0.3535533905932737*m1rOther[9]*mnuOther+0.125*uOther[0]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[6]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[5]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(71,40) = 0.1118033988749895*m0rOther[1]*uOther[7]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[7]*mnuOther+0.125*m0rOther[3]*uOther[5]*mnuOther+0.125*uOther[3]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[4]*mnuOther+0.125*uOther[2]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1rOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(71,41) = 0.125*m0rOther[9]*uOther[9]*mnuOther+0.125*m0rOther[8]*uOther[8]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[7]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[7]*mnuOther-0.3162277660168379*m1rOther[7]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[7]*mnuOther+0.125*m0rOther[6]*uOther[6]*mnuOther+0.225*m0rOther[5]*uOther[5]*mnuOther+0.225*m0rOther[4]*uOther[4]*mnuOther+0.125*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.225*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(71,42) = 0.1118033988749895*m0rOther[4]*uOther[8]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[7]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[7]*mnuOther+0.125*m0rOther[5]*uOther[6]*mnuOther+0.125*uOther[5]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1rOther[4]*mnuOther+0.125*uOther[0]*m0rOther[4]*mnuOther+0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(71,43) = 0.1118033988749895*m0rOther[5]*uOther[9]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[7]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[7]*mnuOther+0.125*m0rOther[4]*uOther[6]*mnuOther+0.125*uOther[4]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1rOther[5]*mnuOther+0.125*uOther[0]*m0rOther[5]*mnuOther+0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(71,44) = 0.1118033988749895*m0rOther[2]*uOther[8]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[7]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[7]*mnuOther+0.125*m0rOther[3]*uOther[6]*mnuOther+0.125*uOther[3]*m0rOther[6]*mnuOther+0.225*m0rOther[1]*uOther[4]*mnuOther+0.225*uOther[1]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1rOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(71,45) = 0.1118033988749895*m0rOther[3]*uOther[9]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[7]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[7]*mnuOther+0.125*m0rOther[2]*uOther[6]*mnuOther+0.125*uOther[2]*m0rOther[6]*mnuOther+0.225*m0rOther[1]*uOther[5]*mnuOther+0.225*uOther[1]*m0rOther[5]*mnuOther+0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1rOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(71,46) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*uOther[1]*m0rOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther+0.125*uOther[2]*m0rOther[5]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther+0.125*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(71,47) = 0.1964285714285714*m0rOther[1]*uOther[7]*mnuOther+0.1964285714285714*uOther[1]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[5]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[4]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[1]*mnuOther-0.3162277660168379*m1rOther[1]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(71,48) = 0.125*m0rOther[1]*uOther[8]*mnuOther+0.125*uOther[1]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[4]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[4]*mnuOther; 
  data->AEM_S(71,49) = 0.125*m0rOther[1]*uOther[9]*mnuOther+0.125*uOther[1]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[5]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[5]*mnuOther; 
  data->AEM_S(72,40) = 0.1118033988749895*m0rOther[2]*uOther[8]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[8]*mnuOther+0.125*m0rOther[3]*uOther[6]*mnuOther+0.125*uOther[3]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[4]*mnuOther+0.125*uOther[1]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1rOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(72,41) = 0.1118033988749895*m0rOther[4]*uOther[8]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[7]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[7]*mnuOther+0.125*m0rOther[5]*uOther[6]*mnuOther+0.125*uOther[5]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1rOther[4]*mnuOther+0.125*uOther[0]*m0rOther[4]*mnuOther+0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(72,42) = 0.125*m0rOther[9]*uOther[9]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[8]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[8]*mnuOther-0.3162277660168379*m1rOther[8]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[8]*mnuOther+0.125*m0rOther[7]*uOther[7]*mnuOther+0.225*m0rOther[6]*uOther[6]*mnuOther+0.125*m0rOther[5]*uOther[5]*mnuOther+0.225*m0rOther[4]*uOther[4]*mnuOther+0.125*m0rOther[3]*uOther[3]*mnuOther+0.225*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(72,43) = 0.1118033988749895*m0rOther[6]*uOther[9]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[8]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[8]*mnuOther+0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1rOther[6]*mnuOther+0.125*uOther[0]*m0rOther[6]*mnuOther+0.125*m0rOther[4]*uOther[5]*mnuOther+0.125*uOther[4]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(72,44) = 0.1118033988749895*m0rOther[1]*uOther[8]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[7]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[7]*mnuOther+0.125*m0rOther[3]*uOther[5]*mnuOther+0.125*uOther[3]*m0rOther[5]*mnuOther+0.225*m0rOther[2]*uOther[4]*mnuOther+0.225*uOther[2]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1rOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(72,45) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*uOther[1]*m0rOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther+0.125*uOther[2]*m0rOther[5]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther+0.125*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(72,46) = 0.1118033988749895*m0rOther[3]*uOther[9]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[8]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[8]*mnuOther+0.225*m0rOther[2]*uOther[6]*mnuOther+0.225*uOther[2]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[5]*mnuOther+0.125*uOther[1]*m0rOther[5]*mnuOther+0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1rOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(72,47) = 0.125*m0rOther[2]*uOther[7]*mnuOther+0.125*uOther[2]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[4]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[4]*mnuOther; 
  data->AEM_S(72,48) = 0.1964285714285714*m0rOther[2]*uOther[8]*mnuOther+0.1964285714285714*uOther[2]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[6]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[4]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[2]*mnuOther-0.3162277660168379*m1rOther[2]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(72,49) = 0.125*m0rOther[2]*uOther[9]*mnuOther+0.125*uOther[2]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[6]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[6]*mnuOther; 
  data->AEM_S(73,40) = 0.1118033988749895*m0rOther[3]*uOther[9]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[9]*mnuOther+0.125*m0rOther[2]*uOther[6]*mnuOther+0.125*uOther[2]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[5]*mnuOther+0.125*uOther[1]*m0rOther[5]*mnuOther+0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1rOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(73,41) = 0.1118033988749895*m0rOther[5]*uOther[9]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[7]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[7]*mnuOther+0.125*m0rOther[4]*uOther[6]*mnuOther+0.125*uOther[4]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1rOther[5]*mnuOther+0.125*uOther[0]*m0rOther[5]*mnuOther+0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(73,42) = 0.1118033988749895*m0rOther[6]*uOther[9]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[8]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[8]*mnuOther+0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1rOther[6]*mnuOther+0.125*uOther[0]*m0rOther[6]*mnuOther+0.125*m0rOther[4]*uOther[5]*mnuOther+0.125*uOther[4]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(73,43) = 0.1964285714285714*m0rOther[9]*uOther[9]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[9]*mnuOther-0.3162277660168379*m1rOther[9]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[9]*mnuOther+0.125*m0rOther[8]*uOther[8]*mnuOther+0.125*m0rOther[7]*uOther[7]*mnuOther+0.225*m0rOther[6]*uOther[6]*mnuOther+0.225*m0rOther[5]*uOther[5]*mnuOther+0.125*m0rOther[4]*uOther[4]*mnuOther+0.225*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(73,44) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*uOther[1]*m0rOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther+0.125*uOther[2]*m0rOther[5]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther+0.125*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(73,45) = 0.1118033988749895*m0rOther[1]*uOther[9]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[7]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[7]*mnuOther+0.225*m0rOther[3]*uOther[5]*mnuOther+0.225*uOther[3]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[4]*mnuOther+0.125*uOther[2]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1rOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(73,46) = 0.1118033988749895*m0rOther[2]*uOther[9]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[8]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[8]*mnuOther+0.225*m0rOther[3]*uOther[6]*mnuOther+0.225*uOther[3]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[4]*mnuOther+0.125*uOther[1]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1rOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(73,47) = 0.125*m0rOther[3]*uOther[7]*mnuOther+0.125*uOther[3]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[5]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[5]*mnuOther; 
  data->AEM_S(73,48) = 0.125*m0rOther[3]*uOther[8]*mnuOther+0.125*uOther[3]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[6]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[6]*mnuOther; 
  data->AEM_S(73,49) = 0.1964285714285714*m0rOther[3]*uOther[9]*mnuOther+0.1964285714285714*uOther[3]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[6]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[5]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[3]*mnuOther-0.3162277660168379*m1rOther[3]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(74,40) = 0.1118033988749895*m0rOther[4]*uOther[8]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[7]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[7]*mnuOther+0.125*m0rOther[5]*uOther[6]*mnuOther+0.125*uOther[5]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1rOther[4]*mnuOther+0.125*uOther[0]*m0rOther[4]*mnuOther+0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(74,41) = 0.1118033988749895*m0rOther[2]*uOther[8]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[7]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[7]*mnuOther+0.125*m0rOther[3]*uOther[6]*mnuOther+0.125*uOther[3]*m0rOther[6]*mnuOther+0.225*m0rOther[1]*uOther[4]*mnuOther+0.225*uOther[1]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1rOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(74,42) = 0.1118033988749895*m0rOther[1]*uOther[8]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[7]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[7]*mnuOther+0.125*m0rOther[3]*uOther[5]*mnuOther+0.125*uOther[3]*m0rOther[5]*mnuOther+0.225*m0rOther[2]*uOther[4]*mnuOther+0.225*uOther[2]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1rOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(74,43) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*uOther[1]*m0rOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther+0.125*uOther[2]*m0rOther[5]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther+0.125*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(74,44) = 0.125*m0rOther[9]*uOther[9]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[8]*mnuOther+0.1*m0rOther[7]*uOther[8]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[8]*mnuOther-0.3162277660168379*m1rOther[8]*mnuOther+0.1*uOther[7]*m0rOther[8]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[8]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[7]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[7]*mnuOther-0.3162277660168379*m1rOther[7]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[7]*mnuOther+0.225*m0rOther[6]*uOther[6]*mnuOther+0.225*m0rOther[5]*uOther[5]*mnuOther+0.405*m0rOther[4]*uOther[4]*mnuOther+0.125*m0rOther[3]*uOther[3]*mnuOther+0.225*m0rOther[2]*uOther[2]*mnuOther+0.225*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(74,45) = 0.1118033988749895*m0rOther[6]*uOther[9]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[8]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[7]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[7]*mnuOther+0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1rOther[6]*mnuOther+0.125*uOther[0]*m0rOther[6]*mnuOther+0.225*m0rOther[4]*uOther[5]*mnuOther+0.225*uOther[4]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(74,46) = 0.1118033988749895*m0rOther[5]*uOther[9]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[8]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[7]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[7]*mnuOther+0.225*m0rOther[4]*uOther[6]*mnuOther+0.225*uOther[4]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1rOther[5]*mnuOther+0.125*uOther[0]*m0rOther[5]*mnuOther+0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(74,47) = 0.1*m0rOther[4]*uOther[8]*mnuOther+0.1*uOther[4]*m0rOther[8]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[7]*mnuOther+0.1964285714285714*uOther[4]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[6]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[4]*mnuOther-0.3162277660168379*m1rOther[4]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[2]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(74,48) = 0.1964285714285714*m0rOther[4]*uOther[8]*mnuOther+0.1964285714285714*uOther[4]*m0rOther[8]*mnuOther+0.1*m0rOther[4]*uOther[7]*mnuOther+0.1*uOther[4]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[6]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[4]*mnuOther-0.3162277660168379*m1rOther[4]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[2]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(74,49) = 0.125*m0rOther[4]*uOther[9]*mnuOther+0.125*uOther[4]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[6]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[6]*mnuOther; 
  data->AEM_S(75,40) = 0.1118033988749895*m0rOther[5]*uOther[9]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[7]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[7]*mnuOther+0.125*m0rOther[4]*uOther[6]*mnuOther+0.125*uOther[4]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1rOther[5]*mnuOther+0.125*uOther[0]*m0rOther[5]*mnuOther+0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(75,41) = 0.1118033988749895*m0rOther[3]*uOther[9]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[7]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[7]*mnuOther+0.125*m0rOther[2]*uOther[6]*mnuOther+0.125*uOther[2]*m0rOther[6]*mnuOther+0.225*m0rOther[1]*uOther[5]*mnuOther+0.225*uOther[1]*m0rOther[5]*mnuOther+0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1rOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(75,42) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*uOther[1]*m0rOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther+0.125*uOther[2]*m0rOther[5]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther+0.125*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(75,43) = 0.1118033988749895*m0rOther[1]*uOther[9]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[7]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[7]*mnuOther+0.225*m0rOther[3]*uOther[5]*mnuOther+0.225*uOther[3]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[4]*mnuOther+0.125*uOther[2]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[1]*mnuOther-0.3535533905932737*m1rOther[1]*mnuOther+0.125*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(75,44) = 0.1118033988749895*m0rOther[6]*uOther[9]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[8]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[7]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[7]*mnuOther+0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1rOther[6]*mnuOther+0.125*uOther[0]*m0rOther[6]*mnuOther+0.225*m0rOther[4]*uOther[5]*mnuOther+0.225*uOther[4]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(75,45) = 0.1964285714285714*m0rOther[9]*uOther[9]*mnuOther+0.1*m0rOther[7]*uOther[9]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[9]*mnuOther-0.3162277660168379*m1rOther[9]*mnuOther+0.1*uOther[7]*m0rOther[9]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[9]*mnuOther+0.125*m0rOther[8]*uOther[8]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[7]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[7]*mnuOther-0.3162277660168379*m1rOther[7]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[7]*mnuOther+0.225*m0rOther[6]*uOther[6]*mnuOther+0.405*m0rOther[5]*uOther[5]*mnuOther+0.225*m0rOther[4]*uOther[4]*mnuOther+0.225*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.225*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(75,46) = 0.1118033988749895*m0rOther[4]*uOther[9]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[8]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[7]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[7]*mnuOther+0.225*m0rOther[5]*uOther[6]*mnuOther+0.225*uOther[5]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1rOther[4]*mnuOther+0.125*uOther[0]*m0rOther[4]*mnuOther+0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(75,47) = 0.1*m0rOther[5]*uOther[9]*mnuOther+0.1*uOther[5]*m0rOther[9]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[7]*mnuOther+0.1964285714285714*uOther[5]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[6]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[5]*mnuOther-0.3162277660168379*m1rOther[5]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[3]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(75,48) = 0.125*m0rOther[5]*uOther[8]*mnuOther+0.125*uOther[5]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[6]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[6]*mnuOther; 
  data->AEM_S(75,49) = 0.1964285714285714*m0rOther[5]*uOther[9]*mnuOther+0.1964285714285714*uOther[5]*m0rOther[9]*mnuOther+0.1*m0rOther[5]*uOther[7]*mnuOther+0.1*uOther[5]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[6]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[5]*mnuOther-0.3162277660168379*m1rOther[5]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[3]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(76,40) = 0.1118033988749895*m0rOther[6]*uOther[9]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[8]*mnuOther+0.1118033988749895*uOther[6]*m0rOther[8]*mnuOther+0.125*m0rOther[0]*uOther[6]*mnuOther-0.3535533905932737*m1rOther[6]*mnuOther+0.125*uOther[0]*m0rOther[6]*mnuOther+0.125*m0rOther[4]*uOther[5]*mnuOther+0.125*uOther[4]*m0rOther[5]*mnuOther+0.125*m0rOther[2]*uOther[3]*mnuOther+0.125*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(76,41) = 0.125*m0rOther[1]*uOther[6]*mnuOther+0.125*uOther[1]*m0rOther[6]*mnuOther+0.125*m0rOther[2]*uOther[5]*mnuOther+0.125*uOther[2]*m0rOther[5]*mnuOther+0.125*m0rOther[3]*uOther[4]*mnuOther+0.125*uOther[3]*m0rOther[4]*mnuOther; 
  data->AEM_S(76,42) = 0.1118033988749895*m0rOther[3]*uOther[9]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[8]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[8]*mnuOther+0.225*m0rOther[2]*uOther[6]*mnuOther+0.225*uOther[2]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[5]*mnuOther+0.125*uOther[1]*m0rOther[5]*mnuOther+0.125*m0rOther[0]*uOther[3]*mnuOther-0.3535533905932737*m1rOther[3]*mnuOther+0.125*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(76,43) = 0.1118033988749895*m0rOther[2]*uOther[9]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[8]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[8]*mnuOther+0.225*m0rOther[3]*uOther[6]*mnuOther+0.225*uOther[3]*m0rOther[6]*mnuOther+0.125*m0rOther[1]*uOther[4]*mnuOther+0.125*uOther[1]*m0rOther[4]*mnuOther+0.125*m0rOther[0]*uOther[2]*mnuOther-0.3535533905932737*m1rOther[2]*mnuOther+0.125*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(76,44) = 0.1118033988749895*m0rOther[5]*uOther[9]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[8]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[7]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[7]*mnuOther+0.225*m0rOther[4]*uOther[6]*mnuOther+0.225*uOther[4]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[5]*mnuOther-0.3535533905932737*m1rOther[5]*mnuOther+0.125*uOther[0]*m0rOther[5]*mnuOther+0.125*m0rOther[1]*uOther[3]*mnuOther+0.125*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(76,45) = 0.1118033988749895*m0rOther[4]*uOther[9]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[8]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[7]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[7]*mnuOther+0.225*m0rOther[5]*uOther[6]*mnuOther+0.225*uOther[5]*m0rOther[6]*mnuOther+0.125*m0rOther[0]*uOther[4]*mnuOther-0.3535533905932737*m1rOther[4]*mnuOther+0.125*uOther[0]*m0rOther[4]*mnuOther+0.125*m0rOther[1]*uOther[2]*mnuOther+0.125*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(76,46) = 0.1964285714285714*m0rOther[9]*uOther[9]*mnuOther+0.1*m0rOther[8]*uOther[9]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[9]*mnuOther-0.3162277660168379*m1rOther[9]*mnuOther+0.1*uOther[8]*m0rOther[9]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[9]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[8]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[8]*mnuOther-0.3162277660168379*m1rOther[8]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[8]*mnuOther+0.125*m0rOther[7]*uOther[7]*mnuOther+0.405*m0rOther[6]*uOther[6]*mnuOther+0.225*m0rOther[5]*uOther[5]*mnuOther+0.225*m0rOther[4]*uOther[4]*mnuOther+0.225*m0rOther[3]*uOther[3]*mnuOther+0.225*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(76,47) = 0.125*m0rOther[6]*uOther[7]*mnuOther+0.125*uOther[6]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[5]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[5]*mnuOther; 
  data->AEM_S(76,48) = 0.1*m0rOther[6]*uOther[9]*mnuOther+0.1*uOther[6]*m0rOther[9]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[8]*mnuOther+0.1964285714285714*uOther[6]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[6]*mnuOther-0.3162277660168379*m1rOther[6]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[5]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[3]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(76,49) = 0.1964285714285714*m0rOther[6]*uOther[9]*mnuOther+0.1964285714285714*uOther[6]*m0rOther[9]*mnuOther+0.1*m0rOther[6]*uOther[8]*mnuOther+0.1*uOther[6]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[6]*mnuOther-0.3162277660168379*m1rOther[6]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[5]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[3]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(77,40) = 0.07985957062499249*m0rOther[7]*uOther[7]*mnuOther+0.125*m0rOther[0]*uOther[7]*mnuOther-0.3535533905932737*m1rOther[7]*mnuOther+0.125*uOther[0]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[5]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[4]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(77,41) = 0.1964285714285714*m0rOther[1]*uOther[7]*mnuOther+0.1964285714285714*uOther[1]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[5]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[4]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[1]*mnuOther-0.3162277660168379*m1rOther[1]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(77,42) = 0.125*m0rOther[2]*uOther[7]*mnuOther+0.125*uOther[2]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[4]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[4]*mnuOther; 
  data->AEM_S(77,43) = 0.125*m0rOther[3]*uOther[7]*mnuOther+0.125*uOther[3]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[5]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[5]*mnuOther; 
  data->AEM_S(77,44) = 0.1*m0rOther[4]*uOther[8]*mnuOther+0.1*uOther[4]*m0rOther[8]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[7]*mnuOther+0.1964285714285714*uOther[4]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[6]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[4]*mnuOther-0.3162277660168379*m1rOther[4]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[2]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(77,45) = 0.1*m0rOther[5]*uOther[9]*mnuOther+0.1*uOther[5]*m0rOther[9]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[7]*mnuOther+0.1964285714285714*uOther[5]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[6]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[5]*mnuOther-0.3162277660168379*m1rOther[5]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[3]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(77,46) = 0.125*m0rOther[6]*uOther[7]*mnuOther+0.125*uOther[6]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[5]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[5]*mnuOther; 
  data->AEM_S(77,47) = 0.125*m0rOther[9]*uOther[9]*mnuOther+0.125*m0rOther[8]*uOther[8]*mnuOther+0.2678571428571428*m0rOther[7]*uOther[7]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[7]*mnuOther-0.2258769757263128*m1rOther[7]*mnuOther+0.07985957062499249*uOther[0]*m0rOther[7]*mnuOther+0.125*m0rOther[6]*uOther[6]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[5]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[4]*mnuOther+0.125*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.1964285714285714*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(77,48) = 0.125*m0rOther[7]*uOther[8]*mnuOther+0.125*uOther[7]*m0rOther[8]*mnuOther+0.1*m0rOther[4]*uOther[4]*mnuOther; 
  data->AEM_S(77,49) = 0.125*m0rOther[7]*uOther[9]*mnuOther+0.125*uOther[7]*m0rOther[9]*mnuOther+0.1*m0rOther[5]*uOther[5]*mnuOther; 
  data->AEM_S(78,40) = 0.07985957062499249*m0rOther[8]*uOther[8]*mnuOther+0.125*m0rOther[0]*uOther[8]*mnuOther-0.3535533905932737*m1rOther[8]*mnuOther+0.125*uOther[0]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[6]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[4]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(78,41) = 0.125*m0rOther[1]*uOther[8]*mnuOther+0.125*uOther[1]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[4]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[4]*mnuOther; 
  data->AEM_S(78,42) = 0.1964285714285714*m0rOther[2]*uOther[8]*mnuOther+0.1964285714285714*uOther[2]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[6]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[4]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[2]*mnuOther-0.3162277660168379*m1rOther[2]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(78,43) = 0.125*m0rOther[3]*uOther[8]*mnuOther+0.125*uOther[3]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[6]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[6]*mnuOther; 
  data->AEM_S(78,44) = 0.1964285714285714*m0rOther[4]*uOther[8]*mnuOther+0.1964285714285714*uOther[4]*m0rOther[8]*mnuOther+0.1*m0rOther[4]*uOther[7]*mnuOther+0.1*uOther[4]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[6]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[4]*mnuOther-0.3162277660168379*m1rOther[4]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[4]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[2]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(78,45) = 0.125*m0rOther[5]*uOther[8]*mnuOther+0.125*uOther[5]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[6]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[6]*mnuOther; 
  data->AEM_S(78,46) = 0.1*m0rOther[6]*uOther[9]*mnuOther+0.1*uOther[6]*m0rOther[9]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[8]*mnuOther+0.1964285714285714*uOther[6]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[6]*mnuOther-0.3162277660168379*m1rOther[6]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[5]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[3]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(78,47) = 0.125*m0rOther[7]*uOther[8]*mnuOther+0.125*uOther[7]*m0rOther[8]*mnuOther+0.1*m0rOther[4]*uOther[4]*mnuOther; 
  data->AEM_S(78,48) = 0.125*m0rOther[9]*uOther[9]*mnuOther+0.2678571428571428*m0rOther[8]*uOther[8]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[8]*mnuOther-0.2258769757263128*m1rOther[8]*mnuOther+0.07985957062499249*uOther[0]*m0rOther[8]*mnuOther+0.125*m0rOther[7]*uOther[7]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[6]*mnuOther+0.125*m0rOther[5]*uOther[5]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[4]*mnuOther+0.125*m0rOther[3]*uOther[3]*mnuOther+0.1964285714285714*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
  data->AEM_S(78,49) = 0.125*m0rOther[8]*uOther[9]*mnuOther+0.125*uOther[8]*m0rOther[9]*mnuOther+0.1*m0rOther[6]*uOther[6]*mnuOther; 
  data->AEM_S(79,40) = 0.07985957062499249*m0rOther[9]*uOther[9]*mnuOther+0.125*m0rOther[0]*uOther[9]*mnuOther-0.3535533905932737*m1rOther[9]*mnuOther+0.125*uOther[0]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[6]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[5]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(79,41) = 0.125*m0rOther[1]*uOther[9]*mnuOther+0.125*uOther[1]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[5]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[5]*mnuOther; 
  data->AEM_S(79,42) = 0.125*m0rOther[2]*uOther[9]*mnuOther+0.125*uOther[2]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[6]*mnuOther+0.1118033988749895*uOther[3]*m0rOther[6]*mnuOther; 
  data->AEM_S(79,43) = 0.1964285714285714*m0rOther[3]*uOther[9]*mnuOther+0.1964285714285714*uOther[3]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[6]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[5]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[3]*mnuOther-0.3162277660168379*m1rOther[3]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[3]*mnuOther; 
  data->AEM_S(79,44) = 0.125*m0rOther[4]*uOther[9]*mnuOther+0.125*uOther[4]*m0rOther[9]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[6]*mnuOther+0.1118033988749895*uOther[5]*m0rOther[6]*mnuOther; 
  data->AEM_S(79,45) = 0.1964285714285714*m0rOther[5]*uOther[9]*mnuOther+0.1964285714285714*uOther[5]*m0rOther[9]*mnuOther+0.1*m0rOther[5]*uOther[7]*mnuOther+0.1*uOther[5]*m0rOther[7]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[6]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[5]*mnuOther-0.3162277660168379*m1rOther[5]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[3]*mnuOther+0.1118033988749895*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(79,46) = 0.1964285714285714*m0rOther[6]*uOther[9]*mnuOther+0.1964285714285714*uOther[6]*m0rOther[9]*mnuOther+0.1*m0rOther[6]*uOther[8]*mnuOther+0.1*uOther[6]*m0rOther[8]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[6]*mnuOther-0.3162277660168379*m1rOther[6]*mnuOther+0.1118033988749895*uOther[0]*m0rOther[6]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[5]*mnuOther+0.1118033988749895*uOther[4]*m0rOther[5]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[3]*mnuOther+0.1118033988749895*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(79,47) = 0.125*m0rOther[7]*uOther[9]*mnuOther+0.125*uOther[7]*m0rOther[9]*mnuOther+0.1*m0rOther[5]*uOther[5]*mnuOther; 
  data->AEM_S(79,48) = 0.125*m0rOther[8]*uOther[9]*mnuOther+0.125*uOther[8]*m0rOther[9]*mnuOther+0.1*m0rOther[6]*uOther[6]*mnuOther; 
  data->AEM_S(79,49) = 0.2678571428571428*m0rOther[9]*uOther[9]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[9]*mnuOther-0.2258769757263128*m1rOther[9]*mnuOther+0.07985957062499249*uOther[0]*m0rOther[9]*mnuOther+0.125*m0rOther[8]*uOther[8]*mnuOther+0.125*m0rOther[7]*uOther[7]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[6]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[5]*mnuOther+0.125*m0rOther[4]*uOther[4]*mnuOther+0.1964285714285714*m0rOther[3]*uOther[3]*mnuOther+0.125*m0rOther[2]*uOther[2]*mnuOther+0.125*m0rOther[1]*uOther[1]*mnuOther+0.125*m0rOther[0]*uOther[0]*mnuOther-0.3535533905932737*m1rOther[0]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfX-m0Self*m1OtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[0] = 0.3535533905932737*m0rOther[9]*m1rSelf[9]-0.3535533905932737*m0rSelf[9]*m1rOther[9]+0.3535533905932737*m0rOther[8]*m1rSelf[8]-0.3535533905932737*m0rSelf[8]*m1rOther[8]+0.3535533905932737*m0rOther[7]*m1rSelf[7]-0.3535533905932737*m0rSelf[7]*m1rOther[7]+0.3535533905932737*m0rOther[6]*m1rSelf[6]-0.3535533905932737*m0rSelf[6]*m1rOther[6]+0.3535533905932737*m0rOther[5]*m1rSelf[5]-0.3535533905932737*m0rSelf[5]*m1rOther[5]+0.3535533905932737*m0rOther[4]*m1rSelf[4]-0.3535533905932737*m0rSelf[4]*m1rOther[4]+0.3535533905932737*m0rOther[3]*m1rSelf[3]-0.3535533905932737*m0rSelf[3]*m1rOther[3]+0.3535533905932737*m0rOther[2]*m1rSelf[2]-0.3535533905932737*m0rSelf[2]*m1rOther[2]+0.3535533905932737*m0rOther[1]*m1rSelf[1]-0.3535533905932737*m0rSelf[1]*m1rOther[1]+0.3535533905932737*m0rOther[0]*m1rSelf[0]-0.3535533905932737*m0rSelf[0]*m1rOther[0]; 
  m1EffD[1] = 0.3162277660168379*m0rOther[1]*m1rSelf[7]-0.3162277660168379*m0rSelf[1]*m1rOther[7]-0.3162277660168379*m1rOther[1]*m0rSelf[7]+0.3162277660168379*m1rSelf[1]*m0rOther[7]+0.3535533905932737*m0rOther[3]*m1rSelf[5]-0.3535533905932737*m0rSelf[3]*m1rOther[5]-0.3535533905932737*m1rOther[3]*m0rSelf[5]+0.3535533905932737*m1rSelf[3]*m0rOther[5]+0.3535533905932737*m0rOther[2]*m1rSelf[4]-0.3535533905932737*m0rSelf[2]*m1rOther[4]-0.3535533905932737*m1rOther[2]*m0rSelf[4]+0.3535533905932737*m1rSelf[2]*m0rOther[4]+0.3535533905932737*m0rOther[0]*m1rSelf[1]-0.3535533905932737*m0rSelf[0]*m1rOther[1]-0.3535533905932737*m1rOther[0]*m0rSelf[1]+0.3535533905932737*m1rSelf[0]*m0rOther[1]; 
  m1EffD[2] = 0.3162277660168379*m0rOther[2]*m1rSelf[8]-0.3162277660168379*m0rSelf[2]*m1rOther[8]-0.3162277660168379*m1rOther[2]*m0rSelf[8]+0.3162277660168379*m1rSelf[2]*m0rOther[8]+0.3535533905932737*m0rOther[3]*m1rSelf[6]-0.3535533905932737*m0rSelf[3]*m1rOther[6]-0.3535533905932737*m1rOther[3]*m0rSelf[6]+0.3535533905932737*m1rSelf[3]*m0rOther[6]+0.3535533905932737*m0rOther[1]*m1rSelf[4]-0.3535533905932737*m0rSelf[1]*m1rOther[4]-0.3535533905932737*m1rOther[1]*m0rSelf[4]+0.3535533905932737*m1rSelf[1]*m0rOther[4]+0.3535533905932737*m0rOther[0]*m1rSelf[2]-0.3535533905932737*m0rSelf[0]*m1rOther[2]-0.3535533905932737*m1rOther[0]*m0rSelf[2]+0.3535533905932737*m1rSelf[0]*m0rOther[2]; 
  m1EffD[3] = 0.3162277660168379*m0rOther[3]*m1rSelf[9]-0.3162277660168379*m0rSelf[3]*m1rOther[9]-0.3162277660168379*m1rOther[3]*m0rSelf[9]+0.3162277660168379*m1rSelf[3]*m0rOther[9]+0.3535533905932737*m0rOther[2]*m1rSelf[6]-0.3535533905932737*m0rSelf[2]*m1rOther[6]-0.3535533905932737*m1rOther[2]*m0rSelf[6]+0.3535533905932737*m1rSelf[2]*m0rOther[6]+0.3535533905932737*m0rOther[1]*m1rSelf[5]-0.3535533905932737*m0rSelf[1]*m1rOther[5]-0.3535533905932737*m1rOther[1]*m0rSelf[5]+0.3535533905932737*m1rSelf[1]*m0rOther[5]+0.3535533905932737*m0rOther[0]*m1rSelf[3]-0.3535533905932737*m0rSelf[0]*m1rOther[3]-0.3535533905932737*m1rOther[0]*m0rSelf[3]+0.3535533905932737*m1rSelf[0]*m0rOther[3]; 
  m1EffD[4] = 0.3162277660168379*m0rOther[4]*m1rSelf[8]-0.3162277660168379*m0rSelf[4]*m1rOther[8]-0.3162277660168379*m1rOther[4]*m0rSelf[8]+0.3162277660168379*m1rSelf[4]*m0rOther[8]+0.3162277660168379*m0rOther[4]*m1rSelf[7]-0.3162277660168379*m0rSelf[4]*m1rOther[7]-0.3162277660168379*m1rOther[4]*m0rSelf[7]+0.3162277660168379*m1rSelf[4]*m0rOther[7]+0.3535533905932737*m0rOther[5]*m1rSelf[6]-0.3535533905932737*m0rSelf[5]*m1rOther[6]-0.3535533905932737*m1rOther[5]*m0rSelf[6]+0.3535533905932737*m1rSelf[5]*m0rOther[6]+0.3535533905932737*m0rOther[0]*m1rSelf[4]-0.3535533905932737*m0rSelf[0]*m1rOther[4]-0.3535533905932737*m1rOther[0]*m0rSelf[4]+0.3535533905932737*m1rSelf[0]*m0rOther[4]+0.3535533905932737*m0rOther[1]*m1rSelf[2]-0.3535533905932737*m0rSelf[1]*m1rOther[2]-0.3535533905932737*m1rOther[1]*m0rSelf[2]+0.3535533905932737*m1rSelf[1]*m0rOther[2]; 
  m1EffD[5] = 0.3162277660168379*m0rOther[5]*m1rSelf[9]-0.3162277660168379*m0rSelf[5]*m1rOther[9]-0.3162277660168379*m1rOther[5]*m0rSelf[9]+0.3162277660168379*m1rSelf[5]*m0rOther[9]+0.3162277660168379*m0rOther[5]*m1rSelf[7]-0.3162277660168379*m0rSelf[5]*m1rOther[7]-0.3162277660168379*m1rOther[5]*m0rSelf[7]+0.3162277660168379*m1rSelf[5]*m0rOther[7]+0.3535533905932737*m0rOther[4]*m1rSelf[6]-0.3535533905932737*m0rSelf[4]*m1rOther[6]-0.3535533905932737*m1rOther[4]*m0rSelf[6]+0.3535533905932737*m1rSelf[4]*m0rOther[6]+0.3535533905932737*m0rOther[0]*m1rSelf[5]-0.3535533905932737*m0rSelf[0]*m1rOther[5]-0.3535533905932737*m1rOther[0]*m0rSelf[5]+0.3535533905932737*m1rSelf[0]*m0rOther[5]+0.3535533905932737*m0rOther[1]*m1rSelf[3]-0.3535533905932737*m0rSelf[1]*m1rOther[3]-0.3535533905932737*m1rOther[1]*m0rSelf[3]+0.3535533905932737*m1rSelf[1]*m0rOther[3]; 
  m1EffD[6] = 0.3162277660168379*m0rOther[6]*m1rSelf[9]-0.3162277660168379*m0rSelf[6]*m1rOther[9]-0.3162277660168379*m1rOther[6]*m0rSelf[9]+0.3162277660168379*m1rSelf[6]*m0rOther[9]+0.3162277660168379*m0rOther[6]*m1rSelf[8]-0.3162277660168379*m0rSelf[6]*m1rOther[8]-0.3162277660168379*m1rOther[6]*m0rSelf[8]+0.3162277660168379*m1rSelf[6]*m0rOther[8]+0.3535533905932737*m0rOther[0]*m1rSelf[6]-0.3535533905932737*m0rSelf[0]*m1rOther[6]-0.3535533905932737*m1rOther[0]*m0rSelf[6]+0.3535533905932737*m1rSelf[0]*m0rOther[6]+0.3535533905932737*m0rOther[4]*m1rSelf[5]-0.3535533905932737*m0rSelf[4]*m1rOther[5]-0.3535533905932737*m1rOther[4]*m0rSelf[5]+0.3535533905932737*m1rSelf[4]*m0rOther[5]+0.3535533905932737*m0rOther[2]*m1rSelf[3]-0.3535533905932737*m0rSelf[2]*m1rOther[3]-0.3535533905932737*m1rOther[2]*m0rSelf[3]+0.3535533905932737*m1rSelf[2]*m0rOther[3]; 
  m1EffD[7] = 0.2258769757263128*m0rOther[7]*m1rSelf[7]+0.3535533905932737*m0rOther[0]*m1rSelf[7]-0.2258769757263128*m0rSelf[7]*m1rOther[7]-0.3535533905932737*m0rSelf[0]*m1rOther[7]-0.3535533905932737*m1rOther[0]*m0rSelf[7]+0.3535533905932737*m1rSelf[0]*m0rOther[7]+0.3162277660168379*m0rOther[5]*m1rSelf[5]-0.3162277660168379*m0rSelf[5]*m1rOther[5]+0.3162277660168379*m0rOther[4]*m1rSelf[4]-0.3162277660168379*m0rSelf[4]*m1rOther[4]+0.3162277660168379*m0rOther[1]*m1rSelf[1]-0.3162277660168379*m0rSelf[1]*m1rOther[1]; 
  m1EffD[8] = 0.2258769757263128*m0rOther[8]*m1rSelf[8]+0.3535533905932737*m0rOther[0]*m1rSelf[8]-0.2258769757263128*m0rSelf[8]*m1rOther[8]-0.3535533905932737*m0rSelf[0]*m1rOther[8]-0.3535533905932737*m1rOther[0]*m0rSelf[8]+0.3535533905932737*m1rSelf[0]*m0rOther[8]+0.3162277660168379*m0rOther[6]*m1rSelf[6]-0.3162277660168379*m0rSelf[6]*m1rOther[6]+0.3162277660168379*m0rOther[4]*m1rSelf[4]-0.3162277660168379*m0rSelf[4]*m1rOther[4]+0.3162277660168379*m0rOther[2]*m1rSelf[2]-0.3162277660168379*m0rSelf[2]*m1rOther[2]; 
  m1EffD[9] = 0.2258769757263128*m0rOther[9]*m1rSelf[9]+0.3535533905932737*m0rOther[0]*m1rSelf[9]-0.2258769757263128*m0rSelf[9]*m1rOther[9]-0.3535533905932737*m0rSelf[0]*m1rOther[9]-0.3535533905932737*m1rOther[0]*m0rSelf[9]+0.3535533905932737*m1rSelf[0]*m0rOther[9]+0.3162277660168379*m0rOther[6]*m1rSelf[6]-0.3162277660168379*m0rSelf[6]*m1rOther[6]+0.3162277660168379*m0rOther[5]*m1rSelf[5]-0.3162277660168379*m0rSelf[5]*m1rOther[5]+0.3162277660168379*m0rOther[3]*m1rSelf[3]-0.3162277660168379*m0rSelf[3]*m1rOther[3]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(10,10); 
  dataDiv->AEM_S(0,0) = 0.3535533905932737*m0rSelf[0]*mnuSelf+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(0,2) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(0,3) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(0,4) = 0.3535533905932737*m0rSelf[4]*mnuSelf+0.3535533905932737*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(0,5) = 0.3535533905932737*m0rSelf[5]*mnuSelf+0.3535533905932737*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(0,6) = 0.3535533905932737*m0rSelf[6]*mnuSelf+0.3535533905932737*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(0,7) = 0.3535533905932737*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(0,8) = 0.3535533905932737*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rOther[8]*mnuOther; 
  dataDiv->AEM_S(0,9) = 0.3535533905932737*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rOther[9]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(1,2) = 0.3535533905932737*m0rSelf[4]*mnuSelf+0.3535533905932737*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(1,3) = 0.3535533905932737*m0rSelf[5]*mnuSelf+0.3535533905932737*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(1,4) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(1,5) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(1,7) = 0.3162277660168379*m0rSelf[1]*mnuSelf+0.3162277660168379*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(2,0) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,1) = 0.3535533905932737*m0rSelf[4]*mnuSelf+0.3535533905932737*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(2,2) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(2,3) = 0.3535533905932737*m0rSelf[6]*mnuSelf+0.3535533905932737*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(2,4) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(2,6) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(2,8) = 0.3162277660168379*m0rSelf[2]*mnuSelf+0.3162277660168379*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(3,0) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,1) = 0.3535533905932737*m0rSelf[5]*mnuSelf+0.3535533905932737*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(3,2) = 0.3535533905932737*m0rSelf[6]*mnuSelf+0.3535533905932737*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(3,3) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.3162277660168379*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(3,5) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(3,6) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(3,9) = 0.3162277660168379*m0rSelf[3]*mnuSelf+0.3162277660168379*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(4,0) = 0.3535533905932737*m0rSelf[4]*mnuSelf+0.3535533905932737*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(4,1) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(4,2) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(4,4) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.3162277660168379*m0rOther[8]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(4,5) = 0.3535533905932737*m0rSelf[6]*mnuSelf+0.3535533905932737*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(4,6) = 0.3535533905932737*m0rSelf[5]*mnuSelf+0.3535533905932737*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(4,7) = 0.3162277660168379*m0rSelf[4]*mnuSelf+0.3162277660168379*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(4,8) = 0.3162277660168379*m0rSelf[4]*mnuSelf+0.3162277660168379*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(5,0) = 0.3535533905932737*m0rSelf[5]*mnuSelf+0.3535533905932737*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(5,1) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(5,3) = 0.3535533905932737*m0rSelf[1]*mnuSelf+0.3535533905932737*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(5,4) = 0.3535533905932737*m0rSelf[6]*mnuSelf+0.3535533905932737*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(5,5) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(5,6) = 0.3535533905932737*m0rSelf[4]*mnuSelf+0.3535533905932737*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(5,7) = 0.3162277660168379*m0rSelf[5]*mnuSelf+0.3162277660168379*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(5,9) = 0.3162277660168379*m0rSelf[5]*mnuSelf+0.3162277660168379*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(6,0) = 0.3535533905932737*m0rSelf[6]*mnuSelf+0.3535533905932737*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(6,2) = 0.3535533905932737*m0rSelf[3]*mnuSelf+0.3535533905932737*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(6,3) = 0.3535533905932737*m0rSelf[2]*mnuSelf+0.3535533905932737*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(6,4) = 0.3535533905932737*m0rSelf[5]*mnuSelf+0.3535533905932737*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(6,5) = 0.3535533905932737*m0rSelf[4]*mnuSelf+0.3535533905932737*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(6,6) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.3162277660168379*m0rOther[9]*mnuOther+0.3162277660168379*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(6,8) = 0.3162277660168379*m0rSelf[6]*mnuSelf+0.3162277660168379*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(6,9) = 0.3162277660168379*m0rSelf[6]*mnuSelf+0.3162277660168379*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(7,0) = 0.3535533905932737*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(7,1) = 0.3162277660168379*m0rSelf[1]*mnuSelf+0.3162277660168379*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(7,4) = 0.3162277660168379*m0rSelf[4]*mnuSelf+0.3162277660168379*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(7,5) = 0.3162277660168379*m0rSelf[5]*mnuSelf+0.3162277660168379*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(7,7) = 0.2258769757263128*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.2258769757263128*m0rOther[7]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(8,0) = 0.3535533905932737*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rOther[8]*mnuOther; 
  dataDiv->AEM_S(8,2) = 0.3162277660168379*m0rSelf[2]*mnuSelf+0.3162277660168379*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(8,4) = 0.3162277660168379*m0rSelf[4]*mnuSelf+0.3162277660168379*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(8,6) = 0.3162277660168379*m0rSelf[6]*mnuSelf+0.3162277660168379*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(8,8) = 0.2258769757263128*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.2258769757263128*m0rOther[8]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(9,0) = 0.3535533905932737*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rOther[9]*mnuOther; 
  dataDiv->AEM_S(9,3) = 0.3162277660168379*m0rSelf[3]*mnuSelf+0.3162277660168379*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(9,5) = 0.3162277660168379*m0rSelf[5]*mnuSelf+0.3162277660168379*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(9,6) = 0.3162277660168379*m0rSelf[6]*mnuSelf+0.3162277660168379*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(9,9) = 0.2258769757263128*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf+0.2258769757263128*m0rOther[9]*mnuOther+0.3535533905932737*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[0],m1EffD[1],m1EffD[2],m1EffD[3],m1EffD[4],m1EffD[5],m1EffD[6],m1EffD[7],m1EffD[8],m1EffD[9]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+0,10,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (-2.0*m1EffD[0]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (-2.0*m1EffD[1]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (-2.0*m1EffD[2]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += (-2.0*m1EffD[3]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (-2.0*m1EffD[4]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (-2.0*m1EffD[5]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
  m1Relax[6] += (-2.0*m1EffD[6]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (-2.0*m1EffD[7]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
  m1Relax[8] += (-2.0*m1EffD[8]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
  m1Relax[9] += (-2.0*m1EffD[9]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[9]*mnuSelf-1.0*m1rOther[9]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(50,10) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(50,11) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(50,12) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(50,13) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(50,14) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(50,15) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(50,16) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(50,17) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(50,18) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(50,19) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(51,10) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(51,11) = 0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(51,12) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(51,13) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(51,14) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(51,15) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(51,17) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(52,10) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(52,11) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(52,12) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(52,13) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(52,14) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(52,16) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(52,18) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,10) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(53,11) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(53,12) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(53,13) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(53,15) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(53,16) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,19) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(54,10) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(54,11) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(54,12) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(54,14) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(54,15) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(54,16) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(54,17) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(54,18) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(55,10) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(55,11) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(55,13) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(55,14) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(55,15) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(55,16) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(55,17) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(55,19) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(56,10) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(56,12) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(56,13) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,14) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(56,15) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(56,16) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(56,18) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(56,19) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(57,10) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(57,11) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(57,14) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(57,15) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(57,17) = 0.2258769757263128*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(58,10) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(58,12) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(58,14) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(58,16) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(58,18) = 0.2258769757263128*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(59,10) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(59,13) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(59,15) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(59,16) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(59,19) = 0.2258769757263128*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(50,30) = -0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(50,31) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(50,32) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(50,33) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(50,34) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(50,35) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(50,36) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(50,37) = -0.3535533905932737*cMSelf[17]*mnuSelf; 
  data->AEM_S(50,38) = -0.3535533905932737*cMSelf[18]*mnuSelf; 
  data->AEM_S(50,39) = -0.3535533905932737*cMSelf[19]*mnuSelf; 
  data->AEM_S(51,30) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(51,31) = (-0.3162277660168379*cMSelf[17]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(51,32) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(51,33) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(51,34) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(51,35) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(51,37) = -0.3162277660168379*cMSelf[11]*mnuSelf; 
  data->AEM_S(52,30) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(52,31) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(52,32) = (-0.3162277660168379*cMSelf[18]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(52,33) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(52,34) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(52,36) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(52,38) = -0.3162277660168379*cMSelf[12]*mnuSelf; 
  data->AEM_S(53,30) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(53,31) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(53,32) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(53,33) = (-0.3162277660168379*cMSelf[19]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(53,35) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(53,36) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(53,39) = -0.3162277660168379*cMSelf[13]*mnuSelf; 
  data->AEM_S(54,30) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(54,31) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(54,32) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(54,34) = (-0.3162277660168379*cMSelf[18]*mnuSelf)-0.3162277660168379*cMSelf[17]*mnuSelf-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(54,35) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(54,36) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(54,37) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(54,38) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(55,30) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(55,31) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(55,33) = -0.3535533905932737*cMSelf[11]*mnuSelf; 
  data->AEM_S(55,34) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(55,35) = (-0.3162277660168379*cMSelf[19]*mnuSelf)-0.3162277660168379*cMSelf[17]*mnuSelf-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(55,36) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(55,37) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(55,39) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(56,30) = -0.3535533905932737*cMSelf[16]*mnuSelf; 
  data->AEM_S(56,32) = -0.3535533905932737*cMSelf[13]*mnuSelf; 
  data->AEM_S(56,33) = -0.3535533905932737*cMSelf[12]*mnuSelf; 
  data->AEM_S(56,34) = -0.3535533905932737*cMSelf[15]*mnuSelf; 
  data->AEM_S(56,35) = -0.3535533905932737*cMSelf[14]*mnuSelf; 
  data->AEM_S(56,36) = (-0.3162277660168379*cMSelf[19]*mnuSelf)-0.3162277660168379*cMSelf[18]*mnuSelf-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(56,38) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(56,39) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(57,30) = -0.3535533905932737*cMSelf[17]*mnuSelf; 
  data->AEM_S(57,31) = -0.3162277660168379*cMSelf[11]*mnuSelf; 
  data->AEM_S(57,34) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(57,35) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(57,37) = (-0.2258769757263128*cMSelf[17]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(58,30) = -0.3535533905932737*cMSelf[18]*mnuSelf; 
  data->AEM_S(58,32) = -0.3162277660168379*cMSelf[12]*mnuSelf; 
  data->AEM_S(58,34) = -0.3162277660168379*cMSelf[14]*mnuSelf; 
  data->AEM_S(58,36) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(58,38) = (-0.2258769757263128*cMSelf[18]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
  data->AEM_S(59,30) = -0.3535533905932737*cMSelf[19]*mnuSelf; 
  data->AEM_S(59,33) = -0.3162277660168379*cMSelf[13]*mnuSelf; 
  data->AEM_S(59,35) = -0.3162277660168379*cMSelf[15]*mnuSelf; 
  data->AEM_S(59,36) = -0.3162277660168379*cMSelf[16]*mnuSelf; 
  data->AEM_S(59,39) = (-0.2258769757263128*cMSelf[19]*mnuSelf)-0.3535533905932737*cMSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(50,50) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(50,51) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(50,52) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(50,53) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(50,54) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(50,55) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(50,56) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(50,57) = -0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(50,58) = -0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(50,59) = -0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(51,50) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(51,51) = (-0.3162277660168379*m0rOther[7]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(51,52) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(51,53) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(51,54) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(51,55) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(51,57) = -0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(52,50) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(52,51) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(52,52) = (-0.3162277660168379*m0rOther[8]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(52,53) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(52,54) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(52,56) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(52,58) = -0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(53,50) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(53,51) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(53,52) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(53,53) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(53,55) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(53,56) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(53,59) = -0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(54,50) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(54,51) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(54,52) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(54,54) = (-0.3162277660168379*m0rOther[8]*mnuOther)-0.3162277660168379*m0rOther[7]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(54,55) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(54,56) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(54,57) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(54,58) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(55,50) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(55,51) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(55,53) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(55,54) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(55,55) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3162277660168379*m0rOther[7]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(55,56) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(55,57) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(55,59) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(56,50) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(56,52) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(56,53) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(56,54) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(56,55) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(56,56) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3162277660168379*m0rOther[8]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(56,58) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(56,59) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(57,50) = -0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(57,51) = -0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(57,54) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(57,55) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(57,57) = (-0.2258769757263128*m0rOther[7]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(58,50) = -0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(58,52) = -0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(58,54) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(58,56) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(58,58) = (-0.2258769757263128*m0rOther[8]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(59,50) = -0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(59,53) = -0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(59,55) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(59,56) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(59,59) = (-0.2258769757263128*m0rOther[9]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(50,70) = 0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(50,71) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(50,72) = 0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(50,73) = 0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(50,74) = 0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(50,75) = 0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(50,76) = 0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(50,77) = 0.3535533905932737*cMOther[17]*mnuOther; 
  data->AEM_S(50,78) = 0.3535533905932737*cMOther[18]*mnuOther; 
  data->AEM_S(50,79) = 0.3535533905932737*cMOther[19]*mnuOther; 
  data->AEM_S(51,70) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(51,71) = 0.3162277660168379*cMOther[17]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(51,72) = 0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(51,73) = 0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(51,74) = 0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(51,75) = 0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(51,77) = 0.3162277660168379*cMOther[11]*mnuOther; 
  data->AEM_S(52,70) = 0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(52,71) = 0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(52,72) = 0.3162277660168379*cMOther[18]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(52,73) = 0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(52,74) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(52,76) = 0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(52,78) = 0.3162277660168379*cMOther[12]*mnuOther; 
  data->AEM_S(53,70) = 0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(53,71) = 0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(53,72) = 0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(53,73) = 0.3162277660168379*cMOther[19]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(53,75) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(53,76) = 0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(53,79) = 0.3162277660168379*cMOther[13]*mnuOther; 
  data->AEM_S(54,70) = 0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(54,71) = 0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(54,72) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(54,74) = 0.3162277660168379*cMOther[18]*mnuOther+0.3162277660168379*cMOther[17]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(54,75) = 0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(54,76) = 0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(54,77) = 0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(54,78) = 0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(55,70) = 0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(55,71) = 0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(55,73) = 0.3535533905932737*cMOther[11]*mnuOther; 
  data->AEM_S(55,74) = 0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(55,75) = 0.3162277660168379*cMOther[19]*mnuOther+0.3162277660168379*cMOther[17]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(55,76) = 0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(55,77) = 0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(55,79) = 0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(56,70) = 0.3535533905932737*cMOther[16]*mnuOther; 
  data->AEM_S(56,72) = 0.3535533905932737*cMOther[13]*mnuOther; 
  data->AEM_S(56,73) = 0.3535533905932737*cMOther[12]*mnuOther; 
  data->AEM_S(56,74) = 0.3535533905932737*cMOther[15]*mnuOther; 
  data->AEM_S(56,75) = 0.3535533905932737*cMOther[14]*mnuOther; 
  data->AEM_S(56,76) = 0.3162277660168379*cMOther[19]*mnuOther+0.3162277660168379*cMOther[18]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(56,78) = 0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(56,79) = 0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(57,70) = 0.3535533905932737*cMOther[17]*mnuOther; 
  data->AEM_S(57,71) = 0.3162277660168379*cMOther[11]*mnuOther; 
  data->AEM_S(57,74) = 0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(57,75) = 0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(57,77) = 0.2258769757263128*cMOther[17]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(58,70) = 0.3535533905932737*cMOther[18]*mnuOther; 
  data->AEM_S(58,72) = 0.3162277660168379*cMOther[12]*mnuOther; 
  data->AEM_S(58,74) = 0.3162277660168379*cMOther[14]*mnuOther; 
  data->AEM_S(58,76) = 0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(58,78) = 0.2258769757263128*cMOther[18]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
  data->AEM_S(59,70) = 0.3535533905932737*cMOther[19]*mnuOther; 
  data->AEM_S(59,73) = 0.3162277660168379*cMOther[13]*mnuOther; 
  data->AEM_S(59,75) = 0.3162277660168379*cMOther[15]*mnuOther; 
  data->AEM_S(59,76) = 0.3162277660168379*cMOther[16]*mnuOther; 
  data->AEM_S(59,79) = 0.2258769757263128*cMOther[19]*mnuOther+0.3535533905932737*cMOther[10]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(70,10) = (-0.125*m0rSelf[9]*uSelf[19]*mnuSelf)-0.125*m0rSelf[8]*uSelf[18]*mnuSelf-0.125*m0rSelf[7]*uSelf[17]*mnuSelf-0.125*m0rSelf[6]*uSelf[16]*mnuSelf-0.125*m0rSelf[5]*uSelf[15]*mnuSelf-0.125*m0rSelf[4]*uSelf[14]*mnuSelf-0.125*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[1]*uSelf[11]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(70,11) = (-0.1118033988749895*m0rSelf[1]*uSelf[17]*mnuSelf)-0.125*m0rSelf[3]*uSelf[15]*mnuSelf-0.125*m0rSelf[2]*uSelf[14]*mnuSelf-0.125*m0rSelf[5]*uSelf[13]*mnuSelf-0.125*m0rSelf[4]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[11]*mnuSelf-0.125*m0rSelf[0]*uSelf[11]*mnuSelf+0.3535533905932737*m1rSelf[11]*mnuSelf-0.125*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,12) = (-0.1118033988749895*m0rSelf[2]*uSelf[18]*mnuSelf)-0.125*m0rSelf[3]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[14]*mnuSelf-0.125*m0rSelf[6]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[12]*mnuSelf-0.125*m0rSelf[0]*uSelf[12]*mnuSelf+0.3535533905932737*m1rSelf[12]*mnuSelf-0.125*m0rSelf[4]*uSelf[11]*mnuSelf-0.125*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,13) = (-0.1118033988749895*m0rSelf[3]*uSelf[19]*mnuSelf)-0.125*m0rSelf[2]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[13]*mnuSelf-0.125*m0rSelf[0]*uSelf[13]*mnuSelf+0.3535533905932737*m1rSelf[13]*mnuSelf-0.125*m0rSelf[6]*uSelf[12]*mnuSelf-0.125*m0rSelf[5]*uSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,14) = (-0.1118033988749895*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[17]*mnuSelf-0.125*m0rSelf[5]*uSelf[16]*mnuSelf-0.125*m0rSelf[6]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[14]*mnuSelf-0.125*m0rSelf[0]*uSelf[14]*mnuSelf+0.3535533905932737*m1rSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[12]*mnuSelf-0.125*m0rSelf[2]*uSelf[11]*mnuSelf-0.125*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,15) = (-0.1118033988749895*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[17]*mnuSelf-0.125*m0rSelf[4]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[15]*mnuSelf-0.125*m0rSelf[0]*uSelf[15]*mnuSelf+0.3535533905932737*m1rSelf[15]*mnuSelf-0.125*m0rSelf[6]*uSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[11]*mnuSelf-0.125*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,16) = (-0.1118033988749895*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[16]*mnuSelf-0.125*m0rSelf[0]*uSelf[16]*mnuSelf+0.3535533905932737*m1rSelf[16]*mnuSelf-0.125*m0rSelf[4]*uSelf[15]*mnuSelf-0.125*m0rSelf[5]*uSelf[14]*mnuSelf-0.125*m0rSelf[2]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,17) = (-0.07985957062499249*m0rSelf[7]*uSelf[17]*mnuSelf)-0.125*m0rSelf[0]*uSelf[17]*mnuSelf+0.3535533905932737*m1rSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[11]*mnuSelf-0.125*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,18) = (-0.07985957062499249*m0rSelf[8]*uSelf[18]*mnuSelf)-0.125*m0rSelf[0]*uSelf[18]*mnuSelf+0.3535533905932737*m1rSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[8]*uSelf[10]*mnuSelf; 
  data->AEM_S(70,19) = (-0.07985957062499249*m0rSelf[9]*uSelf[19]*mnuSelf)-0.125*m0rSelf[0]*uSelf[19]*mnuSelf+0.3535533905932737*m1rSelf[19]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[9]*uSelf[10]*mnuSelf; 
  data->AEM_S(71,10) = (-0.1118033988749895*m0rSelf[1]*uSelf[17]*mnuSelf)-0.125*m0rSelf[3]*uSelf[15]*mnuSelf-0.125*m0rSelf[2]*uSelf[14]*mnuSelf-0.125*m0rSelf[5]*uSelf[13]*mnuSelf-0.125*m0rSelf[4]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[11]*mnuSelf-0.125*m0rSelf[0]*uSelf[11]*mnuSelf+0.3535533905932737*m1rSelf[11]*mnuSelf-0.125*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(71,11) = (-0.125*m0rSelf[9]*uSelf[19]*mnuSelf)-0.125*m0rSelf[8]*uSelf[18]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[17]*mnuSelf+0.3162277660168379*m1rSelf[17]*mnuSelf-0.125*m0rSelf[6]*uSelf[16]*mnuSelf-0.225*m0rSelf[5]*uSelf[15]*mnuSelf-0.225*m0rSelf[4]*uSelf[14]*mnuSelf-0.125*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[2]*uSelf[12]*mnuSelf-0.225*m0rSelf[1]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(71,12) = (-0.1118033988749895*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[17]*mnuSelf-0.125*m0rSelf[5]*uSelf[16]*mnuSelf-0.125*m0rSelf[6]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[14]*mnuSelf-0.125*m0rSelf[0]*uSelf[14]*mnuSelf+0.3535533905932737*m1rSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[12]*mnuSelf-0.125*m0rSelf[2]*uSelf[11]*mnuSelf-0.125*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(71,13) = (-0.1118033988749895*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[17]*mnuSelf-0.125*m0rSelf[4]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[15]*mnuSelf-0.125*m0rSelf[0]*uSelf[15]*mnuSelf+0.3535533905932737*m1rSelf[15]*mnuSelf-0.125*m0rSelf[6]*uSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[11]*mnuSelf-0.125*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(71,14) = (-0.1118033988749895*m0rSelf[2]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[17]*mnuSelf-0.125*m0rSelf[3]*uSelf[16]*mnuSelf-0.225*m0rSelf[1]*uSelf[14]*mnuSelf-0.125*m0rSelf[6]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[12]*mnuSelf-0.125*m0rSelf[0]*uSelf[12]*mnuSelf+0.3535533905932737*m1rSelf[12]*mnuSelf-0.225*m0rSelf[4]*uSelf[11]*mnuSelf-0.125*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(71,15) = (-0.1118033988749895*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[17]*mnuSelf-0.125*m0rSelf[2]*uSelf[16]*mnuSelf-0.225*m0rSelf[1]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[13]*mnuSelf-0.125*m0rSelf[0]*uSelf[13]*mnuSelf+0.3535533905932737*m1rSelf[13]*mnuSelf-0.125*m0rSelf[6]*uSelf[12]*mnuSelf-0.225*m0rSelf[5]*uSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(71,16) = (-0.125*m0rSelf[1]*uSelf[16]*mnuSelf)-0.125*m0rSelf[2]*uSelf[15]*mnuSelf-0.125*m0rSelf[3]*uSelf[14]*mnuSelf-0.125*m0rSelf[4]*uSelf[13]*mnuSelf-0.125*m0rSelf[5]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(71,17) = (-0.1964285714285714*m0rSelf[1]*uSelf[17]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[12]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[11]*mnuSelf+0.3162277660168379*m1rSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(71,18) = (-0.125*m0rSelf[1]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[12]*mnuSelf-0.125*m0rSelf[8]*uSelf[11]*mnuSelf; 
  data->AEM_S(71,19) = (-0.125*m0rSelf[1]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[13]*mnuSelf-0.125*m0rSelf[9]*uSelf[11]*mnuSelf; 
  data->AEM_S(72,10) = (-0.1118033988749895*m0rSelf[2]*uSelf[18]*mnuSelf)-0.125*m0rSelf[3]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[14]*mnuSelf-0.125*m0rSelf[6]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[12]*mnuSelf-0.125*m0rSelf[0]*uSelf[12]*mnuSelf+0.3535533905932737*m1rSelf[12]*mnuSelf-0.125*m0rSelf[4]*uSelf[11]*mnuSelf-0.125*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(72,11) = (-0.1118033988749895*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[17]*mnuSelf-0.125*m0rSelf[5]*uSelf[16]*mnuSelf-0.125*m0rSelf[6]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[14]*mnuSelf-0.125*m0rSelf[0]*uSelf[14]*mnuSelf+0.3535533905932737*m1rSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[12]*mnuSelf-0.125*m0rSelf[2]*uSelf[11]*mnuSelf-0.125*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(72,12) = (-0.125*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1964285714285714*m0rSelf[8]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[18]*mnuSelf+0.3162277660168379*m1rSelf[18]*mnuSelf-0.125*m0rSelf[7]*uSelf[17]*mnuSelf-0.225*m0rSelf[6]*uSelf[16]*mnuSelf-0.125*m0rSelf[5]*uSelf[15]*mnuSelf-0.225*m0rSelf[4]*uSelf[14]*mnuSelf-0.125*m0rSelf[3]*uSelf[13]*mnuSelf-0.225*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[1]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(72,13) = (-0.1118033988749895*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[16]*mnuSelf-0.125*m0rSelf[0]*uSelf[16]*mnuSelf+0.3535533905932737*m1rSelf[16]*mnuSelf-0.125*m0rSelf[4]*uSelf[15]*mnuSelf-0.125*m0rSelf[5]*uSelf[14]*mnuSelf-0.125*m0rSelf[2]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(72,14) = (-0.1118033988749895*m0rSelf[1]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[17]*mnuSelf-0.125*m0rSelf[3]*uSelf[15]*mnuSelf-0.225*m0rSelf[2]*uSelf[14]*mnuSelf-0.125*m0rSelf[5]*uSelf[13]*mnuSelf-0.225*m0rSelf[4]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[11]*mnuSelf-0.125*m0rSelf[0]*uSelf[11]*mnuSelf+0.3535533905932737*m1rSelf[11]*mnuSelf-0.125*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(72,15) = (-0.125*m0rSelf[1]*uSelf[16]*mnuSelf)-0.125*m0rSelf[2]*uSelf[15]*mnuSelf-0.125*m0rSelf[3]*uSelf[14]*mnuSelf-0.125*m0rSelf[4]*uSelf[13]*mnuSelf-0.125*m0rSelf[5]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(72,16) = (-0.1118033988749895*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[18]*mnuSelf-0.225*m0rSelf[2]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[13]*mnuSelf-0.125*m0rSelf[0]*uSelf[13]*mnuSelf+0.3535533905932737*m1rSelf[13]*mnuSelf-0.225*m0rSelf[6]*uSelf[12]*mnuSelf-0.125*m0rSelf[5]*uSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(72,17) = (-0.125*m0rSelf[2]*uSelf[17]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[14]*mnuSelf-0.125*m0rSelf[7]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[11]*mnuSelf; 
  data->AEM_S(72,18) = (-0.1964285714285714*m0rSelf[2]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[13]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[12]*mnuSelf+0.3162277660168379*m1rSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(72,19) = (-0.125*m0rSelf[2]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[13]*mnuSelf-0.125*m0rSelf[9]*uSelf[12]*mnuSelf; 
  data->AEM_S(73,10) = (-0.1118033988749895*m0rSelf[3]*uSelf[19]*mnuSelf)-0.125*m0rSelf[2]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[13]*mnuSelf-0.125*m0rSelf[0]*uSelf[13]*mnuSelf+0.3535533905932737*m1rSelf[13]*mnuSelf-0.125*m0rSelf[6]*uSelf[12]*mnuSelf-0.125*m0rSelf[5]*uSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(73,11) = (-0.1118033988749895*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[17]*mnuSelf-0.125*m0rSelf[4]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[15]*mnuSelf-0.125*m0rSelf[0]*uSelf[15]*mnuSelf+0.3535533905932737*m1rSelf[15]*mnuSelf-0.125*m0rSelf[6]*uSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[11]*mnuSelf-0.125*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(73,12) = (-0.1118033988749895*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[16]*mnuSelf-0.125*m0rSelf[0]*uSelf[16]*mnuSelf+0.3535533905932737*m1rSelf[16]*mnuSelf-0.125*m0rSelf[4]*uSelf[15]*mnuSelf-0.125*m0rSelf[5]*uSelf[14]*mnuSelf-0.125*m0rSelf[2]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(73,13) = (-0.1964285714285714*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[0]*uSelf[19]*mnuSelf+0.3162277660168379*m1rSelf[19]*mnuSelf-0.125*m0rSelf[8]*uSelf[18]*mnuSelf-0.125*m0rSelf[7]*uSelf[17]*mnuSelf-0.225*m0rSelf[6]*uSelf[16]*mnuSelf-0.225*m0rSelf[5]*uSelf[15]*mnuSelf-0.125*m0rSelf[4]*uSelf[14]*mnuSelf-0.225*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[1]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(73,14) = (-0.125*m0rSelf[1]*uSelf[16]*mnuSelf)-0.125*m0rSelf[2]*uSelf[15]*mnuSelf-0.125*m0rSelf[3]*uSelf[14]*mnuSelf-0.125*m0rSelf[4]*uSelf[13]*mnuSelf-0.125*m0rSelf[5]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(73,15) = (-0.1118033988749895*m0rSelf[1]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[17]*mnuSelf-0.225*m0rSelf[3]*uSelf[15]*mnuSelf-0.125*m0rSelf[2]*uSelf[14]*mnuSelf-0.225*m0rSelf[5]*uSelf[13]*mnuSelf-0.125*m0rSelf[4]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[11]*mnuSelf-0.125*m0rSelf[0]*uSelf[11]*mnuSelf+0.3535533905932737*m1rSelf[11]*mnuSelf-0.125*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(73,16) = (-0.1118033988749895*m0rSelf[2]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[18]*mnuSelf-0.225*m0rSelf[3]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[14]*mnuSelf-0.225*m0rSelf[6]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[12]*mnuSelf-0.125*m0rSelf[0]*uSelf[12]*mnuSelf+0.3535533905932737*m1rSelf[12]*mnuSelf-0.125*m0rSelf[4]*uSelf[11]*mnuSelf-0.125*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(73,17) = (-0.125*m0rSelf[3]*uSelf[17]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[15]*mnuSelf-0.125*m0rSelf[7]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[11]*mnuSelf; 
  data->AEM_S(73,18) = (-0.125*m0rSelf[3]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[16]*mnuSelf-0.125*m0rSelf[8]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[12]*mnuSelf; 
  data->AEM_S(73,19) = (-0.1964285714285714*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[13]*mnuSelf+0.3162277660168379*m1rSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,10) = (-0.1118033988749895*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[17]*mnuSelf-0.125*m0rSelf[5]*uSelf[16]*mnuSelf-0.125*m0rSelf[6]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[14]*mnuSelf-0.125*m0rSelf[0]*uSelf[14]*mnuSelf+0.3535533905932737*m1rSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[12]*mnuSelf-0.125*m0rSelf[2]*uSelf[11]*mnuSelf-0.125*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,11) = (-0.1118033988749895*m0rSelf[2]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[17]*mnuSelf-0.125*m0rSelf[3]*uSelf[16]*mnuSelf-0.225*m0rSelf[1]*uSelf[14]*mnuSelf-0.125*m0rSelf[6]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[12]*mnuSelf-0.125*m0rSelf[0]*uSelf[12]*mnuSelf+0.3535533905932737*m1rSelf[12]*mnuSelf-0.225*m0rSelf[4]*uSelf[11]*mnuSelf-0.125*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,12) = (-0.1118033988749895*m0rSelf[1]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[17]*mnuSelf-0.125*m0rSelf[3]*uSelf[15]*mnuSelf-0.225*m0rSelf[2]*uSelf[14]*mnuSelf-0.125*m0rSelf[5]*uSelf[13]*mnuSelf-0.225*m0rSelf[4]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[11]*mnuSelf-0.125*m0rSelf[0]*uSelf[11]*mnuSelf+0.3535533905932737*m1rSelf[11]*mnuSelf-0.125*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,13) = (-0.125*m0rSelf[1]*uSelf[16]*mnuSelf)-0.125*m0rSelf[2]*uSelf[15]*mnuSelf-0.125*m0rSelf[3]*uSelf[14]*mnuSelf-0.125*m0rSelf[4]*uSelf[13]*mnuSelf-0.125*m0rSelf[5]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(74,14) = (-0.125*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1964285714285714*m0rSelf[8]*uSelf[18]*mnuSelf-0.1*m0rSelf[7]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[18]*mnuSelf+0.3162277660168379*m1rSelf[18]*mnuSelf-0.1*m0rSelf[8]*uSelf[17]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[17]*mnuSelf+0.3162277660168379*m1rSelf[17]*mnuSelf-0.225*m0rSelf[6]*uSelf[16]*mnuSelf-0.225*m0rSelf[5]*uSelf[15]*mnuSelf-0.405*m0rSelf[4]*uSelf[14]*mnuSelf-0.125*m0rSelf[3]*uSelf[13]*mnuSelf-0.225*m0rSelf[2]*uSelf[12]*mnuSelf-0.225*m0rSelf[1]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[10]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(74,15) = (-0.1118033988749895*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[16]*mnuSelf-0.125*m0rSelf[0]*uSelf[16]*mnuSelf+0.3535533905932737*m1rSelf[16]*mnuSelf-0.225*m0rSelf[4]*uSelf[15]*mnuSelf-0.225*m0rSelf[5]*uSelf[14]*mnuSelf-0.125*m0rSelf[2]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,16) = (-0.1118033988749895*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[17]*mnuSelf-0.225*m0rSelf[4]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[15]*mnuSelf-0.125*m0rSelf[0]*uSelf[15]*mnuSelf+0.3535533905932737*m1rSelf[15]*mnuSelf-0.225*m0rSelf[6]*uSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[11]*mnuSelf-0.125*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,17) = (-0.1*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1964285714285714*m0rSelf[4]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[15]*mnuSelf-0.1*m0rSelf[8]*uSelf[14]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[14]*mnuSelf+0.3162277660168379*m1rSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,18) = (-0.1964285714285714*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1*m0rSelf[4]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[14]*mnuSelf-0.1*m0rSelf[7]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[14]*mnuSelf+0.3162277660168379*m1rSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(74,19) = (-0.125*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[15]*mnuSelf-0.125*m0rSelf[9]*uSelf[14]*mnuSelf; 
  data->AEM_S(75,10) = (-0.1118033988749895*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[17]*mnuSelf-0.125*m0rSelf[4]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[15]*mnuSelf-0.125*m0rSelf[0]*uSelf[15]*mnuSelf+0.3535533905932737*m1rSelf[15]*mnuSelf-0.125*m0rSelf[6]*uSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[11]*mnuSelf-0.125*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(75,11) = (-0.1118033988749895*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[17]*mnuSelf-0.125*m0rSelf[2]*uSelf[16]*mnuSelf-0.225*m0rSelf[1]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[13]*mnuSelf-0.125*m0rSelf[0]*uSelf[13]*mnuSelf+0.3535533905932737*m1rSelf[13]*mnuSelf-0.125*m0rSelf[6]*uSelf[12]*mnuSelf-0.225*m0rSelf[5]*uSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(75,12) = (-0.125*m0rSelf[1]*uSelf[16]*mnuSelf)-0.125*m0rSelf[2]*uSelf[15]*mnuSelf-0.125*m0rSelf[3]*uSelf[14]*mnuSelf-0.125*m0rSelf[4]*uSelf[13]*mnuSelf-0.125*m0rSelf[5]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(75,13) = (-0.1118033988749895*m0rSelf[1]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[17]*mnuSelf-0.225*m0rSelf[3]*uSelf[15]*mnuSelf-0.125*m0rSelf[2]*uSelf[14]*mnuSelf-0.225*m0rSelf[5]*uSelf[13]*mnuSelf-0.125*m0rSelf[4]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[11]*mnuSelf-0.125*m0rSelf[0]*uSelf[11]*mnuSelf+0.3535533905932737*m1rSelf[11]*mnuSelf-0.125*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(75,14) = (-0.1118033988749895*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[16]*mnuSelf-0.125*m0rSelf[0]*uSelf[16]*mnuSelf+0.3535533905932737*m1rSelf[16]*mnuSelf-0.225*m0rSelf[4]*uSelf[15]*mnuSelf-0.225*m0rSelf[5]*uSelf[14]*mnuSelf-0.125*m0rSelf[2]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(75,15) = (-0.1964285714285714*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1*m0rSelf[7]*uSelf[19]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[19]*mnuSelf+0.3162277660168379*m1rSelf[19]*mnuSelf-0.125*m0rSelf[8]*uSelf[18]*mnuSelf-0.1*m0rSelf[9]*uSelf[17]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[17]*mnuSelf+0.3162277660168379*m1rSelf[17]*mnuSelf-0.225*m0rSelf[6]*uSelf[16]*mnuSelf-0.405*m0rSelf[5]*uSelf[15]*mnuSelf-0.225*m0rSelf[4]*uSelf[14]*mnuSelf-0.225*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[2]*uSelf[12]*mnuSelf-0.225*m0rSelf[1]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[10]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(75,16) = (-0.1118033988749895*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[17]*mnuSelf-0.225*m0rSelf[5]*uSelf[16]*mnuSelf-0.225*m0rSelf[6]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[14]*mnuSelf-0.125*m0rSelf[0]*uSelf[14]*mnuSelf+0.3535533905932737*m1rSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[12]*mnuSelf-0.125*m0rSelf[2]*uSelf[11]*mnuSelf-0.125*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(75,17) = (-0.1*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1964285714285714*m0rSelf[5]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[16]*mnuSelf-0.1*m0rSelf[9]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[15]*mnuSelf+0.3162277660168379*m1rSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(75,18) = (-0.125*m0rSelf[5]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[16]*mnuSelf-0.125*m0rSelf[8]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[14]*mnuSelf; 
  data->AEM_S(75,19) = (-0.1964285714285714*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1*m0rSelf[5]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[16]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[15]*mnuSelf-0.1*m0rSelf[7]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[15]*mnuSelf+0.3162277660168379*m1rSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(76,10) = (-0.1118033988749895*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[16]*mnuSelf-0.125*m0rSelf[0]*uSelf[16]*mnuSelf+0.3535533905932737*m1rSelf[16]*mnuSelf-0.125*m0rSelf[4]*uSelf[15]*mnuSelf-0.125*m0rSelf[5]*uSelf[14]*mnuSelf-0.125*m0rSelf[2]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(76,11) = (-0.125*m0rSelf[1]*uSelf[16]*mnuSelf)-0.125*m0rSelf[2]*uSelf[15]*mnuSelf-0.125*m0rSelf[3]*uSelf[14]*mnuSelf-0.125*m0rSelf[4]*uSelf[13]*mnuSelf-0.125*m0rSelf[5]*uSelf[12]*mnuSelf-0.125*m0rSelf[6]*uSelf[11]*mnuSelf; 
  data->AEM_S(76,12) = (-0.1118033988749895*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[18]*mnuSelf-0.225*m0rSelf[2]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[13]*mnuSelf-0.125*m0rSelf[0]*uSelf[13]*mnuSelf+0.3535533905932737*m1rSelf[13]*mnuSelf-0.225*m0rSelf[6]*uSelf[12]*mnuSelf-0.125*m0rSelf[5]*uSelf[11]*mnuSelf-0.125*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(76,13) = (-0.1118033988749895*m0rSelf[2]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[18]*mnuSelf-0.225*m0rSelf[3]*uSelf[16]*mnuSelf-0.125*m0rSelf[1]*uSelf[14]*mnuSelf-0.225*m0rSelf[6]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[12]*mnuSelf-0.125*m0rSelf[0]*uSelf[12]*mnuSelf+0.3535533905932737*m1rSelf[12]*mnuSelf-0.125*m0rSelf[4]*uSelf[11]*mnuSelf-0.125*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(76,14) = (-0.1118033988749895*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[17]*mnuSelf-0.225*m0rSelf[4]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[15]*mnuSelf-0.125*m0rSelf[0]*uSelf[15]*mnuSelf+0.3535533905932737*m1rSelf[15]*mnuSelf-0.225*m0rSelf[6]*uSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[13]*mnuSelf-0.125*m0rSelf[3]*uSelf[11]*mnuSelf-0.125*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(76,15) = (-0.1118033988749895*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[17]*mnuSelf-0.225*m0rSelf[5]*uSelf[16]*mnuSelf-0.225*m0rSelf[6]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[14]*mnuSelf-0.125*m0rSelf[0]*uSelf[14]*mnuSelf+0.3535533905932737*m1rSelf[14]*mnuSelf-0.125*m0rSelf[1]*uSelf[12]*mnuSelf-0.125*m0rSelf[2]*uSelf[11]*mnuSelf-0.125*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(76,16) = (-0.1964285714285714*m0rSelf[9]*uSelf[19]*mnuSelf)-0.1*m0rSelf[8]*uSelf[19]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[19]*mnuSelf+0.3162277660168379*m1rSelf[19]*mnuSelf-0.1*m0rSelf[9]*uSelf[18]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[18]*mnuSelf+0.3162277660168379*m1rSelf[18]*mnuSelf-0.125*m0rSelf[7]*uSelf[17]*mnuSelf-0.405*m0rSelf[6]*uSelf[16]*mnuSelf-0.225*m0rSelf[5]*uSelf[15]*mnuSelf-0.225*m0rSelf[4]*uSelf[14]*mnuSelf-0.225*m0rSelf[3]*uSelf[13]*mnuSelf-0.225*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[1]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[10]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(76,17) = (-0.125*m0rSelf[6]*uSelf[17]*mnuSelf)-0.125*m0rSelf[7]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[14]*mnuSelf; 
  data->AEM_S(76,18) = (-0.1*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1964285714285714*m0rSelf[6]*uSelf[18]*mnuSelf-0.1*m0rSelf[9]*uSelf[16]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[16]*mnuSelf+0.3162277660168379*m1rSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(76,19) = (-0.1964285714285714*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1*m0rSelf[6]*uSelf[18]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[16]*mnuSelf-0.1*m0rSelf[8]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[16]*mnuSelf+0.3162277660168379*m1rSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(77,10) = (-0.07985957062499249*m0rSelf[7]*uSelf[17]*mnuSelf)-0.125*m0rSelf[0]*uSelf[17]*mnuSelf+0.3535533905932737*m1rSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[11]*mnuSelf-0.125*m0rSelf[7]*uSelf[10]*mnuSelf; 
  data->AEM_S(77,11) = (-0.1964285714285714*m0rSelf[1]*uSelf[17]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[12]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[11]*mnuSelf+0.3162277660168379*m1rSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[10]*mnuSelf; 
  data->AEM_S(77,12) = (-0.125*m0rSelf[2]*uSelf[17]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[14]*mnuSelf-0.125*m0rSelf[7]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[11]*mnuSelf; 
  data->AEM_S(77,13) = (-0.125*m0rSelf[3]*uSelf[17]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[15]*mnuSelf-0.125*m0rSelf[7]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[11]*mnuSelf; 
  data->AEM_S(77,14) = (-0.1*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1964285714285714*m0rSelf[4]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[15]*mnuSelf-0.1*m0rSelf[8]*uSelf[14]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[14]*mnuSelf+0.3162277660168379*m1rSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(77,15) = (-0.1*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1964285714285714*m0rSelf[5]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[16]*mnuSelf-0.1*m0rSelf[9]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[15]*mnuSelf+0.3162277660168379*m1rSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(77,16) = (-0.125*m0rSelf[6]*uSelf[17]*mnuSelf)-0.125*m0rSelf[7]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[14]*mnuSelf; 
  data->AEM_S(77,17) = (-0.125*m0rSelf[9]*uSelf[19]*mnuSelf)-0.125*m0rSelf[8]*uSelf[18]*mnuSelf-0.2678571428571428*m0rSelf[7]*uSelf[17]*mnuSelf-0.07985957062499249*m0rSelf[0]*uSelf[17]*mnuSelf+0.2258769757263128*m1rSelf[17]*mnuSelf-0.125*m0rSelf[6]*uSelf[16]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[14]*mnuSelf-0.125*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[2]*uSelf[12]*mnuSelf-0.1964285714285714*m0rSelf[1]*uSelf[11]*mnuSelf-0.07985957062499249*m0rSelf[7]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(77,18) = (-0.125*m0rSelf[7]*uSelf[18]*mnuSelf)-0.125*m0rSelf[8]*uSelf[17]*mnuSelf-0.1*m0rSelf[4]*uSelf[14]*mnuSelf; 
  data->AEM_S(77,19) = (-0.125*m0rSelf[7]*uSelf[19]*mnuSelf)-0.125*m0rSelf[9]*uSelf[17]*mnuSelf-0.1*m0rSelf[5]*uSelf[15]*mnuSelf; 
  data->AEM_S(78,10) = (-0.07985957062499249*m0rSelf[8]*uSelf[18]*mnuSelf)-0.125*m0rSelf[0]*uSelf[18]*mnuSelf+0.3535533905932737*m1rSelf[18]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[8]*uSelf[10]*mnuSelf; 
  data->AEM_S(78,11) = (-0.125*m0rSelf[1]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[12]*mnuSelf-0.125*m0rSelf[8]*uSelf[11]*mnuSelf; 
  data->AEM_S(78,12) = (-0.1964285714285714*m0rSelf[2]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[13]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[12]*mnuSelf+0.3162277660168379*m1rSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[10]*mnuSelf; 
  data->AEM_S(78,13) = (-0.125*m0rSelf[3]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[16]*mnuSelf-0.125*m0rSelf[8]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[12]*mnuSelf; 
  data->AEM_S(78,14) = (-0.1964285714285714*m0rSelf[4]*uSelf[18]*mnuSelf)-0.1*m0rSelf[4]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[14]*mnuSelf-0.1*m0rSelf[7]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[14]*mnuSelf+0.3162277660168379*m1rSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[10]*mnuSelf; 
  data->AEM_S(78,15) = (-0.125*m0rSelf[5]*uSelf[18]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[16]*mnuSelf-0.125*m0rSelf[8]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[14]*mnuSelf; 
  data->AEM_S(78,16) = (-0.1*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1964285714285714*m0rSelf[6]*uSelf[18]*mnuSelf-0.1*m0rSelf[9]*uSelf[16]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[16]*mnuSelf+0.3162277660168379*m1rSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(78,17) = (-0.125*m0rSelf[7]*uSelf[18]*mnuSelf)-0.125*m0rSelf[8]*uSelf[17]*mnuSelf-0.1*m0rSelf[4]*uSelf[14]*mnuSelf; 
  data->AEM_S(78,18) = (-0.125*m0rSelf[9]*uSelf[19]*mnuSelf)-0.2678571428571428*m0rSelf[8]*uSelf[18]*mnuSelf-0.07985957062499249*m0rSelf[0]*uSelf[18]*mnuSelf+0.2258769757263128*m1rSelf[18]*mnuSelf-0.125*m0rSelf[7]*uSelf[17]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[16]*mnuSelf-0.125*m0rSelf[5]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[14]*mnuSelf-0.125*m0rSelf[3]*uSelf[13]*mnuSelf-0.1964285714285714*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[1]*uSelf[11]*mnuSelf-0.07985957062499249*m0rSelf[8]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
  data->AEM_S(78,19) = (-0.125*m0rSelf[8]*uSelf[19]*mnuSelf)-0.125*m0rSelf[9]*uSelf[18]*mnuSelf-0.1*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(79,10) = (-0.07985957062499249*m0rSelf[9]*uSelf[19]*mnuSelf)-0.125*m0rSelf[0]*uSelf[19]*mnuSelf+0.3535533905932737*m1rSelf[19]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[9]*uSelf[10]*mnuSelf; 
  data->AEM_S(79,11) = (-0.125*m0rSelf[1]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[13]*mnuSelf-0.125*m0rSelf[9]*uSelf[11]*mnuSelf; 
  data->AEM_S(79,12) = (-0.125*m0rSelf[2]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[13]*mnuSelf-0.125*m0rSelf[9]*uSelf[12]*mnuSelf; 
  data->AEM_S(79,13) = (-0.1964285714285714*m0rSelf[3]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[15]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[13]*mnuSelf+0.3162277660168379*m1rSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[10]*mnuSelf; 
  data->AEM_S(79,14) = (-0.125*m0rSelf[4]*uSelf[19]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[15]*mnuSelf-0.125*m0rSelf[9]*uSelf[14]*mnuSelf; 
  data->AEM_S(79,15) = (-0.1964285714285714*m0rSelf[5]*uSelf[19]*mnuSelf)-0.1*m0rSelf[5]*uSelf[17]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[16]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[15]*mnuSelf-0.1*m0rSelf[7]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[15]*mnuSelf+0.3162277660168379*m1rSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[11]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[10]*mnuSelf; 
  data->AEM_S(79,16) = (-0.1964285714285714*m0rSelf[6]*uSelf[19]*mnuSelf)-0.1*m0rSelf[6]*uSelf[18]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[16]*mnuSelf-0.1*m0rSelf[8]*uSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[16]*mnuSelf+0.3162277660168379*m1rSelf[16]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[15]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[14]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[13]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[12]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[10]*mnuSelf; 
  data->AEM_S(79,17) = (-0.125*m0rSelf[7]*uSelf[19]*mnuSelf)-0.125*m0rSelf[9]*uSelf[17]*mnuSelf-0.1*m0rSelf[5]*uSelf[15]*mnuSelf; 
  data->AEM_S(79,18) = (-0.125*m0rSelf[8]*uSelf[19]*mnuSelf)-0.125*m0rSelf[9]*uSelf[18]*mnuSelf-0.1*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(79,19) = (-0.2678571428571428*m0rSelf[9]*uSelf[19]*mnuSelf)-0.07985957062499249*m0rSelf[0]*uSelf[19]*mnuSelf+0.2258769757263128*m1rSelf[19]*mnuSelf-0.125*m0rSelf[8]*uSelf[18]*mnuSelf-0.125*m0rSelf[7]*uSelf[17]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[16]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[15]*mnuSelf-0.125*m0rSelf[4]*uSelf[14]*mnuSelf-0.1964285714285714*m0rSelf[3]*uSelf[13]*mnuSelf-0.125*m0rSelf[2]*uSelf[12]*mnuSelf-0.125*m0rSelf[1]*uSelf[11]*mnuSelf-0.07985957062499249*m0rSelf[9]*uSelf[10]*mnuSelf-0.125*m0rSelf[0]*uSelf[10]*mnuSelf+0.3535533905932737*m1rSelf[10]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(70,50) = 0.125*m0rOther[9]*uOther[19]*mnuOther+0.125*m0rOther[8]*uOther[18]*mnuOther+0.125*m0rOther[7]*uOther[17]*mnuOther+0.125*m0rOther[6]*uOther[16]*mnuOther+0.125*m0rOther[5]*uOther[15]*mnuOther+0.125*m0rOther[4]*uOther[14]*mnuOther+0.125*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[1]*uOther[11]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(70,51) = 0.1118033988749895*m0rOther[1]*uOther[17]*mnuOther+0.125*m0rOther[3]*uOther[15]*mnuOther+0.125*m0rOther[2]*uOther[14]*mnuOther+0.125*m0rOther[5]*uOther[13]*mnuOther+0.125*m0rOther[4]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[11]*mnuOther+0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1rOther[11]*mnuOther+0.125*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(70,52) = 0.1118033988749895*m0rOther[2]*uOther[18]*mnuOther+0.125*m0rOther[3]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[14]*mnuOther+0.125*m0rOther[6]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[12]*mnuOther+0.125*m0rOther[0]*uOther[12]*mnuOther-0.3535533905932737*m1rOther[12]*mnuOther+0.125*m0rOther[4]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(70,53) = 0.1118033988749895*m0rOther[3]*uOther[19]*mnuOther+0.125*m0rOther[2]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[13]*mnuOther+0.125*m0rOther[0]*uOther[13]*mnuOther-0.3535533905932737*m1rOther[13]*mnuOther+0.125*m0rOther[6]*uOther[12]*mnuOther+0.125*m0rOther[5]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(70,54) = 0.1118033988749895*m0rOther[4]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[17]*mnuOther+0.125*m0rOther[5]*uOther[16]*mnuOther+0.125*m0rOther[6]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[14]*mnuOther+0.125*m0rOther[0]*uOther[14]*mnuOther-0.3535533905932737*m1rOther[14]*mnuOther+0.125*m0rOther[1]*uOther[12]*mnuOther+0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(70,55) = 0.1118033988749895*m0rOther[5]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[17]*mnuOther+0.125*m0rOther[4]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[15]*mnuOther+0.125*m0rOther[0]*uOther[15]*mnuOther-0.3535533905932737*m1rOther[15]*mnuOther+0.125*m0rOther[6]*uOther[14]*mnuOther+0.125*m0rOther[1]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(70,56) = 0.1118033988749895*m0rOther[6]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[16]*mnuOther+0.125*m0rOther[0]*uOther[16]*mnuOther-0.3535533905932737*m1rOther[16]*mnuOther+0.125*m0rOther[4]*uOther[15]*mnuOther+0.125*m0rOther[5]*uOther[14]*mnuOther+0.125*m0rOther[2]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(70,57) = 0.07985957062499249*m0rOther[7]*uOther[17]*mnuOther+0.125*m0rOther[0]*uOther[17]*mnuOther-0.3535533905932737*m1rOther[17]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[11]*mnuOther+0.125*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(70,58) = 0.07985957062499249*m0rOther[8]*uOther[18]*mnuOther+0.125*m0rOther[0]*uOther[18]*mnuOther-0.3535533905932737*m1rOther[18]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[8]*uOther[10]*mnuOther; 
  data->AEM_S(70,59) = 0.07985957062499249*m0rOther[9]*uOther[19]*mnuOther+0.125*m0rOther[0]*uOther[19]*mnuOther-0.3535533905932737*m1rOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[9]*uOther[10]*mnuOther; 
  data->AEM_S(71,50) = 0.1118033988749895*m0rOther[1]*uOther[17]*mnuOther+0.125*m0rOther[3]*uOther[15]*mnuOther+0.125*m0rOther[2]*uOther[14]*mnuOther+0.125*m0rOther[5]*uOther[13]*mnuOther+0.125*m0rOther[4]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[11]*mnuOther+0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1rOther[11]*mnuOther+0.125*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(71,51) = 0.125*m0rOther[9]*uOther[19]*mnuOther+0.125*m0rOther[8]*uOther[18]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[17]*mnuOther-0.3162277660168379*m1rOther[17]*mnuOther+0.125*m0rOther[6]*uOther[16]*mnuOther+0.225*m0rOther[5]*uOther[15]*mnuOther+0.225*m0rOther[4]*uOther[14]*mnuOther+0.125*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[2]*uOther[12]*mnuOther+0.225*m0rOther[1]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(71,52) = 0.1118033988749895*m0rOther[4]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[17]*mnuOther+0.125*m0rOther[5]*uOther[16]*mnuOther+0.125*m0rOther[6]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[14]*mnuOther+0.125*m0rOther[0]*uOther[14]*mnuOther-0.3535533905932737*m1rOther[14]*mnuOther+0.125*m0rOther[1]*uOther[12]*mnuOther+0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(71,53) = 0.1118033988749895*m0rOther[5]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[17]*mnuOther+0.125*m0rOther[4]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[15]*mnuOther+0.125*m0rOther[0]*uOther[15]*mnuOther-0.3535533905932737*m1rOther[15]*mnuOther+0.125*m0rOther[6]*uOther[14]*mnuOther+0.125*m0rOther[1]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(71,54) = 0.1118033988749895*m0rOther[2]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[17]*mnuOther+0.125*m0rOther[3]*uOther[16]*mnuOther+0.225*m0rOther[1]*uOther[14]*mnuOther+0.125*m0rOther[6]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[12]*mnuOther+0.125*m0rOther[0]*uOther[12]*mnuOther-0.3535533905932737*m1rOther[12]*mnuOther+0.225*m0rOther[4]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(71,55) = 0.1118033988749895*m0rOther[3]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[17]*mnuOther+0.125*m0rOther[2]*uOther[16]*mnuOther+0.225*m0rOther[1]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[13]*mnuOther+0.125*m0rOther[0]*uOther[13]*mnuOther-0.3535533905932737*m1rOther[13]*mnuOther+0.125*m0rOther[6]*uOther[12]*mnuOther+0.225*m0rOther[5]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(71,56) = 0.125*m0rOther[1]*uOther[16]*mnuOther+0.125*m0rOther[2]*uOther[15]*mnuOther+0.125*m0rOther[3]*uOther[14]*mnuOther+0.125*m0rOther[4]*uOther[13]*mnuOther+0.125*m0rOther[5]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(71,57) = 0.1964285714285714*m0rOther[1]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[12]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[11]*mnuOther-0.3162277660168379*m1rOther[11]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(71,58) = 0.125*m0rOther[1]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[12]*mnuOther+0.125*m0rOther[8]*uOther[11]*mnuOther; 
  data->AEM_S(71,59) = 0.125*m0rOther[1]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[13]*mnuOther+0.125*m0rOther[9]*uOther[11]*mnuOther; 
  data->AEM_S(72,50) = 0.1118033988749895*m0rOther[2]*uOther[18]*mnuOther+0.125*m0rOther[3]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[14]*mnuOther+0.125*m0rOther[6]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[12]*mnuOther+0.125*m0rOther[0]*uOther[12]*mnuOther-0.3535533905932737*m1rOther[12]*mnuOther+0.125*m0rOther[4]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(72,51) = 0.1118033988749895*m0rOther[4]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[17]*mnuOther+0.125*m0rOther[5]*uOther[16]*mnuOther+0.125*m0rOther[6]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[14]*mnuOther+0.125*m0rOther[0]*uOther[14]*mnuOther-0.3535533905932737*m1rOther[14]*mnuOther+0.125*m0rOther[1]*uOther[12]*mnuOther+0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(72,52) = 0.125*m0rOther[9]*uOther[19]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[18]*mnuOther-0.3162277660168379*m1rOther[18]*mnuOther+0.125*m0rOther[7]*uOther[17]*mnuOther+0.225*m0rOther[6]*uOther[16]*mnuOther+0.125*m0rOther[5]*uOther[15]*mnuOther+0.225*m0rOther[4]*uOther[14]*mnuOther+0.125*m0rOther[3]*uOther[13]*mnuOther+0.225*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[1]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(72,53) = 0.1118033988749895*m0rOther[6]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[16]*mnuOther+0.125*m0rOther[0]*uOther[16]*mnuOther-0.3535533905932737*m1rOther[16]*mnuOther+0.125*m0rOther[4]*uOther[15]*mnuOther+0.125*m0rOther[5]*uOther[14]*mnuOther+0.125*m0rOther[2]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(72,54) = 0.1118033988749895*m0rOther[1]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[17]*mnuOther+0.125*m0rOther[3]*uOther[15]*mnuOther+0.225*m0rOther[2]*uOther[14]*mnuOther+0.125*m0rOther[5]*uOther[13]*mnuOther+0.225*m0rOther[4]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[11]*mnuOther+0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1rOther[11]*mnuOther+0.125*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(72,55) = 0.125*m0rOther[1]*uOther[16]*mnuOther+0.125*m0rOther[2]*uOther[15]*mnuOther+0.125*m0rOther[3]*uOther[14]*mnuOther+0.125*m0rOther[4]*uOther[13]*mnuOther+0.125*m0rOther[5]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(72,56) = 0.1118033988749895*m0rOther[3]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[18]*mnuOther+0.225*m0rOther[2]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[13]*mnuOther+0.125*m0rOther[0]*uOther[13]*mnuOther-0.3535533905932737*m1rOther[13]*mnuOther+0.225*m0rOther[6]*uOther[12]*mnuOther+0.125*m0rOther[5]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(72,57) = 0.125*m0rOther[2]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[14]*mnuOther+0.125*m0rOther[7]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[11]*mnuOther; 
  data->AEM_S(72,58) = 0.1964285714285714*m0rOther[2]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[13]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[12]*mnuOther-0.3162277660168379*m1rOther[12]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(72,59) = 0.125*m0rOther[2]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[13]*mnuOther+0.125*m0rOther[9]*uOther[12]*mnuOther; 
  data->AEM_S(73,50) = 0.1118033988749895*m0rOther[3]*uOther[19]*mnuOther+0.125*m0rOther[2]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[13]*mnuOther+0.125*m0rOther[0]*uOther[13]*mnuOther-0.3535533905932737*m1rOther[13]*mnuOther+0.125*m0rOther[6]*uOther[12]*mnuOther+0.125*m0rOther[5]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(73,51) = 0.1118033988749895*m0rOther[5]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[17]*mnuOther+0.125*m0rOther[4]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[15]*mnuOther+0.125*m0rOther[0]*uOther[15]*mnuOther-0.3535533905932737*m1rOther[15]*mnuOther+0.125*m0rOther[6]*uOther[14]*mnuOther+0.125*m0rOther[1]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(73,52) = 0.1118033988749895*m0rOther[6]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[16]*mnuOther+0.125*m0rOther[0]*uOther[16]*mnuOther-0.3535533905932737*m1rOther[16]*mnuOther+0.125*m0rOther[4]*uOther[15]*mnuOther+0.125*m0rOther[5]*uOther[14]*mnuOther+0.125*m0rOther[2]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(73,53) = 0.1964285714285714*m0rOther[9]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[19]*mnuOther-0.3162277660168379*m1rOther[19]*mnuOther+0.125*m0rOther[8]*uOther[18]*mnuOther+0.125*m0rOther[7]*uOther[17]*mnuOther+0.225*m0rOther[6]*uOther[16]*mnuOther+0.225*m0rOther[5]*uOther[15]*mnuOther+0.125*m0rOther[4]*uOther[14]*mnuOther+0.225*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[1]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(73,54) = 0.125*m0rOther[1]*uOther[16]*mnuOther+0.125*m0rOther[2]*uOther[15]*mnuOther+0.125*m0rOther[3]*uOther[14]*mnuOther+0.125*m0rOther[4]*uOther[13]*mnuOther+0.125*m0rOther[5]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(73,55) = 0.1118033988749895*m0rOther[1]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[17]*mnuOther+0.225*m0rOther[3]*uOther[15]*mnuOther+0.125*m0rOther[2]*uOther[14]*mnuOther+0.225*m0rOther[5]*uOther[13]*mnuOther+0.125*m0rOther[4]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[11]*mnuOther+0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1rOther[11]*mnuOther+0.125*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(73,56) = 0.1118033988749895*m0rOther[2]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[18]*mnuOther+0.225*m0rOther[3]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[14]*mnuOther+0.225*m0rOther[6]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[12]*mnuOther+0.125*m0rOther[0]*uOther[12]*mnuOther-0.3535533905932737*m1rOther[12]*mnuOther+0.125*m0rOther[4]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(73,57) = 0.125*m0rOther[3]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[15]*mnuOther+0.125*m0rOther[7]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[11]*mnuOther; 
  data->AEM_S(73,58) = 0.125*m0rOther[3]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[16]*mnuOther+0.125*m0rOther[8]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[12]*mnuOther; 
  data->AEM_S(73,59) = 0.1964285714285714*m0rOther[3]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[13]*mnuOther-0.3162277660168379*m1rOther[13]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(74,50) = 0.1118033988749895*m0rOther[4]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[17]*mnuOther+0.125*m0rOther[5]*uOther[16]*mnuOther+0.125*m0rOther[6]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[14]*mnuOther+0.125*m0rOther[0]*uOther[14]*mnuOther-0.3535533905932737*m1rOther[14]*mnuOther+0.125*m0rOther[1]*uOther[12]*mnuOther+0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(74,51) = 0.1118033988749895*m0rOther[2]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[17]*mnuOther+0.125*m0rOther[3]*uOther[16]*mnuOther+0.225*m0rOther[1]*uOther[14]*mnuOther+0.125*m0rOther[6]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[12]*mnuOther+0.125*m0rOther[0]*uOther[12]*mnuOther-0.3535533905932737*m1rOther[12]*mnuOther+0.225*m0rOther[4]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(74,52) = 0.1118033988749895*m0rOther[1]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[17]*mnuOther+0.125*m0rOther[3]*uOther[15]*mnuOther+0.225*m0rOther[2]*uOther[14]*mnuOther+0.125*m0rOther[5]*uOther[13]*mnuOther+0.225*m0rOther[4]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[11]*mnuOther+0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1rOther[11]*mnuOther+0.125*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(74,53) = 0.125*m0rOther[1]*uOther[16]*mnuOther+0.125*m0rOther[2]*uOther[15]*mnuOther+0.125*m0rOther[3]*uOther[14]*mnuOther+0.125*m0rOther[4]*uOther[13]*mnuOther+0.125*m0rOther[5]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(74,54) = 0.125*m0rOther[9]*uOther[19]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[18]*mnuOther+0.1*m0rOther[7]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[18]*mnuOther-0.3162277660168379*m1rOther[18]*mnuOther+0.1*m0rOther[8]*uOther[17]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[17]*mnuOther-0.3162277660168379*m1rOther[17]*mnuOther+0.225*m0rOther[6]*uOther[16]*mnuOther+0.225*m0rOther[5]*uOther[15]*mnuOther+0.405*m0rOther[4]*uOther[14]*mnuOther+0.125*m0rOther[3]*uOther[13]*mnuOther+0.225*m0rOther[2]*uOther[12]*mnuOther+0.225*m0rOther[1]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[10]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(74,55) = 0.1118033988749895*m0rOther[6]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[16]*mnuOther+0.125*m0rOther[0]*uOther[16]*mnuOther-0.3535533905932737*m1rOther[16]*mnuOther+0.225*m0rOther[4]*uOther[15]*mnuOther+0.225*m0rOther[5]*uOther[14]*mnuOther+0.125*m0rOther[2]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(74,56) = 0.1118033988749895*m0rOther[5]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[17]*mnuOther+0.225*m0rOther[4]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[15]*mnuOther+0.125*m0rOther[0]*uOther[15]*mnuOther-0.3535533905932737*m1rOther[15]*mnuOther+0.225*m0rOther[6]*uOther[14]*mnuOther+0.125*m0rOther[1]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(74,57) = 0.1*m0rOther[4]*uOther[18]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[15]*mnuOther+0.1*m0rOther[8]*uOther[14]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[14]*mnuOther-0.3162277660168379*m1rOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(74,58) = 0.1964285714285714*m0rOther[4]*uOther[18]*mnuOther+0.1*m0rOther[4]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[14]*mnuOther+0.1*m0rOther[7]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[14]*mnuOther-0.3162277660168379*m1rOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(74,59) = 0.125*m0rOther[4]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[15]*mnuOther+0.125*m0rOther[9]*uOther[14]*mnuOther; 
  data->AEM_S(75,50) = 0.1118033988749895*m0rOther[5]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[17]*mnuOther+0.125*m0rOther[4]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[15]*mnuOther+0.125*m0rOther[0]*uOther[15]*mnuOther-0.3535533905932737*m1rOther[15]*mnuOther+0.125*m0rOther[6]*uOther[14]*mnuOther+0.125*m0rOther[1]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(75,51) = 0.1118033988749895*m0rOther[3]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[17]*mnuOther+0.125*m0rOther[2]*uOther[16]*mnuOther+0.225*m0rOther[1]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[13]*mnuOther+0.125*m0rOther[0]*uOther[13]*mnuOther-0.3535533905932737*m1rOther[13]*mnuOther+0.125*m0rOther[6]*uOther[12]*mnuOther+0.225*m0rOther[5]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(75,52) = 0.125*m0rOther[1]*uOther[16]*mnuOther+0.125*m0rOther[2]*uOther[15]*mnuOther+0.125*m0rOther[3]*uOther[14]*mnuOther+0.125*m0rOther[4]*uOther[13]*mnuOther+0.125*m0rOther[5]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(75,53) = 0.1118033988749895*m0rOther[1]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[17]*mnuOther+0.225*m0rOther[3]*uOther[15]*mnuOther+0.125*m0rOther[2]*uOther[14]*mnuOther+0.225*m0rOther[5]*uOther[13]*mnuOther+0.125*m0rOther[4]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[11]*mnuOther+0.125*m0rOther[0]*uOther[11]*mnuOther-0.3535533905932737*m1rOther[11]*mnuOther+0.125*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(75,54) = 0.1118033988749895*m0rOther[6]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[16]*mnuOther+0.125*m0rOther[0]*uOther[16]*mnuOther-0.3535533905932737*m1rOther[16]*mnuOther+0.225*m0rOther[4]*uOther[15]*mnuOther+0.225*m0rOther[5]*uOther[14]*mnuOther+0.125*m0rOther[2]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(75,55) = 0.1964285714285714*m0rOther[9]*uOther[19]*mnuOther+0.1*m0rOther[7]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[19]*mnuOther-0.3162277660168379*m1rOther[19]*mnuOther+0.125*m0rOther[8]*uOther[18]*mnuOther+0.1*m0rOther[9]*uOther[17]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[17]*mnuOther-0.3162277660168379*m1rOther[17]*mnuOther+0.225*m0rOther[6]*uOther[16]*mnuOther+0.405*m0rOther[5]*uOther[15]*mnuOther+0.225*m0rOther[4]*uOther[14]*mnuOther+0.225*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[2]*uOther[12]*mnuOther+0.225*m0rOther[1]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[10]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(75,56) = 0.1118033988749895*m0rOther[4]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[17]*mnuOther+0.225*m0rOther[5]*uOther[16]*mnuOther+0.225*m0rOther[6]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[14]*mnuOther+0.125*m0rOther[0]*uOther[14]*mnuOther-0.3535533905932737*m1rOther[14]*mnuOther+0.125*m0rOther[1]*uOther[12]*mnuOther+0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(75,57) = 0.1*m0rOther[5]*uOther[19]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[16]*mnuOther+0.1*m0rOther[9]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[15]*mnuOther-0.3162277660168379*m1rOther[15]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(75,58) = 0.125*m0rOther[5]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[16]*mnuOther+0.125*m0rOther[8]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[14]*mnuOther; 
  data->AEM_S(75,59) = 0.1964285714285714*m0rOther[5]*uOther[19]*mnuOther+0.1*m0rOther[5]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[16]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[15]*mnuOther+0.1*m0rOther[7]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[15]*mnuOther-0.3162277660168379*m1rOther[15]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(76,50) = 0.1118033988749895*m0rOther[6]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[16]*mnuOther+0.125*m0rOther[0]*uOther[16]*mnuOther-0.3535533905932737*m1rOther[16]*mnuOther+0.125*m0rOther[4]*uOther[15]*mnuOther+0.125*m0rOther[5]*uOther[14]*mnuOther+0.125*m0rOther[2]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(76,51) = 0.125*m0rOther[1]*uOther[16]*mnuOther+0.125*m0rOther[2]*uOther[15]*mnuOther+0.125*m0rOther[3]*uOther[14]*mnuOther+0.125*m0rOther[4]*uOther[13]*mnuOther+0.125*m0rOther[5]*uOther[12]*mnuOther+0.125*m0rOther[6]*uOther[11]*mnuOther; 
  data->AEM_S(76,52) = 0.1118033988749895*m0rOther[3]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[18]*mnuOther+0.225*m0rOther[2]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[13]*mnuOther+0.125*m0rOther[0]*uOther[13]*mnuOther-0.3535533905932737*m1rOther[13]*mnuOther+0.225*m0rOther[6]*uOther[12]*mnuOther+0.125*m0rOther[5]*uOther[11]*mnuOther+0.125*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(76,53) = 0.1118033988749895*m0rOther[2]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[18]*mnuOther+0.225*m0rOther[3]*uOther[16]*mnuOther+0.125*m0rOther[1]*uOther[14]*mnuOther+0.225*m0rOther[6]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[12]*mnuOther+0.125*m0rOther[0]*uOther[12]*mnuOther-0.3535533905932737*m1rOther[12]*mnuOther+0.125*m0rOther[4]*uOther[11]*mnuOther+0.125*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(76,54) = 0.1118033988749895*m0rOther[5]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[17]*mnuOther+0.225*m0rOther[4]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[15]*mnuOther+0.125*m0rOther[0]*uOther[15]*mnuOther-0.3535533905932737*m1rOther[15]*mnuOther+0.225*m0rOther[6]*uOther[14]*mnuOther+0.125*m0rOther[1]*uOther[13]*mnuOther+0.125*m0rOther[3]*uOther[11]*mnuOther+0.125*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(76,55) = 0.1118033988749895*m0rOther[4]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[17]*mnuOther+0.225*m0rOther[5]*uOther[16]*mnuOther+0.225*m0rOther[6]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[14]*mnuOther+0.125*m0rOther[0]*uOther[14]*mnuOther-0.3535533905932737*m1rOther[14]*mnuOther+0.125*m0rOther[1]*uOther[12]*mnuOther+0.125*m0rOther[2]*uOther[11]*mnuOther+0.125*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(76,56) = 0.1964285714285714*m0rOther[9]*uOther[19]*mnuOther+0.1*m0rOther[8]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[19]*mnuOther-0.3162277660168379*m1rOther[19]*mnuOther+0.1*m0rOther[9]*uOther[18]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[18]*mnuOther-0.3162277660168379*m1rOther[18]*mnuOther+0.125*m0rOther[7]*uOther[17]*mnuOther+0.405*m0rOther[6]*uOther[16]*mnuOther+0.225*m0rOther[5]*uOther[15]*mnuOther+0.225*m0rOther[4]*uOther[14]*mnuOther+0.225*m0rOther[3]*uOther[13]*mnuOther+0.225*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[1]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[10]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(76,57) = 0.125*m0rOther[6]*uOther[17]*mnuOther+0.125*m0rOther[7]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[14]*mnuOther; 
  data->AEM_S(76,58) = 0.1*m0rOther[6]*uOther[19]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[18]*mnuOther+0.1*m0rOther[9]*uOther[16]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[16]*mnuOther-0.3162277660168379*m1rOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(76,59) = 0.1964285714285714*m0rOther[6]*uOther[19]*mnuOther+0.1*m0rOther[6]*uOther[18]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[16]*mnuOther+0.1*m0rOther[8]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[16]*mnuOther-0.3162277660168379*m1rOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(77,50) = 0.07985957062499249*m0rOther[7]*uOther[17]*mnuOther+0.125*m0rOther[0]*uOther[17]*mnuOther-0.3535533905932737*m1rOther[17]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[11]*mnuOther+0.125*m0rOther[7]*uOther[10]*mnuOther; 
  data->AEM_S(77,51) = 0.1964285714285714*m0rOther[1]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[12]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[11]*mnuOther-0.3162277660168379*m1rOther[11]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[10]*mnuOther; 
  data->AEM_S(77,52) = 0.125*m0rOther[2]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[14]*mnuOther+0.125*m0rOther[7]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[11]*mnuOther; 
  data->AEM_S(77,53) = 0.125*m0rOther[3]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[15]*mnuOther+0.125*m0rOther[7]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[11]*mnuOther; 
  data->AEM_S(77,54) = 0.1*m0rOther[4]*uOther[18]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[15]*mnuOther+0.1*m0rOther[8]*uOther[14]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[14]*mnuOther-0.3162277660168379*m1rOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(77,55) = 0.1*m0rOther[5]*uOther[19]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[16]*mnuOther+0.1*m0rOther[9]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[15]*mnuOther-0.3162277660168379*m1rOther[15]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(77,56) = 0.125*m0rOther[6]*uOther[17]*mnuOther+0.125*m0rOther[7]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[14]*mnuOther; 
  data->AEM_S(77,57) = 0.125*m0rOther[9]*uOther[19]*mnuOther+0.125*m0rOther[8]*uOther[18]*mnuOther+0.2678571428571428*m0rOther[7]*uOther[17]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[17]*mnuOther-0.2258769757263128*m1rOther[17]*mnuOther+0.125*m0rOther[6]*uOther[16]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[14]*mnuOther+0.125*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[2]*uOther[12]*mnuOther+0.1964285714285714*m0rOther[1]*uOther[11]*mnuOther+0.07985957062499249*m0rOther[7]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(77,58) = 0.125*m0rOther[7]*uOther[18]*mnuOther+0.125*m0rOther[8]*uOther[17]*mnuOther+0.1*m0rOther[4]*uOther[14]*mnuOther; 
  data->AEM_S(77,59) = 0.125*m0rOther[7]*uOther[19]*mnuOther+0.125*m0rOther[9]*uOther[17]*mnuOther+0.1*m0rOther[5]*uOther[15]*mnuOther; 
  data->AEM_S(78,50) = 0.07985957062499249*m0rOther[8]*uOther[18]*mnuOther+0.125*m0rOther[0]*uOther[18]*mnuOther-0.3535533905932737*m1rOther[18]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[8]*uOther[10]*mnuOther; 
  data->AEM_S(78,51) = 0.125*m0rOther[1]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[12]*mnuOther+0.125*m0rOther[8]*uOther[11]*mnuOther; 
  data->AEM_S(78,52) = 0.1964285714285714*m0rOther[2]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[13]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[12]*mnuOther-0.3162277660168379*m1rOther[12]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[10]*mnuOther; 
  data->AEM_S(78,53) = 0.125*m0rOther[3]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[16]*mnuOther+0.125*m0rOther[8]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[12]*mnuOther; 
  data->AEM_S(78,54) = 0.1964285714285714*m0rOther[4]*uOther[18]*mnuOther+0.1*m0rOther[4]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[14]*mnuOther+0.1*m0rOther[7]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[14]*mnuOther-0.3162277660168379*m1rOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[10]*mnuOther; 
  data->AEM_S(78,55) = 0.125*m0rOther[5]*uOther[18]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[16]*mnuOther+0.125*m0rOther[8]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[14]*mnuOther; 
  data->AEM_S(78,56) = 0.1*m0rOther[6]*uOther[19]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[18]*mnuOther+0.1*m0rOther[9]*uOther[16]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[16]*mnuOther-0.3162277660168379*m1rOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(78,57) = 0.125*m0rOther[7]*uOther[18]*mnuOther+0.125*m0rOther[8]*uOther[17]*mnuOther+0.1*m0rOther[4]*uOther[14]*mnuOther; 
  data->AEM_S(78,58) = 0.125*m0rOther[9]*uOther[19]*mnuOther+0.2678571428571428*m0rOther[8]*uOther[18]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[18]*mnuOther-0.2258769757263128*m1rOther[18]*mnuOther+0.125*m0rOther[7]*uOther[17]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[16]*mnuOther+0.125*m0rOther[5]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[14]*mnuOther+0.125*m0rOther[3]*uOther[13]*mnuOther+0.1964285714285714*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[1]*uOther[11]*mnuOther+0.07985957062499249*m0rOther[8]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
  data->AEM_S(78,59) = 0.125*m0rOther[8]*uOther[19]*mnuOther+0.125*m0rOther[9]*uOther[18]*mnuOther+0.1*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(79,50) = 0.07985957062499249*m0rOther[9]*uOther[19]*mnuOther+0.125*m0rOther[0]*uOther[19]*mnuOther-0.3535533905932737*m1rOther[19]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[9]*uOther[10]*mnuOther; 
  data->AEM_S(79,51) = 0.125*m0rOther[1]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[13]*mnuOther+0.125*m0rOther[9]*uOther[11]*mnuOther; 
  data->AEM_S(79,52) = 0.125*m0rOther[2]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[13]*mnuOther+0.125*m0rOther[9]*uOther[12]*mnuOther; 
  data->AEM_S(79,53) = 0.1964285714285714*m0rOther[3]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[15]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[13]*mnuOther-0.3162277660168379*m1rOther[13]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[10]*mnuOther; 
  data->AEM_S(79,54) = 0.125*m0rOther[4]*uOther[19]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[15]*mnuOther+0.125*m0rOther[9]*uOther[14]*mnuOther; 
  data->AEM_S(79,55) = 0.1964285714285714*m0rOther[5]*uOther[19]*mnuOther+0.1*m0rOther[5]*uOther[17]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[16]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[15]*mnuOther+0.1*m0rOther[7]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[15]*mnuOther-0.3162277660168379*m1rOther[15]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[11]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[10]*mnuOther; 
  data->AEM_S(79,56) = 0.1964285714285714*m0rOther[6]*uOther[19]*mnuOther+0.1*m0rOther[6]*uOther[18]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[16]*mnuOther+0.1*m0rOther[8]*uOther[16]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[16]*mnuOther-0.3162277660168379*m1rOther[16]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[15]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[14]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[13]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[12]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[10]*mnuOther; 
  data->AEM_S(79,57) = 0.125*m0rOther[7]*uOther[19]*mnuOther+0.125*m0rOther[9]*uOther[17]*mnuOther+0.1*m0rOther[5]*uOther[15]*mnuOther; 
  data->AEM_S(79,58) = 0.125*m0rOther[8]*uOther[19]*mnuOther+0.125*m0rOther[9]*uOther[18]*mnuOther+0.1*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(79,59) = 0.2678571428571428*m0rOther[9]*uOther[19]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[19]*mnuOther-0.2258769757263128*m1rOther[19]*mnuOther+0.125*m0rOther[8]*uOther[18]*mnuOther+0.125*m0rOther[7]*uOther[17]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[16]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[15]*mnuOther+0.125*m0rOther[4]*uOther[14]*mnuOther+0.1964285714285714*m0rOther[3]*uOther[13]*mnuOther+0.125*m0rOther[2]*uOther[12]*mnuOther+0.125*m0rOther[1]*uOther[11]*mnuOther+0.07985957062499249*m0rOther[9]*uOther[10]*mnuOther+0.125*m0rOther[0]*uOther[10]*mnuOther-0.3535533905932737*m1rOther[10]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfY-m0Self*m1OtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[10] = 0.3535533905932737*m0rOther[9]*m1rSelf[19]-0.3535533905932737*m0rSelf[9]*m1rOther[19]+0.3535533905932737*m0rOther[8]*m1rSelf[18]-0.3535533905932737*m0rSelf[8]*m1rOther[18]+0.3535533905932737*m0rOther[7]*m1rSelf[17]-0.3535533905932737*m0rSelf[7]*m1rOther[17]+0.3535533905932737*m0rOther[6]*m1rSelf[16]-0.3535533905932737*m0rSelf[6]*m1rOther[16]+0.3535533905932737*m0rOther[5]*m1rSelf[15]-0.3535533905932737*m0rSelf[5]*m1rOther[15]+0.3535533905932737*m0rOther[4]*m1rSelf[14]-0.3535533905932737*m0rSelf[4]*m1rOther[14]+0.3535533905932737*m0rOther[3]*m1rSelf[13]-0.3535533905932737*m0rSelf[3]*m1rOther[13]+0.3535533905932737*m0rOther[2]*m1rSelf[12]-0.3535533905932737*m0rSelf[2]*m1rOther[12]+0.3535533905932737*m0rOther[1]*m1rSelf[11]-0.3535533905932737*m0rSelf[1]*m1rOther[11]+0.3535533905932737*m0rOther[0]*m1rSelf[10]-0.3535533905932737*m0rSelf[0]*m1rOther[10]; 
  m1EffD[11] = 0.3162277660168379*m0rOther[1]*m1rSelf[17]-0.3162277660168379*m0rSelf[1]*m1rOther[17]+0.3535533905932737*m0rOther[3]*m1rSelf[15]-0.3535533905932737*m0rSelf[3]*m1rOther[15]+0.3535533905932737*m0rOther[2]*m1rSelf[14]-0.3535533905932737*m0rSelf[2]*m1rOther[14]+0.3535533905932737*m0rOther[5]*m1rSelf[13]-0.3535533905932737*m0rSelf[5]*m1rOther[13]+0.3535533905932737*m0rOther[4]*m1rSelf[12]-0.3535533905932737*m0rSelf[4]*m1rOther[12]+0.3162277660168379*m0rOther[7]*m1rSelf[11]+0.3535533905932737*m0rOther[0]*m1rSelf[11]-0.3162277660168379*m0rSelf[7]*m1rOther[11]-0.3535533905932737*m0rSelf[0]*m1rOther[11]+0.3535533905932737*m0rOther[1]*m1rSelf[10]-0.3535533905932737*m0rSelf[1]*m1rOther[10]; 
  m1EffD[12] = 0.3162277660168379*m0rOther[2]*m1rSelf[18]-0.3162277660168379*m0rSelf[2]*m1rOther[18]+0.3535533905932737*m0rOther[3]*m1rSelf[16]-0.3535533905932737*m0rSelf[3]*m1rOther[16]+0.3535533905932737*m0rOther[1]*m1rSelf[14]-0.3535533905932737*m0rSelf[1]*m1rOther[14]+0.3535533905932737*m0rOther[6]*m1rSelf[13]-0.3535533905932737*m0rSelf[6]*m1rOther[13]+0.3162277660168379*m0rOther[8]*m1rSelf[12]+0.3535533905932737*m0rOther[0]*m1rSelf[12]-0.3162277660168379*m0rSelf[8]*m1rOther[12]-0.3535533905932737*m0rSelf[0]*m1rOther[12]+0.3535533905932737*m0rOther[4]*m1rSelf[11]-0.3535533905932737*m0rSelf[4]*m1rOther[11]+0.3535533905932737*m0rOther[2]*m1rSelf[10]-0.3535533905932737*m0rSelf[2]*m1rOther[10]; 
  m1EffD[13] = 0.3162277660168379*m0rOther[3]*m1rSelf[19]-0.3162277660168379*m0rSelf[3]*m1rOther[19]+0.3535533905932737*m0rOther[2]*m1rSelf[16]-0.3535533905932737*m0rSelf[2]*m1rOther[16]+0.3535533905932737*m0rOther[1]*m1rSelf[15]-0.3535533905932737*m0rSelf[1]*m1rOther[15]+0.3162277660168379*m0rOther[9]*m1rSelf[13]+0.3535533905932737*m0rOther[0]*m1rSelf[13]-0.3162277660168379*m0rSelf[9]*m1rOther[13]-0.3535533905932737*m0rSelf[0]*m1rOther[13]+0.3535533905932737*m0rOther[6]*m1rSelf[12]-0.3535533905932737*m0rSelf[6]*m1rOther[12]+0.3535533905932737*m0rOther[5]*m1rSelf[11]-0.3535533905932737*m0rSelf[5]*m1rOther[11]+0.3535533905932737*m0rOther[3]*m1rSelf[10]-0.3535533905932737*m0rSelf[3]*m1rOther[10]; 
  m1EffD[14] = 0.3162277660168379*m0rOther[4]*m1rSelf[18]-0.3162277660168379*m0rSelf[4]*m1rOther[18]+0.3162277660168379*m0rOther[4]*m1rSelf[17]-0.3162277660168379*m0rSelf[4]*m1rOther[17]+0.3535533905932737*m0rOther[5]*m1rSelf[16]-0.3535533905932737*m0rSelf[5]*m1rOther[16]+0.3535533905932737*m0rOther[6]*m1rSelf[15]-0.3535533905932737*m0rSelf[6]*m1rOther[15]+0.3162277660168379*m0rOther[8]*m1rSelf[14]+0.3162277660168379*m0rOther[7]*m1rSelf[14]+0.3535533905932737*m0rOther[0]*m1rSelf[14]-0.3162277660168379*m0rSelf[8]*m1rOther[14]-0.3162277660168379*m0rSelf[7]*m1rOther[14]-0.3535533905932737*m0rSelf[0]*m1rOther[14]+0.3535533905932737*m0rOther[1]*m1rSelf[12]-0.3535533905932737*m0rSelf[1]*m1rOther[12]+0.3535533905932737*m0rOther[2]*m1rSelf[11]-0.3535533905932737*m0rSelf[2]*m1rOther[11]+0.3535533905932737*m0rOther[4]*m1rSelf[10]-0.3535533905932737*m0rSelf[4]*m1rOther[10]; 
  m1EffD[15] = 0.3162277660168379*m0rOther[5]*m1rSelf[19]-0.3162277660168379*m0rSelf[5]*m1rOther[19]+0.3162277660168379*m0rOther[5]*m1rSelf[17]-0.3162277660168379*m0rSelf[5]*m1rOther[17]+0.3535533905932737*m0rOther[4]*m1rSelf[16]-0.3535533905932737*m0rSelf[4]*m1rOther[16]+0.3162277660168379*m0rOther[9]*m1rSelf[15]+0.3162277660168379*m0rOther[7]*m1rSelf[15]+0.3535533905932737*m0rOther[0]*m1rSelf[15]-0.3162277660168379*m0rSelf[9]*m1rOther[15]-0.3162277660168379*m0rSelf[7]*m1rOther[15]-0.3535533905932737*m0rSelf[0]*m1rOther[15]+0.3535533905932737*m0rOther[6]*m1rSelf[14]-0.3535533905932737*m0rSelf[6]*m1rOther[14]+0.3535533905932737*m0rOther[1]*m1rSelf[13]-0.3535533905932737*m0rSelf[1]*m1rOther[13]+0.3535533905932737*m0rOther[3]*m1rSelf[11]-0.3535533905932737*m0rSelf[3]*m1rOther[11]+0.3535533905932737*m0rOther[5]*m1rSelf[10]-0.3535533905932737*m0rSelf[5]*m1rOther[10]; 
  m1EffD[16] = 0.3162277660168379*m0rOther[6]*m1rSelf[19]-0.3162277660168379*m0rSelf[6]*m1rOther[19]+0.3162277660168379*m0rOther[6]*m1rSelf[18]-0.3162277660168379*m0rSelf[6]*m1rOther[18]+0.3162277660168379*m0rOther[9]*m1rSelf[16]+0.3162277660168379*m0rOther[8]*m1rSelf[16]+0.3535533905932737*m0rOther[0]*m1rSelf[16]-0.3162277660168379*m0rSelf[9]*m1rOther[16]-0.3162277660168379*m0rSelf[8]*m1rOther[16]-0.3535533905932737*m0rSelf[0]*m1rOther[16]+0.3535533905932737*m0rOther[4]*m1rSelf[15]-0.3535533905932737*m0rSelf[4]*m1rOther[15]+0.3535533905932737*m0rOther[5]*m1rSelf[14]-0.3535533905932737*m0rSelf[5]*m1rOther[14]+0.3535533905932737*m0rOther[2]*m1rSelf[13]-0.3535533905932737*m0rSelf[2]*m1rOther[13]+0.3535533905932737*m0rOther[3]*m1rSelf[12]-0.3535533905932737*m0rSelf[3]*m1rOther[12]+0.3535533905932737*m0rOther[6]*m1rSelf[10]-0.3535533905932737*m0rSelf[6]*m1rOther[10]; 
  m1EffD[17] = 0.2258769757263128*m0rOther[7]*m1rSelf[17]+0.3535533905932737*m0rOther[0]*m1rSelf[17]-0.2258769757263128*m0rSelf[7]*m1rOther[17]-0.3535533905932737*m0rSelf[0]*m1rOther[17]+0.3162277660168379*m0rOther[5]*m1rSelf[15]-0.3162277660168379*m0rSelf[5]*m1rOther[15]+0.3162277660168379*m0rOther[4]*m1rSelf[14]-0.3162277660168379*m0rSelf[4]*m1rOther[14]+0.3162277660168379*m0rOther[1]*m1rSelf[11]-0.3162277660168379*m0rSelf[1]*m1rOther[11]+0.3535533905932737*m0rOther[7]*m1rSelf[10]-0.3535533905932737*m0rSelf[7]*m1rOther[10]; 
  m1EffD[18] = 0.2258769757263128*m0rOther[8]*m1rSelf[18]+0.3535533905932737*m0rOther[0]*m1rSelf[18]-0.2258769757263128*m0rSelf[8]*m1rOther[18]-0.3535533905932737*m0rSelf[0]*m1rOther[18]+0.3162277660168379*m0rOther[6]*m1rSelf[16]-0.3162277660168379*m0rSelf[6]*m1rOther[16]+0.3162277660168379*m0rOther[4]*m1rSelf[14]-0.3162277660168379*m0rSelf[4]*m1rOther[14]+0.3162277660168379*m0rOther[2]*m1rSelf[12]-0.3162277660168379*m0rSelf[2]*m1rOther[12]+0.3535533905932737*m0rOther[8]*m1rSelf[10]-0.3535533905932737*m0rSelf[8]*m1rOther[10]; 
  m1EffD[19] = 0.2258769757263128*m0rOther[9]*m1rSelf[19]+0.3535533905932737*m0rOther[0]*m1rSelf[19]-0.2258769757263128*m0rSelf[9]*m1rOther[19]-0.3535533905932737*m0rSelf[0]*m1rOther[19]+0.3162277660168379*m0rOther[6]*m1rSelf[16]-0.3162277660168379*m0rSelf[6]*m1rOther[16]+0.3162277660168379*m0rOther[5]*m1rSelf[15]-0.3162277660168379*m0rSelf[5]*m1rOther[15]+0.3162277660168379*m0rOther[3]*m1rSelf[13]-0.3162277660168379*m0rSelf[3]*m1rOther[13]+0.3535533905932737*m0rOther[9]*m1rSelf[10]-0.3535533905932737*m0rSelf[9]*m1rOther[10]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[10],m1EffD[11],m1EffD[12],m1EffD[13],m1EffD[14],m1EffD[15],m1EffD[16],m1EffD[17],m1EffD[18],m1EffD[19]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+10,10,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[10] += (-2.0*m1EffD[10]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[10]*mnuSelf-1.0*m1rOther[10]*mnuOther; 
  m1Relax[11] += (-2.0*m1EffD[11]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[11]*mnuSelf-1.0*m1rOther[11]*mnuOther; 
  m1Relax[12] += (-2.0*m1EffD[12]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[12]*mnuSelf-1.0*m1rOther[12]*mnuOther; 
  m1Relax[13] += (-2.0*m1EffD[13]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[13]*mnuSelf-1.0*m1rOther[13]*mnuOther; 
  m1Relax[14] += (-2.0*m1EffD[14]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[14]*mnuSelf-1.0*m1rOther[14]*mnuOther; 
  m1Relax[15] += (-2.0*m1EffD[15]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[15]*mnuSelf-1.0*m1rOther[15]*mnuOther; 
  m1Relax[16] += (-2.0*m1EffD[16]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[16]*mnuSelf-1.0*m1rOther[16]*mnuOther; 
  m1Relax[17] += (-2.0*m1EffD[17]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[17]*mnuSelf-1.0*m1rOther[17]*mnuOther; 
  m1Relax[18] += (-2.0*m1EffD[18]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[18]*mnuSelf-1.0*m1rOther[18]*mnuOther; 
  m1Relax[19] += (-2.0*m1EffD[19]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[19]*mnuSelf-1.0*m1rOther[19]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(60,20) = 0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(60,21) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(60,22) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(60,23) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(60,24) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(60,25) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(60,26) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(60,27) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(60,28) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(60,29) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(61,20) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(61,21) = 0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(61,22) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(61,23) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(61,24) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(61,25) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(61,27) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(62,20) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(62,21) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(62,22) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(62,23) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(62,24) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(62,26) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(62,28) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(63,20) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(63,21) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(63,22) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(63,23) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(63,25) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(63,26) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(63,29) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(64,20) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(64,21) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(64,22) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(64,24) = 0.3162277660168379*m0rSelf[8]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(64,25) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(64,26) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(64,27) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(64,28) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(65,20) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(65,21) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(65,23) = 0.3535533905932737*m0rSelf[1]*mnuSelf; 
  data->AEM_S(65,24) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(65,25) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(65,26) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(65,27) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(65,29) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(66,20) = 0.3535533905932737*m0rSelf[6]*mnuSelf; 
  data->AEM_S(66,22) = 0.3535533905932737*m0rSelf[3]*mnuSelf; 
  data->AEM_S(66,23) = 0.3535533905932737*m0rSelf[2]*mnuSelf; 
  data->AEM_S(66,24) = 0.3535533905932737*m0rSelf[5]*mnuSelf; 
  data->AEM_S(66,25) = 0.3535533905932737*m0rSelf[4]*mnuSelf; 
  data->AEM_S(66,26) = 0.3162277660168379*m0rSelf[9]*mnuSelf+0.3162277660168379*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(66,28) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(66,29) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(67,20) = 0.3535533905932737*m0rSelf[7]*mnuSelf; 
  data->AEM_S(67,21) = 0.3162277660168379*m0rSelf[1]*mnuSelf; 
  data->AEM_S(67,24) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(67,25) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(67,27) = 0.2258769757263128*m0rSelf[7]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(68,20) = 0.3535533905932737*m0rSelf[8]*mnuSelf; 
  data->AEM_S(68,22) = 0.3162277660168379*m0rSelf[2]*mnuSelf; 
  data->AEM_S(68,24) = 0.3162277660168379*m0rSelf[4]*mnuSelf; 
  data->AEM_S(68,26) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(68,28) = 0.2258769757263128*m0rSelf[8]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
  data->AEM_S(69,20) = 0.3535533905932737*m0rSelf[9]*mnuSelf; 
  data->AEM_S(69,23) = 0.3162277660168379*m0rSelf[3]*mnuSelf; 
  data->AEM_S(69,25) = 0.3162277660168379*m0rSelf[5]*mnuSelf; 
  data->AEM_S(69,26) = 0.3162277660168379*m0rSelf[6]*mnuSelf; 
  data->AEM_S(69,29) = 0.2258769757263128*m0rSelf[9]*mnuSelf+0.3535533905932737*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(60,30) = -0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(60,31) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(60,32) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(60,33) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(60,34) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(60,35) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(60,36) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(60,37) = -0.3535533905932737*cMSelf[27]*mnuSelf; 
  data->AEM_S(60,38) = -0.3535533905932737*cMSelf[28]*mnuSelf; 
  data->AEM_S(60,39) = -0.3535533905932737*cMSelf[29]*mnuSelf; 
  data->AEM_S(61,30) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(61,31) = (-0.3162277660168379*cMSelf[27]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(61,32) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(61,33) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(61,34) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(61,35) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(61,37) = -0.3162277660168379*cMSelf[21]*mnuSelf; 
  data->AEM_S(62,30) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(62,31) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(62,32) = (-0.3162277660168379*cMSelf[28]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(62,33) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(62,34) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(62,36) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(62,38) = -0.3162277660168379*cMSelf[22]*mnuSelf; 
  data->AEM_S(63,30) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(63,31) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(63,32) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(63,33) = (-0.3162277660168379*cMSelf[29]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(63,35) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(63,36) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(63,39) = -0.3162277660168379*cMSelf[23]*mnuSelf; 
  data->AEM_S(64,30) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(64,31) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(64,32) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(64,34) = (-0.3162277660168379*cMSelf[28]*mnuSelf)-0.3162277660168379*cMSelf[27]*mnuSelf-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(64,35) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(64,36) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(64,37) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(64,38) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(65,30) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(65,31) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(65,33) = -0.3535533905932737*cMSelf[21]*mnuSelf; 
  data->AEM_S(65,34) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(65,35) = (-0.3162277660168379*cMSelf[29]*mnuSelf)-0.3162277660168379*cMSelf[27]*mnuSelf-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(65,36) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(65,37) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(65,39) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(66,30) = -0.3535533905932737*cMSelf[26]*mnuSelf; 
  data->AEM_S(66,32) = -0.3535533905932737*cMSelf[23]*mnuSelf; 
  data->AEM_S(66,33) = -0.3535533905932737*cMSelf[22]*mnuSelf; 
  data->AEM_S(66,34) = -0.3535533905932737*cMSelf[25]*mnuSelf; 
  data->AEM_S(66,35) = -0.3535533905932737*cMSelf[24]*mnuSelf; 
  data->AEM_S(66,36) = (-0.3162277660168379*cMSelf[29]*mnuSelf)-0.3162277660168379*cMSelf[28]*mnuSelf-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(66,38) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(66,39) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(67,30) = -0.3535533905932737*cMSelf[27]*mnuSelf; 
  data->AEM_S(67,31) = -0.3162277660168379*cMSelf[21]*mnuSelf; 
  data->AEM_S(67,34) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(67,35) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(67,37) = (-0.2258769757263128*cMSelf[27]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(68,30) = -0.3535533905932737*cMSelf[28]*mnuSelf; 
  data->AEM_S(68,32) = -0.3162277660168379*cMSelf[22]*mnuSelf; 
  data->AEM_S(68,34) = -0.3162277660168379*cMSelf[24]*mnuSelf; 
  data->AEM_S(68,36) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(68,38) = (-0.2258769757263128*cMSelf[28]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
  data->AEM_S(69,30) = -0.3535533905932737*cMSelf[29]*mnuSelf; 
  data->AEM_S(69,33) = -0.3162277660168379*cMSelf[23]*mnuSelf; 
  data->AEM_S(69,35) = -0.3162277660168379*cMSelf[25]*mnuSelf; 
  data->AEM_S(69,36) = -0.3162277660168379*cMSelf[26]*mnuSelf; 
  data->AEM_S(69,39) = (-0.2258769757263128*cMSelf[29]*mnuSelf)-0.3535533905932737*cMSelf[20]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(60,60) = -0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(60,61) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(60,62) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(60,63) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(60,64) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(60,65) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(60,66) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(60,67) = -0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(60,68) = -0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(60,69) = -0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(61,60) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(61,61) = (-0.3162277660168379*m0rOther[7]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(61,62) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(61,63) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(61,64) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(61,65) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(61,67) = -0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(62,60) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(62,61) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(62,62) = (-0.3162277660168379*m0rOther[8]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(62,63) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(62,64) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(62,66) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(62,68) = -0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(63,60) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(63,61) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(63,62) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(63,63) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(63,65) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(63,66) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(63,69) = -0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(64,60) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(64,61) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(64,62) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(64,64) = (-0.3162277660168379*m0rOther[8]*mnuOther)-0.3162277660168379*m0rOther[7]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(64,65) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(64,66) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(64,67) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(64,68) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(65,60) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(65,61) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(65,63) = -0.3535533905932737*m0rOther[1]*mnuOther; 
  data->AEM_S(65,64) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(65,65) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3162277660168379*m0rOther[7]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(65,66) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(65,67) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(65,69) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(66,60) = -0.3535533905932737*m0rOther[6]*mnuOther; 
  data->AEM_S(66,62) = -0.3535533905932737*m0rOther[3]*mnuOther; 
  data->AEM_S(66,63) = -0.3535533905932737*m0rOther[2]*mnuOther; 
  data->AEM_S(66,64) = -0.3535533905932737*m0rOther[5]*mnuOther; 
  data->AEM_S(66,65) = -0.3535533905932737*m0rOther[4]*mnuOther; 
  data->AEM_S(66,66) = (-0.3162277660168379*m0rOther[9]*mnuOther)-0.3162277660168379*m0rOther[8]*mnuOther-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(66,68) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(66,69) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(67,60) = -0.3535533905932737*m0rOther[7]*mnuOther; 
  data->AEM_S(67,61) = -0.3162277660168379*m0rOther[1]*mnuOther; 
  data->AEM_S(67,64) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(67,65) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(67,67) = (-0.2258769757263128*m0rOther[7]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(68,60) = -0.3535533905932737*m0rOther[8]*mnuOther; 
  data->AEM_S(68,62) = -0.3162277660168379*m0rOther[2]*mnuOther; 
  data->AEM_S(68,64) = -0.3162277660168379*m0rOther[4]*mnuOther; 
  data->AEM_S(68,66) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(68,68) = (-0.2258769757263128*m0rOther[8]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
  data->AEM_S(69,60) = -0.3535533905932737*m0rOther[9]*mnuOther; 
  data->AEM_S(69,63) = -0.3162277660168379*m0rOther[3]*mnuOther; 
  data->AEM_S(69,65) = -0.3162277660168379*m0rOther[5]*mnuOther; 
  data->AEM_S(69,66) = -0.3162277660168379*m0rOther[6]*mnuOther; 
  data->AEM_S(69,69) = (-0.2258769757263128*m0rOther[9]*mnuOther)-0.3535533905932737*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(60,70) = 0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(60,71) = 0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(60,72) = 0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(60,73) = 0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(60,74) = 0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(60,75) = 0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(60,76) = 0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(60,77) = 0.3535533905932737*cMOther[27]*mnuOther; 
  data->AEM_S(60,78) = 0.3535533905932737*cMOther[28]*mnuOther; 
  data->AEM_S(60,79) = 0.3535533905932737*cMOther[29]*mnuOther; 
  data->AEM_S(61,70) = 0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(61,71) = 0.3162277660168379*cMOther[27]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(61,72) = 0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(61,73) = 0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(61,74) = 0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(61,75) = 0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(61,77) = 0.3162277660168379*cMOther[21]*mnuOther; 
  data->AEM_S(62,70) = 0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(62,71) = 0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(62,72) = 0.3162277660168379*cMOther[28]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(62,73) = 0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(62,74) = 0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(62,76) = 0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(62,78) = 0.3162277660168379*cMOther[22]*mnuOther; 
  data->AEM_S(63,70) = 0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(63,71) = 0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(63,72) = 0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(63,73) = 0.3162277660168379*cMOther[29]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(63,75) = 0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(63,76) = 0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(63,79) = 0.3162277660168379*cMOther[23]*mnuOther; 
  data->AEM_S(64,70) = 0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(64,71) = 0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(64,72) = 0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(64,74) = 0.3162277660168379*cMOther[28]*mnuOther+0.3162277660168379*cMOther[27]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(64,75) = 0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(64,76) = 0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(64,77) = 0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(64,78) = 0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(65,70) = 0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(65,71) = 0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(65,73) = 0.3535533905932737*cMOther[21]*mnuOther; 
  data->AEM_S(65,74) = 0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(65,75) = 0.3162277660168379*cMOther[29]*mnuOther+0.3162277660168379*cMOther[27]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(65,76) = 0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(65,77) = 0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(65,79) = 0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(66,70) = 0.3535533905932737*cMOther[26]*mnuOther; 
  data->AEM_S(66,72) = 0.3535533905932737*cMOther[23]*mnuOther; 
  data->AEM_S(66,73) = 0.3535533905932737*cMOther[22]*mnuOther; 
  data->AEM_S(66,74) = 0.3535533905932737*cMOther[25]*mnuOther; 
  data->AEM_S(66,75) = 0.3535533905932737*cMOther[24]*mnuOther; 
  data->AEM_S(66,76) = 0.3162277660168379*cMOther[29]*mnuOther+0.3162277660168379*cMOther[28]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(66,78) = 0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(66,79) = 0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(67,70) = 0.3535533905932737*cMOther[27]*mnuOther; 
  data->AEM_S(67,71) = 0.3162277660168379*cMOther[21]*mnuOther; 
  data->AEM_S(67,74) = 0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(67,75) = 0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(67,77) = 0.2258769757263128*cMOther[27]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(68,70) = 0.3535533905932737*cMOther[28]*mnuOther; 
  data->AEM_S(68,72) = 0.3162277660168379*cMOther[22]*mnuOther; 
  data->AEM_S(68,74) = 0.3162277660168379*cMOther[24]*mnuOther; 
  data->AEM_S(68,76) = 0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(68,78) = 0.2258769757263128*cMOther[28]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
  data->AEM_S(69,70) = 0.3535533905932737*cMOther[29]*mnuOther; 
  data->AEM_S(69,73) = 0.3162277660168379*cMOther[23]*mnuOther; 
  data->AEM_S(69,75) = 0.3162277660168379*cMOther[25]*mnuOther; 
  data->AEM_S(69,76) = 0.3162277660168379*cMOther[26]*mnuOther; 
  data->AEM_S(69,79) = 0.2258769757263128*cMOther[29]*mnuOther+0.3535533905932737*cMOther[20]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfZ-uSelfZ*m0Self) and uCrossSelfZ ... // 
  data->AEM_S(70,20) = (-0.125*m0rSelf[9]*uSelf[29]*mnuSelf)-0.125*m0rSelf[8]*uSelf[28]*mnuSelf-0.125*m0rSelf[7]*uSelf[27]*mnuSelf-0.125*m0rSelf[6]*uSelf[26]*mnuSelf-0.125*m0rSelf[5]*uSelf[25]*mnuSelf-0.125*m0rSelf[4]*uSelf[24]*mnuSelf-0.125*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[1]*uSelf[21]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(70,21) = (-0.1118033988749895*m0rSelf[1]*uSelf[27]*mnuSelf)-0.125*m0rSelf[3]*uSelf[25]*mnuSelf-0.125*m0rSelf[2]*uSelf[24]*mnuSelf-0.125*m0rSelf[5]*uSelf[23]*mnuSelf-0.125*m0rSelf[4]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[21]*mnuSelf-0.125*m0rSelf[0]*uSelf[21]*mnuSelf+0.3535533905932737*m1rSelf[21]*mnuSelf-0.125*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,22) = (-0.1118033988749895*m0rSelf[2]*uSelf[28]*mnuSelf)-0.125*m0rSelf[3]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[24]*mnuSelf-0.125*m0rSelf[6]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[22]*mnuSelf-0.125*m0rSelf[0]*uSelf[22]*mnuSelf+0.3535533905932737*m1rSelf[22]*mnuSelf-0.125*m0rSelf[4]*uSelf[21]*mnuSelf-0.125*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,23) = (-0.1118033988749895*m0rSelf[3]*uSelf[29]*mnuSelf)-0.125*m0rSelf[2]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[23]*mnuSelf-0.125*m0rSelf[0]*uSelf[23]*mnuSelf+0.3535533905932737*m1rSelf[23]*mnuSelf-0.125*m0rSelf[6]*uSelf[22]*mnuSelf-0.125*m0rSelf[5]*uSelf[21]*mnuSelf-0.125*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,24) = (-0.1118033988749895*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[27]*mnuSelf-0.125*m0rSelf[5]*uSelf[26]*mnuSelf-0.125*m0rSelf[6]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[24]*mnuSelf-0.125*m0rSelf[0]*uSelf[24]*mnuSelf+0.3535533905932737*m1rSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[22]*mnuSelf-0.125*m0rSelf[2]*uSelf[21]*mnuSelf-0.125*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,25) = (-0.1118033988749895*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[27]*mnuSelf-0.125*m0rSelf[4]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[25]*mnuSelf-0.125*m0rSelf[0]*uSelf[25]*mnuSelf+0.3535533905932737*m1rSelf[25]*mnuSelf-0.125*m0rSelf[6]*uSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[21]*mnuSelf-0.125*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,26) = (-0.1118033988749895*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[26]*mnuSelf-0.125*m0rSelf[0]*uSelf[26]*mnuSelf+0.3535533905932737*m1rSelf[26]*mnuSelf-0.125*m0rSelf[4]*uSelf[25]*mnuSelf-0.125*m0rSelf[5]*uSelf[24]*mnuSelf-0.125*m0rSelf[2]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,27) = (-0.07985957062499249*m0rSelf[7]*uSelf[27]*mnuSelf)-0.125*m0rSelf[0]*uSelf[27]*mnuSelf+0.3535533905932737*m1rSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[21]*mnuSelf-0.125*m0rSelf[7]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,28) = (-0.07985957062499249*m0rSelf[8]*uSelf[28]*mnuSelf)-0.125*m0rSelf[0]*uSelf[28]*mnuSelf+0.3535533905932737*m1rSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[8]*uSelf[20]*mnuSelf; 
  data->AEM_S(70,29) = (-0.07985957062499249*m0rSelf[9]*uSelf[29]*mnuSelf)-0.125*m0rSelf[0]*uSelf[29]*mnuSelf+0.3535533905932737*m1rSelf[29]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[9]*uSelf[20]*mnuSelf; 
  data->AEM_S(71,20) = (-0.1118033988749895*m0rSelf[1]*uSelf[27]*mnuSelf)-0.125*m0rSelf[3]*uSelf[25]*mnuSelf-0.125*m0rSelf[2]*uSelf[24]*mnuSelf-0.125*m0rSelf[5]*uSelf[23]*mnuSelf-0.125*m0rSelf[4]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[21]*mnuSelf-0.125*m0rSelf[0]*uSelf[21]*mnuSelf+0.3535533905932737*m1rSelf[21]*mnuSelf-0.125*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(71,21) = (-0.125*m0rSelf[9]*uSelf[29]*mnuSelf)-0.125*m0rSelf[8]*uSelf[28]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[27]*mnuSelf+0.3162277660168379*m1rSelf[27]*mnuSelf-0.125*m0rSelf[6]*uSelf[26]*mnuSelf-0.225*m0rSelf[5]*uSelf[25]*mnuSelf-0.225*m0rSelf[4]*uSelf[24]*mnuSelf-0.125*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[2]*uSelf[22]*mnuSelf-0.225*m0rSelf[1]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(71,22) = (-0.1118033988749895*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[27]*mnuSelf-0.125*m0rSelf[5]*uSelf[26]*mnuSelf-0.125*m0rSelf[6]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[24]*mnuSelf-0.125*m0rSelf[0]*uSelf[24]*mnuSelf+0.3535533905932737*m1rSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[22]*mnuSelf-0.125*m0rSelf[2]*uSelf[21]*mnuSelf-0.125*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(71,23) = (-0.1118033988749895*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[27]*mnuSelf-0.125*m0rSelf[4]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[25]*mnuSelf-0.125*m0rSelf[0]*uSelf[25]*mnuSelf+0.3535533905932737*m1rSelf[25]*mnuSelf-0.125*m0rSelf[6]*uSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[21]*mnuSelf-0.125*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(71,24) = (-0.1118033988749895*m0rSelf[2]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[27]*mnuSelf-0.125*m0rSelf[3]*uSelf[26]*mnuSelf-0.225*m0rSelf[1]*uSelf[24]*mnuSelf-0.125*m0rSelf[6]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[22]*mnuSelf-0.125*m0rSelf[0]*uSelf[22]*mnuSelf+0.3535533905932737*m1rSelf[22]*mnuSelf-0.225*m0rSelf[4]*uSelf[21]*mnuSelf-0.125*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(71,25) = (-0.1118033988749895*m0rSelf[3]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[27]*mnuSelf-0.125*m0rSelf[2]*uSelf[26]*mnuSelf-0.225*m0rSelf[1]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[23]*mnuSelf-0.125*m0rSelf[0]*uSelf[23]*mnuSelf+0.3535533905932737*m1rSelf[23]*mnuSelf-0.125*m0rSelf[6]*uSelf[22]*mnuSelf-0.225*m0rSelf[5]*uSelf[21]*mnuSelf-0.125*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(71,26) = (-0.125*m0rSelf[1]*uSelf[26]*mnuSelf)-0.125*m0rSelf[2]*uSelf[25]*mnuSelf-0.125*m0rSelf[3]*uSelf[24]*mnuSelf-0.125*m0rSelf[4]*uSelf[23]*mnuSelf-0.125*m0rSelf[5]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[21]*mnuSelf; 
  data->AEM_S(71,27) = (-0.1964285714285714*m0rSelf[1]*uSelf[27]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[22]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[21]*mnuSelf+0.3162277660168379*m1rSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(71,28) = (-0.125*m0rSelf[1]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[22]*mnuSelf-0.125*m0rSelf[8]*uSelf[21]*mnuSelf; 
  data->AEM_S(71,29) = (-0.125*m0rSelf[1]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[23]*mnuSelf-0.125*m0rSelf[9]*uSelf[21]*mnuSelf; 
  data->AEM_S(72,20) = (-0.1118033988749895*m0rSelf[2]*uSelf[28]*mnuSelf)-0.125*m0rSelf[3]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[24]*mnuSelf-0.125*m0rSelf[6]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[22]*mnuSelf-0.125*m0rSelf[0]*uSelf[22]*mnuSelf+0.3535533905932737*m1rSelf[22]*mnuSelf-0.125*m0rSelf[4]*uSelf[21]*mnuSelf-0.125*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(72,21) = (-0.1118033988749895*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[27]*mnuSelf-0.125*m0rSelf[5]*uSelf[26]*mnuSelf-0.125*m0rSelf[6]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[24]*mnuSelf-0.125*m0rSelf[0]*uSelf[24]*mnuSelf+0.3535533905932737*m1rSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[22]*mnuSelf-0.125*m0rSelf[2]*uSelf[21]*mnuSelf-0.125*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(72,22) = (-0.125*m0rSelf[9]*uSelf[29]*mnuSelf)-0.1964285714285714*m0rSelf[8]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[28]*mnuSelf+0.3162277660168379*m1rSelf[28]*mnuSelf-0.125*m0rSelf[7]*uSelf[27]*mnuSelf-0.225*m0rSelf[6]*uSelf[26]*mnuSelf-0.125*m0rSelf[5]*uSelf[25]*mnuSelf-0.225*m0rSelf[4]*uSelf[24]*mnuSelf-0.125*m0rSelf[3]*uSelf[23]*mnuSelf-0.225*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[1]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(72,23) = (-0.1118033988749895*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[26]*mnuSelf-0.125*m0rSelf[0]*uSelf[26]*mnuSelf+0.3535533905932737*m1rSelf[26]*mnuSelf-0.125*m0rSelf[4]*uSelf[25]*mnuSelf-0.125*m0rSelf[5]*uSelf[24]*mnuSelf-0.125*m0rSelf[2]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(72,24) = (-0.1118033988749895*m0rSelf[1]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[27]*mnuSelf-0.125*m0rSelf[3]*uSelf[25]*mnuSelf-0.225*m0rSelf[2]*uSelf[24]*mnuSelf-0.125*m0rSelf[5]*uSelf[23]*mnuSelf-0.225*m0rSelf[4]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[21]*mnuSelf-0.125*m0rSelf[0]*uSelf[21]*mnuSelf+0.3535533905932737*m1rSelf[21]*mnuSelf-0.125*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(72,25) = (-0.125*m0rSelf[1]*uSelf[26]*mnuSelf)-0.125*m0rSelf[2]*uSelf[25]*mnuSelf-0.125*m0rSelf[3]*uSelf[24]*mnuSelf-0.125*m0rSelf[4]*uSelf[23]*mnuSelf-0.125*m0rSelf[5]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[21]*mnuSelf; 
  data->AEM_S(72,26) = (-0.1118033988749895*m0rSelf[3]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[28]*mnuSelf-0.225*m0rSelf[2]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[23]*mnuSelf-0.125*m0rSelf[0]*uSelf[23]*mnuSelf+0.3535533905932737*m1rSelf[23]*mnuSelf-0.225*m0rSelf[6]*uSelf[22]*mnuSelf-0.125*m0rSelf[5]*uSelf[21]*mnuSelf-0.125*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(72,27) = (-0.125*m0rSelf[2]*uSelf[27]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[24]*mnuSelf-0.125*m0rSelf[7]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[21]*mnuSelf; 
  data->AEM_S(72,28) = (-0.1964285714285714*m0rSelf[2]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[23]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[22]*mnuSelf+0.3162277660168379*m1rSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(72,29) = (-0.125*m0rSelf[2]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[23]*mnuSelf-0.125*m0rSelf[9]*uSelf[22]*mnuSelf; 
  data->AEM_S(73,20) = (-0.1118033988749895*m0rSelf[3]*uSelf[29]*mnuSelf)-0.125*m0rSelf[2]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[23]*mnuSelf-0.125*m0rSelf[0]*uSelf[23]*mnuSelf+0.3535533905932737*m1rSelf[23]*mnuSelf-0.125*m0rSelf[6]*uSelf[22]*mnuSelf-0.125*m0rSelf[5]*uSelf[21]*mnuSelf-0.125*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(73,21) = (-0.1118033988749895*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[27]*mnuSelf-0.125*m0rSelf[4]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[25]*mnuSelf-0.125*m0rSelf[0]*uSelf[25]*mnuSelf+0.3535533905932737*m1rSelf[25]*mnuSelf-0.125*m0rSelf[6]*uSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[21]*mnuSelf-0.125*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(73,22) = (-0.1118033988749895*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[26]*mnuSelf-0.125*m0rSelf[0]*uSelf[26]*mnuSelf+0.3535533905932737*m1rSelf[26]*mnuSelf-0.125*m0rSelf[4]*uSelf[25]*mnuSelf-0.125*m0rSelf[5]*uSelf[24]*mnuSelf-0.125*m0rSelf[2]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(73,23) = (-0.1964285714285714*m0rSelf[9]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[0]*uSelf[29]*mnuSelf+0.3162277660168379*m1rSelf[29]*mnuSelf-0.125*m0rSelf[8]*uSelf[28]*mnuSelf-0.125*m0rSelf[7]*uSelf[27]*mnuSelf-0.225*m0rSelf[6]*uSelf[26]*mnuSelf-0.225*m0rSelf[5]*uSelf[25]*mnuSelf-0.125*m0rSelf[4]*uSelf[24]*mnuSelf-0.225*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[1]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(73,24) = (-0.125*m0rSelf[1]*uSelf[26]*mnuSelf)-0.125*m0rSelf[2]*uSelf[25]*mnuSelf-0.125*m0rSelf[3]*uSelf[24]*mnuSelf-0.125*m0rSelf[4]*uSelf[23]*mnuSelf-0.125*m0rSelf[5]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[21]*mnuSelf; 
  data->AEM_S(73,25) = (-0.1118033988749895*m0rSelf[1]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[27]*mnuSelf-0.225*m0rSelf[3]*uSelf[25]*mnuSelf-0.125*m0rSelf[2]*uSelf[24]*mnuSelf-0.225*m0rSelf[5]*uSelf[23]*mnuSelf-0.125*m0rSelf[4]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[21]*mnuSelf-0.125*m0rSelf[0]*uSelf[21]*mnuSelf+0.3535533905932737*m1rSelf[21]*mnuSelf-0.125*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(73,26) = (-0.1118033988749895*m0rSelf[2]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[28]*mnuSelf-0.225*m0rSelf[3]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[24]*mnuSelf-0.225*m0rSelf[6]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[22]*mnuSelf-0.125*m0rSelf[0]*uSelf[22]*mnuSelf+0.3535533905932737*m1rSelf[22]*mnuSelf-0.125*m0rSelf[4]*uSelf[21]*mnuSelf-0.125*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(73,27) = (-0.125*m0rSelf[3]*uSelf[27]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[25]*mnuSelf-0.125*m0rSelf[7]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[21]*mnuSelf; 
  data->AEM_S(73,28) = (-0.125*m0rSelf[3]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[26]*mnuSelf-0.125*m0rSelf[8]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[22]*mnuSelf; 
  data->AEM_S(73,29) = (-0.1964285714285714*m0rSelf[3]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[23]*mnuSelf+0.3162277660168379*m1rSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,20) = (-0.1118033988749895*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[27]*mnuSelf-0.125*m0rSelf[5]*uSelf[26]*mnuSelf-0.125*m0rSelf[6]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[24]*mnuSelf-0.125*m0rSelf[0]*uSelf[24]*mnuSelf+0.3535533905932737*m1rSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[22]*mnuSelf-0.125*m0rSelf[2]*uSelf[21]*mnuSelf-0.125*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,21) = (-0.1118033988749895*m0rSelf[2]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[27]*mnuSelf-0.125*m0rSelf[3]*uSelf[26]*mnuSelf-0.225*m0rSelf[1]*uSelf[24]*mnuSelf-0.125*m0rSelf[6]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[22]*mnuSelf-0.125*m0rSelf[0]*uSelf[22]*mnuSelf+0.3535533905932737*m1rSelf[22]*mnuSelf-0.225*m0rSelf[4]*uSelf[21]*mnuSelf-0.125*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,22) = (-0.1118033988749895*m0rSelf[1]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[27]*mnuSelf-0.125*m0rSelf[3]*uSelf[25]*mnuSelf-0.225*m0rSelf[2]*uSelf[24]*mnuSelf-0.125*m0rSelf[5]*uSelf[23]*mnuSelf-0.225*m0rSelf[4]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[21]*mnuSelf-0.125*m0rSelf[0]*uSelf[21]*mnuSelf+0.3535533905932737*m1rSelf[21]*mnuSelf-0.125*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,23) = (-0.125*m0rSelf[1]*uSelf[26]*mnuSelf)-0.125*m0rSelf[2]*uSelf[25]*mnuSelf-0.125*m0rSelf[3]*uSelf[24]*mnuSelf-0.125*m0rSelf[4]*uSelf[23]*mnuSelf-0.125*m0rSelf[5]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[21]*mnuSelf; 
  data->AEM_S(74,24) = (-0.125*m0rSelf[9]*uSelf[29]*mnuSelf)-0.1964285714285714*m0rSelf[8]*uSelf[28]*mnuSelf-0.1*m0rSelf[7]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[28]*mnuSelf+0.3162277660168379*m1rSelf[28]*mnuSelf-0.1*m0rSelf[8]*uSelf[27]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[27]*mnuSelf+0.3162277660168379*m1rSelf[27]*mnuSelf-0.225*m0rSelf[6]*uSelf[26]*mnuSelf-0.225*m0rSelf[5]*uSelf[25]*mnuSelf-0.405*m0rSelf[4]*uSelf[24]*mnuSelf-0.125*m0rSelf[3]*uSelf[23]*mnuSelf-0.225*m0rSelf[2]*uSelf[22]*mnuSelf-0.225*m0rSelf[1]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[20]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(74,25) = (-0.1118033988749895*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[26]*mnuSelf-0.125*m0rSelf[0]*uSelf[26]*mnuSelf+0.3535533905932737*m1rSelf[26]*mnuSelf-0.225*m0rSelf[4]*uSelf[25]*mnuSelf-0.225*m0rSelf[5]*uSelf[24]*mnuSelf-0.125*m0rSelf[2]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,26) = (-0.1118033988749895*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[27]*mnuSelf-0.225*m0rSelf[4]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[25]*mnuSelf-0.125*m0rSelf[0]*uSelf[25]*mnuSelf+0.3535533905932737*m1rSelf[25]*mnuSelf-0.225*m0rSelf[6]*uSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[21]*mnuSelf-0.125*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,27) = (-0.1*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1964285714285714*m0rSelf[4]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[25]*mnuSelf-0.1*m0rSelf[8]*uSelf[24]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[24]*mnuSelf+0.3162277660168379*m1rSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,28) = (-0.1964285714285714*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1*m0rSelf[4]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[24]*mnuSelf-0.1*m0rSelf[7]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[24]*mnuSelf+0.3162277660168379*m1rSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(74,29) = (-0.125*m0rSelf[4]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[25]*mnuSelf-0.125*m0rSelf[9]*uSelf[24]*mnuSelf; 
  data->AEM_S(75,20) = (-0.1118033988749895*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[27]*mnuSelf-0.125*m0rSelf[4]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[25]*mnuSelf-0.125*m0rSelf[0]*uSelf[25]*mnuSelf+0.3535533905932737*m1rSelf[25]*mnuSelf-0.125*m0rSelf[6]*uSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[21]*mnuSelf-0.125*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(75,21) = (-0.1118033988749895*m0rSelf[3]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[27]*mnuSelf-0.125*m0rSelf[2]*uSelf[26]*mnuSelf-0.225*m0rSelf[1]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[23]*mnuSelf-0.125*m0rSelf[0]*uSelf[23]*mnuSelf+0.3535533905932737*m1rSelf[23]*mnuSelf-0.125*m0rSelf[6]*uSelf[22]*mnuSelf-0.225*m0rSelf[5]*uSelf[21]*mnuSelf-0.125*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(75,22) = (-0.125*m0rSelf[1]*uSelf[26]*mnuSelf)-0.125*m0rSelf[2]*uSelf[25]*mnuSelf-0.125*m0rSelf[3]*uSelf[24]*mnuSelf-0.125*m0rSelf[4]*uSelf[23]*mnuSelf-0.125*m0rSelf[5]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[21]*mnuSelf; 
  data->AEM_S(75,23) = (-0.1118033988749895*m0rSelf[1]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[27]*mnuSelf-0.225*m0rSelf[3]*uSelf[25]*mnuSelf-0.125*m0rSelf[2]*uSelf[24]*mnuSelf-0.225*m0rSelf[5]*uSelf[23]*mnuSelf-0.125*m0rSelf[4]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[21]*mnuSelf-0.125*m0rSelf[0]*uSelf[21]*mnuSelf+0.3535533905932737*m1rSelf[21]*mnuSelf-0.125*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(75,24) = (-0.1118033988749895*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[26]*mnuSelf-0.125*m0rSelf[0]*uSelf[26]*mnuSelf+0.3535533905932737*m1rSelf[26]*mnuSelf-0.225*m0rSelf[4]*uSelf[25]*mnuSelf-0.225*m0rSelf[5]*uSelf[24]*mnuSelf-0.125*m0rSelf[2]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(75,25) = (-0.1964285714285714*m0rSelf[9]*uSelf[29]*mnuSelf)-0.1*m0rSelf[7]*uSelf[29]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[29]*mnuSelf+0.3162277660168379*m1rSelf[29]*mnuSelf-0.125*m0rSelf[8]*uSelf[28]*mnuSelf-0.1*m0rSelf[9]*uSelf[27]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[27]*mnuSelf+0.3162277660168379*m1rSelf[27]*mnuSelf-0.225*m0rSelf[6]*uSelf[26]*mnuSelf-0.405*m0rSelf[5]*uSelf[25]*mnuSelf-0.225*m0rSelf[4]*uSelf[24]*mnuSelf-0.225*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[2]*uSelf[22]*mnuSelf-0.225*m0rSelf[1]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[20]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(75,26) = (-0.1118033988749895*m0rSelf[4]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[27]*mnuSelf-0.225*m0rSelf[5]*uSelf[26]*mnuSelf-0.225*m0rSelf[6]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[24]*mnuSelf-0.125*m0rSelf[0]*uSelf[24]*mnuSelf+0.3535533905932737*m1rSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[22]*mnuSelf-0.125*m0rSelf[2]*uSelf[21]*mnuSelf-0.125*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(75,27) = (-0.1*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1964285714285714*m0rSelf[5]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[26]*mnuSelf-0.1*m0rSelf[9]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[25]*mnuSelf+0.3162277660168379*m1rSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(75,28) = (-0.125*m0rSelf[5]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[26]*mnuSelf-0.125*m0rSelf[8]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[24]*mnuSelf; 
  data->AEM_S(75,29) = (-0.1964285714285714*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1*m0rSelf[5]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[26]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[25]*mnuSelf-0.1*m0rSelf[7]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[25]*mnuSelf+0.3162277660168379*m1rSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(76,20) = (-0.1118033988749895*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[6]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[26]*mnuSelf-0.125*m0rSelf[0]*uSelf[26]*mnuSelf+0.3535533905932737*m1rSelf[26]*mnuSelf-0.125*m0rSelf[4]*uSelf[25]*mnuSelf-0.125*m0rSelf[5]*uSelf[24]*mnuSelf-0.125*m0rSelf[2]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(76,21) = (-0.125*m0rSelf[1]*uSelf[26]*mnuSelf)-0.125*m0rSelf[2]*uSelf[25]*mnuSelf-0.125*m0rSelf[3]*uSelf[24]*mnuSelf-0.125*m0rSelf[4]*uSelf[23]*mnuSelf-0.125*m0rSelf[5]*uSelf[22]*mnuSelf-0.125*m0rSelf[6]*uSelf[21]*mnuSelf; 
  data->AEM_S(76,22) = (-0.1118033988749895*m0rSelf[3]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[28]*mnuSelf-0.225*m0rSelf[2]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[23]*mnuSelf-0.125*m0rSelf[0]*uSelf[23]*mnuSelf+0.3535533905932737*m1rSelf[23]*mnuSelf-0.225*m0rSelf[6]*uSelf[22]*mnuSelf-0.125*m0rSelf[5]*uSelf[21]*mnuSelf-0.125*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(76,23) = (-0.1118033988749895*m0rSelf[2]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[28]*mnuSelf-0.225*m0rSelf[3]*uSelf[26]*mnuSelf-0.125*m0rSelf[1]*uSelf[24]*mnuSelf-0.225*m0rSelf[6]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[22]*mnuSelf-0.125*m0rSelf[0]*uSelf[22]*mnuSelf+0.3535533905932737*m1rSelf[22]*mnuSelf-0.125*m0rSelf[4]*uSelf[21]*mnuSelf-0.125*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(76,24) = (-0.1118033988749895*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[27]*mnuSelf-0.225*m0rSelf[4]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[25]*mnuSelf-0.125*m0rSelf[0]*uSelf[25]*mnuSelf+0.3535533905932737*m1rSelf[25]*mnuSelf-0.225*m0rSelf[6]*uSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[23]*mnuSelf-0.125*m0rSelf[3]*uSelf[21]*mnuSelf-0.125*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(76,25) = (-0.1118033988749895*m0rSelf[4]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[27]*mnuSelf-0.225*m0rSelf[5]*uSelf[26]*mnuSelf-0.225*m0rSelf[6]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[7]*uSelf[24]*mnuSelf-0.125*m0rSelf[0]*uSelf[24]*mnuSelf+0.3535533905932737*m1rSelf[24]*mnuSelf-0.125*m0rSelf[1]*uSelf[22]*mnuSelf-0.125*m0rSelf[2]*uSelf[21]*mnuSelf-0.125*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(76,26) = (-0.1964285714285714*m0rSelf[9]*uSelf[29]*mnuSelf)-0.1*m0rSelf[8]*uSelf[29]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[29]*mnuSelf+0.3162277660168379*m1rSelf[29]*mnuSelf-0.1*m0rSelf[9]*uSelf[28]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[28]*mnuSelf+0.3162277660168379*m1rSelf[28]*mnuSelf-0.125*m0rSelf[7]*uSelf[27]*mnuSelf-0.405*m0rSelf[6]*uSelf[26]*mnuSelf-0.225*m0rSelf[5]*uSelf[25]*mnuSelf-0.225*m0rSelf[4]*uSelf[24]*mnuSelf-0.225*m0rSelf[3]*uSelf[23]*mnuSelf-0.225*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[1]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[9]*uSelf[20]*mnuSelf-0.1118033988749895*m0rSelf[8]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(76,27) = (-0.125*m0rSelf[6]*uSelf[27]*mnuSelf)-0.125*m0rSelf[7]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[24]*mnuSelf; 
  data->AEM_S(76,28) = (-0.1*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1964285714285714*m0rSelf[6]*uSelf[28]*mnuSelf-0.1*m0rSelf[9]*uSelf[26]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[26]*mnuSelf+0.3162277660168379*m1rSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(76,29) = (-0.1964285714285714*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1*m0rSelf[6]*uSelf[28]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[26]*mnuSelf-0.1*m0rSelf[8]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[26]*mnuSelf+0.3162277660168379*m1rSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(77,20) = (-0.07985957062499249*m0rSelf[7]*uSelf[27]*mnuSelf)-0.125*m0rSelf[0]*uSelf[27]*mnuSelf+0.3535533905932737*m1rSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[21]*mnuSelf-0.125*m0rSelf[7]*uSelf[20]*mnuSelf; 
  data->AEM_S(77,21) = (-0.1964285714285714*m0rSelf[1]*uSelf[27]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[22]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[21]*mnuSelf+0.3162277660168379*m1rSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[20]*mnuSelf; 
  data->AEM_S(77,22) = (-0.125*m0rSelf[2]*uSelf[27]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[24]*mnuSelf-0.125*m0rSelf[7]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[21]*mnuSelf; 
  data->AEM_S(77,23) = (-0.125*m0rSelf[3]*uSelf[27]*mnuSelf)-0.1118033988749895*m0rSelf[1]*uSelf[25]*mnuSelf-0.125*m0rSelf[7]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[21]*mnuSelf; 
  data->AEM_S(77,24) = (-0.1*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1964285714285714*m0rSelf[4]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[25]*mnuSelf-0.1*m0rSelf[8]*uSelf[24]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[24]*mnuSelf+0.3162277660168379*m1rSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(77,25) = (-0.1*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1964285714285714*m0rSelf[5]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[26]*mnuSelf-0.1*m0rSelf[9]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[7]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[25]*mnuSelf+0.3162277660168379*m1rSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(77,26) = (-0.125*m0rSelf[6]*uSelf[27]*mnuSelf)-0.125*m0rSelf[7]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[24]*mnuSelf; 
  data->AEM_S(77,27) = (-0.125*m0rSelf[9]*uSelf[29]*mnuSelf)-0.125*m0rSelf[8]*uSelf[28]*mnuSelf-0.2678571428571428*m0rSelf[7]*uSelf[27]*mnuSelf-0.07985957062499249*m0rSelf[0]*uSelf[27]*mnuSelf+0.2258769757263128*m1rSelf[27]*mnuSelf-0.125*m0rSelf[6]*uSelf[26]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[24]*mnuSelf-0.125*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[2]*uSelf[22]*mnuSelf-0.1964285714285714*m0rSelf[1]*uSelf[21]*mnuSelf-0.07985957062499249*m0rSelf[7]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(77,28) = (-0.125*m0rSelf[7]*uSelf[28]*mnuSelf)-0.125*m0rSelf[8]*uSelf[27]*mnuSelf-0.1*m0rSelf[4]*uSelf[24]*mnuSelf; 
  data->AEM_S(77,29) = (-0.125*m0rSelf[7]*uSelf[29]*mnuSelf)-0.125*m0rSelf[9]*uSelf[27]*mnuSelf-0.1*m0rSelf[5]*uSelf[25]*mnuSelf; 
  data->AEM_S(78,20) = (-0.07985957062499249*m0rSelf[8]*uSelf[28]*mnuSelf)-0.125*m0rSelf[0]*uSelf[28]*mnuSelf+0.3535533905932737*m1rSelf[28]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[8]*uSelf[20]*mnuSelf; 
  data->AEM_S(78,21) = (-0.125*m0rSelf[1]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[22]*mnuSelf-0.125*m0rSelf[8]*uSelf[21]*mnuSelf; 
  data->AEM_S(78,22) = (-0.1964285714285714*m0rSelf[2]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[23]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[22]*mnuSelf+0.3162277660168379*m1rSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[20]*mnuSelf; 
  data->AEM_S(78,23) = (-0.125*m0rSelf[3]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[26]*mnuSelf-0.125*m0rSelf[8]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[22]*mnuSelf; 
  data->AEM_S(78,24) = (-0.1964285714285714*m0rSelf[4]*uSelf[28]*mnuSelf)-0.1*m0rSelf[4]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[24]*mnuSelf-0.1*m0rSelf[7]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[24]*mnuSelf+0.3162277660168379*m1rSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[20]*mnuSelf; 
  data->AEM_S(78,25) = (-0.125*m0rSelf[5]*uSelf[28]*mnuSelf)-0.1118033988749895*m0rSelf[4]*uSelf[26]*mnuSelf-0.125*m0rSelf[8]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[24]*mnuSelf; 
  data->AEM_S(78,26) = (-0.1*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1964285714285714*m0rSelf[6]*uSelf[28]*mnuSelf-0.1*m0rSelf[9]*uSelf[26]*mnuSelf-0.1964285714285714*m0rSelf[8]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[26]*mnuSelf+0.3162277660168379*m1rSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(78,27) = (-0.125*m0rSelf[7]*uSelf[28]*mnuSelf)-0.125*m0rSelf[8]*uSelf[27]*mnuSelf-0.1*m0rSelf[4]*uSelf[24]*mnuSelf; 
  data->AEM_S(78,28) = (-0.125*m0rSelf[9]*uSelf[29]*mnuSelf)-0.2678571428571428*m0rSelf[8]*uSelf[28]*mnuSelf-0.07985957062499249*m0rSelf[0]*uSelf[28]*mnuSelf+0.2258769757263128*m1rSelf[28]*mnuSelf-0.125*m0rSelf[7]*uSelf[27]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[26]*mnuSelf-0.125*m0rSelf[5]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[4]*uSelf[24]*mnuSelf-0.125*m0rSelf[3]*uSelf[23]*mnuSelf-0.1964285714285714*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[1]*uSelf[21]*mnuSelf-0.07985957062499249*m0rSelf[8]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
  data->AEM_S(78,29) = (-0.125*m0rSelf[8]*uSelf[29]*mnuSelf)-0.125*m0rSelf[9]*uSelf[28]*mnuSelf-0.1*m0rSelf[6]*uSelf[26]*mnuSelf; 
  data->AEM_S(79,20) = (-0.07985957062499249*m0rSelf[9]*uSelf[29]*mnuSelf)-0.125*m0rSelf[0]*uSelf[29]*mnuSelf+0.3535533905932737*m1rSelf[29]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[9]*uSelf[20]*mnuSelf; 
  data->AEM_S(79,21) = (-0.125*m0rSelf[1]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[23]*mnuSelf-0.125*m0rSelf[9]*uSelf[21]*mnuSelf; 
  data->AEM_S(79,22) = (-0.125*m0rSelf[2]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[3]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[23]*mnuSelf-0.125*m0rSelf[9]*uSelf[22]*mnuSelf; 
  data->AEM_S(79,23) = (-0.1964285714285714*m0rSelf[3]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[2]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[25]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[23]*mnuSelf+0.3162277660168379*m1rSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[20]*mnuSelf; 
  data->AEM_S(79,24) = (-0.125*m0rSelf[4]*uSelf[29]*mnuSelf)-0.1118033988749895*m0rSelf[5]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[25]*mnuSelf-0.125*m0rSelf[9]*uSelf[24]*mnuSelf; 
  data->AEM_S(79,25) = (-0.1964285714285714*m0rSelf[5]*uSelf[29]*mnuSelf)-0.1*m0rSelf[5]*uSelf[27]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[26]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[25]*mnuSelf-0.1*m0rSelf[7]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[25]*mnuSelf+0.3162277660168379*m1rSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[1]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[21]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[20]*mnuSelf; 
  data->AEM_S(79,26) = (-0.1964285714285714*m0rSelf[6]*uSelf[29]*mnuSelf)-0.1*m0rSelf[6]*uSelf[28]*mnuSelf-0.1964285714285714*m0rSelf[9]*uSelf[26]*mnuSelf-0.1*m0rSelf[8]*uSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[0]*uSelf[26]*mnuSelf+0.3162277660168379*m1rSelf[26]*mnuSelf-0.1118033988749895*m0rSelf[4]*uSelf[25]*mnuSelf-0.1118033988749895*m0rSelf[5]*uSelf[24]*mnuSelf-0.1118033988749895*m0rSelf[2]*uSelf[23]*mnuSelf-0.1118033988749895*m0rSelf[3]*uSelf[22]*mnuSelf-0.1118033988749895*m0rSelf[6]*uSelf[20]*mnuSelf; 
  data->AEM_S(79,27) = (-0.125*m0rSelf[7]*uSelf[29]*mnuSelf)-0.125*m0rSelf[9]*uSelf[27]*mnuSelf-0.1*m0rSelf[5]*uSelf[25]*mnuSelf; 
  data->AEM_S(79,28) = (-0.125*m0rSelf[8]*uSelf[29]*mnuSelf)-0.125*m0rSelf[9]*uSelf[28]*mnuSelf-0.1*m0rSelf[6]*uSelf[26]*mnuSelf; 
  data->AEM_S(79,29) = (-0.2678571428571428*m0rSelf[9]*uSelf[29]*mnuSelf)-0.07985957062499249*m0rSelf[0]*uSelf[29]*mnuSelf+0.2258769757263128*m1rSelf[29]*mnuSelf-0.125*m0rSelf[8]*uSelf[28]*mnuSelf-0.125*m0rSelf[7]*uSelf[27]*mnuSelf-0.1964285714285714*m0rSelf[6]*uSelf[26]*mnuSelf-0.1964285714285714*m0rSelf[5]*uSelf[25]*mnuSelf-0.125*m0rSelf[4]*uSelf[24]*mnuSelf-0.1964285714285714*m0rSelf[3]*uSelf[23]*mnuSelf-0.125*m0rSelf[2]*uSelf[22]*mnuSelf-0.125*m0rSelf[1]*uSelf[21]*mnuSelf-0.07985957062499249*m0rSelf[9]*uSelf[20]*mnuSelf-0.125*m0rSelf[0]*uSelf[20]*mnuSelf+0.3535533905932737*m1rSelf[20]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherZ-uOtherZ*m0Other) and uCrossOtherZ ... // 
  data->AEM_S(70,60) = 0.125*m0rOther[9]*uOther[29]*mnuOther+0.125*m0rOther[8]*uOther[28]*mnuOther+0.125*m0rOther[7]*uOther[27]*mnuOther+0.125*m0rOther[6]*uOther[26]*mnuOther+0.125*m0rOther[5]*uOther[25]*mnuOther+0.125*m0rOther[4]*uOther[24]*mnuOther+0.125*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[1]*uOther[21]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(70,61) = 0.1118033988749895*m0rOther[1]*uOther[27]*mnuOther+0.125*m0rOther[3]*uOther[25]*mnuOther+0.125*m0rOther[2]*uOther[24]*mnuOther+0.125*m0rOther[5]*uOther[23]*mnuOther+0.125*m0rOther[4]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[21]*mnuOther+0.125*m0rOther[0]*uOther[21]*mnuOther-0.3535533905932737*m1rOther[21]*mnuOther+0.125*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(70,62) = 0.1118033988749895*m0rOther[2]*uOther[28]*mnuOther+0.125*m0rOther[3]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[24]*mnuOther+0.125*m0rOther[6]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[22]*mnuOther+0.125*m0rOther[0]*uOther[22]*mnuOther-0.3535533905932737*m1rOther[22]*mnuOther+0.125*m0rOther[4]*uOther[21]*mnuOther+0.125*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(70,63) = 0.1118033988749895*m0rOther[3]*uOther[29]*mnuOther+0.125*m0rOther[2]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[23]*mnuOther+0.125*m0rOther[0]*uOther[23]*mnuOther-0.3535533905932737*m1rOther[23]*mnuOther+0.125*m0rOther[6]*uOther[22]*mnuOther+0.125*m0rOther[5]*uOther[21]*mnuOther+0.125*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(70,64) = 0.1118033988749895*m0rOther[4]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[27]*mnuOther+0.125*m0rOther[5]*uOther[26]*mnuOther+0.125*m0rOther[6]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[24]*mnuOther+0.125*m0rOther[0]*uOther[24]*mnuOther-0.3535533905932737*m1rOther[24]*mnuOther+0.125*m0rOther[1]*uOther[22]*mnuOther+0.125*m0rOther[2]*uOther[21]*mnuOther+0.125*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(70,65) = 0.1118033988749895*m0rOther[5]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[27]*mnuOther+0.125*m0rOther[4]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[25]*mnuOther+0.125*m0rOther[0]*uOther[25]*mnuOther-0.3535533905932737*m1rOther[25]*mnuOther+0.125*m0rOther[6]*uOther[24]*mnuOther+0.125*m0rOther[1]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[21]*mnuOther+0.125*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(70,66) = 0.1118033988749895*m0rOther[6]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[26]*mnuOther+0.125*m0rOther[0]*uOther[26]*mnuOther-0.3535533905932737*m1rOther[26]*mnuOther+0.125*m0rOther[4]*uOther[25]*mnuOther+0.125*m0rOther[5]*uOther[24]*mnuOther+0.125*m0rOther[2]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(70,67) = 0.07985957062499249*m0rOther[7]*uOther[27]*mnuOther+0.125*m0rOther[0]*uOther[27]*mnuOther-0.3535533905932737*m1rOther[27]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[21]*mnuOther+0.125*m0rOther[7]*uOther[20]*mnuOther; 
  data->AEM_S(70,68) = 0.07985957062499249*m0rOther[8]*uOther[28]*mnuOther+0.125*m0rOther[0]*uOther[28]*mnuOther-0.3535533905932737*m1rOther[28]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[8]*uOther[20]*mnuOther; 
  data->AEM_S(70,69) = 0.07985957062499249*m0rOther[9]*uOther[29]*mnuOther+0.125*m0rOther[0]*uOther[29]*mnuOther-0.3535533905932737*m1rOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[9]*uOther[20]*mnuOther; 
  data->AEM_S(71,60) = 0.1118033988749895*m0rOther[1]*uOther[27]*mnuOther+0.125*m0rOther[3]*uOther[25]*mnuOther+0.125*m0rOther[2]*uOther[24]*mnuOther+0.125*m0rOther[5]*uOther[23]*mnuOther+0.125*m0rOther[4]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[21]*mnuOther+0.125*m0rOther[0]*uOther[21]*mnuOther-0.3535533905932737*m1rOther[21]*mnuOther+0.125*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(71,61) = 0.125*m0rOther[9]*uOther[29]*mnuOther+0.125*m0rOther[8]*uOther[28]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[27]*mnuOther-0.3162277660168379*m1rOther[27]*mnuOther+0.125*m0rOther[6]*uOther[26]*mnuOther+0.225*m0rOther[5]*uOther[25]*mnuOther+0.225*m0rOther[4]*uOther[24]*mnuOther+0.125*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[2]*uOther[22]*mnuOther+0.225*m0rOther[1]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(71,62) = 0.1118033988749895*m0rOther[4]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[27]*mnuOther+0.125*m0rOther[5]*uOther[26]*mnuOther+0.125*m0rOther[6]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[24]*mnuOther+0.125*m0rOther[0]*uOther[24]*mnuOther-0.3535533905932737*m1rOther[24]*mnuOther+0.125*m0rOther[1]*uOther[22]*mnuOther+0.125*m0rOther[2]*uOther[21]*mnuOther+0.125*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(71,63) = 0.1118033988749895*m0rOther[5]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[27]*mnuOther+0.125*m0rOther[4]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[25]*mnuOther+0.125*m0rOther[0]*uOther[25]*mnuOther-0.3535533905932737*m1rOther[25]*mnuOther+0.125*m0rOther[6]*uOther[24]*mnuOther+0.125*m0rOther[1]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[21]*mnuOther+0.125*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(71,64) = 0.1118033988749895*m0rOther[2]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[27]*mnuOther+0.125*m0rOther[3]*uOther[26]*mnuOther+0.225*m0rOther[1]*uOther[24]*mnuOther+0.125*m0rOther[6]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[22]*mnuOther+0.125*m0rOther[0]*uOther[22]*mnuOther-0.3535533905932737*m1rOther[22]*mnuOther+0.225*m0rOther[4]*uOther[21]*mnuOther+0.125*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(71,65) = 0.1118033988749895*m0rOther[3]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[27]*mnuOther+0.125*m0rOther[2]*uOther[26]*mnuOther+0.225*m0rOther[1]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[23]*mnuOther+0.125*m0rOther[0]*uOther[23]*mnuOther-0.3535533905932737*m1rOther[23]*mnuOther+0.125*m0rOther[6]*uOther[22]*mnuOther+0.225*m0rOther[5]*uOther[21]*mnuOther+0.125*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(71,66) = 0.125*m0rOther[1]*uOther[26]*mnuOther+0.125*m0rOther[2]*uOther[25]*mnuOther+0.125*m0rOther[3]*uOther[24]*mnuOther+0.125*m0rOther[4]*uOther[23]*mnuOther+0.125*m0rOther[5]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[21]*mnuOther; 
  data->AEM_S(71,67) = 0.1964285714285714*m0rOther[1]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[22]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[21]*mnuOther-0.3162277660168379*m1rOther[21]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(71,68) = 0.125*m0rOther[1]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[22]*mnuOther+0.125*m0rOther[8]*uOther[21]*mnuOther; 
  data->AEM_S(71,69) = 0.125*m0rOther[1]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[23]*mnuOther+0.125*m0rOther[9]*uOther[21]*mnuOther; 
  data->AEM_S(72,60) = 0.1118033988749895*m0rOther[2]*uOther[28]*mnuOther+0.125*m0rOther[3]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[24]*mnuOther+0.125*m0rOther[6]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[22]*mnuOther+0.125*m0rOther[0]*uOther[22]*mnuOther-0.3535533905932737*m1rOther[22]*mnuOther+0.125*m0rOther[4]*uOther[21]*mnuOther+0.125*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(72,61) = 0.1118033988749895*m0rOther[4]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[27]*mnuOther+0.125*m0rOther[5]*uOther[26]*mnuOther+0.125*m0rOther[6]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[24]*mnuOther+0.125*m0rOther[0]*uOther[24]*mnuOther-0.3535533905932737*m1rOther[24]*mnuOther+0.125*m0rOther[1]*uOther[22]*mnuOther+0.125*m0rOther[2]*uOther[21]*mnuOther+0.125*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(72,62) = 0.125*m0rOther[9]*uOther[29]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[28]*mnuOther-0.3162277660168379*m1rOther[28]*mnuOther+0.125*m0rOther[7]*uOther[27]*mnuOther+0.225*m0rOther[6]*uOther[26]*mnuOther+0.125*m0rOther[5]*uOther[25]*mnuOther+0.225*m0rOther[4]*uOther[24]*mnuOther+0.125*m0rOther[3]*uOther[23]*mnuOther+0.225*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[1]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(72,63) = 0.1118033988749895*m0rOther[6]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[26]*mnuOther+0.125*m0rOther[0]*uOther[26]*mnuOther-0.3535533905932737*m1rOther[26]*mnuOther+0.125*m0rOther[4]*uOther[25]*mnuOther+0.125*m0rOther[5]*uOther[24]*mnuOther+0.125*m0rOther[2]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(72,64) = 0.1118033988749895*m0rOther[1]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[27]*mnuOther+0.125*m0rOther[3]*uOther[25]*mnuOther+0.225*m0rOther[2]*uOther[24]*mnuOther+0.125*m0rOther[5]*uOther[23]*mnuOther+0.225*m0rOther[4]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[21]*mnuOther+0.125*m0rOther[0]*uOther[21]*mnuOther-0.3535533905932737*m1rOther[21]*mnuOther+0.125*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(72,65) = 0.125*m0rOther[1]*uOther[26]*mnuOther+0.125*m0rOther[2]*uOther[25]*mnuOther+0.125*m0rOther[3]*uOther[24]*mnuOther+0.125*m0rOther[4]*uOther[23]*mnuOther+0.125*m0rOther[5]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[21]*mnuOther; 
  data->AEM_S(72,66) = 0.1118033988749895*m0rOther[3]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[28]*mnuOther+0.225*m0rOther[2]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[23]*mnuOther+0.125*m0rOther[0]*uOther[23]*mnuOther-0.3535533905932737*m1rOther[23]*mnuOther+0.225*m0rOther[6]*uOther[22]*mnuOther+0.125*m0rOther[5]*uOther[21]*mnuOther+0.125*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(72,67) = 0.125*m0rOther[2]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[24]*mnuOther+0.125*m0rOther[7]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[21]*mnuOther; 
  data->AEM_S(72,68) = 0.1964285714285714*m0rOther[2]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[23]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[22]*mnuOther-0.3162277660168379*m1rOther[22]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(72,69) = 0.125*m0rOther[2]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[23]*mnuOther+0.125*m0rOther[9]*uOther[22]*mnuOther; 
  data->AEM_S(73,60) = 0.1118033988749895*m0rOther[3]*uOther[29]*mnuOther+0.125*m0rOther[2]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[23]*mnuOther+0.125*m0rOther[0]*uOther[23]*mnuOther-0.3535533905932737*m1rOther[23]*mnuOther+0.125*m0rOther[6]*uOther[22]*mnuOther+0.125*m0rOther[5]*uOther[21]*mnuOther+0.125*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(73,61) = 0.1118033988749895*m0rOther[5]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[27]*mnuOther+0.125*m0rOther[4]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[25]*mnuOther+0.125*m0rOther[0]*uOther[25]*mnuOther-0.3535533905932737*m1rOther[25]*mnuOther+0.125*m0rOther[6]*uOther[24]*mnuOther+0.125*m0rOther[1]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[21]*mnuOther+0.125*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(73,62) = 0.1118033988749895*m0rOther[6]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[26]*mnuOther+0.125*m0rOther[0]*uOther[26]*mnuOther-0.3535533905932737*m1rOther[26]*mnuOther+0.125*m0rOther[4]*uOther[25]*mnuOther+0.125*m0rOther[5]*uOther[24]*mnuOther+0.125*m0rOther[2]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(73,63) = 0.1964285714285714*m0rOther[9]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[29]*mnuOther-0.3162277660168379*m1rOther[29]*mnuOther+0.125*m0rOther[8]*uOther[28]*mnuOther+0.125*m0rOther[7]*uOther[27]*mnuOther+0.225*m0rOther[6]*uOther[26]*mnuOther+0.225*m0rOther[5]*uOther[25]*mnuOther+0.125*m0rOther[4]*uOther[24]*mnuOther+0.225*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[1]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(73,64) = 0.125*m0rOther[1]*uOther[26]*mnuOther+0.125*m0rOther[2]*uOther[25]*mnuOther+0.125*m0rOther[3]*uOther[24]*mnuOther+0.125*m0rOther[4]*uOther[23]*mnuOther+0.125*m0rOther[5]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[21]*mnuOther; 
  data->AEM_S(73,65) = 0.1118033988749895*m0rOther[1]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[27]*mnuOther+0.225*m0rOther[3]*uOther[25]*mnuOther+0.125*m0rOther[2]*uOther[24]*mnuOther+0.225*m0rOther[5]*uOther[23]*mnuOther+0.125*m0rOther[4]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[21]*mnuOther+0.125*m0rOther[0]*uOther[21]*mnuOther-0.3535533905932737*m1rOther[21]*mnuOther+0.125*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(73,66) = 0.1118033988749895*m0rOther[2]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[28]*mnuOther+0.225*m0rOther[3]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[24]*mnuOther+0.225*m0rOther[6]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[22]*mnuOther+0.125*m0rOther[0]*uOther[22]*mnuOther-0.3535533905932737*m1rOther[22]*mnuOther+0.125*m0rOther[4]*uOther[21]*mnuOther+0.125*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(73,67) = 0.125*m0rOther[3]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[25]*mnuOther+0.125*m0rOther[7]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[21]*mnuOther; 
  data->AEM_S(73,68) = 0.125*m0rOther[3]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[26]*mnuOther+0.125*m0rOther[8]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[22]*mnuOther; 
  data->AEM_S(73,69) = 0.1964285714285714*m0rOther[3]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[23]*mnuOther-0.3162277660168379*m1rOther[23]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(74,60) = 0.1118033988749895*m0rOther[4]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[27]*mnuOther+0.125*m0rOther[5]*uOther[26]*mnuOther+0.125*m0rOther[6]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[24]*mnuOther+0.125*m0rOther[0]*uOther[24]*mnuOther-0.3535533905932737*m1rOther[24]*mnuOther+0.125*m0rOther[1]*uOther[22]*mnuOther+0.125*m0rOther[2]*uOther[21]*mnuOther+0.125*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(74,61) = 0.1118033988749895*m0rOther[2]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[27]*mnuOther+0.125*m0rOther[3]*uOther[26]*mnuOther+0.225*m0rOther[1]*uOther[24]*mnuOther+0.125*m0rOther[6]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[22]*mnuOther+0.125*m0rOther[0]*uOther[22]*mnuOther-0.3535533905932737*m1rOther[22]*mnuOther+0.225*m0rOther[4]*uOther[21]*mnuOther+0.125*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(74,62) = 0.1118033988749895*m0rOther[1]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[27]*mnuOther+0.125*m0rOther[3]*uOther[25]*mnuOther+0.225*m0rOther[2]*uOther[24]*mnuOther+0.125*m0rOther[5]*uOther[23]*mnuOther+0.225*m0rOther[4]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[21]*mnuOther+0.125*m0rOther[0]*uOther[21]*mnuOther-0.3535533905932737*m1rOther[21]*mnuOther+0.125*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(74,63) = 0.125*m0rOther[1]*uOther[26]*mnuOther+0.125*m0rOther[2]*uOther[25]*mnuOther+0.125*m0rOther[3]*uOther[24]*mnuOther+0.125*m0rOther[4]*uOther[23]*mnuOther+0.125*m0rOther[5]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[21]*mnuOther; 
  data->AEM_S(74,64) = 0.125*m0rOther[9]*uOther[29]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[28]*mnuOther+0.1*m0rOther[7]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[28]*mnuOther-0.3162277660168379*m1rOther[28]*mnuOther+0.1*m0rOther[8]*uOther[27]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[27]*mnuOther-0.3162277660168379*m1rOther[27]*mnuOther+0.225*m0rOther[6]*uOther[26]*mnuOther+0.225*m0rOther[5]*uOther[25]*mnuOther+0.405*m0rOther[4]*uOther[24]*mnuOther+0.125*m0rOther[3]*uOther[23]*mnuOther+0.225*m0rOther[2]*uOther[22]*mnuOther+0.225*m0rOther[1]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[20]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(74,65) = 0.1118033988749895*m0rOther[6]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[26]*mnuOther+0.125*m0rOther[0]*uOther[26]*mnuOther-0.3535533905932737*m1rOther[26]*mnuOther+0.225*m0rOther[4]*uOther[25]*mnuOther+0.225*m0rOther[5]*uOther[24]*mnuOther+0.125*m0rOther[2]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(74,66) = 0.1118033988749895*m0rOther[5]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[27]*mnuOther+0.225*m0rOther[4]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[25]*mnuOther+0.125*m0rOther[0]*uOther[25]*mnuOther-0.3535533905932737*m1rOther[25]*mnuOther+0.225*m0rOther[6]*uOther[24]*mnuOther+0.125*m0rOther[1]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[21]*mnuOther+0.125*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(74,67) = 0.1*m0rOther[4]*uOther[28]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[25]*mnuOther+0.1*m0rOther[8]*uOther[24]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[24]*mnuOther-0.3162277660168379*m1rOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(74,68) = 0.1964285714285714*m0rOther[4]*uOther[28]*mnuOther+0.1*m0rOther[4]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[24]*mnuOther+0.1*m0rOther[7]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[24]*mnuOther-0.3162277660168379*m1rOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(74,69) = 0.125*m0rOther[4]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[25]*mnuOther+0.125*m0rOther[9]*uOther[24]*mnuOther; 
  data->AEM_S(75,60) = 0.1118033988749895*m0rOther[5]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[27]*mnuOther+0.125*m0rOther[4]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[25]*mnuOther+0.125*m0rOther[0]*uOther[25]*mnuOther-0.3535533905932737*m1rOther[25]*mnuOther+0.125*m0rOther[6]*uOther[24]*mnuOther+0.125*m0rOther[1]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[21]*mnuOther+0.125*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(75,61) = 0.1118033988749895*m0rOther[3]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[27]*mnuOther+0.125*m0rOther[2]*uOther[26]*mnuOther+0.225*m0rOther[1]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[23]*mnuOther+0.125*m0rOther[0]*uOther[23]*mnuOther-0.3535533905932737*m1rOther[23]*mnuOther+0.125*m0rOther[6]*uOther[22]*mnuOther+0.225*m0rOther[5]*uOther[21]*mnuOther+0.125*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(75,62) = 0.125*m0rOther[1]*uOther[26]*mnuOther+0.125*m0rOther[2]*uOther[25]*mnuOther+0.125*m0rOther[3]*uOther[24]*mnuOther+0.125*m0rOther[4]*uOther[23]*mnuOther+0.125*m0rOther[5]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[21]*mnuOther; 
  data->AEM_S(75,63) = 0.1118033988749895*m0rOther[1]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[27]*mnuOther+0.225*m0rOther[3]*uOther[25]*mnuOther+0.125*m0rOther[2]*uOther[24]*mnuOther+0.225*m0rOther[5]*uOther[23]*mnuOther+0.125*m0rOther[4]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[21]*mnuOther+0.125*m0rOther[0]*uOther[21]*mnuOther-0.3535533905932737*m1rOther[21]*mnuOther+0.125*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(75,64) = 0.1118033988749895*m0rOther[6]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[26]*mnuOther+0.125*m0rOther[0]*uOther[26]*mnuOther-0.3535533905932737*m1rOther[26]*mnuOther+0.225*m0rOther[4]*uOther[25]*mnuOther+0.225*m0rOther[5]*uOther[24]*mnuOther+0.125*m0rOther[2]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(75,65) = 0.1964285714285714*m0rOther[9]*uOther[29]*mnuOther+0.1*m0rOther[7]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[29]*mnuOther-0.3162277660168379*m1rOther[29]*mnuOther+0.125*m0rOther[8]*uOther[28]*mnuOther+0.1*m0rOther[9]*uOther[27]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[27]*mnuOther-0.3162277660168379*m1rOther[27]*mnuOther+0.225*m0rOther[6]*uOther[26]*mnuOther+0.405*m0rOther[5]*uOther[25]*mnuOther+0.225*m0rOther[4]*uOther[24]*mnuOther+0.225*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[2]*uOther[22]*mnuOther+0.225*m0rOther[1]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[20]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(75,66) = 0.1118033988749895*m0rOther[4]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[27]*mnuOther+0.225*m0rOther[5]*uOther[26]*mnuOther+0.225*m0rOther[6]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[24]*mnuOther+0.125*m0rOther[0]*uOther[24]*mnuOther-0.3535533905932737*m1rOther[24]*mnuOther+0.125*m0rOther[1]*uOther[22]*mnuOther+0.125*m0rOther[2]*uOther[21]*mnuOther+0.125*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(75,67) = 0.1*m0rOther[5]*uOther[29]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[26]*mnuOther+0.1*m0rOther[9]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[25]*mnuOther-0.3162277660168379*m1rOther[25]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(75,68) = 0.125*m0rOther[5]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[26]*mnuOther+0.125*m0rOther[8]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[24]*mnuOther; 
  data->AEM_S(75,69) = 0.1964285714285714*m0rOther[5]*uOther[29]*mnuOther+0.1*m0rOther[5]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[26]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[25]*mnuOther+0.1*m0rOther[7]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[25]*mnuOther-0.3162277660168379*m1rOther[25]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(76,60) = 0.1118033988749895*m0rOther[6]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[26]*mnuOther+0.125*m0rOther[0]*uOther[26]*mnuOther-0.3535533905932737*m1rOther[26]*mnuOther+0.125*m0rOther[4]*uOther[25]*mnuOther+0.125*m0rOther[5]*uOther[24]*mnuOther+0.125*m0rOther[2]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(76,61) = 0.125*m0rOther[1]*uOther[26]*mnuOther+0.125*m0rOther[2]*uOther[25]*mnuOther+0.125*m0rOther[3]*uOther[24]*mnuOther+0.125*m0rOther[4]*uOther[23]*mnuOther+0.125*m0rOther[5]*uOther[22]*mnuOther+0.125*m0rOther[6]*uOther[21]*mnuOther; 
  data->AEM_S(76,62) = 0.1118033988749895*m0rOther[3]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[28]*mnuOther+0.225*m0rOther[2]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[23]*mnuOther+0.125*m0rOther[0]*uOther[23]*mnuOther-0.3535533905932737*m1rOther[23]*mnuOther+0.225*m0rOther[6]*uOther[22]*mnuOther+0.125*m0rOther[5]*uOther[21]*mnuOther+0.125*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(76,63) = 0.1118033988749895*m0rOther[2]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[28]*mnuOther+0.225*m0rOther[3]*uOther[26]*mnuOther+0.125*m0rOther[1]*uOther[24]*mnuOther+0.225*m0rOther[6]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[22]*mnuOther+0.125*m0rOther[0]*uOther[22]*mnuOther-0.3535533905932737*m1rOther[22]*mnuOther+0.125*m0rOther[4]*uOther[21]*mnuOther+0.125*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(76,64) = 0.1118033988749895*m0rOther[5]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[27]*mnuOther+0.225*m0rOther[4]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[25]*mnuOther+0.125*m0rOther[0]*uOther[25]*mnuOther-0.3535533905932737*m1rOther[25]*mnuOther+0.225*m0rOther[6]*uOther[24]*mnuOther+0.125*m0rOther[1]*uOther[23]*mnuOther+0.125*m0rOther[3]*uOther[21]*mnuOther+0.125*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(76,65) = 0.1118033988749895*m0rOther[4]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[27]*mnuOther+0.225*m0rOther[5]*uOther[26]*mnuOther+0.225*m0rOther[6]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[7]*uOther[24]*mnuOther+0.125*m0rOther[0]*uOther[24]*mnuOther-0.3535533905932737*m1rOther[24]*mnuOther+0.125*m0rOther[1]*uOther[22]*mnuOther+0.125*m0rOther[2]*uOther[21]*mnuOther+0.125*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(76,66) = 0.1964285714285714*m0rOther[9]*uOther[29]*mnuOther+0.1*m0rOther[8]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[29]*mnuOther-0.3162277660168379*m1rOther[29]*mnuOther+0.1*m0rOther[9]*uOther[28]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[28]*mnuOther-0.3162277660168379*m1rOther[28]*mnuOther+0.125*m0rOther[7]*uOther[27]*mnuOther+0.405*m0rOther[6]*uOther[26]*mnuOther+0.225*m0rOther[5]*uOther[25]*mnuOther+0.225*m0rOther[4]*uOther[24]*mnuOther+0.225*m0rOther[3]*uOther[23]*mnuOther+0.225*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[1]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[9]*uOther[20]*mnuOther+0.1118033988749895*m0rOther[8]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(76,67) = 0.125*m0rOther[6]*uOther[27]*mnuOther+0.125*m0rOther[7]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[24]*mnuOther; 
  data->AEM_S(76,68) = 0.1*m0rOther[6]*uOther[29]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[28]*mnuOther+0.1*m0rOther[9]*uOther[26]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[26]*mnuOther-0.3162277660168379*m1rOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(76,69) = 0.1964285714285714*m0rOther[6]*uOther[29]*mnuOther+0.1*m0rOther[6]*uOther[28]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[26]*mnuOther+0.1*m0rOther[8]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[26]*mnuOther-0.3162277660168379*m1rOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(77,60) = 0.07985957062499249*m0rOther[7]*uOther[27]*mnuOther+0.125*m0rOther[0]*uOther[27]*mnuOther-0.3535533905932737*m1rOther[27]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[21]*mnuOther+0.125*m0rOther[7]*uOther[20]*mnuOther; 
  data->AEM_S(77,61) = 0.1964285714285714*m0rOther[1]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[22]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[21]*mnuOther-0.3162277660168379*m1rOther[21]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[20]*mnuOther; 
  data->AEM_S(77,62) = 0.125*m0rOther[2]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[24]*mnuOther+0.125*m0rOther[7]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[21]*mnuOther; 
  data->AEM_S(77,63) = 0.125*m0rOther[3]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[25]*mnuOther+0.125*m0rOther[7]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[21]*mnuOther; 
  data->AEM_S(77,64) = 0.1*m0rOther[4]*uOther[28]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[25]*mnuOther+0.1*m0rOther[8]*uOther[24]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[24]*mnuOther-0.3162277660168379*m1rOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(77,65) = 0.1*m0rOther[5]*uOther[29]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[26]*mnuOther+0.1*m0rOther[9]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[7]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[25]*mnuOther-0.3162277660168379*m1rOther[25]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(77,66) = 0.125*m0rOther[6]*uOther[27]*mnuOther+0.125*m0rOther[7]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[24]*mnuOther; 
  data->AEM_S(77,67) = 0.125*m0rOther[9]*uOther[29]*mnuOther+0.125*m0rOther[8]*uOther[28]*mnuOther+0.2678571428571428*m0rOther[7]*uOther[27]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[27]*mnuOther-0.2258769757263128*m1rOther[27]*mnuOther+0.125*m0rOther[6]*uOther[26]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[24]*mnuOther+0.125*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[2]*uOther[22]*mnuOther+0.1964285714285714*m0rOther[1]*uOther[21]*mnuOther+0.07985957062499249*m0rOther[7]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(77,68) = 0.125*m0rOther[7]*uOther[28]*mnuOther+0.125*m0rOther[8]*uOther[27]*mnuOther+0.1*m0rOther[4]*uOther[24]*mnuOther; 
  data->AEM_S(77,69) = 0.125*m0rOther[7]*uOther[29]*mnuOther+0.125*m0rOther[9]*uOther[27]*mnuOther+0.1*m0rOther[5]*uOther[25]*mnuOther; 
  data->AEM_S(78,60) = 0.07985957062499249*m0rOther[8]*uOther[28]*mnuOther+0.125*m0rOther[0]*uOther[28]*mnuOther-0.3535533905932737*m1rOther[28]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[8]*uOther[20]*mnuOther; 
  data->AEM_S(78,61) = 0.125*m0rOther[1]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[22]*mnuOther+0.125*m0rOther[8]*uOther[21]*mnuOther; 
  data->AEM_S(78,62) = 0.1964285714285714*m0rOther[2]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[23]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[22]*mnuOther-0.3162277660168379*m1rOther[22]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[20]*mnuOther; 
  data->AEM_S(78,63) = 0.125*m0rOther[3]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[26]*mnuOther+0.125*m0rOther[8]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[22]*mnuOther; 
  data->AEM_S(78,64) = 0.1964285714285714*m0rOther[4]*uOther[28]*mnuOther+0.1*m0rOther[4]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[24]*mnuOther+0.1*m0rOther[7]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[24]*mnuOther-0.3162277660168379*m1rOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[20]*mnuOther; 
  data->AEM_S(78,65) = 0.125*m0rOther[5]*uOther[28]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[26]*mnuOther+0.125*m0rOther[8]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[24]*mnuOther; 
  data->AEM_S(78,66) = 0.1*m0rOther[6]*uOther[29]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[28]*mnuOther+0.1*m0rOther[9]*uOther[26]*mnuOther+0.1964285714285714*m0rOther[8]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[26]*mnuOther-0.3162277660168379*m1rOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(78,67) = 0.125*m0rOther[7]*uOther[28]*mnuOther+0.125*m0rOther[8]*uOther[27]*mnuOther+0.1*m0rOther[4]*uOther[24]*mnuOther; 
  data->AEM_S(78,68) = 0.125*m0rOther[9]*uOther[29]*mnuOther+0.2678571428571428*m0rOther[8]*uOther[28]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[28]*mnuOther-0.2258769757263128*m1rOther[28]*mnuOther+0.125*m0rOther[7]*uOther[27]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[26]*mnuOther+0.125*m0rOther[5]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[4]*uOther[24]*mnuOther+0.125*m0rOther[3]*uOther[23]*mnuOther+0.1964285714285714*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[1]*uOther[21]*mnuOther+0.07985957062499249*m0rOther[8]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
  data->AEM_S(78,69) = 0.125*m0rOther[8]*uOther[29]*mnuOther+0.125*m0rOther[9]*uOther[28]*mnuOther+0.1*m0rOther[6]*uOther[26]*mnuOther; 
  data->AEM_S(79,60) = 0.07985957062499249*m0rOther[9]*uOther[29]*mnuOther+0.125*m0rOther[0]*uOther[29]*mnuOther-0.3535533905932737*m1rOther[29]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[9]*uOther[20]*mnuOther; 
  data->AEM_S(79,61) = 0.125*m0rOther[1]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[23]*mnuOther+0.125*m0rOther[9]*uOther[21]*mnuOther; 
  data->AEM_S(79,62) = 0.125*m0rOther[2]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[23]*mnuOther+0.125*m0rOther[9]*uOther[22]*mnuOther; 
  data->AEM_S(79,63) = 0.1964285714285714*m0rOther[3]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[25]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[23]*mnuOther-0.3162277660168379*m1rOther[23]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[20]*mnuOther; 
  data->AEM_S(79,64) = 0.125*m0rOther[4]*uOther[29]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[25]*mnuOther+0.125*m0rOther[9]*uOther[24]*mnuOther; 
  data->AEM_S(79,65) = 0.1964285714285714*m0rOther[5]*uOther[29]*mnuOther+0.1*m0rOther[5]*uOther[27]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[26]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[25]*mnuOther+0.1*m0rOther[7]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[25]*mnuOther-0.3162277660168379*m1rOther[25]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[1]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[21]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[20]*mnuOther; 
  data->AEM_S(79,66) = 0.1964285714285714*m0rOther[6]*uOther[29]*mnuOther+0.1*m0rOther[6]*uOther[28]*mnuOther+0.1964285714285714*m0rOther[9]*uOther[26]*mnuOther+0.1*m0rOther[8]*uOther[26]*mnuOther+0.1118033988749895*m0rOther[0]*uOther[26]*mnuOther-0.3162277660168379*m1rOther[26]*mnuOther+0.1118033988749895*m0rOther[4]*uOther[25]*mnuOther+0.1118033988749895*m0rOther[5]*uOther[24]*mnuOther+0.1118033988749895*m0rOther[2]*uOther[23]*mnuOther+0.1118033988749895*m0rOther[3]*uOther[22]*mnuOther+0.1118033988749895*m0rOther[6]*uOther[20]*mnuOther; 
  data->AEM_S(79,67) = 0.125*m0rOther[7]*uOther[29]*mnuOther+0.125*m0rOther[9]*uOther[27]*mnuOther+0.1*m0rOther[5]*uOther[25]*mnuOther; 
  data->AEM_S(79,68) = 0.125*m0rOther[8]*uOther[29]*mnuOther+0.125*m0rOther[9]*uOther[28]*mnuOther+0.1*m0rOther[6]*uOther[26]*mnuOther; 
  data->AEM_S(79,69) = 0.2678571428571428*m0rOther[9]*uOther[29]*mnuOther+0.07985957062499249*m0rOther[0]*uOther[29]*mnuOther-0.2258769757263128*m1rOther[29]*mnuOther+0.125*m0rOther[8]*uOther[28]*mnuOther+0.125*m0rOther[7]*uOther[27]*mnuOther+0.1964285714285714*m0rOther[6]*uOther[26]*mnuOther+0.1964285714285714*m0rOther[5]*uOther[25]*mnuOther+0.125*m0rOther[4]*uOther[24]*mnuOther+0.1964285714285714*m0rOther[3]*uOther[23]*mnuOther+0.125*m0rOther[2]*uOther[22]*mnuOther+0.125*m0rOther[1]*uOther[21]*mnuOther+0.07985957062499249*m0rOther[9]*uOther[20]*mnuOther+0.125*m0rOther[0]*uOther[20]*mnuOther-0.3535533905932737*m1rOther[20]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfZ-m0Self*m1OtherZ)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[20] = 0.3535533905932737*m0rOther[9]*m1rSelf[29]-0.3535533905932737*m0rSelf[9]*m1rOther[29]+0.3535533905932737*m0rOther[8]*m1rSelf[28]-0.3535533905932737*m0rSelf[8]*m1rOther[28]+0.3535533905932737*m0rOther[7]*m1rSelf[27]-0.3535533905932737*m0rSelf[7]*m1rOther[27]+0.3535533905932737*m0rOther[6]*m1rSelf[26]-0.3535533905932737*m0rSelf[6]*m1rOther[26]+0.3535533905932737*m0rOther[5]*m1rSelf[25]-0.3535533905932737*m0rSelf[5]*m1rOther[25]+0.3535533905932737*m0rOther[4]*m1rSelf[24]-0.3535533905932737*m0rSelf[4]*m1rOther[24]+0.3535533905932737*m0rOther[3]*m1rSelf[23]-0.3535533905932737*m0rSelf[3]*m1rOther[23]+0.3535533905932737*m0rOther[2]*m1rSelf[22]-0.3535533905932737*m0rSelf[2]*m1rOther[22]+0.3535533905932737*m0rOther[1]*m1rSelf[21]-0.3535533905932737*m0rSelf[1]*m1rOther[21]+0.3535533905932737*m0rOther[0]*m1rSelf[20]-0.3535533905932737*m0rSelf[0]*m1rOther[20]; 
  m1EffD[21] = 0.3162277660168379*m0rOther[1]*m1rSelf[27]-0.3162277660168379*m0rSelf[1]*m1rOther[27]+0.3535533905932737*m0rOther[3]*m1rSelf[25]-0.3535533905932737*m0rSelf[3]*m1rOther[25]+0.3535533905932737*m0rOther[2]*m1rSelf[24]-0.3535533905932737*m0rSelf[2]*m1rOther[24]+0.3535533905932737*m0rOther[5]*m1rSelf[23]-0.3535533905932737*m0rSelf[5]*m1rOther[23]+0.3535533905932737*m0rOther[4]*m1rSelf[22]-0.3535533905932737*m0rSelf[4]*m1rOther[22]+0.3162277660168379*m0rOther[7]*m1rSelf[21]+0.3535533905932737*m0rOther[0]*m1rSelf[21]-0.3162277660168379*m0rSelf[7]*m1rOther[21]-0.3535533905932737*m0rSelf[0]*m1rOther[21]+0.3535533905932737*m0rOther[1]*m1rSelf[20]-0.3535533905932737*m0rSelf[1]*m1rOther[20]; 
  m1EffD[22] = 0.3162277660168379*m0rOther[2]*m1rSelf[28]-0.3162277660168379*m0rSelf[2]*m1rOther[28]+0.3535533905932737*m0rOther[3]*m1rSelf[26]-0.3535533905932737*m0rSelf[3]*m1rOther[26]+0.3535533905932737*m0rOther[1]*m1rSelf[24]-0.3535533905932737*m0rSelf[1]*m1rOther[24]+0.3535533905932737*m0rOther[6]*m1rSelf[23]-0.3535533905932737*m0rSelf[6]*m1rOther[23]+0.3162277660168379*m0rOther[8]*m1rSelf[22]+0.3535533905932737*m0rOther[0]*m1rSelf[22]-0.3162277660168379*m0rSelf[8]*m1rOther[22]-0.3535533905932737*m0rSelf[0]*m1rOther[22]+0.3535533905932737*m0rOther[4]*m1rSelf[21]-0.3535533905932737*m0rSelf[4]*m1rOther[21]+0.3535533905932737*m0rOther[2]*m1rSelf[20]-0.3535533905932737*m0rSelf[2]*m1rOther[20]; 
  m1EffD[23] = 0.3162277660168379*m0rOther[3]*m1rSelf[29]-0.3162277660168379*m0rSelf[3]*m1rOther[29]+0.3535533905932737*m0rOther[2]*m1rSelf[26]-0.3535533905932737*m0rSelf[2]*m1rOther[26]+0.3535533905932737*m0rOther[1]*m1rSelf[25]-0.3535533905932737*m0rSelf[1]*m1rOther[25]+0.3162277660168379*m0rOther[9]*m1rSelf[23]+0.3535533905932737*m0rOther[0]*m1rSelf[23]-0.3162277660168379*m0rSelf[9]*m1rOther[23]-0.3535533905932737*m0rSelf[0]*m1rOther[23]+0.3535533905932737*m0rOther[6]*m1rSelf[22]-0.3535533905932737*m0rSelf[6]*m1rOther[22]+0.3535533905932737*m0rOther[5]*m1rSelf[21]-0.3535533905932737*m0rSelf[5]*m1rOther[21]+0.3535533905932737*m0rOther[3]*m1rSelf[20]-0.3535533905932737*m0rSelf[3]*m1rOther[20]; 
  m1EffD[24] = 0.3162277660168379*m0rOther[4]*m1rSelf[28]-0.3162277660168379*m0rSelf[4]*m1rOther[28]+0.3162277660168379*m0rOther[4]*m1rSelf[27]-0.3162277660168379*m0rSelf[4]*m1rOther[27]+0.3535533905932737*m0rOther[5]*m1rSelf[26]-0.3535533905932737*m0rSelf[5]*m1rOther[26]+0.3535533905932737*m0rOther[6]*m1rSelf[25]-0.3535533905932737*m0rSelf[6]*m1rOther[25]+0.3162277660168379*m0rOther[8]*m1rSelf[24]+0.3162277660168379*m0rOther[7]*m1rSelf[24]+0.3535533905932737*m0rOther[0]*m1rSelf[24]-0.3162277660168379*m0rSelf[8]*m1rOther[24]-0.3162277660168379*m0rSelf[7]*m1rOther[24]-0.3535533905932737*m0rSelf[0]*m1rOther[24]+0.3535533905932737*m0rOther[1]*m1rSelf[22]-0.3535533905932737*m0rSelf[1]*m1rOther[22]+0.3535533905932737*m0rOther[2]*m1rSelf[21]-0.3535533905932737*m0rSelf[2]*m1rOther[21]+0.3535533905932737*m0rOther[4]*m1rSelf[20]-0.3535533905932737*m0rSelf[4]*m1rOther[20]; 
  m1EffD[25] = 0.3162277660168379*m0rOther[5]*m1rSelf[29]-0.3162277660168379*m0rSelf[5]*m1rOther[29]+0.3162277660168379*m0rOther[5]*m1rSelf[27]-0.3162277660168379*m0rSelf[5]*m1rOther[27]+0.3535533905932737*m0rOther[4]*m1rSelf[26]-0.3535533905932737*m0rSelf[4]*m1rOther[26]+0.3162277660168379*m0rOther[9]*m1rSelf[25]+0.3162277660168379*m0rOther[7]*m1rSelf[25]+0.3535533905932737*m0rOther[0]*m1rSelf[25]-0.3162277660168379*m0rSelf[9]*m1rOther[25]-0.3162277660168379*m0rSelf[7]*m1rOther[25]-0.3535533905932737*m0rSelf[0]*m1rOther[25]+0.3535533905932737*m0rOther[6]*m1rSelf[24]-0.3535533905932737*m0rSelf[6]*m1rOther[24]+0.3535533905932737*m0rOther[1]*m1rSelf[23]-0.3535533905932737*m0rSelf[1]*m1rOther[23]+0.3535533905932737*m0rOther[3]*m1rSelf[21]-0.3535533905932737*m0rSelf[3]*m1rOther[21]+0.3535533905932737*m0rOther[5]*m1rSelf[20]-0.3535533905932737*m0rSelf[5]*m1rOther[20]; 
  m1EffD[26] = 0.3162277660168379*m0rOther[6]*m1rSelf[29]-0.3162277660168379*m0rSelf[6]*m1rOther[29]+0.3162277660168379*m0rOther[6]*m1rSelf[28]-0.3162277660168379*m0rSelf[6]*m1rOther[28]+0.3162277660168379*m0rOther[9]*m1rSelf[26]+0.3162277660168379*m0rOther[8]*m1rSelf[26]+0.3535533905932737*m0rOther[0]*m1rSelf[26]-0.3162277660168379*m0rSelf[9]*m1rOther[26]-0.3162277660168379*m0rSelf[8]*m1rOther[26]-0.3535533905932737*m0rSelf[0]*m1rOther[26]+0.3535533905932737*m0rOther[4]*m1rSelf[25]-0.3535533905932737*m0rSelf[4]*m1rOther[25]+0.3535533905932737*m0rOther[5]*m1rSelf[24]-0.3535533905932737*m0rSelf[5]*m1rOther[24]+0.3535533905932737*m0rOther[2]*m1rSelf[23]-0.3535533905932737*m0rSelf[2]*m1rOther[23]+0.3535533905932737*m0rOther[3]*m1rSelf[22]-0.3535533905932737*m0rSelf[3]*m1rOther[22]+0.3535533905932737*m0rOther[6]*m1rSelf[20]-0.3535533905932737*m0rSelf[6]*m1rOther[20]; 
  m1EffD[27] = 0.2258769757263128*m0rOther[7]*m1rSelf[27]+0.3535533905932737*m0rOther[0]*m1rSelf[27]-0.2258769757263128*m0rSelf[7]*m1rOther[27]-0.3535533905932737*m0rSelf[0]*m1rOther[27]+0.3162277660168379*m0rOther[5]*m1rSelf[25]-0.3162277660168379*m0rSelf[5]*m1rOther[25]+0.3162277660168379*m0rOther[4]*m1rSelf[24]-0.3162277660168379*m0rSelf[4]*m1rOther[24]+0.3162277660168379*m0rOther[1]*m1rSelf[21]-0.3162277660168379*m0rSelf[1]*m1rOther[21]+0.3535533905932737*m0rOther[7]*m1rSelf[20]-0.3535533905932737*m0rSelf[7]*m1rOther[20]; 
  m1EffD[28] = 0.2258769757263128*m0rOther[8]*m1rSelf[28]+0.3535533905932737*m0rOther[0]*m1rSelf[28]-0.2258769757263128*m0rSelf[8]*m1rOther[28]-0.3535533905932737*m0rSelf[0]*m1rOther[28]+0.3162277660168379*m0rOther[6]*m1rSelf[26]-0.3162277660168379*m0rSelf[6]*m1rOther[26]+0.3162277660168379*m0rOther[4]*m1rSelf[24]-0.3162277660168379*m0rSelf[4]*m1rOther[24]+0.3162277660168379*m0rOther[2]*m1rSelf[22]-0.3162277660168379*m0rSelf[2]*m1rOther[22]+0.3535533905932737*m0rOther[8]*m1rSelf[20]-0.3535533905932737*m0rSelf[8]*m1rOther[20]; 
  m1EffD[29] = 0.2258769757263128*m0rOther[9]*m1rSelf[29]+0.3535533905932737*m0rOther[0]*m1rSelf[29]-0.2258769757263128*m0rSelf[9]*m1rOther[29]-0.3535533905932737*m0rSelf[0]*m1rOther[29]+0.3162277660168379*m0rOther[6]*m1rSelf[26]-0.3162277660168379*m0rSelf[6]*m1rOther[26]+0.3162277660168379*m0rOther[5]*m1rSelf[25]-0.3162277660168379*m0rSelf[5]*m1rOther[25]+0.3162277660168379*m0rOther[3]*m1rSelf[23]-0.3162277660168379*m0rSelf[3]*m1rOther[23]+0.3535533905932737*m0rOther[9]*m1rSelf[20]-0.3535533905932737*m0rSelf[9]*m1rOther[20]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[20],m1EffD[21],m1EffD[22],m1EffD[23],m1EffD[24],m1EffD[25],m1EffD[26],m1EffD[27],m1EffD[28],m1EffD[29]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+20,10,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 3 of momentum relaxation. 
  m1Relax[20] += (-2.0*m1EffD[20]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[20]*mnuSelf-1.0*m1rOther[20]*mnuOther; 
  m1Relax[21] += (-2.0*m1EffD[21]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[21]*mnuSelf-1.0*m1rOther[21]*mnuOther; 
  m1Relax[22] += (-2.0*m1EffD[22]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[22]*mnuSelf-1.0*m1rOther[22]*mnuOther; 
  m1Relax[23] += (-2.0*m1EffD[23]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[23]*mnuSelf-1.0*m1rOther[23]*mnuOther; 
  m1Relax[24] += (-2.0*m1EffD[24]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[24]*mnuSelf-1.0*m1rOther[24]*mnuOther; 
  m1Relax[25] += (-2.0*m1EffD[25]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[25]*mnuSelf-1.0*m1rOther[25]*mnuOther; 
  m1Relax[26] += (-2.0*m1EffD[26]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[26]*mnuSelf-1.0*m1rOther[26]*mnuOther; 
  m1Relax[27] += (-2.0*m1EffD[27]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[27]*mnuSelf-1.0*m1rOther[27]*mnuOther; 
  m1Relax[28] += (-2.0*m1EffD[28]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[28]*mnuSelf-1.0*m1rOther[28]*mnuOther; 
  m1Relax[29] += (-2.0*m1EffD[29]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[29]*mnuSelf-1.0*m1rOther[29]*mnuOther; 
 
  double ucMSelf[10]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.3535533905932737*cMSelf[a0+9]*uSelf[a0+9]+0.3535533905932737*cMSelf[a0+8]*uSelf[a0+8]+0.3535533905932737*cMSelf[a0+7]*uSelf[a0+7]+0.3535533905932737*cMSelf[a0+6]*uSelf[a0+6]+0.3535533905932737*cMSelf[a0+5]*uSelf[a0+5]+0.3535533905932737*cMSelf[a0+4]*uSelf[a0+4]+0.3535533905932737*cMSelf[a0+3]*uSelf[a0+3]+0.3535533905932737*cMSelf[a0+2]*uSelf[a0+2]+0.3535533905932737*cMSelf[a0+1]*uSelf[a0+1]+0.3535533905932737*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.3162277660168379*cMSelf[a0+1]*uSelf[a0+7]+0.3162277660168379*uSelf[a0+1]*cMSelf[a0+7]+0.3535533905932737*cMSelf[a0+3]*uSelf[a0+5]+0.3535533905932737*uSelf[a0+3]*cMSelf[a0+5]+0.3535533905932737*cMSelf[a0+2]*uSelf[a0+4]+0.3535533905932737*uSelf[a0+2]*cMSelf[a0+4]+0.3535533905932737*cMSelf[a0]*uSelf[a0+1]+0.3535533905932737*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.3162277660168379*cMSelf[a0+2]*uSelf[a0+8]+0.3162277660168379*uSelf[a0+2]*cMSelf[a0+8]+0.3535533905932737*cMSelf[a0+3]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+3]*cMSelf[a0+6]+0.3535533905932737*cMSelf[a0+1]*uSelf[a0+4]+0.3535533905932737*uSelf[a0+1]*cMSelf[a0+4]+0.3535533905932737*cMSelf[a0]*uSelf[a0+2]+0.3535533905932737*uSelf[a0]*cMSelf[a0+2]; 
    ucMSelf[3] += 0.3162277660168379*cMSelf[a0+3]*uSelf[a0+9]+0.3162277660168379*uSelf[a0+3]*cMSelf[a0+9]+0.3535533905932737*cMSelf[a0+2]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+2]*cMSelf[a0+6]+0.3535533905932737*cMSelf[a0+1]*uSelf[a0+5]+0.3535533905932737*uSelf[a0+1]*cMSelf[a0+5]+0.3535533905932737*cMSelf[a0]*uSelf[a0+3]+0.3535533905932737*uSelf[a0]*cMSelf[a0+3]; 
    ucMSelf[4] += 0.3162277660168379*cMSelf[a0+4]*uSelf[a0+8]+0.3162277660168379*uSelf[a0+4]*cMSelf[a0+8]+0.3162277660168379*cMSelf[a0+4]*uSelf[a0+7]+0.3162277660168379*uSelf[a0+4]*cMSelf[a0+7]+0.3535533905932737*cMSelf[a0+5]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+5]*cMSelf[a0+6]+0.3535533905932737*cMSelf[a0]*uSelf[a0+4]+0.3535533905932737*uSelf[a0]*cMSelf[a0+4]+0.3535533905932737*cMSelf[a0+1]*uSelf[a0+2]+0.3535533905932737*uSelf[a0+1]*cMSelf[a0+2]; 
    ucMSelf[5] += 0.3162277660168379*cMSelf[a0+5]*uSelf[a0+9]+0.3162277660168379*uSelf[a0+5]*cMSelf[a0+9]+0.3162277660168379*cMSelf[a0+5]*uSelf[a0+7]+0.3162277660168379*uSelf[a0+5]*cMSelf[a0+7]+0.3535533905932737*cMSelf[a0+4]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+4]*cMSelf[a0+6]+0.3535533905932737*cMSelf[a0]*uSelf[a0+5]+0.3535533905932737*uSelf[a0]*cMSelf[a0+5]+0.3535533905932737*cMSelf[a0+1]*uSelf[a0+3]+0.3535533905932737*uSelf[a0+1]*cMSelf[a0+3]; 
    ucMSelf[6] += 0.3162277660168379*cMSelf[a0+6]*uSelf[a0+9]+0.3162277660168379*uSelf[a0+6]*cMSelf[a0+9]+0.3162277660168379*cMSelf[a0+6]*uSelf[a0+8]+0.3162277660168379*uSelf[a0+6]*cMSelf[a0+8]+0.3535533905932737*cMSelf[a0]*uSelf[a0+6]+0.3535533905932737*uSelf[a0]*cMSelf[a0+6]+0.3535533905932737*cMSelf[a0+4]*uSelf[a0+5]+0.3535533905932737*uSelf[a0+4]*cMSelf[a0+5]+0.3535533905932737*cMSelf[a0+2]*uSelf[a0+3]+0.3535533905932737*uSelf[a0+2]*cMSelf[a0+3]; 
    ucMSelf[7] += 0.2258769757263128*cMSelf[a0+7]*uSelf[a0+7]+0.3535533905932737*cMSelf[a0]*uSelf[a0+7]+0.3535533905932737*uSelf[a0]*cMSelf[a0+7]+0.3162277660168379*cMSelf[a0+5]*uSelf[a0+5]+0.3162277660168379*cMSelf[a0+4]*uSelf[a0+4]+0.3162277660168379*cMSelf[a0+1]*uSelf[a0+1]; 
    ucMSelf[8] += 0.2258769757263128*cMSelf[a0+8]*uSelf[a0+8]+0.3535533905932737*cMSelf[a0]*uSelf[a0+8]+0.3535533905932737*uSelf[a0]*cMSelf[a0+8]+0.3162277660168379*cMSelf[a0+6]*uSelf[a0+6]+0.3162277660168379*cMSelf[a0+4]*uSelf[a0+4]+0.3162277660168379*cMSelf[a0+2]*uSelf[a0+2]; 
    ucMSelf[9] += 0.2258769757263128*cMSelf[a0+9]*uSelf[a0+9]+0.3535533905932737*cMSelf[a0]*uSelf[a0+9]+0.3535533905932737*uSelf[a0]*cMSelf[a0+9]+0.3162277660168379*cMSelf[a0+6]*uSelf[a0+6]+0.3162277660168379*cMSelf[a0+5]*uSelf[a0+5]+0.3162277660168379*cMSelf[a0+3]*uSelf[a0+3]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(70,30) = 0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(70,31) = 0.3535533905932737*ucMSelf[1]*mnuSelf+1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(70,32) = 0.3535533905932737*ucMSelf[2]*mnuSelf+1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(70,33) = 0.3535533905932737*ucMSelf[3]*mnuSelf+1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(70,34) = 0.3535533905932737*ucMSelf[4]*mnuSelf+1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(70,35) = 0.3535533905932737*ucMSelf[5]*mnuSelf+1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(70,36) = 0.3535533905932737*ucMSelf[6]*mnuSelf+1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(70,37) = 0.3535533905932737*ucMSelf[7]*mnuSelf+1.060660171779821*m0rSelf[7]*mnuSelf-0.3535533905932737*cESelf[7]*mnuSelf; 
  data->AEM_S(70,38) = 0.3535533905932737*ucMSelf[8]*mnuSelf+1.060660171779821*m0rSelf[8]*mnuSelf-0.3535533905932737*cESelf[8]*mnuSelf; 
  data->AEM_S(70,39) = 0.3535533905932737*ucMSelf[9]*mnuSelf+1.060660171779821*m0rSelf[9]*mnuSelf-0.3535533905932737*cESelf[9]*mnuSelf; 
  data->AEM_S(71,30) = 0.3535533905932737*ucMSelf[1]*mnuSelf+1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(71,31) = 0.3162277660168379*ucMSelf[7]*mnuSelf+0.9486832980505137*m0rSelf[7]*mnuSelf-0.3162277660168379*cESelf[7]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(71,32) = 0.3535533905932737*ucMSelf[4]*mnuSelf+1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(71,33) = 0.3535533905932737*ucMSelf[5]*mnuSelf+1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(71,34) = 0.3535533905932737*ucMSelf[2]*mnuSelf+1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(71,35) = 0.3535533905932737*ucMSelf[3]*mnuSelf+1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(71,37) = 0.3162277660168379*ucMSelf[1]*mnuSelf+0.9486832980505137*m0rSelf[1]*mnuSelf-0.3162277660168379*cESelf[1]*mnuSelf; 
  data->AEM_S(72,30) = 0.3535533905932737*ucMSelf[2]*mnuSelf+1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(72,31) = 0.3535533905932737*ucMSelf[4]*mnuSelf+1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(72,32) = 0.3162277660168379*ucMSelf[8]*mnuSelf+0.9486832980505137*m0rSelf[8]*mnuSelf-0.3162277660168379*cESelf[8]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(72,33) = 0.3535533905932737*ucMSelf[6]*mnuSelf+1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(72,34) = 0.3535533905932737*ucMSelf[1]*mnuSelf+1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(72,36) = 0.3535533905932737*ucMSelf[3]*mnuSelf+1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(72,38) = 0.3162277660168379*ucMSelf[2]*mnuSelf+0.9486832980505137*m0rSelf[2]*mnuSelf-0.3162277660168379*cESelf[2]*mnuSelf; 
  data->AEM_S(73,30) = 0.3535533905932737*ucMSelf[3]*mnuSelf+1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(73,31) = 0.3535533905932737*ucMSelf[5]*mnuSelf+1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(73,32) = 0.3535533905932737*ucMSelf[6]*mnuSelf+1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(73,33) = 0.3162277660168379*ucMSelf[9]*mnuSelf+0.9486832980505137*m0rSelf[9]*mnuSelf-0.3162277660168379*cESelf[9]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(73,35) = 0.3535533905932737*ucMSelf[1]*mnuSelf+1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(73,36) = 0.3535533905932737*ucMSelf[2]*mnuSelf+1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(73,39) = 0.3162277660168379*ucMSelf[3]*mnuSelf+0.9486832980505137*m0rSelf[3]*mnuSelf-0.3162277660168379*cESelf[3]*mnuSelf; 
  data->AEM_S(74,30) = 0.3535533905932737*ucMSelf[4]*mnuSelf+1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(74,31) = 0.3535533905932737*ucMSelf[2]*mnuSelf+1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(74,32) = 0.3535533905932737*ucMSelf[1]*mnuSelf+1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(74,34) = 0.3162277660168379*ucMSelf[8]*mnuSelf+0.9486832980505137*m0rSelf[8]*mnuSelf-0.3162277660168379*cESelf[8]*mnuSelf+0.3162277660168379*ucMSelf[7]*mnuSelf+0.9486832980505137*m0rSelf[7]*mnuSelf-0.3162277660168379*cESelf[7]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(74,35) = 0.3535533905932737*ucMSelf[6]*mnuSelf+1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(74,36) = 0.3535533905932737*ucMSelf[5]*mnuSelf+1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(74,37) = 0.3162277660168379*ucMSelf[4]*mnuSelf+0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(74,38) = 0.3162277660168379*ucMSelf[4]*mnuSelf+0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(75,30) = 0.3535533905932737*ucMSelf[5]*mnuSelf+1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(75,31) = 0.3535533905932737*ucMSelf[3]*mnuSelf+1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(75,33) = 0.3535533905932737*ucMSelf[1]*mnuSelf+1.060660171779821*m0rSelf[1]*mnuSelf-0.3535533905932737*cESelf[1]*mnuSelf; 
  data->AEM_S(75,34) = 0.3535533905932737*ucMSelf[6]*mnuSelf+1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(75,35) = 0.3162277660168379*ucMSelf[9]*mnuSelf+0.9486832980505137*m0rSelf[9]*mnuSelf-0.3162277660168379*cESelf[9]*mnuSelf+0.3162277660168379*ucMSelf[7]*mnuSelf+0.9486832980505137*m0rSelf[7]*mnuSelf-0.3162277660168379*cESelf[7]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(75,36) = 0.3535533905932737*ucMSelf[4]*mnuSelf+1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(75,37) = 0.3162277660168379*ucMSelf[5]*mnuSelf+0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(75,39) = 0.3162277660168379*ucMSelf[5]*mnuSelf+0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(76,30) = 0.3535533905932737*ucMSelf[6]*mnuSelf+1.060660171779821*m0rSelf[6]*mnuSelf-0.3535533905932737*cESelf[6]*mnuSelf; 
  data->AEM_S(76,32) = 0.3535533905932737*ucMSelf[3]*mnuSelf+1.060660171779821*m0rSelf[3]*mnuSelf-0.3535533905932737*cESelf[3]*mnuSelf; 
  data->AEM_S(76,33) = 0.3535533905932737*ucMSelf[2]*mnuSelf+1.060660171779821*m0rSelf[2]*mnuSelf-0.3535533905932737*cESelf[2]*mnuSelf; 
  data->AEM_S(76,34) = 0.3535533905932737*ucMSelf[5]*mnuSelf+1.060660171779821*m0rSelf[5]*mnuSelf-0.3535533905932737*cESelf[5]*mnuSelf; 
  data->AEM_S(76,35) = 0.3535533905932737*ucMSelf[4]*mnuSelf+1.060660171779821*m0rSelf[4]*mnuSelf-0.3535533905932737*cESelf[4]*mnuSelf; 
  data->AEM_S(76,36) = 0.3162277660168379*ucMSelf[9]*mnuSelf+0.9486832980505137*m0rSelf[9]*mnuSelf-0.3162277660168379*cESelf[9]*mnuSelf+0.3162277660168379*ucMSelf[8]*mnuSelf+0.9486832980505137*m0rSelf[8]*mnuSelf-0.3162277660168379*cESelf[8]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(76,38) = 0.3162277660168379*ucMSelf[6]*mnuSelf+0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(76,39) = 0.3162277660168379*ucMSelf[6]*mnuSelf+0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(77,30) = 0.3535533905932737*ucMSelf[7]*mnuSelf+1.060660171779821*m0rSelf[7]*mnuSelf-0.3535533905932737*cESelf[7]*mnuSelf; 
  data->AEM_S(77,31) = 0.3162277660168379*ucMSelf[1]*mnuSelf+0.9486832980505137*m0rSelf[1]*mnuSelf-0.3162277660168379*cESelf[1]*mnuSelf; 
  data->AEM_S(77,34) = 0.3162277660168379*ucMSelf[4]*mnuSelf+0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(77,35) = 0.3162277660168379*ucMSelf[5]*mnuSelf+0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(77,37) = 0.2258769757263128*ucMSelf[7]*mnuSelf+0.6776309271789384*m0rSelf[7]*mnuSelf-0.2258769757263128*cESelf[7]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(78,30) = 0.3535533905932737*ucMSelf[8]*mnuSelf+1.060660171779821*m0rSelf[8]*mnuSelf-0.3535533905932737*cESelf[8]*mnuSelf; 
  data->AEM_S(78,32) = 0.3162277660168379*ucMSelf[2]*mnuSelf+0.9486832980505137*m0rSelf[2]*mnuSelf-0.3162277660168379*cESelf[2]*mnuSelf; 
  data->AEM_S(78,34) = 0.3162277660168379*ucMSelf[4]*mnuSelf+0.9486832980505137*m0rSelf[4]*mnuSelf-0.3162277660168379*cESelf[4]*mnuSelf; 
  data->AEM_S(78,36) = 0.3162277660168379*ucMSelf[6]*mnuSelf+0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(78,38) = 0.2258769757263128*ucMSelf[8]*mnuSelf+0.6776309271789384*m0rSelf[8]*mnuSelf-0.2258769757263128*cESelf[8]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
  data->AEM_S(79,30) = 0.3535533905932737*ucMSelf[9]*mnuSelf+1.060660171779821*m0rSelf[9]*mnuSelf-0.3535533905932737*cESelf[9]*mnuSelf; 
  data->AEM_S(79,33) = 0.3162277660168379*ucMSelf[3]*mnuSelf+0.9486832980505137*m0rSelf[3]*mnuSelf-0.3162277660168379*cESelf[3]*mnuSelf; 
  data->AEM_S(79,35) = 0.3162277660168379*ucMSelf[5]*mnuSelf+0.9486832980505137*m0rSelf[5]*mnuSelf-0.3162277660168379*cESelf[5]*mnuSelf; 
  data->AEM_S(79,36) = 0.3162277660168379*ucMSelf[6]*mnuSelf+0.9486832980505137*m0rSelf[6]*mnuSelf-0.3162277660168379*cESelf[6]*mnuSelf; 
  data->AEM_S(79,39) = 0.2258769757263128*ucMSelf[9]*mnuSelf+0.6776309271789384*m0rSelf[9]*mnuSelf-0.2258769757263128*cESelf[9]*mnuSelf+0.3535533905932737*ucMSelf[0]*mnuSelf+1.060660171779821*m0rSelf[0]*mnuSelf-0.3535533905932737*cESelf[0]*mnuSelf; 
 
  double ucMOther[10]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.3535533905932737*cMOther[a0+9]*uOther[a0+9]+0.3535533905932737*cMOther[a0+8]*uOther[a0+8]+0.3535533905932737*cMOther[a0+7]*uOther[a0+7]+0.3535533905932737*cMOther[a0+6]*uOther[a0+6]+0.3535533905932737*cMOther[a0+5]*uOther[a0+5]+0.3535533905932737*cMOther[a0+4]*uOther[a0+4]+0.3535533905932737*cMOther[a0+3]*uOther[a0+3]+0.3535533905932737*cMOther[a0+2]*uOther[a0+2]+0.3535533905932737*cMOther[a0+1]*uOther[a0+1]+0.3535533905932737*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.3162277660168379*cMOther[a0+1]*uOther[a0+7]+0.3162277660168379*uOther[a0+1]*cMOther[a0+7]+0.3535533905932737*cMOther[a0+3]*uOther[a0+5]+0.3535533905932737*uOther[a0+3]*cMOther[a0+5]+0.3535533905932737*cMOther[a0+2]*uOther[a0+4]+0.3535533905932737*uOther[a0+2]*cMOther[a0+4]+0.3535533905932737*cMOther[a0]*uOther[a0+1]+0.3535533905932737*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.3162277660168379*cMOther[a0+2]*uOther[a0+8]+0.3162277660168379*uOther[a0+2]*cMOther[a0+8]+0.3535533905932737*cMOther[a0+3]*uOther[a0+6]+0.3535533905932737*uOther[a0+3]*cMOther[a0+6]+0.3535533905932737*cMOther[a0+1]*uOther[a0+4]+0.3535533905932737*uOther[a0+1]*cMOther[a0+4]+0.3535533905932737*cMOther[a0]*uOther[a0+2]+0.3535533905932737*uOther[a0]*cMOther[a0+2]; 
    ucMOther[3] += 0.3162277660168379*cMOther[a0+3]*uOther[a0+9]+0.3162277660168379*uOther[a0+3]*cMOther[a0+9]+0.3535533905932737*cMOther[a0+2]*uOther[a0+6]+0.3535533905932737*uOther[a0+2]*cMOther[a0+6]+0.3535533905932737*cMOther[a0+1]*uOther[a0+5]+0.3535533905932737*uOther[a0+1]*cMOther[a0+5]+0.3535533905932737*cMOther[a0]*uOther[a0+3]+0.3535533905932737*uOther[a0]*cMOther[a0+3]; 
    ucMOther[4] += 0.3162277660168379*cMOther[a0+4]*uOther[a0+8]+0.3162277660168379*uOther[a0+4]*cMOther[a0+8]+0.3162277660168379*cMOther[a0+4]*uOther[a0+7]+0.3162277660168379*uOther[a0+4]*cMOther[a0+7]+0.3535533905932737*cMOther[a0+5]*uOther[a0+6]+0.3535533905932737*uOther[a0+5]*cMOther[a0+6]+0.3535533905932737*cMOther[a0]*uOther[a0+4]+0.3535533905932737*uOther[a0]*cMOther[a0+4]+0.3535533905932737*cMOther[a0+1]*uOther[a0+2]+0.3535533905932737*uOther[a0+1]*cMOther[a0+2]; 
    ucMOther[5] += 0.3162277660168379*cMOther[a0+5]*uOther[a0+9]+0.3162277660168379*uOther[a0+5]*cMOther[a0+9]+0.3162277660168379*cMOther[a0+5]*uOther[a0+7]+0.3162277660168379*uOther[a0+5]*cMOther[a0+7]+0.3535533905932737*cMOther[a0+4]*uOther[a0+6]+0.3535533905932737*uOther[a0+4]*cMOther[a0+6]+0.3535533905932737*cMOther[a0]*uOther[a0+5]+0.3535533905932737*uOther[a0]*cMOther[a0+5]+0.3535533905932737*cMOther[a0+1]*uOther[a0+3]+0.3535533905932737*uOther[a0+1]*cMOther[a0+3]; 
    ucMOther[6] += 0.3162277660168379*cMOther[a0+6]*uOther[a0+9]+0.3162277660168379*uOther[a0+6]*cMOther[a0+9]+0.3162277660168379*cMOther[a0+6]*uOther[a0+8]+0.3162277660168379*uOther[a0+6]*cMOther[a0+8]+0.3535533905932737*cMOther[a0]*uOther[a0+6]+0.3535533905932737*uOther[a0]*cMOther[a0+6]+0.3535533905932737*cMOther[a0+4]*uOther[a0+5]+0.3535533905932737*uOther[a0+4]*cMOther[a0+5]+0.3535533905932737*cMOther[a0+2]*uOther[a0+3]+0.3535533905932737*uOther[a0+2]*cMOther[a0+3]; 
    ucMOther[7] += 0.2258769757263128*cMOther[a0+7]*uOther[a0+7]+0.3535533905932737*cMOther[a0]*uOther[a0+7]+0.3535533905932737*uOther[a0]*cMOther[a0+7]+0.3162277660168379*cMOther[a0+5]*uOther[a0+5]+0.3162277660168379*cMOther[a0+4]*uOther[a0+4]+0.3162277660168379*cMOther[a0+1]*uOther[a0+1]; 
    ucMOther[8] += 0.2258769757263128*cMOther[a0+8]*uOther[a0+8]+0.3535533905932737*cMOther[a0]*uOther[a0+8]+0.3535533905932737*uOther[a0]*cMOther[a0+8]+0.3162277660168379*cMOther[a0+6]*uOther[a0+6]+0.3162277660168379*cMOther[a0+4]*uOther[a0+4]+0.3162277660168379*cMOther[a0+2]*uOther[a0+2]; 
    ucMOther[9] += 0.2258769757263128*cMOther[a0+9]*uOther[a0+9]+0.3535533905932737*cMOther[a0]*uOther[a0+9]+0.3535533905932737*uOther[a0]*cMOther[a0+9]+0.3162277660168379*cMOther[a0+6]*uOther[a0+6]+0.3162277660168379*cMOther[a0+5]*uOther[a0+5]+0.3162277660168379*cMOther[a0+3]*uOther[a0+3]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(70,70) = (-0.3535533905932737*ucMOther[0]*mnuOther)-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(70,71) = (-0.3535533905932737*ucMOther[1]*mnuOther)-1.060660171779821*m0rOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(70,72) = (-0.3535533905932737*ucMOther[2]*mnuOther)-1.060660171779821*m0rOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(70,73) = (-0.3535533905932737*ucMOther[3]*mnuOther)-1.060660171779821*m0rOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(70,74) = (-0.3535533905932737*ucMOther[4]*mnuOther)-1.060660171779821*m0rOther[4]*mnuOther+0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(70,75) = (-0.3535533905932737*ucMOther[5]*mnuOther)-1.060660171779821*m0rOther[5]*mnuOther+0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(70,76) = (-0.3535533905932737*ucMOther[6]*mnuOther)-1.060660171779821*m0rOther[6]*mnuOther+0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(70,77) = (-0.3535533905932737*ucMOther[7]*mnuOther)-1.060660171779821*m0rOther[7]*mnuOther+0.3535533905932737*cEOther[7]*mnuOther; 
  data->AEM_S(70,78) = (-0.3535533905932737*ucMOther[8]*mnuOther)-1.060660171779821*m0rOther[8]*mnuOther+0.3535533905932737*cEOther[8]*mnuOther; 
  data->AEM_S(70,79) = (-0.3535533905932737*ucMOther[9]*mnuOther)-1.060660171779821*m0rOther[9]*mnuOther+0.3535533905932737*cEOther[9]*mnuOther; 
  data->AEM_S(71,70) = (-0.3535533905932737*ucMOther[1]*mnuOther)-1.060660171779821*m0rOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(71,71) = (-0.3162277660168379*ucMOther[7]*mnuOther)-0.9486832980505137*m0rOther[7]*mnuOther+0.3162277660168379*cEOther[7]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(71,72) = (-0.3535533905932737*ucMOther[4]*mnuOther)-1.060660171779821*m0rOther[4]*mnuOther+0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(71,73) = (-0.3535533905932737*ucMOther[5]*mnuOther)-1.060660171779821*m0rOther[5]*mnuOther+0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(71,74) = (-0.3535533905932737*ucMOther[2]*mnuOther)-1.060660171779821*m0rOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(71,75) = (-0.3535533905932737*ucMOther[3]*mnuOther)-1.060660171779821*m0rOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(71,77) = (-0.3162277660168379*ucMOther[1]*mnuOther)-0.9486832980505137*m0rOther[1]*mnuOther+0.3162277660168379*cEOther[1]*mnuOther; 
  data->AEM_S(72,70) = (-0.3535533905932737*ucMOther[2]*mnuOther)-1.060660171779821*m0rOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(72,71) = (-0.3535533905932737*ucMOther[4]*mnuOther)-1.060660171779821*m0rOther[4]*mnuOther+0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(72,72) = (-0.3162277660168379*ucMOther[8]*mnuOther)-0.9486832980505137*m0rOther[8]*mnuOther+0.3162277660168379*cEOther[8]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(72,73) = (-0.3535533905932737*ucMOther[6]*mnuOther)-1.060660171779821*m0rOther[6]*mnuOther+0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(72,74) = (-0.3535533905932737*ucMOther[1]*mnuOther)-1.060660171779821*m0rOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(72,76) = (-0.3535533905932737*ucMOther[3]*mnuOther)-1.060660171779821*m0rOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(72,78) = (-0.3162277660168379*ucMOther[2]*mnuOther)-0.9486832980505137*m0rOther[2]*mnuOther+0.3162277660168379*cEOther[2]*mnuOther; 
  data->AEM_S(73,70) = (-0.3535533905932737*ucMOther[3]*mnuOther)-1.060660171779821*m0rOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(73,71) = (-0.3535533905932737*ucMOther[5]*mnuOther)-1.060660171779821*m0rOther[5]*mnuOther+0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(73,72) = (-0.3535533905932737*ucMOther[6]*mnuOther)-1.060660171779821*m0rOther[6]*mnuOther+0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(73,73) = (-0.3162277660168379*ucMOther[9]*mnuOther)-0.9486832980505137*m0rOther[9]*mnuOther+0.3162277660168379*cEOther[9]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(73,75) = (-0.3535533905932737*ucMOther[1]*mnuOther)-1.060660171779821*m0rOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(73,76) = (-0.3535533905932737*ucMOther[2]*mnuOther)-1.060660171779821*m0rOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(73,79) = (-0.3162277660168379*ucMOther[3]*mnuOther)-0.9486832980505137*m0rOther[3]*mnuOther+0.3162277660168379*cEOther[3]*mnuOther; 
  data->AEM_S(74,70) = (-0.3535533905932737*ucMOther[4]*mnuOther)-1.060660171779821*m0rOther[4]*mnuOther+0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(74,71) = (-0.3535533905932737*ucMOther[2]*mnuOther)-1.060660171779821*m0rOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(74,72) = (-0.3535533905932737*ucMOther[1]*mnuOther)-1.060660171779821*m0rOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(74,74) = (-0.3162277660168379*ucMOther[8]*mnuOther)-0.9486832980505137*m0rOther[8]*mnuOther+0.3162277660168379*cEOther[8]*mnuOther-0.3162277660168379*ucMOther[7]*mnuOther-0.9486832980505137*m0rOther[7]*mnuOther+0.3162277660168379*cEOther[7]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(74,75) = (-0.3535533905932737*ucMOther[6]*mnuOther)-1.060660171779821*m0rOther[6]*mnuOther+0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(74,76) = (-0.3535533905932737*ucMOther[5]*mnuOther)-1.060660171779821*m0rOther[5]*mnuOther+0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(74,77) = (-0.3162277660168379*ucMOther[4]*mnuOther)-0.9486832980505137*m0rOther[4]*mnuOther+0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(74,78) = (-0.3162277660168379*ucMOther[4]*mnuOther)-0.9486832980505137*m0rOther[4]*mnuOther+0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(75,70) = (-0.3535533905932737*ucMOther[5]*mnuOther)-1.060660171779821*m0rOther[5]*mnuOther+0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(75,71) = (-0.3535533905932737*ucMOther[3]*mnuOther)-1.060660171779821*m0rOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(75,73) = (-0.3535533905932737*ucMOther[1]*mnuOther)-1.060660171779821*m0rOther[1]*mnuOther+0.3535533905932737*cEOther[1]*mnuOther; 
  data->AEM_S(75,74) = (-0.3535533905932737*ucMOther[6]*mnuOther)-1.060660171779821*m0rOther[6]*mnuOther+0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(75,75) = (-0.3162277660168379*ucMOther[9]*mnuOther)-0.9486832980505137*m0rOther[9]*mnuOther+0.3162277660168379*cEOther[9]*mnuOther-0.3162277660168379*ucMOther[7]*mnuOther-0.9486832980505137*m0rOther[7]*mnuOther+0.3162277660168379*cEOther[7]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(75,76) = (-0.3535533905932737*ucMOther[4]*mnuOther)-1.060660171779821*m0rOther[4]*mnuOther+0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(75,77) = (-0.3162277660168379*ucMOther[5]*mnuOther)-0.9486832980505137*m0rOther[5]*mnuOther+0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(75,79) = (-0.3162277660168379*ucMOther[5]*mnuOther)-0.9486832980505137*m0rOther[5]*mnuOther+0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(76,70) = (-0.3535533905932737*ucMOther[6]*mnuOther)-1.060660171779821*m0rOther[6]*mnuOther+0.3535533905932737*cEOther[6]*mnuOther; 
  data->AEM_S(76,72) = (-0.3535533905932737*ucMOther[3]*mnuOther)-1.060660171779821*m0rOther[3]*mnuOther+0.3535533905932737*cEOther[3]*mnuOther; 
  data->AEM_S(76,73) = (-0.3535533905932737*ucMOther[2]*mnuOther)-1.060660171779821*m0rOther[2]*mnuOther+0.3535533905932737*cEOther[2]*mnuOther; 
  data->AEM_S(76,74) = (-0.3535533905932737*ucMOther[5]*mnuOther)-1.060660171779821*m0rOther[5]*mnuOther+0.3535533905932737*cEOther[5]*mnuOther; 
  data->AEM_S(76,75) = (-0.3535533905932737*ucMOther[4]*mnuOther)-1.060660171779821*m0rOther[4]*mnuOther+0.3535533905932737*cEOther[4]*mnuOther; 
  data->AEM_S(76,76) = (-0.3162277660168379*ucMOther[9]*mnuOther)-0.9486832980505137*m0rOther[9]*mnuOther+0.3162277660168379*cEOther[9]*mnuOther-0.3162277660168379*ucMOther[8]*mnuOther-0.9486832980505137*m0rOther[8]*mnuOther+0.3162277660168379*cEOther[8]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(76,78) = (-0.3162277660168379*ucMOther[6]*mnuOther)-0.9486832980505137*m0rOther[6]*mnuOther+0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(76,79) = (-0.3162277660168379*ucMOther[6]*mnuOther)-0.9486832980505137*m0rOther[6]*mnuOther+0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(77,70) = (-0.3535533905932737*ucMOther[7]*mnuOther)-1.060660171779821*m0rOther[7]*mnuOther+0.3535533905932737*cEOther[7]*mnuOther; 
  data->AEM_S(77,71) = (-0.3162277660168379*ucMOther[1]*mnuOther)-0.9486832980505137*m0rOther[1]*mnuOther+0.3162277660168379*cEOther[1]*mnuOther; 
  data->AEM_S(77,74) = (-0.3162277660168379*ucMOther[4]*mnuOther)-0.9486832980505137*m0rOther[4]*mnuOther+0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(77,75) = (-0.3162277660168379*ucMOther[5]*mnuOther)-0.9486832980505137*m0rOther[5]*mnuOther+0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(77,77) = (-0.2258769757263128*ucMOther[7]*mnuOther)-0.6776309271789384*m0rOther[7]*mnuOther+0.2258769757263128*cEOther[7]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(78,70) = (-0.3535533905932737*ucMOther[8]*mnuOther)-1.060660171779821*m0rOther[8]*mnuOther+0.3535533905932737*cEOther[8]*mnuOther; 
  data->AEM_S(78,72) = (-0.3162277660168379*ucMOther[2]*mnuOther)-0.9486832980505137*m0rOther[2]*mnuOther+0.3162277660168379*cEOther[2]*mnuOther; 
  data->AEM_S(78,74) = (-0.3162277660168379*ucMOther[4]*mnuOther)-0.9486832980505137*m0rOther[4]*mnuOther+0.3162277660168379*cEOther[4]*mnuOther; 
  data->AEM_S(78,76) = (-0.3162277660168379*ucMOther[6]*mnuOther)-0.9486832980505137*m0rOther[6]*mnuOther+0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(78,78) = (-0.2258769757263128*ucMOther[8]*mnuOther)-0.6776309271789384*m0rOther[8]*mnuOther+0.2258769757263128*cEOther[8]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
  data->AEM_S(79,70) = (-0.3535533905932737*ucMOther[9]*mnuOther)-1.060660171779821*m0rOther[9]*mnuOther+0.3535533905932737*cEOther[9]*mnuOther; 
  data->AEM_S(79,73) = (-0.3162277660168379*ucMOther[3]*mnuOther)-0.9486832980505137*m0rOther[3]*mnuOther+0.3162277660168379*cEOther[3]*mnuOther; 
  data->AEM_S(79,75) = (-0.3162277660168379*ucMOther[5]*mnuOther)-0.9486832980505137*m0rOther[5]*mnuOther+0.3162277660168379*cEOther[5]*mnuOther; 
  data->AEM_S(79,76) = (-0.3162277660168379*ucMOther[6]*mnuOther)-0.9486832980505137*m0rOther[6]*mnuOther+0.3162277660168379*cEOther[6]*mnuOther; 
  data->AEM_S(79,79) = (-0.2258769757263128*ucMOther[9]*mnuOther)-0.6776309271789384*m0rOther[9]*mnuOther+0.2258769757263128*cEOther[9]*mnuOther-0.3535533905932737*ucMOther[0]*mnuOther-1.060660171779821*m0rOther[0]*mnuOther+0.3535533905932737*cEOther[0]*mnuOther; 
 
  double kinESelf[10]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinESelf[0] += 0.3535533905932737*m1rSelf[a0+9]*uSelf[a0+9]+0.3535533905932737*m1rSelf[a0+8]*uSelf[a0+8]+0.3535533905932737*m1rSelf[a0+7]*uSelf[a0+7]+0.3535533905932737*m1rSelf[a0+6]*uSelf[a0+6]+0.3535533905932737*m1rSelf[a0+5]*uSelf[a0+5]+0.3535533905932737*m1rSelf[a0+4]*uSelf[a0+4]+0.3535533905932737*m1rSelf[a0+3]*uSelf[a0+3]+0.3535533905932737*m1rSelf[a0+2]*uSelf[a0+2]+0.3535533905932737*m1rSelf[a0+1]*uSelf[a0+1]+0.3535533905932737*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.3162277660168379*m1rSelf[a0+1]*uSelf[a0+7]+0.3162277660168379*uSelf[a0+1]*m1rSelf[a0+7]+0.3535533905932737*m1rSelf[a0+3]*uSelf[a0+5]+0.3535533905932737*uSelf[a0+3]*m1rSelf[a0+5]+0.3535533905932737*m1rSelf[a0+2]*uSelf[a0+4]+0.3535533905932737*uSelf[a0+2]*m1rSelf[a0+4]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+1]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.3162277660168379*m1rSelf[a0+2]*uSelf[a0+8]+0.3162277660168379*uSelf[a0+2]*m1rSelf[a0+8]+0.3535533905932737*m1rSelf[a0+3]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+3]*m1rSelf[a0+6]+0.3535533905932737*m1rSelf[a0+1]*uSelf[a0+4]+0.3535533905932737*uSelf[a0+1]*m1rSelf[a0+4]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+2]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+2]; 
    kinESelf[3] += 0.3162277660168379*m1rSelf[a0+3]*uSelf[a0+9]+0.3162277660168379*uSelf[a0+3]*m1rSelf[a0+9]+0.3535533905932737*m1rSelf[a0+2]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+2]*m1rSelf[a0+6]+0.3535533905932737*m1rSelf[a0+1]*uSelf[a0+5]+0.3535533905932737*uSelf[a0+1]*m1rSelf[a0+5]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+3]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+3]; 
    kinESelf[4] += 0.3162277660168379*m1rSelf[a0+4]*uSelf[a0+8]+0.3162277660168379*uSelf[a0+4]*m1rSelf[a0+8]+0.3162277660168379*m1rSelf[a0+4]*uSelf[a0+7]+0.3162277660168379*uSelf[a0+4]*m1rSelf[a0+7]+0.3535533905932737*m1rSelf[a0+5]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+5]*m1rSelf[a0+6]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+4]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+4]+0.3535533905932737*m1rSelf[a0+1]*uSelf[a0+2]+0.3535533905932737*uSelf[a0+1]*m1rSelf[a0+2]; 
    kinESelf[5] += 0.3162277660168379*m1rSelf[a0+5]*uSelf[a0+9]+0.3162277660168379*uSelf[a0+5]*m1rSelf[a0+9]+0.3162277660168379*m1rSelf[a0+5]*uSelf[a0+7]+0.3162277660168379*uSelf[a0+5]*m1rSelf[a0+7]+0.3535533905932737*m1rSelf[a0+4]*uSelf[a0+6]+0.3535533905932737*uSelf[a0+4]*m1rSelf[a0+6]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+5]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+5]+0.3535533905932737*m1rSelf[a0+1]*uSelf[a0+3]+0.3535533905932737*uSelf[a0+1]*m1rSelf[a0+3]; 
    kinESelf[6] += 0.3162277660168379*m1rSelf[a0+6]*uSelf[a0+9]+0.3162277660168379*uSelf[a0+6]*m1rSelf[a0+9]+0.3162277660168379*m1rSelf[a0+6]*uSelf[a0+8]+0.3162277660168379*uSelf[a0+6]*m1rSelf[a0+8]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+6]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+6]+0.3535533905932737*m1rSelf[a0+4]*uSelf[a0+5]+0.3535533905932737*uSelf[a0+4]*m1rSelf[a0+5]+0.3535533905932737*m1rSelf[a0+2]*uSelf[a0+3]+0.3535533905932737*uSelf[a0+2]*m1rSelf[a0+3]; 
    kinESelf[7] += 0.2258769757263128*m1rSelf[a0+7]*uSelf[a0+7]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+7]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+7]+0.3162277660168379*m1rSelf[a0+5]*uSelf[a0+5]+0.3162277660168379*m1rSelf[a0+4]*uSelf[a0+4]+0.3162277660168379*m1rSelf[a0+1]*uSelf[a0+1]; 
    kinESelf[8] += 0.2258769757263128*m1rSelf[a0+8]*uSelf[a0+8]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+8]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+8]+0.3162277660168379*m1rSelf[a0+6]*uSelf[a0+6]+0.3162277660168379*m1rSelf[a0+4]*uSelf[a0+4]+0.3162277660168379*m1rSelf[a0+2]*uSelf[a0+2]; 
    kinESelf[9] += 0.2258769757263128*m1rSelf[a0+9]*uSelf[a0+9]+0.3535533905932737*m1rSelf[a0]*uSelf[a0+9]+0.3535533905932737*uSelf[a0]*m1rSelf[a0+9]+0.3162277660168379*m1rSelf[a0+6]*uSelf[a0+6]+0.3162277660168379*m1rSelf[a0+5]*uSelf[a0+5]+0.3162277660168379*m1rSelf[a0+3]*uSelf[a0+3]; 
  } 
 
  double kinEOther[10]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    kinEOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.3535533905932737*m1rOther[a0+9]*uOther[a0+9]+0.3535533905932737*m1rOther[a0+8]*uOther[a0+8]+0.3535533905932737*m1rOther[a0+7]*uOther[a0+7]+0.3535533905932737*m1rOther[a0+6]*uOther[a0+6]+0.3535533905932737*m1rOther[a0+5]*uOther[a0+5]+0.3535533905932737*m1rOther[a0+4]*uOther[a0+4]+0.3535533905932737*m1rOther[a0+3]*uOther[a0+3]+0.3535533905932737*m1rOther[a0+2]*uOther[a0+2]+0.3535533905932737*m1rOther[a0+1]*uOther[a0+1]+0.3535533905932737*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.3162277660168379*m1rOther[a0+1]*uOther[a0+7]+0.3162277660168379*uOther[a0+1]*m1rOther[a0+7]+0.3535533905932737*m1rOther[a0+3]*uOther[a0+5]+0.3535533905932737*uOther[a0+3]*m1rOther[a0+5]+0.3535533905932737*m1rOther[a0+2]*uOther[a0+4]+0.3535533905932737*uOther[a0+2]*m1rOther[a0+4]+0.3535533905932737*m1rOther[a0]*uOther[a0+1]+0.3535533905932737*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.3162277660168379*m1rOther[a0+2]*uOther[a0+8]+0.3162277660168379*uOther[a0+2]*m1rOther[a0+8]+0.3535533905932737*m1rOther[a0+3]*uOther[a0+6]+0.3535533905932737*uOther[a0+3]*m1rOther[a0+6]+0.3535533905932737*m1rOther[a0+1]*uOther[a0+4]+0.3535533905932737*uOther[a0+1]*m1rOther[a0+4]+0.3535533905932737*m1rOther[a0]*uOther[a0+2]+0.3535533905932737*uOther[a0]*m1rOther[a0+2]; 
    kinEOther[3] += 0.3162277660168379*m1rOther[a0+3]*uOther[a0+9]+0.3162277660168379*uOther[a0+3]*m1rOther[a0+9]+0.3535533905932737*m1rOther[a0+2]*uOther[a0+6]+0.3535533905932737*uOther[a0+2]*m1rOther[a0+6]+0.3535533905932737*m1rOther[a0+1]*uOther[a0+5]+0.3535533905932737*uOther[a0+1]*m1rOther[a0+5]+0.3535533905932737*m1rOther[a0]*uOther[a0+3]+0.3535533905932737*uOther[a0]*m1rOther[a0+3]; 
    kinEOther[4] += 0.3162277660168379*m1rOther[a0+4]*uOther[a0+8]+0.3162277660168379*uOther[a0+4]*m1rOther[a0+8]+0.3162277660168379*m1rOther[a0+4]*uOther[a0+7]+0.3162277660168379*uOther[a0+4]*m1rOther[a0+7]+0.3535533905932737*m1rOther[a0+5]*uOther[a0+6]+0.3535533905932737*uOther[a0+5]*m1rOther[a0+6]+0.3535533905932737*m1rOther[a0]*uOther[a0+4]+0.3535533905932737*uOther[a0]*m1rOther[a0+4]+0.3535533905932737*m1rOther[a0+1]*uOther[a0+2]+0.3535533905932737*uOther[a0+1]*m1rOther[a0+2]; 
    kinEOther[5] += 0.3162277660168379*m1rOther[a0+5]*uOther[a0+9]+0.3162277660168379*uOther[a0+5]*m1rOther[a0+9]+0.3162277660168379*m1rOther[a0+5]*uOther[a0+7]+0.3162277660168379*uOther[a0+5]*m1rOther[a0+7]+0.3535533905932737*m1rOther[a0+4]*uOther[a0+6]+0.3535533905932737*uOther[a0+4]*m1rOther[a0+6]+0.3535533905932737*m1rOther[a0]*uOther[a0+5]+0.3535533905932737*uOther[a0]*m1rOther[a0+5]+0.3535533905932737*m1rOther[a0+1]*uOther[a0+3]+0.3535533905932737*uOther[a0+1]*m1rOther[a0+3]; 
    kinEOther[6] += 0.3162277660168379*m1rOther[a0+6]*uOther[a0+9]+0.3162277660168379*uOther[a0+6]*m1rOther[a0+9]+0.3162277660168379*m1rOther[a0+6]*uOther[a0+8]+0.3162277660168379*uOther[a0+6]*m1rOther[a0+8]+0.3535533905932737*m1rOther[a0]*uOther[a0+6]+0.3535533905932737*uOther[a0]*m1rOther[a0+6]+0.3535533905932737*m1rOther[a0+4]*uOther[a0+5]+0.3535533905932737*uOther[a0+4]*m1rOther[a0+5]+0.3535533905932737*m1rOther[a0+2]*uOther[a0+3]+0.3535533905932737*uOther[a0+2]*m1rOther[a0+3]; 
    kinEOther[7] += 0.2258769757263128*m1rOther[a0+7]*uOther[a0+7]+0.3535533905932737*m1rOther[a0]*uOther[a0+7]+0.3535533905932737*uOther[a0]*m1rOther[a0+7]+0.3162277660168379*m1rOther[a0+5]*uOther[a0+5]+0.3162277660168379*m1rOther[a0+4]*uOther[a0+4]+0.3162277660168379*m1rOther[a0+1]*uOther[a0+1]; 
    kinEOther[8] += 0.2258769757263128*m1rOther[a0+8]*uOther[a0+8]+0.3535533905932737*m1rOther[a0]*uOther[a0+8]+0.3535533905932737*uOther[a0]*m1rOther[a0+8]+0.3162277660168379*m1rOther[a0+6]*uOther[a0+6]+0.3162277660168379*m1rOther[a0+4]*uOther[a0+4]+0.3162277660168379*m1rOther[a0+2]*uOther[a0+2]; 
    kinEOther[9] += 0.2258769757263128*m1rOther[a0+9]*uOther[a0+9]+0.3535533905932737*m1rOther[a0]*uOther[a0+9]+0.3535533905932737*uOther[a0]*m1rOther[a0+9]+0.3162277660168379*m1rOther[a0+6]*uOther[a0+6]+0.3162277660168379*m1rOther[a0+5]*uOther[a0+5]+0.3162277660168379*m1rOther[a0+3]*uOther[a0+3]; 
  } 
 
  double relKinE[10]; 
  // zero out array with dot product of uSelf-uOther and m1EffD. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.3535533905932737*m1EffD[a0+9]*uSelf[a0+9]-0.3535533905932737*m1EffD[a0+9]*uOther[a0+9]+0.3535533905932737*m1EffD[a0+8]*uSelf[a0+8]-0.3535533905932737*m1EffD[a0+8]*uOther[a0+8]+0.3535533905932737*m1EffD[a0+7]*uSelf[a0+7]-0.3535533905932737*m1EffD[a0+7]*uOther[a0+7]+0.3535533905932737*m1EffD[a0+6]*uSelf[a0+6]-0.3535533905932737*m1EffD[a0+6]*uOther[a0+6]+0.3535533905932737*m1EffD[a0+5]*uSelf[a0+5]-0.3535533905932737*m1EffD[a0+5]*uOther[a0+5]+0.3535533905932737*m1EffD[a0+4]*uSelf[a0+4]-0.3535533905932737*m1EffD[a0+4]*uOther[a0+4]+0.3535533905932737*m1EffD[a0+3]*uSelf[a0+3]-0.3535533905932737*m1EffD[a0+3]*uOther[a0+3]+0.3535533905932737*m1EffD[a0+2]*uSelf[a0+2]-0.3535533905932737*m1EffD[a0+2]*uOther[a0+2]+0.3535533905932737*m1EffD[a0+1]*uSelf[a0+1]-0.3535533905932737*m1EffD[a0+1]*uOther[a0+1]+0.3535533905932737*m1EffD[a0]*uSelf[a0]-0.3535533905932737*m1EffD[a0]*uOther[a0]; 
    relKinE[1] += 0.3162277660168379*m1EffD[a0+1]*uSelf[a0+7]-0.3162277660168379*m1EffD[a0+1]*uOther[a0+7]+0.3162277660168379*uSelf[a0+1]*m1EffD[a0+7]-0.3162277660168379*uOther[a0+1]*m1EffD[a0+7]+0.3535533905932737*m1EffD[a0+3]*uSelf[a0+5]-0.3535533905932737*m1EffD[a0+3]*uOther[a0+5]+0.3535533905932737*uSelf[a0+3]*m1EffD[a0+5]-0.3535533905932737*uOther[a0+3]*m1EffD[a0+5]+0.3535533905932737*m1EffD[a0+2]*uSelf[a0+4]-0.3535533905932737*m1EffD[a0+2]*uOther[a0+4]+0.3535533905932737*uSelf[a0+2]*m1EffD[a0+4]-0.3535533905932737*uOther[a0+2]*m1EffD[a0+4]+0.3535533905932737*m1EffD[a0]*uSelf[a0+1]-0.3535533905932737*m1EffD[a0]*uOther[a0+1]+0.3535533905932737*uSelf[a0]*m1EffD[a0+1]-0.3535533905932737*uOther[a0]*m1EffD[a0+1]; 
    relKinE[2] += 0.3162277660168379*m1EffD[a0+2]*uSelf[a0+8]-0.3162277660168379*m1EffD[a0+2]*uOther[a0+8]+0.3162277660168379*uSelf[a0+2]*m1EffD[a0+8]-0.3162277660168379*uOther[a0+2]*m1EffD[a0+8]+0.3535533905932737*m1EffD[a0+3]*uSelf[a0+6]-0.3535533905932737*m1EffD[a0+3]*uOther[a0+6]+0.3535533905932737*uSelf[a0+3]*m1EffD[a0+6]-0.3535533905932737*uOther[a0+3]*m1EffD[a0+6]+0.3535533905932737*m1EffD[a0+1]*uSelf[a0+4]-0.3535533905932737*m1EffD[a0+1]*uOther[a0+4]+0.3535533905932737*uSelf[a0+1]*m1EffD[a0+4]-0.3535533905932737*uOther[a0+1]*m1EffD[a0+4]+0.3535533905932737*m1EffD[a0]*uSelf[a0+2]-0.3535533905932737*m1EffD[a0]*uOther[a0+2]+0.3535533905932737*uSelf[a0]*m1EffD[a0+2]-0.3535533905932737*uOther[a0]*m1EffD[a0+2]; 
    relKinE[3] += 0.3162277660168379*m1EffD[a0+3]*uSelf[a0+9]-0.3162277660168379*m1EffD[a0+3]*uOther[a0+9]+0.3162277660168379*uSelf[a0+3]*m1EffD[a0+9]-0.3162277660168379*uOther[a0+3]*m1EffD[a0+9]+0.3535533905932737*m1EffD[a0+2]*uSelf[a0+6]-0.3535533905932737*m1EffD[a0+2]*uOther[a0+6]+0.3535533905932737*uSelf[a0+2]*m1EffD[a0+6]-0.3535533905932737*uOther[a0+2]*m1EffD[a0+6]+0.3535533905932737*m1EffD[a0+1]*uSelf[a0+5]-0.3535533905932737*m1EffD[a0+1]*uOther[a0+5]+0.3535533905932737*uSelf[a0+1]*m1EffD[a0+5]-0.3535533905932737*uOther[a0+1]*m1EffD[a0+5]+0.3535533905932737*m1EffD[a0]*uSelf[a0+3]-0.3535533905932737*m1EffD[a0]*uOther[a0+3]+0.3535533905932737*uSelf[a0]*m1EffD[a0+3]-0.3535533905932737*uOther[a0]*m1EffD[a0+3]; 
    relKinE[4] += 0.3162277660168379*m1EffD[a0+4]*uSelf[a0+8]-0.3162277660168379*m1EffD[a0+4]*uOther[a0+8]+0.3162277660168379*uSelf[a0+4]*m1EffD[a0+8]-0.3162277660168379*uOther[a0+4]*m1EffD[a0+8]+0.3162277660168379*m1EffD[a0+4]*uSelf[a0+7]-0.3162277660168379*m1EffD[a0+4]*uOther[a0+7]+0.3162277660168379*uSelf[a0+4]*m1EffD[a0+7]-0.3162277660168379*uOther[a0+4]*m1EffD[a0+7]+0.3535533905932737*m1EffD[a0+5]*uSelf[a0+6]-0.3535533905932737*m1EffD[a0+5]*uOther[a0+6]+0.3535533905932737*uSelf[a0+5]*m1EffD[a0+6]-0.3535533905932737*uOther[a0+5]*m1EffD[a0+6]+0.3535533905932737*m1EffD[a0]*uSelf[a0+4]-0.3535533905932737*m1EffD[a0]*uOther[a0+4]+0.3535533905932737*uSelf[a0]*m1EffD[a0+4]-0.3535533905932737*uOther[a0]*m1EffD[a0+4]+0.3535533905932737*m1EffD[a0+1]*uSelf[a0+2]-0.3535533905932737*m1EffD[a0+1]*uOther[a0+2]+0.3535533905932737*uSelf[a0+1]*m1EffD[a0+2]-0.3535533905932737*uOther[a0+1]*m1EffD[a0+2]; 
    relKinE[5] += 0.3162277660168379*m1EffD[a0+5]*uSelf[a0+9]-0.3162277660168379*m1EffD[a0+5]*uOther[a0+9]+0.3162277660168379*uSelf[a0+5]*m1EffD[a0+9]-0.3162277660168379*uOther[a0+5]*m1EffD[a0+9]+0.3162277660168379*m1EffD[a0+5]*uSelf[a0+7]-0.3162277660168379*m1EffD[a0+5]*uOther[a0+7]+0.3162277660168379*uSelf[a0+5]*m1EffD[a0+7]-0.3162277660168379*uOther[a0+5]*m1EffD[a0+7]+0.3535533905932737*m1EffD[a0+4]*uSelf[a0+6]-0.3535533905932737*m1EffD[a0+4]*uOther[a0+6]+0.3535533905932737*uSelf[a0+4]*m1EffD[a0+6]-0.3535533905932737*uOther[a0+4]*m1EffD[a0+6]+0.3535533905932737*m1EffD[a0]*uSelf[a0+5]-0.3535533905932737*m1EffD[a0]*uOther[a0+5]+0.3535533905932737*uSelf[a0]*m1EffD[a0+5]-0.3535533905932737*uOther[a0]*m1EffD[a0+5]+0.3535533905932737*m1EffD[a0+1]*uSelf[a0+3]-0.3535533905932737*m1EffD[a0+1]*uOther[a0+3]+0.3535533905932737*uSelf[a0+1]*m1EffD[a0+3]-0.3535533905932737*uOther[a0+1]*m1EffD[a0+3]; 
    relKinE[6] += 0.3162277660168379*m1EffD[a0+6]*uSelf[a0+9]-0.3162277660168379*m1EffD[a0+6]*uOther[a0+9]+0.3162277660168379*uSelf[a0+6]*m1EffD[a0+9]-0.3162277660168379*uOther[a0+6]*m1EffD[a0+9]+0.3162277660168379*m1EffD[a0+6]*uSelf[a0+8]-0.3162277660168379*m1EffD[a0+6]*uOther[a0+8]+0.3162277660168379*uSelf[a0+6]*m1EffD[a0+8]-0.3162277660168379*uOther[a0+6]*m1EffD[a0+8]+0.3535533905932737*m1EffD[a0]*uSelf[a0+6]-0.3535533905932737*m1EffD[a0]*uOther[a0+6]+0.3535533905932737*uSelf[a0]*m1EffD[a0+6]-0.3535533905932737*uOther[a0]*m1EffD[a0+6]+0.3535533905932737*m1EffD[a0+4]*uSelf[a0+5]-0.3535533905932737*m1EffD[a0+4]*uOther[a0+5]+0.3535533905932737*uSelf[a0+4]*m1EffD[a0+5]-0.3535533905932737*uOther[a0+4]*m1EffD[a0+5]+0.3535533905932737*m1EffD[a0+2]*uSelf[a0+3]-0.3535533905932737*m1EffD[a0+2]*uOther[a0+3]+0.3535533905932737*uSelf[a0+2]*m1EffD[a0+3]-0.3535533905932737*uOther[a0+2]*m1EffD[a0+3]; 
    relKinE[7] += 0.2258769757263128*m1EffD[a0+7]*uSelf[a0+7]+0.3535533905932737*m1EffD[a0]*uSelf[a0+7]-0.2258769757263128*m1EffD[a0+7]*uOther[a0+7]-0.3535533905932737*m1EffD[a0]*uOther[a0+7]+0.3535533905932737*uSelf[a0]*m1EffD[a0+7]-0.3535533905932737*uOther[a0]*m1EffD[a0+7]+0.3162277660168379*m1EffD[a0+5]*uSelf[a0+5]-0.3162277660168379*m1EffD[a0+5]*uOther[a0+5]+0.3162277660168379*m1EffD[a0+4]*uSelf[a0+4]-0.3162277660168379*m1EffD[a0+4]*uOther[a0+4]+0.3162277660168379*m1EffD[a0+1]*uSelf[a0+1]-0.3162277660168379*m1EffD[a0+1]*uOther[a0+1]; 
    relKinE[8] += 0.2258769757263128*m1EffD[a0+8]*uSelf[a0+8]+0.3535533905932737*m1EffD[a0]*uSelf[a0+8]-0.2258769757263128*m1EffD[a0+8]*uOther[a0+8]-0.3535533905932737*m1EffD[a0]*uOther[a0+8]+0.3535533905932737*uSelf[a0]*m1EffD[a0+8]-0.3535533905932737*uOther[a0]*m1EffD[a0+8]+0.3162277660168379*m1EffD[a0+6]*uSelf[a0+6]-0.3162277660168379*m1EffD[a0+6]*uOther[a0+6]+0.3162277660168379*m1EffD[a0+4]*uSelf[a0+4]-0.3162277660168379*m1EffD[a0+4]*uOther[a0+4]+0.3162277660168379*m1EffD[a0+2]*uSelf[a0+2]-0.3162277660168379*m1EffD[a0+2]*uOther[a0+2]; 
    relKinE[9] += 0.2258769757263128*m1EffD[a0+9]*uSelf[a0+9]+0.3535533905932737*m1EffD[a0]*uSelf[a0+9]-0.2258769757263128*m1EffD[a0+9]*uOther[a0+9]-0.3535533905932737*m1EffD[a0]*uOther[a0+9]+0.3535533905932737*uSelf[a0]*m1EffD[a0+9]-0.3535533905932737*uOther[a0]*m1EffD[a0+9]+0.3162277660168379*m1EffD[a0+6]*uSelf[a0+6]-0.3162277660168379*m1EffD[a0+6]*uOther[a0+6]+0.3162277660168379*m1EffD[a0+5]*uSelf[a0+5]-0.3162277660168379*m1EffD[a0+5]*uOther[a0+5]+0.3162277660168379*m1EffD[a0+3]*uSelf[a0+3]-0.3162277660168379*m1EffD[a0+3]*uOther[a0+3]; 
  } 
 
  // Divide m0Other*(m2Self-kinESelf) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Other and m2Self-uSelf.m1Self. 
  double m0OtherThESelf[10]; 
  m0OtherThESelf[0] = 0.3535533905932737*m0rOther[9]*m2rSelf[9]-0.3535533905932737*kinESelf[9]*m0rOther[9]+0.3535533905932737*m0rOther[8]*m2rSelf[8]-0.3535533905932737*kinESelf[8]*m0rOther[8]+0.3535533905932737*m0rOther[7]*m2rSelf[7]-0.3535533905932737*kinESelf[7]*m0rOther[7]+0.3535533905932737*m0rOther[6]*m2rSelf[6]-0.3535533905932737*kinESelf[6]*m0rOther[6]+0.3535533905932737*m0rOther[5]*m2rSelf[5]-0.3535533905932737*kinESelf[5]*m0rOther[5]+0.3535533905932737*m0rOther[4]*m2rSelf[4]-0.3535533905932737*kinESelf[4]*m0rOther[4]+0.3535533905932737*m0rOther[3]*m2rSelf[3]-0.3535533905932737*kinESelf[3]*m0rOther[3]+0.3535533905932737*m0rOther[2]*m2rSelf[2]-0.3535533905932737*kinESelf[2]*m0rOther[2]+0.3535533905932737*m0rOther[1]*m2rSelf[1]-0.3535533905932737*kinESelf[1]*m0rOther[1]+0.3535533905932737*m0rOther[0]*m2rSelf[0]-0.3535533905932737*kinESelf[0]*m0rOther[0]; 
  m0OtherThESelf[1] = 0.3162277660168379*m0rOther[1]*m2rSelf[7]+0.3162277660168379*m2rSelf[1]*m0rOther[7]-0.3162277660168379*kinESelf[1]*m0rOther[7]-0.3162277660168379*m0rOther[1]*kinESelf[7]+0.3535533905932737*m0rOther[3]*m2rSelf[5]+0.3535533905932737*m2rSelf[3]*m0rOther[5]-0.3535533905932737*kinESelf[3]*m0rOther[5]-0.3535533905932737*m0rOther[3]*kinESelf[5]+0.3535533905932737*m0rOther[2]*m2rSelf[4]+0.3535533905932737*m2rSelf[2]*m0rOther[4]-0.3535533905932737*kinESelf[2]*m0rOther[4]-0.3535533905932737*m0rOther[2]*kinESelf[4]+0.3535533905932737*m0rOther[0]*m2rSelf[1]+0.3535533905932737*m2rSelf[0]*m0rOther[1]-0.3535533905932737*kinESelf[0]*m0rOther[1]-0.3535533905932737*m0rOther[0]*kinESelf[1]; 
  m0OtherThESelf[2] = 0.3162277660168379*m0rOther[2]*m2rSelf[8]+0.3162277660168379*m2rSelf[2]*m0rOther[8]-0.3162277660168379*kinESelf[2]*m0rOther[8]-0.3162277660168379*m0rOther[2]*kinESelf[8]+0.3535533905932737*m0rOther[3]*m2rSelf[6]+0.3535533905932737*m2rSelf[3]*m0rOther[6]-0.3535533905932737*kinESelf[3]*m0rOther[6]-0.3535533905932737*m0rOther[3]*kinESelf[6]+0.3535533905932737*m0rOther[1]*m2rSelf[4]+0.3535533905932737*m2rSelf[1]*m0rOther[4]-0.3535533905932737*kinESelf[1]*m0rOther[4]-0.3535533905932737*m0rOther[1]*kinESelf[4]+0.3535533905932737*m0rOther[0]*m2rSelf[2]+0.3535533905932737*m2rSelf[0]*m0rOther[2]-0.3535533905932737*kinESelf[0]*m0rOther[2]-0.3535533905932737*m0rOther[0]*kinESelf[2]; 
  m0OtherThESelf[3] = 0.3162277660168379*m0rOther[3]*m2rSelf[9]+0.3162277660168379*m2rSelf[3]*m0rOther[9]-0.3162277660168379*kinESelf[3]*m0rOther[9]-0.3162277660168379*m0rOther[3]*kinESelf[9]+0.3535533905932737*m0rOther[2]*m2rSelf[6]+0.3535533905932737*m2rSelf[2]*m0rOther[6]-0.3535533905932737*kinESelf[2]*m0rOther[6]-0.3535533905932737*m0rOther[2]*kinESelf[6]+0.3535533905932737*m0rOther[1]*m2rSelf[5]+0.3535533905932737*m2rSelf[1]*m0rOther[5]-0.3535533905932737*kinESelf[1]*m0rOther[5]-0.3535533905932737*m0rOther[1]*kinESelf[5]+0.3535533905932737*m0rOther[0]*m2rSelf[3]+0.3535533905932737*m2rSelf[0]*m0rOther[3]-0.3535533905932737*kinESelf[0]*m0rOther[3]-0.3535533905932737*m0rOther[0]*kinESelf[3]; 
  m0OtherThESelf[4] = 0.3162277660168379*m0rOther[4]*m2rSelf[8]+0.3162277660168379*m2rSelf[4]*m0rOther[8]-0.3162277660168379*kinESelf[4]*m0rOther[8]-0.3162277660168379*m0rOther[4]*kinESelf[8]+0.3162277660168379*m0rOther[4]*m2rSelf[7]+0.3162277660168379*m2rSelf[4]*m0rOther[7]-0.3162277660168379*kinESelf[4]*m0rOther[7]-0.3162277660168379*m0rOther[4]*kinESelf[7]+0.3535533905932737*m0rOther[5]*m2rSelf[6]+0.3535533905932737*m2rSelf[5]*m0rOther[6]-0.3535533905932737*kinESelf[5]*m0rOther[6]-0.3535533905932737*m0rOther[5]*kinESelf[6]+0.3535533905932737*m0rOther[0]*m2rSelf[4]+0.3535533905932737*m2rSelf[0]*m0rOther[4]-0.3535533905932737*kinESelf[0]*m0rOther[4]-0.3535533905932737*m0rOther[0]*kinESelf[4]+0.3535533905932737*m0rOther[1]*m2rSelf[2]+0.3535533905932737*m2rSelf[1]*m0rOther[2]-0.3535533905932737*kinESelf[1]*m0rOther[2]-0.3535533905932737*m0rOther[1]*kinESelf[2]; 
  m0OtherThESelf[5] = 0.3162277660168379*m0rOther[5]*m2rSelf[9]+0.3162277660168379*m2rSelf[5]*m0rOther[9]-0.3162277660168379*kinESelf[5]*m0rOther[9]-0.3162277660168379*m0rOther[5]*kinESelf[9]+0.3162277660168379*m0rOther[5]*m2rSelf[7]+0.3162277660168379*m2rSelf[5]*m0rOther[7]-0.3162277660168379*kinESelf[5]*m0rOther[7]-0.3162277660168379*m0rOther[5]*kinESelf[7]+0.3535533905932737*m0rOther[4]*m2rSelf[6]+0.3535533905932737*m2rSelf[4]*m0rOther[6]-0.3535533905932737*kinESelf[4]*m0rOther[6]-0.3535533905932737*m0rOther[4]*kinESelf[6]+0.3535533905932737*m0rOther[0]*m2rSelf[5]+0.3535533905932737*m2rSelf[0]*m0rOther[5]-0.3535533905932737*kinESelf[0]*m0rOther[5]-0.3535533905932737*m0rOther[0]*kinESelf[5]+0.3535533905932737*m0rOther[1]*m2rSelf[3]+0.3535533905932737*m2rSelf[1]*m0rOther[3]-0.3535533905932737*kinESelf[1]*m0rOther[3]-0.3535533905932737*m0rOther[1]*kinESelf[3]; 
  m0OtherThESelf[6] = 0.3162277660168379*m0rOther[6]*m2rSelf[9]+0.3162277660168379*m2rSelf[6]*m0rOther[9]-0.3162277660168379*kinESelf[6]*m0rOther[9]-0.3162277660168379*m0rOther[6]*kinESelf[9]+0.3162277660168379*m0rOther[6]*m2rSelf[8]+0.3162277660168379*m2rSelf[6]*m0rOther[8]-0.3162277660168379*kinESelf[6]*m0rOther[8]-0.3162277660168379*m0rOther[6]*kinESelf[8]+0.3535533905932737*m0rOther[0]*m2rSelf[6]+0.3535533905932737*m2rSelf[0]*m0rOther[6]-0.3535533905932737*kinESelf[0]*m0rOther[6]-0.3535533905932737*m0rOther[0]*kinESelf[6]+0.3535533905932737*m0rOther[4]*m2rSelf[5]+0.3535533905932737*m2rSelf[4]*m0rOther[5]-0.3535533905932737*kinESelf[4]*m0rOther[5]-0.3535533905932737*m0rOther[4]*kinESelf[5]+0.3535533905932737*m0rOther[2]*m2rSelf[3]+0.3535533905932737*m2rSelf[2]*m0rOther[3]-0.3535533905932737*kinESelf[2]*m0rOther[3]-0.3535533905932737*m0rOther[2]*kinESelf[3]; 
  m0OtherThESelf[7] = 0.2258769757263128*m0rOther[7]*m2rSelf[7]+0.3535533905932737*m0rOther[0]*m2rSelf[7]-0.2258769757263128*kinESelf[7]*m0rOther[7]+0.3535533905932737*m2rSelf[0]*m0rOther[7]-0.3535533905932737*kinESelf[0]*m0rOther[7]-0.3535533905932737*m0rOther[0]*kinESelf[7]+0.3162277660168379*m0rOther[5]*m2rSelf[5]-0.3162277660168379*kinESelf[5]*m0rOther[5]+0.3162277660168379*m0rOther[4]*m2rSelf[4]-0.3162277660168379*kinESelf[4]*m0rOther[4]+0.3162277660168379*m0rOther[1]*m2rSelf[1]-0.3162277660168379*kinESelf[1]*m0rOther[1]; 
  m0OtherThESelf[8] = 0.2258769757263128*m0rOther[8]*m2rSelf[8]+0.3535533905932737*m0rOther[0]*m2rSelf[8]-0.2258769757263128*kinESelf[8]*m0rOther[8]+0.3535533905932737*m2rSelf[0]*m0rOther[8]-0.3535533905932737*kinESelf[0]*m0rOther[8]-0.3535533905932737*m0rOther[0]*kinESelf[8]+0.3162277660168379*m0rOther[6]*m2rSelf[6]-0.3162277660168379*kinESelf[6]*m0rOther[6]+0.3162277660168379*m0rOther[4]*m2rSelf[4]-0.3162277660168379*kinESelf[4]*m0rOther[4]+0.3162277660168379*m0rOther[2]*m2rSelf[2]-0.3162277660168379*kinESelf[2]*m0rOther[2]; 
  m0OtherThESelf[9] = 0.2258769757263128*m0rOther[9]*m2rSelf[9]+0.3535533905932737*m0rOther[0]*m2rSelf[9]-0.2258769757263128*kinESelf[9]*m0rOther[9]+0.3535533905932737*m2rSelf[0]*m0rOther[9]-0.3535533905932737*kinESelf[0]*m0rOther[9]-0.3535533905932737*m0rOther[0]*kinESelf[9]+0.3162277660168379*m0rOther[6]*m2rSelf[6]-0.3162277660168379*kinESelf[6]*m0rOther[6]+0.3162277660168379*m0rOther[5]*m2rSelf[5]-0.3162277660168379*kinESelf[5]*m0rOther[5]+0.3162277660168379*m0rOther[3]*m2rSelf[3]-0.3162277660168379*kinESelf[3]*m0rOther[3]; 
  dataDiv->BEV_S << m0OtherThESelf[0],m0OtherThESelf[1],m0OtherThESelf[2],m0OtherThESelf[3],m0OtherThESelf[4],m0OtherThESelf[5],m0OtherThESelf[6],m0OtherThESelf[7],m0OtherThESelf[8],m0OtherThESelf[9]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthSelf[10]; 
  Eigen::Map<VectorXd>(effEthSelf,10,1) = dataDiv->u_S; 
 
  // Divide m0Self*(m2Other-kinEOther) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Self and m2Other-uOther.m1Other. 
  double m0SelfThEOther[10]; 
  m0SelfThEOther[0] = 0.3535533905932737*m0rSelf[9]*m2rOther[9]-0.3535533905932737*kinEOther[9]*m0rSelf[9]+0.3535533905932737*m0rSelf[8]*m2rOther[8]-0.3535533905932737*kinEOther[8]*m0rSelf[8]+0.3535533905932737*m0rSelf[7]*m2rOther[7]-0.3535533905932737*kinEOther[7]*m0rSelf[7]+0.3535533905932737*m0rSelf[6]*m2rOther[6]-0.3535533905932737*kinEOther[6]*m0rSelf[6]+0.3535533905932737*m0rSelf[5]*m2rOther[5]-0.3535533905932737*kinEOther[5]*m0rSelf[5]+0.3535533905932737*m0rSelf[4]*m2rOther[4]-0.3535533905932737*kinEOther[4]*m0rSelf[4]+0.3535533905932737*m0rSelf[3]*m2rOther[3]-0.3535533905932737*kinEOther[3]*m0rSelf[3]+0.3535533905932737*m0rSelf[2]*m2rOther[2]-0.3535533905932737*kinEOther[2]*m0rSelf[2]+0.3535533905932737*m0rSelf[1]*m2rOther[1]-0.3535533905932737*kinEOther[1]*m0rSelf[1]+0.3535533905932737*m0rSelf[0]*m2rOther[0]-0.3535533905932737*kinEOther[0]*m0rSelf[0]; 
  m0SelfThEOther[1] = 0.3162277660168379*m0rSelf[1]*m2rOther[7]+0.3162277660168379*m2rOther[1]*m0rSelf[7]-0.3162277660168379*kinEOther[1]*m0rSelf[7]-0.3162277660168379*m0rSelf[1]*kinEOther[7]+0.3535533905932737*m0rSelf[3]*m2rOther[5]+0.3535533905932737*m2rOther[3]*m0rSelf[5]-0.3535533905932737*kinEOther[3]*m0rSelf[5]-0.3535533905932737*m0rSelf[3]*kinEOther[5]+0.3535533905932737*m0rSelf[2]*m2rOther[4]+0.3535533905932737*m2rOther[2]*m0rSelf[4]-0.3535533905932737*kinEOther[2]*m0rSelf[4]-0.3535533905932737*m0rSelf[2]*kinEOther[4]+0.3535533905932737*m0rSelf[0]*m2rOther[1]+0.3535533905932737*m2rOther[0]*m0rSelf[1]-0.3535533905932737*kinEOther[0]*m0rSelf[1]-0.3535533905932737*m0rSelf[0]*kinEOther[1]; 
  m0SelfThEOther[2] = 0.3162277660168379*m0rSelf[2]*m2rOther[8]+0.3162277660168379*m2rOther[2]*m0rSelf[8]-0.3162277660168379*kinEOther[2]*m0rSelf[8]-0.3162277660168379*m0rSelf[2]*kinEOther[8]+0.3535533905932737*m0rSelf[3]*m2rOther[6]+0.3535533905932737*m2rOther[3]*m0rSelf[6]-0.3535533905932737*kinEOther[3]*m0rSelf[6]-0.3535533905932737*m0rSelf[3]*kinEOther[6]+0.3535533905932737*m0rSelf[1]*m2rOther[4]+0.3535533905932737*m2rOther[1]*m0rSelf[4]-0.3535533905932737*kinEOther[1]*m0rSelf[4]-0.3535533905932737*m0rSelf[1]*kinEOther[4]+0.3535533905932737*m0rSelf[0]*m2rOther[2]+0.3535533905932737*m2rOther[0]*m0rSelf[2]-0.3535533905932737*kinEOther[0]*m0rSelf[2]-0.3535533905932737*m0rSelf[0]*kinEOther[2]; 
  m0SelfThEOther[3] = 0.3162277660168379*m0rSelf[3]*m2rOther[9]+0.3162277660168379*m2rOther[3]*m0rSelf[9]-0.3162277660168379*kinEOther[3]*m0rSelf[9]-0.3162277660168379*m0rSelf[3]*kinEOther[9]+0.3535533905932737*m0rSelf[2]*m2rOther[6]+0.3535533905932737*m2rOther[2]*m0rSelf[6]-0.3535533905932737*kinEOther[2]*m0rSelf[6]-0.3535533905932737*m0rSelf[2]*kinEOther[6]+0.3535533905932737*m0rSelf[1]*m2rOther[5]+0.3535533905932737*m2rOther[1]*m0rSelf[5]-0.3535533905932737*kinEOther[1]*m0rSelf[5]-0.3535533905932737*m0rSelf[1]*kinEOther[5]+0.3535533905932737*m0rSelf[0]*m2rOther[3]+0.3535533905932737*m2rOther[0]*m0rSelf[3]-0.3535533905932737*kinEOther[0]*m0rSelf[3]-0.3535533905932737*m0rSelf[0]*kinEOther[3]; 
  m0SelfThEOther[4] = 0.3162277660168379*m0rSelf[4]*m2rOther[8]+0.3162277660168379*m2rOther[4]*m0rSelf[8]-0.3162277660168379*kinEOther[4]*m0rSelf[8]-0.3162277660168379*m0rSelf[4]*kinEOther[8]+0.3162277660168379*m0rSelf[4]*m2rOther[7]+0.3162277660168379*m2rOther[4]*m0rSelf[7]-0.3162277660168379*kinEOther[4]*m0rSelf[7]-0.3162277660168379*m0rSelf[4]*kinEOther[7]+0.3535533905932737*m0rSelf[5]*m2rOther[6]+0.3535533905932737*m2rOther[5]*m0rSelf[6]-0.3535533905932737*kinEOther[5]*m0rSelf[6]-0.3535533905932737*m0rSelf[5]*kinEOther[6]+0.3535533905932737*m0rSelf[0]*m2rOther[4]+0.3535533905932737*m2rOther[0]*m0rSelf[4]-0.3535533905932737*kinEOther[0]*m0rSelf[4]-0.3535533905932737*m0rSelf[0]*kinEOther[4]+0.3535533905932737*m0rSelf[1]*m2rOther[2]+0.3535533905932737*m2rOther[1]*m0rSelf[2]-0.3535533905932737*kinEOther[1]*m0rSelf[2]-0.3535533905932737*m0rSelf[1]*kinEOther[2]; 
  m0SelfThEOther[5] = 0.3162277660168379*m0rSelf[5]*m2rOther[9]+0.3162277660168379*m2rOther[5]*m0rSelf[9]-0.3162277660168379*kinEOther[5]*m0rSelf[9]-0.3162277660168379*m0rSelf[5]*kinEOther[9]+0.3162277660168379*m0rSelf[5]*m2rOther[7]+0.3162277660168379*m2rOther[5]*m0rSelf[7]-0.3162277660168379*kinEOther[5]*m0rSelf[7]-0.3162277660168379*m0rSelf[5]*kinEOther[7]+0.3535533905932737*m0rSelf[4]*m2rOther[6]+0.3535533905932737*m2rOther[4]*m0rSelf[6]-0.3535533905932737*kinEOther[4]*m0rSelf[6]-0.3535533905932737*m0rSelf[4]*kinEOther[6]+0.3535533905932737*m0rSelf[0]*m2rOther[5]+0.3535533905932737*m2rOther[0]*m0rSelf[5]-0.3535533905932737*kinEOther[0]*m0rSelf[5]-0.3535533905932737*m0rSelf[0]*kinEOther[5]+0.3535533905932737*m0rSelf[1]*m2rOther[3]+0.3535533905932737*m2rOther[1]*m0rSelf[3]-0.3535533905932737*kinEOther[1]*m0rSelf[3]-0.3535533905932737*m0rSelf[1]*kinEOther[3]; 
  m0SelfThEOther[6] = 0.3162277660168379*m0rSelf[6]*m2rOther[9]+0.3162277660168379*m2rOther[6]*m0rSelf[9]-0.3162277660168379*kinEOther[6]*m0rSelf[9]-0.3162277660168379*m0rSelf[6]*kinEOther[9]+0.3162277660168379*m0rSelf[6]*m2rOther[8]+0.3162277660168379*m2rOther[6]*m0rSelf[8]-0.3162277660168379*kinEOther[6]*m0rSelf[8]-0.3162277660168379*m0rSelf[6]*kinEOther[8]+0.3535533905932737*m0rSelf[0]*m2rOther[6]+0.3535533905932737*m2rOther[0]*m0rSelf[6]-0.3535533905932737*kinEOther[0]*m0rSelf[6]-0.3535533905932737*m0rSelf[0]*kinEOther[6]+0.3535533905932737*m0rSelf[4]*m2rOther[5]+0.3535533905932737*m2rOther[4]*m0rSelf[5]-0.3535533905932737*kinEOther[4]*m0rSelf[5]-0.3535533905932737*m0rSelf[4]*kinEOther[5]+0.3535533905932737*m0rSelf[2]*m2rOther[3]+0.3535533905932737*m2rOther[2]*m0rSelf[3]-0.3535533905932737*kinEOther[2]*m0rSelf[3]-0.3535533905932737*m0rSelf[2]*kinEOther[3]; 
  m0SelfThEOther[7] = 0.2258769757263128*m0rSelf[7]*m2rOther[7]+0.3535533905932737*m0rSelf[0]*m2rOther[7]-0.2258769757263128*kinEOther[7]*m0rSelf[7]+0.3535533905932737*m2rOther[0]*m0rSelf[7]-0.3535533905932737*kinEOther[0]*m0rSelf[7]-0.3535533905932737*m0rSelf[0]*kinEOther[7]+0.3162277660168379*m0rSelf[5]*m2rOther[5]-0.3162277660168379*kinEOther[5]*m0rSelf[5]+0.3162277660168379*m0rSelf[4]*m2rOther[4]-0.3162277660168379*kinEOther[4]*m0rSelf[4]+0.3162277660168379*m0rSelf[1]*m2rOther[1]-0.3162277660168379*kinEOther[1]*m0rSelf[1]; 
  m0SelfThEOther[8] = 0.2258769757263128*m0rSelf[8]*m2rOther[8]+0.3535533905932737*m0rSelf[0]*m2rOther[8]-0.2258769757263128*kinEOther[8]*m0rSelf[8]+0.3535533905932737*m2rOther[0]*m0rSelf[8]-0.3535533905932737*kinEOther[0]*m0rSelf[8]-0.3535533905932737*m0rSelf[0]*kinEOther[8]+0.3162277660168379*m0rSelf[6]*m2rOther[6]-0.3162277660168379*kinEOther[6]*m0rSelf[6]+0.3162277660168379*m0rSelf[4]*m2rOther[4]-0.3162277660168379*kinEOther[4]*m0rSelf[4]+0.3162277660168379*m0rSelf[2]*m2rOther[2]-0.3162277660168379*kinEOther[2]*m0rSelf[2]; 
  m0SelfThEOther[9] = 0.2258769757263128*m0rSelf[9]*m2rOther[9]+0.3535533905932737*m0rSelf[0]*m2rOther[9]-0.2258769757263128*kinEOther[9]*m0rSelf[9]+0.3535533905932737*m2rOther[0]*m0rSelf[9]-0.3535533905932737*kinEOther[0]*m0rSelf[9]-0.3535533905932737*m0rSelf[0]*kinEOther[9]+0.3162277660168379*m0rSelf[6]*m2rOther[6]-0.3162277660168379*kinEOther[6]*m0rSelf[6]+0.3162277660168379*m0rSelf[5]*m2rOther[5]-0.3162277660168379*kinEOther[5]*m0rSelf[5]+0.3162277660168379*m0rSelf[3]*m2rOther[3]-0.3162277660168379*kinEOther[3]*m0rSelf[3]; 
  dataDiv->BEV_S << m0SelfThEOther[0],m0SelfThEOther[1],m0SelfThEOther[2],m0SelfThEOther[3],m0SelfThEOther[4],m0SelfThEOther[5],m0SelfThEOther[6],m0SelfThEOther[7],m0SelfThEOther[8],m0SelfThEOther[9]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthOther[10]; 
  Eigen::Map<VectorXd>(effEthOther,10,1) = dataDiv->u_S; 
 
  double m2Relax[10]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(1.0*relKinE[0]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[0]*mSelf)/(mSelf+mOther)+(1.0*relKinE[0]*mOther)/(mSelf+mOther)+(2.0*effEthOther[0]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2rOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(1.0*relKinE[1]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[1]*mSelf)/(mSelf+mOther)+(1.0*relKinE[1]*mOther)/(mSelf+mOther)+(2.0*effEthOther[1]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2rOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(1.0*relKinE[2]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[2]*mSelf)/(mSelf+mOther)+(1.0*relKinE[2]*mOther)/(mSelf+mOther)+(2.0*effEthOther[2]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2rOther[2])*mnuOther; 
  m2Relax[3] = betaGreenep1*((-(1.0*relKinE[3]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[3]*mSelf)/(mSelf+mOther)+(1.0*relKinE[3]*mOther)/(mSelf+mOther)+(2.0*effEthOther[3]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[3]-1.0*kinESelf[3])*mnuSelf+(kinEOther[3]-1.0*m2rOther[3])*mnuOther; 
  m2Relax[4] = betaGreenep1*((-(1.0*relKinE[4]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[4]*mSelf)/(mSelf+mOther)+(1.0*relKinE[4]*mOther)/(mSelf+mOther)+(2.0*effEthOther[4]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[4]-1.0*kinESelf[4])*mnuSelf+(kinEOther[4]-1.0*m2rOther[4])*mnuOther; 
  m2Relax[5] = betaGreenep1*((-(1.0*relKinE[5]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[5]*mSelf)/(mSelf+mOther)+(1.0*relKinE[5]*mOther)/(mSelf+mOther)+(2.0*effEthOther[5]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[5]-1.0*kinESelf[5])*mnuSelf+(kinEOther[5]-1.0*m2rOther[5])*mnuOther; 
  m2Relax[6] = betaGreenep1*((-(1.0*relKinE[6]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[6]*mSelf)/(mSelf+mOther)+(1.0*relKinE[6]*mOther)/(mSelf+mOther)+(2.0*effEthOther[6]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[6]-1.0*kinESelf[6])*mnuSelf+(kinEOther[6]-1.0*m2rOther[6])*mnuOther; 
  m2Relax[7] = betaGreenep1*((-(1.0*relKinE[7]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[7]*mSelf)/(mSelf+mOther)+(1.0*relKinE[7]*mOther)/(mSelf+mOther)+(2.0*effEthOther[7]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[7]-1.0*kinESelf[7])*mnuSelf+(kinEOther[7]-1.0*m2rOther[7])*mnuOther; 
  m2Relax[8] = betaGreenep1*((-(1.0*relKinE[8]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[8]*mSelf)/(mSelf+mOther)+(1.0*relKinE[8]*mOther)/(mSelf+mOther)+(2.0*effEthOther[8]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[8]-1.0*kinESelf[8])*mnuSelf+(kinEOther[8]-1.0*m2rOther[8])*mnuOther; 
  m2Relax[9] = betaGreenep1*((-(1.0*relKinE[9]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[9]*mSelf)/(mSelf+mOther)+(1.0*relKinE[9]*mOther)/(mSelf+mOther)+(2.0*effEthOther[9]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[9]-1.0*kinESelf[9])*mnuSelf+(kinEOther[9]-1.0*m2rOther[9])*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<10,20>(40,10).setZero(); 
  data->AEM_S.block<20,10>(50,0).setZero(); 
  data->AEM_S.block<10,10>(50,20).setZero(); 
  data->AEM_S.block<10,10>(60,10).setZero(); 
  data->AEM_S.block<10,20>(40,50).setZero(); 
  data->AEM_S.block<20,10>(50,40).setZero(); 
  data->AEM_S.block<10,10>(50,60).setZero(); 
  data->AEM_S.block<10,10>(60,50).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM1sum[10],mnuM1sum[11],mnuM1sum[12],mnuM1sum[13],mnuM1sum[14],mnuM1sum[15],mnuM1sum[16],mnuM1sum[17],mnuM1sum[18],mnuM1sum[19],mnuM1sum[20],mnuM1sum[21],mnuM1sum[22],mnuM1sum[23],mnuM1sum[24],mnuM1sum[25],mnuM1sum[26],mnuM1sum[27],mnuM1sum[28],mnuM1sum[29],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],mnuM2sum[6],mnuM2sum[7],mnuM2sum[8],mnuM2sum[9],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m1Relax[9],m1Relax[10],m1Relax[11],m1Relax[12],m1Relax[13],m1Relax[14],m1Relax[15],m1Relax[16],m1Relax[17],m1Relax[18],m1Relax[19],m1Relax[20],m1Relax[21],m1Relax[22],m1Relax[23],m1Relax[24],m1Relax[25],m1Relax[26],m1Relax[27],m1Relax[28],m1Relax[29],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5],m2Relax[6],m2Relax[7],m2Relax[8],m2Relax[9]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,30,1) = data->u_S.segment<30>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,10,1) = data->u_S.segment<10>(30); 
 
  Eigen::Map<VectorXd>(uCrossOther,30,1) = data->u_S.segment<30>(40); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,10,1) = data->u_S.segment<10>(70); 
 
} 
 
