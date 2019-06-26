#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmLBOCrossPrimMoments2x3vSer_P1(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  if (1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.5*m0Self[3])-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.5*m0Self[3])-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
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
 
  if (1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.5*m0Other[3])-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.5*m0Other[3])-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
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
  data->AEM_S(0,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(0,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(0,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(0,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(1,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(2,2) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,1) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,2) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,3) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,12) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,13) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,14) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,15) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,12) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,13) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,14) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,15) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,12) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,13) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,14) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,15) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,12) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,13) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,14) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,15) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,16) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,17) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,18) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,19) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,16) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,17) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,18) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,19) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,16) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,17) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,18) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,19) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,16) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,17) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,18) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,19) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,28) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,29) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,30) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,31) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,28) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,29) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,30) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,31) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,28) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,29) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,30) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,31) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,28) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,29) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,30) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,31) = -0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(12,0) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(12,1) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(12,3) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(13,0) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(13,1) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(13,2) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(13,3) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(14,0) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(14,1) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(14,2) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(14,3) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(15,0) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(15,1) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(15,2) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(15,3) = 0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(12,16) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(12,17) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(12,18) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(12,19) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(13,16) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(13,17) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(13,18) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(13,19) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(14,16) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(14,17) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(14,18) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(14,19) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(15,16) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(15,17) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(15,18) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(15,19) = 0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(4,4) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,5) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,6) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,7) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(5,4) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,5) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(5,6) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(5,7) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,4) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,5) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(6,6) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,4) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,5) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,6) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,7) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(4,12) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,13) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(4,14) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(4,15) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,12) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,13) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(5,14) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,15) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,12) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,13) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(6,14) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(6,15) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(7,12) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,13) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(7,14) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(7,15) = -0.5*cMSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(4,20) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(4,21) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(4,22) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(4,23) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(5,20) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(5,21) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,22) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(5,23) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(6,20) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(6,21) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(6,22) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,23) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,20) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(7,21) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(7,22) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,23) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(4,28) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,29) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(4,30) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(4,31) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(5,28) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,29) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(5,30) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(5,31) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,28) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,29) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(6,30) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(6,31) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(7,28) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,29) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(7,30) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(7,31) = -0.5*cMOther[4]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(12,4) = 0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(12,5) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(12,6) = 0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(12,7) = 0.5*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(13,4) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(13,5) = 0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(13,6) = 0.5*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(13,7) = 0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(14,4) = 0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(14,5) = 0.5*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(14,6) = 0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(14,7) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(15,4) = 0.5*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(15,5) = 0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(15,6) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(15,7) = 0.5*m1SrSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(12,20) = 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(12,21) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(12,22) = 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(12,23) = 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(13,20) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(13,21) = 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(13,22) = 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(13,23) = 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(14,20) = 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(14,21) = 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(14,22) = 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(14,23) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(15,20) = 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(15,21) = 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(15,22) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(15,23) = 0.5*m1SrOther[4]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(8,8) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,9) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,10) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,11) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,8) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,9) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(9,10) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,11) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,9) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(10,10) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(10,11) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,8) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,9) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,10) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,11) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(8,12) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,13) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(8,14) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(8,15) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(9,12) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,13) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(9,14) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(9,15) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,12) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,13) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(10,14) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(10,15) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(11,12) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,13) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(11,14) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(11,15) = -0.5*cMSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(8,24) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,25) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(8,26) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,27) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,24) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(9,25) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,26) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,27) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,24) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,25) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(10,26) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(10,27) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,24) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(11,25) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(11,26) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,27) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(8,28) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(8,29) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(8,30) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(8,31) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(9,28) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(9,29) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(9,30) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(9,31) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(10,28) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(10,29) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(10,30) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(10,31) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(11,28) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(11,29) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(11,30) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(11,31) = -0.5*cMOther[8]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfZ and uCrossSelfZ ... // 
  data->AEM_S(12,8) = 0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(12,9) = 0.5*m1SrSelf[9]*mnuSelf; 
  data->AEM_S(12,10) = 0.5*m1SrSelf[10]*mnuSelf; 
  data->AEM_S(12,11) = 0.5*m1SrSelf[11]*mnuSelf; 
  data->AEM_S(13,8) = 0.5*m1SrSelf[9]*mnuSelf; 
  data->AEM_S(13,9) = 0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(13,10) = 0.5*m1SrSelf[11]*mnuSelf; 
  data->AEM_S(13,11) = 0.5*m1SrSelf[10]*mnuSelf; 
  data->AEM_S(14,8) = 0.5*m1SrSelf[10]*mnuSelf; 
  data->AEM_S(14,9) = 0.5*m1SrSelf[11]*mnuSelf; 
  data->AEM_S(14,10) = 0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(14,11) = 0.5*m1SrSelf[9]*mnuSelf; 
  data->AEM_S(15,8) = 0.5*m1SrSelf[11]*mnuSelf; 
  data->AEM_S(15,9) = 0.5*m1SrSelf[10]*mnuSelf; 
  data->AEM_S(15,10) = 0.5*m1SrSelf[9]*mnuSelf; 
  data->AEM_S(15,11) = 0.5*m1SrSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherZ and uCrossOtherZ ... // 
  data->AEM_S(12,24) = 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(12,25) = 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(12,26) = 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(12,27) = 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(13,24) = 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(13,25) = 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(13,26) = 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(13,27) = 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(14,24) = 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(14,25) = 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(14,26) = 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(14,27) = 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(15,24) = 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(15,25) = 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(15,26) = 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(15,27) = 0.5*m1SrOther[8]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
  mnuM1sum[9] += m1rSelf[9]*mnuSelf+m1rOther[9]*mnuOther; 
  mnuM1sum[10] += m1rSelf[10]*mnuSelf+m1rOther[10]*mnuOther; 
  mnuM1sum[11] += m1rSelf[11]*mnuSelf+m1rOther[11]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(12,12) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(12,14) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(12,15) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(13,12) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(13,15) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(14,12) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(14,13) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(14,14) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(14,15) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(15,12) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(15,13) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(15,14) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(15,15) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(12,28) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(12,29) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(12,30) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(12,31) = 0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(13,28) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(13,29) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(13,30) = 0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(13,31) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(14,28) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(14,29) = 0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(14,30) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(14,31) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(15,28) = 0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(15,29) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(15,30) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(15,31) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
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
  data->AEM_S(16,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(16,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(17,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(18,2) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,1) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,2) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,3) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(16,12) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(16,13) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(16,14) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(16,15) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(17,12) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(17,13) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(17,14) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(17,15) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(18,12) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(18,13) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(18,14) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(18,15) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(19,12) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(19,13) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(19,14) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(19,15) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(16,16) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(16,17) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(16,18) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(16,19) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(17,16) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(17,17) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(17,18) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(17,19) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(18,16) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(18,17) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(18,18) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(18,19) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(19,16) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(19,17) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(19,18) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(19,19) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(16,28) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(16,29) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(16,30) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(16,31) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(17,28) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(17,29) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(17,30) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(17,31) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(18,28) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(18,29) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(18,30) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(18,31) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(19,28) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(19,29) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(19,30) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(19,31) = 0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(28,0) = (-0.25*m0rSelf[3]*uSelf[3]*mnuSelf)-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(28,1) = (-0.25*m0rSelf[2]*uSelf[3]*mnuSelf)-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(28,2) = (-0.25*m0rSelf[1]*uSelf[3]*mnuSelf)-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(28,3) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,0) = (-0.25*m0rSelf[2]*uSelf[3]*mnuSelf)-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(29,1) = (-0.45*m0rSelf[3]*uSelf[3]*mnuSelf)-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(29,2) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,3) = (-0.45*m0rSelf[1]*uSelf[3]*mnuSelf)-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,0) = (-0.25*m0rSelf[1]*uSelf[3]*mnuSelf)-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,1) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,2) = (-0.45*m0rSelf[3]*uSelf[3]*mnuSelf)-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(30,3) = (-0.45*m0rSelf[2]*uSelf[3]*mnuSelf)-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,0) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,1) = (-0.45*m0rSelf[1]*uSelf[3]*mnuSelf)-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,2) = (-0.45*m0rSelf[2]*uSelf[3]*mnuSelf)-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,3) = (-0.81*m0rSelf[3]*uSelf[3]*mnuSelf)-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(28,16) = 0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(28,17) = 0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(28,18) = 0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(28,19) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(29,16) = 0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(29,17) = 0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(29,18) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(29,19) = 0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,16) = 0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,17) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,18) = 0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(30,19) = 0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(31,16) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,17) = 0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,18) = 0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(31,19) = 0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfX-m0Self*m1OtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[0] = 0.5*m0rOther[3]*m1rSelf[3]-0.5*m0rSelf[3]*m1rOther[3]+0.5*m0rOther[2]*m1rSelf[2]-0.5*m0rSelf[2]*m1rOther[2]+0.5*m0rOther[1]*m1rSelf[1]-0.5*m0rSelf[1]*m1rOther[1]+0.5*m0rOther[0]*m1rSelf[0]-0.5*m0rSelf[0]*m1rOther[0]; 
  m1EffD[1] = 0.5*m0rOther[2]*m1rSelf[3]-0.5*m0rSelf[2]*m1rOther[3]-0.5*m1rOther[2]*m0rSelf[3]+0.5*m1rSelf[2]*m0rOther[3]+0.5*m0rOther[0]*m1rSelf[1]-0.5*m0rSelf[0]*m1rOther[1]-0.5*m1rOther[0]*m0rSelf[1]+0.5*m1rSelf[0]*m0rOther[1]; 
  m1EffD[2] = 0.5*m0rOther[1]*m1rSelf[3]-0.5*m0rSelf[1]*m1rOther[3]-0.5*m1rOther[1]*m0rSelf[3]+0.5*m1rSelf[1]*m0rOther[3]+0.5*m0rOther[0]*m1rSelf[2]-0.5*m0rSelf[0]*m1rOther[2]-0.5*m1rOther[0]*m0rSelf[2]+0.5*m1rSelf[0]*m0rOther[2]; 
  m1EffD[3] = 0.5*m0rOther[0]*m1rSelf[3]-0.5*m0rSelf[0]*m1rOther[3]-0.5*m1rOther[0]*m0rSelf[3]+0.5*m1rSelf[0]*m0rOther[3]+0.5*m0rOther[1]*m1rSelf[2]-0.5*m0rSelf[1]*m1rOther[2]-0.5*m1rOther[1]*m0rSelf[2]+0.5*m1rSelf[1]*m0rOther[2]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(4,4); 
  dataDiv->AEM_S(0,0) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(0,2) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(0,3) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(1,2) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(1,3) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,1) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(2,2) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(2,3) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(3,0) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,1) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(3,2) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(3,3) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
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
  data->AEM_S(20,4) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,5) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,6) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,7) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,4) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,5) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(21,6) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,7) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,4) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,5) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,6) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(22,7) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,4) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,5) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,6) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,7) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(20,12) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(20,13) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(20,14) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(20,15) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(21,12) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(21,13) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(21,14) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(21,15) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(22,12) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(22,13) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(22,14) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(22,15) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(23,12) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(23,13) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(23,14) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(23,15) = -0.5*cMSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(20,20) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(20,21) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(20,22) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(20,23) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(21,20) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(21,21) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(21,22) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(21,23) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(22,20) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(22,21) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(22,22) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(22,23) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(23,20) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(23,21) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(23,22) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(23,23) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(20,28) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(20,29) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(20,30) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(20,31) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(21,28) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(21,29) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(21,30) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(21,31) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(22,28) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(22,29) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(22,30) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(22,31) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(23,28) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(23,29) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(23,30) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(23,31) = 0.5*cMOther[4]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(28,4) = (-0.25*m0rSelf[3]*uSelf[7]*mnuSelf)-0.25*m0rSelf[2]*uSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(28,5) = (-0.25*m0rSelf[2]*uSelf[7]*mnuSelf)-0.25*m0rSelf[3]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[4]*mnuSelf; 
  data->AEM_S(28,6) = (-0.25*m0rSelf[1]*uSelf[7]*mnuSelf)-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1SrSelf[6]*mnuSelf-0.25*m0rSelf[3]*uSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(28,7) = (-0.25*m0rSelf[0]*uSelf[7]*mnuSelf)+0.5*m1SrSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[5]*mnuSelf-0.25*m0rSelf[3]*uSelf[4]*mnuSelf; 
  data->AEM_S(29,4) = (-0.25*m0rSelf[2]*uSelf[7]*mnuSelf)-0.25*m0rSelf[3]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[4]*mnuSelf; 
  data->AEM_S(29,5) = (-0.45*m0rSelf[3]*uSelf[7]*mnuSelf)-0.25*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(29,6) = (-0.25*m0rSelf[0]*uSelf[7]*mnuSelf)+0.5*m1SrSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[5]*mnuSelf-0.25*m0rSelf[3]*uSelf[4]*mnuSelf; 
  data->AEM_S(29,7) = (-0.45*m0rSelf[1]*uSelf[7]*mnuSelf)-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1SrSelf[6]*mnuSelf-0.45*m0rSelf[3]*uSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(30,4) = (-0.25*m0rSelf[1]*uSelf[7]*mnuSelf)-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1SrSelf[6]*mnuSelf-0.25*m0rSelf[3]*uSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(30,5) = (-0.25*m0rSelf[0]*uSelf[7]*mnuSelf)+0.5*m1SrSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[5]*mnuSelf-0.25*m0rSelf[3]*uSelf[4]*mnuSelf; 
  data->AEM_S(30,6) = (-0.45*m0rSelf[3]*uSelf[7]*mnuSelf)-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(30,7) = (-0.45*m0rSelf[2]*uSelf[7]*mnuSelf)-0.45*m0rSelf[3]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[4]*mnuSelf; 
  data->AEM_S(31,4) = (-0.25*m0rSelf[0]*uSelf[7]*mnuSelf)+0.5*m1SrSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[5]*mnuSelf-0.25*m0rSelf[3]*uSelf[4]*mnuSelf; 
  data->AEM_S(31,5) = (-0.45*m0rSelf[1]*uSelf[7]*mnuSelf)-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1SrSelf[6]*mnuSelf-0.45*m0rSelf[3]*uSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(31,6) = (-0.45*m0rSelf[2]*uSelf[7]*mnuSelf)-0.45*m0rSelf[3]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[4]*mnuSelf; 
  data->AEM_S(31,7) = (-0.81*m0rSelf[3]*uSelf[7]*mnuSelf)-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1SrSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(28,20) = 0.25*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(28,21) = 0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther; 
  data->AEM_S(28,22) = 0.25*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1SrOther[6]*mnuOther+0.25*m0rOther[3]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(28,23) = 0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1SrOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther+0.25*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[3]*uOther[4]*mnuOther; 
  data->AEM_S(29,20) = 0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther; 
  data->AEM_S(29,21) = 0.45*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther+0.45*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(29,22) = 0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1SrOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther+0.25*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[3]*uOther[4]*mnuOther; 
  data->AEM_S(29,23) = 0.45*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1SrOther[6]*mnuOther+0.45*m0rOther[3]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(30,20) = 0.25*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1SrOther[6]*mnuOther+0.25*m0rOther[3]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(30,21) = 0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1SrOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther+0.25*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[3]*uOther[4]*mnuOther; 
  data->AEM_S(30,22) = 0.45*m0rOther[3]*uOther[7]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(30,23) = 0.45*m0rOther[2]*uOther[7]*mnuOther+0.45*m0rOther[3]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther; 
  data->AEM_S(31,20) = 0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1SrOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther+0.25*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[3]*uOther[4]*mnuOther; 
  data->AEM_S(31,21) = 0.45*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1SrOther[6]*mnuOther+0.45*m0rOther[3]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(31,22) = 0.45*m0rOther[2]*uOther[7]*mnuOther+0.45*m0rOther[3]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther; 
  data->AEM_S(31,23) = 0.81*m0rOther[3]*uOther[7]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfY-m0Self*m1OtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[4] = 0.5*m0rOther[3]*m1rSelf[7]-0.5*m0rSelf[3]*m1rOther[7]+0.5*m0rOther[2]*m1rSelf[6]-0.5*m0rSelf[2]*m1rOther[6]+0.5*m0rOther[1]*m1rSelf[5]-0.5*m0rSelf[1]*m1rOther[5]+0.5*m0rOther[0]*m1rSelf[4]-0.5*m0rSelf[0]*m1rOther[4]; 
  m1EffD[5] = 0.5*m0rOther[2]*m1rSelf[7]-0.5*m0rSelf[2]*m1rOther[7]+0.5*m0rOther[3]*m1rSelf[6]-0.5*m0rSelf[3]*m1rOther[6]+0.5*m0rOther[0]*m1rSelf[5]-0.5*m0rSelf[0]*m1rOther[5]+0.5*m0rOther[1]*m1rSelf[4]-0.5*m0rSelf[1]*m1rOther[4]; 
  m1EffD[6] = 0.5*m0rOther[1]*m1rSelf[7]-0.5*m0rSelf[1]*m1rOther[7]+0.5*m0rOther[0]*m1rSelf[6]-0.5*m0rSelf[0]*m1rOther[6]+0.5*m0rOther[3]*m1rSelf[5]-0.5*m0rSelf[3]*m1rOther[5]+0.5*m0rOther[2]*m1rSelf[4]-0.5*m0rSelf[2]*m1rOther[4]; 
  m1EffD[7] = 0.5*m0rOther[0]*m1rSelf[7]-0.5*m0rSelf[0]*m1rOther[7]+0.5*m0rOther[1]*m1rSelf[6]-0.5*m0rSelf[1]*m1rOther[6]+0.5*m0rOther[2]*m1rSelf[5]-0.5*m0rSelf[2]*m1rOther[5]+0.5*m0rOther[3]*m1rSelf[4]-0.5*m0rSelf[3]*m1rOther[4]; 
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
  data->AEM_S(24,8) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(24,9) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(24,10) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,11) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,8) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,9) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(25,10) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,11) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,9) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(26,10) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(26,11) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,8) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,9) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,10) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,11) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(24,12) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(24,13) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(24,14) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(24,15) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(25,12) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(25,13) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(25,14) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(25,15) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(26,12) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(26,13) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(26,14) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(26,15) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(27,12) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(27,13) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(27,14) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(27,15) = -0.5*cMSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(24,24) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(24,25) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(24,26) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(24,27) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(25,24) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(25,25) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(25,26) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(25,27) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(26,24) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(26,25) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(26,26) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(26,27) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(27,24) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(27,25) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(27,26) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(27,27) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(24,28) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(24,29) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(24,30) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(24,31) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(25,28) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(25,29) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(25,30) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(25,31) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(26,28) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(26,29) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(26,30) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(26,31) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(27,28) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(27,29) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(27,30) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(27,31) = 0.5*cMOther[8]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfZ-uSelfZ*m0Self) and uCrossSelfZ ... // 
  data->AEM_S(28,8) = (-0.25*m0rSelf[3]*uSelf[11]*mnuSelf)-0.25*m0rSelf[2]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(28,9) = (-0.25*m0rSelf[2]*uSelf[11]*mnuSelf)-0.25*m0rSelf[3]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1SrSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(28,10) = (-0.25*m0rSelf[1]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1SrSelf[10]*mnuSelf-0.25*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(28,11) = (-0.25*m0rSelf[0]*uSelf[11]*mnuSelf)+0.5*m1SrSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(29,8) = (-0.25*m0rSelf[2]*uSelf[11]*mnuSelf)-0.25*m0rSelf[3]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1SrSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(29,9) = (-0.45*m0rSelf[3]*uSelf[11]*mnuSelf)-0.25*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(29,10) = (-0.25*m0rSelf[0]*uSelf[11]*mnuSelf)+0.5*m1SrSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(29,11) = (-0.45*m0rSelf[1]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1SrSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(30,8) = (-0.25*m0rSelf[1]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1SrSelf[10]*mnuSelf-0.25*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(30,9) = (-0.25*m0rSelf[0]*uSelf[11]*mnuSelf)+0.5*m1SrSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(30,10) = (-0.45*m0rSelf[3]*uSelf[11]*mnuSelf)-0.45*m0rSelf[2]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(30,11) = (-0.45*m0rSelf[2]*uSelf[11]*mnuSelf)-0.45*m0rSelf[3]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1SrSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(31,8) = (-0.25*m0rSelf[0]*uSelf[11]*mnuSelf)+0.5*m1SrSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(31,9) = (-0.45*m0rSelf[1]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1SrSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(31,10) = (-0.45*m0rSelf[2]*uSelf[11]*mnuSelf)-0.45*m0rSelf[3]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1SrSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(31,11) = (-0.81*m0rSelf[3]*uSelf[11]*mnuSelf)-0.45*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1SrSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherZ-uOtherZ*m0Other) and uCrossOtherZ ... // 
  data->AEM_S(28,24) = 0.25*m0rOther[3]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(28,25) = 0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1SrOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(28,26) = 0.25*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1SrOther[10]*mnuOther+0.25*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(28,27) = 0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1SrOther[11]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(29,24) = 0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1SrOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(29,25) = 0.45*m0rOther[3]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(29,26) = 0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1SrOther[11]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(29,27) = 0.45*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1SrOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(30,24) = 0.25*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1SrOther[10]*mnuOther+0.25*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(30,25) = 0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1SrOther[11]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(30,26) = 0.45*m0rOther[3]*uOther[11]*mnuOther+0.45*m0rOther[2]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(30,27) = 0.45*m0rOther[2]*uOther[11]*mnuOther+0.45*m0rOther[3]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1SrOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(31,24) = 0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1SrOther[11]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(31,25) = 0.45*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1SrOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(31,26) = 0.45*m0rOther[2]*uOther[11]*mnuOther+0.45*m0rOther[3]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1SrOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(31,27) = 0.81*m0rOther[3]*uOther[11]*mnuOther+0.45*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1SrOther[8]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfZ-m0Self*m1OtherZ)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[8] = 0.5*m0rOther[3]*m1rSelf[11]-0.5*m0rSelf[3]*m1rOther[11]+0.5*m0rOther[2]*m1rSelf[10]-0.5*m0rSelf[2]*m1rOther[10]+0.5*m0rOther[1]*m1rSelf[9]-0.5*m0rSelf[1]*m1rOther[9]+0.5*m0rOther[0]*m1rSelf[8]-0.5*m0rSelf[0]*m1rOther[8]; 
  m1EffD[9] = 0.5*m0rOther[2]*m1rSelf[11]-0.5*m0rSelf[2]*m1rOther[11]+0.5*m0rOther[3]*m1rSelf[10]-0.5*m0rSelf[3]*m1rOther[10]+0.5*m0rOther[0]*m1rSelf[9]-0.5*m0rSelf[0]*m1rOther[9]+0.5*m0rOther[1]*m1rSelf[8]-0.5*m0rSelf[1]*m1rOther[8]; 
  m1EffD[10] = 0.5*m0rOther[1]*m1rSelf[11]-0.5*m0rSelf[1]*m1rOther[11]+0.5*m0rOther[0]*m1rSelf[10]-0.5*m0rSelf[0]*m1rOther[10]+0.5*m0rOther[3]*m1rSelf[9]-0.5*m0rSelf[3]*m1rOther[9]+0.5*m0rOther[2]*m1rSelf[8]-0.5*m0rSelf[2]*m1rOther[8]; 
  m1EffD[11] = 0.5*m0rOther[0]*m1rSelf[11]-0.5*m0rSelf[0]*m1rOther[11]+0.5*m0rOther[1]*m1rSelf[10]-0.5*m0rSelf[1]*m1rOther[10]+0.5*m0rOther[2]*m1rSelf[9]-0.5*m0rSelf[2]*m1rOther[9]+0.5*m0rOther[3]*m1rSelf[8]-0.5*m0rSelf[3]*m1rOther[8]; 
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
    ucMSelf[0] += 0.5*cMSelf[a0+3]*uSelf[a0+3]+0.5*cMSelf[a0+2]*uSelf[a0+2]+0.5*cMSelf[a0+1]*uSelf[a0+1]+0.5*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.5*cMSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.5*cMSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*cMSelf[a0+2]; 
    ucMSelf[3] += 0.5*cMSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*cMSelf[a0+3]+0.5*cMSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*cMSelf[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(28,12) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(28,13) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(28,14) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(28,15) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(29,12) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(29,13) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(29,14) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(29,15) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(30,12) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(30,13) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(30,14) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(30,15) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(31,12) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(31,13) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(31,14) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(31,15) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
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
    ucMOther[0] += 0.5*cMOther[a0+3]*uOther[a0+3]+0.5*cMOther[a0+2]*uOther[a0+2]+0.5*cMOther[a0+1]*uOther[a0+1]+0.5*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.5*cMOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+1]+0.5*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.5*cMOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+2]+0.5*uOther[a0]*cMOther[a0+2]; 
    ucMOther[3] += 0.5*cMOther[a0]*uOther[a0+3]+0.5*uOther[a0]*cMOther[a0+3]+0.5*cMOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*cMOther[a0+2]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(28,28) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(28,29) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(28,30) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(28,31) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(29,28) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(29,29) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(29,30) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(29,31) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(30,28) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(30,29) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(30,30) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(30,31) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(31,28) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(31,29) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(31,30) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(31,31) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
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
    kinESelf[0] += 0.5*m1rSelf[a0+3]*uSelf[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.5*m1rSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.5*m1rSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]; 
    kinESelf[3] += 0.5*m1rSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]; 
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
    kinEOther[0] += 0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    kinEOther[3] += 0.5*m1rOther[a0]*uOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
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
    relKinE[0] += 0.5*m1EffD[a0+3]*uSelf[a0+3]-0.5*m1EffD[a0+3]*uOther[a0+3]+0.5*m1EffD[a0+2]*uSelf[a0+2]-0.5*m1EffD[a0+2]*uOther[a0+2]+0.5*m1EffD[a0+1]*uSelf[a0+1]-0.5*m1EffD[a0+1]*uOther[a0+1]+0.5*m1EffD[a0]*uSelf[a0]-0.5*m1EffD[a0]*uOther[a0]; 
    relKinE[1] += 0.5*m1EffD[a0+2]*uSelf[a0+3]-0.5*m1EffD[a0+2]*uOther[a0+3]+0.5*uSelf[a0+2]*m1EffD[a0+3]-0.5*uOther[a0+2]*m1EffD[a0+3]+0.5*m1EffD[a0]*uSelf[a0+1]-0.5*m1EffD[a0]*uOther[a0+1]+0.5*uSelf[a0]*m1EffD[a0+1]-0.5*uOther[a0]*m1EffD[a0+1]; 
    relKinE[2] += 0.5*m1EffD[a0+1]*uSelf[a0+3]-0.5*m1EffD[a0+1]*uOther[a0+3]+0.5*uSelf[a0+1]*m1EffD[a0+3]-0.5*uOther[a0+1]*m1EffD[a0+3]+0.5*m1EffD[a0]*uSelf[a0+2]-0.5*m1EffD[a0]*uOther[a0+2]+0.5*uSelf[a0]*m1EffD[a0+2]-0.5*uOther[a0]*m1EffD[a0+2]; 
    relKinE[3] += 0.5*m1EffD[a0]*uSelf[a0+3]-0.5*m1EffD[a0]*uOther[a0+3]+0.5*uSelf[a0]*m1EffD[a0+3]-0.5*uOther[a0]*m1EffD[a0+3]+0.5*m1EffD[a0+1]*uSelf[a0+2]-0.5*m1EffD[a0+1]*uOther[a0+2]+0.5*uSelf[a0+1]*m1EffD[a0+2]-0.5*uOther[a0+1]*m1EffD[a0+2]; 
  } 
 
  // Divide m0Other*(m2Self-kinESelf) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Other and m2Self-uSelf.m1Self. 
  double m0OtherThESelf[4]; 
  m0OtherThESelf[0] = 0.5*m0rOther[3]*m2rSelf[3]-0.5*kinESelf[3]*m0rOther[3]+0.5*m0rOther[2]*m2rSelf[2]-0.5*kinESelf[2]*m0rOther[2]+0.5*m0rOther[1]*m2rSelf[1]-0.5*kinESelf[1]*m0rOther[1]+0.5*m0rOther[0]*m2rSelf[0]-0.5*kinESelf[0]*m0rOther[0]; 
  m0OtherThESelf[1] = 0.5*m0rOther[2]*m2rSelf[3]+0.5*m2rSelf[2]*m0rOther[3]-0.5*kinESelf[2]*m0rOther[3]-0.5*m0rOther[2]*kinESelf[3]+0.5*m0rOther[0]*m2rSelf[1]+0.5*m2rSelf[0]*m0rOther[1]-0.5*kinESelf[0]*m0rOther[1]-0.5*m0rOther[0]*kinESelf[1]; 
  m0OtherThESelf[2] = 0.5*m0rOther[1]*m2rSelf[3]+0.5*m2rSelf[1]*m0rOther[3]-0.5*kinESelf[1]*m0rOther[3]-0.5*m0rOther[1]*kinESelf[3]+0.5*m0rOther[0]*m2rSelf[2]+0.5*m2rSelf[0]*m0rOther[2]-0.5*kinESelf[0]*m0rOther[2]-0.5*m0rOther[0]*kinESelf[2]; 
  m0OtherThESelf[3] = 0.5*m0rOther[0]*m2rSelf[3]+0.5*m2rSelf[0]*m0rOther[3]-0.5*kinESelf[0]*m0rOther[3]-0.5*m0rOther[0]*kinESelf[3]+0.5*m0rOther[1]*m2rSelf[2]+0.5*m2rSelf[1]*m0rOther[2]-0.5*kinESelf[1]*m0rOther[2]-0.5*m0rOther[1]*kinESelf[2]; 
  dataDiv->BEV_S << m0OtherThESelf[0],m0OtherThESelf[1],m0OtherThESelf[2],m0OtherThESelf[3]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthSelf[4]; 
  Eigen::Map<VectorXd>(effEthSelf,4,1) = dataDiv->u_S; 
 
  // Divide m0Self*(m2Other-kinEOther) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Self and m2Other-uOther.m1Other. 
  double m0SelfThEOther[4]; 
  m0SelfThEOther[0] = 0.5*m0rSelf[3]*m2rOther[3]-0.5*kinEOther[3]*m0rSelf[3]+0.5*m0rSelf[2]*m2rOther[2]-0.5*kinEOther[2]*m0rSelf[2]+0.5*m0rSelf[1]*m2rOther[1]-0.5*kinEOther[1]*m0rSelf[1]+0.5*m0rSelf[0]*m2rOther[0]-0.5*kinEOther[0]*m0rSelf[0]; 
  m0SelfThEOther[1] = 0.5*m0rSelf[2]*m2rOther[3]+0.5*m2rOther[2]*m0rSelf[3]-0.5*kinEOther[2]*m0rSelf[3]-0.5*m0rSelf[2]*kinEOther[3]+0.5*m0rSelf[0]*m2rOther[1]+0.5*m2rOther[0]*m0rSelf[1]-0.5*kinEOther[0]*m0rSelf[1]-0.5*m0rSelf[0]*kinEOther[1]; 
  m0SelfThEOther[2] = 0.5*m0rSelf[1]*m2rOther[3]+0.5*m2rOther[1]*m0rSelf[3]-0.5*kinEOther[1]*m0rSelf[3]-0.5*m0rSelf[1]*kinEOther[3]+0.5*m0rSelf[0]*m2rOther[2]+0.5*m2rOther[0]*m0rSelf[2]-0.5*kinEOther[0]*m0rSelf[2]-0.5*m0rSelf[0]*kinEOther[2]; 
  m0SelfThEOther[3] = 0.5*m0rSelf[0]*m2rOther[3]+0.5*m2rOther[0]*m0rSelf[3]-0.5*kinEOther[0]*m0rSelf[3]-0.5*m0rSelf[0]*kinEOther[3]+0.5*m0rSelf[1]*m2rOther[2]+0.5*m2rOther[1]*m0rSelf[2]-0.5*kinEOther[1]*m0rSelf[2]-0.5*m0rSelf[1]*kinEOther[2]; 
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
 
void VmLBOCrossPrimMoments2x3vSer_P2(binOpData_t *data, binOpData_t *dataDiv,const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:       free parameter beta+1. This has to be >0. 
  // nu, m:              collisionality and mass. 
  // m0,m1,m2:           moments of the distribution function. 
  // u,vtSq:             self primitive moments: mean flow velocity and thermal speed squared. 
  // cM,cE:              corrections to momentum and energy conservation due to finite velocity space. 
  // uCross,vtSqCross:   cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-1.936491673103709*m0Self[7])-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.936491673103709*m0Self[7])-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.936491673103709*m0Self[7]-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.936491673103709*m0Self[7]-1.936491673103709*m0Self[6]+1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[8]; 
  double m1rSelf[24]; 
  double m2rSelf[8]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m0rSelf[3] = 0.0; 
    m0rSelf[4] = 0.0; 
    m0rSelf[5] = 0.0; 
    m0rSelf[6] = 0.0; 
    m0rSelf[7] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    m1rSelf[2] = 0.0; 
    m1rSelf[3] = 0.0; 
    m1rSelf[4] = 0.0; 
    m1rSelf[5] = 0.0; 
    m1rSelf[6] = 0.0; 
    m1rSelf[7] = 0.0; 
    m1rSelf[8] = m1Self[8]; 
    m1rSelf[9] = 0.0; 
    m1rSelf[10] = 0.0; 
    m1rSelf[11] = 0.0; 
    m1rSelf[12] = 0.0; 
    m1rSelf[13] = 0.0; 
    m1rSelf[14] = 0.0; 
    m1rSelf[15] = 0.0; 
    m1rSelf[16] = m1Self[16]; 
    m1rSelf[17] = 0.0; 
    m1rSelf[18] = 0.0; 
    m1rSelf[19] = 0.0; 
    m1rSelf[20] = 0.0; 
    m1rSelf[21] = 0.0; 
    m1rSelf[22] = 0.0; 
    m1rSelf[23] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    m2rSelf[3] = 0.0; 
    m2rSelf[4] = 0.0; 
    m2rSelf[5] = 0.0; 
    m2rSelf[6] = 0.0; 
    m2rSelf[7] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
    m0rSelf[3] = m0Self[3]; 
    m0rSelf[4] = m0Self[4]; 
    m0rSelf[5] = m0Self[5]; 
    m0rSelf[6] = m0Self[6]; 
    m0rSelf[7] = m0Self[7]; 
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
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = m2Self[1]; 
    m2rSelf[2] = m2Self[2]; 
    m2rSelf[3] = m2Self[3]; 
    m2rSelf[4] = m2Self[4]; 
    m2rSelf[5] = m2Self[5]; 
    m2rSelf[6] = m2Self[6]; 
    m2rSelf[7] = m2Self[7]; 
  } 
 
  if ((-1.936491673103709*m0Other[7])-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-1.936491673103709*m0Other[7])-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.936491673103709*m0Other[7]-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.936491673103709*m0Other[7]-1.936491673103709*m0Other[6]+1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[8]; 
  double m1rOther[24]; 
  double m2rOther[8]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m0rOther[3] = 0.0; 
    m0rOther[4] = 0.0; 
    m0rOther[5] = 0.0; 
    m0rOther[6] = 0.0; 
    m0rOther[7] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    m1rOther[2] = 0.0; 
    m1rOther[3] = 0.0; 
    m1rOther[4] = 0.0; 
    m1rOther[5] = 0.0; 
    m1rOther[6] = 0.0; 
    m1rOther[7] = 0.0; 
    m1rOther[8] = m1Other[8]; 
    m1rOther[9] = 0.0; 
    m1rOther[10] = 0.0; 
    m1rOther[11] = 0.0; 
    m1rOther[12] = 0.0; 
    m1rOther[13] = 0.0; 
    m1rOther[14] = 0.0; 
    m1rOther[15] = 0.0; 
    m1rOther[16] = m1Other[16]; 
    m1rOther[17] = 0.0; 
    m1rOther[18] = 0.0; 
    m1rOther[19] = 0.0; 
    m1rOther[20] = 0.0; 
    m1rOther[21] = 0.0; 
    m1rOther[22] = 0.0; 
    m1rOther[23] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    m2rOther[3] = 0.0; 
    m2rOther[4] = 0.0; 
    m2rOther[5] = 0.0; 
    m2rOther[6] = 0.0; 
    m2rOther[7] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
    m0rOther[3] = m0Other[3]; 
    m0rOther[4] = m0Other[4]; 
    m0rOther[5] = m0Other[5]; 
    m0rOther[6] = m0Other[6]; 
    m0rOther[7] = m0Other[7]; 
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
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
    m2rOther[3] = m2Other[3]; 
    m2rOther[4] = m2Other[4]; 
    m2rOther[5] = m2Other[5]; 
    m2rOther[6] = m2Other[6]; 
    m2rOther[7] = m2Other[7]; 
  } 
 
  // Declare Eigen matrix and vectors for weak system. 
  data->AEM_S = Eigen::MatrixXd::Zero(64,64); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double mnuM1sum[24]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<24; vd++) 
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
  data->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(1,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,3) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(1,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(1,5) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(1,6) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(1,7) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(2,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(2,3) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(2,4) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(2,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(2,6) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(2,7) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,1) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(3,2) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(3,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(3,6) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,7) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(4,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,2) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(4,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(4,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,6) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,7) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(5,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(5,1) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(5,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(5,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(5,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(5,6) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(5,7) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,0) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(6,1) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(6,2) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(6,3) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,4) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(6,5) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(6,6) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,0) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(7,1) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(7,2) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,3) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,4) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(7,5) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,6) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(7,7) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(0,24) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,25) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,26) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,27) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,28) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,29) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(0,30) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(0,31) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(1,24) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,25) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,26) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,27) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,28) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,29) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(1,30) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,31) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(2,24) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,25) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,26) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,27) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,28) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(2,29) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,30) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(2,31) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,24) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,25) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,26) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,27) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,28) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,29) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,30) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,31) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(4,24) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,25) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,26) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(4,27) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,28) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(4,30) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(4,31) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,24) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,25) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,26) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(5,27) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,29) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,30) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(5,31) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,24) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,25) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(6,26) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(6,27) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,28) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(6,29) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,30) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(6,31) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,24) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,25) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(7,26) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,27) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(7,28) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,29) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(7,30) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,31) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,32) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,33) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,34) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,35) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(0,36) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(0,37) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(0,38) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(0,39) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(1,32) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,33) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,34) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,35) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,36) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(1,37) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(1,38) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(1,39) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(2,32) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,33) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,34) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,35) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(2,36) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(2,37) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(2,38) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(2,39) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(3,32) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,33) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,34) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,35) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,36) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,37) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,38) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(3,39) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(4,32) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(4,33) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(4,34) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(4,35) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,36) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(4,38) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(4,39) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(5,32) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(5,33) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(5,34) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(5,35) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(5,37) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,38) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(5,39) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(6,32) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(6,33) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(6,34) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(6,35) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(6,36) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(6,37) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(6,38) = 0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,39) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(7,32) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(7,33) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(7,34) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(7,35) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(7,36) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(7,37) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(7,38) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(7,39) = 0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,56) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,57) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,58) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,59) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(0,60) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(0,61) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(0,62) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(0,63) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(1,56) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,57) = (-0.4472135954999579*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,58) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,59) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,60) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(1,61) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(1,62) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(1,63) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(2,56) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,57) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,58) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,59) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(2,60) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(2,61) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(2,62) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(2,63) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(3,56) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,57) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,58) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,59) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(3,60) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,61) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,62) = (-0.4*cMOther[7]*mnuOther)-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(3,63) = (-0.4*cMOther[6]*mnuOther)-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(4,56) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,57) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(4,58) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(4,59) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,60) = (-0.31943828249997*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(4,62) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(4,63) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(5,56) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,57) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(5,58) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(5,59) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(5,61) = (-0.31943828249997*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(5,62) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(5,63) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(6,56) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,57) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(6,58) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(6,59) = (-0.4*cMOther[7]*mnuOther)-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(6,60) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(6,61) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(6,62) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.31943828249997*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(6,63) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(7,56) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,57) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(7,58) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(7,59) = (-0.4*cMOther[6]*mnuOther)-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(7,60) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(7,61) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(7,62) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(7,63) = (-0.31943828249997*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(24,0) = 0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(24,1) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(24,2) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(24,3) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(24,4) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(24,5) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(24,6) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(24,7) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(25,0) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(25,1) = 0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(25,2) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(25,3) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(25,4) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(25,5) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(25,6) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(25,7) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(26,0) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(26,1) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(26,2) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(26,3) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(26,4) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(26,5) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(26,6) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(26,7) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(27,0) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(27,1) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(27,2) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(27,3) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(27,4) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(27,5) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(27,6) = 0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(27,7) = 0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(28,0) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(28,1) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(28,2) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(28,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(28,4) = 0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(28,6) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(28,7) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(29,0) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(29,1) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(29,2) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(29,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(29,5) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(29,6) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(29,7) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(30,0) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(30,1) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(30,2) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(30,3) = 0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(30,4) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(30,5) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(30,6) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(30,7) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(31,0) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(31,1) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(31,2) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(31,3) = 0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(31,4) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(31,5) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(31,6) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(31,7) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(24,32) = 0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(24,33) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(24,34) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(24,35) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(24,36) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(24,37) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(24,38) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(24,39) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(25,32) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(25,33) = 0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(25,34) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(25,35) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(25,36) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(25,37) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(25,38) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(25,39) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(26,32) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(26,33) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(26,34) = 0.4472135954999579*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(26,35) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(26,36) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(26,37) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(26,38) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(26,39) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(27,32) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(27,33) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(27,34) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(27,35) = 0.4472135954999579*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(27,36) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(27,37) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(27,38) = 0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(27,39) = 0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(28,32) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(28,33) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(28,34) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(28,35) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(28,36) = 0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(28,38) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(28,39) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(29,32) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(29,33) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(29,34) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(29,35) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(29,37) = 0.31943828249997*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(29,38) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(29,39) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(30,32) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(30,33) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(30,34) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(30,35) = 0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(30,36) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(30,37) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(30,38) = 0.4472135954999579*m1rOther[5]*mnuOther+0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(30,39) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(31,32) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(31,33) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(31,34) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(31,35) = 0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(31,36) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(31,37) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(31,38) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(31,39) = 0.31943828249997*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(8,8) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,9) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,10) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,11) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(8,12) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(8,13) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(8,14) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(8,15) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(9,8) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,9) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(9,10) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,11) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(9,12) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,13) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(9,14) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,15) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(10,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,9) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(10,10) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(10,11) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(10,12) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(10,13) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,14) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(10,15) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,8) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,9) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,10) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,11) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(11,12) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,13) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,14) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,15) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,8) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(12,9) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,10) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(12,11) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(12,12) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,14) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,15) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(13,8) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(13,9) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(13,10) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,11) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(13,13) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(13,15) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,8) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(14,9) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(14,10) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(14,11) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,12) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,13) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(14,14) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(14,15) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,8) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(15,9) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(15,10) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,11) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,12) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(15,13) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,14) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,15) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(8,24) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,25) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(8,26) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(8,27) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(8,28) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(8,29) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(8,30) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(8,31) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(9,24) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,25) = (-0.4472135954999579*cMSelf[12]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(9,26) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(9,27) = (-0.447213595499958*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(9,28) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,29) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(9,30) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(9,31) = -0.5000000000000001*cMSelf[13]*mnuSelf; 
  data->AEM_S(10,24) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,25) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(10,26) = (-0.4472135954999579*cMSelf[13]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(10,27) = (-0.447213595499958*cMSelf[15]*mnuSelf)-0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(10,28) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(10,29) = -0.4472135954999579*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,30) = -0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(10,31) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,24) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,25) = (-0.447213595499958*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(11,26) = (-0.447213595499958*cMSelf[15]*mnuSelf)-0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(11,27) = (-0.4472135954999579*cMSelf[13]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(11,28) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,29) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,30) = (-0.4*cMSelf[15]*mnuSelf)-0.447213595499958*cMSelf[9]*mnuSelf; 
  data->AEM_S(11,31) = (-0.4*cMSelf[14]*mnuSelf)-0.447213595499958*cMSelf[10]*mnuSelf; 
  data->AEM_S(12,24) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(12,25) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(12,26) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(12,27) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(12,28) = (-0.31943828249997*cMSelf[12]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(12,30) = (-0.31943828249997*cMSelf[14]*mnuSelf)-0.5000000000000001*cMSelf[10]*mnuSelf; 
  data->AEM_S(12,31) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(13,24) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(13,25) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(13,26) = -0.4472135954999579*cMSelf[10]*mnuSelf; 
  data->AEM_S(13,27) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(13,29) = (-0.31943828249997*cMSelf[13]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(13,30) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(13,31) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.5000000000000001*cMSelf[9]*mnuSelf; 
  data->AEM_S(14,24) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(14,25) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(14,26) = -0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(14,27) = (-0.4*cMSelf[15]*mnuSelf)-0.447213595499958*cMSelf[9]*mnuSelf; 
  data->AEM_S(14,28) = (-0.31943828249997*cMSelf[14]*mnuSelf)-0.5000000000000001*cMSelf[10]*mnuSelf; 
  data->AEM_S(14,29) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(14,30) = (-0.4472135954999579*cMSelf[13]*mnuSelf)-0.31943828249997*cMSelf[12]*mnuSelf-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(14,31) = -0.4*cMSelf[11]*mnuSelf; 
  data->AEM_S(15,24) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(15,25) = -0.5000000000000001*cMSelf[13]*mnuSelf; 
  data->AEM_S(15,26) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(15,27) = (-0.4*cMSelf[14]*mnuSelf)-0.447213595499958*cMSelf[10]*mnuSelf; 
  data->AEM_S(15,28) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(15,29) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.5000000000000001*cMSelf[9]*mnuSelf; 
  data->AEM_S(15,30) = -0.4*cMSelf[11]*mnuSelf; 
  data->AEM_S(15,31) = (-0.31943828249997*cMSelf[13]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf-0.5*cMSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(8,40) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,41) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(8,42) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,43) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(8,44) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(8,45) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(8,46) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(8,47) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(9,40) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(9,41) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,42) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,43) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(9,44) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(9,45) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(9,46) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(9,47) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(10,40) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,41) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(10,42) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(10,43) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(10,44) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(10,45) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(10,46) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(10,47) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(11,40) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(11,41) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(11,42) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,43) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(11,44) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(11,45) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(11,46) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(11,47) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(12,40) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(12,41) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(12,42) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(12,43) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(12,44) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(12,46) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(12,47) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(13,40) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(13,41) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(13,42) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(13,43) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(13,45) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(13,46) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(13,47) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(14,40) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(14,41) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(14,42) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(14,43) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(14,44) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(14,45) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(14,46) = 0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(14,47) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(15,40) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(15,41) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(15,42) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(15,43) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(15,44) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(15,45) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(15,46) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(15,47) = 0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(8,56) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(8,57) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(8,58) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(8,59) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(8,60) = -0.5*cMOther[12]*mnuOther; 
  data->AEM_S(8,61) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(8,62) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(8,63) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(9,56) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(9,57) = (-0.4472135954999579*cMOther[12]*mnuOther)-0.5*cMOther[8]*mnuOther; 
  data->AEM_S(9,58) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(9,59) = (-0.447213595499958*cMOther[14]*mnuOther)-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(9,60) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(9,61) = -0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(9,62) = -0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(9,63) = -0.5000000000000001*cMOther[13]*mnuOther; 
  data->AEM_S(10,56) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(10,57) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(10,58) = (-0.4472135954999579*cMOther[13]*mnuOther)-0.5*cMOther[8]*mnuOther; 
  data->AEM_S(10,59) = (-0.447213595499958*cMOther[15]*mnuOther)-0.5*cMOther[9]*mnuOther; 
  data->AEM_S(10,60) = -0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(10,61) = -0.4472135954999579*cMOther[10]*mnuOther; 
  data->AEM_S(10,62) = -0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(10,63) = -0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(11,56) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(11,57) = (-0.447213595499958*cMOther[14]*mnuOther)-0.5*cMOther[10]*mnuOther; 
  data->AEM_S(11,58) = (-0.447213595499958*cMOther[15]*mnuOther)-0.5*cMOther[9]*mnuOther; 
  data->AEM_S(11,59) = (-0.4472135954999579*cMOther[13]*mnuOther)-0.4472135954999579*cMOther[12]*mnuOther-0.5*cMOther[8]*mnuOther; 
  data->AEM_S(11,60) = -0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(11,61) = -0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(11,62) = (-0.4*cMOther[15]*mnuOther)-0.447213595499958*cMOther[9]*mnuOther; 
  data->AEM_S(11,63) = (-0.4*cMOther[14]*mnuOther)-0.447213595499958*cMOther[10]*mnuOther; 
  data->AEM_S(12,56) = -0.5*cMOther[12]*mnuOther; 
  data->AEM_S(12,57) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(12,58) = -0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(12,59) = -0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(12,60) = (-0.31943828249997*cMOther[12]*mnuOther)-0.5*cMOther[8]*mnuOther; 
  data->AEM_S(12,62) = (-0.31943828249997*cMOther[14]*mnuOther)-0.5000000000000001*cMOther[10]*mnuOther; 
  data->AEM_S(12,63) = -0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(13,56) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(13,57) = -0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(13,58) = -0.4472135954999579*cMOther[10]*mnuOther; 
  data->AEM_S(13,59) = -0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(13,61) = (-0.31943828249997*cMOther[13]*mnuOther)-0.5*cMOther[8]*mnuOther; 
  data->AEM_S(13,62) = -0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(13,63) = (-0.31943828249997*cMOther[15]*mnuOther)-0.5000000000000001*cMOther[9]*mnuOther; 
  data->AEM_S(14,56) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(14,57) = -0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(14,58) = -0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(14,59) = (-0.4*cMOther[15]*mnuOther)-0.447213595499958*cMOther[9]*mnuOther; 
  data->AEM_S(14,60) = (-0.31943828249997*cMOther[14]*mnuOther)-0.5000000000000001*cMOther[10]*mnuOther; 
  data->AEM_S(14,61) = -0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(14,62) = (-0.4472135954999579*cMOther[13]*mnuOther)-0.31943828249997*cMOther[12]*mnuOther-0.5*cMOther[8]*mnuOther; 
  data->AEM_S(14,63) = -0.4*cMOther[11]*mnuOther; 
  data->AEM_S(15,56) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(15,57) = -0.5000000000000001*cMOther[13]*mnuOther; 
  data->AEM_S(15,58) = -0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(15,59) = (-0.4*cMOther[14]*mnuOther)-0.447213595499958*cMOther[10]*mnuOther; 
  data->AEM_S(15,60) = -0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(15,61) = (-0.31943828249997*cMOther[15]*mnuOther)-0.5000000000000001*cMOther[9]*mnuOther; 
  data->AEM_S(15,62) = -0.4*cMOther[11]*mnuOther; 
  data->AEM_S(15,63) = (-0.31943828249997*cMOther[13]*mnuOther)-0.4472135954999579*cMOther[12]*mnuOther-0.5*cMOther[8]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(24,8) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(24,9) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(24,10) = 0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(24,11) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(24,12) = 0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(24,13) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(24,14) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(24,15) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(25,8) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(25,9) = 0.4472135954999579*m1rSelf[12]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(25,10) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(25,11) = 0.447213595499958*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(25,12) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(25,13) = 0.5000000000000001*m1rSelf[15]*mnuSelf; 
  data->AEM_S(25,14) = 0.447213595499958*m1rSelf[11]*mnuSelf; 
  data->AEM_S(25,15) = 0.5000000000000001*m1rSelf[13]*mnuSelf; 
  data->AEM_S(26,8) = 0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(26,9) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(26,10) = 0.4472135954999579*m1rSelf[13]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(26,11) = 0.447213595499958*m1rSelf[15]*mnuSelf+0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(26,12) = 0.5000000000000001*m1rSelf[14]*mnuSelf; 
  data->AEM_S(26,13) = 0.4472135954999579*m1rSelf[10]*mnuSelf; 
  data->AEM_S(26,14) = 0.5000000000000001*m1rSelf[12]*mnuSelf; 
  data->AEM_S(26,15) = 0.447213595499958*m1rSelf[11]*mnuSelf; 
  data->AEM_S(27,8) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(27,9) = 0.447213595499958*m1rSelf[14]*mnuSelf+0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(27,10) = 0.447213595499958*m1rSelf[15]*mnuSelf+0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(27,11) = 0.4472135954999579*m1rSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(27,12) = 0.4472135954999579*m1rSelf[11]*mnuSelf; 
  data->AEM_S(27,13) = 0.4472135954999579*m1rSelf[11]*mnuSelf; 
  data->AEM_S(27,14) = 0.4*m1rSelf[15]*mnuSelf+0.447213595499958*m1rSelf[9]*mnuSelf; 
  data->AEM_S(27,15) = 0.4*m1rSelf[14]*mnuSelf+0.447213595499958*m1rSelf[10]*mnuSelf; 
  data->AEM_S(28,8) = 0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(28,9) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(28,10) = 0.5000000000000001*m1rSelf[14]*mnuSelf; 
  data->AEM_S(28,11) = 0.4472135954999579*m1rSelf[11]*mnuSelf; 
  data->AEM_S(28,12) = 0.31943828249997*m1rSelf[12]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(28,14) = 0.31943828249997*m1rSelf[14]*mnuSelf+0.5000000000000001*m1rSelf[10]*mnuSelf; 
  data->AEM_S(28,15) = 0.4472135954999579*m1rSelf[15]*mnuSelf; 
  data->AEM_S(29,8) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(29,9) = 0.5000000000000001*m1rSelf[15]*mnuSelf; 
  data->AEM_S(29,10) = 0.4472135954999579*m1rSelf[10]*mnuSelf; 
  data->AEM_S(29,11) = 0.4472135954999579*m1rSelf[11]*mnuSelf; 
  data->AEM_S(29,13) = 0.31943828249997*m1rSelf[13]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(29,14) = 0.4472135954999579*m1rSelf[14]*mnuSelf; 
  data->AEM_S(29,15) = 0.31943828249997*m1rSelf[15]*mnuSelf+0.5000000000000001*m1rSelf[9]*mnuSelf; 
  data->AEM_S(30,8) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(30,9) = 0.447213595499958*m1rSelf[11]*mnuSelf; 
  data->AEM_S(30,10) = 0.5000000000000001*m1rSelf[12]*mnuSelf; 
  data->AEM_S(30,11) = 0.4*m1rSelf[15]*mnuSelf+0.447213595499958*m1rSelf[9]*mnuSelf; 
  data->AEM_S(30,12) = 0.31943828249997*m1rSelf[14]*mnuSelf+0.5000000000000001*m1rSelf[10]*mnuSelf; 
  data->AEM_S(30,13) = 0.4472135954999579*m1rSelf[14]*mnuSelf; 
  data->AEM_S(30,14) = 0.4472135954999579*m1rSelf[13]*mnuSelf+0.31943828249997*m1rSelf[12]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(30,15) = 0.4*m1rSelf[11]*mnuSelf; 
  data->AEM_S(31,8) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(31,9) = 0.5000000000000001*m1rSelf[13]*mnuSelf; 
  data->AEM_S(31,10) = 0.447213595499958*m1rSelf[11]*mnuSelf; 
  data->AEM_S(31,11) = 0.4*m1rSelf[14]*mnuSelf+0.447213595499958*m1rSelf[10]*mnuSelf; 
  data->AEM_S(31,12) = 0.4472135954999579*m1rSelf[15]*mnuSelf; 
  data->AEM_S(31,13) = 0.31943828249997*m1rSelf[15]*mnuSelf+0.5000000000000001*m1rSelf[9]*mnuSelf; 
  data->AEM_S(31,14) = 0.4*m1rSelf[11]*mnuSelf; 
  data->AEM_S(31,15) = 0.31943828249997*m1rSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(24,40) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(24,41) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(24,42) = 0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(24,43) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(24,44) = 0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(24,45) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(24,46) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(24,47) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(25,40) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(25,41) = 0.4472135954999579*m1rOther[12]*mnuOther+0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(25,42) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(25,43) = 0.447213595499958*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(25,44) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(25,45) = 0.5000000000000001*m1rOther[15]*mnuOther; 
  data->AEM_S(25,46) = 0.447213595499958*m1rOther[11]*mnuOther; 
  data->AEM_S(25,47) = 0.5000000000000001*m1rOther[13]*mnuOther; 
  data->AEM_S(26,40) = 0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(26,41) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(26,42) = 0.4472135954999579*m1rOther[13]*mnuOther+0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(26,43) = 0.447213595499958*m1rOther[15]*mnuOther+0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(26,44) = 0.5000000000000001*m1rOther[14]*mnuOther; 
  data->AEM_S(26,45) = 0.4472135954999579*m1rOther[10]*mnuOther; 
  data->AEM_S(26,46) = 0.5000000000000001*m1rOther[12]*mnuOther; 
  data->AEM_S(26,47) = 0.447213595499958*m1rOther[11]*mnuOther; 
  data->AEM_S(27,40) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(27,41) = 0.447213595499958*m1rOther[14]*mnuOther+0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(27,42) = 0.447213595499958*m1rOther[15]*mnuOther+0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(27,43) = 0.4472135954999579*m1rOther[13]*mnuOther+0.4472135954999579*m1rOther[12]*mnuOther+0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(27,44) = 0.4472135954999579*m1rOther[11]*mnuOther; 
  data->AEM_S(27,45) = 0.4472135954999579*m1rOther[11]*mnuOther; 
  data->AEM_S(27,46) = 0.4*m1rOther[15]*mnuOther+0.447213595499958*m1rOther[9]*mnuOther; 
  data->AEM_S(27,47) = 0.4*m1rOther[14]*mnuOther+0.447213595499958*m1rOther[10]*mnuOther; 
  data->AEM_S(28,40) = 0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(28,41) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(28,42) = 0.5000000000000001*m1rOther[14]*mnuOther; 
  data->AEM_S(28,43) = 0.4472135954999579*m1rOther[11]*mnuOther; 
  data->AEM_S(28,44) = 0.31943828249997*m1rOther[12]*mnuOther+0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(28,46) = 0.31943828249997*m1rOther[14]*mnuOther+0.5000000000000001*m1rOther[10]*mnuOther; 
  data->AEM_S(28,47) = 0.4472135954999579*m1rOther[15]*mnuOther; 
  data->AEM_S(29,40) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(29,41) = 0.5000000000000001*m1rOther[15]*mnuOther; 
  data->AEM_S(29,42) = 0.4472135954999579*m1rOther[10]*mnuOther; 
  data->AEM_S(29,43) = 0.4472135954999579*m1rOther[11]*mnuOther; 
  data->AEM_S(29,45) = 0.31943828249997*m1rOther[13]*mnuOther+0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(29,46) = 0.4472135954999579*m1rOther[14]*mnuOther; 
  data->AEM_S(29,47) = 0.31943828249997*m1rOther[15]*mnuOther+0.5000000000000001*m1rOther[9]*mnuOther; 
  data->AEM_S(30,40) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(30,41) = 0.447213595499958*m1rOther[11]*mnuOther; 
  data->AEM_S(30,42) = 0.5000000000000001*m1rOther[12]*mnuOther; 
  data->AEM_S(30,43) = 0.4*m1rOther[15]*mnuOther+0.447213595499958*m1rOther[9]*mnuOther; 
  data->AEM_S(30,44) = 0.31943828249997*m1rOther[14]*mnuOther+0.5000000000000001*m1rOther[10]*mnuOther; 
  data->AEM_S(30,45) = 0.4472135954999579*m1rOther[14]*mnuOther; 
  data->AEM_S(30,46) = 0.4472135954999579*m1rOther[13]*mnuOther+0.31943828249997*m1rOther[12]*mnuOther+0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(30,47) = 0.4*m1rOther[11]*mnuOther; 
  data->AEM_S(31,40) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(31,41) = 0.5000000000000001*m1rOther[13]*mnuOther; 
  data->AEM_S(31,42) = 0.447213595499958*m1rOther[11]*mnuOther; 
  data->AEM_S(31,43) = 0.4*m1rOther[14]*mnuOther+0.447213595499958*m1rOther[10]*mnuOther; 
  data->AEM_S(31,44) = 0.4472135954999579*m1rOther[15]*mnuOther; 
  data->AEM_S(31,45) = 0.31943828249997*m1rOther[15]*mnuOther+0.5000000000000001*m1rOther[9]*mnuOther; 
  data->AEM_S(31,46) = 0.4*m1rOther[11]*mnuOther; 
  data->AEM_S(31,47) = 0.31943828249997*m1rOther[13]*mnuOther+0.4472135954999579*m1rOther[12]*mnuOther+0.5*m1rOther[8]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
  mnuM1sum[9] += m1rSelf[9]*mnuSelf+m1rOther[9]*mnuOther; 
  mnuM1sum[10] += m1rSelf[10]*mnuSelf+m1rOther[10]*mnuOther; 
  mnuM1sum[11] += m1rSelf[11]*mnuSelf+m1rOther[11]*mnuOther; 
  mnuM1sum[12] += m1rSelf[12]*mnuSelf+m1rOther[12]*mnuOther; 
  mnuM1sum[13] += m1rSelf[13]*mnuSelf+m1rOther[13]*mnuOther; 
  mnuM1sum[14] += m1rSelf[14]*mnuSelf+m1rOther[14]*mnuOther; 
  mnuM1sum[15] += m1rSelf[15]*mnuSelf+m1rOther[15]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(16,16) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(16,17) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,18) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,19) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,20) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(16,21) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(16,22) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(16,23) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(17,16) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,17) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(17,18) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,19) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,20) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,21) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(17,22) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,23) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(18,16) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,17) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(18,18) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,19) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(18,20) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(18,21) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,22) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(18,23) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,16) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,17) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,18) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,19) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(19,20) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,21) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,22) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,23) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,16) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(20,17) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,18) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(20,19) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(20,20) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,22) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,23) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(21,16) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(21,17) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(21,18) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,19) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,21) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(21,22) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(21,23) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,16) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(22,17) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,18) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(22,19) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,20) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,21) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(22,22) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(22,23) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,16) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(23,17) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(23,18) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,19) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,20) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(23,21) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,22) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,23) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(16,24) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(16,25) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(16,26) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(16,27) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(16,28) = -0.5*cMSelf[20]*mnuSelf; 
  data->AEM_S(16,29) = -0.5*cMSelf[21]*mnuSelf; 
  data->AEM_S(16,30) = -0.5*cMSelf[22]*mnuSelf; 
  data->AEM_S(16,31) = -0.5*cMSelf[23]*mnuSelf; 
  data->AEM_S(17,24) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(17,25) = (-0.4472135954999579*cMSelf[20]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(17,26) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(17,27) = (-0.447213595499958*cMSelf[22]*mnuSelf)-0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(17,28) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(17,29) = -0.5000000000000001*cMSelf[23]*mnuSelf; 
  data->AEM_S(17,30) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(17,31) = -0.5000000000000001*cMSelf[21]*mnuSelf; 
  data->AEM_S(18,24) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(18,25) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(18,26) = (-0.4472135954999579*cMSelf[21]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(18,27) = (-0.447213595499958*cMSelf[23]*mnuSelf)-0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(18,28) = -0.5000000000000001*cMSelf[22]*mnuSelf; 
  data->AEM_S(18,29) = -0.4472135954999579*cMSelf[18]*mnuSelf; 
  data->AEM_S(18,30) = -0.5000000000000001*cMSelf[20]*mnuSelf; 
  data->AEM_S(18,31) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(19,24) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(19,25) = (-0.447213595499958*cMSelf[22]*mnuSelf)-0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(19,26) = (-0.447213595499958*cMSelf[23]*mnuSelf)-0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(19,27) = (-0.4472135954999579*cMSelf[21]*mnuSelf)-0.4472135954999579*cMSelf[20]*mnuSelf-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(19,28) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(19,29) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(19,30) = (-0.4*cMSelf[23]*mnuSelf)-0.447213595499958*cMSelf[17]*mnuSelf; 
  data->AEM_S(19,31) = (-0.4*cMSelf[22]*mnuSelf)-0.447213595499958*cMSelf[18]*mnuSelf; 
  data->AEM_S(20,24) = -0.5*cMSelf[20]*mnuSelf; 
  data->AEM_S(20,25) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(20,26) = -0.5000000000000001*cMSelf[22]*mnuSelf; 
  data->AEM_S(20,27) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(20,28) = (-0.31943828249997*cMSelf[20]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(20,30) = (-0.31943828249997*cMSelf[22]*mnuSelf)-0.5000000000000001*cMSelf[18]*mnuSelf; 
  data->AEM_S(20,31) = -0.4472135954999579*cMSelf[23]*mnuSelf; 
  data->AEM_S(21,24) = -0.5*cMSelf[21]*mnuSelf; 
  data->AEM_S(21,25) = -0.5000000000000001*cMSelf[23]*mnuSelf; 
  data->AEM_S(21,26) = -0.4472135954999579*cMSelf[18]*mnuSelf; 
  data->AEM_S(21,27) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(21,29) = (-0.31943828249997*cMSelf[21]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(21,30) = -0.4472135954999579*cMSelf[22]*mnuSelf; 
  data->AEM_S(21,31) = (-0.31943828249997*cMSelf[23]*mnuSelf)-0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(22,24) = -0.5*cMSelf[22]*mnuSelf; 
  data->AEM_S(22,25) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(22,26) = -0.5000000000000001*cMSelf[20]*mnuSelf; 
  data->AEM_S(22,27) = (-0.4*cMSelf[23]*mnuSelf)-0.447213595499958*cMSelf[17]*mnuSelf; 
  data->AEM_S(22,28) = (-0.31943828249997*cMSelf[22]*mnuSelf)-0.5000000000000001*cMSelf[18]*mnuSelf; 
  data->AEM_S(22,29) = -0.4472135954999579*cMSelf[22]*mnuSelf; 
  data->AEM_S(22,30) = (-0.4472135954999579*cMSelf[21]*mnuSelf)-0.31943828249997*cMSelf[20]*mnuSelf-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(22,31) = -0.4*cMSelf[19]*mnuSelf; 
  data->AEM_S(23,24) = -0.5*cMSelf[23]*mnuSelf; 
  data->AEM_S(23,25) = -0.5000000000000001*cMSelf[21]*mnuSelf; 
  data->AEM_S(23,26) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(23,27) = (-0.4*cMSelf[22]*mnuSelf)-0.447213595499958*cMSelf[18]*mnuSelf; 
  data->AEM_S(23,28) = -0.4472135954999579*cMSelf[23]*mnuSelf; 
  data->AEM_S(23,29) = (-0.31943828249997*cMSelf[23]*mnuSelf)-0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(23,30) = -0.4*cMSelf[19]*mnuSelf; 
  data->AEM_S(23,31) = (-0.31943828249997*cMSelf[21]*mnuSelf)-0.4472135954999579*cMSelf[20]*mnuSelf-0.5*cMSelf[16]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(16,48) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(16,49) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(16,50) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(16,51) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(16,52) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(16,53) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(16,54) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(16,55) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(17,48) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(17,49) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(17,50) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(17,51) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(17,52) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(17,53) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(17,54) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(17,55) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(18,48) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(18,49) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(18,50) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(18,51) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(18,52) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(18,53) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(18,54) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(18,55) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(19,48) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(19,49) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(19,50) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(19,51) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(19,52) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(19,53) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(19,54) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(19,55) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(20,48) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(20,49) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(20,50) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(20,51) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(20,52) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(20,54) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(20,55) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(21,48) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(21,49) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(21,50) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(21,51) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(21,53) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(21,54) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(21,55) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(22,48) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(22,49) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(22,50) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(22,51) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(22,52) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(22,53) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(22,54) = 0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(22,55) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(23,48) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(23,49) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(23,50) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(23,51) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(23,52) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(23,53) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(23,54) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(23,55) = 0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(16,56) = -0.5*cMOther[16]*mnuOther; 
  data->AEM_S(16,57) = -0.5*cMOther[17]*mnuOther; 
  data->AEM_S(16,58) = -0.5*cMOther[18]*mnuOther; 
  data->AEM_S(16,59) = -0.5*cMOther[19]*mnuOther; 
  data->AEM_S(16,60) = -0.5*cMOther[20]*mnuOther; 
  data->AEM_S(16,61) = -0.5*cMOther[21]*mnuOther; 
  data->AEM_S(16,62) = -0.5*cMOther[22]*mnuOther; 
  data->AEM_S(16,63) = -0.5*cMOther[23]*mnuOther; 
  data->AEM_S(17,56) = -0.5*cMOther[17]*mnuOther; 
  data->AEM_S(17,57) = (-0.4472135954999579*cMOther[20]*mnuOther)-0.5*cMOther[16]*mnuOther; 
  data->AEM_S(17,58) = -0.5*cMOther[19]*mnuOther; 
  data->AEM_S(17,59) = (-0.447213595499958*cMOther[22]*mnuOther)-0.5*cMOther[18]*mnuOther; 
  data->AEM_S(17,60) = -0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(17,61) = -0.5000000000000001*cMOther[23]*mnuOther; 
  data->AEM_S(17,62) = -0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(17,63) = -0.5000000000000001*cMOther[21]*mnuOther; 
  data->AEM_S(18,56) = -0.5*cMOther[18]*mnuOther; 
  data->AEM_S(18,57) = -0.5*cMOther[19]*mnuOther; 
  data->AEM_S(18,58) = (-0.4472135954999579*cMOther[21]*mnuOther)-0.5*cMOther[16]*mnuOther; 
  data->AEM_S(18,59) = (-0.447213595499958*cMOther[23]*mnuOther)-0.5*cMOther[17]*mnuOther; 
  data->AEM_S(18,60) = -0.5000000000000001*cMOther[22]*mnuOther; 
  data->AEM_S(18,61) = -0.4472135954999579*cMOther[18]*mnuOther; 
  data->AEM_S(18,62) = -0.5000000000000001*cMOther[20]*mnuOther; 
  data->AEM_S(18,63) = -0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(19,56) = -0.5*cMOther[19]*mnuOther; 
  data->AEM_S(19,57) = (-0.447213595499958*cMOther[22]*mnuOther)-0.5*cMOther[18]*mnuOther; 
  data->AEM_S(19,58) = (-0.447213595499958*cMOther[23]*mnuOther)-0.5*cMOther[17]*mnuOther; 
  data->AEM_S(19,59) = (-0.4472135954999579*cMOther[21]*mnuOther)-0.4472135954999579*cMOther[20]*mnuOther-0.5*cMOther[16]*mnuOther; 
  data->AEM_S(19,60) = -0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(19,61) = -0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(19,62) = (-0.4*cMOther[23]*mnuOther)-0.447213595499958*cMOther[17]*mnuOther; 
  data->AEM_S(19,63) = (-0.4*cMOther[22]*mnuOther)-0.447213595499958*cMOther[18]*mnuOther; 
  data->AEM_S(20,56) = -0.5*cMOther[20]*mnuOther; 
  data->AEM_S(20,57) = -0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(20,58) = -0.5000000000000001*cMOther[22]*mnuOther; 
  data->AEM_S(20,59) = -0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(20,60) = (-0.31943828249997*cMOther[20]*mnuOther)-0.5*cMOther[16]*mnuOther; 
  data->AEM_S(20,62) = (-0.31943828249997*cMOther[22]*mnuOther)-0.5000000000000001*cMOther[18]*mnuOther; 
  data->AEM_S(20,63) = -0.4472135954999579*cMOther[23]*mnuOther; 
  data->AEM_S(21,56) = -0.5*cMOther[21]*mnuOther; 
  data->AEM_S(21,57) = -0.5000000000000001*cMOther[23]*mnuOther; 
  data->AEM_S(21,58) = -0.4472135954999579*cMOther[18]*mnuOther; 
  data->AEM_S(21,59) = -0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(21,61) = (-0.31943828249997*cMOther[21]*mnuOther)-0.5*cMOther[16]*mnuOther; 
  data->AEM_S(21,62) = -0.4472135954999579*cMOther[22]*mnuOther; 
  data->AEM_S(21,63) = (-0.31943828249997*cMOther[23]*mnuOther)-0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(22,56) = -0.5*cMOther[22]*mnuOther; 
  data->AEM_S(22,57) = -0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(22,58) = -0.5000000000000001*cMOther[20]*mnuOther; 
  data->AEM_S(22,59) = (-0.4*cMOther[23]*mnuOther)-0.447213595499958*cMOther[17]*mnuOther; 
  data->AEM_S(22,60) = (-0.31943828249997*cMOther[22]*mnuOther)-0.5000000000000001*cMOther[18]*mnuOther; 
  data->AEM_S(22,61) = -0.4472135954999579*cMOther[22]*mnuOther; 
  data->AEM_S(22,62) = (-0.4472135954999579*cMOther[21]*mnuOther)-0.31943828249997*cMOther[20]*mnuOther-0.5*cMOther[16]*mnuOther; 
  data->AEM_S(22,63) = -0.4*cMOther[19]*mnuOther; 
  data->AEM_S(23,56) = -0.5*cMOther[23]*mnuOther; 
  data->AEM_S(23,57) = -0.5000000000000001*cMOther[21]*mnuOther; 
  data->AEM_S(23,58) = -0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(23,59) = (-0.4*cMOther[22]*mnuOther)-0.447213595499958*cMOther[18]*mnuOther; 
  data->AEM_S(23,60) = -0.4472135954999579*cMOther[23]*mnuOther; 
  data->AEM_S(23,61) = (-0.31943828249997*cMOther[23]*mnuOther)-0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(23,62) = -0.4*cMOther[19]*mnuOther; 
  data->AEM_S(23,63) = (-0.31943828249997*cMOther[21]*mnuOther)-0.4472135954999579*cMOther[20]*mnuOther-0.5*cMOther[16]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfZ and uCrossSelfZ ... // 
  data->AEM_S(24,16) = 0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(24,17) = 0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(24,18) = 0.5*m1rSelf[18]*mnuSelf; 
  data->AEM_S(24,19) = 0.5*m1rSelf[19]*mnuSelf; 
  data->AEM_S(24,20) = 0.5*m1rSelf[20]*mnuSelf; 
  data->AEM_S(24,21) = 0.5*m1rSelf[21]*mnuSelf; 
  data->AEM_S(24,22) = 0.5*m1rSelf[22]*mnuSelf; 
  data->AEM_S(24,23) = 0.5*m1rSelf[23]*mnuSelf; 
  data->AEM_S(25,16) = 0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(25,17) = 0.4472135954999579*m1rSelf[20]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(25,18) = 0.5*m1rSelf[19]*mnuSelf; 
  data->AEM_S(25,19) = 0.447213595499958*m1rSelf[22]*mnuSelf+0.5*m1rSelf[18]*mnuSelf; 
  data->AEM_S(25,20) = 0.4472135954999579*m1rSelf[17]*mnuSelf; 
  data->AEM_S(25,21) = 0.5000000000000001*m1rSelf[23]*mnuSelf; 
  data->AEM_S(25,22) = 0.447213595499958*m1rSelf[19]*mnuSelf; 
  data->AEM_S(25,23) = 0.5000000000000001*m1rSelf[21]*mnuSelf; 
  data->AEM_S(26,16) = 0.5*m1rSelf[18]*mnuSelf; 
  data->AEM_S(26,17) = 0.5*m1rSelf[19]*mnuSelf; 
  data->AEM_S(26,18) = 0.4472135954999579*m1rSelf[21]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(26,19) = 0.447213595499958*m1rSelf[23]*mnuSelf+0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(26,20) = 0.5000000000000001*m1rSelf[22]*mnuSelf; 
  data->AEM_S(26,21) = 0.4472135954999579*m1rSelf[18]*mnuSelf; 
  data->AEM_S(26,22) = 0.5000000000000001*m1rSelf[20]*mnuSelf; 
  data->AEM_S(26,23) = 0.447213595499958*m1rSelf[19]*mnuSelf; 
  data->AEM_S(27,16) = 0.5*m1rSelf[19]*mnuSelf; 
  data->AEM_S(27,17) = 0.447213595499958*m1rSelf[22]*mnuSelf+0.5*m1rSelf[18]*mnuSelf; 
  data->AEM_S(27,18) = 0.447213595499958*m1rSelf[23]*mnuSelf+0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(27,19) = 0.4472135954999579*m1rSelf[21]*mnuSelf+0.4472135954999579*m1rSelf[20]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(27,20) = 0.4472135954999579*m1rSelf[19]*mnuSelf; 
  data->AEM_S(27,21) = 0.4472135954999579*m1rSelf[19]*mnuSelf; 
  data->AEM_S(27,22) = 0.4*m1rSelf[23]*mnuSelf+0.447213595499958*m1rSelf[17]*mnuSelf; 
  data->AEM_S(27,23) = 0.4*m1rSelf[22]*mnuSelf+0.447213595499958*m1rSelf[18]*mnuSelf; 
  data->AEM_S(28,16) = 0.5*m1rSelf[20]*mnuSelf; 
  data->AEM_S(28,17) = 0.4472135954999579*m1rSelf[17]*mnuSelf; 
  data->AEM_S(28,18) = 0.5000000000000001*m1rSelf[22]*mnuSelf; 
  data->AEM_S(28,19) = 0.4472135954999579*m1rSelf[19]*mnuSelf; 
  data->AEM_S(28,20) = 0.31943828249997*m1rSelf[20]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(28,22) = 0.31943828249997*m1rSelf[22]*mnuSelf+0.5000000000000001*m1rSelf[18]*mnuSelf; 
  data->AEM_S(28,23) = 0.4472135954999579*m1rSelf[23]*mnuSelf; 
  data->AEM_S(29,16) = 0.5*m1rSelf[21]*mnuSelf; 
  data->AEM_S(29,17) = 0.5000000000000001*m1rSelf[23]*mnuSelf; 
  data->AEM_S(29,18) = 0.4472135954999579*m1rSelf[18]*mnuSelf; 
  data->AEM_S(29,19) = 0.4472135954999579*m1rSelf[19]*mnuSelf; 
  data->AEM_S(29,21) = 0.31943828249997*m1rSelf[21]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(29,22) = 0.4472135954999579*m1rSelf[22]*mnuSelf; 
  data->AEM_S(29,23) = 0.31943828249997*m1rSelf[23]*mnuSelf+0.5000000000000001*m1rSelf[17]*mnuSelf; 
  data->AEM_S(30,16) = 0.5*m1rSelf[22]*mnuSelf; 
  data->AEM_S(30,17) = 0.447213595499958*m1rSelf[19]*mnuSelf; 
  data->AEM_S(30,18) = 0.5000000000000001*m1rSelf[20]*mnuSelf; 
  data->AEM_S(30,19) = 0.4*m1rSelf[23]*mnuSelf+0.447213595499958*m1rSelf[17]*mnuSelf; 
  data->AEM_S(30,20) = 0.31943828249997*m1rSelf[22]*mnuSelf+0.5000000000000001*m1rSelf[18]*mnuSelf; 
  data->AEM_S(30,21) = 0.4472135954999579*m1rSelf[22]*mnuSelf; 
  data->AEM_S(30,22) = 0.4472135954999579*m1rSelf[21]*mnuSelf+0.31943828249997*m1rSelf[20]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(30,23) = 0.4*m1rSelf[19]*mnuSelf; 
  data->AEM_S(31,16) = 0.5*m1rSelf[23]*mnuSelf; 
  data->AEM_S(31,17) = 0.5000000000000001*m1rSelf[21]*mnuSelf; 
  data->AEM_S(31,18) = 0.447213595499958*m1rSelf[19]*mnuSelf; 
  data->AEM_S(31,19) = 0.4*m1rSelf[22]*mnuSelf+0.447213595499958*m1rSelf[18]*mnuSelf; 
  data->AEM_S(31,20) = 0.4472135954999579*m1rSelf[23]*mnuSelf; 
  data->AEM_S(31,21) = 0.31943828249997*m1rSelf[23]*mnuSelf+0.5000000000000001*m1rSelf[17]*mnuSelf; 
  data->AEM_S(31,22) = 0.4*m1rSelf[19]*mnuSelf; 
  data->AEM_S(31,23) = 0.31943828249997*m1rSelf[21]*mnuSelf+0.4472135954999579*m1rSelf[20]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherZ and uCrossOtherZ ... // 
  data->AEM_S(24,48) = 0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(24,49) = 0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(24,50) = 0.5*m1rOther[18]*mnuOther; 
  data->AEM_S(24,51) = 0.5*m1rOther[19]*mnuOther; 
  data->AEM_S(24,52) = 0.5*m1rOther[20]*mnuOther; 
  data->AEM_S(24,53) = 0.5*m1rOther[21]*mnuOther; 
  data->AEM_S(24,54) = 0.5*m1rOther[22]*mnuOther; 
  data->AEM_S(24,55) = 0.5*m1rOther[23]*mnuOther; 
  data->AEM_S(25,48) = 0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(25,49) = 0.4472135954999579*m1rOther[20]*mnuOther+0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(25,50) = 0.5*m1rOther[19]*mnuOther; 
  data->AEM_S(25,51) = 0.447213595499958*m1rOther[22]*mnuOther+0.5*m1rOther[18]*mnuOther; 
  data->AEM_S(25,52) = 0.4472135954999579*m1rOther[17]*mnuOther; 
  data->AEM_S(25,53) = 0.5000000000000001*m1rOther[23]*mnuOther; 
  data->AEM_S(25,54) = 0.447213595499958*m1rOther[19]*mnuOther; 
  data->AEM_S(25,55) = 0.5000000000000001*m1rOther[21]*mnuOther; 
  data->AEM_S(26,48) = 0.5*m1rOther[18]*mnuOther; 
  data->AEM_S(26,49) = 0.5*m1rOther[19]*mnuOther; 
  data->AEM_S(26,50) = 0.4472135954999579*m1rOther[21]*mnuOther+0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(26,51) = 0.447213595499958*m1rOther[23]*mnuOther+0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(26,52) = 0.5000000000000001*m1rOther[22]*mnuOther; 
  data->AEM_S(26,53) = 0.4472135954999579*m1rOther[18]*mnuOther; 
  data->AEM_S(26,54) = 0.5000000000000001*m1rOther[20]*mnuOther; 
  data->AEM_S(26,55) = 0.447213595499958*m1rOther[19]*mnuOther; 
  data->AEM_S(27,48) = 0.5*m1rOther[19]*mnuOther; 
  data->AEM_S(27,49) = 0.447213595499958*m1rOther[22]*mnuOther+0.5*m1rOther[18]*mnuOther; 
  data->AEM_S(27,50) = 0.447213595499958*m1rOther[23]*mnuOther+0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(27,51) = 0.4472135954999579*m1rOther[21]*mnuOther+0.4472135954999579*m1rOther[20]*mnuOther+0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(27,52) = 0.4472135954999579*m1rOther[19]*mnuOther; 
  data->AEM_S(27,53) = 0.4472135954999579*m1rOther[19]*mnuOther; 
  data->AEM_S(27,54) = 0.4*m1rOther[23]*mnuOther+0.447213595499958*m1rOther[17]*mnuOther; 
  data->AEM_S(27,55) = 0.4*m1rOther[22]*mnuOther+0.447213595499958*m1rOther[18]*mnuOther; 
  data->AEM_S(28,48) = 0.5*m1rOther[20]*mnuOther; 
  data->AEM_S(28,49) = 0.4472135954999579*m1rOther[17]*mnuOther; 
  data->AEM_S(28,50) = 0.5000000000000001*m1rOther[22]*mnuOther; 
  data->AEM_S(28,51) = 0.4472135954999579*m1rOther[19]*mnuOther; 
  data->AEM_S(28,52) = 0.31943828249997*m1rOther[20]*mnuOther+0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(28,54) = 0.31943828249997*m1rOther[22]*mnuOther+0.5000000000000001*m1rOther[18]*mnuOther; 
  data->AEM_S(28,55) = 0.4472135954999579*m1rOther[23]*mnuOther; 
  data->AEM_S(29,48) = 0.5*m1rOther[21]*mnuOther; 
  data->AEM_S(29,49) = 0.5000000000000001*m1rOther[23]*mnuOther; 
  data->AEM_S(29,50) = 0.4472135954999579*m1rOther[18]*mnuOther; 
  data->AEM_S(29,51) = 0.4472135954999579*m1rOther[19]*mnuOther; 
  data->AEM_S(29,53) = 0.31943828249997*m1rOther[21]*mnuOther+0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(29,54) = 0.4472135954999579*m1rOther[22]*mnuOther; 
  data->AEM_S(29,55) = 0.31943828249997*m1rOther[23]*mnuOther+0.5000000000000001*m1rOther[17]*mnuOther; 
  data->AEM_S(30,48) = 0.5*m1rOther[22]*mnuOther; 
  data->AEM_S(30,49) = 0.447213595499958*m1rOther[19]*mnuOther; 
  data->AEM_S(30,50) = 0.5000000000000001*m1rOther[20]*mnuOther; 
  data->AEM_S(30,51) = 0.4*m1rOther[23]*mnuOther+0.447213595499958*m1rOther[17]*mnuOther; 
  data->AEM_S(30,52) = 0.31943828249997*m1rOther[22]*mnuOther+0.5000000000000001*m1rOther[18]*mnuOther; 
  data->AEM_S(30,53) = 0.4472135954999579*m1rOther[22]*mnuOther; 
  data->AEM_S(30,54) = 0.4472135954999579*m1rOther[21]*mnuOther+0.31943828249997*m1rOther[20]*mnuOther+0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(30,55) = 0.4*m1rOther[19]*mnuOther; 
  data->AEM_S(31,48) = 0.5*m1rOther[23]*mnuOther; 
  data->AEM_S(31,49) = 0.5000000000000001*m1rOther[21]*mnuOther; 
  data->AEM_S(31,50) = 0.447213595499958*m1rOther[19]*mnuOther; 
  data->AEM_S(31,51) = 0.4*m1rOther[22]*mnuOther+0.447213595499958*m1rOther[18]*mnuOther; 
  data->AEM_S(31,52) = 0.4472135954999579*m1rOther[23]*mnuOther; 
  data->AEM_S(31,53) = 0.31943828249997*m1rOther[23]*mnuOther+0.5000000000000001*m1rOther[17]*mnuOther; 
  data->AEM_S(31,54) = 0.4*m1rOther[19]*mnuOther; 
  data->AEM_S(31,55) = 0.31943828249997*m1rOther[21]*mnuOther+0.4472135954999579*m1rOther[20]*mnuOther+0.5*m1rOther[16]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[16] += m1rSelf[16]*mnuSelf+m1rOther[16]*mnuOther; 
  mnuM1sum[17] += m1rSelf[17]*mnuSelf+m1rOther[17]*mnuOther; 
  mnuM1sum[18] += m1rSelf[18]*mnuSelf+m1rOther[18]*mnuOther; 
  mnuM1sum[19] += m1rSelf[19]*mnuSelf+m1rOther[19]*mnuOther; 
  mnuM1sum[20] += m1rSelf[20]*mnuSelf+m1rOther[20]*mnuOther; 
  mnuM1sum[21] += m1rSelf[21]*mnuSelf+m1rOther[21]*mnuOther; 
  mnuM1sum[22] += m1rSelf[22]*mnuSelf+m1rOther[22]*mnuOther; 
  mnuM1sum[23] += m1rSelf[23]*mnuSelf+m1rOther[23]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(24,24) = 1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(24,25) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(24,26) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(24,27) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(24,28) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(24,29) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(24,30) = 1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(24,31) = 1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(25,24) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(25,25) = 1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(25,26) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(25,27) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(25,28) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(25,29) = 1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(25,30) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(25,31) = 1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(26,24) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(26,25) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(26,26) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(26,27) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(26,28) = 1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(26,29) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(26,30) = 1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(26,31) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(27,24) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(27,25) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(27,26) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(27,27) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(27,28) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(27,29) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(27,30) = 1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(27,31) = 1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(28,24) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(28,25) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(28,26) = 1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(28,27) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(28,28) = 0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(28,30) = 0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(28,31) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(29,24) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(29,25) = 1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(29,26) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(29,27) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(29,29) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(29,30) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(29,31) = 0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(30,24) = 1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(30,25) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(30,26) = 1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(30,27) = 1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(30,28) = 0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(30,29) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(30,30) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(30,31) = 1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(31,24) = 1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(31,25) = 1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(31,26) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(31,27) = 1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(31,28) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(31,29) = 0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(31,30) = 1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(31,31) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(24,56) = 1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(24,57) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(24,58) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(24,59) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(24,60) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(24,61) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(24,62) = 1.5*m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(24,63) = 1.5*m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(25,56) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(25,57) = 1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(25,58) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(25,59) = 1.341640786499874*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(25,60) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(25,61) = 1.5*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(25,62) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(25,63) = 1.5*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(26,56) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(26,57) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(26,58) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(26,59) = 1.341640786499874*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(26,60) = 1.5*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(26,61) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(26,62) = 1.5*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(26,63) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(27,56) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(27,57) = 1.341640786499874*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(27,58) = 1.341640786499874*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(27,59) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(27,60) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(27,61) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(27,62) = 1.2*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(27,63) = 1.2*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(28,56) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(28,57) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(28,58) = 1.5*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(28,59) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(28,60) = 0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(28,62) = 0.9583148474999099*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(28,63) = 1.341640786499874*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(29,56) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(29,57) = 1.5*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(29,58) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(29,59) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(29,61) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(29,62) = 1.341640786499874*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(29,63) = 0.9583148474999099*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(30,56) = 1.5*m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(30,57) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(30,58) = 1.5*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(30,59) = 1.2*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(30,60) = 0.9583148474999099*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(30,61) = 1.341640786499874*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(30,62) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(30,63) = 1.2*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(31,56) = 1.5*m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(31,57) = 1.5*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(31,58) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(31,59) = 1.2*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(31,60) = 1.341640786499874*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(31,61) = 0.9583148474999099*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(31,62) = 1.2*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(31,63) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[8]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2rSelf[0]*mnuSelf+m2rOther[0]*mnuOther; 
  mnuM2sum[1] = m2rSelf[1]*mnuSelf+m2rOther[1]*mnuOther; 
  mnuM2sum[2] = m2rSelf[2]*mnuSelf+m2rOther[2]*mnuOther; 
  mnuM2sum[3] = m2rSelf[3]*mnuSelf+m2rOther[3]*mnuOther; 
  mnuM2sum[4] = m2rSelf[4]*mnuSelf+m2rOther[4]*mnuOther; 
  mnuM2sum[5] = m2rSelf[5]*mnuSelf+m2rOther[5]*mnuOther; 
  mnuM2sum[6] = m2rSelf[6]*mnuSelf+m2rOther[6]*mnuOther; 
  mnuM2sum[7] = m2rSelf[7]*mnuSelf+m2rOther[7]*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<8,16>(0,8).setZero(); 
  data->AEM_S.block<16,8>(8,0).setZero(); 
  data->AEM_S.block<8,8>(8,16).setZero(); 
  data->AEM_S.block<8,8>(16,8).setZero(); 
  data->AEM_S.block<8,16>(0,40).setZero(); 
  data->AEM_S.block<16,8>(8,32).setZero(); 
  data->AEM_S.block<8,8>(8,48).setZero(); 
  data->AEM_S.block<8,8>(16,40).setZero(); 
 
  double m1Relax[24]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<24; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  double m1EffD[24]; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(32,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(32,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(32,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(32,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(32,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(32,6) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(32,7) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(33,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(33,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,3) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,5) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(33,6) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,7) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(34,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(34,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(34,3) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(34,4) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(34,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(34,6) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(34,7) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,1) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(35,2) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(35,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(35,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,6) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(35,7) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(36,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(36,2) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(36,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(36,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(36,6) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,7) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(37,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(37,1) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(37,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(37,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(37,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(37,6) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(37,7) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(38,0) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(38,1) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(38,2) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(38,3) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(38,4) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(38,5) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(38,6) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(38,7) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(39,0) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(39,1) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(39,2) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(39,3) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(39,4) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(39,5) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(39,6) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(39,7) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(32,24) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(32,25) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(32,26) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(32,27) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(32,28) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(32,29) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(32,30) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(32,31) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(33,24) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(33,25) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(33,26) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(33,27) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(33,28) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(33,29) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(33,30) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(33,31) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(34,24) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(34,25) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(34,26) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(34,27) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(34,28) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(34,29) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(34,30) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(34,31) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(35,24) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(35,25) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(35,26) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(35,27) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(35,28) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(35,29) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(35,30) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(35,31) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(36,24) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(36,25) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(36,26) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(36,27) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(36,28) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(36,30) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(36,31) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(37,24) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(37,25) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(37,26) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(37,27) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(37,29) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(37,30) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(37,31) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(38,24) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(38,25) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(38,26) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(38,27) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(38,28) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(38,29) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(38,30) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(38,31) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(39,24) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(39,25) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(39,26) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(39,27) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(39,28) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(39,29) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(39,30) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(39,31) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(32,32) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(32,33) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(32,34) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(32,35) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(32,36) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(32,37) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(32,38) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(32,39) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(33,32) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(33,33) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(33,34) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(33,35) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(33,36) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(33,37) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(33,38) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(33,39) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(34,32) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(34,33) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(34,34) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(34,35) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(34,36) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(34,37) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(34,38) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(34,39) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(35,32) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(35,33) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(35,34) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(35,35) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(35,36) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(35,37) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(35,38) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(35,39) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(36,32) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(36,33) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(36,34) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(36,35) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(36,36) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(36,38) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(36,39) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(37,32) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(37,33) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(37,34) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(37,35) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(37,37) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(37,38) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(37,39) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(38,32) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(38,33) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(38,34) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(38,35) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(38,36) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(38,37) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(38,38) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.31943828249997*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(38,39) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(39,32) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(39,33) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(39,34) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(39,35) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(39,36) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(39,37) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(39,38) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(39,39) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(32,56) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(32,57) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(32,58) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(32,59) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(32,60) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(32,61) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(32,62) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(32,63) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(33,56) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(33,57) = 0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(33,58) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(33,59) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(33,60) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(33,61) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(33,62) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(33,63) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(34,56) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(34,57) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(34,58) = 0.4472135954999579*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(34,59) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(34,60) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(34,61) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(34,62) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(34,63) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(35,56) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(35,57) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(35,58) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(35,59) = 0.4472135954999579*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(35,60) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(35,61) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(35,62) = 0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(35,63) = 0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(36,56) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(36,57) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(36,58) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(36,59) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(36,60) = 0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(36,62) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(36,63) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(37,56) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(37,57) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(37,58) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(37,59) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(37,61) = 0.31943828249997*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(37,62) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(37,63) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(38,56) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(38,57) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(38,58) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(38,59) = 0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(38,60) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(38,61) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(38,62) = 0.4472135954999579*cMOther[5]*mnuOther+0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(38,63) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(39,56) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(39,57) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(39,58) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(39,59) = 0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(39,60) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(39,61) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(39,62) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(39,63) = 0.31943828249997*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(56,0) = (-0.25*m0rSelf[7]*uSelf[7]*mnuSelf)-0.25*m0rSelf[6]*uSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(56,1) = (-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(56,2) = (-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,3) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(56,4) = (-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(56,5) = (-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(56,6) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(56,7) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(57,0) = (-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(57,1) = (-0.45*m0rSelf[7]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(57,2) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(57,3) = (-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf)-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(57,4) = (-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(57,5) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(57,6) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(57,7) = (-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf)-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(58,0) = (-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(58,1) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(58,2) = (-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.45*m0rSelf[6]*uSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(58,3) = (-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(58,4) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(58,5) = (-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(58,6) = (-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(58,7) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(59,0) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(59,1) = (-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf)-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(59,2) = (-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(59,3) = (-0.7071428571428572*m0rSelf[7]*uSelf[7]*mnuSelf)-0.4024922359499621*m0rSelf[1]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[6]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.2*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.2*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.81*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(59,4) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(59,5) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(59,6) = (-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(59,7) = (-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf)-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(60,0) = (-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(60,1) = (-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(60,2) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(60,3) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(60,4) = (-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.5357142857142857*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[6]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(60,5) = (-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(60,6) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(60,7) = (-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(61,0) = (-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(61,1) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(61,2) = (-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(61,3) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(61,4) = (-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(61,5) = (-0.5357142857142857*m0rSelf[7]*uSelf[7]*mnuSelf)-0.159719141249985*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(61,6) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(61,7) = (-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(62,0) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(62,1) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(62,2) = (-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(62,3) = (-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(62,4) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(62,5) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(62,6) = (-0.6173469387755102*m0rSelf[7]*uSelf[7]*mnuSelf)-0.351382110749967*m0rSelf[1]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[7]*mnuSelf-0.9642857142857143*m0rSelf[6]*uSelf[6]*mnuSelf-0.2874944542499729*m0rSelf[2]*uSelf[6]*mnuSelf-0.2874944542499729*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(62,7) = (-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(63,0) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(63,1) = (-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf)-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(63,2) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(63,3) = (-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf)-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(63,4) = (-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(63,5) = (-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(63,6) = (-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(63,7) = (-0.9642857142857143*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2874944542499729*m0rSelf[1]*uSelf[7]*mnuSelf-0.2874944542499729*uSelf[1]*m0rSelf[7]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[6]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(56,32) = 0.25*m0rOther[7]*uOther[7]*mnuOther+0.25*m0rOther[6]*uOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(56,33) = 0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(56,34) = 0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(56,35) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(56,36) = 0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(56,37) = 0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(56,38) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(56,39) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(57,32) = 0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(57,33) = 0.45*m0rOther[7]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(57,34) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(57,35) = 0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(57,36) = 0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(57,37) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(57,38) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(57,39) = 0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(58,32) = 0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(58,33) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(58,34) = 0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.45*m0rOther[6]*uOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(58,35) = 0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(58,36) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(58,37) = 0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(58,38) = 0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(58,39) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(59,32) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(59,33) = 0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(59,34) = 0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(59,35) = 0.7071428571428572*m0rOther[7]*uOther[7]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[7]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[6]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[6]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.2*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.2*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(59,36) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(59,37) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(59,38) = 0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(59,39) = 0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(60,32) = 0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(60,33) = 0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(60,34) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(60,35) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(60,36) = 0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[2]*uOther[6]*mnuOther+0.159719141249985*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(60,37) = 0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(60,38) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(60,39) = 0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(61,32) = 0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(61,33) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(61,34) = 0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(61,35) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(61,36) = 0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(61,37) = 0.5357142857142857*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*uOther[1]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(61,38) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(61,39) = 0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(62,32) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(62,33) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(62,34) = 0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(62,35) = 0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(62,36) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(62,37) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(62,38) = 0.6173469387755102*m0rOther[7]*uOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[7]*mnuOther+0.351382110749967*uOther[1]*m0rOther[7]*mnuOther+0.9642857142857143*m0rOther[6]*uOther[6]*mnuOther+0.2874944542499729*m0rOther[2]*uOther[6]*mnuOther+0.2874944542499729*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(62,39) = 0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(63,32) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(63,33) = 0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(63,34) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(63,35) = 0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(63,36) = 0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(63,37) = 0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(63,38) = 0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(63,39) = 0.9642857142857143*m0rOther[7]*uOther[7]*mnuOther+0.2874944542499729*m0rOther[1]*uOther[7]*mnuOther+0.2874944542499729*uOther[1]*m0rOther[7]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[6]*mnuOther+0.351382110749967*m0rOther[2]*uOther[6]*mnuOther+0.351382110749967*uOther[2]*m0rOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfX-m0Self*m1OtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[0] = 0.5*m0rOther[7]*m1rSelf[7]-0.5*m0rSelf[7]*m1rOther[7]+0.5*m0rOther[6]*m1rSelf[6]-0.5*m0rSelf[6]*m1rOther[6]+0.5*m0rOther[5]*m1rSelf[5]-0.5*m0rSelf[5]*m1rOther[5]+0.5*m0rOther[4]*m1rSelf[4]-0.5*m0rSelf[4]*m1rOther[4]+0.5*m0rOther[3]*m1rSelf[3]-0.5*m0rSelf[3]*m1rOther[3]+0.5*m0rOther[2]*m1rSelf[2]-0.5*m0rSelf[2]*m1rOther[2]+0.5*m0rOther[1]*m1rSelf[1]-0.5*m0rSelf[1]*m1rOther[1]+0.5*m0rOther[0]*m1rSelf[0]-0.5*m0rSelf[0]*m1rOther[0]; 
  m1EffD[1] = 0.5000000000000001*m0rOther[5]*m1rSelf[7]-0.5000000000000001*m0rSelf[5]*m1rOther[7]-0.5000000000000001*m1rOther[5]*m0rSelf[7]+0.5000000000000001*m1rSelf[5]*m0rOther[7]+0.447213595499958*m0rOther[3]*m1rSelf[6]-0.447213595499958*m0rSelf[3]*m1rOther[6]-0.447213595499958*m1rOther[3]*m0rSelf[6]+0.447213595499958*m1rSelf[3]*m0rOther[6]+0.4472135954999579*m0rOther[1]*m1rSelf[4]-0.4472135954999579*m0rSelf[1]*m1rOther[4]-0.4472135954999579*m1rOther[1]*m0rSelf[4]+0.4472135954999579*m1rSelf[1]*m0rOther[4]+0.5*m0rOther[2]*m1rSelf[3]-0.5*m0rSelf[2]*m1rOther[3]-0.5*m1rOther[2]*m0rSelf[3]+0.5*m1rSelf[2]*m0rOther[3]+0.5*m0rOther[0]*m1rSelf[1]-0.5*m0rSelf[0]*m1rOther[1]-0.5*m1rOther[0]*m0rSelf[1]+0.5*m1rSelf[0]*m0rOther[1]; 
  m1EffD[2] = 0.447213595499958*m0rOther[3]*m1rSelf[7]-0.447213595499958*m0rSelf[3]*m1rOther[7]-0.447213595499958*m1rOther[3]*m0rSelf[7]+0.447213595499958*m1rSelf[3]*m0rOther[7]+0.5000000000000001*m0rOther[4]*m1rSelf[6]-0.5000000000000001*m0rSelf[4]*m1rOther[6]-0.5000000000000001*m1rOther[4]*m0rSelf[6]+0.5000000000000001*m1rSelf[4]*m0rOther[6]+0.4472135954999579*m0rOther[2]*m1rSelf[5]-0.4472135954999579*m0rSelf[2]*m1rOther[5]-0.4472135954999579*m1rOther[2]*m0rSelf[5]+0.4472135954999579*m1rSelf[2]*m0rOther[5]+0.5*m0rOther[1]*m1rSelf[3]-0.5*m0rSelf[1]*m1rOther[3]-0.5*m1rOther[1]*m0rSelf[3]+0.5*m1rSelf[1]*m0rOther[3]+0.5*m0rOther[0]*m1rSelf[2]-0.5*m0rSelf[0]*m1rOther[2]-0.5*m1rOther[0]*m0rSelf[2]+0.5*m1rSelf[0]*m0rOther[2]; 
  m1EffD[3] = 0.4*m0rOther[6]*m1rSelf[7]+0.447213595499958*m0rOther[2]*m1rSelf[7]-0.4*m0rSelf[6]*m1rOther[7]-0.447213595499958*m0rSelf[2]*m1rOther[7]-0.4*m1rOther[6]*m0rSelf[7]-0.447213595499958*m1rOther[2]*m0rSelf[7]+0.4*m1rSelf[6]*m0rOther[7]+0.447213595499958*m1rSelf[2]*m0rOther[7]+0.447213595499958*m0rOther[1]*m1rSelf[6]-0.447213595499958*m0rSelf[1]*m1rOther[6]-0.447213595499958*m1rOther[1]*m0rSelf[6]+0.447213595499958*m1rSelf[1]*m0rOther[6]+0.4472135954999579*m0rOther[3]*m1rSelf[5]-0.4472135954999579*m0rSelf[3]*m1rOther[5]-0.4472135954999579*m1rOther[3]*m0rSelf[5]+0.4472135954999579*m1rSelf[3]*m0rOther[5]+0.4472135954999579*m0rOther[3]*m1rSelf[4]-0.4472135954999579*m0rSelf[3]*m1rOther[4]-0.4472135954999579*m1rOther[3]*m0rSelf[4]+0.4472135954999579*m1rSelf[3]*m0rOther[4]+0.5*m0rOther[0]*m1rSelf[3]-0.5*m0rSelf[0]*m1rOther[3]-0.5*m1rOther[0]*m0rSelf[3]+0.5*m1rSelf[0]*m0rOther[3]+0.5*m0rOther[1]*m1rSelf[2]-0.5*m0rSelf[1]*m1rOther[2]-0.5*m1rOther[1]*m0rSelf[2]+0.5*m1rSelf[1]*m0rOther[2]; 
  m1EffD[4] = 0.4472135954999579*m0rOther[7]*m1rSelf[7]-0.4472135954999579*m0rSelf[7]*m1rOther[7]+0.31943828249997*m0rOther[6]*m1rSelf[6]+0.5000000000000001*m0rOther[2]*m1rSelf[6]-0.31943828249997*m0rSelf[6]*m1rOther[6]-0.5000000000000001*m0rSelf[2]*m1rOther[6]-0.5000000000000001*m1rOther[2]*m0rSelf[6]+0.5000000000000001*m1rSelf[2]*m0rOther[6]+0.31943828249997*m0rOther[4]*m1rSelf[4]+0.5*m0rOther[0]*m1rSelf[4]-0.31943828249997*m0rSelf[4]*m1rOther[4]-0.5*m0rSelf[0]*m1rOther[4]-0.5*m1rOther[0]*m0rSelf[4]+0.5*m1rSelf[0]*m0rOther[4]+0.4472135954999579*m0rOther[3]*m1rSelf[3]-0.4472135954999579*m0rSelf[3]*m1rOther[3]+0.4472135954999579*m0rOther[1]*m1rSelf[1]-0.4472135954999579*m0rSelf[1]*m1rOther[1]; 
  m1EffD[5] = 0.31943828249997*m0rOther[7]*m1rSelf[7]+0.5000000000000001*m0rOther[1]*m1rSelf[7]-0.31943828249997*m0rSelf[7]*m1rOther[7]-0.5000000000000001*m0rSelf[1]*m1rOther[7]-0.5000000000000001*m1rOther[1]*m0rSelf[7]+0.5000000000000001*m1rSelf[1]*m0rOther[7]+0.4472135954999579*m0rOther[6]*m1rSelf[6]-0.4472135954999579*m0rSelf[6]*m1rOther[6]+0.31943828249997*m0rOther[5]*m1rSelf[5]+0.5*m0rOther[0]*m1rSelf[5]-0.31943828249997*m0rSelf[5]*m1rOther[5]-0.5*m0rSelf[0]*m1rOther[5]-0.5*m1rOther[0]*m0rSelf[5]+0.5*m1rSelf[0]*m0rOther[5]+0.4472135954999579*m0rOther[3]*m1rSelf[3]-0.4472135954999579*m0rSelf[3]*m1rOther[3]+0.4472135954999579*m0rOther[2]*m1rSelf[2]-0.4472135954999579*m0rSelf[2]*m1rOther[2]; 
  m1EffD[6] = 0.4*m0rOther[3]*m1rSelf[7]-0.4*m0rSelf[3]*m1rOther[7]-0.4*m1rOther[3]*m0rSelf[7]+0.4*m1rSelf[3]*m0rOther[7]+0.4472135954999579*m0rOther[5]*m1rSelf[6]+0.31943828249997*m0rOther[4]*m1rSelf[6]+0.5*m0rOther[0]*m1rSelf[6]-0.4472135954999579*m0rSelf[5]*m1rOther[6]-0.31943828249997*m0rSelf[4]*m1rOther[6]-0.5*m0rSelf[0]*m1rOther[6]-0.4472135954999579*m1rOther[5]*m0rSelf[6]-0.31943828249997*m1rOther[4]*m0rSelf[6]-0.5*m1rOther[0]*m0rSelf[6]+0.4472135954999579*m1rSelf[5]*m0rOther[6]+0.31943828249997*m1rSelf[4]*m0rOther[6]+0.5*m1rSelf[0]*m0rOther[6]+0.5000000000000001*m0rOther[2]*m1rSelf[4]-0.5000000000000001*m0rSelf[2]*m1rOther[4]-0.5000000000000001*m1rOther[2]*m0rSelf[4]+0.5000000000000001*m1rSelf[2]*m0rOther[4]+0.447213595499958*m0rOther[1]*m1rSelf[3]-0.447213595499958*m0rSelf[1]*m1rOther[3]-0.447213595499958*m1rOther[1]*m0rSelf[3]+0.447213595499958*m1rSelf[1]*m0rOther[3]; 
  m1EffD[7] = 0.31943828249997*m0rOther[5]*m1rSelf[7]+0.4472135954999579*m0rOther[4]*m1rSelf[7]+0.5*m0rOther[0]*m1rSelf[7]-0.31943828249997*m0rSelf[5]*m1rOther[7]-0.4472135954999579*m0rSelf[4]*m1rOther[7]-0.5*m0rSelf[0]*m1rOther[7]-0.31943828249997*m1rOther[5]*m0rSelf[7]-0.4472135954999579*m1rOther[4]*m0rSelf[7]-0.5*m1rOther[0]*m0rSelf[7]+0.31943828249997*m1rSelf[5]*m0rOther[7]+0.4472135954999579*m1rSelf[4]*m0rOther[7]+0.5*m1rSelf[0]*m0rOther[7]+0.4*m0rOther[3]*m1rSelf[6]-0.4*m0rSelf[3]*m1rOther[6]-0.4*m1rOther[3]*m0rSelf[6]+0.4*m1rSelf[3]*m0rOther[6]+0.5000000000000001*m0rOther[1]*m1rSelf[5]-0.5000000000000001*m0rSelf[1]*m1rOther[5]-0.5000000000000001*m1rOther[1]*m0rSelf[5]+0.5000000000000001*m1rSelf[1]*m0rOther[5]+0.447213595499958*m0rOther[2]*m1rSelf[3]-0.447213595499958*m0rSelf[2]*m1rOther[3]-0.447213595499958*m1rOther[2]*m0rSelf[3]+0.447213595499958*m1rSelf[2]*m0rOther[3]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(8,8); 
  dataDiv->AEM_S(0,0) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(0,2) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(0,3) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(0,4) = 0.5*m0rSelf[4]*mnuSelf+0.5*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(0,5) = 0.5*m0rSelf[5]*mnuSelf+0.5*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(0,6) = 0.5*m0rSelf[6]*mnuSelf+0.5*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(0,7) = 0.5*m0rSelf[7]*mnuSelf+0.5*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(1,2) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(1,3) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf+0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(1,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf+0.4472135954999579*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,5) = 0.5000000000000001*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(1,6) = 0.447213595499958*m0rSelf[3]*mnuSelf+0.447213595499958*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(1,7) = 0.5000000000000001*m0rSelf[5]*mnuSelf+0.5000000000000001*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,1) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(2,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(2,3) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf+0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(2,4) = 0.5000000000000001*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(2,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf+0.4472135954999579*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,6) = 0.5000000000000001*m0rSelf[4]*mnuSelf+0.5000000000000001*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(2,7) = 0.447213595499958*m0rSelf[3]*mnuSelf+0.447213595499958*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,0) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,1) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf+0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(3,2) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf+0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(3,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(3,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,6) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(3,7) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(4,0) = 0.5*m0rSelf[4]*mnuSelf+0.5*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(4,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf+0.4472135954999579*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(4,2) = 0.5000000000000001*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(4,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(4,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(4,6) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf+0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(4,7) = 0.4472135954999579*m0rSelf[7]*mnuSelf+0.4472135954999579*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(5,0) = 0.5*m0rSelf[5]*mnuSelf+0.5*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(5,1) = 0.5000000000000001*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(5,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf+0.4472135954999579*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(5,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(5,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(5,6) = 0.4472135954999579*m0rSelf[6]*mnuSelf+0.4472135954999579*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(5,7) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf+0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(6,0) = 0.5*m0rSelf[6]*mnuSelf+0.5*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(6,1) = 0.447213595499958*m0rSelf[3]*mnuSelf+0.447213595499958*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(6,2) = 0.5000000000000001*m0rSelf[4]*mnuSelf+0.5000000000000001*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(6,3) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf+0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(6,4) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf+0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(6,5) = 0.4472135954999579*m0rSelf[6]*mnuSelf+0.4472135954999579*m0rOther[6]*mnuOther; 
  dataDiv->AEM_S(6,6) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(6,7) = 0.4*m0rSelf[3]*mnuSelf+0.4*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(7,0) = 0.5*m0rSelf[7]*mnuSelf+0.5*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(7,1) = 0.5000000000000001*m0rSelf[5]*mnuSelf+0.5000000000000001*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(7,2) = 0.447213595499958*m0rSelf[3]*mnuSelf+0.447213595499958*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(7,3) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf+0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(7,4) = 0.4472135954999579*m0rSelf[7]*mnuSelf+0.4472135954999579*m0rOther[7]*mnuOther; 
  dataDiv->AEM_S(7,5) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf+0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(7,6) = 0.4*m0rSelf[3]*mnuSelf+0.4*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(7,7) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[0],m1EffD[1],m1EffD[2],m1EffD[3],m1EffD[4],m1EffD[5],m1EffD[6],m1EffD[7]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+0,8,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (-2.0*m1EffD[0]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (-2.0*m1EffD[1]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (-2.0*m1EffD[2]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += (-2.0*m1EffD[3]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (-2.0*m1EffD[4]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (-2.0*m1EffD[5]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
  m1Relax[6] += (-2.0*m1EffD[6]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (-2.0*m1EffD[7]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(40,8) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(40,9) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(40,10) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(40,11) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(40,12) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(40,13) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(40,14) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(40,15) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(41,8) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(41,9) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(41,10) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(41,11) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(41,12) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(41,13) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(41,14) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(41,15) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(42,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(42,9) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(42,10) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(42,11) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(42,12) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(42,13) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(42,14) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(42,15) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,8) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,9) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(43,10) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(43,11) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(43,12) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,13) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(43,14) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(43,15) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(44,8) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(44,9) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(44,10) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(44,11) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(44,12) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(44,14) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(44,15) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(45,8) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(45,9) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(45,10) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(45,11) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(45,13) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(45,14) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(45,15) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(46,8) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(46,9) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(46,10) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(46,11) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(46,12) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(46,13) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(46,14) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(46,15) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(47,8) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(47,9) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(47,10) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(47,11) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(47,12) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(47,13) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(47,14) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(47,15) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(40,24) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(40,25) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(40,26) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(40,27) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(40,28) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(40,29) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(40,30) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(40,31) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(41,24) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(41,25) = (-0.4472135954999579*cMSelf[12]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(41,26) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(41,27) = (-0.447213595499958*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(41,28) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(41,29) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(41,30) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(41,31) = -0.5000000000000001*cMSelf[13]*mnuSelf; 
  data->AEM_S(42,24) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(42,25) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(42,26) = (-0.4472135954999579*cMSelf[13]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(42,27) = (-0.447213595499958*cMSelf[15]*mnuSelf)-0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(42,28) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(42,29) = -0.4472135954999579*cMSelf[10]*mnuSelf; 
  data->AEM_S(42,30) = -0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(42,31) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(43,24) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(43,25) = (-0.447213595499958*cMSelf[14]*mnuSelf)-0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(43,26) = (-0.447213595499958*cMSelf[15]*mnuSelf)-0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(43,27) = (-0.4472135954999579*cMSelf[13]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(43,28) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(43,29) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(43,30) = (-0.4*cMSelf[15]*mnuSelf)-0.447213595499958*cMSelf[9]*mnuSelf; 
  data->AEM_S(43,31) = (-0.4*cMSelf[14]*mnuSelf)-0.447213595499958*cMSelf[10]*mnuSelf; 
  data->AEM_S(44,24) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(44,25) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(44,26) = -0.5000000000000001*cMSelf[14]*mnuSelf; 
  data->AEM_S(44,27) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(44,28) = (-0.31943828249997*cMSelf[12]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(44,30) = (-0.31943828249997*cMSelf[14]*mnuSelf)-0.5000000000000001*cMSelf[10]*mnuSelf; 
  data->AEM_S(44,31) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(45,24) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(45,25) = -0.5000000000000001*cMSelf[15]*mnuSelf; 
  data->AEM_S(45,26) = -0.4472135954999579*cMSelf[10]*mnuSelf; 
  data->AEM_S(45,27) = -0.4472135954999579*cMSelf[11]*mnuSelf; 
  data->AEM_S(45,29) = (-0.31943828249997*cMSelf[13]*mnuSelf)-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(45,30) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(45,31) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.5000000000000001*cMSelf[9]*mnuSelf; 
  data->AEM_S(46,24) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(46,25) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(46,26) = -0.5000000000000001*cMSelf[12]*mnuSelf; 
  data->AEM_S(46,27) = (-0.4*cMSelf[15]*mnuSelf)-0.447213595499958*cMSelf[9]*mnuSelf; 
  data->AEM_S(46,28) = (-0.31943828249997*cMSelf[14]*mnuSelf)-0.5000000000000001*cMSelf[10]*mnuSelf; 
  data->AEM_S(46,29) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(46,30) = (-0.4472135954999579*cMSelf[13]*mnuSelf)-0.31943828249997*cMSelf[12]*mnuSelf-0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(46,31) = -0.4*cMSelf[11]*mnuSelf; 
  data->AEM_S(47,24) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(47,25) = -0.5000000000000001*cMSelf[13]*mnuSelf; 
  data->AEM_S(47,26) = -0.447213595499958*cMSelf[11]*mnuSelf; 
  data->AEM_S(47,27) = (-0.4*cMSelf[14]*mnuSelf)-0.447213595499958*cMSelf[10]*mnuSelf; 
  data->AEM_S(47,28) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(47,29) = (-0.31943828249997*cMSelf[15]*mnuSelf)-0.5000000000000001*cMSelf[9]*mnuSelf; 
  data->AEM_S(47,30) = -0.4*cMSelf[11]*mnuSelf; 
  data->AEM_S(47,31) = (-0.31943828249997*cMSelf[13]*mnuSelf)-0.4472135954999579*cMSelf[12]*mnuSelf-0.5*cMSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(40,40) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(40,41) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(40,42) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(40,43) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(40,44) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(40,45) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(40,46) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(40,47) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(41,40) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(41,41) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(41,42) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(41,43) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(41,44) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(41,45) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(41,46) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(41,47) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(42,40) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(42,41) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(42,42) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(42,43) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(42,44) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(42,45) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(42,46) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(42,47) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(43,40) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(43,41) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(43,42) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(43,43) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(43,44) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(43,45) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(43,46) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(43,47) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(44,40) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(44,41) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(44,42) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(44,43) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(44,44) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(44,46) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(44,47) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(45,40) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(45,41) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(45,42) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(45,43) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(45,45) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(45,46) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(45,47) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(46,40) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(46,41) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(46,42) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(46,43) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(46,44) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(46,45) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(46,46) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.31943828249997*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(46,47) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(47,40) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(47,41) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(47,42) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(47,43) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(47,44) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(47,45) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(47,46) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(47,47) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(40,56) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(40,57) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(40,58) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(40,59) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(40,60) = 0.5*cMOther[12]*mnuOther; 
  data->AEM_S(40,61) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(40,62) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(40,63) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(41,56) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(41,57) = 0.4472135954999579*cMOther[12]*mnuOther+0.5*cMOther[8]*mnuOther; 
  data->AEM_S(41,58) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(41,59) = 0.447213595499958*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(41,60) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(41,61) = 0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(41,62) = 0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(41,63) = 0.5000000000000001*cMOther[13]*mnuOther; 
  data->AEM_S(42,56) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(42,57) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(42,58) = 0.4472135954999579*cMOther[13]*mnuOther+0.5*cMOther[8]*mnuOther; 
  data->AEM_S(42,59) = 0.447213595499958*cMOther[15]*mnuOther+0.5*cMOther[9]*mnuOther; 
  data->AEM_S(42,60) = 0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(42,61) = 0.4472135954999579*cMOther[10]*mnuOther; 
  data->AEM_S(42,62) = 0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(42,63) = 0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(43,56) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(43,57) = 0.447213595499958*cMOther[14]*mnuOther+0.5*cMOther[10]*mnuOther; 
  data->AEM_S(43,58) = 0.447213595499958*cMOther[15]*mnuOther+0.5*cMOther[9]*mnuOther; 
  data->AEM_S(43,59) = 0.4472135954999579*cMOther[13]*mnuOther+0.4472135954999579*cMOther[12]*mnuOther+0.5*cMOther[8]*mnuOther; 
  data->AEM_S(43,60) = 0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(43,61) = 0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(43,62) = 0.4*cMOther[15]*mnuOther+0.447213595499958*cMOther[9]*mnuOther; 
  data->AEM_S(43,63) = 0.4*cMOther[14]*mnuOther+0.447213595499958*cMOther[10]*mnuOther; 
  data->AEM_S(44,56) = 0.5*cMOther[12]*mnuOther; 
  data->AEM_S(44,57) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(44,58) = 0.5000000000000001*cMOther[14]*mnuOther; 
  data->AEM_S(44,59) = 0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(44,60) = 0.31943828249997*cMOther[12]*mnuOther+0.5*cMOther[8]*mnuOther; 
  data->AEM_S(44,62) = 0.31943828249997*cMOther[14]*mnuOther+0.5000000000000001*cMOther[10]*mnuOther; 
  data->AEM_S(44,63) = 0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(45,56) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(45,57) = 0.5000000000000001*cMOther[15]*mnuOther; 
  data->AEM_S(45,58) = 0.4472135954999579*cMOther[10]*mnuOther; 
  data->AEM_S(45,59) = 0.4472135954999579*cMOther[11]*mnuOther; 
  data->AEM_S(45,61) = 0.31943828249997*cMOther[13]*mnuOther+0.5*cMOther[8]*mnuOther; 
  data->AEM_S(45,62) = 0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(45,63) = 0.31943828249997*cMOther[15]*mnuOther+0.5000000000000001*cMOther[9]*mnuOther; 
  data->AEM_S(46,56) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(46,57) = 0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(46,58) = 0.5000000000000001*cMOther[12]*mnuOther; 
  data->AEM_S(46,59) = 0.4*cMOther[15]*mnuOther+0.447213595499958*cMOther[9]*mnuOther; 
  data->AEM_S(46,60) = 0.31943828249997*cMOther[14]*mnuOther+0.5000000000000001*cMOther[10]*mnuOther; 
  data->AEM_S(46,61) = 0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(46,62) = 0.4472135954999579*cMOther[13]*mnuOther+0.31943828249997*cMOther[12]*mnuOther+0.5*cMOther[8]*mnuOther; 
  data->AEM_S(46,63) = 0.4*cMOther[11]*mnuOther; 
  data->AEM_S(47,56) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(47,57) = 0.5000000000000001*cMOther[13]*mnuOther; 
  data->AEM_S(47,58) = 0.447213595499958*cMOther[11]*mnuOther; 
  data->AEM_S(47,59) = 0.4*cMOther[14]*mnuOther+0.447213595499958*cMOther[10]*mnuOther; 
  data->AEM_S(47,60) = 0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(47,61) = 0.31943828249997*cMOther[15]*mnuOther+0.5000000000000001*cMOther[9]*mnuOther; 
  data->AEM_S(47,62) = 0.4*cMOther[11]*mnuOther; 
  data->AEM_S(47,63) = 0.31943828249997*cMOther[13]*mnuOther+0.4472135954999579*cMOther[12]*mnuOther+0.5*cMOther[8]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(56,8) = (-0.25*m0rSelf[7]*uSelf[15]*mnuSelf)-0.25*m0rSelf[6]*uSelf[14]*mnuSelf-0.25*m0rSelf[5]*uSelf[13]*mnuSelf-0.25*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[3]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(56,9) = (-0.2500000000000001*m0rSelf[5]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[11]*mnuSelf-0.25*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(56,10) = (-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf)-0.2500000000000001*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.25*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(56,11) = (-0.2*m0rSelf[6]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.2*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(56,12) = (-0.223606797749979*m0rSelf[7]*uSelf[15]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[4]*uSelf[8]*mnuSelf; 
  data->AEM_S(56,13) = (-0.159719141249985*m0rSelf[7]*uSelf[15]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[9]*mnuSelf-0.25*m0rSelf[5]*uSelf[8]*mnuSelf; 
  data->AEM_S(56,14) = (-0.2*m0rSelf[3]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[6]*uSelf[8]*mnuSelf; 
  data->AEM_S(56,15) = (-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[9]*mnuSelf-0.25*m0rSelf[7]*uSelf[8]*mnuSelf; 
  data->AEM_S(57,8) = (-0.2500000000000001*m0rSelf[5]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[11]*mnuSelf-0.25*m0rSelf[2]*uSelf[11]*mnuSelf-0.25*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(57,9) = (-0.45*m0rSelf[7]*uSelf[15]*mnuSelf)-0.3928571428571428*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[5]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[12]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf-0.45*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(57,10) = (-0.2*m0rSelf[6]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.2*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(57,11) = (-0.4024922359499621*m0rSelf[3]*uSelf[15]*mnuSelf)-0.2*m0rSelf[5]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.447213595499958*m1rSelf[14]*mnuSelf-0.2*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[11]*mnuSelf-0.45*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(57,12) = (-0.223606797749979*m0rSelf[5]*uSelf[15]*mnuSelf)-0.3928571428571429*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(57,13) = (-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[15]*mnuSelf+0.5000000000000001*m1rSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.25*m0rSelf[5]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[8]*mnuSelf; 
  data->AEM_S(57,14) = (-0.351382110749967*m0rSelf[6]*uSelf[15]*mnuSelf)-0.2*m0rSelf[2]*uSelf[15]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[14]*mnuSelf-0.2*m0rSelf[3]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[12]*mnuSelf-0.2*m0rSelf[5]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf-0.2*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(57,15) = (-0.2874944542499729*m0rSelf[7]*uSelf[15]*mnuSelf)-0.45*m0rSelf[1]*uSelf[15]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[14]*mnuSelf-0.2*m0rSelf[2]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[13]*mnuSelf+0.5000000000000001*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[11]*mnuSelf-0.2*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[7]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[8]*mnuSelf; 
  data->AEM_S(58,8) = (-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf)-0.2500000000000001*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[11]*mnuSelf-0.25*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.25*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(58,9) = (-0.2*m0rSelf[6]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.2*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(58,10) = (-0.3928571428571428*m0rSelf[7]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.45*m0rSelf[6]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.25*m0rSelf[4]*uSelf[12]*mnuSelf-0.45*m0rSelf[3]*uSelf[11]*mnuSelf-0.45*m0rSelf[2]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(58,11) = (-0.3928571428571429*m0rSelf[5]*uSelf[15]*mnuSelf)-0.2*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.447213595499958*m1rSelf[15]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.2*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[11]*mnuSelf-0.45*m0rSelf[2]*uSelf[11]*mnuSelf-0.45*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(58,12) = (-0.2*m0rSelf[3]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[14]*mnuSelf+0.5000000000000001*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[12]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[8]*mnuSelf; 
  data->AEM_S(58,13) = (-0.3928571428571429*m0rSelf[3]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(58,14) = (-0.351382110749967*m0rSelf[7]*uSelf[15]*mnuSelf)-0.2*m0rSelf[1]*uSelf[15]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[14]*mnuSelf-0.45*m0rSelf[2]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[12]*mnuSelf+0.5000000000000001*m1rSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[11]*mnuSelf-0.45*m0rSelf[6]*uSelf[10]*mnuSelf-0.2*m0rSelf[7]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[8]*mnuSelf; 
  data->AEM_S(58,15) = (-0.351382110749967*m0rSelf[6]*uSelf[15]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[15]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[14]*mnuSelf-0.2*m0rSelf[1]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[13]*mnuSelf-0.2*m0rSelf[3]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[11]*mnuSelf-0.2*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.2*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(59,8) = (-0.2*m0rSelf[6]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.2*m0rSelf[7]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[11]*mnuSelf-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[10]*mnuSelf-0.25*m0rSelf[1]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(59,9) = (-0.4024922359499621*m0rSelf[3]*uSelf[15]*mnuSelf)-0.2*m0rSelf[5]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.447213595499958*m1rSelf[14]*mnuSelf-0.2*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[11]*mnuSelf-0.45*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(59,10) = (-0.3928571428571429*m0rSelf[5]*uSelf[15]*mnuSelf)-0.2*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.447213595499958*m1rSelf[15]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.2*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[11]*mnuSelf-0.45*m0rSelf[2]*uSelf[11]*mnuSelf-0.45*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(59,11) = (-0.7071428571428572*m0rSelf[7]*uSelf[15]*mnuSelf)-0.4024922359499621*m0rSelf[1]*uSelf[15]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[14]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[13]*mnuSelf-0.2*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.2*m0rSelf[5]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[12]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf-0.81*m0rSelf[3]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[10]*mnuSelf-0.45*m0rSelf[2]*uSelf[10]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[9]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(59,12) = (-0.3513821107499669*m0rSelf[6]*uSelf[15]*mnuSelf)-0.2*m0rSelf[2]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[14]*mnuSelf-0.2*m0rSelf[3]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[12]*mnuSelf-0.2*m0rSelf[5]*uSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.2*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(59,13) = (-0.3513821107499669*m0rSelf[6]*uSelf[15]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[14]*mnuSelf-0.2*m0rSelf[1]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[13]*mnuSelf-0.2*m0rSelf[3]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[11]*mnuSelf-0.2*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.2*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(59,14) = (-0.3513821107499669*m0rSelf[5]*uSelf[15]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[15]*mnuSelf-0.2*m0rSelf[0]*uSelf[15]*mnuSelf+0.4*m1rSelf[15]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[13]*mnuSelf-0.2*m0rSelf[1]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[12]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[10]*mnuSelf-0.2*m0rSelf[5]*uSelf[9]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.447213595499958*m1rSelf[9]*mnuSelf-0.2*m0rSelf[7]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(59,15) = (-0.7071428571428572*m0rSelf[3]*uSelf[15]*mnuSelf)-0.3513821107499669*m0rSelf[5]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[14]*mnuSelf-0.2*m0rSelf[0]*uSelf[14]*mnuSelf+0.4*m1rSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[12]*mnuSelf-0.2*m0rSelf[2]*uSelf[12]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[10]*mnuSelf-0.2*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.447213595499958*m1rSelf[10]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[9]*mnuSelf-0.2*m0rSelf[6]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(60,8) = (-0.223606797749979*m0rSelf[7]*uSelf[15]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[4]*uSelf[8]*mnuSelf; 
  data->AEM_S(60,9) = (-0.223606797749979*m0rSelf[5]*uSelf[15]*mnuSelf)-0.3928571428571429*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(60,10) = (-0.2*m0rSelf[3]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[14]*mnuSelf+0.5000000000000001*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[12]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[8]*mnuSelf; 
  data->AEM_S(60,11) = (-0.3513821107499669*m0rSelf[6]*uSelf[15]*mnuSelf)-0.2*m0rSelf[2]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[14]*mnuSelf-0.2*m0rSelf[3]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[12]*mnuSelf-0.2*m0rSelf[5]*uSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.2*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(60,12) = (-0.3928571428571428*m0rSelf[7]*uSelf[15]*mnuSelf)-0.5357142857142857*m0rSelf[6]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[5]*uSelf[13]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[12]*mnuSelf+0.31943828249997*m1rSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[11]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[10]*mnuSelf-0.25*m0rSelf[2]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[9]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(60,13) = (-0.1428571428571428*m0rSelf[7]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[5]*uSelf[12]*mnuSelf-0.2*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[9]*mnuSelf; 
  data->AEM_S(60,14) = (-0.3513821107499669*m0rSelf[3]*uSelf[15]*mnuSelf)-0.1428571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[14]*mnuSelf+0.31943828249997*m1rSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[10]*mnuSelf+0.5000000000000001*m1rSelf[10]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[9]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(60,15) = (-0.1428571428571428*m0rSelf[5]*uSelf[15]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[11]*mnuSelf-0.2*m0rSelf[2]*uSelf[11]*mnuSelf-0.2*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[8]*mnuSelf; 
  data->AEM_S(61,8) = (-0.159719141249985*m0rSelf[7]*uSelf[15]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[9]*mnuSelf-0.25*m0rSelf[5]*uSelf[8]*mnuSelf; 
  data->AEM_S(61,9) = (-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[15]*mnuSelf+0.5000000000000001*m1rSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.25*m0rSelf[5]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[8]*mnuSelf; 
  data->AEM_S(61,10) = (-0.3928571428571429*m0rSelf[3]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(61,11) = (-0.3513821107499669*m0rSelf[6]*uSelf[15]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[14]*mnuSelf-0.2*m0rSelf[1]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[13]*mnuSelf-0.2*m0rSelf[3]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[11]*mnuSelf-0.2*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.2*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(61,12) = (-0.1428571428571428*m0rSelf[7]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[5]*uSelf[12]*mnuSelf-0.2*m0rSelf[3]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[9]*mnuSelf; 
  data->AEM_S(61,13) = (-0.5357142857142857*m0rSelf[7]*uSelf[15]*mnuSelf)-0.159719141249985*m0rSelf[1]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[14]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[13]*mnuSelf+0.31943828249997*m1rSelf[13]*mnuSelf-0.25*m0rSelf[4]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(61,14) = (-0.3513821107499669*m0rSelf[3]*uSelf[15]*mnuSelf)-0.3928571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[11]*mnuSelf-0.2*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[10]*mnuSelf-0.2*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[8]*mnuSelf; 
  data->AEM_S(61,15) = (-0.5357142857142857*m0rSelf[5]*uSelf[15]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[15]*mnuSelf+0.31943828249997*m1rSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[14]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[9]*mnuSelf+0.5000000000000001*m1rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(62,8) = (-0.2*m0rSelf[3]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[12]*mnuSelf-0.2*m0rSelf[7]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[6]*uSelf[8]*mnuSelf; 
  data->AEM_S(62,9) = (-0.351382110749967*m0rSelf[6]*uSelf[15]*mnuSelf)-0.2*m0rSelf[2]*uSelf[15]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[14]*mnuSelf-0.2*m0rSelf[3]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[12]*mnuSelf-0.2*m0rSelf[5]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf-0.2*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(62,10) = (-0.351382110749967*m0rSelf[7]*uSelf[15]*mnuSelf)-0.2*m0rSelf[1]*uSelf[15]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[14]*mnuSelf-0.45*m0rSelf[2]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[12]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[12]*mnuSelf+0.5000000000000001*m1rSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[11]*mnuSelf-0.45*m0rSelf[6]*uSelf[10]*mnuSelf-0.2*m0rSelf[7]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[8]*mnuSelf; 
  data->AEM_S(62,11) = (-0.3513821107499669*m0rSelf[5]*uSelf[15]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[15]*mnuSelf-0.2*m0rSelf[0]*uSelf[15]*mnuSelf+0.4*m1rSelf[15]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[13]*mnuSelf-0.2*m0rSelf[1]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[12]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[10]*mnuSelf-0.2*m0rSelf[5]*uSelf[9]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.447213595499958*m1rSelf[9]*mnuSelf-0.2*m0rSelf[7]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(62,12) = (-0.3513821107499669*m0rSelf[3]*uSelf[15]*mnuSelf)-0.1428571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[14]*mnuSelf+0.31943828249997*m1rSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[10]*mnuSelf+0.5000000000000001*m1rSelf[10]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[9]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(62,13) = (-0.3513821107499669*m0rSelf[3]*uSelf[15]*mnuSelf)-0.3928571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[11]*mnuSelf-0.2*m0rSelf[1]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[10]*mnuSelf-0.2*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[8]*mnuSelf; 
  data->AEM_S(62,14) = (-0.6173469387755102*m0rSelf[7]*uSelf[15]*mnuSelf)-0.351382110749967*m0rSelf[1]*uSelf[15]*mnuSelf-0.9642857142857143*m0rSelf[6]*uSelf[14]*mnuSelf-0.2874944542499729*m0rSelf[2]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[12]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[12]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[12]*mnuSelf+0.31943828249997*m1rSelf[12]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[11]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[10]*mnuSelf-0.45*m0rSelf[2]*uSelf[10]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(62,15) = (-0.6173469387755102*m0rSelf[6]*uSelf[15]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[15]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[14]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[11]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[11]*mnuSelf-0.2*m0rSelf[0]*uSelf[11]*mnuSelf+0.4*m1rSelf[11]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[10]*mnuSelf-0.2*m0rSelf[1]*uSelf[10]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[9]*mnuSelf-0.2*m0rSelf[2]*uSelf[9]*mnuSelf-0.2*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,8) = (-0.159719141249985*m0rSelf[5]*uSelf[15]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.2*m0rSelf[3]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[12]*mnuSelf-0.2*m0rSelf[6]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[9]*mnuSelf-0.25*m0rSelf[7]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,9) = (-0.2874944542499729*m0rSelf[7]*uSelf[15]*mnuSelf)-0.45*m0rSelf[1]*uSelf[15]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[14]*mnuSelf-0.2*m0rSelf[2]*uSelf[14]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[13]*mnuSelf+0.5000000000000001*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[11]*mnuSelf-0.2*m0rSelf[6]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[7]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,10) = (-0.351382110749967*m0rSelf[6]*uSelf[15]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[15]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[14]*mnuSelf-0.2*m0rSelf[1]*uSelf[14]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[13]*mnuSelf-0.2*m0rSelf[3]*uSelf[12]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[11]*mnuSelf-0.2*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.447213595499958*m1rSelf[11]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.2*m0rSelf[6]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,11) = (-0.7071428571428572*m0rSelf[3]*uSelf[15]*mnuSelf)-0.3513821107499669*m0rSelf[5]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[14]*mnuSelf-0.2*m0rSelf[0]*uSelf[14]*mnuSelf+0.4*m1rSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[13]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[12]*mnuSelf-0.2*m0rSelf[2]*uSelf[12]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[11]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[10]*mnuSelf-0.2*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.447213595499958*m1rSelf[10]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[9]*mnuSelf-0.2*m0rSelf[6]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,12) = (-0.1428571428571428*m0rSelf[5]*uSelf[15]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[14]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[11]*mnuSelf-0.2*m0rSelf[2]*uSelf[11]*mnuSelf-0.2*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,13) = (-0.5357142857142857*m0rSelf[5]*uSelf[15]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[15]*mnuSelf+0.31943828249997*m1rSelf[15]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[14]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[11]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[9]*mnuSelf+0.5000000000000001*m1rSelf[9]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[8]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,14) = (-0.6173469387755102*m0rSelf[6]*uSelf[15]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[15]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[14]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[14]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[13]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[12]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[11]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[11]*mnuSelf-0.2*m0rSelf[0]*uSelf[11]*mnuSelf+0.4*m1rSelf[11]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[10]*mnuSelf-0.2*m0rSelf[1]*uSelf[10]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[9]*mnuSelf-0.2*m0rSelf[2]*uSelf[9]*mnuSelf-0.2*m0rSelf[3]*uSelf[8]*mnuSelf; 
  data->AEM_S(63,15) = (-0.9642857142857143*m0rSelf[7]*uSelf[15]*mnuSelf)-0.2874944542499729*m0rSelf[1]*uSelf[15]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[14]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[14]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[13]*mnuSelf+0.31943828249997*m1rSelf[13]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[12]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[12]*mnuSelf+0.4472135954999579*m1rSelf[12]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[11]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[10]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[9]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(56,40) = 0.25*m0rOther[7]*uOther[15]*mnuOther+0.25*m0rOther[6]*uOther[14]*mnuOther+0.25*m0rOther[5]*uOther[13]*mnuOther+0.25*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[3]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(56,41) = 0.2500000000000001*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.223606797749979*m0rOther[6]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(56,42) = 0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[7]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.25*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(56,43) = 0.2*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.2*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[6]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(56,44) = 0.223606797749979*m0rOther[7]*uOther[15]*mnuOther+0.159719141249985*m0rOther[6]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[4]*uOther[8]*mnuOther; 
  data->AEM_S(56,45) = 0.159719141249985*m0rOther[7]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[6]*uOther[14]*mnuOther+0.159719141249985*m0rOther[5]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[5]*uOther[8]*mnuOther; 
  data->AEM_S(56,46) = 0.2*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[6]*uOther[13]*mnuOther+0.159719141249985*m0rOther[6]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[6]*uOther[8]*mnuOther; 
  data->AEM_S(56,47) = 0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.159719141249985*m0rOther[7]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[9]*mnuOther+0.25*m0rOther[7]*uOther[8]*mnuOther; 
  data->AEM_S(57,40) = 0.2500000000000001*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.223606797749979*m0rOther[6]*uOther[11]*mnuOther+0.25*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(57,41) = 0.45*m0rOther[7]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[5]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[0]*uOther[12]*mnuOther-0.4472135954999579*m1rOther[12]*mnuOther+0.45*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[6]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(57,42) = 0.2*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.2*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[6]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(57,43) = 0.4024922359499621*m0rOther[3]*uOther[15]*mnuOther+0.2*m0rOther[5]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.447213595499958*m1rOther[14]*mnuOther+0.2*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[11]*mnuOther+0.45*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.223606797749979*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[6]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(57,44) = 0.223606797749979*m0rOther[5]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[7]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(57,45) = 0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[15]*mnuOther-0.5000000000000001*m1rOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.159719141249985*m0rOther[7]*uOther[13]*mnuOther+0.25*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.25*m0rOther[5]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[8]*mnuOther; 
  data->AEM_S(57,46) = 0.351382110749967*m0rOther[6]*uOther[15]*mnuOther+0.2*m0rOther[2]*uOther[15]*mnuOther+0.351382110749967*m0rOther[7]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[14]*mnuOther+0.2*m0rOther[3]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[12]*mnuOther+0.2*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.447213595499958*m1rOther[11]*mnuOther+0.2*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(57,47) = 0.2874944542499729*m0rOther[7]*uOther[15]*mnuOther+0.45*m0rOther[1]*uOther[15]*mnuOther+0.351382110749967*m0rOther[6]*uOther[14]*mnuOther+0.2*m0rOther[2]*uOther[14]*mnuOther+0.159719141249985*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[13]*mnuOther-0.5000000000000001*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[11]*mnuOther+0.2*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[7]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[8]*mnuOther; 
  data->AEM_S(58,40) = 0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[7]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.25*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(58,41) = 0.2*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.2*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[6]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(58,42) = 0.3928571428571428*m0rOther[7]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.45*m0rOther[6]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.25*m0rOther[4]*uOther[12]*mnuOther+0.45*m0rOther[3]*uOther[11]*mnuOther+0.45*m0rOther[2]*uOther[10]*mnuOther+0.223606797749979*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(58,43) = 0.3928571428571429*m0rOther[5]*uOther[15]*mnuOther+0.2*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.447213595499958*m1rOther[15]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.2*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[11]*mnuOther+0.45*m0rOther[2]*uOther[11]*mnuOther+0.45*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[7]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(58,44) = 0.2*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[14]*mnuOther-0.5000000000000001*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[6]*uOther[13]*mnuOther+0.159719141249985*m0rOther[6]*uOther[12]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[8]*mnuOther; 
  data->AEM_S(58,45) = 0.3928571428571429*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[6]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.4472135954999579*m1rOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(58,46) = 0.351382110749967*m0rOther[7]*uOther[15]*mnuOther+0.2*m0rOther[1]*uOther[15]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[14]*mnuOther+0.45*m0rOther[2]*uOther[14]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.159719141249985*m0rOther[4]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[12]*mnuOther-0.5000000000000001*m1rOther[12]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[11]*mnuOther+0.45*m0rOther[6]*uOther[10]*mnuOther+0.2*m0rOther[7]*uOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[8]*mnuOther; 
  data->AEM_S(58,47) = 0.351382110749967*m0rOther[6]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[15]*mnuOther+0.351382110749967*m0rOther[7]*uOther[14]*mnuOther+0.2*m0rOther[1]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[13]*mnuOther+0.2*m0rOther[3]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[11]*mnuOther+0.2*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.447213595499958*m1rOther[11]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.2*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(59,40) = 0.2*m0rOther[6]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.2*m0rOther[7]*uOther[14]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther+0.223606797749979*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[7]*uOther[10]*mnuOther+0.25*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[6]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(59,41) = 0.4024922359499621*m0rOther[3]*uOther[15]*mnuOther+0.2*m0rOther[5]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.447213595499958*m1rOther[14]*mnuOther+0.2*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[11]*mnuOther+0.45*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.223606797749979*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[6]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(59,42) = 0.3928571428571429*m0rOther[5]*uOther[15]*mnuOther+0.2*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.447213595499958*m1rOther[15]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.2*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[11]*mnuOther+0.45*m0rOther[2]*uOther[11]*mnuOther+0.45*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[7]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(59,43) = 0.7071428571428572*m0rOther[7]*uOther[15]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[15]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[14]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[13]*mnuOther+0.2*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.2*m0rOther[5]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[0]*uOther[12]*mnuOther-0.4472135954999579*m1rOther[12]*mnuOther+0.81*m0rOther[3]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[10]*mnuOther+0.45*m0rOther[2]*uOther[10]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[9]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(59,44) = 0.3513821107499669*m0rOther[6]*uOther[15]*mnuOther+0.2*m0rOther[2]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[14]*mnuOther+0.2*m0rOther[3]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[12]*mnuOther+0.2*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.2*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(59,45) = 0.3513821107499669*m0rOther[6]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[14]*mnuOther+0.2*m0rOther[1]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[13]*mnuOther+0.2*m0rOther[3]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[11]*mnuOther+0.2*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.2*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(59,46) = 0.3513821107499669*m0rOther[5]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[15]*mnuOther+0.2*m0rOther[0]*uOther[15]*mnuOther-0.4*m1rOther[15]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[13]*mnuOther+0.2*m0rOther[1]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[12]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[10]*mnuOther+0.2*m0rOther[5]*uOther[9]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.447213595499958*m1rOther[9]*mnuOther+0.2*m0rOther[7]*uOther[8]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(59,47) = 0.7071428571428572*m0rOther[3]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[14]*mnuOther+0.2*m0rOther[0]*uOther[14]*mnuOther-0.4*m1rOther[14]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[12]*mnuOther+0.2*m0rOther[2]*uOther[12]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[10]*mnuOther+0.2*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.447213595499958*m1rOther[10]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[9]*mnuOther+0.2*m0rOther[6]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(60,40) = 0.223606797749979*m0rOther[7]*uOther[15]*mnuOther+0.159719141249985*m0rOther[6]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[4]*uOther[8]*mnuOther; 
  data->AEM_S(60,41) = 0.223606797749979*m0rOther[5]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[7]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(60,42) = 0.2*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[14]*mnuOther-0.5000000000000001*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[6]*uOther[13]*mnuOther+0.159719141249985*m0rOther[6]*uOther[12]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[8]*mnuOther; 
  data->AEM_S(60,43) = 0.3513821107499669*m0rOther[6]*uOther[15]*mnuOther+0.2*m0rOther[2]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[14]*mnuOther+0.2*m0rOther[3]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[12]*mnuOther+0.2*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.2*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(60,44) = 0.3928571428571428*m0rOther[7]*uOther[15]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[14]*mnuOther+0.159719141249985*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[5]*uOther[13]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[12]*mnuOther+0.159719141249985*m0rOther[0]*uOther[12]*mnuOther-0.31943828249997*m1rOther[12]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[11]*mnuOther+0.159719141249985*m0rOther[6]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[9]*mnuOther+0.159719141249985*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(60,45) = 0.1428571428571428*m0rOther[7]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[5]*uOther[12]*mnuOther+0.2*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[7]*uOther[9]*mnuOther; 
  data->AEM_S(60,46) = 0.3513821107499669*m0rOther[3]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[14]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[14]*mnuOther+0.159719141249985*m0rOther[0]*uOther[14]*mnuOther-0.31943828249997*m1rOther[14]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[12]*mnuOther+0.159719141249985*m0rOther[2]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[10]*mnuOther-0.5000000000000001*m1rOther[10]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[9]*mnuOther+0.159719141249985*m0rOther[6]*uOther[8]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(60,47) = 0.1428571428571428*m0rOther[5]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[11]*mnuOther+0.2*m0rOther[2]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[7]*uOther[8]*mnuOther; 
  data->AEM_S(61,40) = 0.159719141249985*m0rOther[7]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[6]*uOther[14]*mnuOther+0.159719141249985*m0rOther[5]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[5]*uOther[8]*mnuOther; 
  data->AEM_S(61,41) = 0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[15]*mnuOther-0.5000000000000001*m1rOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.159719141249985*m0rOther[7]*uOther[13]*mnuOther+0.25*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.25*m0rOther[5]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[8]*mnuOther; 
  data->AEM_S(61,42) = 0.3928571428571429*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[6]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.4472135954999579*m1rOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(61,43) = 0.3513821107499669*m0rOther[6]*uOther[15]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[14]*mnuOther+0.2*m0rOther[1]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[13]*mnuOther+0.2*m0rOther[3]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[11]*mnuOther+0.2*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.2*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(61,44) = 0.1428571428571428*m0rOther[7]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[5]*uOther[12]*mnuOther+0.2*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[7]*uOther[9]*mnuOther; 
  data->AEM_S(61,45) = 0.5357142857142857*m0rOther[7]*uOther[15]*mnuOther+0.159719141249985*m0rOther[1]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[14]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[13]*mnuOther+0.159719141249985*m0rOther[0]*uOther[13]*mnuOther-0.31943828249997*m1rOther[13]*mnuOther+0.25*m0rOther[4]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[10]*mnuOther+0.159719141249985*m0rOther[7]*uOther[9]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(61,46) = 0.3513821107499669*m0rOther[3]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.4472135954999579*m1rOther[14]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[13]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[11]*mnuOther+0.2*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[10]*mnuOther+0.2*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[6]*uOther[8]*mnuOther; 
  data->AEM_S(61,47) = 0.5357142857142857*m0rOther[5]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[15]*mnuOther+0.159719141249985*m0rOther[0]*uOther[15]*mnuOther-0.31943828249997*m1rOther[15]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[14]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[13]*mnuOther+0.159719141249985*m0rOther[1]*uOther[13]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[10]*mnuOther+0.159719141249985*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[9]*mnuOther-0.5000000000000001*m1rOther[9]*mnuOther+0.159719141249985*m0rOther[7]*uOther[8]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(62,40) = 0.2*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.159719141249985*m0rOther[4]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[6]*uOther[13]*mnuOther+0.159719141249985*m0rOther[6]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[12]*mnuOther+0.2*m0rOther[7]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[6]*uOther[8]*mnuOther; 
  data->AEM_S(62,41) = 0.351382110749967*m0rOther[6]*uOther[15]*mnuOther+0.2*m0rOther[2]*uOther[15]*mnuOther+0.351382110749967*m0rOther[7]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[14]*mnuOther+0.2*m0rOther[3]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[12]*mnuOther+0.2*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.447213595499958*m1rOther[11]*mnuOther+0.2*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(62,42) = 0.351382110749967*m0rOther[7]*uOther[15]*mnuOther+0.2*m0rOther[1]*uOther[15]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[14]*mnuOther+0.45*m0rOther[2]*uOther[14]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.159719141249985*m0rOther[4]*uOther[12]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[12]*mnuOther-0.5000000000000001*m1rOther[12]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[11]*mnuOther+0.45*m0rOther[6]*uOther[10]*mnuOther+0.2*m0rOther[7]*uOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[8]*mnuOther; 
  data->AEM_S(62,43) = 0.3513821107499669*m0rOther[5]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[15]*mnuOther+0.2*m0rOther[0]*uOther[15]*mnuOther-0.4*m1rOther[15]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[13]*mnuOther+0.2*m0rOther[1]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[12]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[10]*mnuOther+0.2*m0rOther[5]*uOther[9]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.447213595499958*m1rOther[9]*mnuOther+0.2*m0rOther[7]*uOther[8]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(62,44) = 0.3513821107499669*m0rOther[3]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[14]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[14]*mnuOther+0.159719141249985*m0rOther[0]*uOther[14]*mnuOther-0.31943828249997*m1rOther[14]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[12]*mnuOther+0.159719141249985*m0rOther[2]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[5]*uOther[10]*mnuOther+0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[10]*mnuOther-0.5000000000000001*m1rOther[10]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[9]*mnuOther+0.159719141249985*m0rOther[6]*uOther[8]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(62,45) = 0.3513821107499669*m0rOther[3]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.4472135954999579*m1rOther[14]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[13]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[12]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[11]*mnuOther+0.2*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[4]*uOther[10]*mnuOther+0.2*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[6]*uOther[8]*mnuOther; 
  data->AEM_S(62,46) = 0.6173469387755102*m0rOther[7]*uOther[15]*mnuOther+0.351382110749967*m0rOther[1]*uOther[15]*mnuOther+0.9642857142857143*m0rOther[6]*uOther[14]*mnuOther+0.2874944542499729*m0rOther[2]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[13]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[12]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[12]*mnuOther+0.159719141249985*m0rOther[0]*uOther[12]*mnuOther-0.31943828249997*m1rOther[12]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[11]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[10]*mnuOther+0.45*m0rOther[2]*uOther[10]*mnuOther+0.351382110749967*m0rOther[7]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.159719141249985*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(62,47) = 0.6173469387755102*m0rOther[6]*uOther[15]*mnuOther+0.351382110749967*m0rOther[2]*uOther[15]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[14]*mnuOther+0.351382110749967*m0rOther[1]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[11]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[11]*mnuOther+0.2*m0rOther[0]*uOther[11]*mnuOther-0.4*m1rOther[11]*mnuOther+0.351382110749967*m0rOther[7]*uOther[10]*mnuOther+0.2*m0rOther[1]*uOther[10]*mnuOther+0.351382110749967*m0rOther[6]*uOther[9]*mnuOther+0.2*m0rOther[2]*uOther[9]*mnuOther+0.2*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(63,40) = 0.159719141249985*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.2*m0rOther[3]*uOther[14]*mnuOther+0.159719141249985*m0rOther[7]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[7]*uOther[12]*mnuOther+0.2*m0rOther[6]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[9]*mnuOther+0.25*m0rOther[7]*uOther[8]*mnuOther; 
  data->AEM_S(63,41) = 0.2874944542499729*m0rOther[7]*uOther[15]*mnuOther+0.45*m0rOther[1]*uOther[15]*mnuOther+0.351382110749967*m0rOther[6]*uOther[14]*mnuOther+0.2*m0rOther[2]*uOther[14]*mnuOther+0.159719141249985*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[13]*mnuOther-0.5000000000000001*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[11]*mnuOther+0.2*m0rOther[6]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[7]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[8]*mnuOther; 
  data->AEM_S(63,42) = 0.351382110749967*m0rOther[6]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[15]*mnuOther+0.351382110749967*m0rOther[7]*uOther[14]*mnuOther+0.2*m0rOther[1]*uOther[14]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[13]*mnuOther+0.2*m0rOther[3]*uOther[12]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[11]*mnuOther+0.2*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.447213595499958*m1rOther[11]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.2*m0rOther[6]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(63,43) = 0.7071428571428572*m0rOther[3]*uOther[15]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[14]*mnuOther+0.2*m0rOther[0]*uOther[14]*mnuOther-0.4*m1rOther[14]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[13]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[12]*mnuOther+0.2*m0rOther[2]*uOther[12]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[11]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[10]*mnuOther+0.2*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.447213595499958*m1rOther[10]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[9]*mnuOther+0.2*m0rOther[6]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther; 
  data->AEM_S(63,44) = 0.1428571428571428*m0rOther[5]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[14]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[11]*mnuOther+0.2*m0rOther[2]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[7]*uOther[8]*mnuOther; 
  data->AEM_S(63,45) = 0.5357142857142857*m0rOther[5]*uOther[15]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[15]*mnuOther+0.159719141249985*m0rOther[0]*uOther[15]*mnuOther-0.31943828249997*m1rOther[15]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[14]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[13]*mnuOther+0.159719141249985*m0rOther[1]*uOther[13]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[12]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[11]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[10]*mnuOther+0.159719141249985*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[9]*mnuOther-0.5000000000000001*m1rOther[9]*mnuOther+0.159719141249985*m0rOther[7]*uOther[8]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[8]*mnuOther; 
  data->AEM_S(63,46) = 0.6173469387755102*m0rOther[6]*uOther[15]*mnuOther+0.351382110749967*m0rOther[2]*uOther[15]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[14]*mnuOther+0.351382110749967*m0rOther[1]*uOther[14]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[13]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[12]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[11]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[11]*mnuOther+0.2*m0rOther[0]*uOther[11]*mnuOther-0.4*m1rOther[11]*mnuOther+0.351382110749967*m0rOther[7]*uOther[10]*mnuOther+0.2*m0rOther[1]*uOther[10]*mnuOther+0.351382110749967*m0rOther[6]*uOther[9]*mnuOther+0.2*m0rOther[2]*uOther[9]*mnuOther+0.2*m0rOther[3]*uOther[8]*mnuOther; 
  data->AEM_S(63,47) = 0.9642857142857143*m0rOther[7]*uOther[15]*mnuOther+0.2874944542499729*m0rOther[1]*uOther[15]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[14]*mnuOther+0.351382110749967*m0rOther[2]*uOther[14]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[13]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[13]*mnuOther+0.159719141249985*m0rOther[0]*uOther[13]*mnuOther-0.31943828249997*m1rOther[13]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[12]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[12]*mnuOther+0.223606797749979*m0rOther[0]*uOther[12]*mnuOther-0.4472135954999579*m1rOther[12]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[11]*mnuOther+0.351382110749967*m0rOther[6]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[10]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[9]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.159719141249985*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfY-m0Self*m1OtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[8] = 0.5*m0rOther[7]*m1rSelf[15]-0.5*m0rSelf[7]*m1rOther[15]+0.5*m0rOther[6]*m1rSelf[14]-0.5*m0rSelf[6]*m1rOther[14]+0.5*m0rOther[5]*m1rSelf[13]-0.5*m0rSelf[5]*m1rOther[13]+0.5*m0rOther[4]*m1rSelf[12]-0.5*m0rSelf[4]*m1rOther[12]+0.5*m0rOther[3]*m1rSelf[11]-0.5*m0rSelf[3]*m1rOther[11]+0.5*m0rOther[2]*m1rSelf[10]-0.5*m0rSelf[2]*m1rOther[10]+0.5*m0rOther[1]*m1rSelf[9]-0.5*m0rSelf[1]*m1rOther[9]+0.5*m0rOther[0]*m1rSelf[8]-0.5*m0rSelf[0]*m1rOther[8]; 
  m1EffD[9] = 0.5000000000000001*m0rOther[5]*m1rSelf[15]-0.5000000000000001*m0rSelf[5]*m1rOther[15]+0.447213595499958*m0rOther[3]*m1rSelf[14]-0.447213595499958*m0rSelf[3]*m1rOther[14]+0.5000000000000001*m0rOther[7]*m1rSelf[13]-0.5000000000000001*m0rSelf[7]*m1rOther[13]+0.4472135954999579*m0rOther[1]*m1rSelf[12]-0.4472135954999579*m0rSelf[1]*m1rOther[12]+0.447213595499958*m0rOther[6]*m1rSelf[11]+0.5*m0rOther[2]*m1rSelf[11]-0.447213595499958*m0rSelf[6]*m1rOther[11]-0.5*m0rSelf[2]*m1rOther[11]+0.5*m0rOther[3]*m1rSelf[10]-0.5*m0rSelf[3]*m1rOther[10]+0.4472135954999579*m0rOther[4]*m1rSelf[9]+0.5*m0rOther[0]*m1rSelf[9]-0.4472135954999579*m0rSelf[4]*m1rOther[9]-0.5*m0rSelf[0]*m1rOther[9]+0.5*m0rOther[1]*m1rSelf[8]-0.5*m0rSelf[1]*m1rOther[8]; 
  m1EffD[10] = 0.447213595499958*m0rOther[3]*m1rSelf[15]-0.447213595499958*m0rSelf[3]*m1rOther[15]+0.5000000000000001*m0rOther[4]*m1rSelf[14]-0.5000000000000001*m0rSelf[4]*m1rOther[14]+0.4472135954999579*m0rOther[2]*m1rSelf[13]-0.4472135954999579*m0rSelf[2]*m1rOther[13]+0.5000000000000001*m0rOther[6]*m1rSelf[12]-0.5000000000000001*m0rSelf[6]*m1rOther[12]+0.447213595499958*m0rOther[7]*m1rSelf[11]+0.5*m0rOther[1]*m1rSelf[11]-0.447213595499958*m0rSelf[7]*m1rOther[11]-0.5*m0rSelf[1]*m1rOther[11]+0.4472135954999579*m0rOther[5]*m1rSelf[10]+0.5*m0rOther[0]*m1rSelf[10]-0.4472135954999579*m0rSelf[5]*m1rOther[10]-0.5*m0rSelf[0]*m1rOther[10]+0.5*m0rOther[3]*m1rSelf[9]-0.5*m0rSelf[3]*m1rOther[9]+0.5*m0rOther[2]*m1rSelf[8]-0.5*m0rSelf[2]*m1rOther[8]; 
  m1EffD[11] = 0.4*m0rOther[6]*m1rSelf[15]+0.447213595499958*m0rOther[2]*m1rSelf[15]-0.4*m0rSelf[6]*m1rOther[15]-0.447213595499958*m0rSelf[2]*m1rOther[15]+0.4*m0rOther[7]*m1rSelf[14]+0.447213595499958*m0rOther[1]*m1rSelf[14]-0.4*m0rSelf[7]*m1rOther[14]-0.447213595499958*m0rSelf[1]*m1rOther[14]+0.4472135954999579*m0rOther[3]*m1rSelf[13]-0.4472135954999579*m0rSelf[3]*m1rOther[13]+0.4472135954999579*m0rOther[3]*m1rSelf[12]-0.4472135954999579*m0rSelf[3]*m1rOther[12]+0.4472135954999579*m0rOther[5]*m1rSelf[11]+0.4472135954999579*m0rOther[4]*m1rSelf[11]+0.5*m0rOther[0]*m1rSelf[11]-0.4472135954999579*m0rSelf[5]*m1rOther[11]-0.4472135954999579*m0rSelf[4]*m1rOther[11]-0.5*m0rSelf[0]*m1rOther[11]+0.447213595499958*m0rOther[7]*m1rSelf[10]+0.5*m0rOther[1]*m1rSelf[10]-0.447213595499958*m0rSelf[7]*m1rOther[10]-0.5*m0rSelf[1]*m1rOther[10]+0.447213595499958*m0rOther[6]*m1rSelf[9]+0.5*m0rOther[2]*m1rSelf[9]-0.447213595499958*m0rSelf[6]*m1rOther[9]-0.5*m0rSelf[2]*m1rOther[9]+0.5*m0rOther[3]*m1rSelf[8]-0.5*m0rSelf[3]*m1rOther[8]; 
  m1EffD[12] = 0.4472135954999579*m0rOther[7]*m1rSelf[15]-0.4472135954999579*m0rSelf[7]*m1rOther[15]+0.31943828249997*m0rOther[6]*m1rSelf[14]+0.5000000000000001*m0rOther[2]*m1rSelf[14]-0.31943828249997*m0rSelf[6]*m1rOther[14]-0.5000000000000001*m0rSelf[2]*m1rOther[14]+0.31943828249997*m0rOther[4]*m1rSelf[12]+0.5*m0rOther[0]*m1rSelf[12]-0.31943828249997*m0rSelf[4]*m1rOther[12]-0.5*m0rSelf[0]*m1rOther[12]+0.4472135954999579*m0rOther[3]*m1rSelf[11]-0.4472135954999579*m0rSelf[3]*m1rOther[11]+0.5000000000000001*m0rOther[6]*m1rSelf[10]-0.5000000000000001*m0rSelf[6]*m1rOther[10]+0.4472135954999579*m0rOther[1]*m1rSelf[9]-0.4472135954999579*m0rSelf[1]*m1rOther[9]+0.5*m0rOther[4]*m1rSelf[8]-0.5*m0rSelf[4]*m1rOther[8]; 
  m1EffD[13] = 0.31943828249997*m0rOther[7]*m1rSelf[15]+0.5000000000000001*m0rOther[1]*m1rSelf[15]-0.31943828249997*m0rSelf[7]*m1rOther[15]-0.5000000000000001*m0rSelf[1]*m1rOther[15]+0.4472135954999579*m0rOther[6]*m1rSelf[14]-0.4472135954999579*m0rSelf[6]*m1rOther[14]+0.31943828249997*m0rOther[5]*m1rSelf[13]+0.5*m0rOther[0]*m1rSelf[13]-0.31943828249997*m0rSelf[5]*m1rOther[13]-0.5*m0rSelf[0]*m1rOther[13]+0.4472135954999579*m0rOther[3]*m1rSelf[11]-0.4472135954999579*m0rSelf[3]*m1rOther[11]+0.4472135954999579*m0rOther[2]*m1rSelf[10]-0.4472135954999579*m0rSelf[2]*m1rOther[10]+0.5000000000000001*m0rOther[7]*m1rSelf[9]-0.5000000000000001*m0rSelf[7]*m1rOther[9]+0.5*m0rOther[5]*m1rSelf[8]-0.5*m0rSelf[5]*m1rOther[8]; 
  m1EffD[14] = 0.4*m0rOther[3]*m1rSelf[15]-0.4*m0rSelf[3]*m1rOther[15]+0.4472135954999579*m0rOther[5]*m1rSelf[14]+0.31943828249997*m0rOther[4]*m1rSelf[14]+0.5*m0rOther[0]*m1rSelf[14]-0.4472135954999579*m0rSelf[5]*m1rOther[14]-0.31943828249997*m0rSelf[4]*m1rOther[14]-0.5*m0rSelf[0]*m1rOther[14]+0.4472135954999579*m0rOther[6]*m1rSelf[13]-0.4472135954999579*m0rSelf[6]*m1rOther[13]+0.31943828249997*m0rOther[6]*m1rSelf[12]+0.5000000000000001*m0rOther[2]*m1rSelf[12]-0.31943828249997*m0rSelf[6]*m1rOther[12]-0.5000000000000001*m0rSelf[2]*m1rOther[12]+0.4*m0rOther[7]*m1rSelf[11]+0.447213595499958*m0rOther[1]*m1rSelf[11]-0.4*m0rSelf[7]*m1rOther[11]-0.447213595499958*m0rSelf[1]*m1rOther[11]+0.5000000000000001*m0rOther[4]*m1rSelf[10]-0.5000000000000001*m0rSelf[4]*m1rOther[10]+0.447213595499958*m0rOther[3]*m1rSelf[9]-0.447213595499958*m0rSelf[3]*m1rOther[9]+0.5*m0rOther[6]*m1rSelf[8]-0.5*m0rSelf[6]*m1rOther[8]; 
  m1EffD[15] = 0.31943828249997*m0rOther[5]*m1rSelf[15]+0.4472135954999579*m0rOther[4]*m1rSelf[15]+0.5*m0rOther[0]*m1rSelf[15]-0.31943828249997*m0rSelf[5]*m1rOther[15]-0.4472135954999579*m0rSelf[4]*m1rOther[15]-0.5*m0rSelf[0]*m1rOther[15]+0.4*m0rOther[3]*m1rSelf[14]-0.4*m0rSelf[3]*m1rOther[14]+0.31943828249997*m0rOther[7]*m1rSelf[13]+0.5000000000000001*m0rOther[1]*m1rSelf[13]-0.31943828249997*m0rSelf[7]*m1rOther[13]-0.5000000000000001*m0rSelf[1]*m1rOther[13]+0.4472135954999579*m0rOther[7]*m1rSelf[12]-0.4472135954999579*m0rSelf[7]*m1rOther[12]+0.4*m0rOther[6]*m1rSelf[11]+0.447213595499958*m0rOther[2]*m1rSelf[11]-0.4*m0rSelf[6]*m1rOther[11]-0.447213595499958*m0rSelf[2]*m1rOther[11]+0.447213595499958*m0rOther[3]*m1rSelf[10]-0.447213595499958*m0rSelf[3]*m1rOther[10]+0.5000000000000001*m0rOther[5]*m1rSelf[9]-0.5000000000000001*m0rSelf[5]*m1rOther[9]+0.5*m0rOther[7]*m1rSelf[8]-0.5*m0rSelf[7]*m1rOther[8]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[8],m1EffD[9],m1EffD[10],m1EffD[11],m1EffD[12],m1EffD[13],m1EffD[14],m1EffD[15]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+8,8,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[8] += (-2.0*m1EffD[8]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
  m1Relax[9] += (-2.0*m1EffD[9]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[9]*mnuSelf-1.0*m1rOther[9]*mnuOther; 
  m1Relax[10] += (-2.0*m1EffD[10]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[10]*mnuSelf-1.0*m1rOther[10]*mnuOther; 
  m1Relax[11] += (-2.0*m1EffD[11]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[11]*mnuSelf-1.0*m1rOther[11]*mnuOther; 
  m1Relax[12] += (-2.0*m1EffD[12]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[12]*mnuSelf-1.0*m1rOther[12]*mnuOther; 
  m1Relax[13] += (-2.0*m1EffD[13]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[13]*mnuSelf-1.0*m1rOther[13]*mnuOther; 
  m1Relax[14] += (-2.0*m1EffD[14]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[14]*mnuSelf-1.0*m1rOther[14]*mnuOther; 
  m1Relax[15] += (-2.0*m1EffD[15]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[15]*mnuSelf-1.0*m1rOther[15]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(48,16) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(48,17) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(48,18) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(48,19) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(48,20) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(48,21) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(48,22) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(48,23) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(49,16) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(49,17) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(49,18) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(49,19) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(49,20) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(49,21) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(49,22) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(49,23) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(50,16) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(50,17) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(50,18) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(50,19) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(50,20) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(50,21) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(50,22) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(50,23) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(51,16) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(51,17) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(51,18) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(51,19) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(51,20) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(51,21) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(51,22) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(51,23) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(52,16) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(52,17) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(52,18) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(52,19) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(52,20) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(52,22) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(52,23) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(53,16) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(53,17) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(53,18) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(53,19) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(53,21) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(53,22) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(53,23) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(54,16) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(54,17) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(54,18) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(54,19) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(54,20) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(54,21) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(54,22) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(54,23) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(55,16) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(55,17) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(55,18) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(55,19) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(55,20) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(55,21) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(55,22) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(55,23) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(48,24) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(48,25) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(48,26) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(48,27) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(48,28) = -0.5*cMSelf[20]*mnuSelf; 
  data->AEM_S(48,29) = -0.5*cMSelf[21]*mnuSelf; 
  data->AEM_S(48,30) = -0.5*cMSelf[22]*mnuSelf; 
  data->AEM_S(48,31) = -0.5*cMSelf[23]*mnuSelf; 
  data->AEM_S(49,24) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(49,25) = (-0.4472135954999579*cMSelf[20]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(49,26) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(49,27) = (-0.447213595499958*cMSelf[22]*mnuSelf)-0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(49,28) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(49,29) = -0.5000000000000001*cMSelf[23]*mnuSelf; 
  data->AEM_S(49,30) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(49,31) = -0.5000000000000001*cMSelf[21]*mnuSelf; 
  data->AEM_S(50,24) = -0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(50,25) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(50,26) = (-0.4472135954999579*cMSelf[21]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(50,27) = (-0.447213595499958*cMSelf[23]*mnuSelf)-0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(50,28) = -0.5000000000000001*cMSelf[22]*mnuSelf; 
  data->AEM_S(50,29) = -0.4472135954999579*cMSelf[18]*mnuSelf; 
  data->AEM_S(50,30) = -0.5000000000000001*cMSelf[20]*mnuSelf; 
  data->AEM_S(50,31) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(51,24) = -0.5*cMSelf[19]*mnuSelf; 
  data->AEM_S(51,25) = (-0.447213595499958*cMSelf[22]*mnuSelf)-0.5*cMSelf[18]*mnuSelf; 
  data->AEM_S(51,26) = (-0.447213595499958*cMSelf[23]*mnuSelf)-0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(51,27) = (-0.4472135954999579*cMSelf[21]*mnuSelf)-0.4472135954999579*cMSelf[20]*mnuSelf-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(51,28) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(51,29) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(51,30) = (-0.4*cMSelf[23]*mnuSelf)-0.447213595499958*cMSelf[17]*mnuSelf; 
  data->AEM_S(51,31) = (-0.4*cMSelf[22]*mnuSelf)-0.447213595499958*cMSelf[18]*mnuSelf; 
  data->AEM_S(52,24) = -0.5*cMSelf[20]*mnuSelf; 
  data->AEM_S(52,25) = -0.4472135954999579*cMSelf[17]*mnuSelf; 
  data->AEM_S(52,26) = -0.5000000000000001*cMSelf[22]*mnuSelf; 
  data->AEM_S(52,27) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(52,28) = (-0.31943828249997*cMSelf[20]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(52,30) = (-0.31943828249997*cMSelf[22]*mnuSelf)-0.5000000000000001*cMSelf[18]*mnuSelf; 
  data->AEM_S(52,31) = -0.4472135954999579*cMSelf[23]*mnuSelf; 
  data->AEM_S(53,24) = -0.5*cMSelf[21]*mnuSelf; 
  data->AEM_S(53,25) = -0.5000000000000001*cMSelf[23]*mnuSelf; 
  data->AEM_S(53,26) = -0.4472135954999579*cMSelf[18]*mnuSelf; 
  data->AEM_S(53,27) = -0.4472135954999579*cMSelf[19]*mnuSelf; 
  data->AEM_S(53,29) = (-0.31943828249997*cMSelf[21]*mnuSelf)-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(53,30) = -0.4472135954999579*cMSelf[22]*mnuSelf; 
  data->AEM_S(53,31) = (-0.31943828249997*cMSelf[23]*mnuSelf)-0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(54,24) = -0.5*cMSelf[22]*mnuSelf; 
  data->AEM_S(54,25) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(54,26) = -0.5000000000000001*cMSelf[20]*mnuSelf; 
  data->AEM_S(54,27) = (-0.4*cMSelf[23]*mnuSelf)-0.447213595499958*cMSelf[17]*mnuSelf; 
  data->AEM_S(54,28) = (-0.31943828249997*cMSelf[22]*mnuSelf)-0.5000000000000001*cMSelf[18]*mnuSelf; 
  data->AEM_S(54,29) = -0.4472135954999579*cMSelf[22]*mnuSelf; 
  data->AEM_S(54,30) = (-0.4472135954999579*cMSelf[21]*mnuSelf)-0.31943828249997*cMSelf[20]*mnuSelf-0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(54,31) = -0.4*cMSelf[19]*mnuSelf; 
  data->AEM_S(55,24) = -0.5*cMSelf[23]*mnuSelf; 
  data->AEM_S(55,25) = -0.5000000000000001*cMSelf[21]*mnuSelf; 
  data->AEM_S(55,26) = -0.447213595499958*cMSelf[19]*mnuSelf; 
  data->AEM_S(55,27) = (-0.4*cMSelf[22]*mnuSelf)-0.447213595499958*cMSelf[18]*mnuSelf; 
  data->AEM_S(55,28) = -0.4472135954999579*cMSelf[23]*mnuSelf; 
  data->AEM_S(55,29) = (-0.31943828249997*cMSelf[23]*mnuSelf)-0.5000000000000001*cMSelf[17]*mnuSelf; 
  data->AEM_S(55,30) = -0.4*cMSelf[19]*mnuSelf; 
  data->AEM_S(55,31) = (-0.31943828249997*cMSelf[21]*mnuSelf)-0.4472135954999579*cMSelf[20]*mnuSelf-0.5*cMSelf[16]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(48,48) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(48,49) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(48,50) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(48,51) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(48,52) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(48,53) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(48,54) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(48,55) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(49,48) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(49,49) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(49,50) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(49,51) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(49,52) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(49,53) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(49,54) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(49,55) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(50,48) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(50,49) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(50,50) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(50,51) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(50,52) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(50,53) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(50,54) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(50,55) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(51,48) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(51,49) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(51,50) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(51,51) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(51,52) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(51,53) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(51,54) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(51,55) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(52,48) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(52,49) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(52,50) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(52,51) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(52,52) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(52,54) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(52,55) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(53,48) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(53,49) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(53,50) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(53,51) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(53,53) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(53,54) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(53,55) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(54,48) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(54,49) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(54,50) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(54,51) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(54,52) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(54,53) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(54,54) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.31943828249997*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(54,55) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(55,48) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(55,49) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(55,50) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(55,51) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(55,52) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(55,53) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(55,54) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(55,55) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(48,56) = 0.5*cMOther[16]*mnuOther; 
  data->AEM_S(48,57) = 0.5*cMOther[17]*mnuOther; 
  data->AEM_S(48,58) = 0.5*cMOther[18]*mnuOther; 
  data->AEM_S(48,59) = 0.5*cMOther[19]*mnuOther; 
  data->AEM_S(48,60) = 0.5*cMOther[20]*mnuOther; 
  data->AEM_S(48,61) = 0.5*cMOther[21]*mnuOther; 
  data->AEM_S(48,62) = 0.5*cMOther[22]*mnuOther; 
  data->AEM_S(48,63) = 0.5*cMOther[23]*mnuOther; 
  data->AEM_S(49,56) = 0.5*cMOther[17]*mnuOther; 
  data->AEM_S(49,57) = 0.4472135954999579*cMOther[20]*mnuOther+0.5*cMOther[16]*mnuOther; 
  data->AEM_S(49,58) = 0.5*cMOther[19]*mnuOther; 
  data->AEM_S(49,59) = 0.447213595499958*cMOther[22]*mnuOther+0.5*cMOther[18]*mnuOther; 
  data->AEM_S(49,60) = 0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(49,61) = 0.5000000000000001*cMOther[23]*mnuOther; 
  data->AEM_S(49,62) = 0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(49,63) = 0.5000000000000001*cMOther[21]*mnuOther; 
  data->AEM_S(50,56) = 0.5*cMOther[18]*mnuOther; 
  data->AEM_S(50,57) = 0.5*cMOther[19]*mnuOther; 
  data->AEM_S(50,58) = 0.4472135954999579*cMOther[21]*mnuOther+0.5*cMOther[16]*mnuOther; 
  data->AEM_S(50,59) = 0.447213595499958*cMOther[23]*mnuOther+0.5*cMOther[17]*mnuOther; 
  data->AEM_S(50,60) = 0.5000000000000001*cMOther[22]*mnuOther; 
  data->AEM_S(50,61) = 0.4472135954999579*cMOther[18]*mnuOther; 
  data->AEM_S(50,62) = 0.5000000000000001*cMOther[20]*mnuOther; 
  data->AEM_S(50,63) = 0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(51,56) = 0.5*cMOther[19]*mnuOther; 
  data->AEM_S(51,57) = 0.447213595499958*cMOther[22]*mnuOther+0.5*cMOther[18]*mnuOther; 
  data->AEM_S(51,58) = 0.447213595499958*cMOther[23]*mnuOther+0.5*cMOther[17]*mnuOther; 
  data->AEM_S(51,59) = 0.4472135954999579*cMOther[21]*mnuOther+0.4472135954999579*cMOther[20]*mnuOther+0.5*cMOther[16]*mnuOther; 
  data->AEM_S(51,60) = 0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(51,61) = 0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(51,62) = 0.4*cMOther[23]*mnuOther+0.447213595499958*cMOther[17]*mnuOther; 
  data->AEM_S(51,63) = 0.4*cMOther[22]*mnuOther+0.447213595499958*cMOther[18]*mnuOther; 
  data->AEM_S(52,56) = 0.5*cMOther[20]*mnuOther; 
  data->AEM_S(52,57) = 0.4472135954999579*cMOther[17]*mnuOther; 
  data->AEM_S(52,58) = 0.5000000000000001*cMOther[22]*mnuOther; 
  data->AEM_S(52,59) = 0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(52,60) = 0.31943828249997*cMOther[20]*mnuOther+0.5*cMOther[16]*mnuOther; 
  data->AEM_S(52,62) = 0.31943828249997*cMOther[22]*mnuOther+0.5000000000000001*cMOther[18]*mnuOther; 
  data->AEM_S(52,63) = 0.4472135954999579*cMOther[23]*mnuOther; 
  data->AEM_S(53,56) = 0.5*cMOther[21]*mnuOther; 
  data->AEM_S(53,57) = 0.5000000000000001*cMOther[23]*mnuOther; 
  data->AEM_S(53,58) = 0.4472135954999579*cMOther[18]*mnuOther; 
  data->AEM_S(53,59) = 0.4472135954999579*cMOther[19]*mnuOther; 
  data->AEM_S(53,61) = 0.31943828249997*cMOther[21]*mnuOther+0.5*cMOther[16]*mnuOther; 
  data->AEM_S(53,62) = 0.4472135954999579*cMOther[22]*mnuOther; 
  data->AEM_S(53,63) = 0.31943828249997*cMOther[23]*mnuOther+0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(54,56) = 0.5*cMOther[22]*mnuOther; 
  data->AEM_S(54,57) = 0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(54,58) = 0.5000000000000001*cMOther[20]*mnuOther; 
  data->AEM_S(54,59) = 0.4*cMOther[23]*mnuOther+0.447213595499958*cMOther[17]*mnuOther; 
  data->AEM_S(54,60) = 0.31943828249997*cMOther[22]*mnuOther+0.5000000000000001*cMOther[18]*mnuOther; 
  data->AEM_S(54,61) = 0.4472135954999579*cMOther[22]*mnuOther; 
  data->AEM_S(54,62) = 0.4472135954999579*cMOther[21]*mnuOther+0.31943828249997*cMOther[20]*mnuOther+0.5*cMOther[16]*mnuOther; 
  data->AEM_S(54,63) = 0.4*cMOther[19]*mnuOther; 
  data->AEM_S(55,56) = 0.5*cMOther[23]*mnuOther; 
  data->AEM_S(55,57) = 0.5000000000000001*cMOther[21]*mnuOther; 
  data->AEM_S(55,58) = 0.447213595499958*cMOther[19]*mnuOther; 
  data->AEM_S(55,59) = 0.4*cMOther[22]*mnuOther+0.447213595499958*cMOther[18]*mnuOther; 
  data->AEM_S(55,60) = 0.4472135954999579*cMOther[23]*mnuOther; 
  data->AEM_S(55,61) = 0.31943828249997*cMOther[23]*mnuOther+0.5000000000000001*cMOther[17]*mnuOther; 
  data->AEM_S(55,62) = 0.4*cMOther[19]*mnuOther; 
  data->AEM_S(55,63) = 0.31943828249997*cMOther[21]*mnuOther+0.4472135954999579*cMOther[20]*mnuOther+0.5*cMOther[16]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfZ-uSelfZ*m0Self) and uCrossSelfZ ... // 
  data->AEM_S(56,16) = (-0.25*m0rSelf[7]*uSelf[23]*mnuSelf)-0.25*m0rSelf[6]*uSelf[22]*mnuSelf-0.25*m0rSelf[5]*uSelf[21]*mnuSelf-0.25*m0rSelf[4]*uSelf[20]*mnuSelf-0.25*m0rSelf[3]*uSelf[19]*mnuSelf-0.25*m0rSelf[2]*uSelf[18]*mnuSelf-0.25*m0rSelf[1]*uSelf[17]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(56,17) = (-0.2500000000000001*m0rSelf[5]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[22]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[19]*mnuSelf-0.25*m0rSelf[2]*uSelf[19]*mnuSelf-0.25*m0rSelf[3]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.25*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(56,18) = (-0.223606797749979*m0rSelf[3]*uSelf[23]*mnuSelf)-0.2500000000000001*m0rSelf[4]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[21]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[19]*mnuSelf-0.25*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[18]*mnuSelf-0.25*m0rSelf[0]*uSelf[18]*mnuSelf+0.5*m1rSelf[18]*mnuSelf-0.25*m0rSelf[3]*uSelf[17]*mnuSelf-0.25*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(56,19) = (-0.2*m0rSelf[6]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[23]*mnuSelf-0.2*m0rSelf[7]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[19]*mnuSelf-0.25*m0rSelf[0]*uSelf[19]*mnuSelf+0.5*m1rSelf[19]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[18]*mnuSelf-0.25*m0rSelf[1]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[17]*mnuSelf-0.25*m0rSelf[2]*uSelf[17]*mnuSelf-0.25*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(56,20) = (-0.223606797749979*m0rSelf[7]*uSelf[23]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[22]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[20]*mnuSelf-0.25*m0rSelf[0]*uSelf[20]*mnuSelf+0.5*m1rSelf[20]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[19]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf-0.25*m0rSelf[4]*uSelf[16]*mnuSelf; 
  data->AEM_S(56,21) = (-0.159719141249985*m0rSelf[7]*uSelf[23]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[23]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[21]*mnuSelf-0.25*m0rSelf[0]*uSelf[21]*mnuSelf+0.5*m1rSelf[21]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[18]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[17]*mnuSelf-0.25*m0rSelf[5]*uSelf[16]*mnuSelf; 
  data->AEM_S(56,22) = (-0.2*m0rSelf[3]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[22]*mnuSelf-0.25*m0rSelf[0]*uSelf[22]*mnuSelf+0.5*m1rSelf[22]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[20]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[20]*mnuSelf-0.2*m0rSelf[7]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[19]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.25*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(56,23) = (-0.159719141249985*m0rSelf[5]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[23]*mnuSelf-0.25*m0rSelf[0]*uSelf[23]*mnuSelf+0.5*m1rSelf[23]*mnuSelf-0.2*m0rSelf[3]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[21]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[20]*mnuSelf-0.2*m0rSelf[6]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[18]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[17]*mnuSelf-0.25*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(57,16) = (-0.2500000000000001*m0rSelf[5]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[22]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[19]*mnuSelf-0.25*m0rSelf[2]*uSelf[19]*mnuSelf-0.25*m0rSelf[3]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.25*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(57,17) = (-0.45*m0rSelf[7]*uSelf[23]*mnuSelf)-0.3928571428571428*m0rSelf[6]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[22]*mnuSelf-0.25*m0rSelf[5]*uSelf[21]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[20]*mnuSelf+0.4472135954999579*m1rSelf[20]*mnuSelf-0.45*m0rSelf[3]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[18]*mnuSelf-0.25*m0rSelf[2]*uSelf[18]*mnuSelf-0.45*m0rSelf[1]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(57,18) = (-0.2*m0rSelf[6]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[23]*mnuSelf-0.2*m0rSelf[7]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[19]*mnuSelf-0.25*m0rSelf[0]*uSelf[19]*mnuSelf+0.5*m1rSelf[19]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[18]*mnuSelf-0.25*m0rSelf[1]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[17]*mnuSelf-0.25*m0rSelf[2]*uSelf[17]*mnuSelf-0.25*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(57,19) = (-0.4024922359499621*m0rSelf[3]*uSelf[23]*mnuSelf)-0.2*m0rSelf[5]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[22]*mnuSelf+0.447213595499958*m1rSelf[22]*mnuSelf-0.2*m0rSelf[6]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[21]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[19]*mnuSelf-0.45*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[18]*mnuSelf-0.25*m0rSelf[0]*uSelf[18]*mnuSelf+0.5*m1rSelf[18]*mnuSelf-0.45*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[16]*mnuSelf-0.25*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(57,20) = (-0.223606797749979*m0rSelf[5]*uSelf[23]*mnuSelf)-0.3928571428571429*m0rSelf[3]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[21]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.4472135954999579*m1rSelf[17]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(57,21) = (-0.159719141249985*m0rSelf[5]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[23]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[23]*mnuSelf+0.5000000000000001*m1rSelf[23]*mnuSelf-0.2*m0rSelf[3]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[21]*mnuSelf-0.25*m0rSelf[1]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[20]*mnuSelf-0.2*m0rSelf[6]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[18]*mnuSelf-0.25*m0rSelf[5]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(57,22) = (-0.351382110749967*m0rSelf[6]*uSelf[23]*mnuSelf)-0.2*m0rSelf[2]*uSelf[23]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[22]*mnuSelf-0.2*m0rSelf[3]*uSelf[21]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[20]*mnuSelf-0.2*m0rSelf[5]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.447213595499958*m1rSelf[19]*mnuSelf-0.2*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(57,23) = (-0.2874944542499729*m0rSelf[7]*uSelf[23]*mnuSelf)-0.45*m0rSelf[1]*uSelf[23]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[22]*mnuSelf-0.2*m0rSelf[2]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[21]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[21]*mnuSelf+0.5000000000000001*m1rSelf[21]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[19]*mnuSelf-0.2*m0rSelf[6]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[18]*mnuSelf-0.45*m0rSelf[7]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[16]*mnuSelf; 
  data->AEM_S(58,16) = (-0.223606797749979*m0rSelf[3]*uSelf[23]*mnuSelf)-0.2500000000000001*m0rSelf[4]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[21]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[19]*mnuSelf-0.25*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[18]*mnuSelf-0.25*m0rSelf[0]*uSelf[18]*mnuSelf+0.5*m1rSelf[18]*mnuSelf-0.25*m0rSelf[3]*uSelf[17]*mnuSelf-0.25*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(58,17) = (-0.2*m0rSelf[6]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[23]*mnuSelf-0.2*m0rSelf[7]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[19]*mnuSelf-0.25*m0rSelf[0]*uSelf[19]*mnuSelf+0.5*m1rSelf[19]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[18]*mnuSelf-0.25*m0rSelf[1]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[17]*mnuSelf-0.25*m0rSelf[2]*uSelf[17]*mnuSelf-0.25*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(58,18) = (-0.3928571428571428*m0rSelf[7]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[23]*mnuSelf-0.45*m0rSelf[6]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[21]*mnuSelf+0.4472135954999579*m1rSelf[21]*mnuSelf-0.25*m0rSelf[4]*uSelf[20]*mnuSelf-0.45*m0rSelf[3]*uSelf[19]*mnuSelf-0.45*m0rSelf[2]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[17]*mnuSelf-0.25*m0rSelf[1]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(58,19) = (-0.3928571428571429*m0rSelf[5]*uSelf[23]*mnuSelf)-0.2*m0rSelf[4]*uSelf[23]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[23]*mnuSelf+0.447213595499958*m1rSelf[23]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[21]*mnuSelf-0.2*m0rSelf[7]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[19]*mnuSelf-0.45*m0rSelf[2]*uSelf[19]*mnuSelf-0.45*m0rSelf[3]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[16]*mnuSelf-0.25*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(58,20) = (-0.2*m0rSelf[3]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[22]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[22]*mnuSelf+0.5000000000000001*m1rSelf[22]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[20]*mnuSelf-0.25*m0rSelf[2]*uSelf[20]*mnuSelf-0.2*m0rSelf[7]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[19]*mnuSelf-0.25*m0rSelf[4]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(58,21) = (-0.3928571428571429*m0rSelf[3]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[19]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[18]*mnuSelf+0.4472135954999579*m1rSelf[18]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(58,22) = (-0.351382110749967*m0rSelf[7]*uSelf[23]*mnuSelf)-0.2*m0rSelf[1]*uSelf[23]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[22]*mnuSelf-0.45*m0rSelf[2]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[20]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[20]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[20]*mnuSelf+0.5000000000000001*m1rSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[19]*mnuSelf-0.45*m0rSelf[6]*uSelf[18]*mnuSelf-0.2*m0rSelf[7]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[16]*mnuSelf; 
  data->AEM_S(58,23) = (-0.351382110749967*m0rSelf[6]*uSelf[23]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[23]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[22]*mnuSelf-0.2*m0rSelf[1]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[21]*mnuSelf-0.2*m0rSelf[3]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[19]*mnuSelf-0.2*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.447213595499958*m1rSelf[19]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,16) = (-0.2*m0rSelf[6]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[23]*mnuSelf-0.2*m0rSelf[7]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[19]*mnuSelf-0.25*m0rSelf[0]*uSelf[19]*mnuSelf+0.5*m1rSelf[19]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[18]*mnuSelf-0.25*m0rSelf[1]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[17]*mnuSelf-0.25*m0rSelf[2]*uSelf[17]*mnuSelf-0.25*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,17) = (-0.4024922359499621*m0rSelf[3]*uSelf[23]*mnuSelf)-0.2*m0rSelf[5]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[22]*mnuSelf+0.447213595499958*m1rSelf[22]*mnuSelf-0.2*m0rSelf[6]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[21]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[19]*mnuSelf-0.45*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[18]*mnuSelf-0.25*m0rSelf[0]*uSelf[18]*mnuSelf+0.5*m1rSelf[18]*mnuSelf-0.45*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[16]*mnuSelf-0.25*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,18) = (-0.3928571428571429*m0rSelf[5]*uSelf[23]*mnuSelf)-0.2*m0rSelf[4]*uSelf[23]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[23]*mnuSelf+0.447213595499958*m1rSelf[23]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[21]*mnuSelf-0.2*m0rSelf[7]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[19]*mnuSelf-0.45*m0rSelf[2]*uSelf[19]*mnuSelf-0.45*m0rSelf[3]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[16]*mnuSelf-0.25*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,19) = (-0.7071428571428572*m0rSelf[7]*uSelf[23]*mnuSelf)-0.4024922359499621*m0rSelf[1]*uSelf[23]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[22]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[21]*mnuSelf-0.2*m0rSelf[4]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[21]*mnuSelf+0.4472135954999579*m1rSelf[21]*mnuSelf-0.2*m0rSelf[5]*uSelf[20]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[20]*mnuSelf+0.4472135954999579*m1rSelf[20]*mnuSelf-0.81*m0rSelf[3]*uSelf[19]*mnuSelf-0.4024922359499621*m0rSelf[6]*uSelf[18]*mnuSelf-0.45*m0rSelf[2]*uSelf[18]*mnuSelf-0.4024922359499621*m0rSelf[7]*uSelf[17]*mnuSelf-0.45*m0rSelf[1]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(59,20) = (-0.3513821107499669*m0rSelf[6]*uSelf[23]*mnuSelf)-0.2*m0rSelf[2]*uSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[22]*mnuSelf-0.2*m0rSelf[3]*uSelf[21]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[20]*mnuSelf-0.2*m0rSelf[5]*uSelf[19]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.4472135954999579*m1rSelf[19]*mnuSelf-0.2*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,21) = (-0.3513821107499669*m0rSelf[6]*uSelf[23]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[22]*mnuSelf-0.2*m0rSelf[1]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[21]*mnuSelf-0.2*m0rSelf[3]*uSelf[20]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[19]*mnuSelf-0.2*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.4472135954999579*m1rSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,22) = (-0.3513821107499669*m0rSelf[5]*uSelf[23]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[23]*mnuSelf-0.2*m0rSelf[0]*uSelf[23]*mnuSelf+0.4*m1rSelf[23]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[21]*mnuSelf-0.2*m0rSelf[1]*uSelf[21]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[20]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[19]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[19]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[18]*mnuSelf-0.2*m0rSelf[5]*uSelf[17]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.447213595499958*m1rSelf[17]*mnuSelf-0.2*m0rSelf[7]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(59,23) = (-0.7071428571428572*m0rSelf[3]*uSelf[23]*mnuSelf)-0.3513821107499669*m0rSelf[5]*uSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[22]*mnuSelf-0.2*m0rSelf[0]*uSelf[22]*mnuSelf+0.4*m1rSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[21]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[21]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[20]*mnuSelf-0.2*m0rSelf[2]*uSelf[20]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[19]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[18]*mnuSelf-0.2*m0rSelf[4]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[18]*mnuSelf+0.447213595499958*m1rSelf[18]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[17]*mnuSelf-0.2*m0rSelf[6]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(60,16) = (-0.223606797749979*m0rSelf[7]*uSelf[23]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[22]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[20]*mnuSelf-0.25*m0rSelf[0]*uSelf[20]*mnuSelf+0.5*m1rSelf[20]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[19]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf-0.25*m0rSelf[4]*uSelf[16]*mnuSelf; 
  data->AEM_S(60,17) = (-0.223606797749979*m0rSelf[5]*uSelf[23]*mnuSelf)-0.3928571428571429*m0rSelf[3]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[21]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.4472135954999579*m1rSelf[17]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(60,18) = (-0.2*m0rSelf[3]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[22]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[22]*mnuSelf+0.5000000000000001*m1rSelf[22]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[20]*mnuSelf-0.25*m0rSelf[2]*uSelf[20]*mnuSelf-0.2*m0rSelf[7]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[19]*mnuSelf-0.25*m0rSelf[4]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(60,19) = (-0.3513821107499669*m0rSelf[6]*uSelf[23]*mnuSelf)-0.2*m0rSelf[2]*uSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[22]*mnuSelf-0.2*m0rSelf[3]*uSelf[21]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[20]*mnuSelf-0.2*m0rSelf[5]*uSelf[19]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.4472135954999579*m1rSelf[19]*mnuSelf-0.2*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.3928571428571429*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(60,20) = (-0.3928571428571428*m0rSelf[7]*uSelf[23]*mnuSelf)-0.5357142857142857*m0rSelf[6]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[22]*mnuSelf-0.25*m0rSelf[5]*uSelf[21]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[20]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[20]*mnuSelf+0.31943828249997*m1rSelf[20]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[19]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[18]*mnuSelf-0.25*m0rSelf[2]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(60,21) = (-0.1428571428571428*m0rSelf[7]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[23]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[22]*mnuSelf-0.25*m0rSelf[4]*uSelf[21]*mnuSelf-0.25*m0rSelf[5]*uSelf[20]*mnuSelf-0.2*m0rSelf[3]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[17]*mnuSelf; 
  data->AEM_S(60,22) = (-0.3513821107499669*m0rSelf[3]*uSelf[23]*mnuSelf)-0.1428571428571428*m0rSelf[5]*uSelf[22]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[22]*mnuSelf+0.31943828249997*m1rSelf[22]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[21]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[20]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[18]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[18]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[18]*mnuSelf+0.5000000000000001*m1rSelf[18]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(60,23) = (-0.1428571428571428*m0rSelf[5]*uSelf[23]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[23]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[23]*mnuSelf+0.4472135954999579*m1rSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[22]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[21]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[19]*mnuSelf-0.2*m0rSelf[2]*uSelf[19]*mnuSelf-0.2*m0rSelf[3]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(61,16) = (-0.159719141249985*m0rSelf[7]*uSelf[23]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[23]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[21]*mnuSelf-0.25*m0rSelf[0]*uSelf[21]*mnuSelf+0.5*m1rSelf[21]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[18]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[17]*mnuSelf-0.25*m0rSelf[5]*uSelf[16]*mnuSelf; 
  data->AEM_S(61,17) = (-0.159719141249985*m0rSelf[5]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[23]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[23]*mnuSelf+0.5000000000000001*m1rSelf[23]*mnuSelf-0.2*m0rSelf[3]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[21]*mnuSelf-0.25*m0rSelf[1]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[20]*mnuSelf-0.2*m0rSelf[6]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[18]*mnuSelf-0.25*m0rSelf[5]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(61,18) = (-0.3928571428571429*m0rSelf[3]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[19]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[18]*mnuSelf+0.4472135954999579*m1rSelf[18]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(61,19) = (-0.3513821107499669*m0rSelf[6]*uSelf[23]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[22]*mnuSelf-0.2*m0rSelf[1]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[21]*mnuSelf-0.2*m0rSelf[3]*uSelf[20]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[19]*mnuSelf-0.2*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.4472135954999579*m1rSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(61,20) = (-0.1428571428571428*m0rSelf[7]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[23]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[22]*mnuSelf-0.25*m0rSelf[4]*uSelf[21]*mnuSelf-0.25*m0rSelf[5]*uSelf[20]*mnuSelf-0.2*m0rSelf[3]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[17]*mnuSelf; 
  data->AEM_S(61,21) = (-0.5357142857142857*m0rSelf[7]*uSelf[23]*mnuSelf)-0.159719141249985*m0rSelf[1]*uSelf[23]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[22]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[21]*mnuSelf+0.31943828249997*m1rSelf[21]*mnuSelf-0.25*m0rSelf[4]*uSelf[20]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[19]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[18]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[17]*mnuSelf-0.25*m0rSelf[1]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(61,22) = (-0.3513821107499669*m0rSelf[3]*uSelf[23]*mnuSelf)-0.3928571428571428*m0rSelf[5]*uSelf[22]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[22]*mnuSelf+0.4472135954999579*m1rSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[19]*mnuSelf-0.2*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[18]*mnuSelf-0.2*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(61,23) = (-0.5357142857142857*m0rSelf[5]*uSelf[23]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[23]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[23]*mnuSelf+0.31943828249997*m1rSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[22]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[18]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[17]*mnuSelf+0.5000000000000001*m1rSelf[17]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(62,16) = (-0.2*m0rSelf[3]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[5]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[22]*mnuSelf-0.25*m0rSelf[0]*uSelf[22]*mnuSelf+0.5*m1rSelf[22]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[20]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[20]*mnuSelf-0.2*m0rSelf[7]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[19]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf-0.25*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(62,17) = (-0.351382110749967*m0rSelf[6]*uSelf[23]*mnuSelf)-0.2*m0rSelf[2]*uSelf[23]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[22]*mnuSelf-0.2*m0rSelf[3]*uSelf[21]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[20]*mnuSelf-0.2*m0rSelf[5]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.447213595499958*m1rSelf[19]*mnuSelf-0.2*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(62,18) = (-0.351382110749967*m0rSelf[7]*uSelf[23]*mnuSelf)-0.2*m0rSelf[1]*uSelf[23]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[22]*mnuSelf-0.45*m0rSelf[2]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[20]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[20]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[20]*mnuSelf+0.5000000000000001*m1rSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[19]*mnuSelf-0.45*m0rSelf[6]*uSelf[18]*mnuSelf-0.2*m0rSelf[7]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[16]*mnuSelf; 
  data->AEM_S(62,19) = (-0.3513821107499669*m0rSelf[5]*uSelf[23]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[23]*mnuSelf-0.2*m0rSelf[0]*uSelf[23]*mnuSelf+0.4*m1rSelf[23]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[21]*mnuSelf-0.2*m0rSelf[1]*uSelf[21]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[20]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[19]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[19]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[18]*mnuSelf-0.2*m0rSelf[5]*uSelf[17]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.447213595499958*m1rSelf[17]*mnuSelf-0.2*m0rSelf[7]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(62,20) = (-0.3513821107499669*m0rSelf[3]*uSelf[23]*mnuSelf)-0.1428571428571428*m0rSelf[5]*uSelf[22]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[22]*mnuSelf+0.31943828249997*m1rSelf[22]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[21]*mnuSelf-0.5357142857142857*m0rSelf[6]*uSelf[20]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[18]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[18]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[18]*mnuSelf+0.5000000000000001*m1rSelf[18]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[6]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(62,21) = (-0.3513821107499669*m0rSelf[3]*uSelf[23]*mnuSelf)-0.3928571428571428*m0rSelf[5]*uSelf[22]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[22]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[22]*mnuSelf+0.4472135954999579*m1rSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[7]*uSelf[19]*mnuSelf-0.2*m0rSelf[1]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[18]*mnuSelf-0.2*m0rSelf[3]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[16]*mnuSelf; 
  data->AEM_S(62,22) = (-0.6173469387755102*m0rSelf[7]*uSelf[23]*mnuSelf)-0.351382110749967*m0rSelf[1]*uSelf[23]*mnuSelf-0.9642857142857143*m0rSelf[6]*uSelf[22]*mnuSelf-0.2874944542499729*m0rSelf[2]*uSelf[22]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[21]*mnuSelf+0.4472135954999579*m1rSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[20]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[20]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[20]*mnuSelf+0.31943828249997*m1rSelf[20]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[19]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[18]*mnuSelf-0.45*m0rSelf[2]*uSelf[18]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[17]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(62,23) = (-0.6173469387755102*m0rSelf[6]*uSelf[23]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[23]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[22]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[21]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[19]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[19]*mnuSelf-0.2*m0rSelf[0]*uSelf[19]*mnuSelf+0.4*m1rSelf[19]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[18]*mnuSelf-0.2*m0rSelf[1]*uSelf[18]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[17]*mnuSelf-0.2*m0rSelf[2]*uSelf[17]*mnuSelf-0.2*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,16) = (-0.159719141249985*m0rSelf[5]*uSelf[23]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[23]*mnuSelf-0.25*m0rSelf[0]*uSelf[23]*mnuSelf+0.5*m1rSelf[23]*mnuSelf-0.2*m0rSelf[3]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[21]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[20]*mnuSelf-0.2*m0rSelf[6]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[18]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[17]*mnuSelf-0.25*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,17) = (-0.2874944542499729*m0rSelf[7]*uSelf[23]*mnuSelf)-0.45*m0rSelf[1]*uSelf[23]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[22]*mnuSelf-0.2*m0rSelf[2]*uSelf[22]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[21]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[21]*mnuSelf+0.5000000000000001*m1rSelf[21]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[20]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[19]*mnuSelf-0.2*m0rSelf[6]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[18]*mnuSelf-0.45*m0rSelf[7]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[5]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,18) = (-0.351382110749967*m0rSelf[6]*uSelf[23]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[23]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[22]*mnuSelf-0.2*m0rSelf[1]*uSelf[22]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[21]*mnuSelf-0.2*m0rSelf[3]*uSelf[20]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[19]*mnuSelf-0.2*m0rSelf[4]*uSelf[19]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[19]*mnuSelf+0.447213595499958*m1rSelf[19]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[18]*mnuSelf-0.2*m0rSelf[6]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,19) = (-0.7071428571428572*m0rSelf[3]*uSelf[23]*mnuSelf)-0.3513821107499669*m0rSelf[5]*uSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[22]*mnuSelf-0.2*m0rSelf[0]*uSelf[22]*mnuSelf+0.4*m1rSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[21]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[21]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[20]*mnuSelf-0.2*m0rSelf[2]*uSelf[20]*mnuSelf-0.7071428571428572*m0rSelf[7]*uSelf[19]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[5]*uSelf[18]*mnuSelf-0.2*m0rSelf[4]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[18]*mnuSelf+0.447213595499958*m1rSelf[18]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[17]*mnuSelf-0.2*m0rSelf[6]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,20) = (-0.1428571428571428*m0rSelf[5]*uSelf[23]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[23]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[23]*mnuSelf+0.4472135954999579*m1rSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[22]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[21]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[21]*mnuSelf-0.3928571428571428*m0rSelf[7]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[19]*mnuSelf-0.2*m0rSelf[2]*uSelf[19]*mnuSelf-0.2*m0rSelf[3]*uSelf[18]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[7]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,21) = (-0.5357142857142857*m0rSelf[5]*uSelf[23]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[23]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[23]*mnuSelf+0.31943828249997*m1rSelf[23]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[22]*mnuSelf-0.5357142857142857*m0rSelf[7]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[7]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[6]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[19]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[18]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[17]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[17]*mnuSelf+0.5000000000000001*m1rSelf[17]*mnuSelf-0.159719141249985*m0rSelf[7]*uSelf[16]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,22) = (-0.6173469387755102*m0rSelf[6]*uSelf[23]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[23]*mnuSelf-0.6173469387755102*m0rSelf[7]*uSelf[22]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[22]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[21]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[20]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[19]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[19]*mnuSelf-0.2*m0rSelf[0]*uSelf[19]*mnuSelf+0.4*m1rSelf[19]*mnuSelf-0.351382110749967*m0rSelf[7]*uSelf[18]*mnuSelf-0.2*m0rSelf[1]*uSelf[18]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[17]*mnuSelf-0.2*m0rSelf[2]*uSelf[17]*mnuSelf-0.2*m0rSelf[3]*uSelf[16]*mnuSelf; 
  data->AEM_S(63,23) = (-0.9642857142857143*m0rSelf[7]*uSelf[23]*mnuSelf)-0.2874944542499729*m0rSelf[1]*uSelf[23]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[22]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[22]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[21]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[21]*mnuSelf+0.31943828249997*m1rSelf[21]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[20]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[20]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[20]*mnuSelf+0.4472135954999579*m1rSelf[20]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[19]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[18]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[18]*mnuSelf-0.2874944542499729*m0rSelf[7]*uSelf[17]*mnuSelf-0.45*m0rSelf[1]*uSelf[17]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherZ-uOtherZ*m0Other) and uCrossOtherZ ... // 
  data->AEM_S(56,48) = 0.25*m0rOther[7]*uOther[23]*mnuOther+0.25*m0rOther[6]*uOther[22]*mnuOther+0.25*m0rOther[5]*uOther[21]*mnuOther+0.25*m0rOther[4]*uOther[20]*mnuOther+0.25*m0rOther[3]*uOther[19]*mnuOther+0.25*m0rOther[2]*uOther[18]*mnuOther+0.25*m0rOther[1]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(56,49) = 0.2500000000000001*m0rOther[5]*uOther[23]*mnuOther+0.223606797749979*m0rOther[3]*uOther[22]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[21]*mnuOther+0.223606797749979*m0rOther[1]*uOther[20]*mnuOther+0.223606797749979*m0rOther[6]*uOther[19]*mnuOther+0.25*m0rOther[2]*uOther[19]*mnuOther+0.25*m0rOther[3]*uOther[18]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.25*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(56,50) = 0.223606797749979*m0rOther[3]*uOther[23]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[22]*mnuOther+0.223606797749979*m0rOther[2]*uOther[21]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[20]*mnuOther+0.223606797749979*m0rOther[7]*uOther[19]*mnuOther+0.25*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[5]*uOther[18]*mnuOther+0.25*m0rOther[0]*uOther[18]*mnuOther-0.5*m1rOther[18]*mnuOther+0.25*m0rOther[3]*uOther[17]*mnuOther+0.25*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(56,51) = 0.2*m0rOther[6]*uOther[23]*mnuOther+0.223606797749979*m0rOther[2]*uOther[23]*mnuOther+0.2*m0rOther[7]*uOther[22]*mnuOther+0.223606797749979*m0rOther[1]*uOther[22]*mnuOther+0.223606797749979*m0rOther[3]*uOther[21]*mnuOther+0.223606797749979*m0rOther[3]*uOther[20]*mnuOther+0.223606797749979*m0rOther[5]*uOther[19]*mnuOther+0.223606797749979*m0rOther[4]*uOther[19]*mnuOther+0.25*m0rOther[0]*uOther[19]*mnuOther-0.5*m1rOther[19]*mnuOther+0.223606797749979*m0rOther[7]*uOther[18]*mnuOther+0.25*m0rOther[1]*uOther[18]*mnuOther+0.223606797749979*m0rOther[6]*uOther[17]*mnuOther+0.25*m0rOther[2]*uOther[17]*mnuOther+0.25*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(56,52) = 0.223606797749979*m0rOther[7]*uOther[23]*mnuOther+0.159719141249985*m0rOther[6]*uOther[22]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[22]*mnuOther+0.159719141249985*m0rOther[4]*uOther[20]*mnuOther+0.25*m0rOther[0]*uOther[20]*mnuOther-0.5*m1rOther[20]*mnuOther+0.223606797749979*m0rOther[3]*uOther[19]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.25*m0rOther[4]*uOther[16]*mnuOther; 
  data->AEM_S(56,53) = 0.159719141249985*m0rOther[7]*uOther[23]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[23]*mnuOther+0.223606797749979*m0rOther[6]*uOther[22]*mnuOther+0.159719141249985*m0rOther[5]*uOther[21]*mnuOther+0.25*m0rOther[0]*uOther[21]*mnuOther-0.5*m1rOther[21]*mnuOther+0.223606797749979*m0rOther[3]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[17]*mnuOther+0.25*m0rOther[5]*uOther[16]*mnuOther; 
  data->AEM_S(56,54) = 0.2*m0rOther[3]*uOther[23]*mnuOther+0.223606797749979*m0rOther[5]*uOther[22]*mnuOther+0.159719141249985*m0rOther[4]*uOther[22]*mnuOther+0.25*m0rOther[0]*uOther[22]*mnuOther-0.5*m1rOther[22]*mnuOther+0.223606797749979*m0rOther[6]*uOther[21]*mnuOther+0.159719141249985*m0rOther[6]*uOther[20]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[20]*mnuOther+0.2*m0rOther[7]*uOther[19]*mnuOther+0.223606797749979*m0rOther[1]*uOther[19]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[18]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.25*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(56,55) = 0.159719141249985*m0rOther[5]*uOther[23]*mnuOther+0.223606797749979*m0rOther[4]*uOther[23]*mnuOther+0.25*m0rOther[0]*uOther[23]*mnuOther-0.5*m1rOther[23]*mnuOther+0.2*m0rOther[3]*uOther[22]*mnuOther+0.159719141249985*m0rOther[7]*uOther[21]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[21]*mnuOther+0.223606797749979*m0rOther[7]*uOther[20]*mnuOther+0.2*m0rOther[6]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[17]*mnuOther+0.25*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(57,48) = 0.2500000000000001*m0rOther[5]*uOther[23]*mnuOther+0.223606797749979*m0rOther[3]*uOther[22]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[21]*mnuOther+0.223606797749979*m0rOther[1]*uOther[20]*mnuOther+0.223606797749979*m0rOther[6]*uOther[19]*mnuOther+0.25*m0rOther[2]*uOther[19]*mnuOther+0.25*m0rOther[3]*uOther[18]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.25*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(57,49) = 0.45*m0rOther[7]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[22]*mnuOther+0.223606797749979*m0rOther[2]*uOther[22]*mnuOther+0.25*m0rOther[5]*uOther[21]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[20]*mnuOther+0.223606797749979*m0rOther[0]*uOther[20]*mnuOther-0.4472135954999579*m1rOther[20]*mnuOther+0.45*m0rOther[3]*uOther[19]*mnuOther+0.223606797749979*m0rOther[6]*uOther[18]*mnuOther+0.25*m0rOther[2]*uOther[18]*mnuOther+0.45*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(57,50) = 0.2*m0rOther[6]*uOther[23]*mnuOther+0.223606797749979*m0rOther[2]*uOther[23]*mnuOther+0.2*m0rOther[7]*uOther[22]*mnuOther+0.223606797749979*m0rOther[1]*uOther[22]*mnuOther+0.223606797749979*m0rOther[3]*uOther[21]*mnuOther+0.223606797749979*m0rOther[3]*uOther[20]*mnuOther+0.223606797749979*m0rOther[5]*uOther[19]*mnuOther+0.223606797749979*m0rOther[4]*uOther[19]*mnuOther+0.25*m0rOther[0]*uOther[19]*mnuOther-0.5*m1rOther[19]*mnuOther+0.223606797749979*m0rOther[7]*uOther[18]*mnuOther+0.25*m0rOther[1]*uOther[18]*mnuOther+0.223606797749979*m0rOther[6]*uOther[17]*mnuOther+0.25*m0rOther[2]*uOther[17]*mnuOther+0.25*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(57,51) = 0.4024922359499621*m0rOther[3]*uOther[23]*mnuOther+0.2*m0rOther[5]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[22]*mnuOther+0.223606797749979*m0rOther[0]*uOther[22]*mnuOther-0.447213595499958*m1rOther[22]*mnuOther+0.2*m0rOther[6]*uOther[21]*mnuOther+0.223606797749979*m0rOther[2]*uOther[21]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[20]*mnuOther+0.223606797749979*m0rOther[2]*uOther[20]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[19]*mnuOther+0.45*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[5]*uOther[18]*mnuOther+0.223606797749979*m0rOther[4]*uOther[18]*mnuOther+0.25*m0rOther[0]*uOther[18]*mnuOther-0.5*m1rOther[18]*mnuOther+0.45*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[6]*uOther[16]*mnuOther+0.25*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(57,52) = 0.223606797749979*m0rOther[5]*uOther[23]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[22]*mnuOther+0.223606797749979*m0rOther[7]*uOther[21]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.4472135954999579*m1rOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(57,53) = 0.159719141249985*m0rOther[5]*uOther[23]*mnuOther+0.223606797749979*m0rOther[4]*uOther[23]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[23]*mnuOther-0.5000000000000001*m1rOther[23]*mnuOther+0.2*m0rOther[3]*uOther[22]*mnuOther+0.159719141249985*m0rOther[7]*uOther[21]*mnuOther+0.25*m0rOther[1]*uOther[21]*mnuOther+0.223606797749979*m0rOther[7]*uOther[20]*mnuOther+0.2*m0rOther[6]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[18]*mnuOther+0.25*m0rOther[5]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(57,54) = 0.351382110749967*m0rOther[6]*uOther[23]*mnuOther+0.2*m0rOther[2]*uOther[23]*mnuOther+0.351382110749967*m0rOther[7]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[22]*mnuOther+0.2*m0rOther[3]*uOther[21]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[20]*mnuOther+0.2*m0rOther[5]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.447213595499958*m1rOther[19]*mnuOther+0.2*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(57,55) = 0.2874944542499729*m0rOther[7]*uOther[23]*mnuOther+0.45*m0rOther[1]*uOther[23]*mnuOther+0.351382110749967*m0rOther[6]*uOther[22]*mnuOther+0.2*m0rOther[2]*uOther[22]*mnuOther+0.159719141249985*m0rOther[5]*uOther[21]*mnuOther+0.223606797749979*m0rOther[4]*uOther[21]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[21]*mnuOther-0.5000000000000001*m1rOther[21]*mnuOther+0.223606797749979*m0rOther[5]*uOther[20]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[19]*mnuOther+0.2*m0rOther[6]*uOther[18]*mnuOther+0.223606797749979*m0rOther[2]*uOther[18]*mnuOther+0.45*m0rOther[7]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[16]*mnuOther; 
  data->AEM_S(58,48) = 0.223606797749979*m0rOther[3]*uOther[23]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[22]*mnuOther+0.223606797749979*m0rOther[2]*uOther[21]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[20]*mnuOther+0.223606797749979*m0rOther[7]*uOther[19]*mnuOther+0.25*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[5]*uOther[18]*mnuOther+0.25*m0rOther[0]*uOther[18]*mnuOther-0.5*m1rOther[18]*mnuOther+0.25*m0rOther[3]*uOther[17]*mnuOther+0.25*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(58,49) = 0.2*m0rOther[6]*uOther[23]*mnuOther+0.223606797749979*m0rOther[2]*uOther[23]*mnuOther+0.2*m0rOther[7]*uOther[22]*mnuOther+0.223606797749979*m0rOther[1]*uOther[22]*mnuOther+0.223606797749979*m0rOther[3]*uOther[21]*mnuOther+0.223606797749979*m0rOther[3]*uOther[20]*mnuOther+0.223606797749979*m0rOther[5]*uOther[19]*mnuOther+0.223606797749979*m0rOther[4]*uOther[19]*mnuOther+0.25*m0rOther[0]*uOther[19]*mnuOther-0.5*m1rOther[19]*mnuOther+0.223606797749979*m0rOther[7]*uOther[18]*mnuOther+0.25*m0rOther[1]*uOther[18]*mnuOther+0.223606797749979*m0rOther[6]*uOther[17]*mnuOther+0.25*m0rOther[2]*uOther[17]*mnuOther+0.25*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(58,50) = 0.3928571428571428*m0rOther[7]*uOther[23]*mnuOther+0.223606797749979*m0rOther[1]*uOther[23]*mnuOther+0.45*m0rOther[6]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[21]*mnuOther+0.223606797749979*m0rOther[0]*uOther[21]*mnuOther-0.4472135954999579*m1rOther[21]*mnuOther+0.25*m0rOther[4]*uOther[20]*mnuOther+0.45*m0rOther[3]*uOther[19]*mnuOther+0.45*m0rOther[2]*uOther[18]*mnuOther+0.223606797749979*m0rOther[7]*uOther[17]*mnuOther+0.25*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[5]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(58,51) = 0.3928571428571429*m0rOther[5]*uOther[23]*mnuOther+0.2*m0rOther[4]*uOther[23]*mnuOther+0.223606797749979*m0rOther[0]*uOther[23]*mnuOther-0.447213595499958*m1rOther[23]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[21]*mnuOther+0.223606797749979*m0rOther[1]*uOther[21]*mnuOther+0.2*m0rOther[7]*uOther[20]*mnuOther+0.223606797749979*m0rOther[1]*uOther[20]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[19]*mnuOther+0.45*m0rOther[2]*uOther[19]*mnuOther+0.45*m0rOther[3]*uOther[18]*mnuOther+0.223606797749979*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.223606797749979*m0rOther[7]*uOther[16]*mnuOther+0.25*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(58,52) = 0.2*m0rOther[3]*uOther[23]*mnuOther+0.223606797749979*m0rOther[5]*uOther[22]*mnuOther+0.159719141249985*m0rOther[4]*uOther[22]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[22]*mnuOther-0.5000000000000001*m1rOther[22]*mnuOther+0.223606797749979*m0rOther[6]*uOther[21]*mnuOther+0.159719141249985*m0rOther[6]*uOther[20]*mnuOther+0.25*m0rOther[2]*uOther[20]*mnuOther+0.2*m0rOther[7]*uOther[19]*mnuOther+0.223606797749979*m0rOther[1]*uOther[19]*mnuOther+0.25*m0rOther[4]*uOther[18]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(58,53) = 0.3928571428571429*m0rOther[3]*uOther[23]*mnuOther+0.223606797749979*m0rOther[4]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[21]*mnuOther+0.223606797749979*m0rOther[6]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[19]*mnuOther+0.223606797749979*m0rOther[1]*uOther[19]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[18]*mnuOther+0.223606797749979*m0rOther[0]*uOther[18]*mnuOther-0.4472135954999579*m1rOther[18]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(58,54) = 0.351382110749967*m0rOther[7]*uOther[23]*mnuOther+0.2*m0rOther[1]*uOther[23]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[22]*mnuOther+0.45*m0rOther[2]*uOther[22]*mnuOther+0.223606797749979*m0rOther[4]*uOther[21]*mnuOther+0.223606797749979*m0rOther[5]*uOther[20]*mnuOther+0.159719141249985*m0rOther[4]*uOther[20]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[20]*mnuOther-0.5000000000000001*m1rOther[20]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[19]*mnuOther+0.45*m0rOther[6]*uOther[18]*mnuOther+0.2*m0rOther[7]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[16]*mnuOther; 
  data->AEM_S(58,55) = 0.351382110749967*m0rOther[6]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[23]*mnuOther+0.351382110749967*m0rOther[7]*uOther[22]*mnuOther+0.2*m0rOther[1]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[21]*mnuOther+0.2*m0rOther[3]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[19]*mnuOther+0.2*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.447213595499958*m1rOther[19]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(59,48) = 0.2*m0rOther[6]*uOther[23]*mnuOther+0.223606797749979*m0rOther[2]*uOther[23]*mnuOther+0.2*m0rOther[7]*uOther[22]*mnuOther+0.223606797749979*m0rOther[1]*uOther[22]*mnuOther+0.223606797749979*m0rOther[3]*uOther[21]*mnuOther+0.223606797749979*m0rOther[3]*uOther[20]*mnuOther+0.223606797749979*m0rOther[5]*uOther[19]*mnuOther+0.223606797749979*m0rOther[4]*uOther[19]*mnuOther+0.25*m0rOther[0]*uOther[19]*mnuOther-0.5*m1rOther[19]*mnuOther+0.223606797749979*m0rOther[7]*uOther[18]*mnuOther+0.25*m0rOther[1]*uOther[18]*mnuOther+0.223606797749979*m0rOther[6]*uOther[17]*mnuOther+0.25*m0rOther[2]*uOther[17]*mnuOther+0.25*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(59,49) = 0.4024922359499621*m0rOther[3]*uOther[23]*mnuOther+0.2*m0rOther[5]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[22]*mnuOther+0.223606797749979*m0rOther[0]*uOther[22]*mnuOther-0.447213595499958*m1rOther[22]*mnuOther+0.2*m0rOther[6]*uOther[21]*mnuOther+0.223606797749979*m0rOther[2]*uOther[21]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[20]*mnuOther+0.223606797749979*m0rOther[2]*uOther[20]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[19]*mnuOther+0.45*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[5]*uOther[18]*mnuOther+0.223606797749979*m0rOther[4]*uOther[18]*mnuOther+0.25*m0rOther[0]*uOther[18]*mnuOther-0.5*m1rOther[18]*mnuOther+0.45*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[6]*uOther[16]*mnuOther+0.25*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(59,50) = 0.3928571428571429*m0rOther[5]*uOther[23]*mnuOther+0.2*m0rOther[4]*uOther[23]*mnuOther+0.223606797749979*m0rOther[0]*uOther[23]*mnuOther-0.447213595499958*m1rOther[23]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[21]*mnuOther+0.223606797749979*m0rOther[1]*uOther[21]*mnuOther+0.2*m0rOther[7]*uOther[20]*mnuOther+0.223606797749979*m0rOther[1]*uOther[20]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[19]*mnuOther+0.45*m0rOther[2]*uOther[19]*mnuOther+0.45*m0rOther[3]*uOther[18]*mnuOther+0.223606797749979*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.223606797749979*m0rOther[7]*uOther[16]*mnuOther+0.25*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(59,51) = 0.7071428571428572*m0rOther[7]*uOther[23]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[23]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[22]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[21]*mnuOther+0.2*m0rOther[4]*uOther[21]*mnuOther+0.223606797749979*m0rOther[0]*uOther[21]*mnuOther-0.4472135954999579*m1rOther[21]*mnuOther+0.2*m0rOther[5]*uOther[20]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[20]*mnuOther+0.223606797749979*m0rOther[0]*uOther[20]*mnuOther-0.4472135954999579*m1rOther[20]*mnuOther+0.81*m0rOther[3]*uOther[19]*mnuOther+0.4024922359499621*m0rOther[6]*uOther[18]*mnuOther+0.45*m0rOther[2]*uOther[18]*mnuOther+0.4024922359499621*m0rOther[7]*uOther[17]*mnuOther+0.45*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[5]*uOther[16]*mnuOther+0.223606797749979*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(59,52) = 0.3513821107499669*m0rOther[6]*uOther[23]*mnuOther+0.2*m0rOther[2]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[22]*mnuOther+0.2*m0rOther[3]*uOther[21]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[20]*mnuOther+0.2*m0rOther[5]*uOther[19]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.4472135954999579*m1rOther[19]*mnuOther+0.2*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(59,53) = 0.3513821107499669*m0rOther[6]*uOther[23]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[22]*mnuOther+0.2*m0rOther[1]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[21]*mnuOther+0.2*m0rOther[3]*uOther[20]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[19]*mnuOther+0.2*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.4472135954999579*m1rOther[19]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(59,54) = 0.3513821107499669*m0rOther[5]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[23]*mnuOther+0.2*m0rOther[0]*uOther[23]*mnuOther-0.4*m1rOther[23]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[22]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[21]*mnuOther+0.2*m0rOther[1]*uOther[21]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[20]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[19]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[19]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[18]*mnuOther+0.2*m0rOther[5]*uOther[17]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.447213595499958*m1rOther[17]*mnuOther+0.2*m0rOther[7]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(59,55) = 0.7071428571428572*m0rOther[3]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[22]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[22]*mnuOther+0.2*m0rOther[0]*uOther[22]*mnuOther-0.4*m1rOther[22]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[21]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[21]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[20]*mnuOther+0.2*m0rOther[2]*uOther[20]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[19]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[18]*mnuOther+0.2*m0rOther[4]*uOther[18]*mnuOther+0.223606797749979*m0rOther[0]*uOther[18]*mnuOther-0.447213595499958*m1rOther[18]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[17]*mnuOther+0.2*m0rOther[6]*uOther[16]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(60,48) = 0.223606797749979*m0rOther[7]*uOther[23]*mnuOther+0.159719141249985*m0rOther[6]*uOther[22]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[22]*mnuOther+0.159719141249985*m0rOther[4]*uOther[20]*mnuOther+0.25*m0rOther[0]*uOther[20]*mnuOther-0.5*m1rOther[20]*mnuOther+0.223606797749979*m0rOther[3]*uOther[19]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.25*m0rOther[4]*uOther[16]*mnuOther; 
  data->AEM_S(60,49) = 0.223606797749979*m0rOther[5]*uOther[23]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[22]*mnuOther+0.223606797749979*m0rOther[7]*uOther[21]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.4472135954999579*m1rOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(60,50) = 0.2*m0rOther[3]*uOther[23]*mnuOther+0.223606797749979*m0rOther[5]*uOther[22]*mnuOther+0.159719141249985*m0rOther[4]*uOther[22]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[22]*mnuOther-0.5000000000000001*m1rOther[22]*mnuOther+0.223606797749979*m0rOther[6]*uOther[21]*mnuOther+0.159719141249985*m0rOther[6]*uOther[20]*mnuOther+0.25*m0rOther[2]*uOther[20]*mnuOther+0.2*m0rOther[7]*uOther[19]*mnuOther+0.223606797749979*m0rOther[1]*uOther[19]*mnuOther+0.25*m0rOther[4]*uOther[18]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(60,51) = 0.3513821107499669*m0rOther[6]*uOther[23]*mnuOther+0.2*m0rOther[2]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[22]*mnuOther+0.2*m0rOther[3]*uOther[21]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[20]*mnuOther+0.2*m0rOther[5]*uOther[19]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.4472135954999579*m1rOther[19]*mnuOther+0.2*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.3928571428571429*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(60,52) = 0.3928571428571428*m0rOther[7]*uOther[23]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[22]*mnuOther+0.159719141249985*m0rOther[2]*uOther[22]*mnuOther+0.25*m0rOther[5]*uOther[21]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[20]*mnuOther+0.159719141249985*m0rOther[0]*uOther[20]*mnuOther-0.31943828249997*m1rOther[20]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[19]*mnuOther+0.159719141249985*m0rOther[6]*uOther[18]*mnuOther+0.25*m0rOther[2]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[17]*mnuOther+0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(60,53) = 0.1428571428571428*m0rOther[7]*uOther[23]*mnuOther+0.223606797749979*m0rOther[1]*uOther[23]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[22]*mnuOther+0.223606797749979*m0rOther[2]*uOther[22]*mnuOther+0.25*m0rOther[4]*uOther[21]*mnuOther+0.25*m0rOther[5]*uOther[20]*mnuOther+0.2*m0rOther[3]*uOther[19]*mnuOther+0.223606797749979*m0rOther[6]*uOther[18]*mnuOther+0.223606797749979*m0rOther[7]*uOther[17]*mnuOther; 
  data->AEM_S(60,54) = 0.3513821107499669*m0rOther[3]*uOther[23]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[22]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[22]*mnuOther+0.159719141249985*m0rOther[0]*uOther[22]*mnuOther-0.31943828249997*m1rOther[22]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[21]*mnuOther+0.223606797749979*m0rOther[2]*uOther[21]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[20]*mnuOther+0.159719141249985*m0rOther[2]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[5]*uOther[18]*mnuOther+0.159719141249985*m0rOther[4]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[18]*mnuOther-0.5000000000000001*m1rOther[18]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[17]*mnuOther+0.159719141249985*m0rOther[6]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(60,55) = 0.1428571428571428*m0rOther[5]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[23]*mnuOther+0.223606797749979*m0rOther[0]*uOther[23]*mnuOther-0.4472135954999579*m1rOther[23]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[22]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[21]*mnuOther+0.223606797749979*m0rOther[1]*uOther[21]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[19]*mnuOther+0.2*m0rOther[2]*uOther[19]*mnuOther+0.2*m0rOther[3]*uOther[18]*mnuOther+0.223606797749979*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(61,48) = 0.159719141249985*m0rOther[7]*uOther[23]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[23]*mnuOther+0.223606797749979*m0rOther[6]*uOther[22]*mnuOther+0.159719141249985*m0rOther[5]*uOther[21]*mnuOther+0.25*m0rOther[0]*uOther[21]*mnuOther-0.5*m1rOther[21]*mnuOther+0.223606797749979*m0rOther[3]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[17]*mnuOther+0.25*m0rOther[5]*uOther[16]*mnuOther; 
  data->AEM_S(61,49) = 0.159719141249985*m0rOther[5]*uOther[23]*mnuOther+0.223606797749979*m0rOther[4]*uOther[23]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[23]*mnuOther-0.5000000000000001*m1rOther[23]*mnuOther+0.2*m0rOther[3]*uOther[22]*mnuOther+0.159719141249985*m0rOther[7]*uOther[21]*mnuOther+0.25*m0rOther[1]*uOther[21]*mnuOther+0.223606797749979*m0rOther[7]*uOther[20]*mnuOther+0.2*m0rOther[6]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[18]*mnuOther+0.25*m0rOther[5]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(61,50) = 0.3928571428571429*m0rOther[3]*uOther[23]*mnuOther+0.223606797749979*m0rOther[4]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[21]*mnuOther+0.223606797749979*m0rOther[6]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[19]*mnuOther+0.223606797749979*m0rOther[1]*uOther[19]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[18]*mnuOther+0.223606797749979*m0rOther[0]*uOther[18]*mnuOther-0.4472135954999579*m1rOther[18]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(61,51) = 0.3513821107499669*m0rOther[6]*uOther[23]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[22]*mnuOther+0.2*m0rOther[1]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[21]*mnuOther+0.2*m0rOther[3]*uOther[20]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[19]*mnuOther+0.2*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.4472135954999579*m1rOther[19]*mnuOther+0.3928571428571429*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(61,52) = 0.1428571428571428*m0rOther[7]*uOther[23]*mnuOther+0.223606797749979*m0rOther[1]*uOther[23]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[22]*mnuOther+0.223606797749979*m0rOther[2]*uOther[22]*mnuOther+0.25*m0rOther[4]*uOther[21]*mnuOther+0.25*m0rOther[5]*uOther[20]*mnuOther+0.2*m0rOther[3]*uOther[19]*mnuOther+0.223606797749979*m0rOther[6]*uOther[18]*mnuOther+0.223606797749979*m0rOther[7]*uOther[17]*mnuOther; 
  data->AEM_S(61,53) = 0.5357142857142857*m0rOther[7]*uOther[23]*mnuOther+0.159719141249985*m0rOther[1]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[22]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[21]*mnuOther+0.159719141249985*m0rOther[0]*uOther[21]*mnuOther-0.31943828249997*m1rOther[21]*mnuOther+0.25*m0rOther[4]*uOther[20]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[19]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[18]*mnuOther+0.159719141249985*m0rOther[7]*uOther[17]*mnuOther+0.25*m0rOther[1]*uOther[17]*mnuOther+0.159719141249985*m0rOther[5]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(61,54) = 0.3513821107499669*m0rOther[3]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[22]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[22]*mnuOther+0.223606797749979*m0rOther[0]*uOther[22]*mnuOther-0.4472135954999579*m1rOther[22]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[21]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[20]*mnuOther+0.223606797749979*m0rOther[2]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[19]*mnuOther+0.2*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[4]*uOther[18]*mnuOther+0.2*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(61,55) = 0.5357142857142857*m0rOther[5]*uOther[23]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[23]*mnuOther+0.159719141249985*m0rOther[0]*uOther[23]*mnuOther-0.31943828249997*m1rOther[23]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[22]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[21]*mnuOther+0.159719141249985*m0rOther[1]*uOther[21]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[20]*mnuOther+0.223606797749979*m0rOther[1]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[18]*mnuOther+0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[17]*mnuOther-0.5000000000000001*m1rOther[17]*mnuOther+0.159719141249985*m0rOther[7]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(62,48) = 0.2*m0rOther[3]*uOther[23]*mnuOther+0.223606797749979*m0rOther[5]*uOther[22]*mnuOther+0.159719141249985*m0rOther[4]*uOther[22]*mnuOther+0.25*m0rOther[0]*uOther[22]*mnuOther-0.5*m1rOther[22]*mnuOther+0.223606797749979*m0rOther[6]*uOther[21]*mnuOther+0.159719141249985*m0rOther[6]*uOther[20]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[20]*mnuOther+0.2*m0rOther[7]*uOther[19]*mnuOther+0.223606797749979*m0rOther[1]*uOther[19]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[18]*mnuOther+0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.25*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(62,49) = 0.351382110749967*m0rOther[6]*uOther[23]*mnuOther+0.2*m0rOther[2]*uOther[23]*mnuOther+0.351382110749967*m0rOther[7]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[22]*mnuOther+0.2*m0rOther[3]*uOther[21]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[20]*mnuOther+0.2*m0rOther[5]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.447213595499958*m1rOther[19]*mnuOther+0.2*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(62,50) = 0.351382110749967*m0rOther[7]*uOther[23]*mnuOther+0.2*m0rOther[1]*uOther[23]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[22]*mnuOther+0.45*m0rOther[2]*uOther[22]*mnuOther+0.223606797749979*m0rOther[4]*uOther[21]*mnuOther+0.223606797749979*m0rOther[5]*uOther[20]*mnuOther+0.159719141249985*m0rOther[4]*uOther[20]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[20]*mnuOther-0.5000000000000001*m1rOther[20]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[19]*mnuOther+0.45*m0rOther[6]*uOther[18]*mnuOther+0.2*m0rOther[7]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[16]*mnuOther; 
  data->AEM_S(62,51) = 0.3513821107499669*m0rOther[5]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[23]*mnuOther+0.2*m0rOther[0]*uOther[23]*mnuOther-0.4*m1rOther[23]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[22]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[21]*mnuOther+0.2*m0rOther[1]*uOther[21]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[20]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[19]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[19]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[18]*mnuOther+0.2*m0rOther[5]*uOther[17]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.447213595499958*m1rOther[17]*mnuOther+0.2*m0rOther[7]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(62,52) = 0.3513821107499669*m0rOther[3]*uOther[23]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[22]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[22]*mnuOther+0.159719141249985*m0rOther[0]*uOther[22]*mnuOther-0.31943828249997*m1rOther[22]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[21]*mnuOther+0.223606797749979*m0rOther[2]*uOther[21]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[20]*mnuOther+0.159719141249985*m0rOther[2]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[5]*uOther[18]*mnuOther+0.159719141249985*m0rOther[4]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[18]*mnuOther-0.5000000000000001*m1rOther[18]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[17]*mnuOther+0.159719141249985*m0rOther[6]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(62,53) = 0.3513821107499669*m0rOther[3]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[22]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[22]*mnuOther+0.223606797749979*m0rOther[0]*uOther[22]*mnuOther-0.4472135954999579*m1rOther[22]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[21]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[20]*mnuOther+0.223606797749979*m0rOther[2]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[7]*uOther[19]*mnuOther+0.2*m0rOther[1]*uOther[19]*mnuOther+0.223606797749979*m0rOther[4]*uOther[18]*mnuOther+0.2*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[6]*uOther[16]*mnuOther; 
  data->AEM_S(62,54) = 0.6173469387755102*m0rOther[7]*uOther[23]*mnuOther+0.351382110749967*m0rOther[1]*uOther[23]*mnuOther+0.9642857142857143*m0rOther[6]*uOther[22]*mnuOther+0.2874944542499729*m0rOther[2]*uOther[22]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[21]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[21]*mnuOther+0.223606797749979*m0rOther[0]*uOther[21]*mnuOther-0.4472135954999579*m1rOther[21]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[20]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[20]*mnuOther+0.159719141249985*m0rOther[0]*uOther[20]*mnuOther-0.31943828249997*m1rOther[20]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[19]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[18]*mnuOther+0.45*m0rOther[2]*uOther[18]*mnuOther+0.351382110749967*m0rOther[7]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[5]*uOther[16]*mnuOther+0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(62,55) = 0.6173469387755102*m0rOther[6]*uOther[23]*mnuOther+0.351382110749967*m0rOther[2]*uOther[23]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[22]*mnuOther+0.351382110749967*m0rOther[1]*uOther[22]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[21]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[19]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[19]*mnuOther+0.2*m0rOther[0]*uOther[19]*mnuOther-0.4*m1rOther[19]*mnuOther+0.351382110749967*m0rOther[7]*uOther[18]*mnuOther+0.2*m0rOther[1]*uOther[18]*mnuOther+0.351382110749967*m0rOther[6]*uOther[17]*mnuOther+0.2*m0rOther[2]*uOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(63,48) = 0.159719141249985*m0rOther[5]*uOther[23]*mnuOther+0.223606797749979*m0rOther[4]*uOther[23]*mnuOther+0.25*m0rOther[0]*uOther[23]*mnuOther-0.5*m1rOther[23]*mnuOther+0.2*m0rOther[3]*uOther[22]*mnuOther+0.159719141249985*m0rOther[7]*uOther[21]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[21]*mnuOther+0.223606797749979*m0rOther[7]*uOther[20]*mnuOther+0.2*m0rOther[6]*uOther[19]*mnuOther+0.223606797749979*m0rOther[2]*uOther[19]*mnuOther+0.223606797749979*m0rOther[3]*uOther[18]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[17]*mnuOther+0.25*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(63,49) = 0.2874944542499729*m0rOther[7]*uOther[23]*mnuOther+0.45*m0rOther[1]*uOther[23]*mnuOther+0.351382110749967*m0rOther[6]*uOther[22]*mnuOther+0.2*m0rOther[2]*uOther[22]*mnuOther+0.159719141249985*m0rOther[5]*uOther[21]*mnuOther+0.223606797749979*m0rOther[4]*uOther[21]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[21]*mnuOther-0.5000000000000001*m1rOther[21]*mnuOther+0.223606797749979*m0rOther[5]*uOther[20]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[19]*mnuOther+0.2*m0rOther[6]*uOther[18]*mnuOther+0.223606797749979*m0rOther[2]*uOther[18]*mnuOther+0.45*m0rOther[7]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[5]*uOther[16]*mnuOther; 
  data->AEM_S(63,50) = 0.351382110749967*m0rOther[6]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[23]*mnuOther+0.351382110749967*m0rOther[7]*uOther[22]*mnuOther+0.2*m0rOther[1]*uOther[22]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[21]*mnuOther+0.2*m0rOther[3]*uOther[20]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[19]*mnuOther+0.2*m0rOther[4]*uOther[19]*mnuOther+0.223606797749979*m0rOther[0]*uOther[19]*mnuOther-0.447213595499958*m1rOther[19]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[18]*mnuOther+0.223606797749979*m0rOther[1]*uOther[18]*mnuOther+0.2*m0rOther[6]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(63,51) = 0.7071428571428572*m0rOther[3]*uOther[23]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[22]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[22]*mnuOther+0.2*m0rOther[0]*uOther[22]*mnuOther-0.4*m1rOther[22]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[21]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[21]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[20]*mnuOther+0.2*m0rOther[2]*uOther[20]*mnuOther+0.7071428571428572*m0rOther[7]*uOther[19]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[5]*uOther[18]*mnuOther+0.2*m0rOther[4]*uOther[18]*mnuOther+0.223606797749979*m0rOther[0]*uOther[18]*mnuOther-0.447213595499958*m1rOther[18]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[17]*mnuOther+0.2*m0rOther[6]*uOther[16]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther; 
  data->AEM_S(63,52) = 0.1428571428571428*m0rOther[5]*uOther[23]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[23]*mnuOther+0.223606797749979*m0rOther[0]*uOther[23]*mnuOther-0.4472135954999579*m1rOther[23]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[22]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[21]*mnuOther+0.223606797749979*m0rOther[1]*uOther[21]*mnuOther+0.3928571428571428*m0rOther[7]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[19]*mnuOther+0.2*m0rOther[2]*uOther[19]*mnuOther+0.2*m0rOther[3]*uOther[18]*mnuOther+0.223606797749979*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[7]*uOther[16]*mnuOther; 
  data->AEM_S(63,53) = 0.5357142857142857*m0rOther[5]*uOther[23]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[23]*mnuOther+0.159719141249985*m0rOther[0]*uOther[23]*mnuOther-0.31943828249997*m1rOther[23]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[22]*mnuOther+0.5357142857142857*m0rOther[7]*uOther[21]*mnuOther+0.159719141249985*m0rOther[1]*uOther[21]*mnuOther+0.1428571428571428*m0rOther[7]*uOther[20]*mnuOther+0.223606797749979*m0rOther[1]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[6]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[19]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[18]*mnuOther+0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[4]*uOther[17]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[17]*mnuOther-0.5000000000000001*m1rOther[17]*mnuOther+0.159719141249985*m0rOther[7]*uOther[16]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[16]*mnuOther; 
  data->AEM_S(63,54) = 0.6173469387755102*m0rOther[6]*uOther[23]*mnuOther+0.351382110749967*m0rOther[2]*uOther[23]*mnuOther+0.6173469387755102*m0rOther[7]*uOther[22]*mnuOther+0.351382110749967*m0rOther[1]*uOther[22]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[21]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[20]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[19]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[19]*mnuOther+0.2*m0rOther[0]*uOther[19]*mnuOther-0.4*m1rOther[19]*mnuOther+0.351382110749967*m0rOther[7]*uOther[18]*mnuOther+0.2*m0rOther[1]*uOther[18]*mnuOther+0.351382110749967*m0rOther[6]*uOther[17]*mnuOther+0.2*m0rOther[2]*uOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther; 
  data->AEM_S(63,55) = 0.9642857142857143*m0rOther[7]*uOther[23]*mnuOther+0.2874944542499729*m0rOther[1]*uOther[23]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[22]*mnuOther+0.351382110749967*m0rOther[2]*uOther[22]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[21]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[21]*mnuOther+0.159719141249985*m0rOther[0]*uOther[21]*mnuOther-0.31943828249997*m1rOther[21]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[20]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[20]*mnuOther+0.223606797749979*m0rOther[0]*uOther[20]*mnuOther-0.4472135954999579*m1rOther[20]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[19]*mnuOther+0.351382110749967*m0rOther[6]*uOther[18]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[18]*mnuOther+0.2874944542499729*m0rOther[7]*uOther[17]*mnuOther+0.45*m0rOther[1]*uOther[17]*mnuOther+0.159719141249985*m0rOther[5]*uOther[16]*mnuOther+0.223606797749979*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfZ-m0Self*m1OtherZ)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[16] = 0.5*m0rOther[7]*m1rSelf[23]-0.5*m0rSelf[7]*m1rOther[23]+0.5*m0rOther[6]*m1rSelf[22]-0.5*m0rSelf[6]*m1rOther[22]+0.5*m0rOther[5]*m1rSelf[21]-0.5*m0rSelf[5]*m1rOther[21]+0.5*m0rOther[4]*m1rSelf[20]-0.5*m0rSelf[4]*m1rOther[20]+0.5*m0rOther[3]*m1rSelf[19]-0.5*m0rSelf[3]*m1rOther[19]+0.5*m0rOther[2]*m1rSelf[18]-0.5*m0rSelf[2]*m1rOther[18]+0.5*m0rOther[1]*m1rSelf[17]-0.5*m0rSelf[1]*m1rOther[17]+0.5*m0rOther[0]*m1rSelf[16]-0.5*m0rSelf[0]*m1rOther[16]; 
  m1EffD[17] = 0.5000000000000001*m0rOther[5]*m1rSelf[23]-0.5000000000000001*m0rSelf[5]*m1rOther[23]+0.447213595499958*m0rOther[3]*m1rSelf[22]-0.447213595499958*m0rSelf[3]*m1rOther[22]+0.5000000000000001*m0rOther[7]*m1rSelf[21]-0.5000000000000001*m0rSelf[7]*m1rOther[21]+0.4472135954999579*m0rOther[1]*m1rSelf[20]-0.4472135954999579*m0rSelf[1]*m1rOther[20]+0.447213595499958*m0rOther[6]*m1rSelf[19]+0.5*m0rOther[2]*m1rSelf[19]-0.447213595499958*m0rSelf[6]*m1rOther[19]-0.5*m0rSelf[2]*m1rOther[19]+0.5*m0rOther[3]*m1rSelf[18]-0.5*m0rSelf[3]*m1rOther[18]+0.4472135954999579*m0rOther[4]*m1rSelf[17]+0.5*m0rOther[0]*m1rSelf[17]-0.4472135954999579*m0rSelf[4]*m1rOther[17]-0.5*m0rSelf[0]*m1rOther[17]+0.5*m0rOther[1]*m1rSelf[16]-0.5*m0rSelf[1]*m1rOther[16]; 
  m1EffD[18] = 0.447213595499958*m0rOther[3]*m1rSelf[23]-0.447213595499958*m0rSelf[3]*m1rOther[23]+0.5000000000000001*m0rOther[4]*m1rSelf[22]-0.5000000000000001*m0rSelf[4]*m1rOther[22]+0.4472135954999579*m0rOther[2]*m1rSelf[21]-0.4472135954999579*m0rSelf[2]*m1rOther[21]+0.5000000000000001*m0rOther[6]*m1rSelf[20]-0.5000000000000001*m0rSelf[6]*m1rOther[20]+0.447213595499958*m0rOther[7]*m1rSelf[19]+0.5*m0rOther[1]*m1rSelf[19]-0.447213595499958*m0rSelf[7]*m1rOther[19]-0.5*m0rSelf[1]*m1rOther[19]+0.4472135954999579*m0rOther[5]*m1rSelf[18]+0.5*m0rOther[0]*m1rSelf[18]-0.4472135954999579*m0rSelf[5]*m1rOther[18]-0.5*m0rSelf[0]*m1rOther[18]+0.5*m0rOther[3]*m1rSelf[17]-0.5*m0rSelf[3]*m1rOther[17]+0.5*m0rOther[2]*m1rSelf[16]-0.5*m0rSelf[2]*m1rOther[16]; 
  m1EffD[19] = 0.4*m0rOther[6]*m1rSelf[23]+0.447213595499958*m0rOther[2]*m1rSelf[23]-0.4*m0rSelf[6]*m1rOther[23]-0.447213595499958*m0rSelf[2]*m1rOther[23]+0.4*m0rOther[7]*m1rSelf[22]+0.447213595499958*m0rOther[1]*m1rSelf[22]-0.4*m0rSelf[7]*m1rOther[22]-0.447213595499958*m0rSelf[1]*m1rOther[22]+0.4472135954999579*m0rOther[3]*m1rSelf[21]-0.4472135954999579*m0rSelf[3]*m1rOther[21]+0.4472135954999579*m0rOther[3]*m1rSelf[20]-0.4472135954999579*m0rSelf[3]*m1rOther[20]+0.4472135954999579*m0rOther[5]*m1rSelf[19]+0.4472135954999579*m0rOther[4]*m1rSelf[19]+0.5*m0rOther[0]*m1rSelf[19]-0.4472135954999579*m0rSelf[5]*m1rOther[19]-0.4472135954999579*m0rSelf[4]*m1rOther[19]-0.5*m0rSelf[0]*m1rOther[19]+0.447213595499958*m0rOther[7]*m1rSelf[18]+0.5*m0rOther[1]*m1rSelf[18]-0.447213595499958*m0rSelf[7]*m1rOther[18]-0.5*m0rSelf[1]*m1rOther[18]+0.447213595499958*m0rOther[6]*m1rSelf[17]+0.5*m0rOther[2]*m1rSelf[17]-0.447213595499958*m0rSelf[6]*m1rOther[17]-0.5*m0rSelf[2]*m1rOther[17]+0.5*m0rOther[3]*m1rSelf[16]-0.5*m0rSelf[3]*m1rOther[16]; 
  m1EffD[20] = 0.4472135954999579*m0rOther[7]*m1rSelf[23]-0.4472135954999579*m0rSelf[7]*m1rOther[23]+0.31943828249997*m0rOther[6]*m1rSelf[22]+0.5000000000000001*m0rOther[2]*m1rSelf[22]-0.31943828249997*m0rSelf[6]*m1rOther[22]-0.5000000000000001*m0rSelf[2]*m1rOther[22]+0.31943828249997*m0rOther[4]*m1rSelf[20]+0.5*m0rOther[0]*m1rSelf[20]-0.31943828249997*m0rSelf[4]*m1rOther[20]-0.5*m0rSelf[0]*m1rOther[20]+0.4472135954999579*m0rOther[3]*m1rSelf[19]-0.4472135954999579*m0rSelf[3]*m1rOther[19]+0.5000000000000001*m0rOther[6]*m1rSelf[18]-0.5000000000000001*m0rSelf[6]*m1rOther[18]+0.4472135954999579*m0rOther[1]*m1rSelf[17]-0.4472135954999579*m0rSelf[1]*m1rOther[17]+0.5*m0rOther[4]*m1rSelf[16]-0.5*m0rSelf[4]*m1rOther[16]; 
  m1EffD[21] = 0.31943828249997*m0rOther[7]*m1rSelf[23]+0.5000000000000001*m0rOther[1]*m1rSelf[23]-0.31943828249997*m0rSelf[7]*m1rOther[23]-0.5000000000000001*m0rSelf[1]*m1rOther[23]+0.4472135954999579*m0rOther[6]*m1rSelf[22]-0.4472135954999579*m0rSelf[6]*m1rOther[22]+0.31943828249997*m0rOther[5]*m1rSelf[21]+0.5*m0rOther[0]*m1rSelf[21]-0.31943828249997*m0rSelf[5]*m1rOther[21]-0.5*m0rSelf[0]*m1rOther[21]+0.4472135954999579*m0rOther[3]*m1rSelf[19]-0.4472135954999579*m0rSelf[3]*m1rOther[19]+0.4472135954999579*m0rOther[2]*m1rSelf[18]-0.4472135954999579*m0rSelf[2]*m1rOther[18]+0.5000000000000001*m0rOther[7]*m1rSelf[17]-0.5000000000000001*m0rSelf[7]*m1rOther[17]+0.5*m0rOther[5]*m1rSelf[16]-0.5*m0rSelf[5]*m1rOther[16]; 
  m1EffD[22] = 0.4*m0rOther[3]*m1rSelf[23]-0.4*m0rSelf[3]*m1rOther[23]+0.4472135954999579*m0rOther[5]*m1rSelf[22]+0.31943828249997*m0rOther[4]*m1rSelf[22]+0.5*m0rOther[0]*m1rSelf[22]-0.4472135954999579*m0rSelf[5]*m1rOther[22]-0.31943828249997*m0rSelf[4]*m1rOther[22]-0.5*m0rSelf[0]*m1rOther[22]+0.4472135954999579*m0rOther[6]*m1rSelf[21]-0.4472135954999579*m0rSelf[6]*m1rOther[21]+0.31943828249997*m0rOther[6]*m1rSelf[20]+0.5000000000000001*m0rOther[2]*m1rSelf[20]-0.31943828249997*m0rSelf[6]*m1rOther[20]-0.5000000000000001*m0rSelf[2]*m1rOther[20]+0.4*m0rOther[7]*m1rSelf[19]+0.447213595499958*m0rOther[1]*m1rSelf[19]-0.4*m0rSelf[7]*m1rOther[19]-0.447213595499958*m0rSelf[1]*m1rOther[19]+0.5000000000000001*m0rOther[4]*m1rSelf[18]-0.5000000000000001*m0rSelf[4]*m1rOther[18]+0.447213595499958*m0rOther[3]*m1rSelf[17]-0.447213595499958*m0rSelf[3]*m1rOther[17]+0.5*m0rOther[6]*m1rSelf[16]-0.5*m0rSelf[6]*m1rOther[16]; 
  m1EffD[23] = 0.31943828249997*m0rOther[5]*m1rSelf[23]+0.4472135954999579*m0rOther[4]*m1rSelf[23]+0.5*m0rOther[0]*m1rSelf[23]-0.31943828249997*m0rSelf[5]*m1rOther[23]-0.4472135954999579*m0rSelf[4]*m1rOther[23]-0.5*m0rSelf[0]*m1rOther[23]+0.4*m0rOther[3]*m1rSelf[22]-0.4*m0rSelf[3]*m1rOther[22]+0.31943828249997*m0rOther[7]*m1rSelf[21]+0.5000000000000001*m0rOther[1]*m1rSelf[21]-0.31943828249997*m0rSelf[7]*m1rOther[21]-0.5000000000000001*m0rSelf[1]*m1rOther[21]+0.4472135954999579*m0rOther[7]*m1rSelf[20]-0.4472135954999579*m0rSelf[7]*m1rOther[20]+0.4*m0rOther[6]*m1rSelf[19]+0.447213595499958*m0rOther[2]*m1rSelf[19]-0.4*m0rSelf[6]*m1rOther[19]-0.447213595499958*m0rSelf[2]*m1rOther[19]+0.447213595499958*m0rOther[3]*m1rSelf[18]-0.447213595499958*m0rSelf[3]*m1rOther[18]+0.5000000000000001*m0rOther[5]*m1rSelf[17]-0.5000000000000001*m0rSelf[5]*m1rOther[17]+0.5*m0rOther[7]*m1rSelf[16]-0.5*m0rSelf[7]*m1rOther[16]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[16],m1EffD[17],m1EffD[18],m1EffD[19],m1EffD[20],m1EffD[21],m1EffD[22],m1EffD[23]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+16,8,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 3 of momentum relaxation. 
  m1Relax[16] += (-2.0*m1EffD[16]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[16]*mnuSelf-1.0*m1rOther[16]*mnuOther; 
  m1Relax[17] += (-2.0*m1EffD[17]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[17]*mnuSelf-1.0*m1rOther[17]*mnuOther; 
  m1Relax[18] += (-2.0*m1EffD[18]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[18]*mnuSelf-1.0*m1rOther[18]*mnuOther; 
  m1Relax[19] += (-2.0*m1EffD[19]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[19]*mnuSelf-1.0*m1rOther[19]*mnuOther; 
  m1Relax[20] += (-2.0*m1EffD[20]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[20]*mnuSelf-1.0*m1rOther[20]*mnuOther; 
  m1Relax[21] += (-2.0*m1EffD[21]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[21]*mnuSelf-1.0*m1rOther[21]*mnuOther; 
  m1Relax[22] += (-2.0*m1EffD[22]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[22]*mnuSelf-1.0*m1rOther[22]*mnuOther; 
  m1Relax[23] += (-2.0*m1EffD[23]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[23]*mnuSelf-1.0*m1rOther[23]*mnuOther; 
 
  double ucMSelf[8]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.5*cMSelf[a0+7]*uSelf[a0+7]+0.5*cMSelf[a0+6]*uSelf[a0+6]+0.5*cMSelf[a0+5]*uSelf[a0+5]+0.5*cMSelf[a0+4]*uSelf[a0+4]+0.5*cMSelf[a0+3]*uSelf[a0+3]+0.5*cMSelf[a0+2]*uSelf[a0+2]+0.5*cMSelf[a0+1]*uSelf[a0+1]+0.5*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.5000000000000001*cMSelf[a0+5]*uSelf[a0+7]+0.5000000000000001*uSelf[a0+5]*cMSelf[a0+7]+0.447213595499958*cMSelf[a0+3]*uSelf[a0+6]+0.447213595499958*uSelf[a0+3]*cMSelf[a0+6]+0.4472135954999579*cMSelf[a0+1]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+1]*cMSelf[a0+4]+0.5*cMSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.447213595499958*cMSelf[a0+3]*uSelf[a0+7]+0.447213595499958*uSelf[a0+3]*cMSelf[a0+7]+0.5000000000000001*cMSelf[a0+4]*uSelf[a0+6]+0.5000000000000001*uSelf[a0+4]*cMSelf[a0+6]+0.4472135954999579*cMSelf[a0+2]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+2]*cMSelf[a0+5]+0.5*cMSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*cMSelf[a0+2]; 
    ucMSelf[3] += 0.4*cMSelf[a0+6]*uSelf[a0+7]+0.447213595499958*cMSelf[a0+2]*uSelf[a0+7]+0.4*uSelf[a0+6]*cMSelf[a0+7]+0.447213595499958*uSelf[a0+2]*cMSelf[a0+7]+0.447213595499958*cMSelf[a0+1]*uSelf[a0+6]+0.447213595499958*uSelf[a0+1]*cMSelf[a0+6]+0.4472135954999579*cMSelf[a0+3]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+3]*cMSelf[a0+5]+0.4472135954999579*cMSelf[a0+3]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+3]*cMSelf[a0+4]+0.5*cMSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*cMSelf[a0+3]+0.5*cMSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*cMSelf[a0+2]; 
    ucMSelf[4] += 0.4472135954999579*cMSelf[a0+7]*uSelf[a0+7]+0.31943828249997*cMSelf[a0+6]*uSelf[a0+6]+0.5000000000000001*cMSelf[a0+2]*uSelf[a0+6]+0.5000000000000001*uSelf[a0+2]*cMSelf[a0+6]+0.31943828249997*cMSelf[a0+4]*uSelf[a0+4]+0.5*cMSelf[a0]*uSelf[a0+4]+0.5*uSelf[a0]*cMSelf[a0+4]+0.4472135954999579*cMSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*cMSelf[a0+1]*uSelf[a0+1]; 
    ucMSelf[5] += 0.31943828249997*cMSelf[a0+7]*uSelf[a0+7]+0.5000000000000001*cMSelf[a0+1]*uSelf[a0+7]+0.5000000000000001*uSelf[a0+1]*cMSelf[a0+7]+0.4472135954999579*cMSelf[a0+6]*uSelf[a0+6]+0.31943828249997*cMSelf[a0+5]*uSelf[a0+5]+0.5*cMSelf[a0]*uSelf[a0+5]+0.5*uSelf[a0]*cMSelf[a0+5]+0.4472135954999579*cMSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*cMSelf[a0+2]*uSelf[a0+2]; 
    ucMSelf[6] += 0.4*cMSelf[a0+3]*uSelf[a0+7]+0.4*uSelf[a0+3]*cMSelf[a0+7]+0.4472135954999579*cMSelf[a0+5]*uSelf[a0+6]+0.31943828249997*cMSelf[a0+4]*uSelf[a0+6]+0.5*cMSelf[a0]*uSelf[a0+6]+0.4472135954999579*uSelf[a0+5]*cMSelf[a0+6]+0.31943828249997*uSelf[a0+4]*cMSelf[a0+6]+0.5*uSelf[a0]*cMSelf[a0+6]+0.5000000000000001*cMSelf[a0+2]*uSelf[a0+4]+0.5000000000000001*uSelf[a0+2]*cMSelf[a0+4]+0.447213595499958*cMSelf[a0+1]*uSelf[a0+3]+0.447213595499958*uSelf[a0+1]*cMSelf[a0+3]; 
    ucMSelf[7] += 0.31943828249997*cMSelf[a0+5]*uSelf[a0+7]+0.4472135954999579*cMSelf[a0+4]*uSelf[a0+7]+0.5*cMSelf[a0]*uSelf[a0+7]+0.31943828249997*uSelf[a0+5]*cMSelf[a0+7]+0.4472135954999579*uSelf[a0+4]*cMSelf[a0+7]+0.5*uSelf[a0]*cMSelf[a0+7]+0.4*cMSelf[a0+3]*uSelf[a0+6]+0.4*uSelf[a0+3]*cMSelf[a0+6]+0.5000000000000001*cMSelf[a0+1]*uSelf[a0+5]+0.5000000000000001*uSelf[a0+1]*cMSelf[a0+5]+0.447213595499958*cMSelf[a0+2]*uSelf[a0+3]+0.447213595499958*uSelf[a0+2]*cMSelf[a0+3]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(56,24) = 0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(56,25) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(56,26) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(56,27) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(56,28) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(56,29) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(56,30) = 0.5*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(56,31) = 0.5*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(57,24) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(57,25) = 0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(57,26) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(57,27) = 0.447213595499958*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(57,28) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(57,29) = 0.5000000000000001*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(57,30) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(57,31) = 0.5000000000000001*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(58,24) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(58,25) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(58,26) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(58,27) = 0.447213595499958*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(58,28) = 0.5000000000000001*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(58,29) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(58,30) = 0.5000000000000001*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(58,31) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(59,24) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(59,25) = 0.447213595499958*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(59,26) = 0.447213595499958*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(59,27) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(59,28) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(59,29) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(59,30) = 0.4*ucMSelf[7]*mnuSelf+1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.447213595499958*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(59,31) = 0.4*ucMSelf[6]*mnuSelf+1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.447213595499958*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(60,24) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(60,25) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(60,26) = 0.5000000000000001*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(60,27) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(60,28) = 0.31943828249997*ucMSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(60,30) = 0.31943828249997*ucMSelf[6]*mnuSelf+0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+0.5000000000000001*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(60,31) = 0.4472135954999579*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(61,24) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(61,25) = 0.5000000000000001*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(61,26) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(61,27) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(61,29) = 0.31943828249997*ucMSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(61,30) = 0.4472135954999579*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(61,31) = 0.31943828249997*ucMSelf[7]*mnuSelf+0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+0.5000000000000001*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(62,24) = 0.5*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(62,25) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(62,26) = 0.5000000000000001*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(62,27) = 0.4*ucMSelf[7]*mnuSelf+1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.447213595499958*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(62,28) = 0.31943828249997*ucMSelf[6]*mnuSelf+0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+0.5000000000000001*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(62,29) = 0.4472135954999579*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(62,30) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.31943828249997*ucMSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(62,31) = 0.4*ucMSelf[3]*mnuSelf+1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(63,24) = 0.5*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(63,25) = 0.5000000000000001*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(63,26) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(63,27) = 0.4*ucMSelf[6]*mnuSelf+1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.447213595499958*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(63,28) = 0.4472135954999579*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(63,29) = 0.31943828249997*ucMSelf[7]*mnuSelf+0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+0.5000000000000001*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(63,30) = 0.4*ucMSelf[3]*mnuSelf+1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(63,31) = 0.31943828249997*ucMSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  double ucMOther[8]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.5*cMOther[a0+7]*uOther[a0+7]+0.5*cMOther[a0+6]*uOther[a0+6]+0.5*cMOther[a0+5]*uOther[a0+5]+0.5*cMOther[a0+4]*uOther[a0+4]+0.5*cMOther[a0+3]*uOther[a0+3]+0.5*cMOther[a0+2]*uOther[a0+2]+0.5*cMOther[a0+1]*uOther[a0+1]+0.5*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.5000000000000001*cMOther[a0+5]*uOther[a0+7]+0.5000000000000001*uOther[a0+5]*cMOther[a0+7]+0.447213595499958*cMOther[a0+3]*uOther[a0+6]+0.447213595499958*uOther[a0+3]*cMOther[a0+6]+0.4472135954999579*cMOther[a0+1]*uOther[a0+4]+0.4472135954999579*uOther[a0+1]*cMOther[a0+4]+0.5*cMOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+1]+0.5*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.447213595499958*cMOther[a0+3]*uOther[a0+7]+0.447213595499958*uOther[a0+3]*cMOther[a0+7]+0.5000000000000001*cMOther[a0+4]*uOther[a0+6]+0.5000000000000001*uOther[a0+4]*cMOther[a0+6]+0.4472135954999579*cMOther[a0+2]*uOther[a0+5]+0.4472135954999579*uOther[a0+2]*cMOther[a0+5]+0.5*cMOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+2]+0.5*uOther[a0]*cMOther[a0+2]; 
    ucMOther[3] += 0.4*cMOther[a0+6]*uOther[a0+7]+0.447213595499958*cMOther[a0+2]*uOther[a0+7]+0.4*uOther[a0+6]*cMOther[a0+7]+0.447213595499958*uOther[a0+2]*cMOther[a0+7]+0.447213595499958*cMOther[a0+1]*uOther[a0+6]+0.447213595499958*uOther[a0+1]*cMOther[a0+6]+0.4472135954999579*cMOther[a0+3]*uOther[a0+5]+0.4472135954999579*uOther[a0+3]*cMOther[a0+5]+0.4472135954999579*cMOther[a0+3]*uOther[a0+4]+0.4472135954999579*uOther[a0+3]*cMOther[a0+4]+0.5*cMOther[a0]*uOther[a0+3]+0.5*uOther[a0]*cMOther[a0+3]+0.5*cMOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*cMOther[a0+2]; 
    ucMOther[4] += 0.4472135954999579*cMOther[a0+7]*uOther[a0+7]+0.31943828249997*cMOther[a0+6]*uOther[a0+6]+0.5000000000000001*cMOther[a0+2]*uOther[a0+6]+0.5000000000000001*uOther[a0+2]*cMOther[a0+6]+0.31943828249997*cMOther[a0+4]*uOther[a0+4]+0.5*cMOther[a0]*uOther[a0+4]+0.5*uOther[a0]*cMOther[a0+4]+0.4472135954999579*cMOther[a0+3]*uOther[a0+3]+0.4472135954999579*cMOther[a0+1]*uOther[a0+1]; 
    ucMOther[5] += 0.31943828249997*cMOther[a0+7]*uOther[a0+7]+0.5000000000000001*cMOther[a0+1]*uOther[a0+7]+0.5000000000000001*uOther[a0+1]*cMOther[a0+7]+0.4472135954999579*cMOther[a0+6]*uOther[a0+6]+0.31943828249997*cMOther[a0+5]*uOther[a0+5]+0.5*cMOther[a0]*uOther[a0+5]+0.5*uOther[a0]*cMOther[a0+5]+0.4472135954999579*cMOther[a0+3]*uOther[a0+3]+0.4472135954999579*cMOther[a0+2]*uOther[a0+2]; 
    ucMOther[6] += 0.4*cMOther[a0+3]*uOther[a0+7]+0.4*uOther[a0+3]*cMOther[a0+7]+0.4472135954999579*cMOther[a0+5]*uOther[a0+6]+0.31943828249997*cMOther[a0+4]*uOther[a0+6]+0.5*cMOther[a0]*uOther[a0+6]+0.4472135954999579*uOther[a0+5]*cMOther[a0+6]+0.31943828249997*uOther[a0+4]*cMOther[a0+6]+0.5*uOther[a0]*cMOther[a0+6]+0.5000000000000001*cMOther[a0+2]*uOther[a0+4]+0.5000000000000001*uOther[a0+2]*cMOther[a0+4]+0.447213595499958*cMOther[a0+1]*uOther[a0+3]+0.447213595499958*uOther[a0+1]*cMOther[a0+3]; 
    ucMOther[7] += 0.31943828249997*cMOther[a0+5]*uOther[a0+7]+0.4472135954999579*cMOther[a0+4]*uOther[a0+7]+0.5*cMOther[a0]*uOther[a0+7]+0.31943828249997*uOther[a0+5]*cMOther[a0+7]+0.4472135954999579*uOther[a0+4]*cMOther[a0+7]+0.5*uOther[a0]*cMOther[a0+7]+0.4*cMOther[a0+3]*uOther[a0+6]+0.4*uOther[a0+3]*cMOther[a0+6]+0.5000000000000001*cMOther[a0+1]*uOther[a0+5]+0.5000000000000001*uOther[a0+1]*cMOther[a0+5]+0.447213595499958*cMOther[a0+2]*uOther[a0+3]+0.447213595499958*uOther[a0+2]*cMOther[a0+3]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(56,56) = (-0.5*ucMOther[0]*mnuOther)-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(56,57) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(56,58) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(56,59) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(56,60) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(56,61) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(56,62) = (-0.5*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5*cEOther[6]*mnuOther; 
  data->AEM_S(56,63) = (-0.5*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5*cEOther[7]*mnuOther; 
  data->AEM_S(57,56) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(57,57) = (-0.4472135954999579*ucMOther[4]*mnuOther)-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(57,58) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(57,59) = (-0.447213595499958*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-0.5*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(57,60) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(57,61) = (-0.5000000000000001*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(57,62) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(57,63) = (-0.5000000000000001*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(58,56) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(58,57) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(58,58) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(58,59) = (-0.447213595499958*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-0.5*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(58,60) = (-0.5000000000000001*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(58,61) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(58,62) = (-0.5000000000000001*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(58,63) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(59,56) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(59,57) = (-0.447213595499958*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-0.5*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(59,58) = (-0.447213595499958*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-0.5*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(59,59) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.4472135954999579*ucMOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(59,60) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(59,61) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(59,62) = (-0.4*ucMOther[7]*mnuOther)-1.2*m0rOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.447213595499958*ucMOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(59,63) = (-0.4*ucMOther[6]*mnuOther)-1.2*m0rOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.447213595499958*ucMOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(60,56) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(60,57) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(60,58) = (-0.5000000000000001*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(60,59) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(60,60) = (-0.31943828249997*ucMOther[4]*mnuOther)-0.9583148474999099*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(60,62) = (-0.31943828249997*ucMOther[6]*mnuOther)-0.9583148474999099*m0rOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-0.5000000000000001*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(60,63) = (-0.4472135954999579*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(61,56) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(61,57) = (-0.5000000000000001*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(61,58) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(61,59) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(61,61) = (-0.31943828249997*ucMOther[5]*mnuOther)-0.9583148474999099*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(61,62) = (-0.4472135954999579*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(61,63) = (-0.31943828249997*ucMOther[7]*mnuOther)-0.9583148474999099*m0rOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-0.5000000000000001*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(62,56) = (-0.5*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5*cEOther[6]*mnuOther; 
  data->AEM_S(62,57) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(62,58) = (-0.5000000000000001*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(62,59) = (-0.4*ucMOther[7]*mnuOther)-1.2*m0rOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.447213595499958*ucMOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(62,60) = (-0.31943828249997*ucMOther[6]*mnuOther)-0.9583148474999099*m0rOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-0.5000000000000001*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(62,61) = (-0.4472135954999579*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(62,62) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.31943828249997*ucMOther[4]*mnuOther-0.9583148474999099*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(62,63) = (-0.4*ucMOther[3]*mnuOther)-1.2*m0rOther[3]*mnuOther+0.4*cEOther[3]*mnuOther; 
  data->AEM_S(63,56) = (-0.5*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5*cEOther[7]*mnuOther; 
  data->AEM_S(63,57) = (-0.5000000000000001*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(63,58) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(63,59) = (-0.4*ucMOther[6]*mnuOther)-1.2*m0rOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.447213595499958*ucMOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(63,60) = (-0.4472135954999579*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(63,61) = (-0.31943828249997*ucMOther[7]*mnuOther)-0.9583148474999099*m0rOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-0.5000000000000001*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(63,62) = (-0.4*ucMOther[3]*mnuOther)-1.2*m0rOther[3]*mnuOther+0.4*cEOther[3]*mnuOther; 
  data->AEM_S(63,63) = (-0.31943828249997*ucMOther[5]*mnuOther)-0.9583148474999099*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.4472135954999579*ucMOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[8]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinESelf[0] += 0.5*m1rSelf[a0+7]*uSelf[a0+7]+0.5*m1rSelf[a0+6]*uSelf[a0+6]+0.5*m1rSelf[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0+3]*uSelf[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]; 
    kinESelf[1] += 0.5000000000000001*m1rSelf[a0+5]*uSelf[a0+7]+0.5000000000000001*uSelf[a0+5]*m1rSelf[a0+7]+0.447213595499958*m1rSelf[a0+3]*uSelf[a0+6]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+6]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+1]*m1rSelf[a0+4]+0.5*m1rSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]; 
    kinESelf[2] += 0.447213595499958*m1rSelf[a0+3]*uSelf[a0+7]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+7]+0.5000000000000001*m1rSelf[a0+4]*uSelf[a0+6]+0.5000000000000001*uSelf[a0+4]*m1rSelf[a0+6]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+2]*m1rSelf[a0+5]+0.5*m1rSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]; 
    kinESelf[3] += 0.4*m1rSelf[a0+6]*uSelf[a0+7]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+7]+0.4*uSelf[a0+6]*m1rSelf[a0+7]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+7]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+6]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+6]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+5]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+4]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]; 
    kinESelf[4] += 0.4472135954999579*m1rSelf[a0+7]*uSelf[a0+7]+0.31943828249997*m1rSelf[a0+6]*uSelf[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+6]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+4]+0.5*uSelf[a0]*m1rSelf[a0+4]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+1]; 
    kinESelf[5] += 0.31943828249997*m1rSelf[a0+7]*uSelf[a0+7]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+7]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+7]+0.4472135954999579*m1rSelf[a0+6]*uSelf[a0+6]+0.31943828249997*m1rSelf[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0]*uSelf[a0+5]+0.5*uSelf[a0]*m1rSelf[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+2]; 
    kinESelf[6] += 0.4*m1rSelf[a0+3]*uSelf[a0+7]+0.4*uSelf[a0+3]*m1rSelf[a0+7]+0.4472135954999579*m1rSelf[a0+5]*uSelf[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+6]+0.5*m1rSelf[a0]*uSelf[a0+6]+0.4472135954999579*uSelf[a0+5]*m1rSelf[a0+6]+0.31943828249997*uSelf[a0+4]*m1rSelf[a0+6]+0.5*uSelf[a0]*m1rSelf[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+4]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+4]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+3]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+3]; 
    kinESelf[7] += 0.31943828249997*m1rSelf[a0+5]*uSelf[a0+7]+0.4472135954999579*m1rSelf[a0+4]*uSelf[a0+7]+0.5*m1rSelf[a0]*uSelf[a0+7]+0.31943828249997*uSelf[a0+5]*m1rSelf[a0+7]+0.4472135954999579*uSelf[a0+4]*m1rSelf[a0+7]+0.5*uSelf[a0]*m1rSelf[a0+7]+0.4*m1rSelf[a0+3]*uSelf[a0+6]+0.4*uSelf[a0+3]*m1rSelf[a0+6]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+5]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+5]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+3]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+3]; 
  } 
 
  double kinEOther[8]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    kinEOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.5*m1rOther[a0+7]*uOther[a0+7]+0.5*m1rOther[a0+6]*uOther[a0+6]+0.5*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.5000000000000001*m1rOther[a0+5]*uOther[a0+7]+0.5000000000000001*uOther[a0+5]*m1rOther[a0+7]+0.447213595499958*m1rOther[a0+3]*uOther[a0+6]+0.447213595499958*uOther[a0+3]*m1rOther[a0+6]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+4]+0.4472135954999579*uOther[a0+1]*m1rOther[a0+4]+0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.447213595499958*m1rOther[a0+3]*uOther[a0+7]+0.447213595499958*uOther[a0+3]*m1rOther[a0+7]+0.5000000000000001*m1rOther[a0+4]*uOther[a0+6]+0.5000000000000001*uOther[a0+4]*m1rOther[a0+6]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+5]+0.4472135954999579*uOther[a0+2]*m1rOther[a0+5]+0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    kinEOther[3] += 0.4*m1rOther[a0+6]*uOther[a0+7]+0.447213595499958*m1rOther[a0+2]*uOther[a0+7]+0.4*uOther[a0+6]*m1rOther[a0+7]+0.447213595499958*uOther[a0+2]*m1rOther[a0+7]+0.447213595499958*m1rOther[a0+1]*uOther[a0+6]+0.447213595499958*uOther[a0+1]*m1rOther[a0+6]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+5]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+4]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
    kinEOther[4] += 0.4472135954999579*m1rOther[a0+7]*uOther[a0+7]+0.31943828249997*m1rOther[a0+6]*uOther[a0+6]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+6]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+6]+0.31943828249997*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+4]+0.5*uOther[a0]*m1rOther[a0+4]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+1]; 
    kinEOther[5] += 0.31943828249997*m1rOther[a0+7]*uOther[a0+7]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+7]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+7]+0.4472135954999579*m1rOther[a0+6]*uOther[a0+6]+0.31943828249997*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rOther[a0]*uOther[a0+5]+0.5*uOther[a0]*m1rOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+2]; 
    kinEOther[6] += 0.4*m1rOther[a0+3]*uOther[a0+7]+0.4*uOther[a0+3]*m1rOther[a0+7]+0.4472135954999579*m1rOther[a0+5]*uOther[a0+6]+0.31943828249997*m1rOther[a0+4]*uOther[a0+6]+0.5*m1rOther[a0]*uOther[a0+6]+0.4472135954999579*uOther[a0+5]*m1rOther[a0+6]+0.31943828249997*uOther[a0+4]*m1rOther[a0+6]+0.5*uOther[a0]*m1rOther[a0+6]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+4]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+4]+0.447213595499958*m1rOther[a0+1]*uOther[a0+3]+0.447213595499958*uOther[a0+1]*m1rOther[a0+3]; 
    kinEOther[7] += 0.31943828249997*m1rOther[a0+5]*uOther[a0+7]+0.4472135954999579*m1rOther[a0+4]*uOther[a0+7]+0.5*m1rOther[a0]*uOther[a0+7]+0.31943828249997*uOther[a0+5]*m1rOther[a0+7]+0.4472135954999579*uOther[a0+4]*m1rOther[a0+7]+0.5*uOther[a0]*m1rOther[a0+7]+0.4*m1rOther[a0+3]*uOther[a0+6]+0.4*uOther[a0+3]*m1rOther[a0+6]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+5]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+5]+0.447213595499958*m1rOther[a0+2]*uOther[a0+3]+0.447213595499958*uOther[a0+2]*m1rOther[a0+3]; 
  } 
 
  double relKinE[8]; 
  // zero out array with dot product of uSelf-uOther and m1EffD. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.5*m1EffD[a0+7]*uSelf[a0+7]-0.5*m1EffD[a0+7]*uOther[a0+7]+0.5*m1EffD[a0+6]*uSelf[a0+6]-0.5*m1EffD[a0+6]*uOther[a0+6]+0.5*m1EffD[a0+5]*uSelf[a0+5]-0.5*m1EffD[a0+5]*uOther[a0+5]+0.5*m1EffD[a0+4]*uSelf[a0+4]-0.5*m1EffD[a0+4]*uOther[a0+4]+0.5*m1EffD[a0+3]*uSelf[a0+3]-0.5*m1EffD[a0+3]*uOther[a0+3]+0.5*m1EffD[a0+2]*uSelf[a0+2]-0.5*m1EffD[a0+2]*uOther[a0+2]+0.5*m1EffD[a0+1]*uSelf[a0+1]-0.5*m1EffD[a0+1]*uOther[a0+1]+0.5*m1EffD[a0]*uSelf[a0]-0.5*m1EffD[a0]*uOther[a0]; 
    relKinE[1] += 0.5000000000000001*m1EffD[a0+5]*uSelf[a0+7]-0.5000000000000001*m1EffD[a0+5]*uOther[a0+7]+0.5000000000000001*uSelf[a0+5]*m1EffD[a0+7]-0.5000000000000001*uOther[a0+5]*m1EffD[a0+7]+0.447213595499958*m1EffD[a0+3]*uSelf[a0+6]-0.447213595499958*m1EffD[a0+3]*uOther[a0+6]+0.447213595499958*uSelf[a0+3]*m1EffD[a0+6]-0.447213595499958*uOther[a0+3]*m1EffD[a0+6]+0.4472135954999579*m1EffD[a0+1]*uSelf[a0+4]-0.4472135954999579*m1EffD[a0+1]*uOther[a0+4]+0.4472135954999579*uSelf[a0+1]*m1EffD[a0+4]-0.4472135954999579*uOther[a0+1]*m1EffD[a0+4]+0.5*m1EffD[a0+2]*uSelf[a0+3]-0.5*m1EffD[a0+2]*uOther[a0+3]+0.5*uSelf[a0+2]*m1EffD[a0+3]-0.5*uOther[a0+2]*m1EffD[a0+3]+0.5*m1EffD[a0]*uSelf[a0+1]-0.5*m1EffD[a0]*uOther[a0+1]+0.5*uSelf[a0]*m1EffD[a0+1]-0.5*uOther[a0]*m1EffD[a0+1]; 
    relKinE[2] += 0.447213595499958*m1EffD[a0+3]*uSelf[a0+7]-0.447213595499958*m1EffD[a0+3]*uOther[a0+7]+0.447213595499958*uSelf[a0+3]*m1EffD[a0+7]-0.447213595499958*uOther[a0+3]*m1EffD[a0+7]+0.5000000000000001*m1EffD[a0+4]*uSelf[a0+6]-0.5000000000000001*m1EffD[a0+4]*uOther[a0+6]+0.5000000000000001*uSelf[a0+4]*m1EffD[a0+6]-0.5000000000000001*uOther[a0+4]*m1EffD[a0+6]+0.4472135954999579*m1EffD[a0+2]*uSelf[a0+5]-0.4472135954999579*m1EffD[a0+2]*uOther[a0+5]+0.4472135954999579*uSelf[a0+2]*m1EffD[a0+5]-0.4472135954999579*uOther[a0+2]*m1EffD[a0+5]+0.5*m1EffD[a0+1]*uSelf[a0+3]-0.5*m1EffD[a0+1]*uOther[a0+3]+0.5*uSelf[a0+1]*m1EffD[a0+3]-0.5*uOther[a0+1]*m1EffD[a0+3]+0.5*m1EffD[a0]*uSelf[a0+2]-0.5*m1EffD[a0]*uOther[a0+2]+0.5*uSelf[a0]*m1EffD[a0+2]-0.5*uOther[a0]*m1EffD[a0+2]; 
    relKinE[3] += 0.4*m1EffD[a0+6]*uSelf[a0+7]+0.447213595499958*m1EffD[a0+2]*uSelf[a0+7]-0.4*m1EffD[a0+6]*uOther[a0+7]-0.447213595499958*m1EffD[a0+2]*uOther[a0+7]+0.4*uSelf[a0+6]*m1EffD[a0+7]-0.4*uOther[a0+6]*m1EffD[a0+7]+0.447213595499958*uSelf[a0+2]*m1EffD[a0+7]-0.447213595499958*uOther[a0+2]*m1EffD[a0+7]+0.447213595499958*m1EffD[a0+1]*uSelf[a0+6]-0.447213595499958*m1EffD[a0+1]*uOther[a0+6]+0.447213595499958*uSelf[a0+1]*m1EffD[a0+6]-0.447213595499958*uOther[a0+1]*m1EffD[a0+6]+0.4472135954999579*m1EffD[a0+3]*uSelf[a0+5]-0.4472135954999579*m1EffD[a0+3]*uOther[a0+5]+0.4472135954999579*uSelf[a0+3]*m1EffD[a0+5]-0.4472135954999579*uOther[a0+3]*m1EffD[a0+5]+0.4472135954999579*m1EffD[a0+3]*uSelf[a0+4]-0.4472135954999579*m1EffD[a0+3]*uOther[a0+4]+0.4472135954999579*uSelf[a0+3]*m1EffD[a0+4]-0.4472135954999579*uOther[a0+3]*m1EffD[a0+4]+0.5*m1EffD[a0]*uSelf[a0+3]-0.5*m1EffD[a0]*uOther[a0+3]+0.5*uSelf[a0]*m1EffD[a0+3]-0.5*uOther[a0]*m1EffD[a0+3]+0.5*m1EffD[a0+1]*uSelf[a0+2]-0.5*m1EffD[a0+1]*uOther[a0+2]+0.5*uSelf[a0+1]*m1EffD[a0+2]-0.5*uOther[a0+1]*m1EffD[a0+2]; 
    relKinE[4] += 0.4472135954999579*m1EffD[a0+7]*uSelf[a0+7]-0.4472135954999579*m1EffD[a0+7]*uOther[a0+7]+0.31943828249997*m1EffD[a0+6]*uSelf[a0+6]+0.5000000000000001*m1EffD[a0+2]*uSelf[a0+6]-0.31943828249997*m1EffD[a0+6]*uOther[a0+6]-0.5000000000000001*m1EffD[a0+2]*uOther[a0+6]+0.5000000000000001*uSelf[a0+2]*m1EffD[a0+6]-0.5000000000000001*uOther[a0+2]*m1EffD[a0+6]+0.31943828249997*m1EffD[a0+4]*uSelf[a0+4]+0.5*m1EffD[a0]*uSelf[a0+4]-0.31943828249997*m1EffD[a0+4]*uOther[a0+4]-0.5*m1EffD[a0]*uOther[a0+4]+0.5*uSelf[a0]*m1EffD[a0+4]-0.5*uOther[a0]*m1EffD[a0+4]+0.4472135954999579*m1EffD[a0+3]*uSelf[a0+3]-0.4472135954999579*m1EffD[a0+3]*uOther[a0+3]+0.4472135954999579*m1EffD[a0+1]*uSelf[a0+1]-0.4472135954999579*m1EffD[a0+1]*uOther[a0+1]; 
    relKinE[5] += 0.31943828249997*m1EffD[a0+7]*uSelf[a0+7]+0.5000000000000001*m1EffD[a0+1]*uSelf[a0+7]-0.31943828249997*m1EffD[a0+7]*uOther[a0+7]-0.5000000000000001*m1EffD[a0+1]*uOther[a0+7]+0.5000000000000001*uSelf[a0+1]*m1EffD[a0+7]-0.5000000000000001*uOther[a0+1]*m1EffD[a0+7]+0.4472135954999579*m1EffD[a0+6]*uSelf[a0+6]-0.4472135954999579*m1EffD[a0+6]*uOther[a0+6]+0.31943828249997*m1EffD[a0+5]*uSelf[a0+5]+0.5*m1EffD[a0]*uSelf[a0+5]-0.31943828249997*m1EffD[a0+5]*uOther[a0+5]-0.5*m1EffD[a0]*uOther[a0+5]+0.5*uSelf[a0]*m1EffD[a0+5]-0.5*uOther[a0]*m1EffD[a0+5]+0.4472135954999579*m1EffD[a0+3]*uSelf[a0+3]-0.4472135954999579*m1EffD[a0+3]*uOther[a0+3]+0.4472135954999579*m1EffD[a0+2]*uSelf[a0+2]-0.4472135954999579*m1EffD[a0+2]*uOther[a0+2]; 
    relKinE[6] += 0.4*m1EffD[a0+3]*uSelf[a0+7]-0.4*m1EffD[a0+3]*uOther[a0+7]+0.4*uSelf[a0+3]*m1EffD[a0+7]-0.4*uOther[a0+3]*m1EffD[a0+7]+0.4472135954999579*m1EffD[a0+5]*uSelf[a0+6]+0.31943828249997*m1EffD[a0+4]*uSelf[a0+6]+0.5*m1EffD[a0]*uSelf[a0+6]-0.4472135954999579*m1EffD[a0+5]*uOther[a0+6]-0.31943828249997*m1EffD[a0+4]*uOther[a0+6]-0.5*m1EffD[a0]*uOther[a0+6]+0.4472135954999579*uSelf[a0+5]*m1EffD[a0+6]-0.4472135954999579*uOther[a0+5]*m1EffD[a0+6]+0.31943828249997*uSelf[a0+4]*m1EffD[a0+6]-0.31943828249997*uOther[a0+4]*m1EffD[a0+6]+0.5*uSelf[a0]*m1EffD[a0+6]-0.5*uOther[a0]*m1EffD[a0+6]+0.5000000000000001*m1EffD[a0+2]*uSelf[a0+4]-0.5000000000000001*m1EffD[a0+2]*uOther[a0+4]+0.5000000000000001*uSelf[a0+2]*m1EffD[a0+4]-0.5000000000000001*uOther[a0+2]*m1EffD[a0+4]+0.447213595499958*m1EffD[a0+1]*uSelf[a0+3]-0.447213595499958*m1EffD[a0+1]*uOther[a0+3]+0.447213595499958*uSelf[a0+1]*m1EffD[a0+3]-0.447213595499958*uOther[a0+1]*m1EffD[a0+3]; 
    relKinE[7] += 0.31943828249997*m1EffD[a0+5]*uSelf[a0+7]+0.4472135954999579*m1EffD[a0+4]*uSelf[a0+7]+0.5*m1EffD[a0]*uSelf[a0+7]-0.31943828249997*m1EffD[a0+5]*uOther[a0+7]-0.4472135954999579*m1EffD[a0+4]*uOther[a0+7]-0.5*m1EffD[a0]*uOther[a0+7]+0.31943828249997*uSelf[a0+5]*m1EffD[a0+7]-0.31943828249997*uOther[a0+5]*m1EffD[a0+7]+0.4472135954999579*uSelf[a0+4]*m1EffD[a0+7]-0.4472135954999579*uOther[a0+4]*m1EffD[a0+7]+0.5*uSelf[a0]*m1EffD[a0+7]-0.5*uOther[a0]*m1EffD[a0+7]+0.4*m1EffD[a0+3]*uSelf[a0+6]-0.4*m1EffD[a0+3]*uOther[a0+6]+0.4*uSelf[a0+3]*m1EffD[a0+6]-0.4*uOther[a0+3]*m1EffD[a0+6]+0.5000000000000001*m1EffD[a0+1]*uSelf[a0+5]-0.5000000000000001*m1EffD[a0+1]*uOther[a0+5]+0.5000000000000001*uSelf[a0+1]*m1EffD[a0+5]-0.5000000000000001*uOther[a0+1]*m1EffD[a0+5]+0.447213595499958*m1EffD[a0+2]*uSelf[a0+3]-0.447213595499958*m1EffD[a0+2]*uOther[a0+3]+0.447213595499958*uSelf[a0+2]*m1EffD[a0+3]-0.447213595499958*uOther[a0+2]*m1EffD[a0+3]; 
  } 
 
  // Divide m0Other*(m2Self-kinESelf) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Other and m2Self-uSelf.m1Self. 
  double m0OtherThESelf[8]; 
  m0OtherThESelf[0] = 0.5*m0rOther[7]*m2rSelf[7]-0.5*kinESelf[7]*m0rOther[7]+0.5*m0rOther[6]*m2rSelf[6]-0.5*kinESelf[6]*m0rOther[6]+0.5*m0rOther[5]*m2rSelf[5]-0.5*kinESelf[5]*m0rOther[5]+0.5*m0rOther[4]*m2rSelf[4]-0.5*kinESelf[4]*m0rOther[4]+0.5*m0rOther[3]*m2rSelf[3]-0.5*kinESelf[3]*m0rOther[3]+0.5*m0rOther[2]*m2rSelf[2]-0.5*kinESelf[2]*m0rOther[2]+0.5*m0rOther[1]*m2rSelf[1]-0.5*kinESelf[1]*m0rOther[1]+0.5*m0rOther[0]*m2rSelf[0]-0.5*kinESelf[0]*m0rOther[0]; 
  m0OtherThESelf[1] = 0.5000000000000001*m0rOther[5]*m2rSelf[7]+0.5000000000000001*m2rSelf[5]*m0rOther[7]-0.5000000000000001*kinESelf[5]*m0rOther[7]-0.5000000000000001*m0rOther[5]*kinESelf[7]+0.447213595499958*m0rOther[3]*m2rSelf[6]+0.447213595499958*m2rSelf[3]*m0rOther[6]-0.447213595499958*kinESelf[3]*m0rOther[6]-0.447213595499958*m0rOther[3]*kinESelf[6]+0.4472135954999579*m0rOther[1]*m2rSelf[4]+0.4472135954999579*m2rSelf[1]*m0rOther[4]-0.4472135954999579*kinESelf[1]*m0rOther[4]-0.4472135954999579*m0rOther[1]*kinESelf[4]+0.5*m0rOther[2]*m2rSelf[3]+0.5*m2rSelf[2]*m0rOther[3]-0.5*kinESelf[2]*m0rOther[3]-0.5*m0rOther[2]*kinESelf[3]+0.5*m0rOther[0]*m2rSelf[1]+0.5*m2rSelf[0]*m0rOther[1]-0.5*kinESelf[0]*m0rOther[1]-0.5*m0rOther[0]*kinESelf[1]; 
  m0OtherThESelf[2] = 0.447213595499958*m0rOther[3]*m2rSelf[7]+0.447213595499958*m2rSelf[3]*m0rOther[7]-0.447213595499958*kinESelf[3]*m0rOther[7]-0.447213595499958*m0rOther[3]*kinESelf[7]+0.5000000000000001*m0rOther[4]*m2rSelf[6]+0.5000000000000001*m2rSelf[4]*m0rOther[6]-0.5000000000000001*kinESelf[4]*m0rOther[6]-0.5000000000000001*m0rOther[4]*kinESelf[6]+0.4472135954999579*m0rOther[2]*m2rSelf[5]+0.4472135954999579*m2rSelf[2]*m0rOther[5]-0.4472135954999579*kinESelf[2]*m0rOther[5]-0.4472135954999579*m0rOther[2]*kinESelf[5]+0.5*m0rOther[1]*m2rSelf[3]+0.5*m2rSelf[1]*m0rOther[3]-0.5*kinESelf[1]*m0rOther[3]-0.5*m0rOther[1]*kinESelf[3]+0.5*m0rOther[0]*m2rSelf[2]+0.5*m2rSelf[0]*m0rOther[2]-0.5*kinESelf[0]*m0rOther[2]-0.5*m0rOther[0]*kinESelf[2]; 
  m0OtherThESelf[3] = 0.4*m0rOther[6]*m2rSelf[7]+0.447213595499958*m0rOther[2]*m2rSelf[7]+0.4*m2rSelf[6]*m0rOther[7]-0.4*kinESelf[6]*m0rOther[7]+0.447213595499958*m2rSelf[2]*m0rOther[7]-0.447213595499958*kinESelf[2]*m0rOther[7]-0.4*m0rOther[6]*kinESelf[7]-0.447213595499958*m0rOther[2]*kinESelf[7]+0.447213595499958*m0rOther[1]*m2rSelf[6]+0.447213595499958*m2rSelf[1]*m0rOther[6]-0.447213595499958*kinESelf[1]*m0rOther[6]-0.447213595499958*m0rOther[1]*kinESelf[6]+0.4472135954999579*m0rOther[3]*m2rSelf[5]+0.4472135954999579*m2rSelf[3]*m0rOther[5]-0.4472135954999579*kinESelf[3]*m0rOther[5]-0.4472135954999579*m0rOther[3]*kinESelf[5]+0.4472135954999579*m0rOther[3]*m2rSelf[4]+0.4472135954999579*m2rSelf[3]*m0rOther[4]-0.4472135954999579*kinESelf[3]*m0rOther[4]-0.4472135954999579*m0rOther[3]*kinESelf[4]+0.5*m0rOther[0]*m2rSelf[3]+0.5*m2rSelf[0]*m0rOther[3]-0.5*kinESelf[0]*m0rOther[3]-0.5*m0rOther[0]*kinESelf[3]+0.5*m0rOther[1]*m2rSelf[2]+0.5*m2rSelf[1]*m0rOther[2]-0.5*kinESelf[1]*m0rOther[2]-0.5*m0rOther[1]*kinESelf[2]; 
  m0OtherThESelf[4] = 0.4472135954999579*m0rOther[7]*m2rSelf[7]-0.4472135954999579*kinESelf[7]*m0rOther[7]+0.31943828249997*m0rOther[6]*m2rSelf[6]+0.5000000000000001*m0rOther[2]*m2rSelf[6]-0.31943828249997*kinESelf[6]*m0rOther[6]+0.5000000000000001*m2rSelf[2]*m0rOther[6]-0.5000000000000001*kinESelf[2]*m0rOther[6]-0.5000000000000001*m0rOther[2]*kinESelf[6]+0.31943828249997*m0rOther[4]*m2rSelf[4]+0.5*m0rOther[0]*m2rSelf[4]-0.31943828249997*kinESelf[4]*m0rOther[4]+0.5*m2rSelf[0]*m0rOther[4]-0.5*kinESelf[0]*m0rOther[4]-0.5*m0rOther[0]*kinESelf[4]+0.4472135954999579*m0rOther[3]*m2rSelf[3]-0.4472135954999579*kinESelf[3]*m0rOther[3]+0.4472135954999579*m0rOther[1]*m2rSelf[1]-0.4472135954999579*kinESelf[1]*m0rOther[1]; 
  m0OtherThESelf[5] = 0.31943828249997*m0rOther[7]*m2rSelf[7]+0.5000000000000001*m0rOther[1]*m2rSelf[7]-0.31943828249997*kinESelf[7]*m0rOther[7]+0.5000000000000001*m2rSelf[1]*m0rOther[7]-0.5000000000000001*kinESelf[1]*m0rOther[7]-0.5000000000000001*m0rOther[1]*kinESelf[7]+0.4472135954999579*m0rOther[6]*m2rSelf[6]-0.4472135954999579*kinESelf[6]*m0rOther[6]+0.31943828249997*m0rOther[5]*m2rSelf[5]+0.5*m0rOther[0]*m2rSelf[5]-0.31943828249997*kinESelf[5]*m0rOther[5]+0.5*m2rSelf[0]*m0rOther[5]-0.5*kinESelf[0]*m0rOther[5]-0.5*m0rOther[0]*kinESelf[5]+0.4472135954999579*m0rOther[3]*m2rSelf[3]-0.4472135954999579*kinESelf[3]*m0rOther[3]+0.4472135954999579*m0rOther[2]*m2rSelf[2]-0.4472135954999579*kinESelf[2]*m0rOther[2]; 
  m0OtherThESelf[6] = 0.4*m0rOther[3]*m2rSelf[7]+0.4*m2rSelf[3]*m0rOther[7]-0.4*kinESelf[3]*m0rOther[7]-0.4*m0rOther[3]*kinESelf[7]+0.4472135954999579*m0rOther[5]*m2rSelf[6]+0.31943828249997*m0rOther[4]*m2rSelf[6]+0.5*m0rOther[0]*m2rSelf[6]+0.4472135954999579*m2rSelf[5]*m0rOther[6]-0.4472135954999579*kinESelf[5]*m0rOther[6]+0.31943828249997*m2rSelf[4]*m0rOther[6]-0.31943828249997*kinESelf[4]*m0rOther[6]+0.5*m2rSelf[0]*m0rOther[6]-0.5*kinESelf[0]*m0rOther[6]-0.4472135954999579*m0rOther[5]*kinESelf[6]-0.31943828249997*m0rOther[4]*kinESelf[6]-0.5*m0rOther[0]*kinESelf[6]+0.5000000000000001*m0rOther[2]*m2rSelf[4]+0.5000000000000001*m2rSelf[2]*m0rOther[4]-0.5000000000000001*kinESelf[2]*m0rOther[4]-0.5000000000000001*m0rOther[2]*kinESelf[4]+0.447213595499958*m0rOther[1]*m2rSelf[3]+0.447213595499958*m2rSelf[1]*m0rOther[3]-0.447213595499958*kinESelf[1]*m0rOther[3]-0.447213595499958*m0rOther[1]*kinESelf[3]; 
  m0OtherThESelf[7] = 0.31943828249997*m0rOther[5]*m2rSelf[7]+0.4472135954999579*m0rOther[4]*m2rSelf[7]+0.5*m0rOther[0]*m2rSelf[7]+0.31943828249997*m2rSelf[5]*m0rOther[7]-0.31943828249997*kinESelf[5]*m0rOther[7]+0.4472135954999579*m2rSelf[4]*m0rOther[7]-0.4472135954999579*kinESelf[4]*m0rOther[7]+0.5*m2rSelf[0]*m0rOther[7]-0.5*kinESelf[0]*m0rOther[7]-0.31943828249997*m0rOther[5]*kinESelf[7]-0.4472135954999579*m0rOther[4]*kinESelf[7]-0.5*m0rOther[0]*kinESelf[7]+0.4*m0rOther[3]*m2rSelf[6]+0.4*m2rSelf[3]*m0rOther[6]-0.4*kinESelf[3]*m0rOther[6]-0.4*m0rOther[3]*kinESelf[6]+0.5000000000000001*m0rOther[1]*m2rSelf[5]+0.5000000000000001*m2rSelf[1]*m0rOther[5]-0.5000000000000001*kinESelf[1]*m0rOther[5]-0.5000000000000001*m0rOther[1]*kinESelf[5]+0.447213595499958*m0rOther[2]*m2rSelf[3]+0.447213595499958*m2rSelf[2]*m0rOther[3]-0.447213595499958*kinESelf[2]*m0rOther[3]-0.447213595499958*m0rOther[2]*kinESelf[3]; 
  dataDiv->BEV_S << m0OtherThESelf[0],m0OtherThESelf[1],m0OtherThESelf[2],m0OtherThESelf[3],m0OtherThESelf[4],m0OtherThESelf[5],m0OtherThESelf[6],m0OtherThESelf[7]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthSelf[8]; 
  Eigen::Map<VectorXd>(effEthSelf,8,1) = dataDiv->u_S; 
 
  // Divide m0Self*(m2Other-kinEOther) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Self and m2Other-uOther.m1Other. 
  double m0SelfThEOther[8]; 
  m0SelfThEOther[0] = 0.5*m0rSelf[7]*m2rOther[7]-0.5*kinEOther[7]*m0rSelf[7]+0.5*m0rSelf[6]*m2rOther[6]-0.5*kinEOther[6]*m0rSelf[6]+0.5*m0rSelf[5]*m2rOther[5]-0.5*kinEOther[5]*m0rSelf[5]+0.5*m0rSelf[4]*m2rOther[4]-0.5*kinEOther[4]*m0rSelf[4]+0.5*m0rSelf[3]*m2rOther[3]-0.5*kinEOther[3]*m0rSelf[3]+0.5*m0rSelf[2]*m2rOther[2]-0.5*kinEOther[2]*m0rSelf[2]+0.5*m0rSelf[1]*m2rOther[1]-0.5*kinEOther[1]*m0rSelf[1]+0.5*m0rSelf[0]*m2rOther[0]-0.5*kinEOther[0]*m0rSelf[0]; 
  m0SelfThEOther[1] = 0.5000000000000001*m0rSelf[5]*m2rOther[7]+0.5000000000000001*m2rOther[5]*m0rSelf[7]-0.5000000000000001*kinEOther[5]*m0rSelf[7]-0.5000000000000001*m0rSelf[5]*kinEOther[7]+0.447213595499958*m0rSelf[3]*m2rOther[6]+0.447213595499958*m2rOther[3]*m0rSelf[6]-0.447213595499958*kinEOther[3]*m0rSelf[6]-0.447213595499958*m0rSelf[3]*kinEOther[6]+0.4472135954999579*m0rSelf[1]*m2rOther[4]+0.4472135954999579*m2rOther[1]*m0rSelf[4]-0.4472135954999579*kinEOther[1]*m0rSelf[4]-0.4472135954999579*m0rSelf[1]*kinEOther[4]+0.5*m0rSelf[2]*m2rOther[3]+0.5*m2rOther[2]*m0rSelf[3]-0.5*kinEOther[2]*m0rSelf[3]-0.5*m0rSelf[2]*kinEOther[3]+0.5*m0rSelf[0]*m2rOther[1]+0.5*m2rOther[0]*m0rSelf[1]-0.5*kinEOther[0]*m0rSelf[1]-0.5*m0rSelf[0]*kinEOther[1]; 
  m0SelfThEOther[2] = 0.447213595499958*m0rSelf[3]*m2rOther[7]+0.447213595499958*m2rOther[3]*m0rSelf[7]-0.447213595499958*kinEOther[3]*m0rSelf[7]-0.447213595499958*m0rSelf[3]*kinEOther[7]+0.5000000000000001*m0rSelf[4]*m2rOther[6]+0.5000000000000001*m2rOther[4]*m0rSelf[6]-0.5000000000000001*kinEOther[4]*m0rSelf[6]-0.5000000000000001*m0rSelf[4]*kinEOther[6]+0.4472135954999579*m0rSelf[2]*m2rOther[5]+0.4472135954999579*m2rOther[2]*m0rSelf[5]-0.4472135954999579*kinEOther[2]*m0rSelf[5]-0.4472135954999579*m0rSelf[2]*kinEOther[5]+0.5*m0rSelf[1]*m2rOther[3]+0.5*m2rOther[1]*m0rSelf[3]-0.5*kinEOther[1]*m0rSelf[3]-0.5*m0rSelf[1]*kinEOther[3]+0.5*m0rSelf[0]*m2rOther[2]+0.5*m2rOther[0]*m0rSelf[2]-0.5*kinEOther[0]*m0rSelf[2]-0.5*m0rSelf[0]*kinEOther[2]; 
  m0SelfThEOther[3] = 0.4*m0rSelf[6]*m2rOther[7]+0.447213595499958*m0rSelf[2]*m2rOther[7]+0.4*m2rOther[6]*m0rSelf[7]-0.4*kinEOther[6]*m0rSelf[7]+0.447213595499958*m2rOther[2]*m0rSelf[7]-0.447213595499958*kinEOther[2]*m0rSelf[7]-0.4*m0rSelf[6]*kinEOther[7]-0.447213595499958*m0rSelf[2]*kinEOther[7]+0.447213595499958*m0rSelf[1]*m2rOther[6]+0.447213595499958*m2rOther[1]*m0rSelf[6]-0.447213595499958*kinEOther[1]*m0rSelf[6]-0.447213595499958*m0rSelf[1]*kinEOther[6]+0.4472135954999579*m0rSelf[3]*m2rOther[5]+0.4472135954999579*m2rOther[3]*m0rSelf[5]-0.4472135954999579*kinEOther[3]*m0rSelf[5]-0.4472135954999579*m0rSelf[3]*kinEOther[5]+0.4472135954999579*m0rSelf[3]*m2rOther[4]+0.4472135954999579*m2rOther[3]*m0rSelf[4]-0.4472135954999579*kinEOther[3]*m0rSelf[4]-0.4472135954999579*m0rSelf[3]*kinEOther[4]+0.5*m0rSelf[0]*m2rOther[3]+0.5*m2rOther[0]*m0rSelf[3]-0.5*kinEOther[0]*m0rSelf[3]-0.5*m0rSelf[0]*kinEOther[3]+0.5*m0rSelf[1]*m2rOther[2]+0.5*m2rOther[1]*m0rSelf[2]-0.5*kinEOther[1]*m0rSelf[2]-0.5*m0rSelf[1]*kinEOther[2]; 
  m0SelfThEOther[4] = 0.4472135954999579*m0rSelf[7]*m2rOther[7]-0.4472135954999579*kinEOther[7]*m0rSelf[7]+0.31943828249997*m0rSelf[6]*m2rOther[6]+0.5000000000000001*m0rSelf[2]*m2rOther[6]-0.31943828249997*kinEOther[6]*m0rSelf[6]+0.5000000000000001*m2rOther[2]*m0rSelf[6]-0.5000000000000001*kinEOther[2]*m0rSelf[6]-0.5000000000000001*m0rSelf[2]*kinEOther[6]+0.31943828249997*m0rSelf[4]*m2rOther[4]+0.5*m0rSelf[0]*m2rOther[4]-0.31943828249997*kinEOther[4]*m0rSelf[4]+0.5*m2rOther[0]*m0rSelf[4]-0.5*kinEOther[0]*m0rSelf[4]-0.5*m0rSelf[0]*kinEOther[4]+0.4472135954999579*m0rSelf[3]*m2rOther[3]-0.4472135954999579*kinEOther[3]*m0rSelf[3]+0.4472135954999579*m0rSelf[1]*m2rOther[1]-0.4472135954999579*kinEOther[1]*m0rSelf[1]; 
  m0SelfThEOther[5] = 0.31943828249997*m0rSelf[7]*m2rOther[7]+0.5000000000000001*m0rSelf[1]*m2rOther[7]-0.31943828249997*kinEOther[7]*m0rSelf[7]+0.5000000000000001*m2rOther[1]*m0rSelf[7]-0.5000000000000001*kinEOther[1]*m0rSelf[7]-0.5000000000000001*m0rSelf[1]*kinEOther[7]+0.4472135954999579*m0rSelf[6]*m2rOther[6]-0.4472135954999579*kinEOther[6]*m0rSelf[6]+0.31943828249997*m0rSelf[5]*m2rOther[5]+0.5*m0rSelf[0]*m2rOther[5]-0.31943828249997*kinEOther[5]*m0rSelf[5]+0.5*m2rOther[0]*m0rSelf[5]-0.5*kinEOther[0]*m0rSelf[5]-0.5*m0rSelf[0]*kinEOther[5]+0.4472135954999579*m0rSelf[3]*m2rOther[3]-0.4472135954999579*kinEOther[3]*m0rSelf[3]+0.4472135954999579*m0rSelf[2]*m2rOther[2]-0.4472135954999579*kinEOther[2]*m0rSelf[2]; 
  m0SelfThEOther[6] = 0.4*m0rSelf[3]*m2rOther[7]+0.4*m2rOther[3]*m0rSelf[7]-0.4*kinEOther[3]*m0rSelf[7]-0.4*m0rSelf[3]*kinEOther[7]+0.4472135954999579*m0rSelf[5]*m2rOther[6]+0.31943828249997*m0rSelf[4]*m2rOther[6]+0.5*m0rSelf[0]*m2rOther[6]+0.4472135954999579*m2rOther[5]*m0rSelf[6]-0.4472135954999579*kinEOther[5]*m0rSelf[6]+0.31943828249997*m2rOther[4]*m0rSelf[6]-0.31943828249997*kinEOther[4]*m0rSelf[6]+0.5*m2rOther[0]*m0rSelf[6]-0.5*kinEOther[0]*m0rSelf[6]-0.4472135954999579*m0rSelf[5]*kinEOther[6]-0.31943828249997*m0rSelf[4]*kinEOther[6]-0.5*m0rSelf[0]*kinEOther[6]+0.5000000000000001*m0rSelf[2]*m2rOther[4]+0.5000000000000001*m2rOther[2]*m0rSelf[4]-0.5000000000000001*kinEOther[2]*m0rSelf[4]-0.5000000000000001*m0rSelf[2]*kinEOther[4]+0.447213595499958*m0rSelf[1]*m2rOther[3]+0.447213595499958*m2rOther[1]*m0rSelf[3]-0.447213595499958*kinEOther[1]*m0rSelf[3]-0.447213595499958*m0rSelf[1]*kinEOther[3]; 
  m0SelfThEOther[7] = 0.31943828249997*m0rSelf[5]*m2rOther[7]+0.4472135954999579*m0rSelf[4]*m2rOther[7]+0.5*m0rSelf[0]*m2rOther[7]+0.31943828249997*m2rOther[5]*m0rSelf[7]-0.31943828249997*kinEOther[5]*m0rSelf[7]+0.4472135954999579*m2rOther[4]*m0rSelf[7]-0.4472135954999579*kinEOther[4]*m0rSelf[7]+0.5*m2rOther[0]*m0rSelf[7]-0.5*kinEOther[0]*m0rSelf[7]-0.31943828249997*m0rSelf[5]*kinEOther[7]-0.4472135954999579*m0rSelf[4]*kinEOther[7]-0.5*m0rSelf[0]*kinEOther[7]+0.4*m0rSelf[3]*m2rOther[6]+0.4*m2rOther[3]*m0rSelf[6]-0.4*kinEOther[3]*m0rSelf[6]-0.4*m0rSelf[3]*kinEOther[6]+0.5000000000000001*m0rSelf[1]*m2rOther[5]+0.5000000000000001*m2rOther[1]*m0rSelf[5]-0.5000000000000001*kinEOther[1]*m0rSelf[5]-0.5000000000000001*m0rSelf[1]*kinEOther[5]+0.447213595499958*m0rSelf[2]*m2rOther[3]+0.447213595499958*m2rOther[2]*m0rSelf[3]-0.447213595499958*kinEOther[2]*m0rSelf[3]-0.447213595499958*m0rSelf[2]*kinEOther[3]; 
  dataDiv->BEV_S << m0SelfThEOther[0],m0SelfThEOther[1],m0SelfThEOther[2],m0SelfThEOther[3],m0SelfThEOther[4],m0SelfThEOther[5],m0SelfThEOther[6],m0SelfThEOther[7]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthOther[8]; 
  Eigen::Map<VectorXd>(effEthOther,8,1) = dataDiv->u_S; 
 
  double m2Relax[8]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(1.0*relKinE[0]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[0]*mSelf)/(mSelf+mOther)+(1.0*relKinE[0]*mOther)/(mSelf+mOther)+(2.0*effEthOther[0]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2rOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(1.0*relKinE[1]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[1]*mSelf)/(mSelf+mOther)+(1.0*relKinE[1]*mOther)/(mSelf+mOther)+(2.0*effEthOther[1]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2rOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(1.0*relKinE[2]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[2]*mSelf)/(mSelf+mOther)+(1.0*relKinE[2]*mOther)/(mSelf+mOther)+(2.0*effEthOther[2]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2rOther[2])*mnuOther; 
  m2Relax[3] = betaGreenep1*((-(1.0*relKinE[3]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[3]*mSelf)/(mSelf+mOther)+(1.0*relKinE[3]*mOther)/(mSelf+mOther)+(2.0*effEthOther[3]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[3]-1.0*kinESelf[3])*mnuSelf+(kinEOther[3]-1.0*m2rOther[3])*mnuOther; 
  m2Relax[4] = betaGreenep1*((-(1.0*relKinE[4]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[4]*mSelf)/(mSelf+mOther)+(1.0*relKinE[4]*mOther)/(mSelf+mOther)+(2.0*effEthOther[4]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[4]-1.0*kinESelf[4])*mnuSelf+(kinEOther[4]-1.0*m2rOther[4])*mnuOther; 
  m2Relax[5] = betaGreenep1*((-(1.0*relKinE[5]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[5]*mSelf)/(mSelf+mOther)+(1.0*relKinE[5]*mOther)/(mSelf+mOther)+(2.0*effEthOther[5]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[5]-1.0*kinESelf[5])*mnuSelf+(kinEOther[5]-1.0*m2rOther[5])*mnuOther; 
  m2Relax[6] = betaGreenep1*((-(1.0*relKinE[6]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[6]*mSelf)/(mSelf+mOther)+(1.0*relKinE[6]*mOther)/(mSelf+mOther)+(2.0*effEthOther[6]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[6]-1.0*kinESelf[6])*mnuSelf+(kinEOther[6]-1.0*m2rOther[6])*mnuOther; 
  m2Relax[7] = betaGreenep1*((-(1.0*relKinE[7]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[7]*mSelf)/(mSelf+mOther)+(1.0*relKinE[7]*mOther)/(mSelf+mOther)+(2.0*effEthOther[7]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[7]-1.0*kinESelf[7])*mnuSelf+(kinEOther[7]-1.0*m2rOther[7])*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<8,16>(32,8).setZero(); 
  data->AEM_S.block<16,8>(40,0).setZero(); 
  data->AEM_S.block<8,8>(40,16).setZero(); 
  data->AEM_S.block<8,8>(48,8).setZero(); 
  data->AEM_S.block<8,16>(32,40).setZero(); 
  data->AEM_S.block<16,8>(40,32).setZero(); 
  data->AEM_S.block<8,8>(40,48).setZero(); 
  data->AEM_S.block<8,8>(48,40).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM1sum[10],mnuM1sum[11],mnuM1sum[12],mnuM1sum[13],mnuM1sum[14],mnuM1sum[15],mnuM1sum[16],mnuM1sum[17],mnuM1sum[18],mnuM1sum[19],mnuM1sum[20],mnuM1sum[21],mnuM1sum[22],mnuM1sum[23],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],mnuM2sum[6],mnuM2sum[7],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m1Relax[9],m1Relax[10],m1Relax[11],m1Relax[12],m1Relax[13],m1Relax[14],m1Relax[15],m1Relax[16],m1Relax[17],m1Relax[18],m1Relax[19],m1Relax[20],m1Relax[21],m1Relax[22],m1Relax[23],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5],m2Relax[6],m2Relax[7]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,24,1) = data->u_S.segment<24>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,8,1) = data->u_S.segment<8>(24); 
 
  Eigen::Map<VectorXd>(uCrossOther,24,1) = data->u_S.segment<24>(32); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,8,1) = data->u_S.segment<8>(56); 
 
} 
 
