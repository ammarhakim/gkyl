#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmCrossPrimMoments2x3vMax_P1(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[9]; 
  double m2rSelf[3]; 
  double m0SrSelf[3]; 
  double m1SrSelf[9]; 
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
    m1rSelf[6] = m1Self[6]; 
    m1rSelf[7] = 0.0; 
    m1rSelf[8] = 0.0; 
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
    m1SrSelf[6] = m1SSelf[6]; 
    m1SrSelf[7] = 0.0; 
    m1SrSelf[8] = 0.0; 
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
    m1rSelf[6] = m1Self[6]; 
    m1rSelf[7] = m1Self[7]; 
    m1rSelf[8] = m1Self[8]; 
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
    m1SrSelf[6] = m1SSelf[6]; 
    m1SrSelf[7] = m1SSelf[7]; 
    m1SrSelf[8] = m1SSelf[8]; 
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
  double m1rOther[9]; 
  double m2rOther[3]; 
  double m0SrOther[3]; 
  double m1SrOther[9]; 
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
    m1rOther[6] = m1Other[6]; 
    m1rOther[7] = 0.0; 
    m1rOther[8] = 0.0; 
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
    m1SrOther[6] = m1SOther[6]; 
    m1SrOther[7] = 0.0; 
    m1SrOther[8] = 0.0; 
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
    m1rOther[6] = m1Other[6]; 
    m1rOther[7] = m1Other[7]; 
    m1rOther[8] = m1Other[8]; 
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
    m1SrOther[6] = m1SOther[6]; 
    m1SrOther[7] = m1SOther[7]; 
    m1SrOther[8] = m1SOther[8]; 
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
    m2SrOther[2] = m2SOther[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(24,24); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
  double mnuM1sum[9]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<9; vd++) 
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
  data->AEM_S(0,9) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,10) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,11) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,9) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,10) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,9) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,11) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,12) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,13) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,14) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,12) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,13) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,12) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,14) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,21) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,22) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,23) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,21) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,22) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,21) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,23) = -0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(9,0) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(9,1) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(9,2) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(10,0) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(10,1) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(11,0) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(11,2) = 0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(9,12) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(9,13) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(9,14) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(10,12) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(10,13) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(11,12) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(11,14) = 0.5*m1SrOther[0]*mnuOther; 
 
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
  data->AEM_S(3,9) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,10) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(3,11) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(4,9) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,10) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,9) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,11) = -0.5*cMSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(3,15) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,16) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,17) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(4,15) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(4,16) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,15) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(5,17) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(3,21) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,22) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(3,23) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(4,21) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,22) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(5,21) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,23) = -0.5*cMOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(9,3) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(9,4) = 0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(9,5) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(10,3) = 0.5*m1SrSelf[4]*mnuSelf; 
  data->AEM_S(10,4) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(11,3) = 0.5*m1SrSelf[5]*mnuSelf; 
  data->AEM_S(11,5) = 0.5*m1SrSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(9,15) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(9,16) = 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(9,17) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(10,15) = 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(10,16) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(11,15) = 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(11,17) = 0.5*m1SrOther[3]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(6,6) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,6) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,7) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,6) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,8) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(6,9) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,10) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(6,11) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(7,9) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,10) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(8,9) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,11) = -0.5*cMSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(6,18) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,19) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(6,20) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(7,18) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,19) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,18) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,20) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(6,21) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,22) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(6,23) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(7,21) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,22) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(8,21) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(8,23) = -0.5*cMOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfZ and uCrossSelfZ ... // 
  data->AEM_S(9,6) = 0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(9,7) = 0.5*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(9,8) = 0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(10,6) = 0.5*m1SrSelf[7]*mnuSelf; 
  data->AEM_S(10,7) = 0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(11,6) = 0.5*m1SrSelf[8]*mnuSelf; 
  data->AEM_S(11,8) = 0.5*m1SrSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherZ and uCrossOtherZ ... // 
  data->AEM_S(9,18) = 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(9,19) = 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(9,20) = 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(10,18) = 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(10,19) = 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(11,18) = 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(11,20) = 0.5*m1SrOther[6]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(9,9) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(9,10) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(9,11) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(10,9) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(10,10) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(11,9) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(11,11) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(9,21) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(9,22) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(9,23) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(10,21) = 0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(10,22) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(11,21) = 0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(11,23) = 0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[3]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2SrSelf[0]*mnuSelf+m2SrOther[0]*mnuOther; 
  mnuM2sum[1] = m2SrSelf[1]*mnuSelf+m2SrOther[1]*mnuOther; 
  mnuM2sum[2] = m2SrSelf[2]*mnuSelf+m2SrOther[2]*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<3,6>(0,3).setZero(); 
  data->AEM_S.block<6,3>(3,0).setZero(); 
  data->AEM_S.block<3,3>(3,6).setZero(); 
  data->AEM_S.block<3,3>(6,3).setZero(); 
  data->AEM_S.block<3,6>(0,15).setZero(); 
  data->AEM_S.block<6,3>(3,12).setZero(); 
  data->AEM_S.block<3,3>(3,18).setZero(); 
  data->AEM_S.block<3,3>(6,15).setZero(); 
 
  double m1Relax[9]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<9; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(12,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(14,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,2) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(12,9) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(12,10) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(12,11) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(13,9) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(13,10) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(14,9) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(14,11) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(12,12) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(12,13) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(12,14) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(13,12) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(13,13) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(14,12) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(14,14) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(12,21) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(12,22) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(12,23) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(13,21) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(13,22) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(14,21) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(14,23) = 0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(21,0) = (-0.25*m0rSelf[2]*uSelf[2]*mnuSelf)-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(21,1) = (-0.25*m0rSelf[0]*uSelf[1]*mnuSelf)+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,2) = (-0.25*m0rSelf[0]*uSelf[2]*mnuSelf)+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,0) = (-0.25*m0rSelf[0]*uSelf[1]*mnuSelf)+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,1) = (-0.25*m0rSelf[2]*uSelf[2]*mnuSelf)-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(22,2) = (-0.25*m0rSelf[1]*uSelf[2]*mnuSelf)-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,0) = (-0.25*m0rSelf[0]*uSelf[2]*mnuSelf)+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,1) = (-0.25*m0rSelf[1]*uSelf[2]*mnuSelf)-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,2) = (-0.45*m0rSelf[2]*uSelf[2]*mnuSelf)-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(21,12) = 0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(21,13) = 0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(21,14) = 0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(22,12) = 0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(22,13) = 0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(22,14) = 0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(23,12) = 0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(23,13) = 0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(23,14) = 0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(15,3) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(15,4) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,5) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,4) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(17,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,5) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(15,9) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(15,10) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(15,11) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(16,9) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(16,10) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(17,9) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(17,11) = -0.5*cMSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(15,15) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(15,16) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(15,17) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(16,15) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(16,16) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(17,15) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(17,17) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(15,21) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(15,22) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(15,23) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(16,21) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(16,22) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(17,21) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(17,23) = 0.5*cMOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(21,3) = (-0.25*m0rSelf[2]*uSelf[5]*mnuSelf)-0.25*m0rSelf[1]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(21,4) = (-0.25*m0rSelf[0]*uSelf[4]*mnuSelf)+0.5*m1SrSelf[4]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(21,5) = (-0.25*m0rSelf[0]*uSelf[5]*mnuSelf)+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf; 
  data->AEM_S(22,3) = (-0.25*m0rSelf[0]*uSelf[4]*mnuSelf)+0.5*m1SrSelf[4]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(22,4) = (-0.25*m0rSelf[2]*uSelf[5]*mnuSelf)-0.45*m0rSelf[1]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(22,5) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(23,3) = (-0.25*m0rSelf[0]*uSelf[5]*mnuSelf)+0.5*m1SrSelf[5]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf; 
  data->AEM_S(23,4) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*m0rSelf[2]*uSelf[4]*mnuSelf; 
  data->AEM_S(23,5) = (-0.45*m0rSelf[2]*uSelf[5]*mnuSelf)-0.25*m0rSelf[1]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1SrSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(21,15) = 0.25*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(21,16) = 0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(21,17) = 0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther; 
  data->AEM_S(22,15) = 0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1SrOther[4]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(22,16) = 0.25*m0rOther[2]*uOther[5]*mnuOther+0.45*m0rOther[1]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(22,17) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(23,15) = 0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1SrOther[5]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther; 
  data->AEM_S(23,16) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther; 
  data->AEM_S(23,17) = 0.45*m0rOther[2]*uOther[5]*mnuOther+0.25*m0rOther[1]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[3] += (m1rOther[3]-1.0*m1rSelf[3])*betaGreenep1*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (m1rOther[4]-1.0*m1rSelf[4])*betaGreenep1*mnuSelf+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (m1rOther[5]-1.0*m1rSelf[5])*betaGreenep1*mnuSelf+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(18,6) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,7) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(18,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,6) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,7) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,6) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,8) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(18,9) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(18,10) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(18,11) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(19,9) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(19,10) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(20,9) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(20,11) = -0.5*cMSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(18,18) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(18,19) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(18,20) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(19,18) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(19,19) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(20,18) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(20,20) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(18,21) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(18,22) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(18,23) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(19,21) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(19,22) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(20,21) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(20,23) = 0.5*cMOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfZ-uSelfZ*m0Self) and uCrossSelfZ ... // 
  data->AEM_S(21,6) = (-0.25*m0rSelf[2]*uSelf[8]*mnuSelf)-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(21,7) = (-0.25*m0rSelf[0]*uSelf[7]*mnuSelf)+0.5*m1SrSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(21,8) = (-0.25*m0rSelf[0]*uSelf[8]*mnuSelf)+0.5*m1SrSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(22,6) = (-0.25*m0rSelf[0]*uSelf[7]*mnuSelf)+0.5*m1SrSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(22,7) = (-0.25*m0rSelf[2]*uSelf[8]*mnuSelf)-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1SrSelf[6]*mnuSelf; 
  data->AEM_S(22,8) = (-0.25*m0rSelf[1]*uSelf[8]*mnuSelf)-0.25*m0rSelf[2]*uSelf[7]*mnuSelf; 
  data->AEM_S(23,6) = (-0.25*m0rSelf[0]*uSelf[8]*mnuSelf)+0.5*m1SrSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(23,7) = (-0.25*m0rSelf[1]*uSelf[8]*mnuSelf)-0.25*m0rSelf[2]*uSelf[7]*mnuSelf; 
  data->AEM_S(23,8) = (-0.45*m0rSelf[2]*uSelf[8]*mnuSelf)-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1SrSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherZ-uOtherZ*m0Other) and uCrossOtherZ ... // 
  data->AEM_S(21,18) = 0.25*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(21,19) = 0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1SrOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(21,20) = 0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1SrOther[8]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(22,18) = 0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1SrOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(22,19) = 0.25*m0rOther[2]*uOther[8]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(22,20) = 0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther; 
  data->AEM_S(23,18) = 0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1SrOther[8]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(23,19) = 0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther; 
  data->AEM_S(23,20) = 0.45*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1SrOther[6]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of momentum relaxation. 
  m1Relax[6] += (m1rOther[6]-1.0*m1rSelf[6])*betaGreenep1*mnuSelf+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (m1rOther[7]-1.0*m1rSelf[7])*betaGreenep1*mnuSelf+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
  m1Relax[8] += (m1rOther[8]-1.0*m1rSelf[8])*betaGreenep1*mnuSelf+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
 
  double ucMSelf[3]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.5*cMSelf[a0+2]*uSelf[a0+2]+0.5*cMSelf[a0+1]*uSelf[a0+1]+0.5*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.5*cMSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.5*cMSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*cMSelf[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(21,9) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(21,10) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(21,11) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(22,9) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(22,10) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(23,9) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(23,11) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  double ucMOther[3]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.5*cMOther[a0+2]*uOther[a0+2]+0.5*cMOther[a0+1]*uOther[a0+1]+0.5*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.5*cMOther[a0]*uOther[a0+1]+0.5*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.5*cMOther[a0]*uOther[a0+2]+0.5*uOther[a0]*cMOther[a0+2]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(21,21) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(21,22) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(21,23) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(22,21) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(22,22) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(23,21) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(23,23) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[3]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
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
  for (unsigned short int vd=0; vd<3; vd++) 
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
  for (unsigned short int vd=0; vd<3; vd++) 
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
  data->AEM_S.block<3,6>(12,3).setZero(); 
  data->AEM_S.block<6,3>(15,0).setZero(); 
  data->AEM_S.block<3,3>(15,6).setZero(); 
  data->AEM_S.block<3,3>(18,3).setZero(); 
  data->AEM_S.block<3,6>(12,15).setZero(); 
  data->AEM_S.block<6,3>(15,12).setZero(); 
  data->AEM_S.block<3,3>(15,18).setZero(); 
  data->AEM_S.block<3,3>(18,15).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m2Relax[0],m2Relax[1],m2Relax[2]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,9,1) = data->u_S.segment<9>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,3,1) = data->u_S.segment<3>(9); 
 
  Eigen::Map<VectorXd>(uCrossOther,9,1) = data->u_S.segment<9>(12); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,3,1) = data->u_S.segment<3>(21); 
 
} 
 
void VmCrossPrimMoments2x3vMax_P2(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[18]; 
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
    m1rSelf[12] = m1Self[12]; 
    m1rSelf[13] = 0.0; 
    m1rSelf[14] = 0.0; 
    m1rSelf[15] = 0.0; 
    m1rSelf[16] = 0.0; 
    m1rSelf[17] = 0.0; 
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
    m1rSelf[12] = m1Self[12]; 
    m1rSelf[13] = m1Self[13]; 
    m1rSelf[14] = m1Self[14]; 
    m1rSelf[15] = m1Self[15]; 
    m1rSelf[16] = m1Self[16]; 
    m1rSelf[17] = m1Self[17]; 
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
  double m1rOther[18]; 
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
    m1rOther[12] = m1Other[12]; 
    m1rOther[13] = 0.0; 
    m1rOther[14] = 0.0; 
    m1rOther[15] = 0.0; 
    m1rOther[16] = 0.0; 
    m1rOther[17] = 0.0; 
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
    m1rOther[12] = m1Other[12]; 
    m1rOther[13] = m1Other[13]; 
    m1rOther[14] = m1Other[14]; 
    m1rOther[15] = m1Other[15]; 
    m1rOther[16] = m1Other[16]; 
    m1rOther[17] = m1Other[17]; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
    m2rOther[3] = m2Other[3]; 
    m2rOther[4] = m2Other[4]; 
    m2rOther[5] = m2Other[5]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(48,48); 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
  double mnuM1sum[18]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<18; vd++) 
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
  data->AEM_S(0,18) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,19) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,20) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,21) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,22) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,23) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(1,18) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,19) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,20) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,21) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,22) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,18) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,19) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,20) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,21) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,23) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,18) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,19) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,20) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,21) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,22) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,23) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,18) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,19) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,21) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,22) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,18) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,20) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(5,21) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,23) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,24) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,25) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,26) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,27) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(0,28) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(0,29) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(1,24) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,25) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,26) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,27) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,28) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(2,24) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,25) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,26) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,27) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(2,29) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(3,24) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,25) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,26) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,27) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,28) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,29) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,24) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(4,25) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(4,27) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,28) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,24) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(5,26) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(5,27) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(5,29) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,42) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,43) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,44) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,45) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(0,46) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(0,47) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(1,42) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,43) = (-0.4472135954999579*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,44) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,45) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,46) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(2,42) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,43) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,44) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,45) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(2,47) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(3,42) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,43) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,44) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,45) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(3,46) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,47) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,42) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,43) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(4,45) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,46) = (-0.31943828249997*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(5,42) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,44) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(5,45) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(5,47) = (-0.31943828249997*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(18,0) = 0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(18,1) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(18,2) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(18,3) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(18,4) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(18,5) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(19,0) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(19,1) = 0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(19,2) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(19,3) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(19,4) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(20,0) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(20,1) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(20,2) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(20,3) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(20,5) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(21,0) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(21,1) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(21,2) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(21,3) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(21,4) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(21,5) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(22,0) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(22,1) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(22,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(22,4) = 0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(23,0) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(23,2) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(23,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(23,5) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(18,24) = 0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(18,25) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(18,26) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(18,27) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(18,28) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(18,29) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(19,24) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(19,25) = 0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(19,26) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(19,27) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(19,28) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(20,24) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(20,25) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(20,26) = 0.4472135954999579*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(20,27) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(20,29) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(21,24) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(21,25) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(21,26) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(21,27) = 0.4472135954999579*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(21,28) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(21,29) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(22,24) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(22,25) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(22,27) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(22,28) = 0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(23,24) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(23,26) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(23,27) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(23,29) = 0.31943828249997*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
 
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
  data->AEM_S(6,18) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,19) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(6,20) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(6,21) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(6,22) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(6,23) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(7,18) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,19) = (-0.4472135954999579*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(7,20) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(7,21) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(7,22) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(8,18) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(8,19) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(8,20) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(8,21) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(8,23) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(9,18) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,19) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(9,20) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(9,21) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.4472135954999579*cMSelf[10]*mnuSelf-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(9,22) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(9,23) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(10,18) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(10,19) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(10,21) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(10,22) = (-0.31943828249997*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(11,18) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(11,20) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(11,21) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(11,23) = (-0.31943828249997*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(6,30) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,31) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(6,32) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(6,33) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(6,34) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(6,35) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(7,30) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,31) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(7,32) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(7,33) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(7,34) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(8,30) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,31) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(8,32) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,33) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(8,35) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(9,30) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,31) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(9,32) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(9,33) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,34) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(9,35) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(10,30) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(10,31) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(10,33) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(10,34) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(11,30) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(11,32) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(11,33) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(11,35) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(6,42) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,43) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(6,44) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(6,45) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(6,46) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(6,47) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(7,42) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,43) = (-0.4472135954999579*cMOther[10]*mnuOther)-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(7,44) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(7,45) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(7,46) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(8,42) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(8,43) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(8,44) = (-0.4472135954999579*cMOther[11]*mnuOther)-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(8,45) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(8,47) = -0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(9,42) = -0.5*cMOther[9]*mnuOther; 
  data->AEM_S(9,43) = -0.5*cMOther[8]*mnuOther; 
  data->AEM_S(9,44) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(9,45) = (-0.4472135954999579*cMOther[11]*mnuOther)-0.4472135954999579*cMOther[10]*mnuOther-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(9,46) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(9,47) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(10,42) = -0.5*cMOther[10]*mnuOther; 
  data->AEM_S(10,43) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(10,45) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(10,46) = (-0.31943828249997*cMOther[10]*mnuOther)-0.5*cMOther[6]*mnuOther; 
  data->AEM_S(11,42) = -0.5*cMOther[11]*mnuOther; 
  data->AEM_S(11,44) = -0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(11,45) = -0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(11,47) = (-0.31943828249997*cMOther[11]*mnuOther)-0.5*cMOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(18,6) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(18,7) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(18,8) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(18,9) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(18,10) = 0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(18,11) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(19,6) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(19,7) = 0.4472135954999579*m1rSelf[10]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(19,8) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(19,9) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(19,10) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(20,6) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(20,7) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(20,8) = 0.4472135954999579*m1rSelf[11]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(20,9) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(20,11) = 0.4472135954999579*m1rSelf[8]*mnuSelf; 
  data->AEM_S(21,6) = 0.5*m1rSelf[9]*mnuSelf; 
  data->AEM_S(21,7) = 0.5*m1rSelf[8]*mnuSelf; 
  data->AEM_S(21,8) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(21,9) = 0.4472135954999579*m1rSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(21,10) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(21,11) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(22,6) = 0.5*m1rSelf[10]*mnuSelf; 
  data->AEM_S(22,7) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(22,9) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(22,10) = 0.31943828249997*m1rSelf[10]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(23,6) = 0.5*m1rSelf[11]*mnuSelf; 
  data->AEM_S(23,8) = 0.4472135954999579*m1rSelf[8]*mnuSelf; 
  data->AEM_S(23,9) = 0.4472135954999579*m1rSelf[9]*mnuSelf; 
  data->AEM_S(23,11) = 0.31943828249997*m1rSelf[11]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(18,30) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(18,31) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(18,32) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(18,33) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(18,34) = 0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(18,35) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(19,30) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(19,31) = 0.4472135954999579*m1rOther[10]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(19,32) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(19,33) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(19,34) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(20,30) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(20,31) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(20,32) = 0.4472135954999579*m1rOther[11]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(20,33) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(20,35) = 0.4472135954999579*m1rOther[8]*mnuOther; 
  data->AEM_S(21,30) = 0.5*m1rOther[9]*mnuOther; 
  data->AEM_S(21,31) = 0.5*m1rOther[8]*mnuOther; 
  data->AEM_S(21,32) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(21,33) = 0.4472135954999579*m1rOther[11]*mnuOther+0.4472135954999579*m1rOther[10]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(21,34) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(21,35) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(22,30) = 0.5*m1rOther[10]*mnuOther; 
  data->AEM_S(22,31) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(22,33) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(22,34) = 0.31943828249997*m1rOther[10]*mnuOther+0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(23,30) = 0.5*m1rOther[11]*mnuOther; 
  data->AEM_S(23,32) = 0.4472135954999579*m1rOther[8]*mnuOther; 
  data->AEM_S(23,33) = 0.4472135954999579*m1rOther[9]*mnuOther; 
  data->AEM_S(23,35) = 0.31943828249997*m1rOther[11]*mnuOther+0.5*m1rOther[6]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
  mnuM1sum[9] += m1rSelf[9]*mnuSelf+m1rOther[9]*mnuOther; 
  mnuM1sum[10] += m1rSelf[10]*mnuSelf+m1rOther[10]*mnuOther; 
  mnuM1sum[11] += m1rSelf[11]*mnuSelf+m1rOther[11]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(12,12) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,14) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,15) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(12,16) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(12,17) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(13,12) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(13,15) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,16) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,12) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,13) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(14,14) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(14,15) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,17) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,12) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,13) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,14) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,15) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(15,16) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(15,17) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,12) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(16,13) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,15) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,16) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(17,12) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(17,14) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,15) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,17) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(12,18) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(12,19) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(12,20) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(12,21) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(12,22) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(12,23) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(13,18) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(13,19) = (-0.4472135954999579*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(13,20) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(13,21) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(13,22) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(14,18) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(14,19) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(14,20) = (-0.4472135954999579*cMSelf[17]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(14,21) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(14,23) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(15,18) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(15,19) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(15,20) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(15,21) = (-0.4472135954999579*cMSelf[17]*mnuSelf)-0.4472135954999579*cMSelf[16]*mnuSelf-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(15,22) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(15,23) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(16,18) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(16,19) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(16,21) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(16,22) = (-0.31943828249997*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(17,18) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(17,20) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(17,21) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(17,23) = (-0.31943828249997*cMSelf[17]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(12,36) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(12,37) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(12,38) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(12,39) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(12,40) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(12,41) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(13,36) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(13,37) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(13,38) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(13,39) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(13,40) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(14,36) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(14,37) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(14,38) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(14,39) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(14,41) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(15,36) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(15,37) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(15,38) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(15,39) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(15,40) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(15,41) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(16,36) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(16,37) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(16,39) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(16,40) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(17,36) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(17,38) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(17,39) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(17,41) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(12,42) = -0.5*cMOther[12]*mnuOther; 
  data->AEM_S(12,43) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(12,44) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(12,45) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(12,46) = -0.5*cMOther[16]*mnuOther; 
  data->AEM_S(12,47) = -0.5*cMOther[17]*mnuOther; 
  data->AEM_S(13,42) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(13,43) = (-0.4472135954999579*cMOther[16]*mnuOther)-0.5*cMOther[12]*mnuOther; 
  data->AEM_S(13,44) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(13,45) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(13,46) = -0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(14,42) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(14,43) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(14,44) = (-0.4472135954999579*cMOther[17]*mnuOther)-0.5*cMOther[12]*mnuOther; 
  data->AEM_S(14,45) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(14,47) = -0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(15,42) = -0.5*cMOther[15]*mnuOther; 
  data->AEM_S(15,43) = -0.5*cMOther[14]*mnuOther; 
  data->AEM_S(15,44) = -0.5*cMOther[13]*mnuOther; 
  data->AEM_S(15,45) = (-0.4472135954999579*cMOther[17]*mnuOther)-0.4472135954999579*cMOther[16]*mnuOther-0.5*cMOther[12]*mnuOther; 
  data->AEM_S(15,46) = -0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(15,47) = -0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(16,42) = -0.5*cMOther[16]*mnuOther; 
  data->AEM_S(16,43) = -0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(16,45) = -0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(16,46) = (-0.31943828249997*cMOther[16]*mnuOther)-0.5*cMOther[12]*mnuOther; 
  data->AEM_S(17,42) = -0.5*cMOther[17]*mnuOther; 
  data->AEM_S(17,44) = -0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(17,45) = -0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(17,47) = (-0.31943828249997*cMOther[17]*mnuOther)-0.5*cMOther[12]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfZ and uCrossSelfZ ... // 
  data->AEM_S(18,12) = 0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(18,13) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(18,14) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(18,15) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(18,16) = 0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(18,17) = 0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(19,12) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(19,13) = 0.4472135954999579*m1rSelf[16]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(19,14) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(19,15) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(19,16) = 0.4472135954999579*m1rSelf[13]*mnuSelf; 
  data->AEM_S(20,12) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(20,13) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(20,14) = 0.4472135954999579*m1rSelf[17]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(20,15) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(20,17) = 0.4472135954999579*m1rSelf[14]*mnuSelf; 
  data->AEM_S(21,12) = 0.5*m1rSelf[15]*mnuSelf; 
  data->AEM_S(21,13) = 0.5*m1rSelf[14]*mnuSelf; 
  data->AEM_S(21,14) = 0.5*m1rSelf[13]*mnuSelf; 
  data->AEM_S(21,15) = 0.4472135954999579*m1rSelf[17]*mnuSelf+0.4472135954999579*m1rSelf[16]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(21,16) = 0.4472135954999579*m1rSelf[15]*mnuSelf; 
  data->AEM_S(21,17) = 0.4472135954999579*m1rSelf[15]*mnuSelf; 
  data->AEM_S(22,12) = 0.5*m1rSelf[16]*mnuSelf; 
  data->AEM_S(22,13) = 0.4472135954999579*m1rSelf[13]*mnuSelf; 
  data->AEM_S(22,15) = 0.4472135954999579*m1rSelf[15]*mnuSelf; 
  data->AEM_S(22,16) = 0.31943828249997*m1rSelf[16]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(23,12) = 0.5*m1rSelf[17]*mnuSelf; 
  data->AEM_S(23,14) = 0.4472135954999579*m1rSelf[14]*mnuSelf; 
  data->AEM_S(23,15) = 0.4472135954999579*m1rSelf[15]*mnuSelf; 
  data->AEM_S(23,17) = 0.31943828249997*m1rSelf[17]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherZ and uCrossOtherZ ... // 
  data->AEM_S(18,36) = 0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(18,37) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(18,38) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(18,39) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(18,40) = 0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(18,41) = 0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(19,36) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(19,37) = 0.4472135954999579*m1rOther[16]*mnuOther+0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(19,38) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(19,39) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(19,40) = 0.4472135954999579*m1rOther[13]*mnuOther; 
  data->AEM_S(20,36) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(20,37) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(20,38) = 0.4472135954999579*m1rOther[17]*mnuOther+0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(20,39) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(20,41) = 0.4472135954999579*m1rOther[14]*mnuOther; 
  data->AEM_S(21,36) = 0.5*m1rOther[15]*mnuOther; 
  data->AEM_S(21,37) = 0.5*m1rOther[14]*mnuOther; 
  data->AEM_S(21,38) = 0.5*m1rOther[13]*mnuOther; 
  data->AEM_S(21,39) = 0.4472135954999579*m1rOther[17]*mnuOther+0.4472135954999579*m1rOther[16]*mnuOther+0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(21,40) = 0.4472135954999579*m1rOther[15]*mnuOther; 
  data->AEM_S(21,41) = 0.4472135954999579*m1rOther[15]*mnuOther; 
  data->AEM_S(22,36) = 0.5*m1rOther[16]*mnuOther; 
  data->AEM_S(22,37) = 0.4472135954999579*m1rOther[13]*mnuOther; 
  data->AEM_S(22,39) = 0.4472135954999579*m1rOther[15]*mnuOther; 
  data->AEM_S(22,40) = 0.31943828249997*m1rOther[16]*mnuOther+0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(23,36) = 0.5*m1rOther[17]*mnuOther; 
  data->AEM_S(23,38) = 0.4472135954999579*m1rOther[14]*mnuOther; 
  data->AEM_S(23,39) = 0.4472135954999579*m1rOther[15]*mnuOther; 
  data->AEM_S(23,41) = 0.31943828249997*m1rOther[17]*mnuOther+0.5*m1rOther[12]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[12] += m1rSelf[12]*mnuSelf+m1rOther[12]*mnuOther; 
  mnuM1sum[13] += m1rSelf[13]*mnuSelf+m1rOther[13]*mnuOther; 
  mnuM1sum[14] += m1rSelf[14]*mnuSelf+m1rOther[14]*mnuOther; 
  mnuM1sum[15] += m1rSelf[15]*mnuSelf+m1rOther[15]*mnuOther; 
  mnuM1sum[16] += m1rSelf[16]*mnuSelf+m1rOther[16]*mnuOther; 
  mnuM1sum[17] += m1rSelf[17]*mnuSelf+m1rOther[17]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(18,18) = 1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(18,19) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(18,20) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(18,21) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(18,22) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(18,23) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(19,18) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(19,19) = 1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(19,20) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(19,21) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(19,22) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(20,18) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(20,19) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(20,20) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(20,21) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(20,23) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(21,18) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(21,19) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(21,20) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(21,21) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(21,22) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(21,23) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(22,18) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(22,19) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(22,21) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(22,22) = 0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(23,18) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(23,20) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(23,21) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(23,23) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(18,42) = 1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(18,43) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(18,44) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(18,45) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(18,46) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(18,47) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(19,42) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(19,43) = 1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(19,44) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(19,45) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(19,46) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(20,42) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(20,43) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(20,44) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(20,45) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(20,47) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(21,42) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(21,43) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(21,44) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(21,45) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(21,46) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(21,47) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(22,42) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(22,43) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(22,45) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(22,46) = 0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(23,42) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(23,44) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(23,45) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(23,47) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[6]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2rSelf[0]*mnuSelf+m2rOther[0]*mnuOther; 
  mnuM2sum[1] = m2rSelf[1]*mnuSelf+m2rOther[1]*mnuOther; 
  mnuM2sum[2] = m2rSelf[2]*mnuSelf+m2rOther[2]*mnuOther; 
  mnuM2sum[3] = m2rSelf[3]*mnuSelf+m2rOther[3]*mnuOther; 
  mnuM2sum[4] = m2rSelf[4]*mnuSelf+m2rOther[4]*mnuOther; 
  mnuM2sum[5] = m2rSelf[5]*mnuSelf+m2rOther[5]*mnuOther; 
 
  // Set other entries to 0. // 
  data->AEM_S.block<6,12>(0,6).setZero(); 
  data->AEM_S.block<12,6>(6,0).setZero(); 
  data->AEM_S.block<6,6>(6,12).setZero(); 
  data->AEM_S.block<6,6>(12,6).setZero(); 
  data->AEM_S.block<6,12>(0,30).setZero(); 
  data->AEM_S.block<12,6>(6,24).setZero(); 
  data->AEM_S.block<6,6>(6,36).setZero(); 
  data->AEM_S.block<6,6>(12,30).setZero(); 
 
  double m1Relax[18]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<18; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(24,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(24,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(24,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(24,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(24,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(25,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(25,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(25,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(26,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(26,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(26,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(26,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,1) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,2) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(27,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(27,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(28,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(28,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(28,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(28,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(29,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(29,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(24,18) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(24,19) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(24,20) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(24,21) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(24,22) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(24,23) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(25,18) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(25,19) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(25,20) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(25,21) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(25,22) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(26,18) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(26,19) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(26,20) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(26,21) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(26,23) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(27,18) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(27,19) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(27,20) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(27,21) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(27,22) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(27,23) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(28,18) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(28,19) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(28,21) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(28,22) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(29,18) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(29,20) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(29,21) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(29,23) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
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
  data->AEM_S(24,42) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(24,43) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(24,44) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(24,45) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(24,46) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(24,47) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(25,42) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(25,43) = 0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(25,44) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(25,45) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(25,46) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(26,42) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(26,43) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(26,44) = 0.4472135954999579*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(26,45) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(26,47) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(27,42) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(27,43) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(27,44) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(27,45) = 0.4472135954999579*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(27,46) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(27,47) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(28,42) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(28,43) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(28,45) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(28,46) = 0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(29,42) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(29,44) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(29,45) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(29,47) = 0.31943828249997*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(42,0) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(42,1) = (-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(42,2) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(42,3) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(42,4) = (-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf)-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(42,5) = (-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(43,0) = (-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(43,1) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(43,2) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(43,3) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(43,4) = (-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf)-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(43,5) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(44,0) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(44,1) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(44,2) = (-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf)-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(44,3) = (-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(44,4) = (-0.25*m0rSelf[2]*uSelf[4]*mnuSelf)-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(44,5) = (-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(45,0) = (-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(45,1) = (-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(45,2) = (-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf)-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(45,3) = (-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf)-0.2*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.2*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.81*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(45,4) = (-0.2*m0rSelf[3]*uSelf[5]*mnuSelf)-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(45,5) = (-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(46,0) = (-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf)-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(46,1) = (-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf)-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(46,2) = (-0.25*m0rSelf[2]*uSelf[4]*mnuSelf)-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(46,3) = (-0.2*m0rSelf[3]*uSelf[5]*mnuSelf)-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(46,4) = (-0.25*m0rSelf[5]*uSelf[5]*mnuSelf)-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(46,5) = (-0.25*m0rSelf[4]*uSelf[5]*mnuSelf)-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(47,0) = (-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf)-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(47,1) = (-0.25*m0rSelf[1]*uSelf[5]*mnuSelf)-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(47,2) = (-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(47,3) = (-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf)-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(47,4) = (-0.25*m0rSelf[4]*uSelf[5]*mnuSelf)-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(47,5) = (-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf)-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(42,24) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(42,25) = 0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(42,26) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(42,27) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(42,28) = 0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(42,29) = 0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(43,24) = 0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(43,25) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(43,26) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(43,27) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(43,28) = 0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(43,29) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(44,24) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(44,25) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(44,26) = 0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(44,27) = 0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(44,28) = 0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(44,29) = 0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(45,24) = 0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(45,25) = 0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(45,26) = 0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(45,27) = 0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.2*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.2*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(45,28) = 0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(45,29) = 0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(46,24) = 0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(46,25) = 0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(46,26) = 0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(46,27) = 0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(46,28) = 0.25*m0rOther[5]*uOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(46,29) = 0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(47,24) = 0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(47,25) = 0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(47,26) = 0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(47,27) = 0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(47,28) = 0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(47,29) = 0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (m1rOther[0]-1.0*m1rSelf[0])*betaGreenep1*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (m1rOther[1]-1.0*m1rSelf[1])*betaGreenep1*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (m1rOther[2]-1.0*m1rSelf[2])*betaGreenep1*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += (m1rOther[3]-1.0*m1rSelf[3])*betaGreenep1*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (m1rOther[4]-1.0*m1rSelf[4])*betaGreenep1*mnuSelf+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (m1rOther[5]-1.0*m1rSelf[5])*betaGreenep1*mnuSelf+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(30,6) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(30,7) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(30,8) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,9) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(30,10) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(30,11) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(31,6) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,7) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(31,8) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,9) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,10) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(32,6) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(32,7) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(32,8) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(32,9) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(32,11) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,6) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,7) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(33,8) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(33,9) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(33,10) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(33,11) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,6) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(34,7) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(34,9) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(34,10) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(35,6) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(35,8) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(35,9) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(35,11) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(30,18) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(30,19) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(30,20) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(30,21) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(30,22) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(30,23) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(31,18) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(31,19) = (-0.4472135954999579*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(31,20) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(31,21) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(31,22) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(32,18) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(32,19) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(32,20) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(32,21) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(32,23) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(33,18) = -0.5*cMSelf[9]*mnuSelf; 
  data->AEM_S(33,19) = -0.5*cMSelf[8]*mnuSelf; 
  data->AEM_S(33,20) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(33,21) = (-0.4472135954999579*cMSelf[11]*mnuSelf)-0.4472135954999579*cMSelf[10]*mnuSelf-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(33,22) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(33,23) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(34,18) = -0.5*cMSelf[10]*mnuSelf; 
  data->AEM_S(34,19) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(34,21) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(34,22) = (-0.31943828249997*cMSelf[10]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(35,18) = -0.5*cMSelf[11]*mnuSelf; 
  data->AEM_S(35,20) = -0.4472135954999579*cMSelf[8]*mnuSelf; 
  data->AEM_S(35,21) = -0.4472135954999579*cMSelf[9]*mnuSelf; 
  data->AEM_S(35,23) = (-0.31943828249997*cMSelf[11]*mnuSelf)-0.5*cMSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(30,30) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(30,31) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(30,32) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(30,33) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(30,34) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(30,35) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(31,30) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(31,31) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(31,32) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(31,33) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(31,34) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(32,30) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(32,31) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(32,32) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(32,33) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(32,35) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(33,30) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(33,31) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(33,32) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(33,33) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(33,34) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(33,35) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(34,30) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(34,31) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(34,33) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(34,34) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(35,30) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(35,32) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(35,33) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(35,35) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(30,42) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(30,43) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(30,44) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(30,45) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(30,46) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(30,47) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(31,42) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(31,43) = 0.4472135954999579*cMOther[10]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(31,44) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(31,45) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(31,46) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(32,42) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(32,43) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(32,44) = 0.4472135954999579*cMOther[11]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(32,45) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(32,47) = 0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(33,42) = 0.5*cMOther[9]*mnuOther; 
  data->AEM_S(33,43) = 0.5*cMOther[8]*mnuOther; 
  data->AEM_S(33,44) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(33,45) = 0.4472135954999579*cMOther[11]*mnuOther+0.4472135954999579*cMOther[10]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(33,46) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(33,47) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(34,42) = 0.5*cMOther[10]*mnuOther; 
  data->AEM_S(34,43) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(34,45) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(34,46) = 0.31943828249997*cMOther[10]*mnuOther+0.5*cMOther[6]*mnuOther; 
  data->AEM_S(35,42) = 0.5*cMOther[11]*mnuOther; 
  data->AEM_S(35,44) = 0.4472135954999579*cMOther[8]*mnuOther; 
  data->AEM_S(35,45) = 0.4472135954999579*cMOther[9]*mnuOther; 
  data->AEM_S(35,47) = 0.31943828249997*cMOther[11]*mnuOther+0.5*cMOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(42,6) = (-0.25*m0rSelf[5]*uSelf[11]*mnuSelf)-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.25*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(42,7) = (-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf)-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(42,8) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.25*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(42,9) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(42,10) = (-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[4]*uSelf[6]*mnuSelf; 
  data->AEM_S(42,11) = (-0.159719141249985*m0rSelf[5]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[6]*mnuSelf; 
  data->AEM_S(43,6) = (-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf)-0.25*m0rSelf[2]*uSelf[9]*mnuSelf-0.25*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(43,7) = (-0.25*m0rSelf[5]*uSelf[11]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(43,8) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(43,9) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.45*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(43,10) = (-0.3928571428571428*m0rSelf[1]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(43,11) = (-0.25*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[7]*mnuSelf; 
  data->AEM_S(44,6) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.25*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.25*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(44,7) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(44,8) = (-0.3928571428571428*m0rSelf[5]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.45*m0rSelf[3]*uSelf[9]*mnuSelf-0.45*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(44,9) = (-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.45*m0rSelf[2]*uSelf[9]*mnuSelf-0.45*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(44,10) = (-0.25*m0rSelf[2]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[4]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf; 
  data->AEM_S(44,11) = (-0.3928571428571428*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(45,6) = (-0.223606797749979*m0rSelf[3]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[9]*mnuSelf-0.25*m0rSelf[0]*uSelf[9]*mnuSelf+0.5*m1rSelf[9]*mnuSelf-0.25*m0rSelf[1]*uSelf[8]*mnuSelf-0.25*m0rSelf[2]*uSelf[7]*mnuSelf-0.25*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(45,7) = (-0.223606797749979*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[10]*mnuSelf-0.45*m0rSelf[1]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[8]*mnuSelf-0.25*m0rSelf[0]*uSelf[8]*mnuSelf+0.5*m1rSelf[8]*mnuSelf-0.45*m0rSelf[3]*uSelf[7]*mnuSelf-0.25*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(45,8) = (-0.223606797749979*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[10]*mnuSelf-0.45*m0rSelf[2]*uSelf[9]*mnuSelf-0.45*m0rSelf[3]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.25*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(45,9) = (-0.3928571428571428*m0rSelf[5]*uSelf[11]*mnuSelf)-0.2*m0rSelf[4]*uSelf[11]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[11]*mnuSelf+0.4472135954999579*m1rSelf[11]*mnuSelf-0.2*m0rSelf[5]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[10]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[10]*mnuSelf+0.4472135954999579*m1rSelf[10]*mnuSelf-0.81*m0rSelf[3]*uSelf[9]*mnuSelf-0.45*m0rSelf[2]*uSelf[8]*mnuSelf-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(45,10) = (-0.2*m0rSelf[3]*uSelf[11]*mnuSelf)-0.3928571428571428*m0rSelf[3]*uSelf[10]*mnuSelf-0.2*m0rSelf[5]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(45,11) = (-0.3928571428571428*m0rSelf[3]*uSelf[11]*mnuSelf)-0.2*m0rSelf[3]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[9]*mnuSelf-0.2*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(46,6) = (-0.159719141249985*m0rSelf[4]*uSelf[10]*mnuSelf)-0.25*m0rSelf[0]*uSelf[10]*mnuSelf+0.5*m1rSelf[10]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.25*m0rSelf[4]*uSelf[6]*mnuSelf; 
  data->AEM_S(46,7) = (-0.3928571428571428*m0rSelf[1]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(46,8) = (-0.25*m0rSelf[2]*uSelf[10]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.25*m0rSelf[4]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf; 
  data->AEM_S(46,9) = (-0.2*m0rSelf[3]*uSelf[11]*mnuSelf)-0.3928571428571428*m0rSelf[3]*uSelf[10]*mnuSelf-0.2*m0rSelf[5]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(46,10) = (-0.25*m0rSelf[5]*uSelf[11]*mnuSelf)-0.5357142857142857*m0rSelf[4]*uSelf[10]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[10]*mnuSelf+0.31943828249997*m1rSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[9]*mnuSelf-0.25*m0rSelf[2]*uSelf[8]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(46,11) = (-0.25*m0rSelf[4]*uSelf[11]*mnuSelf)-0.25*m0rSelf[5]*uSelf[10]*mnuSelf-0.2*m0rSelf[3]*uSelf[9]*mnuSelf; 
  data->AEM_S(47,6) = (-0.159719141249985*m0rSelf[5]*uSelf[11]*mnuSelf)-0.25*m0rSelf[0]*uSelf[11]*mnuSelf+0.5*m1rSelf[11]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[6]*mnuSelf; 
  data->AEM_S(47,7) = (-0.25*m0rSelf[1]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[8]*mnuSelf-0.25*m0rSelf[5]*uSelf[7]*mnuSelf; 
  data->AEM_S(47,8) = (-0.3928571428571428*m0rSelf[2]*uSelf[11]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[8]*mnuSelf+0.4472135954999579*m1rSelf[8]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(47,9) = (-0.3928571428571428*m0rSelf[3]*uSelf[11]*mnuSelf)-0.2*m0rSelf[3]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[9]*mnuSelf-0.2*m0rSelf[4]*uSelf[9]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[9]*mnuSelf+0.4472135954999579*m1rSelf[9]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[8]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf; 
  data->AEM_S(47,10) = (-0.25*m0rSelf[4]*uSelf[11]*mnuSelf)-0.25*m0rSelf[5]*uSelf[10]*mnuSelf-0.2*m0rSelf[3]*uSelf[9]*mnuSelf; 
  data->AEM_S(47,11) = (-0.5357142857142857*m0rSelf[5]*uSelf[11]*mnuSelf)-0.159719141249985*m0rSelf[0]*uSelf[11]*mnuSelf+0.31943828249997*m1rSelf[11]*mnuSelf-0.25*m0rSelf[4]*uSelf[10]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[9]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[8]*mnuSelf-0.25*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(42,30) = 0.25*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(42,31) = 0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(42,32) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.25*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(42,33) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(42,34) = 0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[4]*uOther[6]*mnuOther; 
  data->AEM_S(42,35) = 0.159719141249985*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[6]*mnuOther; 
  data->AEM_S(43,30) = 0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.25*m0rOther[2]*uOther[9]*mnuOther+0.25*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(43,31) = 0.25*m0rOther[5]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.4472135954999579*m1rOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(43,32) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(43,33) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.45*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(43,34) = 0.3928571428571428*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(43,35) = 0.25*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[7]*mnuOther; 
  data->AEM_S(44,30) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.25*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.25*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(44,31) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(44,32) = 0.3928571428571428*m0rOther[5]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.45*m0rOther[3]*uOther[9]*mnuOther+0.45*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(44,33) = 0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.45*m0rOther[2]*uOther[9]*mnuOther+0.45*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(44,34) = 0.25*m0rOther[2]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[4]*uOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther; 
  data->AEM_S(44,35) = 0.3928571428571428*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[0]*uOther[8]*mnuOther-0.4472135954999579*m1rOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(45,30) = 0.223606797749979*m0rOther[3]*uOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[10]*mnuOther+0.223606797749979*m0rOther[5]*uOther[9]*mnuOther+0.223606797749979*m0rOther[4]*uOther[9]*mnuOther+0.25*m0rOther[0]*uOther[9]*mnuOther-0.5*m1rOther[9]*mnuOther+0.25*m0rOther[1]*uOther[8]*mnuOther+0.25*m0rOther[2]*uOther[7]*mnuOther+0.25*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(45,31) = 0.223606797749979*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[10]*mnuOther+0.45*m0rOther[1]*uOther[9]*mnuOther+0.223606797749979*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[4]*uOther[8]*mnuOther+0.25*m0rOther[0]*uOther[8]*mnuOther-0.5*m1rOther[8]*mnuOther+0.45*m0rOther[3]*uOther[7]*mnuOther+0.25*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(45,32) = 0.223606797749979*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[10]*mnuOther+0.45*m0rOther[2]*uOther[9]*mnuOther+0.45*m0rOther[3]*uOther[8]*mnuOther+0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.25*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(45,33) = 0.3928571428571428*m0rOther[5]*uOther[11]*mnuOther+0.2*m0rOther[4]*uOther[11]*mnuOther+0.223606797749979*m0rOther[0]*uOther[11]*mnuOther-0.4472135954999579*m1rOther[11]*mnuOther+0.2*m0rOther[5]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[10]*mnuOther+0.223606797749979*m0rOther[0]*uOther[10]*mnuOther-0.4472135954999579*m1rOther[10]*mnuOther+0.81*m0rOther[3]*uOther[9]*mnuOther+0.45*m0rOther[2]*uOther[8]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(45,34) = 0.2*m0rOther[3]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[10]*mnuOther+0.2*m0rOther[5]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(45,35) = 0.3928571428571428*m0rOther[3]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[9]*mnuOther+0.2*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(46,30) = 0.159719141249985*m0rOther[4]*uOther[10]*mnuOther+0.25*m0rOther[0]*uOther[10]*mnuOther-0.5*m1rOther[10]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.25*m0rOther[4]*uOther[6]*mnuOther; 
  data->AEM_S(46,31) = 0.3928571428571428*m0rOther[1]*uOther[10]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(46,32) = 0.25*m0rOther[2]*uOther[10]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.25*m0rOther[4]*uOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther; 
  data->AEM_S(46,33) = 0.2*m0rOther[3]*uOther[11]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[10]*mnuOther+0.2*m0rOther[5]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(46,34) = 0.25*m0rOther[5]*uOther[11]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[10]*mnuOther+0.159719141249985*m0rOther[0]*uOther[10]*mnuOther-0.31943828249997*m1rOther[10]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[9]*mnuOther+0.25*m0rOther[2]*uOther[8]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(46,35) = 0.25*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[5]*uOther[10]*mnuOther+0.2*m0rOther[3]*uOther[9]*mnuOther; 
  data->AEM_S(47,30) = 0.159719141249985*m0rOther[5]*uOther[11]*mnuOther+0.25*m0rOther[0]*uOther[11]*mnuOther-0.5*m1rOther[11]*mnuOther+0.223606797749979*m0rOther[3]*uOther[9]*mnuOther+0.223606797749979*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[6]*mnuOther; 
  data->AEM_S(47,31) = 0.25*m0rOther[1]*uOther[11]*mnuOther+0.223606797749979*m0rOther[2]*uOther[9]*mnuOther+0.223606797749979*m0rOther[3]*uOther[8]*mnuOther+0.25*m0rOther[5]*uOther[7]*mnuOther; 
  data->AEM_S(47,32) = 0.3928571428571428*m0rOther[2]*uOther[11]*mnuOther+0.223606797749979*m0rOther[1]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[8]*mnuOther+0.223606797749979*m0rOther[0]*uOther[8]*mnuOther-0.4472135954999579*m1rOther[8]*mnuOther+0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(47,33) = 0.3928571428571428*m0rOther[3]*uOther[11]*mnuOther+0.2*m0rOther[3]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[9]*mnuOther+0.2*m0rOther[4]*uOther[9]*mnuOther+0.223606797749979*m0rOther[0]*uOther[9]*mnuOther-0.4472135954999579*m1rOther[9]*mnuOther+0.223606797749979*m0rOther[1]*uOther[8]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther; 
  data->AEM_S(47,34) = 0.25*m0rOther[4]*uOther[11]*mnuOther+0.25*m0rOther[5]*uOther[10]*mnuOther+0.2*m0rOther[3]*uOther[9]*mnuOther; 
  data->AEM_S(47,35) = 0.5357142857142857*m0rOther[5]*uOther[11]*mnuOther+0.159719141249985*m0rOther[0]*uOther[11]*mnuOther-0.31943828249997*m1rOther[11]*mnuOther+0.25*m0rOther[4]*uOther[10]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[9]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[8]*mnuOther+0.25*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*m0rOther[5]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[6] += (m1rOther[6]-1.0*m1rSelf[6])*betaGreenep1*mnuSelf+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (m1rOther[7]-1.0*m1rSelf[7])*betaGreenep1*mnuSelf+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
  m1Relax[8] += (m1rOther[8]-1.0*m1rSelf[8])*betaGreenep1*mnuSelf+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
  m1Relax[9] += (m1rOther[9]-1.0*m1rSelf[9])*betaGreenep1*mnuSelf+m1rSelf[9]*mnuSelf-1.0*m1rOther[9]*mnuOther; 
  m1Relax[10] += (m1rOther[10]-1.0*m1rSelf[10])*betaGreenep1*mnuSelf+m1rSelf[10]*mnuSelf-1.0*m1rOther[10]*mnuOther; 
  m1Relax[11] += (m1rOther[11]-1.0*m1rSelf[11])*betaGreenep1*mnuSelf+m1rSelf[11]*mnuSelf-1.0*m1rOther[11]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(36,12) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(36,13) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(36,14) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(36,15) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(36,16) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(36,17) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(37,12) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(37,13) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(37,14) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(37,15) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(37,16) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(38,12) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(38,13) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(38,14) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(38,15) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(38,17) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(39,12) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(39,13) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(39,14) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(39,15) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(39,16) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(39,17) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(40,12) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(40,13) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(40,15) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(40,16) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(41,12) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(41,14) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(41,15) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(41,17) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(36,18) = -0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(36,19) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(36,20) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(36,21) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(36,22) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(36,23) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(37,18) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(37,19) = (-0.4472135954999579*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(37,20) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(37,21) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(37,22) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(38,18) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(38,19) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(38,20) = (-0.4472135954999579*cMSelf[17]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(38,21) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(38,23) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(39,18) = -0.5*cMSelf[15]*mnuSelf; 
  data->AEM_S(39,19) = -0.5*cMSelf[14]*mnuSelf; 
  data->AEM_S(39,20) = -0.5*cMSelf[13]*mnuSelf; 
  data->AEM_S(39,21) = (-0.4472135954999579*cMSelf[17]*mnuSelf)-0.4472135954999579*cMSelf[16]*mnuSelf-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(39,22) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(39,23) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(40,18) = -0.5*cMSelf[16]*mnuSelf; 
  data->AEM_S(40,19) = -0.4472135954999579*cMSelf[13]*mnuSelf; 
  data->AEM_S(40,21) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(40,22) = (-0.31943828249997*cMSelf[16]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
  data->AEM_S(41,18) = -0.5*cMSelf[17]*mnuSelf; 
  data->AEM_S(41,20) = -0.4472135954999579*cMSelf[14]*mnuSelf; 
  data->AEM_S(41,21) = -0.4472135954999579*cMSelf[15]*mnuSelf; 
  data->AEM_S(41,23) = (-0.31943828249997*cMSelf[17]*mnuSelf)-0.5*cMSelf[12]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(36,36) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(36,37) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(36,38) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(36,39) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(36,40) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(36,41) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(37,36) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(37,37) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(37,38) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(37,39) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(37,40) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(38,36) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(38,37) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(38,38) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(38,39) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(38,41) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(39,36) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(39,37) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(39,38) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(39,39) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(39,40) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(39,41) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(40,36) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(40,37) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(40,39) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(40,40) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(41,36) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(41,38) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(41,39) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(41,41) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(36,42) = 0.5*cMOther[12]*mnuOther; 
  data->AEM_S(36,43) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(36,44) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(36,45) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(36,46) = 0.5*cMOther[16]*mnuOther; 
  data->AEM_S(36,47) = 0.5*cMOther[17]*mnuOther; 
  data->AEM_S(37,42) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(37,43) = 0.4472135954999579*cMOther[16]*mnuOther+0.5*cMOther[12]*mnuOther; 
  data->AEM_S(37,44) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(37,45) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(37,46) = 0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(38,42) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(38,43) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(38,44) = 0.4472135954999579*cMOther[17]*mnuOther+0.5*cMOther[12]*mnuOther; 
  data->AEM_S(38,45) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(38,47) = 0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(39,42) = 0.5*cMOther[15]*mnuOther; 
  data->AEM_S(39,43) = 0.5*cMOther[14]*mnuOther; 
  data->AEM_S(39,44) = 0.5*cMOther[13]*mnuOther; 
  data->AEM_S(39,45) = 0.4472135954999579*cMOther[17]*mnuOther+0.4472135954999579*cMOther[16]*mnuOther+0.5*cMOther[12]*mnuOther; 
  data->AEM_S(39,46) = 0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(39,47) = 0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(40,42) = 0.5*cMOther[16]*mnuOther; 
  data->AEM_S(40,43) = 0.4472135954999579*cMOther[13]*mnuOther; 
  data->AEM_S(40,45) = 0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(40,46) = 0.31943828249997*cMOther[16]*mnuOther+0.5*cMOther[12]*mnuOther; 
  data->AEM_S(41,42) = 0.5*cMOther[17]*mnuOther; 
  data->AEM_S(41,44) = 0.4472135954999579*cMOther[14]*mnuOther; 
  data->AEM_S(41,45) = 0.4472135954999579*cMOther[15]*mnuOther; 
  data->AEM_S(41,47) = 0.31943828249997*cMOther[17]*mnuOther+0.5*cMOther[12]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfZ-uSelfZ*m0Self) and uCrossSelfZ ... // 
  data->AEM_S(42,12) = (-0.25*m0rSelf[5]*uSelf[17]*mnuSelf)-0.25*m0rSelf[4]*uSelf[16]*mnuSelf-0.25*m0rSelf[3]*uSelf[15]*mnuSelf-0.25*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[1]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(42,13) = (-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf)-0.25*m0rSelf[2]*uSelf[15]*mnuSelf-0.25*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf; 
  data->AEM_S(42,14) = (-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf)-0.25*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.25*m0rSelf[3]*uSelf[13]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf; 
  data->AEM_S(42,15) = (-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.25*m0rSelf[1]*uSelf[14]*mnuSelf-0.25*m0rSelf[2]*uSelf[13]*mnuSelf-0.25*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(42,16) = (-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf)-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.25*m0rSelf[4]*uSelf[12]*mnuSelf; 
  data->AEM_S(42,17) = (-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf)-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[5]*uSelf[12]*mnuSelf; 
  data->AEM_S(43,12) = (-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf)-0.25*m0rSelf[2]*uSelf[15]*mnuSelf-0.25*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf; 
  data->AEM_S(43,13) = (-0.25*m0rSelf[5]*uSelf[17]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[16]*mnuSelf+0.4472135954999579*m1rSelf[16]*mnuSelf-0.45*m0rSelf[3]*uSelf[15]*mnuSelf-0.25*m0rSelf[2]*uSelf[14]*mnuSelf-0.45*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(43,14) = (-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.25*m0rSelf[1]*uSelf[14]*mnuSelf-0.25*m0rSelf[2]*uSelf[13]*mnuSelf-0.25*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(43,15) = (-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf-0.45*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.45*m0rSelf[3]*uSelf[13]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf; 
  data->AEM_S(43,16) = (-0.3928571428571428*m0rSelf[1]*uSelf[16]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf; 
  data->AEM_S(43,17) = (-0.25*m0rSelf[1]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.25*m0rSelf[5]*uSelf[13]*mnuSelf; 
  data->AEM_S(44,12) = (-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf)-0.25*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.25*m0rSelf[3]*uSelf[13]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf; 
  data->AEM_S(44,13) = (-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.25*m0rSelf[1]*uSelf[14]*mnuSelf-0.25*m0rSelf[2]*uSelf[13]*mnuSelf-0.25*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(44,14) = (-0.3928571428571428*m0rSelf[5]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.4472135954999579*m1rSelf[17]*mnuSelf-0.25*m0rSelf[4]*uSelf[16]*mnuSelf-0.45*m0rSelf[3]*uSelf[15]*mnuSelf-0.45*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(44,15) = (-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf-0.45*m0rSelf[2]*uSelf[15]*mnuSelf-0.45*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf; 
  data->AEM_S(44,16) = (-0.25*m0rSelf[2]*uSelf[16]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.25*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf; 
  data->AEM_S(44,17) = (-0.3928571428571428*m0rSelf[2]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf; 
  data->AEM_S(45,12) = (-0.223606797749979*m0rSelf[3]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[3]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[15]*mnuSelf-0.25*m0rSelf[0]*uSelf[15]*mnuSelf+0.5*m1rSelf[15]*mnuSelf-0.25*m0rSelf[1]*uSelf[14]*mnuSelf-0.25*m0rSelf[2]*uSelf[13]*mnuSelf-0.25*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(45,13) = (-0.223606797749979*m0rSelf[2]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[16]*mnuSelf-0.45*m0rSelf[1]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[14]*mnuSelf-0.25*m0rSelf[0]*uSelf[14]*mnuSelf+0.5*m1rSelf[14]*mnuSelf-0.45*m0rSelf[3]*uSelf[13]*mnuSelf-0.25*m0rSelf[2]*uSelf[12]*mnuSelf; 
  data->AEM_S(45,14) = (-0.223606797749979*m0rSelf[1]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[16]*mnuSelf-0.45*m0rSelf[2]*uSelf[15]*mnuSelf-0.45*m0rSelf[3]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[13]*mnuSelf-0.25*m0rSelf[0]*uSelf[13]*mnuSelf+0.5*m1rSelf[13]*mnuSelf-0.25*m0rSelf[1]*uSelf[12]*mnuSelf; 
  data->AEM_S(45,15) = (-0.3928571428571428*m0rSelf[5]*uSelf[17]*mnuSelf)-0.2*m0rSelf[4]*uSelf[17]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[17]*mnuSelf+0.4472135954999579*m1rSelf[17]*mnuSelf-0.2*m0rSelf[5]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[16]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[16]*mnuSelf+0.4472135954999579*m1rSelf[16]*mnuSelf-0.81*m0rSelf[3]*uSelf[15]*mnuSelf-0.45*m0rSelf[2]*uSelf[14]*mnuSelf-0.45*m0rSelf[1]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[12]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(45,16) = (-0.2*m0rSelf[3]*uSelf[17]*mnuSelf)-0.3928571428571428*m0rSelf[3]*uSelf[16]*mnuSelf-0.2*m0rSelf[5]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(45,17) = (-0.3928571428571428*m0rSelf[3]*uSelf[17]*mnuSelf)-0.2*m0rSelf[3]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[15]*mnuSelf-0.2*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(46,12) = (-0.159719141249985*m0rSelf[4]*uSelf[16]*mnuSelf)-0.25*m0rSelf[0]*uSelf[16]*mnuSelf+0.5*m1rSelf[16]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[13]*mnuSelf-0.25*m0rSelf[4]*uSelf[12]*mnuSelf; 
  data->AEM_S(46,13) = (-0.3928571428571428*m0rSelf[1]*uSelf[16]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[13]*mnuSelf+0.4472135954999579*m1rSelf[13]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[12]*mnuSelf; 
  data->AEM_S(46,14) = (-0.25*m0rSelf[2]*uSelf[16]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.25*m0rSelf[4]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf; 
  data->AEM_S(46,15) = (-0.2*m0rSelf[3]*uSelf[17]*mnuSelf)-0.3928571428571428*m0rSelf[3]*uSelf[16]*mnuSelf-0.2*m0rSelf[5]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(46,16) = (-0.25*m0rSelf[5]*uSelf[17]*mnuSelf)-0.5357142857142857*m0rSelf[4]*uSelf[16]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[16]*mnuSelf+0.31943828249997*m1rSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[15]*mnuSelf-0.25*m0rSelf[2]*uSelf[14]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
  data->AEM_S(46,17) = (-0.25*m0rSelf[4]*uSelf[17]*mnuSelf)-0.25*m0rSelf[5]*uSelf[16]*mnuSelf-0.2*m0rSelf[3]*uSelf[15]*mnuSelf; 
  data->AEM_S(47,12) = (-0.159719141249985*m0rSelf[5]*uSelf[17]*mnuSelf)-0.25*m0rSelf[0]*uSelf[17]*mnuSelf+0.5*m1rSelf[17]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[5]*uSelf[12]*mnuSelf; 
  data->AEM_S(47,13) = (-0.25*m0rSelf[1]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[14]*mnuSelf-0.25*m0rSelf[5]*uSelf[13]*mnuSelf; 
  data->AEM_S(47,14) = (-0.3928571428571428*m0rSelf[2]*uSelf[17]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[14]*mnuSelf+0.4472135954999579*m1rSelf[14]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[12]*mnuSelf; 
  data->AEM_S(47,15) = (-0.3928571428571428*m0rSelf[3]*uSelf[17]*mnuSelf)-0.2*m0rSelf[3]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[15]*mnuSelf-0.2*m0rSelf[4]*uSelf[15]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[15]*mnuSelf+0.4472135954999579*m1rSelf[15]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[14]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[13]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[12]*mnuSelf; 
  data->AEM_S(47,16) = (-0.25*m0rSelf[4]*uSelf[17]*mnuSelf)-0.25*m0rSelf[5]*uSelf[16]*mnuSelf-0.2*m0rSelf[3]*uSelf[15]*mnuSelf; 
  data->AEM_S(47,17) = (-0.5357142857142857*m0rSelf[5]*uSelf[17]*mnuSelf)-0.159719141249985*m0rSelf[0]*uSelf[17]*mnuSelf+0.31943828249997*m1rSelf[17]*mnuSelf-0.25*m0rSelf[4]*uSelf[16]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[15]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[14]*mnuSelf-0.25*m0rSelf[1]*uSelf[13]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[12]*mnuSelf-0.25*m0rSelf[0]*uSelf[12]*mnuSelf+0.5*m1rSelf[12]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherZ-uOtherZ*m0Other) and uCrossOtherZ ... // 
  data->AEM_S(42,36) = 0.25*m0rOther[5]*uOther[17]*mnuOther+0.25*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[3]*uOther[15]*mnuOther+0.25*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[1]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(42,37) = 0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.25*m0rOther[2]*uOther[15]*mnuOther+0.25*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther; 
  data->AEM_S(42,38) = 0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.25*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.25*m0rOther[3]*uOther[13]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther; 
  data->AEM_S(42,39) = 0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther+0.223606797749979*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.25*m0rOther[1]*uOther[14]*mnuOther+0.25*m0rOther[2]*uOther[13]*mnuOther+0.25*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(42,40) = 0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.25*m0rOther[4]*uOther[12]*mnuOther; 
  data->AEM_S(42,41) = 0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[5]*uOther[12]*mnuOther; 
  data->AEM_S(43,36) = 0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.25*m0rOther[2]*uOther[15]*mnuOther+0.25*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther; 
  data->AEM_S(43,37) = 0.25*m0rOther[5]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[16]*mnuOther+0.223606797749979*m0rOther[0]*uOther[16]*mnuOther-0.4472135954999579*m1rOther[16]*mnuOther+0.45*m0rOther[3]*uOther[15]*mnuOther+0.25*m0rOther[2]*uOther[14]*mnuOther+0.45*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(43,38) = 0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther+0.223606797749979*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.25*m0rOther[1]*uOther[14]*mnuOther+0.25*m0rOther[2]*uOther[13]*mnuOther+0.25*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(43,39) = 0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther+0.45*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.223606797749979*m0rOther[4]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.45*m0rOther[3]*uOther[13]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther; 
  data->AEM_S(43,40) = 0.3928571428571428*m0rOther[1]*uOther[16]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther; 
  data->AEM_S(43,41) = 0.25*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.25*m0rOther[5]*uOther[13]*mnuOther; 
  data->AEM_S(44,36) = 0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.25*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.25*m0rOther[3]*uOther[13]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther; 
  data->AEM_S(44,37) = 0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther+0.223606797749979*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.25*m0rOther[1]*uOther[14]*mnuOther+0.25*m0rOther[2]*uOther[13]*mnuOther+0.25*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(44,38) = 0.3928571428571428*m0rOther[5]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.4472135954999579*m1rOther[17]*mnuOther+0.25*m0rOther[4]*uOther[16]*mnuOther+0.45*m0rOther[3]*uOther[15]*mnuOther+0.45*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(44,39) = 0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.45*m0rOther[2]*uOther[15]*mnuOther+0.45*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther; 
  data->AEM_S(44,40) = 0.25*m0rOther[2]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.25*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther; 
  data->AEM_S(44,41) = 0.3928571428571428*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.4472135954999579*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther; 
  data->AEM_S(45,36) = 0.223606797749979*m0rOther[3]*uOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[16]*mnuOther+0.223606797749979*m0rOther[5]*uOther[15]*mnuOther+0.223606797749979*m0rOther[4]*uOther[15]*mnuOther+0.25*m0rOther[0]*uOther[15]*mnuOther-0.5*m1rOther[15]*mnuOther+0.25*m0rOther[1]*uOther[14]*mnuOther+0.25*m0rOther[2]*uOther[13]*mnuOther+0.25*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(45,37) = 0.223606797749979*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[16]*mnuOther+0.45*m0rOther[1]*uOther[15]*mnuOther+0.223606797749979*m0rOther[5]*uOther[14]*mnuOther+0.223606797749979*m0rOther[4]*uOther[14]*mnuOther+0.25*m0rOther[0]*uOther[14]*mnuOther-0.5*m1rOther[14]*mnuOther+0.45*m0rOther[3]*uOther[13]*mnuOther+0.25*m0rOther[2]*uOther[12]*mnuOther; 
  data->AEM_S(45,38) = 0.223606797749979*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[16]*mnuOther+0.45*m0rOther[2]*uOther[15]*mnuOther+0.45*m0rOther[3]*uOther[14]*mnuOther+0.223606797749979*m0rOther[5]*uOther[13]*mnuOther+0.223606797749979*m0rOther[4]*uOther[13]*mnuOther+0.25*m0rOther[0]*uOther[13]*mnuOther-0.5*m1rOther[13]*mnuOther+0.25*m0rOther[1]*uOther[12]*mnuOther; 
  data->AEM_S(45,39) = 0.3928571428571428*m0rOther[5]*uOther[17]*mnuOther+0.2*m0rOther[4]*uOther[17]*mnuOther+0.223606797749979*m0rOther[0]*uOther[17]*mnuOther-0.4472135954999579*m1rOther[17]*mnuOther+0.2*m0rOther[5]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[16]*mnuOther+0.223606797749979*m0rOther[0]*uOther[16]*mnuOther-0.4472135954999579*m1rOther[16]*mnuOther+0.81*m0rOther[3]*uOther[15]*mnuOther+0.45*m0rOther[2]*uOther[14]*mnuOther+0.45*m0rOther[1]*uOther[13]*mnuOther+0.223606797749979*m0rOther[5]*uOther[12]*mnuOther+0.223606797749979*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(45,40) = 0.2*m0rOther[3]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[16]*mnuOther+0.2*m0rOther[5]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(45,41) = 0.3928571428571428*m0rOther[3]*uOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[15]*mnuOther+0.2*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(46,36) = 0.159719141249985*m0rOther[4]*uOther[16]*mnuOther+0.25*m0rOther[0]*uOther[16]*mnuOther-0.5*m1rOther[16]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[13]*mnuOther+0.25*m0rOther[4]*uOther[12]*mnuOther; 
  data->AEM_S(46,37) = 0.3928571428571428*m0rOther[1]*uOther[16]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[13]*mnuOther+0.223606797749979*m0rOther[0]*uOther[13]*mnuOther-0.4472135954999579*m1rOther[13]*mnuOther+0.223606797749979*m0rOther[1]*uOther[12]*mnuOther; 
  data->AEM_S(46,38) = 0.25*m0rOther[2]*uOther[16]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.25*m0rOther[4]*uOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther; 
  data->AEM_S(46,39) = 0.2*m0rOther[3]*uOther[17]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[16]*mnuOther+0.2*m0rOther[5]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(46,40) = 0.25*m0rOther[5]*uOther[17]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[16]*mnuOther+0.159719141249985*m0rOther[0]*uOther[16]*mnuOther-0.31943828249997*m1rOther[16]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[15]*mnuOther+0.25*m0rOther[2]*uOther[14]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[13]*mnuOther+0.159719141249985*m0rOther[4]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther; 
  data->AEM_S(46,41) = 0.25*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[5]*uOther[16]*mnuOther+0.2*m0rOther[3]*uOther[15]*mnuOther; 
  data->AEM_S(47,36) = 0.159719141249985*m0rOther[5]*uOther[17]*mnuOther+0.25*m0rOther[0]*uOther[17]*mnuOther-0.5*m1rOther[17]*mnuOther+0.223606797749979*m0rOther[3]*uOther[15]*mnuOther+0.223606797749979*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[5]*uOther[12]*mnuOther; 
  data->AEM_S(47,37) = 0.25*m0rOther[1]*uOther[17]*mnuOther+0.223606797749979*m0rOther[2]*uOther[15]*mnuOther+0.223606797749979*m0rOther[3]*uOther[14]*mnuOther+0.25*m0rOther[5]*uOther[13]*mnuOther; 
  data->AEM_S(47,38) = 0.3928571428571428*m0rOther[2]*uOther[17]*mnuOther+0.223606797749979*m0rOther[1]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[14]*mnuOther+0.223606797749979*m0rOther[0]*uOther[14]*mnuOther-0.4472135954999579*m1rOther[14]*mnuOther+0.223606797749979*m0rOther[3]*uOther[13]*mnuOther+0.223606797749979*m0rOther[2]*uOther[12]*mnuOther; 
  data->AEM_S(47,39) = 0.3928571428571428*m0rOther[3]*uOther[17]*mnuOther+0.2*m0rOther[3]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[15]*mnuOther+0.2*m0rOther[4]*uOther[15]*mnuOther+0.223606797749979*m0rOther[0]*uOther[15]*mnuOther-0.4472135954999579*m1rOther[15]*mnuOther+0.223606797749979*m0rOther[1]*uOther[14]*mnuOther+0.223606797749979*m0rOther[2]*uOther[13]*mnuOther+0.223606797749979*m0rOther[3]*uOther[12]*mnuOther; 
  data->AEM_S(47,40) = 0.25*m0rOther[4]*uOther[17]*mnuOther+0.25*m0rOther[5]*uOther[16]*mnuOther+0.2*m0rOther[3]*uOther[15]*mnuOther; 
  data->AEM_S(47,41) = 0.5357142857142857*m0rOther[5]*uOther[17]*mnuOther+0.159719141249985*m0rOther[0]*uOther[17]*mnuOther-0.31943828249997*m1rOther[17]*mnuOther+0.25*m0rOther[4]*uOther[16]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[15]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[14]*mnuOther+0.25*m0rOther[1]*uOther[13]*mnuOther+0.159719141249985*m0rOther[5]*uOther[12]*mnuOther+0.25*m0rOther[0]*uOther[12]*mnuOther-0.5*m1rOther[12]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of momentum relaxation. 
  m1Relax[12] += (m1rOther[12]-1.0*m1rSelf[12])*betaGreenep1*mnuSelf+m1rSelf[12]*mnuSelf-1.0*m1rOther[12]*mnuOther; 
  m1Relax[13] += (m1rOther[13]-1.0*m1rSelf[13])*betaGreenep1*mnuSelf+m1rSelf[13]*mnuSelf-1.0*m1rOther[13]*mnuOther; 
  m1Relax[14] += (m1rOther[14]-1.0*m1rSelf[14])*betaGreenep1*mnuSelf+m1rSelf[14]*mnuSelf-1.0*m1rOther[14]*mnuOther; 
  m1Relax[15] += (m1rOther[15]-1.0*m1rSelf[15])*betaGreenep1*mnuSelf+m1rSelf[15]*mnuSelf-1.0*m1rOther[15]*mnuOther; 
  m1Relax[16] += (m1rOther[16]-1.0*m1rSelf[16])*betaGreenep1*mnuSelf+m1rSelf[16]*mnuSelf-1.0*m1rOther[16]*mnuOther; 
  m1Relax[17] += (m1rOther[17]-1.0*m1rSelf[17])*betaGreenep1*mnuSelf+m1rSelf[17]*mnuSelf-1.0*m1rOther[17]*mnuOther; 
 
  double ucMSelf[6]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
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
  data->AEM_S(42,18) = 0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(42,19) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(42,20) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(42,21) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(42,22) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(42,23) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(43,18) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(43,19) = 0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(43,20) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(43,21) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(43,22) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(44,18) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(44,19) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(44,20) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(44,21) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(44,23) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(45,18) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(45,19) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(45,20) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(45,21) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(45,22) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(45,23) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(46,18) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(46,19) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(46,21) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(46,22) = 0.31943828249997*ucMSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(47,18) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(47,20) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(47,21) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(47,23) = 0.31943828249997*ucMSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  double ucMOther[6]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
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
  data->AEM_S(42,42) = (-0.5*ucMOther[0]*mnuOther)-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(42,43) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(42,44) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(42,45) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(42,46) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(42,47) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(43,42) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(43,43) = (-0.4472135954999579*ucMOther[4]*mnuOther)-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(43,44) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(43,45) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(43,46) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(44,42) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(44,43) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(44,44) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(44,45) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(44,47) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(45,42) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(45,43) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(45,44) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(45,45) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.4472135954999579*ucMOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(45,46) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(45,47) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(46,42) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(46,43) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(46,45) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(46,46) = (-0.31943828249997*ucMOther[4]*mnuOther)-0.9583148474999099*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(47,42) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(47,44) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(47,45) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(47,47) = (-0.31943828249997*ucMOther[5]*mnuOther)-0.9583148474999099*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[6]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
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
  for (unsigned short int vd=0; vd<3; vd++) 
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
  for (unsigned short int vd=0; vd<3; vd++) 
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
  data->AEM_S.block<6,12>(24,6).setZero(); 
  data->AEM_S.block<12,6>(30,0).setZero(); 
  data->AEM_S.block<6,6>(30,12).setZero(); 
  data->AEM_S.block<6,6>(36,6).setZero(); 
  data->AEM_S.block<6,12>(24,30).setZero(); 
  data->AEM_S.block<12,6>(30,24).setZero(); 
  data->AEM_S.block<6,6>(30,36).setZero(); 
  data->AEM_S.block<6,6>(36,30).setZero(); 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM1sum[10],mnuM1sum[11],mnuM1sum[12],mnuM1sum[13],mnuM1sum[14],mnuM1sum[15],mnuM1sum[16],mnuM1sum[17],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m1Relax[8],m1Relax[9],m1Relax[10],m1Relax[11],m1Relax[12],m1Relax[13],m1Relax[14],m1Relax[15],m1Relax[16],m1Relax[17],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,18,1) = data->u_S.segment<18>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,6,1) = data->u_S.segment<6>(18); 
 
  Eigen::Map<VectorXd>(uCrossOther,18,1) = data->u_S.segment<18>(24); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,6,1) = data->u_S.segment<6>(42); 
 
} 
 
