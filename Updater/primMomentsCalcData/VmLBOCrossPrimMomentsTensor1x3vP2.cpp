#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmLBOCrossPrimMoments1x3vTensor_P2(binOpData_t *data, binOpData_t *dataDiv,const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[9]; 
  double m2rSelf[3]; 
  double cMrSelf[9]; 
  double cErSelf[3]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m1rSelf[0] = m1Self[0]; 
    m1rSelf[1] = 0.0; 
    m1rSelf[2] = 0.0; 
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = 0.0; 
    cMrSelf[2] = 0.0; 
    m1rSelf[3] = m1Self[3]; 
    m1rSelf[4] = 0.0; 
    m1rSelf[5] = 0.0; 
    cMrSelf[3] = cMSelf[3]; 
    cMrSelf[4] = 0.0; 
    cMrSelf[5] = 0.0; 
    m1rSelf[6] = m1Self[6]; 
    m1rSelf[7] = 0.0; 
    m1rSelf[8] = 0.0; 
    cMrSelf[6] = cMSelf[6]; 
    cMrSelf[7] = 0.0; 
    cMrSelf[8] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = 0.0; 
    cErSelf[2] = 0.0; 
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
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = cMSelf[1]; 
    cMrSelf[2] = cMSelf[2]; 
    cMrSelf[3] = cMSelf[3]; 
    cMrSelf[4] = cMSelf[4]; 
    cMrSelf[5] = cMSelf[5]; 
    cMrSelf[6] = cMSelf[6]; 
    cMrSelf[7] = cMSelf[7]; 
    cMrSelf[8] = cMSelf[8]; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = cESelf[1]; 
    cErSelf[2] = cESelf[2]; 
  } 
 
  if (1.58113883008419*m0Other[2]-1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.58113883008419*m0Other[2]+1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[3]; 
  double m1rOther[9]; 
  double m2rOther[3]; 
  double cMrOther[9]; 
  double cErOther[3]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m1rOther[0] = m1Other[0]; 
    m1rOther[1] = 0.0; 
    m1rOther[2] = 0.0; 
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = 0.0; 
    cMrOther[2] = 0.0; 
    m1rOther[3] = m1Other[3]; 
    m1rOther[4] = 0.0; 
    m1rOther[5] = 0.0; 
    cMrOther[3] = cMOther[3]; 
    cMrOther[4] = 0.0; 
    cMrOther[5] = 0.0; 
    m1rOther[6] = m1Other[6]; 
    m1rOther[7] = 0.0; 
    m1rOther[8] = 0.0; 
    cMrOther[6] = cMOther[6]; 
    cMrOther[7] = 0.0; 
    cMrOther[8] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = 0.0; 
    cErOther[2] = 0.0; 
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
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = cMOther[1]; 
    cMrOther[2] = cMOther[2]; 
    cMrOther[3] = cMOther[3]; 
    cMrOther[4] = cMOther[4]; 
    cMrOther[5] = cMOther[5]; 
    cMrOther[6] = cMOther[6]; 
    cMrOther[7] = cMOther[7]; 
    cMrOther[8] = cMOther[8]; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = cEOther[1]; 
    cErOther[2] = cEOther[2]; 
  } 
 
  // Declare Eigen matrix and vectors for weak system. 
  data->AEM_S = Eigen::MatrixXd::Zero(24,24); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double mnuM1sum[9]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<9; vd++) 
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
  data->AEM_S(0,9) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
  data->AEM_S(0,10) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(0,11) = -0.7071067811865475*cMrSelf[2]*mnuSelf; 
  data->AEM_S(1,9) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(1,10) = (-0.6324555320336759*cMrSelf[2]*mnuSelf)-0.7071067811865475*cMrSelf[0]*mnuSelf; 
  data->AEM_S(1,11) = -0.6324555320336759*cMrSelf[1]*mnuSelf; 
  data->AEM_S(2,9) = -0.7071067811865475*cMrSelf[2]*mnuSelf; 
  data->AEM_S(2,10) = -0.6324555320336759*cMrSelf[1]*mnuSelf; 
  data->AEM_S(2,11) = (-0.4517539514526256*cMrSelf[2]*mnuSelf)-0.7071067811865475*cMrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,12) = 0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(0,13) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(0,14) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(1,12) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(1,13) = 0.6324555320336759*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(1,14) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(2,12) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(2,13) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(2,14) = 0.4517539514526256*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,21) = -0.7071067811865475*cMrOther[0]*mnuOther; 
  data->AEM_S(0,22) = -0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(0,23) = -0.7071067811865475*cMrOther[2]*mnuOther; 
  data->AEM_S(1,21) = -0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(1,22) = (-0.6324555320336759*cMrOther[2]*mnuOther)-0.7071067811865475*cMrOther[0]*mnuOther; 
  data->AEM_S(1,23) = -0.6324555320336759*cMrOther[1]*mnuOther; 
  data->AEM_S(2,21) = -0.7071067811865475*cMrOther[2]*mnuOther; 
  data->AEM_S(2,22) = -0.6324555320336759*cMrOther[1]*mnuOther; 
  data->AEM_S(2,23) = (-0.4517539514526256*cMrOther[2]*mnuOther)-0.7071067811865475*cMrOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(9,0) = 0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(9,1) = 0.7071067811865475*m1rSelf[1]*mnuSelf; 
  data->AEM_S(9,2) = 0.7071067811865475*m1rSelf[2]*mnuSelf; 
  data->AEM_S(10,0) = 0.7071067811865475*m1rSelf[1]*mnuSelf; 
  data->AEM_S(10,1) = 0.6324555320336759*m1rSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(10,2) = 0.6324555320336759*m1rSelf[1]*mnuSelf; 
  data->AEM_S(11,0) = 0.7071067811865475*m1rSelf[2]*mnuSelf; 
  data->AEM_S(11,1) = 0.6324555320336759*m1rSelf[1]*mnuSelf; 
  data->AEM_S(11,2) = 0.4517539514526256*m1rSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(9,12) = 0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(9,13) = 0.7071067811865475*m1rOther[1]*mnuOther; 
  data->AEM_S(9,14) = 0.7071067811865475*m1rOther[2]*mnuOther; 
  data->AEM_S(10,12) = 0.7071067811865475*m1rOther[1]*mnuOther; 
  data->AEM_S(10,13) = 0.6324555320336759*m1rOther[2]*mnuOther+0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(10,14) = 0.6324555320336759*m1rOther[1]*mnuOther; 
  data->AEM_S(11,12) = 0.7071067811865475*m1rOther[2]*mnuOther; 
  data->AEM_S(11,13) = 0.6324555320336759*m1rOther[1]*mnuOther; 
  data->AEM_S(11,14) = 0.4517539514526256*m1rOther[2]*mnuOther+0.7071067811865475*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(3,3) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(3,4) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(3,5) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(4,3) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(4,4) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(4,5) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,3) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(5,4) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(5,5) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(3,9) = -0.7071067811865475*cMrSelf[3]*mnuSelf; 
  data->AEM_S(3,10) = -0.7071067811865475*cMrSelf[4]*mnuSelf; 
  data->AEM_S(3,11) = -0.7071067811865475*cMrSelf[5]*mnuSelf; 
  data->AEM_S(4,9) = -0.7071067811865475*cMrSelf[4]*mnuSelf; 
  data->AEM_S(4,10) = (-0.6324555320336759*cMrSelf[5]*mnuSelf)-0.7071067811865475*cMrSelf[3]*mnuSelf; 
  data->AEM_S(4,11) = -0.6324555320336759*cMrSelf[4]*mnuSelf; 
  data->AEM_S(5,9) = -0.7071067811865475*cMrSelf[5]*mnuSelf; 
  data->AEM_S(5,10) = -0.6324555320336759*cMrSelf[4]*mnuSelf; 
  data->AEM_S(5,11) = (-0.4517539514526256*cMrSelf[5]*mnuSelf)-0.7071067811865475*cMrSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(3,15) = 0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(3,16) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(3,17) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(4,15) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(4,16) = 0.6324555320336759*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(4,17) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(5,15) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(5,16) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(5,17) = 0.4517539514526256*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(3,21) = -0.7071067811865475*cMrOther[3]*mnuOther; 
  data->AEM_S(3,22) = -0.7071067811865475*cMrOther[4]*mnuOther; 
  data->AEM_S(3,23) = -0.7071067811865475*cMrOther[5]*mnuOther; 
  data->AEM_S(4,21) = -0.7071067811865475*cMrOther[4]*mnuOther; 
  data->AEM_S(4,22) = (-0.6324555320336759*cMrOther[5]*mnuOther)-0.7071067811865475*cMrOther[3]*mnuOther; 
  data->AEM_S(4,23) = -0.6324555320336759*cMrOther[4]*mnuOther; 
  data->AEM_S(5,21) = -0.7071067811865475*cMrOther[5]*mnuOther; 
  data->AEM_S(5,22) = -0.6324555320336759*cMrOther[4]*mnuOther; 
  data->AEM_S(5,23) = (-0.4517539514526256*cMrOther[5]*mnuOther)-0.7071067811865475*cMrOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfY and uCrossSelfY ... // 
  data->AEM_S(9,3) = 0.7071067811865475*m1rSelf[3]*mnuSelf; 
  data->AEM_S(9,4) = 0.7071067811865475*m1rSelf[4]*mnuSelf; 
  data->AEM_S(9,5) = 0.7071067811865475*m1rSelf[5]*mnuSelf; 
  data->AEM_S(10,3) = 0.7071067811865475*m1rSelf[4]*mnuSelf; 
  data->AEM_S(10,4) = 0.6324555320336759*m1rSelf[5]*mnuSelf+0.7071067811865475*m1rSelf[3]*mnuSelf; 
  data->AEM_S(10,5) = 0.6324555320336759*m1rSelf[4]*mnuSelf; 
  data->AEM_S(11,3) = 0.7071067811865475*m1rSelf[5]*mnuSelf; 
  data->AEM_S(11,4) = 0.6324555320336759*m1rSelf[4]*mnuSelf; 
  data->AEM_S(11,5) = 0.4517539514526256*m1rSelf[5]*mnuSelf+0.7071067811865475*m1rSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherY and uCrossOtherY ... // 
  data->AEM_S(9,15) = 0.7071067811865475*m1rOther[3]*mnuOther; 
  data->AEM_S(9,16) = 0.7071067811865475*m1rOther[4]*mnuOther; 
  data->AEM_S(9,17) = 0.7071067811865475*m1rOther[5]*mnuOther; 
  data->AEM_S(10,15) = 0.7071067811865475*m1rOther[4]*mnuOther; 
  data->AEM_S(10,16) = 0.6324555320336759*m1rOther[5]*mnuOther+0.7071067811865475*m1rOther[3]*mnuOther; 
  data->AEM_S(10,17) = 0.6324555320336759*m1rOther[4]*mnuOther; 
  data->AEM_S(11,15) = 0.7071067811865475*m1rOther[5]*mnuOther; 
  data->AEM_S(11,16) = 0.6324555320336759*m1rOther[4]*mnuOther; 
  data->AEM_S(11,17) = 0.4517539514526256*m1rOther[5]*mnuOther+0.7071067811865475*m1rOther[3]*mnuOther; 
 
  // ... Contribution to RHS vector from component 2 of mnuM1Self+mnuM1Other. 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(6,6) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(6,7) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(6,8) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(7,6) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(7,7) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(7,8) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,6) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,7) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,8) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to momentum conservation (self) ... // 
  data->AEM_S(6,9) = -0.7071067811865475*cMrSelf[6]*mnuSelf; 
  data->AEM_S(6,10) = -0.7071067811865475*cMrSelf[7]*mnuSelf; 
  data->AEM_S(6,11) = -0.7071067811865475*cMrSelf[8]*mnuSelf; 
  data->AEM_S(7,9) = -0.7071067811865475*cMrSelf[7]*mnuSelf; 
  data->AEM_S(7,10) = (-0.6324555320336759*cMrSelf[8]*mnuSelf)-0.7071067811865475*cMrSelf[6]*mnuSelf; 
  data->AEM_S(7,11) = -0.6324555320336759*cMrSelf[7]*mnuSelf; 
  data->AEM_S(8,9) = -0.7071067811865475*cMrSelf[8]*mnuSelf; 
  data->AEM_S(8,10) = -0.6324555320336759*cMrSelf[7]*mnuSelf; 
  data->AEM_S(8,11) = (-0.4517539514526256*cMrSelf[8]*mnuSelf)-0.7071067811865475*cMrSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(6,18) = 0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(6,19) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(6,20) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(7,18) = 0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(7,19) = 0.6324555320336759*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(7,20) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(8,18) = 0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(8,19) = 0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(8,20) = 0.4517539514526256*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(6,21) = -0.7071067811865475*cMrOther[6]*mnuOther; 
  data->AEM_S(6,22) = -0.7071067811865475*cMrOther[7]*mnuOther; 
  data->AEM_S(6,23) = -0.7071067811865475*cMrOther[8]*mnuOther; 
  data->AEM_S(7,21) = -0.7071067811865475*cMrOther[7]*mnuOther; 
  data->AEM_S(7,22) = (-0.6324555320336759*cMrOther[8]*mnuOther)-0.7071067811865475*cMrOther[6]*mnuOther; 
  data->AEM_S(7,23) = -0.6324555320336759*cMrOther[7]*mnuOther; 
  data->AEM_S(8,21) = -0.7071067811865475*cMrOther[8]*mnuOther; 
  data->AEM_S(8,22) = -0.6324555320336759*cMrOther[7]*mnuOther; 
  data->AEM_S(8,23) = (-0.4517539514526256*cMrOther[8]*mnuOther)-0.7071067811865475*cMrOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfZ and uCrossSelfZ ... // 
  data->AEM_S(9,6) = 0.7071067811865475*m1rSelf[6]*mnuSelf; 
  data->AEM_S(9,7) = 0.7071067811865475*m1rSelf[7]*mnuSelf; 
  data->AEM_S(9,8) = 0.7071067811865475*m1rSelf[8]*mnuSelf; 
  data->AEM_S(10,6) = 0.7071067811865475*m1rSelf[7]*mnuSelf; 
  data->AEM_S(10,7) = 0.6324555320336759*m1rSelf[8]*mnuSelf+0.7071067811865475*m1rSelf[6]*mnuSelf; 
  data->AEM_S(10,8) = 0.6324555320336759*m1rSelf[7]*mnuSelf; 
  data->AEM_S(11,6) = 0.7071067811865475*m1rSelf[8]*mnuSelf; 
  data->AEM_S(11,7) = 0.6324555320336759*m1rSelf[7]*mnuSelf; 
  data->AEM_S(11,8) = 0.4517539514526256*m1rSelf[8]*mnuSelf+0.7071067811865475*m1rSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherZ and uCrossOtherZ ... // 
  data->AEM_S(9,18) = 0.7071067811865475*m1rOther[6]*mnuOther; 
  data->AEM_S(9,19) = 0.7071067811865475*m1rOther[7]*mnuOther; 
  data->AEM_S(9,20) = 0.7071067811865475*m1rOther[8]*mnuOther; 
  data->AEM_S(10,18) = 0.7071067811865475*m1rOther[7]*mnuOther; 
  data->AEM_S(10,19) = 0.6324555320336759*m1rOther[8]*mnuOther+0.7071067811865475*m1rOther[6]*mnuOther; 
  data->AEM_S(10,20) = 0.6324555320336759*m1rOther[7]*mnuOther; 
  data->AEM_S(11,18) = 0.7071067811865475*m1rOther[8]*mnuOther; 
  data->AEM_S(11,19) = 0.6324555320336759*m1rOther[7]*mnuOther; 
  data->AEM_S(11,20) = 0.4517539514526256*m1rOther[8]*mnuOther+0.7071067811865475*m1rOther[6]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(9,9) = 2.121320343559642*m0rSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
  data->AEM_S(9,10) = 2.121320343559642*m0rSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(9,11) = 2.121320343559642*m0rSelf[2]*mnuSelf-0.7071067811865475*cErSelf[2]*mnuSelf; 
  data->AEM_S(10,9) = 2.121320343559642*m0rSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(10,10) = 1.897366596101028*m0rSelf[2]*mnuSelf-0.6324555320336759*cErSelf[2]*mnuSelf+2.121320343559642*m0rSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
  data->AEM_S(10,11) = 1.897366596101028*m0rSelf[1]*mnuSelf-0.6324555320336759*cErSelf[1]*mnuSelf; 
  data->AEM_S(11,9) = 2.121320343559642*m0rSelf[2]*mnuSelf-0.7071067811865475*cErSelf[2]*mnuSelf; 
  data->AEM_S(11,10) = 1.897366596101028*m0rSelf[1]*mnuSelf-0.6324555320336759*cErSelf[1]*mnuSelf; 
  data->AEM_S(11,11) = 1.355261854357877*m0rSelf[2]*mnuSelf-0.4517539514526256*cErSelf[2]*mnuSelf+2.121320343559642*m0rSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(9,21) = 2.121320343559642*m0rOther[0]*mnuOther-0.7071067811865475*cErOther[0]*mnuOther; 
  data->AEM_S(9,22) = 2.121320343559642*m0rOther[1]*mnuOther-0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(9,23) = 2.121320343559642*m0rOther[2]*mnuOther-0.7071067811865475*cErOther[2]*mnuOther; 
  data->AEM_S(10,21) = 2.121320343559642*m0rOther[1]*mnuOther-0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(10,22) = 1.897366596101028*m0rOther[2]*mnuOther-0.6324555320336759*cErOther[2]*mnuOther+2.121320343559642*m0rOther[0]*mnuOther-0.7071067811865475*cErOther[0]*mnuOther; 
  data->AEM_S(10,23) = 1.897366596101028*m0rOther[1]*mnuOther-0.6324555320336759*cErOther[1]*mnuOther; 
  data->AEM_S(11,21) = 2.121320343559642*m0rOther[2]*mnuOther-0.7071067811865475*cErOther[2]*mnuOther; 
  data->AEM_S(11,22) = 1.897366596101028*m0rOther[1]*mnuOther-0.6324555320336759*cErOther[1]*mnuOther; 
  data->AEM_S(11,23) = 1.355261854357877*m0rOther[2]*mnuOther-0.4517539514526256*cErOther[2]*mnuOther+2.121320343559642*m0rOther[0]*mnuOther-0.7071067811865475*cErOther[0]*mnuOther; 
 
  double mnuM2sum[3]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2rSelf[0]*mnuSelf+m2rOther[0]*mnuOther; 
  mnuM2sum[1] = m2rSelf[1]*mnuSelf+m2rOther[1]*mnuOther; 
  mnuM2sum[2] = m2rSelf[2]*mnuSelf+m2rOther[2]*mnuOther; 
 
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
 
  double m1EffD[9]; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(12,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(12,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,1) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(13,2) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,0) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,1) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(14,2) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(12,9) = -0.7071067811865475*cMrSelf[0]*mnuSelf; 
  data->AEM_S(12,10) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(12,11) = -0.7071067811865475*cMrSelf[2]*mnuSelf; 
  data->AEM_S(13,9) = -0.7071067811865475*cMrSelf[1]*mnuSelf; 
  data->AEM_S(13,10) = (-0.6324555320336759*cMrSelf[2]*mnuSelf)-0.7071067811865475*cMrSelf[0]*mnuSelf; 
  data->AEM_S(13,11) = -0.6324555320336759*cMrSelf[1]*mnuSelf; 
  data->AEM_S(14,9) = -0.7071067811865475*cMrSelf[2]*mnuSelf; 
  data->AEM_S(14,10) = -0.6324555320336759*cMrSelf[1]*mnuSelf; 
  data->AEM_S(14,11) = (-0.4517539514526256*cMrSelf[2]*mnuSelf)-0.7071067811865475*cMrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(12,12) = -0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(12,13) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(12,14) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(13,12) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(13,13) = (-0.6324555320336759*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(13,14) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(14,12) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(14,13) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(14,14) = (-0.4517539514526256*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(12,21) = 0.7071067811865475*cMrOther[0]*mnuOther; 
  data->AEM_S(12,22) = 0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(12,23) = 0.7071067811865475*cMrOther[2]*mnuOther; 
  data->AEM_S(13,21) = 0.7071067811865475*cMrOther[1]*mnuOther; 
  data->AEM_S(13,22) = 0.6324555320336759*cMrOther[2]*mnuOther+0.7071067811865475*cMrOther[0]*mnuOther; 
  data->AEM_S(13,23) = 0.6324555320336759*cMrOther[1]*mnuOther; 
  data->AEM_S(14,21) = 0.7071067811865475*cMrOther[2]*mnuOther; 
  data->AEM_S(14,22) = 0.6324555320336759*cMrOther[1]*mnuOther; 
  data->AEM_S(14,23) = 0.4517539514526256*cMrOther[2]*mnuOther+0.7071067811865475*cMrOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(21,0) = (-0.5*m0rSelf[2]*uSelf[2]*mnuSelf)-0.5*m0rSelf[1]*uSelf[1]*mnuSelf-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(21,1) = (-0.4472135954999579*m0rSelf[1]*uSelf[2]*mnuSelf)-0.4472135954999579*uSelf[1]*m0rSelf[2]*mnuSelf-0.5*m0rSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m1rSelf[1]*mnuSelf-0.5*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(21,2) = (-0.31943828249997*m0rSelf[2]*uSelf[2]*mnuSelf)-0.5*m0rSelf[0]*uSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[2]*mnuSelf-0.5*uSelf[0]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(22,0) = (-0.4472135954999579*m0rSelf[1]*uSelf[2]*mnuSelf)-0.4472135954999579*uSelf[1]*m0rSelf[2]*mnuSelf-0.5*m0rSelf[0]*uSelf[1]*mnuSelf+0.7071067811865475*m1rSelf[1]*mnuSelf-0.5*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,1) = (-0.7857142857142857*m0rSelf[2]*uSelf[2]*mnuSelf)-0.4472135954999579*m0rSelf[0]*uSelf[2]*mnuSelf+0.6324555320336759*m1rSelf[2]*mnuSelf-0.4472135954999579*uSelf[0]*m0rSelf[2]*mnuSelf-0.9*m0rSelf[1]*uSelf[1]*mnuSelf-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
  data->AEM_S(22,2) = (-0.7857142857142857*m0rSelf[1]*uSelf[2]*mnuSelf)-0.7857142857142857*uSelf[1]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[1]*mnuSelf+0.6324555320336759*m1rSelf[1]*mnuSelf-0.4472135954999579*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,0) = (-0.31943828249997*m0rSelf[2]*uSelf[2]*mnuSelf)-0.5*m0rSelf[0]*uSelf[2]*mnuSelf+0.7071067811865475*m1rSelf[2]*mnuSelf-0.5*uSelf[0]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(23,1) = (-0.7857142857142857*m0rSelf[1]*uSelf[2]*mnuSelf)-0.7857142857142857*uSelf[1]*m0rSelf[2]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[1]*mnuSelf+0.6324555320336759*m1rSelf[1]*mnuSelf-0.4472135954999579*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,2) = (-1.071428571428571*m0rSelf[2]*uSelf[2]*mnuSelf)-0.31943828249997*m0rSelf[0]*uSelf[2]*mnuSelf+0.4517539514526256*m1rSelf[2]*mnuSelf-0.31943828249997*uSelf[0]*m0rSelf[2]*mnuSelf-0.7857142857142857*m0rSelf[1]*uSelf[1]*mnuSelf-0.5*m0rSelf[0]*uSelf[0]*mnuSelf+0.7071067811865475*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(21,12) = 0.5*m0rOther[2]*uOther[2]*mnuOther+0.5*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(21,13) = 0.4472135954999579*m0rOther[1]*uOther[2]*mnuOther+0.4472135954999579*uOther[1]*m0rOther[2]*mnuOther+0.5*m0rOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1rOther[1]*mnuOther+0.5*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(21,14) = 0.31943828249997*m0rOther[2]*uOther[2]*mnuOther+0.5*m0rOther[0]*uOther[2]*mnuOther-0.7071067811865475*m1rOther[2]*mnuOther+0.5*uOther[0]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(22,12) = 0.4472135954999579*m0rOther[1]*uOther[2]*mnuOther+0.4472135954999579*uOther[1]*m0rOther[2]*mnuOther+0.5*m0rOther[0]*uOther[1]*mnuOther-0.7071067811865475*m1rOther[1]*mnuOther+0.5*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(22,13) = 0.7857142857142857*m0rOther[2]*uOther[2]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[2]*mnuOther-0.6324555320336759*m1rOther[2]*mnuOther+0.4472135954999579*uOther[0]*m0rOther[2]*mnuOther+0.9*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1rOther[0]*mnuOther; 
  data->AEM_S(22,14) = 0.7857142857142857*m0rOther[1]*uOther[2]*mnuOther+0.7857142857142857*uOther[1]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[1]*mnuOther-0.6324555320336759*m1rOther[1]*mnuOther+0.4472135954999579*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(23,12) = 0.31943828249997*m0rOther[2]*uOther[2]*mnuOther+0.5*m0rOther[0]*uOther[2]*mnuOther-0.7071067811865475*m1rOther[2]*mnuOther+0.5*uOther[0]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(23,13) = 0.7857142857142857*m0rOther[1]*uOther[2]*mnuOther+0.7857142857142857*uOther[1]*m0rOther[2]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[1]*mnuOther-0.6324555320336759*m1rOther[1]*mnuOther+0.4472135954999579*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(23,14) = 1.071428571428571*m0rOther[2]*uOther[2]*mnuOther+0.31943828249997*m0rOther[0]*uOther[2]*mnuOther-0.4517539514526256*m1rOther[2]*mnuOther+0.31943828249997*uOther[0]*m0rOther[2]*mnuOther+0.7857142857142857*m0rOther[1]*uOther[1]*mnuOther+0.5*m0rOther[0]*uOther[0]*mnuOther-0.7071067811865475*m1rOther[0]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfX-m0Self*m1OtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[0] = 0.7071067811865475*m0rOther[2]*m1rSelf[2]-0.7071067811865475*m0rSelf[2]*m1rOther[2]+0.7071067811865475*m0rOther[1]*m1rSelf[1]-0.7071067811865475*m0rSelf[1]*m1rOther[1]+0.7071067811865475*m0rOther[0]*m1rSelf[0]-0.7071067811865475*m0rSelf[0]*m1rOther[0]; 
  m1EffD[1] = 0.6324555320336759*m0rOther[1]*m1rSelf[2]-0.6324555320336759*m0rSelf[1]*m1rOther[2]-0.6324555320336759*m1rOther[1]*m0rSelf[2]+0.6324555320336759*m1rSelf[1]*m0rOther[2]+0.7071067811865475*m0rOther[0]*m1rSelf[1]-0.7071067811865475*m0rSelf[0]*m1rOther[1]-0.7071067811865475*m1rOther[0]*m0rSelf[1]+0.7071067811865475*m1rSelf[0]*m0rOther[1]; 
  m1EffD[2] = 0.4517539514526256*m0rOther[2]*m1rSelf[2]+0.7071067811865475*m0rOther[0]*m1rSelf[2]-0.4517539514526256*m0rSelf[2]*m1rOther[2]-0.7071067811865475*m0rSelf[0]*m1rOther[2]-0.7071067811865475*m1rOther[0]*m0rSelf[2]+0.7071067811865475*m1rSelf[0]*m0rOther[2]+0.6324555320336759*m0rOther[1]*m1rSelf[1]-0.6324555320336759*m0rSelf[1]*m1rOther[1]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(3,3); 
  dataDiv->AEM_S(0,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf+0.7071067811865475*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf+0.7071067811865475*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(0,2) = 0.7071067811865475*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf+0.7071067811865475*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf+0.6324555320336759*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(1,2) = 0.6324555320336759*m0rSelf[1]*mnuSelf+0.6324555320336759*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(2,0) = 0.7071067811865475*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,1) = 0.6324555320336759*m0rSelf[1]*mnuSelf+0.6324555320336759*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(2,2) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf+0.4517539514526256*m0rOther[2]*mnuOther+0.7071067811865475*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[0],m1EffD[1],m1EffD[2]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+0,3,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += (-2.0*m1EffD[0]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += (-2.0*m1EffD[1]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += (-2.0*m1EffD[2]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfY ... // 
  data->AEM_S(15,3) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(15,4) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,5) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,3) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,4) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(16,5) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,3) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,4) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,5) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(15,9) = -0.7071067811865475*cMrSelf[3]*mnuSelf; 
  data->AEM_S(15,10) = -0.7071067811865475*cMrSelf[4]*mnuSelf; 
  data->AEM_S(15,11) = -0.7071067811865475*cMrSelf[5]*mnuSelf; 
  data->AEM_S(16,9) = -0.7071067811865475*cMrSelf[4]*mnuSelf; 
  data->AEM_S(16,10) = (-0.6324555320336759*cMrSelf[5]*mnuSelf)-0.7071067811865475*cMrSelf[3]*mnuSelf; 
  data->AEM_S(16,11) = -0.6324555320336759*cMrSelf[4]*mnuSelf; 
  data->AEM_S(17,9) = -0.7071067811865475*cMrSelf[5]*mnuSelf; 
  data->AEM_S(17,10) = -0.6324555320336759*cMrSelf[4]*mnuSelf; 
  data->AEM_S(17,11) = (-0.4517539514526256*cMrSelf[5]*mnuSelf)-0.7071067811865475*cMrSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(15,15) = -0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(15,16) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(15,17) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(16,15) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(16,16) = (-0.6324555320336759*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(16,17) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(17,15) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(17,16) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(17,17) = (-0.4517539514526256*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(15,21) = 0.7071067811865475*cMrOther[3]*mnuOther; 
  data->AEM_S(15,22) = 0.7071067811865475*cMrOther[4]*mnuOther; 
  data->AEM_S(15,23) = 0.7071067811865475*cMrOther[5]*mnuOther; 
  data->AEM_S(16,21) = 0.7071067811865475*cMrOther[4]*mnuOther; 
  data->AEM_S(16,22) = 0.6324555320336759*cMrOther[5]*mnuOther+0.7071067811865475*cMrOther[3]*mnuOther; 
  data->AEM_S(16,23) = 0.6324555320336759*cMrOther[4]*mnuOther; 
  data->AEM_S(17,21) = 0.7071067811865475*cMrOther[5]*mnuOther; 
  data->AEM_S(17,22) = 0.6324555320336759*cMrOther[4]*mnuOther; 
  data->AEM_S(17,23) = 0.4517539514526256*cMrOther[5]*mnuOther+0.7071067811865475*cMrOther[3]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfY-uSelfY*m0Self) and uCrossSelfY ... // 
  data->AEM_S(21,3) = (-0.5*m0rSelf[2]*uSelf[5]*mnuSelf)-0.5*m0rSelf[1]*uSelf[4]*mnuSelf-0.5*m0rSelf[0]*uSelf[3]*mnuSelf+0.7071067811865475*m1rSelf[3]*mnuSelf; 
  data->AEM_S(21,4) = (-0.4472135954999579*m0rSelf[1]*uSelf[5]*mnuSelf)-0.4472135954999579*m0rSelf[2]*uSelf[4]*mnuSelf-0.5*m0rSelf[0]*uSelf[4]*mnuSelf+0.7071067811865475*m1rSelf[4]*mnuSelf-0.5*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(21,5) = (-0.31943828249997*m0rSelf[2]*uSelf[5]*mnuSelf)-0.5*m0rSelf[0]*uSelf[5]*mnuSelf+0.7071067811865475*m1rSelf[5]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[4]*mnuSelf-0.5*m0rSelf[2]*uSelf[3]*mnuSelf; 
  data->AEM_S(22,3) = (-0.4472135954999579*m0rSelf[1]*uSelf[5]*mnuSelf)-0.4472135954999579*m0rSelf[2]*uSelf[4]*mnuSelf-0.5*m0rSelf[0]*uSelf[4]*mnuSelf+0.7071067811865475*m1rSelf[4]*mnuSelf-0.5*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(22,4) = (-0.7857142857142857*m0rSelf[2]*uSelf[5]*mnuSelf)-0.4472135954999579*m0rSelf[0]*uSelf[5]*mnuSelf+0.6324555320336759*m1rSelf[5]*mnuSelf-0.9*m0rSelf[1]*uSelf[4]*mnuSelf-0.4472135954999579*m0rSelf[2]*uSelf[3]*mnuSelf-0.5*m0rSelf[0]*uSelf[3]*mnuSelf+0.7071067811865475*m1rSelf[3]*mnuSelf; 
  data->AEM_S(22,5) = (-0.7857142857142857*m0rSelf[1]*uSelf[5]*mnuSelf)-0.7857142857142857*m0rSelf[2]*uSelf[4]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[4]*mnuSelf+0.6324555320336759*m1rSelf[4]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(23,3) = (-0.31943828249997*m0rSelf[2]*uSelf[5]*mnuSelf)-0.5*m0rSelf[0]*uSelf[5]*mnuSelf+0.7071067811865475*m1rSelf[5]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[4]*mnuSelf-0.5*m0rSelf[2]*uSelf[3]*mnuSelf; 
  data->AEM_S(23,4) = (-0.7857142857142857*m0rSelf[1]*uSelf[5]*mnuSelf)-0.7857142857142857*m0rSelf[2]*uSelf[4]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[4]*mnuSelf+0.6324555320336759*m1rSelf[4]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[3]*mnuSelf; 
  data->AEM_S(23,5) = (-1.071428571428571*m0rSelf[2]*uSelf[5]*mnuSelf)-0.31943828249997*m0rSelf[0]*uSelf[5]*mnuSelf+0.4517539514526256*m1rSelf[5]*mnuSelf-0.7857142857142857*m0rSelf[1]*uSelf[4]*mnuSelf-0.31943828249997*m0rSelf[2]*uSelf[3]*mnuSelf-0.5*m0rSelf[0]*uSelf[3]*mnuSelf+0.7071067811865475*m1rSelf[3]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherY-uOtherY*m0Other) and uCrossOtherY ... // 
  data->AEM_S(21,15) = 0.5*m0rOther[2]*uOther[5]*mnuOther+0.5*m0rOther[1]*uOther[4]*mnuOther+0.5*m0rOther[0]*uOther[3]*mnuOther-0.7071067811865475*m1rOther[3]*mnuOther; 
  data->AEM_S(21,16) = 0.4472135954999579*m0rOther[1]*uOther[5]*mnuOther+0.4472135954999579*m0rOther[2]*uOther[4]*mnuOther+0.5*m0rOther[0]*uOther[4]*mnuOther-0.7071067811865475*m1rOther[4]*mnuOther+0.5*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(21,17) = 0.31943828249997*m0rOther[2]*uOther[5]*mnuOther+0.5*m0rOther[0]*uOther[5]*mnuOther-0.7071067811865475*m1rOther[5]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[4]*mnuOther+0.5*m0rOther[2]*uOther[3]*mnuOther; 
  data->AEM_S(22,15) = 0.4472135954999579*m0rOther[1]*uOther[5]*mnuOther+0.4472135954999579*m0rOther[2]*uOther[4]*mnuOther+0.5*m0rOther[0]*uOther[4]*mnuOther-0.7071067811865475*m1rOther[4]*mnuOther+0.5*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(22,16) = 0.7857142857142857*m0rOther[2]*uOther[5]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[5]*mnuOther-0.6324555320336759*m1rOther[5]*mnuOther+0.9*m0rOther[1]*uOther[4]*mnuOther+0.4472135954999579*m0rOther[2]*uOther[3]*mnuOther+0.5*m0rOther[0]*uOther[3]*mnuOther-0.7071067811865475*m1rOther[3]*mnuOther; 
  data->AEM_S(22,17) = 0.7857142857142857*m0rOther[1]*uOther[5]*mnuOther+0.7857142857142857*m0rOther[2]*uOther[4]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[4]*mnuOther-0.6324555320336759*m1rOther[4]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(23,15) = 0.31943828249997*m0rOther[2]*uOther[5]*mnuOther+0.5*m0rOther[0]*uOther[5]*mnuOther-0.7071067811865475*m1rOther[5]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[4]*mnuOther+0.5*m0rOther[2]*uOther[3]*mnuOther; 
  data->AEM_S(23,16) = 0.7857142857142857*m0rOther[1]*uOther[5]*mnuOther+0.7857142857142857*m0rOther[2]*uOther[4]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[4]*mnuOther-0.6324555320336759*m1rOther[4]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[3]*mnuOther; 
  data->AEM_S(23,17) = 1.071428571428571*m0rOther[2]*uOther[5]*mnuOther+0.31943828249997*m0rOther[0]*uOther[5]*mnuOther-0.4517539514526256*m1rOther[5]*mnuOther+0.7857142857142857*m0rOther[1]*uOther[4]*mnuOther+0.31943828249997*m0rOther[2]*uOther[3]*mnuOther+0.5*m0rOther[0]*uOther[3]*mnuOther-0.7071067811865475*m1rOther[3]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfY-m0Self*m1OtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[3] = 0.7071067811865475*m0rOther[2]*m1rSelf[5]-0.7071067811865475*m0rSelf[2]*m1rOther[5]+0.7071067811865475*m0rOther[1]*m1rSelf[4]-0.7071067811865475*m0rSelf[1]*m1rOther[4]+0.7071067811865475*m0rOther[0]*m1rSelf[3]-0.7071067811865475*m0rSelf[0]*m1rOther[3]; 
  m1EffD[4] = 0.6324555320336759*m0rOther[1]*m1rSelf[5]-0.6324555320336759*m0rSelf[1]*m1rOther[5]+0.6324555320336759*m0rOther[2]*m1rSelf[4]+0.7071067811865475*m0rOther[0]*m1rSelf[4]-0.6324555320336759*m0rSelf[2]*m1rOther[4]-0.7071067811865475*m0rSelf[0]*m1rOther[4]+0.7071067811865475*m0rOther[1]*m1rSelf[3]-0.7071067811865475*m0rSelf[1]*m1rOther[3]; 
  m1EffD[5] = 0.4517539514526256*m0rOther[2]*m1rSelf[5]+0.7071067811865475*m0rOther[0]*m1rSelf[5]-0.4517539514526256*m0rSelf[2]*m1rOther[5]-0.7071067811865475*m0rSelf[0]*m1rOther[5]+0.6324555320336759*m0rOther[1]*m1rSelf[4]-0.6324555320336759*m0rSelf[1]*m1rOther[4]+0.7071067811865475*m0rOther[2]*m1rSelf[3]-0.7071067811865475*m0rSelf[2]*m1rOther[3]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[3],m1EffD[4],m1EffD[5]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+3,3,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 2 of momentum relaxation. 
  m1Relax[3] += (-2.0*m1EffD[3]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += (-2.0*m1EffD[4]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += (-2.0*m1EffD[5]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfZ ... // 
  data->AEM_S(18,6) = 0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,7) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(18,8) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,6) = 0.7071067811865475*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,7) = 0.6324555320336759*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
  data->AEM_S(19,8) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,6) = 0.7071067811865475*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,7) = 0.6324555320336759*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,8) = 0.4517539514526256*m0rSelf[2]*mnuSelf+0.7071067811865475*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(18,9) = -0.7071067811865475*cMrSelf[6]*mnuSelf; 
  data->AEM_S(18,10) = -0.7071067811865475*cMrSelf[7]*mnuSelf; 
  data->AEM_S(18,11) = -0.7071067811865475*cMrSelf[8]*mnuSelf; 
  data->AEM_S(19,9) = -0.7071067811865475*cMrSelf[7]*mnuSelf; 
  data->AEM_S(19,10) = (-0.6324555320336759*cMrSelf[8]*mnuSelf)-0.7071067811865475*cMrSelf[6]*mnuSelf; 
  data->AEM_S(19,11) = -0.6324555320336759*cMrSelf[7]*mnuSelf; 
  data->AEM_S(20,9) = -0.7071067811865475*cMrSelf[8]*mnuSelf; 
  data->AEM_S(20,10) = -0.6324555320336759*cMrSelf[7]*mnuSelf; 
  data->AEM_S(20,11) = (-0.4517539514526256*cMrSelf[8]*mnuSelf)-0.7071067811865475*cMrSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(18,18) = -0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(18,19) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(18,20) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(19,18) = -0.7071067811865475*m0rOther[1]*mnuOther; 
  data->AEM_S(19,19) = (-0.6324555320336759*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
  data->AEM_S(19,20) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(20,18) = -0.7071067811865475*m0rOther[2]*mnuOther; 
  data->AEM_S(20,19) = -0.6324555320336759*m0rOther[1]*mnuOther; 
  data->AEM_S(20,20) = (-0.4517539514526256*m0rOther[2]*mnuOther)-0.7071067811865475*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(18,21) = 0.7071067811865475*cMrOther[6]*mnuOther; 
  data->AEM_S(18,22) = 0.7071067811865475*cMrOther[7]*mnuOther; 
  data->AEM_S(18,23) = 0.7071067811865475*cMrOther[8]*mnuOther; 
  data->AEM_S(19,21) = 0.7071067811865475*cMrOther[7]*mnuOther; 
  data->AEM_S(19,22) = 0.6324555320336759*cMrOther[8]*mnuOther+0.7071067811865475*cMrOther[6]*mnuOther; 
  data->AEM_S(19,23) = 0.6324555320336759*cMrOther[7]*mnuOther; 
  data->AEM_S(20,21) = 0.7071067811865475*cMrOther[8]*mnuOther; 
  data->AEM_S(20,22) = 0.6324555320336759*cMrOther[7]*mnuOther; 
  data->AEM_S(20,23) = 0.4517539514526256*cMrOther[8]*mnuOther+0.7071067811865475*cMrOther[6]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfZ-uSelfZ*m0Self) and uCrossSelfZ ... // 
  data->AEM_S(21,6) = (-0.5*m0rSelf[2]*uSelf[8]*mnuSelf)-0.5*m0rSelf[1]*uSelf[7]*mnuSelf-0.5*m0rSelf[0]*uSelf[6]*mnuSelf+0.7071067811865475*m1rSelf[6]*mnuSelf; 
  data->AEM_S(21,7) = (-0.4472135954999579*m0rSelf[1]*uSelf[8]*mnuSelf)-0.4472135954999579*m0rSelf[2]*uSelf[7]*mnuSelf-0.5*m0rSelf[0]*uSelf[7]*mnuSelf+0.7071067811865475*m1rSelf[7]*mnuSelf-0.5*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(21,8) = (-0.31943828249997*m0rSelf[2]*uSelf[8]*mnuSelf)-0.5*m0rSelf[0]*uSelf[8]*mnuSelf+0.7071067811865475*m1rSelf[8]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[7]*mnuSelf-0.5*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(22,6) = (-0.4472135954999579*m0rSelf[1]*uSelf[8]*mnuSelf)-0.4472135954999579*m0rSelf[2]*uSelf[7]*mnuSelf-0.5*m0rSelf[0]*uSelf[7]*mnuSelf+0.7071067811865475*m1rSelf[7]*mnuSelf-0.5*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(22,7) = (-0.7857142857142857*m0rSelf[2]*uSelf[8]*mnuSelf)-0.4472135954999579*m0rSelf[0]*uSelf[8]*mnuSelf+0.6324555320336759*m1rSelf[8]*mnuSelf-0.9*m0rSelf[1]*uSelf[7]*mnuSelf-0.4472135954999579*m0rSelf[2]*uSelf[6]*mnuSelf-0.5*m0rSelf[0]*uSelf[6]*mnuSelf+0.7071067811865475*m1rSelf[6]*mnuSelf; 
  data->AEM_S(22,8) = (-0.7857142857142857*m0rSelf[1]*uSelf[8]*mnuSelf)-0.7857142857142857*m0rSelf[2]*uSelf[7]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[7]*mnuSelf+0.6324555320336759*m1rSelf[7]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(23,6) = (-0.31943828249997*m0rSelf[2]*uSelf[8]*mnuSelf)-0.5*m0rSelf[0]*uSelf[8]*mnuSelf+0.7071067811865475*m1rSelf[8]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[7]*mnuSelf-0.5*m0rSelf[2]*uSelf[6]*mnuSelf; 
  data->AEM_S(23,7) = (-0.7857142857142857*m0rSelf[1]*uSelf[8]*mnuSelf)-0.7857142857142857*m0rSelf[2]*uSelf[7]*mnuSelf-0.4472135954999579*m0rSelf[0]*uSelf[7]*mnuSelf+0.6324555320336759*m1rSelf[7]*mnuSelf-0.4472135954999579*m0rSelf[1]*uSelf[6]*mnuSelf; 
  data->AEM_S(23,8) = (-1.071428571428571*m0rSelf[2]*uSelf[8]*mnuSelf)-0.31943828249997*m0rSelf[0]*uSelf[8]*mnuSelf+0.4517539514526256*m1rSelf[8]*mnuSelf-0.7857142857142857*m0rSelf[1]*uSelf[7]*mnuSelf-0.31943828249997*m0rSelf[2]*uSelf[6]*mnuSelf-0.5*m0rSelf[0]*uSelf[6]*mnuSelf+0.7071067811865475*m1rSelf[6]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherZ-uOtherZ*m0Other) and uCrossOtherZ ... // 
  data->AEM_S(21,18) = 0.5*m0rOther[2]*uOther[8]*mnuOther+0.5*m0rOther[1]*uOther[7]*mnuOther+0.5*m0rOther[0]*uOther[6]*mnuOther-0.7071067811865475*m1rOther[6]*mnuOther; 
  data->AEM_S(21,19) = 0.4472135954999579*m0rOther[1]*uOther[8]*mnuOther+0.4472135954999579*m0rOther[2]*uOther[7]*mnuOther+0.5*m0rOther[0]*uOther[7]*mnuOther-0.7071067811865475*m1rOther[7]*mnuOther+0.5*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(21,20) = 0.31943828249997*m0rOther[2]*uOther[8]*mnuOther+0.5*m0rOther[0]*uOther[8]*mnuOther-0.7071067811865475*m1rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[7]*mnuOther+0.5*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(22,18) = 0.4472135954999579*m0rOther[1]*uOther[8]*mnuOther+0.4472135954999579*m0rOther[2]*uOther[7]*mnuOther+0.5*m0rOther[0]*uOther[7]*mnuOther-0.7071067811865475*m1rOther[7]*mnuOther+0.5*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(22,19) = 0.7857142857142857*m0rOther[2]*uOther[8]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[8]*mnuOther-0.6324555320336759*m1rOther[8]*mnuOther+0.9*m0rOther[1]*uOther[7]*mnuOther+0.4472135954999579*m0rOther[2]*uOther[6]*mnuOther+0.5*m0rOther[0]*uOther[6]*mnuOther-0.7071067811865475*m1rOther[6]*mnuOther; 
  data->AEM_S(22,20) = 0.7857142857142857*m0rOther[1]*uOther[8]*mnuOther+0.7857142857142857*m0rOther[2]*uOther[7]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[7]*mnuOther-0.6324555320336759*m1rOther[7]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(23,18) = 0.31943828249997*m0rOther[2]*uOther[8]*mnuOther+0.5*m0rOther[0]*uOther[8]*mnuOther-0.7071067811865475*m1rOther[8]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[7]*mnuOther+0.5*m0rOther[2]*uOther[6]*mnuOther; 
  data->AEM_S(23,19) = 0.7857142857142857*m0rOther[1]*uOther[8]*mnuOther+0.7857142857142857*m0rOther[2]*uOther[7]*mnuOther+0.4472135954999579*m0rOther[0]*uOther[7]*mnuOther-0.6324555320336759*m1rOther[7]*mnuOther+0.4472135954999579*m0rOther[1]*uOther[6]*mnuOther; 
  data->AEM_S(23,20) = 1.071428571428571*m0rOther[2]*uOther[8]*mnuOther+0.31943828249997*m0rOther[0]*uOther[8]*mnuOther-0.4517539514526256*m1rOther[8]*mnuOther+0.7857142857142857*m0rOther[1]*uOther[7]*mnuOther+0.31943828249997*m0rOther[2]*uOther[6]*mnuOther+0.5*m0rOther[0]*uOther[6]*mnuOther-0.7071067811865475*m1rOther[6]*mnuOther; 
 
  // ... Divide (m0Other*m1SelfZ-m0Self*m1OtherZ)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute m0Other*m1Self-m0Self*m1Other. 
  m1EffD[6] = 0.7071067811865475*m0rOther[2]*m1rSelf[8]-0.7071067811865475*m0rSelf[2]*m1rOther[8]+0.7071067811865475*m0rOther[1]*m1rSelf[7]-0.7071067811865475*m0rSelf[1]*m1rOther[7]+0.7071067811865475*m0rOther[0]*m1rSelf[6]-0.7071067811865475*m0rSelf[0]*m1rOther[6]; 
  m1EffD[7] = 0.6324555320336759*m0rOther[1]*m1rSelf[8]-0.6324555320336759*m0rSelf[1]*m1rOther[8]+0.6324555320336759*m0rOther[2]*m1rSelf[7]+0.7071067811865475*m0rOther[0]*m1rSelf[7]-0.6324555320336759*m0rSelf[2]*m1rOther[7]-0.7071067811865475*m0rSelf[0]*m1rOther[7]+0.7071067811865475*m0rOther[1]*m1rSelf[6]-0.7071067811865475*m0rSelf[1]*m1rOther[6]; 
  m1EffD[8] = 0.4517539514526256*m0rOther[2]*m1rSelf[8]+0.7071067811865475*m0rOther[0]*m1rSelf[8]-0.4517539514526256*m0rSelf[2]*m1rOther[8]-0.7071067811865475*m0rSelf[0]*m1rOther[8]+0.6324555320336759*m0rOther[1]*m1rSelf[7]-0.6324555320336759*m0rSelf[1]*m1rOther[7]+0.7071067811865475*m0rOther[2]*m1rSelf[6]-0.7071067811865475*m0rSelf[2]*m1rOther[6]; 
  // Fill BEV. 
  dataDiv->BEV_S << m1EffD[6],m1EffD[7],m1EffD[8]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(m1EffD+6,3,1) = dataDiv->u_S; 
 
  // ... Contribution to RHS vector from component 3 of momentum relaxation. 
  m1Relax[6] += (-2.0*m1EffD[6]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += (-2.0*m1EffD[7]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
  m1Relax[8] += (-2.0*m1EffD[8]*betaGreenep1*mnuOther*mnuSelf)+m1rSelf[8]*mnuSelf-1.0*m1rOther[8]*mnuOther; 
 
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
    ucMSelf[0] += 0.7071067811865475*cMrSelf[a0+2]*uSelf[a0+2]+0.7071067811865475*cMrSelf[a0+1]*uSelf[a0+1]+0.7071067811865475*cMrSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.6324555320336759*cMrSelf[a0+1]*uSelf[a0+2]+0.6324555320336759*uSelf[a0+1]*cMrSelf[a0+2]+0.7071067811865475*cMrSelf[a0]*uSelf[a0+1]+0.7071067811865475*uSelf[a0]*cMrSelf[a0+1]; 
    ucMSelf[2] += 0.4517539514526256*cMrSelf[a0+2]*uSelf[a0+2]+0.7071067811865475*cMrSelf[a0]*uSelf[a0+2]+0.7071067811865475*uSelf[a0]*cMrSelf[a0+2]+0.6324555320336759*cMrSelf[a0+1]*uSelf[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(21,9) = 0.7071067811865475*ucMSelf[0]*mnuSelf+2.121320343559642*m0rSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
  data->AEM_S(21,10) = 0.7071067811865475*ucMSelf[1]*mnuSelf+2.121320343559642*m0rSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(21,11) = 0.7071067811865475*ucMSelf[2]*mnuSelf+2.121320343559642*m0rSelf[2]*mnuSelf-0.7071067811865475*cErSelf[2]*mnuSelf; 
  data->AEM_S(22,9) = 0.7071067811865475*ucMSelf[1]*mnuSelf+2.121320343559642*m0rSelf[1]*mnuSelf-0.7071067811865475*cErSelf[1]*mnuSelf; 
  data->AEM_S(22,10) = 0.6324555320336759*ucMSelf[2]*mnuSelf+1.897366596101028*m0rSelf[2]*mnuSelf-0.6324555320336759*cErSelf[2]*mnuSelf+0.7071067811865475*ucMSelf[0]*mnuSelf+2.121320343559642*m0rSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
  data->AEM_S(22,11) = 0.6324555320336759*ucMSelf[1]*mnuSelf+1.897366596101028*m0rSelf[1]*mnuSelf-0.6324555320336759*cErSelf[1]*mnuSelf; 
  data->AEM_S(23,9) = 0.7071067811865475*ucMSelf[2]*mnuSelf+2.121320343559642*m0rSelf[2]*mnuSelf-0.7071067811865475*cErSelf[2]*mnuSelf; 
  data->AEM_S(23,10) = 0.6324555320336759*ucMSelf[1]*mnuSelf+1.897366596101028*m0rSelf[1]*mnuSelf-0.6324555320336759*cErSelf[1]*mnuSelf; 
  data->AEM_S(23,11) = 0.4517539514526256*ucMSelf[2]*mnuSelf+1.355261854357877*m0rSelf[2]*mnuSelf-0.4517539514526256*cErSelf[2]*mnuSelf+0.7071067811865475*ucMSelf[0]*mnuSelf+2.121320343559642*m0rSelf[0]*mnuSelf-0.7071067811865475*cErSelf[0]*mnuSelf; 
 
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
    ucMOther[0] += 0.7071067811865475*cMrOther[a0+2]*uOther[a0+2]+0.7071067811865475*cMrOther[a0+1]*uOther[a0+1]+0.7071067811865475*cMrOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.6324555320336759*cMrOther[a0+1]*uOther[a0+2]+0.6324555320336759*uOther[a0+1]*cMrOther[a0+2]+0.7071067811865475*cMrOther[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*cMrOther[a0+1]; 
    ucMOther[2] += 0.4517539514526256*cMrOther[a0+2]*uOther[a0+2]+0.7071067811865475*cMrOther[a0]*uOther[a0+2]+0.7071067811865475*uOther[a0]*cMrOther[a0+2]+0.6324555320336759*cMrOther[a0+1]*uOther[a0+1]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(21,21) = (-0.7071067811865475*ucMOther[0]*mnuOther)-2.121320343559642*m0rOther[0]*mnuOther+0.7071067811865475*cErOther[0]*mnuOther; 
  data->AEM_S(21,22) = (-0.7071067811865475*ucMOther[1]*mnuOther)-2.121320343559642*m0rOther[1]*mnuOther+0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(21,23) = (-0.7071067811865475*ucMOther[2]*mnuOther)-2.121320343559642*m0rOther[2]*mnuOther+0.7071067811865475*cErOther[2]*mnuOther; 
  data->AEM_S(22,21) = (-0.7071067811865475*ucMOther[1]*mnuOther)-2.121320343559642*m0rOther[1]*mnuOther+0.7071067811865475*cErOther[1]*mnuOther; 
  data->AEM_S(22,22) = (-0.6324555320336759*ucMOther[2]*mnuOther)-1.897366596101028*m0rOther[2]*mnuOther+0.6324555320336759*cErOther[2]*mnuOther-0.7071067811865475*ucMOther[0]*mnuOther-2.121320343559642*m0rOther[0]*mnuOther+0.7071067811865475*cErOther[0]*mnuOther; 
  data->AEM_S(22,23) = (-0.6324555320336759*ucMOther[1]*mnuOther)-1.897366596101028*m0rOther[1]*mnuOther+0.6324555320336759*cErOther[1]*mnuOther; 
  data->AEM_S(23,21) = (-0.7071067811865475*ucMOther[2]*mnuOther)-2.121320343559642*m0rOther[2]*mnuOther+0.7071067811865475*cErOther[2]*mnuOther; 
  data->AEM_S(23,22) = (-0.6324555320336759*ucMOther[1]*mnuOther)-1.897366596101028*m0rOther[1]*mnuOther+0.6324555320336759*cErOther[1]*mnuOther; 
  data->AEM_S(23,23) = (-0.4517539514526256*ucMOther[2]*mnuOther)-1.355261854357877*m0rOther[2]*mnuOther+0.4517539514526256*cErOther[2]*mnuOther-0.7071067811865475*ucMOther[0]*mnuOther-2.121320343559642*m0rOther[0]*mnuOther+0.7071067811865475*cErOther[0]*mnuOther; 
 
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
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.7071067811865475*m1rOther[a0+2]*uOther[a0+2]+0.7071067811865475*m1rOther[a0+1]*uOther[a0+1]+0.7071067811865475*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.6324555320336759*m1rOther[a0+1]*uOther[a0+2]+0.6324555320336759*uOther[a0+1]*m1rOther[a0+2]+0.7071067811865475*m1rOther[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.4517539514526256*m1rOther[a0+2]*uOther[a0+2]+0.7071067811865475*m1rOther[a0]*uOther[a0+2]+0.7071067811865475*uOther[a0]*m1rOther[a0+2]+0.6324555320336759*m1rOther[a0+1]*uOther[a0+1]; 
  } 
 
  double relKinE[3]; 
  // zero out array with dot product of uSelf-uOther and m1EffD. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.7071067811865475*m1EffD[a0+2]*uSelf[a0+2]-0.7071067811865475*m1EffD[a0+2]*uOther[a0+2]+0.7071067811865475*m1EffD[a0+1]*uSelf[a0+1]-0.7071067811865475*m1EffD[a0+1]*uOther[a0+1]+0.7071067811865475*m1EffD[a0]*uSelf[a0]-0.7071067811865475*m1EffD[a0]*uOther[a0]; 
    relKinE[1] += 0.6324555320336759*m1EffD[a0+1]*uSelf[a0+2]-0.6324555320336759*m1EffD[a0+1]*uOther[a0+2]+0.6324555320336759*uSelf[a0+1]*m1EffD[a0+2]-0.6324555320336759*uOther[a0+1]*m1EffD[a0+2]+0.7071067811865475*m1EffD[a0]*uSelf[a0+1]-0.7071067811865475*m1EffD[a0]*uOther[a0+1]+0.7071067811865475*uSelf[a0]*m1EffD[a0+1]-0.7071067811865475*uOther[a0]*m1EffD[a0+1]; 
    relKinE[2] += 0.4517539514526256*m1EffD[a0+2]*uSelf[a0+2]+0.7071067811865475*m1EffD[a0]*uSelf[a0+2]-0.4517539514526256*m1EffD[a0+2]*uOther[a0+2]-0.7071067811865475*m1EffD[a0]*uOther[a0+2]+0.7071067811865475*uSelf[a0]*m1EffD[a0+2]-0.7071067811865475*uOther[a0]*m1EffD[a0+2]+0.6324555320336759*m1EffD[a0+1]*uSelf[a0+1]-0.6324555320336759*m1EffD[a0+1]*uOther[a0+1]; 
  } 
 
  // Divide m0Other*(m2Self-kinESelf) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Other and m2Self-uSelf.m1Self. 
  double m0OtherThESelf[3]; 
  m0OtherThESelf[0] = 0.7071067811865475*m0rOther[2]*m2rSelf[2]-0.7071067811865475*kinESelf[2]*m0rOther[2]+0.7071067811865475*m0rOther[1]*m2rSelf[1]-0.7071067811865475*kinESelf[1]*m0rOther[1]+0.7071067811865475*m0rOther[0]*m2rSelf[0]-0.7071067811865475*kinESelf[0]*m0rOther[0]; 
  m0OtherThESelf[1] = 0.6324555320336759*m0rOther[1]*m2rSelf[2]+0.6324555320336759*m2rSelf[1]*m0rOther[2]-0.6324555320336759*kinESelf[1]*m0rOther[2]-0.6324555320336759*m0rOther[1]*kinESelf[2]+0.7071067811865475*m0rOther[0]*m2rSelf[1]+0.7071067811865475*m2rSelf[0]*m0rOther[1]-0.7071067811865475*kinESelf[0]*m0rOther[1]-0.7071067811865475*m0rOther[0]*kinESelf[1]; 
  m0OtherThESelf[2] = 0.4517539514526256*m0rOther[2]*m2rSelf[2]+0.7071067811865475*m0rOther[0]*m2rSelf[2]-0.4517539514526256*kinESelf[2]*m0rOther[2]+0.7071067811865475*m2rSelf[0]*m0rOther[2]-0.7071067811865475*kinESelf[0]*m0rOther[2]-0.7071067811865475*m0rOther[0]*kinESelf[2]+0.6324555320336759*m0rOther[1]*m2rSelf[1]-0.6324555320336759*kinESelf[1]*m0rOther[1]; 
  dataDiv->BEV_S << m0OtherThESelf[0],m0OtherThESelf[1],m0OtherThESelf[2]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthSelf[3]; 
  Eigen::Map<VectorXd>(effEthSelf,3,1) = dataDiv->u_S; 
 
  // Divide m0Self*(m2Other-kinEOther) by mnuSelf*m0Self+mnuOther*m0Other. 
  // Product of m0Self and m2Other-uOther.m1Other. 
  double m0SelfThEOther[3]; 
  m0SelfThEOther[0] = 0.7071067811865475*m0rSelf[2]*m2rOther[2]-0.7071067811865475*kinEOther[2]*m0rSelf[2]+0.7071067811865475*m0rSelf[1]*m2rOther[1]-0.7071067811865475*kinEOther[1]*m0rSelf[1]+0.7071067811865475*m0rSelf[0]*m2rOther[0]-0.7071067811865475*kinEOther[0]*m0rSelf[0]; 
  m0SelfThEOther[1] = 0.6324555320336759*m0rSelf[1]*m2rOther[2]+0.6324555320336759*m2rOther[1]*m0rSelf[2]-0.6324555320336759*kinEOther[1]*m0rSelf[2]-0.6324555320336759*m0rSelf[1]*kinEOther[2]+0.7071067811865475*m0rSelf[0]*m2rOther[1]+0.7071067811865475*m2rOther[0]*m0rSelf[1]-0.7071067811865475*kinEOther[0]*m0rSelf[1]-0.7071067811865475*m0rSelf[0]*kinEOther[1]; 
  m0SelfThEOther[2] = 0.4517539514526256*m0rSelf[2]*m2rOther[2]+0.7071067811865475*m0rSelf[0]*m2rOther[2]-0.4517539514526256*kinEOther[2]*m0rSelf[2]+0.7071067811865475*m2rOther[0]*m0rSelf[2]-0.7071067811865475*kinEOther[0]*m0rSelf[2]-0.7071067811865475*m0rSelf[0]*kinEOther[2]+0.6324555320336759*m0rSelf[1]*m2rOther[1]-0.6324555320336759*kinEOther[1]*m0rSelf[1]; 
  dataDiv->BEV_S << m0SelfThEOther[0],m0SelfThEOther[1],m0SelfThEOther[2]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double effEthOther[3]; 
  Eigen::Map<VectorXd>(effEthOther,3,1) = dataDiv->u_S; 
 
  double m2Relax[3]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(1.0*relKinE[0]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[0]*mSelf)/(mSelf+mOther)+(1.0*relKinE[0]*mOther)/(mSelf+mOther)+(2.0*effEthOther[0]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2rOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(1.0*relKinE[1]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[1]*mSelf)/(mSelf+mOther)+(1.0*relKinE[1]*mOther)/(mSelf+mOther)+(2.0*effEthOther[1]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2rOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(1.0*relKinE[2]*mSelf)/(mSelf+mOther))-(2.0*effEthSelf[2]*mSelf)/(mSelf+mOther)+(1.0*relKinE[2]*mOther)/(mSelf+mOther)+(2.0*effEthOther[2]*mOther)/(mSelf+mOther))*mnuOther*mnuSelf+(m2rSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2rOther[2])*mnuOther; 
 
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
 
