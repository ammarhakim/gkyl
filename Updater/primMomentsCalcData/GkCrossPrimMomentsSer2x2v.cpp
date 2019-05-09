#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkCrossPrimMoments2x2vSer_P1(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[4]; 
  double m2rSelf[4]; 
  double m0SrSelf[4]; 
  double m1SrSelf[4]; 
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
  double m1rOther[4]; 
  double m2rOther[4]; 
  double m0SrOther[4]; 
  double m1SrOther[4]; 
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
    m2SrOther[0] = m2SOther[0]; 
    m2SrOther[1] = m2SOther[1]; 
    m2SrOther[2] = m2SOther[2]; 
    m2SrOther[3] = m2SOther[3]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(16,16); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double deltaSelf = sqrt(mnuOther/mnuSelf); 
  double mnuM1sum[4]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<4; vd++) 
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
  data->AEM_S(0,4) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,5) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,6) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,7) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,4) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,5) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,6) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,7) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,4) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,5) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,6) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,7) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,4) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,5) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,6) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,7) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,8) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,9) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,10) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,11) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,8) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,9) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,10) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,11) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,8) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,9) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,10) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,11) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,8) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,9) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,10) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,11) = 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,12) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,13) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,14) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,15) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,12) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,13) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,14) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,15) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,12) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,13) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,14) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,15) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,12) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,13) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,14) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,15) = -0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(4,0) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(4,1) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(4,2) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(4,3) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(5,0) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(5,1) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(5,2) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(5,3) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(6,0) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(6,1) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(6,2) = 0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(6,3) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(7,0) = 0.5*m1SrSelf[3]*mnuSelf; 
  data->AEM_S(7,1) = 0.5*m1SrSelf[2]*mnuSelf; 
  data->AEM_S(7,2) = 0.5*m1SrSelf[1]*mnuSelf; 
  data->AEM_S(7,3) = 0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(4,8) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(4,9) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(4,10) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(4,11) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(5,8) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(5,9) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(5,10) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(5,11) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(6,8) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(6,9) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(6,10) = 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(6,11) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(7,8) = 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(7,9) = 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(7,10) = 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(7,11) = 0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(4,4) = m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(4,5) = m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(4,6) = m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(4,7) = m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(5,4) = m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(5,5) = m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(5,6) = m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(5,7) = m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(6,4) = m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(6,5) = m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(6,6) = m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(6,7) = m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(7,4) = m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(7,5) = m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(7,6) = m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(7,7) = m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(4,12) = m0rOther[0]*mnuOther+0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(4,13) = m0rOther[1]*mnuOther+0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(4,14) = m0rOther[2]*mnuOther+0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(4,15) = m0rOther[3]*mnuOther+0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(5,12) = m0rOther[1]*mnuOther+0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(5,13) = m0rOther[0]*mnuOther+0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(5,14) = m0rOther[3]*mnuOther+0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(5,15) = m0rOther[2]*mnuOther+0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(6,12) = m0rOther[2]*mnuOther+0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(6,13) = m0rOther[3]*mnuOther+0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(6,14) = m0rOther[0]*mnuOther+0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(6,15) = m0rOther[1]*mnuOther+0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(7,12) = m0rOther[3]*mnuOther+0.5*m0SrOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(7,13) = m0rOther[2]*mnuOther+0.5*m0SrOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(7,14) = m0rOther[1]*mnuOther+0.5*m0SrOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(7,15) = m0rOther[0]*mnuOther+0.5*m0SrOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
  double mnuM2sum[4]; 
  // ... Contribution to RHS vector from mnuM2Self+mnuM2Other. 
  mnuM2sum[0] = m2SrSelf[0]*mnuSelf+m2SrOther[0]*mnuOther; 
  mnuM2sum[1] = m2SrSelf[1]*mnuSelf+m2SrOther[1]*mnuOther; 
  mnuM2sum[2] = m2SrSelf[2]*mnuSelf+m2SrOther[2]*mnuOther; 
  mnuM2sum[3] = m2SrSelf[3]*mnuSelf+m2SrOther[3]*mnuOther; 
 
  double m1Relax[4]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(8,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(8,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(8,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(8,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(9,1) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(9,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(9,3) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(10,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(10,2) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(10,3) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(11,1) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(11,2) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(11,3) = 0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(8,4) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(8,5) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(8,6) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(8,7) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(9,4) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(9,5) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(9,6) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(9,7) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(10,4) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(10,5) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(10,6) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(10,7) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(11,4) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(11,5) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(11,6) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(11,7) = -0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(8,8) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,9) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(8,10) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,11) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,8) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(9,9) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,10) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,11) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,8) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,9) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(10,10) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(10,11) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,8) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(11,9) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(11,10) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,11) = -0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(8,12) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(8,13) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(8,14) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(8,15) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(9,12) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(9,13) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(9,14) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(9,15) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(10,12) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(10,13) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(10,14) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(10,15) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(11,12) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(11,13) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(11,14) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(11,15) = 0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(12,0) = (-0.25*m0rSelf[3]*uSelf[3]*mnuSelf)-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(12,1) = (-0.25*m0rSelf[2]*uSelf[3]*mnuSelf)-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = (-0.25*m0rSelf[1]*uSelf[3]*mnuSelf)-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(12,3) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,0) = (-0.25*m0rSelf[2]*uSelf[3]*mnuSelf)-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(13,1) = (-0.45*m0rSelf[3]*uSelf[3]*mnuSelf)-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(13,2) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(13,3) = (-0.45*m0rSelf[1]*uSelf[3]*mnuSelf)-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,0) = (-0.25*m0rSelf[1]*uSelf[3]*mnuSelf)-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,1) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(14,2) = (-0.45*m0rSelf[3]*uSelf[3]*mnuSelf)-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
  data->AEM_S(14,3) = (-0.45*m0rSelf[2]*uSelf[3]*mnuSelf)-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,0) = (-0.25*m0rSelf[0]*uSelf[3]*mnuSelf)+0.5*m1SrSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,1) = (-0.45*m0rSelf[1]*uSelf[3]*mnuSelf)-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1SrSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(15,2) = (-0.45*m0rSelf[2]*uSelf[3]*mnuSelf)-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1SrSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(15,3) = (-0.81*m0rSelf[3]*uSelf[3]*mnuSelf)-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1SrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(12,8) = 0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(12,9) = 0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(12,10) = 0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(12,11) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(13,8) = 0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(13,9) = 0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(13,10) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(13,11) = 0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(14,8) = 0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(14,9) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(14,10) = 0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(14,11) = 0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(15,8) = 0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1SrOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(15,9) = 0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1SrOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(15,10) = 0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1SrOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(15,11) = 0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1SrOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += betaGreenep1*(m1rOther[0]*deltaSelf-1.0*m1rSelf[0]*deltaSelf)*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += betaGreenep1*(m1rOther[1]*deltaSelf-1.0*m1rSelf[1]*deltaSelf)*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += betaGreenep1*(m1rOther[2]*deltaSelf-1.0*m1rSelf[2]*deltaSelf)*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += betaGreenep1*(m1rOther[3]*deltaSelf-1.0*m1rSelf[3]*deltaSelf)*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
 
  double ucMSelf[4]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.5*cMSelf[a0+3]*uSelf[a0+3]+0.5*cMSelf[a0+2]*uSelf[a0+2]+0.5*cMSelf[a0+1]*uSelf[a0+1]+0.5*cMSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.5*cMSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*cMSelf[a0+1]; 
    ucMSelf[2] += 0.5*cMSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*cMSelf[a0+3]+0.5*cMSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*cMSelf[a0+2]; 
    ucMSelf[3] += 0.5*cMSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*cMSelf[a0+3]+0.5*cMSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*cMSelf[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(12,4) = 0.5*ucMSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(12,5) = 0.5*ucMSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(12,6) = 0.5*ucMSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(12,7) = 0.5*ucMSelf[3]*mnuSelf+m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(13,4) = 0.5*ucMSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(13,5) = 0.5*ucMSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(13,6) = 0.5*ucMSelf[3]*mnuSelf+m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(13,7) = 0.5*ucMSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(14,4) = 0.5*ucMSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(14,5) = 0.5*ucMSelf[3]*mnuSelf+m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(14,6) = 0.5*ucMSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(14,7) = 0.5*ucMSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(15,4) = 0.5*ucMSelf[3]*mnuSelf+m0rSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(15,5) = 0.5*ucMSelf[2]*mnuSelf+m0rSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(15,6) = 0.5*ucMSelf[1]*mnuSelf+m0rSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(15,7) = 0.5*ucMSelf[0]*mnuSelf+m0rSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  double ucMOther[4]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMOther[0] += 0.5*cMOther[a0+3]*uOther[a0+3]+0.5*cMOther[a0+2]*uOther[a0+2]+0.5*cMOther[a0+1]*uOther[a0+1]+0.5*cMOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.5*cMOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+1]+0.5*uOther[a0]*cMOther[a0+1]; 
    ucMOther[2] += 0.5*cMOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*cMOther[a0+3]+0.5*cMOther[a0]*uOther[a0+2]+0.5*uOther[a0]*cMOther[a0+2]; 
    ucMOther[3] += 0.5*cMOther[a0]*uOther[a0+3]+0.5*uOther[a0]*cMOther[a0+3]+0.5*cMOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*cMOther[a0+2]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(12,12) = (-0.5*ucMOther[0]*mnuOther)-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(12,13) = (-0.5*ucMOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(12,14) = (-0.5*ucMOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(12,15) = (-0.5*ucMOther[3]*mnuOther)-1.0*m0rOther[3]*mnuOther-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(13,12) = (-0.5*ucMOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(13,13) = (-0.5*ucMOther[0]*mnuOther)-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(13,14) = (-0.5*ucMOther[3]*mnuOther)-1.0*m0rOther[3]*mnuOther-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(13,15) = (-0.5*ucMOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(14,12) = (-0.5*ucMOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(14,13) = (-0.5*ucMOther[3]*mnuOther)-1.0*m0rOther[3]*mnuOther-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(14,14) = (-0.5*ucMOther[0]*mnuOther)-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(14,15) = (-0.5*ucMOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(15,12) = (-0.5*ucMOther[3]*mnuOther)-1.0*m0rOther[3]*mnuOther-0.5*m0SrOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(15,13) = (-0.5*ucMOther[2]*mnuOther)-1.0*m0rOther[2]*mnuOther-0.5*m0SrOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(15,14) = (-0.5*ucMOther[1]*mnuOther)-1.0*m0rOther[1]*mnuOther-0.5*m0SrOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(15,15) = (-0.5*ucMOther[0]*mnuOther)-1.0*m0rOther[0]*mnuOther-0.5*m0SrOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[4]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
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
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    kinEOther[0] += 0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0]; 
    kinEOther[1] += 0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    kinEOther[2] += 0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    kinEOther[3] += 0.5*m1rOther[a0]*uOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
  } 
 
  double relKinE[4]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.5*m1rSelf[a0+3]*uSelf[a0+3]-0.5*m1rOther[a0+3]*uSelf[a0+3]-0.5*m1rSelf[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]-0.5*m1rOther[a0+2]*uSelf[a0+2]-0.5*m1rSelf[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]-0.5*m1rOther[a0+1]*uSelf[a0+1]-0.5*m1rSelf[a0+1]*uOther[a0+1]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]-0.5*m1rOther[a0]*uSelf[a0]-0.5*m1rSelf[a0]*uOther[a0]+0.5*m1rOther[a0]*uOther[a0]; 
    relKinE[1] += 0.5*m1rSelf[a0+2]*uSelf[a0+3]-0.5*m1rOther[a0+2]*uSelf[a0+3]-0.5*m1rSelf[a0+2]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]-0.5*uOther[a0+2]*m1rSelf[a0+3]-0.5*uSelf[a0+2]*m1rOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]-0.5*m1rOther[a0]*uSelf[a0+1]-0.5*m1rSelf[a0]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]-0.5*uOther[a0]*m1rSelf[a0+1]-0.5*uSelf[a0]*m1rOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    relKinE[2] += 0.5*m1rSelf[a0+1]*uSelf[a0+3]-0.5*m1rOther[a0+1]*uSelf[a0+3]-0.5*m1rSelf[a0+1]*uOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]-0.5*uOther[a0+1]*m1rSelf[a0+3]-0.5*uSelf[a0+1]*m1rOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]-0.5*m1rOther[a0]*uSelf[a0+2]-0.5*m1rSelf[a0]*uOther[a0+2]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]-0.5*uOther[a0]*m1rSelf[a0+2]-0.5*uSelf[a0]*m1rOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    relKinE[3] += 0.5*m1rSelf[a0]*uSelf[a0+3]-0.5*m1rOther[a0]*uSelf[a0+3]-0.5*m1rSelf[a0]*uOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]-0.5*uOther[a0]*m1rSelf[a0+3]-0.5*uSelf[a0]*m1rOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]-0.5*m1rOther[a0+1]*uSelf[a0+2]-0.5*m1rSelf[a0+1]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]-0.5*uOther[a0+1]*m1rSelf[a0+2]-0.5*uSelf[a0+1]*m1rOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
  } 
 
  double m2Relax[4]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(0.5*relKinE[0]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[0]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[0]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[0]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[0]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[0]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2SrOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(0.5*relKinE[1]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[1]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[1]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[1]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[1]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[1]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2SrOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(0.5*relKinE[2]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[2]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[2]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[2]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[2]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[2]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2SrOther[2])*mnuOther; 
  m2Relax[3] = betaGreenep1*((-(0.5*relKinE[3]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[3]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[3]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[3]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[3]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[3]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2SrSelf[3]-1.0*kinESelf[3])*mnuSelf+(kinEOther[3]-1.0*m2SrOther[3])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,4,1) = data->u_S.segment<4>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,4,1) = data->u_S.segment<4>(4); 
 
  Eigen::Map<VectorXd>(uCrossOther,4,1) = data->u_S.segment<4>(8); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,4,1) = data->u_S.segment<4>(12); 
 
} 
 
void GkCrossPrimMoments2x2vSer_P2(binOpData_t *data, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  double m1rSelf[8]; 
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
  double m1rOther[8]; 
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
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = m2Other[1]; 
    m2rOther[2] = m2Other[2]; 
    m2rOther[3] = m2Other[3]; 
    m2rOther[4] = m2Other[4]; 
    m2rOther[5] = m2Other[5]; 
    m2rOther[6] = m2Other[6]; 
    m2rOther[7] = m2Other[7]; 
  } 
 
  // Declare Eigen matrix and vectors for weak division. 
  data->AEM_S = Eigen::MatrixXd::Zero(32,32); 
 
  double mnuSelf   = mSelf*nuSelf; 
  double mnuOther  = mOther*nuOther; 
  double deltaSelf = sqrt(mnuOther/mnuSelf); 
  double mnuM1sum[8]; 
  // zero out array with sum of m*nu*m1. 
  for (unsigned short int vd=0; vd<8; vd++) 
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
  data->AEM_S(0,8) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(0,9) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(0,10) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(0,11) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(0,12) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(0,13) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(0,14) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(0,15) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(1,8) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,9) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(1,10) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,11) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(1,12) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(1,13) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(1,14) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(1,15) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(2,8) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,9) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(2,10) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(2,11) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(2,12) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(2,13) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(2,14) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(2,15) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,8) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,9) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(3,10) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,11) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(3,12) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,13) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(3,14) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(3,15) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(4,8) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(4,9) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(4,10) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(4,11) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(4,12) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(4,14) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(4,15) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,8) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(5,9) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(5,10) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(5,11) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(5,13) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(5,14) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(5,15) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,8) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,9) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(6,10) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(6,11) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(6,12) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(6,13) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(6,14) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(6,15) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,8) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,9) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(7,10) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,11) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(7,12) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(7,13) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(7,14) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(7,15) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,16) = 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,17) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,18) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,19) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(0,20) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(0,21) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(0,22) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(0,23) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(1,16) = 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,17) = 0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,18) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,19) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(1,20) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(1,21) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(1,22) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(1,23) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(2,16) = 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,17) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,18) = 0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,19) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(2,20) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(2,21) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(2,22) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(2,23) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(3,16) = 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,17) = 0.447213595499958*m0rOther[6]*mnuOther+0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,18) = 0.447213595499958*m0rOther[7]*mnuOther+0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,19) = 0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(3,20) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,21) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(3,22) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(3,23) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(4,16) = 0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(4,17) = 0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(4,18) = 0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(4,19) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(4,20) = 0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(4,22) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(4,23) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(5,16) = 0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(5,17) = 0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(5,18) = 0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(5,19) = 0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(5,21) = 0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,22) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(5,23) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(6,16) = 0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(6,17) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(6,18) = 0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(6,19) = 0.4*m0rOther[7]*mnuOther+0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(6,20) = 0.31943828249997*m0rOther[6]*mnuOther+0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(6,21) = 0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(6,22) = 0.4472135954999579*m0rOther[5]*mnuOther+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,23) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(7,16) = 0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(7,17) = 0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(7,18) = 0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(7,19) = 0.4*m0rOther[6]*mnuOther+0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(7,20) = 0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(7,21) = 0.31943828249997*m0rOther[7]*mnuOther+0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(7,22) = 0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(7,23) = 0.31943828249997*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,24) = -0.5*cMOther[0]*mnuOther; 
  data->AEM_S(0,25) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(0,26) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(0,27) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(0,28) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(0,29) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(0,30) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(0,31) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(1,24) = -0.5*cMOther[1]*mnuOther; 
  data->AEM_S(1,25) = (-0.4472135954999579*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(1,26) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(1,27) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(1,28) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(1,29) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(1,30) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(1,31) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(2,24) = -0.5*cMOther[2]*mnuOther; 
  data->AEM_S(2,25) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(2,26) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(2,27) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(2,28) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(2,29) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(2,30) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(2,31) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(3,24) = -0.5*cMOther[3]*mnuOther; 
  data->AEM_S(3,25) = (-0.447213595499958*cMOther[6]*mnuOther)-0.5*cMOther[2]*mnuOther; 
  data->AEM_S(3,26) = (-0.447213595499958*cMOther[7]*mnuOther)-0.5*cMOther[1]*mnuOther; 
  data->AEM_S(3,27) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(3,28) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,29) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(3,30) = (-0.4*cMOther[7]*mnuOther)-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(3,31) = (-0.4*cMOther[6]*mnuOther)-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(4,24) = -0.5*cMOther[4]*mnuOther; 
  data->AEM_S(4,25) = -0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(4,26) = -0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(4,27) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(4,28) = (-0.31943828249997*cMOther[4]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(4,30) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(4,31) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(5,24) = -0.5*cMOther[5]*mnuOther; 
  data->AEM_S(5,25) = -0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(5,26) = -0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(5,27) = -0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(5,29) = (-0.31943828249997*cMOther[5]*mnuOther)-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(5,30) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(5,31) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(6,24) = -0.5*cMOther[6]*mnuOther; 
  data->AEM_S(6,25) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(6,26) = -0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(6,27) = (-0.4*cMOther[7]*mnuOther)-0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(6,28) = (-0.31943828249997*cMOther[6]*mnuOther)-0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(6,29) = -0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(6,30) = (-0.4472135954999579*cMOther[5]*mnuOther)-0.31943828249997*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
  data->AEM_S(6,31) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(7,24) = -0.5*cMOther[7]*mnuOther; 
  data->AEM_S(7,25) = -0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(7,26) = -0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(7,27) = (-0.4*cMOther[6]*mnuOther)-0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(7,28) = -0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(7,29) = (-0.31943828249997*cMOther[7]*mnuOther)-0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(7,30) = -0.4*cMOther[3]*mnuOther; 
  data->AEM_S(7,31) = (-0.31943828249997*cMOther[5]*mnuOther)-0.4472135954999579*cMOther[4]*mnuOther-0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, m1SelfX and uCrossSelfX ... // 
  data->AEM_S(8,0) = 0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(8,1) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(8,2) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(8,3) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(8,4) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(8,5) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(8,6) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(8,7) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(9,0) = 0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(9,1) = 0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(9,2) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(9,3) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(9,4) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(9,5) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(9,6) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(9,7) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(10,0) = 0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(10,1) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(10,2) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(10,3) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(10,4) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(10,5) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(10,6) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(10,7) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(11,0) = 0.5*m1rSelf[3]*mnuSelf; 
  data->AEM_S(11,1) = 0.447213595499958*m1rSelf[6]*mnuSelf+0.5*m1rSelf[2]*mnuSelf; 
  data->AEM_S(11,2) = 0.447213595499958*m1rSelf[7]*mnuSelf+0.5*m1rSelf[1]*mnuSelf; 
  data->AEM_S(11,3) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(11,4) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(11,5) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(11,6) = 0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(11,7) = 0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(12,0) = 0.5*m1rSelf[4]*mnuSelf; 
  data->AEM_S(12,1) = 0.4472135954999579*m1rSelf[1]*mnuSelf; 
  data->AEM_S(12,2) = 0.5000000000000001*m1rSelf[6]*mnuSelf; 
  data->AEM_S(12,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(12,4) = 0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(12,6) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(12,7) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(13,0) = 0.5*m1rSelf[5]*mnuSelf; 
  data->AEM_S(13,1) = 0.5000000000000001*m1rSelf[7]*mnuSelf; 
  data->AEM_S(13,2) = 0.4472135954999579*m1rSelf[2]*mnuSelf; 
  data->AEM_S(13,3) = 0.4472135954999579*m1rSelf[3]*mnuSelf; 
  data->AEM_S(13,5) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(13,6) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(13,7) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(14,0) = 0.5*m1rSelf[6]*mnuSelf; 
  data->AEM_S(14,1) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(14,2) = 0.5000000000000001*m1rSelf[4]*mnuSelf; 
  data->AEM_S(14,3) = 0.4*m1rSelf[7]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf; 
  data->AEM_S(14,4) = 0.31943828249997*m1rSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf; 
  data->AEM_S(14,5) = 0.4472135954999579*m1rSelf[6]*mnuSelf; 
  data->AEM_S(14,6) = 0.4472135954999579*m1rSelf[5]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(14,7) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(15,0) = 0.5*m1rSelf[7]*mnuSelf; 
  data->AEM_S(15,1) = 0.5000000000000001*m1rSelf[5]*mnuSelf; 
  data->AEM_S(15,2) = 0.447213595499958*m1rSelf[3]*mnuSelf; 
  data->AEM_S(15,3) = 0.4*m1rSelf[6]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf; 
  data->AEM_S(15,4) = 0.4472135954999579*m1rSelf[7]*mnuSelf; 
  data->AEM_S(15,5) = 0.31943828249997*m1rSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf; 
  data->AEM_S(15,6) = 0.4*m1rSelf[3]*mnuSelf; 
  data->AEM_S(15,7) = 0.31943828249997*m1rSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, m1OtherX and uCrossOtherX ... // 
  data->AEM_S(8,16) = 0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(8,17) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(8,18) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(8,19) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(8,20) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(8,21) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(8,22) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(8,23) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(9,16) = 0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(9,17) = 0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(9,18) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(9,19) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(9,20) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(9,21) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(9,22) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(9,23) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(10,16) = 0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(10,17) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(10,18) = 0.4472135954999579*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(10,19) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(10,20) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(10,21) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(10,22) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(10,23) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(11,16) = 0.5*m1rOther[3]*mnuOther; 
  data->AEM_S(11,17) = 0.447213595499958*m1rOther[6]*mnuOther+0.5*m1rOther[2]*mnuOther; 
  data->AEM_S(11,18) = 0.447213595499958*m1rOther[7]*mnuOther+0.5*m1rOther[1]*mnuOther; 
  data->AEM_S(11,19) = 0.4472135954999579*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(11,20) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(11,21) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(11,22) = 0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(11,23) = 0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(12,16) = 0.5*m1rOther[4]*mnuOther; 
  data->AEM_S(12,17) = 0.4472135954999579*m1rOther[1]*mnuOther; 
  data->AEM_S(12,18) = 0.5000000000000001*m1rOther[6]*mnuOther; 
  data->AEM_S(12,19) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(12,20) = 0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(12,22) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(12,23) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(13,16) = 0.5*m1rOther[5]*mnuOther; 
  data->AEM_S(13,17) = 0.5000000000000001*m1rOther[7]*mnuOther; 
  data->AEM_S(13,18) = 0.4472135954999579*m1rOther[2]*mnuOther; 
  data->AEM_S(13,19) = 0.4472135954999579*m1rOther[3]*mnuOther; 
  data->AEM_S(13,21) = 0.31943828249997*m1rOther[5]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(13,22) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(13,23) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(14,16) = 0.5*m1rOther[6]*mnuOther; 
  data->AEM_S(14,17) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(14,18) = 0.5000000000000001*m1rOther[4]*mnuOther; 
  data->AEM_S(14,19) = 0.4*m1rOther[7]*mnuOther+0.447213595499958*m1rOther[1]*mnuOther; 
  data->AEM_S(14,20) = 0.31943828249997*m1rOther[6]*mnuOther+0.5000000000000001*m1rOther[2]*mnuOther; 
  data->AEM_S(14,21) = 0.4472135954999579*m1rOther[6]*mnuOther; 
  data->AEM_S(14,22) = 0.4472135954999579*m1rOther[5]*mnuOther+0.31943828249997*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(14,23) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(15,16) = 0.5*m1rOther[7]*mnuOther; 
  data->AEM_S(15,17) = 0.5000000000000001*m1rOther[5]*mnuOther; 
  data->AEM_S(15,18) = 0.447213595499958*m1rOther[3]*mnuOther; 
  data->AEM_S(15,19) = 0.4*m1rOther[6]*mnuOther+0.447213595499958*m1rOther[2]*mnuOther; 
  data->AEM_S(15,20) = 0.4472135954999579*m1rOther[7]*mnuOther; 
  data->AEM_S(15,21) = 0.31943828249997*m1rOther[7]*mnuOther+0.5000000000000001*m1rOther[1]*mnuOther; 
  data->AEM_S(15,22) = 0.4*m1rOther[3]*mnuOther; 
  data->AEM_S(15,23) = 0.31943828249997*m1rOther[5]*mnuOther+0.4472135954999579*m1rOther[4]*mnuOther+0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of mnuM1Self+mnuM1Other. 
  mnuM1sum[0] += m1rSelf[0]*mnuSelf+m1rOther[0]*mnuOther; 
  mnuM1sum[1] += m1rSelf[1]*mnuSelf+m1rOther[1]*mnuOther; 
  mnuM1sum[2] += m1rSelf[2]*mnuSelf+m1rOther[2]*mnuOther; 
  mnuM1sum[3] += m1rSelf[3]*mnuSelf+m1rOther[3]*mnuOther; 
  mnuM1sum[4] += m1rSelf[4]*mnuSelf+m1rOther[4]*mnuOther; 
  mnuM1sum[5] += m1rSelf[5]*mnuSelf+m1rOther[5]*mnuOther; 
  mnuM1sum[6] += m1rSelf[6]*mnuSelf+m1rOther[6]*mnuOther; 
  mnuM1sum[7] += m1rSelf[7]*mnuSelf+m1rOther[7]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(8,8) = 1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(8,9) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(8,10) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(8,11) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(8,12) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(8,13) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(8,14) = 1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(8,15) = 1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(9,8) = 1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(9,9) = 1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(9,10) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(9,11) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(9,12) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(9,13) = 1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(9,14) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(9,15) = 1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(10,8) = 1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(10,9) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(10,10) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(10,11) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(10,12) = 1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(10,13) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(10,14) = 1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(10,15) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(11,8) = 1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(11,9) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(11,10) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(11,11) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(11,12) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(11,13) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(11,14) = 1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(11,15) = 1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(12,8) = 1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(12,9) = 1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(12,10) = 1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(12,11) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(12,12) = 0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(12,14) = 0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(12,15) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(13,8) = 1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(13,9) = 1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(13,10) = 1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(13,11) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(13,13) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(13,15) = 0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(14,8) = 1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(14,9) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(14,10) = 1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(14,11) = 1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(14,12) = 0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(14,13) = 1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(14,14) = 1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(14,15) = 1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(15,8) = 1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(15,9) = 1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(15,10) = 1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(15,11) = 1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(15,12) = 1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(15,13) = 0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(15,14) = 1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(15,15) = 0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(8,24) = 1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(8,25) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(8,26) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(8,27) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(8,28) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(8,29) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(8,30) = 1.5*m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(8,31) = 1.5*m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(9,24) = 1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(9,25) = 1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(9,26) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(9,27) = 1.341640786499874*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(9,28) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(9,29) = 1.5*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(9,30) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(9,31) = 1.5*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(10,24) = 1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(10,25) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(10,26) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(10,27) = 1.341640786499874*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(10,28) = 1.5*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(10,29) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(10,30) = 1.5*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(10,31) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(11,24) = 1.5*m0rOther[3]*mnuOther-0.5*cEOther[3]*mnuOther; 
  data->AEM_S(11,25) = 1.341640786499874*m0rOther[6]*mnuOther-0.447213595499958*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5*cEOther[2]*mnuOther; 
  data->AEM_S(11,26) = 1.341640786499874*m0rOther[7]*mnuOther-0.447213595499958*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5*cEOther[1]*mnuOther; 
  data->AEM_S(11,27) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(11,28) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(11,29) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(11,30) = 1.2*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(11,31) = 1.2*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(12,24) = 1.5*m0rOther[4]*mnuOther-0.5*cEOther[4]*mnuOther; 
  data->AEM_S(12,25) = 1.341640786499874*m0rOther[1]*mnuOther-0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(12,26) = 1.5*m0rOther[6]*mnuOther-0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(12,27) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(12,28) = 0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(12,30) = 0.9583148474999099*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(12,31) = 1.341640786499874*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(13,24) = 1.5*m0rOther[5]*mnuOther-0.5*cEOther[5]*mnuOther; 
  data->AEM_S(13,25) = 1.5*m0rOther[7]*mnuOther-0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(13,26) = 1.341640786499874*m0rOther[2]*mnuOther-0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(13,27) = 1.341640786499874*m0rOther[3]*mnuOther-0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(13,29) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(13,30) = 1.341640786499874*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(13,31) = 0.9583148474999099*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(14,24) = 1.5*m0rOther[6]*mnuOther-0.5*cEOther[6]*mnuOther; 
  data->AEM_S(14,25) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(14,26) = 1.5*m0rOther[4]*mnuOther-0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(14,27) = 1.2*m0rOther[7]*mnuOther-0.4*cEOther[7]*mnuOther+1.341640786499874*m0rOther[1]*mnuOther-0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(14,28) = 0.9583148474999099*m0rOther[6]*mnuOther-0.31943828249997*cEOther[6]*mnuOther+1.5*m0rOther[2]*mnuOther-0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(14,29) = 1.341640786499874*m0rOther[6]*mnuOther-0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(14,30) = 1.341640786499874*m0rOther[5]*mnuOther-0.4472135954999579*cEOther[5]*mnuOther+0.9583148474999099*m0rOther[4]*mnuOther-0.31943828249997*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
  data->AEM_S(14,31) = 1.2*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(15,24) = 1.5*m0rOther[7]*mnuOther-0.5*cEOther[7]*mnuOther; 
  data->AEM_S(15,25) = 1.5*m0rOther[5]*mnuOther-0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(15,26) = 1.341640786499874*m0rOther[3]*mnuOther-0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(15,27) = 1.2*m0rOther[6]*mnuOther-0.4*cEOther[6]*mnuOther+1.341640786499874*m0rOther[2]*mnuOther-0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(15,28) = 1.341640786499874*m0rOther[7]*mnuOther-0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(15,29) = 0.9583148474999099*m0rOther[7]*mnuOther-0.31943828249997*cEOther[7]*mnuOther+1.5*m0rOther[1]*mnuOther-0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(15,30) = 1.2*m0rOther[3]*mnuOther-0.4*cEOther[3]*mnuOther; 
  data->AEM_S(15,31) = 0.9583148474999099*m0rOther[5]*mnuOther-0.31943828249997*cEOther[5]*mnuOther+1.341640786499874*m0rOther[4]*mnuOther-0.4472135954999579*cEOther[4]*mnuOther+1.5*m0rOther[0]*mnuOther-0.5*cEOther[0]*mnuOther; 
 
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
 
  double m1Relax[8]; 
  // zero out array with sum of momentum relaxation terms. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    m1Relax[vd] = 0.0; 
  } 
 
  // ... Relaxation block from weak multiply of mSelf, nuSelf, M0Self and uCrossSelfX ... // 
  data->AEM_S(16,0) = 0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(16,1) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(16,2) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(16,3) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(16,4) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(16,5) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(16,6) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(16,7) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(17,0) = 0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(17,2) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,3) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(17,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(17,5) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(17,6) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(17,7) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(18,0) = 0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,1) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(18,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(18,3) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(18,4) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(18,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(18,6) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(18,7) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,0) = 0.5*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,1) = 0.447213595499958*m0rSelf[6]*mnuSelf+0.5*m0rSelf[2]*mnuSelf; 
  data->AEM_S(19,2) = 0.447213595499958*m0rSelf[7]*mnuSelf+0.5*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(19,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(19,6) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(19,7) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,0) = 0.5*m0rSelf[4]*mnuSelf; 
  data->AEM_S(20,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf; 
  data->AEM_S(20,2) = 0.5000000000000001*m0rSelf[6]*mnuSelf; 
  data->AEM_S(20,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(20,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(20,6) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(20,7) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(21,0) = 0.5*m0rSelf[5]*mnuSelf; 
  data->AEM_S(21,1) = 0.5000000000000001*m0rSelf[7]*mnuSelf; 
  data->AEM_S(21,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf; 
  data->AEM_S(21,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf; 
  data->AEM_S(21,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(21,6) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(21,7) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,0) = 0.5*m0rSelf[6]*mnuSelf; 
  data->AEM_S(22,1) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(22,2) = 0.5000000000000001*m0rSelf[4]*mnuSelf; 
  data->AEM_S(22,3) = 0.4*m0rSelf[7]*mnuSelf+0.447213595499958*m0rSelf[1]*mnuSelf; 
  data->AEM_S(22,4) = 0.31943828249997*m0rSelf[6]*mnuSelf+0.5000000000000001*m0rSelf[2]*mnuSelf; 
  data->AEM_S(22,5) = 0.4472135954999579*m0rSelf[6]*mnuSelf; 
  data->AEM_S(22,6) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
  data->AEM_S(22,7) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,0) = 0.5*m0rSelf[7]*mnuSelf; 
  data->AEM_S(23,1) = 0.5000000000000001*m0rSelf[5]*mnuSelf; 
  data->AEM_S(23,2) = 0.447213595499958*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,3) = 0.4*m0rSelf[6]*mnuSelf+0.447213595499958*m0rSelf[2]*mnuSelf; 
  data->AEM_S(23,4) = 0.4472135954999579*m0rSelf[7]*mnuSelf; 
  data->AEM_S(23,5) = 0.31943828249997*m0rSelf[7]*mnuSelf+0.5000000000000001*m0rSelf[1]*mnuSelf; 
  data->AEM_S(23,6) = 0.4*m0rSelf[3]*mnuSelf; 
  data->AEM_S(23,7) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf; 
 
  // ... Block from correction to (self) 1st moment of collision operator ... // 
  data->AEM_S(16,8) = -0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(16,9) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(16,10) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(16,11) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(16,12) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(16,13) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(16,14) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(16,15) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(17,8) = -0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(17,9) = (-0.4472135954999579*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(17,10) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(17,11) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(17,12) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(17,13) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(17,14) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(17,15) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(18,8) = -0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(18,9) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(18,10) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(18,11) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(18,12) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(18,13) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(18,14) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(18,15) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(19,8) = -0.5*cMSelf[3]*mnuSelf; 
  data->AEM_S(19,9) = (-0.447213595499958*cMSelf[6]*mnuSelf)-0.5*cMSelf[2]*mnuSelf; 
  data->AEM_S(19,10) = (-0.447213595499958*cMSelf[7]*mnuSelf)-0.5*cMSelf[1]*mnuSelf; 
  data->AEM_S(19,11) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(19,12) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(19,13) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(19,14) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(19,15) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(20,8) = -0.5*cMSelf[4]*mnuSelf; 
  data->AEM_S(20,9) = -0.4472135954999579*cMSelf[1]*mnuSelf; 
  data->AEM_S(20,10) = -0.5000000000000001*cMSelf[6]*mnuSelf; 
  data->AEM_S(20,11) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(20,12) = (-0.31943828249997*cMSelf[4]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(20,14) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(20,15) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(21,8) = -0.5*cMSelf[5]*mnuSelf; 
  data->AEM_S(21,9) = -0.5000000000000001*cMSelf[7]*mnuSelf; 
  data->AEM_S(21,10) = -0.4472135954999579*cMSelf[2]*mnuSelf; 
  data->AEM_S(21,11) = -0.4472135954999579*cMSelf[3]*mnuSelf; 
  data->AEM_S(21,13) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(21,14) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(21,15) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(22,8) = -0.5*cMSelf[6]*mnuSelf; 
  data->AEM_S(22,9) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(22,10) = -0.5000000000000001*cMSelf[4]*mnuSelf; 
  data->AEM_S(22,11) = (-0.4*cMSelf[7]*mnuSelf)-0.447213595499958*cMSelf[1]*mnuSelf; 
  data->AEM_S(22,12) = (-0.31943828249997*cMSelf[6]*mnuSelf)-0.5000000000000001*cMSelf[2]*mnuSelf; 
  data->AEM_S(22,13) = -0.4472135954999579*cMSelf[6]*mnuSelf; 
  data->AEM_S(22,14) = (-0.4472135954999579*cMSelf[5]*mnuSelf)-0.31943828249997*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
  data->AEM_S(22,15) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(23,8) = -0.5*cMSelf[7]*mnuSelf; 
  data->AEM_S(23,9) = -0.5000000000000001*cMSelf[5]*mnuSelf; 
  data->AEM_S(23,10) = -0.447213595499958*cMSelf[3]*mnuSelf; 
  data->AEM_S(23,11) = (-0.4*cMSelf[6]*mnuSelf)-0.447213595499958*cMSelf[2]*mnuSelf; 
  data->AEM_S(23,12) = -0.4472135954999579*cMSelf[7]*mnuSelf; 
  data->AEM_S(23,13) = (-0.31943828249997*cMSelf[7]*mnuSelf)-0.5000000000000001*cMSelf[1]*mnuSelf; 
  data->AEM_S(23,14) = -0.4*cMSelf[3]*mnuSelf; 
  data->AEM_S(23,15) = (-0.31943828249997*cMSelf[5]*mnuSelf)-0.4472135954999579*cMSelf[4]*mnuSelf-0.5*cMSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(16,16) = -0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(16,17) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(16,18) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(16,19) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(16,20) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(16,21) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(16,22) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(16,23) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(17,16) = -0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(17,17) = (-0.4472135954999579*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(17,18) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(17,19) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(17,20) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(17,21) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(17,22) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(17,23) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(18,16) = -0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(18,17) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(18,18) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(18,19) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(18,20) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(18,21) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(18,22) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(18,23) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(19,16) = -0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(19,17) = (-0.447213595499958*m0rOther[6]*mnuOther)-0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(19,18) = (-0.447213595499958*m0rOther[7]*mnuOther)-0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(19,19) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(19,20) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(19,21) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(19,22) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(19,23) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(20,16) = -0.5*m0rOther[4]*mnuOther; 
  data->AEM_S(20,17) = -0.4472135954999579*m0rOther[1]*mnuOther; 
  data->AEM_S(20,18) = -0.5000000000000001*m0rOther[6]*mnuOther; 
  data->AEM_S(20,19) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(20,20) = (-0.31943828249997*m0rOther[4]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(20,22) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(20,23) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(21,16) = -0.5*m0rOther[5]*mnuOther; 
  data->AEM_S(21,17) = -0.5000000000000001*m0rOther[7]*mnuOther; 
  data->AEM_S(21,18) = -0.4472135954999579*m0rOther[2]*mnuOther; 
  data->AEM_S(21,19) = -0.4472135954999579*m0rOther[3]*mnuOther; 
  data->AEM_S(21,21) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(21,22) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(21,23) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(22,16) = -0.5*m0rOther[6]*mnuOther; 
  data->AEM_S(22,17) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(22,18) = -0.5000000000000001*m0rOther[4]*mnuOther; 
  data->AEM_S(22,19) = (-0.4*m0rOther[7]*mnuOther)-0.447213595499958*m0rOther[1]*mnuOther; 
  data->AEM_S(22,20) = (-0.31943828249997*m0rOther[6]*mnuOther)-0.5000000000000001*m0rOther[2]*mnuOther; 
  data->AEM_S(22,21) = -0.4472135954999579*m0rOther[6]*mnuOther; 
  data->AEM_S(22,22) = (-0.4472135954999579*m0rOther[5]*mnuOther)-0.31943828249997*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(22,23) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(23,16) = -0.5*m0rOther[7]*mnuOther; 
  data->AEM_S(23,17) = -0.5000000000000001*m0rOther[5]*mnuOther; 
  data->AEM_S(23,18) = -0.447213595499958*m0rOther[3]*mnuOther; 
  data->AEM_S(23,19) = (-0.4*m0rOther[6]*mnuOther)-0.447213595499958*m0rOther[2]*mnuOther; 
  data->AEM_S(23,20) = -0.4472135954999579*m0rOther[7]*mnuOther; 
  data->AEM_S(23,21) = (-0.31943828249997*m0rOther[7]*mnuOther)-0.5000000000000001*m0rOther[1]*mnuOther; 
  data->AEM_S(23,22) = -0.4*m0rOther[3]*mnuOther; 
  data->AEM_S(23,23) = (-0.31943828249997*m0rOther[5]*mnuOther)-0.4472135954999579*m0rOther[4]*mnuOther-0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to (other) 1st moment of collision operator ... // 
  data->AEM_S(16,24) = 0.5*cMOther[0]*mnuOther; 
  data->AEM_S(16,25) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(16,26) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(16,27) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(16,28) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(16,29) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(16,30) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(16,31) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(17,24) = 0.5*cMOther[1]*mnuOther; 
  data->AEM_S(17,25) = 0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(17,26) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(17,27) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(17,28) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(17,29) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(17,30) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(17,31) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(18,24) = 0.5*cMOther[2]*mnuOther; 
  data->AEM_S(18,25) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(18,26) = 0.4472135954999579*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(18,27) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(18,28) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(18,29) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(18,30) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(18,31) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(19,24) = 0.5*cMOther[3]*mnuOther; 
  data->AEM_S(19,25) = 0.447213595499958*cMOther[6]*mnuOther+0.5*cMOther[2]*mnuOther; 
  data->AEM_S(19,26) = 0.447213595499958*cMOther[7]*mnuOther+0.5*cMOther[1]*mnuOther; 
  data->AEM_S(19,27) = 0.4472135954999579*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(19,28) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(19,29) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(19,30) = 0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(19,31) = 0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(20,24) = 0.5*cMOther[4]*mnuOther; 
  data->AEM_S(20,25) = 0.4472135954999579*cMOther[1]*mnuOther; 
  data->AEM_S(20,26) = 0.5000000000000001*cMOther[6]*mnuOther; 
  data->AEM_S(20,27) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(20,28) = 0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(20,30) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(20,31) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(21,24) = 0.5*cMOther[5]*mnuOther; 
  data->AEM_S(21,25) = 0.5000000000000001*cMOther[7]*mnuOther; 
  data->AEM_S(21,26) = 0.4472135954999579*cMOther[2]*mnuOther; 
  data->AEM_S(21,27) = 0.4472135954999579*cMOther[3]*mnuOther; 
  data->AEM_S(21,29) = 0.31943828249997*cMOther[5]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(21,30) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(21,31) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(22,24) = 0.5*cMOther[6]*mnuOther; 
  data->AEM_S(22,25) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(22,26) = 0.5000000000000001*cMOther[4]*mnuOther; 
  data->AEM_S(22,27) = 0.4*cMOther[7]*mnuOther+0.447213595499958*cMOther[1]*mnuOther; 
  data->AEM_S(22,28) = 0.31943828249997*cMOther[6]*mnuOther+0.5000000000000001*cMOther[2]*mnuOther; 
  data->AEM_S(22,29) = 0.4472135954999579*cMOther[6]*mnuOther; 
  data->AEM_S(22,30) = 0.4472135954999579*cMOther[5]*mnuOther+0.31943828249997*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
  data->AEM_S(22,31) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(23,24) = 0.5*cMOther[7]*mnuOther; 
  data->AEM_S(23,25) = 0.5000000000000001*cMOther[5]*mnuOther; 
  data->AEM_S(23,26) = 0.447213595499958*cMOther[3]*mnuOther; 
  data->AEM_S(23,27) = 0.4*cMOther[6]*mnuOther+0.447213595499958*cMOther[2]*mnuOther; 
  data->AEM_S(23,28) = 0.4472135954999579*cMOther[7]*mnuOther; 
  data->AEM_S(23,29) = 0.31943828249997*cMOther[7]*mnuOther+0.5000000000000001*cMOther[1]*mnuOther; 
  data->AEM_S(23,30) = 0.4*cMOther[3]*mnuOther; 
  data->AEM_S(23,31) = 0.31943828249997*cMOther[5]*mnuOther+0.4472135954999579*cMOther[4]*mnuOther+0.5*cMOther[0]*mnuOther; 
 
  // ... Block from weak multiply of mSelf, nuSelf, (m1SelfX-uSelfX*m0Self) and uCrossSelfX ... // 
  data->AEM_S(24,0) = (-0.25*m0rSelf[7]*uSelf[7]*mnuSelf)-0.25*m0rSelf[6]*uSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(24,1) = (-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(24,2) = (-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,3) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(24,4) = (-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(24,5) = (-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(24,6) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(24,7) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,0) = (-0.2500000000000001*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2500000000000001*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[2]*uSelf[3]*mnuSelf-0.25*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,1) = (-0.45*m0rSelf[7]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(25,2) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(25,3) = (-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf)-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(25,4) = (-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(25,5) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(25,6) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(25,7) = (-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf)-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(26,0) = (-0.223606797749979*m0rSelf[3]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[3]*m0rSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[1]*uSelf[3]*mnuSelf-0.25*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,1) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,2) = (-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.45*m0rSelf[6]*uSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.45*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(26,3) = (-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(26,4) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(26,5) = (-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(26,6) = (-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(26,7) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,0) = (-0.2*m0rSelf[6]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[2]*uSelf[7]*mnuSelf-0.2*uSelf[6]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[3]*m0rSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[3]*mnuSelf+0.5*m1rSelf[3]*mnuSelf-0.25*uSelf[0]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[1]*uSelf[2]*mnuSelf-0.25*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,1) = (-0.4024922359499621*m0rSelf[3]*uSelf[7]*mnuSelf)-0.4024922359499621*uSelf[3]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[5]*uSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.447213595499958*m1rSelf[6]*mnuSelf-0.2*uSelf[5]*m0rSelf[6]*mnuSelf-0.3928571428571429*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[1]*uSelf[3]*mnuSelf-0.45*uSelf[1]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[2]*mnuSelf+0.5*m1rSelf[2]*mnuSelf-0.25*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,2) = (-0.3928571428571429*m0rSelf[5]*uSelf[7]*mnuSelf)-0.2*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.447213595499958*m1rSelf[7]*mnuSelf-0.3928571428571429*uSelf[5]*m0rSelf[7]*mnuSelf-0.2*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.45*m0rSelf[2]*uSelf[3]*mnuSelf-0.45*uSelf[2]*m0rSelf[3]*mnuSelf-0.25*m0rSelf[0]*uSelf[1]*mnuSelf+0.5*m1rSelf[1]*mnuSelf-0.25*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,3) = (-0.7071428571428572*m0rSelf[7]*uSelf[7]*mnuSelf)-0.4024922359499621*m0rSelf[1]*uSelf[7]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[6]*uSelf[6]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[6]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.2*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.2*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.81*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(27,4) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,5) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(27,6) = (-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(27,7) = (-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf)-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(28,0) = (-0.223606797749979*m0rSelf[7]*uSelf[7]*mnuSelf)-0.159719141249985*m0rSelf[6]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[6]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.25*m0rSelf[0]*uSelf[4]*mnuSelf+0.5*m1rSelf[4]*mnuSelf-0.25*uSelf[0]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(28,1) = (-0.223606797749979*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.4472135954999579*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(28,2) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[6]*mnuSelf+0.5000000000000001*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[2]*uSelf[4]*mnuSelf-0.25*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(28,3) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(28,4) = (-0.3928571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.5357142857142857*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[6]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[5]*uSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.25*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(28,5) = (-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(28,6) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(28,7) = (-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,0) = (-0.159719141249985*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2500000000000001*m0rSelf[1]*uSelf[7]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[6]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.25*m0rSelf[0]*uSelf[5]*mnuSelf+0.5*m1rSelf[5]*mnuSelf-0.25*uSelf[0]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(29,1) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[7]*mnuSelf+0.5000000000000001*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[1]*uSelf[5]*mnuSelf-0.25*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,2) = (-0.3928571428571429*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3928571428571429*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.4472135954999579*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,3) = (-0.3513821107499669*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571429*m0rSelf[2]*uSelf[7]*mnuSelf-0.3513821107499669*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571428*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.4472135954999579*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(29,4) = (-0.1428571428571428*m0rSelf[7]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[1]*uSelf[7]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[6]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[6]*mnuSelf-0.25*m0rSelf[4]*uSelf[5]*mnuSelf-0.25*uSelf[4]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[3]*mnuSelf; 
  data->AEM_S(29,5) = (-0.5357142857142857*m0rSelf[7]*uSelf[7]*mnuSelf)-0.159719141249985*m0rSelf[1]*uSelf[7]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[6]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.25*m0rSelf[4]*uSelf[4]*mnuSelf-0.3928571428571428*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.25*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(29,6) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(29,7) = (-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(30,0) = (-0.2*m0rSelf[3]*uSelf[7]*mnuSelf)-0.2*uSelf[3]*m0rSelf[7]*mnuSelf-0.223606797749979*m0rSelf[5]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[6]*mnuSelf-0.25*m0rSelf[0]*uSelf[6]*mnuSelf+0.5*m1rSelf[6]*mnuSelf-0.223606797749979*uSelf[5]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[4]*m0rSelf[6]*mnuSelf-0.25*uSelf[0]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[2]*uSelf[4]*mnuSelf-0.2500000000000001*uSelf[2]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(30,1) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.2*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.2*uSelf[2]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[6]*mnuSelf-0.3928571428571428*uSelf[1]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[3]*uSelf[5]*mnuSelf-0.2*uSelf[3]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,2) = (-0.351382110749967*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2*m0rSelf[1]*uSelf[7]*mnuSelf-0.2*uSelf[1]*m0rSelf[7]*mnuSelf-0.2874944542499729*m0rSelf[6]*uSelf[6]*mnuSelf-0.45*m0rSelf[2]*uSelf[6]*mnuSelf-0.45*uSelf[2]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[4]*uSelf[4]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[4]*mnuSelf+0.5000000000000001*m1rSelf[4]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[1]*mnuSelf; 
  data->AEM_S(30,3) = (-0.3513821107499669*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3513821107499669*m0rSelf[4]*uSelf[7]*mnuSelf-0.2*m0rSelf[0]*uSelf[7]*mnuSelf+0.4*m1rSelf[7]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[7]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[7]*mnuSelf-0.2*uSelf[0]*m0rSelf[7]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[6]*mnuSelf-0.7071428571428572*uSelf[3]*m0rSelf[6]*mnuSelf-0.2*m0rSelf[1]*uSelf[5]*mnuSelf-0.2*uSelf[1]*m0rSelf[5]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[4]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[2]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[2]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[1]*mnuSelf+0.447213595499958*m1rSelf[1]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(30,4) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.1428571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[6]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[6]*mnuSelf+0.31943828249997*m1rSelf[6]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.5357142857142857*uSelf[4]*m0rSelf[6]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[5]*mnuSelf-0.159719141249985*m0rSelf[2]*uSelf[4]*mnuSelf-0.159719141249985*uSelf[2]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[1]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[1]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[2]*mnuSelf+0.5000000000000001*m1rSelf[2]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(30,5) = (-0.3513821107499669*m0rSelf[3]*uSelf[7]*mnuSelf)-0.3513821107499669*uSelf[3]*m0rSelf[7]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[6]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[6]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[6]*mnuSelf+0.4472135954999579*m1rSelf[6]*mnuSelf-0.3928571428571428*uSelf[5]*m0rSelf[6]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[6]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[1]*uSelf[3]*mnuSelf-0.2*uSelf[1]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(30,6) = (-0.6173469387755102*m0rSelf[7]*uSelf[7]*mnuSelf)-0.351382110749967*m0rSelf[1]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[7]*mnuSelf-0.9642857142857143*m0rSelf[6]*uSelf[6]*mnuSelf-0.2874944542499729*m0rSelf[2]*uSelf[6]*mnuSelf-0.2874944542499729*uSelf[2]*m0rSelf[6]*mnuSelf-0.3928571428571428*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[5]*mnuSelf+0.4472135954999579*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[5]*mnuSelf-0.5357142857142857*m0rSelf[4]*uSelf[4]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[4]*mnuSelf+0.31943828249997*m1rSelf[4]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.45*m0rSelf[2]*uSelf[2]*mnuSelf-0.3928571428571428*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
  data->AEM_S(30,7) = (-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,0) = (-0.159719141249985*m0rSelf[5]*uSelf[7]*mnuSelf)-0.223606797749979*m0rSelf[4]*uSelf[7]*mnuSelf-0.25*m0rSelf[0]*uSelf[7]*mnuSelf+0.5*m1rSelf[7]*mnuSelf-0.159719141249985*uSelf[5]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[7]*mnuSelf-0.25*uSelf[0]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[3]*uSelf[6]*mnuSelf-0.2*uSelf[3]*m0rSelf[6]*mnuSelf-0.2500000000000001*m0rSelf[1]*uSelf[5]*mnuSelf-0.2500000000000001*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[3]*mnuSelf-0.223606797749979*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,1) = (-0.2874944542499729*m0rSelf[7]*uSelf[7]*mnuSelf)-0.45*m0rSelf[1]*uSelf[7]*mnuSelf-0.45*uSelf[1]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[6]*uSelf[6]*mnuSelf-0.2*m0rSelf[2]*uSelf[6]*mnuSelf-0.2*uSelf[2]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[5]*uSelf[5]*mnuSelf-0.223606797749979*m0rSelf[4]*uSelf[5]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[5]*mnuSelf+0.5000000000000001*m1rSelf[5]*mnuSelf-0.223606797749979*uSelf[4]*m0rSelf[5]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[5]*mnuSelf-0.4024922359499621*m0rSelf[3]*uSelf[3]*mnuSelf-0.223606797749979*m0rSelf[2]*uSelf[2]*mnuSelf; 
  data->AEM_S(31,2) = (-0.351382110749967*m0rSelf[6]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[2]*uSelf[7]*mnuSelf-0.351382110749967*uSelf[6]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[2]*m0rSelf[7]*mnuSelf-0.2*m0rSelf[1]*uSelf[6]*mnuSelf-0.2*uSelf[1]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[3]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[3]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[3]*uSelf[4]*mnuSelf-0.2*uSelf[3]*m0rSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[3]*mnuSelf+0.447213595499958*m1rSelf[3]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[2]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,3) = (-0.7071428571428572*m0rSelf[3]*uSelf[7]*mnuSelf)-0.7071428571428572*uSelf[3]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[5]*uSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[4]*uSelf[6]*mnuSelf-0.2*m0rSelf[0]*uSelf[6]*mnuSelf+0.4*m1rSelf[6]*mnuSelf-0.3513821107499669*uSelf[5]*m0rSelf[6]*mnuSelf-0.3513821107499669*uSelf[4]*m0rSelf[6]*mnuSelf-0.2*uSelf[0]*m0rSelf[6]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[5]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[4]*mnuSelf-0.2*uSelf[2]*m0rSelf[4]*mnuSelf-0.4024922359499621*m0rSelf[1]*uSelf[3]*mnuSelf-0.4024922359499621*uSelf[1]*m0rSelf[3]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[2]*mnuSelf+0.447213595499958*m1rSelf[2]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,4) = (-0.1428571428571428*m0rSelf[5]*uSelf[7]*mnuSelf)-0.3928571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[7]*mnuSelf+0.4472135954999579*m1rSelf[7]*mnuSelf-0.1428571428571428*uSelf[5]*m0rSelf[7]*mnuSelf-0.3928571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[5]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[5]*mnuSelf-0.2*m0rSelf[2]*uSelf[3]*mnuSelf-0.2*uSelf[2]*m0rSelf[3]*mnuSelf; 
  data->AEM_S(31,5) = (-0.5357142857142857*m0rSelf[5]*uSelf[7]*mnuSelf)-0.1428571428571428*m0rSelf[4]*uSelf[7]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[7]*mnuSelf+0.31943828249997*m1rSelf[7]*mnuSelf-0.5357142857142857*uSelf[5]*m0rSelf[7]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[7]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[7]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[6]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[6]*mnuSelf-0.159719141249985*m0rSelf[1]*uSelf[5]*mnuSelf-0.159719141249985*uSelf[1]*m0rSelf[5]*mnuSelf-0.223606797749979*m0rSelf[1]*uSelf[4]*mnuSelf-0.223606797749979*uSelf[1]*m0rSelf[4]*mnuSelf-0.3928571428571429*m0rSelf[2]*uSelf[3]*mnuSelf-0.3928571428571429*uSelf[2]*m0rSelf[3]*mnuSelf-0.2500000000000001*m0rSelf[0]*uSelf[1]*mnuSelf+0.5000000000000001*m1rSelf[1]*mnuSelf-0.2500000000000001*uSelf[0]*m0rSelf[1]*mnuSelf; 
  data->AEM_S(31,6) = (-0.6173469387755102*m0rSelf[6]*uSelf[7]*mnuSelf)-0.351382110749967*m0rSelf[2]*uSelf[7]*mnuSelf-0.6173469387755102*uSelf[6]*m0rSelf[7]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[7]*mnuSelf-0.351382110749967*m0rSelf[1]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[1]*m0rSelf[6]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[5]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[5]*mnuSelf-0.3513821107499669*m0rSelf[3]*uSelf[4]*mnuSelf-0.3513821107499669*uSelf[3]*m0rSelf[4]*mnuSelf-0.2*m0rSelf[0]*uSelf[3]*mnuSelf+0.4*m1rSelf[3]*mnuSelf-0.2*uSelf[0]*m0rSelf[3]*mnuSelf-0.2*m0rSelf[1]*uSelf[2]*mnuSelf-0.2*uSelf[1]*m0rSelf[2]*mnuSelf; 
  data->AEM_S(31,7) = (-0.9642857142857143*m0rSelf[7]*uSelf[7]*mnuSelf)-0.2874944542499729*m0rSelf[1]*uSelf[7]*mnuSelf-0.2874944542499729*uSelf[1]*m0rSelf[7]*mnuSelf-0.6173469387755102*m0rSelf[6]*uSelf[6]*mnuSelf-0.351382110749967*m0rSelf[2]*uSelf[6]*mnuSelf-0.351382110749967*uSelf[2]*m0rSelf[6]*mnuSelf-0.5357142857142857*m0rSelf[5]*uSelf[5]*mnuSelf-0.1428571428571428*m0rSelf[4]*uSelf[5]*mnuSelf-0.159719141249985*m0rSelf[0]*uSelf[5]*mnuSelf+0.31943828249997*m1rSelf[5]*mnuSelf-0.1428571428571428*uSelf[4]*m0rSelf[5]*mnuSelf-0.159719141249985*uSelf[0]*m0rSelf[5]*mnuSelf-0.3928571428571428*m0rSelf[4]*uSelf[4]*mnuSelf-0.223606797749979*m0rSelf[0]*uSelf[4]*mnuSelf+0.4472135954999579*m1rSelf[4]*mnuSelf-0.223606797749979*uSelf[0]*m0rSelf[4]*mnuSelf-0.7071428571428572*m0rSelf[3]*uSelf[3]*mnuSelf-0.3928571428571428*m0rSelf[2]*uSelf[2]*mnuSelf-0.45*m0rSelf[1]*uSelf[1]*mnuSelf-0.25*m0rSelf[0]*uSelf[0]*mnuSelf+0.5*m1rSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, (m1OtherX-uOtherX*m0Other) and uCrossOtherX ... // 
  data->AEM_S(24,16) = 0.25*m0rOther[7]*uOther[7]*mnuOther+0.25*m0rOther[6]*uOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(24,17) = 0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(24,18) = 0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(24,19) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(24,20) = 0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(24,21) = 0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(24,22) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(24,23) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(25,16) = 0.2500000000000001*m0rOther[5]*uOther[7]*mnuOther+0.2500000000000001*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[3]*uOther[6]*mnuOther+0.223606797749979*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.25*m0rOther[2]*uOther[3]*mnuOther+0.25*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(25,17) = 0.45*m0rOther[7]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(25,18) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(25,19) = 0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(25,20) = 0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(25,21) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(25,22) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(25,23) = 0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(26,16) = 0.223606797749979*m0rOther[3]*uOther[7]*mnuOther+0.223606797749979*uOther[3]*m0rOther[7]*mnuOther+0.2500000000000001*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.25*m0rOther[1]*uOther[3]*mnuOther+0.25*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(26,17) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(26,18) = 0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.45*m0rOther[6]*uOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.45*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(26,19) = 0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(26,20) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(26,21) = 0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(26,22) = 0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(26,23) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(27,16) = 0.2*m0rOther[6]*uOther[7]*mnuOther+0.223606797749979*m0rOther[2]*uOther[7]*mnuOther+0.2*uOther[6]*m0rOther[7]*mnuOther+0.223606797749979*uOther[2]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[6]*mnuOther+0.223606797749979*uOther[1]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[3]*uOther[5]*mnuOther+0.223606797749979*uOther[3]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[4]*mnuOther+0.223606797749979*uOther[3]*m0rOther[4]*mnuOther+0.25*m0rOther[0]*uOther[3]*mnuOther-0.5*m1rOther[3]*mnuOther+0.25*uOther[0]*m0rOther[3]*mnuOther+0.25*m0rOther[1]*uOther[2]*mnuOther+0.25*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(27,17) = 0.4024922359499621*m0rOther[3]*uOther[7]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[7]*mnuOther+0.2*m0rOther[5]*uOther[6]*mnuOther+0.3928571428571429*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.447213595499958*m1rOther[6]*mnuOther+0.2*uOther[5]*m0rOther[6]*mnuOther+0.3928571428571429*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.45*m0rOther[1]*uOther[3]*mnuOther+0.45*uOther[1]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[2]*mnuOther-0.5*m1rOther[2]*mnuOther+0.25*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(27,18) = 0.3928571428571429*m0rOther[5]*uOther[7]*mnuOther+0.2*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.447213595499958*m1rOther[7]*mnuOther+0.3928571428571429*uOther[5]*m0rOther[7]*mnuOther+0.2*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[6]*mnuOther+0.4024922359499621*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.45*m0rOther[2]*uOther[3]*mnuOther+0.45*uOther[2]*m0rOther[3]*mnuOther+0.25*m0rOther[0]*uOther[1]*mnuOther-0.5*m1rOther[1]*mnuOther+0.25*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(27,19) = 0.7071428571428572*m0rOther[7]*uOther[7]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[7]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[6]*uOther[6]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[6]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.2*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.2*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.81*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(27,20) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(27,21) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(27,22) = 0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(27,23) = 0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(28,16) = 0.223606797749979*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[6]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[6]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.25*m0rOther[0]*uOther[4]*mnuOther-0.5*m1rOther[4]*mnuOther+0.25*uOther[0]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(28,17) = 0.223606797749979*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[6]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.4472135954999579*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(28,18) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[6]*mnuOther-0.5000000000000001*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[6]*mnuOther+0.25*m0rOther[2]*uOther[4]*mnuOther+0.25*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(28,19) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(28,20) = 0.3928571428571428*m0rOther[7]*uOther[7]*mnuOther+0.5357142857142857*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[2]*uOther[6]*mnuOther+0.159719141249985*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[5]*uOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.25*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(28,21) = 0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(28,22) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(28,23) = 0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(29,16) = 0.159719141249985*m0rOther[7]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[7]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[6]*uOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.25*m0rOther[0]*uOther[5]*mnuOther-0.5*m1rOther[5]*mnuOther+0.25*uOther[0]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(29,17) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[7]*mnuOther-0.5000000000000001*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.25*m0rOther[1]*uOther[5]*mnuOther+0.25*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(29,18) = 0.3928571428571429*m0rOther[3]*uOther[7]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*uOther[4]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.4472135954999579*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(29,19) = 0.3513821107499669*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[7]*mnuOther+0.3513821107499669*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571428*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.4472135954999579*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(29,20) = 0.1428571428571428*m0rOther[7]*uOther[7]*mnuOther+0.223606797749979*m0rOther[1]*uOther[7]*mnuOther+0.223606797749979*uOther[1]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[6]*uOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[6]*mnuOther+0.223606797749979*uOther[2]*m0rOther[6]*mnuOther+0.25*m0rOther[4]*uOther[5]*mnuOther+0.25*uOther[4]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[3]*mnuOther; 
  data->AEM_S(29,21) = 0.5357142857142857*m0rOther[7]*uOther[7]*mnuOther+0.159719141249985*m0rOther[1]*uOther[7]*mnuOther+0.159719141249985*uOther[1]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[6]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.25*m0rOther[4]*uOther[4]*mnuOther+0.3928571428571428*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.25*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(29,22) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(29,23) = 0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(30,16) = 0.2*m0rOther[3]*uOther[7]*mnuOther+0.2*uOther[3]*m0rOther[7]*mnuOther+0.223606797749979*m0rOther[5]*uOther[6]*mnuOther+0.159719141249985*m0rOther[4]*uOther[6]*mnuOther+0.25*m0rOther[0]*uOther[6]*mnuOther-0.5*m1rOther[6]*mnuOther+0.223606797749979*uOther[5]*m0rOther[6]*mnuOther+0.159719141249985*uOther[4]*m0rOther[6]*mnuOther+0.25*uOther[0]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[2]*uOther[4]*mnuOther+0.2500000000000001*uOther[2]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[1]*uOther[3]*mnuOther+0.223606797749979*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(30,17) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.2*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.2*uOther[2]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[6]*mnuOther+0.3928571428571428*uOther[1]*m0rOther[6]*mnuOther+0.2*m0rOther[3]*uOther[5]*mnuOther+0.2*uOther[3]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[4]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,18) = 0.351382110749967*m0rOther[7]*uOther[7]*mnuOther+0.2*m0rOther[1]*uOther[7]*mnuOther+0.2*uOther[1]*m0rOther[7]*mnuOther+0.2874944542499729*m0rOther[6]*uOther[6]*mnuOther+0.45*m0rOther[2]*uOther[6]*mnuOther+0.45*uOther[2]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[4]*uOther[4]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[4]*mnuOther-0.5000000000000001*m1rOther[4]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[1]*mnuOther; 
  data->AEM_S(30,19) = 0.3513821107499669*m0rOther[5]*uOther[7]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[7]*mnuOther+0.2*m0rOther[0]*uOther[7]*mnuOther-0.4*m1rOther[7]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[7]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[7]*mnuOther+0.2*uOther[0]*m0rOther[7]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[6]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[6]*mnuOther+0.2*m0rOther[1]*uOther[5]*mnuOther+0.2*uOther[1]*m0rOther[5]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[4]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[2]*uOther[3]*mnuOther+0.4024922359499621*uOther[2]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[1]*mnuOther-0.447213595499958*m1rOther[1]*mnuOther+0.223606797749979*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(30,20) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.1428571428571428*m0rOther[5]*uOther[6]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[6]*mnuOther+0.159719141249985*m0rOther[0]*uOther[6]*mnuOther-0.31943828249997*m1rOther[6]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[6]*mnuOther+0.5357142857142857*uOther[4]*m0rOther[6]*mnuOther+0.159719141249985*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[5]*mnuOther+0.223606797749979*uOther[2]*m0rOther[5]*mnuOther+0.159719141249985*m0rOther[2]*uOther[4]*mnuOther+0.159719141249985*uOther[2]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[1]*uOther[3]*mnuOther+0.3928571428571429*uOther[1]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[2]*mnuOther-0.5000000000000001*m1rOther[2]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(30,21) = 0.3513821107499669*m0rOther[3]*uOther[7]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[7]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[6]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[6]*mnuOther+0.223606797749979*m0rOther[0]*uOther[6]*mnuOther-0.4472135954999579*m1rOther[6]*mnuOther+0.3928571428571428*uOther[5]*m0rOther[6]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[6]*mnuOther+0.223606797749979*uOther[0]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[2]*uOther[4]*mnuOther+0.223606797749979*uOther[2]*m0rOther[4]*mnuOther+0.2*m0rOther[1]*uOther[3]*mnuOther+0.2*uOther[1]*m0rOther[3]*mnuOther; 
  data->AEM_S(30,22) = 0.6173469387755102*m0rOther[7]*uOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[7]*mnuOther+0.351382110749967*uOther[1]*m0rOther[7]*mnuOther+0.9642857142857143*m0rOther[6]*uOther[6]*mnuOther+0.2874944542499729*m0rOther[2]*uOther[6]*mnuOther+0.2874944542499729*uOther[2]*m0rOther[6]*mnuOther+0.3928571428571428*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.223606797749979*m0rOther[0]*uOther[5]*mnuOther-0.4472135954999579*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.223606797749979*uOther[0]*m0rOther[5]*mnuOther+0.5357142857142857*m0rOther[4]*uOther[4]*mnuOther+0.159719141249985*m0rOther[0]*uOther[4]*mnuOther-0.31943828249997*m1rOther[4]*mnuOther+0.159719141249985*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.45*m0rOther[2]*uOther[2]*mnuOther+0.3928571428571428*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
  data->AEM_S(30,23) = 0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,16) = 0.159719141249985*m0rOther[5]*uOther[7]*mnuOther+0.223606797749979*m0rOther[4]*uOther[7]*mnuOther+0.25*m0rOther[0]*uOther[7]*mnuOther-0.5*m1rOther[7]*mnuOther+0.159719141249985*uOther[5]*m0rOther[7]*mnuOther+0.223606797749979*uOther[4]*m0rOther[7]*mnuOther+0.25*uOther[0]*m0rOther[7]*mnuOther+0.2*m0rOther[3]*uOther[6]*mnuOther+0.2*uOther[3]*m0rOther[6]*mnuOther+0.2500000000000001*m0rOther[1]*uOther[5]*mnuOther+0.2500000000000001*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[2]*uOther[3]*mnuOther+0.223606797749979*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(31,17) = 0.2874944542499729*m0rOther[7]*uOther[7]*mnuOther+0.45*m0rOther[1]*uOther[7]*mnuOther+0.45*uOther[1]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[6]*uOther[6]*mnuOther+0.2*m0rOther[2]*uOther[6]*mnuOther+0.2*uOther[2]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[5]*uOther[5]*mnuOther+0.223606797749979*m0rOther[4]*uOther[5]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[5]*mnuOther-0.5000000000000001*m1rOther[5]*mnuOther+0.223606797749979*uOther[4]*m0rOther[5]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[5]*mnuOther+0.4024922359499621*m0rOther[3]*uOther[3]*mnuOther+0.223606797749979*m0rOther[2]*uOther[2]*mnuOther; 
  data->AEM_S(31,18) = 0.351382110749967*m0rOther[6]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[7]*mnuOther+0.351382110749967*uOther[6]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[2]*m0rOther[7]*mnuOther+0.2*m0rOther[1]*uOther[6]*mnuOther+0.2*uOther[1]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[3]*uOther[5]*mnuOther+0.3928571428571429*uOther[3]*m0rOther[5]*mnuOther+0.2*m0rOther[3]*uOther[4]*mnuOther+0.2*uOther[3]*m0rOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[3]*mnuOther-0.447213595499958*m1rOther[3]*mnuOther+0.223606797749979*uOther[0]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[1]*uOther[2]*mnuOther+0.223606797749979*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,19) = 0.7071428571428572*m0rOther[3]*uOther[7]*mnuOther+0.7071428571428572*uOther[3]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[5]*uOther[6]*mnuOther+0.3513821107499669*m0rOther[4]*uOther[6]*mnuOther+0.2*m0rOther[0]*uOther[6]*mnuOther-0.4*m1rOther[6]*mnuOther+0.3513821107499669*uOther[5]*m0rOther[6]*mnuOther+0.3513821107499669*uOther[4]*m0rOther[6]*mnuOther+0.2*uOther[0]*m0rOther[6]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[5]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[4]*mnuOther+0.2*uOther[2]*m0rOther[4]*mnuOther+0.4024922359499621*m0rOther[1]*uOther[3]*mnuOther+0.4024922359499621*uOther[1]*m0rOther[3]*mnuOther+0.223606797749979*m0rOther[0]*uOther[2]*mnuOther-0.447213595499958*m1rOther[2]*mnuOther+0.223606797749979*uOther[0]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,20) = 0.1428571428571428*m0rOther[5]*uOther[7]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[7]*mnuOther+0.223606797749979*m0rOther[0]*uOther[7]*mnuOther-0.4472135954999579*m1rOther[7]*mnuOther+0.1428571428571428*uOther[5]*m0rOther[7]*mnuOther+0.3928571428571428*uOther[4]*m0rOther[7]*mnuOther+0.223606797749979*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.223606797749979*m0rOther[1]*uOther[5]*mnuOther+0.223606797749979*uOther[1]*m0rOther[5]*mnuOther+0.2*m0rOther[2]*uOther[3]*mnuOther+0.2*uOther[2]*m0rOther[3]*mnuOther; 
  data->AEM_S(31,21) = 0.5357142857142857*m0rOther[5]*uOther[7]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[7]*mnuOther+0.159719141249985*m0rOther[0]*uOther[7]*mnuOther-0.31943828249997*m1rOther[7]*mnuOther+0.5357142857142857*uOther[5]*m0rOther[7]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[7]*mnuOther+0.159719141249985*uOther[0]*m0rOther[7]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[6]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[6]*mnuOther+0.159719141249985*m0rOther[1]*uOther[5]*mnuOther+0.159719141249985*uOther[1]*m0rOther[5]*mnuOther+0.223606797749979*m0rOther[1]*uOther[4]*mnuOther+0.223606797749979*uOther[1]*m0rOther[4]*mnuOther+0.3928571428571429*m0rOther[2]*uOther[3]*mnuOther+0.3928571428571429*uOther[2]*m0rOther[3]*mnuOther+0.2500000000000001*m0rOther[0]*uOther[1]*mnuOther-0.5000000000000001*m1rOther[1]*mnuOther+0.2500000000000001*uOther[0]*m0rOther[1]*mnuOther; 
  data->AEM_S(31,22) = 0.6173469387755102*m0rOther[6]*uOther[7]*mnuOther+0.351382110749967*m0rOther[2]*uOther[7]*mnuOther+0.6173469387755102*uOther[6]*m0rOther[7]*mnuOther+0.351382110749967*uOther[2]*m0rOther[7]*mnuOther+0.351382110749967*m0rOther[1]*uOther[6]*mnuOther+0.351382110749967*uOther[1]*m0rOther[6]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[5]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[5]*mnuOther+0.3513821107499669*m0rOther[3]*uOther[4]*mnuOther+0.3513821107499669*uOther[3]*m0rOther[4]*mnuOther+0.2*m0rOther[0]*uOther[3]*mnuOther-0.4*m1rOther[3]*mnuOther+0.2*uOther[0]*m0rOther[3]*mnuOther+0.2*m0rOther[1]*uOther[2]*mnuOther+0.2*uOther[1]*m0rOther[2]*mnuOther; 
  data->AEM_S(31,23) = 0.9642857142857143*m0rOther[7]*uOther[7]*mnuOther+0.2874944542499729*m0rOther[1]*uOther[7]*mnuOther+0.2874944542499729*uOther[1]*m0rOther[7]*mnuOther+0.6173469387755102*m0rOther[6]*uOther[6]*mnuOther+0.351382110749967*m0rOther[2]*uOther[6]*mnuOther+0.351382110749967*uOther[2]*m0rOther[6]*mnuOther+0.5357142857142857*m0rOther[5]*uOther[5]*mnuOther+0.1428571428571428*m0rOther[4]*uOther[5]*mnuOther+0.159719141249985*m0rOther[0]*uOther[5]*mnuOther-0.31943828249997*m1rOther[5]*mnuOther+0.1428571428571428*uOther[4]*m0rOther[5]*mnuOther+0.159719141249985*uOther[0]*m0rOther[5]*mnuOther+0.3928571428571428*m0rOther[4]*uOther[4]*mnuOther+0.223606797749979*m0rOther[0]*uOther[4]*mnuOther-0.4472135954999579*m1rOther[4]*mnuOther+0.223606797749979*uOther[0]*m0rOther[4]*mnuOther+0.7071428571428572*m0rOther[3]*uOther[3]*mnuOther+0.3928571428571428*m0rOther[2]*uOther[2]*mnuOther+0.45*m0rOther[1]*uOther[1]*mnuOther+0.25*m0rOther[0]*uOther[0]*mnuOther-0.5*m1rOther[0]*mnuOther; 
 
  // ... Contribution to RHS vector from component 1 of momentum relaxation. 
  m1Relax[0] += betaGreenep1*(m1rOther[0]*deltaSelf-1.0*m1rSelf[0]*deltaSelf)*mnuSelf+m1rSelf[0]*mnuSelf-1.0*m1rOther[0]*mnuOther; 
  m1Relax[1] += betaGreenep1*(m1rOther[1]*deltaSelf-1.0*m1rSelf[1]*deltaSelf)*mnuSelf+m1rSelf[1]*mnuSelf-1.0*m1rOther[1]*mnuOther; 
  m1Relax[2] += betaGreenep1*(m1rOther[2]*deltaSelf-1.0*m1rSelf[2]*deltaSelf)*mnuSelf+m1rSelf[2]*mnuSelf-1.0*m1rOther[2]*mnuOther; 
  m1Relax[3] += betaGreenep1*(m1rOther[3]*deltaSelf-1.0*m1rSelf[3]*deltaSelf)*mnuSelf+m1rSelf[3]*mnuSelf-1.0*m1rOther[3]*mnuOther; 
  m1Relax[4] += betaGreenep1*(m1rOther[4]*deltaSelf-1.0*m1rSelf[4]*deltaSelf)*mnuSelf+m1rSelf[4]*mnuSelf-1.0*m1rOther[4]*mnuOther; 
  m1Relax[5] += betaGreenep1*(m1rOther[5]*deltaSelf-1.0*m1rSelf[5]*deltaSelf)*mnuSelf+m1rSelf[5]*mnuSelf-1.0*m1rOther[5]*mnuOther; 
  m1Relax[6] += betaGreenep1*(m1rOther[6]*deltaSelf-1.0*m1rSelf[6]*deltaSelf)*mnuSelf+m1rSelf[6]*mnuSelf-1.0*m1rOther[6]*mnuOther; 
  m1Relax[7] += betaGreenep1*(m1rOther[7]*deltaSelf-1.0*m1rSelf[7]*deltaSelf)*mnuSelf+m1rSelf[7]*mnuSelf-1.0*m1rOther[7]*mnuOther; 
 
  double ucMSelf[8]; 
  // Zero out array with dot product of uSelf and cMSelf. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    ucMSelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
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
  data->AEM_S(24,8) = 0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(24,9) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(24,10) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(24,11) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(24,12) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(24,13) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(24,14) = 0.5*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(24,15) = 0.5*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(25,8) = 0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(25,9) = 0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(25,10) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(25,11) = 0.447213595499958*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(25,12) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(25,13) = 0.5000000000000001*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(25,14) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(25,15) = 0.5000000000000001*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(26,8) = 0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(26,9) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(26,10) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(26,11) = 0.447213595499958*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(26,12) = 0.5000000000000001*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(26,13) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(26,14) = 0.5000000000000001*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(26,15) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(27,8) = 0.5*ucMSelf[3]*mnuSelf+1.5*m0rSelf[3]*mnuSelf-0.5*cESelf[3]*mnuSelf; 
  data->AEM_S(27,9) = 0.447213595499958*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.447213595499958*cESelf[6]*mnuSelf+0.5*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5*cESelf[2]*mnuSelf; 
  data->AEM_S(27,10) = 0.447213595499958*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.447213595499958*cESelf[7]*mnuSelf+0.5*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5*cESelf[1]*mnuSelf; 
  data->AEM_S(27,11) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(27,12) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(27,13) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(27,14) = 0.4*ucMSelf[7]*mnuSelf+1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.447213595499958*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(27,15) = 0.4*ucMSelf[6]*mnuSelf+1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.447213595499958*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(28,8) = 0.5*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5*cESelf[4]*mnuSelf; 
  data->AEM_S(28,9) = 0.4472135954999579*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.4472135954999579*cESelf[1]*mnuSelf; 
  data->AEM_S(28,10) = 0.5000000000000001*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5000000000000001*cESelf[6]*mnuSelf; 
  data->AEM_S(28,11) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(28,12) = 0.31943828249997*ucMSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(28,14) = 0.31943828249997*ucMSelf[6]*mnuSelf+0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+0.5000000000000001*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(28,15) = 0.4472135954999579*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(29,8) = 0.5*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5*cESelf[5]*mnuSelf; 
  data->AEM_S(29,9) = 0.5000000000000001*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5000000000000001*cESelf[7]*mnuSelf; 
  data->AEM_S(29,10) = 0.4472135954999579*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.4472135954999579*cESelf[2]*mnuSelf; 
  data->AEM_S(29,11) = 0.4472135954999579*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.4472135954999579*cESelf[3]*mnuSelf; 
  data->AEM_S(29,13) = 0.31943828249997*ucMSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(29,14) = 0.4472135954999579*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(29,15) = 0.31943828249997*ucMSelf[7]*mnuSelf+0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+0.5000000000000001*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(30,8) = 0.5*ucMSelf[6]*mnuSelf+1.5*m0rSelf[6]*mnuSelf-0.5*cESelf[6]*mnuSelf; 
  data->AEM_S(30,9) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(30,10) = 0.5000000000000001*ucMSelf[4]*mnuSelf+1.5*m0rSelf[4]*mnuSelf-0.5000000000000001*cESelf[4]*mnuSelf; 
  data->AEM_S(30,11) = 0.4*ucMSelf[7]*mnuSelf+1.2*m0rSelf[7]*mnuSelf-0.4*cESelf[7]*mnuSelf+0.447213595499958*ucMSelf[1]*mnuSelf+1.341640786499874*m0rSelf[1]*mnuSelf-0.447213595499958*cESelf[1]*mnuSelf; 
  data->AEM_S(30,12) = 0.31943828249997*ucMSelf[6]*mnuSelf+0.9583148474999099*m0rSelf[6]*mnuSelf-0.31943828249997*cESelf[6]*mnuSelf+0.5000000000000001*ucMSelf[2]*mnuSelf+1.5*m0rSelf[2]*mnuSelf-0.5000000000000001*cESelf[2]*mnuSelf; 
  data->AEM_S(30,13) = 0.4472135954999579*ucMSelf[6]*mnuSelf+1.341640786499874*m0rSelf[6]*mnuSelf-0.4472135954999579*cESelf[6]*mnuSelf; 
  data->AEM_S(30,14) = 0.4472135954999579*ucMSelf[5]*mnuSelf+1.341640786499874*m0rSelf[5]*mnuSelf-0.4472135954999579*cESelf[5]*mnuSelf+0.31943828249997*ucMSelf[4]*mnuSelf+0.9583148474999099*m0rSelf[4]*mnuSelf-0.31943828249997*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
  data->AEM_S(30,15) = 0.4*ucMSelf[3]*mnuSelf+1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(31,8) = 0.5*ucMSelf[7]*mnuSelf+1.5*m0rSelf[7]*mnuSelf-0.5*cESelf[7]*mnuSelf; 
  data->AEM_S(31,9) = 0.5000000000000001*ucMSelf[5]*mnuSelf+1.5*m0rSelf[5]*mnuSelf-0.5000000000000001*cESelf[5]*mnuSelf; 
  data->AEM_S(31,10) = 0.447213595499958*ucMSelf[3]*mnuSelf+1.341640786499874*m0rSelf[3]*mnuSelf-0.447213595499958*cESelf[3]*mnuSelf; 
  data->AEM_S(31,11) = 0.4*ucMSelf[6]*mnuSelf+1.2*m0rSelf[6]*mnuSelf-0.4*cESelf[6]*mnuSelf+0.447213595499958*ucMSelf[2]*mnuSelf+1.341640786499874*m0rSelf[2]*mnuSelf-0.447213595499958*cESelf[2]*mnuSelf; 
  data->AEM_S(31,12) = 0.4472135954999579*ucMSelf[7]*mnuSelf+1.341640786499874*m0rSelf[7]*mnuSelf-0.4472135954999579*cESelf[7]*mnuSelf; 
  data->AEM_S(31,13) = 0.31943828249997*ucMSelf[7]*mnuSelf+0.9583148474999099*m0rSelf[7]*mnuSelf-0.31943828249997*cESelf[7]*mnuSelf+0.5000000000000001*ucMSelf[1]*mnuSelf+1.5*m0rSelf[1]*mnuSelf-0.5000000000000001*cESelf[1]*mnuSelf; 
  data->AEM_S(31,14) = 0.4*ucMSelf[3]*mnuSelf+1.2*m0rSelf[3]*mnuSelf-0.4*cESelf[3]*mnuSelf; 
  data->AEM_S(31,15) = 0.31943828249997*ucMSelf[5]*mnuSelf+0.9583148474999099*m0rSelf[5]*mnuSelf-0.31943828249997*cESelf[5]*mnuSelf+0.4472135954999579*ucMSelf[4]*mnuSelf+1.341640786499874*m0rSelf[4]*mnuSelf-0.4472135954999579*cESelf[4]*mnuSelf+0.5*ucMSelf[0]*mnuSelf+1.5*m0rSelf[0]*mnuSelf-0.5*cESelf[0]*mnuSelf; 
 
  double ucMOther[8]; 
  // Zero out array with dot product of uOther and cMOther. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    ucMOther[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
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
  data->AEM_S(24,24) = (-0.5*ucMOther[0]*mnuOther)-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(24,25) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(24,26) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(24,27) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(24,28) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(24,29) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(24,30) = (-0.5*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5*cEOther[6]*mnuOther; 
  data->AEM_S(24,31) = (-0.5*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5*cEOther[7]*mnuOther; 
  data->AEM_S(25,24) = (-0.5*ucMOther[1]*mnuOther)-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(25,25) = (-0.4472135954999579*ucMOther[4]*mnuOther)-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(25,26) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(25,27) = (-0.447213595499958*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-0.5*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(25,28) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(25,29) = (-0.5000000000000001*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(25,30) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(25,31) = (-0.5000000000000001*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(26,24) = (-0.5*ucMOther[2]*mnuOther)-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(26,25) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(26,26) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(26,27) = (-0.447213595499958*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-0.5*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(26,28) = (-0.5000000000000001*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(26,29) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(26,30) = (-0.5000000000000001*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(26,31) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(27,24) = (-0.5*ucMOther[3]*mnuOther)-1.5*m0rOther[3]*mnuOther+0.5*cEOther[3]*mnuOther; 
  data->AEM_S(27,25) = (-0.447213595499958*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.447213595499958*cEOther[6]*mnuOther-0.5*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5*cEOther[2]*mnuOther; 
  data->AEM_S(27,26) = (-0.447213595499958*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.447213595499958*cEOther[7]*mnuOther-0.5*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5*cEOther[1]*mnuOther; 
  data->AEM_S(27,27) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.4472135954999579*ucMOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(27,28) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(27,29) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(27,30) = (-0.4*ucMOther[7]*mnuOther)-1.2*m0rOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.447213595499958*ucMOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(27,31) = (-0.4*ucMOther[6]*mnuOther)-1.2*m0rOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.447213595499958*ucMOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(28,24) = (-0.5*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5*cEOther[4]*mnuOther; 
  data->AEM_S(28,25) = (-0.4472135954999579*ucMOther[1]*mnuOther)-1.341640786499874*m0rOther[1]*mnuOther+0.4472135954999579*cEOther[1]*mnuOther; 
  data->AEM_S(28,26) = (-0.5000000000000001*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5000000000000001*cEOther[6]*mnuOther; 
  data->AEM_S(28,27) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(28,28) = (-0.31943828249997*ucMOther[4]*mnuOther)-0.9583148474999099*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(28,30) = (-0.31943828249997*ucMOther[6]*mnuOther)-0.9583148474999099*m0rOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-0.5000000000000001*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(28,31) = (-0.4472135954999579*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(29,24) = (-0.5*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5*cEOther[5]*mnuOther; 
  data->AEM_S(29,25) = (-0.5000000000000001*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5000000000000001*cEOther[7]*mnuOther; 
  data->AEM_S(29,26) = (-0.4472135954999579*ucMOther[2]*mnuOther)-1.341640786499874*m0rOther[2]*mnuOther+0.4472135954999579*cEOther[2]*mnuOther; 
  data->AEM_S(29,27) = (-0.4472135954999579*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.4472135954999579*cEOther[3]*mnuOther; 
  data->AEM_S(29,29) = (-0.31943828249997*ucMOther[5]*mnuOther)-0.9583148474999099*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(29,30) = (-0.4472135954999579*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(29,31) = (-0.31943828249997*ucMOther[7]*mnuOther)-0.9583148474999099*m0rOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-0.5000000000000001*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(30,24) = (-0.5*ucMOther[6]*mnuOther)-1.5*m0rOther[6]*mnuOther+0.5*cEOther[6]*mnuOther; 
  data->AEM_S(30,25) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(30,26) = (-0.5000000000000001*ucMOther[4]*mnuOther)-1.5*m0rOther[4]*mnuOther+0.5000000000000001*cEOther[4]*mnuOther; 
  data->AEM_S(30,27) = (-0.4*ucMOther[7]*mnuOther)-1.2*m0rOther[7]*mnuOther+0.4*cEOther[7]*mnuOther-0.447213595499958*ucMOther[1]*mnuOther-1.341640786499874*m0rOther[1]*mnuOther+0.447213595499958*cEOther[1]*mnuOther; 
  data->AEM_S(30,28) = (-0.31943828249997*ucMOther[6]*mnuOther)-0.9583148474999099*m0rOther[6]*mnuOther+0.31943828249997*cEOther[6]*mnuOther-0.5000000000000001*ucMOther[2]*mnuOther-1.5*m0rOther[2]*mnuOther+0.5000000000000001*cEOther[2]*mnuOther; 
  data->AEM_S(30,29) = (-0.4472135954999579*ucMOther[6]*mnuOther)-1.341640786499874*m0rOther[6]*mnuOther+0.4472135954999579*cEOther[6]*mnuOther; 
  data->AEM_S(30,30) = (-0.4472135954999579*ucMOther[5]*mnuOther)-1.341640786499874*m0rOther[5]*mnuOther+0.4472135954999579*cEOther[5]*mnuOther-0.31943828249997*ucMOther[4]*mnuOther-0.9583148474999099*m0rOther[4]*mnuOther+0.31943828249997*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
  data->AEM_S(30,31) = (-0.4*ucMOther[3]*mnuOther)-1.2*m0rOther[3]*mnuOther+0.4*cEOther[3]*mnuOther; 
  data->AEM_S(31,24) = (-0.5*ucMOther[7]*mnuOther)-1.5*m0rOther[7]*mnuOther+0.5*cEOther[7]*mnuOther; 
  data->AEM_S(31,25) = (-0.5000000000000001*ucMOther[5]*mnuOther)-1.5*m0rOther[5]*mnuOther+0.5000000000000001*cEOther[5]*mnuOther; 
  data->AEM_S(31,26) = (-0.447213595499958*ucMOther[3]*mnuOther)-1.341640786499874*m0rOther[3]*mnuOther+0.447213595499958*cEOther[3]*mnuOther; 
  data->AEM_S(31,27) = (-0.4*ucMOther[6]*mnuOther)-1.2*m0rOther[6]*mnuOther+0.4*cEOther[6]*mnuOther-0.447213595499958*ucMOther[2]*mnuOther-1.341640786499874*m0rOther[2]*mnuOther+0.447213595499958*cEOther[2]*mnuOther; 
  data->AEM_S(31,28) = (-0.4472135954999579*ucMOther[7]*mnuOther)-1.341640786499874*m0rOther[7]*mnuOther+0.4472135954999579*cEOther[7]*mnuOther; 
  data->AEM_S(31,29) = (-0.31943828249997*ucMOther[7]*mnuOther)-0.9583148474999099*m0rOther[7]*mnuOther+0.31943828249997*cEOther[7]*mnuOther-0.5000000000000001*ucMOther[1]*mnuOther-1.5*m0rOther[1]*mnuOther+0.5000000000000001*cEOther[1]*mnuOther; 
  data->AEM_S(31,30) = (-0.4*ucMOther[3]*mnuOther)-1.2*m0rOther[3]*mnuOther+0.4*cEOther[3]*mnuOther; 
  data->AEM_S(31,31) = (-0.31943828249997*ucMOther[5]*mnuOther)-0.9583148474999099*m0rOther[5]*mnuOther+0.31943828249997*cEOther[5]*mnuOther-0.4472135954999579*ucMOther[4]*mnuOther-1.341640786499874*m0rOther[4]*mnuOther+0.4472135954999579*cEOther[4]*mnuOther-0.5*ucMOther[0]*mnuOther-1.5*m0rOther[0]*mnuOther+0.5*cEOther[0]*mnuOther; 
 
  double kinESelf[8]; 
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    kinESelf[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
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
  for (unsigned short int vd=0; vd<1; vd++) 
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
  // zero out array with dot product of u and m1. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += 0.5*m1rSelf[a0+7]*uSelf[a0+7]-0.5*m1rOther[a0+7]*uSelf[a0+7]-0.5*m1rSelf[a0+7]*uOther[a0+7]+0.5*m1rOther[a0+7]*uOther[a0+7]+0.5*m1rSelf[a0+6]*uSelf[a0+6]-0.5*m1rOther[a0+6]*uSelf[a0+6]-0.5*m1rSelf[a0+6]*uOther[a0+6]+0.5*m1rOther[a0+6]*uOther[a0+6]+0.5*m1rSelf[a0+5]*uSelf[a0+5]-0.5*m1rOther[a0+5]*uSelf[a0+5]-0.5*m1rSelf[a0+5]*uOther[a0+5]+0.5*m1rOther[a0+5]*uOther[a0+5]+0.5*m1rSelf[a0+4]*uSelf[a0+4]-0.5*m1rOther[a0+4]*uSelf[a0+4]-0.5*m1rSelf[a0+4]*uOther[a0+4]+0.5*m1rOther[a0+4]*uOther[a0+4]+0.5*m1rSelf[a0+3]*uSelf[a0+3]-0.5*m1rOther[a0+3]*uSelf[a0+3]-0.5*m1rSelf[a0+3]*uOther[a0+3]+0.5*m1rOther[a0+3]*uOther[a0+3]+0.5*m1rSelf[a0+2]*uSelf[a0+2]-0.5*m1rOther[a0+2]*uSelf[a0+2]-0.5*m1rSelf[a0+2]*uOther[a0+2]+0.5*m1rOther[a0+2]*uOther[a0+2]+0.5*m1rSelf[a0+1]*uSelf[a0+1]-0.5*m1rOther[a0+1]*uSelf[a0+1]-0.5*m1rSelf[a0+1]*uOther[a0+1]+0.5*m1rOther[a0+1]*uOther[a0+1]+0.5*m1rSelf[a0]*uSelf[a0]-0.5*m1rOther[a0]*uSelf[a0]-0.5*m1rSelf[a0]*uOther[a0]+0.5*m1rOther[a0]*uOther[a0]; 
    relKinE[1] += 0.5000000000000001*m1rSelf[a0+5]*uSelf[a0+7]-0.5000000000000001*m1rOther[a0+5]*uSelf[a0+7]-0.5000000000000001*m1rSelf[a0+5]*uOther[a0+7]+0.5000000000000001*m1rOther[a0+5]*uOther[a0+7]+0.5000000000000001*uSelf[a0+5]*m1rSelf[a0+7]-0.5000000000000001*uOther[a0+5]*m1rSelf[a0+7]-0.5000000000000001*uSelf[a0+5]*m1rOther[a0+7]+0.5000000000000001*uOther[a0+5]*m1rOther[a0+7]+0.447213595499958*m1rSelf[a0+3]*uSelf[a0+6]-0.447213595499958*m1rOther[a0+3]*uSelf[a0+6]-0.447213595499958*m1rSelf[a0+3]*uOther[a0+6]+0.447213595499958*m1rOther[a0+3]*uOther[a0+6]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+6]-0.447213595499958*uOther[a0+3]*m1rSelf[a0+6]-0.447213595499958*uSelf[a0+3]*m1rOther[a0+6]+0.447213595499958*uOther[a0+3]*m1rOther[a0+6]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+4]-0.4472135954999579*m1rOther[a0+1]*uSelf[a0+4]-0.4472135954999579*m1rSelf[a0+1]*uOther[a0+4]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+4]+0.4472135954999579*uSelf[a0+1]*m1rSelf[a0+4]-0.4472135954999579*uOther[a0+1]*m1rSelf[a0+4]-0.4472135954999579*uSelf[a0+1]*m1rOther[a0+4]+0.4472135954999579*uOther[a0+1]*m1rOther[a0+4]+0.5*m1rSelf[a0+2]*uSelf[a0+3]-0.5*m1rOther[a0+2]*uSelf[a0+3]-0.5*m1rSelf[a0+2]*uOther[a0+3]+0.5*m1rOther[a0+2]*uOther[a0+3]+0.5*uSelf[a0+2]*m1rSelf[a0+3]-0.5*uOther[a0+2]*m1rSelf[a0+3]-0.5*uSelf[a0+2]*m1rOther[a0+3]+0.5*uOther[a0+2]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+1]-0.5*m1rOther[a0]*uSelf[a0+1]-0.5*m1rSelf[a0]*uOther[a0+1]+0.5*m1rOther[a0]*uOther[a0+1]+0.5*uSelf[a0]*m1rSelf[a0+1]-0.5*uOther[a0]*m1rSelf[a0+1]-0.5*uSelf[a0]*m1rOther[a0+1]+0.5*uOther[a0]*m1rOther[a0+1]; 
    relKinE[2] += 0.447213595499958*m1rSelf[a0+3]*uSelf[a0+7]-0.447213595499958*m1rOther[a0+3]*uSelf[a0+7]-0.447213595499958*m1rSelf[a0+3]*uOther[a0+7]+0.447213595499958*m1rOther[a0+3]*uOther[a0+7]+0.447213595499958*uSelf[a0+3]*m1rSelf[a0+7]-0.447213595499958*uOther[a0+3]*m1rSelf[a0+7]-0.447213595499958*uSelf[a0+3]*m1rOther[a0+7]+0.447213595499958*uOther[a0+3]*m1rOther[a0+7]+0.5000000000000001*m1rSelf[a0+4]*uSelf[a0+6]-0.5000000000000001*m1rOther[a0+4]*uSelf[a0+6]-0.5000000000000001*m1rSelf[a0+4]*uOther[a0+6]+0.5000000000000001*m1rOther[a0+4]*uOther[a0+6]+0.5000000000000001*uSelf[a0+4]*m1rSelf[a0+6]-0.5000000000000001*uOther[a0+4]*m1rSelf[a0+6]-0.5000000000000001*uSelf[a0+4]*m1rOther[a0+6]+0.5000000000000001*uOther[a0+4]*m1rOther[a0+6]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+5]-0.4472135954999579*m1rOther[a0+2]*uSelf[a0+5]-0.4472135954999579*m1rSelf[a0+2]*uOther[a0+5]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+5]+0.4472135954999579*uSelf[a0+2]*m1rSelf[a0+5]-0.4472135954999579*uOther[a0+2]*m1rSelf[a0+5]-0.4472135954999579*uSelf[a0+2]*m1rOther[a0+5]+0.4472135954999579*uOther[a0+2]*m1rOther[a0+5]+0.5*m1rSelf[a0+1]*uSelf[a0+3]-0.5*m1rOther[a0+1]*uSelf[a0+3]-0.5*m1rSelf[a0+1]*uOther[a0+3]+0.5*m1rOther[a0+1]*uOther[a0+3]+0.5*uSelf[a0+1]*m1rSelf[a0+3]-0.5*uOther[a0+1]*m1rSelf[a0+3]-0.5*uSelf[a0+1]*m1rOther[a0+3]+0.5*uOther[a0+1]*m1rOther[a0+3]+0.5*m1rSelf[a0]*uSelf[a0+2]-0.5*m1rOther[a0]*uSelf[a0+2]-0.5*m1rSelf[a0]*uOther[a0+2]+0.5*m1rOther[a0]*uOther[a0+2]+0.5*uSelf[a0]*m1rSelf[a0+2]-0.5*uOther[a0]*m1rSelf[a0+2]-0.5*uSelf[a0]*m1rOther[a0+2]+0.5*uOther[a0]*m1rOther[a0+2]; 
    relKinE[3] += 0.4*m1rSelf[a0+6]*uSelf[a0+7]-0.4*m1rOther[a0+6]*uSelf[a0+7]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+7]-0.447213595499958*m1rOther[a0+2]*uSelf[a0+7]-0.4*m1rSelf[a0+6]*uOther[a0+7]+0.4*m1rOther[a0+6]*uOther[a0+7]-0.447213595499958*m1rSelf[a0+2]*uOther[a0+7]+0.447213595499958*m1rOther[a0+2]*uOther[a0+7]+0.4*uSelf[a0+6]*m1rSelf[a0+7]-0.4*uOther[a0+6]*m1rSelf[a0+7]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+7]-0.447213595499958*uOther[a0+2]*m1rSelf[a0+7]-0.4*uSelf[a0+6]*m1rOther[a0+7]+0.4*uOther[a0+6]*m1rOther[a0+7]-0.447213595499958*uSelf[a0+2]*m1rOther[a0+7]+0.447213595499958*uOther[a0+2]*m1rOther[a0+7]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+6]-0.447213595499958*m1rOther[a0+1]*uSelf[a0+6]-0.447213595499958*m1rSelf[a0+1]*uOther[a0+6]+0.447213595499958*m1rOther[a0+1]*uOther[a0+6]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+6]-0.447213595499958*uOther[a0+1]*m1rSelf[a0+6]-0.447213595499958*uSelf[a0+1]*m1rOther[a0+6]+0.447213595499958*uOther[a0+1]*m1rOther[a0+6]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+5]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+5]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+5]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+5]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+5]-0.4472135954999579*uOther[a0+3]*m1rSelf[a0+5]-0.4472135954999579*uSelf[a0+3]*m1rOther[a0+5]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+4]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+4]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+4]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+4]+0.4472135954999579*uSelf[a0+3]*m1rSelf[a0+4]-0.4472135954999579*uOther[a0+3]*m1rSelf[a0+4]-0.4472135954999579*uSelf[a0+3]*m1rOther[a0+4]+0.4472135954999579*uOther[a0+3]*m1rOther[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+3]-0.5*m1rOther[a0]*uSelf[a0+3]-0.5*m1rSelf[a0]*uOther[a0+3]+0.5*m1rOther[a0]*uOther[a0+3]+0.5*uSelf[a0]*m1rSelf[a0+3]-0.5*uOther[a0]*m1rSelf[a0+3]-0.5*uSelf[a0]*m1rOther[a0+3]+0.5*uOther[a0]*m1rOther[a0+3]+0.5*m1rSelf[a0+1]*uSelf[a0+2]-0.5*m1rOther[a0+1]*uSelf[a0+2]-0.5*m1rSelf[a0+1]*uOther[a0+2]+0.5*m1rOther[a0+1]*uOther[a0+2]+0.5*uSelf[a0+1]*m1rSelf[a0+2]-0.5*uOther[a0+1]*m1rSelf[a0+2]-0.5*uSelf[a0+1]*m1rOther[a0+2]+0.5*uOther[a0+1]*m1rOther[a0+2]; 
    relKinE[4] += 0.4472135954999579*m1rSelf[a0+7]*uSelf[a0+7]-0.4472135954999579*m1rOther[a0+7]*uSelf[a0+7]-0.4472135954999579*m1rSelf[a0+7]*uOther[a0+7]+0.4472135954999579*m1rOther[a0+7]*uOther[a0+7]+0.31943828249997*m1rSelf[a0+6]*uSelf[a0+6]-0.31943828249997*m1rOther[a0+6]*uSelf[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+6]-0.5000000000000001*m1rOther[a0+2]*uSelf[a0+6]-0.31943828249997*m1rSelf[a0+6]*uOther[a0+6]+0.31943828249997*m1rOther[a0+6]*uOther[a0+6]-0.5000000000000001*m1rSelf[a0+2]*uOther[a0+6]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+6]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+6]-0.5000000000000001*uOther[a0+2]*m1rSelf[a0+6]-0.5000000000000001*uSelf[a0+2]*m1rOther[a0+6]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+4]-0.31943828249997*m1rOther[a0+4]*uSelf[a0+4]+0.5*m1rSelf[a0]*uSelf[a0+4]-0.5*m1rOther[a0]*uSelf[a0+4]-0.31943828249997*m1rSelf[a0+4]*uOther[a0+4]+0.31943828249997*m1rOther[a0+4]*uOther[a0+4]-0.5*m1rSelf[a0]*uOther[a0+4]+0.5*m1rOther[a0]*uOther[a0+4]+0.5*uSelf[a0]*m1rSelf[a0+4]-0.5*uOther[a0]*m1rSelf[a0+4]-0.5*uSelf[a0]*m1rOther[a0+4]+0.5*uOther[a0]*m1rOther[a0+4]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rSelf[a0+1]*uSelf[a0+1]-0.4472135954999579*m1rOther[a0+1]*uSelf[a0+1]-0.4472135954999579*m1rSelf[a0+1]*uOther[a0+1]+0.4472135954999579*m1rOther[a0+1]*uOther[a0+1]; 
    relKinE[5] += 0.31943828249997*m1rSelf[a0+7]*uSelf[a0+7]-0.31943828249997*m1rOther[a0+7]*uSelf[a0+7]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+7]-0.5000000000000001*m1rOther[a0+1]*uSelf[a0+7]-0.31943828249997*m1rSelf[a0+7]*uOther[a0+7]+0.31943828249997*m1rOther[a0+7]*uOther[a0+7]-0.5000000000000001*m1rSelf[a0+1]*uOther[a0+7]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+7]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+7]-0.5000000000000001*uOther[a0+1]*m1rSelf[a0+7]-0.5000000000000001*uSelf[a0+1]*m1rOther[a0+7]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+7]+0.4472135954999579*m1rSelf[a0+6]*uSelf[a0+6]-0.4472135954999579*m1rOther[a0+6]*uSelf[a0+6]-0.4472135954999579*m1rSelf[a0+6]*uOther[a0+6]+0.4472135954999579*m1rOther[a0+6]*uOther[a0+6]+0.31943828249997*m1rSelf[a0+5]*uSelf[a0+5]-0.31943828249997*m1rOther[a0+5]*uSelf[a0+5]+0.5*m1rSelf[a0]*uSelf[a0+5]-0.5*m1rOther[a0]*uSelf[a0+5]-0.31943828249997*m1rSelf[a0+5]*uOther[a0+5]+0.31943828249997*m1rOther[a0+5]*uOther[a0+5]-0.5*m1rSelf[a0]*uOther[a0+5]+0.5*m1rOther[a0]*uOther[a0+5]+0.5*uSelf[a0]*m1rSelf[a0+5]-0.5*uOther[a0]*m1rSelf[a0+5]-0.5*uSelf[a0]*m1rOther[a0+5]+0.5*uOther[a0]*m1rOther[a0+5]+0.4472135954999579*m1rSelf[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rOther[a0+3]*uSelf[a0+3]-0.4472135954999579*m1rSelf[a0+3]*uOther[a0+3]+0.4472135954999579*m1rOther[a0+3]*uOther[a0+3]+0.4472135954999579*m1rSelf[a0+2]*uSelf[a0+2]-0.4472135954999579*m1rOther[a0+2]*uSelf[a0+2]-0.4472135954999579*m1rSelf[a0+2]*uOther[a0+2]+0.4472135954999579*m1rOther[a0+2]*uOther[a0+2]; 
    relKinE[6] += 0.4*m1rSelf[a0+3]*uSelf[a0+7]-0.4*m1rOther[a0+3]*uSelf[a0+7]-0.4*m1rSelf[a0+3]*uOther[a0+7]+0.4*m1rOther[a0+3]*uOther[a0+7]+0.4*uSelf[a0+3]*m1rSelf[a0+7]-0.4*uOther[a0+3]*m1rSelf[a0+7]-0.4*uSelf[a0+3]*m1rOther[a0+7]+0.4*uOther[a0+3]*m1rOther[a0+7]+0.4472135954999579*m1rSelf[a0+5]*uSelf[a0+6]-0.4472135954999579*m1rOther[a0+5]*uSelf[a0+6]+0.31943828249997*m1rSelf[a0+4]*uSelf[a0+6]-0.31943828249997*m1rOther[a0+4]*uSelf[a0+6]+0.5*m1rSelf[a0]*uSelf[a0+6]-0.5*m1rOther[a0]*uSelf[a0+6]-0.4472135954999579*m1rSelf[a0+5]*uOther[a0+6]+0.4472135954999579*m1rOther[a0+5]*uOther[a0+6]-0.31943828249997*m1rSelf[a0+4]*uOther[a0+6]+0.31943828249997*m1rOther[a0+4]*uOther[a0+6]-0.5*m1rSelf[a0]*uOther[a0+6]+0.5*m1rOther[a0]*uOther[a0+6]+0.4472135954999579*uSelf[a0+5]*m1rSelf[a0+6]-0.4472135954999579*uOther[a0+5]*m1rSelf[a0+6]+0.31943828249997*uSelf[a0+4]*m1rSelf[a0+6]-0.31943828249997*uOther[a0+4]*m1rSelf[a0+6]+0.5*uSelf[a0]*m1rSelf[a0+6]-0.5*uOther[a0]*m1rSelf[a0+6]-0.4472135954999579*uSelf[a0+5]*m1rOther[a0+6]+0.4472135954999579*uOther[a0+5]*m1rOther[a0+6]-0.31943828249997*uSelf[a0+4]*m1rOther[a0+6]+0.31943828249997*uOther[a0+4]*m1rOther[a0+6]-0.5*uSelf[a0]*m1rOther[a0+6]+0.5*uOther[a0]*m1rOther[a0+6]+0.5000000000000001*m1rSelf[a0+2]*uSelf[a0+4]-0.5000000000000001*m1rOther[a0+2]*uSelf[a0+4]-0.5000000000000001*m1rSelf[a0+2]*uOther[a0+4]+0.5000000000000001*m1rOther[a0+2]*uOther[a0+4]+0.5000000000000001*uSelf[a0+2]*m1rSelf[a0+4]-0.5000000000000001*uOther[a0+2]*m1rSelf[a0+4]-0.5000000000000001*uSelf[a0+2]*m1rOther[a0+4]+0.5000000000000001*uOther[a0+2]*m1rOther[a0+4]+0.447213595499958*m1rSelf[a0+1]*uSelf[a0+3]-0.447213595499958*m1rOther[a0+1]*uSelf[a0+3]-0.447213595499958*m1rSelf[a0+1]*uOther[a0+3]+0.447213595499958*m1rOther[a0+1]*uOther[a0+3]+0.447213595499958*uSelf[a0+1]*m1rSelf[a0+3]-0.447213595499958*uOther[a0+1]*m1rSelf[a0+3]-0.447213595499958*uSelf[a0+1]*m1rOther[a0+3]+0.447213595499958*uOther[a0+1]*m1rOther[a0+3]; 
    relKinE[7] += 0.31943828249997*m1rSelf[a0+5]*uSelf[a0+7]-0.31943828249997*m1rOther[a0+5]*uSelf[a0+7]+0.4472135954999579*m1rSelf[a0+4]*uSelf[a0+7]-0.4472135954999579*m1rOther[a0+4]*uSelf[a0+7]+0.5*m1rSelf[a0]*uSelf[a0+7]-0.5*m1rOther[a0]*uSelf[a0+7]-0.31943828249997*m1rSelf[a0+5]*uOther[a0+7]+0.31943828249997*m1rOther[a0+5]*uOther[a0+7]-0.4472135954999579*m1rSelf[a0+4]*uOther[a0+7]+0.4472135954999579*m1rOther[a0+4]*uOther[a0+7]-0.5*m1rSelf[a0]*uOther[a0+7]+0.5*m1rOther[a0]*uOther[a0+7]+0.31943828249997*uSelf[a0+5]*m1rSelf[a0+7]-0.31943828249997*uOther[a0+5]*m1rSelf[a0+7]+0.4472135954999579*uSelf[a0+4]*m1rSelf[a0+7]-0.4472135954999579*uOther[a0+4]*m1rSelf[a0+7]+0.5*uSelf[a0]*m1rSelf[a0+7]-0.5*uOther[a0]*m1rSelf[a0+7]-0.31943828249997*uSelf[a0+5]*m1rOther[a0+7]+0.31943828249997*uOther[a0+5]*m1rOther[a0+7]-0.4472135954999579*uSelf[a0+4]*m1rOther[a0+7]+0.4472135954999579*uOther[a0+4]*m1rOther[a0+7]-0.5*uSelf[a0]*m1rOther[a0+7]+0.5*uOther[a0]*m1rOther[a0+7]+0.4*m1rSelf[a0+3]*uSelf[a0+6]-0.4*m1rOther[a0+3]*uSelf[a0+6]-0.4*m1rSelf[a0+3]*uOther[a0+6]+0.4*m1rOther[a0+3]*uOther[a0+6]+0.4*uSelf[a0+3]*m1rSelf[a0+6]-0.4*uOther[a0+3]*m1rSelf[a0+6]-0.4*uSelf[a0+3]*m1rOther[a0+6]+0.4*uOther[a0+3]*m1rOther[a0+6]+0.5000000000000001*m1rSelf[a0+1]*uSelf[a0+5]-0.5000000000000001*m1rOther[a0+1]*uSelf[a0+5]-0.5000000000000001*m1rSelf[a0+1]*uOther[a0+5]+0.5000000000000001*m1rOther[a0+1]*uOther[a0+5]+0.5000000000000001*uSelf[a0+1]*m1rSelf[a0+5]-0.5000000000000001*uOther[a0+1]*m1rSelf[a0+5]-0.5000000000000001*uSelf[a0+1]*m1rOther[a0+5]+0.5000000000000001*uOther[a0+1]*m1rOther[a0+5]+0.447213595499958*m1rSelf[a0+2]*uSelf[a0+3]-0.447213595499958*m1rOther[a0+2]*uSelf[a0+3]-0.447213595499958*m1rSelf[a0+2]*uOther[a0+3]+0.447213595499958*m1rOther[a0+2]*uOther[a0+3]+0.447213595499958*uSelf[a0+2]*m1rSelf[a0+3]-0.447213595499958*uOther[a0+2]*m1rSelf[a0+3]-0.447213595499958*uSelf[a0+2]*m1rOther[a0+3]+0.447213595499958*uOther[a0+2]*m1rOther[a0+3]; 
  } 
 
  double m2Relax[8]; 
  // ... Contribution to RHS vector from energy relaxation. 
  m2Relax[0] = betaGreenep1*((-(0.5*relKinE[0]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[0]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[0]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[0]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[0]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[0]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[0]-1.0*kinESelf[0])*mnuSelf+(kinEOther[0]-1.0*m2rOther[0])*mnuOther; 
  m2Relax[1] = betaGreenep1*((-(0.5*relKinE[1]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[1]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[1]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[1]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[1]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[1]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[1]-1.0*kinESelf[1])*mnuSelf+(kinEOther[1]-1.0*m2rOther[1])*mnuOther; 
  m2Relax[2] = betaGreenep1*((-(0.5*relKinE[2]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[2]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[2]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[2]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[2]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[2]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[2]-1.0*kinESelf[2])*mnuSelf+(kinEOther[2]-1.0*m2rOther[2])*mnuOther; 
  m2Relax[3] = betaGreenep1*((-(0.5*relKinE[3]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[3]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[3]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[3]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[3]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[3]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[3]-1.0*kinESelf[3])*mnuSelf+(kinEOther[3]-1.0*m2rOther[3])*mnuOther; 
  m2Relax[4] = betaGreenep1*((-(0.5*relKinE[4]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[4]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[4]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[4]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[4]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[4]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[4]-1.0*kinESelf[4])*mnuSelf+(kinEOther[4]-1.0*m2rOther[4])*mnuOther; 
  m2Relax[5] = betaGreenep1*((-(0.5*relKinE[5]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[5]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[5]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[5]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[5]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[5]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[5]-1.0*kinESelf[5])*mnuSelf+(kinEOther[5]-1.0*m2rOther[5])*mnuOther; 
  m2Relax[6] = betaGreenep1*((-(0.5*relKinE[6]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[6]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[6]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[6]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[6]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[6]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[6]-1.0*kinESelf[6])*mnuSelf+(kinEOther[6]-1.0*m2rOther[6])*mnuOther; 
  m2Relax[7] = betaGreenep1*((-(0.5*relKinE[7]*deltaSelf*mSelf)/(mSelf+mOther))-(1.0*m2rSelf[7]*deltaSelf*mSelf)/(mSelf+mOther)+(kinESelf[7]*deltaSelf*mSelf)/(mSelf+mOther)+(0.5*relKinE[7]*deltaSelf*mOther)/(mSelf+mOther)+(m2rOther[7]*deltaSelf*mOther)/(mSelf+mOther)-(1.0*kinEOther[7]*deltaSelf*mOther)/(mSelf+mOther))*mnuSelf+(m2rSelf[7]-1.0*kinESelf[7])*mnuSelf+(kinEOther[7]-1.0*m2rOther[7])*mnuOther; 
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3],mnuM2sum[4],mnuM2sum[5],mnuM2sum[6],mnuM2sum[7],m1Relax[0],m1Relax[1],m1Relax[2],m1Relax[3],m1Relax[4],m1Relax[5],m1Relax[6],m1Relax[7],m2Relax[0],m2Relax[1],m2Relax[2],m2Relax[3],m2Relax[4],m2Relax[5],m2Relax[6],m2Relax[7]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,8,1) = data->u_S.segment<8>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,8,1) = data->u_S.segment<8>(8); 
 
  Eigen::Map<VectorXd>(uCrossOther,8,1) = data->u_S.segment<8>(16); 
 
  Eigen::Map<VectorXd>(vtSqCrossOther,8,1) = data->u_S.segment<8>(24); 
 
} 
 
