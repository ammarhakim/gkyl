#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmLBOGCrossPrimMoments2x3vSer_P1(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  if (1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if (1.5*m2SSelf[3]-0.8660254037844386*m2SSelf[2]-0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
  if ((-1.5*m0Self[3])-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if ((-1.5*m2SSelf[3])-0.8660254037844386*m2SSelf[2]+0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
  if ((-1.5*m0Self[3])+0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if ((-1.5*m2SSelf[3])+0.8660254037844386*m2SSelf[2]-0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
  if (1.5*m0Self[3]+0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if (1.5*m2SSelf[3]+0.8660254037844386*m2SSelf[2]+0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
 
  double m0rSelf[4]; 
  double m1rSelf[12]; 
  double m2rSelf[4]; 
  double cMrSelf[12]; 
  double cErSelf[4]; 
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
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = 0.0; 
    cMrSelf[2] = 0.0; 
    cMrSelf[3] = 0.0; 
    m1rSelf[4] = m1Self[4]; 
    m1rSelf[5] = 0.0; 
    m1rSelf[6] = 0.0; 
    m1rSelf[7] = 0.0; 
    cMrSelf[4] = cMSelf[4]; 
    cMrSelf[5] = 0.0; 
    cMrSelf[6] = 0.0; 
    cMrSelf[7] = 0.0; 
    m1rSelf[8] = m1Self[8]; 
    m1rSelf[9] = 0.0; 
    m1rSelf[10] = 0.0; 
    m1rSelf[11] = 0.0; 
    cMrSelf[8] = cMSelf[8]; 
    cMrSelf[9] = 0.0; 
    cMrSelf[10] = 0.0; 
    cMrSelf[11] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    m2rSelf[3] = 0.0; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = 0.0; 
    cErSelf[2] = 0.0; 
    cErSelf[3] = 0.0; 
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
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = cMSelf[1]; 
    cMrSelf[2] = cMSelf[2]; 
    cMrSelf[3] = cMSelf[3]; 
    cMrSelf[4] = cMSelf[4]; 
    cMrSelf[5] = cMSelf[5]; 
    cMrSelf[6] = cMSelf[6]; 
    cMrSelf[7] = cMSelf[7]; 
    cMrSelf[8] = cMSelf[8]; 
    cMrSelf[9] = cMSelf[9]; 
    cMrSelf[10] = cMSelf[10]; 
    cMrSelf[11] = cMSelf[11]; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = cESelf[1]; 
    cErSelf[2] = cESelf[2]; 
    cErSelf[3] = cESelf[3]; 
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
 
  cellAvg = false;
  if (1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if (1.5*m2SOther[3]-0.8660254037844386*m2SOther[2]-0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
  if ((-1.5*m0Other[3])-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if ((-1.5*m2SOther[3])-0.8660254037844386*m2SOther[2]+0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
  if ((-1.5*m0Other[3])+0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if ((-1.5*m2SOther[3])+0.8660254037844386*m2SOther[2]-0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
  if (1.5*m0Other[3]+0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if (1.5*m2SOther[3]+0.8660254037844386*m2SOther[2]+0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
 
  double m0rOther[4]; 
  double m1rOther[12]; 
  double m2rOther[4]; 
  double cMrOther[12]; 
  double cErOther[4]; 
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
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = 0.0; 
    cMrOther[2] = 0.0; 
    cMrOther[3] = 0.0; 
    m1rOther[4] = m1Other[4]; 
    m1rOther[5] = 0.0; 
    m1rOther[6] = 0.0; 
    m1rOther[7] = 0.0; 
    cMrOther[4] = cMOther[4]; 
    cMrOther[5] = 0.0; 
    cMrOther[6] = 0.0; 
    cMrOther[7] = 0.0; 
    m1rOther[8] = m1Other[8]; 
    m1rOther[9] = 0.0; 
    m1rOther[10] = 0.0; 
    m1rOther[11] = 0.0; 
    cMrOther[8] = cMOther[8]; 
    cMrOther[9] = 0.0; 
    cMrOther[10] = 0.0; 
    cMrOther[11] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    m2rOther[3] = 0.0; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = 0.0; 
    cErOther[2] = 0.0; 
    cErOther[3] = 0.0; 
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
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = cMOther[1]; 
    cMrOther[2] = cMOther[2]; 
    cMrOther[3] = cMOther[3]; 
    cMrOther[4] = cMOther[4]; 
    cMrOther[5] = cMOther[5]; 
    cMrOther[6] = cMOther[6]; 
    cMrOther[7] = cMOther[7]; 
    cMrOther[8] = cMOther[8]; 
    cMrOther[9] = cMOther[9]; 
    cMrOther[10] = cMOther[10]; 
    cMrOther[11] = cMOther[11]; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = cEOther[1]; 
    cErOther[2] = cEOther[2]; 
    cErOther[3] = cEOther[3]; 
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
  data->AEM_S.setZero(); 
 
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
  data->AEM_S(0,12) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(0,13) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(0,14) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(0,15) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(1,12) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(1,13) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(1,14) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(1,15) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(2,12) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(2,13) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(2,14) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(2,15) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(3,12) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(3,13) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(3,14) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(3,15) = -0.5*cMrSelf[0]*mnuSelf; 
 
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
  data->AEM_S(0,28) = -0.5*cMrOther[0]*mnuOther; 
  data->AEM_S(0,29) = -0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(0,30) = -0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(0,31) = -0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(1,28) = -0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(1,29) = -0.5*cMrOther[0]*mnuOther; 
  data->AEM_S(1,30) = -0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(1,31) = -0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(2,28) = -0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(2,29) = -0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(2,30) = -0.5*cMrOther[0]*mnuOther; 
  data->AEM_S(2,31) = -0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(3,28) = -0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(3,29) = -0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(3,30) = -0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(3,31) = -0.5*cMrOther[0]*mnuOther; 
 
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
  data->AEM_S(4,12) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(4,13) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(4,14) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(4,15) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(5,12) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(5,13) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(5,14) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(5,15) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(6,12) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(6,13) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(6,14) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(6,15) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(7,12) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(7,13) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(7,14) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(7,15) = -0.5*cMrSelf[4]*mnuSelf; 
 
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
  data->AEM_S(4,28) = -0.5*cMrOther[4]*mnuOther; 
  data->AEM_S(4,29) = -0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(4,30) = -0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(4,31) = -0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(5,28) = -0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(5,29) = -0.5*cMrOther[4]*mnuOther; 
  data->AEM_S(5,30) = -0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(5,31) = -0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(6,28) = -0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(6,29) = -0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(6,30) = -0.5*cMrOther[4]*mnuOther; 
  data->AEM_S(6,31) = -0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(7,28) = -0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(7,29) = -0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(7,30) = -0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(7,31) = -0.5*cMrOther[4]*mnuOther; 
 
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
  data->AEM_S(8,12) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(8,13) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(8,14) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(8,15) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(9,12) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(9,13) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(9,14) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(9,15) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(10,12) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(10,13) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(10,14) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(10,15) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(11,12) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(11,13) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(11,14) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(11,15) = -0.5*cMrSelf[8]*mnuSelf; 
 
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
  data->AEM_S(8,28) = -0.5*cMrOther[8]*mnuOther; 
  data->AEM_S(8,29) = -0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(8,30) = -0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(8,31) = -0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(9,28) = -0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(9,29) = -0.5*cMrOther[8]*mnuOther; 
  data->AEM_S(9,30) = -0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(9,31) = -0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(10,28) = -0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(10,29) = -0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(10,30) = -0.5*cMrOther[8]*mnuOther; 
  data->AEM_S(10,31) = -0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(11,28) = -0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(11,29) = -0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(11,30) = -0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(11,31) = -0.5*cMrOther[8]*mnuOther; 
 
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
  data->AEM_S(12,12) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(12,14) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(12,15) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(13,12) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(13,15) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(14,12) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(14,13) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(14,14) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(14,15) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(15,12) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(15,13) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(15,14) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(15,15) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(12,28) = 0.5*m0SrOther[0]*mnuOther-0.5*cErOther[0]*mnuOther; 
  data->AEM_S(12,29) = 0.5*m0SrOther[1]*mnuOther-0.5*cErOther[1]*mnuOther; 
  data->AEM_S(12,30) = 0.5*m0SrOther[2]*mnuOther-0.5*cErOther[2]*mnuOther; 
  data->AEM_S(12,31) = 0.5*m0SrOther[3]*mnuOther-0.5*cErOther[3]*mnuOther; 
  data->AEM_S(13,28) = 0.5*m0SrOther[1]*mnuOther-0.5*cErOther[1]*mnuOther; 
  data->AEM_S(13,29) = 0.5*m0SrOther[0]*mnuOther-0.5*cErOther[0]*mnuOther; 
  data->AEM_S(13,30) = 0.5*m0SrOther[3]*mnuOther-0.5*cErOther[3]*mnuOther; 
  data->AEM_S(13,31) = 0.5*m0SrOther[2]*mnuOther-0.5*cErOther[2]*mnuOther; 
  data->AEM_S(14,28) = 0.5*m0SrOther[2]*mnuOther-0.5*cErOther[2]*mnuOther; 
  data->AEM_S(14,29) = 0.5*m0SrOther[3]*mnuOther-0.5*cErOther[3]*mnuOther; 
  data->AEM_S(14,30) = 0.5*m0SrOther[0]*mnuOther-0.5*cErOther[0]*mnuOther; 
  data->AEM_S(14,31) = 0.5*m0SrOther[1]*mnuOther-0.5*cErOther[1]*mnuOther; 
  data->AEM_S(15,28) = 0.5*m0SrOther[3]*mnuOther-0.5*cErOther[3]*mnuOther; 
  data->AEM_S(15,29) = 0.5*m0SrOther[2]*mnuOther-0.5*cErOther[2]*mnuOther; 
  data->AEM_S(15,30) = 0.5*m0SrOther[1]*mnuOther-0.5*cErOther[1]*mnuOther; 
  data->AEM_S(15,31) = 0.5*m0SrOther[0]*mnuOther-0.5*cErOther[0]*mnuOther; 
 
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
  data->AEM_S(16,12) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(16,13) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(16,14) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(16,15) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(17,12) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(17,13) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(17,14) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(17,15) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(18,12) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(18,13) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(18,14) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(18,15) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(19,12) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(19,13) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(19,14) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(19,15) = -0.5*cMrSelf[0]*mnuSelf; 
 
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
  data->AEM_S(16,28) = 0.5*cMrOther[0]*mnuOther; 
  data->AEM_S(16,29) = 0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(16,30) = 0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(16,31) = 0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(17,28) = 0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(17,29) = 0.5*cMrOther[0]*mnuOther; 
  data->AEM_S(17,30) = 0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(17,31) = 0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(18,28) = 0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(18,29) = 0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(18,30) = 0.5*cMrOther[0]*mnuOther; 
  data->AEM_S(18,31) = 0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(19,28) = 0.5*cMrOther[3]*mnuOther; 
  data->AEM_S(19,29) = 0.5*cMrOther[2]*mnuOther; 
  data->AEM_S(19,30) = 0.5*cMrOther[1]*mnuOther; 
  data->AEM_S(19,31) = 0.5*cMrOther[0]*mnuOther; 
 
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
  data->AEM_S(20,12) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(20,13) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(20,14) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(20,15) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(21,12) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(21,13) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(21,14) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(21,15) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(22,12) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(22,13) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(22,14) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(22,15) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(23,12) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(23,13) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(23,14) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(23,15) = -0.5*cMrSelf[4]*mnuSelf; 
 
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
  data->AEM_S(20,28) = 0.5*cMrOther[4]*mnuOther; 
  data->AEM_S(20,29) = 0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(20,30) = 0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(20,31) = 0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(21,28) = 0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(21,29) = 0.5*cMrOther[4]*mnuOther; 
  data->AEM_S(21,30) = 0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(21,31) = 0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(22,28) = 0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(22,29) = 0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(22,30) = 0.5*cMrOther[4]*mnuOther; 
  data->AEM_S(22,31) = 0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(23,28) = 0.5*cMrOther[7]*mnuOther; 
  data->AEM_S(23,29) = 0.5*cMrOther[6]*mnuOther; 
  data->AEM_S(23,30) = 0.5*cMrOther[5]*mnuOther; 
  data->AEM_S(23,31) = 0.5*cMrOther[4]*mnuOther; 
 
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
  data->AEM_S(24,12) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(24,13) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(24,14) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(24,15) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(25,12) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(25,13) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(25,14) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(25,15) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(26,12) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(26,13) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(26,14) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(26,15) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(27,12) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(27,13) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(27,14) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(27,15) = -0.5*cMrSelf[8]*mnuSelf; 
 
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
  data->AEM_S(24,28) = 0.5*cMrOther[8]*mnuOther; 
  data->AEM_S(24,29) = 0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(24,30) = 0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(24,31) = 0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(25,28) = 0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(25,29) = 0.5*cMrOther[8]*mnuOther; 
  data->AEM_S(25,30) = 0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(25,31) = 0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(26,28) = 0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(26,29) = 0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(26,30) = 0.5*cMrOther[8]*mnuOther; 
  data->AEM_S(26,31) = 0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(27,28) = 0.5*cMrOther[11]*mnuOther; 
  data->AEM_S(27,29) = 0.5*cMrOther[10]*mnuOther; 
  data->AEM_S(27,30) = 0.5*cMrOther[9]*mnuOther; 
  data->AEM_S(27,31) = 0.5*cMrOther[8]*mnuOther; 
 
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
    ucMSelf[0] += 0.5*cMrSelf[a0+3]*uSelf[a0+3]+0.5*cMrSelf[a0+2]*uSelf[a0+2]+0.5*cMrSelf[a0+1]*uSelf[a0+1]+0.5*cMrSelf[a0]*uSelf[a0]; 
    ucMSelf[1] += 0.5*cMrSelf[a0+2]*uSelf[a0+3]+0.5*uSelf[a0+2]*cMrSelf[a0+3]+0.5*cMrSelf[a0]*uSelf[a0+1]+0.5*uSelf[a0]*cMrSelf[a0+1]; 
    ucMSelf[2] += 0.5*cMrSelf[a0+1]*uSelf[a0+3]+0.5*uSelf[a0+1]*cMrSelf[a0+3]+0.5*cMrSelf[a0]*uSelf[a0+2]+0.5*uSelf[a0]*cMrSelf[a0+2]; 
    ucMSelf[3] += 0.5*cMrSelf[a0]*uSelf[a0+3]+0.5*uSelf[a0]*cMrSelf[a0+3]+0.5*cMrSelf[a0+1]*uSelf[a0+2]+0.5*uSelf[a0+1]*cMrSelf[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  data->AEM_S(28,12) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(28,13) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(28,14) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(28,15) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(29,12) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(29,13) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(29,14) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(29,15) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(30,12) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(30,13) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(30,14) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(30,15) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(31,12) = 0.5*ucMSelf[3]*mnuSelf+0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(31,13) = 0.5*ucMSelf[2]*mnuSelf+0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(31,14) = 0.5*ucMSelf[1]*mnuSelf+0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(31,15) = 0.5*ucMSelf[0]*mnuSelf+0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
 
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
    ucMOther[0] += 0.5*cMrOther[a0+3]*uOther[a0+3]+0.5*cMrOther[a0+2]*uOther[a0+2]+0.5*cMrOther[a0+1]*uOther[a0+1]+0.5*cMrOther[a0]*uOther[a0]; 
    ucMOther[1] += 0.5*cMrOther[a0+2]*uOther[a0+3]+0.5*uOther[a0+2]*cMrOther[a0+3]+0.5*cMrOther[a0]*uOther[a0+1]+0.5*uOther[a0]*cMrOther[a0+1]; 
    ucMOther[2] += 0.5*cMrOther[a0+1]*uOther[a0+3]+0.5*uOther[a0+1]*cMrOther[a0+3]+0.5*cMrOther[a0]*uOther[a0+2]+0.5*uOther[a0]*cMrOther[a0+2]; 
    ucMOther[3] += 0.5*cMrOther[a0]*uOther[a0+3]+0.5*uOther[a0]*cMrOther[a0+3]+0.5*cMrOther[a0+1]*uOther[a0+2]+0.5*uOther[a0+1]*cMrOther[a0+2]; 
  } 
 
  // ... Block from correction to (other) 2nd moment of collision operator ... // 
  data->AEM_S(28,28) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cErOther[0]*mnuOther; 
  data->AEM_S(28,29) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cErOther[1]*mnuOther; 
  data->AEM_S(28,30) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cErOther[2]*mnuOther; 
  data->AEM_S(28,31) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cErOther[3]*mnuOther; 
  data->AEM_S(29,28) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cErOther[1]*mnuOther; 
  data->AEM_S(29,29) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cErOther[0]*mnuOther; 
  data->AEM_S(29,30) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cErOther[3]*mnuOther; 
  data->AEM_S(29,31) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cErOther[2]*mnuOther; 
  data->AEM_S(30,28) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cErOther[2]*mnuOther; 
  data->AEM_S(30,29) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cErOther[3]*mnuOther; 
  data->AEM_S(30,30) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cErOther[0]*mnuOther; 
  data->AEM_S(30,31) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cErOther[1]*mnuOther; 
  data->AEM_S(31,28) = (-0.5*ucMOther[3]*mnuOther)-0.5*m0SrOther[3]*mnuOther+0.5*cErOther[3]*mnuOther; 
  data->AEM_S(31,29) = (-0.5*ucMOther[2]*mnuOther)-0.5*m0SrOther[2]*mnuOther+0.5*cErOther[2]*mnuOther; 
  data->AEM_S(31,30) = (-0.5*ucMOther[1]*mnuOther)-0.5*m0SrOther[1]*mnuOther+0.5*cErOther[1]*mnuOther; 
  data->AEM_S(31,31) = (-0.5*ucMOther[0]*mnuOther)-0.5*m0SrOther[0]*mnuOther+0.5*cErOther[0]*mnuOther; 
 
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
 
void VmLBOECrossPrimMoments2x3vSer_P1(binOpData_t *data, binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *m1Self, const double *m2Self, const double *uSelf, const double *vtSqSelf, const double *cMSelf, const double *cESelf, const double *m0SSelf, const double *m1SSelf, const double *m2SSelf, const double mOther, const double nuOther, const double *m0Other, const double *m1Other, const double *m2Other, const double *uOther, const double *vtSqOther, const double *cMOther, const double *cEOther, const double *m0SOther, const double *m1SOther, const double *m2SOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
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
  if (1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if (1.5*m2SSelf[3]-0.8660254037844386*m2SSelf[2]-0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
  if ((-1.5*m0Self[3])-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if ((-1.5*m2SSelf[3])-0.8660254037844386*m2SSelf[2]+0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
  if ((-1.5*m0Self[3])+0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if ((-1.5*m2SSelf[3])+0.8660254037844386*m2SSelf[2]-0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
  if (1.5*m0Self[3]+0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) cellAvg = true; 
  if (1.5*m2SSelf[3]+0.8660254037844386*m2SSelf[2]+0.8660254037844386*m2SSelf[1]+0.5*m2SSelf[0] < 0) cellAvg = true; 
 
  double m0rSelf[4]; 
  double m1rSelf[12]; 
  double m2rSelf[4]; 
  double cMrSelf[12]; 
  double cErSelf[4]; 
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
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = 0.0; 
    cMrSelf[2] = 0.0; 
    cMrSelf[3] = 0.0; 
    m1rSelf[4] = m1Self[4]; 
    m1rSelf[5] = 0.0; 
    m1rSelf[6] = 0.0; 
    m1rSelf[7] = 0.0; 
    cMrSelf[4] = cMSelf[4]; 
    cMrSelf[5] = 0.0; 
    cMrSelf[6] = 0.0; 
    cMrSelf[7] = 0.0; 
    m1rSelf[8] = m1Self[8]; 
    m1rSelf[9] = 0.0; 
    m1rSelf[10] = 0.0; 
    m1rSelf[11] = 0.0; 
    cMrSelf[8] = cMSelf[8]; 
    cMrSelf[9] = 0.0; 
    cMrSelf[10] = 0.0; 
    cMrSelf[11] = 0.0; 
    m2rSelf[0] = m2Self[0]; 
    m2rSelf[1] = 0.0; 
    m2rSelf[2] = 0.0; 
    m2rSelf[3] = 0.0; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = 0.0; 
    cErSelf[2] = 0.0; 
    cErSelf[3] = 0.0; 
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
    cMrSelf[0] = cMSelf[0]; 
    cMrSelf[1] = cMSelf[1]; 
    cMrSelf[2] = cMSelf[2]; 
    cMrSelf[3] = cMSelf[3]; 
    cMrSelf[4] = cMSelf[4]; 
    cMrSelf[5] = cMSelf[5]; 
    cMrSelf[6] = cMSelf[6]; 
    cMrSelf[7] = cMSelf[7]; 
    cMrSelf[8] = cMSelf[8]; 
    cMrSelf[9] = cMSelf[9]; 
    cMrSelf[10] = cMSelf[10]; 
    cMrSelf[11] = cMSelf[11]; 
    cErSelf[0] = cESelf[0]; 
    cErSelf[1] = cESelf[1]; 
    cErSelf[2] = cESelf[2]; 
    cErSelf[3] = cESelf[3]; 
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
 
  cellAvg = false;
  if (1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if (1.5*m2SOther[3]-0.8660254037844386*m2SOther[2]-0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
  if ((-1.5*m0Other[3])-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if ((-1.5*m2SOther[3])-0.8660254037844386*m2SOther[2]+0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
  if ((-1.5*m0Other[3])+0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if ((-1.5*m2SOther[3])+0.8660254037844386*m2SOther[2]-0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
  if (1.5*m0Other[3]+0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) cellAvg = true; 
  if (1.5*m2SOther[3]+0.8660254037844386*m2SOther[2]+0.8660254037844386*m2SOther[1]+0.5*m2SOther[0] < 0) cellAvg = true; 
 
  double m0rOther[4]; 
  double m1rOther[12]; 
  double m2rOther[4]; 
  double cMrOther[12]; 
  double cErOther[4]; 
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
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = 0.0; 
    cMrOther[2] = 0.0; 
    cMrOther[3] = 0.0; 
    m1rOther[4] = m1Other[4]; 
    m1rOther[5] = 0.0; 
    m1rOther[6] = 0.0; 
    m1rOther[7] = 0.0; 
    cMrOther[4] = cMOther[4]; 
    cMrOther[5] = 0.0; 
    cMrOther[6] = 0.0; 
    cMrOther[7] = 0.0; 
    m1rOther[8] = m1Other[8]; 
    m1rOther[9] = 0.0; 
    m1rOther[10] = 0.0; 
    m1rOther[11] = 0.0; 
    cMrOther[8] = cMOther[8]; 
    cMrOther[9] = 0.0; 
    cMrOther[10] = 0.0; 
    cMrOther[11] = 0.0; 
    m2rOther[0] = m2Other[0]; 
    m2rOther[1] = 0.0; 
    m2rOther[2] = 0.0; 
    m2rOther[3] = 0.0; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = 0.0; 
    cErOther[2] = 0.0; 
    cErOther[3] = 0.0; 
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
    cMrOther[0] = cMOther[0]; 
    cMrOther[1] = cMOther[1]; 
    cMrOther[2] = cMOther[2]; 
    cMrOther[3] = cMOther[3]; 
    cMrOther[4] = cMOther[4]; 
    cMrOther[5] = cMOther[5]; 
    cMrOther[6] = cMOther[6]; 
    cMrOther[7] = cMOther[7]; 
    cMrOther[8] = cMOther[8]; 
    cMrOther[9] = cMOther[9]; 
    cMrOther[10] = cMOther[10]; 
    cMrOther[11] = cMOther[11]; 
    cErOther[0] = cEOther[0]; 
    cErOther[1] = cEOther[1]; 
    cErOther[2] = cEOther[2]; 
    cErOther[3] = cEOther[3]; 
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
  data->AEM_S.setZero(); 
 
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
  data->AEM_S(0,12) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(0,13) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(0,14) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(0,15) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(1,12) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(1,13) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(1,14) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(1,15) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(2,12) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(2,13) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(2,14) = -0.5*cMrSelf[0]*mnuSelf; 
  data->AEM_S(2,15) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(3,12) = -0.5*cMrSelf[3]*mnuSelf; 
  data->AEM_S(3,13) = -0.5*cMrSelf[2]*mnuSelf; 
  data->AEM_S(3,14) = -0.5*cMrSelf[1]*mnuSelf; 
  data->AEM_S(3,15) = -0.5*cMrSelf[0]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherX ... // 
  data->AEM_S(0,0) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(0,1) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(0,2) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(0,3) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,0) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(1,1) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(1,2) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(1,3) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,0) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(2,1) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(2,2) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(2,3) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,0) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(3,1) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(3,2) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(3,3) += 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(0,12) += -(0.5*cMrOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(0,13) += -(0.5*cMrOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(0,14) += -(0.5*cMrOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(0,15) += -(0.5*cMrOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(1,12) += -(0.5*cMrOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(1,13) += -(0.5*cMrOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(1,14) += -(0.5*cMrOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(1,15) += -(0.5*cMrOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(2,12) += -(0.5*cMrOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(2,13) += -(0.5*cMrOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(2,14) += -(0.5*cMrOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(2,15) += -(0.5*cMrOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(3,12) += -(0.5*cMrOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(3,13) += -(0.5*cMrOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(3,14) += -(0.5*cMrOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(3,15) += -(0.5*cMrOther[0]*mSelf*mnuOther)/mOther; 
 
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
  data->AEM_S(12,0) += 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(12,1) += 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(12,2) += 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(12,3) += 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(13,0) += 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(13,1) += 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(13,2) += 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(13,3) += 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(14,0) += 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(14,1) += 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(14,2) += 0.5*m1SrOther[0]*mnuOther; 
  data->AEM_S(14,3) += 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(15,0) += 0.5*m1SrOther[3]*mnuOther; 
  data->AEM_S(15,1) += 0.5*m1SrOther[2]*mnuOther; 
  data->AEM_S(15,2) += 0.5*m1SrOther[1]*mnuOther; 
  data->AEM_S(15,3) += 0.5*m1SrOther[0]*mnuOther; 
 
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
  data->AEM_S(4,12) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(4,13) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(4,14) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(4,15) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(5,12) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(5,13) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(5,14) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(5,15) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(6,12) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(6,13) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(6,14) = -0.5*cMrSelf[4]*mnuSelf; 
  data->AEM_S(6,15) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(7,12) = -0.5*cMrSelf[7]*mnuSelf; 
  data->AEM_S(7,13) = -0.5*cMrSelf[6]*mnuSelf; 
  data->AEM_S(7,14) = -0.5*cMrSelf[5]*mnuSelf; 
  data->AEM_S(7,15) = -0.5*cMrSelf[4]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherY ... // 
  data->AEM_S(4,4) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(4,5) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(4,6) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(4,7) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(5,4) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(5,5) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(5,6) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(5,7) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(6,4) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(6,5) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(6,6) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(6,7) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,4) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(7,5) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(7,6) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(7,7) += 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(4,12) += -(0.5*cMrOther[4]*mSelf*mnuOther)/mOther; 
  data->AEM_S(4,13) += -(0.5*cMrOther[5]*mSelf*mnuOther)/mOther; 
  data->AEM_S(4,14) += -(0.5*cMrOther[6]*mSelf*mnuOther)/mOther; 
  data->AEM_S(4,15) += -(0.5*cMrOther[7]*mSelf*mnuOther)/mOther; 
  data->AEM_S(5,12) += -(0.5*cMrOther[5]*mSelf*mnuOther)/mOther; 
  data->AEM_S(5,13) += -(0.5*cMrOther[4]*mSelf*mnuOther)/mOther; 
  data->AEM_S(5,14) += -(0.5*cMrOther[7]*mSelf*mnuOther)/mOther; 
  data->AEM_S(5,15) += -(0.5*cMrOther[6]*mSelf*mnuOther)/mOther; 
  data->AEM_S(6,12) += -(0.5*cMrOther[6]*mSelf*mnuOther)/mOther; 
  data->AEM_S(6,13) += -(0.5*cMrOther[7]*mSelf*mnuOther)/mOther; 
  data->AEM_S(6,14) += -(0.5*cMrOther[4]*mSelf*mnuOther)/mOther; 
  data->AEM_S(6,15) += -(0.5*cMrOther[5]*mSelf*mnuOther)/mOther; 
  data->AEM_S(7,12) += -(0.5*cMrOther[7]*mSelf*mnuOther)/mOther; 
  data->AEM_S(7,13) += -(0.5*cMrOther[6]*mSelf*mnuOther)/mOther; 
  data->AEM_S(7,14) += -(0.5*cMrOther[5]*mSelf*mnuOther)/mOther; 
  data->AEM_S(7,15) += -(0.5*cMrOther[4]*mSelf*mnuOther)/mOther; 
 
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
  data->AEM_S(12,4) += 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(12,5) += 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(12,6) += 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(12,7) += 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(13,4) += 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(13,5) += 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(13,6) += 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(13,7) += 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(14,4) += 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(14,5) += 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(14,6) += 0.5*m1SrOther[4]*mnuOther; 
  data->AEM_S(14,7) += 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(15,4) += 0.5*m1SrOther[7]*mnuOther; 
  data->AEM_S(15,5) += 0.5*m1SrOther[6]*mnuOther; 
  data->AEM_S(15,6) += 0.5*m1SrOther[5]*mnuOther; 
  data->AEM_S(15,7) += 0.5*m1SrOther[4]*mnuOther; 
 
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
  data->AEM_S(8,12) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(8,13) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(8,14) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(8,15) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(9,12) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(9,13) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(9,14) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(9,15) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(10,12) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(10,13) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(10,14) = -0.5*cMrSelf[8]*mnuSelf; 
  data->AEM_S(10,15) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(11,12) = -0.5*cMrSelf[11]*mnuSelf; 
  data->AEM_S(11,13) = -0.5*cMrSelf[10]*mnuSelf; 
  data->AEM_S(11,14) = -0.5*cMrSelf[9]*mnuSelf; 
  data->AEM_S(11,15) = -0.5*cMrSelf[8]*mnuSelf; 
 
  // ... Block from weak multiply of mOther, nuOther, M0Other and uCrossOtherZ ... // 
  data->AEM_S(8,8) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(8,9) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(8,10) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(8,11) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,8) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(9,9) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(9,10) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(9,11) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,8) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(10,9) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(10,10) += 0.5*m0rOther[0]*mnuOther; 
  data->AEM_S(10,11) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,8) += 0.5*m0rOther[3]*mnuOther; 
  data->AEM_S(11,9) += 0.5*m0rOther[2]*mnuOther; 
  data->AEM_S(11,10) += 0.5*m0rOther[1]*mnuOther; 
  data->AEM_S(11,11) += 0.5*m0rOther[0]*mnuOther; 
 
  // ... Block from correction to momentum conservation (other) ... // 
  data->AEM_S(8,12) += -(0.5*cMrOther[8]*mSelf*mnuOther)/mOther; 
  data->AEM_S(8,13) += -(0.5*cMrOther[9]*mSelf*mnuOther)/mOther; 
  data->AEM_S(8,14) += -(0.5*cMrOther[10]*mSelf*mnuOther)/mOther; 
  data->AEM_S(8,15) += -(0.5*cMrOther[11]*mSelf*mnuOther)/mOther; 
  data->AEM_S(9,12) += -(0.5*cMrOther[9]*mSelf*mnuOther)/mOther; 
  data->AEM_S(9,13) += -(0.5*cMrOther[8]*mSelf*mnuOther)/mOther; 
  data->AEM_S(9,14) += -(0.5*cMrOther[11]*mSelf*mnuOther)/mOther; 
  data->AEM_S(9,15) += -(0.5*cMrOther[10]*mSelf*mnuOther)/mOther; 
  data->AEM_S(10,12) += -(0.5*cMrOther[10]*mSelf*mnuOther)/mOther; 
  data->AEM_S(10,13) += -(0.5*cMrOther[11]*mSelf*mnuOther)/mOther; 
  data->AEM_S(10,14) += -(0.5*cMrOther[8]*mSelf*mnuOther)/mOther; 
  data->AEM_S(10,15) += -(0.5*cMrOther[9]*mSelf*mnuOther)/mOther; 
  data->AEM_S(11,12) += -(0.5*cMrOther[11]*mSelf*mnuOther)/mOther; 
  data->AEM_S(11,13) += -(0.5*cMrOther[10]*mSelf*mnuOther)/mOther; 
  data->AEM_S(11,14) += -(0.5*cMrOther[9]*mSelf*mnuOther)/mOther; 
  data->AEM_S(11,15) += -(0.5*cMrOther[8]*mSelf*mnuOther)/mOther; 
 
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
  data->AEM_S(12,8) += 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(12,9) += 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(12,10) += 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(12,11) += 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(13,8) += 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(13,9) += 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(13,10) += 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(13,11) += 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(14,8) += 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(14,9) += 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(14,10) += 0.5*m1SrOther[8]*mnuOther; 
  data->AEM_S(14,11) += 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(15,8) += 0.5*m1SrOther[11]*mnuOther; 
  data->AEM_S(15,9) += 0.5*m1SrOther[10]*mnuOther; 
  data->AEM_S(15,10) += 0.5*m1SrOther[9]*mnuOther; 
  data->AEM_S(15,11) += 0.5*m1SrOther[8]*mnuOther; 
 
  // ... Contribution to RHS vector from component 3 of mnuM1Self+mnuM1Other. 
  mnuM1sum[8] += m1rSelf[8]*mnuSelf+m1rOther[8]*mnuOther; 
  mnuM1sum[9] += m1rSelf[9]*mnuSelf+m1rOther[9]*mnuOther; 
  mnuM1sum[10] += m1rSelf[10]*mnuSelf+m1rOther[10]*mnuOther; 
  mnuM1sum[11] += m1rSelf[11]*mnuSelf+m1rOther[11]*mnuOther; 
 
  // ... Block from correction to energy conservation (self) ... // 
  data->AEM_S(12,12) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(12,13) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(12,14) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(12,15) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(13,12) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(13,13) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(13,14) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(13,15) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(14,12) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(14,13) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(14,14) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
  data->AEM_S(14,15) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(15,12) = 0.5*m0SrSelf[3]*mnuSelf-0.5*cErSelf[3]*mnuSelf; 
  data->AEM_S(15,13) = 0.5*m0SrSelf[2]*mnuSelf-0.5*cErSelf[2]*mnuSelf; 
  data->AEM_S(15,14) = 0.5*m0SrSelf[1]*mnuSelf-0.5*cErSelf[1]*mnuSelf; 
  data->AEM_S(15,15) = 0.5*m0SrSelf[0]*mnuSelf-0.5*cErSelf[0]*mnuSelf; 
 
  // ... Block from correction to energy conservation (other) ... // 
  data->AEM_S(12,12) += (0.5*m0SrOther[0]*mSelf*mnuOther)/mOther-(0.5*cErOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(12,13) += (0.5*m0SrOther[1]*mSelf*mnuOther)/mOther-(0.5*cErOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(12,14) += (0.5*m0SrOther[2]*mSelf*mnuOther)/mOther-(0.5*cErOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(12,15) += (0.5*m0SrOther[3]*mSelf*mnuOther)/mOther-(0.5*cErOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(13,12) += (0.5*m0SrOther[1]*mSelf*mnuOther)/mOther-(0.5*cErOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(13,13) += (0.5*m0SrOther[0]*mSelf*mnuOther)/mOther-(0.5*cErOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(13,14) += (0.5*m0SrOther[3]*mSelf*mnuOther)/mOther-(0.5*cErOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(13,15) += (0.5*m0SrOther[2]*mSelf*mnuOther)/mOther-(0.5*cErOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(14,12) += (0.5*m0SrOther[2]*mSelf*mnuOther)/mOther-(0.5*cErOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(14,13) += (0.5*m0SrOther[3]*mSelf*mnuOther)/mOther-(0.5*cErOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(14,14) += (0.5*m0SrOther[0]*mSelf*mnuOther)/mOther-(0.5*cErOther[0]*mSelf*mnuOther)/mOther; 
  data->AEM_S(14,15) += (0.5*m0SrOther[1]*mSelf*mnuOther)/mOther-(0.5*cErOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(15,12) += (0.5*m0SrOther[3]*mSelf*mnuOther)/mOther-(0.5*cErOther[3]*mSelf*mnuOther)/mOther; 
  data->AEM_S(15,13) += (0.5*m0SrOther[2]*mSelf*mnuOther)/mOther-(0.5*cErOther[2]*mSelf*mnuOther)/mOther; 
  data->AEM_S(15,14) += (0.5*m0SrOther[1]*mSelf*mnuOther)/mOther-(0.5*cErOther[1]*mSelf*mnuOther)/mOther; 
  data->AEM_S(15,15) += (0.5*m0SrOther[0]*mSelf*mnuOther)/mOther-(0.5*cErOther[0]*mSelf*mnuOther)/mOther; 
 
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
 
  // ....... RHS vector is composed of m1 and m2 .......... // 
  data->BEV_S << mnuM1sum[0],mnuM1sum[1],mnuM1sum[2],mnuM1sum[3],mnuM1sum[4],mnuM1sum[5],mnuM1sum[6],mnuM1sum[7],mnuM1sum[8],mnuM1sum[9],mnuM1sum[10],mnuM1sum[11],mnuM2sum[0],mnuM2sum[1],mnuM2sum[2],mnuM2sum[3]; 
 
  data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
  Eigen::Map<VectorXd>(uCrossSelf,12,1) = data->u_S.segment<12>(0); 
 
  Eigen::Map<VectorXd>(vtSqCrossSelf,4,1) = data->u_S.segment<4>(12); 
 
  uCrossOther[0] = uCrossSelf[0]; 
  uCrossOther[1] = uCrossSelf[1]; 
  uCrossOther[2] = uCrossSelf[2]; 
  uCrossOther[3] = uCrossSelf[3]; 
  uCrossOther[4] = uCrossSelf[4]; 
  uCrossOther[5] = uCrossSelf[5]; 
  uCrossOther[6] = uCrossSelf[6]; 
  uCrossOther[7] = uCrossSelf[7]; 
  uCrossOther[8] = uCrossSelf[8]; 
  uCrossOther[9] = uCrossSelf[9]; 
  uCrossOther[10] = uCrossSelf[10]; 
  uCrossOther[11] = uCrossSelf[11]; 
  vtSqCrossOther[0] = (vtSqCrossSelf[0]*mSelf)/mOther; 
  vtSqCrossOther[1] = (vtSqCrossSelf[1]*mSelf)/mOther; 
  vtSqCrossOther[2] = (vtSqCrossSelf[2]*mSelf)/mOther; 
  vtSqCrossOther[3] = (vtSqCrossSelf[3]*mSelf)/mOther; 
} 
 
