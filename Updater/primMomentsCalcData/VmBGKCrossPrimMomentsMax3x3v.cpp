#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmBGKCrossPrimMoments3x3vMax_P1(binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *m0Other, const double *uOther, const double *vtSqOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:     free parameter beta+1. This has to be >0. 
  // m, nu:            mass and collisionality. 
  // m0:               zeroth moment of the distribution function. 
  // u,vtSq:           self primitive moments: mean flow velocity and thermal speed squared. 
  // uCross,vtSqCross: cross primitive moments: mean flow velocity and thermal speed squared. 
 
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
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m0rSelf[3] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
    m0rSelf[3] = m0Self[3]; 
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
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m0rOther[3] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
    m0rOther[3] = m0Other[3]; 
  } 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
 
  double uRelDmnu[12]; 
 
  // ... Divide (uSelfX-uOtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[0] = uSelf[0]-1.0*uOther[0]; 
  uRelDmnu[1] = uSelf[1]-1.0*uOther[1]; 
  uRelDmnu[2] = uSelf[2]-1.0*uOther[2]; 
  uRelDmnu[3] = uSelf[3]-1.0*uOther[3]; 
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
  dataDiv->BEV_S << uRelDmnu[0],uRelDmnu[1],uRelDmnu[2],uRelDmnu[3]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+0,4,1) = dataDiv->u_S; 
 
  // ... Component 1 of cross-velocity of this species ... // 
  uCrossSelf[0] = ((-0.3535533905932737*m0rOther[3]*uRelDmnu[3])-0.3535533905932737*m0rOther[2]*uRelDmnu[2]-0.3535533905932737*m0rOther[1]*uRelDmnu[1]-0.3535533905932737*m0rOther[0]*uRelDmnu[0])*betaGreenep1*mnuOther+uSelf[0]; 
  uCrossSelf[1] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[1])-0.3535533905932737*uRelDmnu[0]*m0rOther[1])*betaGreenep1*mnuOther+uSelf[1]; 
  uCrossSelf[2] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[2])-0.3535533905932737*uRelDmnu[0]*m0rOther[2])*betaGreenep1*mnuOther+uSelf[2]; 
  uCrossSelf[3] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[3])-0.3535533905932737*uRelDmnu[0]*m0rOther[3])*betaGreenep1*mnuOther+uSelf[3]; 
 
  // ... Component 1 of cross-velocity of the other species ... // 
  uCrossOther[0] = (0.3535533905932737*m0rSelf[3]*uRelDmnu[3]+0.3535533905932737*m0rSelf[2]*uRelDmnu[2]+0.3535533905932737*m0rSelf[1]*uRelDmnu[1]+0.3535533905932737*m0rSelf[0]*uRelDmnu[0])*betaGreenep1*mnuSelf+uOther[0]; 
  uCrossOther[1] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[1]+0.3535533905932737*uRelDmnu[0]*m0rSelf[1])*betaGreenep1*mnuSelf+uOther[1]; 
  uCrossOther[2] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[2]+0.3535533905932737*uRelDmnu[0]*m0rSelf[2])*betaGreenep1*mnuSelf+uOther[2]; 
  uCrossOther[3] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[3]+0.3535533905932737*uRelDmnu[0]*m0rSelf[3])*betaGreenep1*mnuSelf+uOther[3]; 
 
  // ... Divide (uSelfY-uOtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[4] = uSelf[4]-1.0*uOther[4]; 
  uRelDmnu[5] = uSelf[5]-1.0*uOther[5]; 
  uRelDmnu[6] = uSelf[6]-1.0*uOther[6]; 
  uRelDmnu[7] = uSelf[7]-1.0*uOther[7]; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[4],uRelDmnu[5],uRelDmnu[6],uRelDmnu[7]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+4,4,1) = dataDiv->u_S; 
 
  // ... Component 2 of cross-velocity of this species ... // 
  uCrossSelf[4] = ((-0.3535533905932737*m0rOther[3]*uRelDmnu[7])-0.3535533905932737*m0rOther[2]*uRelDmnu[6]-0.3535533905932737*m0rOther[1]*uRelDmnu[5]-0.3535533905932737*m0rOther[0]*uRelDmnu[4])*betaGreenep1*mnuOther+uSelf[4]; 
  uCrossSelf[5] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[5])-0.3535533905932737*m0rOther[1]*uRelDmnu[4])*betaGreenep1*mnuOther+uSelf[5]; 
  uCrossSelf[6] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[6])-0.3535533905932737*m0rOther[2]*uRelDmnu[4])*betaGreenep1*mnuOther+uSelf[6]; 
  uCrossSelf[7] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[7])-0.3535533905932737*m0rOther[3]*uRelDmnu[4])*betaGreenep1*mnuOther+uSelf[7]; 
 
  // ... Component 2 of cross-velocity of the other species ... // 
  uCrossOther[4] = (0.3535533905932737*m0rSelf[3]*uRelDmnu[7]+0.3535533905932737*m0rSelf[2]*uRelDmnu[6]+0.3535533905932737*m0rSelf[1]*uRelDmnu[5]+0.3535533905932737*m0rSelf[0]*uRelDmnu[4])*betaGreenep1*mnuSelf+uOther[4]; 
  uCrossOther[5] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[5]+0.3535533905932737*m0rSelf[1]*uRelDmnu[4])*betaGreenep1*mnuSelf+uOther[5]; 
  uCrossOther[6] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[6]+0.3535533905932737*m0rSelf[2]*uRelDmnu[4])*betaGreenep1*mnuSelf+uOther[6]; 
  uCrossOther[7] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[7]+0.3535533905932737*m0rSelf[3]*uRelDmnu[4])*betaGreenep1*mnuSelf+uOther[7]; 
 
  // ... Divide (uSelfZ-uOtherZ)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[8] = uSelf[8]-1.0*uOther[8]; 
  uRelDmnu[9] = uSelf[9]-1.0*uOther[9]; 
  uRelDmnu[10] = uSelf[10]-1.0*uOther[10]; 
  uRelDmnu[11] = uSelf[11]-1.0*uOther[11]; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[8],uRelDmnu[9],uRelDmnu[10],uRelDmnu[11]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+8,4,1) = dataDiv->u_S; 
 
  // ... Component 3 of cross-velocity of this species ... // 
  uCrossSelf[8] = ((-0.3535533905932737*m0rOther[3]*uRelDmnu[11])-0.3535533905932737*m0rOther[2]*uRelDmnu[10]-0.3535533905932737*m0rOther[1]*uRelDmnu[9]-0.3535533905932737*m0rOther[0]*uRelDmnu[8])*betaGreenep1*mnuOther+uSelf[8]; 
  uCrossSelf[9] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[9])-0.3535533905932737*m0rOther[1]*uRelDmnu[8])*betaGreenep1*mnuOther+uSelf[9]; 
  uCrossSelf[10] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[10])-0.3535533905932737*m0rOther[2]*uRelDmnu[8])*betaGreenep1*mnuOther+uSelf[10]; 
  uCrossSelf[11] = ((-0.3535533905932737*m0rOther[0]*uRelDmnu[11])-0.3535533905932737*m0rOther[3]*uRelDmnu[8])*betaGreenep1*mnuOther+uSelf[11]; 
 
  // ... Component 3 of cross-velocity of the other species ... // 
  uCrossOther[8] = (0.3535533905932737*m0rSelf[3]*uRelDmnu[11]+0.3535533905932737*m0rSelf[2]*uRelDmnu[10]+0.3535533905932737*m0rSelf[1]*uRelDmnu[9]+0.3535533905932737*m0rSelf[0]*uRelDmnu[8])*betaGreenep1*mnuSelf+uOther[8]; 
  uCrossOther[9] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[9]+0.3535533905932737*m0rSelf[1]*uRelDmnu[8])*betaGreenep1*mnuSelf+uOther[9]; 
  uCrossOther[10] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[10]+0.3535533905932737*m0rSelf[2]*uRelDmnu[8])*betaGreenep1*mnuSelf+uOther[10]; 
  uCrossOther[11] = (0.3535533905932737*m0rSelf[0]*uRelDmnu[11]+0.3535533905932737*m0rSelf[3]*uRelDmnu[8])*betaGreenep1*mnuSelf+uOther[11]; 
 
  double uRelSq[4]; 
  // Zero out array with dot product of uSelf-uOther with itself. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    uRelSq[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    uRelSq[0] += 0.3535533905932737*uSelf[a0+3]*uSelf[a0+3]-0.7071067811865475*uOther[a0+3]*uSelf[a0+3]+0.3535533905932737*uOther[a0+3]*uOther[a0+3]+0.3535533905932737*uSelf[a0+2]*uSelf[a0+2]-0.7071067811865475*uOther[a0+2]*uSelf[a0+2]+0.3535533905932737*uOther[a0+2]*uOther[a0+2]+0.3535533905932737*uSelf[a0+1]*uSelf[a0+1]-0.7071067811865475*uOther[a0+1]*uSelf[a0+1]+0.3535533905932737*uOther[a0+1]*uOther[a0+1]+0.3535533905932737*uSelf[a0]*uSelf[a0]-0.7071067811865475*uOther[a0]*uSelf[a0]+0.3535533905932737*uOther[a0]*uOther[a0]; 
    uRelSq[1] += 0.7071067811865475*uSelf[a0]*uSelf[a0+1]-0.7071067811865475*uOther[a0]*uSelf[a0+1]-0.7071067811865475*uSelf[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*uOther[a0+1]; 
    uRelSq[2] += 0.7071067811865475*uSelf[a0]*uSelf[a0+2]-0.7071067811865475*uOther[a0]*uSelf[a0+2]-0.7071067811865475*uSelf[a0]*uOther[a0+2]+0.7071067811865475*uOther[a0]*uOther[a0+2]; 
    uRelSq[3] += 0.7071067811865475*uSelf[a0]*uSelf[a0+3]-0.7071067811865475*uOther[a0]*uSelf[a0+3]-0.7071067811865475*uSelf[a0]*uOther[a0+3]+0.7071067811865475*uOther[a0]*uOther[a0+3]; 
  } 
 
  double relKinE[4]; 
  // Zero out array with ((beta+1)/2)*(mSelf+mOther)*(uSelf-uOther) . uRelDmnu. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += betaGreenep1*(0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+3]*mSelf-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+3]*mSelf+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+2]*mSelf-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+2]*mSelf+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+1]*mSelf-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+1]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0]*mSelf+0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+3]*mOther-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+3]*mOther+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+2]*mOther-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+2]*mOther+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+1]*mOther-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+1]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0]*mOther); 
    relKinE[1] += betaGreenep1*(0.1767766952966368*uRelDmnu[a0]*uSelf[a0+1]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+1]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+1]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+1]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+1]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+1]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+1]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+1]*mOther); 
    relKinE[2] += betaGreenep1*(0.1767766952966368*uRelDmnu[a0]*uSelf[a0+2]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+2]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+2]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+2]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+2]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+2]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+2]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+2]*mOther); 
    relKinE[3] += betaGreenep1*(0.1767766952966368*uRelDmnu[a0]*uSelf[a0+3]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+3]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+3]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+3]*mOther); 
  } 
 
  double Tdiff[4]; 
  Tdiff[0] = vtSqSelf[0]*mSelf-1.0*vtSqOther[0]*mOther; 
  Tdiff[1] = vtSqSelf[1]*mSelf-1.0*vtSqOther[1]*mOther; 
  Tdiff[2] = vtSqSelf[2]*mSelf-1.0*vtSqOther[2]*mOther; 
  Tdiff[3] = vtSqSelf[3]*mSelf-1.0*vtSqOther[3]*mOther; 
 
  double diffSelf[4]; 
  diffSelf[0] = (0.3535533905932737*m0rOther[3]*relKinE[3]+0.3535533905932737*m0rOther[2]*relKinE[2]+0.3535533905932737*m0rOther[1]*relKinE[1]+0.3535533905932737*m0rOther[0]*relKinE[0])*mnuOther-1.0*uRelSq[0]*mOther+3.0*Tdiff[0]; 
  diffSelf[1] = (0.3535533905932737*m0rOther[0]*relKinE[1]+0.3535533905932737*relKinE[0]*m0rOther[1])*mnuOther-1.0*uRelSq[1]*mOther+3.0*Tdiff[1]; 
  diffSelf[2] = (0.3535533905932737*m0rOther[0]*relKinE[2]+0.3535533905932737*relKinE[0]*m0rOther[2])*mnuOther-1.0*uRelSq[2]*mOther+3.0*Tdiff[2]; 
  diffSelf[3] = (0.3535533905932737*m0rOther[0]*relKinE[3]+0.3535533905932737*relKinE[0]*m0rOther[3])*mnuOther-1.0*uRelSq[3]*mOther+3.0*Tdiff[3]; 
 
  // Divide diffSelf by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffSelf[0],diffSelf[1],diffSelf[2],diffSelf[3]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaSelf[4]; 
  Eigen::Map<VectorXd>(vtSqDeltaSelf,4,1) = dataDiv->u_S; 
 
  double diffOther[4]; 
  diffOther[0] = (0.3535533905932737*m0rSelf[3]*relKinE[3]+0.3535533905932737*m0rSelf[2]*relKinE[2]+0.3535533905932737*m0rSelf[1]*relKinE[1]+0.3535533905932737*m0rSelf[0]*relKinE[0])*mnuSelf-1.0*uRelSq[0]*mSelf-3.0*Tdiff[0]; 
  diffOther[1] = (0.3535533905932737*m0rSelf[0]*relKinE[1]+0.3535533905932737*relKinE[0]*m0rSelf[1])*mnuSelf-1.0*uRelSq[1]*mSelf-3.0*Tdiff[1]; 
  diffOther[2] = (0.3535533905932737*m0rSelf[0]*relKinE[2]+0.3535533905932737*relKinE[0]*m0rSelf[2])*mnuSelf-1.0*uRelSq[2]*mSelf-3.0*Tdiff[2]; 
  diffOther[3] = (0.3535533905932737*m0rSelf[0]*relKinE[3]+0.3535533905932737*relKinE[0]*m0rSelf[3])*mnuSelf-1.0*uRelSq[3]*mSelf-3.0*Tdiff[3]; 
 
  // Divide diffOther by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffOther[0],diffOther[1],diffOther[2],diffOther[3]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaOther[4]; 
  Eigen::Map<VectorXd>(vtSqDeltaOther,4,1) = dataDiv->u_S; 
 
  // ... Cross-thermal speeds (squared) ... // 
  double deltaFacOther = (0.6666666666666666*mnuOther)/(mSelf+mOther); 
  vtSqCrossSelf[0] = (-0.3535533905932737*m0rOther[3]*vtSqDeltaSelf[3]*deltaFacOther)-0.3535533905932737*m0rOther[2]*vtSqDeltaSelf[2]*deltaFacOther-0.3535533905932737*m0rOther[1]*vtSqDeltaSelf[1]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[0]*deltaFacOther+vtSqSelf[0]; 
  vtSqCrossSelf[1] = (-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[1]*deltaFacOther)-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[1]*deltaFacOther+vtSqSelf[1]; 
  vtSqCrossSelf[2] = (-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[2]*deltaFacOther)-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[2]*deltaFacOther+vtSqSelf[2]; 
  vtSqCrossSelf[3] = (-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[3]*deltaFacOther)-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[3]*deltaFacOther+vtSqSelf[3]; 
 
  double deltaFacSelf = (0.6666666666666666*mnuSelf)/(mSelf+mOther); 
  vtSqCrossOther[0] = (-0.3535533905932737*m0rSelf[3]*vtSqDeltaOther[3]*deltaFacSelf)-0.3535533905932737*m0rSelf[2]*vtSqDeltaOther[2]*deltaFacSelf-0.3535533905932737*m0rSelf[1]*vtSqDeltaOther[1]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[0]*deltaFacSelf+vtSqOther[0]; 
  vtSqCrossOther[1] = (-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[1]*deltaFacSelf)-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[1]*deltaFacSelf+vtSqOther[1]; 
  vtSqCrossOther[2] = (-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[2]*deltaFacSelf)-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[2]*deltaFacSelf+vtSqOther[2]; 
  vtSqCrossOther[3] = (-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[3]*deltaFacSelf)-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[3]*deltaFacSelf+vtSqOther[3]; 
 
} 
 
void VmBGKCrossPrimMoments3x3vMax_P2(binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *m0Other, const double *uOther, const double *vtSqOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:     free parameter beta+1. This has to be >0. 
  // m, nu:            mass and collisionality. 
  // m0:               zeroth moment of the distribution function. 
  // u,vtSq:           self primitive moments: mean flow velocity and thermal speed squared. 
  // uCross,vtSqCross: cross primitive moments: mean flow velocity and thermal speed squared. 
 
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
  } 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
 
  double uRelDmnu[30]; 
 
  // ... Divide (uSelfX-uOtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[0] = uSelf[0]-1.0*uOther[0]; 
  uRelDmnu[1] = uSelf[1]-1.0*uOther[1]; 
  uRelDmnu[2] = uSelf[2]-1.0*uOther[2]; 
  uRelDmnu[3] = uSelf[3]-1.0*uOther[3]; 
  uRelDmnu[4] = uSelf[4]-1.0*uOther[4]; 
  uRelDmnu[5] = uSelf[5]-1.0*uOther[5]; 
  uRelDmnu[6] = uSelf[6]-1.0*uOther[6]; 
  uRelDmnu[7] = uSelf[7]-1.0*uOther[7]; 
  uRelDmnu[8] = uSelf[8]-1.0*uOther[8]; 
  uRelDmnu[9] = uSelf[9]-1.0*uOther[9]; 
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
  dataDiv->BEV_S << uRelDmnu[0],uRelDmnu[1],uRelDmnu[2],uRelDmnu[3],uRelDmnu[4],uRelDmnu[5],uRelDmnu[6],uRelDmnu[7],uRelDmnu[8],uRelDmnu[9]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+0,10,1) = dataDiv->u_S; 
 
  // ... Component 1 of cross-velocity of this species ... // 
  uCrossSelf[0] = ((-0.3535533905932737*m0rOther[9]*uRelDmnu[9])-0.3535533905932737*m0rOther[8]*uRelDmnu[8]-0.3535533905932737*m0rOther[7]*uRelDmnu[7]-0.3535533905932737*m0rOther[6]*uRelDmnu[6]-0.3535533905932737*m0rOther[5]*uRelDmnu[5]-0.3535533905932737*m0rOther[4]*uRelDmnu[4]-0.3535533905932737*m0rOther[3]*uRelDmnu[3]-0.3535533905932737*m0rOther[2]*uRelDmnu[2]-0.3535533905932737*m0rOther[1]*uRelDmnu[1]-0.3535533905932737*m0rOther[0]*uRelDmnu[0])*betaGreenep1*mnuOther+uSelf[0]; 
  uCrossSelf[1] = ((-0.3162277660168379*m0rOther[1]*uRelDmnu[7])-0.3162277660168379*uRelDmnu[1]*m0rOther[7]-0.3535533905932737*m0rOther[3]*uRelDmnu[5]-0.3535533905932737*uRelDmnu[3]*m0rOther[5]-0.3535533905932737*m0rOther[2]*uRelDmnu[4]-0.3535533905932737*uRelDmnu[2]*m0rOther[4]-0.3535533905932737*m0rOther[0]*uRelDmnu[1]-0.3535533905932737*uRelDmnu[0]*m0rOther[1])*betaGreenep1*mnuOther+uSelf[1]; 
  uCrossSelf[2] = ((-0.3162277660168379*m0rOther[2]*uRelDmnu[8])-0.3162277660168379*uRelDmnu[2]*m0rOther[8]-0.3535533905932737*m0rOther[3]*uRelDmnu[6]-0.3535533905932737*uRelDmnu[3]*m0rOther[6]-0.3535533905932737*m0rOther[1]*uRelDmnu[4]-0.3535533905932737*uRelDmnu[1]*m0rOther[4]-0.3535533905932737*m0rOther[0]*uRelDmnu[2]-0.3535533905932737*uRelDmnu[0]*m0rOther[2])*betaGreenep1*mnuOther+uSelf[2]; 
  uCrossSelf[3] = ((-0.3162277660168379*m0rOther[3]*uRelDmnu[9])-0.3162277660168379*uRelDmnu[3]*m0rOther[9]-0.3535533905932737*m0rOther[2]*uRelDmnu[6]-0.3535533905932737*uRelDmnu[2]*m0rOther[6]-0.3535533905932737*m0rOther[1]*uRelDmnu[5]-0.3535533905932737*uRelDmnu[1]*m0rOther[5]-0.3535533905932737*m0rOther[0]*uRelDmnu[3]-0.3535533905932737*uRelDmnu[0]*m0rOther[3])*betaGreenep1*mnuOther+uSelf[3]; 
  uCrossSelf[4] = ((-0.3162277660168379*m0rOther[4]*uRelDmnu[8])-0.3162277660168379*uRelDmnu[4]*m0rOther[8]-0.3162277660168379*m0rOther[4]*uRelDmnu[7]-0.3162277660168379*uRelDmnu[4]*m0rOther[7]-0.3535533905932737*m0rOther[5]*uRelDmnu[6]-0.3535533905932737*uRelDmnu[5]*m0rOther[6]-0.3535533905932737*m0rOther[0]*uRelDmnu[4]-0.3535533905932737*uRelDmnu[0]*m0rOther[4]-0.3535533905932737*m0rOther[1]*uRelDmnu[2]-0.3535533905932737*uRelDmnu[1]*m0rOther[2])*betaGreenep1*mnuOther+uSelf[4]; 
  uCrossSelf[5] = ((-0.3162277660168379*m0rOther[5]*uRelDmnu[9])-0.3162277660168379*uRelDmnu[5]*m0rOther[9]-0.3162277660168379*m0rOther[5]*uRelDmnu[7]-0.3162277660168379*uRelDmnu[5]*m0rOther[7]-0.3535533905932737*m0rOther[4]*uRelDmnu[6]-0.3535533905932737*uRelDmnu[4]*m0rOther[6]-0.3535533905932737*m0rOther[0]*uRelDmnu[5]-0.3535533905932737*uRelDmnu[0]*m0rOther[5]-0.3535533905932737*m0rOther[1]*uRelDmnu[3]-0.3535533905932737*uRelDmnu[1]*m0rOther[3])*betaGreenep1*mnuOther+uSelf[5]; 
  uCrossSelf[6] = ((-0.3162277660168379*m0rOther[6]*uRelDmnu[9])-0.3162277660168379*uRelDmnu[6]*m0rOther[9]-0.3162277660168379*m0rOther[6]*uRelDmnu[8]-0.3162277660168379*uRelDmnu[6]*m0rOther[8]-0.3535533905932737*m0rOther[0]*uRelDmnu[6]-0.3535533905932737*uRelDmnu[0]*m0rOther[6]-0.3535533905932737*m0rOther[4]*uRelDmnu[5]-0.3535533905932737*uRelDmnu[4]*m0rOther[5]-0.3535533905932737*m0rOther[2]*uRelDmnu[3]-0.3535533905932737*uRelDmnu[2]*m0rOther[3])*betaGreenep1*mnuOther+uSelf[6]; 
  uCrossSelf[7] = ((-0.2258769757263128*m0rOther[7]*uRelDmnu[7])-0.3535533905932737*m0rOther[0]*uRelDmnu[7]-0.3535533905932737*uRelDmnu[0]*m0rOther[7]-0.3162277660168379*m0rOther[5]*uRelDmnu[5]-0.3162277660168379*m0rOther[4]*uRelDmnu[4]-0.3162277660168379*m0rOther[1]*uRelDmnu[1])*betaGreenep1*mnuOther+uSelf[7]; 
  uCrossSelf[8] = ((-0.2258769757263128*m0rOther[8]*uRelDmnu[8])-0.3535533905932737*m0rOther[0]*uRelDmnu[8]-0.3535533905932737*uRelDmnu[0]*m0rOther[8]-0.3162277660168379*m0rOther[6]*uRelDmnu[6]-0.3162277660168379*m0rOther[4]*uRelDmnu[4]-0.3162277660168379*m0rOther[2]*uRelDmnu[2])*betaGreenep1*mnuOther+uSelf[8]; 
  uCrossSelf[9] = ((-0.2258769757263128*m0rOther[9]*uRelDmnu[9])-0.3535533905932737*m0rOther[0]*uRelDmnu[9]-0.3535533905932737*uRelDmnu[0]*m0rOther[9]-0.3162277660168379*m0rOther[6]*uRelDmnu[6]-0.3162277660168379*m0rOther[5]*uRelDmnu[5]-0.3162277660168379*m0rOther[3]*uRelDmnu[3])*betaGreenep1*mnuOther+uSelf[9]; 
 
  // ... Component 1 of cross-velocity of the other species ... // 
  uCrossOther[0] = (0.3535533905932737*m0rSelf[9]*uRelDmnu[9]+0.3535533905932737*m0rSelf[8]*uRelDmnu[8]+0.3535533905932737*m0rSelf[7]*uRelDmnu[7]+0.3535533905932737*m0rSelf[6]*uRelDmnu[6]+0.3535533905932737*m0rSelf[5]*uRelDmnu[5]+0.3535533905932737*m0rSelf[4]*uRelDmnu[4]+0.3535533905932737*m0rSelf[3]*uRelDmnu[3]+0.3535533905932737*m0rSelf[2]*uRelDmnu[2]+0.3535533905932737*m0rSelf[1]*uRelDmnu[1]+0.3535533905932737*m0rSelf[0]*uRelDmnu[0])*betaGreenep1*mnuSelf+uOther[0]; 
  uCrossOther[1] = (0.3162277660168379*m0rSelf[1]*uRelDmnu[7]+0.3162277660168379*uRelDmnu[1]*m0rSelf[7]+0.3535533905932737*m0rSelf[3]*uRelDmnu[5]+0.3535533905932737*uRelDmnu[3]*m0rSelf[5]+0.3535533905932737*m0rSelf[2]*uRelDmnu[4]+0.3535533905932737*uRelDmnu[2]*m0rSelf[4]+0.3535533905932737*m0rSelf[0]*uRelDmnu[1]+0.3535533905932737*uRelDmnu[0]*m0rSelf[1])*betaGreenep1*mnuSelf+uOther[1]; 
  uCrossOther[2] = (0.3162277660168379*m0rSelf[2]*uRelDmnu[8]+0.3162277660168379*uRelDmnu[2]*m0rSelf[8]+0.3535533905932737*m0rSelf[3]*uRelDmnu[6]+0.3535533905932737*uRelDmnu[3]*m0rSelf[6]+0.3535533905932737*m0rSelf[1]*uRelDmnu[4]+0.3535533905932737*uRelDmnu[1]*m0rSelf[4]+0.3535533905932737*m0rSelf[0]*uRelDmnu[2]+0.3535533905932737*uRelDmnu[0]*m0rSelf[2])*betaGreenep1*mnuSelf+uOther[2]; 
  uCrossOther[3] = (0.3162277660168379*m0rSelf[3]*uRelDmnu[9]+0.3162277660168379*uRelDmnu[3]*m0rSelf[9]+0.3535533905932737*m0rSelf[2]*uRelDmnu[6]+0.3535533905932737*uRelDmnu[2]*m0rSelf[6]+0.3535533905932737*m0rSelf[1]*uRelDmnu[5]+0.3535533905932737*uRelDmnu[1]*m0rSelf[5]+0.3535533905932737*m0rSelf[0]*uRelDmnu[3]+0.3535533905932737*uRelDmnu[0]*m0rSelf[3])*betaGreenep1*mnuSelf+uOther[3]; 
  uCrossOther[4] = (0.3162277660168379*m0rSelf[4]*uRelDmnu[8]+0.3162277660168379*uRelDmnu[4]*m0rSelf[8]+0.3162277660168379*m0rSelf[4]*uRelDmnu[7]+0.3162277660168379*uRelDmnu[4]*m0rSelf[7]+0.3535533905932737*m0rSelf[5]*uRelDmnu[6]+0.3535533905932737*uRelDmnu[5]*m0rSelf[6]+0.3535533905932737*m0rSelf[0]*uRelDmnu[4]+0.3535533905932737*uRelDmnu[0]*m0rSelf[4]+0.3535533905932737*m0rSelf[1]*uRelDmnu[2]+0.3535533905932737*uRelDmnu[1]*m0rSelf[2])*betaGreenep1*mnuSelf+uOther[4]; 
  uCrossOther[5] = (0.3162277660168379*m0rSelf[5]*uRelDmnu[9]+0.3162277660168379*uRelDmnu[5]*m0rSelf[9]+0.3162277660168379*m0rSelf[5]*uRelDmnu[7]+0.3162277660168379*uRelDmnu[5]*m0rSelf[7]+0.3535533905932737*m0rSelf[4]*uRelDmnu[6]+0.3535533905932737*uRelDmnu[4]*m0rSelf[6]+0.3535533905932737*m0rSelf[0]*uRelDmnu[5]+0.3535533905932737*uRelDmnu[0]*m0rSelf[5]+0.3535533905932737*m0rSelf[1]*uRelDmnu[3]+0.3535533905932737*uRelDmnu[1]*m0rSelf[3])*betaGreenep1*mnuSelf+uOther[5]; 
  uCrossOther[6] = (0.3162277660168379*m0rSelf[6]*uRelDmnu[9]+0.3162277660168379*uRelDmnu[6]*m0rSelf[9]+0.3162277660168379*m0rSelf[6]*uRelDmnu[8]+0.3162277660168379*uRelDmnu[6]*m0rSelf[8]+0.3535533905932737*m0rSelf[0]*uRelDmnu[6]+0.3535533905932737*uRelDmnu[0]*m0rSelf[6]+0.3535533905932737*m0rSelf[4]*uRelDmnu[5]+0.3535533905932737*uRelDmnu[4]*m0rSelf[5]+0.3535533905932737*m0rSelf[2]*uRelDmnu[3]+0.3535533905932737*uRelDmnu[2]*m0rSelf[3])*betaGreenep1*mnuSelf+uOther[6]; 
  uCrossOther[7] = (0.2258769757263128*m0rSelf[7]*uRelDmnu[7]+0.3535533905932737*m0rSelf[0]*uRelDmnu[7]+0.3535533905932737*uRelDmnu[0]*m0rSelf[7]+0.3162277660168379*m0rSelf[5]*uRelDmnu[5]+0.3162277660168379*m0rSelf[4]*uRelDmnu[4]+0.3162277660168379*m0rSelf[1]*uRelDmnu[1])*betaGreenep1*mnuSelf+uOther[7]; 
  uCrossOther[8] = (0.2258769757263128*m0rSelf[8]*uRelDmnu[8]+0.3535533905932737*m0rSelf[0]*uRelDmnu[8]+0.3535533905932737*uRelDmnu[0]*m0rSelf[8]+0.3162277660168379*m0rSelf[6]*uRelDmnu[6]+0.3162277660168379*m0rSelf[4]*uRelDmnu[4]+0.3162277660168379*m0rSelf[2]*uRelDmnu[2])*betaGreenep1*mnuSelf+uOther[8]; 
  uCrossOther[9] = (0.2258769757263128*m0rSelf[9]*uRelDmnu[9]+0.3535533905932737*m0rSelf[0]*uRelDmnu[9]+0.3535533905932737*uRelDmnu[0]*m0rSelf[9]+0.3162277660168379*m0rSelf[6]*uRelDmnu[6]+0.3162277660168379*m0rSelf[5]*uRelDmnu[5]+0.3162277660168379*m0rSelf[3]*uRelDmnu[3])*betaGreenep1*mnuSelf+uOther[9]; 
 
  // ... Divide (uSelfY-uOtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[10] = uSelf[10]-1.0*uOther[10]; 
  uRelDmnu[11] = uSelf[11]-1.0*uOther[11]; 
  uRelDmnu[12] = uSelf[12]-1.0*uOther[12]; 
  uRelDmnu[13] = uSelf[13]-1.0*uOther[13]; 
  uRelDmnu[14] = uSelf[14]-1.0*uOther[14]; 
  uRelDmnu[15] = uSelf[15]-1.0*uOther[15]; 
  uRelDmnu[16] = uSelf[16]-1.0*uOther[16]; 
  uRelDmnu[17] = uSelf[17]-1.0*uOther[17]; 
  uRelDmnu[18] = uSelf[18]-1.0*uOther[18]; 
  uRelDmnu[19] = uSelf[19]-1.0*uOther[19]; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[10],uRelDmnu[11],uRelDmnu[12],uRelDmnu[13],uRelDmnu[14],uRelDmnu[15],uRelDmnu[16],uRelDmnu[17],uRelDmnu[18],uRelDmnu[19]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+10,10,1) = dataDiv->u_S; 
 
  // ... Component 2 of cross-velocity of this species ... // 
  uCrossSelf[10] = ((-0.3535533905932737*m0rOther[9]*uRelDmnu[19])-0.3535533905932737*m0rOther[8]*uRelDmnu[18]-0.3535533905932737*m0rOther[7]*uRelDmnu[17]-0.3535533905932737*m0rOther[6]*uRelDmnu[16]-0.3535533905932737*m0rOther[5]*uRelDmnu[15]-0.3535533905932737*m0rOther[4]*uRelDmnu[14]-0.3535533905932737*m0rOther[3]*uRelDmnu[13]-0.3535533905932737*m0rOther[2]*uRelDmnu[12]-0.3535533905932737*m0rOther[1]*uRelDmnu[11]-0.3535533905932737*m0rOther[0]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[10]; 
  uCrossSelf[11] = ((-0.3162277660168379*m0rOther[1]*uRelDmnu[17])-0.3535533905932737*m0rOther[3]*uRelDmnu[15]-0.3535533905932737*m0rOther[2]*uRelDmnu[14]-0.3535533905932737*m0rOther[5]*uRelDmnu[13]-0.3535533905932737*m0rOther[4]*uRelDmnu[12]-0.3162277660168379*m0rOther[7]*uRelDmnu[11]-0.3535533905932737*m0rOther[0]*uRelDmnu[11]-0.3535533905932737*m0rOther[1]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[11]; 
  uCrossSelf[12] = ((-0.3162277660168379*m0rOther[2]*uRelDmnu[18])-0.3535533905932737*m0rOther[3]*uRelDmnu[16]-0.3535533905932737*m0rOther[1]*uRelDmnu[14]-0.3535533905932737*m0rOther[6]*uRelDmnu[13]-0.3162277660168379*m0rOther[8]*uRelDmnu[12]-0.3535533905932737*m0rOther[0]*uRelDmnu[12]-0.3535533905932737*m0rOther[4]*uRelDmnu[11]-0.3535533905932737*m0rOther[2]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[12]; 
  uCrossSelf[13] = ((-0.3162277660168379*m0rOther[3]*uRelDmnu[19])-0.3535533905932737*m0rOther[2]*uRelDmnu[16]-0.3535533905932737*m0rOther[1]*uRelDmnu[15]-0.3162277660168379*m0rOther[9]*uRelDmnu[13]-0.3535533905932737*m0rOther[0]*uRelDmnu[13]-0.3535533905932737*m0rOther[6]*uRelDmnu[12]-0.3535533905932737*m0rOther[5]*uRelDmnu[11]-0.3535533905932737*m0rOther[3]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[13]; 
  uCrossSelf[14] = ((-0.3162277660168379*m0rOther[4]*uRelDmnu[18])-0.3162277660168379*m0rOther[4]*uRelDmnu[17]-0.3535533905932737*m0rOther[5]*uRelDmnu[16]-0.3535533905932737*m0rOther[6]*uRelDmnu[15]-0.3162277660168379*m0rOther[8]*uRelDmnu[14]-0.3162277660168379*m0rOther[7]*uRelDmnu[14]-0.3535533905932737*m0rOther[0]*uRelDmnu[14]-0.3535533905932737*m0rOther[1]*uRelDmnu[12]-0.3535533905932737*m0rOther[2]*uRelDmnu[11]-0.3535533905932737*m0rOther[4]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[14]; 
  uCrossSelf[15] = ((-0.3162277660168379*m0rOther[5]*uRelDmnu[19])-0.3162277660168379*m0rOther[5]*uRelDmnu[17]-0.3535533905932737*m0rOther[4]*uRelDmnu[16]-0.3162277660168379*m0rOther[9]*uRelDmnu[15]-0.3162277660168379*m0rOther[7]*uRelDmnu[15]-0.3535533905932737*m0rOther[0]*uRelDmnu[15]-0.3535533905932737*m0rOther[6]*uRelDmnu[14]-0.3535533905932737*m0rOther[1]*uRelDmnu[13]-0.3535533905932737*m0rOther[3]*uRelDmnu[11]-0.3535533905932737*m0rOther[5]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[15]; 
  uCrossSelf[16] = ((-0.3162277660168379*m0rOther[6]*uRelDmnu[19])-0.3162277660168379*m0rOther[6]*uRelDmnu[18]-0.3162277660168379*m0rOther[9]*uRelDmnu[16]-0.3162277660168379*m0rOther[8]*uRelDmnu[16]-0.3535533905932737*m0rOther[0]*uRelDmnu[16]-0.3535533905932737*m0rOther[4]*uRelDmnu[15]-0.3535533905932737*m0rOther[5]*uRelDmnu[14]-0.3535533905932737*m0rOther[2]*uRelDmnu[13]-0.3535533905932737*m0rOther[3]*uRelDmnu[12]-0.3535533905932737*m0rOther[6]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[16]; 
  uCrossSelf[17] = ((-0.2258769757263128*m0rOther[7]*uRelDmnu[17])-0.3535533905932737*m0rOther[0]*uRelDmnu[17]-0.3162277660168379*m0rOther[5]*uRelDmnu[15]-0.3162277660168379*m0rOther[4]*uRelDmnu[14]-0.3162277660168379*m0rOther[1]*uRelDmnu[11]-0.3535533905932737*m0rOther[7]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[17]; 
  uCrossSelf[18] = ((-0.2258769757263128*m0rOther[8]*uRelDmnu[18])-0.3535533905932737*m0rOther[0]*uRelDmnu[18]-0.3162277660168379*m0rOther[6]*uRelDmnu[16]-0.3162277660168379*m0rOther[4]*uRelDmnu[14]-0.3162277660168379*m0rOther[2]*uRelDmnu[12]-0.3535533905932737*m0rOther[8]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[18]; 
  uCrossSelf[19] = ((-0.2258769757263128*m0rOther[9]*uRelDmnu[19])-0.3535533905932737*m0rOther[0]*uRelDmnu[19]-0.3162277660168379*m0rOther[6]*uRelDmnu[16]-0.3162277660168379*m0rOther[5]*uRelDmnu[15]-0.3162277660168379*m0rOther[3]*uRelDmnu[13]-0.3535533905932737*m0rOther[9]*uRelDmnu[10])*betaGreenep1*mnuOther+uSelf[19]; 
 
  // ... Component 2 of cross-velocity of the other species ... // 
  uCrossOther[10] = (0.3535533905932737*m0rSelf[9]*uRelDmnu[19]+0.3535533905932737*m0rSelf[8]*uRelDmnu[18]+0.3535533905932737*m0rSelf[7]*uRelDmnu[17]+0.3535533905932737*m0rSelf[6]*uRelDmnu[16]+0.3535533905932737*m0rSelf[5]*uRelDmnu[15]+0.3535533905932737*m0rSelf[4]*uRelDmnu[14]+0.3535533905932737*m0rSelf[3]*uRelDmnu[13]+0.3535533905932737*m0rSelf[2]*uRelDmnu[12]+0.3535533905932737*m0rSelf[1]*uRelDmnu[11]+0.3535533905932737*m0rSelf[0]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[10]; 
  uCrossOther[11] = (0.3162277660168379*m0rSelf[1]*uRelDmnu[17]+0.3535533905932737*m0rSelf[3]*uRelDmnu[15]+0.3535533905932737*m0rSelf[2]*uRelDmnu[14]+0.3535533905932737*m0rSelf[5]*uRelDmnu[13]+0.3535533905932737*m0rSelf[4]*uRelDmnu[12]+0.3162277660168379*m0rSelf[7]*uRelDmnu[11]+0.3535533905932737*m0rSelf[0]*uRelDmnu[11]+0.3535533905932737*m0rSelf[1]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[11]; 
  uCrossOther[12] = (0.3162277660168379*m0rSelf[2]*uRelDmnu[18]+0.3535533905932737*m0rSelf[3]*uRelDmnu[16]+0.3535533905932737*m0rSelf[1]*uRelDmnu[14]+0.3535533905932737*m0rSelf[6]*uRelDmnu[13]+0.3162277660168379*m0rSelf[8]*uRelDmnu[12]+0.3535533905932737*m0rSelf[0]*uRelDmnu[12]+0.3535533905932737*m0rSelf[4]*uRelDmnu[11]+0.3535533905932737*m0rSelf[2]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[12]; 
  uCrossOther[13] = (0.3162277660168379*m0rSelf[3]*uRelDmnu[19]+0.3535533905932737*m0rSelf[2]*uRelDmnu[16]+0.3535533905932737*m0rSelf[1]*uRelDmnu[15]+0.3162277660168379*m0rSelf[9]*uRelDmnu[13]+0.3535533905932737*m0rSelf[0]*uRelDmnu[13]+0.3535533905932737*m0rSelf[6]*uRelDmnu[12]+0.3535533905932737*m0rSelf[5]*uRelDmnu[11]+0.3535533905932737*m0rSelf[3]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[13]; 
  uCrossOther[14] = (0.3162277660168379*m0rSelf[4]*uRelDmnu[18]+0.3162277660168379*m0rSelf[4]*uRelDmnu[17]+0.3535533905932737*m0rSelf[5]*uRelDmnu[16]+0.3535533905932737*m0rSelf[6]*uRelDmnu[15]+0.3162277660168379*m0rSelf[8]*uRelDmnu[14]+0.3162277660168379*m0rSelf[7]*uRelDmnu[14]+0.3535533905932737*m0rSelf[0]*uRelDmnu[14]+0.3535533905932737*m0rSelf[1]*uRelDmnu[12]+0.3535533905932737*m0rSelf[2]*uRelDmnu[11]+0.3535533905932737*m0rSelf[4]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[14]; 
  uCrossOther[15] = (0.3162277660168379*m0rSelf[5]*uRelDmnu[19]+0.3162277660168379*m0rSelf[5]*uRelDmnu[17]+0.3535533905932737*m0rSelf[4]*uRelDmnu[16]+0.3162277660168379*m0rSelf[9]*uRelDmnu[15]+0.3162277660168379*m0rSelf[7]*uRelDmnu[15]+0.3535533905932737*m0rSelf[0]*uRelDmnu[15]+0.3535533905932737*m0rSelf[6]*uRelDmnu[14]+0.3535533905932737*m0rSelf[1]*uRelDmnu[13]+0.3535533905932737*m0rSelf[3]*uRelDmnu[11]+0.3535533905932737*m0rSelf[5]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[15]; 
  uCrossOther[16] = (0.3162277660168379*m0rSelf[6]*uRelDmnu[19]+0.3162277660168379*m0rSelf[6]*uRelDmnu[18]+0.3162277660168379*m0rSelf[9]*uRelDmnu[16]+0.3162277660168379*m0rSelf[8]*uRelDmnu[16]+0.3535533905932737*m0rSelf[0]*uRelDmnu[16]+0.3535533905932737*m0rSelf[4]*uRelDmnu[15]+0.3535533905932737*m0rSelf[5]*uRelDmnu[14]+0.3535533905932737*m0rSelf[2]*uRelDmnu[13]+0.3535533905932737*m0rSelf[3]*uRelDmnu[12]+0.3535533905932737*m0rSelf[6]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[16]; 
  uCrossOther[17] = (0.2258769757263128*m0rSelf[7]*uRelDmnu[17]+0.3535533905932737*m0rSelf[0]*uRelDmnu[17]+0.3162277660168379*m0rSelf[5]*uRelDmnu[15]+0.3162277660168379*m0rSelf[4]*uRelDmnu[14]+0.3162277660168379*m0rSelf[1]*uRelDmnu[11]+0.3535533905932737*m0rSelf[7]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[17]; 
  uCrossOther[18] = (0.2258769757263128*m0rSelf[8]*uRelDmnu[18]+0.3535533905932737*m0rSelf[0]*uRelDmnu[18]+0.3162277660168379*m0rSelf[6]*uRelDmnu[16]+0.3162277660168379*m0rSelf[4]*uRelDmnu[14]+0.3162277660168379*m0rSelf[2]*uRelDmnu[12]+0.3535533905932737*m0rSelf[8]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[18]; 
  uCrossOther[19] = (0.2258769757263128*m0rSelf[9]*uRelDmnu[19]+0.3535533905932737*m0rSelf[0]*uRelDmnu[19]+0.3162277660168379*m0rSelf[6]*uRelDmnu[16]+0.3162277660168379*m0rSelf[5]*uRelDmnu[15]+0.3162277660168379*m0rSelf[3]*uRelDmnu[13]+0.3535533905932737*m0rSelf[9]*uRelDmnu[10])*betaGreenep1*mnuSelf+uOther[19]; 
 
  // ... Divide (uSelfZ-uOtherZ)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[20] = uSelf[20]-1.0*uOther[20]; 
  uRelDmnu[21] = uSelf[21]-1.0*uOther[21]; 
  uRelDmnu[22] = uSelf[22]-1.0*uOther[22]; 
  uRelDmnu[23] = uSelf[23]-1.0*uOther[23]; 
  uRelDmnu[24] = uSelf[24]-1.0*uOther[24]; 
  uRelDmnu[25] = uSelf[25]-1.0*uOther[25]; 
  uRelDmnu[26] = uSelf[26]-1.0*uOther[26]; 
  uRelDmnu[27] = uSelf[27]-1.0*uOther[27]; 
  uRelDmnu[28] = uSelf[28]-1.0*uOther[28]; 
  uRelDmnu[29] = uSelf[29]-1.0*uOther[29]; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[20],uRelDmnu[21],uRelDmnu[22],uRelDmnu[23],uRelDmnu[24],uRelDmnu[25],uRelDmnu[26],uRelDmnu[27],uRelDmnu[28],uRelDmnu[29]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+20,10,1) = dataDiv->u_S; 
 
  // ... Component 3 of cross-velocity of this species ... // 
  uCrossSelf[20] = ((-0.3535533905932737*m0rOther[9]*uRelDmnu[29])-0.3535533905932737*m0rOther[8]*uRelDmnu[28]-0.3535533905932737*m0rOther[7]*uRelDmnu[27]-0.3535533905932737*m0rOther[6]*uRelDmnu[26]-0.3535533905932737*m0rOther[5]*uRelDmnu[25]-0.3535533905932737*m0rOther[4]*uRelDmnu[24]-0.3535533905932737*m0rOther[3]*uRelDmnu[23]-0.3535533905932737*m0rOther[2]*uRelDmnu[22]-0.3535533905932737*m0rOther[1]*uRelDmnu[21]-0.3535533905932737*m0rOther[0]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[20]; 
  uCrossSelf[21] = ((-0.3162277660168379*m0rOther[1]*uRelDmnu[27])-0.3535533905932737*m0rOther[3]*uRelDmnu[25]-0.3535533905932737*m0rOther[2]*uRelDmnu[24]-0.3535533905932737*m0rOther[5]*uRelDmnu[23]-0.3535533905932737*m0rOther[4]*uRelDmnu[22]-0.3162277660168379*m0rOther[7]*uRelDmnu[21]-0.3535533905932737*m0rOther[0]*uRelDmnu[21]-0.3535533905932737*m0rOther[1]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[21]; 
  uCrossSelf[22] = ((-0.3162277660168379*m0rOther[2]*uRelDmnu[28])-0.3535533905932737*m0rOther[3]*uRelDmnu[26]-0.3535533905932737*m0rOther[1]*uRelDmnu[24]-0.3535533905932737*m0rOther[6]*uRelDmnu[23]-0.3162277660168379*m0rOther[8]*uRelDmnu[22]-0.3535533905932737*m0rOther[0]*uRelDmnu[22]-0.3535533905932737*m0rOther[4]*uRelDmnu[21]-0.3535533905932737*m0rOther[2]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[22]; 
  uCrossSelf[23] = ((-0.3162277660168379*m0rOther[3]*uRelDmnu[29])-0.3535533905932737*m0rOther[2]*uRelDmnu[26]-0.3535533905932737*m0rOther[1]*uRelDmnu[25]-0.3162277660168379*m0rOther[9]*uRelDmnu[23]-0.3535533905932737*m0rOther[0]*uRelDmnu[23]-0.3535533905932737*m0rOther[6]*uRelDmnu[22]-0.3535533905932737*m0rOther[5]*uRelDmnu[21]-0.3535533905932737*m0rOther[3]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[23]; 
  uCrossSelf[24] = ((-0.3162277660168379*m0rOther[4]*uRelDmnu[28])-0.3162277660168379*m0rOther[4]*uRelDmnu[27]-0.3535533905932737*m0rOther[5]*uRelDmnu[26]-0.3535533905932737*m0rOther[6]*uRelDmnu[25]-0.3162277660168379*m0rOther[8]*uRelDmnu[24]-0.3162277660168379*m0rOther[7]*uRelDmnu[24]-0.3535533905932737*m0rOther[0]*uRelDmnu[24]-0.3535533905932737*m0rOther[1]*uRelDmnu[22]-0.3535533905932737*m0rOther[2]*uRelDmnu[21]-0.3535533905932737*m0rOther[4]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[24]; 
  uCrossSelf[25] = ((-0.3162277660168379*m0rOther[5]*uRelDmnu[29])-0.3162277660168379*m0rOther[5]*uRelDmnu[27]-0.3535533905932737*m0rOther[4]*uRelDmnu[26]-0.3162277660168379*m0rOther[9]*uRelDmnu[25]-0.3162277660168379*m0rOther[7]*uRelDmnu[25]-0.3535533905932737*m0rOther[0]*uRelDmnu[25]-0.3535533905932737*m0rOther[6]*uRelDmnu[24]-0.3535533905932737*m0rOther[1]*uRelDmnu[23]-0.3535533905932737*m0rOther[3]*uRelDmnu[21]-0.3535533905932737*m0rOther[5]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[25]; 
  uCrossSelf[26] = ((-0.3162277660168379*m0rOther[6]*uRelDmnu[29])-0.3162277660168379*m0rOther[6]*uRelDmnu[28]-0.3162277660168379*m0rOther[9]*uRelDmnu[26]-0.3162277660168379*m0rOther[8]*uRelDmnu[26]-0.3535533905932737*m0rOther[0]*uRelDmnu[26]-0.3535533905932737*m0rOther[4]*uRelDmnu[25]-0.3535533905932737*m0rOther[5]*uRelDmnu[24]-0.3535533905932737*m0rOther[2]*uRelDmnu[23]-0.3535533905932737*m0rOther[3]*uRelDmnu[22]-0.3535533905932737*m0rOther[6]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[26]; 
  uCrossSelf[27] = ((-0.2258769757263128*m0rOther[7]*uRelDmnu[27])-0.3535533905932737*m0rOther[0]*uRelDmnu[27]-0.3162277660168379*m0rOther[5]*uRelDmnu[25]-0.3162277660168379*m0rOther[4]*uRelDmnu[24]-0.3162277660168379*m0rOther[1]*uRelDmnu[21]-0.3535533905932737*m0rOther[7]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[27]; 
  uCrossSelf[28] = ((-0.2258769757263128*m0rOther[8]*uRelDmnu[28])-0.3535533905932737*m0rOther[0]*uRelDmnu[28]-0.3162277660168379*m0rOther[6]*uRelDmnu[26]-0.3162277660168379*m0rOther[4]*uRelDmnu[24]-0.3162277660168379*m0rOther[2]*uRelDmnu[22]-0.3535533905932737*m0rOther[8]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[28]; 
  uCrossSelf[29] = ((-0.2258769757263128*m0rOther[9]*uRelDmnu[29])-0.3535533905932737*m0rOther[0]*uRelDmnu[29]-0.3162277660168379*m0rOther[6]*uRelDmnu[26]-0.3162277660168379*m0rOther[5]*uRelDmnu[25]-0.3162277660168379*m0rOther[3]*uRelDmnu[23]-0.3535533905932737*m0rOther[9]*uRelDmnu[20])*betaGreenep1*mnuOther+uSelf[29]; 
 
  // ... Component 3 of cross-velocity of the other species ... // 
  uCrossOther[20] = (0.3535533905932737*m0rSelf[9]*uRelDmnu[29]+0.3535533905932737*m0rSelf[8]*uRelDmnu[28]+0.3535533905932737*m0rSelf[7]*uRelDmnu[27]+0.3535533905932737*m0rSelf[6]*uRelDmnu[26]+0.3535533905932737*m0rSelf[5]*uRelDmnu[25]+0.3535533905932737*m0rSelf[4]*uRelDmnu[24]+0.3535533905932737*m0rSelf[3]*uRelDmnu[23]+0.3535533905932737*m0rSelf[2]*uRelDmnu[22]+0.3535533905932737*m0rSelf[1]*uRelDmnu[21]+0.3535533905932737*m0rSelf[0]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[20]; 
  uCrossOther[21] = (0.3162277660168379*m0rSelf[1]*uRelDmnu[27]+0.3535533905932737*m0rSelf[3]*uRelDmnu[25]+0.3535533905932737*m0rSelf[2]*uRelDmnu[24]+0.3535533905932737*m0rSelf[5]*uRelDmnu[23]+0.3535533905932737*m0rSelf[4]*uRelDmnu[22]+0.3162277660168379*m0rSelf[7]*uRelDmnu[21]+0.3535533905932737*m0rSelf[0]*uRelDmnu[21]+0.3535533905932737*m0rSelf[1]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[21]; 
  uCrossOther[22] = (0.3162277660168379*m0rSelf[2]*uRelDmnu[28]+0.3535533905932737*m0rSelf[3]*uRelDmnu[26]+0.3535533905932737*m0rSelf[1]*uRelDmnu[24]+0.3535533905932737*m0rSelf[6]*uRelDmnu[23]+0.3162277660168379*m0rSelf[8]*uRelDmnu[22]+0.3535533905932737*m0rSelf[0]*uRelDmnu[22]+0.3535533905932737*m0rSelf[4]*uRelDmnu[21]+0.3535533905932737*m0rSelf[2]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[22]; 
  uCrossOther[23] = (0.3162277660168379*m0rSelf[3]*uRelDmnu[29]+0.3535533905932737*m0rSelf[2]*uRelDmnu[26]+0.3535533905932737*m0rSelf[1]*uRelDmnu[25]+0.3162277660168379*m0rSelf[9]*uRelDmnu[23]+0.3535533905932737*m0rSelf[0]*uRelDmnu[23]+0.3535533905932737*m0rSelf[6]*uRelDmnu[22]+0.3535533905932737*m0rSelf[5]*uRelDmnu[21]+0.3535533905932737*m0rSelf[3]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[23]; 
  uCrossOther[24] = (0.3162277660168379*m0rSelf[4]*uRelDmnu[28]+0.3162277660168379*m0rSelf[4]*uRelDmnu[27]+0.3535533905932737*m0rSelf[5]*uRelDmnu[26]+0.3535533905932737*m0rSelf[6]*uRelDmnu[25]+0.3162277660168379*m0rSelf[8]*uRelDmnu[24]+0.3162277660168379*m0rSelf[7]*uRelDmnu[24]+0.3535533905932737*m0rSelf[0]*uRelDmnu[24]+0.3535533905932737*m0rSelf[1]*uRelDmnu[22]+0.3535533905932737*m0rSelf[2]*uRelDmnu[21]+0.3535533905932737*m0rSelf[4]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[24]; 
  uCrossOther[25] = (0.3162277660168379*m0rSelf[5]*uRelDmnu[29]+0.3162277660168379*m0rSelf[5]*uRelDmnu[27]+0.3535533905932737*m0rSelf[4]*uRelDmnu[26]+0.3162277660168379*m0rSelf[9]*uRelDmnu[25]+0.3162277660168379*m0rSelf[7]*uRelDmnu[25]+0.3535533905932737*m0rSelf[0]*uRelDmnu[25]+0.3535533905932737*m0rSelf[6]*uRelDmnu[24]+0.3535533905932737*m0rSelf[1]*uRelDmnu[23]+0.3535533905932737*m0rSelf[3]*uRelDmnu[21]+0.3535533905932737*m0rSelf[5]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[25]; 
  uCrossOther[26] = (0.3162277660168379*m0rSelf[6]*uRelDmnu[29]+0.3162277660168379*m0rSelf[6]*uRelDmnu[28]+0.3162277660168379*m0rSelf[9]*uRelDmnu[26]+0.3162277660168379*m0rSelf[8]*uRelDmnu[26]+0.3535533905932737*m0rSelf[0]*uRelDmnu[26]+0.3535533905932737*m0rSelf[4]*uRelDmnu[25]+0.3535533905932737*m0rSelf[5]*uRelDmnu[24]+0.3535533905932737*m0rSelf[2]*uRelDmnu[23]+0.3535533905932737*m0rSelf[3]*uRelDmnu[22]+0.3535533905932737*m0rSelf[6]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[26]; 
  uCrossOther[27] = (0.2258769757263128*m0rSelf[7]*uRelDmnu[27]+0.3535533905932737*m0rSelf[0]*uRelDmnu[27]+0.3162277660168379*m0rSelf[5]*uRelDmnu[25]+0.3162277660168379*m0rSelf[4]*uRelDmnu[24]+0.3162277660168379*m0rSelf[1]*uRelDmnu[21]+0.3535533905932737*m0rSelf[7]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[27]; 
  uCrossOther[28] = (0.2258769757263128*m0rSelf[8]*uRelDmnu[28]+0.3535533905932737*m0rSelf[0]*uRelDmnu[28]+0.3162277660168379*m0rSelf[6]*uRelDmnu[26]+0.3162277660168379*m0rSelf[4]*uRelDmnu[24]+0.3162277660168379*m0rSelf[2]*uRelDmnu[22]+0.3535533905932737*m0rSelf[8]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[28]; 
  uCrossOther[29] = (0.2258769757263128*m0rSelf[9]*uRelDmnu[29]+0.3535533905932737*m0rSelf[0]*uRelDmnu[29]+0.3162277660168379*m0rSelf[6]*uRelDmnu[26]+0.3162277660168379*m0rSelf[5]*uRelDmnu[25]+0.3162277660168379*m0rSelf[3]*uRelDmnu[23]+0.3535533905932737*m0rSelf[9]*uRelDmnu[20])*betaGreenep1*mnuSelf+uOther[29]; 
 
  double uRelSq[10]; 
  // Zero out array with dot product of uSelf-uOther with itself. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    uRelSq[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    uRelSq[0] += 0.3535533905932737*uSelf[a0+9]*uSelf[a0+9]-0.7071067811865475*uOther[a0+9]*uSelf[a0+9]+0.3535533905932737*uOther[a0+9]*uOther[a0+9]+0.3535533905932737*uSelf[a0+8]*uSelf[a0+8]-0.7071067811865475*uOther[a0+8]*uSelf[a0+8]+0.3535533905932737*uOther[a0+8]*uOther[a0+8]+0.3535533905932737*uSelf[a0+7]*uSelf[a0+7]-0.7071067811865475*uOther[a0+7]*uSelf[a0+7]+0.3535533905932737*uOther[a0+7]*uOther[a0+7]+0.3535533905932737*uSelf[a0+6]*uSelf[a0+6]-0.7071067811865475*uOther[a0+6]*uSelf[a0+6]+0.3535533905932737*uOther[a0+6]*uOther[a0+6]+0.3535533905932737*uSelf[a0+5]*uSelf[a0+5]-0.7071067811865475*uOther[a0+5]*uSelf[a0+5]+0.3535533905932737*uOther[a0+5]*uOther[a0+5]+0.3535533905932737*uSelf[a0+4]*uSelf[a0+4]-0.7071067811865475*uOther[a0+4]*uSelf[a0+4]+0.3535533905932737*uOther[a0+4]*uOther[a0+4]+0.3535533905932737*uSelf[a0+3]*uSelf[a0+3]-0.7071067811865475*uOther[a0+3]*uSelf[a0+3]+0.3535533905932737*uOther[a0+3]*uOther[a0+3]+0.3535533905932737*uSelf[a0+2]*uSelf[a0+2]-0.7071067811865475*uOther[a0+2]*uSelf[a0+2]+0.3535533905932737*uOther[a0+2]*uOther[a0+2]+0.3535533905932737*uSelf[a0+1]*uSelf[a0+1]-0.7071067811865475*uOther[a0+1]*uSelf[a0+1]+0.3535533905932737*uOther[a0+1]*uOther[a0+1]+0.3535533905932737*uSelf[a0]*uSelf[a0]-0.7071067811865475*uOther[a0]*uSelf[a0]+0.3535533905932737*uOther[a0]*uOther[a0]; 
    uRelSq[1] += 0.6324555320336759*uSelf[a0+1]*uSelf[a0+7]-0.6324555320336759*uOther[a0+1]*uSelf[a0+7]-0.6324555320336759*uSelf[a0+1]*uOther[a0+7]+0.6324555320336759*uOther[a0+1]*uOther[a0+7]+0.7071067811865475*uSelf[a0+3]*uSelf[a0+5]-0.7071067811865475*uOther[a0+3]*uSelf[a0+5]-0.7071067811865475*uSelf[a0+3]*uOther[a0+5]+0.7071067811865475*uOther[a0+3]*uOther[a0+5]+0.7071067811865475*uSelf[a0+2]*uSelf[a0+4]-0.7071067811865475*uOther[a0+2]*uSelf[a0+4]-0.7071067811865475*uSelf[a0+2]*uOther[a0+4]+0.7071067811865475*uOther[a0+2]*uOther[a0+4]+0.7071067811865475*uSelf[a0]*uSelf[a0+1]-0.7071067811865475*uOther[a0]*uSelf[a0+1]-0.7071067811865475*uSelf[a0]*uOther[a0+1]+0.7071067811865475*uOther[a0]*uOther[a0+1]; 
    uRelSq[2] += 0.6324555320336759*uSelf[a0+2]*uSelf[a0+8]-0.6324555320336759*uOther[a0+2]*uSelf[a0+8]-0.6324555320336759*uSelf[a0+2]*uOther[a0+8]+0.6324555320336759*uOther[a0+2]*uOther[a0+8]+0.7071067811865475*uSelf[a0+3]*uSelf[a0+6]-0.7071067811865475*uOther[a0+3]*uSelf[a0+6]-0.7071067811865475*uSelf[a0+3]*uOther[a0+6]+0.7071067811865475*uOther[a0+3]*uOther[a0+6]+0.7071067811865475*uSelf[a0+1]*uSelf[a0+4]-0.7071067811865475*uOther[a0+1]*uSelf[a0+4]-0.7071067811865475*uSelf[a0+1]*uOther[a0+4]+0.7071067811865475*uOther[a0+1]*uOther[a0+4]+0.7071067811865475*uSelf[a0]*uSelf[a0+2]-0.7071067811865475*uOther[a0]*uSelf[a0+2]-0.7071067811865475*uSelf[a0]*uOther[a0+2]+0.7071067811865475*uOther[a0]*uOther[a0+2]; 
    uRelSq[3] += 0.6324555320336759*uSelf[a0+3]*uSelf[a0+9]-0.6324555320336759*uOther[a0+3]*uSelf[a0+9]-0.6324555320336759*uSelf[a0+3]*uOther[a0+9]+0.6324555320336759*uOther[a0+3]*uOther[a0+9]+0.7071067811865475*uSelf[a0+2]*uSelf[a0+6]-0.7071067811865475*uOther[a0+2]*uSelf[a0+6]-0.7071067811865475*uSelf[a0+2]*uOther[a0+6]+0.7071067811865475*uOther[a0+2]*uOther[a0+6]+0.7071067811865475*uSelf[a0+1]*uSelf[a0+5]-0.7071067811865475*uOther[a0+1]*uSelf[a0+5]-0.7071067811865475*uSelf[a0+1]*uOther[a0+5]+0.7071067811865475*uOther[a0+1]*uOther[a0+5]+0.7071067811865475*uSelf[a0]*uSelf[a0+3]-0.7071067811865475*uOther[a0]*uSelf[a0+3]-0.7071067811865475*uSelf[a0]*uOther[a0+3]+0.7071067811865475*uOther[a0]*uOther[a0+3]; 
    uRelSq[4] += 0.6324555320336759*uSelf[a0+4]*uSelf[a0+8]-0.6324555320336759*uOther[a0+4]*uSelf[a0+8]-0.6324555320336759*uSelf[a0+4]*uOther[a0+8]+0.6324555320336759*uOther[a0+4]*uOther[a0+8]+0.6324555320336759*uSelf[a0+4]*uSelf[a0+7]-0.6324555320336759*uOther[a0+4]*uSelf[a0+7]-0.6324555320336759*uSelf[a0+4]*uOther[a0+7]+0.6324555320336759*uOther[a0+4]*uOther[a0+7]+0.7071067811865475*uSelf[a0+5]*uSelf[a0+6]-0.7071067811865475*uOther[a0+5]*uSelf[a0+6]-0.7071067811865475*uSelf[a0+5]*uOther[a0+6]+0.7071067811865475*uOther[a0+5]*uOther[a0+6]+0.7071067811865475*uSelf[a0]*uSelf[a0+4]-0.7071067811865475*uOther[a0]*uSelf[a0+4]-0.7071067811865475*uSelf[a0]*uOther[a0+4]+0.7071067811865475*uOther[a0]*uOther[a0+4]+0.7071067811865475*uSelf[a0+1]*uSelf[a0+2]-0.7071067811865475*uOther[a0+1]*uSelf[a0+2]-0.7071067811865475*uSelf[a0+1]*uOther[a0+2]+0.7071067811865475*uOther[a0+1]*uOther[a0+2]; 
    uRelSq[5] += 0.6324555320336759*uSelf[a0+5]*uSelf[a0+9]-0.6324555320336759*uOther[a0+5]*uSelf[a0+9]-0.6324555320336759*uSelf[a0+5]*uOther[a0+9]+0.6324555320336759*uOther[a0+5]*uOther[a0+9]+0.6324555320336759*uSelf[a0+5]*uSelf[a0+7]-0.6324555320336759*uOther[a0+5]*uSelf[a0+7]-0.6324555320336759*uSelf[a0+5]*uOther[a0+7]+0.6324555320336759*uOther[a0+5]*uOther[a0+7]+0.7071067811865475*uSelf[a0+4]*uSelf[a0+6]-0.7071067811865475*uOther[a0+4]*uSelf[a0+6]-0.7071067811865475*uSelf[a0+4]*uOther[a0+6]+0.7071067811865475*uOther[a0+4]*uOther[a0+6]+0.7071067811865475*uSelf[a0]*uSelf[a0+5]-0.7071067811865475*uOther[a0]*uSelf[a0+5]-0.7071067811865475*uSelf[a0]*uOther[a0+5]+0.7071067811865475*uOther[a0]*uOther[a0+5]+0.7071067811865475*uSelf[a0+1]*uSelf[a0+3]-0.7071067811865475*uOther[a0+1]*uSelf[a0+3]-0.7071067811865475*uSelf[a0+1]*uOther[a0+3]+0.7071067811865475*uOther[a0+1]*uOther[a0+3]; 
    uRelSq[6] += 0.6324555320336759*uSelf[a0+6]*uSelf[a0+9]-0.6324555320336759*uOther[a0+6]*uSelf[a0+9]-0.6324555320336759*uSelf[a0+6]*uOther[a0+9]+0.6324555320336759*uOther[a0+6]*uOther[a0+9]+0.6324555320336759*uSelf[a0+6]*uSelf[a0+8]-0.6324555320336759*uOther[a0+6]*uSelf[a0+8]-0.6324555320336759*uSelf[a0+6]*uOther[a0+8]+0.6324555320336759*uOther[a0+6]*uOther[a0+8]+0.7071067811865475*uSelf[a0]*uSelf[a0+6]-0.7071067811865475*uOther[a0]*uSelf[a0+6]-0.7071067811865475*uSelf[a0]*uOther[a0+6]+0.7071067811865475*uOther[a0]*uOther[a0+6]+0.7071067811865475*uSelf[a0+4]*uSelf[a0+5]-0.7071067811865475*uOther[a0+4]*uSelf[a0+5]-0.7071067811865475*uSelf[a0+4]*uOther[a0+5]+0.7071067811865475*uOther[a0+4]*uOther[a0+5]+0.7071067811865475*uSelf[a0+2]*uSelf[a0+3]-0.7071067811865475*uOther[a0+2]*uSelf[a0+3]-0.7071067811865475*uSelf[a0+2]*uOther[a0+3]+0.7071067811865475*uOther[a0+2]*uOther[a0+3]; 
    uRelSq[7] += 0.2258769757263128*uSelf[a0+7]*uSelf[a0+7]-0.4517539514526256*uOther[a0+7]*uSelf[a0+7]+0.7071067811865475*uSelf[a0]*uSelf[a0+7]-0.7071067811865475*uOther[a0]*uSelf[a0+7]+0.2258769757263128*uOther[a0+7]*uOther[a0+7]-0.7071067811865475*uSelf[a0]*uOther[a0+7]+0.7071067811865475*uOther[a0]*uOther[a0+7]+0.3162277660168379*uSelf[a0+5]*uSelf[a0+5]-0.6324555320336759*uOther[a0+5]*uSelf[a0+5]+0.3162277660168379*uOther[a0+5]*uOther[a0+5]+0.3162277660168379*uSelf[a0+4]*uSelf[a0+4]-0.6324555320336759*uOther[a0+4]*uSelf[a0+4]+0.3162277660168379*uOther[a0+4]*uOther[a0+4]+0.3162277660168379*uSelf[a0+1]*uSelf[a0+1]-0.6324555320336759*uOther[a0+1]*uSelf[a0+1]+0.3162277660168379*uOther[a0+1]*uOther[a0+1]; 
    uRelSq[8] += 0.2258769757263128*uSelf[a0+8]*uSelf[a0+8]-0.4517539514526256*uOther[a0+8]*uSelf[a0+8]+0.7071067811865475*uSelf[a0]*uSelf[a0+8]-0.7071067811865475*uOther[a0]*uSelf[a0+8]+0.2258769757263128*uOther[a0+8]*uOther[a0+8]-0.7071067811865475*uSelf[a0]*uOther[a0+8]+0.7071067811865475*uOther[a0]*uOther[a0+8]+0.3162277660168379*uSelf[a0+6]*uSelf[a0+6]-0.6324555320336759*uOther[a0+6]*uSelf[a0+6]+0.3162277660168379*uOther[a0+6]*uOther[a0+6]+0.3162277660168379*uSelf[a0+4]*uSelf[a0+4]-0.6324555320336759*uOther[a0+4]*uSelf[a0+4]+0.3162277660168379*uOther[a0+4]*uOther[a0+4]+0.3162277660168379*uSelf[a0+2]*uSelf[a0+2]-0.6324555320336759*uOther[a0+2]*uSelf[a0+2]+0.3162277660168379*uOther[a0+2]*uOther[a0+2]; 
    uRelSq[9] += 0.2258769757263128*uSelf[a0+9]*uSelf[a0+9]-0.4517539514526256*uOther[a0+9]*uSelf[a0+9]+0.7071067811865475*uSelf[a0]*uSelf[a0+9]-0.7071067811865475*uOther[a0]*uSelf[a0+9]+0.2258769757263128*uOther[a0+9]*uOther[a0+9]-0.7071067811865475*uSelf[a0]*uOther[a0+9]+0.7071067811865475*uOther[a0]*uOther[a0+9]+0.3162277660168379*uSelf[a0+6]*uSelf[a0+6]-0.6324555320336759*uOther[a0+6]*uSelf[a0+6]+0.3162277660168379*uOther[a0+6]*uOther[a0+6]+0.3162277660168379*uSelf[a0+5]*uSelf[a0+5]-0.6324555320336759*uOther[a0+5]*uSelf[a0+5]+0.3162277660168379*uOther[a0+5]*uOther[a0+5]+0.3162277660168379*uSelf[a0+3]*uSelf[a0+3]-0.6324555320336759*uOther[a0+3]*uSelf[a0+3]+0.3162277660168379*uOther[a0+3]*uOther[a0+3]; 
  } 
 
  double relKinE[10]; 
  // Zero out array with ((beta+1)/2)*(mSelf+mOther)*(uSelf-uOther) . uRelDmnu. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += betaGreenep1*(0.1767766952966368*uRelDmnu[a0+9]*uSelf[a0+9]*mSelf-0.1767766952966368*uOther[a0+9]*uRelDmnu[a0+9]*mSelf+0.1767766952966368*uRelDmnu[a0+8]*uSelf[a0+8]*mSelf-0.1767766952966368*uOther[a0+8]*uRelDmnu[a0+8]*mSelf+0.1767766952966368*uRelDmnu[a0+7]*uSelf[a0+7]*mSelf-0.1767766952966368*uOther[a0+7]*uRelDmnu[a0+7]*mSelf+0.1767766952966368*uRelDmnu[a0+6]*uSelf[a0+6]*mSelf-0.1767766952966368*uOther[a0+6]*uRelDmnu[a0+6]*mSelf+0.1767766952966368*uRelDmnu[a0+5]*uSelf[a0+5]*mSelf-0.1767766952966368*uOther[a0+5]*uRelDmnu[a0+5]*mSelf+0.1767766952966368*uRelDmnu[a0+4]*uSelf[a0+4]*mSelf-0.1767766952966368*uOther[a0+4]*uRelDmnu[a0+4]*mSelf+0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+3]*mSelf-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+3]*mSelf+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+2]*mSelf-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+2]*mSelf+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+1]*mSelf-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+1]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0]*mSelf+0.1767766952966368*uRelDmnu[a0+9]*uSelf[a0+9]*mOther-0.1767766952966368*uOther[a0+9]*uRelDmnu[a0+9]*mOther+0.1767766952966368*uRelDmnu[a0+8]*uSelf[a0+8]*mOther-0.1767766952966368*uOther[a0+8]*uRelDmnu[a0+8]*mOther+0.1767766952966368*uRelDmnu[a0+7]*uSelf[a0+7]*mOther-0.1767766952966368*uOther[a0+7]*uRelDmnu[a0+7]*mOther+0.1767766952966368*uRelDmnu[a0+6]*uSelf[a0+6]*mOther-0.1767766952966368*uOther[a0+6]*uRelDmnu[a0+6]*mOther+0.1767766952966368*uRelDmnu[a0+5]*uSelf[a0+5]*mOther-0.1767766952966368*uOther[a0+5]*uRelDmnu[a0+5]*mOther+0.1767766952966368*uRelDmnu[a0+4]*uSelf[a0+4]*mOther-0.1767766952966368*uOther[a0+4]*uRelDmnu[a0+4]*mOther+0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+3]*mOther-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+3]*mOther+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+2]*mOther-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+2]*mOther+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+1]*mOther-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+1]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0]*mOther); 
    relKinE[1] += betaGreenep1*(0.1581138830084189*uRelDmnu[a0+1]*uSelf[a0+7]*mSelf+0.1581138830084189*uSelf[a0+1]*uRelDmnu[a0+7]*mSelf-0.1581138830084189*uOther[a0+1]*uRelDmnu[a0+7]*mSelf-0.1581138830084189*uRelDmnu[a0+1]*uOther[a0+7]*mSelf+0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+5]*mSelf+0.1767766952966368*uSelf[a0+3]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uRelDmnu[a0+3]*uOther[a0+5]*mSelf+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+4]*mSelf+0.1767766952966368*uSelf[a0+2]*uRelDmnu[a0+4]*mSelf-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+4]*mSelf-0.1767766952966368*uRelDmnu[a0+2]*uOther[a0+4]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+1]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+1]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+1]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+1]*mSelf+0.1581138830084189*uRelDmnu[a0+1]*uSelf[a0+7]*mOther+0.1581138830084189*uSelf[a0+1]*uRelDmnu[a0+7]*mOther-0.1581138830084189*uOther[a0+1]*uRelDmnu[a0+7]*mOther-0.1581138830084189*uRelDmnu[a0+1]*uOther[a0+7]*mOther+0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+5]*mOther+0.1767766952966368*uSelf[a0+3]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uRelDmnu[a0+3]*uOther[a0+5]*mOther+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+4]*mOther+0.1767766952966368*uSelf[a0+2]*uRelDmnu[a0+4]*mOther-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+4]*mOther-0.1767766952966368*uRelDmnu[a0+2]*uOther[a0+4]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+1]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+1]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+1]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+1]*mOther); 
    relKinE[2] += betaGreenep1*(0.1581138830084189*uRelDmnu[a0+2]*uSelf[a0+8]*mSelf+0.1581138830084189*uSelf[a0+2]*uRelDmnu[a0+8]*mSelf-0.1581138830084189*uOther[a0+2]*uRelDmnu[a0+8]*mSelf-0.1581138830084189*uRelDmnu[a0+2]*uOther[a0+8]*mSelf+0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+6]*mSelf+0.1767766952966368*uSelf[a0+3]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uRelDmnu[a0+3]*uOther[a0+6]*mSelf+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+4]*mSelf+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+4]*mSelf-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+4]*mSelf-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+4]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+2]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+2]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+2]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+2]*mSelf+0.1581138830084189*uRelDmnu[a0+2]*uSelf[a0+8]*mOther+0.1581138830084189*uSelf[a0+2]*uRelDmnu[a0+8]*mOther-0.1581138830084189*uOther[a0+2]*uRelDmnu[a0+8]*mOther-0.1581138830084189*uRelDmnu[a0+2]*uOther[a0+8]*mOther+0.1767766952966368*uRelDmnu[a0+3]*uSelf[a0+6]*mOther+0.1767766952966368*uSelf[a0+3]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uOther[a0+3]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uRelDmnu[a0+3]*uOther[a0+6]*mOther+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+4]*mOther+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+4]*mOther-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+4]*mOther-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+4]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+2]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+2]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+2]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+2]*mOther); 
    relKinE[3] += betaGreenep1*(0.1581138830084189*uRelDmnu[a0+3]*uSelf[a0+9]*mSelf+0.1581138830084189*uSelf[a0+3]*uRelDmnu[a0+9]*mSelf-0.1581138830084189*uOther[a0+3]*uRelDmnu[a0+9]*mSelf-0.1581138830084189*uRelDmnu[a0+3]*uOther[a0+9]*mSelf+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+6]*mSelf+0.1767766952966368*uSelf[a0+2]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uRelDmnu[a0+2]*uOther[a0+6]*mSelf+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+5]*mSelf+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+5]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+3]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+3]*mSelf+0.1581138830084189*uRelDmnu[a0+3]*uSelf[a0+9]*mOther+0.1581138830084189*uSelf[a0+3]*uRelDmnu[a0+9]*mOther-0.1581138830084189*uOther[a0+3]*uRelDmnu[a0+9]*mOther-0.1581138830084189*uRelDmnu[a0+3]*uOther[a0+9]*mOther+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+6]*mOther+0.1767766952966368*uSelf[a0+2]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uRelDmnu[a0+2]*uOther[a0+6]*mOther+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+5]*mOther+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+5]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+3]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+3]*mOther); 
    relKinE[4] += betaGreenep1*(0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+8]*mSelf+0.1581138830084189*uSelf[a0+4]*uRelDmnu[a0+8]*mSelf-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+8]*mSelf-0.1581138830084189*uRelDmnu[a0+4]*uOther[a0+8]*mSelf+0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+7]*mSelf+0.1581138830084189*uSelf[a0+4]*uRelDmnu[a0+7]*mSelf-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+7]*mSelf-0.1581138830084189*uRelDmnu[a0+4]*uOther[a0+7]*mSelf+0.1767766952966368*uRelDmnu[a0+5]*uSelf[a0+6]*mSelf+0.1767766952966368*uSelf[a0+5]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uOther[a0+5]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uRelDmnu[a0+5]*uOther[a0+6]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+4]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+4]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+4]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+4]*mSelf+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+2]*mSelf+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+2]*mSelf-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+2]*mSelf-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+2]*mSelf+0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+8]*mOther+0.1581138830084189*uSelf[a0+4]*uRelDmnu[a0+8]*mOther-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+8]*mOther-0.1581138830084189*uRelDmnu[a0+4]*uOther[a0+8]*mOther+0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+7]*mOther+0.1581138830084189*uSelf[a0+4]*uRelDmnu[a0+7]*mOther-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+7]*mOther-0.1581138830084189*uRelDmnu[a0+4]*uOther[a0+7]*mOther+0.1767766952966368*uRelDmnu[a0+5]*uSelf[a0+6]*mOther+0.1767766952966368*uSelf[a0+5]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uOther[a0+5]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uRelDmnu[a0+5]*uOther[a0+6]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+4]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+4]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+4]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+4]*mOther+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+2]*mOther+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+2]*mOther-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+2]*mOther-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+2]*mOther); 
    relKinE[5] += betaGreenep1*(0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+9]*mSelf+0.1581138830084189*uSelf[a0+5]*uRelDmnu[a0+9]*mSelf-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+9]*mSelf-0.1581138830084189*uRelDmnu[a0+5]*uOther[a0+9]*mSelf+0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+7]*mSelf+0.1581138830084189*uSelf[a0+5]*uRelDmnu[a0+7]*mSelf-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+7]*mSelf-0.1581138830084189*uRelDmnu[a0+5]*uOther[a0+7]*mSelf+0.1767766952966368*uRelDmnu[a0+4]*uSelf[a0+6]*mSelf+0.1767766952966368*uSelf[a0+4]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uOther[a0+4]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uRelDmnu[a0+4]*uOther[a0+6]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+5]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+5]*mSelf+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+3]*mSelf+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+3]*mSelf+0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+9]*mOther+0.1581138830084189*uSelf[a0+5]*uRelDmnu[a0+9]*mOther-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+9]*mOther-0.1581138830084189*uRelDmnu[a0+5]*uOther[a0+9]*mOther+0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+7]*mOther+0.1581138830084189*uSelf[a0+5]*uRelDmnu[a0+7]*mOther-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+7]*mOther-0.1581138830084189*uRelDmnu[a0+5]*uOther[a0+7]*mOther+0.1767766952966368*uRelDmnu[a0+4]*uSelf[a0+6]*mOther+0.1767766952966368*uSelf[a0+4]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uOther[a0+4]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uRelDmnu[a0+4]*uOther[a0+6]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+5]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+5]*mOther+0.1767766952966368*uRelDmnu[a0+1]*uSelf[a0+3]*mOther+0.1767766952966368*uSelf[a0+1]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uOther[a0+1]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uRelDmnu[a0+1]*uOther[a0+3]*mOther); 
    relKinE[6] += betaGreenep1*(0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+9]*mSelf+0.1581138830084189*uSelf[a0+6]*uRelDmnu[a0+9]*mSelf-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+9]*mSelf-0.1581138830084189*uRelDmnu[a0+6]*uOther[a0+9]*mSelf+0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+8]*mSelf+0.1581138830084189*uSelf[a0+6]*uRelDmnu[a0+8]*mSelf-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+8]*mSelf-0.1581138830084189*uRelDmnu[a0+6]*uOther[a0+8]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+6]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+6]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+6]*mSelf+0.1767766952966368*uRelDmnu[a0+4]*uSelf[a0+5]*mSelf+0.1767766952966368*uSelf[a0+4]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uOther[a0+4]*uRelDmnu[a0+5]*mSelf-0.1767766952966368*uRelDmnu[a0+4]*uOther[a0+5]*mSelf+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+3]*mSelf+0.1767766952966368*uSelf[a0+2]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+3]*mSelf-0.1767766952966368*uRelDmnu[a0+2]*uOther[a0+3]*mSelf+0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+9]*mOther+0.1581138830084189*uSelf[a0+6]*uRelDmnu[a0+9]*mOther-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+9]*mOther-0.1581138830084189*uRelDmnu[a0+6]*uOther[a0+9]*mOther+0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+8]*mOther+0.1581138830084189*uSelf[a0+6]*uRelDmnu[a0+8]*mOther-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+8]*mOther-0.1581138830084189*uRelDmnu[a0+6]*uOther[a0+8]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+6]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+6]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+6]*mOther+0.1767766952966368*uRelDmnu[a0+4]*uSelf[a0+5]*mOther+0.1767766952966368*uSelf[a0+4]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uOther[a0+4]*uRelDmnu[a0+5]*mOther-0.1767766952966368*uRelDmnu[a0+4]*uOther[a0+5]*mOther+0.1767766952966368*uRelDmnu[a0+2]*uSelf[a0+3]*mOther+0.1767766952966368*uSelf[a0+2]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uOther[a0+2]*uRelDmnu[a0+3]*mOther-0.1767766952966368*uRelDmnu[a0+2]*uOther[a0+3]*mOther); 
    relKinE[7] += betaGreenep1*(0.1129384878631564*uRelDmnu[a0+7]*uSelf[a0+7]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+7]*mSelf-0.1129384878631564*uOther[a0+7]*uRelDmnu[a0+7]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+7]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+7]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+7]*mSelf+0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+5]*mSelf-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+5]*mSelf+0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+4]*mSelf-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+4]*mSelf+0.1581138830084189*uRelDmnu[a0+1]*uSelf[a0+1]*mSelf-0.1581138830084189*uOther[a0+1]*uRelDmnu[a0+1]*mSelf+0.1129384878631564*uRelDmnu[a0+7]*uSelf[a0+7]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+7]*mOther-0.1129384878631564*uOther[a0+7]*uRelDmnu[a0+7]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+7]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+7]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+7]*mOther+0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+5]*mOther-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+5]*mOther+0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+4]*mOther-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+4]*mOther+0.1581138830084189*uRelDmnu[a0+1]*uSelf[a0+1]*mOther-0.1581138830084189*uOther[a0+1]*uRelDmnu[a0+1]*mOther); 
    relKinE[8] += betaGreenep1*(0.1129384878631564*uRelDmnu[a0+8]*uSelf[a0+8]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+8]*mSelf-0.1129384878631564*uOther[a0+8]*uRelDmnu[a0+8]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+8]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+8]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+8]*mSelf+0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+6]*mSelf-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+6]*mSelf+0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+4]*mSelf-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+4]*mSelf+0.1581138830084189*uRelDmnu[a0+2]*uSelf[a0+2]*mSelf-0.1581138830084189*uOther[a0+2]*uRelDmnu[a0+2]*mSelf+0.1129384878631564*uRelDmnu[a0+8]*uSelf[a0+8]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+8]*mOther-0.1129384878631564*uOther[a0+8]*uRelDmnu[a0+8]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+8]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+8]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+8]*mOther+0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+6]*mOther-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+6]*mOther+0.1581138830084189*uRelDmnu[a0+4]*uSelf[a0+4]*mOther-0.1581138830084189*uOther[a0+4]*uRelDmnu[a0+4]*mOther+0.1581138830084189*uRelDmnu[a0+2]*uSelf[a0+2]*mOther-0.1581138830084189*uOther[a0+2]*uRelDmnu[a0+2]*mOther); 
    relKinE[9] += betaGreenep1*(0.1129384878631564*uRelDmnu[a0+9]*uSelf[a0+9]*mSelf+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+9]*mSelf-0.1129384878631564*uOther[a0+9]*uRelDmnu[a0+9]*mSelf+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+9]*mSelf-0.1767766952966368*uOther[a0]*uRelDmnu[a0+9]*mSelf-0.1767766952966368*uRelDmnu[a0]*uOther[a0+9]*mSelf+0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+6]*mSelf-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+6]*mSelf+0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+5]*mSelf-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+5]*mSelf+0.1581138830084189*uRelDmnu[a0+3]*uSelf[a0+3]*mSelf-0.1581138830084189*uOther[a0+3]*uRelDmnu[a0+3]*mSelf+0.1129384878631564*uRelDmnu[a0+9]*uSelf[a0+9]*mOther+0.1767766952966368*uRelDmnu[a0]*uSelf[a0+9]*mOther-0.1129384878631564*uOther[a0+9]*uRelDmnu[a0+9]*mOther+0.1767766952966368*uSelf[a0]*uRelDmnu[a0+9]*mOther-0.1767766952966368*uOther[a0]*uRelDmnu[a0+9]*mOther-0.1767766952966368*uRelDmnu[a0]*uOther[a0+9]*mOther+0.1581138830084189*uRelDmnu[a0+6]*uSelf[a0+6]*mOther-0.1581138830084189*uOther[a0+6]*uRelDmnu[a0+6]*mOther+0.1581138830084189*uRelDmnu[a0+5]*uSelf[a0+5]*mOther-0.1581138830084189*uOther[a0+5]*uRelDmnu[a0+5]*mOther+0.1581138830084189*uRelDmnu[a0+3]*uSelf[a0+3]*mOther-0.1581138830084189*uOther[a0+3]*uRelDmnu[a0+3]*mOther); 
  } 
 
  double Tdiff[10]; 
  Tdiff[0] = vtSqSelf[0]*mSelf-1.0*vtSqOther[0]*mOther; 
  Tdiff[1] = vtSqSelf[1]*mSelf-1.0*vtSqOther[1]*mOther; 
  Tdiff[2] = vtSqSelf[2]*mSelf-1.0*vtSqOther[2]*mOther; 
  Tdiff[3] = vtSqSelf[3]*mSelf-1.0*vtSqOther[3]*mOther; 
  Tdiff[4] = vtSqSelf[4]*mSelf-1.0*vtSqOther[4]*mOther; 
  Tdiff[5] = vtSqSelf[5]*mSelf-1.0*vtSqOther[5]*mOther; 
  Tdiff[6] = vtSqSelf[6]*mSelf-1.0*vtSqOther[6]*mOther; 
  Tdiff[7] = vtSqSelf[7]*mSelf-1.0*vtSqOther[7]*mOther; 
  Tdiff[8] = vtSqSelf[8]*mSelf-1.0*vtSqOther[8]*mOther; 
  Tdiff[9] = vtSqSelf[9]*mSelf-1.0*vtSqOther[9]*mOther; 
 
  double diffSelf[10]; 
  diffSelf[0] = (0.3535533905932737*m0rOther[9]*relKinE[9]+0.3535533905932737*m0rOther[8]*relKinE[8]+0.3535533905932737*m0rOther[7]*relKinE[7]+0.3535533905932737*m0rOther[6]*relKinE[6]+0.3535533905932737*m0rOther[5]*relKinE[5]+0.3535533905932737*m0rOther[4]*relKinE[4]+0.3535533905932737*m0rOther[3]*relKinE[3]+0.3535533905932737*m0rOther[2]*relKinE[2]+0.3535533905932737*m0rOther[1]*relKinE[1]+0.3535533905932737*m0rOther[0]*relKinE[0])*mnuOther-1.0*uRelSq[0]*mOther+3.0*Tdiff[0]; 
  diffSelf[1] = (0.3162277660168379*m0rOther[1]*relKinE[7]+0.3162277660168379*relKinE[1]*m0rOther[7]+0.3535533905932737*m0rOther[3]*relKinE[5]+0.3535533905932737*relKinE[3]*m0rOther[5]+0.3535533905932737*m0rOther[2]*relKinE[4]+0.3535533905932737*relKinE[2]*m0rOther[4]+0.3535533905932737*m0rOther[0]*relKinE[1]+0.3535533905932737*relKinE[0]*m0rOther[1])*mnuOther-1.0*uRelSq[1]*mOther+3.0*Tdiff[1]; 
  diffSelf[2] = (0.3162277660168379*m0rOther[2]*relKinE[8]+0.3162277660168379*relKinE[2]*m0rOther[8]+0.3535533905932737*m0rOther[3]*relKinE[6]+0.3535533905932737*relKinE[3]*m0rOther[6]+0.3535533905932737*m0rOther[1]*relKinE[4]+0.3535533905932737*relKinE[1]*m0rOther[4]+0.3535533905932737*m0rOther[0]*relKinE[2]+0.3535533905932737*relKinE[0]*m0rOther[2])*mnuOther-1.0*uRelSq[2]*mOther+3.0*Tdiff[2]; 
  diffSelf[3] = (0.3162277660168379*m0rOther[3]*relKinE[9]+0.3162277660168379*relKinE[3]*m0rOther[9]+0.3535533905932737*m0rOther[2]*relKinE[6]+0.3535533905932737*relKinE[2]*m0rOther[6]+0.3535533905932737*m0rOther[1]*relKinE[5]+0.3535533905932737*relKinE[1]*m0rOther[5]+0.3535533905932737*m0rOther[0]*relKinE[3]+0.3535533905932737*relKinE[0]*m0rOther[3])*mnuOther-1.0*uRelSq[3]*mOther+3.0*Tdiff[3]; 
  diffSelf[4] = (0.3162277660168379*m0rOther[4]*relKinE[8]+0.3162277660168379*relKinE[4]*m0rOther[8]+0.3162277660168379*m0rOther[4]*relKinE[7]+0.3162277660168379*relKinE[4]*m0rOther[7]+0.3535533905932737*m0rOther[5]*relKinE[6]+0.3535533905932737*relKinE[5]*m0rOther[6]+0.3535533905932737*m0rOther[0]*relKinE[4]+0.3535533905932737*relKinE[0]*m0rOther[4]+0.3535533905932737*m0rOther[1]*relKinE[2]+0.3535533905932737*relKinE[1]*m0rOther[2])*mnuOther-1.0*uRelSq[4]*mOther+3.0*Tdiff[4]; 
  diffSelf[5] = (0.3162277660168379*m0rOther[5]*relKinE[9]+0.3162277660168379*relKinE[5]*m0rOther[9]+0.3162277660168379*m0rOther[5]*relKinE[7]+0.3162277660168379*relKinE[5]*m0rOther[7]+0.3535533905932737*m0rOther[4]*relKinE[6]+0.3535533905932737*relKinE[4]*m0rOther[6]+0.3535533905932737*m0rOther[0]*relKinE[5]+0.3535533905932737*relKinE[0]*m0rOther[5]+0.3535533905932737*m0rOther[1]*relKinE[3]+0.3535533905932737*relKinE[1]*m0rOther[3])*mnuOther-1.0*uRelSq[5]*mOther+3.0*Tdiff[5]; 
  diffSelf[6] = (0.3162277660168379*m0rOther[6]*relKinE[9]+0.3162277660168379*relKinE[6]*m0rOther[9]+0.3162277660168379*m0rOther[6]*relKinE[8]+0.3162277660168379*relKinE[6]*m0rOther[8]+0.3535533905932737*m0rOther[0]*relKinE[6]+0.3535533905932737*relKinE[0]*m0rOther[6]+0.3535533905932737*m0rOther[4]*relKinE[5]+0.3535533905932737*relKinE[4]*m0rOther[5]+0.3535533905932737*m0rOther[2]*relKinE[3]+0.3535533905932737*relKinE[2]*m0rOther[3])*mnuOther-1.0*uRelSq[6]*mOther+3.0*Tdiff[6]; 
  diffSelf[7] = (0.2258769757263128*m0rOther[7]*relKinE[7]+0.3535533905932737*m0rOther[0]*relKinE[7]+0.3535533905932737*relKinE[0]*m0rOther[7]+0.3162277660168379*m0rOther[5]*relKinE[5]+0.3162277660168379*m0rOther[4]*relKinE[4]+0.3162277660168379*m0rOther[1]*relKinE[1])*mnuOther-1.0*uRelSq[7]*mOther+3.0*Tdiff[7]; 
  diffSelf[8] = (0.2258769757263128*m0rOther[8]*relKinE[8]+0.3535533905932737*m0rOther[0]*relKinE[8]+0.3535533905932737*relKinE[0]*m0rOther[8]+0.3162277660168379*m0rOther[6]*relKinE[6]+0.3162277660168379*m0rOther[4]*relKinE[4]+0.3162277660168379*m0rOther[2]*relKinE[2])*mnuOther-1.0*uRelSq[8]*mOther+3.0*Tdiff[8]; 
  diffSelf[9] = (0.2258769757263128*m0rOther[9]*relKinE[9]+0.3535533905932737*m0rOther[0]*relKinE[9]+0.3535533905932737*relKinE[0]*m0rOther[9]+0.3162277660168379*m0rOther[6]*relKinE[6]+0.3162277660168379*m0rOther[5]*relKinE[5]+0.3162277660168379*m0rOther[3]*relKinE[3])*mnuOther-1.0*uRelSq[9]*mOther+3.0*Tdiff[9]; 
 
  // Divide diffSelf by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffSelf[0],diffSelf[1],diffSelf[2],diffSelf[3],diffSelf[4],diffSelf[5],diffSelf[6],diffSelf[7],diffSelf[8],diffSelf[9]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaSelf[10]; 
  Eigen::Map<VectorXd>(vtSqDeltaSelf,10,1) = dataDiv->u_S; 
 
  double diffOther[10]; 
  diffOther[0] = (0.3535533905932737*m0rSelf[9]*relKinE[9]+0.3535533905932737*m0rSelf[8]*relKinE[8]+0.3535533905932737*m0rSelf[7]*relKinE[7]+0.3535533905932737*m0rSelf[6]*relKinE[6]+0.3535533905932737*m0rSelf[5]*relKinE[5]+0.3535533905932737*m0rSelf[4]*relKinE[4]+0.3535533905932737*m0rSelf[3]*relKinE[3]+0.3535533905932737*m0rSelf[2]*relKinE[2]+0.3535533905932737*m0rSelf[1]*relKinE[1]+0.3535533905932737*m0rSelf[0]*relKinE[0])*mnuSelf-1.0*uRelSq[0]*mSelf-3.0*Tdiff[0]; 
  diffOther[1] = (0.3162277660168379*m0rSelf[1]*relKinE[7]+0.3162277660168379*relKinE[1]*m0rSelf[7]+0.3535533905932737*m0rSelf[3]*relKinE[5]+0.3535533905932737*relKinE[3]*m0rSelf[5]+0.3535533905932737*m0rSelf[2]*relKinE[4]+0.3535533905932737*relKinE[2]*m0rSelf[4]+0.3535533905932737*m0rSelf[0]*relKinE[1]+0.3535533905932737*relKinE[0]*m0rSelf[1])*mnuSelf-1.0*uRelSq[1]*mSelf-3.0*Tdiff[1]; 
  diffOther[2] = (0.3162277660168379*m0rSelf[2]*relKinE[8]+0.3162277660168379*relKinE[2]*m0rSelf[8]+0.3535533905932737*m0rSelf[3]*relKinE[6]+0.3535533905932737*relKinE[3]*m0rSelf[6]+0.3535533905932737*m0rSelf[1]*relKinE[4]+0.3535533905932737*relKinE[1]*m0rSelf[4]+0.3535533905932737*m0rSelf[0]*relKinE[2]+0.3535533905932737*relKinE[0]*m0rSelf[2])*mnuSelf-1.0*uRelSq[2]*mSelf-3.0*Tdiff[2]; 
  diffOther[3] = (0.3162277660168379*m0rSelf[3]*relKinE[9]+0.3162277660168379*relKinE[3]*m0rSelf[9]+0.3535533905932737*m0rSelf[2]*relKinE[6]+0.3535533905932737*relKinE[2]*m0rSelf[6]+0.3535533905932737*m0rSelf[1]*relKinE[5]+0.3535533905932737*relKinE[1]*m0rSelf[5]+0.3535533905932737*m0rSelf[0]*relKinE[3]+0.3535533905932737*relKinE[0]*m0rSelf[3])*mnuSelf-1.0*uRelSq[3]*mSelf-3.0*Tdiff[3]; 
  diffOther[4] = (0.3162277660168379*m0rSelf[4]*relKinE[8]+0.3162277660168379*relKinE[4]*m0rSelf[8]+0.3162277660168379*m0rSelf[4]*relKinE[7]+0.3162277660168379*relKinE[4]*m0rSelf[7]+0.3535533905932737*m0rSelf[5]*relKinE[6]+0.3535533905932737*relKinE[5]*m0rSelf[6]+0.3535533905932737*m0rSelf[0]*relKinE[4]+0.3535533905932737*relKinE[0]*m0rSelf[4]+0.3535533905932737*m0rSelf[1]*relKinE[2]+0.3535533905932737*relKinE[1]*m0rSelf[2])*mnuSelf-1.0*uRelSq[4]*mSelf-3.0*Tdiff[4]; 
  diffOther[5] = (0.3162277660168379*m0rSelf[5]*relKinE[9]+0.3162277660168379*relKinE[5]*m0rSelf[9]+0.3162277660168379*m0rSelf[5]*relKinE[7]+0.3162277660168379*relKinE[5]*m0rSelf[7]+0.3535533905932737*m0rSelf[4]*relKinE[6]+0.3535533905932737*relKinE[4]*m0rSelf[6]+0.3535533905932737*m0rSelf[0]*relKinE[5]+0.3535533905932737*relKinE[0]*m0rSelf[5]+0.3535533905932737*m0rSelf[1]*relKinE[3]+0.3535533905932737*relKinE[1]*m0rSelf[3])*mnuSelf-1.0*uRelSq[5]*mSelf-3.0*Tdiff[5]; 
  diffOther[6] = (0.3162277660168379*m0rSelf[6]*relKinE[9]+0.3162277660168379*relKinE[6]*m0rSelf[9]+0.3162277660168379*m0rSelf[6]*relKinE[8]+0.3162277660168379*relKinE[6]*m0rSelf[8]+0.3535533905932737*m0rSelf[0]*relKinE[6]+0.3535533905932737*relKinE[0]*m0rSelf[6]+0.3535533905932737*m0rSelf[4]*relKinE[5]+0.3535533905932737*relKinE[4]*m0rSelf[5]+0.3535533905932737*m0rSelf[2]*relKinE[3]+0.3535533905932737*relKinE[2]*m0rSelf[3])*mnuSelf-1.0*uRelSq[6]*mSelf-3.0*Tdiff[6]; 
  diffOther[7] = (0.2258769757263128*m0rSelf[7]*relKinE[7]+0.3535533905932737*m0rSelf[0]*relKinE[7]+0.3535533905932737*relKinE[0]*m0rSelf[7]+0.3162277660168379*m0rSelf[5]*relKinE[5]+0.3162277660168379*m0rSelf[4]*relKinE[4]+0.3162277660168379*m0rSelf[1]*relKinE[1])*mnuSelf-1.0*uRelSq[7]*mSelf-3.0*Tdiff[7]; 
  diffOther[8] = (0.2258769757263128*m0rSelf[8]*relKinE[8]+0.3535533905932737*m0rSelf[0]*relKinE[8]+0.3535533905932737*relKinE[0]*m0rSelf[8]+0.3162277660168379*m0rSelf[6]*relKinE[6]+0.3162277660168379*m0rSelf[4]*relKinE[4]+0.3162277660168379*m0rSelf[2]*relKinE[2])*mnuSelf-1.0*uRelSq[8]*mSelf-3.0*Tdiff[8]; 
  diffOther[9] = (0.2258769757263128*m0rSelf[9]*relKinE[9]+0.3535533905932737*m0rSelf[0]*relKinE[9]+0.3535533905932737*relKinE[0]*m0rSelf[9]+0.3162277660168379*m0rSelf[6]*relKinE[6]+0.3162277660168379*m0rSelf[5]*relKinE[5]+0.3162277660168379*m0rSelf[3]*relKinE[3])*mnuSelf-1.0*uRelSq[9]*mSelf-3.0*Tdiff[9]; 
 
  // Divide diffOther by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffOther[0],diffOther[1],diffOther[2],diffOther[3],diffOther[4],diffOther[5],diffOther[6],diffOther[7],diffOther[8],diffOther[9]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaOther[10]; 
  Eigen::Map<VectorXd>(vtSqDeltaOther,10,1) = dataDiv->u_S; 
 
  // ... Cross-thermal speeds (squared) ... // 
  double deltaFacOther = (0.6666666666666666*mnuOther)/(mSelf+mOther); 
  vtSqCrossSelf[0] = (-0.3535533905932737*m0rOther[9]*vtSqDeltaSelf[9]*deltaFacOther)-0.3535533905932737*m0rOther[8]*vtSqDeltaSelf[8]*deltaFacOther-0.3535533905932737*m0rOther[7]*vtSqDeltaSelf[7]*deltaFacOther-0.3535533905932737*m0rOther[6]*vtSqDeltaSelf[6]*deltaFacOther-0.3535533905932737*m0rOther[5]*vtSqDeltaSelf[5]*deltaFacOther-0.3535533905932737*m0rOther[4]*vtSqDeltaSelf[4]*deltaFacOther-0.3535533905932737*m0rOther[3]*vtSqDeltaSelf[3]*deltaFacOther-0.3535533905932737*m0rOther[2]*vtSqDeltaSelf[2]*deltaFacOther-0.3535533905932737*m0rOther[1]*vtSqDeltaSelf[1]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[0]*deltaFacOther+vtSqSelf[0]; 
  vtSqCrossSelf[1] = (-0.3162277660168379*m0rOther[1]*vtSqDeltaSelf[7]*deltaFacOther)-0.3162277660168379*vtSqDeltaSelf[1]*m0rOther[7]*deltaFacOther-0.3535533905932737*m0rOther[3]*vtSqDeltaSelf[5]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[3]*m0rOther[5]*deltaFacOther-0.3535533905932737*m0rOther[2]*vtSqDeltaSelf[4]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[2]*m0rOther[4]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[1]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[1]*deltaFacOther+vtSqSelf[1]; 
  vtSqCrossSelf[2] = (-0.3162277660168379*m0rOther[2]*vtSqDeltaSelf[8]*deltaFacOther)-0.3162277660168379*vtSqDeltaSelf[2]*m0rOther[8]*deltaFacOther-0.3535533905932737*m0rOther[3]*vtSqDeltaSelf[6]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[3]*m0rOther[6]*deltaFacOther-0.3535533905932737*m0rOther[1]*vtSqDeltaSelf[4]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[1]*m0rOther[4]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[2]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[2]*deltaFacOther+vtSqSelf[2]; 
  vtSqCrossSelf[3] = (-0.3162277660168379*m0rOther[3]*vtSqDeltaSelf[9]*deltaFacOther)-0.3162277660168379*vtSqDeltaSelf[3]*m0rOther[9]*deltaFacOther-0.3535533905932737*m0rOther[2]*vtSqDeltaSelf[6]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[2]*m0rOther[6]*deltaFacOther-0.3535533905932737*m0rOther[1]*vtSqDeltaSelf[5]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[1]*m0rOther[5]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[3]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[3]*deltaFacOther+vtSqSelf[3]; 
  vtSqCrossSelf[4] = (-0.3162277660168379*m0rOther[4]*vtSqDeltaSelf[8]*deltaFacOther)-0.3162277660168379*vtSqDeltaSelf[4]*m0rOther[8]*deltaFacOther-0.3162277660168379*m0rOther[4]*vtSqDeltaSelf[7]*deltaFacOther-0.3162277660168379*vtSqDeltaSelf[4]*m0rOther[7]*deltaFacOther-0.3535533905932737*m0rOther[5]*vtSqDeltaSelf[6]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[5]*m0rOther[6]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[4]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[4]*deltaFacOther-0.3535533905932737*m0rOther[1]*vtSqDeltaSelf[2]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[1]*m0rOther[2]*deltaFacOther+vtSqSelf[4]; 
  vtSqCrossSelf[5] = (-0.3162277660168379*m0rOther[5]*vtSqDeltaSelf[9]*deltaFacOther)-0.3162277660168379*vtSqDeltaSelf[5]*m0rOther[9]*deltaFacOther-0.3162277660168379*m0rOther[5]*vtSqDeltaSelf[7]*deltaFacOther-0.3162277660168379*vtSqDeltaSelf[5]*m0rOther[7]*deltaFacOther-0.3535533905932737*m0rOther[4]*vtSqDeltaSelf[6]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[4]*m0rOther[6]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[5]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[5]*deltaFacOther-0.3535533905932737*m0rOther[1]*vtSqDeltaSelf[3]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[1]*m0rOther[3]*deltaFacOther+vtSqSelf[5]; 
  vtSqCrossSelf[6] = (-0.3162277660168379*m0rOther[6]*vtSqDeltaSelf[9]*deltaFacOther)-0.3162277660168379*vtSqDeltaSelf[6]*m0rOther[9]*deltaFacOther-0.3162277660168379*m0rOther[6]*vtSqDeltaSelf[8]*deltaFacOther-0.3162277660168379*vtSqDeltaSelf[6]*m0rOther[8]*deltaFacOther-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[6]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[6]*deltaFacOther-0.3535533905932737*m0rOther[4]*vtSqDeltaSelf[5]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[4]*m0rOther[5]*deltaFacOther-0.3535533905932737*m0rOther[2]*vtSqDeltaSelf[3]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[2]*m0rOther[3]*deltaFacOther+vtSqSelf[6]; 
  vtSqCrossSelf[7] = (-0.2258769757263128*m0rOther[7]*vtSqDeltaSelf[7]*deltaFacOther)-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[7]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[7]*deltaFacOther-0.3162277660168379*m0rOther[5]*vtSqDeltaSelf[5]*deltaFacOther-0.3162277660168379*m0rOther[4]*vtSqDeltaSelf[4]*deltaFacOther-0.3162277660168379*m0rOther[1]*vtSqDeltaSelf[1]*deltaFacOther+vtSqSelf[7]; 
  vtSqCrossSelf[8] = (-0.2258769757263128*m0rOther[8]*vtSqDeltaSelf[8]*deltaFacOther)-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[8]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[8]*deltaFacOther-0.3162277660168379*m0rOther[6]*vtSqDeltaSelf[6]*deltaFacOther-0.3162277660168379*m0rOther[4]*vtSqDeltaSelf[4]*deltaFacOther-0.3162277660168379*m0rOther[2]*vtSqDeltaSelf[2]*deltaFacOther+vtSqSelf[8]; 
  vtSqCrossSelf[9] = (-0.2258769757263128*m0rOther[9]*vtSqDeltaSelf[9]*deltaFacOther)-0.3535533905932737*m0rOther[0]*vtSqDeltaSelf[9]*deltaFacOther-0.3535533905932737*vtSqDeltaSelf[0]*m0rOther[9]*deltaFacOther-0.3162277660168379*m0rOther[6]*vtSqDeltaSelf[6]*deltaFacOther-0.3162277660168379*m0rOther[5]*vtSqDeltaSelf[5]*deltaFacOther-0.3162277660168379*m0rOther[3]*vtSqDeltaSelf[3]*deltaFacOther+vtSqSelf[9]; 
 
  double deltaFacSelf = (0.6666666666666666*mnuSelf)/(mSelf+mOther); 
  vtSqCrossOther[0] = (-0.3535533905932737*m0rSelf[9]*vtSqDeltaOther[9]*deltaFacSelf)-0.3535533905932737*m0rSelf[8]*vtSqDeltaOther[8]*deltaFacSelf-0.3535533905932737*m0rSelf[7]*vtSqDeltaOther[7]*deltaFacSelf-0.3535533905932737*m0rSelf[6]*vtSqDeltaOther[6]*deltaFacSelf-0.3535533905932737*m0rSelf[5]*vtSqDeltaOther[5]*deltaFacSelf-0.3535533905932737*m0rSelf[4]*vtSqDeltaOther[4]*deltaFacSelf-0.3535533905932737*m0rSelf[3]*vtSqDeltaOther[3]*deltaFacSelf-0.3535533905932737*m0rSelf[2]*vtSqDeltaOther[2]*deltaFacSelf-0.3535533905932737*m0rSelf[1]*vtSqDeltaOther[1]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[0]*deltaFacSelf+vtSqOther[0]; 
  vtSqCrossOther[1] = (-0.3162277660168379*m0rSelf[1]*vtSqDeltaOther[7]*deltaFacSelf)-0.3162277660168379*vtSqDeltaOther[1]*m0rSelf[7]*deltaFacSelf-0.3535533905932737*m0rSelf[3]*vtSqDeltaOther[5]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[3]*m0rSelf[5]*deltaFacSelf-0.3535533905932737*m0rSelf[2]*vtSqDeltaOther[4]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[2]*m0rSelf[4]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[1]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[1]*deltaFacSelf+vtSqOther[1]; 
  vtSqCrossOther[2] = (-0.3162277660168379*m0rSelf[2]*vtSqDeltaOther[8]*deltaFacSelf)-0.3162277660168379*vtSqDeltaOther[2]*m0rSelf[8]*deltaFacSelf-0.3535533905932737*m0rSelf[3]*vtSqDeltaOther[6]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[3]*m0rSelf[6]*deltaFacSelf-0.3535533905932737*m0rSelf[1]*vtSqDeltaOther[4]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[1]*m0rSelf[4]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[2]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[2]*deltaFacSelf+vtSqOther[2]; 
  vtSqCrossOther[3] = (-0.3162277660168379*m0rSelf[3]*vtSqDeltaOther[9]*deltaFacSelf)-0.3162277660168379*vtSqDeltaOther[3]*m0rSelf[9]*deltaFacSelf-0.3535533905932737*m0rSelf[2]*vtSqDeltaOther[6]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[2]*m0rSelf[6]*deltaFacSelf-0.3535533905932737*m0rSelf[1]*vtSqDeltaOther[5]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[1]*m0rSelf[5]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[3]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[3]*deltaFacSelf+vtSqOther[3]; 
  vtSqCrossOther[4] = (-0.3162277660168379*m0rSelf[4]*vtSqDeltaOther[8]*deltaFacSelf)-0.3162277660168379*vtSqDeltaOther[4]*m0rSelf[8]*deltaFacSelf-0.3162277660168379*m0rSelf[4]*vtSqDeltaOther[7]*deltaFacSelf-0.3162277660168379*vtSqDeltaOther[4]*m0rSelf[7]*deltaFacSelf-0.3535533905932737*m0rSelf[5]*vtSqDeltaOther[6]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[5]*m0rSelf[6]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[4]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[4]*deltaFacSelf-0.3535533905932737*m0rSelf[1]*vtSqDeltaOther[2]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[1]*m0rSelf[2]*deltaFacSelf+vtSqOther[4]; 
  vtSqCrossOther[5] = (-0.3162277660168379*m0rSelf[5]*vtSqDeltaOther[9]*deltaFacSelf)-0.3162277660168379*vtSqDeltaOther[5]*m0rSelf[9]*deltaFacSelf-0.3162277660168379*m0rSelf[5]*vtSqDeltaOther[7]*deltaFacSelf-0.3162277660168379*vtSqDeltaOther[5]*m0rSelf[7]*deltaFacSelf-0.3535533905932737*m0rSelf[4]*vtSqDeltaOther[6]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[4]*m0rSelf[6]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[5]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[5]*deltaFacSelf-0.3535533905932737*m0rSelf[1]*vtSqDeltaOther[3]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[1]*m0rSelf[3]*deltaFacSelf+vtSqOther[5]; 
  vtSqCrossOther[6] = (-0.3162277660168379*m0rSelf[6]*vtSqDeltaOther[9]*deltaFacSelf)-0.3162277660168379*vtSqDeltaOther[6]*m0rSelf[9]*deltaFacSelf-0.3162277660168379*m0rSelf[6]*vtSqDeltaOther[8]*deltaFacSelf-0.3162277660168379*vtSqDeltaOther[6]*m0rSelf[8]*deltaFacSelf-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[6]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[6]*deltaFacSelf-0.3535533905932737*m0rSelf[4]*vtSqDeltaOther[5]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[4]*m0rSelf[5]*deltaFacSelf-0.3535533905932737*m0rSelf[2]*vtSqDeltaOther[3]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[2]*m0rSelf[3]*deltaFacSelf+vtSqOther[6]; 
  vtSqCrossOther[7] = (-0.2258769757263128*m0rSelf[7]*vtSqDeltaOther[7]*deltaFacSelf)-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[7]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[7]*deltaFacSelf-0.3162277660168379*m0rSelf[5]*vtSqDeltaOther[5]*deltaFacSelf-0.3162277660168379*m0rSelf[4]*vtSqDeltaOther[4]*deltaFacSelf-0.3162277660168379*m0rSelf[1]*vtSqDeltaOther[1]*deltaFacSelf+vtSqOther[7]; 
  vtSqCrossOther[8] = (-0.2258769757263128*m0rSelf[8]*vtSqDeltaOther[8]*deltaFacSelf)-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[8]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[8]*deltaFacSelf-0.3162277660168379*m0rSelf[6]*vtSqDeltaOther[6]*deltaFacSelf-0.3162277660168379*m0rSelf[4]*vtSqDeltaOther[4]*deltaFacSelf-0.3162277660168379*m0rSelf[2]*vtSqDeltaOther[2]*deltaFacSelf+vtSqOther[8]; 
  vtSqCrossOther[9] = (-0.2258769757263128*m0rSelf[9]*vtSqDeltaOther[9]*deltaFacSelf)-0.3535533905932737*m0rSelf[0]*vtSqDeltaOther[9]*deltaFacSelf-0.3535533905932737*vtSqDeltaOther[0]*m0rSelf[9]*deltaFacSelf-0.3162277660168379*m0rSelf[6]*vtSqDeltaOther[6]*deltaFacSelf-0.3162277660168379*m0rSelf[5]*vtSqDeltaOther[5]*deltaFacSelf-0.3162277660168379*m0rSelf[3]*vtSqDeltaOther[3]*deltaFacSelf+vtSqOther[9]; 
 
} 
 
