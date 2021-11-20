#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmBGKCrossPrimMoments1x2vSer_P1(binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *m0Other, const double *uOther, const double *vtSqOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:     free parameter beta+1. This has to be >0. 
  // m, nu:            mass and collisionality. 
  // m0:               zeroth moment of the distribution function. 
  // u,vtSq:           self primitive moments: mean flow velocity and thermal speed squared. 
  // uCross,vtSqCross: cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (0.7071067811865475*m0Self[0]-1.224744871391589*m0Self[1] < 0) { 
    cellAvg = true;
  }
  if (1.224744871391589*m0Self[1]+0.7071067811865475*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[2]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
  } 
 
  if (0.7071067811865475*m0Other[0]-1.224744871391589*m0Other[1] < 0) { 
    cellAvg = true;
  }
  if (1.224744871391589*m0Other[1]+0.7071067811865475*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[2]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
  } 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
 
  double uRelDmnu[4]; 
 
  // ... Divide (uSelfX-uOtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[0] = uSelf[0]-1.0*uOther[0]; 
  uRelDmnu[1] = uSelf[1]-1.0*uOther[1]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S.setZero(); 
  dataDiv->AEM_S(0,0) = 0.7071067811865475*m0rSelf[0]*mnuSelf+0.7071067811865475*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.7071067811865475*m0rSelf[1]*mnuSelf+0.7071067811865475*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.7071067811865475*m0rSelf[1]*mnuSelf+0.7071067811865475*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.7071067811865475*m0rSelf[0]*mnuSelf+0.7071067811865475*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[0],uRelDmnu[1]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+0,2,1) = dataDiv->u_S; 
 
  // ... Component 1 of cross-velocity of this species ... // 
  uCrossSelf[0] = ((-0.7071067811865475*m0rOther[1]*uRelDmnu[1])-0.7071067811865475*m0rOther[0]*uRelDmnu[0])*betaGreenep1*mnuOther+uSelf[0]; 
  uCrossSelf[1] = ((-0.7071067811865475*m0rOther[0]*uRelDmnu[1])-0.7071067811865475*uRelDmnu[0]*m0rOther[1])*betaGreenep1*mnuOther+uSelf[1]; 
 
  // ... Component 1 of cross-velocity of the other species ... // 
  uCrossOther[0] = (0.7071067811865475*m0rSelf[1]*uRelDmnu[1]+0.7071067811865475*m0rSelf[0]*uRelDmnu[0])*betaGreenep1*mnuSelf+uOther[0]; 
  uCrossOther[1] = (0.7071067811865475*m0rSelf[0]*uRelDmnu[1]+0.7071067811865475*uRelDmnu[0]*m0rSelf[1])*betaGreenep1*mnuSelf+uOther[1]; 
 
  // ... Divide (uSelfY-uOtherY)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[2] = uSelf[2]-1.0*uOther[2]; 
  uRelDmnu[3] = uSelf[3]-1.0*uOther[3]; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[2],uRelDmnu[3]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+2,2,1) = dataDiv->u_S; 
 
  // ... Component 2 of cross-velocity of this species ... // 
  uCrossSelf[2] = ((-0.7071067811865475*m0rOther[1]*uRelDmnu[3])-0.7071067811865475*m0rOther[0]*uRelDmnu[2])*betaGreenep1*mnuOther+uSelf[2]; 
  uCrossSelf[3] = ((-0.7071067811865475*m0rOther[0]*uRelDmnu[3])-0.7071067811865475*m0rOther[1]*uRelDmnu[2])*betaGreenep1*mnuOther+uSelf[3]; 
 
  // ... Component 2 of cross-velocity of the other species ... // 
  uCrossOther[2] = (0.7071067811865475*m0rSelf[1]*uRelDmnu[3]+0.7071067811865475*m0rSelf[0]*uRelDmnu[2])*betaGreenep1*mnuSelf+uOther[2]; 
  uCrossOther[3] = (0.7071067811865475*m0rSelf[0]*uRelDmnu[3]+0.7071067811865475*m0rSelf[1]*uRelDmnu[2])*betaGreenep1*mnuSelf+uOther[3]; 
 
  double uRelSq[2]; 
  // Zero out array with dot product of uSelf-uOther with itself. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    uRelSq[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    uRelSq[0] += 0.7071067811865475*uSelf[a0+1]*uSelf[a0+1]-1.414213562373095*uOther[a0+1]*uSelf[a0+1]+0.7071067811865475*uOther[a0+1]*uOther[a0+1]+0.7071067811865475*uSelf[a0]*uSelf[a0]-1.414213562373095*uOther[a0]*uSelf[a0]+0.7071067811865475*uOther[a0]*uOther[a0]; 
    uRelSq[1] += 1.414213562373095*uSelf[a0]*uSelf[a0+1]-1.414213562373095*uOther[a0]*uSelf[a0+1]-1.414213562373095*uSelf[a0]*uOther[a0+1]+1.414213562373095*uOther[a0]*uOther[a0+1]; 
  } 
 
  double relKinE[2]; 
  // Zero out array with ((beta+1)/2)*(mSelf+mOther)*(uSelf-uOther) . uRelDmnu. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += betaGreenep1*(0.3535533905932737*uRelDmnu[a0+1]*uSelf[a0+1]*mSelf-0.3535533905932737*uOther[a0+1]*uRelDmnu[a0+1]*mSelf+0.3535533905932737*uRelDmnu[a0]*uSelf[a0]*mSelf-0.3535533905932737*uOther[a0]*uRelDmnu[a0]*mSelf+0.3535533905932737*uRelDmnu[a0+1]*uSelf[a0+1]*mOther-0.3535533905932737*uOther[a0+1]*uRelDmnu[a0+1]*mOther+0.3535533905932737*uRelDmnu[a0]*uSelf[a0]*mOther-0.3535533905932737*uOther[a0]*uRelDmnu[a0]*mOther); 
    relKinE[1] += betaGreenep1*(0.3535533905932737*uRelDmnu[a0]*uSelf[a0+1]*mSelf+0.3535533905932737*uSelf[a0]*uRelDmnu[a0+1]*mSelf-0.3535533905932737*uOther[a0]*uRelDmnu[a0+1]*mSelf-0.3535533905932737*uRelDmnu[a0]*uOther[a0+1]*mSelf+0.3535533905932737*uRelDmnu[a0]*uSelf[a0+1]*mOther+0.3535533905932737*uSelf[a0]*uRelDmnu[a0+1]*mOther-0.3535533905932737*uOther[a0]*uRelDmnu[a0+1]*mOther-0.3535533905932737*uRelDmnu[a0]*uOther[a0+1]*mOther); 
  } 
 
  double Tdiff[2]; 
  Tdiff[0] = vtSqSelf[0]*mSelf-1.0*vtSqOther[0]*mOther; 
  Tdiff[1] = vtSqSelf[1]*mSelf-1.0*vtSqOther[1]*mOther; 
 
  double diffSelf[2]; 
  diffSelf[0] = (0.7071067811865475*m0rOther[1]*relKinE[1]+0.7071067811865475*m0rOther[0]*relKinE[0])*mnuOther-1.0*uRelSq[0]*mOther+2.0*Tdiff[0]; 
  diffSelf[1] = (0.7071067811865475*m0rOther[0]*relKinE[1]+0.7071067811865475*relKinE[0]*m0rOther[1])*mnuOther-1.0*uRelSq[1]*mOther+2.0*Tdiff[1]; 
 
  // Divide diffSelf by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffSelf[0],diffSelf[1]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaSelf[2]; 
  Eigen::Map<VectorXd>(vtSqDeltaSelf,2,1) = dataDiv->u_S; 
 
  double diffOther[2]; 
  diffOther[0] = (0.7071067811865475*m0rSelf[1]*relKinE[1]+0.7071067811865475*m0rSelf[0]*relKinE[0])*mnuSelf-1.0*uRelSq[0]*mSelf-2.0*Tdiff[0]; 
  diffOther[1] = (0.7071067811865475*m0rSelf[0]*relKinE[1]+0.7071067811865475*relKinE[0]*m0rSelf[1])*mnuSelf-1.0*uRelSq[1]*mSelf-2.0*Tdiff[1]; 
 
  // Divide diffOther by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffOther[0],diffOther[1]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaOther[2]; 
  Eigen::Map<VectorXd>(vtSqDeltaOther,2,1) = dataDiv->u_S; 
 
  // ... Cross-thermal speeds (squared) ... // 
  double deltaFacOther = mnuOther/(mSelf+mOther); 
  vtSqCrossSelf[0] = (-0.7071067811865475*m0rOther[1]*vtSqDeltaSelf[1]*deltaFacOther)-0.7071067811865475*m0rOther[0]*vtSqDeltaSelf[0]*deltaFacOther+vtSqSelf[0]; 
  vtSqCrossSelf[1] = (-0.7071067811865475*m0rOther[0]*vtSqDeltaSelf[1]*deltaFacOther)-0.7071067811865475*vtSqDeltaSelf[0]*m0rOther[1]*deltaFacOther+vtSqSelf[1]; 
 
  double deltaFacSelf = mnuSelf/(mSelf+mOther); 
  vtSqCrossOther[0] = (-0.7071067811865475*m0rSelf[1]*vtSqDeltaOther[1]*deltaFacSelf)-0.7071067811865475*m0rSelf[0]*vtSqDeltaOther[0]*deltaFacSelf+vtSqOther[0]; 
  vtSqCrossOther[1] = (-0.7071067811865475*m0rSelf[0]*vtSqDeltaOther[1]*deltaFacSelf)-0.7071067811865475*vtSqDeltaOther[0]*m0rSelf[1]*deltaFacSelf+vtSqOther[1]; 
 
} 
 
