#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkBGKCrossPrimMoments2x2vMax_P1(binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *m0Other, const double *uOther, const double *vtSqOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:     free parameter beta+1. This has to be >0. 
  // m, nu:            mass and collisionality. 
  // m0:               zeroth moment of the distribution function. 
  // u,vtSq:           self primitive moments: mean flow velocity and thermal speed squared. 
  // uCross,vtSqCross: cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if ((-0.8660254037844386*m0Self[2])-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Self[2])+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[3]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
  } 
 
  if ((-0.8660254037844386*m0Other[2])-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if ((-0.8660254037844386*m0Other[2])+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[3]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
  } 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
 
  double uRelDmnu[3]; 
 
  // ... Divide (uSelfX-uOtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[0] = uSelf[0]-1.0*uOther[0]; 
  uRelDmnu[1] = uSelf[1]-1.0*uOther[1]; 
  uRelDmnu[2] = uSelf[2]-1.0*uOther[2]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(3,3); 
  dataDiv->AEM_S(0,0) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(0,2) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,2) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[0],uRelDmnu[1],uRelDmnu[2]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+0,3,1) = dataDiv->u_S; 
 
  // ... Component 1 of cross-velocity of this species ... // 
  uCrossSelf[0] = ((-0.5*m0rOther[2]*uRelDmnu[2])-0.5*m0rOther[1]*uRelDmnu[1]-0.5*m0rOther[0]*uRelDmnu[0])*betaGreenep1*mnuOther+uSelf[0]; 
  uCrossSelf[1] = ((-0.5*m0rOther[0]*uRelDmnu[1])-0.5*uRelDmnu[0]*m0rOther[1])*betaGreenep1*mnuOther+uSelf[1]; 
  uCrossSelf[2] = ((-0.5*m0rOther[0]*uRelDmnu[2])-0.5*uRelDmnu[0]*m0rOther[2])*betaGreenep1*mnuOther+uSelf[2]; 
 
  // ... Component 1 of cross-velocity of the other species ... // 
  uCrossOther[0] = (0.5*m0rSelf[2]*uRelDmnu[2]+0.5*m0rSelf[1]*uRelDmnu[1]+0.5*m0rSelf[0]*uRelDmnu[0])*betaGreenep1*mnuSelf+uOther[0]; 
  uCrossOther[1] = (0.5*m0rSelf[0]*uRelDmnu[1]+0.5*uRelDmnu[0]*m0rSelf[1])*betaGreenep1*mnuSelf+uOther[1]; 
  uCrossOther[2] = (0.5*m0rSelf[0]*uRelDmnu[2]+0.5*uRelDmnu[0]*m0rSelf[2])*betaGreenep1*mnuSelf+uOther[2]; 
 
  double uRelSq[3]; 
  // Zero out array with dot product of uSelf-uOther with itself. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    uRelSq[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    uRelSq[0] += 0.5*uSelf[a0+2]*uSelf[a0+2]-1.0*uOther[a0+2]*uSelf[a0+2]+0.5*uOther[a0+2]*uOther[a0+2]+0.5*uSelf[a0+1]*uSelf[a0+1]-1.0*uOther[a0+1]*uSelf[a0+1]+0.5*uOther[a0+1]*uOther[a0+1]+0.5*uSelf[a0]*uSelf[a0]-1.0*uOther[a0]*uSelf[a0]+0.5*uOther[a0]*uOther[a0]; 
    uRelSq[1] += uSelf[a0]*uSelf[a0+1]-1.0*uOther[a0]*uSelf[a0+1]-1.0*uSelf[a0]*uOther[a0+1]+uOther[a0]*uOther[a0+1]; 
    uRelSq[2] += uSelf[a0]*uSelf[a0+2]-1.0*uOther[a0]*uSelf[a0+2]-1.0*uSelf[a0]*uOther[a0+2]+uOther[a0]*uOther[a0+2]; 
  } 
 
  double relKinE[3]; 
  // Zero out array with ((beta+1)/2)*(mSelf+mOther)*(uSelf-uOther) . uRelDmnu. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += betaGreenep1*(0.25*uRelDmnu[a0+2]*uSelf[a0+2]*mSelf-0.25*uOther[a0+2]*uRelDmnu[a0+2]*mSelf+0.25*uRelDmnu[a0+1]*uSelf[a0+1]*mSelf-0.25*uOther[a0+1]*uRelDmnu[a0+1]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0]*mSelf-0.25*uOther[a0]*uRelDmnu[a0]*mSelf+0.25*uRelDmnu[a0+2]*uSelf[a0+2]*mOther-0.25*uOther[a0+2]*uRelDmnu[a0+2]*mOther+0.25*uRelDmnu[a0+1]*uSelf[a0+1]*mOther-0.25*uOther[a0+1]*uRelDmnu[a0+1]*mOther+0.25*uRelDmnu[a0]*uSelf[a0]*mOther-0.25*uOther[a0]*uRelDmnu[a0]*mOther); 
    relKinE[1] += betaGreenep1*(0.25*uRelDmnu[a0]*uSelf[a0+1]*mSelf+0.25*uSelf[a0]*uRelDmnu[a0+1]*mSelf-0.25*uOther[a0]*uRelDmnu[a0+1]*mSelf-0.25*uRelDmnu[a0]*uOther[a0+1]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0+1]*mOther+0.25*uSelf[a0]*uRelDmnu[a0+1]*mOther-0.25*uOther[a0]*uRelDmnu[a0+1]*mOther-0.25*uRelDmnu[a0]*uOther[a0+1]*mOther); 
    relKinE[2] += betaGreenep1*(0.25*uRelDmnu[a0]*uSelf[a0+2]*mSelf+0.25*uSelf[a0]*uRelDmnu[a0+2]*mSelf-0.25*uOther[a0]*uRelDmnu[a0+2]*mSelf-0.25*uRelDmnu[a0]*uOther[a0+2]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0+2]*mOther+0.25*uSelf[a0]*uRelDmnu[a0+2]*mOther-0.25*uOther[a0]*uRelDmnu[a0+2]*mOther-0.25*uRelDmnu[a0]*uOther[a0+2]*mOther); 
  } 
 
  double Tdiff[3]; 
  Tdiff[0] = vtSqSelf[0]*mSelf-1.0*vtSqOther[0]*mOther; 
  Tdiff[1] = vtSqSelf[1]*mSelf-1.0*vtSqOther[1]*mOther; 
  Tdiff[2] = vtSqSelf[2]*mSelf-1.0*vtSqOther[2]*mOther; 
 
  double diffSelf[3]; 
  diffSelf[0] = (0.5*m0rOther[2]*relKinE[2]+0.5*m0rOther[1]*relKinE[1]+0.5*m0rOther[0]*relKinE[0])*mnuOther-1.0*uRelSq[0]*mOther+3.0*Tdiff[0]; 
  diffSelf[1] = (0.5*m0rOther[0]*relKinE[1]+0.5*relKinE[0]*m0rOther[1])*mnuOther-1.0*uRelSq[1]*mOther+3.0*Tdiff[1]; 
  diffSelf[2] = (0.5*m0rOther[0]*relKinE[2]+0.5*relKinE[0]*m0rOther[2])*mnuOther-1.0*uRelSq[2]*mOther+3.0*Tdiff[2]; 
 
  // Divide diffSelf by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffSelf[0],diffSelf[1],diffSelf[2]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaSelf[3]; 
  Eigen::Map<VectorXd>(vtSqDeltaSelf,3,1) = dataDiv->u_S; 
 
  double diffOther[3]; 
  diffOther[0] = (0.5*m0rSelf[2]*relKinE[2]+0.5*m0rSelf[1]*relKinE[1]+0.5*m0rSelf[0]*relKinE[0])*mnuSelf-1.0*uRelSq[0]*mSelf-3.0*Tdiff[0]; 
  diffOther[1] = (0.5*m0rSelf[0]*relKinE[1]+0.5*relKinE[0]*m0rSelf[1])*mnuSelf-1.0*uRelSq[1]*mSelf-3.0*Tdiff[1]; 
  diffOther[2] = (0.5*m0rSelf[0]*relKinE[2]+0.5*relKinE[0]*m0rSelf[2])*mnuSelf-1.0*uRelSq[2]*mSelf-3.0*Tdiff[2]; 
 
  // Divide diffOther by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffOther[0],diffOther[1],diffOther[2]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaOther[3]; 
  Eigen::Map<VectorXd>(vtSqDeltaOther,3,1) = dataDiv->u_S; 
 
  // ... Cross-thermal speeds (squared) ... // 
  double deltaFacOther = (0.6666666666666666*mnuOther)/(mSelf+mOther); 
  vtSqCrossSelf[0] = (-0.5*m0rOther[2]*vtSqDeltaSelf[2]*deltaFacOther)-0.5*m0rOther[1]*vtSqDeltaSelf[1]*deltaFacOther-0.5*m0rOther[0]*vtSqDeltaSelf[0]*deltaFacOther+vtSqSelf[0]; 
  vtSqCrossSelf[1] = (-0.5*m0rOther[0]*vtSqDeltaSelf[1]*deltaFacOther)-0.5*vtSqDeltaSelf[0]*m0rOther[1]*deltaFacOther+vtSqSelf[1]; 
  vtSqCrossSelf[2] = (-0.5*m0rOther[0]*vtSqDeltaSelf[2]*deltaFacOther)-0.5*vtSqDeltaSelf[0]*m0rOther[2]*deltaFacOther+vtSqSelf[2]; 
 
  double deltaFacSelf = (0.6666666666666666*mnuSelf)/(mSelf+mOther); 
  vtSqCrossOther[0] = (-0.5*m0rSelf[2]*vtSqDeltaOther[2]*deltaFacSelf)-0.5*m0rSelf[1]*vtSqDeltaOther[1]*deltaFacSelf-0.5*m0rSelf[0]*vtSqDeltaOther[0]*deltaFacSelf+vtSqOther[0]; 
  vtSqCrossOther[1] = (-0.5*m0rSelf[0]*vtSqDeltaOther[1]*deltaFacSelf)-0.5*vtSqDeltaOther[0]*m0rSelf[1]*deltaFacSelf+vtSqOther[1]; 
  vtSqCrossOther[2] = (-0.5*m0rSelf[0]*vtSqDeltaOther[2]*deltaFacSelf)-0.5*vtSqDeltaOther[0]*m0rSelf[2]*deltaFacSelf+vtSqOther[2]; 
 
} 
 
void GkBGKCrossPrimMoments2x2vMax_P2(binOpData_t *dataDiv, const double betaGreenep1, const double mSelf, const double nuSelf, const double *m0Self, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *m0Other, const double *uOther, const double *vtSqOther, double *uCrossSelf, double *vtSqCrossSelf, double *uCrossOther, double *vtSqCrossOther) 
{ 
  // betaGreenep1:     free parameter beta+1. This has to be >0. 
  // m, nu:            mass and collisionality. 
  // m0:               zeroth moment of the distribution function. 
  // u,vtSq:           self primitive moments: mean flow velocity and thermal speed squared. 
  // uCross,vtSqCross: cross primitive moments: mean flow velocity and thermal speed squared. 
 
  // If a corner value is below zero, use cell average m0.
  bool cellAvg = false;
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]-0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]+0.8660254037844386*m0Self[2]-0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]-1.5*m0Self[3]-0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Self[5]+1.118033988749895*m0Self[4]+1.5*m0Self[3]+0.8660254037844386*m0Self[2]+0.8660254037844386*m0Self[1]+0.5*m0Self[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rSelf[6]; 
  if (cellAvg) { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = 0.0; 
    m0rSelf[2] = 0.0; 
    m0rSelf[3] = 0.0; 
    m0rSelf[4] = 0.0; 
    m0rSelf[5] = 0.0; 
  } else { 
    m0rSelf[0] = m0Self[0]; 
    m0rSelf[1] = m0Self[1]; 
    m0rSelf[2] = m0Self[2]; 
    m0rSelf[3] = m0Self[3]; 
    m0rSelf[4] = m0Self[4]; 
    m0rSelf[5] = m0Self[5]; 
  } 
 
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]-0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]+0.8660254037844386*m0Other[2]-0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]-1.5*m0Other[3]-0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
  if (1.118033988749895*m0Other[5]+1.118033988749895*m0Other[4]+1.5*m0Other[3]+0.8660254037844386*m0Other[2]+0.8660254037844386*m0Other[1]+0.5*m0Other[0] < 0) { 
    cellAvg = true;
  }
 
  double m0rOther[6]; 
  if (cellAvg) { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = 0.0; 
    m0rOther[2] = 0.0; 
    m0rOther[3] = 0.0; 
    m0rOther[4] = 0.0; 
    m0rOther[5] = 0.0; 
  } else { 
    m0rOther[0] = m0Other[0]; 
    m0rOther[1] = m0Other[1]; 
    m0rOther[2] = m0Other[2]; 
    m0rOther[3] = m0Other[3]; 
    m0rOther[4] = m0Other[4]; 
    m0rOther[5] = m0Other[5]; 
  } 
 
  double mnuSelf  = mSelf*nuSelf; 
  double mnuOther = mOther*nuOther; 
 
  double uRelDmnu[6]; 
 
  // ... Divide (uSelfX-uOtherX)/(mnuSelf*m0Self+mnuOther*m0Other) ... // 
  // Compute (uSelf-uOther). 
  uRelDmnu[0] = uSelf[0]-1.0*uOther[0]; 
  uRelDmnu[1] = uSelf[1]-1.0*uOther[1]; 
  uRelDmnu[2] = uSelf[2]-1.0*uOther[2]; 
  uRelDmnu[3] = uSelf[3]-1.0*uOther[3]; 
  uRelDmnu[4] = uSelf[4]-1.0*uOther[4]; 
  uRelDmnu[5] = uSelf[5]-1.0*uOther[5]; 
  // Fill AEM matrix. 
  dataDiv->AEM_S = Eigen::MatrixXd::Zero(6,6); 
  dataDiv->AEM_S(0,0) = 0.5*m0rSelf[0]*mnuSelf+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(0,1) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(0,2) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(0,3) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(0,4) = 0.5*m0rSelf[4]*mnuSelf+0.5*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(0,5) = 0.5*m0rSelf[5]*mnuSelf+0.5*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(1,0) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(1,1) = 0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(1,2) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(1,3) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(1,4) = 0.4472135954999579*m0rSelf[1]*mnuSelf+0.4472135954999579*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(2,0) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(2,1) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(2,2) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.4472135954999579*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(2,3) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(2,5) = 0.4472135954999579*m0rSelf[2]*mnuSelf+0.4472135954999579*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(3,0) = 0.5*m0rSelf[3]*mnuSelf+0.5*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,1) = 0.5*m0rSelf[2]*mnuSelf+0.5*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(3,2) = 0.5*m0rSelf[1]*mnuSelf+0.5*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(3,3) = 0.4472135954999579*m0rSelf[5]*mnuSelf+0.4472135954999579*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.4472135954999579*m0rOther[5]*mnuOther+0.4472135954999579*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(3,4) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(3,5) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(4,0) = 0.5*m0rSelf[4]*mnuSelf+0.5*m0rOther[4]*mnuOther; 
  dataDiv->AEM_S(4,1) = 0.4472135954999579*m0rSelf[1]*mnuSelf+0.4472135954999579*m0rOther[1]*mnuOther; 
  dataDiv->AEM_S(4,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(4,4) = 0.31943828249997*m0rSelf[4]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.31943828249997*m0rOther[4]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  dataDiv->AEM_S(5,0) = 0.5*m0rSelf[5]*mnuSelf+0.5*m0rOther[5]*mnuOther; 
  dataDiv->AEM_S(5,2) = 0.4472135954999579*m0rSelf[2]*mnuSelf+0.4472135954999579*m0rOther[2]*mnuOther; 
  dataDiv->AEM_S(5,3) = 0.4472135954999579*m0rSelf[3]*mnuSelf+0.4472135954999579*m0rOther[3]*mnuOther; 
  dataDiv->AEM_S(5,5) = 0.31943828249997*m0rSelf[5]*mnuSelf+0.5*m0rSelf[0]*mnuSelf+0.31943828249997*m0rOther[5]*mnuOther+0.5*m0rOther[0]*mnuOther; 
  // Fill BEV. 
  dataDiv->BEV_S << uRelDmnu[0],uRelDmnu[1],uRelDmnu[2],uRelDmnu[3],uRelDmnu[4],uRelDmnu[5]; 
  // Invert system of equations from weak division. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(uRelDmnu+0,6,1) = dataDiv->u_S; 
 
  // ... Component 1 of cross-velocity of this species ... // 
  uCrossSelf[0] = ((-0.5*m0rOther[5]*uRelDmnu[5])-0.5*m0rOther[4]*uRelDmnu[4]-0.5*m0rOther[3]*uRelDmnu[3]-0.5*m0rOther[2]*uRelDmnu[2]-0.5*m0rOther[1]*uRelDmnu[1]-0.5*m0rOther[0]*uRelDmnu[0])*betaGreenep1*mnuOther+uSelf[0]; 
  uCrossSelf[1] = ((-0.4472135954999579*m0rOther[1]*uRelDmnu[4])-0.4472135954999579*uRelDmnu[1]*m0rOther[4]-0.5*m0rOther[2]*uRelDmnu[3]-0.5*uRelDmnu[2]*m0rOther[3]-0.5*m0rOther[0]*uRelDmnu[1]-0.5*uRelDmnu[0]*m0rOther[1])*betaGreenep1*mnuOther+uSelf[1]; 
  uCrossSelf[2] = ((-0.4472135954999579*m0rOther[2]*uRelDmnu[5])-0.4472135954999579*uRelDmnu[2]*m0rOther[5]-0.5*m0rOther[1]*uRelDmnu[3]-0.5*uRelDmnu[1]*m0rOther[3]-0.5*m0rOther[0]*uRelDmnu[2]-0.5*uRelDmnu[0]*m0rOther[2])*betaGreenep1*mnuOther+uSelf[2]; 
  uCrossSelf[3] = ((-0.4472135954999579*m0rOther[3]*uRelDmnu[5])-0.4472135954999579*uRelDmnu[3]*m0rOther[5]-0.4472135954999579*m0rOther[3]*uRelDmnu[4]-0.4472135954999579*uRelDmnu[3]*m0rOther[4]-0.5*m0rOther[0]*uRelDmnu[3]-0.5*uRelDmnu[0]*m0rOther[3]-0.5*m0rOther[1]*uRelDmnu[2]-0.5*uRelDmnu[1]*m0rOther[2])*betaGreenep1*mnuOther+uSelf[3]; 
  uCrossSelf[4] = ((-0.31943828249997*m0rOther[4]*uRelDmnu[4])-0.5*m0rOther[0]*uRelDmnu[4]-0.5*uRelDmnu[0]*m0rOther[4]-0.4472135954999579*m0rOther[3]*uRelDmnu[3]-0.4472135954999579*m0rOther[1]*uRelDmnu[1])*betaGreenep1*mnuOther+uSelf[4]; 
  uCrossSelf[5] = ((-0.31943828249997*m0rOther[5]*uRelDmnu[5])-0.5*m0rOther[0]*uRelDmnu[5]-0.5*uRelDmnu[0]*m0rOther[5]-0.4472135954999579*m0rOther[3]*uRelDmnu[3]-0.4472135954999579*m0rOther[2]*uRelDmnu[2])*betaGreenep1*mnuOther+uSelf[5]; 
 
  // ... Component 1 of cross-velocity of the other species ... // 
  uCrossOther[0] = (0.5*m0rSelf[5]*uRelDmnu[5]+0.5*m0rSelf[4]*uRelDmnu[4]+0.5*m0rSelf[3]*uRelDmnu[3]+0.5*m0rSelf[2]*uRelDmnu[2]+0.5*m0rSelf[1]*uRelDmnu[1]+0.5*m0rSelf[0]*uRelDmnu[0])*betaGreenep1*mnuSelf+uOther[0]; 
  uCrossOther[1] = (0.4472135954999579*m0rSelf[1]*uRelDmnu[4]+0.4472135954999579*uRelDmnu[1]*m0rSelf[4]+0.5*m0rSelf[2]*uRelDmnu[3]+0.5*uRelDmnu[2]*m0rSelf[3]+0.5*m0rSelf[0]*uRelDmnu[1]+0.5*uRelDmnu[0]*m0rSelf[1])*betaGreenep1*mnuSelf+uOther[1]; 
  uCrossOther[2] = (0.4472135954999579*m0rSelf[2]*uRelDmnu[5]+0.4472135954999579*uRelDmnu[2]*m0rSelf[5]+0.5*m0rSelf[1]*uRelDmnu[3]+0.5*uRelDmnu[1]*m0rSelf[3]+0.5*m0rSelf[0]*uRelDmnu[2]+0.5*uRelDmnu[0]*m0rSelf[2])*betaGreenep1*mnuSelf+uOther[2]; 
  uCrossOther[3] = (0.4472135954999579*m0rSelf[3]*uRelDmnu[5]+0.4472135954999579*uRelDmnu[3]*m0rSelf[5]+0.4472135954999579*m0rSelf[3]*uRelDmnu[4]+0.4472135954999579*uRelDmnu[3]*m0rSelf[4]+0.5*m0rSelf[0]*uRelDmnu[3]+0.5*uRelDmnu[0]*m0rSelf[3]+0.5*m0rSelf[1]*uRelDmnu[2]+0.5*uRelDmnu[1]*m0rSelf[2])*betaGreenep1*mnuSelf+uOther[3]; 
  uCrossOther[4] = (0.31943828249997*m0rSelf[4]*uRelDmnu[4]+0.5*m0rSelf[0]*uRelDmnu[4]+0.5*uRelDmnu[0]*m0rSelf[4]+0.4472135954999579*m0rSelf[3]*uRelDmnu[3]+0.4472135954999579*m0rSelf[1]*uRelDmnu[1])*betaGreenep1*mnuSelf+uOther[4]; 
  uCrossOther[5] = (0.31943828249997*m0rSelf[5]*uRelDmnu[5]+0.5*m0rSelf[0]*uRelDmnu[5]+0.5*uRelDmnu[0]*m0rSelf[5]+0.4472135954999579*m0rSelf[3]*uRelDmnu[3]+0.4472135954999579*m0rSelf[2]*uRelDmnu[2])*betaGreenep1*mnuSelf+uOther[5]; 
 
  double uRelSq[6]; 
  // Zero out array with dot product of uSelf-uOther with itself. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    uRelSq[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    uRelSq[0] += 0.5*uSelf[a0+5]*uSelf[a0+5]-1.0*uOther[a0+5]*uSelf[a0+5]+0.5*uOther[a0+5]*uOther[a0+5]+0.5*uSelf[a0+4]*uSelf[a0+4]-1.0*uOther[a0+4]*uSelf[a0+4]+0.5*uOther[a0+4]*uOther[a0+4]+0.5*uSelf[a0+3]*uSelf[a0+3]-1.0*uOther[a0+3]*uSelf[a0+3]+0.5*uOther[a0+3]*uOther[a0+3]+0.5*uSelf[a0+2]*uSelf[a0+2]-1.0*uOther[a0+2]*uSelf[a0+2]+0.5*uOther[a0+2]*uOther[a0+2]+0.5*uSelf[a0+1]*uSelf[a0+1]-1.0*uOther[a0+1]*uSelf[a0+1]+0.5*uOther[a0+1]*uOther[a0+1]+0.5*uSelf[a0]*uSelf[a0]-1.0*uOther[a0]*uSelf[a0]+0.5*uOther[a0]*uOther[a0]; 
    uRelSq[1] += 0.8944271909999159*uSelf[a0+1]*uSelf[a0+4]-0.8944271909999159*uOther[a0+1]*uSelf[a0+4]-0.8944271909999159*uSelf[a0+1]*uOther[a0+4]+0.8944271909999159*uOther[a0+1]*uOther[a0+4]+uSelf[a0+2]*uSelf[a0+3]-1.0*uOther[a0+2]*uSelf[a0+3]-1.0*uSelf[a0+2]*uOther[a0+3]+uOther[a0+2]*uOther[a0+3]+uSelf[a0]*uSelf[a0+1]-1.0*uOther[a0]*uSelf[a0+1]-1.0*uSelf[a0]*uOther[a0+1]+uOther[a0]*uOther[a0+1]; 
    uRelSq[2] += 0.8944271909999159*uSelf[a0+2]*uSelf[a0+5]-0.8944271909999159*uOther[a0+2]*uSelf[a0+5]-0.8944271909999159*uSelf[a0+2]*uOther[a0+5]+0.8944271909999159*uOther[a0+2]*uOther[a0+5]+uSelf[a0+1]*uSelf[a0+3]-1.0*uOther[a0+1]*uSelf[a0+3]-1.0*uSelf[a0+1]*uOther[a0+3]+uOther[a0+1]*uOther[a0+3]+uSelf[a0]*uSelf[a0+2]-1.0*uOther[a0]*uSelf[a0+2]-1.0*uSelf[a0]*uOther[a0+2]+uOther[a0]*uOther[a0+2]; 
    uRelSq[3] += 0.8944271909999159*uSelf[a0+3]*uSelf[a0+5]-0.8944271909999159*uOther[a0+3]*uSelf[a0+5]-0.8944271909999159*uSelf[a0+3]*uOther[a0+5]+0.8944271909999159*uOther[a0+3]*uOther[a0+5]+0.8944271909999159*uSelf[a0+3]*uSelf[a0+4]-0.8944271909999159*uOther[a0+3]*uSelf[a0+4]-0.8944271909999159*uSelf[a0+3]*uOther[a0+4]+0.8944271909999159*uOther[a0+3]*uOther[a0+4]+uSelf[a0]*uSelf[a0+3]-1.0*uOther[a0]*uSelf[a0+3]-1.0*uSelf[a0]*uOther[a0+3]+uOther[a0]*uOther[a0+3]+uSelf[a0+1]*uSelf[a0+2]-1.0*uOther[a0+1]*uSelf[a0+2]-1.0*uSelf[a0+1]*uOther[a0+2]+uOther[a0+1]*uOther[a0+2]; 
    uRelSq[4] += 0.31943828249997*uSelf[a0+4]*uSelf[a0+4]-0.6388765649999399*uOther[a0+4]*uSelf[a0+4]+uSelf[a0]*uSelf[a0+4]-1.0*uOther[a0]*uSelf[a0+4]+0.31943828249997*uOther[a0+4]*uOther[a0+4]-1.0*uSelf[a0]*uOther[a0+4]+uOther[a0]*uOther[a0+4]+0.4472135954999579*uSelf[a0+3]*uSelf[a0+3]-0.8944271909999159*uOther[a0+3]*uSelf[a0+3]+0.4472135954999579*uOther[a0+3]*uOther[a0+3]+0.4472135954999579*uSelf[a0+1]*uSelf[a0+1]-0.8944271909999159*uOther[a0+1]*uSelf[a0+1]+0.4472135954999579*uOther[a0+1]*uOther[a0+1]; 
    uRelSq[5] += 0.31943828249997*uSelf[a0+5]*uSelf[a0+5]-0.6388765649999399*uOther[a0+5]*uSelf[a0+5]+uSelf[a0]*uSelf[a0+5]-1.0*uOther[a0]*uSelf[a0+5]+0.31943828249997*uOther[a0+5]*uOther[a0+5]-1.0*uSelf[a0]*uOther[a0+5]+uOther[a0]*uOther[a0+5]+0.4472135954999579*uSelf[a0+3]*uSelf[a0+3]-0.8944271909999159*uOther[a0+3]*uSelf[a0+3]+0.4472135954999579*uOther[a0+3]*uOther[a0+3]+0.4472135954999579*uSelf[a0+2]*uSelf[a0+2]-0.8944271909999159*uOther[a0+2]*uSelf[a0+2]+0.4472135954999579*uOther[a0+2]*uOther[a0+2]; 
  } 
 
  double relKinE[6]; 
  // Zero out array with ((beta+1)/2)*(mSelf+mOther)*(uSelf-uOther) . uRelDmnu. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    relKinE[vd] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    relKinE[0] += betaGreenep1*(0.25*uRelDmnu[a0+5]*uSelf[a0+5]*mSelf-0.25*uOther[a0+5]*uRelDmnu[a0+5]*mSelf+0.25*uRelDmnu[a0+4]*uSelf[a0+4]*mSelf-0.25*uOther[a0+4]*uRelDmnu[a0+4]*mSelf+0.25*uRelDmnu[a0+3]*uSelf[a0+3]*mSelf-0.25*uOther[a0+3]*uRelDmnu[a0+3]*mSelf+0.25*uRelDmnu[a0+2]*uSelf[a0+2]*mSelf-0.25*uOther[a0+2]*uRelDmnu[a0+2]*mSelf+0.25*uRelDmnu[a0+1]*uSelf[a0+1]*mSelf-0.25*uOther[a0+1]*uRelDmnu[a0+1]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0]*mSelf-0.25*uOther[a0]*uRelDmnu[a0]*mSelf+0.25*uRelDmnu[a0+5]*uSelf[a0+5]*mOther-0.25*uOther[a0+5]*uRelDmnu[a0+5]*mOther+0.25*uRelDmnu[a0+4]*uSelf[a0+4]*mOther-0.25*uOther[a0+4]*uRelDmnu[a0+4]*mOther+0.25*uRelDmnu[a0+3]*uSelf[a0+3]*mOther-0.25*uOther[a0+3]*uRelDmnu[a0+3]*mOther+0.25*uRelDmnu[a0+2]*uSelf[a0+2]*mOther-0.25*uOther[a0+2]*uRelDmnu[a0+2]*mOther+0.25*uRelDmnu[a0+1]*uSelf[a0+1]*mOther-0.25*uOther[a0+1]*uRelDmnu[a0+1]*mOther+0.25*uRelDmnu[a0]*uSelf[a0]*mOther-0.25*uOther[a0]*uRelDmnu[a0]*mOther); 
    relKinE[1] += betaGreenep1*(0.223606797749979*uRelDmnu[a0+1]*uSelf[a0+4]*mSelf+0.223606797749979*uSelf[a0+1]*uRelDmnu[a0+4]*mSelf-0.223606797749979*uOther[a0+1]*uRelDmnu[a0+4]*mSelf-0.223606797749979*uRelDmnu[a0+1]*uOther[a0+4]*mSelf+0.25*uRelDmnu[a0+2]*uSelf[a0+3]*mSelf+0.25*uSelf[a0+2]*uRelDmnu[a0+3]*mSelf-0.25*uOther[a0+2]*uRelDmnu[a0+3]*mSelf-0.25*uRelDmnu[a0+2]*uOther[a0+3]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0+1]*mSelf+0.25*uSelf[a0]*uRelDmnu[a0+1]*mSelf-0.25*uOther[a0]*uRelDmnu[a0+1]*mSelf-0.25*uRelDmnu[a0]*uOther[a0+1]*mSelf+0.223606797749979*uRelDmnu[a0+1]*uSelf[a0+4]*mOther+0.223606797749979*uSelf[a0+1]*uRelDmnu[a0+4]*mOther-0.223606797749979*uOther[a0+1]*uRelDmnu[a0+4]*mOther-0.223606797749979*uRelDmnu[a0+1]*uOther[a0+4]*mOther+0.25*uRelDmnu[a0+2]*uSelf[a0+3]*mOther+0.25*uSelf[a0+2]*uRelDmnu[a0+3]*mOther-0.25*uOther[a0+2]*uRelDmnu[a0+3]*mOther-0.25*uRelDmnu[a0+2]*uOther[a0+3]*mOther+0.25*uRelDmnu[a0]*uSelf[a0+1]*mOther+0.25*uSelf[a0]*uRelDmnu[a0+1]*mOther-0.25*uOther[a0]*uRelDmnu[a0+1]*mOther-0.25*uRelDmnu[a0]*uOther[a0+1]*mOther); 
    relKinE[2] += betaGreenep1*(0.223606797749979*uRelDmnu[a0+2]*uSelf[a0+5]*mSelf+0.223606797749979*uSelf[a0+2]*uRelDmnu[a0+5]*mSelf-0.223606797749979*uOther[a0+2]*uRelDmnu[a0+5]*mSelf-0.223606797749979*uRelDmnu[a0+2]*uOther[a0+5]*mSelf+0.25*uRelDmnu[a0+1]*uSelf[a0+3]*mSelf+0.25*uSelf[a0+1]*uRelDmnu[a0+3]*mSelf-0.25*uOther[a0+1]*uRelDmnu[a0+3]*mSelf-0.25*uRelDmnu[a0+1]*uOther[a0+3]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0+2]*mSelf+0.25*uSelf[a0]*uRelDmnu[a0+2]*mSelf-0.25*uOther[a0]*uRelDmnu[a0+2]*mSelf-0.25*uRelDmnu[a0]*uOther[a0+2]*mSelf+0.223606797749979*uRelDmnu[a0+2]*uSelf[a0+5]*mOther+0.223606797749979*uSelf[a0+2]*uRelDmnu[a0+5]*mOther-0.223606797749979*uOther[a0+2]*uRelDmnu[a0+5]*mOther-0.223606797749979*uRelDmnu[a0+2]*uOther[a0+5]*mOther+0.25*uRelDmnu[a0+1]*uSelf[a0+3]*mOther+0.25*uSelf[a0+1]*uRelDmnu[a0+3]*mOther-0.25*uOther[a0+1]*uRelDmnu[a0+3]*mOther-0.25*uRelDmnu[a0+1]*uOther[a0+3]*mOther+0.25*uRelDmnu[a0]*uSelf[a0+2]*mOther+0.25*uSelf[a0]*uRelDmnu[a0+2]*mOther-0.25*uOther[a0]*uRelDmnu[a0+2]*mOther-0.25*uRelDmnu[a0]*uOther[a0+2]*mOther); 
    relKinE[3] += betaGreenep1*(0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+5]*mSelf+0.223606797749979*uSelf[a0+3]*uRelDmnu[a0+5]*mSelf-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+5]*mSelf-0.223606797749979*uRelDmnu[a0+3]*uOther[a0+5]*mSelf+0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+4]*mSelf+0.223606797749979*uSelf[a0+3]*uRelDmnu[a0+4]*mSelf-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+4]*mSelf-0.223606797749979*uRelDmnu[a0+3]*uOther[a0+4]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0+3]*mSelf+0.25*uSelf[a0]*uRelDmnu[a0+3]*mSelf-0.25*uOther[a0]*uRelDmnu[a0+3]*mSelf-0.25*uRelDmnu[a0]*uOther[a0+3]*mSelf+0.25*uRelDmnu[a0+1]*uSelf[a0+2]*mSelf+0.25*uSelf[a0+1]*uRelDmnu[a0+2]*mSelf-0.25*uOther[a0+1]*uRelDmnu[a0+2]*mSelf-0.25*uRelDmnu[a0+1]*uOther[a0+2]*mSelf+0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+5]*mOther+0.223606797749979*uSelf[a0+3]*uRelDmnu[a0+5]*mOther-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+5]*mOther-0.223606797749979*uRelDmnu[a0+3]*uOther[a0+5]*mOther+0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+4]*mOther+0.223606797749979*uSelf[a0+3]*uRelDmnu[a0+4]*mOther-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+4]*mOther-0.223606797749979*uRelDmnu[a0+3]*uOther[a0+4]*mOther+0.25*uRelDmnu[a0]*uSelf[a0+3]*mOther+0.25*uSelf[a0]*uRelDmnu[a0+3]*mOther-0.25*uOther[a0]*uRelDmnu[a0+3]*mOther-0.25*uRelDmnu[a0]*uOther[a0+3]*mOther+0.25*uRelDmnu[a0+1]*uSelf[a0+2]*mOther+0.25*uSelf[a0+1]*uRelDmnu[a0+2]*mOther-0.25*uOther[a0+1]*uRelDmnu[a0+2]*mOther-0.25*uRelDmnu[a0+1]*uOther[a0+2]*mOther); 
    relKinE[4] += betaGreenep1*(0.159719141249985*uRelDmnu[a0+4]*uSelf[a0+4]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0+4]*mSelf-0.159719141249985*uOther[a0+4]*uRelDmnu[a0+4]*mSelf+0.25*uSelf[a0]*uRelDmnu[a0+4]*mSelf-0.25*uOther[a0]*uRelDmnu[a0+4]*mSelf-0.25*uRelDmnu[a0]*uOther[a0+4]*mSelf+0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+3]*mSelf-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+3]*mSelf+0.223606797749979*uRelDmnu[a0+1]*uSelf[a0+1]*mSelf-0.223606797749979*uOther[a0+1]*uRelDmnu[a0+1]*mSelf+0.159719141249985*uRelDmnu[a0+4]*uSelf[a0+4]*mOther+0.25*uRelDmnu[a0]*uSelf[a0+4]*mOther-0.159719141249985*uOther[a0+4]*uRelDmnu[a0+4]*mOther+0.25*uSelf[a0]*uRelDmnu[a0+4]*mOther-0.25*uOther[a0]*uRelDmnu[a0+4]*mOther-0.25*uRelDmnu[a0]*uOther[a0+4]*mOther+0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+3]*mOther-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+3]*mOther+0.223606797749979*uRelDmnu[a0+1]*uSelf[a0+1]*mOther-0.223606797749979*uOther[a0+1]*uRelDmnu[a0+1]*mOther); 
    relKinE[5] += betaGreenep1*(0.159719141249985*uRelDmnu[a0+5]*uSelf[a0+5]*mSelf+0.25*uRelDmnu[a0]*uSelf[a0+5]*mSelf-0.159719141249985*uOther[a0+5]*uRelDmnu[a0+5]*mSelf+0.25*uSelf[a0]*uRelDmnu[a0+5]*mSelf-0.25*uOther[a0]*uRelDmnu[a0+5]*mSelf-0.25*uRelDmnu[a0]*uOther[a0+5]*mSelf+0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+3]*mSelf-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+3]*mSelf+0.223606797749979*uRelDmnu[a0+2]*uSelf[a0+2]*mSelf-0.223606797749979*uOther[a0+2]*uRelDmnu[a0+2]*mSelf+0.159719141249985*uRelDmnu[a0+5]*uSelf[a0+5]*mOther+0.25*uRelDmnu[a0]*uSelf[a0+5]*mOther-0.159719141249985*uOther[a0+5]*uRelDmnu[a0+5]*mOther+0.25*uSelf[a0]*uRelDmnu[a0+5]*mOther-0.25*uOther[a0]*uRelDmnu[a0+5]*mOther-0.25*uRelDmnu[a0]*uOther[a0+5]*mOther+0.223606797749979*uRelDmnu[a0+3]*uSelf[a0+3]*mOther-0.223606797749979*uOther[a0+3]*uRelDmnu[a0+3]*mOther+0.223606797749979*uRelDmnu[a0+2]*uSelf[a0+2]*mOther-0.223606797749979*uOther[a0+2]*uRelDmnu[a0+2]*mOther); 
  } 
 
  double Tdiff[6]; 
  Tdiff[0] = vtSqSelf[0]*mSelf-1.0*vtSqOther[0]*mOther; 
  Tdiff[1] = vtSqSelf[1]*mSelf-1.0*vtSqOther[1]*mOther; 
  Tdiff[2] = vtSqSelf[2]*mSelf-1.0*vtSqOther[2]*mOther; 
  Tdiff[3] = vtSqSelf[3]*mSelf-1.0*vtSqOther[3]*mOther; 
  Tdiff[4] = vtSqSelf[4]*mSelf-1.0*vtSqOther[4]*mOther; 
  Tdiff[5] = vtSqSelf[5]*mSelf-1.0*vtSqOther[5]*mOther; 
 
  double diffSelf[6]; 
  diffSelf[0] = (0.5*m0rOther[5]*relKinE[5]+0.5*m0rOther[4]*relKinE[4]+0.5*m0rOther[3]*relKinE[3]+0.5*m0rOther[2]*relKinE[2]+0.5*m0rOther[1]*relKinE[1]+0.5*m0rOther[0]*relKinE[0])*mnuOther-1.0*uRelSq[0]*mOther+3.0*Tdiff[0]; 
  diffSelf[1] = (0.4472135954999579*m0rOther[1]*relKinE[4]+0.4472135954999579*relKinE[1]*m0rOther[4]+0.5*m0rOther[2]*relKinE[3]+0.5*relKinE[2]*m0rOther[3]+0.5*m0rOther[0]*relKinE[1]+0.5*relKinE[0]*m0rOther[1])*mnuOther-1.0*uRelSq[1]*mOther+3.0*Tdiff[1]; 
  diffSelf[2] = (0.4472135954999579*m0rOther[2]*relKinE[5]+0.4472135954999579*relKinE[2]*m0rOther[5]+0.5*m0rOther[1]*relKinE[3]+0.5*relKinE[1]*m0rOther[3]+0.5*m0rOther[0]*relKinE[2]+0.5*relKinE[0]*m0rOther[2])*mnuOther-1.0*uRelSq[2]*mOther+3.0*Tdiff[2]; 
  diffSelf[3] = (0.4472135954999579*m0rOther[3]*relKinE[5]+0.4472135954999579*relKinE[3]*m0rOther[5]+0.4472135954999579*m0rOther[3]*relKinE[4]+0.4472135954999579*relKinE[3]*m0rOther[4]+0.5*m0rOther[0]*relKinE[3]+0.5*relKinE[0]*m0rOther[3]+0.5*m0rOther[1]*relKinE[2]+0.5*relKinE[1]*m0rOther[2])*mnuOther-1.0*uRelSq[3]*mOther+3.0*Tdiff[3]; 
  diffSelf[4] = (0.31943828249997*m0rOther[4]*relKinE[4]+0.5*m0rOther[0]*relKinE[4]+0.5*relKinE[0]*m0rOther[4]+0.4472135954999579*m0rOther[3]*relKinE[3]+0.4472135954999579*m0rOther[1]*relKinE[1])*mnuOther-1.0*uRelSq[4]*mOther+3.0*Tdiff[4]; 
  diffSelf[5] = (0.31943828249997*m0rOther[5]*relKinE[5]+0.5*m0rOther[0]*relKinE[5]+0.5*relKinE[0]*m0rOther[5]+0.4472135954999579*m0rOther[3]*relKinE[3]+0.4472135954999579*m0rOther[2]*relKinE[2])*mnuOther-1.0*uRelSq[5]*mOther+3.0*Tdiff[5]; 
 
  // Divide diffSelf by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffSelf[0],diffSelf[1],diffSelf[2],diffSelf[3],diffSelf[4],diffSelf[5]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaSelf[6]; 
  Eigen::Map<VectorXd>(vtSqDeltaSelf,6,1) = dataDiv->u_S; 
 
  double diffOther[6]; 
  diffOther[0] = (0.5*m0rSelf[5]*relKinE[5]+0.5*m0rSelf[4]*relKinE[4]+0.5*m0rSelf[3]*relKinE[3]+0.5*m0rSelf[2]*relKinE[2]+0.5*m0rSelf[1]*relKinE[1]+0.5*m0rSelf[0]*relKinE[0])*mnuSelf-1.0*uRelSq[0]*mSelf-3.0*Tdiff[0]; 
  diffOther[1] = (0.4472135954999579*m0rSelf[1]*relKinE[4]+0.4472135954999579*relKinE[1]*m0rSelf[4]+0.5*m0rSelf[2]*relKinE[3]+0.5*relKinE[2]*m0rSelf[3]+0.5*m0rSelf[0]*relKinE[1]+0.5*relKinE[0]*m0rSelf[1])*mnuSelf-1.0*uRelSq[1]*mSelf-3.0*Tdiff[1]; 
  diffOther[2] = (0.4472135954999579*m0rSelf[2]*relKinE[5]+0.4472135954999579*relKinE[2]*m0rSelf[5]+0.5*m0rSelf[1]*relKinE[3]+0.5*relKinE[1]*m0rSelf[3]+0.5*m0rSelf[0]*relKinE[2]+0.5*relKinE[0]*m0rSelf[2])*mnuSelf-1.0*uRelSq[2]*mSelf-3.0*Tdiff[2]; 
  diffOther[3] = (0.4472135954999579*m0rSelf[3]*relKinE[5]+0.4472135954999579*relKinE[3]*m0rSelf[5]+0.4472135954999579*m0rSelf[3]*relKinE[4]+0.4472135954999579*relKinE[3]*m0rSelf[4]+0.5*m0rSelf[0]*relKinE[3]+0.5*relKinE[0]*m0rSelf[3]+0.5*m0rSelf[1]*relKinE[2]+0.5*relKinE[1]*m0rSelf[2])*mnuSelf-1.0*uRelSq[3]*mSelf-3.0*Tdiff[3]; 
  diffOther[4] = (0.31943828249997*m0rSelf[4]*relKinE[4]+0.5*m0rSelf[0]*relKinE[4]+0.5*relKinE[0]*m0rSelf[4]+0.4472135954999579*m0rSelf[3]*relKinE[3]+0.4472135954999579*m0rSelf[1]*relKinE[1])*mnuSelf-1.0*uRelSq[4]*mSelf-3.0*Tdiff[4]; 
  diffOther[5] = (0.31943828249997*m0rSelf[5]*relKinE[5]+0.5*m0rSelf[0]*relKinE[5]+0.5*relKinE[0]*m0rSelf[5]+0.4472135954999579*m0rSelf[3]*relKinE[3]+0.4472135954999579*m0rSelf[2]*relKinE[2])*mnuSelf-1.0*uRelSq[5]*mSelf-3.0*Tdiff[5]; 
 
  // Divide diffOther by mnuSelf*m0Self+mnuOther*m0Other. 
  dataDiv->BEV_S << diffOther[0],diffOther[1],diffOther[2],diffOther[3],diffOther[4],diffOther[5]; 
  // Invert system of equations from weak division. dataDiv.AEM was filled earlier. 
  dataDiv->u_S = dataDiv->AEM_S.colPivHouseholderQr().solve(dataDiv->BEV_S); 
  // Copy data from Eigen vector. 
  double vtSqDeltaOther[6]; 
  Eigen::Map<VectorXd>(vtSqDeltaOther,6,1) = dataDiv->u_S; 
 
  // ... Cross-thermal speeds (squared) ... // 
  double deltaFacOther = (0.6666666666666666*mnuOther)/(mSelf+mOther); 
  vtSqCrossSelf[0] = (-0.5*m0rOther[5]*vtSqDeltaSelf[5]*deltaFacOther)-0.5*m0rOther[4]*vtSqDeltaSelf[4]*deltaFacOther-0.5*m0rOther[3]*vtSqDeltaSelf[3]*deltaFacOther-0.5*m0rOther[2]*vtSqDeltaSelf[2]*deltaFacOther-0.5*m0rOther[1]*vtSqDeltaSelf[1]*deltaFacOther-0.5*m0rOther[0]*vtSqDeltaSelf[0]*deltaFacOther+vtSqSelf[0]; 
  vtSqCrossSelf[1] = (-0.4472135954999579*m0rOther[1]*vtSqDeltaSelf[4]*deltaFacOther)-0.4472135954999579*vtSqDeltaSelf[1]*m0rOther[4]*deltaFacOther-0.5*m0rOther[2]*vtSqDeltaSelf[3]*deltaFacOther-0.5*vtSqDeltaSelf[2]*m0rOther[3]*deltaFacOther-0.5*m0rOther[0]*vtSqDeltaSelf[1]*deltaFacOther-0.5*vtSqDeltaSelf[0]*m0rOther[1]*deltaFacOther+vtSqSelf[1]; 
  vtSqCrossSelf[2] = (-0.4472135954999579*m0rOther[2]*vtSqDeltaSelf[5]*deltaFacOther)-0.4472135954999579*vtSqDeltaSelf[2]*m0rOther[5]*deltaFacOther-0.5*m0rOther[1]*vtSqDeltaSelf[3]*deltaFacOther-0.5*vtSqDeltaSelf[1]*m0rOther[3]*deltaFacOther-0.5*m0rOther[0]*vtSqDeltaSelf[2]*deltaFacOther-0.5*vtSqDeltaSelf[0]*m0rOther[2]*deltaFacOther+vtSqSelf[2]; 
  vtSqCrossSelf[3] = (-0.4472135954999579*m0rOther[3]*vtSqDeltaSelf[5]*deltaFacOther)-0.4472135954999579*vtSqDeltaSelf[3]*m0rOther[5]*deltaFacOther-0.4472135954999579*m0rOther[3]*vtSqDeltaSelf[4]*deltaFacOther-0.4472135954999579*vtSqDeltaSelf[3]*m0rOther[4]*deltaFacOther-0.5*m0rOther[0]*vtSqDeltaSelf[3]*deltaFacOther-0.5*vtSqDeltaSelf[0]*m0rOther[3]*deltaFacOther-0.5*m0rOther[1]*vtSqDeltaSelf[2]*deltaFacOther-0.5*vtSqDeltaSelf[1]*m0rOther[2]*deltaFacOther+vtSqSelf[3]; 
  vtSqCrossSelf[4] = (-0.31943828249997*m0rOther[4]*vtSqDeltaSelf[4]*deltaFacOther)-0.5*m0rOther[0]*vtSqDeltaSelf[4]*deltaFacOther-0.5*vtSqDeltaSelf[0]*m0rOther[4]*deltaFacOther-0.4472135954999579*m0rOther[3]*vtSqDeltaSelf[3]*deltaFacOther-0.4472135954999579*m0rOther[1]*vtSqDeltaSelf[1]*deltaFacOther+vtSqSelf[4]; 
  vtSqCrossSelf[5] = (-0.31943828249997*m0rOther[5]*vtSqDeltaSelf[5]*deltaFacOther)-0.5*m0rOther[0]*vtSqDeltaSelf[5]*deltaFacOther-0.5*vtSqDeltaSelf[0]*m0rOther[5]*deltaFacOther-0.4472135954999579*m0rOther[3]*vtSqDeltaSelf[3]*deltaFacOther-0.4472135954999579*m0rOther[2]*vtSqDeltaSelf[2]*deltaFacOther+vtSqSelf[5]; 
 
  double deltaFacSelf = (0.6666666666666666*mnuSelf)/(mSelf+mOther); 
  vtSqCrossOther[0] = (-0.5*m0rSelf[5]*vtSqDeltaOther[5]*deltaFacSelf)-0.5*m0rSelf[4]*vtSqDeltaOther[4]*deltaFacSelf-0.5*m0rSelf[3]*vtSqDeltaOther[3]*deltaFacSelf-0.5*m0rSelf[2]*vtSqDeltaOther[2]*deltaFacSelf-0.5*m0rSelf[1]*vtSqDeltaOther[1]*deltaFacSelf-0.5*m0rSelf[0]*vtSqDeltaOther[0]*deltaFacSelf+vtSqOther[0]; 
  vtSqCrossOther[1] = (-0.4472135954999579*m0rSelf[1]*vtSqDeltaOther[4]*deltaFacSelf)-0.4472135954999579*vtSqDeltaOther[1]*m0rSelf[4]*deltaFacSelf-0.5*m0rSelf[2]*vtSqDeltaOther[3]*deltaFacSelf-0.5*vtSqDeltaOther[2]*m0rSelf[3]*deltaFacSelf-0.5*m0rSelf[0]*vtSqDeltaOther[1]*deltaFacSelf-0.5*vtSqDeltaOther[0]*m0rSelf[1]*deltaFacSelf+vtSqOther[1]; 
  vtSqCrossOther[2] = (-0.4472135954999579*m0rSelf[2]*vtSqDeltaOther[5]*deltaFacSelf)-0.4472135954999579*vtSqDeltaOther[2]*m0rSelf[5]*deltaFacSelf-0.5*m0rSelf[1]*vtSqDeltaOther[3]*deltaFacSelf-0.5*vtSqDeltaOther[1]*m0rSelf[3]*deltaFacSelf-0.5*m0rSelf[0]*vtSqDeltaOther[2]*deltaFacSelf-0.5*vtSqDeltaOther[0]*m0rSelf[2]*deltaFacSelf+vtSqOther[2]; 
  vtSqCrossOther[3] = (-0.4472135954999579*m0rSelf[3]*vtSqDeltaOther[5]*deltaFacSelf)-0.4472135954999579*vtSqDeltaOther[3]*m0rSelf[5]*deltaFacSelf-0.4472135954999579*m0rSelf[3]*vtSqDeltaOther[4]*deltaFacSelf-0.4472135954999579*vtSqDeltaOther[3]*m0rSelf[4]*deltaFacSelf-0.5*m0rSelf[0]*vtSqDeltaOther[3]*deltaFacSelf-0.5*vtSqDeltaOther[0]*m0rSelf[3]*deltaFacSelf-0.5*m0rSelf[1]*vtSqDeltaOther[2]*deltaFacSelf-0.5*vtSqDeltaOther[1]*m0rSelf[2]*deltaFacSelf+vtSqOther[3]; 
  vtSqCrossOther[4] = (-0.31943828249997*m0rSelf[4]*vtSqDeltaOther[4]*deltaFacSelf)-0.5*m0rSelf[0]*vtSqDeltaOther[4]*deltaFacSelf-0.5*vtSqDeltaOther[0]*m0rSelf[4]*deltaFacSelf-0.4472135954999579*m0rSelf[3]*vtSqDeltaOther[3]*deltaFacSelf-0.4472135954999579*m0rSelf[1]*vtSqDeltaOther[1]*deltaFacSelf+vtSqOther[4]; 
  vtSqCrossOther[5] = (-0.31943828249997*m0rSelf[5]*vtSqDeltaOther[5]*deltaFacSelf)-0.5*m0rSelf[0]*vtSqDeltaOther[5]*deltaFacSelf-0.5*vtSqDeltaOther[0]*m0rSelf[5]*deltaFacSelf-0.4472135954999579*m0rSelf[3]*vtSqDeltaOther[3]*deltaFacSelf-0.4472135954999579*m0rSelf[2]*vtSqDeltaOther[2]*deltaFacSelf+vtSqOther[5]; 
 
} 
 
