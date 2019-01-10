#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmCrossPrimMomentsGreene1x1vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[2]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[2]; 
  for (unsigned short int k=0; k<2; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uSMuOSq[0] += 0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+1]+0.7071067811865475*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.414213562373095*uSMuO[a0]*uSMuO[a0+1]; 
  } 
 
  // ..... Get the cross flow velocity uSelf2 ..... // 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uCross[a0] = 0.5*(uSelf[a0]+uOther[a0])-0.5*uSMuO[a0]*beta; 
    uCross[a0+1] = 0.5*(uSelf[a0+1]+uOther[a0+1])-0.5*uSMuO[a0+1]*beta; 
 
  } 
 
  double mBetaFrac = (0.5*(beta+1.0))/(mRat+1.0); 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = (vtSqOther[0]+uSMuOSq[0])*mBetaFrac*mRat-1.0*vtSqSelf[0]*mBetaFrac+vtSqSelf[0]; 
  vtSqCross[1] = (vtSqOther[1]+uSMuOSq[1])*mBetaFrac*mRat-1.0*vtSqSelf[1]*mBetaFrac+vtSqSelf[1]; 
 
} 
 
void VmCrossPrimMomentsGreene1x1vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[3]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[3]; 
  for (unsigned short int k=0; k<3; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uSMuOSq[0] += 0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+2]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+1]+0.7071067811865475*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.264911064067352*uSMuO[a0+1]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.4517539514526256*uSMuO[a0+2]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+2]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+1]; 
  } 
 
  // ..... Get the cross flow velocity uSelf2 ..... // 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uCross[a0] = 0.5*(uSelf[a0]+uOther[a0])-0.5*uSMuO[a0]*beta; 
    uCross[a0+1] = 0.5*(uSelf[a0+1]+uOther[a0+1])-0.5*uSMuO[a0+1]*beta; 
    uCross[a0+2] = 0.5*(uSelf[a0+2]+uOther[a0+2])-0.5*uSMuO[a0+2]*beta; 
 
  } 
 
  double mBetaFrac = (0.5*(beta+1.0))/(mRat+1.0); 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = (vtSqOther[0]+uSMuOSq[0])*mBetaFrac*mRat-1.0*vtSqSelf[0]*mBetaFrac+vtSqSelf[0]; 
  vtSqCross[1] = (vtSqOther[1]+uSMuOSq[1])*mBetaFrac*mRat-1.0*vtSqSelf[1]*mBetaFrac+vtSqSelf[1]; 
  vtSqCross[2] = (vtSqOther[2]+uSMuOSq[2])*mBetaFrac*mRat-1.0*vtSqSelf[2]*mBetaFrac+vtSqSelf[2]; 
 
} 
 
void VmCrossPrimMomentsGreene1x1vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[4]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uSMuOSq[0] += 0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+3]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+2]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+1]+0.7071067811865475*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.242118006816237*uSMuO[a0+2]*uSMuO[a0+3]+1.264911064067352*uSMuO[a0+1]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.421637021355784*uSMuO[a0+3]*uSMuO[a0+3]+1.242118006816237*uSMuO[a0+1]*uSMuO[a0+3]+0.4517539514526256*uSMuO[a0+2]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+2]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+1]; 
    uSMuOSq[3] += 0.8432740427115681*uSMuO[a0+2]*uSMuO[a0+3]+1.414213562373095*uSMuO[a0]*uSMuO[a0+3]+1.242118006816237*uSMuO[a0+1]*uSMuO[a0+2]; 
  } 
 
  // ..... Get the cross flow velocity uSelf2 ..... // 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uCross[a0] = 0.5*(uSelf[a0]+uOther[a0])-0.5*uSMuO[a0]*beta; 
    uCross[a0+1] = 0.5*(uSelf[a0+1]+uOther[a0+1])-0.5*uSMuO[a0+1]*beta; 
    uCross[a0+2] = 0.5*(uSelf[a0+2]+uOther[a0+2])-0.5*uSMuO[a0+2]*beta; 
    uCross[a0+3] = 0.5*(uSelf[a0+3]+uOther[a0+3])-0.5*uSMuO[a0+3]*beta; 
 
  } 
 
  double mBetaFrac = (0.5*(beta+1.0))/(mRat+1.0); 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = (vtSqOther[0]+uSMuOSq[0])*mBetaFrac*mRat-1.0*vtSqSelf[0]*mBetaFrac+vtSqSelf[0]; 
  vtSqCross[1] = (vtSqOther[1]+uSMuOSq[1])*mBetaFrac*mRat-1.0*vtSqSelf[1]*mBetaFrac+vtSqSelf[1]; 
  vtSqCross[2] = (vtSqOther[2]+uSMuOSq[2])*mBetaFrac*mRat-1.0*vtSqSelf[2]*mBetaFrac+vtSqSelf[2]; 
  vtSqCross[3] = (vtSqOther[3]+uSMuOSq[3])*mBetaFrac*mRat-1.0*vtSqSelf[3]*mBetaFrac+vtSqSelf[3]; 
 
} 
 
void VmCrossPrimMomentsHeavyIon1x1vMax_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[2]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[2]; 
  for (unsigned short int k=0; k<2; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uSMuOSq[0] += 0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+1]+0.7071067811865475*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.414213562373095*uSMuO[a0]*uSMuO[a0+1]; 
  } 
 
  // ..... Get the cross flow velocity uCross ..... // 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uCross[a0] = uOther[a0]; 
    uCross[a0+1] = uOther[a0+1]; 
 
  } 
 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = uSMuOSq[0]*mRat+vtSqSelf[0]; 
  vtSqCross[1] = uSMuOSq[1]*mRat+vtSqSelf[1]; 
 
} 
 
void VmCrossPrimMomentsHeavyIon1x1vMax_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[3]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[3]; 
  for (unsigned short int k=0; k<3; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uSMuOSq[0] += 0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+2]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+1]+0.7071067811865475*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.264911064067352*uSMuO[a0+1]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.4517539514526256*uSMuO[a0+2]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+2]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+1]; 
  } 
 
  // ..... Get the cross flow velocity uCross ..... // 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uCross[a0] = uOther[a0]; 
    uCross[a0+1] = uOther[a0+1]; 
    uCross[a0+2] = uOther[a0+2]; 
 
  } 
 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = uSMuOSq[0]*mRat+vtSqSelf[0]; 
  vtSqCross[1] = uSMuOSq[1]*mRat+vtSqSelf[1]; 
  vtSqCross[2] = uSMuOSq[2]*mRat+vtSqSelf[2]; 
 
} 
 
void VmCrossPrimMomentsHeavyIon1x1vMax_P3(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[4]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uSMuOSq[0] += 0.7071067811865475*uSMuO[a0+3]*uSMuO[a0+3]+0.7071067811865475*uSMuO[a0+2]*uSMuO[a0+2]+0.7071067811865475*uSMuO[a0+1]*uSMuO[a0+1]+0.7071067811865475*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.242118006816237*uSMuO[a0+2]*uSMuO[a0+3]+1.264911064067352*uSMuO[a0+1]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.421637021355784*uSMuO[a0+3]*uSMuO[a0+3]+1.242118006816237*uSMuO[a0+1]*uSMuO[a0+3]+0.4517539514526256*uSMuO[a0+2]*uSMuO[a0+2]+1.414213562373095*uSMuO[a0]*uSMuO[a0+2]+0.6324555320336759*uSMuO[a0+1]*uSMuO[a0+1]; 
    uSMuOSq[3] += 0.8432740427115681*uSMuO[a0+2]*uSMuO[a0+3]+1.414213562373095*uSMuO[a0]*uSMuO[a0+3]+1.242118006816237*uSMuO[a0+1]*uSMuO[a0+2]; 
  } 
 
  // ..... Get the cross flow velocity uCross ..... // 
  for (unsigned short int vd=0; vd<1; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uCross[a0] = uOther[a0]; 
    uCross[a0+1] = uOther[a0+1]; 
    uCross[a0+2] = uOther[a0+2]; 
    uCross[a0+3] = uOther[a0+3]; 
 
  } 
 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = uSMuOSq[0]*mRat+vtSqSelf[0]; 
  vtSqCross[1] = uSMuOSq[1]*mRat+vtSqSelf[1]; 
  vtSqCross[2] = uSMuOSq[2]*mRat+vtSqSelf[2]; 
  vtSqCross[3] = uSMuOSq[3]*mRat+vtSqSelf[3]; 
 
} 
 
