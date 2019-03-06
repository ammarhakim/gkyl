#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void GkCrossPrimMomentsGreene2x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[8]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uSMuOSq[0] += 0.5*uSMuO[a0+3]*uSMuO[a0+3]+0.5*uSMuO[a0+2]*uSMuO[a0+2]+0.5*uSMuO[a0+1]*uSMuO[a0+1]+0.5*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += uSMuO[a0+2]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += uSMuO[a0+1]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += uSMuO[a0]*uSMuO[a0+3]+uSMuO[a0+1]*uSMuO[a0+2]; 
  } 
 
  // ..... Get the cross flow velocity uSelf2 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uCross[a0] = 0.5*(uSelf[a0]+uOther[a0])-0.5*uSMuO[a0]*beta; 
    uCross[a0+1] = 0.5*(uSelf[a0+1]+uOther[a0+1])-0.5*uSMuO[a0+1]*beta; 
    uCross[a0+2] = 0.5*(uSelf[a0+2]+uOther[a0+2])-0.5*uSMuO[a0+2]*beta; 
    uCross[a0+3] = 0.5*(uSelf[a0+3]+uOther[a0+3])-0.5*uSMuO[a0+3]*beta; 
 
  } 
 
  double mBetaFrac = (0.5*(beta+1.0))/(mRat+1.0); 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = (vtSqOther[0]+0.3333333333333333*uSMuOSq[0])*mBetaFrac*mRat-1.0*vtSqSelf[0]*mBetaFrac+vtSqSelf[0]; 
  vtSqCross[1] = (vtSqOther[1]+0.3333333333333333*uSMuOSq[1])*mBetaFrac*mRat-1.0*vtSqSelf[1]*mBetaFrac+vtSqSelf[1]; 
  vtSqCross[2] = (vtSqOther[2]+0.3333333333333333*uSMuOSq[2])*mBetaFrac*mRat-1.0*vtSqSelf[2]*mBetaFrac+vtSqSelf[2]; 
  vtSqCross[3] = (vtSqOther[3]+0.3333333333333333*uSMuOSq[3])*mBetaFrac*mRat-1.0*vtSqSelf[3]*mBetaFrac+vtSqSelf[3]; 
 
} 
 
void GkCrossPrimMomentsGreene2x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[16]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
  uSMuO[8] = uSelf[8]-1.0*uOther[8]; 
  uSMuO[9] = uSelf[9]-1.0*uOther[9]; 
  uSMuO[10] = uSelf[10]-1.0*uOther[10]; 
  uSMuO[11] = uSelf[11]-1.0*uOther[11]; 
  uSMuO[12] = uSelf[12]-1.0*uOther[12]; 
  uSMuO[13] = uSelf[13]-1.0*uOther[13]; 
  uSMuO[14] = uSelf[14]-1.0*uOther[14]; 
  uSMuO[15] = uSelf[15]-1.0*uOther[15]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[8]; 
  for (unsigned short int k=0; k<8; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uSMuOSq[0] += 0.5*uSMuO[a0+7]*uSMuO[a0+7]+0.5*uSMuO[a0+6]*uSMuO[a0+6]+0.5*uSMuO[a0+5]*uSMuO[a0+5]+0.5*uSMuO[a0+4]*uSMuO[a0+4]+0.5*uSMuO[a0+3]*uSMuO[a0+3]+0.5*uSMuO[a0+2]*uSMuO[a0+2]+0.5*uSMuO[a0+1]*uSMuO[a0+1]+0.5*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.0*uSMuO[a0+5]*uSMuO[a0+7]+0.8944271909999161*uSMuO[a0+3]*uSMuO[a0+6]+0.8944271909999159*uSMuO[a0+1]*uSMuO[a0+4]+uSMuO[a0+2]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.8944271909999161*uSMuO[a0+3]*uSMuO[a0+7]+1.0*uSMuO[a0+4]*uSMuO[a0+6]+0.8944271909999159*uSMuO[a0+2]*uSMuO[a0+5]+uSMuO[a0+1]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += 0.8*uSMuO[a0+6]*uSMuO[a0+7]+0.8944271909999161*uSMuO[a0+2]*uSMuO[a0+7]+0.8944271909999161*uSMuO[a0+1]*uSMuO[a0+6]+0.8944271909999159*uSMuO[a0+3]*uSMuO[a0+5]+0.8944271909999159*uSMuO[a0+3]*uSMuO[a0+4]+uSMuO[a0]*uSMuO[a0+3]+uSMuO[a0+1]*uSMuO[a0+2]; 
    uSMuOSq[4] += 0.4472135954999579*uSMuO[a0+7]*uSMuO[a0+7]+0.31943828249997*uSMuO[a0+6]*uSMuO[a0+6]+1.0*uSMuO[a0+2]*uSMuO[a0+6]+0.31943828249997*uSMuO[a0+4]*uSMuO[a0+4]+uSMuO[a0]*uSMuO[a0+4]+0.4472135954999579*uSMuO[a0+3]*uSMuO[a0+3]+0.4472135954999579*uSMuO[a0+1]*uSMuO[a0+1]; 
    uSMuOSq[5] += 0.31943828249997*uSMuO[a0+7]*uSMuO[a0+7]+1.0*uSMuO[a0+1]*uSMuO[a0+7]+0.4472135954999579*uSMuO[a0+6]*uSMuO[a0+6]+0.31943828249997*uSMuO[a0+5]*uSMuO[a0+5]+uSMuO[a0]*uSMuO[a0+5]+0.4472135954999579*uSMuO[a0+3]*uSMuO[a0+3]+0.4472135954999579*uSMuO[a0+2]*uSMuO[a0+2]; 
    uSMuOSq[6] += 0.8*uSMuO[a0+3]*uSMuO[a0+7]+0.8944271909999159*uSMuO[a0+5]*uSMuO[a0+6]+0.6388765649999399*uSMuO[a0+4]*uSMuO[a0+6]+uSMuO[a0]*uSMuO[a0+6]+1.0*uSMuO[a0+2]*uSMuO[a0+4]+0.8944271909999161*uSMuO[a0+1]*uSMuO[a0+3]; 
    uSMuOSq[7] += 0.6388765649999399*uSMuO[a0+5]*uSMuO[a0+7]+0.8944271909999159*uSMuO[a0+4]*uSMuO[a0+7]+uSMuO[a0]*uSMuO[a0+7]+0.8*uSMuO[a0+3]*uSMuO[a0+6]+1.0*uSMuO[a0+1]*uSMuO[a0+5]+0.8944271909999161*uSMuO[a0+2]*uSMuO[a0+3]; 
  } 
 
  // ..... Get the cross flow velocity uSelf2 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uCross[a0] = 0.5*(uSelf[a0]+uOther[a0])-0.5*uSMuO[a0]*beta; 
    uCross[a0+1] = 0.5*(uSelf[a0+1]+uOther[a0+1])-0.5*uSMuO[a0+1]*beta; 
    uCross[a0+2] = 0.5*(uSelf[a0+2]+uOther[a0+2])-0.5*uSMuO[a0+2]*beta; 
    uCross[a0+3] = 0.5*(uSelf[a0+3]+uOther[a0+3])-0.5*uSMuO[a0+3]*beta; 
    uCross[a0+4] = 0.5*(uSelf[a0+4]+uOther[a0+4])-0.5*uSMuO[a0+4]*beta; 
    uCross[a0+5] = 0.5*(uSelf[a0+5]+uOther[a0+5])-0.5*uSMuO[a0+5]*beta; 
    uCross[a0+6] = 0.5*(uSelf[a0+6]+uOther[a0+6])-0.5*uSMuO[a0+6]*beta; 
    uCross[a0+7] = 0.5*(uSelf[a0+7]+uOther[a0+7])-0.5*uSMuO[a0+7]*beta; 
 
  } 
 
  double mBetaFrac = (0.5*(beta+1.0))/(mRat+1.0); 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = (vtSqOther[0]+0.3333333333333333*uSMuOSq[0])*mBetaFrac*mRat-1.0*vtSqSelf[0]*mBetaFrac+vtSqSelf[0]; 
  vtSqCross[1] = (vtSqOther[1]+0.3333333333333333*uSMuOSq[1])*mBetaFrac*mRat-1.0*vtSqSelf[1]*mBetaFrac+vtSqSelf[1]; 
  vtSqCross[2] = (vtSqOther[2]+0.3333333333333333*uSMuOSq[2])*mBetaFrac*mRat-1.0*vtSqSelf[2]*mBetaFrac+vtSqSelf[2]; 
  vtSqCross[3] = (vtSqOther[3]+0.3333333333333333*uSMuOSq[3])*mBetaFrac*mRat-1.0*vtSqSelf[3]*mBetaFrac+vtSqSelf[3]; 
  vtSqCross[4] = (vtSqOther[4]+0.3333333333333333*uSMuOSq[4])*mBetaFrac*mRat-1.0*vtSqSelf[4]*mBetaFrac+vtSqSelf[4]; 
  vtSqCross[5] = (vtSqOther[5]+0.3333333333333333*uSMuOSq[5])*mBetaFrac*mRat-1.0*vtSqSelf[5]*mBetaFrac+vtSqSelf[5]; 
  vtSqCross[6] = (vtSqOther[6]+0.3333333333333333*uSMuOSq[6])*mBetaFrac*mRat-1.0*vtSqSelf[6]*mBetaFrac+vtSqSelf[6]; 
  vtSqCross[7] = (vtSqOther[7]+0.3333333333333333*uSMuOSq[7])*mBetaFrac*mRat-1.0*vtSqSelf[7]*mBetaFrac+vtSqSelf[7]; 
 
} 
 
void GkCrossPrimMomentsHeavyIon2x2vSer_P1(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[8]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uSMuOSq[0] += 0.5*uSMuO[a0+3]*uSMuO[a0+3]+0.5*uSMuO[a0+2]*uSMuO[a0+2]+0.5*uSMuO[a0+1]*uSMuO[a0+1]+0.5*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += uSMuO[a0+2]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += uSMuO[a0+1]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += uSMuO[a0]*uSMuO[a0+3]+uSMuO[a0+1]*uSMuO[a0+2]; 
  } 
 
  // ..... Get the cross flow velocity uCross ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uCross[a0] = uOther[a0]; 
    uCross[a0+1] = uOther[a0+1]; 
    uCross[a0+2] = uOther[a0+2]; 
    uCross[a0+3] = uOther[a0+3]; 
 
  } 
 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = 0.3333333333333333*uSMuOSq[0]*mRat+vtSqSelf[0]; 
  vtSqCross[1] = 0.3333333333333333*uSMuOSq[1]*mRat+vtSqSelf[1]; 
  vtSqCross[2] = 0.3333333333333333*uSMuOSq[2]*mRat+vtSqSelf[2]; 
  vtSqCross[3] = 0.3333333333333333*uSMuOSq[3]*mRat+vtSqSelf[3]; 
 
} 
 
void GkCrossPrimMomentsHeavyIon2x2vSer_P2(const double mRat, const double beta, const double *uSelf, const double *vtSqSelf, const double *uOther, const double *vtSqOther, double *uCross, double *vtSqCross) 
{ 
  // mRat:              mass ratio = m_other/m_self. 
  // uSelf, vtSqSelf:   bulk flow velocity and T/m of self species. 
  // uOther, vtSqOther: bulk flow velocity and T/m of other species. 
  // uCross:            bulk flow velocity for cross-species collision term. 
  // vtSqCross:         squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity uSelf-uOther ..... // 
  double uSMuO[16]; 
  uSMuO[0] = uSelf[0]-1.0*uOther[0]; 
  uSMuO[1] = uSelf[1]-1.0*uOther[1]; 
  uSMuO[2] = uSelf[2]-1.0*uOther[2]; 
  uSMuO[3] = uSelf[3]-1.0*uOther[3]; 
  uSMuO[4] = uSelf[4]-1.0*uOther[4]; 
  uSMuO[5] = uSelf[5]-1.0*uOther[5]; 
  uSMuO[6] = uSelf[6]-1.0*uOther[6]; 
  uSMuO[7] = uSelf[7]-1.0*uOther[7]; 
  uSMuO[8] = uSelf[8]-1.0*uOther[8]; 
  uSMuO[9] = uSelf[9]-1.0*uOther[9]; 
  uSMuO[10] = uSelf[10]-1.0*uOther[10]; 
  uSMuO[11] = uSelf[11]-1.0*uOther[11]; 
  uSMuO[12] = uSelf[12]-1.0*uOther[12]; 
  uSMuO[13] = uSelf[13]-1.0*uOther[13]; 
  uSMuO[14] = uSelf[14]-1.0*uOther[14]; 
  uSMuO[15] = uSelf[15]-1.0*uOther[15]; 
 
  // ..... Get the relative speed squared (uSelf-uOther)^2 ..... // 
  double uSMuOSq[8]; 
  for (unsigned short int k=0; k<8; k++) 
  { 
    uSMuOSq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uSMuOSq[0] += 0.5*uSMuO[a0+7]*uSMuO[a0+7]+0.5*uSMuO[a0+6]*uSMuO[a0+6]+0.5*uSMuO[a0+5]*uSMuO[a0+5]+0.5*uSMuO[a0+4]*uSMuO[a0+4]+0.5*uSMuO[a0+3]*uSMuO[a0+3]+0.5*uSMuO[a0+2]*uSMuO[a0+2]+0.5*uSMuO[a0+1]*uSMuO[a0+1]+0.5*uSMuO[a0]*uSMuO[a0]; 
    uSMuOSq[1] += 1.0*uSMuO[a0+5]*uSMuO[a0+7]+0.8944271909999161*uSMuO[a0+3]*uSMuO[a0+6]+0.8944271909999159*uSMuO[a0+1]*uSMuO[a0+4]+uSMuO[a0+2]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+1]; 
    uSMuOSq[2] += 0.8944271909999161*uSMuO[a0+3]*uSMuO[a0+7]+1.0*uSMuO[a0+4]*uSMuO[a0+6]+0.8944271909999159*uSMuO[a0+2]*uSMuO[a0+5]+uSMuO[a0+1]*uSMuO[a0+3]+uSMuO[a0]*uSMuO[a0+2]; 
    uSMuOSq[3] += 0.8*uSMuO[a0+6]*uSMuO[a0+7]+0.8944271909999161*uSMuO[a0+2]*uSMuO[a0+7]+0.8944271909999161*uSMuO[a0+1]*uSMuO[a0+6]+0.8944271909999159*uSMuO[a0+3]*uSMuO[a0+5]+0.8944271909999159*uSMuO[a0+3]*uSMuO[a0+4]+uSMuO[a0]*uSMuO[a0+3]+uSMuO[a0+1]*uSMuO[a0+2]; 
    uSMuOSq[4] += 0.4472135954999579*uSMuO[a0+7]*uSMuO[a0+7]+0.31943828249997*uSMuO[a0+6]*uSMuO[a0+6]+1.0*uSMuO[a0+2]*uSMuO[a0+6]+0.31943828249997*uSMuO[a0+4]*uSMuO[a0+4]+uSMuO[a0]*uSMuO[a0+4]+0.4472135954999579*uSMuO[a0+3]*uSMuO[a0+3]+0.4472135954999579*uSMuO[a0+1]*uSMuO[a0+1]; 
    uSMuOSq[5] += 0.31943828249997*uSMuO[a0+7]*uSMuO[a0+7]+1.0*uSMuO[a0+1]*uSMuO[a0+7]+0.4472135954999579*uSMuO[a0+6]*uSMuO[a0+6]+0.31943828249997*uSMuO[a0+5]*uSMuO[a0+5]+uSMuO[a0]*uSMuO[a0+5]+0.4472135954999579*uSMuO[a0+3]*uSMuO[a0+3]+0.4472135954999579*uSMuO[a0+2]*uSMuO[a0+2]; 
    uSMuOSq[6] += 0.8*uSMuO[a0+3]*uSMuO[a0+7]+0.8944271909999159*uSMuO[a0+5]*uSMuO[a0+6]+0.6388765649999399*uSMuO[a0+4]*uSMuO[a0+6]+uSMuO[a0]*uSMuO[a0+6]+1.0*uSMuO[a0+2]*uSMuO[a0+4]+0.8944271909999161*uSMuO[a0+1]*uSMuO[a0+3]; 
    uSMuOSq[7] += 0.6388765649999399*uSMuO[a0+5]*uSMuO[a0+7]+0.8944271909999159*uSMuO[a0+4]*uSMuO[a0+7]+uSMuO[a0]*uSMuO[a0+7]+0.8*uSMuO[a0+3]*uSMuO[a0+6]+1.0*uSMuO[a0+1]*uSMuO[a0+5]+0.8944271909999161*uSMuO[a0+2]*uSMuO[a0+3]; 
  } 
 
  // ..... Get the cross flow velocity uCross ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uCross[a0] = uOther[a0]; 
    uCross[a0+1] = uOther[a0+1]; 
    uCross[a0+2] = uOther[a0+2]; 
    uCross[a0+3] = uOther[a0+3]; 
    uCross[a0+4] = uOther[a0+4]; 
    uCross[a0+5] = uOther[a0+5]; 
    uCross[a0+6] = uOther[a0+6]; 
    uCross[a0+7] = uOther[a0+7]; 
 
  } 
 
  // ..... Get the cross thermal speed squared vtSqCross ..... // 
  vtSqCross[0] = 0.3333333333333333*uSMuOSq[0]*mRat+vtSqSelf[0]; 
  vtSqCross[1] = 0.3333333333333333*uSMuOSq[1]*mRat+vtSqSelf[1]; 
  vtSqCross[2] = 0.3333333333333333*uSMuOSq[2]*mRat+vtSqSelf[2]; 
  vtSqCross[3] = 0.3333333333333333*uSMuOSq[3]*mRat+vtSqSelf[3]; 
  vtSqCross[4] = 0.3333333333333333*uSMuOSq[4]*mRat+vtSqSelf[4]; 
  vtSqCross[5] = 0.3333333333333333*uSMuOSq[5]*mRat+vtSqSelf[5]; 
  vtSqCross[6] = 0.3333333333333333*uSMuOSq[6]*mRat+vtSqSelf[6]; 
  vtSqCross[7] = 0.3333333333333333*uSMuOSq[7]*mRat+vtSqSelf[7]; 
 
} 
 
