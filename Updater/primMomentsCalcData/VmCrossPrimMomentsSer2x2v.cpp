#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmCrossPrimMomentsGreene2x2vSer12_P1(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[8]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
  u1Mu2[4] = u1[4]-1.0*u2[4]; 
  u1Mu2[5] = u1[5]-1.0*u2[5]; 
  u1Mu2[6] = u1[6]-1.0*u2[6]; 
  u1Mu2[7] = u1[7]-1.0*u2[7]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    u1Mu2Sq[0] += 0.5*u1Mu2[a0+3]*u1Mu2[a0+3]+0.5*u1Mu2[a0+2]*u1Mu2[a0+2]+0.5*u1Mu2[a0+1]*u1Mu2[a0+1]+0.5*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += u1Mu2[a0+2]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+1]; 
    u1Mu2Sq[2] += u1Mu2[a0+1]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+2]; 
    u1Mu2Sq[3] += u1Mu2[a0]*u1Mu2[a0+3]+u1Mu2[a0+1]*u1Mu2[a0+2]; 
  } 
 
  // ..... Get the relative flow velocity u12 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uCross[a0] = 0.5*(u2[a0]+u1[a0])-0.5*u1Mu2[a0]*beta; 
    uCross[a0+1] = 0.5*(u2[a0+1]+u1[a0+1])-0.5*u1Mu2[a0+1]*beta; 
    uCross[a0+2] = 0.5*(u2[a0+2]+u1[a0+2])-0.5*u1Mu2[a0+2]*beta; 
    uCross[a0+3] = 0.5*(u2[a0+3]+u1[a0+3])-0.5*u1Mu2[a0+3]*beta; 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq12 ..... // 
  vtSqCross[0] = ((-1.0*vtSq1[0])-0.08333333333333333*u1Mu2Sq[0])*m1Dm2*mBetaFrac+(vtSq2[0]+0.25*u1Mu2Sq[0])*mBetaFrac+vtSq1[0]; 
  vtSqCross[1] = ((-1.0*vtSq1[1])-0.08333333333333333*u1Mu2Sq[1])*m1Dm2*mBetaFrac+(vtSq2[1]+0.25*u1Mu2Sq[1])*mBetaFrac+vtSq1[1]; 
  vtSqCross[2] = ((-1.0*vtSq1[2])-0.08333333333333333*u1Mu2Sq[2])*m1Dm2*mBetaFrac+(vtSq2[2]+0.25*u1Mu2Sq[2])*mBetaFrac+vtSq1[2]; 
  vtSqCross[3] = ((-1.0*vtSq1[3])-0.08333333333333333*u1Mu2Sq[3])*m1Dm2*mBetaFrac+(vtSq2[3]+0.25*u1Mu2Sq[3])*mBetaFrac+vtSq1[3]; 
 
} 
 
void VmCrossPrimMomentsGreene2x2vSer21_P1(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[8]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
  u1Mu2[4] = u1[4]-1.0*u2[4]; 
  u1Mu2[5] = u1[5]-1.0*u2[5]; 
  u1Mu2[6] = u1[6]-1.0*u2[6]; 
  u1Mu2[7] = u1[7]-1.0*u2[7]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[4]; 
  for (unsigned short int k=0; k<4; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    u1Mu2Sq[0] += 0.5*u1Mu2[a0+3]*u1Mu2[a0+3]+0.5*u1Mu2[a0+2]*u1Mu2[a0+2]+0.5*u1Mu2[a0+1]*u1Mu2[a0+1]+0.5*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += u1Mu2[a0+2]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+1]; 
    u1Mu2Sq[2] += u1Mu2[a0+1]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+2]; 
    u1Mu2Sq[3] += u1Mu2[a0]*u1Mu2[a0+3]+u1Mu2[a0+1]*u1Mu2[a0+2]; 
  } 
 
  // ..... Get the relative flow velocity u21 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    uCross[a0] = 0.5*(u1Mu2[a0]*beta+u2[a0]+u1[a0]); 
    uCross[a0+1] = 0.5*(u1Mu2[a0+1]*beta+u2[a0+1]+u1[a0+1]); 
    uCross[a0+2] = 0.5*(u1Mu2[a0+2]*beta+u2[a0+2]+u1[a0+2]); 
    uCross[a0+3] = 0.5*(u1Mu2[a0+3]*beta+u2[a0+3]+u1[a0+3]); 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq21 ..... // 
  vtSqCross[0] = (vtSq1[0]+0.25*u1Mu2Sq[0])*m1Dm2*mBetaFrac+((-1.0*vtSq2[0])-0.08333333333333333*u1Mu2Sq[0])*mBetaFrac+vtSq2[0]; 
  vtSqCross[1] = (vtSq1[1]+0.25*u1Mu2Sq[1])*m1Dm2*mBetaFrac+((-1.0*vtSq2[1])-0.08333333333333333*u1Mu2Sq[1])*mBetaFrac+vtSq2[1]; 
  vtSqCross[2] = (vtSq1[2]+0.25*u1Mu2Sq[2])*m1Dm2*mBetaFrac+((-1.0*vtSq2[2])-0.08333333333333333*u1Mu2Sq[2])*mBetaFrac+vtSq2[2]; 
  vtSqCross[3] = (vtSq1[3]+0.25*u1Mu2Sq[3])*m1Dm2*mBetaFrac+((-1.0*vtSq2[3])-0.08333333333333333*u1Mu2Sq[3])*mBetaFrac+vtSq2[3]; 
 
} 
 
void VmCrossPrimMomentsGreene2x2vSer12_P2(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[16]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
  u1Mu2[4] = u1[4]-1.0*u2[4]; 
  u1Mu2[5] = u1[5]-1.0*u2[5]; 
  u1Mu2[6] = u1[6]-1.0*u2[6]; 
  u1Mu2[7] = u1[7]-1.0*u2[7]; 
  u1Mu2[8] = u1[8]-1.0*u2[8]; 
  u1Mu2[9] = u1[9]-1.0*u2[9]; 
  u1Mu2[10] = u1[10]-1.0*u2[10]; 
  u1Mu2[11] = u1[11]-1.0*u2[11]; 
  u1Mu2[12] = u1[12]-1.0*u2[12]; 
  u1Mu2[13] = u1[13]-1.0*u2[13]; 
  u1Mu2[14] = u1[14]-1.0*u2[14]; 
  u1Mu2[15] = u1[15]-1.0*u2[15]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[8]; 
  for (unsigned short int k=0; k<8; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    u1Mu2Sq[0] += 0.5*u1Mu2[a0+7]*u1Mu2[a0+7]+0.5*u1Mu2[a0+6]*u1Mu2[a0+6]+0.5*u1Mu2[a0+5]*u1Mu2[a0+5]+0.5*u1Mu2[a0+4]*u1Mu2[a0+4]+0.5*u1Mu2[a0+3]*u1Mu2[a0+3]+0.5*u1Mu2[a0+2]*u1Mu2[a0+2]+0.5*u1Mu2[a0+1]*u1Mu2[a0+1]+0.5*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += 1.0*u1Mu2[a0+5]*u1Mu2[a0+7]+0.8944271909999161*u1Mu2[a0+3]*u1Mu2[a0+6]+0.8944271909999159*u1Mu2[a0+1]*u1Mu2[a0+4]+u1Mu2[a0+2]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+1]; 
    u1Mu2Sq[2] += 0.8944271909999161*u1Mu2[a0+3]*u1Mu2[a0+7]+1.0*u1Mu2[a0+4]*u1Mu2[a0+6]+0.8944271909999159*u1Mu2[a0+2]*u1Mu2[a0+5]+u1Mu2[a0+1]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+2]; 
    u1Mu2Sq[3] += 0.8*u1Mu2[a0+6]*u1Mu2[a0+7]+0.8944271909999161*u1Mu2[a0+2]*u1Mu2[a0+7]+0.8944271909999161*u1Mu2[a0+1]*u1Mu2[a0+6]+0.8944271909999159*u1Mu2[a0+3]*u1Mu2[a0+5]+0.8944271909999159*u1Mu2[a0+3]*u1Mu2[a0+4]+u1Mu2[a0]*u1Mu2[a0+3]+u1Mu2[a0+1]*u1Mu2[a0+2]; 
    u1Mu2Sq[4] += 0.4472135954999579*u1Mu2[a0+7]*u1Mu2[a0+7]+0.31943828249997*u1Mu2[a0+6]*u1Mu2[a0+6]+1.0*u1Mu2[a0+2]*u1Mu2[a0+6]+0.31943828249997*u1Mu2[a0+4]*u1Mu2[a0+4]+u1Mu2[a0]*u1Mu2[a0+4]+0.4472135954999579*u1Mu2[a0+3]*u1Mu2[a0+3]+0.4472135954999579*u1Mu2[a0+1]*u1Mu2[a0+1]; 
    u1Mu2Sq[5] += 0.31943828249997*u1Mu2[a0+7]*u1Mu2[a0+7]+1.0*u1Mu2[a0+1]*u1Mu2[a0+7]+0.4472135954999579*u1Mu2[a0+6]*u1Mu2[a0+6]+0.31943828249997*u1Mu2[a0+5]*u1Mu2[a0+5]+u1Mu2[a0]*u1Mu2[a0+5]+0.4472135954999579*u1Mu2[a0+3]*u1Mu2[a0+3]+0.4472135954999579*u1Mu2[a0+2]*u1Mu2[a0+2]; 
    u1Mu2Sq[6] += 0.8*u1Mu2[a0+3]*u1Mu2[a0+7]+0.8944271909999159*u1Mu2[a0+5]*u1Mu2[a0+6]+0.6388765649999399*u1Mu2[a0+4]*u1Mu2[a0+6]+u1Mu2[a0]*u1Mu2[a0+6]+1.0*u1Mu2[a0+2]*u1Mu2[a0+4]+0.8944271909999161*u1Mu2[a0+1]*u1Mu2[a0+3]; 
    u1Mu2Sq[7] += 0.6388765649999399*u1Mu2[a0+5]*u1Mu2[a0+7]+0.8944271909999159*u1Mu2[a0+4]*u1Mu2[a0+7]+u1Mu2[a0]*u1Mu2[a0+7]+0.8*u1Mu2[a0+3]*u1Mu2[a0+6]+1.0*u1Mu2[a0+1]*u1Mu2[a0+5]+0.8944271909999161*u1Mu2[a0+2]*u1Mu2[a0+3]; 
  } 
 
  // ..... Get the relative flow velocity u12 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uCross[a0] = 0.5*(u2[a0]+u1[a0])-0.5*u1Mu2[a0]*beta; 
    uCross[a0+1] = 0.5*(u2[a0+1]+u1[a0+1])-0.5*u1Mu2[a0+1]*beta; 
    uCross[a0+2] = 0.5*(u2[a0+2]+u1[a0+2])-0.5*u1Mu2[a0+2]*beta; 
    uCross[a0+3] = 0.5*(u2[a0+3]+u1[a0+3])-0.5*u1Mu2[a0+3]*beta; 
    uCross[a0+4] = 0.5*(u2[a0+4]+u1[a0+4])-0.5*u1Mu2[a0+4]*beta; 
    uCross[a0+5] = 0.5*(u2[a0+5]+u1[a0+5])-0.5*u1Mu2[a0+5]*beta; 
    uCross[a0+6] = 0.5*(u2[a0+6]+u1[a0+6])-0.5*u1Mu2[a0+6]*beta; 
    uCross[a0+7] = 0.5*(u2[a0+7]+u1[a0+7])-0.5*u1Mu2[a0+7]*beta; 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq12 ..... // 
  vtSqCross[0] = ((-1.0*vtSq1[0])-0.08333333333333333*u1Mu2Sq[0])*m1Dm2*mBetaFrac+(vtSq2[0]+0.25*u1Mu2Sq[0])*mBetaFrac+vtSq1[0]; 
  vtSqCross[1] = ((-1.0*vtSq1[1])-0.08333333333333333*u1Mu2Sq[1])*m1Dm2*mBetaFrac+(vtSq2[1]+0.25*u1Mu2Sq[1])*mBetaFrac+vtSq1[1]; 
  vtSqCross[2] = ((-1.0*vtSq1[2])-0.08333333333333333*u1Mu2Sq[2])*m1Dm2*mBetaFrac+(vtSq2[2]+0.25*u1Mu2Sq[2])*mBetaFrac+vtSq1[2]; 
  vtSqCross[3] = ((-1.0*vtSq1[3])-0.08333333333333333*u1Mu2Sq[3])*m1Dm2*mBetaFrac+(vtSq2[3]+0.25*u1Mu2Sq[3])*mBetaFrac+vtSq1[3]; 
  vtSqCross[4] = ((-1.0*vtSq1[4])-0.08333333333333333*u1Mu2Sq[4])*m1Dm2*mBetaFrac+(vtSq2[4]+0.25*u1Mu2Sq[4])*mBetaFrac+vtSq1[4]; 
  vtSqCross[5] = ((-1.0*vtSq1[5])-0.08333333333333333*u1Mu2Sq[5])*m1Dm2*mBetaFrac+(vtSq2[5]+0.25*u1Mu2Sq[5])*mBetaFrac+vtSq1[5]; 
  vtSqCross[6] = ((-1.0*vtSq1[6])-0.08333333333333333*u1Mu2Sq[6])*m1Dm2*mBetaFrac+(vtSq2[6]+0.25*u1Mu2Sq[6])*mBetaFrac+vtSq1[6]; 
  vtSqCross[7] = ((-1.0*vtSq1[7])-0.08333333333333333*u1Mu2Sq[7])*m1Dm2*mBetaFrac+(vtSq2[7]+0.25*u1Mu2Sq[7])*mBetaFrac+vtSq1[7]; 
 
} 
 
void VmCrossPrimMomentsGreene2x2vSer21_P2(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[16]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
  u1Mu2[4] = u1[4]-1.0*u2[4]; 
  u1Mu2[5] = u1[5]-1.0*u2[5]; 
  u1Mu2[6] = u1[6]-1.0*u2[6]; 
  u1Mu2[7] = u1[7]-1.0*u2[7]; 
  u1Mu2[8] = u1[8]-1.0*u2[8]; 
  u1Mu2[9] = u1[9]-1.0*u2[9]; 
  u1Mu2[10] = u1[10]-1.0*u2[10]; 
  u1Mu2[11] = u1[11]-1.0*u2[11]; 
  u1Mu2[12] = u1[12]-1.0*u2[12]; 
  u1Mu2[13] = u1[13]-1.0*u2[13]; 
  u1Mu2[14] = u1[14]-1.0*u2[14]; 
  u1Mu2[15] = u1[15]-1.0*u2[15]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[8]; 
  for (unsigned short int k=0; k<8; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    u1Mu2Sq[0] += 0.5*u1Mu2[a0+7]*u1Mu2[a0+7]+0.5*u1Mu2[a0+6]*u1Mu2[a0+6]+0.5*u1Mu2[a0+5]*u1Mu2[a0+5]+0.5*u1Mu2[a0+4]*u1Mu2[a0+4]+0.5*u1Mu2[a0+3]*u1Mu2[a0+3]+0.5*u1Mu2[a0+2]*u1Mu2[a0+2]+0.5*u1Mu2[a0+1]*u1Mu2[a0+1]+0.5*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += 1.0*u1Mu2[a0+5]*u1Mu2[a0+7]+0.8944271909999161*u1Mu2[a0+3]*u1Mu2[a0+6]+0.8944271909999159*u1Mu2[a0+1]*u1Mu2[a0+4]+u1Mu2[a0+2]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+1]; 
    u1Mu2Sq[2] += 0.8944271909999161*u1Mu2[a0+3]*u1Mu2[a0+7]+1.0*u1Mu2[a0+4]*u1Mu2[a0+6]+0.8944271909999159*u1Mu2[a0+2]*u1Mu2[a0+5]+u1Mu2[a0+1]*u1Mu2[a0+3]+u1Mu2[a0]*u1Mu2[a0+2]; 
    u1Mu2Sq[3] += 0.8*u1Mu2[a0+6]*u1Mu2[a0+7]+0.8944271909999161*u1Mu2[a0+2]*u1Mu2[a0+7]+0.8944271909999161*u1Mu2[a0+1]*u1Mu2[a0+6]+0.8944271909999159*u1Mu2[a0+3]*u1Mu2[a0+5]+0.8944271909999159*u1Mu2[a0+3]*u1Mu2[a0+4]+u1Mu2[a0]*u1Mu2[a0+3]+u1Mu2[a0+1]*u1Mu2[a0+2]; 
    u1Mu2Sq[4] += 0.4472135954999579*u1Mu2[a0+7]*u1Mu2[a0+7]+0.31943828249997*u1Mu2[a0+6]*u1Mu2[a0+6]+1.0*u1Mu2[a0+2]*u1Mu2[a0+6]+0.31943828249997*u1Mu2[a0+4]*u1Mu2[a0+4]+u1Mu2[a0]*u1Mu2[a0+4]+0.4472135954999579*u1Mu2[a0+3]*u1Mu2[a0+3]+0.4472135954999579*u1Mu2[a0+1]*u1Mu2[a0+1]; 
    u1Mu2Sq[5] += 0.31943828249997*u1Mu2[a0+7]*u1Mu2[a0+7]+1.0*u1Mu2[a0+1]*u1Mu2[a0+7]+0.4472135954999579*u1Mu2[a0+6]*u1Mu2[a0+6]+0.31943828249997*u1Mu2[a0+5]*u1Mu2[a0+5]+u1Mu2[a0]*u1Mu2[a0+5]+0.4472135954999579*u1Mu2[a0+3]*u1Mu2[a0+3]+0.4472135954999579*u1Mu2[a0+2]*u1Mu2[a0+2]; 
    u1Mu2Sq[6] += 0.8*u1Mu2[a0+3]*u1Mu2[a0+7]+0.8944271909999159*u1Mu2[a0+5]*u1Mu2[a0+6]+0.6388765649999399*u1Mu2[a0+4]*u1Mu2[a0+6]+u1Mu2[a0]*u1Mu2[a0+6]+1.0*u1Mu2[a0+2]*u1Mu2[a0+4]+0.8944271909999161*u1Mu2[a0+1]*u1Mu2[a0+3]; 
    u1Mu2Sq[7] += 0.6388765649999399*u1Mu2[a0+5]*u1Mu2[a0+7]+0.8944271909999159*u1Mu2[a0+4]*u1Mu2[a0+7]+u1Mu2[a0]*u1Mu2[a0+7]+0.8*u1Mu2[a0+3]*u1Mu2[a0+6]+1.0*u1Mu2[a0+1]*u1Mu2[a0+5]+0.8944271909999161*u1Mu2[a0+2]*u1Mu2[a0+3]; 
  } 
 
  // ..... Get the relative flow velocity u21 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    uCross[a0] = 0.5*(u1Mu2[a0]*beta+u2[a0]+u1[a0]); 
    uCross[a0+1] = 0.5*(u1Mu2[a0+1]*beta+u2[a0+1]+u1[a0+1]); 
    uCross[a0+2] = 0.5*(u1Mu2[a0+2]*beta+u2[a0+2]+u1[a0+2]); 
    uCross[a0+3] = 0.5*(u1Mu2[a0+3]*beta+u2[a0+3]+u1[a0+3]); 
    uCross[a0+4] = 0.5*(u1Mu2[a0+4]*beta+u2[a0+4]+u1[a0+4]); 
    uCross[a0+5] = 0.5*(u1Mu2[a0+5]*beta+u2[a0+5]+u1[a0+5]); 
    uCross[a0+6] = 0.5*(u1Mu2[a0+6]*beta+u2[a0+6]+u1[a0+6]); 
    uCross[a0+7] = 0.5*(u1Mu2[a0+7]*beta+u2[a0+7]+u1[a0+7]); 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq21 ..... // 
  vtSqCross[0] = (vtSq1[0]+0.25*u1Mu2Sq[0])*m1Dm2*mBetaFrac+((-1.0*vtSq2[0])-0.08333333333333333*u1Mu2Sq[0])*mBetaFrac+vtSq2[0]; 
  vtSqCross[1] = (vtSq1[1]+0.25*u1Mu2Sq[1])*m1Dm2*mBetaFrac+((-1.0*vtSq2[1])-0.08333333333333333*u1Mu2Sq[1])*mBetaFrac+vtSq2[1]; 
  vtSqCross[2] = (vtSq1[2]+0.25*u1Mu2Sq[2])*m1Dm2*mBetaFrac+((-1.0*vtSq2[2])-0.08333333333333333*u1Mu2Sq[2])*mBetaFrac+vtSq2[2]; 
  vtSqCross[3] = (vtSq1[3]+0.25*u1Mu2Sq[3])*m1Dm2*mBetaFrac+((-1.0*vtSq2[3])-0.08333333333333333*u1Mu2Sq[3])*mBetaFrac+vtSq2[3]; 
  vtSqCross[4] = (vtSq1[4]+0.25*u1Mu2Sq[4])*m1Dm2*mBetaFrac+((-1.0*vtSq2[4])-0.08333333333333333*u1Mu2Sq[4])*mBetaFrac+vtSq2[4]; 
  vtSqCross[5] = (vtSq1[5]+0.25*u1Mu2Sq[5])*m1Dm2*mBetaFrac+((-1.0*vtSq2[5])-0.08333333333333333*u1Mu2Sq[5])*mBetaFrac+vtSq2[5]; 
  vtSqCross[6] = (vtSq1[6]+0.25*u1Mu2Sq[6])*m1Dm2*mBetaFrac+((-1.0*vtSq2[6])-0.08333333333333333*u1Mu2Sq[6])*mBetaFrac+vtSq2[6]; 
  vtSqCross[7] = (vtSq1[7]+0.25*u1Mu2Sq[7])*m1Dm2*mBetaFrac+((-1.0*vtSq2[7])-0.08333333333333333*u1Mu2Sq[7])*mBetaFrac+vtSq2[7]; 
 
} 
 
