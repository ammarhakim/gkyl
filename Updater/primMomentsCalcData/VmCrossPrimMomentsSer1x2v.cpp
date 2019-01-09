#include <PrimMomentsModDecl.h> 
 
using namespace Eigen; 
 
void VmCrossPrimMomentsGreene1x2vSer12_P1(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[4]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[2]; 
  for (unsigned short int k=0; k<2; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    u1Mu2Sq[0] += 0.7071067811865475*u1Mu2[a0+1]*u1Mu2[a0+1]+0.7071067811865475*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += 1.414213562373095*u1Mu2[a0]*u1Mu2[a0+1]; 
  } 
 
  // ..... Get the relative flow velocity u12 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uCross[a0] = 0.5*(u2[a0]+u1[a0])-0.5*u1Mu2[a0]*beta; 
    uCross[a0+1] = 0.5*(u2[a0+1]+u1[a0+1])-0.5*u1Mu2[a0+1]*beta; 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq12 ..... // 
  vtSqCross[0] = ((-1.0*vtSq1[0])-0.08333333333333333*u1Mu2Sq[0])*m1Dm2*mBetaFrac+(vtSq2[0]+0.25*u1Mu2Sq[0])*mBetaFrac+vtSq1[0]; 
  vtSqCross[1] = ((-1.0*vtSq1[1])-0.08333333333333333*u1Mu2Sq[1])*m1Dm2*mBetaFrac+(vtSq2[1]+0.25*u1Mu2Sq[1])*mBetaFrac+vtSq1[1]; 
 
} 
 
void VmCrossPrimMomentsGreene1x2vSer21_P1(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[4]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[2]; 
  for (unsigned short int k=0; k<2; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    u1Mu2Sq[0] += 0.7071067811865475*u1Mu2[a0+1]*u1Mu2[a0+1]+0.7071067811865475*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += 1.414213562373095*u1Mu2[a0]*u1Mu2[a0+1]; 
  } 
 
  // ..... Get the relative flow velocity u21 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    uCross[a0] = 0.5*(u1Mu2[a0]*beta+u2[a0]+u1[a0]); 
    uCross[a0+1] = 0.5*(u1Mu2[a0+1]*beta+u2[a0+1]+u1[a0+1]); 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq21 ..... // 
  vtSqCross[0] = (vtSq1[0]+0.25*u1Mu2Sq[0])*m1Dm2*mBetaFrac+((-1.0*vtSq2[0])-0.08333333333333333*u1Mu2Sq[0])*mBetaFrac+vtSq2[0]; 
  vtSqCross[1] = (vtSq1[1]+0.25*u1Mu2Sq[1])*m1Dm2*mBetaFrac+((-1.0*vtSq2[1])-0.08333333333333333*u1Mu2Sq[1])*mBetaFrac+vtSq2[1]; 
 
} 
 
void VmCrossPrimMomentsGreene1x2vSer12_P2(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[6]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
  u1Mu2[4] = u1[4]-1.0*u2[4]; 
  u1Mu2[5] = u1[5]-1.0*u2[5]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[3]; 
  for (unsigned short int k=0; k<3; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    u1Mu2Sq[0] += 0.7071067811865475*u1Mu2[a0+2]*u1Mu2[a0+2]+0.7071067811865475*u1Mu2[a0+1]*u1Mu2[a0+1]+0.7071067811865475*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += 1.264911064067352*u1Mu2[a0+1]*u1Mu2[a0+2]+1.414213562373095*u1Mu2[a0]*u1Mu2[a0+1]; 
    u1Mu2Sq[2] += 0.4517539514526256*u1Mu2[a0+2]*u1Mu2[a0+2]+1.414213562373095*u1Mu2[a0]*u1Mu2[a0+2]+0.6324555320336759*u1Mu2[a0+1]*u1Mu2[a0+1]; 
  } 
 
  // ..... Get the relative flow velocity u12 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uCross[a0] = 0.5*(u2[a0]+u1[a0])-0.5*u1Mu2[a0]*beta; 
    uCross[a0+1] = 0.5*(u2[a0+1]+u1[a0+1])-0.5*u1Mu2[a0+1]*beta; 
    uCross[a0+2] = 0.5*(u2[a0+2]+u1[a0+2])-0.5*u1Mu2[a0+2]*beta; 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq12 ..... // 
  vtSqCross[0] = ((-1.0*vtSq1[0])-0.08333333333333333*u1Mu2Sq[0])*m1Dm2*mBetaFrac+(vtSq2[0]+0.25*u1Mu2Sq[0])*mBetaFrac+vtSq1[0]; 
  vtSqCross[1] = ((-1.0*vtSq1[1])-0.08333333333333333*u1Mu2Sq[1])*m1Dm2*mBetaFrac+(vtSq2[1]+0.25*u1Mu2Sq[1])*mBetaFrac+vtSq1[1]; 
  vtSqCross[2] = ((-1.0*vtSq1[2])-0.08333333333333333*u1Mu2Sq[2])*m1Dm2*mBetaFrac+(vtSq2[2]+0.25*u1Mu2Sq[2])*mBetaFrac+vtSq1[2]; 
 
} 
 
void VmCrossPrimMomentsGreene1x2vSer21_P2(const double m1Dm2, const double beta, const double *n1, const double *u1, const double *vtSq1, const double *n2, const double *u2, const double *vtSq2, double *uCross, double *vtSqCross) 
{ 
  // mRat:          mass ratio = m_1/m_2. 
  // n1, u1, vtSq1: number density, bulk flow velocity and T_1/m_1 of first species. 
  // n2, u2, vtSq2: number density, bulk flow velocity and T_1/m_1 of second species. 
  // uCross:        bulk flow velocity for cross-species collision term. 
  // vtSqCross:     squared thermal speed for cross-species collision term. 
 
  // ..... Compute and save the relative velocity u1-u2 ..... // 
  double u1Mu2[6]; 
  u1Mu2[0] = u1[0]-1.0*u2[0]; 
  u1Mu2[1] = u1[1]-1.0*u2[1]; 
  u1Mu2[2] = u1[2]-1.0*u2[2]; 
  u1Mu2[3] = u1[3]-1.0*u2[3]; 
  u1Mu2[4] = u1[4]-1.0*u2[4]; 
  u1Mu2[5] = u1[5]-1.0*u2[5]; 
 
  // ..... Get the relative speed squared (u1-u2)^2 ..... // 
  double u1Mu2Sq[3]; 
  for (unsigned short int k=0; k<3; k++) 
  { 
    u1Mu2Sq[k] = 0.0; 
  } 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    u1Mu2Sq[0] += 0.7071067811865475*u1Mu2[a0+2]*u1Mu2[a0+2]+0.7071067811865475*u1Mu2[a0+1]*u1Mu2[a0+1]+0.7071067811865475*u1Mu2[a0]*u1Mu2[a0]; 
    u1Mu2Sq[1] += 1.264911064067352*u1Mu2[a0+1]*u1Mu2[a0+2]+1.414213562373095*u1Mu2[a0]*u1Mu2[a0+1]; 
    u1Mu2Sq[2] += 0.4517539514526256*u1Mu2[a0+2]*u1Mu2[a0+2]+1.414213562373095*u1Mu2[a0]*u1Mu2[a0+2]+0.6324555320336759*u1Mu2[a0+1]*u1Mu2[a0+1]; 
  } 
 
  // ..... Get the relative flow velocity u21 ..... // 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    uCross[a0] = 0.5*(u1Mu2[a0]*beta+u2[a0]+u1[a0]); 
    uCross[a0+1] = 0.5*(u1Mu2[a0+1]*beta+u2[a0+1]+u1[a0+1]); 
    uCross[a0+2] = 0.5*(u1Mu2[a0+2]*beta+u2[a0+2]+u1[a0+2]); 
 
  } 
 
  double mBetaFrac = (beta+1.0)/(m1Dm2+1.0); 
  // ..... Get the relative thermal speed squared vtSq21 ..... // 
  vtSqCross[0] = (vtSq1[0]+0.25*u1Mu2Sq[0])*m1Dm2*mBetaFrac+((-1.0*vtSq2[0])-0.08333333333333333*u1Mu2Sq[0])*mBetaFrac+vtSq2[0]; 
  vtSqCross[1] = (vtSq1[1]+0.25*u1Mu2Sq[1])*m1Dm2*mBetaFrac+((-1.0*vtSq2[1])-0.08333333333333333*u1Mu2Sq[1])*mBetaFrac+vtSq2[1]; 
  vtSqCross[2] = (vtSq1[2]+0.25*u1Mu2Sq[2])*m1Dm2*mBetaFrac+((-1.0*vtSq2[2])-0.08333333333333333*u1Mu2Sq[2])*mBetaFrac+vtSq2[2]; 
 
} 
 
