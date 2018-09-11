#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x2vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[4]; 
  tmp[0] = 0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[0]*B[3]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<4; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x2vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[10]; 
  tmp[0] = 0.7071067811865475*A[2]*B[7]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6324555320336759*A[1]*B[7]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[1]*B[4]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[1]*B[5]+0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.6324555320336759*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.7071067811865475*A[1]*B[2]; 
  tmp[5] = 0.6324555320336759*A[2]*B[5]+0.7071067811865475*A[0]*B[5]+0.7071067811865475*A[1]*B[3]; 
  tmp[6] = 0.7071067811865475*A[0]*B[6]; 
  tmp[7] = 0.4517539514526256*A[2]*B[7]+0.7071067811865475*A[0]*B[7]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[8] = 0.7071067811865475*A[0]*B[8]; 
  tmp[9] = 0.7071067811865475*A[0]*B[9]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<10; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x2vMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[20]; 
  tmp[0] = 0.7071067811865475*A[3]*B[17]+0.7071067811865475*A[2]*B[7]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6210590034081186*A[2]*B[17]+0.6210590034081186*A[3]*B[7]+0.6324555320336759*A[1]*B[7]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[2]*B[11]+0.7071067811865475*A[1]*B[4]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[2]*B[13]+0.7071067811865475*A[1]*B[5]+0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.6210590034081187*A[3]*B[11]+0.632455532033676*A[1]*B[11]+0.6324555320336759*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.7071067811865475*A[1]*B[2]; 
  tmp[5] = 0.6210590034081187*A[3]*B[13]+0.632455532033676*A[1]*B[13]+0.6324555320336759*A[2]*B[5]+0.7071067811865475*A[0]*B[5]+0.7071067811865475*A[1]*B[3]; 
  tmp[6] = 0.7071067811865475*A[1]*B[10]+0.7071067811865475*A[0]*B[6]; 
  tmp[7] = 0.421637021355784*A[3]*B[17]+0.6210590034081186*A[1]*B[17]+0.4517539514526256*A[2]*B[7]+0.7071067811865475*A[0]*B[7]+0.6210590034081186*B[1]*A[3]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[8] = 0.7071067811865475*A[1]*B[12]+0.7071067811865475*A[0]*B[8]; 
  tmp[9] = 0.7071067811865475*A[1]*B[15]+0.7071067811865475*A[0]*B[9]; 
  tmp[10] = 0.6324555320336759*A[2]*B[10]+0.7071067811865475*A[0]*B[10]+0.7071067811865475*A[1]*B[6]; 
  tmp[11] = 0.4517539514526256*A[2]*B[11]+0.7071067811865475*A[0]*B[11]+0.6210590034081187*A[3]*B[4]+0.632455532033676*A[1]*B[4]+0.7071067811865475*A[2]*B[2]; 
  tmp[12] = 0.6324555320336759*A[2]*B[12]+0.7071067811865475*A[0]*B[12]+0.7071067811865475*A[1]*B[8]; 
  tmp[13] = 0.4517539514526256*A[2]*B[13]+0.7071067811865475*A[0]*B[13]+0.6210590034081187*A[3]*B[5]+0.632455532033676*A[1]*B[5]+0.7071067811865475*A[2]*B[3]; 
  tmp[14] = 0.7071067811865475*A[0]*B[14]; 
  tmp[15] = 0.6324555320336759*A[2]*B[15]+0.7071067811865475*A[0]*B[15]+0.7071067811865475*A[1]*B[9]; 
  tmp[16] = 0.7071067811865475*A[0]*B[16]; 
  tmp[17] = 0.421637021355784*A[2]*B[17]+0.7071067811865475*A[0]*B[17]+0.421637021355784*A[3]*B[7]+0.6210590034081186*A[1]*B[7]+0.7071067811865475*B[0]*A[3]+0.6210590034081186*B[1]*A[2]; 
  tmp[18] = 0.7071067811865475*A[0]*B[18]; 
  tmp[19] = 0.7071067811865475*A[0]*B[19]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<20; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide1x2vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.7071067811865475*A[0]-1.224744871391589*A[1] < 0) { 
    avgA = true;
  }
  if (1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  double As[2]; 
  double Bs[4]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(4,4); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(4);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(4);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*As[0]; 
  AEM(0,1) = 0.7071067811865475*As[1]; 
  AEM(0,2) = 0.7071067811865475*As[1]; 
  AEM(0,3) = 0.7071067811865475*As[0]; 
 
  // Fill BEV. 
  BEV << Bs[0],Bs[1],Bs[2],Bs[3]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,4,1) = u; 
 
} 
 
void CartFieldBinOpDivide1x2vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
  if (1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  double As[3]; 
  double Bs[10]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = 0.0; 
    Bs[5] = 0.0; 
    Bs[6] = B[6]; 
    Bs[7] = 0.0; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
    Bs[7] = B[7]; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(10,10); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(10);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(10);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*As[0]; 
  AEM(0,1) = 0.7071067811865475*As[1]; 
  AEM(0,3) = 0.7071067811865475*As[1]; 
  AEM(0,4) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  AEM(0,8) = 0.7071067811865475*As[0]; 
  AEM(1,4) = 0.7071067811865475*As[1]; 
  AEM(2,1) = 0.7071067811865475*As[2]; 
  AEM(2,2) = 0.6324555320336759*As[1]; 
 
  // Fill BEV. 
  BEV << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,10,1) = u; 
 
} 
 
void CartFieldBinOpDivide1x2vMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.870828693386971*A[3])+1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
  if (1.870828693386971*A[3]+1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  double As[4]; 
  double Bs[20]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = 0.0; 
    Bs[5] = 0.0; 
    Bs[6] = B[6]; 
    Bs[7] = 0.0; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
    Bs[10] = 0.0; 
    Bs[11] = 0.0; 
    Bs[12] = 0.0; 
    Bs[13] = 0.0; 
    Bs[14] = B[14]; 
    Bs[15] = 0.0; 
    Bs[16] = B[16]; 
    Bs[17] = 0.0; 
    Bs[18] = B[18]; 
    Bs[19] = B[19]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
    Bs[7] = B[7]; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
    Bs[10] = B[10]; 
    Bs[11] = B[11]; 
    Bs[12] = B[12]; 
    Bs[13] = B[13]; 
    Bs[14] = B[14]; 
    Bs[15] = B[15]; 
    Bs[16] = B[16]; 
    Bs[17] = B[17]; 
    Bs[18] = B[18]; 
    Bs[19] = B[19]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(20,20); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(20);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(20);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*As[0]; 
  AEM(0,1) = 0.7071067811865475*As[1]; 
  AEM(0,4) = 0.7071067811865475*As[1]; 
  AEM(0,5) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  AEM(0,10) = 0.7071067811865475*As[0]; 
  AEM(0,15) = 0.7071067811865475*As[0]; 
  AEM(0,18) = 0.7071067811865475*As[1]; 
  AEM(1,3) = 0.7071067811865475*As[1]; 
  AEM(1,8) = 0.7071067811865475*As[2]; 
  AEM(1,9) = 0.6210590034081186*As[3]+0.6324555320336759*As[1]; 
  AEM(2,6) = 0.7071067811865475*As[2]; 
  AEM(2,15) = 0.7071067811865475*As[2]; 
  AEM(3,8) = 0.7071067811865475*As[3]; 
  AEM(3,9) = 0.6210590034081186*As[2]; 
 
  // Fill BEV. 
  BEV << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9],Bs[10],Bs[11],Bs[12],Bs[13],Bs[14],Bs[15],Bs[16],Bs[17],Bs[18],Bs[19]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,20,1) = u; 
 
} 
 
