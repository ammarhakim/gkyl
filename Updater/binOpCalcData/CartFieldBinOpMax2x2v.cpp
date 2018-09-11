#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2x2vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[5]; 
  tmp[0] = 0.5*A[2]*B[2]+0.5*A[1]*B[1]+0.5*A[0]*B[0]; 
  tmp[1] = 0.5*A[0]*B[1]+0.5*B[0]*A[1]; 
  tmp[2] = 0.5*A[0]*B[2]+0.5*B[0]*A[2]; 
  tmp[3] = 0.5*A[0]*B[3]; 
  tmp[4] = 0.5*A[0]*B[4]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<5; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply2x2vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[15]; 
  tmp[0] = 0.5*A[5]*B[12]+0.5*A[4]*B[11]+0.5*A[3]*B[5]+0.5*A[2]*B[2]+0.5*A[1]*B[1]+0.5*A[0]*B[0]; 
  tmp[1] = 0.4472135954999579*A[1]*B[11]+0.5*A[2]*B[5]+0.4472135954999579*B[1]*A[4]+0.5*B[2]*A[3]+0.5*A[0]*B[1]+0.5*B[0]*A[1]; 
  tmp[2] = 0.4472135954999579*A[2]*B[12]+0.5*A[1]*B[5]+0.4472135954999579*B[2]*A[5]+0.5*B[1]*A[3]+0.5*A[0]*B[2]+0.5*B[0]*A[2]; 
  tmp[3] = 0.5*A[2]*B[7]+0.5*A[1]*B[6]+0.5*A[0]*B[3]; 
  tmp[4] = 0.5*A[2]*B[9]+0.5*A[1]*B[8]+0.5*A[0]*B[4]; 
  tmp[5] = 0.4472135954999579*A[3]*B[12]+0.4472135954999579*A[3]*B[11]+0.4472135954999579*A[5]*B[5]+0.4472135954999579*A[4]*B[5]+0.5*A[0]*B[5]+0.5*B[0]*A[3]+0.5*A[1]*B[2]+0.5*B[1]*A[2]; 
  tmp[6] = 0.5*A[3]*B[7]+0.4472135954999579*A[4]*B[6]+0.5*A[0]*B[6]+0.5*A[1]*B[3]; 
  tmp[7] = 0.4472135954999579*A[5]*B[7]+0.5*A[0]*B[7]+0.5*A[3]*B[6]+0.5*A[2]*B[3]; 
  tmp[8] = 0.5*A[3]*B[9]+0.4472135954999579*A[4]*B[8]+0.5*A[0]*B[8]+0.5*A[1]*B[4]; 
  tmp[9] = 0.4472135954999579*A[5]*B[9]+0.5*A[0]*B[9]+0.5*A[3]*B[8]+0.5*A[2]*B[4]; 
  tmp[10] = 0.5*A[0]*B[10]; 
  tmp[11] = 0.31943828249997*A[4]*B[11]+0.5*A[0]*B[11]+0.4472135954999579*A[3]*B[5]+0.5*B[0]*A[4]+0.4472135954999579*A[1]*B[1]; 
  tmp[12] = 0.31943828249997*A[5]*B[12]+0.5*A[0]*B[12]+0.4472135954999579*A[3]*B[5]+0.5*B[0]*A[5]+0.4472135954999579*A[2]*B[2]; 
  tmp[13] = 0.5*A[0]*B[13]; 
  tmp[14] = 0.5*A[0]*B[14]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<15; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide2x2vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-0.8660254037844386*A[2])-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.8660254037844386*A[2])-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.8660254037844386*A[2])+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.8660254037844386*A[2])+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
 
  double As[3]; 
  double Bs[5]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = 0.0; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(5,5); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(5);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(5);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*As[0]; 
  AEM(0,1) = 0.5*As[1]; 
  AEM(0,2) = 0.5*As[2]; 
  AEM(0,3) = 0.5*As[1]; 
  AEM(0,4) = 0.5*As[0]; 
  AEM(1,1) = 0.5*As[2]; 
  AEM(1,3) = 0.5*As[0]; 
 
  // Fill BEV. 
  BEV << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,5,1) = u; 
 
} 
 
void CartFieldBinOpDivide2x2vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
 
  double As[6]; 
  double Bs[15]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = 0.0; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = 0.0; 
    Bs[6] = 0.0; 
    Bs[7] = 0.0; 
    Bs[8] = 0.0; 
    Bs[9] = 0.0; 
    Bs[10] = B[10]; 
    Bs[11] = 0.0; 
    Bs[12] = 0.0; 
    Bs[13] = B[13]; 
    Bs[14] = B[14]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
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
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(15,15); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(15);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(15);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*As[0]; 
  AEM(0,1) = 0.5*As[1]; 
  AEM(0,2) = 0.5*As[2]; 
  AEM(0,5) = 0.5*As[3]; 
  AEM(0,6) = 0.5*As[1]; 
  AEM(0,7) = 0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(0,8) = 0.5*As[3]; 
  AEM(0,11) = 0.5*As[2]; 
  AEM(0,12) = 0.5*As[2]; 
  AEM(0,13) = 0.5*As[3]; 
  AEM(0,14) = 0.4472135954999579*As[5]+0.5*As[0]; 
  AEM(1,2) = 0.5*As[1]; 
  AEM(1,6) = 0.5*As[0]; 
  AEM(1,13) = 0.5*As[0]; 
  AEM(2,0) = 0.5*As[3]; 
  AEM(2,1) = 0.5*As[2]; 
  AEM(2,2) = 0.5*As[1]; 
  AEM(2,5) = 0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(2,9) = 0.5*As[1]; 
  AEM(3,0) = 0.5*As[2]; 
  AEM(3,7) = 0.5*As[1]; 
  AEM(3,13) = 0.5*As[2]; 
  AEM(4,6) = 0.5*As[4]; 
  AEM(4,7) = 0.4472135954999579*As[1]; 
  AEM(4,11) = 0.4472135954999579*As[3]; 
  AEM(4,12) = 0.5*As[5]; 
  AEM(4,14) = 0.4472135954999579*As[2]; 
  AEM(5,2) = 0.4472135954999579*As[3]; 
 
  // Fill BEV. 
  BEV << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9],Bs[10],Bs[11],Bs[12],Bs[13],Bs[14]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,15,1) = u; 
 
} 
 
