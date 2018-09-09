#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2x3vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[6]; 
  tmp[0] = 0.5*A[2]*B[2]+0.5*A[1]*B[1]+0.5*A[0]*B[0]; 
  tmp[1] = 0.5*A[0]*B[1]+0.5*B[0]*A[1]; 
  tmp[2] = 0.5*A[0]*B[2]+0.5*B[0]*A[2]; 
  tmp[3] = 0.5*A[0]*B[3]; 
  tmp[4] = 0.5*A[0]*B[4]; 
  tmp[5] = 0.5*A[0]*B[5]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<6; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply2x3vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[21]; 
  tmp[0] = 0.5*A[5]*B[17]+0.5*A[4]*B[16]+0.5*A[3]*B[6]+0.5*A[2]*B[2]+0.5*A[1]*B[1]+0.5*A[0]*B[0]; 
  tmp[1] = 0.4472135954999579*A[1]*B[16]+0.5*A[2]*B[6]+0.4472135954999579*B[1]*A[4]+0.5*B[2]*A[3]+0.5*A[0]*B[1]+0.5*B[0]*A[1]; 
  tmp[2] = 0.4472135954999579*A[2]*B[17]+0.5*A[1]*B[6]+0.4472135954999579*B[2]*A[5]+0.5*B[1]*A[3]+0.5*A[0]*B[2]+0.5*B[0]*A[2]; 
  tmp[3] = 0.5*A[2]*B[8]+0.5*A[1]*B[7]+0.5*A[0]*B[3]; 
  tmp[4] = 0.5*A[2]*B[10]+0.5*A[1]*B[9]+0.5*A[0]*B[4]; 
  tmp[5] = 0.5*A[2]*B[13]+0.5*A[1]*B[12]+0.5*A[0]*B[5]; 
  tmp[6] = 0.4472135954999579*A[3]*B[17]+0.4472135954999579*A[3]*B[16]+0.4472135954999579*A[5]*B[6]+0.4472135954999579*A[4]*B[6]+0.5*A[0]*B[6]+0.5*B[0]*A[3]+0.5*A[1]*B[2]+0.5*B[1]*A[2]; 
  tmp[7] = 0.5*A[3]*B[8]+0.4472135954999579*A[4]*B[7]+0.5*A[0]*B[7]+0.5*A[1]*B[3]; 
  tmp[8] = 0.4472135954999579*A[5]*B[8]+0.5*A[0]*B[8]+0.5*A[3]*B[7]+0.5*A[2]*B[3]; 
  tmp[9] = 0.5*A[3]*B[10]+0.4472135954999579*A[4]*B[9]+0.5*A[0]*B[9]+0.5*A[1]*B[4]; 
  tmp[10] = 0.4472135954999579*A[5]*B[10]+0.5*A[0]*B[10]+0.5*A[3]*B[9]+0.5*A[2]*B[4]; 
  tmp[11] = 0.5*A[0]*B[11]; 
  tmp[12] = 0.5*A[3]*B[13]+0.4472135954999579*A[4]*B[12]+0.5*A[0]*B[12]+0.5*A[1]*B[5]; 
  tmp[13] = 0.4472135954999579*A[5]*B[13]+0.5*A[0]*B[13]+0.5*A[3]*B[12]+0.5*A[2]*B[5]; 
  tmp[14] = 0.5*A[0]*B[14]; 
  tmp[15] = 0.5*A[0]*B[15]; 
  tmp[16] = 0.31943828249997*A[4]*B[16]+0.5*A[0]*B[16]+0.4472135954999579*A[3]*B[6]+0.5*B[0]*A[4]+0.4472135954999579*A[1]*B[1]; 
  tmp[17] = 0.31943828249997*A[5]*B[17]+0.5*A[0]*B[17]+0.4472135954999579*A[3]*B[6]+0.5*B[0]*A[5]+0.4472135954999579*A[2]*B[2]; 
  tmp[18] = 0.5*A[0]*B[18]; 
  tmp[19] = 0.5*A[0]*B[19]; 
  tmp[20] = 0.5*A[0]*B[20]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<21; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide2x3vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(6,6); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(6);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(6);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*As[0]; 
  AEM(0,1) = 0.5*As[1]; 
  AEM(0,2) = 0.5*As[2]; 
  AEM(0,3) = 0.5*As[1]; 
  AEM(0,4) = 0.5*As[0]; 
  AEM(1,0) = 0.5*As[2]; 
  AEM(1,2) = 0.5*As[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,6,1) = u; 
 
} 
 
void CartFieldBinOpDivide2x3vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(21,21); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(21);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(21);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*As[0]; 
  AEM(0,1) = 0.5*As[1]; 
  AEM(0,2) = 0.5*As[2]; 
  AEM(0,6) = 0.5*As[1]; 
  AEM(0,7) = 0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(0,8) = 0.5*As[3]; 
  AEM(0,12) = 0.5*As[2]; 
  AEM(0,13) = 0.5*As[3]; 
  AEM(0,14) = 0.4472135954999579*As[5]+0.5*As[0]; 
  AEM(1,0) = 0.5*As[0]; 
  AEM(1,7) = 0.5*As[0]; 
  AEM(1,14) = 0.5*As[0]; 
  AEM(1,15) = 0.5*As[3]; 
  AEM(1,16) = 0.5*As[2]; 
  AEM(1,17) = 0.5*As[1]; 
  AEM(2,3) = 0.5*As[1]; 
  AEM(2,9) = 0.5*As[2]; 
  AEM(2,16) = 0.5*As[1]; 
  AEM(3,1) = 0.5*As[2]; 
  AEM(3,14) = 0.5*As[1]; 
  AEM(3,20) = 0.5*As[2]; 
  AEM(4,12) = 0.5*As[4]; 
  AEM(4,13) = 0.4472135954999579*As[1]; 
  AEM(4,18) = 0.5*As[5]; 
  AEM(4,20) = 0.4472135954999579*As[2]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10],B[11],B[12],B[13],B[14],B[15],B[16],B[17],B[18],B[19],B[20]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,21,1) = u; 
 
} 
 
