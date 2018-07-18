#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x1vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[3]; 
  tmp[0] = 0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[0]*B[2]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<3; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x1vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[6]; 
  tmp[0] = 0.7071067811865475*A[2]*B[4]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6324555320336759*A[1]*B[4]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[1]*B[3]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.6324555320336759*A[2]*B[3]+0.7071067811865475*A[0]*B[3]+0.7071067811865475*A[1]*B[2]; 
  tmp[4] = 0.4517539514526256*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[5] = 0.7071067811865475*A[0]*B[5]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<6; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x1vMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[10]; 
  tmp[0] = 0.7071067811865475*A[3]*B[8]+0.7071067811865475*A[2]*B[4]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6210590034081186*A[2]*B[8]+0.6210590034081186*A[3]*B[4]+0.6324555320336759*A[1]*B[4]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[2]*B[6]+0.7071067811865475*A[1]*B[3]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.6210590034081187*A[3]*B[6]+0.632455532033676*A[1]*B[6]+0.6324555320336759*A[2]*B[3]+0.7071067811865475*A[0]*B[3]+0.7071067811865475*A[1]*B[2]; 
  tmp[4] = 0.421637021355784*A[3]*B[8]+0.6210590034081186*A[1]*B[8]+0.4517539514526256*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.6210590034081186*B[1]*A[3]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[5] = 0.7071067811865475*A[1]*B[7]+0.7071067811865475*A[0]*B[5]; 
  tmp[6] = 0.4517539514526256*A[2]*B[6]+0.7071067811865475*A[0]*B[6]+0.6210590034081187*A[3]*B[3]+0.632455532033676*A[1]*B[3]+0.7071067811865475*A[2]*B[2]; 
  tmp[7] = 0.6324555320336759*A[2]*B[7]+0.7071067811865475*A[0]*B[7]+0.7071067811865475*A[1]*B[5]; 
  tmp[8] = 0.421637021355784*A[2]*B[8]+0.7071067811865475*A[0]*B[8]+0.421637021355784*A[3]*B[4]+0.6210590034081186*A[1]*B[4]+0.7071067811865475*B[0]*A[3]+0.6210590034081186*B[1]*A[2]; 
  tmp[9] = 0.7071067811865475*A[0]*B[9]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<10; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide1x1vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(3,3); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(3);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(3);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.7071067811865475*A[0]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,3,1) = u; 
 
} 
 
void CartFieldBinOpDivide1x1vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(6,6); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(6);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(6);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(0,4) = 0.7071067811865475*A[2]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(1,4) = 0.6324555320336759*A[1]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
  AEM(2,3) = 0.7071067811865475*A[1]; 
  AEM(3,2) = 0.7071067811865475*A[1]; 
  AEM(3,3) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(4,0) = 0.7071067811865475*A[2]; 
  AEM(4,1) = 0.6324555320336759*A[1]; 
  AEM(4,4) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(5,5) = 0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,6,1) = u; 
 
} 
 
void CartFieldBinOpDivide1x1vMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(10,10); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(10);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(10);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(0,4) = 0.7071067811865475*A[2]; 
  AEM(0,8) = 0.7071067811865475*A[3]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(1,4) = 0.6210590034081186*A[3]+0.6324555320336759*A[1]; 
  AEM(1,8) = 0.6210590034081186*A[2]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
  AEM(2,3) = 0.7071067811865475*A[1]; 
  AEM(2,6) = 0.7071067811865475*A[2]; 
  AEM(3,2) = 0.7071067811865475*A[1]; 
  AEM(3,3) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(3,6) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(4,0) = 0.7071067811865475*A[2]; 
  AEM(4,1) = 0.6210590034081186*A[3]+0.6324555320336759*A[1]; 
  AEM(4,4) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(4,8) = 0.421637021355784*A[3]+0.6210590034081186*A[1]; 
  AEM(5,5) = 0.7071067811865475*A[0]; 
  AEM(5,7) = 0.7071067811865475*A[1]; 
  AEM(6,2) = 0.7071067811865475*A[2]; 
  AEM(6,3) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(6,6) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(7,5) = 0.7071067811865475*A[1]; 
  AEM(7,7) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(8,0) = 0.7071067811865475*A[3]; 
  AEM(8,1) = 0.6210590034081186*A[2]; 
  AEM(8,4) = 0.421637021355784*A[3]+0.6210590034081186*A[1]; 
  AEM(8,8) = 0.421637021355784*A[2]+0.7071067811865475*A[0]; 
  AEM(9,9) = 0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,10,1) = u; 
 
} 
 
