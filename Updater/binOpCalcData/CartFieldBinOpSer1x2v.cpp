#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[8]; 
 
  tmp[0] = 0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[1]*B[4]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[1]*B[5]+0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.7071067811865475*A[0]*B[4]+0.7071067811865475*A[1]*B[2]; 
  tmp[5] = 0.7071067811865475*A[0]*B[5]+0.7071067811865475*A[1]*B[3]; 
  tmp[6] = 0.7071067811865475*A[1]*B[7]+0.7071067811865475*A[0]*B[6]; 
  tmp[7] = 0.7071067811865475*A[0]*B[7]+0.7071067811865475*A[1]*B[6]; 
 
  // This tmp allows for in-place multiplication. 
  for (short int i=0; i<8; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[20]; 
 
  tmp[0] = 0.7071067811865475*A[2]*B[7]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6324555320336759*A[1]*B[7]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[2]*B[11]+0.7071067811865475*A[1]*B[4]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[2]*B[13]+0.7071067811865475*A[1]*B[5]+0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.632455532033676*A[1]*B[11]+0.6324555320336759*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.7071067811865475*A[1]*B[2]; 
  tmp[5] = 0.632455532033676*A[1]*B[13]+0.6324555320336759*A[2]*B[5]+0.7071067811865475*A[0]*B[5]+0.7071067811865475*A[1]*B[3]; 
  tmp[6] = 0.7071067811865475*A[2]*B[17]+0.7071067811865475*A[1]*B[10]+0.7071067811865475*A[0]*B[6]; 
  tmp[7] = 0.4517539514526256*A[2]*B[7]+0.7071067811865475*A[0]*B[7]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[8] = 0.7071067811865475*A[1]*B[12]+0.7071067811865475*A[0]*B[8]; 
  tmp[9] = 0.7071067811865475*A[1]*B[15]+0.7071067811865475*A[0]*B[9]; 
  tmp[10] = 0.6324555320336759*A[1]*B[17]+0.6324555320336759*A[2]*B[10]+0.7071067811865475*A[0]*B[10]+0.7071067811865475*A[1]*B[6]; 
  tmp[11] = 0.4517539514526256*A[2]*B[11]+0.7071067811865475*A[0]*B[11]+0.632455532033676*A[1]*B[4]+0.7071067811865475*A[2]*B[2]; 
  tmp[12] = 0.6324555320336759*A[2]*B[12]+0.7071067811865475*A[0]*B[12]+0.7071067811865475*A[1]*B[8]; 
  tmp[13] = 0.4517539514526256*A[2]*B[13]+0.7071067811865475*A[0]*B[13]+0.632455532033676*A[1]*B[5]+0.7071067811865475*A[2]*B[3]; 
  tmp[14] = 0.7071067811865475*A[1]*B[18]+0.7071067811865475*A[0]*B[14]; 
  tmp[15] = 0.6324555320336759*A[2]*B[15]+0.7071067811865475*A[0]*B[15]+0.7071067811865475*A[1]*B[9]; 
  tmp[16] = 0.7071067811865475*A[1]*B[19]+0.7071067811865475*A[0]*B[16]; 
  tmp[17] = 0.4517539514526256*A[2]*B[17]+0.7071067811865475*A[0]*B[17]+0.6324555320336759*A[1]*B[10]+0.7071067811865475*A[2]*B[6]; 
  tmp[18] = 0.6324555320336759*A[2]*B[18]+0.7071067811865475*A[0]*B[18]+0.7071067811865475*A[1]*B[14]; 
  tmp[19] = 0.6324555320336759*A[2]*B[19]+0.7071067811865475*A[0]*B[19]+0.7071067811865475*A[1]*B[16]; 
 
  // This tmp allows for in-place multiplication. 
  for (short int i=0; i<20; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide1x2vSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(8,8); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(8); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(8); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.7071067811865475*A[0]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
  AEM(2,4) = 0.7071067811865475*A[1]; 
  AEM(3,3) = 0.7071067811865475*A[0]; 
  AEM(3,5) = 0.7071067811865475*A[1]; 
  AEM(4,2) = 0.7071067811865475*A[1]; 
  AEM(4,4) = 0.7071067811865475*A[0]; 
  AEM(5,3) = 0.7071067811865475*A[1]; 
  AEM(5,5) = 0.7071067811865475*A[0]; 
  AEM(6,6) = 0.7071067811865475*A[0]; 
  AEM(6,7) = 0.7071067811865475*A[1]; 
  AEM(7,6) = 0.7071067811865475*A[1]; 
  AEM(7,7) = 0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,8,1) = u; 
 
} 
 
void CartFieldBinOpDivide1x2vSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(20,20); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(20); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(20); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(0,7) = 0.7071067811865475*A[2]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(1,7) = 0.6324555320336759*A[1]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
  AEM(2,4) = 0.7071067811865475*A[1]; 
  AEM(2,11) = 0.7071067811865475*A[2]; 
  AEM(3,3) = 0.7071067811865475*A[0]; 
  AEM(3,5) = 0.7071067811865475*A[1]; 
  AEM(3,13) = 0.7071067811865475*A[2]; 
  AEM(4,2) = 0.7071067811865475*A[1]; 
  AEM(4,4) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(4,11) = 0.632455532033676*A[1]; 
  AEM(5,3) = 0.7071067811865475*A[1]; 
  AEM(5,5) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(5,13) = 0.632455532033676*A[1]; 
  AEM(6,6) = 0.7071067811865475*A[0]; 
  AEM(6,10) = 0.7071067811865475*A[1]; 
  AEM(6,17) = 0.7071067811865475*A[2]; 
  AEM(7,0) = 0.7071067811865475*A[2]; 
  AEM(7,1) = 0.6324555320336759*A[1]; 
  AEM(7,7) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(8,8) = 0.7071067811865475*A[0]; 
  AEM(8,12) = 0.7071067811865475*A[1]; 
  AEM(9,9) = 0.7071067811865475*A[0]; 
  AEM(9,15) = 0.7071067811865475*A[1]; 
  AEM(10,6) = 0.7071067811865475*A[1]; 
  AEM(10,10) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(10,17) = 0.6324555320336759*A[1]; 
  AEM(11,2) = 0.7071067811865475*A[2]; 
  AEM(11,4) = 0.632455532033676*A[1]; 
  AEM(11,11) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(12,8) = 0.7071067811865475*A[1]; 
  AEM(12,12) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(13,3) = 0.7071067811865475*A[2]; 
  AEM(13,5) = 0.632455532033676*A[1]; 
  AEM(13,13) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(14,14) = 0.7071067811865475*A[0]; 
  AEM(14,18) = 0.7071067811865475*A[1]; 
  AEM(15,9) = 0.7071067811865475*A[1]; 
  AEM(15,15) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(16,16) = 0.7071067811865475*A[0]; 
  AEM(16,19) = 0.7071067811865475*A[1]; 
  AEM(17,6) = 0.7071067811865475*A[2]; 
  AEM(17,10) = 0.6324555320336759*A[1]; 
  AEM(17,17) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(18,14) = 0.7071067811865475*A[1]; 
  AEM(18,18) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(19,16) = 0.7071067811865475*A[1]; 
  AEM(19,19) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10],B[11],B[12],B[13],B[14],B[15],B[16],B[17],B[18],B[19]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,20,1) = u; 
 
} 
 
