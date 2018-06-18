#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x3vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[5]; 
  tmp[0] = 0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.7071067811865475*A[0]*B[4]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<5; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x3vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[15]; 
  tmp[0] = 0.7071067811865475*A[2]*B[11]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6324555320336759*A[1]*B[11]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[1]*B[5]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[1]*B[6]+0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.7071067811865475*A[1]*B[8]+0.7071067811865475*A[0]*B[4]; 
  tmp[5] = 0.6324555320336759*A[2]*B[5]+0.7071067811865475*A[0]*B[5]+0.7071067811865475*A[1]*B[2]; 
  tmp[6] = 0.6324555320336759*A[2]*B[6]+0.7071067811865475*A[0]*B[6]+0.7071067811865475*A[1]*B[3]; 
  tmp[7] = 0.7071067811865475*A[0]*B[7]; 
  tmp[8] = 0.6324555320336759*A[2]*B[8]+0.7071067811865475*A[0]*B[8]+0.7071067811865475*A[1]*B[4]; 
  tmp[9] = 0.7071067811865475*A[0]*B[9]; 
  tmp[10] = 0.7071067811865475*A[0]*B[10]; 
  tmp[11] = 0.4517539514526256*A[2]*B[11]+0.7071067811865475*A[0]*B[11]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[12] = 0.7071067811865475*A[0]*B[12]; 
  tmp[13] = 0.7071067811865475*A[0]*B[13]; 
  tmp[14] = 0.7071067811865475*A[0]*B[14]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<15; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x3vMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[35]; 
  tmp[0] = 0.7071067811865475*A[3]*B[31]+0.7071067811865475*A[2]*B[11]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6210590034081186*A[2]*B[31]+0.6210590034081186*A[3]*B[11]+0.6324555320336759*A[1]*B[11]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[2]*B[19]+0.7071067811865475*A[1]*B[5]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[2]*B[21]+0.7071067811865475*A[1]*B[6]+0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.7071067811865475*A[2]*B[25]+0.7071067811865475*A[1]*B[8]+0.7071067811865475*A[0]*B[4]; 
  tmp[5] = 0.6210590034081187*A[3]*B[19]+0.632455532033676*A[1]*B[19]+0.6324555320336759*A[2]*B[5]+0.7071067811865475*A[0]*B[5]+0.7071067811865475*A[1]*B[2]; 
  tmp[6] = 0.6210590034081187*A[3]*B[21]+0.632455532033676*A[1]*B[21]+0.6324555320336759*A[2]*B[6]+0.7071067811865475*A[0]*B[6]+0.7071067811865475*A[1]*B[3]; 
  tmp[7] = 0.7071067811865475*A[1]*B[15]+0.7071067811865475*A[0]*B[7]; 
  tmp[8] = 0.6210590034081187*A[3]*B[25]+0.632455532033676*A[1]*B[25]+0.6324555320336759*A[2]*B[8]+0.7071067811865475*A[0]*B[8]+0.7071067811865475*A[1]*B[4]; 
  tmp[9] = 0.7071067811865475*A[1]*B[16]+0.7071067811865475*A[0]*B[9]; 
  tmp[10] = 0.7071067811865475*A[1]*B[17]+0.7071067811865475*A[0]*B[10]; 
  tmp[11] = 0.421637021355784*A[3]*B[31]+0.6210590034081186*A[1]*B[31]+0.4517539514526256*A[2]*B[11]+0.7071067811865475*A[0]*B[11]+0.6210590034081186*B[1]*A[3]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[12] = 0.7071067811865475*A[1]*B[20]+0.7071067811865475*A[0]*B[12]; 
  tmp[13] = 0.7071067811865475*A[1]*B[23]+0.7071067811865475*A[0]*B[13]; 
  tmp[14] = 0.7071067811865475*A[1]*B[28]+0.7071067811865475*A[0]*B[14]; 
  tmp[15] = 0.6324555320336759*A[2]*B[15]+0.7071067811865475*A[0]*B[15]+0.7071067811865475*A[1]*B[7]; 
  tmp[16] = 0.6324555320336759*A[2]*B[16]+0.7071067811865475*A[0]*B[16]+0.7071067811865475*A[1]*B[9]; 
  tmp[17] = 0.6324555320336759*A[2]*B[17]+0.7071067811865475*A[0]*B[17]+0.7071067811865475*A[1]*B[10]; 
  tmp[18] = 0.7071067811865475*A[0]*B[18]; 
  tmp[19] = 0.4517539514526256*A[2]*B[19]+0.7071067811865475*A[0]*B[19]+0.6210590034081187*A[3]*B[5]+0.632455532033676*A[1]*B[5]+0.7071067811865475*A[2]*B[2]; 
  tmp[20] = 0.6324555320336759*A[2]*B[20]+0.7071067811865475*A[0]*B[20]+0.7071067811865475*A[1]*B[12]; 
  tmp[21] = 0.4517539514526256*A[2]*B[21]+0.7071067811865475*A[0]*B[21]+0.6210590034081187*A[3]*B[6]+0.632455532033676*A[1]*B[6]+0.7071067811865475*A[2]*B[3]; 
  tmp[22] = 0.7071067811865475*A[0]*B[22]; 
  tmp[23] = 0.6324555320336759*A[2]*B[23]+0.7071067811865475*A[0]*B[23]+0.7071067811865475*A[1]*B[13]; 
  tmp[24] = 0.7071067811865475*A[0]*B[24]; 
  tmp[25] = 0.4517539514526256*A[2]*B[25]+0.7071067811865475*A[0]*B[25]+0.6210590034081187*A[3]*B[8]+0.632455532033676*A[1]*B[8]+0.7071067811865475*A[2]*B[4]; 
  tmp[26] = 0.7071067811865475*A[0]*B[26]; 
  tmp[27] = 0.7071067811865475*A[0]*B[27]; 
  tmp[28] = 0.6324555320336759*A[2]*B[28]+0.7071067811865475*A[0]*B[28]+0.7071067811865475*A[1]*B[14]; 
  tmp[29] = 0.7071067811865475*A[0]*B[29]; 
  tmp[30] = 0.7071067811865475*A[0]*B[30]; 
  tmp[31] = 0.421637021355784*A[2]*B[31]+0.7071067811865475*A[0]*B[31]+0.421637021355784*A[3]*B[11]+0.6210590034081186*A[1]*B[11]+0.7071067811865475*B[0]*A[3]+0.6210590034081186*B[1]*A[2]; 
  tmp[32] = 0.7071067811865475*A[0]*B[32]; 
  tmp[33] = 0.7071067811865475*A[0]*B[33]; 
  tmp[34] = 0.7071067811865475*A[0]*B[34]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<35; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide1x3vMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(5,5); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(5); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(5); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.7071067811865475*A[0]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
  AEM(3,3) = 0.7071067811865475*A[0]; 
  AEM(4,4) = 0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,5,1) = u; 
 
} 
 
void CartFieldBinOpDivide1x3vMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(15,15); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(15); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(15); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(0,11) = 0.7071067811865475*A[2]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(1,11) = 0.6324555320336759*A[1]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
  AEM(2,5) = 0.7071067811865475*A[1]; 
  AEM(3,3) = 0.7071067811865475*A[0]; 
  AEM(3,6) = 0.7071067811865475*A[1]; 
  AEM(4,4) = 0.7071067811865475*A[0]; 
  AEM(4,8) = 0.7071067811865475*A[1]; 
  AEM(5,2) = 0.7071067811865475*A[1]; 
  AEM(5,5) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(6,3) = 0.7071067811865475*A[1]; 
  AEM(6,6) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(7,7) = 0.7071067811865475*A[0]; 
  AEM(8,4) = 0.7071067811865475*A[1]; 
  AEM(8,8) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(9,9) = 0.7071067811865475*A[0]; 
  AEM(10,10) = 0.7071067811865475*A[0]; 
  AEM(11,0) = 0.7071067811865475*A[2]; 
  AEM(11,1) = 0.6324555320336759*A[1]; 
  AEM(11,11) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(12,12) = 0.7071067811865475*A[0]; 
  AEM(13,13) = 0.7071067811865475*A[0]; 
  AEM(14,14) = 0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10],B[11],B[12],B[13],B[14]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,15,1) = u; 
 
} 
 
void CartFieldBinOpDivide1x3vMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(35,35); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(35); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(35); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(0,11) = 0.7071067811865475*A[2]; 
  AEM(0,31) = 0.7071067811865475*A[3]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(1,11) = 0.6210590034081186*A[3]+0.6324555320336759*A[1]; 
  AEM(1,31) = 0.6210590034081186*A[2]; 
  AEM(2,2) = 0.7071067811865475*A[0]; 
  AEM(2,5) = 0.7071067811865475*A[1]; 
  AEM(2,19) = 0.7071067811865475*A[2]; 
  AEM(3,3) = 0.7071067811865475*A[0]; 
  AEM(3,6) = 0.7071067811865475*A[1]; 
  AEM(3,21) = 0.7071067811865475*A[2]; 
  AEM(4,4) = 0.7071067811865475*A[0]; 
  AEM(4,8) = 0.7071067811865475*A[1]; 
  AEM(4,25) = 0.7071067811865475*A[2]; 
  AEM(5,2) = 0.7071067811865475*A[1]; 
  AEM(5,5) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(5,19) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(6,3) = 0.7071067811865475*A[1]; 
  AEM(6,6) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(6,21) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(7,7) = 0.7071067811865475*A[0]; 
  AEM(7,15) = 0.7071067811865475*A[1]; 
  AEM(8,4) = 0.7071067811865475*A[1]; 
  AEM(8,8) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(8,25) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(9,9) = 0.7071067811865475*A[0]; 
  AEM(9,16) = 0.7071067811865475*A[1]; 
  AEM(10,10) = 0.7071067811865475*A[0]; 
  AEM(10,17) = 0.7071067811865475*A[1]; 
  AEM(11,0) = 0.7071067811865475*A[2]; 
  AEM(11,1) = 0.6210590034081186*A[3]+0.6324555320336759*A[1]; 
  AEM(11,11) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(11,31) = 0.421637021355784*A[3]+0.6210590034081186*A[1]; 
  AEM(12,12) = 0.7071067811865475*A[0]; 
  AEM(12,20) = 0.7071067811865475*A[1]; 
  AEM(13,13) = 0.7071067811865475*A[0]; 
  AEM(13,23) = 0.7071067811865475*A[1]; 
  AEM(14,14) = 0.7071067811865475*A[0]; 
  AEM(14,28) = 0.7071067811865475*A[1]; 
  AEM(15,7) = 0.7071067811865475*A[1]; 
  AEM(15,15) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(16,9) = 0.7071067811865475*A[1]; 
  AEM(16,16) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(17,10) = 0.7071067811865475*A[1]; 
  AEM(17,17) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(18,18) = 0.7071067811865475*A[0]; 
  AEM(19,2) = 0.7071067811865475*A[2]; 
  AEM(19,5) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(19,19) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(20,12) = 0.7071067811865475*A[1]; 
  AEM(20,20) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(21,3) = 0.7071067811865475*A[2]; 
  AEM(21,6) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(21,21) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(22,22) = 0.7071067811865475*A[0]; 
  AEM(23,13) = 0.7071067811865475*A[1]; 
  AEM(23,23) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(24,24) = 0.7071067811865475*A[0]; 
  AEM(25,4) = 0.7071067811865475*A[2]; 
  AEM(25,8) = 0.6210590034081187*A[3]+0.632455532033676*A[1]; 
  AEM(25,25) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
  AEM(26,26) = 0.7071067811865475*A[0]; 
  AEM(27,27) = 0.7071067811865475*A[0]; 
  AEM(28,14) = 0.7071067811865475*A[1]; 
  AEM(28,28) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(29,29) = 0.7071067811865475*A[0]; 
  AEM(30,30) = 0.7071067811865475*A[0]; 
  AEM(31,0) = 0.7071067811865475*A[3]; 
  AEM(31,1) = 0.6210590034081186*A[2]; 
  AEM(31,11) = 0.421637021355784*A[3]+0.6210590034081186*A[1]; 
  AEM(31,31) = 0.421637021355784*A[2]+0.7071067811865475*A[0]; 
  AEM(32,32) = 0.7071067811865475*A[0]; 
  AEM(33,33) = 0.7071067811865475*A[0]; 
  AEM(34,34) = 0.7071067811865475*A[0]; 
 
  // Fill BEV. 
  BEV << B[0],B[1],B[2],B[3],B[4],B[5],B[6],B[7],B[8],B[9],B[10],B[11],B[12],B[13],B[14],B[15],B[16],B[17],B[18],B[19],B[20],B[21],B[22],B[23],B[24],B[25],B[26],B[27],B[28],B[29],B[30],B[31],B[32],B[33],B[34]; 
 
  // Solve the system of equations. 
  u = AEM.colPivHouseholderQr().solve(BEV); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,35,1) = u; 
 
} 
 
