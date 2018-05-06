#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[2*Ncomp]; 
 
  for (short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 2*vd; 
    short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[b0] = 0.7071067811865475*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0]*B[b0]; 
    tmp[b0+1] = 0.7071067811865475*A[a0]*B[b0+1]+0.7071067811865475*A[a0+1]*B[b0]; 
  } 
 
  // This tmp allows for in-place multiplication. 
  for (short int i=0; i<2*Ncomp; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[3*Ncomp]; 
 
  for (short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 3*vd; 
    short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[b0] = 0.7071067811865475*A[a0+2]*B[b0+2]+0.7071067811865475*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0]*B[b0]; 
    tmp[b0+1] = 0.6324555320336759*A[a0+1]*B[b0+2]+0.6324555320336759*A[a0+2]*B[b0+1]+0.7071067811865475*A[a0]*B[b0+1]+0.7071067811865475*A[a0+1]*B[b0]; 
    tmp[b0+2] = 0.4517539514526256*A[a0+2]*B[b0+2]+0.7071067811865475*A[a0]*B[b0+2]+0.6324555320336759*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0+2]*B[b0]; 
  } 
 
  // This tmp allows for in-place multiplication. 
  for (short int i=0; i<3*Ncomp; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide1xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(2,2); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(2); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(2); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.7071067811865475*A[0]; 
 
  for(short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 2*vd; 
    // Fill BEV. 
    BEV << B[b0],B[b0+1]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*2,2,1) = u; 
 
  } 
} 
 
void CartFieldBinOpDivide1xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(3,3); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(3); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(3); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*A[0]; 
  AEM(0,1) = 0.7071067811865475*A[1]; 
  AEM(0,2) = 0.7071067811865475*A[2]; 
  AEM(1,0) = 0.7071067811865475*A[1]; 
  AEM(1,1) = 0.6324555320336759*A[2]+0.7071067811865475*A[0]; 
  AEM(1,2) = 0.6324555320336759*A[1]; 
  AEM(2,0) = 0.7071067811865475*A[2]; 
  AEM(2,1) = 0.6324555320336759*A[1]; 
  AEM(2,2) = 0.4517539514526256*A[2]+0.7071067811865475*A[0]; 
 
  for(short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 3*vd; 
    // Fill BEV. 
    BEV << B[b0],B[b0+1],B[b0+2]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*3,3,1) = u; 
 
  } 
} 
 
