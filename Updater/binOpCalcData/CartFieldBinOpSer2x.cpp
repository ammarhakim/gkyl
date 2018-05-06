#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[4*Ncomp]; 
 
  for (short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 4*vd; 
    short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[b0] = 0.5*A[a0+3]*B[b0+3]+0.5*A[a0+2]*B[b0+2]+0.5*A[a0+1]*B[b0+1]+0.5*A[a0]*B[b0]; 
    tmp[b0+1] = 0.5*A[a0+2]*B[b0+3]+0.5*A[a0+3]*B[b0+2]+0.5*A[a0]*B[b0+1]+0.5*A[a0+1]*B[b0]; 
    tmp[b0+2] = 0.5*A[a0+1]*B[b0+3]+0.5*A[a0]*B[b0+2]+0.5*A[a0+3]*B[b0+1]+0.5*A[a0+2]*B[b0]; 
    tmp[b0+3] = 0.5*A[a0]*B[b0+3]+0.5*A[a0+1]*B[b0+2]+0.5*A[a0+2]*B[b0+1]+0.5*A[a0+3]*B[b0]; 
  } 
 
  // This tmp allows for in-place multiplication. 
  for (short int i=0; i<4*Ncomp; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply2xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[8*Ncomp]; 
 
  for (short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 8*vd; 
    short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[b0] = 0.5*A[a0+7]*B[b0+7]+0.5*A[a0+6]*B[b0+6]+0.5*A[a0+5]*B[b0+5]+0.5*A[a0+4]*B[b0+4]+0.5*A[a0+3]*B[b0+3]+0.5*A[a0+2]*B[b0+2]+0.5*A[a0+1]*B[b0+1]+0.5*A[a0]*B[b0]; 
    tmp[b0+1] = 0.5000000000000001*A[a0+5]*B[b0+7]+0.447213595499958*A[a0+3]*B[b0+6]+0.5000000000000001*A[a0+7]*B[b0+5]+0.4472135954999579*A[a0+1]*B[b0+4]+0.447213595499958*A[a0+6]*B[b0+3]+0.5*A[a0+2]*B[b0+3]+0.5*A[a0+3]*B[b0+2]+0.4472135954999579*A[a0+4]*B[b0+1]+0.5*A[a0]*B[b0+1]+0.5*A[a0+1]*B[b0]; 
    tmp[b0+2] = 0.447213595499958*A[a0+3]*B[b0+7]+0.5000000000000001*A[a0+4]*B[b0+6]+0.4472135954999579*A[a0+2]*B[b0+5]+0.5000000000000001*A[a0+6]*B[b0+4]+0.447213595499958*A[a0+7]*B[b0+3]+0.5*A[a0+1]*B[b0+3]+0.4472135954999579*A[a0+5]*B[b0+2]+0.5*A[a0]*B[b0+2]+0.5*A[a0+3]*B[b0+1]+0.5*A[a0+2]*B[b0]; 
    tmp[b0+3] = 0.4*A[a0+6]*B[b0+7]+0.447213595499958*A[a0+2]*B[b0+7]+0.4*A[a0+7]*B[b0+6]+0.447213595499958*A[a0+1]*B[b0+6]+0.4472135954999579*A[a0+3]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+4]+0.4472135954999579*A[a0+5]*B[b0+3]+0.4472135954999579*A[a0+4]*B[b0+3]+0.5*A[a0]*B[b0+3]+0.447213595499958*A[a0+7]*B[b0+2]+0.5*A[a0+1]*B[b0+2]+0.447213595499958*A[a0+6]*B[b0+1]+0.5*A[a0+2]*B[b0+1]+0.5*A[a0+3]*B[b0]; 
    tmp[b0+4] = 0.4472135954999579*A[a0+7]*B[b0+7]+0.31943828249997*A[a0+6]*B[b0+6]+0.5000000000000001*A[a0+2]*B[b0+6]+0.31943828249997*A[a0+4]*B[b0+4]+0.5*A[a0]*B[b0+4]+0.4472135954999579*A[a0+3]*B[b0+3]+0.5000000000000001*A[a0+6]*B[b0+2]+0.4472135954999579*A[a0+1]*B[b0+1]+0.5*A[a0+4]*B[b0]; 
    tmp[b0+5] = 0.31943828249997*A[a0+7]*B[b0+7]+0.5000000000000001*A[a0+1]*B[b0+7]+0.4472135954999579*A[a0+6]*B[b0+6]+0.31943828249997*A[a0+5]*B[b0+5]+0.5*A[a0]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+3]+0.4472135954999579*A[a0+2]*B[b0+2]+0.5000000000000001*A[a0+7]*B[b0+1]+0.5*A[a0+5]*B[b0]; 
    tmp[b0+6] = 0.4*A[a0+3]*B[b0+7]+0.4472135954999579*A[a0+5]*B[b0+6]+0.31943828249997*A[a0+4]*B[b0+6]+0.5*A[a0]*B[b0+6]+0.4472135954999579*A[a0+6]*B[b0+5]+0.31943828249997*A[a0+6]*B[b0+4]+0.5000000000000001*A[a0+2]*B[b0+4]+0.4*A[a0+7]*B[b0+3]+0.447213595499958*A[a0+1]*B[b0+3]+0.5000000000000001*A[a0+4]*B[b0+2]+0.447213595499958*A[a0+3]*B[b0+1]+0.5*A[a0+6]*B[b0]; 
    tmp[b0+7] = 0.31943828249997*A[a0+5]*B[b0+7]+0.4472135954999579*A[a0+4]*B[b0+7]+0.5*A[a0]*B[b0+7]+0.4*A[a0+3]*B[b0+6]+0.31943828249997*A[a0+7]*B[b0+5]+0.5000000000000001*A[a0+1]*B[b0+5]+0.4472135954999579*A[a0+7]*B[b0+4]+0.4*A[a0+6]*B[b0+3]+0.447213595499958*A[a0+2]*B[b0+3]+0.447213595499958*A[a0+3]*B[b0+2]+0.5000000000000001*A[a0+5]*B[b0+1]+0.5*A[a0+7]*B[b0]; 
  } 
 
  // This tmp allows for in-place multiplication. 
  for (short int i=0; i<8*Ncomp; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide2xSer_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(4,4); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(4); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(4); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*A[0]; 
  AEM(0,1) = 0.5*A[1]; 
  AEM(0,2) = 0.5*A[2]; 
  AEM(0,3) = 0.5*A[3]; 
  AEM(1,0) = 0.5*A[1]; 
  AEM(1,1) = 0.5*A[0]; 
  AEM(1,2) = 0.5*A[3]; 
  AEM(1,3) = 0.5*A[2]; 
  AEM(2,0) = 0.5*A[2]; 
  AEM(2,1) = 0.5*A[3]; 
  AEM(2,2) = 0.5*A[0]; 
  AEM(2,3) = 0.5*A[1]; 
  AEM(3,0) = 0.5*A[3]; 
  AEM(3,1) = 0.5*A[2]; 
  AEM(3,2) = 0.5*A[1]; 
  AEM(3,3) = 0.5*A[0]; 
 
  for(short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 4*vd; 
    // Fill BEV. 
    BEV << B[b0],B[b0+1],B[b0+2],B[b0+3]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*4,4,1) = u; 
 
  } 
} 
 
void CartFieldBinOpDivide2xSer_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM(8,8); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV(8); 
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u(8); 
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*A[0]; 
  AEM(0,1) = 0.5*A[1]; 
  AEM(0,2) = 0.5*A[2]; 
  AEM(0,3) = 0.5*A[3]; 
  AEM(0,4) = 0.5*A[4]; 
  AEM(0,5) = 0.5*A[5]; 
  AEM(0,6) = 0.5*A[6]; 
  AEM(0,7) = 0.5*A[7]; 
  AEM(1,0) = 0.5*A[1]; 
  AEM(1,1) = 0.4472135954999579*A[4]+0.5*A[0]; 
  AEM(1,2) = 0.5*A[3]; 
  AEM(1,3) = 0.447213595499958*A[6]+0.5*A[2]; 
  AEM(1,4) = 0.4472135954999579*A[1]; 
  AEM(1,5) = 0.5000000000000001*A[7]; 
  AEM(1,6) = 0.447213595499958*A[3]; 
  AEM(1,7) = 0.5000000000000001*A[5]; 
  AEM(2,0) = 0.5*A[2]; 
  AEM(2,1) = 0.5*A[3]; 
  AEM(2,2) = 0.4472135954999579*A[5]+0.5*A[0]; 
  AEM(2,3) = 0.447213595499958*A[7]+0.5*A[1]; 
  AEM(2,4) = 0.5000000000000001*A[6]; 
  AEM(2,5) = 0.4472135954999579*A[2]; 
  AEM(2,6) = 0.5000000000000001*A[4]; 
  AEM(2,7) = 0.447213595499958*A[3]; 
  AEM(3,0) = 0.5*A[3]; 
  AEM(3,1) = 0.447213595499958*A[6]+0.5*A[2]; 
  AEM(3,2) = 0.447213595499958*A[7]+0.5*A[1]; 
  AEM(3,3) = 0.4472135954999579*A[5]+0.4472135954999579*A[4]+0.5*A[0]; 
  AEM(3,4) = 0.4472135954999579*A[3]; 
  AEM(3,5) = 0.4472135954999579*A[3]; 
  AEM(3,6) = 0.4*A[7]+0.447213595499958*A[1]; 
  AEM(3,7) = 0.4*A[6]+0.447213595499958*A[2]; 
  AEM(4,0) = 0.5*A[4]; 
  AEM(4,1) = 0.4472135954999579*A[1]; 
  AEM(4,2) = 0.5000000000000001*A[6]; 
  AEM(4,3) = 0.4472135954999579*A[3]; 
  AEM(4,4) = 0.31943828249997*A[4]+0.5*A[0]; 
  AEM(4,6) = 0.31943828249997*A[6]+0.5000000000000001*A[2]; 
  AEM(4,7) = 0.4472135954999579*A[7]; 
  AEM(5,0) = 0.5*A[5]; 
  AEM(5,1) = 0.5000000000000001*A[7]; 
  AEM(5,2) = 0.4472135954999579*A[2]; 
  AEM(5,3) = 0.4472135954999579*A[3]; 
  AEM(5,5) = 0.31943828249997*A[5]+0.5*A[0]; 
  AEM(5,6) = 0.4472135954999579*A[6]; 
  AEM(5,7) = 0.31943828249997*A[7]+0.5000000000000001*A[1]; 
  AEM(6,0) = 0.5*A[6]; 
  AEM(6,1) = 0.447213595499958*A[3]; 
  AEM(6,2) = 0.5000000000000001*A[4]; 
  AEM(6,3) = 0.4*A[7]+0.447213595499958*A[1]; 
  AEM(6,4) = 0.31943828249997*A[6]+0.5000000000000001*A[2]; 
  AEM(6,5) = 0.4472135954999579*A[6]; 
  AEM(6,6) = 0.4472135954999579*A[5]+0.31943828249997*A[4]+0.5*A[0]; 
  AEM(6,7) = 0.4*A[3]; 
  AEM(7,0) = 0.5*A[7]; 
  AEM(7,1) = 0.5000000000000001*A[5]; 
  AEM(7,2) = 0.447213595499958*A[3]; 
  AEM(7,3) = 0.4*A[6]+0.447213595499958*A[2]; 
  AEM(7,4) = 0.4472135954999579*A[7]; 
  AEM(7,5) = 0.31943828249997*A[7]+0.5000000000000001*A[1]; 
  AEM(7,6) = 0.4*A[3]; 
  AEM(7,7) = 0.31943828249997*A[5]+0.4472135954999579*A[4]+0.5*A[0]; 
 
  for(short int vd=0; vd<Ncomp; vd++) 
  { 
    short int b0 = 8*vd; 
    // Fill BEV. 
    BEV << B[b0],B[b0+1],B[b0+2],B[b0+3],B[b0+4],B[b0+5],B[b0+6],B[b0+7]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*8,8,1) = u; 
 
  } 
} 
 
