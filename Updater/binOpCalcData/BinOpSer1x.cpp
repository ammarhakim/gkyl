#include <math.h> 
#include <BinOpModDecl.h> 
 
using namespace Eigen; 
 
void BinOpMultiply1xSer_P1(const double *A, const double *B, const short int vDim, double *out) 
{ 
  for(short int vd=0; vd<vDim; vd++) 
  { 
    short int c0 = 2*vd; 
    // Component-wise (of the vectors) multiplication. 
    out[c0] = 0.7071067811865475*A[c0+1]*B[c0+1]+0.7071067811865475*A[c0]*B[c0]; 
    out[c0+1] = 0.7071067811865475*A[c0]*B[c0+1]+0.7071067811865475*B[c0]*A[c0+1]; 
  } 
} 
 
void BinOpMultiply1xSer_P2(const double *A, const double *B, const short int vDim, double *out) 
{ 
  for(short int vd=0; vd<vDim; vd++) 
  { 
    short int c0 = 3*vd; 
    // Component-wise (of the vectors) multiplication. 
    out[c0] = 0.7071067811865475*A[c0+2]*B[c0+2]+0.7071067811865475*A[c0+1]*B[c0+1]+0.7071067811865475*A[c0]*B[c0]; 
    out[c0+1] = 0.6324555320336759*A[c0+1]*B[c0+2]+0.6324555320336759*B[c0+1]*A[c0+2]+0.7071067811865475*A[c0]*B[c0+1]+0.7071067811865475*B[c0]*A[c0+1]; 
    out[c0+2] = 0.4517539514526256*A[c0+2]*B[c0+2]+0.7071067811865475*A[c0]*B[c0+2]+0.7071067811865475*B[c0]*A[c0+2]+0.6324555320336759*A[c0+1]*B[c0+1]; 
  } 
} 
 
void BinOpDivide1xSer_P1(const double *A, const double *B, const short int vDim, double *out) 
{ 
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
 
  for(short int vd=0; vd<vDim; vd++) 
  { 
    // Fill BEV. 
    BEV << B[2*vd],B[2*vd+1]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*2,2,1) = u; 
 
  } 
} 
 
void BinOpDivide1xSer_P2(const double *A, const double *B, const short int vDim, double *out) 
{ 
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
 
  for(short int vd=0; vd<vDim; vd++) 
  { 
    // Fill BEV. 
    BEV << B[3*vd],B[3*vd+1],B[3*vd+2]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*3,3,1) = u; 
 
  } 
} 
 
