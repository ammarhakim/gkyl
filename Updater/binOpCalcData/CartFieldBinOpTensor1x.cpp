#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1xTensor_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[2]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 2*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.7071067811865475*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0]*B[b0]; 
    tmp[1] = 0.7071067811865475*A[a0]*B[b0+1]+0.7071067811865475*A[a0+1]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<2; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply1xTensor_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[3]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 3*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.7071067811865475*A[a0+2]*B[b0+2]+0.7071067811865475*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0]*B[b0]; 
    tmp[1] = 0.6324555320336759*A[a0+1]*B[b0+2]+0.6324555320336759*A[a0+2]*B[b0+1]+0.7071067811865475*A[a0]*B[b0+1]+0.7071067811865475*A[a0+1]*B[b0]; 
    tmp[2] = 0.4517539514526256*A[a0+2]*B[b0+2]+0.7071067811865475*A[a0]*B[b0+2]+0.6324555320336759*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0+2]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<3; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply1xTensor_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[4]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 4*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.7071067811865475*A[a0+3]*B[b0+3]+0.7071067811865475*A[a0+2]*B[b0+2]+0.7071067811865475*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0]*B[b0]; 
    tmp[1] = 0.6210590034081186*A[a0+2]*B[b0+3]+0.6210590034081186*A[a0+3]*B[b0+2]+0.6324555320336759*A[a0+1]*B[b0+2]+0.6324555320336759*A[a0+2]*B[b0+1]+0.7071067811865475*A[a0]*B[b0+1]+0.7071067811865475*A[a0+1]*B[b0]; 
    tmp[2] = 0.421637021355784*A[a0+3]*B[b0+3]+0.6210590034081186*A[a0+1]*B[b0+3]+0.4517539514526256*A[a0+2]*B[b0+2]+0.7071067811865475*A[a0]*B[b0+2]+0.6210590034081186*A[a0+3]*B[b0+1]+0.6324555320336759*A[a0+1]*B[b0+1]+0.7071067811865475*A[a0+2]*B[b0]; 
    tmp[3] = 0.421637021355784*A[a0+2]*B[b0+3]+0.7071067811865475*A[a0]*B[b0+3]+0.421637021355784*A[a0+3]*B[b0+2]+0.6210590034081186*A[a0+1]*B[b0+2]+0.6210590034081186*A[a0+2]*B[b0+1]+0.7071067811865475*A[a0+3]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<4; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpDivide1xTensor_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.7071067811865475*A[0]-1.224744871391589*A[1] < 0) { 
    avgA = true;
  }
  if (1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  unsigned short int b0; 
  double As[2]; 
  double Bs[2*Ncomp]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 2*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = 0.0; 
    } 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 2*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = B[b0+1]; 
    } 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(2,2); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(2);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(2);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*As[0]; 
  AEM(0,1) = 0.7071067811865475*As[1]; 
  AEM(1,0) = 0.7071067811865475*As[1]; 
  AEM(1,1) = 0.7071067811865475*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 2*vd; 
    // Fill BEV. 
    BEV << Bs[b0],Bs[b0+1]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*2,2,1) = u; 
  } 
} 
 
void CartFieldBinOpDivide1xTensor_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
  if (1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  unsigned short int b0; 
  double As[3]; 
  double Bs[3*Ncomp]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 3*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = 0.0; 
      Bs[b0+2] = 0.0; 
    } 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 3*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = B[b0+1]; 
      Bs[b0+2] = B[b0+2]; 
    } 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(3,3); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(3);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(3);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.7071067811865475*As[0]; 
  AEM(0,1) = 0.7071067811865475*As[1]; 
  AEM(0,2) = 0.7071067811865475*As[2]; 
  AEM(1,0) = 0.7071067811865475*As[1]; 
  AEM(1,1) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  AEM(1,2) = 0.6324555320336759*As[1]; 
  AEM(2,0) = 0.7071067811865475*As[2]; 
  AEM(2,1) = 0.6324555320336759*As[1]; 
  AEM(2,2) = 0.4517539514526256*As[2]+0.7071067811865475*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 3*vd; 
    // Fill BEV. 
    BEV << Bs[b0],Bs[b0+1],Bs[b0+2]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*3,3,1) = u; 
  } 
} 
 
void CartFieldBinOpDivide1xTensor_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.870828693386971*A[3])+1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
  if (1.870828693386971*A[3]+1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  unsigned short int b0; 
  double As[4]; 
  double Bs[4*Ncomp]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 4*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = 0.0; 
      Bs[b0+2] = 0.0; 
      Bs[b0+3] = 0.0; 
    } 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 4*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = B[b0+1]; 
      Bs[b0+2] = B[b0+2]; 
      Bs[b0+3] = B[b0+3]; 
    } 
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
  AEM(0,2) = 0.7071067811865475*As[2]; 
  AEM(0,3) = 0.7071067811865475*As[3]; 
  AEM(1,0) = 0.7071067811865475*As[1]; 
  AEM(1,1) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  AEM(1,2) = 0.6210590034081186*As[3]+0.6324555320336759*As[1]; 
  AEM(1,3) = 0.6210590034081186*As[2]; 
  AEM(2,0) = 0.7071067811865475*As[2]; 
  AEM(2,1) = 0.6210590034081186*As[3]+0.6324555320336759*As[1]; 
  AEM(2,2) = 0.4517539514526256*As[2]+0.7071067811865475*As[0]; 
  AEM(2,3) = 0.421637021355784*As[3]+0.6210590034081186*As[1]; 
  AEM(3,0) = 0.7071067811865475*As[3]; 
  AEM(3,1) = 0.6210590034081186*As[2]; 
  AEM(3,2) = 0.421637021355784*As[3]+0.6210590034081186*As[1]; 
  AEM(3,3) = 0.421637021355784*As[2]+0.7071067811865475*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 4*vd; 
    // Fill BEV. 
    BEV << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*4,4,1) = u; 
  } 
} 
 
void CartFieldBinOpDotProduct1xTensor_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<2; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.7071067811865475*A[a0+1]*B[a0+1]+0.7071067811865475*A[a0]*B[a0]; 
    out[1] += 0.7071067811865475*A[a0]*B[a0+1]+0.7071067811865475*B[a0]*A[a0+1]; 
  } 
 
} 
 
void CartFieldBinOpDotProduct1xTensor_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<3; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 3*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.7071067811865475*A[a0+2]*B[a0+2]+0.7071067811865475*A[a0+1]*B[a0+1]+0.7071067811865475*A[a0]*B[a0]; 
    out[1] += 0.6324555320336759*A[a0+1]*B[a0+2]+0.6324555320336759*B[a0+1]*A[a0+2]+0.7071067811865475*A[a0]*B[a0+1]+0.7071067811865475*B[a0]*A[a0+1]; 
    out[2] += 0.4517539514526256*A[a0+2]*B[a0+2]+0.7071067811865475*A[a0]*B[a0+2]+0.7071067811865475*B[a0]*A[a0+2]+0.6324555320336759*A[a0+1]*B[a0+1]; 
  } 
 
} 
 
void CartFieldBinOpDotProduct1xTensor_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<4; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.7071067811865475*A[a0+3]*B[a0+3]+0.7071067811865475*A[a0+2]*B[a0+2]+0.7071067811865475*A[a0+1]*B[a0+1]+0.7071067811865475*A[a0]*B[a0]; 
    out[1] += 0.6210590034081186*A[a0+2]*B[a0+3]+0.6210590034081186*B[a0+2]*A[a0+3]+0.6324555320336759*A[a0+1]*B[a0+2]+0.6324555320336759*B[a0+1]*A[a0+2]+0.7071067811865475*A[a0]*B[a0+1]+0.7071067811865475*B[a0]*A[a0+1]; 
    out[2] += 0.421637021355784*A[a0+3]*B[a0+3]+0.6210590034081186*A[a0+1]*B[a0+3]+0.6210590034081186*B[a0+1]*A[a0+3]+0.4517539514526256*A[a0+2]*B[a0+2]+0.7071067811865475*A[a0]*B[a0+2]+0.7071067811865475*B[a0]*A[a0+2]+0.6324555320336759*A[a0+1]*B[a0+1]; 
    out[3] += 0.421637021355784*A[a0+2]*B[a0+3]+0.7071067811865475*A[a0]*B[a0+3]+0.421637021355784*B[a0+2]*A[a0+3]+0.7071067811865475*B[a0]*A[a0+3]+0.6210590034081186*A[a0+1]*B[a0+2]+0.6210590034081186*B[a0+1]*A[a0+2]; 
  } 
 
} 
 
