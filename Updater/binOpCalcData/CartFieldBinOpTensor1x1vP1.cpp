#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x1vTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[4]; 
 
  tmp[0] = 0.5*(A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.5*(A[2]*B[3]+B[2]*A[3]+A[0]*B[1]+B[0]*A[1]); 
  tmp[2] = 0.5*(A[1]*B[3]+B[1]*A[3]+A[0]*B[2]+B[0]*A[2]); 
  tmp[3] = 0.5*(A[0]*B[3]+B[0]*A[3]+A[1]*B[2]+B[1]*A[2]); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<4; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x1vTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[4]; 
  tmp[0] = 0.7071067811865476*A[1]*B[1]+0.7071067811865476*A[0]*B[0]; 
  tmp[1] = 0.7071067811865476*A[0]*B[1]+0.7071067811865476*B[0]*A[1]; 
  tmp[2] = 0.7071067811865476*A[1]*B[3]+0.7071067811865476*A[0]*B[2]; 
  tmp[3] = 0.7071067811865476*A[0]*B[3]+0.7071067811865476*A[1]*B[2]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<4; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x1vTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.7071067811865475*A[0]-1.224744871391589*A[1] < 0.0) { 
    avgA = true;
  }
  if (1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
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
    Bs[3] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(4,4); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,2) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,3) = 0.7071067811865475*As[0]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,4,1) = data->u_D; 
 
} 
 
