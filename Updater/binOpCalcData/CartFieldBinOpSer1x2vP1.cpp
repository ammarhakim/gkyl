#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x2vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[8]; 
 
  tmp[0] = 0.3535533905932737*(A[7]*B[7]+A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.3535533905932737*(A[6]*B[7]+B[6]*A[7]+A[3]*B[5]+B[3]*A[5]+A[2]*B[4]+B[2]*A[4]+A[0]*B[1]+B[0]*A[1]); 
  tmp[2] = 0.3535533905932737*(A[5]*B[7]+B[5]*A[7]+A[3]*B[6]+B[3]*A[6]+A[1]*B[4]+B[1]*A[4]+A[0]*B[2]+B[0]*A[2]); 
  tmp[3] = 0.3535533905932737*(A[4]*B[7]+B[4]*A[7]+A[2]*B[6]+B[2]*A[6]+A[1]*B[5]+B[1]*A[5]+A[0]*B[3]+B[0]*A[3]); 
  tmp[4] = 0.3535533905932737*(A[3]*B[7]+B[3]*A[7]+A[5]*B[6]+B[5]*A[6]+A[0]*B[4]+B[0]*A[4]+A[1]*B[2]+B[1]*A[2]); 
  tmp[5] = 0.3535533905932737*(A[2]*B[7]+B[2]*A[7]+A[4]*B[6]+B[4]*A[6]+A[0]*B[5]+B[0]*A[5]+A[1]*B[3]+B[1]*A[3]); 
  tmp[6] = 0.3535533905932737*(A[1]*B[7]+B[1]*A[7]+A[0]*B[6]+B[0]*A[6]+A[4]*B[5]+B[4]*A[5]+A[2]*B[3]+B[2]*A[3]); 
  tmp[7] = 0.3535533905932737*(A[0]*B[7]+B[0]*A[7]+A[1]*B[6]+B[1]*A[6]+A[2]*B[5]+B[2]*A[5]+A[3]*B[4]+B[3]*A[4]); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<8; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x2vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
  for (unsigned short int i=0; i<8; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x2vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
  double Bs[8]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = 0.0; 
    Bs[5] = 0.0; 
    Bs[6] = B[6]; 
    Bs[7] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
    Bs[7] = B[7]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(8,8); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,2) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,3) = 0.7071067811865475*As[0]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,8,1) = data->u_D; 
 
} 
 
