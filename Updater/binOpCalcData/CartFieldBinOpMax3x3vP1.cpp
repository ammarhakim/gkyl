#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply3x3vMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[7]; 
 
  tmp[0] = 0.125*(A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.125*(A[0]*B[1]+B[0]*A[1]); 
  tmp[2] = 0.125*(A[0]*B[2]+B[0]*A[2]); 
  tmp[3] = 0.125*(A[0]*B[3]+B[0]*A[3]); 
  tmp[4] = 0.125*(A[0]*B[4]+B[0]*A[4]); 
  tmp[5] = 0.125*(A[0]*B[5]+B[0]*A[5]); 
  tmp[6] = 0.125*(A[0]*B[6]+B[0]*A[6]); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<7; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply3x3vMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[7]; 
  tmp[0] = 0.3535533905932738*A[3]*B[3]+0.3535533905932738*A[2]*B[2]+0.3535533905932738*A[1]*B[1]+0.3535533905932738*A[0]*B[0]; 
  tmp[1] = 0.3535533905932738*A[0]*B[1]+0.3535533905932738*B[0]*A[1]; 
  tmp[2] = 0.3535533905932738*A[0]*B[2]+0.3535533905932738*B[0]*A[2]; 
  tmp[3] = 0.3535533905932738*A[0]*B[3]+0.3535533905932738*B[0]*A[3]; 
  tmp[4] = 0.3535533905932737*A[0]*B[4]; 
  tmp[5] = 0.3535533905932737*A[0]*B[5]; 
  tmp[6] = 0.3535533905932737*A[0]*B[6]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<7; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide3x3vMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.3535533905932737*A[0]-0.6123724356957944*(A[3]+A[2]+A[1]) < 0.0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*(A[3]+A[2]))+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])+0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])+0.6123724356957944*(A[2]+A[1])+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (0.6123724356957944*A[3]-0.6123724356957944*(A[2]+A[1])+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (0.6123724356957944*(A[3]+A[2])-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
  if (0.6123724356957944*(A[3]+A[2]+A[1])+0.3535533905932737*A[0] < 0.0) { 
    avgA = true;
  }
 
  double As[4]; 
  double Bs[7]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = 0.0; 
    Bs[3] = 0.0; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(7,7); 
  data->AEM_D(0,0) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,1) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,2) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,3) = 0.3535533905932737*As[3]; 
  data->AEM_D(0,4) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,5) = 0.3535533905932737*As[0]; 
  data->AEM_D(1,1) = 0.3535533905932737*As[2]; 
  data->AEM_D(1,3) = 0.3535533905932737*As[0]; 
  data->AEM_D(1,5) = 0.3535533905932737*As[3]; 
  data->AEM_D(2,1) = 0.3535533905932737*As[0]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,7,1) = data->u_D; 
 
} 
 
