#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x3vTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[16]; 
 
  tmp[0] = 0.25*(A[15]*B[15]+A[14]*B[14]+A[13]*B[13]+A[12]*B[12]+A[11]*B[11]+A[10]*B[10]+A[9]*B[9]+A[8]*B[8]+A[7]*B[7]+A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.25*(A[14]*B[15]+B[14]*A[15]+A[10]*B[13]+B[10]*A[13]+A[9]*B[12]+B[9]*A[12]+A[7]*B[11]+B[7]*A[11]+A[4]*B[8]+B[4]*A[8]+A[3]*B[6]+B[3]*A[6]+A[2]*B[5]+B[2]*A[5]+A[0]*B[1]+B[0]*A[1]); 
  tmp[2] = 0.25*(A[13]*B[15]+B[13]*A[15]+A[10]*B[14]+B[10]*A[14]+A[8]*B[12]+B[8]*A[12]+A[6]*B[11]+B[6]*A[11]+A[4]*B[9]+B[4]*A[9]+A[3]*B[7]+B[3]*A[7]+A[1]*B[5]+B[1]*A[5]+A[0]*B[2]+B[0]*A[2]); 
  tmp[3] = 0.25*(A[12]*B[15]+B[12]*A[15]+A[9]*B[14]+B[9]*A[14]+A[8]*B[13]+B[8]*A[13]+A[5]*B[11]+B[5]*A[11]+A[4]*B[10]+B[4]*A[10]+A[2]*B[7]+B[2]*A[7]+A[1]*B[6]+B[1]*A[6]+A[0]*B[3]+B[0]*A[3]); 
  tmp[4] = 0.25*(A[11]*B[15]+B[11]*A[15]+A[7]*B[14]+B[7]*A[14]+A[6]*B[13]+B[6]*A[13]+A[5]*B[12]+B[5]*A[12]+A[3]*B[10]+B[3]*A[10]+A[2]*B[9]+B[2]*A[9]+A[1]*B[8]+B[1]*A[8]+A[0]*B[4]+B[0]*A[4]); 
  tmp[5] = 0.25*(A[10]*B[15]+B[10]*A[15]+A[13]*B[14]+B[13]*A[14]+A[4]*B[12]+B[4]*A[12]+A[3]*B[11]+B[3]*A[11]+A[8]*B[9]+B[8]*A[9]+A[6]*B[7]+B[6]*A[7]+A[0]*B[5]+B[0]*A[5]+A[1]*B[2]+B[1]*A[2]); 
  tmp[6] = 0.25*(A[9]*B[15]+B[9]*A[15]+A[12]*B[14]+B[12]*A[14]+A[4]*B[13]+B[4]*A[13]+A[2]*B[11]+B[2]*A[11]+A[8]*B[10]+B[8]*A[10]+A[5]*B[7]+B[5]*A[7]+A[0]*B[6]+B[0]*A[6]+A[1]*B[3]+B[1]*A[3]); 
  tmp[7] = 0.25*(A[8]*B[15]+B[8]*A[15]+A[4]*B[14]+B[4]*A[14]+A[12]*B[13]+B[12]*A[13]+A[1]*B[11]+B[1]*A[11]+A[9]*B[10]+B[9]*A[10]+A[0]*B[7]+B[0]*A[7]+A[5]*B[6]+B[5]*A[6]+A[2]*B[3]+B[2]*A[3]); 
  tmp[8] = 0.25*(A[7]*B[15]+B[7]*A[15]+A[11]*B[14]+B[11]*A[14]+A[3]*B[13]+B[3]*A[13]+A[2]*B[12]+B[2]*A[12]+A[6]*B[10]+B[6]*A[10]+A[5]*B[9]+B[5]*A[9]+A[0]*B[8]+B[0]*A[8]+A[1]*B[4]+B[1]*A[4]); 
  tmp[9] = 0.25*(A[6]*B[15]+B[6]*A[15]+A[3]*B[14]+B[3]*A[14]+A[11]*B[13]+B[11]*A[13]+A[1]*B[12]+B[1]*A[12]+A[7]*B[10]+B[7]*A[10]+A[0]*B[9]+B[0]*A[9]+A[5]*B[8]+B[5]*A[8]+A[2]*B[4]+B[2]*A[4]); 
  tmp[10] = 0.25*(A[5]*B[15]+B[5]*A[15]+A[2]*B[14]+B[2]*A[14]+A[1]*B[13]+B[1]*A[13]+A[11]*B[12]+B[11]*A[12]+A[0]*B[10]+B[0]*A[10]+A[7]*B[9]+B[7]*A[9]+A[6]*B[8]+B[6]*A[8]+A[3]*B[4]+B[3]*A[4]); 
  tmp[11] = 0.25*(A[4]*B[15]+B[4]*A[15]+A[8]*B[14]+B[8]*A[14]+A[9]*B[13]+B[9]*A[13]+A[10]*B[12]+B[10]*A[12]+A[0]*B[11]+B[0]*A[11]+A[1]*B[7]+B[1]*A[7]+A[2]*B[6]+B[2]*A[6]+A[3]*B[5]+B[3]*A[5]); 
  tmp[12] = 0.25*(A[3]*B[15]+B[3]*A[15]+A[6]*B[14]+B[6]*A[14]+A[7]*B[13]+B[7]*A[13]+A[0]*B[12]+B[0]*A[12]+A[10]*B[11]+B[10]*A[11]+A[1]*B[9]+B[1]*A[9]+A[2]*B[8]+B[2]*A[8]+A[4]*B[5]+B[4]*A[5]); 
  tmp[13] = 0.25*(A[2]*B[15]+B[2]*A[15]+A[5]*B[14]+B[5]*A[14]+A[0]*B[13]+B[0]*A[13]+A[7]*B[12]+B[7]*A[12]+A[9]*B[11]+B[9]*A[11]+A[1]*B[10]+B[1]*A[10]+A[3]*B[8]+B[3]*A[8]+A[4]*B[6]+B[4]*A[6]); 
  tmp[14] = 0.25*(A[1]*B[15]+B[1]*A[15]+A[0]*B[14]+B[0]*A[14]+A[5]*B[13]+B[5]*A[13]+A[6]*B[12]+B[6]*A[12]+A[8]*B[11]+B[8]*A[11]+A[2]*B[10]+B[2]*A[10]+A[3]*B[9]+B[3]*A[9]+A[4]*B[7]+B[4]*A[7]); 
  tmp[15] = 0.25*(A[0]*B[15]+B[0]*A[15]+A[1]*B[14]+B[1]*A[14]+A[2]*B[13]+B[2]*A[13]+A[3]*B[12]+B[3]*A[12]+A[4]*B[11]+B[4]*A[11]+A[5]*B[10]+B[5]*A[10]+A[6]*B[9]+B[6]*A[9]+A[7]*B[8]+B[7]*A[8]); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<16; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x3vTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[16]; 
  tmp[0] = 0.7071067811865476*A[1]*B[1]+0.7071067811865476*A[0]*B[0]; 
  tmp[1] = 0.7071067811865476*A[0]*B[1]+0.7071067811865476*B[0]*A[1]; 
  tmp[2] = 0.7071067811865476*A[1]*B[5]+0.7071067811865476*A[0]*B[2]; 
  tmp[3] = 0.7071067811865476*A[1]*B[6]+0.7071067811865476*A[0]*B[3]; 
  tmp[4] = 0.7071067811865476*A[1]*B[8]+0.7071067811865476*A[0]*B[4]; 
  tmp[5] = 0.7071067811865476*A[0]*B[5]+0.7071067811865476*A[1]*B[2]; 
  tmp[6] = 0.7071067811865476*A[0]*B[6]+0.7071067811865476*A[1]*B[3]; 
  tmp[7] = 0.7071067811865476*A[1]*B[11]+0.7071067811865476*A[0]*B[7]; 
  tmp[8] = 0.7071067811865476*A[0]*B[8]+0.7071067811865476*A[1]*B[4]; 
  tmp[9] = 0.7071067811865476*A[1]*B[12]+0.7071067811865476*A[0]*B[9]; 
  tmp[10] = 0.7071067811865476*A[1]*B[13]+0.7071067811865476*A[0]*B[10]; 
  tmp[11] = 0.7071067811865476*A[0]*B[11]+0.7071067811865476*A[1]*B[7]; 
  tmp[12] = 0.7071067811865476*A[0]*B[12]+0.7071067811865476*A[1]*B[9]; 
  tmp[13] = 0.7071067811865476*A[0]*B[13]+0.7071067811865476*A[1]*B[10]; 
  tmp[14] = 0.7071067811865476*A[1]*B[15]+0.7071067811865476*A[0]*B[14]; 
  tmp[15] = 0.7071067811865476*A[0]*B[15]+0.7071067811865476*A[1]*B[14]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<16; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x3vTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
  double Bs[16]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = 0.0; 
    Bs[6] = 0.0; 
    Bs[7] = B[7]; 
    Bs[8] = 0.0; 
    Bs[9] = B[9]; 
    Bs[10] = B[10]; 
    Bs[11] = 0.0; 
    Bs[12] = 0.0; 
    Bs[13] = 0.0; 
    Bs[14] = B[14]; 
    Bs[15] = 0.0; 
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
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
    Bs[10] = B[10]; 
    Bs[11] = B[11]; 
    Bs[12] = B[12]; 
    Bs[13] = B[13]; 
    Bs[14] = B[14]; 
    Bs[15] = B[15]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(16,16); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,2) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,3) = 0.7071067811865475*As[0]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9],Bs[10],Bs[11],Bs[12],Bs[13],Bs[14],Bs[15]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,16,1) = data->u_D; 
 
} 
 
