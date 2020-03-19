#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x2vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[10]; 
 
  tmp[0] = 0.3535533905932737*(A[9]*B[9]+A[8]*B[8]+A[7]*B[7]+A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.07071067811865474*(4.47213595499958*(A[1]*B[7]+B[1]*A[7])+5.0*(A[3]*B[5]+B[3]*A[5]+A[2]*B[4]+B[2]*A[4]+A[0]*B[1]+B[0]*A[1])); 
  tmp[2] = 0.07071067811865474*(4.47213595499958*(A[2]*B[8]+B[2]*A[8])+5.0*(A[3]*B[6]+B[3]*A[6]+A[1]*B[4]+B[1]*A[4]+A[0]*B[2]+B[0]*A[2])); 
  tmp[3] = 0.07071067811865474*(4.47213595499958*(A[3]*B[9]+B[3]*A[9])+5.0*(A[2]*B[6]+B[2]*A[6]+A[1]*B[5]+B[1]*A[5]+A[0]*B[3]+B[0]*A[3])); 
  tmp[4] = 0.07071067811865474*(4.47213595499958*(A[4]*B[8]+B[4]*A[8]+A[4]*B[7]+B[4]*A[7])+5.0*(A[5]*B[6]+B[5]*A[6]+A[0]*B[4]+B[0]*A[4]+A[1]*B[2]+B[1]*A[2])); 
  tmp[5] = 0.07071067811865474*(4.47213595499958*(A[5]*B[9]+B[5]*A[9]+A[5]*B[7]+B[5]*A[7])+5.0*(A[4]*B[6]+B[4]*A[6]+A[0]*B[5]+B[0]*A[5]+A[1]*B[3]+B[1]*A[3])); 
  tmp[6] = 0.07071067811865474*(4.47213595499958*(A[6]*B[9]+B[6]*A[9]+A[6]*B[8]+B[6]*A[8])+5.0*(A[0]*B[6]+B[0]*A[6]+A[4]*B[5]+B[4]*A[5]+A[2]*B[3]+B[2]*A[3])); 
  tmp[7] = 0.01010152544552211*((22.3606797749979*A[7]+35.0*A[0])*B[7]+35.0*B[0]*A[7]+31.30495168499706*(A[5]*B[5]+A[4]*B[4]+A[1]*B[1])); 
  tmp[8] = 0.01010152544552211*((22.3606797749979*A[8]+35.0*A[0])*B[8]+35.0*B[0]*A[8]+31.30495168499706*(A[6]*B[6]+A[4]*B[4]+A[2]*B[2])); 
  tmp[9] = 0.01010152544552211*((22.3606797749979*A[9]+35.0*A[0])*B[9]+35.0*B[0]*A[9]+31.30495168499706*(A[6]*B[6]+A[5]*B[5]+A[3]*B[3])); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<10; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x2vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[10]; 
  tmp[0] = 0.7071067811865475*A[2]*B[7]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6324555320336759*A[1]*B[7]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[1]*B[4]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[1]*B[5]+0.7071067811865475*A[0]*B[3]; 
  tmp[4] = 0.6324555320336759*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.7071067811865475*A[1]*B[2]; 
  tmp[5] = 0.6324555320336759*A[2]*B[5]+0.7071067811865475*A[0]*B[5]+0.7071067811865475*A[1]*B[3]; 
  tmp[6] = 0.7071067811865475*A[0]*B[6]; 
  tmp[7] = 0.4517539514526257*A[2]*B[7]+0.7071067811865475*A[0]*B[7]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[8] = 0.7071067811865475*A[0]*B[8]; 
  tmp[9] = 0.7071067811865475*A[0]*B[9]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<10; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x2vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
    avgA = true;
  }
 
  double As[3]; 
  double Bs[10]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
    Bs[4] = 0.0; 
    Bs[5] = 0.0; 
    Bs[6] = B[6]; 
    Bs[7] = 0.0; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
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
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(10,10); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,3) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,4) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_D(0,8) = 0.7071067811865475*As[0]; 
  data->AEM_D(1,4) = 0.7071067811865475*As[1]; 
  data->AEM_D(2,1) = 0.7071067811865475*As[2]; 
  data->AEM_D(2,2) = 0.6324555320336759*As[1]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,10,1) = data->u_D; 
 
} 
 
