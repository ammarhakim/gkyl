#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2x2vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[15]; 
 
  tmp[0] = 0.25*(A[14]*B[14]+A[13]*B[13]+A[12]*B[12]+A[11]*B[11]+A[10]*B[10]+A[9]*B[9]+A[8]*B[8]+A[7]*B[7]+A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.05*(4.47213595499958*(A[1]*B[11]+B[1]*A[11])+5.0*(A[4]*B[8]+B[4]*A[8]+A[3]*B[6]+B[3]*A[6]+A[2]*B[5]+B[2]*A[5]+A[0]*B[1]+B[0]*A[1])); 
  tmp[2] = 0.05*(4.47213595499958*(A[2]*B[12]+B[2]*A[12])+5.0*(A[4]*B[9]+B[4]*A[9]+A[3]*B[7]+B[3]*A[7]+A[1]*B[5]+B[1]*A[5]+A[0]*B[2]+B[0]*A[2])); 
  tmp[3] = 0.05*(4.47213595499958*(A[3]*B[13]+B[3]*A[13])+5.0*(A[4]*B[10]+B[4]*A[10]+A[2]*B[7]+B[2]*A[7]+A[1]*B[6]+B[1]*A[6]+A[0]*B[3]+B[0]*A[3])); 
  tmp[4] = 0.05*(4.47213595499958*(A[4]*B[14]+B[4]*A[14])+5.0*(A[3]*B[10]+B[3]*A[10]+A[2]*B[9]+B[2]*A[9]+A[1]*B[8]+B[1]*A[8]+A[0]*B[4]+B[0]*A[4])); 
  tmp[5] = 0.05*(4.47213595499958*(A[5]*B[12]+B[5]*A[12]+A[5]*B[11]+B[5]*A[11])+5.0*(A[8]*B[9]+B[8]*A[9]+A[6]*B[7]+B[6]*A[7]+A[0]*B[5]+B[0]*A[5]+A[1]*B[2]+B[1]*A[2])); 
  tmp[6] = 0.05*(4.47213595499958*(A[6]*B[13]+B[6]*A[13]+A[6]*B[11]+B[6]*A[11])+5.0*(A[8]*B[10]+B[8]*A[10]+A[5]*B[7]+B[5]*A[7]+A[0]*B[6]+B[0]*A[6]+A[1]*B[3]+B[1]*A[3])); 
  tmp[7] = 0.05*(4.47213595499958*(A[7]*B[13]+B[7]*A[13]+A[7]*B[12]+B[7]*A[12])+5.0*(A[9]*B[10]+B[9]*A[10]+A[0]*B[7]+B[0]*A[7]+A[5]*B[6]+B[5]*A[6]+A[2]*B[3]+B[2]*A[3])); 
  tmp[8] = 0.05*(4.47213595499958*(A[8]*B[14]+B[8]*A[14]+A[8]*B[11]+B[8]*A[11])+5.0*(A[6]*B[10]+B[6]*A[10]+A[5]*B[9]+B[5]*A[9]+A[0]*B[8]+B[0]*A[8]+A[1]*B[4]+B[1]*A[4])); 
  tmp[9] = 0.05*(4.47213595499958*(A[9]*B[14]+B[9]*A[14]+A[9]*B[12]+B[9]*A[12])+5.0*(A[7]*B[10]+B[7]*A[10]+A[0]*B[9]+B[0]*A[9]+A[5]*B[8]+B[5]*A[8]+A[2]*B[4]+B[2]*A[4])); 
  tmp[10] = 0.05*(4.47213595499958*(A[10]*B[14]+B[10]*A[14]+A[10]*B[13]+B[10]*A[13])+5.0*(A[0]*B[10]+B[0]*A[10]+A[7]*B[9]+B[7]*A[9]+A[6]*B[8]+B[6]*A[8]+A[3]*B[4]+B[3]*A[4])); 
  tmp[11] = 0.007142857142857143*((22.3606797749979*A[11]+35.0*A[0])*B[11]+35.0*B[0]*A[11]+31.30495168499706*(A[8]*B[8]+A[6]*B[6]+A[5]*B[5]+A[1]*B[1])); 
  tmp[12] = 0.007142857142857143*((22.3606797749979*A[12]+35.0*A[0])*B[12]+35.0*B[0]*A[12]+31.30495168499706*(A[9]*B[9]+A[7]*B[7]+A[5]*B[5]+A[2]*B[2])); 
  tmp[13] = 0.007142857142857143*((22.3606797749979*A[13]+35.0*A[0])*B[13]+35.0*B[0]*A[13]+31.30495168499706*(A[10]*B[10]+A[7]*B[7]+A[6]*B[6]+A[3]*B[3])); 
  tmp[14] = 0.007142857142857143*((22.3606797749979*A[14]+35.0*A[0])*B[14]+35.0*B[0]*A[14]+31.30495168499706*(A[10]*B[10]+A[9]*B[9]+A[8]*B[8]+A[4]*B[4])); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<15; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply2x2vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[15]; 
  tmp[0] = 0.5*A[5]*B[12]+0.5*A[4]*B[11]+0.5*A[3]*B[5]+0.5*A[2]*B[2]+0.5*A[1]*B[1]+0.5*A[0]*B[0]; 
  tmp[1] = 0.447213595499958*A[1]*B[11]+0.5*A[2]*B[5]+0.447213595499958*B[1]*A[4]+0.5*B[2]*A[3]+0.5*A[0]*B[1]+0.5*B[0]*A[1]; 
  tmp[2] = 0.447213595499958*A[2]*B[12]+0.5*A[1]*B[5]+0.447213595499958*B[2]*A[5]+0.5*B[1]*A[3]+0.5*A[0]*B[2]+0.5*B[0]*A[2]; 
  tmp[3] = 0.5*A[2]*B[7]+0.5*A[1]*B[6]+0.5*A[0]*B[3]; 
  tmp[4] = 0.5*A[2]*B[9]+0.5*A[1]*B[8]+0.5*A[0]*B[4]; 
  tmp[5] = 0.447213595499958*A[3]*B[12]+0.447213595499958*A[3]*B[11]+0.447213595499958*A[5]*B[5]+0.447213595499958*A[4]*B[5]+0.5*A[0]*B[5]+0.5*B[0]*A[3]+0.5*A[1]*B[2]+0.5*B[1]*A[2]; 
  tmp[6] = 0.5*A[3]*B[7]+0.447213595499958*A[4]*B[6]+0.5*A[0]*B[6]+0.5*A[1]*B[3]; 
  tmp[7] = 0.447213595499958*A[5]*B[7]+0.5*A[0]*B[7]+0.5*A[3]*B[6]+0.5*A[2]*B[3]; 
  tmp[8] = 0.5*A[3]*B[9]+0.447213595499958*A[4]*B[8]+0.5*A[0]*B[8]+0.5*A[1]*B[4]; 
  tmp[9] = 0.447213595499958*A[5]*B[9]+0.5*A[0]*B[9]+0.5*A[3]*B[8]+0.5*A[2]*B[4]; 
  tmp[10] = 0.5*A[0]*B[10]; 
  tmp[11] = 0.31943828249997*A[4]*B[11]+0.5*A[0]*B[11]+0.4472135954999579*A[3]*B[5]+0.5*B[0]*A[4]+0.4472135954999579*A[1]*B[1]; 
  tmp[12] = 0.31943828249997*A[5]*B[12]+0.5*A[0]*B[12]+0.4472135954999579*A[3]*B[5]+0.5*B[0]*A[5]+0.4472135954999579*A[2]*B[2]; 
  tmp[13] = 0.5*A[0]*B[13]; 
  tmp[14] = 0.5*A[0]*B[14]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<15; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide2x2vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.118033988749895*(A[5]+A[4])+1.5*A[3]-0.8660254037844386*(A[2]+A[1])+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.118033988749895*(A[5]+A[4])-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.118033988749895*(A[5]+A[4])-1.5*A[3]+0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.118033988749895*(A[5]+A[4])+1.5*A[3]+0.8660254037844386*(A[2]+A[1])+0.5*A[0] < 0.0) { 
    avgA = true;
  }
 
  double As[6]; 
  double Bs[15]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = 0.0; 
    Bs[3] = B[3]; 
    Bs[4] = B[4]; 
    Bs[5] = 0.0; 
    Bs[6] = 0.0; 
    Bs[7] = 0.0; 
    Bs[8] = 0.0; 
    Bs[9] = 0.0; 
    Bs[10] = B[10]; 
    Bs[11] = 0.0; 
    Bs[12] = 0.0; 
    Bs[13] = B[13]; 
    Bs[14] = B[14]; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
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
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(15,15); 
  data->AEM_D(0,0) = 0.5*As[0]; 
  data->AEM_D(0,1) = 0.5*As[1]; 
  data->AEM_D(0,2) = 0.5*As[2]; 
  data->AEM_D(0,5) = 0.5*As[3]; 
  data->AEM_D(0,6) = 0.5*As[1]; 
  data->AEM_D(0,7) = 0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_D(0,8) = 0.5*As[3]; 
  data->AEM_D(0,11) = 0.5*As[2]; 
  data->AEM_D(0,12) = 0.5*As[2]; 
  data->AEM_D(0,13) = 0.5*As[3]; 
  data->AEM_D(0,14) = 0.4472135954999579*As[5]+0.5*As[0]; 
  data->AEM_D(1,2) = 0.5*As[1]; 
  data->AEM_D(1,6) = 0.5*As[0]; 
  data->AEM_D(1,13) = 0.5*As[0]; 
  data->AEM_D(2,0) = 0.5*As[3]; 
  data->AEM_D(2,1) = 0.5*As[2]; 
  data->AEM_D(2,2) = 0.5*As[1]; 
  data->AEM_D(2,5) = 0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_D(2,9) = 0.5*As[1]; 
  data->AEM_D(3,0) = 0.5*As[2]; 
  data->AEM_D(3,7) = 0.5*As[1]; 
  data->AEM_D(3,13) = 0.5*As[2]; 
  data->AEM_D(4,6) = 0.5*As[4]; 
  data->AEM_D(4,7) = 0.4472135954999579*As[1]; 
  data->AEM_D(4,11) = 0.4472135954999579*As[3]; 
  data->AEM_D(4,12) = 0.5*As[5]; 
  data->AEM_D(4,14) = 0.4472135954999579*As[2]; 
  data->AEM_D(5,2) = 0.4472135954999579*As[3]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9],Bs[10],Bs[11],Bs[12],Bs[13],Bs[14]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,15,1) = data->u_D; 
 
} 
 
