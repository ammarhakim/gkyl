#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2xMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    tmp[0] = 0.5*A[a0+2]*B[b0+2]+0.5*A[a0+1]*B[b0+1]+0.5*A[a0]*B[b0]; 
    tmp[1] = 0.5*A[a0]*B[b0+1]+0.5*A[a0+1]*B[b0]; 
    tmp[2] = 0.5*A[a0]*B[b0+2]+0.5*A[a0+2]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<3; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply2xMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[6]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 6*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.5*A[a0+5]*B[b0+5]+0.5*A[a0+4]*B[b0+4]+0.5*A[a0+3]*B[b0+3]+0.5*A[a0+2]*B[b0+2]+0.5*A[a0+1]*B[b0+1]+0.5*A[a0]*B[b0]; 
    tmp[1] = 0.4472135954999579*A[a0+1]*B[b0+4]+0.5*A[a0+2]*B[b0+3]+0.5*A[a0+3]*B[b0+2]+0.4472135954999579*A[a0+4]*B[b0+1]+0.5*A[a0]*B[b0+1]+0.5*A[a0+1]*B[b0]; 
    tmp[2] = 0.4472135954999579*A[a0+2]*B[b0+5]+0.5*A[a0+1]*B[b0+3]+0.4472135954999579*A[a0+5]*B[b0+2]+0.5*A[a0]*B[b0+2]+0.5*A[a0+3]*B[b0+1]+0.5*A[a0+2]*B[b0]; 
    tmp[3] = 0.4472135954999579*A[a0+3]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+4]+0.4472135954999579*A[a0+5]*B[b0+3]+0.4472135954999579*A[a0+4]*B[b0+3]+0.5*A[a0]*B[b0+3]+0.5*A[a0+1]*B[b0+2]+0.5*A[a0+2]*B[b0+1]+0.5*A[a0+3]*B[b0]; 
    tmp[4] = 0.31943828249997*A[a0+4]*B[b0+4]+0.5*A[a0]*B[b0+4]+0.4472135954999579*A[a0+3]*B[b0+3]+0.4472135954999579*A[a0+1]*B[b0+1]+0.5*A[a0+4]*B[b0]; 
    tmp[5] = 0.31943828249997*A[a0+5]*B[b0+5]+0.5*A[a0]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+3]+0.4472135954999579*A[a0+2]*B[b0+2]+0.5*A[a0+5]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<6; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpDivide2xMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-0.8660254037844386*A[2])-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.8660254037844386*A[2])-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.8660254037844386*A[2])+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.8660254037844386*A[2])+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
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
 
  // Fill AEM_S matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(3,3);
  data->AEM_S(0,0) = 0.5*As[0]; 
  data->AEM_S(0,1) = 0.5*As[1]; 
  data->AEM_S(0,2) = 0.5*As[2]; 
  data->AEM_S(1,0) = 0.5*As[1]; 
  data->AEM_S(1,1) = 0.5*As[0]; 
  data->AEM_S(2,0) = 0.5*As[2]; 
  data->AEM_S(2,2) = 0.5*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 3*vd; 
    // Fill BEV_S. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2]; 
 
    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*3,3,1) = data->u_S; 
  } 
} 
 
void CartFieldBinOpDivide2xMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
 
  unsigned short int b0; 
  double As[6]; 
  double Bs[6*Ncomp]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 6*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = 0.0; 
      Bs[b0+2] = 0.0; 
      Bs[b0+3] = 0.0; 
      Bs[b0+4] = 0.0; 
      Bs[b0+5] = 0.0; 
    } 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 6*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = B[b0+1]; 
      Bs[b0+2] = B[b0+2]; 
      Bs[b0+3] = B[b0+3]; 
      Bs[b0+4] = B[b0+4]; 
      Bs[b0+5] = B[b0+5]; 
    } 
  } 
 
  // Fill AEM_S matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(6,6);
  data->AEM_S(0,0) = 0.5*As[0]; 
  data->AEM_S(0,1) = 0.5*As[1]; 
  data->AEM_S(0,2) = 0.5*As[2]; 
  data->AEM_S(0,3) = 0.5*As[3]; 
  data->AEM_S(0,4) = 0.5*As[4]; 
  data->AEM_S(0,5) = 0.5*As[5]; 
  data->AEM_S(1,0) = 0.5*As[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(1,2) = 0.5*As[3]; 
  data->AEM_S(1,3) = 0.5*As[2]; 
  data->AEM_S(1,4) = 0.4472135954999579*As[1]; 
  data->AEM_S(2,0) = 0.5*As[2]; 
  data->AEM_S(2,1) = 0.5*As[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*As[5]+0.5*As[0]; 
  data->AEM_S(2,3) = 0.5*As[1]; 
  data->AEM_S(2,5) = 0.4472135954999579*As[2]; 
  data->AEM_S(3,0) = 0.5*As[3]; 
  data->AEM_S(3,1) = 0.5*As[2]; 
  data->AEM_S(3,2) = 0.5*As[1]; 
  data->AEM_S(3,3) = 0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(3,4) = 0.4472135954999579*As[3]; 
  data->AEM_S(3,5) = 0.4472135954999579*As[3]; 
  data->AEM_S(4,0) = 0.5*As[4]; 
  data->AEM_S(4,1) = 0.4472135954999579*As[1]; 
  data->AEM_S(4,3) = 0.4472135954999579*As[3]; 
  data->AEM_S(4,4) = 0.31943828249997*As[4]+0.5*As[0]; 
  data->AEM_S(5,0) = 0.5*As[5]; 
  data->AEM_S(5,2) = 0.4472135954999579*As[2]; 
  data->AEM_S(5,3) = 0.4472135954999579*As[3]; 
  data->AEM_S(5,5) = 0.31943828249997*As[5]+0.5*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 6*vd; 
    // Fill BEV_S. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3],Bs[b0+4],Bs[b0+5]; 
 
    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*6,6,1) = data->u_S; 
  } 
} 
 
void CartFieldBinOpDotProduct2xMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    out[0] += 0.5*A[a0+2]*B[a0+2]+0.5*A[a0+1]*B[a0+1]+0.5*A[a0]*B[a0]; 
    out[1] += 0.5*A[a0]*B[a0+1]+0.5*B[a0]*A[a0+1]; 
    out[2] += 0.5*A[a0]*B[a0+2]+0.5*B[a0]*A[a0+2]; 
  } 
 
} 
 
void CartFieldBinOpDotProduct2xMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<6; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 6*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.5*A[a0+5]*B[a0+5]+0.5*A[a0+4]*B[a0+4]+0.5*A[a0+3]*B[a0+3]+0.5*A[a0+2]*B[a0+2]+0.5*A[a0+1]*B[a0+1]+0.5*A[a0]*B[a0]; 
    out[1] += 0.4472135954999579*A[a0+1]*B[a0+4]+0.4472135954999579*B[a0+1]*A[a0+4]+0.5*A[a0+2]*B[a0+3]+0.5*B[a0+2]*A[a0+3]+0.5*A[a0]*B[a0+1]+0.5*B[a0]*A[a0+1]; 
    out[2] += 0.4472135954999579*A[a0+2]*B[a0+5]+0.4472135954999579*B[a0+2]*A[a0+5]+0.5*A[a0+1]*B[a0+3]+0.5*B[a0+1]*A[a0+3]+0.5*A[a0]*B[a0+2]+0.5*B[a0]*A[a0+2]; 
    out[3] += 0.4472135954999579*A[a0+3]*B[a0+5]+0.4472135954999579*B[a0+3]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+4]+0.4472135954999579*B[a0+3]*A[a0+4]+0.5*A[a0]*B[a0+3]+0.5*B[a0]*A[a0+3]+0.5*A[a0+1]*B[a0+2]+0.5*B[a0+1]*A[a0+2]; 
    out[4] += 0.31943828249997*A[a0+4]*B[a0+4]+0.5*A[a0]*B[a0+4]+0.5*B[a0]*A[a0+4]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+1]*B[a0+1]; 
    out[5] += 0.31943828249997*A[a0+5]*B[a0+5]+0.5*A[a0]*B[a0+5]+0.5*B[a0]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+2]*B[a0+2]; 
  } 
 
} 
 
