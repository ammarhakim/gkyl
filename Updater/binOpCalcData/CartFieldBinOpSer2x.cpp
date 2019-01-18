#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    tmp[0] = 0.5*A[a0+3]*B[b0+3]+0.5*A[a0+2]*B[b0+2]+0.5*A[a0+1]*B[b0+1]+0.5*A[a0]*B[b0]; 
    tmp[1] = 0.5*A[a0+2]*B[b0+3]+0.5*A[a0+3]*B[b0+2]+0.5*A[a0]*B[b0+1]+0.5*A[a0+1]*B[b0]; 
    tmp[2] = 0.5*A[a0+1]*B[b0+3]+0.5*A[a0]*B[b0+2]+0.5*A[a0+3]*B[b0+1]+0.5*A[a0+2]*B[b0]; 
    tmp[3] = 0.5*A[a0]*B[b0+3]+0.5*A[a0+1]*B[b0+2]+0.5*A[a0+2]*B[b0+1]+0.5*A[a0+3]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<4; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply2xSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[8]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 8*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.5*A[a0+7]*B[b0+7]+0.5*A[a0+6]*B[b0+6]+0.5*A[a0+5]*B[b0+5]+0.5*A[a0+4]*B[b0+4]+0.5*A[a0+3]*B[b0+3]+0.5*A[a0+2]*B[b0+2]+0.5*A[a0+1]*B[b0+1]+0.5*A[a0]*B[b0]; 
    tmp[1] = 0.5000000000000001*A[a0+5]*B[b0+7]+0.447213595499958*A[a0+3]*B[b0+6]+0.5000000000000001*A[a0+7]*B[b0+5]+0.4472135954999579*A[a0+1]*B[b0+4]+0.447213595499958*A[a0+6]*B[b0+3]+0.5*A[a0+2]*B[b0+3]+0.5*A[a0+3]*B[b0+2]+0.4472135954999579*A[a0+4]*B[b0+1]+0.5*A[a0]*B[b0+1]+0.5*A[a0+1]*B[b0]; 
    tmp[2] = 0.447213595499958*A[a0+3]*B[b0+7]+0.5000000000000001*A[a0+4]*B[b0+6]+0.4472135954999579*A[a0+2]*B[b0+5]+0.5000000000000001*A[a0+6]*B[b0+4]+0.447213595499958*A[a0+7]*B[b0+3]+0.5*A[a0+1]*B[b0+3]+0.4472135954999579*A[a0+5]*B[b0+2]+0.5*A[a0]*B[b0+2]+0.5*A[a0+3]*B[b0+1]+0.5*A[a0+2]*B[b0]; 
    tmp[3] = 0.4*A[a0+6]*B[b0+7]+0.447213595499958*A[a0+2]*B[b0+7]+0.4*A[a0+7]*B[b0+6]+0.447213595499958*A[a0+1]*B[b0+6]+0.4472135954999579*A[a0+3]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+4]+0.4472135954999579*A[a0+5]*B[b0+3]+0.4472135954999579*A[a0+4]*B[b0+3]+0.5*A[a0]*B[b0+3]+0.447213595499958*A[a0+7]*B[b0+2]+0.5*A[a0+1]*B[b0+2]+0.447213595499958*A[a0+6]*B[b0+1]+0.5*A[a0+2]*B[b0+1]+0.5*A[a0+3]*B[b0]; 
    tmp[4] = 0.4472135954999579*A[a0+7]*B[b0+7]+0.31943828249997*A[a0+6]*B[b0+6]+0.5000000000000001*A[a0+2]*B[b0+6]+0.31943828249997*A[a0+4]*B[b0+4]+0.5*A[a0]*B[b0+4]+0.4472135954999579*A[a0+3]*B[b0+3]+0.5000000000000001*A[a0+6]*B[b0+2]+0.4472135954999579*A[a0+1]*B[b0+1]+0.5*A[a0+4]*B[b0]; 
    tmp[5] = 0.31943828249997*A[a0+7]*B[b0+7]+0.5000000000000001*A[a0+1]*B[b0+7]+0.4472135954999579*A[a0+6]*B[b0+6]+0.31943828249997*A[a0+5]*B[b0+5]+0.5*A[a0]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+3]+0.4472135954999579*A[a0+2]*B[b0+2]+0.5000000000000001*A[a0+7]*B[b0+1]+0.5*A[a0+5]*B[b0]; 
    tmp[6] = 0.4*A[a0+3]*B[b0+7]+0.4472135954999579*A[a0+5]*B[b0+6]+0.31943828249997*A[a0+4]*B[b0+6]+0.5*A[a0]*B[b0+6]+0.4472135954999579*A[a0+6]*B[b0+5]+0.31943828249997*A[a0+6]*B[b0+4]+0.5000000000000001*A[a0+2]*B[b0+4]+0.4*A[a0+7]*B[b0+3]+0.447213595499958*A[a0+1]*B[b0+3]+0.5000000000000001*A[a0+4]*B[b0+2]+0.447213595499958*A[a0+3]*B[b0+1]+0.5*A[a0+6]*B[b0]; 
    tmp[7] = 0.31943828249997*A[a0+5]*B[b0+7]+0.4472135954999579*A[a0+4]*B[b0+7]+0.5*A[a0]*B[b0+7]+0.4*A[a0+3]*B[b0+6]+0.31943828249997*A[a0+7]*B[b0+5]+0.5000000000000001*A[a0+1]*B[b0+5]+0.4472135954999579*A[a0+7]*B[b0+4]+0.4*A[a0+6]*B[b0+3]+0.447213595499958*A[a0+2]*B[b0+3]+0.447213595499958*A[a0+3]*B[b0+2]+0.5000000000000001*A[a0+5]*B[b0+1]+0.5*A[a0+7]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<8; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpDivide2xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-1.5*A[3])-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-1.5*A[3])-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
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
 
  // Fill AEM_S matrix. 
  data->AEM_S(0,0) = 0.5*As[0]; 
  data->AEM_S(0,1) = 0.5*As[1]; 
  data->AEM_S(0,2) = 0.5*As[2]; 
  data->AEM_S(0,3) = 0.5*As[3]; 
  data->AEM_S(1,0) = 0.5*As[1]; 
  data->AEM_S(1,1) = 0.5*As[0]; 
  data->AEM_S(1,2) = 0.5*As[3]; 
  data->AEM_S(1,3) = 0.5*As[2]; 
  data->AEM_S(2,0) = 0.5*As[2]; 
  data->AEM_S(2,1) = 0.5*As[3]; 
  data->AEM_S(2,2) = 0.5*As[0]; 
  data->AEM_S(2,3) = 0.5*As[1]; 
  data->AEM_S(3,0) = 0.5*As[3]; 
  data->AEM_S(3,1) = 0.5*As[2]; 
  data->AEM_S(3,2) = 0.5*As[1]; 
  data->AEM_S(3,3) = 0.5*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 4*vd; 
    // Fill BEV_S. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3]; 
 
    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*4,4,1) = data->u_S; 
  } 
} 
 
void CartFieldBinOpDivide2xSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.936491673103709*A[7])-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-1.936491673103709*A[7])-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if (1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
 
  unsigned short int b0; 
  double As[8]; 
  double Bs[8*Ncomp]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
    As[6] = 0.0; 
    As[7] = 0.0; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 8*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = 0.0; 
      Bs[b0+2] = 0.0; 
      Bs[b0+3] = 0.0; 
      Bs[b0+4] = 0.0; 
      Bs[b0+5] = 0.0; 
      Bs[b0+6] = 0.0; 
      Bs[b0+7] = 0.0; 
    } 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
    As[6] = A[6]; 
    As[7] = A[7]; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 8*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = B[b0+1]; 
      Bs[b0+2] = B[b0+2]; 
      Bs[b0+3] = B[b0+3]; 
      Bs[b0+4] = B[b0+4]; 
      Bs[b0+5] = B[b0+5]; 
      Bs[b0+6] = B[b0+6]; 
      Bs[b0+7] = B[b0+7]; 
    } 
  } 
 
  // Fill AEM_S matrix. 
  data->AEM_S(0,0) = 0.5*As[0]; 
  data->AEM_S(0,1) = 0.5*As[1]; 
  data->AEM_S(0,2) = 0.5*As[2]; 
  data->AEM_S(0,3) = 0.5*As[3]; 
  data->AEM_S(0,4) = 0.5*As[4]; 
  data->AEM_S(0,5) = 0.5*As[5]; 
  data->AEM_S(0,6) = 0.5*As[6]; 
  data->AEM_S(0,7) = 0.5*As[7]; 
  data->AEM_S(1,0) = 0.5*As[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(1,2) = 0.5*As[3]; 
  data->AEM_S(1,3) = 0.447213595499958*As[6]+0.5*As[2]; 
  data->AEM_S(1,4) = 0.4472135954999579*As[1]; 
  data->AEM_S(1,5) = 0.5000000000000001*As[7]; 
  data->AEM_S(1,6) = 0.447213595499958*As[3]; 
  data->AEM_S(1,7) = 0.5000000000000001*As[5]; 
  data->AEM_S(2,0) = 0.5*As[2]; 
  data->AEM_S(2,1) = 0.5*As[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*As[5]+0.5*As[0]; 
  data->AEM_S(2,3) = 0.447213595499958*As[7]+0.5*As[1]; 
  data->AEM_S(2,4) = 0.5000000000000001*As[6]; 
  data->AEM_S(2,5) = 0.4472135954999579*As[2]; 
  data->AEM_S(2,6) = 0.5000000000000001*As[4]; 
  data->AEM_S(2,7) = 0.447213595499958*As[3]; 
  data->AEM_S(3,0) = 0.5*As[3]; 
  data->AEM_S(3,1) = 0.447213595499958*As[6]+0.5*As[2]; 
  data->AEM_S(3,2) = 0.447213595499958*As[7]+0.5*As[1]; 
  data->AEM_S(3,3) = 0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(3,4) = 0.4472135954999579*As[3]; 
  data->AEM_S(3,5) = 0.4472135954999579*As[3]; 
  data->AEM_S(3,6) = 0.4*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(3,7) = 0.4*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(4,0) = 0.5*As[4]; 
  data->AEM_S(4,1) = 0.4472135954999579*As[1]; 
  data->AEM_S(4,2) = 0.5000000000000001*As[6]; 
  data->AEM_S(4,3) = 0.4472135954999579*As[3]; 
  data->AEM_S(4,4) = 0.31943828249997*As[4]+0.5*As[0]; 
  data->AEM_S(4,6) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  data->AEM_S(4,7) = 0.4472135954999579*As[7]; 
  data->AEM_S(5,0) = 0.5*As[5]; 
  data->AEM_S(5,1) = 0.5000000000000001*As[7]; 
  data->AEM_S(5,2) = 0.4472135954999579*As[2]; 
  data->AEM_S(5,3) = 0.4472135954999579*As[3]; 
  data->AEM_S(5,5) = 0.31943828249997*As[5]+0.5*As[0]; 
  data->AEM_S(5,6) = 0.4472135954999579*As[6]; 
  data->AEM_S(5,7) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  data->AEM_S(6,0) = 0.5*As[6]; 
  data->AEM_S(6,1) = 0.447213595499958*As[3]; 
  data->AEM_S(6,2) = 0.5000000000000001*As[4]; 
  data->AEM_S(6,3) = 0.4*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(6,4) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  data->AEM_S(6,5) = 0.4472135954999579*As[6]; 
  data->AEM_S(6,6) = 0.4472135954999579*As[5]+0.31943828249997*As[4]+0.5*As[0]; 
  data->AEM_S(6,7) = 0.4*As[3]; 
  data->AEM_S(7,0) = 0.5*As[7]; 
  data->AEM_S(7,1) = 0.5000000000000001*As[5]; 
  data->AEM_S(7,2) = 0.447213595499958*As[3]; 
  data->AEM_S(7,3) = 0.4*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(7,4) = 0.4472135954999579*As[7]; 
  data->AEM_S(7,5) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  data->AEM_S(7,6) = 0.4*As[3]; 
  data->AEM_S(7,7) = 0.31943828249997*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 8*vd; 
    // Fill BEV_S. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3],Bs[b0+4],Bs[b0+5],Bs[b0+6],Bs[b0+7]; 
 
    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*8,8,1) = data->u_S; 
  } 
} 
 
void CartFieldBinOpDotProduct2xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    out[0] += 0.5*A[a0+3]*B[a0+3]+0.5*A[a0+2]*B[a0+2]+0.5*A[a0+1]*B[a0+1]+0.5*A[a0]*B[a0]; 
    out[1] += 0.5*A[a0+2]*B[a0+3]+0.5*B[a0+2]*A[a0+3]+0.5*A[a0]*B[a0+1]+0.5*B[a0]*A[a0+1]; 
    out[2] += 0.5*A[a0+1]*B[a0+3]+0.5*B[a0+1]*A[a0+3]+0.5*A[a0]*B[a0+2]+0.5*B[a0]*A[a0+2]; 
    out[3] += 0.5*A[a0]*B[a0+3]+0.5*B[a0]*A[a0+3]+0.5*A[a0+1]*B[a0+2]+0.5*B[a0+1]*A[a0+2]; 
  } 
 
} 
 
void CartFieldBinOpDotProduct2xSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<8; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 8*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.5*A[a0+7]*B[a0+7]+0.5*A[a0+6]*B[a0+6]+0.5*A[a0+5]*B[a0+5]+0.5*A[a0+4]*B[a0+4]+0.5*A[a0+3]*B[a0+3]+0.5*A[a0+2]*B[a0+2]+0.5*A[a0+1]*B[a0+1]+0.5*A[a0]*B[a0]; 
    out[1] += 0.5000000000000001*A[a0+5]*B[a0+7]+0.5000000000000001*B[a0+5]*A[a0+7]+0.447213595499958*A[a0+3]*B[a0+6]+0.447213595499958*B[a0+3]*A[a0+6]+0.4472135954999579*A[a0+1]*B[a0+4]+0.4472135954999579*B[a0+1]*A[a0+4]+0.5*A[a0+2]*B[a0+3]+0.5*B[a0+2]*A[a0+3]+0.5*A[a0]*B[a0+1]+0.5*B[a0]*A[a0+1]; 
    out[2] += 0.447213595499958*A[a0+3]*B[a0+7]+0.447213595499958*B[a0+3]*A[a0+7]+0.5000000000000001*A[a0+4]*B[a0+6]+0.5000000000000001*B[a0+4]*A[a0+6]+0.4472135954999579*A[a0+2]*B[a0+5]+0.4472135954999579*B[a0+2]*A[a0+5]+0.5*A[a0+1]*B[a0+3]+0.5*B[a0+1]*A[a0+3]+0.5*A[a0]*B[a0+2]+0.5*B[a0]*A[a0+2]; 
    out[3] += 0.4*A[a0+6]*B[a0+7]+0.447213595499958*A[a0+2]*B[a0+7]+0.4*B[a0+6]*A[a0+7]+0.447213595499958*B[a0+2]*A[a0+7]+0.447213595499958*A[a0+1]*B[a0+6]+0.447213595499958*B[a0+1]*A[a0+6]+0.4472135954999579*A[a0+3]*B[a0+5]+0.4472135954999579*B[a0+3]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+4]+0.4472135954999579*B[a0+3]*A[a0+4]+0.5*A[a0]*B[a0+3]+0.5*B[a0]*A[a0+3]+0.5*A[a0+1]*B[a0+2]+0.5*B[a0+1]*A[a0+2]; 
    out[4] += 0.4472135954999579*A[a0+7]*B[a0+7]+0.31943828249997*A[a0+6]*B[a0+6]+0.5000000000000001*A[a0+2]*B[a0+6]+0.5000000000000001*B[a0+2]*A[a0+6]+0.31943828249997*A[a0+4]*B[a0+4]+0.5*A[a0]*B[a0+4]+0.5*B[a0]*A[a0+4]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+1]*B[a0+1]; 
    out[5] += 0.31943828249997*A[a0+7]*B[a0+7]+0.5000000000000001*A[a0+1]*B[a0+7]+0.5000000000000001*B[a0+1]*A[a0+7]+0.4472135954999579*A[a0+6]*B[a0+6]+0.31943828249997*A[a0+5]*B[a0+5]+0.5*A[a0]*B[a0+5]+0.5*B[a0]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+2]*B[a0+2]; 
    out[6] += 0.4*A[a0+3]*B[a0+7]+0.4*B[a0+3]*A[a0+7]+0.4472135954999579*A[a0+5]*B[a0+6]+0.31943828249997*A[a0+4]*B[a0+6]+0.5*A[a0]*B[a0+6]+0.4472135954999579*B[a0+5]*A[a0+6]+0.31943828249997*B[a0+4]*A[a0+6]+0.5*B[a0]*A[a0+6]+0.5000000000000001*A[a0+2]*B[a0+4]+0.5000000000000001*B[a0+2]*A[a0+4]+0.447213595499958*A[a0+1]*B[a0+3]+0.447213595499958*B[a0+1]*A[a0+3]; 
    out[7] += 0.31943828249997*A[a0+5]*B[a0+7]+0.4472135954999579*A[a0+4]*B[a0+7]+0.5*A[a0]*B[a0+7]+0.31943828249997*B[a0+5]*A[a0+7]+0.4472135954999579*B[a0+4]*A[a0+7]+0.5*B[a0]*A[a0+7]+0.4*A[a0+3]*B[a0+6]+0.4*B[a0+3]*A[a0+6]+0.5000000000000001*A[a0+1]*B[a0+5]+0.5000000000000001*B[a0+1]*A[a0+5]+0.447213595499958*A[a0+2]*B[a0+3]+0.447213595499958*B[a0+2]*A[a0+3]; 
  } 
 
} 
 
