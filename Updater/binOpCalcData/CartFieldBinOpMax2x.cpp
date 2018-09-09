#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2xMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
void CartFieldBinOpMultiply2xMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
void CartFieldBinOpMultiply2xMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[10]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 10*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.5*A[a0+9]*B[b0+9]+0.5*A[a0+8]*B[b0+8]+0.5*A[a0+7]*B[b0+7]+0.5*A[a0+6]*B[b0+6]+0.5*A[a0+5]*B[b0+5]+0.5*A[a0+4]*B[b0+4]+0.5*A[a0+3]*B[b0+3]+0.5*A[a0+2]*B[b0+2]+0.5*A[a0+1]*B[b0+1]+0.5*A[a0]*B[b0]; 
    tmp[1] = 0.4391550328268398*A[a0+4]*B[b0+8]+0.5000000000000001*A[a0+5]*B[b0+7]+0.447213595499958*A[a0+3]*B[b0+6]+0.5000000000000001*A[a0+7]*B[b0+5]+0.4391550328268398*A[a0+8]*B[b0+4]+0.4472135954999579*A[a0+1]*B[b0+4]+0.447213595499958*A[a0+6]*B[b0+3]+0.5*A[a0+2]*B[b0+3]+0.5*A[a0+3]*B[b0+2]+0.4472135954999579*A[a0+4]*B[b0+1]+0.5*A[a0]*B[b0+1]+0.5*A[a0+1]*B[b0]; 
    tmp[2] = 0.4391550328268398*A[a0+5]*B[b0+9]+0.447213595499958*A[a0+3]*B[b0+7]+0.5000000000000001*A[a0+4]*B[b0+6]+0.4391550328268398*A[a0+9]*B[b0+5]+0.4472135954999579*A[a0+2]*B[b0+5]+0.5000000000000001*A[a0+6]*B[b0+4]+0.447213595499958*A[a0+7]*B[b0+3]+0.5*A[a0+1]*B[b0+3]+0.4472135954999579*A[a0+5]*B[b0+2]+0.5*A[a0]*B[b0+2]+0.5*A[a0+3]*B[b0+1]+0.5*A[a0+2]*B[b0]; 
    tmp[3] = 0.4391550328268399*A[a0+7]*B[b0+9]+0.4391550328268399*A[a0+6]*B[b0+8]+0.4391550328268399*A[a0+9]*B[b0+7]+0.4*A[a0+6]*B[b0+7]+0.447213595499958*A[a0+2]*B[b0+7]+0.4391550328268399*A[a0+8]*B[b0+6]+0.4*A[a0+7]*B[b0+6]+0.447213595499958*A[a0+1]*B[b0+6]+0.4472135954999579*A[a0+3]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+4]+0.4472135954999579*A[a0+5]*B[b0+3]+0.4472135954999579*A[a0+4]*B[b0+3]+0.5*A[a0]*B[b0+3]+0.447213595499958*A[a0+7]*B[b0+2]+0.5*A[a0+1]*B[b0+2]+0.447213595499958*A[a0+6]*B[b0+1]+0.5*A[a0+2]*B[b0+1]+0.5*A[a0+3]*B[b0]; 
    tmp[4] = 0.2981423969999719*A[a0+8]*B[b0+8]+0.4391550328268398*A[a0+1]*B[b0+8]+0.4472135954999579*A[a0+7]*B[b0+7]+0.31943828249997*A[a0+6]*B[b0+6]+0.5000000000000001*A[a0+2]*B[b0+6]+0.31943828249997*A[a0+4]*B[b0+4]+0.5*A[a0]*B[b0+4]+0.4472135954999579*A[a0+3]*B[b0+3]+0.5000000000000001*A[a0+6]*B[b0+2]+0.4391550328268398*A[a0+8]*B[b0+1]+0.4472135954999579*A[a0+1]*B[b0+1]+0.5*A[a0+4]*B[b0]; 
    tmp[5] = 0.2981423969999719*A[a0+9]*B[b0+9]+0.4391550328268398*A[a0+2]*B[b0+9]+0.31943828249997*A[a0+7]*B[b0+7]+0.5000000000000001*A[a0+1]*B[b0+7]+0.4472135954999579*A[a0+6]*B[b0+6]+0.31943828249997*A[a0+5]*B[b0+5]+0.5*A[a0]*B[b0+5]+0.4472135954999579*A[a0+3]*B[b0+3]+0.4391550328268398*A[a0+9]*B[b0+2]+0.4472135954999579*A[a0+2]*B[b0+2]+0.5000000000000001*A[a0+7]*B[b0+1]+0.5*A[a0+5]*B[b0]; 
    tmp[6] = 0.4391550328268399*A[a0+3]*B[b0+8]+0.4*A[a0+3]*B[b0+7]+0.4472135954999579*A[a0+5]*B[b0+6]+0.31943828249997*A[a0+4]*B[b0+6]+0.5*A[a0]*B[b0+6]+0.4472135954999579*A[a0+6]*B[b0+5]+0.31943828249997*A[a0+6]*B[b0+4]+0.5000000000000001*A[a0+2]*B[b0+4]+0.4391550328268399*A[a0+8]*B[b0+3]+0.4*A[a0+7]*B[b0+3]+0.447213595499958*A[a0+1]*B[b0+3]+0.5000000000000001*A[a0+4]*B[b0+2]+0.447213595499958*A[a0+3]*B[b0+1]+0.5*A[a0+6]*B[b0]; 
    tmp[7] = 0.4391550328268399*A[a0+3]*B[b0+9]+0.31943828249997*A[a0+5]*B[b0+7]+0.4472135954999579*A[a0+4]*B[b0+7]+0.5*A[a0]*B[b0+7]+0.4*A[a0+3]*B[b0+6]+0.31943828249997*A[a0+7]*B[b0+5]+0.5000000000000001*A[a0+1]*B[b0+5]+0.4472135954999579*A[a0+7]*B[b0+4]+0.4391550328268399*A[a0+9]*B[b0+3]+0.4*A[a0+6]*B[b0+3]+0.447213595499958*A[a0+2]*B[b0+3]+0.447213595499958*A[a0+3]*B[b0+2]+0.5000000000000001*A[a0+5]*B[b0+1]+0.5*A[a0+7]*B[b0]; 
    tmp[8] = 0.2981423969999719*A[a0+4]*B[b0+8]+0.5*A[a0]*B[b0+8]+0.4391550328268399*A[a0+3]*B[b0+6]+0.2981423969999719*A[a0+8]*B[b0+4]+0.4391550328268398*A[a0+1]*B[b0+4]+0.4391550328268399*A[a0+6]*B[b0+3]+0.4391550328268398*A[a0+4]*B[b0+1]+0.5*A[a0+8]*B[b0]; 
    tmp[9] = 0.2981423969999719*A[a0+5]*B[b0+9]+0.5*A[a0]*B[b0+9]+0.4391550328268399*A[a0+3]*B[b0+7]+0.2981423969999719*A[a0+9]*B[b0+5]+0.4391550328268398*A[a0+2]*B[b0+5]+0.4391550328268399*A[a0+7]*B[b0+3]+0.4391550328268398*A[a0+5]*B[b0+2]+0.5*A[a0+9]*B[b0]; 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<10; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpDivide2xMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
  double As[3]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(3,3); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(3);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(3);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*As[0]; 
  AEM(0,1) = 0.5*As[1]; 
  AEM(0,2) = 0.5*As[2]; 
  AEM(1,0) = 0.5*As[1]; 
  AEM(1,1) = 0.5*As[0]; 
  AEM(2,0) = 0.5*As[2]; 
  AEM(2,2) = 0.5*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 3*vd; 
    // Fill BEV. 
    BEV << B[b0],B[b0+1],B[b0+2]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*3,3,1) = u; 
  } 
} 
 
void CartFieldBinOpDivide2xMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
  double As[6]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(6,6); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(6);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(6);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*As[0]; 
  AEM(0,1) = 0.5*As[1]; 
  AEM(0,2) = 0.5*As[2]; 
  AEM(0,3) = 0.5*As[3]; 
  AEM(0,4) = 0.5*As[4]; 
  AEM(0,5) = 0.5*As[5]; 
  AEM(1,0) = 0.5*As[1]; 
  AEM(1,1) = 0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(1,2) = 0.5*As[3]; 
  AEM(1,3) = 0.5*As[2]; 
  AEM(1,4) = 0.4472135954999579*As[1]; 
  AEM(2,0) = 0.5*As[2]; 
  AEM(2,1) = 0.5*As[3]; 
  AEM(2,2) = 0.4472135954999579*As[5]+0.5*As[0]; 
  AEM(2,3) = 0.5*As[1]; 
  AEM(2,5) = 0.4472135954999579*As[2]; 
  AEM(3,0) = 0.5*As[3]; 
  AEM(3,1) = 0.5*As[2]; 
  AEM(3,2) = 0.5*As[1]; 
  AEM(3,3) = 0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(3,4) = 0.4472135954999579*As[3]; 
  AEM(3,5) = 0.4472135954999579*As[3]; 
  AEM(4,0) = 0.5*As[4]; 
  AEM(4,1) = 0.4472135954999579*As[1]; 
  AEM(4,3) = 0.4472135954999579*As[3]; 
  AEM(4,4) = 0.31943828249997*As[4]+0.5*As[0]; 
  AEM(5,0) = 0.5*As[5]; 
  AEM(5,2) = 0.4472135954999579*As[2]; 
  AEM(5,3) = 0.4472135954999579*As[3]; 
  AEM(5,5) = 0.31943828249997*As[5]+0.5*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 6*vd; 
    // Fill BEV. 
    BEV << B[b0],B[b0+1],B[b0+2],B[b0+3],B[b0+4],B[b0+5]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*6,6,1) = u; 
  } 
} 
 
void CartFieldBinOpDivide2xMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.322875655532295*A[9])-1.322875655532295*A[8]-1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-1.322875655532295*A[9])-1.322875655532295*A[8]-1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-1.322875655532295*A[9])+1.322875655532295*A[8]+1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
  if ((-1.322875655532295*A[9])+1.322875655532295*A[8]+1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0) { 
    avgA = true;
  }
 
  double As[10]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    As[4] = 0.0; 
    As[5] = 0.0; 
    As[6] = 0.0; 
    As[7] = 0.0; 
    As[8] = 0.0; 
    As[9] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    As[2] = A[2]; 
    As[3] = A[3]; 
    As[4] = A[4]; 
    As[5] = A[5]; 
    As[6] = A[6]; 
    As[7] = A[7]; 
    As[8] = A[8]; 
    As[9] = A[9]; 
  } 
 
  // Declare Eigen Matrix with triple basis tensor dotted with B vector. 
  Eigen::MatrixXd AEM = Eigen::MatrixXd::Zero(10,10); 
  // Declare Eigen Vector with coefficients of B. 
  Eigen::VectorXd BEV = Eigen::VectorXd::Zero(10);  
  // Declare vector with solution to system of equations. 
  Eigen::VectorXd u = Eigen::VectorXd::Zero(10);  
 
  // Fill AEM matrix. 
  AEM(0,0) = 0.5*As[0]; 
  AEM(0,1) = 0.5*As[1]; 
  AEM(0,2) = 0.5*As[2]; 
  AEM(0,3) = 0.5*As[3]; 
  AEM(0,4) = 0.5*As[4]; 
  AEM(0,5) = 0.5*As[5]; 
  AEM(0,6) = 0.5*As[6]; 
  AEM(0,7) = 0.5*As[7]; 
  AEM(0,8) = 0.5*As[8]; 
  AEM(0,9) = 0.5*As[9]; 
  AEM(1,0) = 0.5*As[1]; 
  AEM(1,1) = 0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(1,2) = 0.5*As[3]; 
  AEM(1,3) = 0.447213595499958*As[6]+0.5*As[2]; 
  AEM(1,4) = 0.4391550328268398*As[8]+0.4472135954999579*As[1]; 
  AEM(1,5) = 0.5000000000000001*As[7]; 
  AEM(1,6) = 0.447213595499958*As[3]; 
  AEM(1,7) = 0.5000000000000001*As[5]; 
  AEM(1,8) = 0.4391550328268398*As[4]; 
  AEM(2,0) = 0.5*As[2]; 
  AEM(2,1) = 0.5*As[3]; 
  AEM(2,2) = 0.4472135954999579*As[5]+0.5*As[0]; 
  AEM(2,3) = 0.447213595499958*As[7]+0.5*As[1]; 
  AEM(2,4) = 0.5000000000000001*As[6]; 
  AEM(2,5) = 0.4391550328268398*As[9]+0.4472135954999579*As[2]; 
  AEM(2,6) = 0.5000000000000001*As[4]; 
  AEM(2,7) = 0.447213595499958*As[3]; 
  AEM(2,9) = 0.4391550328268398*As[5]; 
  AEM(3,0) = 0.5*As[3]; 
  AEM(3,1) = 0.447213595499958*As[6]+0.5*As[2]; 
  AEM(3,2) = 0.447213595499958*As[7]+0.5*As[1]; 
  AEM(3,3) = 0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(3,4) = 0.4472135954999579*As[3]; 
  AEM(3,5) = 0.4472135954999579*As[3]; 
  AEM(3,6) = 0.4391550328268399*As[8]+0.4*As[7]+0.447213595499958*As[1]; 
  AEM(3,7) = 0.4391550328268399*As[9]+0.4*As[6]+0.447213595499958*As[2]; 
  AEM(3,8) = 0.4391550328268399*As[6]; 
  AEM(3,9) = 0.4391550328268399*As[7]; 
  AEM(4,0) = 0.5*As[4]; 
  AEM(4,1) = 0.4391550328268398*As[8]+0.4472135954999579*As[1]; 
  AEM(4,2) = 0.5000000000000001*As[6]; 
  AEM(4,3) = 0.4472135954999579*As[3]; 
  AEM(4,4) = 0.31943828249997*As[4]+0.5*As[0]; 
  AEM(4,6) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  AEM(4,7) = 0.4472135954999579*As[7]; 
  AEM(4,8) = 0.2981423969999719*As[8]+0.4391550328268398*As[1]; 
  AEM(5,0) = 0.5*As[5]; 
  AEM(5,1) = 0.5000000000000001*As[7]; 
  AEM(5,2) = 0.4391550328268398*As[9]+0.4472135954999579*As[2]; 
  AEM(5,3) = 0.4472135954999579*As[3]; 
  AEM(5,5) = 0.31943828249997*As[5]+0.5*As[0]; 
  AEM(5,6) = 0.4472135954999579*As[6]; 
  AEM(5,7) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  AEM(5,9) = 0.2981423969999719*As[9]+0.4391550328268398*As[2]; 
  AEM(6,0) = 0.5*As[6]; 
  AEM(6,1) = 0.447213595499958*As[3]; 
  AEM(6,2) = 0.5000000000000001*As[4]; 
  AEM(6,3) = 0.4391550328268399*As[8]+0.4*As[7]+0.447213595499958*As[1]; 
  AEM(6,4) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  AEM(6,5) = 0.4472135954999579*As[6]; 
  AEM(6,6) = 0.4472135954999579*As[5]+0.31943828249997*As[4]+0.5*As[0]; 
  AEM(6,7) = 0.4*As[3]; 
  AEM(6,8) = 0.4391550328268399*As[3]; 
  AEM(7,0) = 0.5*As[7]; 
  AEM(7,1) = 0.5000000000000001*As[5]; 
  AEM(7,2) = 0.447213595499958*As[3]; 
  AEM(7,3) = 0.4391550328268399*As[9]+0.4*As[6]+0.447213595499958*As[2]; 
  AEM(7,4) = 0.4472135954999579*As[7]; 
  AEM(7,5) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  AEM(7,6) = 0.4*As[3]; 
  AEM(7,7) = 0.31943828249997*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  AEM(7,9) = 0.4391550328268399*As[3]; 
  AEM(8,0) = 0.5*As[8]; 
  AEM(8,1) = 0.4391550328268398*As[4]; 
  AEM(8,3) = 0.4391550328268399*As[6]; 
  AEM(8,4) = 0.2981423969999719*As[8]+0.4391550328268398*As[1]; 
  AEM(8,6) = 0.4391550328268399*As[3]; 
  AEM(8,8) = 0.2981423969999719*As[4]+0.5*As[0]; 
  AEM(9,0) = 0.5*As[9]; 
  AEM(9,2) = 0.4391550328268398*As[5]; 
  AEM(9,3) = 0.4391550328268399*As[7]; 
  AEM(9,5) = 0.2981423969999719*As[9]+0.4391550328268398*As[2]; 
  AEM(9,7) = 0.4391550328268399*As[3]; 
  AEM(9,9) = 0.2981423969999719*As[5]+0.5*As[0]; 
 
  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 10*vd; 
    // Fill BEV. 
    BEV << B[b0],B[b0+1],B[b0+2],B[b0+3],B[b0+4],B[b0+5],B[b0+6],B[b0+7],B[b0+8],B[b0+9]; 
 
    // Solve the system of equations. 
    u = AEM.colPivHouseholderQr().solve(BEV); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*10,10,1) = u; 
  } 
} 
 
void CartFieldBinOpDotProduct2xMax_P1(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
void CartFieldBinOpDotProduct2xMax_P2(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
void CartFieldBinOpDotProduct2xMax_P3(const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<10; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 10*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.5*A[a0+9]*B[a0+9]+0.5*A[a0+8]*B[a0+8]+0.5*A[a0+7]*B[a0+7]+0.5*A[a0+6]*B[a0+6]+0.5*A[a0+5]*B[a0+5]+0.5*A[a0+4]*B[a0+4]+0.5*A[a0+3]*B[a0+3]+0.5*A[a0+2]*B[a0+2]+0.5*A[a0+1]*B[a0+1]+0.5*A[a0]*B[a0]; 
    out[1] += 0.4391550328268398*A[a0+4]*B[a0+8]+0.4391550328268398*B[a0+4]*A[a0+8]+0.5000000000000001*A[a0+5]*B[a0+7]+0.5000000000000001*B[a0+5]*A[a0+7]+0.447213595499958*A[a0+3]*B[a0+6]+0.447213595499958*B[a0+3]*A[a0+6]+0.4472135954999579*A[a0+1]*B[a0+4]+0.4472135954999579*B[a0+1]*A[a0+4]+0.5*A[a0+2]*B[a0+3]+0.5*B[a0+2]*A[a0+3]+0.5*A[a0]*B[a0+1]+0.5*B[a0]*A[a0+1]; 
    out[2] += 0.4391550328268398*A[a0+5]*B[a0+9]+0.4391550328268398*B[a0+5]*A[a0+9]+0.447213595499958*A[a0+3]*B[a0+7]+0.447213595499958*B[a0+3]*A[a0+7]+0.5000000000000001*A[a0+4]*B[a0+6]+0.5000000000000001*B[a0+4]*A[a0+6]+0.4472135954999579*A[a0+2]*B[a0+5]+0.4472135954999579*B[a0+2]*A[a0+5]+0.5*A[a0+1]*B[a0+3]+0.5*B[a0+1]*A[a0+3]+0.5*A[a0]*B[a0+2]+0.5*B[a0]*A[a0+2]; 
    out[3] += 0.4391550328268399*A[a0+7]*B[a0+9]+0.4391550328268399*B[a0+7]*A[a0+9]+0.4391550328268399*A[a0+6]*B[a0+8]+0.4391550328268399*B[a0+6]*A[a0+8]+0.4*A[a0+6]*B[a0+7]+0.447213595499958*A[a0+2]*B[a0+7]+0.4*B[a0+6]*A[a0+7]+0.447213595499958*B[a0+2]*A[a0+7]+0.447213595499958*A[a0+1]*B[a0+6]+0.447213595499958*B[a0+1]*A[a0+6]+0.4472135954999579*A[a0+3]*B[a0+5]+0.4472135954999579*B[a0+3]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+4]+0.4472135954999579*B[a0+3]*A[a0+4]+0.5*A[a0]*B[a0+3]+0.5*B[a0]*A[a0+3]+0.5*A[a0+1]*B[a0+2]+0.5*B[a0+1]*A[a0+2]; 
    out[4] += 0.2981423969999719*A[a0+8]*B[a0+8]+0.4391550328268398*A[a0+1]*B[a0+8]+0.4391550328268398*B[a0+1]*A[a0+8]+0.4472135954999579*A[a0+7]*B[a0+7]+0.31943828249997*A[a0+6]*B[a0+6]+0.5000000000000001*A[a0+2]*B[a0+6]+0.5000000000000001*B[a0+2]*A[a0+6]+0.31943828249997*A[a0+4]*B[a0+4]+0.5*A[a0]*B[a0+4]+0.5*B[a0]*A[a0+4]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+1]*B[a0+1]; 
    out[5] += 0.2981423969999719*A[a0+9]*B[a0+9]+0.4391550328268398*A[a0+2]*B[a0+9]+0.4391550328268398*B[a0+2]*A[a0+9]+0.31943828249997*A[a0+7]*B[a0+7]+0.5000000000000001*A[a0+1]*B[a0+7]+0.5000000000000001*B[a0+1]*A[a0+7]+0.4472135954999579*A[a0+6]*B[a0+6]+0.31943828249997*A[a0+5]*B[a0+5]+0.5*A[a0]*B[a0+5]+0.5*B[a0]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+2]*B[a0+2]; 
    out[6] += 0.4391550328268399*A[a0+3]*B[a0+8]+0.4391550328268399*B[a0+3]*A[a0+8]+0.4*A[a0+3]*B[a0+7]+0.4*B[a0+3]*A[a0+7]+0.4472135954999579*A[a0+5]*B[a0+6]+0.31943828249997*A[a0+4]*B[a0+6]+0.5*A[a0]*B[a0+6]+0.4472135954999579*B[a0+5]*A[a0+6]+0.31943828249997*B[a0+4]*A[a0+6]+0.5*B[a0]*A[a0+6]+0.5000000000000001*A[a0+2]*B[a0+4]+0.5000000000000001*B[a0+2]*A[a0+4]+0.447213595499958*A[a0+1]*B[a0+3]+0.447213595499958*B[a0+1]*A[a0+3]; 
    out[7] += 0.4391550328268399*A[a0+3]*B[a0+9]+0.4391550328268399*B[a0+3]*A[a0+9]+0.31943828249997*A[a0+5]*B[a0+7]+0.4472135954999579*A[a0+4]*B[a0+7]+0.5*A[a0]*B[a0+7]+0.31943828249997*B[a0+5]*A[a0+7]+0.4472135954999579*B[a0+4]*A[a0+7]+0.5*B[a0]*A[a0+7]+0.4*A[a0+3]*B[a0+6]+0.4*B[a0+3]*A[a0+6]+0.5000000000000001*A[a0+1]*B[a0+5]+0.5000000000000001*B[a0+1]*A[a0+5]+0.447213595499958*A[a0+2]*B[a0+3]+0.447213595499958*B[a0+2]*A[a0+3]; 
    out[8] += 0.2981423969999719*A[a0+4]*B[a0+8]+0.5*A[a0]*B[a0+8]+0.2981423969999719*B[a0+4]*A[a0+8]+0.5*B[a0]*A[a0+8]+0.4391550328268399*A[a0+3]*B[a0+6]+0.4391550328268399*B[a0+3]*A[a0+6]+0.4391550328268398*A[a0+1]*B[a0+4]+0.4391550328268398*B[a0+1]*A[a0+4]; 
    out[9] += 0.2981423969999719*A[a0+5]*B[a0+9]+0.5*A[a0]*B[a0+9]+0.2981423969999719*B[a0+5]*A[a0+9]+0.5*B[a0]*A[a0+9]+0.4391550328268399*A[a0+3]*B[a0+7]+0.4391550328268399*B[a0+3]*A[a0+7]+0.4391550328268398*A[a0+2]*B[a0+5]+0.4391550328268398*B[a0+2]*A[a0+5]; 
  } 
 
} 
 
