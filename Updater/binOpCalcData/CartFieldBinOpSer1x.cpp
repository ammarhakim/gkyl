#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    tmp[0] = 0.7071067811865475*(A[a0+1]*B[b0+1]+A[a0]*B[b0]); 
    tmp[1] = 0.7071067811865475*(A[a0]*B[b0+1]+A[a0+1]*B[b0]); 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<2; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply1xSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    tmp[0] = 0.7071067811865475*(A[a0+2]*B[b0+2]+A[a0+1]*B[b0+1]+A[a0]*B[b0]); 
    tmp[1] = 0.1414213562373095*(4.47213595499958*A[a0+1]*B[b0+2]+(4.47213595499958*A[a0+2]+5.0*A[a0])*B[b0+1]+5.0*A[a0+1]*B[b0]); 
    tmp[2] = 0.02020305089104421*((22.3606797749979*A[a0+2]+35.0*A[a0])*B[b0+2]+31.30495168499706*A[a0+1]*B[b0+1]+35.0*A[a0+2]*B[b0]); 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<3; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply1xSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    tmp[0] = 0.7071067811865475*(A[a0+3]*B[b0+3]+A[a0+2]*B[b0+2]+A[a0+1]*B[b0+1]+A[a0]*B[b0]); 
    tmp[1] = 0.02020305089104421*(30.7408522978788*A[a0+2]*B[b0+3]+(30.7408522978788*A[a0+3]+31.30495168499706*A[a0+1])*B[b0+2]+(31.30495168499706*A[a0+2]+35.0*A[a0])*B[b0+1]+35.0*A[a0+1]*B[b0]); 
    tmp[2] = 0.006734350297014738*((62.60990336999411*A[a0+3]+92.22255689363637*A[a0+1])*B[b0+3]+(67.0820393249937*A[a0+2]+105.0*A[a0])*B[b0+2]+(92.22255689363637*A[a0+3]+93.91485505499116*A[a0+1])*B[b0+1]+105.0*A[a0+2]*B[b0]); 
    tmp[3] = 0.006734350297014738*((62.60990336999411*A[a0+2]+105.0*A[a0])*B[b0+3]+(62.60990336999411*A[a0+3]+92.22255689363637*A[a0+1])*B[b0+2]+92.22255689363637*A[a0+2]*B[b0+1]+105.0*A[a0+3]*B[b0]); 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<4; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpDivide1xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.7071067811865475*A[0]-1.224744871391589*A[1] < 0.0) { 
    avgA = true;
  }
  if (1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
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
 
  // Fill AEM matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(2,2); 
  data->AEM_S(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_S(1,0) = 0.7071067811865475*As[1]; 
  data->AEM_S(1,1) = 0.7071067811865475*As[0]; 

  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 2*vd; 
    // Fill BEV. 
    data->BEV_S << Bs[b0],Bs[b0+1]; 

    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*2,2,1) = data->u_S; 
  } 
}

void CartFieldBinOpDivide1xSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
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
 
  // Fill AEM matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(3,3); 
  data->AEM_S(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_S(0,2) = 0.7071067811865475*As[2]; 
  data->AEM_S(1,0) = 0.7071067811865475*As[1]; 
  data->AEM_S(1,1) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_S(1,2) = 0.6324555320336759*As[1]; 
  data->AEM_S(2,0) = 0.7071067811865475*As[2]; 
  data->AEM_S(2,1) = 0.6324555320336759*As[1]; 
  data->AEM_S(2,2) = 0.4517539514526256*As[2]+0.7071067811865475*As[0]; 

  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 3*vd; 
    // Fill BEV. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2]; 

    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*3,3,1) = data->u_S; 
  } 
}

void CartFieldBinOpDivide1xSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.870828693386971*A[3])+1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.870828693386971*A[3]+1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) { 
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
 
  // Fill AEM matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(4,4); 
  data->AEM_S(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_S(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_S(0,2) = 0.7071067811865475*As[2]; 
  data->AEM_S(0,3) = 0.7071067811865475*As[3]; 
  data->AEM_S(1,0) = 0.7071067811865475*As[1]; 
  data->AEM_S(1,1) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_S(1,2) = 0.6210590034081186*As[3]+0.6324555320336759*As[1]; 
  data->AEM_S(1,3) = 0.6210590034081186*As[2]; 
  data->AEM_S(2,0) = 0.7071067811865475*As[2]; 
  data->AEM_S(2,1) = 0.6210590034081186*As[3]+0.6324555320336759*As[1]; 
  data->AEM_S(2,2) = 0.4517539514526256*As[2]+0.7071067811865475*As[0]; 
  data->AEM_S(2,3) = 0.421637021355784*As[3]+0.6210590034081186*As[1]; 
  data->AEM_S(3,0) = 0.7071067811865475*As[3]; 
  data->AEM_S(3,1) = 0.6210590034081186*As[2]; 
  data->AEM_S(3,2) = 0.421637021355784*As[3]+0.6210590034081186*As[1]; 
  data->AEM_S(3,3) = 0.421637021355784*As[2]+0.7071067811865475*As[0]; 

  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 4*vd; 
    // Fill BEV. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3]; 

    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*4,4,1) = data->u_S; 
  } 
}

void CartFieldBinOpDividePositivity1xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If A<0 at corners, but A>0 near positivity control points, use cell average A.
  bool expA = false;
  bool avgA = false;
  if ((0.7071067811865475*A[0]-1.224744871391589*A[1] < 0.0) && (0.7071067811865475*A[0]-0.4898979485566357*A[1] > 0.0)) { 
    expA = true;
  }
  if ((1.224744871391589*A[1]+0.7071067811865475*A[0] < 0.0) && (0.4898979485566357*A[1]+0.7071067811865475*A[0] > 0.0)) { 
    expA = true;
  }
  // If A is zero near positivity control points, use cell average A.
  if (0.7071067811865475*A[0]-0.4898979485566357*A[1] < 0.0) { 
    avgA = true;
  }
  if (0.4898979485566357*A[1]+0.7071067811865475*A[0] < 0.0) { 
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
 
  if ((avgA) || (!expA)) {
    // Fill AEM matrix. 
    data->AEM_S = Eigen::MatrixXd::Zero(2,2); 
    data->AEM_S(0,0) = 0.7071067811865475*As[0]; 
    data->AEM_S(0,1) = 0.7071067811865475*As[1]; 
    data->AEM_S(1,0) = 0.7071067811865475*As[1]; 
    data->AEM_S(1,1) = 0.7071067811865475*As[0]; 

    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 2*vd; 
      // Fill BEV. 
      data->BEV_S << Bs[b0],Bs[b0+1]; 

      // Solve the system of equations. 
      data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
      // Copy data from Eigen vector. 
      Eigen::Map<VectorXd>(out+vd*2,2,1) = data->u_S; 
    } 
  } else {
    double xBar[1];
    double g1[1];
    xBar[0] = (0.5773502691896258*As[1])/As[0]; 

    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBar[0]*xBar[0])-(1.0*xBar[0]*xBar[0]*xBar[0])/(1.0-1.0*xBar[0]*xBar[0]); 

    // Fill AEM matrix. 
    data->AEM_S = Eigen::MatrixXd::Zero(2,2); 
    data->AEM_S(0,0) = 0.7071067811865475*As[0]; 
    data->AEM_S(0,1) = 0.7071067811865475*As[1]; 
    data->AEM_S(1,0) = 0.7071067811865475*As[1]; 
    data->AEM_S(1,1) = 2.121320343559642*As[0]-(2.449489742783178*As[1])/g1[0]; 

    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 2*vd; 
      // Fill BEV. 
      data->BEV_S << Bs[b0],Bs[b0+1]; 

      // Solve the system of equations. 
      data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
      // Copy data from Eigen vector. 
      Eigen::Map<VectorXd>(out+vd*2,2,1) = data->u_S; 
    } 
  };
}

void CartFieldBinOpDotProduct1xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
void CartFieldBinOpDotProduct1xSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
void CartFieldBinOpDotProduct1xSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
