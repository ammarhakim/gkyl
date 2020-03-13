#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply2xTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
    tmp[0] = 0.5*(A[a0+3]*B[b0+3]+A[a0+2]*B[b0+2]+A[a0+1]*B[b0+1]+A[a0]*B[b0]); 
    tmp[1] = 0.5*(A[a0+2]*B[b0+3]+A[a0+3]*B[b0+2]+A[a0]*B[b0+1]+A[a0+1]*B[b0]); 
    tmp[2] = 0.5*(A[a0+1]*B[b0+3]+A[a0]*B[b0+2]+A[a0+3]*B[b0+1]+A[a0+2]*B[b0]); 
    tmp[3] = 0.5*(A[a0]*B[b0+3]+A[a0+1]*B[b0+2]+A[a0+2]*B[b0+1]+A[a0+3]*B[b0]); 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<4; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply2xTensor_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[9]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 9*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.5*(A[a0+8]*B[b0+8]+A[a0+7]*B[b0+7]+A[a0+6]*B[b0+6]+A[a0+5]*B[b0+5]+A[a0+4]*B[b0+4]+A[a0+3]*B[b0+3]+A[a0+2]*B[b0+2]+A[a0+1]*B[b0+1]+A[a0]*B[b0]); 
    tmp[1] = 0.03333333333333333*(13.41640786499874*A[a0+7]*B[b0+8]+(13.41640786499874*A[a0+8]+15.0*A[a0+5])*B[b0+7]+13.41640786499874*A[a0+3]*B[b0+6]+15.0*A[a0+7]*B[b0+5]+13.41640786499874*A[a0+1]*B[b0+4]+(13.41640786499874*A[a0+6]+15.0*A[a0+2])*B[b0+3]+15.0*A[a0+3]*B[b0+2]+(13.41640786499874*A[a0+4]+15.0*A[a0])*B[b0+1]+15.0*A[a0+1]*B[b0]); 
    tmp[2] = 0.03333333333333333*(13.41640786499874*A[a0+6]*B[b0+8]+13.41640786499874*A[a0+3]*B[b0+7]+(13.41640786499874*A[a0+8]+15.0*A[a0+4])*B[b0+6]+13.41640786499874*A[a0+2]*B[b0+5]+15.0*A[a0+6]*B[b0+4]+(13.41640786499874*A[a0+7]+15.0*A[a0+1])*B[b0+3]+(13.41640786499874*A[a0+5]+15.0*A[a0])*B[b0+2]+15.0*A[a0+3]*B[b0+1]+15.0*A[a0+2]*B[b0]); 
    tmp[3] = 0.03333333333333333*(12.0*A[a0+3]*B[b0+8]+(12.0*A[a0+6]+13.41640786499874*A[a0+2])*B[b0+7]+(12.0*A[a0+7]+13.41640786499874*A[a0+1])*B[b0+6]+13.41640786499874*A[a0+3]*B[b0+5]+13.41640786499874*A[a0+3]*B[b0+4]+(12.0*A[a0+8]+13.41640786499874*A[a0+5]+13.41640786499874*A[a0+4]+15.0*A[a0])*B[b0+3]+(13.41640786499874*A[a0+7]+15.0*A[a0+1])*B[b0+2]+(13.41640786499874*A[a0+6]+15.0*A[a0+2])*B[b0+1]+15.0*A[a0+3]*B[b0]); 
    tmp[4] = 0.004761904761904762*((67.0820393249937*A[a0+8]+105.0*A[a0+5])*B[b0+8]+93.91485505499116*A[a0+7]*B[b0+7]+(67.0820393249937*A[a0+6]+105.0*A[a0+2])*B[b0+6]+105.0*A[a0+8]*B[b0+5]+(67.0820393249937*A[a0+4]+105.0*A[a0])*B[b0+4]+93.91485505499116*A[a0+3]*B[b0+3]+105.0*A[a0+6]*B[b0+2]+93.91485505499116*A[a0+1]*B[b0+1]+105.0*A[a0+4]*B[b0]); 
    tmp[5] = 0.004761904761904762*((67.0820393249937*A[a0+8]+105.0*A[a0+4])*B[b0+8]+(67.0820393249937*A[a0+7]+105.0*A[a0+1])*B[b0+7]+93.91485505499116*A[a0+6]*B[b0+6]+(67.0820393249937*A[a0+5]+105.0*A[a0])*B[b0+5]+105.0*A[a0+8]*B[b0+4]+93.91485505499116*A[a0+3]*B[b0+3]+93.91485505499116*A[a0+2]*B[b0+2]+105.0*A[a0+7]*B[b0+1]+105.0*A[a0+5]*B[b0]); 
    tmp[6] = 0.004761904761904762*((60.0*A[a0+6]+93.91485505499116*A[a0+2])*B[b0+8]+84.0*A[a0+3]*B[b0+7]+(60.0*A[a0+8]+93.91485505499116*A[a0+5]+67.0820393249937*A[a0+4]+105.0*A[a0])*B[b0+6]+93.91485505499116*A[a0+6]*B[b0+5]+(67.0820393249937*A[a0+6]+105.0*A[a0+2])*B[b0+4]+(84.0*A[a0+7]+93.91485505499116*A[a0+1])*B[b0+3]+(93.91485505499116*A[a0+8]+105.0*A[a0+4])*B[b0+2]+93.91485505499116*A[a0+3]*B[b0+1]+105.0*A[a0+6]*B[b0]); 
    tmp[7] = 0.004761904761904762*((60.0*A[a0+7]+93.91485505499116*A[a0+1])*B[b0+8]+(60.0*A[a0+8]+67.0820393249937*A[a0+5]+93.91485505499116*A[a0+4]+105.0*A[a0])*B[b0+7]+84.0*A[a0+3]*B[b0+6]+(67.0820393249937*A[a0+7]+105.0*A[a0+1])*B[b0+5]+93.91485505499116*A[a0+7]*B[b0+4]+(84.0*A[a0+6]+93.91485505499116*A[a0+2])*B[b0+3]+93.91485505499116*A[a0+3]*B[b0+2]+(93.91485505499116*A[a0+8]+105.0*A[a0+5])*B[b0+1]+105.0*A[a0+7]*B[b0]); 
    tmp[8] = 6.802721088435374e-4*((300.0*A[a0+8]+469.5742752749559*A[a0+5]+469.5742752749559*A[a0+4]+735.0*A[a0])*B[b0+8]+(420.0*A[a0+7]+657.4039853849381*A[a0+1])*B[b0+7]+(420.0*A[a0+6]+657.4039853849381*A[a0+2])*B[b0+6]+(469.5742752749559*A[a0+8]+735.0*A[a0+4])*B[b0+5]+(469.5742752749559*A[a0+8]+735.0*A[a0+5])*B[b0+4]+588.0*A[a0+3]*B[b0+3]+657.4039853849381*A[a0+6]*B[b0+2]+657.4039853849381*A[a0+7]*B[b0+1]+735.0*A[a0+8]*B[b0]); 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<9; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpDivide2xTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-1.5*A[3])+0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-1.5*A[3])-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.5*A[3]+0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
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
    // Fill BEV. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3]; 

    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*4,4,1) = data->u_S; 
  } 
}

void CartFieldBinOpDivide2xTensor_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (2.5*A[8]-1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (2.5*A[8]-1.936491673103709*A[7]+1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]+0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (2.5*A[8]+1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (2.5*A[8]+1.936491673103709*A[7]+1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]+0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
 
  unsigned short int b0; 
  double As[9]; 
  double Bs[9*Ncomp]; 
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
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 9*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = 0.0; 
      Bs[b0+2] = 0.0; 
      Bs[b0+3] = 0.0; 
      Bs[b0+4] = 0.0; 
      Bs[b0+5] = 0.0; 
      Bs[b0+6] = 0.0; 
      Bs[b0+7] = 0.0; 
      Bs[b0+8] = 0.0; 
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
    As[8] = A[8]; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 9*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = B[b0+1]; 
      Bs[b0+2] = B[b0+2]; 
      Bs[b0+3] = B[b0+3]; 
      Bs[b0+4] = B[b0+4]; 
      Bs[b0+5] = B[b0+5]; 
      Bs[b0+6] = B[b0+6]; 
      Bs[b0+7] = B[b0+7]; 
      Bs[b0+8] = B[b0+8]; 
    } 
  } 
 
  // Fill AEM matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(9,9); 
  data->AEM_S(0,0) = 0.5*As[0]; 
  data->AEM_S(0,1) = 0.5*As[1]; 
  data->AEM_S(0,2) = 0.5*As[2]; 
  data->AEM_S(0,3) = 0.5*As[3]; 
  data->AEM_S(0,4) = 0.5*As[4]; 
  data->AEM_S(0,5) = 0.5*As[5]; 
  data->AEM_S(0,6) = 0.5*As[6]; 
  data->AEM_S(0,7) = 0.5*As[7]; 
  data->AEM_S(0,8) = 0.5*As[8]; 
  data->AEM_S(1,0) = 0.5*As[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(1,2) = 0.5*As[3]; 
  data->AEM_S(1,3) = 0.447213595499958*As[6]+0.5*As[2]; 
  data->AEM_S(1,4) = 0.4472135954999579*As[1]; 
  data->AEM_S(1,5) = 0.5000000000000001*As[7]; 
  data->AEM_S(1,6) = 0.447213595499958*As[3]; 
  data->AEM_S(1,7) = 0.447213595499958*As[8]+0.5000000000000001*As[5]; 
  data->AEM_S(1,8) = 0.447213595499958*As[7]; 
  data->AEM_S(2,0) = 0.5*As[2]; 
  data->AEM_S(2,1) = 0.5*As[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*As[5]+0.5*As[0]; 
  data->AEM_S(2,3) = 0.447213595499958*As[7]+0.5*As[1]; 
  data->AEM_S(2,4) = 0.5000000000000001*As[6]; 
  data->AEM_S(2,5) = 0.4472135954999579*As[2]; 
  data->AEM_S(2,6) = 0.447213595499958*As[8]+0.5000000000000001*As[4]; 
  data->AEM_S(2,7) = 0.447213595499958*As[3]; 
  data->AEM_S(2,8) = 0.447213595499958*As[6]; 
  data->AEM_S(3,0) = 0.5*As[3]; 
  data->AEM_S(3,1) = 0.447213595499958*As[6]+0.5*As[2]; 
  data->AEM_S(3,2) = 0.447213595499958*As[7]+0.5*As[1]; 
  data->AEM_S(3,3) = 0.4*As[8]+0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(3,4) = 0.4472135954999579*As[3]; 
  data->AEM_S(3,5) = 0.4472135954999579*As[3]; 
  data->AEM_S(3,6) = 0.4*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(3,7) = 0.4*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(3,8) = 0.4*As[3]; 
  data->AEM_S(4,0) = 0.5*As[4]; 
  data->AEM_S(4,1) = 0.4472135954999579*As[1]; 
  data->AEM_S(4,2) = 0.5000000000000001*As[6]; 
  data->AEM_S(4,3) = 0.4472135954999579*As[3]; 
  data->AEM_S(4,4) = 0.31943828249997*As[4]+0.5*As[0]; 
  data->AEM_S(4,5) = 0.5*As[8]; 
  data->AEM_S(4,6) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  data->AEM_S(4,7) = 0.4472135954999579*As[7]; 
  data->AEM_S(4,8) = 0.31943828249997*As[8]+0.5*As[5]; 
  data->AEM_S(5,0) = 0.5*As[5]; 
  data->AEM_S(5,1) = 0.5000000000000001*As[7]; 
  data->AEM_S(5,2) = 0.4472135954999579*As[2]; 
  data->AEM_S(5,3) = 0.4472135954999579*As[3]; 
  data->AEM_S(5,4) = 0.5*As[8]; 
  data->AEM_S(5,5) = 0.31943828249997*As[5]+0.5*As[0]; 
  data->AEM_S(5,6) = 0.4472135954999579*As[6]; 
  data->AEM_S(5,7) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  data->AEM_S(5,8) = 0.31943828249997*As[8]+0.5*As[4]; 
  data->AEM_S(6,0) = 0.5*As[6]; 
  data->AEM_S(6,1) = 0.447213595499958*As[3]; 
  data->AEM_S(6,2) = 0.447213595499958*As[8]+0.5000000000000001*As[4]; 
  data->AEM_S(6,3) = 0.4*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(6,4) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  data->AEM_S(6,5) = 0.4472135954999579*As[6]; 
  data->AEM_S(6,6) = 0.2857142857142857*As[8]+0.4472135954999579*As[5]+0.31943828249997*As[4]+0.5*As[0]; 
  data->AEM_S(6,7) = 0.4*As[3]; 
  data->AEM_S(6,8) = 0.2857142857142857*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(7,0) = 0.5*As[7]; 
  data->AEM_S(7,1) = 0.447213595499958*As[8]+0.5000000000000001*As[5]; 
  data->AEM_S(7,2) = 0.447213595499958*As[3]; 
  data->AEM_S(7,3) = 0.4*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(7,4) = 0.4472135954999579*As[7]; 
  data->AEM_S(7,5) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  data->AEM_S(7,6) = 0.4*As[3]; 
  data->AEM_S(7,7) = 0.2857142857142857*As[8]+0.31943828249997*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(7,8) = 0.2857142857142857*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(8,0) = 0.5*As[8]; 
  data->AEM_S(8,1) = 0.447213595499958*As[7]; 
  data->AEM_S(8,2) = 0.447213595499958*As[6]; 
  data->AEM_S(8,3) = 0.4*As[3]; 
  data->AEM_S(8,4) = 0.31943828249997*As[8]+0.5*As[5]; 
  data->AEM_S(8,5) = 0.31943828249997*As[8]+0.5*As[4]; 
  data->AEM_S(8,6) = 0.2857142857142857*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(8,7) = 0.2857142857142857*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(8,8) = 0.2040816326530612*As[8]+0.31943828249997*As[5]+0.31943828249997*As[4]+0.5*As[0]; 

  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 9*vd; 
    // Fill BEV. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3],Bs[b0+4],Bs[b0+5],Bs[b0+6],Bs[b0+7],Bs[b0+8]; 

    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*9,9,1) = data->u_S; 
  } 
}

void CartFieldBinOpDotProduct2xTensor_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
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
 
void CartFieldBinOpDotProduct2xTensor_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<9; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 9*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.5*A[a0+8]*B[a0+8]+0.5*A[a0+7]*B[a0+7]+0.5*A[a0+6]*B[a0+6]+0.5*A[a0+5]*B[a0+5]+0.5*A[a0+4]*B[a0+4]+0.5*A[a0+3]*B[a0+3]+0.5*A[a0+2]*B[a0+2]+0.5*A[a0+1]*B[a0+1]+0.5*A[a0]*B[a0]; 
    out[1] += 0.447213595499958*A[a0+7]*B[a0+8]+0.447213595499958*B[a0+7]*A[a0+8]+0.5000000000000001*A[a0+5]*B[a0+7]+0.5000000000000001*B[a0+5]*A[a0+7]+0.447213595499958*A[a0+3]*B[a0+6]+0.447213595499958*B[a0+3]*A[a0+6]+0.4472135954999579*A[a0+1]*B[a0+4]+0.4472135954999579*B[a0+1]*A[a0+4]+0.5*A[a0+2]*B[a0+3]+0.5*B[a0+2]*A[a0+3]+0.5*A[a0]*B[a0+1]+0.5*B[a0]*A[a0+1]; 
    out[2] += 0.447213595499958*A[a0+6]*B[a0+8]+0.447213595499958*B[a0+6]*A[a0+8]+0.447213595499958*A[a0+3]*B[a0+7]+0.447213595499958*B[a0+3]*A[a0+7]+0.5000000000000001*A[a0+4]*B[a0+6]+0.5000000000000001*B[a0+4]*A[a0+6]+0.4472135954999579*A[a0+2]*B[a0+5]+0.4472135954999579*B[a0+2]*A[a0+5]+0.5*A[a0+1]*B[a0+3]+0.5*B[a0+1]*A[a0+3]+0.5*A[a0]*B[a0+2]+0.5*B[a0]*A[a0+2]; 
    out[3] += 0.4*A[a0+3]*B[a0+8]+0.4*B[a0+3]*A[a0+8]+0.4*A[a0+6]*B[a0+7]+0.447213595499958*A[a0+2]*B[a0+7]+0.4*B[a0+6]*A[a0+7]+0.447213595499958*B[a0+2]*A[a0+7]+0.447213595499958*A[a0+1]*B[a0+6]+0.447213595499958*B[a0+1]*A[a0+6]+0.4472135954999579*A[a0+3]*B[a0+5]+0.4472135954999579*B[a0+3]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+4]+0.4472135954999579*B[a0+3]*A[a0+4]+0.5*A[a0]*B[a0+3]+0.5*B[a0]*A[a0+3]+0.5*A[a0+1]*B[a0+2]+0.5*B[a0+1]*A[a0+2]; 
    out[4] += 0.31943828249997*A[a0+8]*B[a0+8]+0.5*A[a0+5]*B[a0+8]+0.5*B[a0+5]*A[a0+8]+0.4472135954999579*A[a0+7]*B[a0+7]+0.31943828249997*A[a0+6]*B[a0+6]+0.5000000000000001*A[a0+2]*B[a0+6]+0.5000000000000001*B[a0+2]*A[a0+6]+0.31943828249997*A[a0+4]*B[a0+4]+0.5*A[a0]*B[a0+4]+0.5*B[a0]*A[a0+4]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+1]*B[a0+1]; 
    out[5] += 0.31943828249997*A[a0+8]*B[a0+8]+0.5*A[a0+4]*B[a0+8]+0.5*B[a0+4]*A[a0+8]+0.31943828249997*A[a0+7]*B[a0+7]+0.5000000000000001*A[a0+1]*B[a0+7]+0.5000000000000001*B[a0+1]*A[a0+7]+0.4472135954999579*A[a0+6]*B[a0+6]+0.31943828249997*A[a0+5]*B[a0+5]+0.5*A[a0]*B[a0+5]+0.5*B[a0]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+2]*B[a0+2]; 
    out[6] += 0.2857142857142857*A[a0+6]*B[a0+8]+0.447213595499958*A[a0+2]*B[a0+8]+0.2857142857142857*B[a0+6]*A[a0+8]+0.447213595499958*B[a0+2]*A[a0+8]+0.4*A[a0+3]*B[a0+7]+0.4*B[a0+3]*A[a0+7]+0.4472135954999579*A[a0+5]*B[a0+6]+0.31943828249997*A[a0+4]*B[a0+6]+0.5*A[a0]*B[a0+6]+0.4472135954999579*B[a0+5]*A[a0+6]+0.31943828249997*B[a0+4]*A[a0+6]+0.5*B[a0]*A[a0+6]+0.5000000000000001*A[a0+2]*B[a0+4]+0.5000000000000001*B[a0+2]*A[a0+4]+0.447213595499958*A[a0+1]*B[a0+3]+0.447213595499958*B[a0+1]*A[a0+3]; 
    out[7] += 0.2857142857142857*A[a0+7]*B[a0+8]+0.447213595499958*A[a0+1]*B[a0+8]+0.2857142857142857*B[a0+7]*A[a0+8]+0.447213595499958*B[a0+1]*A[a0+8]+0.31943828249997*A[a0+5]*B[a0+7]+0.4472135954999579*A[a0+4]*B[a0+7]+0.5*A[a0]*B[a0+7]+0.31943828249997*B[a0+5]*A[a0+7]+0.4472135954999579*B[a0+4]*A[a0+7]+0.5*B[a0]*A[a0+7]+0.4*A[a0+3]*B[a0+6]+0.4*B[a0+3]*A[a0+6]+0.5000000000000001*A[a0+1]*B[a0+5]+0.5000000000000001*B[a0+1]*A[a0+5]+0.447213595499958*A[a0+2]*B[a0+3]+0.447213595499958*B[a0+2]*A[a0+3]; 
    out[8] += 0.2040816326530612*A[a0+8]*B[a0+8]+0.31943828249997*A[a0+5]*B[a0+8]+0.31943828249997*A[a0+4]*B[a0+8]+0.5*A[a0]*B[a0+8]+0.31943828249997*B[a0+5]*A[a0+8]+0.31943828249997*B[a0+4]*A[a0+8]+0.5*B[a0]*A[a0+8]+0.2857142857142857*A[a0+7]*B[a0+7]+0.447213595499958*A[a0+1]*B[a0+7]+0.447213595499958*B[a0+1]*A[a0+7]+0.2857142857142857*A[a0+6]*B[a0+6]+0.447213595499958*A[a0+2]*B[a0+6]+0.447213595499958*B[a0+2]*A[a0+6]+0.5*A[a0+4]*B[a0+5]+0.5*B[a0+4]*A[a0+5]+0.4*A[a0+3]*B[a0+3]; 
  } 
 
} 
 
