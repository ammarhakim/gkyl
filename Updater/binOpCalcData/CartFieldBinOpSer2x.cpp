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
    tmp[0] = 0.5*(A[a0+7]*B[b0+7]+A[a0+6]*B[b0+6]+A[a0+5]*B[b0+5]+A[a0+4]*B[b0+4]+A[a0+3]*B[b0+3]+A[a0+2]*B[b0+2]+A[a0+1]*B[b0+1]+A[a0]*B[b0]); 
    tmp[1] = 0.03333333333333333*(15.0*A[a0+5]*B[b0+7]+13.41640786499874*A[a0+3]*B[b0+6]+15.0*A[a0+7]*B[b0+5]+13.41640786499874*A[a0+1]*B[b0+4]+(13.41640786499874*A[a0+6]+15.0*A[a0+2])*B[b0+3]+15.0*A[a0+3]*B[b0+2]+(13.41640786499874*A[a0+4]+15.0*A[a0])*B[b0+1]+15.0*A[a0+1]*B[b0]); 
    tmp[2] = 0.03333333333333333*(13.41640786499874*A[a0+3]*B[b0+7]+15.0*A[a0+4]*B[b0+6]+13.41640786499874*A[a0+2]*B[b0+5]+15.0*A[a0+6]*B[b0+4]+(13.41640786499874*A[a0+7]+15.0*A[a0+1])*B[b0+3]+(13.41640786499874*A[a0+5]+15.0*A[a0])*B[b0+2]+15.0*A[a0+3]*B[b0+1]+15.0*A[a0+2]*B[b0]); 
    tmp[3] = 0.03333333333333333*((12.0*A[a0+6]+13.41640786499874*A[a0+2])*B[b0+7]+(12.0*A[a0+7]+13.41640786499874*A[a0+1])*B[b0+6]+13.41640786499874*A[a0+3]*B[b0+5]+13.41640786499874*A[a0+3]*B[b0+4]+(13.41640786499874*A[a0+5]+13.41640786499874*A[a0+4]+15.0*A[a0])*B[b0+3]+(13.41640786499874*A[a0+7]+15.0*A[a0+1])*B[b0+2]+(13.41640786499874*A[a0+6]+15.0*A[a0+2])*B[b0+1]+15.0*A[a0+3]*B[b0]); 
    tmp[4] = 0.004761904761904762*(93.91485505499116*A[a0+7]*B[b0+7]+(67.0820393249937*A[a0+6]+105.0*A[a0+2])*B[b0+6]+(67.0820393249937*A[a0+4]+105.0*A[a0])*B[b0+4]+93.91485505499116*A[a0+3]*B[b0+3]+105.0*A[a0+6]*B[b0+2]+93.91485505499116*A[a0+1]*B[b0+1]+105.0*A[a0+4]*B[b0]); 
    tmp[5] = 0.004761904761904762*((67.0820393249937*A[a0+7]+105.0*A[a0+1])*B[b0+7]+93.91485505499116*A[a0+6]*B[b0+6]+(67.0820393249937*A[a0+5]+105.0*A[a0])*B[b0+5]+93.91485505499116*A[a0+3]*B[b0+3]+93.91485505499116*A[a0+2]*B[b0+2]+105.0*A[a0+7]*B[b0+1]+105.0*A[a0+5]*B[b0]); 
    tmp[6] = 0.004761904761904762*(84.0*A[a0+3]*B[b0+7]+(93.91485505499116*A[a0+5]+67.0820393249937*A[a0+4]+105.0*A[a0])*B[b0+6]+93.91485505499116*A[a0+6]*B[b0+5]+(67.0820393249937*A[a0+6]+105.0*A[a0+2])*B[b0+4]+(84.0*A[a0+7]+93.91485505499116*A[a0+1])*B[b0+3]+105.0*A[a0+4]*B[b0+2]+93.91485505499116*A[a0+3]*B[b0+1]+105.0*A[a0+6]*B[b0]); 
    tmp[7] = 0.004761904761904762*((67.0820393249937*A[a0+5]+93.91485505499116*A[a0+4]+105.0*A[a0])*B[b0+7]+84.0*A[a0+3]*B[b0+6]+(67.0820393249937*A[a0+7]+105.0*A[a0+1])*B[b0+5]+93.91485505499116*A[a0+7]*B[b0+4]+(84.0*A[a0+6]+93.91485505499116*A[a0+2])*B[b0+3]+93.91485505499116*A[a0+3]*B[b0+2]+105.0*A[a0+5]*B[b0+1]+105.0*A[a0+7]*B[b0]); 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<8; i++) 
    { 
      out[b0+i] = tmp[i]; 
    } 
  } 
 
} 
 
void CartFieldBinOpMultiply2xSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[12]; 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int b0 = 12*vd; 
    unsigned short int a0 = b0*eqNcomp; 
    // Component-wise (of the vectors) multiplication. 
    tmp[0] = 0.5*(A[a0+11]*B[b0+11]+A[a0+10]*B[b0+10]+A[a0+9]*B[b0+9]+A[a0+8]*B[b0+8]+A[a0+7]*B[b0+7]+A[a0+6]*B[b0+6]+A[a0+5]*B[b0+5]+A[a0+4]*B[b0+4]+A[a0+3]*B[b0+3]+A[a0+2]*B[b0+2]+A[a0+1]*B[b0+1]+A[a0]*B[b0]); 
    tmp[1] = 0.004761904761904762*(105.0*A[a0+9]*B[b0+11]+92.22255689363637*A[a0+6]*B[b0+10]+105.0*A[a0+11]*B[b0+9]+92.22255689363637*A[a0+4]*B[b0+8]+105.0*A[a0+5]*B[b0+7]+(92.22255689363637*A[a0+10]+93.91485505499116*A[a0+3])*B[b0+6]+105.0*A[a0+7]*B[b0+5]+(92.22255689363637*A[a0+8]+93.91485505499116*A[a0+1])*B[b0+4]+(93.91485505499116*A[a0+6]+105.0*A[a0+2])*B[b0+3]+105.0*A[a0+3]*B[b0+2]+(93.91485505499116*A[a0+4]+105.0*A[a0])*B[b0+1]+105.0*A[a0+1]*B[b0]); 
    tmp[2] = 0.004761904761904762*(92.22255689363637*A[a0+7]*B[b0+11]+105.0*A[a0+8]*B[b0+10]+92.22255689363637*A[a0+5]*B[b0+9]+105.0*A[a0+10]*B[b0+8]+(92.22255689363637*A[a0+11]+93.91485505499116*A[a0+3])*B[b0+7]+105.0*A[a0+4]*B[b0+6]+(92.22255689363637*A[a0+9]+93.91485505499116*A[a0+2])*B[b0+5]+105.0*A[a0+6]*B[b0+4]+(93.91485505499116*A[a0+7]+105.0*A[a0+1])*B[b0+3]+(93.91485505499116*A[a0+5]+105.0*A[a0])*B[b0+2]+105.0*A[a0+3]*B[b0+1]+105.0*A[a0+2]*B[b0]); 
    tmp[3] = 0.004761904761904762*(92.22255689363638*A[a0+5]*B[b0+11]+92.22255689363638*A[a0+4]*B[b0+10]+92.2225568936364*A[a0+7]*B[b0+9]+92.2225568936364*A[a0+6]*B[b0+8]+(92.2225568936364*A[a0+9]+84.0*A[a0+6]+93.91485505499116*A[a0+2])*B[b0+7]+(92.2225568936364*A[a0+8]+84.0*A[a0+7]+93.91485505499116*A[a0+1])*B[b0+6]+(92.22255689363638*A[a0+11]+93.91485505499116*A[a0+3])*B[b0+5]+(92.22255689363638*A[a0+10]+93.91485505499116*A[a0+3])*B[b0+4]+(93.91485505499116*A[a0+5]+93.91485505499116*A[a0+4]+105.0*A[a0])*B[b0+3]+(93.91485505499116*A[a0+7]+105.0*A[a0+1])*B[b0+2]+(93.91485505499116*A[a0+6]+105.0*A[a0+2])*B[b0+1]+105.0*A[a0+3]*B[b0]); 
    tmp[4] = 0.004761904761904762*(93.91485505499116*A[a0+11]*B[b0+11]+(62.60990336999411*A[a0+10]+92.22255689363638*A[a0+3])*B[b0+10]+(62.60990336999411*A[a0+8]+92.22255689363637*A[a0+1])*B[b0+8]+93.91485505499116*A[a0+7]*B[b0+7]+(67.0820393249937*A[a0+6]+105.0*A[a0+2])*B[b0+6]+(67.0820393249937*A[a0+4]+105.0*A[a0])*B[b0+4]+(92.22255689363638*A[a0+10]+93.91485505499116*A[a0+3])*B[b0+3]+105.0*A[a0+6]*B[b0+2]+(92.22255689363637*A[a0+8]+93.91485505499116*A[a0+1])*B[b0+1]+105.0*A[a0+4]*B[b0]); 
    tmp[5] = 0.004761904761904762*((62.60990336999411*A[a0+11]+92.22255689363638*A[a0+3])*B[b0+11]+93.91485505499116*A[a0+10]*B[b0+10]+(62.60990336999411*A[a0+9]+92.22255689363637*A[a0+2])*B[b0+9]+(67.0820393249937*A[a0+7]+105.0*A[a0+1])*B[b0+7]+93.91485505499116*A[a0+6]*B[b0+6]+(67.0820393249937*A[a0+5]+105.0*A[a0])*B[b0+5]+(92.22255689363638*A[a0+11]+93.91485505499116*A[a0+3])*B[b0+3]+(92.22255689363637*A[a0+9]+93.91485505499116*A[a0+2])*B[b0+2]+105.0*A[a0+7]*B[b0+1]+105.0*A[a0+5]*B[b0]); 
    tmp[6] = 0.001587301587301587*(247.4590875276153*A[a0+7]*B[b0+11]+(187.8297101099824*A[a0+8]+247.4590875276153*A[a0+7]+276.6676706809091*A[a0+1])*B[b0+10]+(187.8297101099824*A[a0+10]+276.6676706809092*A[a0+3])*B[b0+8]+(247.4590875276153*A[a0+11]+247.4590875276153*A[a0+10]+252.0*A[a0+3])*B[b0+7]+(281.7445651649735*A[a0+5]+201.2461179749811*A[a0+4]+315.0*A[a0])*B[b0+6]+281.7445651649735*A[a0+6]*B[b0+5]+(201.2461179749811*A[a0+6]+314.9999999999999*A[a0+2])*B[b0+4]+(276.6676706809092*A[a0+8]+252.0*A[a0+7]+281.7445651649734*A[a0+1])*B[b0+3]+314.9999999999999*A[a0+4]*B[b0+2]+(276.6676706809091*A[a0+10]+281.7445651649734*A[a0+3])*B[b0+1]+315.0*A[a0+6]*B[b0]); 
    tmp[7] = 0.001587301587301587*((187.8297101099824*A[a0+9]+247.4590875276153*A[a0+6]+276.6676706809091*A[a0+2])*B[b0+11]+247.4590875276153*A[a0+6]*B[b0+10]+(187.8297101099824*A[a0+11]+276.6676706809092*A[a0+3])*B[b0+9]+(201.2461179749811*A[a0+5]+281.7445651649735*A[a0+4]+315.0*A[a0])*B[b0+7]+(247.4590875276153*A[a0+11]+247.4590875276153*A[a0+10]+252.0*A[a0+3])*B[b0+6]+(201.2461179749811*A[a0+7]+314.9999999999999*A[a0+1])*B[b0+5]+281.7445651649735*A[a0+7]*B[b0+4]+(276.6676706809092*A[a0+9]+252.0*A[a0+6]+281.7445651649734*A[a0+2])*B[b0+3]+(276.6676706809091*A[a0+11]+281.7445651649734*A[a0+3])*B[b0+2]+314.9999999999999*A[a0+5]*B[b0+1]+315.0*A[a0+7]*B[b0]); 
    tmp[8] = 0.001587301587301587*((187.8297101099824*A[a0+6]+315.0*A[a0+2])*B[b0+10]+(187.8297101099823*A[a0+4]+315.0*A[a0])*B[b0+8]+(187.8297101099824*A[a0+10]+276.6676706809092*A[a0+3])*B[b0+6]+(187.8297101099823*A[a0+8]+276.6676706809091*A[a0+1])*B[b0+4]+276.6676706809092*A[a0+6]*B[b0+3]+315.0*A[a0+10]*B[b0+2]+276.6676706809091*A[a0+4]*B[b0+1]+315.0*A[a0+8]*B[b0]); 
    tmp[9] = 0.001587301587301587*((187.8297101099824*A[a0+7]+315.0*A[a0+1])*B[b0+11]+(187.8297101099823*A[a0+5]+315.0*A[a0])*B[b0+9]+(187.8297101099824*A[a0+11]+276.6676706809092*A[a0+3])*B[b0+7]+(187.8297101099823*A[a0+9]+276.6676706809091*A[a0+2])*B[b0+5]+276.6676706809092*A[a0+7]*B[b0+3]+276.6676706809091*A[a0+5]*B[b0+2]+315.0*A[a0+11]*B[b0+1]+315.0*A[a0+9]*B[b0]); 
    tmp[10] = 0.001587301587301587*((281.7445651649735*A[a0+5]+187.8297101099823*A[a0+4]+315.0*A[a0])*B[b0+10]+(187.8297101099824*A[a0+6]+315.0*A[a0+2])*B[b0+8]+247.4590875276153*A[a0+6]*B[b0+7]+(187.8297101099824*A[a0+8]+247.4590875276153*A[a0+7]+276.6676706809091*A[a0+1])*B[b0+6]+281.7445651649735*A[a0+10]*B[b0+5]+(187.8297101099823*A[a0+10]+276.6676706809092*A[a0+3])*B[b0+4]+276.6676706809092*A[a0+4]*B[b0+3]+315.0*A[a0+8]*B[b0+2]+276.6676706809091*A[a0+6]*B[b0+1]+315.0*A[a0+10]*B[b0]); 
    tmp[11] = 0.001587301587301587*((187.8297101099823*A[a0+5]+281.7445651649735*A[a0+4]+315.0*A[a0])*B[b0+11]+(187.8297101099824*A[a0+7]+315.0*A[a0+1])*B[b0+9]+(187.8297101099824*A[a0+9]+247.4590875276153*A[a0+6]+276.6676706809091*A[a0+2])*B[b0+7]+247.4590875276153*A[a0+7]*B[b0+6]+(187.8297101099823*A[a0+11]+276.6676706809092*A[a0+3])*B[b0+5]+281.7445651649735*A[a0+11]*B[b0+4]+276.6676706809092*A[a0+5]*B[b0+3]+276.6676706809091*A[a0+7]*B[b0+2]+315.0*A[a0+9]*B[b0+1]+315.0*A[a0+11]*B[b0]); 
 
    // This tmp allows for in-place multiplication. 
    for (unsigned short int i=0; i<12; i++) 
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

void CartFieldBinOpDivide2xSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.936491673103709*A[7])-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-1.936491673103709*A[7])+1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]+0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (1.936491673103709*A[7]+1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]+0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
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
 
  // Fill AEM matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(8,8); 
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
    // Fill BEV. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3],Bs[b0+4],Bs[b0+5],Bs[b0+6],Bs[b0+7]; 

    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*8,8,1) = data->u_S; 
  } 
}

void CartFieldBinOpDivide2xSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (2.29128784747792*A[11]+2.29128784747792*A[10]-1.322875655532295*A[9]-1.322875655532295*A[8]-1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-2.29128784747792*A[11])-2.29128784747792*A[10]+1.322875655532295*A[9]-1.322875655532295*A[8]-1.936491673103709*A[7]+1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]+0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-2.29128784747792*A[11])-2.29128784747792*A[10]-1.322875655532295*A[9]+1.322875655532295*A[8]+1.936491673103709*A[7]-1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]-1.5*A[3]-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (2.29128784747792*A[11]+2.29128784747792*A[10]+1.322875655532295*A[9]+1.322875655532295*A[8]+1.936491673103709*A[7]+1.936491673103709*A[6]+1.118033988749895*A[5]+1.118033988749895*A[4]+1.5*A[3]+0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
 
  unsigned short int b0; 
  double As[12]; 
  double Bs[12*Ncomp]; 
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
    As[10] = 0.0; 
    As[11] = 0.0; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 12*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = 0.0; 
      Bs[b0+2] = 0.0; 
      Bs[b0+3] = 0.0; 
      Bs[b0+4] = 0.0; 
      Bs[b0+5] = 0.0; 
      Bs[b0+6] = 0.0; 
      Bs[b0+7] = 0.0; 
      Bs[b0+8] = 0.0; 
      Bs[b0+9] = 0.0; 
      Bs[b0+10] = 0.0; 
      Bs[b0+11] = 0.0; 
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
    As[9] = A[9]; 
    As[10] = A[10]; 
    As[11] = A[11]; 
    for(unsigned short int vd=0; vd<Ncomp; vd++) 
    { 
      b0 = 12*vd; 
      Bs[b0] = B[b0]; 
      Bs[b0+1] = B[b0+1]; 
      Bs[b0+2] = B[b0+2]; 
      Bs[b0+3] = B[b0+3]; 
      Bs[b0+4] = B[b0+4]; 
      Bs[b0+5] = B[b0+5]; 
      Bs[b0+6] = B[b0+6]; 
      Bs[b0+7] = B[b0+7]; 
      Bs[b0+8] = B[b0+8]; 
      Bs[b0+9] = B[b0+9]; 
      Bs[b0+10] = B[b0+10]; 
      Bs[b0+11] = B[b0+11]; 
    } 
  } 
 
  // Fill AEM matrix. 
  data->AEM_S = Eigen::MatrixXd::Zero(12,12); 
  data->AEM_S(0,0) = 0.5*As[0]; 
  data->AEM_S(0,1) = 0.5*As[1]; 
  data->AEM_S(0,2) = 0.5*As[2]; 
  data->AEM_S(0,3) = 0.5*As[3]; 
  data->AEM_S(0,4) = 0.5*As[4]; 
  data->AEM_S(0,5) = 0.5*As[5]; 
  data->AEM_S(0,6) = 0.5*As[6]; 
  data->AEM_S(0,7) = 0.5*As[7]; 
  data->AEM_S(0,8) = 0.5*As[8]; 
  data->AEM_S(0,9) = 0.5*As[9]; 
  data->AEM_S(0,10) = 0.5*As[10]; 
  data->AEM_S(0,11) = 0.5*As[11]; 
  data->AEM_S(1,0) = 0.5*As[1]; 
  data->AEM_S(1,1) = 0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(1,2) = 0.5*As[3]; 
  data->AEM_S(1,3) = 0.447213595499958*As[6]+0.5*As[2]; 
  data->AEM_S(1,4) = 0.4391550328268398*As[8]+0.4472135954999579*As[1]; 
  data->AEM_S(1,5) = 0.5000000000000001*As[7]; 
  data->AEM_S(1,6) = 0.4391550328268399*As[10]+0.447213595499958*As[3]; 
  data->AEM_S(1,7) = 0.5000000000000001*As[5]; 
  data->AEM_S(1,8) = 0.4391550328268398*As[4]; 
  data->AEM_S(1,9) = 0.5*As[11]; 
  data->AEM_S(1,10) = 0.4391550328268399*As[6]; 
  data->AEM_S(1,11) = 0.5*As[9]; 
  data->AEM_S(2,0) = 0.5*As[2]; 
  data->AEM_S(2,1) = 0.5*As[3]; 
  data->AEM_S(2,2) = 0.4472135954999579*As[5]+0.5*As[0]; 
  data->AEM_S(2,3) = 0.447213595499958*As[7]+0.5*As[1]; 
  data->AEM_S(2,4) = 0.5000000000000001*As[6]; 
  data->AEM_S(2,5) = 0.4391550328268398*As[9]+0.4472135954999579*As[2]; 
  data->AEM_S(2,6) = 0.5000000000000001*As[4]; 
  data->AEM_S(2,7) = 0.4391550328268399*As[11]+0.447213595499958*As[3]; 
  data->AEM_S(2,8) = 0.5*As[10]; 
  data->AEM_S(2,9) = 0.4391550328268398*As[5]; 
  data->AEM_S(2,10) = 0.5*As[8]; 
  data->AEM_S(2,11) = 0.4391550328268399*As[7]; 
  data->AEM_S(3,0) = 0.5*As[3]; 
  data->AEM_S(3,1) = 0.447213595499958*As[6]+0.5*As[2]; 
  data->AEM_S(3,2) = 0.447213595499958*As[7]+0.5*As[1]; 
  data->AEM_S(3,3) = 0.4472135954999579*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(3,4) = 0.4391550328268399*As[10]+0.4472135954999579*As[3]; 
  data->AEM_S(3,5) = 0.4391550328268399*As[11]+0.4472135954999579*As[3]; 
  data->AEM_S(3,6) = 0.4391550328268399*As[8]+0.4*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(3,7) = 0.4391550328268399*As[9]+0.4*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(3,8) = 0.4391550328268399*As[6]; 
  data->AEM_S(3,9) = 0.4391550328268399*As[7]; 
  data->AEM_S(3,10) = 0.4391550328268399*As[4]; 
  data->AEM_S(3,11) = 0.4391550328268399*As[5]; 
  data->AEM_S(4,0) = 0.5*As[4]; 
  data->AEM_S(4,1) = 0.4391550328268398*As[8]+0.4472135954999579*As[1]; 
  data->AEM_S(4,2) = 0.5000000000000001*As[6]; 
  data->AEM_S(4,3) = 0.4391550328268399*As[10]+0.4472135954999579*As[3]; 
  data->AEM_S(4,4) = 0.31943828249997*As[4]+0.5*As[0]; 
  data->AEM_S(4,6) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  data->AEM_S(4,7) = 0.4472135954999579*As[7]; 
  data->AEM_S(4,8) = 0.2981423969999719*As[8]+0.4391550328268398*As[1]; 
  data->AEM_S(4,10) = 0.2981423969999719*As[10]+0.4391550328268399*As[3]; 
  data->AEM_S(4,11) = 0.4472135954999579*As[11]; 
  data->AEM_S(5,0) = 0.5*As[5]; 
  data->AEM_S(5,1) = 0.5000000000000001*As[7]; 
  data->AEM_S(5,2) = 0.4391550328268398*As[9]+0.4472135954999579*As[2]; 
  data->AEM_S(5,3) = 0.4391550328268399*As[11]+0.4472135954999579*As[3]; 
  data->AEM_S(5,5) = 0.31943828249997*As[5]+0.5*As[0]; 
  data->AEM_S(5,6) = 0.4472135954999579*As[6]; 
  data->AEM_S(5,7) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  data->AEM_S(5,9) = 0.2981423969999719*As[9]+0.4391550328268398*As[2]; 
  data->AEM_S(5,10) = 0.4472135954999579*As[10]; 
  data->AEM_S(5,11) = 0.2981423969999719*As[11]+0.4391550328268399*As[3]; 
  data->AEM_S(6,0) = 0.5*As[6]; 
  data->AEM_S(6,1) = 0.4391550328268399*As[10]+0.447213595499958*As[3]; 
  data->AEM_S(6,2) = 0.5000000000000001*As[4]; 
  data->AEM_S(6,3) = 0.4391550328268399*As[8]+0.4*As[7]+0.447213595499958*As[1]; 
  data->AEM_S(6,4) = 0.31943828249997*As[6]+0.5000000000000001*As[2]; 
  data->AEM_S(6,5) = 0.4472135954999579*As[6]; 
  data->AEM_S(6,6) = 0.4472135954999579*As[5]+0.31943828249997*As[4]+0.5*As[0]; 
  data->AEM_S(6,7) = 0.3927922024247863*As[11]+0.3927922024247863*As[10]+0.4*As[3]; 
  data->AEM_S(6,8) = 0.2981423969999719*As[10]+0.4391550328268399*As[3]; 
  data->AEM_S(6,10) = 0.2981423969999719*As[8]+0.3927922024247863*As[7]+0.4391550328268399*As[1]; 
  data->AEM_S(6,11) = 0.3927922024247863*As[7]; 
  data->AEM_S(7,0) = 0.5*As[7]; 
  data->AEM_S(7,1) = 0.5000000000000001*As[5]; 
  data->AEM_S(7,2) = 0.4391550328268399*As[11]+0.447213595499958*As[3]; 
  data->AEM_S(7,3) = 0.4391550328268399*As[9]+0.4*As[6]+0.447213595499958*As[2]; 
  data->AEM_S(7,4) = 0.4472135954999579*As[7]; 
  data->AEM_S(7,5) = 0.31943828249997*As[7]+0.5000000000000001*As[1]; 
  data->AEM_S(7,6) = 0.3927922024247863*As[11]+0.3927922024247863*As[10]+0.4*As[3]; 
  data->AEM_S(7,7) = 0.31943828249997*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 
  data->AEM_S(7,9) = 0.2981423969999719*As[11]+0.4391550328268399*As[3]; 
  data->AEM_S(7,10) = 0.3927922024247863*As[6]; 
  data->AEM_S(7,11) = 0.2981423969999719*As[9]+0.3927922024247863*As[6]+0.4391550328268399*As[2]; 
  data->AEM_S(8,0) = 0.5*As[8]; 
  data->AEM_S(8,1) = 0.4391550328268398*As[4]; 
  data->AEM_S(8,2) = 0.5*As[10]; 
  data->AEM_S(8,3) = 0.4391550328268399*As[6]; 
  data->AEM_S(8,4) = 0.2981423969999719*As[8]+0.4391550328268398*As[1]; 
  data->AEM_S(8,6) = 0.2981423969999719*As[10]+0.4391550328268399*As[3]; 
  data->AEM_S(8,8) = 0.2981423969999719*As[4]+0.5*As[0]; 
  data->AEM_S(8,10) = 0.2981423969999719*As[6]+0.5*As[2]; 
  data->AEM_S(9,0) = 0.5*As[9]; 
  data->AEM_S(9,1) = 0.5*As[11]; 
  data->AEM_S(9,2) = 0.4391550328268398*As[5]; 
  data->AEM_S(9,3) = 0.4391550328268399*As[7]; 
  data->AEM_S(9,5) = 0.2981423969999719*As[9]+0.4391550328268398*As[2]; 
  data->AEM_S(9,7) = 0.2981423969999719*As[11]+0.4391550328268399*As[3]; 
  data->AEM_S(9,9) = 0.2981423969999719*As[5]+0.5*As[0]; 
  data->AEM_S(9,11) = 0.2981423969999719*As[7]+0.5*As[1]; 
  data->AEM_S(10,0) = 0.5*As[10]; 
  data->AEM_S(10,1) = 0.4391550328268399*As[6]; 
  data->AEM_S(10,2) = 0.5*As[8]; 
  data->AEM_S(10,3) = 0.4391550328268399*As[4]; 
  data->AEM_S(10,4) = 0.2981423969999719*As[10]+0.4391550328268399*As[3]; 
  data->AEM_S(10,5) = 0.4472135954999579*As[10]; 
  data->AEM_S(10,6) = 0.2981423969999719*As[8]+0.3927922024247863*As[7]+0.4391550328268399*As[1]; 
  data->AEM_S(10,7) = 0.3927922024247863*As[6]; 
  data->AEM_S(10,8) = 0.2981423969999719*As[6]+0.5*As[2]; 
  data->AEM_S(10,10) = 0.4472135954999579*As[5]+0.2981423969999719*As[4]+0.5*As[0]; 
  data->AEM_S(11,0) = 0.5*As[11]; 
  data->AEM_S(11,1) = 0.5*As[9]; 
  data->AEM_S(11,2) = 0.4391550328268399*As[7]; 
  data->AEM_S(11,3) = 0.4391550328268399*As[5]; 
  data->AEM_S(11,4) = 0.4472135954999579*As[11]; 
  data->AEM_S(11,5) = 0.2981423969999719*As[11]+0.4391550328268399*As[3]; 
  data->AEM_S(11,6) = 0.3927922024247863*As[7]; 
  data->AEM_S(11,7) = 0.2981423969999719*As[9]+0.3927922024247863*As[6]+0.4391550328268399*As[2]; 
  data->AEM_S(11,9) = 0.2981423969999719*As[7]+0.5*As[1]; 
  data->AEM_S(11,11) = 0.2981423969999719*As[5]+0.4472135954999579*As[4]+0.5*As[0]; 

  for(unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    b0 = 12*vd; 
    // Fill BEV. 
    data->BEV_S << Bs[b0],Bs[b0+1],Bs[b0+2],Bs[b0+3],Bs[b0+4],Bs[b0+5],Bs[b0+6],Bs[b0+7],Bs[b0+8],Bs[b0+9],Bs[b0+10],Bs[b0+11]; 

    // Solve the system of equations. 
    data->u_S = data->AEM_S.colPivHouseholderQr().solve(data->BEV_S); 
 
    // Copy data from Eigen vector. 
    Eigen::Map<VectorXd>(out+vd*12,12,1) = data->u_S; 
  } 
}

void CartFieldBinOpDividePositivity2xSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       denominator field (must be a scalar field). 
  // B:       numerator field (can be scalar or vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field (same number of components as B). 
 
  // If A<0 at corners, but A>0 near positivity control points, use cell average A.
  bool expA = false;
  bool avgA = false;
  if ((1.5*A[3]-0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) && (0.24*A[3]-0.3464101615137755*A[2]-0.3464101615137755*A[1]+0.5*A[0] > 0.0)) { 
    expA = true;
  }
  if (((-1.5*A[3])+0.8660254037844386*A[2]-0.8660254037844386*A[1]+0.5*A[0] < 0.0) && ((-0.24*A[3])+0.3464101615137755*A[2]-0.3464101615137755*A[1]+0.5*A[0] > 0.0)) { 
    expA = true;
  }
  if (((-1.5*A[3])-0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) && ((-0.24*A[3])-0.3464101615137755*A[2]+0.3464101615137755*A[1]+0.5*A[0] > 0.0)) { 
    expA = true;
  }
  if ((1.5*A[3]+0.8660254037844386*A[2]+0.8660254037844386*A[1]+0.5*A[0] < 0.0) && (0.24*A[3]+0.3464101615137755*A[2]+0.3464101615137755*A[1]+0.5*A[0] > 0.0)) { 
    expA = true;
  }
  // If A is zero near positivity control points, use cell average A.
  if (0.24*A[3]-0.3464101615137755*A[2]-0.3464101615137755*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-0.24*A[3])+0.3464101615137755*A[2]-0.3464101615137755*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if ((-0.24*A[3])-0.3464101615137755*A[2]+0.3464101615137755*A[1]+0.5*A[0] < 0.0) { 
    avgA = true;
  }
  if (0.24*A[3]+0.3464101615137755*A[2]+0.3464101615137755*A[1]+0.5*A[0] < 0.0) { 
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
 
  if ((avgA) || (!expA)) {
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
  } else {
    double xBar[4];
    double g1[4];
    xBar[0] = (1.0*As[3])/(1.732050807568877*As[2]-3.0*As[0])-(1.732050807568878*As[1])/(1.732050807568877*As[2]-3.0*As[0]); 
    xBar[1] = (1.0*As[3])/(1.732050807568877*As[2]+3.0*As[0])+(1.732050807568878*As[1])/(1.732050807568877*As[2]+3.0*As[0]); 

    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBar[0]*xBar[0])-(1.0*xBar[0]*xBar[0]*xBar[0])/(1.0-1.0*xBar[0]*xBar[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBar[1]*xBar[1])-(1.0*xBar[1]*xBar[1]*xBar[1])/(1.0-1.0*xBar[1]*xBar[1]); 

    xBar[2] = (1.0*As[3])/(1.732050807568877*As[1]-3.0*As[0])-(1.732050807568878*As[2])/(1.732050807568877*As[1]-3.0*As[0]); 
    xBar[3] = (1.0*As[3])/(1.732050807568877*As[1]+3.0*As[0])+(1.732050807568878*As[2])/(1.732050807568877*As[1]+3.0*As[0]); 

    g1[2] = (3.0*xBar[2])/(1.0-1.0*xBar[2]*xBar[2])-(1.0*xBar[2]*xBar[2]*xBar[2])/(1.0-1.0*xBar[2]*xBar[2]); 
    g1[3] = (3.0*xBar[3])/(1.0-1.0*xBar[3]*xBar[3])-(1.0*xBar[3]*xBar[3]*xBar[3])/(1.0-1.0*xBar[3]*xBar[3]); 

    // Fill AEM matrix. 
    data->AEM_S = Eigen::MatrixXd::Zero(4,4); 
    data->AEM_S(0,0) = 0.5*As[0]; 
    data->AEM_S(0,1) = 0.5*As[1]; 
    data->AEM_S(0,2) = 0.5*As[2]; 
    data->AEM_S(0,3) = 0.5*As[3]; 
    data->AEM_S(1,0) = 0.5*As[1]; 
    data->AEM_S(1,1) = (-(0.25*As[3])/g1[1])+(0.25*As[3])/g1[0]-(0.4330127018922193*As[1])/g1[1]-(0.4330127018922193*As[1])/g1[0]+As[0]; 
    data->AEM_S(1,2) = 0.5*As[3]; 
    data->AEM_S(1,3) = (-(0.4330127018922193*As[3])/g1[1])-(0.4330127018922193*As[3])/g1[0]+As[2]-(0.75*As[1])/g1[1]+(0.75*As[1])/g1[0]; 
    data->AEM_S(2,0) = 0.5*As[2]; 
    data->AEM_S(2,1) = 0.5*As[3]; 
    data->AEM_S(2,2) = (-(0.25*As[3])/g1[3])-(0.4330127018922193*As[2])/g1[3]+(0.25*As[3])/g1[2]-(0.4330127018922193*As[2])/g1[2]+As[0]; 
    data->AEM_S(2,3) = (-(0.4330127018922193*As[3])/g1[3])-(0.75*As[2])/g1[3]-(0.4330127018922193*As[3])/g1[2]+(0.75*As[2])/g1[2]+As[1]; 
    data->AEM_S(3,0) = 0.5*As[3]; 
    data->AEM_S(3,1) = (-(0.4330127018922193*As[3])/g1[1])-(0.4330127018922193*As[3])/g1[0]+As[2]-(0.75*As[1])/g1[1]+(0.75*As[1])/g1[0]; 
    data->AEM_S(3,2) = (-(0.4330127018922193*As[3])/g1[3])-(0.75*As[2])/g1[3]-(0.4330127018922193*As[3])/g1[2]+(0.75*As[2])/g1[2]+As[1]; 
    data->AEM_S(3,3) = (-(0.25*As[3])/g1[3])-(0.4330127018922193*As[2])/g1[3]+(0.25*As[3])/g1[2]-(0.25*As[3])/g1[1]+(0.25*As[3])/g1[0]-(0.4330127018922193*As[2])/g1[2]-(0.4330127018922193*As[1])/g1[1]-(0.4330127018922193*As[1])/g1[0]+1.5*As[0]; 

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
  };
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
 
void CartFieldBinOpDotProduct2xSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should be 1 here). 
  // out:     output field A.B (out only has one component). 
 
  // zero out. This is ok in this operator because there is no in-place dot-product. 
  for (unsigned short int vd=0; vd<12; vd++) 
  { 
    out[vd] = 0.0; 
  } 
 
  for (unsigned short int vd=0; vd<Ncomp; vd++) 
  { 
    unsigned short int a0 = 12*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    out[0] += 0.5*A[a0+11]*B[a0+11]+0.5*A[a0+10]*B[a0+10]+0.5*A[a0+9]*B[a0+9]+0.5*A[a0+8]*B[a0+8]+0.5*A[a0+7]*B[a0+7]+0.5*A[a0+6]*B[a0+6]+0.5*A[a0+5]*B[a0+5]+0.5*A[a0+4]*B[a0+4]+0.5*A[a0+3]*B[a0+3]+0.5*A[a0+2]*B[a0+2]+0.5*A[a0+1]*B[a0+1]+0.5*A[a0]*B[a0]; 
    out[1] += 0.5*A[a0+9]*B[a0+11]+0.5*B[a0+9]*A[a0+11]+0.4391550328268399*A[a0+6]*B[a0+10]+0.4391550328268399*B[a0+6]*A[a0+10]+0.4391550328268398*A[a0+4]*B[a0+8]+0.4391550328268398*B[a0+4]*A[a0+8]+0.5000000000000001*A[a0+5]*B[a0+7]+0.5000000000000001*B[a0+5]*A[a0+7]+0.447213595499958*A[a0+3]*B[a0+6]+0.447213595499958*B[a0+3]*A[a0+6]+0.4472135954999579*A[a0+1]*B[a0+4]+0.4472135954999579*B[a0+1]*A[a0+4]+0.5*A[a0+2]*B[a0+3]+0.5*B[a0+2]*A[a0+3]+0.5*A[a0]*B[a0+1]+0.5*B[a0]*A[a0+1]; 
    out[2] += 0.4391550328268399*A[a0+7]*B[a0+11]+0.4391550328268399*B[a0+7]*A[a0+11]+0.5*A[a0+8]*B[a0+10]+0.5*B[a0+8]*A[a0+10]+0.4391550328268398*A[a0+5]*B[a0+9]+0.4391550328268398*B[a0+5]*A[a0+9]+0.447213595499958*A[a0+3]*B[a0+7]+0.447213595499958*B[a0+3]*A[a0+7]+0.5000000000000001*A[a0+4]*B[a0+6]+0.5000000000000001*B[a0+4]*A[a0+6]+0.4472135954999579*A[a0+2]*B[a0+5]+0.4472135954999579*B[a0+2]*A[a0+5]+0.5*A[a0+1]*B[a0+3]+0.5*B[a0+1]*A[a0+3]+0.5*A[a0]*B[a0+2]+0.5*B[a0]*A[a0+2]; 
    out[3] += 0.4391550328268399*A[a0+5]*B[a0+11]+0.4391550328268399*B[a0+5]*A[a0+11]+0.4391550328268399*A[a0+4]*B[a0+10]+0.4391550328268399*B[a0+4]*A[a0+10]+0.4391550328268399*A[a0+7]*B[a0+9]+0.4391550328268399*B[a0+7]*A[a0+9]+0.4391550328268399*A[a0+6]*B[a0+8]+0.4391550328268399*B[a0+6]*A[a0+8]+0.4*A[a0+6]*B[a0+7]+0.447213595499958*A[a0+2]*B[a0+7]+0.4*B[a0+6]*A[a0+7]+0.447213595499958*B[a0+2]*A[a0+7]+0.447213595499958*A[a0+1]*B[a0+6]+0.447213595499958*B[a0+1]*A[a0+6]+0.4472135954999579*A[a0+3]*B[a0+5]+0.4472135954999579*B[a0+3]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+4]+0.4472135954999579*B[a0+3]*A[a0+4]+0.5*A[a0]*B[a0+3]+0.5*B[a0]*A[a0+3]+0.5*A[a0+1]*B[a0+2]+0.5*B[a0+1]*A[a0+2]; 
    out[4] += 0.4472135954999579*A[a0+11]*B[a0+11]+0.2981423969999719*A[a0+10]*B[a0+10]+0.4391550328268399*A[a0+3]*B[a0+10]+0.4391550328268399*B[a0+3]*A[a0+10]+0.2981423969999719*A[a0+8]*B[a0+8]+0.4391550328268398*A[a0+1]*B[a0+8]+0.4391550328268398*B[a0+1]*A[a0+8]+0.4472135954999579*A[a0+7]*B[a0+7]+0.31943828249997*A[a0+6]*B[a0+6]+0.5000000000000001*A[a0+2]*B[a0+6]+0.5000000000000001*B[a0+2]*A[a0+6]+0.31943828249997*A[a0+4]*B[a0+4]+0.5*A[a0]*B[a0+4]+0.5*B[a0]*A[a0+4]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+1]*B[a0+1]; 
    out[5] += 0.2981423969999719*A[a0+11]*B[a0+11]+0.4391550328268399*A[a0+3]*B[a0+11]+0.4391550328268399*B[a0+3]*A[a0+11]+0.4472135954999579*A[a0+10]*B[a0+10]+0.2981423969999719*A[a0+9]*B[a0+9]+0.4391550328268398*A[a0+2]*B[a0+9]+0.4391550328268398*B[a0+2]*A[a0+9]+0.31943828249997*A[a0+7]*B[a0+7]+0.5000000000000001*A[a0+1]*B[a0+7]+0.5000000000000001*B[a0+1]*A[a0+7]+0.4472135954999579*A[a0+6]*B[a0+6]+0.31943828249997*A[a0+5]*B[a0+5]+0.5*A[a0]*B[a0+5]+0.5*B[a0]*A[a0+5]+0.4472135954999579*A[a0+3]*B[a0+3]+0.4472135954999579*A[a0+2]*B[a0+2]; 
    out[6] += 0.3927922024247863*A[a0+7]*B[a0+11]+0.3927922024247863*B[a0+7]*A[a0+11]+0.2981423969999719*A[a0+8]*B[a0+10]+0.3927922024247863*A[a0+7]*B[a0+10]+0.4391550328268399*A[a0+1]*B[a0+10]+0.2981423969999719*B[a0+8]*A[a0+10]+0.3927922024247863*B[a0+7]*A[a0+10]+0.4391550328268399*B[a0+1]*A[a0+10]+0.4391550328268399*A[a0+3]*B[a0+8]+0.4391550328268399*B[a0+3]*A[a0+8]+0.4*A[a0+3]*B[a0+7]+0.4*B[a0+3]*A[a0+7]+0.4472135954999579*A[a0+5]*B[a0+6]+0.31943828249997*A[a0+4]*B[a0+6]+0.5*A[a0]*B[a0+6]+0.4472135954999579*B[a0+5]*A[a0+6]+0.31943828249997*B[a0+4]*A[a0+6]+0.5*B[a0]*A[a0+6]+0.5000000000000001*A[a0+2]*B[a0+4]+0.5000000000000001*B[a0+2]*A[a0+4]+0.447213595499958*A[a0+1]*B[a0+3]+0.447213595499958*B[a0+1]*A[a0+3]; 
    out[7] += 0.2981423969999719*A[a0+9]*B[a0+11]+0.3927922024247863*A[a0+6]*B[a0+11]+0.4391550328268399*A[a0+2]*B[a0+11]+0.2981423969999719*B[a0+9]*A[a0+11]+0.3927922024247863*B[a0+6]*A[a0+11]+0.4391550328268399*B[a0+2]*A[a0+11]+0.3927922024247863*A[a0+6]*B[a0+10]+0.3927922024247863*B[a0+6]*A[a0+10]+0.4391550328268399*A[a0+3]*B[a0+9]+0.4391550328268399*B[a0+3]*A[a0+9]+0.31943828249997*A[a0+5]*B[a0+7]+0.4472135954999579*A[a0+4]*B[a0+7]+0.5*A[a0]*B[a0+7]+0.31943828249997*B[a0+5]*A[a0+7]+0.4472135954999579*B[a0+4]*A[a0+7]+0.5*B[a0]*A[a0+7]+0.4*A[a0+3]*B[a0+6]+0.4*B[a0+3]*A[a0+6]+0.5000000000000001*A[a0+1]*B[a0+5]+0.5000000000000001*B[a0+1]*A[a0+5]+0.447213595499958*A[a0+2]*B[a0+3]+0.447213595499958*B[a0+2]*A[a0+3]; 
    out[8] += 0.2981423969999719*A[a0+6]*B[a0+10]+0.5*A[a0+2]*B[a0+10]+0.2981423969999719*B[a0+6]*A[a0+10]+0.5*B[a0+2]*A[a0+10]+0.2981423969999719*A[a0+4]*B[a0+8]+0.5*A[a0]*B[a0+8]+0.2981423969999719*B[a0+4]*A[a0+8]+0.5*B[a0]*A[a0+8]+0.4391550328268399*A[a0+3]*B[a0+6]+0.4391550328268399*B[a0+3]*A[a0+6]+0.4391550328268398*A[a0+1]*B[a0+4]+0.4391550328268398*B[a0+1]*A[a0+4]; 
    out[9] += 0.2981423969999719*A[a0+7]*B[a0+11]+0.5*A[a0+1]*B[a0+11]+0.2981423969999719*B[a0+7]*A[a0+11]+0.5*B[a0+1]*A[a0+11]+0.2981423969999719*A[a0+5]*B[a0+9]+0.5*A[a0]*B[a0+9]+0.2981423969999719*B[a0+5]*A[a0+9]+0.5*B[a0]*A[a0+9]+0.4391550328268399*A[a0+3]*B[a0+7]+0.4391550328268399*B[a0+3]*A[a0+7]+0.4391550328268398*A[a0+2]*B[a0+5]+0.4391550328268398*B[a0+2]*A[a0+5]; 
    out[10] += 0.4472135954999579*A[a0+5]*B[a0+10]+0.2981423969999719*A[a0+4]*B[a0+10]+0.5*A[a0]*B[a0+10]+0.4472135954999579*B[a0+5]*A[a0+10]+0.2981423969999719*B[a0+4]*A[a0+10]+0.5*B[a0]*A[a0+10]+0.2981423969999719*A[a0+6]*B[a0+8]+0.5*A[a0+2]*B[a0+8]+0.2981423969999719*B[a0+6]*A[a0+8]+0.5*B[a0+2]*A[a0+8]+0.3927922024247863*A[a0+6]*B[a0+7]+0.3927922024247863*B[a0+6]*A[a0+7]+0.4391550328268399*A[a0+1]*B[a0+6]+0.4391550328268399*B[a0+1]*A[a0+6]+0.4391550328268399*A[a0+3]*B[a0+4]+0.4391550328268399*B[a0+3]*A[a0+4]; 
    out[11] += 0.2981423969999719*A[a0+5]*B[a0+11]+0.4472135954999579*A[a0+4]*B[a0+11]+0.5*A[a0]*B[a0+11]+0.2981423969999719*B[a0+5]*A[a0+11]+0.4472135954999579*B[a0+4]*A[a0+11]+0.5*B[a0]*A[a0+11]+0.2981423969999719*A[a0+7]*B[a0+9]+0.5*A[a0+1]*B[a0+9]+0.2981423969999719*B[a0+7]*A[a0+9]+0.5*B[a0+1]*A[a0+9]+0.3927922024247863*A[a0+6]*B[a0+7]+0.4391550328268399*A[a0+2]*B[a0+7]+0.3927922024247863*B[a0+6]*A[a0+7]+0.4391550328268399*B[a0+2]*A[a0+7]+0.4391550328268399*A[a0+3]*B[a0+5]+0.4391550328268399*B[a0+3]*A[a0+5]; 
  } 
 
} 
 
