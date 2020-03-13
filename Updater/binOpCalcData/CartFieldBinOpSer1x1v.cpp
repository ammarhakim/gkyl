#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply1x1vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[4]; 
 
  tmp[0] = 0.5*(A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.5*(A[2]*B[3]+B[2]*A[3]+A[0]*B[1]+B[0]*A[1]); 
  tmp[2] = 0.5*(A[1]*B[3]+B[1]*A[3]+A[0]*B[2]+B[0]*A[2]); 
  tmp[3] = 0.5*(A[0]*B[3]+B[0]*A[3]+A[1]*B[2]+B[1]*A[2]); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<4; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x1vSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[8]; 
 
  tmp[0] = 0.5*(A[7]*B[7]+A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.03333333333333333*(15.0*A[5]*B[7]+15.0*B[5]*A[7]+13.41640786499874*A[3]*B[6]+13.41640786499874*B[3]*A[6]+13.41640786499874*A[1]*B[4]+13.41640786499874*B[1]*A[4]+15.0*A[2]*B[3]+15.0*B[2]*A[3]+15.0*A[0]*B[1]+15.0*B[0]*A[1]); 
  tmp[2] = 0.03333333333333333*(13.41640786499874*A[3]*B[7]+13.41640786499874*B[3]*A[7]+15.0*A[4]*B[6]+15.0*B[4]*A[6]+13.41640786499874*A[2]*B[5]+13.41640786499874*B[2]*A[5]+15.0*A[1]*B[3]+15.0*B[1]*A[3]+15.0*A[0]*B[2]+15.0*B[0]*A[2]); 
  tmp[3] = 0.03333333333333333*((12.0*A[6]+13.41640786499874*A[2])*B[7]+(12.0*B[6]+13.41640786499874*B[2])*A[7]+13.41640786499874*A[1]*B[6]+13.41640786499874*B[1]*A[6]+13.41640786499874*A[3]*B[5]+13.41640786499874*B[3]*A[5]+13.41640786499874*A[3]*B[4]+13.41640786499874*B[3]*A[4]+15.0*A[0]*B[3]+15.0*B[0]*A[3]+15.0*A[1]*B[2]+15.0*B[1]*A[2]); 
  tmp[4] = 0.004761904761904762*(93.91485505499116*A[7]*B[7]+(67.0820393249937*A[6]+105.0*A[2])*B[6]+105.0*B[2]*A[6]+(67.0820393249937*A[4]+105.0*A[0])*B[4]+105.0*B[0]*A[4]+93.91485505499116*A[3]*B[3]+93.91485505499116*A[1]*B[1]); 
  tmp[5] = 0.004761904761904762*((67.0820393249937*A[7]+105.0*A[1])*B[7]+105.0*B[1]*A[7]+93.91485505499116*A[6]*B[6]+(67.0820393249937*A[5]+105.0*A[0])*B[5]+105.0*B[0]*A[5]+93.91485505499116*A[3]*B[3]+93.91485505499116*A[2]*B[2]); 
  tmp[6] = 0.004761904761904762*(84.0*A[3]*B[7]+84.0*B[3]*A[7]+(93.91485505499116*A[5]+67.0820393249937*A[4]+105.0*A[0])*B[6]+(93.91485505499116*B[5]+67.0820393249937*B[4]+105.0*B[0])*A[6]+105.0*A[2]*B[4]+105.0*B[2]*A[4]+93.91485505499116*A[1]*B[3]+93.91485505499116*B[1]*A[3]); 
  tmp[7] = 0.004761904761904762*((67.0820393249937*A[5]+93.91485505499116*A[4]+105.0*A[0])*B[7]+(67.0820393249937*B[5]+93.91485505499116*B[4]+105.0*B[0])*A[7]+84.0*A[3]*B[6]+84.0*B[3]*A[6]+105.0*A[1]*B[5]+105.0*B[1]*A[5]+93.91485505499116*A[2]*B[3]+93.91485505499116*B[2]*A[3]); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<8; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply1x1vSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field. 
  // B:       scalar/vector field (must be vector if A is vector). 
  // Ncomp:   number of components of B (could be 1D, 2D, 3D, vector). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else. 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[12]; 
 
  tmp[0] = 0.5*(A[11]*B[11]+A[10]*B[10]+A[9]*B[9]+A[8]*B[8]+A[7]*B[7]+A[6]*B[6]+A[5]*B[5]+A[4]*B[4]+A[3]*B[3]+A[2]*B[2]+A[1]*B[1]+A[0]*B[0]); 
  tmp[1] = 0.004761904761904762*(105.0*A[9]*B[11]+105.0*B[9]*A[11]+92.22255689363637*A[6]*B[10]+92.22255689363637*B[6]*A[10]+92.22255689363637*A[4]*B[8]+92.22255689363637*B[4]*A[8]+105.0*A[5]*B[7]+105.0*B[5]*A[7]+93.91485505499116*A[3]*B[6]+93.91485505499116*B[3]*A[6]+93.91485505499116*A[1]*B[4]+93.91485505499116*B[1]*A[4]+105.0*A[2]*B[3]+105.0*B[2]*A[3]+105.0*A[0]*B[1]+105.0*B[0]*A[1]); 
  tmp[2] = 0.004761904761904762*(92.22255689363637*A[7]*B[11]+92.22255689363637*B[7]*A[11]+105.0*A[8]*B[10]+105.0*B[8]*A[10]+92.22255689363637*A[5]*B[9]+92.22255689363637*B[5]*A[9]+93.91485505499116*A[3]*B[7]+93.91485505499116*B[3]*A[7]+105.0*A[4]*B[6]+105.0*B[4]*A[6]+93.91485505499116*A[2]*B[5]+93.91485505499116*B[2]*A[5]+105.0*A[1]*B[3]+105.0*B[1]*A[3]+105.0*A[0]*B[2]+105.0*B[0]*A[2]); 
  tmp[3] = 0.004761904761904762*(92.22255689363638*A[5]*B[11]+92.22255689363638*B[5]*A[11]+92.22255689363638*A[4]*B[10]+92.22255689363638*B[4]*A[10]+92.2225568936364*A[7]*B[9]+92.2225568936364*B[7]*A[9]+92.2225568936364*A[6]*B[8]+92.2225568936364*B[6]*A[8]+(84.0*A[6]+93.91485505499116*A[2])*B[7]+(84.0*B[6]+93.91485505499116*B[2])*A[7]+93.91485505499116*A[1]*B[6]+93.91485505499116*B[1]*A[6]+93.91485505499116*A[3]*B[5]+93.91485505499116*B[3]*A[5]+93.91485505499116*A[3]*B[4]+93.91485505499116*B[3]*A[4]+105.0*A[0]*B[3]+105.0*B[0]*A[3]+105.0*A[1]*B[2]+105.0*B[1]*A[2]); 
  tmp[4] = 0.004761904761904762*(93.91485505499116*A[11]*B[11]+(62.60990336999411*A[10]+92.22255689363638*A[3])*B[10]+92.22255689363638*B[3]*A[10]+(62.60990336999411*A[8]+92.22255689363637*A[1])*B[8]+92.22255689363637*B[1]*A[8]+93.91485505499116*A[7]*B[7]+(67.0820393249937*A[6]+105.0*A[2])*B[6]+105.0*B[2]*A[6]+(67.0820393249937*A[4]+105.0*A[0])*B[4]+105.0*B[0]*A[4]+93.91485505499116*A[3]*B[3]+93.91485505499116*A[1]*B[1]); 
  tmp[5] = 0.004761904761904762*((62.60990336999411*A[11]+92.22255689363638*A[3])*B[11]+92.22255689363638*B[3]*A[11]+93.91485505499116*A[10]*B[10]+(62.60990336999411*A[9]+92.22255689363637*A[2])*B[9]+92.22255689363637*B[2]*A[9]+(67.0820393249937*A[7]+105.0*A[1])*B[7]+105.0*B[1]*A[7]+93.91485505499116*A[6]*B[6]+(67.0820393249937*A[5]+105.0*A[0])*B[5]+105.0*B[0]*A[5]+93.91485505499116*A[3]*B[3]+93.91485505499116*A[2]*B[2]); 
  tmp[6] = 0.001587301587301587*(247.4590875276153*A[7]*B[11]+247.4590875276153*B[7]*A[11]+(187.8297101099824*A[8]+247.4590875276153*A[7]+276.6676706809091*A[1])*B[10]+(187.8297101099824*B[8]+247.4590875276153*B[7]+276.6676706809091*B[1])*A[10]+276.6676706809092*A[3]*B[8]+276.6676706809092*B[3]*A[8]+252.0*A[3]*B[7]+252.0*B[3]*A[7]+(281.7445651649735*A[5]+201.2461179749811*A[4]+315.0*A[0])*B[6]+(281.7445651649735*B[5]+201.2461179749811*B[4]+315.0*B[0])*A[6]+314.9999999999999*A[2]*B[4]+314.9999999999999*B[2]*A[4]+281.7445651649734*A[1]*B[3]+281.7445651649734*B[1]*A[3]); 
  tmp[7] = 0.001587301587301587*((187.8297101099824*A[9]+247.4590875276153*A[6]+276.6676706809091*A[2])*B[11]+(187.8297101099824*B[9]+247.4590875276153*B[6]+276.6676706809091*B[2])*A[11]+247.4590875276153*A[6]*B[10]+247.4590875276153*B[6]*A[10]+276.6676706809092*A[3]*B[9]+276.6676706809092*B[3]*A[9]+(201.2461179749811*A[5]+281.7445651649735*A[4]+315.0*A[0])*B[7]+(201.2461179749811*B[5]+281.7445651649735*B[4]+315.0*B[0])*A[7]+252.0*A[3]*B[6]+252.0*B[3]*A[6]+314.9999999999999*A[1]*B[5]+314.9999999999999*B[1]*A[5]+281.7445651649734*A[2]*B[3]+281.7445651649734*B[2]*A[3]); 
  tmp[8] = 0.001587301587301587*((187.8297101099824*A[6]+315.0*A[2])*B[10]+(187.8297101099824*B[6]+315.0*B[2])*A[10]+(187.8297101099823*A[4]+315.0*A[0])*B[8]+(187.8297101099823*B[4]+315.0*B[0])*A[8]+276.6676706809092*A[3]*B[6]+276.6676706809092*B[3]*A[6]+276.6676706809091*A[1]*B[4]+276.6676706809091*B[1]*A[4]); 
  tmp[9] = 0.001587301587301587*((187.8297101099824*A[7]+315.0*A[1])*B[11]+(187.8297101099824*B[7]+315.0*B[1])*A[11]+(187.8297101099823*A[5]+315.0*A[0])*B[9]+(187.8297101099823*B[5]+315.0*B[0])*A[9]+276.6676706809092*A[3]*B[7]+276.6676706809092*B[3]*A[7]+276.6676706809091*A[2]*B[5]+276.6676706809091*B[2]*A[5]); 
  tmp[10] = 0.001587301587301587*((281.7445651649735*A[5]+187.8297101099823*A[4]+315.0*A[0])*B[10]+(281.7445651649735*B[5]+187.8297101099823*B[4]+315.0*B[0])*A[10]+(187.8297101099824*A[6]+315.0*A[2])*B[8]+(187.8297101099824*B[6]+315.0*B[2])*A[8]+247.4590875276153*A[6]*B[7]+247.4590875276153*B[6]*A[7]+276.6676706809091*A[1]*B[6]+276.6676706809091*B[1]*A[6]+276.6676706809092*A[3]*B[4]+276.6676706809092*B[3]*A[4]); 
  tmp[11] = 0.001587301587301587*((187.8297101099823*A[5]+281.7445651649735*A[4]+315.0*A[0])*B[11]+(187.8297101099823*B[5]+281.7445651649735*B[4]+315.0*B[0])*A[11]+(187.8297101099824*A[7]+315.0*A[1])*B[9]+(187.8297101099824*B[7]+315.0*B[1])*A[9]+(247.4590875276153*A[6]+276.6676706809091*A[2])*B[7]+(247.4590875276153*B[6]+276.6676706809091*B[2])*A[7]+276.6676706809092*A[3]*B[5]+276.6676706809092*B[3]*A[5]); 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<12; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x1vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[4]; 
  tmp[0] = 0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[1]*B[3]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.7071067811865475*A[0]*B[3]+0.7071067811865475*A[1]*B[2]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<4; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x1vSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[8]; 
  tmp[0] = 0.7071067811865475*A[2]*B[4]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6324555320336759*A[1]*B[4]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865475*A[2]*B[6]+0.7071067811865475*A[1]*B[3]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.632455532033676*A[1]*B[6]+0.6324555320336759*A[2]*B[3]+0.7071067811865475*A[0]*B[3]+0.7071067811865475*A[1]*B[2]; 
  tmp[4] = 0.4517539514526256*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[5] = 0.7071067811865475*A[1]*B[7]+0.7071067811865475*A[0]*B[5]; 
  tmp[6] = 0.4517539514526256*A[2]*B[6]+0.7071067811865475*A[0]*B[6]+0.632455532033676*A[1]*B[3]+0.7071067811865475*A[2]*B[2]; 
  tmp[7] = 0.6324555320336759*A[2]*B[7]+0.7071067811865475*A[0]*B[7]+0.7071067811865475*A[1]*B[5]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<8; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseMultiply1x1vSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[12]; 
  tmp[0] = 0.7071067811865475*A[3]*B[8]+0.7071067811865475*A[2]*B[4]+0.7071067811865475*A[1]*B[1]+0.7071067811865475*A[0]*B[0]; 
  tmp[1] = 0.6210590034081186*A[2]*B[8]+0.6210590034081186*A[3]*B[4]+0.6324555320336759*A[1]*B[4]+0.6324555320336759*B[1]*A[2]+0.7071067811865475*A[0]*B[1]+0.7071067811865475*B[0]*A[1]; 
  tmp[2] = 0.7071067811865474*A[3]*B[10]+0.7071067811865475*A[2]*B[6]+0.7071067811865475*A[1]*B[3]+0.7071067811865475*A[0]*B[2]; 
  tmp[3] = 0.6210590034081187*A[2]*B[10]+0.6210590034081187*A[3]*B[6]+0.632455532033676*A[1]*B[6]+0.6324555320336759*A[2]*B[3]+0.7071067811865475*A[0]*B[3]+0.7071067811865475*A[1]*B[2]; 
  tmp[4] = 0.421637021355784*A[3]*B[8]+0.6210590034081186*A[1]*B[8]+0.4517539514526256*A[2]*B[4]+0.7071067811865475*A[0]*B[4]+0.6210590034081186*B[1]*A[3]+0.7071067811865475*B[0]*A[2]+0.6324555320336759*A[1]*B[1]; 
  tmp[5] = 0.7071067811865475*A[1]*B[7]+0.7071067811865475*A[0]*B[5]; 
  tmp[6] = 0.4216370213557839*A[3]*B[10]+0.6210590034081187*A[1]*B[10]+0.4517539514526256*A[2]*B[6]+0.7071067811865475*A[0]*B[6]+0.6210590034081187*A[3]*B[3]+0.632455532033676*A[1]*B[3]+0.7071067811865475*A[2]*B[2]; 
  tmp[7] = 0.6324555320336759*A[2]*B[7]+0.7071067811865475*A[0]*B[7]+0.7071067811865475*A[1]*B[5]; 
  tmp[8] = 0.421637021355784*A[2]*B[8]+0.7071067811865475*A[0]*B[8]+0.421637021355784*A[3]*B[4]+0.6210590034081186*A[1]*B[4]+0.7071067811865475*B[0]*A[3]+0.6210590034081186*B[1]*A[2]; 
  tmp[9] = 0.7071067811865474*A[1]*B[11]+0.7071067811865475*A[0]*B[9]; 
  tmp[10] = 0.421637021355784*A[2]*B[10]+0.7071067811865475*A[0]*B[10]+0.4216370213557839*A[3]*B[6]+0.6210590034081187*A[1]*B[6]+0.6210590034081187*A[2]*B[3]+0.7071067811865474*B[2]*A[3]; 
  tmp[11] = 0.6324555320336759*A[2]*B[11]+0.7071067811865475*A[0]*B[11]+0.7071067811865474*A[1]*B[9]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<12; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x1vSer_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.7071067811865475*A[0]-1.224744871391589*A[1] < 0) { 
    avgA = true;
  }
  if (1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  double As[2]; 
  double Bs[4]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = 0.0; 
  } else { 
    As[0] = A[0]; 
    As[1] = A[1]; 
    Bs[0] = B[0]; 
    Bs[1] = B[1]; 
    Bs[2] = B[2]; 
    Bs[3] = B[3]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(4,4); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,2) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,3) = 0.7071067811865475*As[0]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,4,1) = data->u_D; 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x1vSer_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
  if (1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  double As[3]; 
  double Bs[8]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = 0.0; 
    Bs[4] = 0.0; 
    Bs[5] = B[5]; 
    Bs[6] = 0.0; 
    Bs[7] = 0.0; 
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
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(8,8); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,3) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,4) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_D(1,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(1,3) = 0.7071067811865475*As[1]; 
  data->AEM_D(1,4) = 0.7071067811865475*As[2]; 
  data->AEM_D(1,5) = 0.6324555320336759*As[1]; 
  data->AEM_D(2,4) = 0.7071067811865475*As[2]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,8,1) = data->u_D; 
 
} 
 
void CartFieldBinOpConfPhaseDivide1x1vSer_P3(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-1.870828693386971*A[3])+1.58113883008419*A[2]-1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
  if (1.870828693386971*A[3]+1.58113883008419*A[2]+1.224744871391589*A[1]+0.7071067811865475*A[0] < 0) { 
    avgA = true;
  }
 
  double As[4]; 
  double Bs[12]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = B[2]; 
    Bs[3] = 0.0; 
    Bs[4] = 0.0; 
    Bs[5] = B[5]; 
    Bs[6] = 0.0; 
    Bs[7] = 0.0; 
    Bs[8] = 0.0; 
    Bs[9] = B[9]; 
    Bs[10] = 0.0; 
    Bs[11] = 0.0; 
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
    Bs[7] = B[7]; 
    Bs[8] = B[8]; 
    Bs[9] = B[9]; 
    Bs[10] = B[10]; 
    Bs[11] = B[11]; 
  } 
 
  // Fill AEM matrix. 
  data->AEM_D = Eigen::MatrixXd::Zero(12,12); 
  data->AEM_D(0,0) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,1) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,4) = 0.7071067811865475*As[1]; 
  data->AEM_D(0,5) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_D(0,10) = 0.7071067811865475*As[0]; 
  data->AEM_D(0,11) = 0.7071067811865475*As[1]; 
  data->AEM_D(1,2) = 0.7071067811865475*As[1]; 
  data->AEM_D(1,3) = 0.6324555320336759*As[2]+0.7071067811865475*As[0]; 
  data->AEM_D(1,4) = 0.7071067811865475*As[2]; 
  data->AEM_D(1,5) = 0.6210590034081186*As[3]+0.6324555320336759*As[1]; 
  data->AEM_D(2,2) = 0.7071067811865475*As[2]; 
  data->AEM_D(2,3) = 0.6210590034081187*As[3]+0.632455532033676*As[1]; 
  data->AEM_D(2,8) = 0.7071067811865475*As[3]; 
  data->AEM_D(2,9) = 0.6210590034081186*As[2]; 
  data->AEM_D(3,6) = 0.7071067811865474*As[3]; 
  data->AEM_D(3,7) = 0.6210590034081187*As[2]; 
 
  // Fill BEV. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9],Bs[10],Bs[11]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,12,1) = data->u_D; 
 
} 
 
