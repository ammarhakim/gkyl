#include <math.h> 
#include <CartFieldBinOpModDecl.h> 
 
using namespace Eigen; 
 
void CartFieldBinOpMultiply3x3vMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[7]; 
  tmp[0] = 0.3535533905932737*A[3]*B[3]+0.3535533905932737*A[2]*B[2]+0.3535533905932737*A[1]*B[1]+0.3535533905932737*A[0]*B[0]; 
  tmp[1] = 0.3535533905932737*A[0]*B[1]+0.3535533905932737*B[0]*A[1]; 
  tmp[2] = 0.3535533905932737*A[0]*B[2]+0.3535533905932737*B[0]*A[2]; 
  tmp[3] = 0.3535533905932737*A[0]*B[3]+0.3535533905932737*B[0]*A[3]; 
  tmp[4] = 0.3535533905932737*A[0]*B[4]; 
  tmp[5] = 0.3535533905932737*A[0]*B[5]; 
  tmp[6] = 0.3535533905932737*A[0]*B[6]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<7; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpMultiply3x3vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       scalar/vector field in configuration space. 
  // B:       scalar field in phase space. 
  // Ncomp:   number of components of B (should =1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (should=1 here). 
  // out:     output field A*B (same number of components as B). 
 
  double tmp[28]; 
  tmp[0] = 0.3535533905932737*A[9]*B[24]+0.3535533905932737*A[8]*B[23]+0.3535533905932737*A[7]*B[22]+0.3535533905932737*A[6]*B[9]+0.3535533905932737*A[5]*B[8]+0.3535533905932737*A[4]*B[7]+0.3535533905932737*A[3]*B[3]+0.3535533905932737*A[2]*B[2]+0.3535533905932737*A[1]*B[1]+0.3535533905932737*A[0]*B[0]; 
  tmp[1] = 0.3162277660168379*A[1]*B[22]+0.3535533905932737*A[3]*B[8]+0.3535533905932737*A[2]*B[7]+0.3162277660168379*B[1]*A[7]+0.3535533905932737*B[3]*A[5]+0.3535533905932737*B[2]*A[4]+0.3535533905932737*A[0]*B[1]+0.3535533905932737*B[0]*A[1]; 
  tmp[2] = 0.3162277660168379*A[2]*B[23]+0.3535533905932737*A[3]*B[9]+0.3162277660168379*B[2]*A[8]+0.3535533905932737*A[1]*B[7]+0.3535533905932737*B[3]*A[6]+0.3535533905932737*B[1]*A[4]+0.3535533905932737*A[0]*B[2]+0.3535533905932737*B[0]*A[2]; 
  tmp[3] = 0.3162277660168379*A[3]*B[24]+0.3535533905932737*A[2]*B[9]+0.3162277660168379*B[3]*A[9]+0.3535533905932737*A[1]*B[8]+0.3535533905932737*B[2]*A[6]+0.3535533905932737*B[1]*A[5]+0.3535533905932737*A[0]*B[3]+0.3535533905932737*B[0]*A[3]; 
  tmp[4] = 0.3535533905932737*A[3]*B[12]+0.3535533905932737*A[2]*B[11]+0.3535533905932737*A[1]*B[10]+0.3535533905932737*A[0]*B[4]; 
  tmp[5] = 0.3535533905932737*A[3]*B[15]+0.3535533905932737*A[2]*B[14]+0.3535533905932737*A[1]*B[13]+0.3535533905932737*A[0]*B[5]; 
  tmp[6] = 0.3535533905932737*A[3]*B[19]+0.3535533905932737*A[2]*B[18]+0.3535533905932737*A[1]*B[17]+0.3535533905932737*A[0]*B[6]; 
  tmp[7] = 0.3162277660168379*A[4]*B[23]+0.3162277660168379*A[4]*B[22]+0.3535533905932737*A[5]*B[9]+0.3535533905932737*A[6]*B[8]+0.3162277660168379*B[7]*A[8]+0.3162277660168379*A[7]*B[7]+0.3535533905932737*A[0]*B[7]+0.3535533905932737*B[0]*A[4]+0.3535533905932737*A[1]*B[2]+0.3535533905932737*B[1]*A[2]; 
  tmp[8] = 0.3162277660168379*A[5]*B[24]+0.3162277660168379*A[5]*B[22]+0.3535533905932737*A[4]*B[9]+0.3162277660168379*B[8]*A[9]+0.3162277660168379*A[7]*B[8]+0.3535533905932737*A[0]*B[8]+0.3535533905932737*A[6]*B[7]+0.3535533905932737*B[0]*A[5]+0.3535533905932737*A[1]*B[3]+0.3535533905932737*B[1]*A[3]; 
  tmp[9] = 0.3162277660168379*A[6]*B[24]+0.3162277660168379*A[6]*B[23]+0.3162277660168379*A[9]*B[9]+0.3162277660168379*A[8]*B[9]+0.3535533905932737*A[0]*B[9]+0.3535533905932737*A[4]*B[8]+0.3535533905932737*A[5]*B[7]+0.3535533905932737*B[0]*A[6]+0.3535533905932737*A[2]*B[3]+0.3535533905932737*B[2]*A[3]; 
  tmp[10] = 0.3535533905932737*A[5]*B[12]+0.3535533905932737*A[4]*B[11]+0.3162277660168379*A[7]*B[10]+0.3535533905932737*A[0]*B[10]+0.3535533905932737*A[1]*B[4]; 
  tmp[11] = 0.3535533905932737*A[6]*B[12]+0.3162277660168379*A[8]*B[11]+0.3535533905932737*A[0]*B[11]+0.3535533905932737*A[4]*B[10]+0.3535533905932737*A[2]*B[4]; 
  tmp[12] = 0.3162277660168379*A[9]*B[12]+0.3535533905932737*A[0]*B[12]+0.3535533905932737*A[6]*B[11]+0.3535533905932737*A[5]*B[10]+0.3535533905932737*A[3]*B[4]; 
  tmp[13] = 0.3535533905932737*A[5]*B[15]+0.3535533905932737*A[4]*B[14]+0.3162277660168379*A[7]*B[13]+0.3535533905932737*A[0]*B[13]+0.3535533905932737*A[1]*B[5]; 
  tmp[14] = 0.3535533905932737*A[6]*B[15]+0.3162277660168379*A[8]*B[14]+0.3535533905932737*A[0]*B[14]+0.3535533905932737*A[4]*B[13]+0.3535533905932737*A[2]*B[5]; 
  tmp[15] = 0.3162277660168379*A[9]*B[15]+0.3535533905932737*A[0]*B[15]+0.3535533905932737*A[6]*B[14]+0.3535533905932737*A[5]*B[13]+0.3535533905932737*A[3]*B[5]; 
  tmp[16] = 0.3535533905932737*A[0]*B[16]; 
  tmp[17] = 0.3535533905932737*A[5]*B[19]+0.3535533905932737*A[4]*B[18]+0.3162277660168379*A[7]*B[17]+0.3535533905932737*A[0]*B[17]+0.3535533905932737*A[1]*B[6]; 
  tmp[18] = 0.3535533905932737*A[6]*B[19]+0.3162277660168379*A[8]*B[18]+0.3535533905932737*A[0]*B[18]+0.3535533905932737*A[4]*B[17]+0.3535533905932737*A[2]*B[6]; 
  tmp[19] = 0.3162277660168379*A[9]*B[19]+0.3535533905932737*A[0]*B[19]+0.3535533905932737*A[6]*B[18]+0.3535533905932737*A[5]*B[17]+0.3535533905932737*A[3]*B[6]; 
  tmp[20] = 0.3535533905932737*A[0]*B[20]; 
  tmp[21] = 0.3535533905932737*A[0]*B[21]; 
  tmp[22] = 0.2258769757263128*A[7]*B[22]+0.3535533905932737*A[0]*B[22]+0.3162277660168379*A[5]*B[8]+0.3162277660168379*A[4]*B[7]+0.3535533905932737*B[0]*A[7]+0.3162277660168379*A[1]*B[1]; 
  tmp[23] = 0.2258769757263128*A[8]*B[23]+0.3535533905932737*A[0]*B[23]+0.3162277660168379*A[6]*B[9]+0.3535533905932737*B[0]*A[8]+0.3162277660168379*A[4]*B[7]+0.3162277660168379*A[2]*B[2]; 
  tmp[24] = 0.2258769757263128*A[9]*B[24]+0.3535533905932737*A[0]*B[24]+0.3162277660168379*A[6]*B[9]+0.3535533905932737*B[0]*A[9]+0.3162277660168379*A[5]*B[8]+0.3162277660168379*A[3]*B[3]; 
  tmp[25] = 0.3535533905932737*A[0]*B[25]; 
  tmp[26] = 0.3535533905932737*A[0]*B[26]; 
  tmp[27] = 0.3535533905932737*A[0]*B[27]; 
 
  // This tmp allows for in-place multiplication. 
  for (unsigned short int i=0; i<28; i++) 
  { 
    out[i] = tmp[i]; 
  } 
 
} 
 
void CartFieldBinOpDivide3x3vMax_P1(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if ((-0.6123724356957944*A[3])-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
 
  double As[4]; 
  double Bs[7]; 
  if (avgA) { 
    As[0] = A[0]; 
    As[1] = 0.0; 
    As[2] = 0.0; 
    As[3] = 0.0; 
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = 0.0; 
    Bs[3] = 0.0; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
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
  } 
 
  // Fill AEM_D matrix. 
  data->AEM_D(0,0) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,1) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,2) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,3) = 0.3535533905932737*As[3]; 
  data->AEM_D(0,4) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,5) = 0.3535533905932737*As[0]; 
  data->AEM_D(1,1) = 0.3535533905932737*As[2]; 
  data->AEM_D(1,3) = 0.3535533905932737*As[0]; 
  data->AEM_D(1,5) = 0.3535533905932737*As[3]; 
  data->AEM_D(2,1) = 0.3535533905932737*As[0]; 
 
  // Fill BEV_D. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,7,1) = data->u_D; 
 
} 
 
void CartFieldBinOpDivide3x3vMax_P2(binOpData_t* data, const double *A, const double *B, const short int Ncomp, const short int eqNcomp, double *out) 
{ 
  // A:       configuration space denominator field (must be a scalar field). 
  // B:       phase space numerator field (must be a scalar field). 
  // Ncomp:   number of components of B (=1 here). 
  // eqNcomp: =1 if A:numComponents=B:numComponents, =0 else (=1 here). 
  // out:     output field (same number of components as B). 
 
  // If a corner value is below zero, use cell average A.
  bool avgA = false;
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]+1.060660171779821*A[5]+1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]-0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
  if (0.7905694150420947*A[9]+0.7905694150420947*A[8]+0.7905694150420947*A[7]+1.060660171779821*A[6]-1.060660171779821*A[5]-1.060660171779821*A[4]-0.6123724356957944*A[3]-0.6123724356957944*A[2]+0.6123724356957944*A[1]+0.3535533905932737*A[0] < 0) { 
    avgA = true;
  }
 
  double As[10]; 
  double Bs[28]; 
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
    Bs[0] = B[0]; 
    Bs[1] = 0.0; 
    Bs[2] = 0.0; 
    Bs[3] = 0.0; 
    Bs[4] = B[4]; 
    Bs[5] = B[5]; 
    Bs[6] = B[6]; 
    Bs[7] = 0.0; 
    Bs[8] = 0.0; 
    Bs[9] = 0.0; 
    Bs[10] = 0.0; 
    Bs[11] = 0.0; 
    Bs[12] = 0.0; 
    Bs[13] = 0.0; 
    Bs[14] = 0.0; 
    Bs[15] = 0.0; 
    Bs[16] = B[16]; 
    Bs[17] = 0.0; 
    Bs[18] = 0.0; 
    Bs[19] = 0.0; 
    Bs[20] = B[20]; 
    Bs[21] = B[21]; 
    Bs[22] = 0.0; 
    Bs[23] = 0.0; 
    Bs[24] = 0.0; 
    Bs[25] = B[25]; 
    Bs[26] = B[26]; 
    Bs[27] = B[27]; 
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
    Bs[16] = B[16]; 
    Bs[17] = B[17]; 
    Bs[18] = B[18]; 
    Bs[19] = B[19]; 
    Bs[20] = B[20]; 
    Bs[21] = B[21]; 
    Bs[22] = B[22]; 
    Bs[23] = B[23]; 
    Bs[24] = B[24]; 
    Bs[25] = B[25]; 
    Bs[26] = B[26]; 
    Bs[27] = B[27]; 
  } 
 
  // Fill AEM_D matrix. 
  data->AEM_D(0,0) = 0.3535533905932737*As[0]; 
  data->AEM_D(0,1) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,2) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,3) = 0.3535533905932737*As[3]; 
  data->AEM_D(0,7) = 0.3535533905932737*As[4]; 
  data->AEM_D(0,8) = 0.3535533905932737*As[5]; 
  data->AEM_D(0,9) = 0.3535533905932737*As[6]; 
  data->AEM_D(0,10) = 0.3535533905932737*As[1]; 
  data->AEM_D(0,11) = 0.3162277660168379*As[7]+0.3535533905932737*As[0]; 
  data->AEM_D(0,12) = 0.3535533905932737*As[4]; 
  data->AEM_D(0,13) = 0.3535533905932737*As[5]; 
  data->AEM_D(0,17) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,18) = 0.3535533905932737*As[3]; 
  data->AEM_D(0,20) = 0.3535533905932737*As[2]; 
  data->AEM_D(0,21) = 0.3535533905932737*As[4]; 
  data->AEM_D(0,22) = 0.3162277660168379*As[8]+0.3535533905932737*As[0]; 
  data->AEM_D(0,23) = 0.3535533905932737*As[6]; 
  data->AEM_D(0,27) = 0.3535533905932737*As[1]; 
  data->AEM_D(1,1) = 0.3535533905932737*As[3]; 
  data->AEM_D(1,2) = 0.3535533905932737*As[3]; 
  data->AEM_D(1,3) = 0.3535533905932737*As[5]; 
  data->AEM_D(1,4) = 0.3535533905932737*As[6]; 
  data->AEM_D(1,5) = 0.3162277660168379*As[9]+0.3535533905932737*As[0]; 
  data->AEM_D(1,10) = 0.3535533905932737*As[1]; 
  data->AEM_D(1,11) = 0.3535533905932737*As[2]; 
  data->AEM_D(1,16) = 0.3535533905932737*As[0]; 
  data->AEM_D(1,27) = 0.3535533905932737*As[0]; 
  data->AEM_D(2,10) = 0.3535533905932737*As[0]; 
  data->AEM_D(2,14) = 0.3535533905932737*As[4]; 
  data->AEM_D(2,15) = 0.3535533905932737*As[2]; 
  data->AEM_D(2,16) = 0.3535533905932737*As[1]; 
  data->AEM_D(2,21) = 0.3162277660168379*As[8]+0.3162277660168379*As[7]+0.3535533905932737*As[0]; 
  data->AEM_D(2,22) = 0.3535533905932737*As[6]; 
  data->AEM_D(2,23) = 0.3535533905932737*As[5]; 
  data->AEM_D(2,24) = 0.3535533905932737*As[5]; 
  data->AEM_D(2,25) = 0.3535533905932737*As[3]; 
  data->AEM_D(2,27) = 0.3535533905932737*As[1]; 
  data->AEM_D(3,3) = 0.3535533905932737*As[6]; 
  data->AEM_D(3,4) = 0.3162277660168379*As[9]+0.3162277660168379*As[7]+0.3535533905932737*As[0]; 
  data->AEM_D(3,5) = 0.3535533905932737*As[4]; 
  data->AEM_D(3,6) = 0.3535533905932737*As[6]; 
  data->AEM_D(3,8) = 0.3535533905932737*As[3]; 
  data->AEM_D(3,9) = 0.3535533905932737*As[2]; 
  data->AEM_D(3,13) = 0.3535533905932737*As[5]; 
  data->AEM_D(3,14) = 0.3535533905932737*As[4]; 
  data->AEM_D(3,15) = 0.3162277660168379*As[9]+0.3162277660168379*As[8]+0.3535533905932737*As[0]; 
  data->AEM_D(3,20) = 0.3535533905932737*As[1]; 
  data->AEM_D(4,2) = 0.3535533905932737*As[2]; 
  data->AEM_D(4,12) = 0.3535533905932737*As[3]; 
  data->AEM_D(4,23) = 0.3535533905932737*As[1]; 
  data->AEM_D(5,5) = 0.3535533905932737*As[2]; 
  data->AEM_D(5,15) = 0.3535533905932737*As[3]; 
  data->AEM_D(6,8) = 0.3535533905932737*As[1]; 
  data->AEM_D(6,18) = 0.3535533905932737*As[2]; 
  data->AEM_D(7,0) = 0.3535533905932737*As[3]; 
  data->AEM_D(7,24) = 0.3535533905932737*As[7]; 
  data->AEM_D(7,25) = 0.3162277660168379*As[1]; 
  data->AEM_D(8,3) = 0.3162277660168379*As[4]; 
  data->AEM_D(8,4) = 0.3162277660168379*As[5]; 
  data->AEM_D(8,6) = 0.3535533905932737*As[8]; 
  data->AEM_D(8,8) = 0.3162277660168379*As[2]; 
  data->AEM_D(8,13) = 0.3162277660168379*As[4]; 
  data->AEM_D(8,15) = 0.3162277660168379*As[6]; 
  data->AEM_D(8,16) = 0.3535533905932737*As[9]; 
  data->AEM_D(8,19) = 0.3162277660168379*As[3]; 
  data->AEM_D(8,24) = 0.3162277660168379*As[5]; 
  data->AEM_D(8,25) = 0.3162277660168379*As[6]; 
 
  // Fill BEV_D. 
  data->BEV_D << Bs[0],Bs[1],Bs[2],Bs[3],Bs[4],Bs[5],Bs[6],Bs[7],Bs[8],Bs[9],Bs[10],Bs[11],Bs[12],Bs[13],Bs[14],Bs[15],Bs[16],Bs[17],Bs[18],Bs[19],Bs[20],Bs[21],Bs[22],Bs[23],Bs[24],Bs[25],Bs[26],Bs[27]; 
 
  // Solve the system of equations. 
  data->u_D = data->AEM_D.colPivHouseholderQr().solve(data->BEV_D); 
 
  // Copy data from Eigen vector. 
  Eigen::Map<VectorXd>(out,28,1) = data->u_D; 
 
} 
 
