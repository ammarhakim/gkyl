#include <LinearHyperModDecl.h>

void LinearHyperMax3DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  }
} 
void LinearHyperMax3DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[2]*vfact; 
  out[5] += 1.732050807568877*f[3]*vfact; 
  out[7] += 3.872983346207417*f[1]*vfact; 
  }
} 
void LinearHyperMax3DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[2]*vfact; 
  out[5] += 1.732050807568877*f[3]*vfact; 
  out[7] += 3.872983346207417*f[1]*vfact; 
  out[10] += 1.732050807568877*f[6]*vfact; 
  out[11] += 3.872983346207417*f[4]*vfact; 
  out[12] += 1.732050807568877*f[8]*vfact; 
  out[13] += 3.872983346207417*f[5]*vfact; 
  out[15] += 1.732050807568877*f[9]*vfact; 
  out[17] += 5.916079783099617*f[7]*vfact+2.645751311064591*f[0]*vfact; 
  }
} 
void LinearHyperMax3DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[2]*vfact; 
  out[5] += 1.732050807568877*f[3]*vfact; 
  out[7] += 3.872983346207417*f[1]*vfact; 
  out[10] += 1.732050807568877*f[6]*vfact; 
  out[11] += 3.872983346207417*f[4]*vfact; 
  out[12] += 1.732050807568877*f[8]*vfact; 
  out[13] += 3.872983346207417*f[5]*vfact; 
  out[15] += 1.732050807568877*f[9]*vfact; 
  out[17] += 5.916079783099617*f[7]*vfact+2.645751311064591*f[0]*vfact; 
  out[20] += 3.872983346207417*f[10]*vfact; 
  out[21] += 1.732050807568877*f[14]*vfact; 
  out[22] += 1.732050807568877*f[16]*vfact; 
  out[23] += 3.872983346207417*f[12]*vfact; 
  out[24] += 3.872983346207417*f[15]*vfact; 
  out[26] += 5.916079783099617*f[11]*vfact+2.645751311064591*f[2]*vfact; 
  out[27] += 1.732050807568877*f[18]*vfact; 
  out[28] += 5.916079783099617*f[13]*vfact+2.645751311064591*f[3]*vfact; 
  out[30] += 1.732050807568877*f[19]*vfact; 
  out[32] += 7.937253933193772*f[17]*vfact+5.196152422706631*f[1]*vfact; 
  }
} 
void LinearHyperMax3DP1_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  }
} 
void LinearHyperMax3DP2_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[3]*vfact; 
  out[8] += 3.872983346207417*f[2]*vfact; 
  }
} 
void LinearHyperMax3DP3_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[3]*vfact; 
  out[8] += 3.872983346207417*f[2]*vfact; 
  out[10] += 1.732050807568877*f[5]*vfact; 
  out[11] += 1.732050807568877*f[7]*vfact; 
  out[12] += 3.872983346207417*f[4]*vfact; 
  out[14] += 3.872983346207417*f[6]*vfact; 
  out[16] += 1.732050807568877*f[9]*vfact; 
  out[18] += 5.916079783099617*f[8]*vfact+2.645751311064591*f[0]*vfact; 
  }
} 
void LinearHyperMax3DP4_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[3]*vfact; 
  out[8] += 3.872983346207417*f[2]*vfact; 
  out[10] += 1.732050807568877*f[5]*vfact; 
  out[11] += 1.732050807568877*f[7]*vfact; 
  out[12] += 3.872983346207417*f[4]*vfact; 
  out[14] += 3.872983346207417*f[6]*vfact; 
  out[16] += 1.732050807568877*f[9]*vfact; 
  out[18] += 5.916079783099617*f[8]*vfact+2.645751311064591*f[0]*vfact; 
  out[20] += 1.732050807568877*f[13]*vfact; 
  out[21] += 3.872983346207417*f[10]*vfact; 
  out[22] += 1.732050807568877*f[15]*vfact; 
  out[23] += 3.872983346207417*f[11]*vfact; 
  out[25] += 3.872983346207417*f[16]*vfact; 
  out[26] += 1.732050807568877*f[17]*vfact; 
  out[27] += 5.916079783099617*f[12]*vfact+2.645751311064591*f[1]*vfact; 
  out[29] += 5.916079783099617*f[14]*vfact+2.645751311064591*f[3]*vfact; 
  out[31] += 1.732050807568877*f[19]*vfact; 
  out[33] += 7.937253933193772*f[18]*vfact+5.196152422706631*f[2]*vfact; 
  }
} 
void LinearHyperMax3DP1_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[3] += 1.732050807568877*f[0]*vfact; 
  }
} 
void LinearHyperMax3DP2_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[3] += 1.732050807568877*f[0]*vfact; 
  out[5] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[2]*vfact; 
  out[9] += 3.872983346207417*f[3]*vfact; 
  }
} 
void LinearHyperMax3DP3_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[3] += 1.732050807568877*f[0]*vfact; 
  out[5] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[2]*vfact; 
  out[9] += 3.872983346207417*f[3]*vfact; 
  out[10] += 1.732050807568877*f[4]*vfact; 
  out[13] += 1.732050807568877*f[7]*vfact; 
  out[14] += 1.732050807568877*f[8]*vfact; 
  out[15] += 3.872983346207417*f[5]*vfact; 
  out[16] += 3.872983346207417*f[6]*vfact; 
  out[19] += 5.916079783099617*f[9]*vfact+2.645751311064591*f[0]*vfact; 
  }
} 
void LinearHyperMax3DP4_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[3] += 1.732050807568877*f[0]*vfact; 
  out[5] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[2]*vfact; 
  out[9] += 3.872983346207417*f[3]*vfact; 
  out[10] += 1.732050807568877*f[4]*vfact; 
  out[13] += 1.732050807568877*f[7]*vfact; 
  out[14] += 1.732050807568877*f[8]*vfact; 
  out[15] += 3.872983346207417*f[5]*vfact; 
  out[16] += 3.872983346207417*f[6]*vfact; 
  out[19] += 5.916079783099617*f[9]*vfact+2.645751311064591*f[0]*vfact; 
  out[20] += 1.732050807568877*f[11]*vfact; 
  out[21] += 1.732050807568877*f[12]*vfact; 
  out[22] += 3.872983346207417*f[10]*vfact; 
  out[24] += 3.872983346207417*f[13]*vfact; 
  out[25] += 3.872983346207417*f[14]*vfact; 
  out[28] += 1.732050807568877*f[17]*vfact; 
  out[29] += 1.732050807568877*f[18]*vfact; 
  out[30] += 5.916079783099617*f[15]*vfact+2.645751311064591*f[1]*vfact; 
  out[31] += 5.916079783099617*f[16]*vfact+2.645751311064591*f[2]*vfact; 
  out[34] += 7.937253933193772*f[19]*vfact+5.196152422706631*f[3]*vfact; 
  }
} 
