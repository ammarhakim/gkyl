void LinearHyperSer3DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[2]*vfact; 
  out[5] += 1.732050807568877*f[3]*vfact; 
  out[7] += 1.732050807568877*f[6]*vfact; 
  }
} 
void LinearHyperSer3DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[17] += 3.872983346207417*f[10]*vfact; 
  out[18] += 1.732050807568877*f[14]*vfact; 
  out[19] += 1.732050807568877*f[16]*vfact; 
  }
} 
void LinearHyperSer3DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[23] += 5.916079783099617*f[11]*vfact+2.645751311064591*f[2]*vfact; 
  out[24] += 1.732050807568877*f[18]*vfact; 
  out[25] += 5.916079783099617*f[13]*vfact+2.645751311064591*f[3]*vfact; 
  out[27] += 1.732050807568877*f[19]*vfact; 
  out[29] += 5.916079783099617*f[20]*vfact+2.645751311064591*f[6]*vfact; 
  out[30] += 1.732050807568877*f[26]*vfact; 
  out[31] += 1.732050807568877*f[28]*vfact; 
  }
} 
void LinearHyperSer3DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[35] += 3.872983346207417*f[21]*vfact; 
  out[36] += 3.872983346207417*f[22]*vfact; 
  out[37] += 1.732050807568877*f[25]*vfact; 
  out[38] += 5.916079783099617*f[20]*vfact+2.645751311064591*f[6]*vfact; 
  out[39] += 1.732050807568877*f[29]*vfact; 
  out[40] += 1.732050807568877*f[31]*vfact; 
  out[41] += 7.937253933193772*f[26]*vfact+5.196152422706631*f[4]*vfact; 
  out[42] += 1.732050807568877*f[33]*vfact; 
  out[43] += 7.937253933193772*f[28]*vfact+5.196152422706631*f[5]*vfact; 
  out[45] += 1.732050807568877*f[34]*vfact; 
  out[47] += 7.937253933193772*f[38]*vfact+5.196152422706631*f[10]*vfact; 
  out[48] += 1.732050807568877*f[44]*vfact; 
  out[49] += 1.732050807568877*f[46]*vfact; 
  }
} 
void LinearHyperSer3DP1_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  out[4] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[3]*vfact; 
  out[7] += 1.732050807568877*f[5]*vfact; 
  }
} 
void LinearHyperSer3DP2_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[17] += 1.732050807568877*f[13]*vfact; 
  out[18] += 3.872983346207417*f[10]*vfact; 
  out[19] += 1.732050807568877*f[15]*vfact; 
  }
} 
void LinearHyperSer3DP3_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[23] += 1.732050807568877*f[17]*vfact; 
  out[24] += 5.916079783099617*f[12]*vfact+2.645751311064591*f[1]*vfact; 
  out[26] += 5.916079783099617*f[14]*vfact+2.645751311064591*f[3]*vfact; 
  out[28] += 1.732050807568877*f[19]*vfact; 
  out[29] += 1.732050807568877*f[25]*vfact; 
  out[30] += 5.916079783099617*f[21]*vfact+2.645751311064591*f[5]*vfact; 
  out[31] += 1.732050807568877*f[27]*vfact; 
  }
} 
void LinearHyperSer3DP4_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[35] += 3.872983346207417*f[20]*vfact; 
  out[36] += 1.732050807568877*f[24]*vfact; 
  out[37] += 3.872983346207417*f[22]*vfact; 
  out[38] += 1.732050807568877*f[28]*vfact; 
  out[39] += 5.916079783099617*f[21]*vfact+2.645751311064591*f[5]*vfact; 
  out[40] += 1.732050807568877*f[30]*vfact; 
  out[41] += 1.732050807568877*f[32]*vfact; 
  out[42] += 7.937253933193772*f[27]*vfact+5.196152422706631*f[4]*vfact; 
  out[44] += 7.937253933193772*f[29]*vfact+5.196152422706631*f[6]*vfact; 
  out[46] += 1.732050807568877*f[34]*vfact; 
  out[47] += 1.732050807568877*f[43]*vfact; 
  out[48] += 7.937253933193772*f[39]*vfact+5.196152422706631*f[10]*vfact; 
  out[49] += 1.732050807568877*f[45]*vfact; 
  }
} 
void LinearHyperSer3DP1_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[3] += 1.732050807568877*f[0]*vfact; 
  out[5] += 1.732050807568877*f[1]*vfact; 
  out[6] += 1.732050807568877*f[2]*vfact; 
  out[7] += 1.732050807568877*f[4]*vfact; 
  }
} 
void LinearHyperSer3DP2_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[17] += 1.732050807568877*f[11]*vfact; 
  out[18] += 1.732050807568877*f[12]*vfact; 
  out[19] += 3.872983346207417*f[10]*vfact; 
  }
} 
void LinearHyperSer3DP3_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[25] += 1.732050807568877*f[17]*vfact; 
  out[26] += 1.732050807568877*f[18]*vfact; 
  out[27] += 5.916079783099617*f[15]*vfact+2.645751311064591*f[1]*vfact; 
  out[28] += 5.916079783099617*f[16]*vfact+2.645751311064591*f[2]*vfact; 
  out[29] += 1.732050807568877*f[23]*vfact; 
  out[30] += 1.732050807568877*f[24]*vfact; 
  out[31] += 5.916079783099617*f[22]*vfact+2.645751311064591*f[4]*vfact; 
  }
} 
void LinearHyperSer3DP4_Vol3(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
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
  out[35] += 1.732050807568877*f[23]*vfact; 
  out[36] += 3.872983346207417*f[20]*vfact; 
  out[37] += 3.872983346207417*f[21]*vfact; 
  out[38] += 1.732050807568877*f[26]*vfact; 
  out[39] += 1.732050807568877*f[27]*vfact; 
  out[40] += 5.916079783099617*f[22]*vfact+2.645751311064591*f[4]*vfact; 
  out[43] += 1.732050807568877*f[32]*vfact; 
  out[44] += 1.732050807568877*f[33]*vfact; 
  out[45] += 7.937253933193772*f[30]*vfact+5.196152422706631*f[5]*vfact; 
  out[46] += 7.937253933193772*f[31]*vfact+5.196152422706631*f[6]*vfact; 
  out[47] += 1.732050807568877*f[41]*vfact; 
  out[48] += 1.732050807568877*f[42]*vfact; 
  out[49] += 7.937253933193772*f[40]*vfact+5.196152422706631*f[10]*vfact; 
  }
} 
