void LinearHyperMax2DP1_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  }
} 
void LinearHyperMax2DP2_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[3] += 1.732050807568877*f[2]*vfact; 
  out[4] += 3.872983346207417*f[1]*vfact; 
  }
} 
void LinearHyperMax2DP3_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[3] += 1.732050807568877*f[2]*vfact; 
  out[4] += 3.872983346207417*f[1]*vfact; 
  out[6] += 3.872983346207417*f[3]*vfact; 
  out[7] += 1.732050807568877*f[5]*vfact; 
  out[8] += 5.916079783099617*f[4]*vfact+2.645751311064591*f[0]*vfact; 
  }
} 
void LinearHyperMax2DP4_Vol1(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[1] += 1.732050807568877*f[0]*vfact; 
  out[3] += 1.732050807568877*f[2]*vfact; 
  out[4] += 3.872983346207417*f[1]*vfact; 
  out[6] += 3.872983346207417*f[3]*vfact; 
  out[7] += 1.732050807568877*f[5]*vfact; 
  out[8] += 5.916079783099617*f[4]*vfact+2.645751311064591*f[0]*vfact; 
  out[10] += 3.872983346207417*f[7]*vfact; 
  out[11] += 5.916079783099617*f[6]*vfact+2.645751311064591*f[2]*vfact; 
  out[12] += 1.732050807568877*f[9]*vfact; 
  out[13] += 7.937253933193772*f[8]*vfact+5.196152422706631*f[1]*vfact; 
  }
} 
void LinearHyperMax2DP1_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  }
} 
void LinearHyperMax2DP2_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  out[3] += 1.732050807568877*f[1]*vfact; 
  out[5] += 3.872983346207417*f[2]*vfact; 
  }
} 
void LinearHyperMax2DP3_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  out[3] += 1.732050807568877*f[1]*vfact; 
  out[5] += 3.872983346207417*f[2]*vfact; 
  out[6] += 1.732050807568877*f[4]*vfact; 
  out[7] += 3.872983346207417*f[3]*vfact; 
  out[9] += 5.916079783099617*f[5]*vfact+2.645751311064591*f[0]*vfact; 
  }
} 
void LinearHyperMax2DP4_Vol2(int meqn, int mbasis, double vfact, const double *fIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &fIn[mbasis*m]; 
  double *out = &volOut[mbasis*m]; 
  out[2] += 1.732050807568877*f[0]*vfact; 
  out[3] += 1.732050807568877*f[1]*vfact; 
  out[5] += 3.872983346207417*f[2]*vfact; 
  out[6] += 1.732050807568877*f[4]*vfact; 
  out[7] += 3.872983346207417*f[3]*vfact; 
  out[9] += 5.916079783099617*f[5]*vfact+2.645751311064591*f[0]*vfact; 
  out[10] += 3.872983346207417*f[6]*vfact; 
  out[11] += 1.732050807568877*f[8]*vfact; 
  out[12] += 5.916079783099617*f[7]*vfact+2.645751311064591*f[1]*vfact; 
  out[14] += 7.937253933193772*f[9]*vfact+5.196152422706631*f[2]*vfact; 
  }
} 
