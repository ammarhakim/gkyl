/* Polyorder 1 */ 
void ModalMax3DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  }
}
void ModalMax3DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  }
}
void ModalMax3DP1_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  }
}
void ModalMax3DP1_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  }
}
void ModalMax3DP1_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += 1.224744871391589*f[0]*sfact; 
  }
}
void ModalMax3DP1_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += -1.224744871391589*f[0]*sfact; 
  }
}
/* Polyorder 2 */ 
void ModalMax3DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += 1.224744871391589*f[1]*sfact; 
  out[5] += 1.224744871391589*f[2]*sfact; 
  out[6] += 0.7071067811865475*f[3]*sfact; 
  out[7] += 1.58113883008419*f[0]*sfact; 
  out[8] += 0.7071067811865475*f[4]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  }
}
void ModalMax3DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += -1.224744871391589*f[1]*sfact; 
  out[5] += -1.224744871391589*f[2]*sfact; 
  out[6] += 0.7071067811865475*f[3]*sfact; 
  out[7] += 1.58113883008419*f[0]*sfact; 
  out[8] += 0.7071067811865475*f[4]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  }
}
void ModalMax3DP2_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += 1.224744871391589*f[1]*sfact; 
  out[5] += 0.7071067811865475*f[3]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 1.58113883008419*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  }
}
void ModalMax3DP2_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += -1.224744871391589*f[1]*sfact; 
  out[5] += 0.7071067811865475*f[3]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 1.58113883008419*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  }
}
void ModalMax3DP2_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += 1.224744871391589*f[0]*sfact; 
  out[4] += 0.7071067811865475*f[3]*sfact; 
  out[5] += 1.224744871391589*f[1]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 0.7071067811865475*f[5]*sfact; 
  out[9] += 1.58113883008419*f[0]*sfact; 
  }
}
void ModalMax3DP2_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += -1.224744871391589*f[0]*sfact; 
  out[4] += 0.7071067811865475*f[3]*sfact; 
  out[5] += -1.224744871391589*f[1]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 0.7071067811865475*f[5]*sfact; 
  out[9] += 1.58113883008419*f[0]*sfact; 
  }
}
/* Polyorder 3 */ 
void ModalMax3DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += 1.224744871391589*f[1]*sfact; 
  out[5] += 1.224744871391589*f[2]*sfact; 
  out[6] += 0.7071067811865475*f[3]*sfact; 
  out[7] += 1.58113883008419*f[0]*sfact; 
  out[8] += 0.7071067811865475*f[4]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += 1.224744871391589*f[3]*sfact; 
  out[11] += 1.58113883008419*f[1]*sfact; 
  out[12] += 1.224744871391589*f[4]*sfact; 
  out[13] += 1.58113883008419*f[2]*sfact; 
  out[14] += 0.7071067811865475*f[6]*sfact; 
  out[15] += 1.224744871391589*f[5]*sfact; 
  out[16] += 0.7071067811865475*f[7]*sfact; 
  out[17] += 1.870828693386971*f[0]*sfact; 
  out[18] += 0.7071067811865475*f[8]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  }
}
void ModalMax3DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += -1.224744871391589*f[1]*sfact; 
  out[5] += -1.224744871391589*f[2]*sfact; 
  out[6] += 0.7071067811865475*f[3]*sfact; 
  out[7] += 1.58113883008419*f[0]*sfact; 
  out[8] += 0.7071067811865475*f[4]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += -1.224744871391589*f[3]*sfact; 
  out[11] += 1.58113883008419*f[1]*sfact; 
  out[12] += -1.224744871391589*f[4]*sfact; 
  out[13] += 1.58113883008419*f[2]*sfact; 
  out[14] += 0.7071067811865475*f[6]*sfact; 
  out[15] += -1.224744871391589*f[5]*sfact; 
  out[16] += 0.7071067811865475*f[7]*sfact; 
  out[17] += -1.870828693386971*f[0]*sfact; 
  out[18] += 0.7071067811865475*f[8]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  }
}
void ModalMax3DP3_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += 1.224744871391589*f[1]*sfact; 
  out[5] += 0.7071067811865475*f[3]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 1.58113883008419*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += 1.224744871391589*f[3]*sfact; 
  out[11] += 1.224744871391589*f[4]*sfact; 
  out[12] += 1.58113883008419*f[1]*sfact; 
  out[13] += 0.7071067811865475*f[6]*sfact; 
  out[14] += 1.58113883008419*f[2]*sfact; 
  out[15] += 0.7071067811865475*f[7]*sfact; 
  out[16] += 1.224744871391589*f[5]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += 1.870828693386971*f[0]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  }
}
void ModalMax3DP3_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += -1.224744871391589*f[1]*sfact; 
  out[5] += 0.7071067811865475*f[3]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 1.58113883008419*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += -1.224744871391589*f[3]*sfact; 
  out[11] += -1.224744871391589*f[4]*sfact; 
  out[12] += 1.58113883008419*f[1]*sfact; 
  out[13] += 0.7071067811865475*f[6]*sfact; 
  out[14] += 1.58113883008419*f[2]*sfact; 
  out[15] += 0.7071067811865475*f[7]*sfact; 
  out[16] += -1.224744871391589*f[5]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += -1.870828693386971*f[0]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  }
}
void ModalMax3DP3_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += 1.224744871391589*f[0]*sfact; 
  out[4] += 0.7071067811865475*f[3]*sfact; 
  out[5] += 1.224744871391589*f[1]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 0.7071067811865475*f[5]*sfact; 
  out[9] += 1.58113883008419*f[0]*sfact; 
  out[10] += 1.224744871391589*f[3]*sfact; 
  out[11] += 0.7071067811865475*f[6]*sfact; 
  out[12] += 0.7071067811865475*f[7]*sfact; 
  out[13] += 1.224744871391589*f[4]*sfact; 
  out[14] += 1.224744871391589*f[5]*sfact; 
  out[15] += 1.58113883008419*f[1]*sfact; 
  out[16] += 1.58113883008419*f[2]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += 0.7071067811865475*f[9]*sfact; 
  out[19] += 1.870828693386971*f[0]*sfact; 
  }
}
void ModalMax3DP3_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += -1.224744871391589*f[0]*sfact; 
  out[4] += 0.7071067811865475*f[3]*sfact; 
  out[5] += -1.224744871391589*f[1]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 0.7071067811865475*f[5]*sfact; 
  out[9] += 1.58113883008419*f[0]*sfact; 
  out[10] += -1.224744871391589*f[3]*sfact; 
  out[11] += 0.7071067811865475*f[6]*sfact; 
  out[12] += 0.7071067811865475*f[7]*sfact; 
  out[13] += -1.224744871391589*f[4]*sfact; 
  out[14] += -1.224744871391589*f[5]*sfact; 
  out[15] += 1.58113883008419*f[1]*sfact; 
  out[16] += 1.58113883008419*f[2]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += 0.7071067811865475*f[9]*sfact; 
  out[19] += -1.870828693386971*f[0]*sfact; 
  }
}
/* Polyorder 4 */ 
void ModalMax3DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += 1.224744871391589*f[1]*sfact; 
  out[5] += 1.224744871391589*f[2]*sfact; 
  out[6] += 0.7071067811865475*f[3]*sfact; 
  out[7] += 1.58113883008419*f[0]*sfact; 
  out[8] += 0.7071067811865475*f[4]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += 1.224744871391589*f[3]*sfact; 
  out[11] += 1.58113883008419*f[1]*sfact; 
  out[12] += 1.224744871391589*f[4]*sfact; 
  out[13] += 1.58113883008419*f[2]*sfact; 
  out[14] += 0.7071067811865475*f[6]*sfact; 
  out[15] += 1.224744871391589*f[5]*sfact; 
  out[16] += 0.7071067811865475*f[7]*sfact; 
  out[17] += 1.870828693386971*f[0]*sfact; 
  out[18] += 0.7071067811865475*f[8]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  out[20] += 1.58113883008419*f[3]*sfact; 
  out[21] += 1.224744871391589*f[6]*sfact; 
  out[22] += 1.224744871391589*f[7]*sfact; 
  out[23] += 1.58113883008419*f[4]*sfact; 
  out[24] += 1.58113883008419*f[5]*sfact; 
  out[25] += 0.7071067811865475*f[10]*sfact; 
  out[26] += 1.870828693386971*f[1]*sfact; 
  out[27] += 1.224744871391589*f[8]*sfact; 
  out[28] += 1.870828693386971*f[2]*sfact; 
  out[29] += 0.7071067811865475*f[11]*sfact; 
  out[30] += 1.224744871391589*f[9]*sfact; 
  out[31] += 0.7071067811865475*f[12]*sfact; 
  out[32] += 2.121320343559642*f[0]*sfact; 
  out[33] += 0.7071067811865475*f[13]*sfact; 
  out[34] += 0.7071067811865475*f[14]*sfact; 
  }
}
void ModalMax3DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += -1.224744871391589*f[1]*sfact; 
  out[5] += -1.224744871391589*f[2]*sfact; 
  out[6] += 0.7071067811865475*f[3]*sfact; 
  out[7] += 1.58113883008419*f[0]*sfact; 
  out[8] += 0.7071067811865475*f[4]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += -1.224744871391589*f[3]*sfact; 
  out[11] += 1.58113883008419*f[1]*sfact; 
  out[12] += -1.224744871391589*f[4]*sfact; 
  out[13] += 1.58113883008419*f[2]*sfact; 
  out[14] += 0.7071067811865475*f[6]*sfact; 
  out[15] += -1.224744871391589*f[5]*sfact; 
  out[16] += 0.7071067811865475*f[7]*sfact; 
  out[17] += -1.870828693386971*f[0]*sfact; 
  out[18] += 0.7071067811865475*f[8]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  out[20] += 1.58113883008419*f[3]*sfact; 
  out[21] += -1.224744871391589*f[6]*sfact; 
  out[22] += -1.224744871391589*f[7]*sfact; 
  out[23] += 1.58113883008419*f[4]*sfact; 
  out[24] += 1.58113883008419*f[5]*sfact; 
  out[25] += 0.7071067811865475*f[10]*sfact; 
  out[26] += -1.870828693386971*f[1]*sfact; 
  out[27] += -1.224744871391589*f[8]*sfact; 
  out[28] += -1.870828693386971*f[2]*sfact; 
  out[29] += 0.7071067811865475*f[11]*sfact; 
  out[30] += -1.224744871391589*f[9]*sfact; 
  out[31] += 0.7071067811865475*f[12]*sfact; 
  out[32] += 2.121320343559642*f[0]*sfact; 
  out[33] += 0.7071067811865475*f[13]*sfact; 
  out[34] += 0.7071067811865475*f[14]*sfact; 
  }
}
void ModalMax3DP4_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += 1.224744871391589*f[1]*sfact; 
  out[5] += 0.7071067811865475*f[3]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 1.58113883008419*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += 1.224744871391589*f[3]*sfact; 
  out[11] += 1.224744871391589*f[4]*sfact; 
  out[12] += 1.58113883008419*f[1]*sfact; 
  out[13] += 0.7071067811865475*f[6]*sfact; 
  out[14] += 1.58113883008419*f[2]*sfact; 
  out[15] += 0.7071067811865475*f[7]*sfact; 
  out[16] += 1.224744871391589*f[5]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += 1.870828693386971*f[0]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  out[20] += 1.224744871391589*f[6]*sfact; 
  out[21] += 1.58113883008419*f[3]*sfact; 
  out[22] += 1.224744871391589*f[7]*sfact; 
  out[23] += 1.58113883008419*f[4]*sfact; 
  out[24] += 0.7071067811865475*f[10]*sfact; 
  out[25] += 1.58113883008419*f[5]*sfact; 
  out[26] += 1.224744871391589*f[8]*sfact; 
  out[27] += 1.870828693386971*f[1]*sfact; 
  out[28] += 0.7071067811865475*f[11]*sfact; 
  out[29] += 1.870828693386971*f[2]*sfact; 
  out[30] += 0.7071067811865475*f[12]*sfact; 
  out[31] += 1.224744871391589*f[9]*sfact; 
  out[32] += 0.7071067811865475*f[13]*sfact; 
  out[33] += 2.121320343559642*f[0]*sfact; 
  out[34] += 0.7071067811865475*f[14]*sfact; 
  }
}
void ModalMax3DP4_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += 0.7071067811865475*f[2]*sfact; 
  out[4] += -1.224744871391589*f[1]*sfact; 
  out[5] += 0.7071067811865475*f[3]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 1.58113883008419*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[5]*sfact; 
  out[10] += -1.224744871391589*f[3]*sfact; 
  out[11] += -1.224744871391589*f[4]*sfact; 
  out[12] += 1.58113883008419*f[1]*sfact; 
  out[13] += 0.7071067811865475*f[6]*sfact; 
  out[14] += 1.58113883008419*f[2]*sfact; 
  out[15] += 0.7071067811865475*f[7]*sfact; 
  out[16] += -1.224744871391589*f[5]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += -1.870828693386971*f[0]*sfact; 
  out[19] += 0.7071067811865475*f[9]*sfact; 
  out[20] += -1.224744871391589*f[6]*sfact; 
  out[21] += 1.58113883008419*f[3]*sfact; 
  out[22] += -1.224744871391589*f[7]*sfact; 
  out[23] += 1.58113883008419*f[4]*sfact; 
  out[24] += 0.7071067811865475*f[10]*sfact; 
  out[25] += 1.58113883008419*f[5]*sfact; 
  out[26] += -1.224744871391589*f[8]*sfact; 
  out[27] += -1.870828693386971*f[1]*sfact; 
  out[28] += 0.7071067811865475*f[11]*sfact; 
  out[29] += -1.870828693386971*f[2]*sfact; 
  out[30] += 0.7071067811865475*f[12]*sfact; 
  out[31] += -1.224744871391589*f[9]*sfact; 
  out[32] += 0.7071067811865475*f[13]*sfact; 
  out[33] += 2.121320343559642*f[0]*sfact; 
  out[34] += 0.7071067811865475*f[14]*sfact; 
  }
}
void ModalMax3DP4_SurfToVol3_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += 1.224744871391589*f[0]*sfact; 
  out[4] += 0.7071067811865475*f[3]*sfact; 
  out[5] += 1.224744871391589*f[1]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 0.7071067811865475*f[5]*sfact; 
  out[9] += 1.58113883008419*f[0]*sfact; 
  out[10] += 1.224744871391589*f[3]*sfact; 
  out[11] += 0.7071067811865475*f[6]*sfact; 
  out[12] += 0.7071067811865475*f[7]*sfact; 
  out[13] += 1.224744871391589*f[4]*sfact; 
  out[14] += 1.224744871391589*f[5]*sfact; 
  out[15] += 1.58113883008419*f[1]*sfact; 
  out[16] += 1.58113883008419*f[2]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += 0.7071067811865475*f[9]*sfact; 
  out[19] += 1.870828693386971*f[0]*sfact; 
  out[20] += 1.224744871391589*f[6]*sfact; 
  out[21] += 1.224744871391589*f[7]*sfact; 
  out[22] += 1.58113883008419*f[3]*sfact; 
  out[23] += 0.7071067811865475*f[10]*sfact; 
  out[24] += 1.58113883008419*f[4]*sfact; 
  out[25] += 1.58113883008419*f[5]*sfact; 
  out[26] += 0.7071067811865475*f[11]*sfact; 
  out[27] += 0.7071067811865475*f[12]*sfact; 
  out[28] += 1.224744871391589*f[8]*sfact; 
  out[29] += 1.224744871391589*f[9]*sfact; 
  out[30] += 1.870828693386971*f[1]*sfact; 
  out[31] += 1.870828693386971*f[2]*sfact; 
  out[32] += 0.7071067811865475*f[13]*sfact; 
  out[33] += 0.7071067811865475*f[14]*sfact; 
  out[34] += 2.121320343559642*f[0]*sfact; 
  }
}
void ModalMax3DP4_SurfToVol3_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 0.7071067811865475*f[2]*sfact; 
  out[3] += -1.224744871391589*f[0]*sfact; 
  out[4] += 0.7071067811865475*f[3]*sfact; 
  out[5] += -1.224744871391589*f[1]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 0.7071067811865475*f[4]*sfact; 
  out[8] += 0.7071067811865475*f[5]*sfact; 
  out[9] += 1.58113883008419*f[0]*sfact; 
  out[10] += -1.224744871391589*f[3]*sfact; 
  out[11] += 0.7071067811865475*f[6]*sfact; 
  out[12] += 0.7071067811865475*f[7]*sfact; 
  out[13] += -1.224744871391589*f[4]*sfact; 
  out[14] += -1.224744871391589*f[5]*sfact; 
  out[15] += 1.58113883008419*f[1]*sfact; 
  out[16] += 1.58113883008419*f[2]*sfact; 
  out[17] += 0.7071067811865475*f[8]*sfact; 
  out[18] += 0.7071067811865475*f[9]*sfact; 
  out[19] += -1.870828693386971*f[0]*sfact; 
  out[20] += -1.224744871391589*f[6]*sfact; 
  out[21] += -1.224744871391589*f[7]*sfact; 
  out[22] += 1.58113883008419*f[3]*sfact; 
  out[23] += 0.7071067811865475*f[10]*sfact; 
  out[24] += 1.58113883008419*f[4]*sfact; 
  out[25] += 1.58113883008419*f[5]*sfact; 
  out[26] += 0.7071067811865475*f[11]*sfact; 
  out[27] += 0.7071067811865475*f[12]*sfact; 
  out[28] += -1.224744871391589*f[8]*sfact; 
  out[29] += -1.224744871391589*f[9]*sfact; 
  out[30] += -1.870828693386971*f[1]*sfact; 
  out[31] += -1.870828693386971*f[2]*sfact; 
  out[32] += 0.7071067811865475*f[13]*sfact; 
  out[33] += 0.7071067811865475*f[14]*sfact; 
  out[34] += 2.121320343559642*f[0]*sfact; 
  }
}
