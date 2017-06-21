#include <BasisSurfToVolModDecl.h>

/* Polyorder 1 */ 
void ModalSer2DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  }
}
void ModalSer2DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  }
}
void ModalSer2DP1_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  }
}
void ModalSer2DP1_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  }
}
/* Polyorder 2 */ 
void ModalSer2DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  out[4] += 1.58113883008419*f[0]*sfact; 
  out[5] += 0.7071067811865475*f[2]*sfact; 
  out[6] += 1.58113883008419*f[1]*sfact; 
  out[7] += 1.224744871391589*f[2]*sfact; 
  }
}
void ModalSer2DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  out[4] += 1.58113883008419*f[0]*sfact; 
  out[5] += 0.7071067811865475*f[2]*sfact; 
  out[6] += 1.58113883008419*f[1]*sfact; 
  out[7] += -1.224744871391589*f[2]*sfact; 
  }
}
void ModalSer2DP2_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  out[4] += 0.7071067811865475*f[2]*sfact; 
  out[5] += 1.58113883008419*f[0]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 1.58113883008419*f[1]*sfact; 
  }
}
void ModalSer2DP2_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  out[4] += 0.7071067811865475*f[2]*sfact; 
  out[5] += 1.58113883008419*f[0]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 1.58113883008419*f[1]*sfact; 
  }
}
/* Polyorder 3 */ 
void ModalSer2DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  out[4] += 1.58113883008419*f[0]*sfact; 
  out[5] += 0.7071067811865475*f[2]*sfact; 
  out[6] += 1.58113883008419*f[1]*sfact; 
  out[7] += 1.224744871391589*f[2]*sfact; 
  out[8] += 1.870828693386971*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[3]*sfact; 
  out[10] += 1.870828693386971*f[1]*sfact; 
  out[11] += 1.224744871391589*f[3]*sfact; 
  }
}
void ModalSer2DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  out[4] += 1.58113883008419*f[0]*sfact; 
  out[5] += 0.7071067811865475*f[2]*sfact; 
  out[6] += 1.58113883008419*f[1]*sfact; 
  out[7] += -1.224744871391589*f[2]*sfact; 
  out[8] += -1.870828693386971*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[3]*sfact; 
  out[10] += -1.870828693386971*f[1]*sfact; 
  out[11] += -1.224744871391589*f[3]*sfact; 
  }
}
void ModalSer2DP3_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  out[4] += 0.7071067811865475*f[2]*sfact; 
  out[5] += 1.58113883008419*f[0]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 1.58113883008419*f[1]*sfact; 
  out[8] += 0.7071067811865475*f[3]*sfact; 
  out[9] += 1.870828693386971*f[0]*sfact; 
  out[10] += 1.224744871391589*f[3]*sfact; 
  out[11] += 1.870828693386971*f[1]*sfact; 
  }
}
void ModalSer2DP3_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  out[4] += 0.7071067811865475*f[2]*sfact; 
  out[5] += 1.58113883008419*f[0]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 1.58113883008419*f[1]*sfact; 
  out[8] += 0.7071067811865475*f[3]*sfact; 
  out[9] += -1.870828693386971*f[0]*sfact; 
  out[10] += -1.224744871391589*f[3]*sfact; 
  out[11] += -1.870828693386971*f[1]*sfact; 
  }
}
/* Polyorder 4 */ 
void ModalSer2DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  out[4] += 1.58113883008419*f[0]*sfact; 
  out[5] += 0.7071067811865475*f[2]*sfact; 
  out[6] += 1.58113883008419*f[1]*sfact; 
  out[7] += 1.224744871391589*f[2]*sfact; 
  out[8] += 1.870828693386971*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[3]*sfact; 
  out[10] += 1.58113883008419*f[2]*sfact; 
  out[11] += 1.870828693386971*f[1]*sfact; 
  out[12] += 1.224744871391589*f[3]*sfact; 
  out[13] += 2.121320343559642*f[0]*sfact; 
  out[14] += 0.7071067811865475*f[4]*sfact; 
  out[15] += 2.121320343559642*f[1]*sfact; 
  out[16] += 1.224744871391589*f[4]*sfact; 
  }
}
void ModalSer2DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 0.7071067811865475*f[1]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  out[4] += 1.58113883008419*f[0]*sfact; 
  out[5] += 0.7071067811865475*f[2]*sfact; 
  out[6] += 1.58113883008419*f[1]*sfact; 
  out[7] += -1.224744871391589*f[2]*sfact; 
  out[8] += -1.870828693386971*f[0]*sfact; 
  out[9] += 0.7071067811865475*f[3]*sfact; 
  out[10] += 1.58113883008419*f[2]*sfact; 
  out[11] += -1.870828693386971*f[1]*sfact; 
  out[12] += -1.224744871391589*f[3]*sfact; 
  out[13] += 2.121320343559642*f[0]*sfact; 
  out[14] += 0.7071067811865475*f[4]*sfact; 
  out[15] += 2.121320343559642*f[1]*sfact; 
  out[16] += -1.224744871391589*f[4]*sfact; 
  }
}
void ModalSer2DP4_SurfToVol2_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += 1.224744871391589*f[0]*sfact; 
  out[3] += 1.224744871391589*f[1]*sfact; 
  out[4] += 0.7071067811865475*f[2]*sfact; 
  out[5] += 1.58113883008419*f[0]*sfact; 
  out[6] += 1.224744871391589*f[2]*sfact; 
  out[7] += 1.58113883008419*f[1]*sfact; 
  out[8] += 0.7071067811865475*f[3]*sfact; 
  out[9] += 1.870828693386971*f[0]*sfact; 
  out[10] += 1.58113883008419*f[2]*sfact; 
  out[11] += 1.224744871391589*f[3]*sfact; 
  out[12] += 1.870828693386971*f[1]*sfact; 
  out[13] += 0.7071067811865475*f[4]*sfact; 
  out[14] += 2.121320343559642*f[0]*sfact; 
  out[15] += 1.224744871391589*f[4]*sfact; 
  out[16] += 2.121320343559642*f[1]*sfact; 
  }
}
void ModalSer2DP4_SurfToVol2_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 0.7071067811865475*f[1]*sfact; 
  out[2] += -1.224744871391589*f[0]*sfact; 
  out[3] += -1.224744871391589*f[1]*sfact; 
  out[4] += 0.7071067811865475*f[2]*sfact; 
  out[5] += 1.58113883008419*f[0]*sfact; 
  out[6] += -1.224744871391589*f[2]*sfact; 
  out[7] += 1.58113883008419*f[1]*sfact; 
  out[8] += 0.7071067811865475*f[3]*sfact; 
  out[9] += -1.870828693386971*f[0]*sfact; 
  out[10] += 1.58113883008419*f[2]*sfact; 
  out[11] += -1.224744871391589*f[3]*sfact; 
  out[12] += -1.870828693386971*f[1]*sfact; 
  out[13] += 0.7071067811865475*f[4]*sfact; 
  out[14] += 2.121320343559642*f[0]*sfact; 
  out[15] += -1.224744871391589*f[4]*sfact; 
  out[16] += 2.121320343559642*f[1]*sfact; 
  }
}
