#include <BasisSurfToVolModDecl.h>

/* Polyorder 1 */ 
void ModalMax1DP1_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  }
}
void ModalMax1DP1_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  }
}
/* Polyorder 2 */ 
void ModalMax1DP2_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 1.58113883008419*f[0]*sfact; 
  }
}
void ModalMax1DP2_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 1.58113883008419*f[0]*sfact; 
  }
}
/* Polyorder 3 */ 
void ModalMax1DP3_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 1.58113883008419*f[0]*sfact; 
  out[3] += 1.870828693386971*f[0]*sfact; 
  }
}
void ModalMax1DP3_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 1.58113883008419*f[0]*sfact; 
  out[3] += -1.870828693386971*f[0]*sfact; 
  }
}
/* Polyorder 4 */ 
void ModalMax1DP4_SurfToVol1_Left(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += 1.224744871391589*f[0]*sfact; 
  out[2] += 1.58113883008419*f[0]*sfact; 
  out[3] += 1.870828693386971*f[0]*sfact; 
  out[4] += 2.121320343559642*f[0]*sfact; 
  }
}
void ModalMax1DP4_SurfToVol1_Right(int meqn, int mbasis, int msurf, double sfact, const double *surfIn, double *volOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  {
  const double *f = &surfIn[msurf*m]; 
  double *out = &volOut[mbasis*m]; 
  out[0] += 0.7071067811865475*f[0]*sfact; 
  out[1] += -1.224744871391589*f[0]*sfact; 
  out[2] += 1.58113883008419*f[0]*sfact; 
  out[3] += -1.870828693386971*f[0]*sfact; 
  out[4] += 2.121320343559642*f[0]*sfact; 
  }
}
