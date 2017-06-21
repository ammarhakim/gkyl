#include <BasisVolToSurfModDecl.h> 
// polyOrder 1 
//    dir 1 
void ModalSerendipBasisSurf1DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
  } 
} 
void ModalSerendipBasisSurf1DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 0.7071067811865475*_in[0]-1.224744871391589*_in[1]; 
  } 
} 


// polyOrder 2 
//    dir 1 
void ModalSerendipBasisSurf1DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.58113883008419*_in[2]+1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
  } 
} 
void ModalSerendipBasisSurf1DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.58113883008419*_in[2]-1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
  } 
} 


// polyOrder 3 
//    dir 1 
void ModalSerendipBasisSurf1DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.870828693386971*_in[3]+1.58113883008419*_in[2]+1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
  } 
} 
void ModalSerendipBasisSurf1DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = (-1.870828693386971*_in[3])+1.58113883008419*_in[2]-1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
  } 
} 


// polyOrder 4 
//    dir 1 
void ModalSerendipBasisSurf1DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 2.121320343559642*_in[4]+1.870828693386971*_in[3]+1.58113883008419*_in[2]+1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
  } 
} 
void ModalSerendipBasisSurf1DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 2.121320343559642*_in[4]-1.870828693386971*_in[3]+1.58113883008419*_in[2]-1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
  } 
} 


