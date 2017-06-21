#include <BasisVolToSurfModDecl.h> 
// polyOrder 1 
//    dir 1 
void ModalSerendipBasisSurf2DP1_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
   _out[1] = 1.224744871391589*_in[3]+0.7071067811865475*_in[2]; 
  } 
} 
void ModalSerendipBasisSurf2DP1_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 0.7071067811865475*_in[0]-1.224744871391589*_in[1]; 
   _out[1] = 0.7071067811865475*_in[2]-1.224744871391589*_in[3]; 
  } 
} 

//    dir 2 
void ModalSerendipBasisSurf2DP1_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.224744871391589*_in[2]+0.7071067811865475*_in[0]; 
   _out[1] = 1.224744871391589*_in[3]+0.7071067811865475*_in[1]; 
  } 
} 
void ModalSerendipBasisSurf2DP1_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 0.7071067811865475*_in[0]-1.224744871391589*_in[2]; 
   _out[1] = 0.7071067811865475*_in[1]-1.224744871391589*_in[3]; 
  } 
} 


// polyOrder 2 
//    dir 1 
void ModalSerendipBasisSurf2DP2_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.58113883008419*_in[4]+1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
   _out[1] = 1.58113883008419*_in[6]+1.224744871391589*_in[3]+0.7071067811865475*_in[2]; 
   _out[2] = 1.224744871391589*_in[7]+0.7071067811865475*_in[5]; 
  } 
} 
void ModalSerendipBasisSurf2DP2_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.58113883008419*_in[4]-1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
   _out[1] = 1.58113883008419*_in[6]-1.224744871391589*_in[3]+0.7071067811865475*_in[2]; 
   _out[2] = 0.7071067811865475*_in[5]-1.224744871391589*_in[7]; 
  } 
} 

//    dir 2 
void ModalSerendipBasisSurf2DP2_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.58113883008419*_in[5]+1.224744871391589*_in[2]+0.7071067811865475*_in[0]; 
   _out[1] = 1.58113883008419*_in[7]+1.224744871391589*_in[3]+0.7071067811865475*_in[1]; 
   _out[2] = 1.224744871391589*_in[6]+0.7071067811865475*_in[4]; 
  } 
} 
void ModalSerendipBasisSurf2DP2_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.58113883008419*_in[5]-1.224744871391589*_in[2]+0.7071067811865475*_in[0]; 
   _out[1] = 1.58113883008419*_in[7]-1.224744871391589*_in[3]+0.7071067811865475*_in[1]; 
   _out[2] = 0.7071067811865475*_in[4]-1.224744871391589*_in[6]; 
  } 
} 


// polyOrder 3 
//    dir 1 
void ModalSerendipBasisSurf2DP3_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.870828693386971*_in[8]+1.58113883008419*_in[4]+1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
   _out[1] = 1.870828693386971*_in[10]+1.58113883008419*_in[6]+1.224744871391589*_in[3]+0.7071067811865475*_in[2]; 
   _out[2] = 1.224744871391589*_in[7]+0.7071067811865475*_in[5]; 
   _out[3] = 1.224744871391589*_in[11]+0.7071067811865475*_in[9]; 
  } 
} 
void ModalSerendipBasisSurf2DP3_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = (-1.870828693386971*_in[8])+1.58113883008419*_in[4]-1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
   _out[1] = (-1.870828693386971*_in[10])+1.58113883008419*_in[6]-1.224744871391589*_in[3]+0.7071067811865475*_in[2]; 
   _out[2] = 0.7071067811865475*_in[5]-1.224744871391589*_in[7]; 
   _out[3] = 0.7071067811865475*_in[9]-1.224744871391589*_in[11]; 
  } 
} 

//    dir 2 
void ModalSerendipBasisSurf2DP3_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 1.870828693386971*_in[9]+1.58113883008419*_in[5]+1.224744871391589*_in[2]+0.7071067811865475*_in[0]; 
   _out[1] = 1.870828693386971*_in[11]+1.58113883008419*_in[7]+1.224744871391589*_in[3]+0.7071067811865475*_in[1]; 
   _out[2] = 1.224744871391589*_in[6]+0.7071067811865475*_in[4]; 
   _out[3] = 1.224744871391589*_in[10]+0.7071067811865475*_in[8]; 
  } 
} 
void ModalSerendipBasisSurf2DP3_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = (-1.870828693386971*_in[9])+1.58113883008419*_in[5]-1.224744871391589*_in[2]+0.7071067811865475*_in[0]; 
   _out[1] = (-1.870828693386971*_in[11])+1.58113883008419*_in[7]-1.224744871391589*_in[3]+0.7071067811865475*_in[1]; 
   _out[2] = 0.7071067811865475*_in[4]-1.224744871391589*_in[6]; 
   _out[3] = 0.7071067811865475*_in[8]-1.224744871391589*_in[10]; 
  } 
} 


// polyOrder 4 
//    dir 1 
void ModalSerendipBasisSurf2DP4_Upper1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 2.121320343559642*_in[13]+1.870828693386971*_in[8]+1.58113883008419*_in[4]+1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
   _out[1] = 2.121320343559642*_in[15]+1.870828693386971*_in[11]+1.58113883008419*_in[6]+1.224744871391589*_in[3]+0.7071067811865475*_in[2]; 
   _out[2] = 1.58113883008419*_in[10]+1.224744871391589*_in[7]+0.7071067811865475*_in[5]; 
   _out[3] = 1.224744871391589*_in[12]+0.7071067811865475*_in[9]; 
   _out[4] = 1.224744871391589*_in[16]+0.7071067811865475*_in[14]; 
  } 
} 
void ModalSerendipBasisSurf2DP4_Lower1(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 2.121320343559642*_in[13]-1.870828693386971*_in[8]+1.58113883008419*_in[4]-1.224744871391589*_in[1]+0.7071067811865475*_in[0]; 
   _out[1] = 2.121320343559642*_in[15]-1.870828693386971*_in[11]+1.58113883008419*_in[6]-1.224744871391589*_in[3]+0.7071067811865475*_in[2]; 
   _out[2] = 1.58113883008419*_in[10]-1.224744871391589*_in[7]+0.7071067811865475*_in[5]; 
   _out[3] = 0.7071067811865475*_in[9]-1.224744871391589*_in[12]; 
   _out[4] = 0.7071067811865475*_in[14]-1.224744871391589*_in[16]; 
  } 
} 

//    dir 2 
void ModalSerendipBasisSurf2DP4_Upper2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 2.121320343559642*_in[14]+1.870828693386971*_in[9]+1.58113883008419*_in[5]+1.224744871391589*_in[2]+0.7071067811865475*_in[0]; 
   _out[1] = 2.121320343559642*_in[16]+1.870828693386971*_in[12]+1.58113883008419*_in[7]+1.224744871391589*_in[3]+0.7071067811865475*_in[1]; 
   _out[2] = 1.58113883008419*_in[10]+1.224744871391589*_in[6]+0.7071067811865475*_in[4]; 
   _out[3] = 1.224744871391589*_in[11]+0.7071067811865475*_in[8]; 
   _out[4] = 1.224744871391589*_in[15]+0.7071067811865475*_in[13]; 
  } 
} 
void ModalSerendipBasisSurf2DP4_Lower2(int meqn, int mbasis, int msurf, const double *volIn, double *surfOut) 
{ 
  for (unsigned m=0; m<meqn; ++m) 
  { 
  double *_out = &surfOut[msurf*m]; 
  const double *_in = &volIn[mbasis*m]; 
   _out[0] = 2.121320343559642*_in[14]-1.870828693386971*_in[9]+1.58113883008419*_in[5]-1.224744871391589*_in[2]+0.7071067811865475*_in[0]; 
   _out[1] = 2.121320343559642*_in[16]-1.870828693386971*_in[12]+1.58113883008419*_in[7]-1.224744871391589*_in[3]+0.7071067811865475*_in[1]; 
   _out[2] = 1.58113883008419*_in[10]-1.224744871391589*_in[6]+0.7071067811865475*_in[4]; 
   _out[3] = 0.7071067811865475*_in[8]-1.224744871391589*_in[11]; 
   _out[4] = 0.7071067811865475*_in[13]-1.224744871391589*_in[15]; 
  } 
} 


