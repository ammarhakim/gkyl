#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[17] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[23] += 22.9128784747792*rdxFnu[0]*f[4]; 
  out[25] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[29] += 22.9128784747792*rdxFnu[0]*f[10]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVol3xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[17] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[18] += 22.9128784747792*rdxFnu[1]*f[2]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[21] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[23] += 22.9128784747792*rdxFnu[0]*f[4]; 
  out[24] += 22.9128784747792*rdxFnu[1]*f[4]; 
  out[25] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[26] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[29] += 22.9128784747792*rdxFnu[0]*f[10]; 
  out[30] += 22.9128784747792*rdxFnu[1]*f[10]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol3xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[17] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[18] += 22.9128784747792*rdxFnu[1]*f[2]; 
  out[19] += 22.9128784747792*rdxFnu[2]*f[3]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[21] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[22] += 6.708203932499369*rdxFnu[2]*f[4]; 
  out[23] += 22.9128784747792*rdxFnu[0]*f[4]; 
  out[24] += 22.9128784747792*rdxFnu[1]*f[4]; 
  out[25] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[26] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[27] += 22.9128784747792*rdxFnu[2]*f[5]; 
  out[28] += 22.9128784747792*rdxFnu[2]*f[6]; 
  out[29] += 22.9128784747792*rdxFnu[0]*f[10]; 
  out[30] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[31] += 22.9128784747792*rdxFnu[2]*f[10]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstDiffusionVol3xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[16] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[17] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[19] += 22.9128784747792*rdxFnu[1]*f[3]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[23] += 22.9128784747792*rdxFnu[0]*f[4]; 
  out[25] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[27] += 22.9128784747792*rdxFnu[1]*f[5]; 
  out[28] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[29] += 22.9128784747792*rdxFnu[0]*f[10]; 
  out[31] += 22.9128784747792*rdxFnu[1]*f[10]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol3xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[8] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[12] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[18] += 22.9128784747792*rdxFnu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[24] += 22.9128784747792*rdxFnu[0]*f[4]; 
  out[26] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[30] += 22.9128784747792*rdxFnu[0]*f[10]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVol3xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[8] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[12] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[16] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[18] += 22.9128784747792*rdxFnu[0]*f[2]; 
  out[19] += 22.9128784747792*rdxFnu[1]*f[3]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[24] += 22.9128784747792*rdxFnu[0]*f[4]; 
  out[26] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[27] += 22.9128784747792*rdxFnu[1]*f[5]; 
  out[28] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[30] += 22.9128784747792*rdxFnu[0]*f[10]; 
  out[31] += 22.9128784747792*rdxFnu[1]*f[10]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol3xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[9] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[15] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[16] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[19] += 22.9128784747792*rdxFnu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[27] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[28] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[31] += 22.9128784747792*rdxFnu[0]*f[10]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol3xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol3xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol3xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol3xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol3xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol3xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol3xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol3xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol3xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol3xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol3xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol3xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol3xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol3xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol3xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[96]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 3.622844186547359*rdxF[0]*f[7]*nu[17]+1.620185174601965*f[0]*rdxF[0]*nu[17]+2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[4] += 3.622844186547359*rdxF[0]*f[11]*nu[17]+1.620185174601965*rdxF[0]*f[2]*nu[17]+2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 3.622844186547359*rdxF[0]*f[13]*nu[17]+1.620185174601965*rdxF[0]*f[3]*nu[17]+2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[7] += 9.48683298050514*rdxF[0]*f[17]*nu[17]+10.86853255964208*rdxF[0]*f[1]*nu[17]+7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[10] += 3.622844186547359*rdxF[0]*nu[17]*f[20]+1.620185174601965*rdxF[0]*f[6]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 9.486832980505138*rdxF[0]*nu[17]*f[23]+10.86853255964208*rdxF[0]*f[4]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 1.620185174601965*rdxF[0]*f[8]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 9.486832980505138*rdxF[0]*nu[17]*f[25]+10.86853255964208*rdxF[0]*f[5]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[15] += 1.620185174601965*rdxF[0]*f[9]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[17] += 26.08879069638912*rdxF[0]*f[7]*nu[17]+14.8492424049175*f[0]*rdxF[0]*nu[17]+14.23024947075771*rdxF[0]*nu[7]*f[17]+18.1142209327368*rdxF[0]*f[1]*nu[7]+10.86853255964208*rdxF[0]*nu[1]*f[7]+9.721111047611789*f[0]*rdxF[0]*nu[1]+8.100925873009823*nu[0]*rdxF[0]*f[1]; 
  out[20] += 9.48683298050514*rdxF[0]*nu[17]*f[29]+7.115124735378852*rdxF[0]*nu[7]*f[20]+10.86853255964208*rdxF[0]*f[10]*nu[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[21] += 2.371708245126284*rdxF[0]*nu[7]*f[21]+1.620185174601965*rdxF[0]*f[14]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[22] += 2.371708245126284*rdxF[0]*nu[7]*f[22]+1.620185174601965*rdxF[0]*f[16]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 
  out[23] += 14.23024947075771*rdxF[0]*nu[7]*f[23]+26.08879069638912*rdxF[0]*f[11]*nu[17]+14.8492424049175*rdxF[0]*f[2]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[11]+18.1142209327368*rdxF[0]*f[4]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[4]+9.721111047611787*rdxF[0]*nu[1]*f[2]; 
  out[24] += 2.371708245126284*rdxF[0]*nu[7]*f[24]+1.620185174601965*rdxF[0]*nu[17]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[18]; 
  out[25] += 14.23024947075771*rdxF[0]*nu[7]*f[25]+26.08879069638912*rdxF[0]*f[13]*nu[17]+14.8492424049175*rdxF[0]*f[3]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[13]+18.1142209327368*rdxF[0]*f[5]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[5]+9.721111047611787*rdxF[0]*nu[1]*f[3]; 
  out[27] += 2.371708245126284*rdxF[0]*nu[7]*f[27]+1.620185174601965*rdxF[0]*nu[17]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[19]; 
  out[29] += 14.23024947075771*rdxF[0]*nu[7]*f[29]+26.08879069638912*rdxF[0]*nu[17]*f[20]+10.86853255964208*rdxF[0]*nu[1]*f[20]+14.8492424049175*rdxF[0]*f[6]*nu[17]+18.1142209327368*rdxF[0]*nu[7]*f[10]+8.100925873009823*nu[0]*rdxF[0]*f[10]+9.721111047611789*rdxF[0]*nu[1]*f[6]; 
  out[30] += 2.371708245126284*rdxF[0]*nu[7]*f[30]+1.620185174601965*rdxF[0]*nu[17]*f[26]+1.060660171779821*rdxF[0]*nu[1]*f[26]; 
  out[31] += 2.371708245126284*rdxF[0]*nu[7]*f[31]+1.620185174601965*rdxF[0]*nu[17]*f[28]+1.060660171779821*rdxF[0]*nu[1]*f[28]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol3xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[96]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 3.622844186547359*rdxF[0]*f[7]*nu[17]+1.620185174601965*f[0]*rdxF[0]*nu[17]+2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[2] += 3.622844186547359*rdxF[1]*f[8]*nu[50]+1.620185174601965*f[0]*rdxF[1]*nu[50]+2.371708245126284*rdxF[1]*f[2]*nu[40]+1.060660171779821*f[0]*rdxF[1]*nu[34]; 
  out[4] += 3.622844186547359*rdxF[1]*f[12]*nu[50]+1.620185174601965*f[1]*rdxF[1]*nu[50]+2.371708245126284*rdxF[1]*f[4]*nu[40]+1.060660171779821*f[1]*rdxF[1]*nu[34]+3.622844186547359*rdxF[0]*f[11]*nu[17]+1.620185174601965*rdxF[0]*f[2]*nu[17]+2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 3.622844186547359*rdxF[0]*f[13]*nu[17]+1.620185174601965*rdxF[0]*f[3]*nu[17]+2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 3.622844186547359*rdxF[1]*f[14]*nu[50]+1.620185174601965*rdxF[1]*f[3]*nu[50]+2.371708245126284*rdxF[1]*f[6]*nu[40]+1.060660171779821*rdxF[1]*f[3]*nu[34]; 
  out[7] += 9.48683298050514*rdxF[0]*f[17]*nu[17]+10.86853255964208*rdxF[0]*f[1]*nu[17]+7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[8] += 9.48683298050514*rdxF[1]*f[18]*nu[50]+10.86853255964208*rdxF[1]*f[2]*nu[50]+7.115124735378852*rdxF[1]*f[8]*nu[40]+5.303300858899105*f[0]*rdxF[1]*nu[40]+4.743416490252569*rdxF[1]*f[2]*nu[34]+2.371708245126284*f[0]*rdxF[1]*nu[32]; 
  out[10] += 3.622844186547359*rdxF[1]*f[21]*nu[50]+1.620185174601965*rdxF[1]*f[5]*nu[50]+2.371708245126284*rdxF[1]*f[10]*nu[40]+1.060660171779821*rdxF[1]*f[5]*nu[34]+3.622844186547359*rdxF[0]*nu[17]*f[20]+1.620185174601965*rdxF[0]*f[6]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 1.620185174601965*rdxF[1]*f[7]*nu[50]+2.371708245126284*rdxF[1]*f[11]*nu[40]+1.060660171779821*rdxF[1]*f[7]*nu[34]+9.486832980505138*rdxF[0]*nu[17]*f[23]+10.86853255964208*rdxF[0]*f[4]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 9.486832980505138*rdxF[1]*f[24]*nu[50]+10.86853255964208*rdxF[1]*f[4]*nu[50]+7.115124735378852*rdxF[1]*f[12]*nu[40]+5.303300858899106*f[1]*rdxF[1]*nu[40]+4.743416490252569*rdxF[1]*f[4]*nu[34]+2.371708245126284*f[1]*rdxF[1]*nu[32]+1.620185174601965*rdxF[0]*f[8]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 9.486832980505138*rdxF[0]*nu[17]*f[25]+10.86853255964208*rdxF[0]*f[5]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[14] += 9.486832980505138*rdxF[1]*f[26]*nu[50]+10.86853255964208*rdxF[1]*f[6]*nu[50]+7.115124735378852*rdxF[1]*f[14]*nu[40]+5.303300858899106*rdxF[1]*f[3]*nu[40]+4.743416490252569*rdxF[1]*f[6]*nu[34]+2.371708245126284*rdxF[1]*f[3]*nu[32]; 
  out[15] += 1.620185174601965*rdxF[0]*f[9]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[16] += 1.620185174601965*rdxF[1]*f[9]*nu[50]+2.371708245126284*rdxF[1]*f[16]*nu[40]+1.060660171779821*rdxF[1]*f[9]*nu[34]; 
  out[17] += 26.08879069638912*rdxF[0]*f[7]*nu[17]+14.8492424049175*f[0]*rdxF[0]*nu[17]+14.23024947075771*rdxF[0]*nu[7]*f[17]+18.1142209327368*rdxF[0]*f[1]*nu[7]+10.86853255964208*rdxF[0]*nu[1]*f[7]+9.721111047611789*f[0]*rdxF[0]*nu[1]+8.100925873009823*nu[0]*rdxF[0]*f[1]; 
  out[18] += 26.08879069638912*rdxF[1]*f[8]*nu[50]+14.8492424049175*f[0]*rdxF[1]*nu[50]+14.23024947075771*rdxF[1]*f[18]*nu[40]+18.1142209327368*rdxF[1]*f[2]*nu[40]+10.86853255964208*rdxF[1]*f[8]*nu[34]+9.721111047611789*f[0]*rdxF[1]*nu[34]+8.100925873009823*rdxF[1]*f[2]*nu[32]; 
  out[20] += 1.620185174601965*rdxF[1]*f[13]*nu[50]+2.371708245126284*rdxF[1]*f[20]*nu[40]+1.060660171779821*rdxF[1]*f[13]*nu[34]+9.48683298050514*rdxF[0]*nu[17]*f[29]+7.115124735378852*rdxF[0]*nu[7]*f[20]+10.86853255964208*rdxF[0]*f[10]*nu[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[21] += 9.48683298050514*rdxF[1]*f[30]*nu[50]+10.86853255964208*rdxF[1]*f[10]*nu[50]+7.115124735378852*rdxF[1]*f[21]*nu[40]+5.303300858899105*rdxF[1]*f[5]*nu[40]+4.743416490252569*rdxF[1]*f[10]*nu[34]+2.371708245126284*rdxF[1]*f[5]*nu[32]+2.371708245126284*rdxF[0]*nu[7]*f[21]+1.620185174601965*rdxF[0]*f[14]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[22] += 1.620185174601965*rdxF[1]*f[15]*nu[50]+2.371708245126284*rdxF[1]*f[22]*nu[40]+1.060660171779821*rdxF[1]*f[15]*nu[34]+2.371708245126284*rdxF[0]*nu[7]*f[22]+1.620185174601965*rdxF[0]*f[16]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 
  out[23] += 1.620185174601965*rdxF[1]*f[17]*nu[50]+2.371708245126284*rdxF[1]*f[23]*nu[40]+1.060660171779821*rdxF[1]*f[17]*nu[34]+14.23024947075771*rdxF[0]*nu[7]*f[23]+26.08879069638912*rdxF[0]*f[11]*nu[17]+14.8492424049175*rdxF[0]*f[2]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[11]+18.1142209327368*rdxF[0]*f[4]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[4]+9.721111047611787*rdxF[0]*nu[1]*f[2]; 
  out[24] += 26.08879069638912*rdxF[1]*f[12]*nu[50]+14.8492424049175*f[1]*rdxF[1]*nu[50]+14.23024947075771*rdxF[1]*f[24]*nu[40]+18.1142209327368*rdxF[1]*f[4]*nu[40]+10.86853255964208*rdxF[1]*f[12]*nu[34]+9.721111047611787*f[1]*rdxF[1]*nu[34]+8.100925873009823*rdxF[1]*f[4]*nu[32]+2.371708245126284*rdxF[0]*nu[7]*f[24]+1.620185174601965*rdxF[0]*nu[17]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[18]; 
  out[25] += 14.23024947075771*rdxF[0]*nu[7]*f[25]+26.08879069638912*rdxF[0]*f[13]*nu[17]+14.8492424049175*rdxF[0]*f[3]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[13]+18.1142209327368*rdxF[0]*f[5]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[5]+9.721111047611787*rdxF[0]*nu[1]*f[3]; 
  out[26] += 26.08879069638912*rdxF[1]*f[14]*nu[50]+14.8492424049175*rdxF[1]*f[3]*nu[50]+14.23024947075771*rdxF[1]*f[26]*nu[40]+18.1142209327368*rdxF[1]*f[6]*nu[40]+10.86853255964208*rdxF[1]*f[14]*nu[34]+9.721111047611787*rdxF[1]*f[3]*nu[34]+8.100925873009823*rdxF[1]*f[6]*nu[32]; 
  out[27] += 2.371708245126284*rdxF[0]*nu[7]*f[27]+1.620185174601965*rdxF[0]*nu[17]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[19]; 
  out[28] += 1.620185174601965*rdxF[1]*f[19]*nu[50]+2.371708245126284*rdxF[1]*f[28]*nu[40]+1.060660171779821*rdxF[1]*f[19]*nu[34]; 
  out[29] += 1.620185174601965*rdxF[1]*f[25]*nu[50]+2.371708245126284*rdxF[1]*f[29]*nu[40]+1.060660171779821*rdxF[1]*f[25]*nu[34]+14.23024947075771*rdxF[0]*nu[7]*f[29]+26.08879069638912*rdxF[0]*nu[17]*f[20]+10.86853255964208*rdxF[0]*nu[1]*f[20]+14.8492424049175*rdxF[0]*f[6]*nu[17]+18.1142209327368*rdxF[0]*nu[7]*f[10]+8.100925873009823*nu[0]*rdxF[0]*f[10]+9.721111047611789*rdxF[0]*nu[1]*f[6]; 
  out[30] += 26.08879069638912*rdxF[1]*f[21]*nu[50]+14.8492424049175*rdxF[1]*f[5]*nu[50]+14.23024947075771*rdxF[1]*f[30]*nu[40]+18.1142209327368*rdxF[1]*f[10]*nu[40]+10.86853255964208*rdxF[1]*f[21]*nu[34]+9.721111047611789*rdxF[1]*f[5]*nu[34]+8.100925873009823*rdxF[1]*f[10]*nu[32]+2.371708245126284*rdxF[0]*nu[7]*f[30]+1.620185174601965*rdxF[0]*nu[17]*f[26]+1.060660171779821*rdxF[0]*nu[1]*f[26]; 
  out[31] += 1.620185174601965*rdxF[1]*f[27]*nu[50]+2.371708245126284*rdxF[1]*f[31]*nu[40]+1.060660171779821*rdxF[1]*f[27]*nu[34]+2.371708245126284*rdxF[0]*nu[7]*f[31]+1.620185174601965*rdxF[0]*nu[17]*f[28]+1.060660171779821*rdxF[0]*nu[1]*f[28]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+rdxF[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol3xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[96]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 

  out[1] += 3.622844186547359*rdxF[0]*f[7]*nu[17]+1.620185174601965*f[0]*rdxF[0]*nu[17]+2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[2] += 3.622844186547359*rdxF[1]*f[8]*nu[50]+1.620185174601965*f[0]*rdxF[1]*nu[50]+2.371708245126284*rdxF[1]*f[2]*nu[40]+1.060660171779821*f[0]*rdxF[1]*nu[34]; 
  out[3] += 3.622844186547359*rdxF[2]*f[9]*nu[83]+1.620185174601965*f[0]*rdxF[2]*nu[83]+2.371708245126284*rdxF[2]*f[3]*nu[73]+1.060660171779821*f[0]*rdxF[2]*nu[67]; 
  out[4] += 3.622844186547359*rdxF[1]*f[12]*nu[50]+1.620185174601965*f[1]*rdxF[1]*nu[50]+2.371708245126284*rdxF[1]*f[4]*nu[40]+1.060660171779821*f[1]*rdxF[1]*nu[34]+3.622844186547359*rdxF[0]*f[11]*nu[17]+1.620185174601965*rdxF[0]*f[2]*nu[17]+2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 3.622844186547359*rdxF[2]*f[15]*nu[83]+1.620185174601965*f[1]*rdxF[2]*nu[83]+2.371708245126284*rdxF[2]*f[5]*nu[73]+1.060660171779821*f[1]*rdxF[2]*nu[67]+3.622844186547359*rdxF[0]*f[13]*nu[17]+1.620185174601965*rdxF[0]*f[3]*nu[17]+2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 3.622844186547359*rdxF[2]*f[16]*nu[83]+1.620185174601965*f[2]*rdxF[2]*nu[83]+2.371708245126284*rdxF[2]*f[6]*nu[73]+1.060660171779821*f[2]*rdxF[2]*nu[67]+3.622844186547359*rdxF[1]*f[14]*nu[50]+1.620185174601965*rdxF[1]*f[3]*nu[50]+2.371708245126284*rdxF[1]*f[6]*nu[40]+1.060660171779821*rdxF[1]*f[3]*nu[34]; 
  out[7] += 9.48683298050514*rdxF[0]*f[17]*nu[17]+10.86853255964208*rdxF[0]*f[1]*nu[17]+7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[8] += 9.48683298050514*rdxF[1]*f[18]*nu[50]+10.86853255964208*rdxF[1]*f[2]*nu[50]+7.115124735378852*rdxF[1]*f[8]*nu[40]+5.303300858899105*f[0]*rdxF[1]*nu[40]+4.743416490252569*rdxF[1]*f[2]*nu[34]+2.371708245126284*f[0]*rdxF[1]*nu[32]; 
  out[9] += 9.48683298050514*rdxF[2]*f[19]*nu[83]+10.86853255964208*rdxF[2]*f[3]*nu[83]+7.115124735378852*rdxF[2]*f[9]*nu[73]+5.303300858899105*f[0]*rdxF[2]*nu[73]+4.743416490252569*rdxF[2]*f[3]*nu[67]+2.371708245126284*f[0]*rdxF[2]*nu[64]; 
  out[10] += 3.622844186547359*rdxF[2]*f[22]*nu[83]+1.620185174601965*rdxF[2]*f[4]*nu[83]+2.371708245126284*rdxF[2]*f[10]*nu[73]+1.060660171779821*rdxF[2]*f[4]*nu[67]+3.622844186547359*rdxF[1]*f[21]*nu[50]+1.620185174601965*rdxF[1]*f[5]*nu[50]+2.371708245126284*rdxF[1]*f[10]*nu[40]+1.060660171779821*rdxF[1]*f[5]*nu[34]+3.622844186547359*rdxF[0]*nu[17]*f[20]+1.620185174601965*rdxF[0]*f[6]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 1.620185174601965*rdxF[1]*f[7]*nu[50]+2.371708245126284*rdxF[1]*f[11]*nu[40]+1.060660171779821*rdxF[1]*f[7]*nu[34]+9.486832980505138*rdxF[0]*nu[17]*f[23]+10.86853255964208*rdxF[0]*f[4]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 9.486832980505138*rdxF[1]*f[24]*nu[50]+10.86853255964208*rdxF[1]*f[4]*nu[50]+7.115124735378852*rdxF[1]*f[12]*nu[40]+5.303300858899106*f[1]*rdxF[1]*nu[40]+4.743416490252569*rdxF[1]*f[4]*nu[34]+2.371708245126284*f[1]*rdxF[1]*nu[32]+1.620185174601965*rdxF[0]*f[8]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 1.620185174601965*rdxF[2]*f[7]*nu[83]+2.371708245126284*rdxF[2]*f[13]*nu[73]+1.060660171779821*rdxF[2]*f[7]*nu[67]+9.486832980505138*rdxF[0]*nu[17]*f[25]+10.86853255964208*rdxF[0]*f[5]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[14] += 1.620185174601965*rdxF[2]*f[8]*nu[83]+2.371708245126284*rdxF[2]*f[14]*nu[73]+1.060660171779821*rdxF[2]*f[8]*nu[67]+9.486832980505138*rdxF[1]*f[26]*nu[50]+10.86853255964208*rdxF[1]*f[6]*nu[50]+7.115124735378852*rdxF[1]*f[14]*nu[40]+5.303300858899106*rdxF[1]*f[3]*nu[40]+4.743416490252569*rdxF[1]*f[6]*nu[34]+2.371708245126284*rdxF[1]*f[3]*nu[32]; 
  out[15] += 9.486832980505138*rdxF[2]*f[27]*nu[83]+10.86853255964208*rdxF[2]*f[5]*nu[83]+7.115124735378852*rdxF[2]*f[15]*nu[73]+5.303300858899106*f[1]*rdxF[2]*nu[73]+4.743416490252569*rdxF[2]*f[5]*nu[67]+2.371708245126284*f[1]*rdxF[2]*nu[64]+1.620185174601965*rdxF[0]*f[9]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[16] += 9.486832980505138*rdxF[2]*f[28]*nu[83]+10.86853255964208*rdxF[2]*f[6]*nu[83]+7.115124735378852*rdxF[2]*f[16]*nu[73]+5.303300858899106*f[2]*rdxF[2]*nu[73]+4.743416490252569*rdxF[2]*f[6]*nu[67]+2.371708245126284*f[2]*rdxF[2]*nu[64]+1.620185174601965*rdxF[1]*f[9]*nu[50]+2.371708245126284*rdxF[1]*f[16]*nu[40]+1.060660171779821*rdxF[1]*f[9]*nu[34]; 
  out[17] += 26.08879069638912*rdxF[0]*f[7]*nu[17]+14.8492424049175*f[0]*rdxF[0]*nu[17]+14.23024947075771*rdxF[0]*nu[7]*f[17]+18.1142209327368*rdxF[0]*f[1]*nu[7]+10.86853255964208*rdxF[0]*nu[1]*f[7]+9.721111047611789*f[0]*rdxF[0]*nu[1]+8.100925873009823*nu[0]*rdxF[0]*f[1]; 
  out[18] += 26.08879069638912*rdxF[1]*f[8]*nu[50]+14.8492424049175*f[0]*rdxF[1]*nu[50]+14.23024947075771*rdxF[1]*f[18]*nu[40]+18.1142209327368*rdxF[1]*f[2]*nu[40]+10.86853255964208*rdxF[1]*f[8]*nu[34]+9.721111047611789*f[0]*rdxF[1]*nu[34]+8.100925873009823*rdxF[1]*f[2]*nu[32]; 
  out[19] += 26.08879069638912*rdxF[2]*f[9]*nu[83]+14.8492424049175*f[0]*rdxF[2]*nu[83]+14.23024947075771*rdxF[2]*f[19]*nu[73]+18.1142209327368*rdxF[2]*f[3]*nu[73]+10.86853255964208*rdxF[2]*f[9]*nu[67]+9.721111047611789*f[0]*rdxF[2]*nu[67]+8.100925873009823*rdxF[2]*f[3]*nu[64]; 
  out[20] += 1.620185174601965*rdxF[2]*f[11]*nu[83]+2.371708245126284*rdxF[2]*f[20]*nu[73]+1.060660171779821*rdxF[2]*f[11]*nu[67]+1.620185174601965*rdxF[1]*f[13]*nu[50]+2.371708245126284*rdxF[1]*f[20]*nu[40]+1.060660171779821*rdxF[1]*f[13]*nu[34]+9.48683298050514*rdxF[0]*nu[17]*f[29]+7.115124735378852*rdxF[0]*nu[7]*f[20]+10.86853255964208*rdxF[0]*f[10]*nu[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[21] += 1.620185174601965*rdxF[2]*f[12]*nu[83]+2.371708245126284*rdxF[2]*f[21]*nu[73]+1.060660171779821*rdxF[2]*f[12]*nu[67]+9.48683298050514*rdxF[1]*f[30]*nu[50]+10.86853255964208*rdxF[1]*f[10]*nu[50]+7.115124735378852*rdxF[1]*f[21]*nu[40]+5.303300858899105*rdxF[1]*f[5]*nu[40]+4.743416490252569*rdxF[1]*f[10]*nu[34]+2.371708245126284*rdxF[1]*f[5]*nu[32]+2.371708245126284*rdxF[0]*nu[7]*f[21]+1.620185174601965*rdxF[0]*f[14]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[22] += 9.48683298050514*rdxF[2]*f[31]*nu[83]+10.86853255964208*rdxF[2]*f[10]*nu[83]+7.115124735378852*rdxF[2]*f[22]*nu[73]+5.303300858899105*rdxF[2]*f[4]*nu[73]+4.743416490252569*rdxF[2]*f[10]*nu[67]+2.371708245126284*rdxF[2]*f[4]*nu[64]+1.620185174601965*rdxF[1]*f[15]*nu[50]+2.371708245126284*rdxF[1]*f[22]*nu[40]+1.060660171779821*rdxF[1]*f[15]*nu[34]+2.371708245126284*rdxF[0]*nu[7]*f[22]+1.620185174601965*rdxF[0]*f[16]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 
  out[23] += 1.620185174601965*rdxF[1]*f[17]*nu[50]+2.371708245126284*rdxF[1]*f[23]*nu[40]+1.060660171779821*rdxF[1]*f[17]*nu[34]+14.23024947075771*rdxF[0]*nu[7]*f[23]+26.08879069638912*rdxF[0]*f[11]*nu[17]+14.8492424049175*rdxF[0]*f[2]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[11]+18.1142209327368*rdxF[0]*f[4]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[4]+9.721111047611787*rdxF[0]*nu[1]*f[2]; 
  out[24] += 26.08879069638912*rdxF[1]*f[12]*nu[50]+14.8492424049175*f[1]*rdxF[1]*nu[50]+14.23024947075771*rdxF[1]*f[24]*nu[40]+18.1142209327368*rdxF[1]*f[4]*nu[40]+10.86853255964208*rdxF[1]*f[12]*nu[34]+9.721111047611787*f[1]*rdxF[1]*nu[34]+8.100925873009823*rdxF[1]*f[4]*nu[32]+2.371708245126284*rdxF[0]*nu[7]*f[24]+1.620185174601965*rdxF[0]*nu[17]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[18]; 
  out[25] += 1.620185174601965*rdxF[2]*f[17]*nu[83]+2.371708245126284*rdxF[2]*f[25]*nu[73]+1.060660171779821*rdxF[2]*f[17]*nu[67]+14.23024947075771*rdxF[0]*nu[7]*f[25]+26.08879069638912*rdxF[0]*f[13]*nu[17]+14.8492424049175*rdxF[0]*f[3]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[13]+18.1142209327368*rdxF[0]*f[5]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[5]+9.721111047611787*rdxF[0]*nu[1]*f[3]; 
  out[26] += 1.620185174601965*rdxF[2]*f[18]*nu[83]+2.371708245126284*rdxF[2]*f[26]*nu[73]+1.060660171779821*rdxF[2]*f[18]*nu[67]+26.08879069638912*rdxF[1]*f[14]*nu[50]+14.8492424049175*rdxF[1]*f[3]*nu[50]+14.23024947075771*rdxF[1]*f[26]*nu[40]+18.1142209327368*rdxF[1]*f[6]*nu[40]+10.86853255964208*rdxF[1]*f[14]*nu[34]+9.721111047611787*rdxF[1]*f[3]*nu[34]+8.100925873009823*rdxF[1]*f[6]*nu[32]; 
  out[27] += 26.08879069638912*rdxF[2]*f[15]*nu[83]+14.8492424049175*f[1]*rdxF[2]*nu[83]+14.23024947075771*rdxF[2]*f[27]*nu[73]+18.1142209327368*rdxF[2]*f[5]*nu[73]+10.86853255964208*rdxF[2]*f[15]*nu[67]+9.721111047611787*f[1]*rdxF[2]*nu[67]+8.100925873009823*rdxF[2]*f[5]*nu[64]+2.371708245126284*rdxF[0]*nu[7]*f[27]+1.620185174601965*rdxF[0]*nu[17]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[19]; 
  out[28] += 26.08879069638912*rdxF[2]*f[16]*nu[83]+14.8492424049175*f[2]*rdxF[2]*nu[83]+14.23024947075771*rdxF[2]*f[28]*nu[73]+18.1142209327368*rdxF[2]*f[6]*nu[73]+10.86853255964208*rdxF[2]*f[16]*nu[67]+9.721111047611787*f[2]*rdxF[2]*nu[67]+8.100925873009823*rdxF[2]*f[6]*nu[64]+1.620185174601965*rdxF[1]*f[19]*nu[50]+2.371708245126284*rdxF[1]*f[28]*nu[40]+1.060660171779821*rdxF[1]*f[19]*nu[34]; 
  out[29] += 1.620185174601965*rdxF[2]*f[23]*nu[83]+2.371708245126284*rdxF[2]*f[29]*nu[73]+1.060660171779821*rdxF[2]*f[23]*nu[67]+1.620185174601965*rdxF[1]*f[25]*nu[50]+2.371708245126284*rdxF[1]*f[29]*nu[40]+1.060660171779821*rdxF[1]*f[25]*nu[34]+14.23024947075771*rdxF[0]*nu[7]*f[29]+26.08879069638912*rdxF[0]*nu[17]*f[20]+10.86853255964208*rdxF[0]*nu[1]*f[20]+14.8492424049175*rdxF[0]*f[6]*nu[17]+18.1142209327368*rdxF[0]*nu[7]*f[10]+8.100925873009823*nu[0]*rdxF[0]*f[10]+9.721111047611789*rdxF[0]*nu[1]*f[6]; 
  out[30] += 1.620185174601965*rdxF[2]*f[24]*nu[83]+2.371708245126284*rdxF[2]*f[30]*nu[73]+1.060660171779821*rdxF[2]*f[24]*nu[67]+26.08879069638912*rdxF[1]*f[21]*nu[50]+14.8492424049175*rdxF[1]*f[5]*nu[50]+14.23024947075771*rdxF[1]*f[30]*nu[40]+18.1142209327368*rdxF[1]*f[10]*nu[40]+10.86853255964208*rdxF[1]*f[21]*nu[34]+9.721111047611789*rdxF[1]*f[5]*nu[34]+8.100925873009823*rdxF[1]*f[10]*nu[32]+2.371708245126284*rdxF[0]*nu[7]*f[30]+1.620185174601965*rdxF[0]*nu[17]*f[26]+1.060660171779821*rdxF[0]*nu[1]*f[26]; 
  out[31] += 26.08879069638912*rdxF[2]*f[22]*nu[83]+14.8492424049175*rdxF[2]*f[4]*nu[83]+14.23024947075771*rdxF[2]*f[31]*nu[73]+18.1142209327368*rdxF[2]*f[10]*nu[73]+10.86853255964208*rdxF[2]*f[22]*nu[67]+9.721111047611789*rdxF[2]*f[4]*nu[67]+8.100925873009823*rdxF[2]*f[10]*nu[64]+1.620185174601965*rdxF[1]*f[27]*nu[50]+2.371708245126284*rdxF[1]*f[31]*nu[40]+1.060660171779821*rdxF[1]*f[27]*nu[34]+2.371708245126284*rdxF[0]*nu[7]*f[31]+1.620185174601965*rdxF[0]*nu[17]*f[28]+1.060660171779821*rdxF[0]*nu[1]*f[28]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+rdxF[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+rdxF[2]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol3xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[96]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[1] += 3.622844186547359*rdxF[0]*f[7]*nu[17]+1.620185174601965*f[0]*rdxF[0]*nu[17]+2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[3] += 3.622844186547359*rdxF[1]*f[9]*nu[83]+1.620185174601965*f[0]*rdxF[1]*nu[83]+2.371708245126284*rdxF[1]*f[3]*nu[73]+1.060660171779821*f[0]*rdxF[1]*nu[67]; 
  out[4] += 3.622844186547359*rdxF[0]*f[11]*nu[17]+1.620185174601965*rdxF[0]*f[2]*nu[17]+2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 3.622844186547359*rdxF[1]*f[15]*nu[83]+1.620185174601965*f[1]*rdxF[1]*nu[83]+2.371708245126284*rdxF[1]*f[5]*nu[73]+1.060660171779821*f[1]*rdxF[1]*nu[67]+3.622844186547359*rdxF[0]*f[13]*nu[17]+1.620185174601965*rdxF[0]*f[3]*nu[17]+2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 3.622844186547359*rdxF[1]*f[16]*nu[83]+1.620185174601965*rdxF[1]*f[2]*nu[83]+2.371708245126284*rdxF[1]*f[6]*nu[73]+1.060660171779821*rdxF[1]*f[2]*nu[67]; 
  out[7] += 9.48683298050514*rdxF[0]*f[17]*nu[17]+10.86853255964208*rdxF[0]*f[1]*nu[17]+7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[9] += 9.48683298050514*rdxF[1]*f[19]*nu[83]+10.86853255964208*rdxF[1]*f[3]*nu[83]+7.115124735378852*rdxF[1]*f[9]*nu[73]+5.303300858899105*f[0]*rdxF[1]*nu[73]+4.743416490252569*rdxF[1]*f[3]*nu[67]+2.371708245126284*f[0]*rdxF[1]*nu[64]; 
  out[10] += 3.622844186547359*rdxF[1]*f[22]*nu[83]+1.620185174601965*rdxF[1]*f[4]*nu[83]+2.371708245126284*rdxF[1]*f[10]*nu[73]+1.060660171779821*rdxF[1]*f[4]*nu[67]+3.622844186547359*rdxF[0]*nu[17]*f[20]+1.620185174601965*rdxF[0]*f[6]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 9.486832980505138*rdxF[0]*nu[17]*f[23]+10.86853255964208*rdxF[0]*f[4]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 1.620185174601965*rdxF[0]*f[8]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 1.620185174601965*rdxF[1]*f[7]*nu[83]+2.371708245126284*rdxF[1]*f[13]*nu[73]+1.060660171779821*rdxF[1]*f[7]*nu[67]+9.486832980505138*rdxF[0]*nu[17]*f[25]+10.86853255964208*rdxF[0]*f[5]*nu[17]+7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[14] += 1.620185174601965*rdxF[1]*f[8]*nu[83]+2.371708245126284*rdxF[1]*f[14]*nu[73]+1.060660171779821*rdxF[1]*f[8]*nu[67]; 
  out[15] += 9.486832980505138*rdxF[1]*f[27]*nu[83]+10.86853255964208*rdxF[1]*f[5]*nu[83]+7.115124735378852*rdxF[1]*f[15]*nu[73]+5.303300858899106*f[1]*rdxF[1]*nu[73]+4.743416490252569*rdxF[1]*f[5]*nu[67]+2.371708245126284*f[1]*rdxF[1]*nu[64]+1.620185174601965*rdxF[0]*f[9]*nu[17]+2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[16] += 9.486832980505138*rdxF[1]*f[28]*nu[83]+10.86853255964208*rdxF[1]*f[6]*nu[83]+7.115124735378852*rdxF[1]*f[16]*nu[73]+5.303300858899106*rdxF[1]*f[2]*nu[73]+4.743416490252569*rdxF[1]*f[6]*nu[67]+2.371708245126284*rdxF[1]*f[2]*nu[64]; 
  out[17] += 26.08879069638912*rdxF[0]*f[7]*nu[17]+14.8492424049175*f[0]*rdxF[0]*nu[17]+14.23024947075771*rdxF[0]*nu[7]*f[17]+18.1142209327368*rdxF[0]*f[1]*nu[7]+10.86853255964208*rdxF[0]*nu[1]*f[7]+9.721111047611789*f[0]*rdxF[0]*nu[1]+8.100925873009823*nu[0]*rdxF[0]*f[1]; 
  out[19] += 26.08879069638912*rdxF[1]*f[9]*nu[83]+14.8492424049175*f[0]*rdxF[1]*nu[83]+14.23024947075771*rdxF[1]*f[19]*nu[73]+18.1142209327368*rdxF[1]*f[3]*nu[73]+10.86853255964208*rdxF[1]*f[9]*nu[67]+9.721111047611789*f[0]*rdxF[1]*nu[67]+8.100925873009823*rdxF[1]*f[3]*nu[64]; 
  out[20] += 1.620185174601965*rdxF[1]*f[11]*nu[83]+2.371708245126284*rdxF[1]*f[20]*nu[73]+1.060660171779821*rdxF[1]*f[11]*nu[67]+9.48683298050514*rdxF[0]*nu[17]*f[29]+7.115124735378852*rdxF[0]*nu[7]*f[20]+10.86853255964208*rdxF[0]*f[10]*nu[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[21] += 1.620185174601965*rdxF[1]*f[12]*nu[83]+2.371708245126284*rdxF[1]*f[21]*nu[73]+1.060660171779821*rdxF[1]*f[12]*nu[67]+2.371708245126284*rdxF[0]*nu[7]*f[21]+1.620185174601965*rdxF[0]*f[14]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[22] += 9.48683298050514*rdxF[1]*f[31]*nu[83]+10.86853255964208*rdxF[1]*f[10]*nu[83]+7.115124735378852*rdxF[1]*f[22]*nu[73]+5.303300858899105*rdxF[1]*f[4]*nu[73]+4.743416490252569*rdxF[1]*f[10]*nu[67]+2.371708245126284*rdxF[1]*f[4]*nu[64]+2.371708245126284*rdxF[0]*nu[7]*f[22]+1.620185174601965*rdxF[0]*f[16]*nu[17]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 
  out[23] += 14.23024947075771*rdxF[0]*nu[7]*f[23]+26.08879069638912*rdxF[0]*f[11]*nu[17]+14.8492424049175*rdxF[0]*f[2]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[11]+18.1142209327368*rdxF[0]*f[4]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[4]+9.721111047611787*rdxF[0]*nu[1]*f[2]; 
  out[24] += 2.371708245126284*rdxF[0]*nu[7]*f[24]+1.620185174601965*rdxF[0]*nu[17]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[18]; 
  out[25] += 1.620185174601965*rdxF[1]*f[17]*nu[83]+2.371708245126284*rdxF[1]*f[25]*nu[73]+1.060660171779821*rdxF[1]*f[17]*nu[67]+14.23024947075771*rdxF[0]*nu[7]*f[25]+26.08879069638912*rdxF[0]*f[13]*nu[17]+14.8492424049175*rdxF[0]*f[3]*nu[17]+10.86853255964208*rdxF[0]*nu[1]*f[13]+18.1142209327368*rdxF[0]*f[5]*nu[7]+8.100925873009823*nu[0]*rdxF[0]*f[5]+9.721111047611787*rdxF[0]*nu[1]*f[3]; 
  out[26] += 1.620185174601965*rdxF[1]*f[18]*nu[83]+2.371708245126284*rdxF[1]*f[26]*nu[73]+1.060660171779821*rdxF[1]*f[18]*nu[67]; 
  out[27] += 26.08879069638912*rdxF[1]*f[15]*nu[83]+14.8492424049175*f[1]*rdxF[1]*nu[83]+14.23024947075771*rdxF[1]*f[27]*nu[73]+18.1142209327368*rdxF[1]*f[5]*nu[73]+10.86853255964208*rdxF[1]*f[15]*nu[67]+9.721111047611787*f[1]*rdxF[1]*nu[67]+8.100925873009823*rdxF[1]*f[5]*nu[64]+2.371708245126284*rdxF[0]*nu[7]*f[27]+1.620185174601965*rdxF[0]*nu[17]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[19]; 
  out[28] += 26.08879069638912*rdxF[1]*f[16]*nu[83]+14.8492424049175*rdxF[1]*f[2]*nu[83]+14.23024947075771*rdxF[1]*f[28]*nu[73]+18.1142209327368*rdxF[1]*f[6]*nu[73]+10.86853255964208*rdxF[1]*f[16]*nu[67]+9.721111047611787*rdxF[1]*f[2]*nu[67]+8.100925873009823*rdxF[1]*f[6]*nu[64]; 
  out[29] += 1.620185174601965*rdxF[1]*f[23]*nu[83]+2.371708245126284*rdxF[1]*f[29]*nu[73]+1.060660171779821*rdxF[1]*f[23]*nu[67]+14.23024947075771*rdxF[0]*nu[7]*f[29]+26.08879069638912*rdxF[0]*nu[17]*f[20]+10.86853255964208*rdxF[0]*nu[1]*f[20]+14.8492424049175*rdxF[0]*f[6]*nu[17]+18.1142209327368*rdxF[0]*nu[7]*f[10]+8.100925873009823*nu[0]*rdxF[0]*f[10]+9.721111047611789*rdxF[0]*nu[1]*f[6]; 
  out[30] += 1.620185174601965*rdxF[1]*f[24]*nu[83]+2.371708245126284*rdxF[1]*f[30]*nu[73]+1.060660171779821*rdxF[1]*f[24]*nu[67]+2.371708245126284*rdxF[0]*nu[7]*f[30]+1.620185174601965*rdxF[0]*nu[17]*f[26]+1.060660171779821*rdxF[0]*nu[1]*f[26]; 
  out[31] += 26.08879069638912*rdxF[1]*f[22]*nu[83]+14.8492424049175*rdxF[1]*f[4]*nu[83]+14.23024947075771*rdxF[1]*f[31]*nu[73]+18.1142209327368*rdxF[1]*f[10]*nu[73]+10.86853255964208*rdxF[1]*f[22]*nu[67]+9.721111047611789*rdxF[1]*f[4]*nu[67]+8.100925873009823*rdxF[1]*f[10]*nu[64]+2.371708245126284*rdxF[0]*nu[7]*f[31]+1.620185174601965*rdxF[0]*nu[17]*f[28]+1.060660171779821*rdxF[0]*nu[1]*f[28]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+rdxF[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol3xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[96]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 3.622844186547359*rdxF[0]*f[8]*nu[50]+1.620185174601965*f[0]*rdxF[0]*nu[50]+2.371708245126284*rdxF[0]*f[2]*nu[40]+1.060660171779821*f[0]*rdxF[0]*nu[34]; 
  out[4] += 3.622844186547359*rdxF[0]*f[12]*nu[50]+1.620185174601965*rdxF[0]*f[1]*nu[50]+2.371708245126284*rdxF[0]*f[4]*nu[40]+1.060660171779821*rdxF[0]*f[1]*nu[34]; 
  out[6] += 3.622844186547359*rdxF[0]*f[14]*nu[50]+1.620185174601965*rdxF[0]*f[3]*nu[50]+2.371708245126284*rdxF[0]*f[6]*nu[40]+1.060660171779821*rdxF[0]*f[3]*nu[34]; 
  out[8] += 9.48683298050514*rdxF[0]*f[18]*nu[50]+10.86853255964208*rdxF[0]*f[2]*nu[50]+7.115124735378852*rdxF[0]*f[8]*nu[40]+5.303300858899105*f[0]*rdxF[0]*nu[40]+4.743416490252569*rdxF[0]*f[2]*nu[34]+2.371708245126284*f[0]*rdxF[0]*nu[32]; 
  out[10] += 3.622844186547359*rdxF[0]*f[21]*nu[50]+1.620185174601965*rdxF[0]*f[5]*nu[50]+2.371708245126284*rdxF[0]*f[10]*nu[40]+1.060660171779821*rdxF[0]*f[5]*nu[34]; 
  out[11] += 1.620185174601965*rdxF[0]*f[7]*nu[50]+2.371708245126284*rdxF[0]*f[11]*nu[40]+1.060660171779821*rdxF[0]*f[7]*nu[34]; 
  out[12] += 9.486832980505138*rdxF[0]*f[24]*nu[50]+10.86853255964208*rdxF[0]*f[4]*nu[50]+7.115124735378852*rdxF[0]*f[12]*nu[40]+5.303300858899106*rdxF[0]*f[1]*nu[40]+4.743416490252569*rdxF[0]*f[4]*nu[34]+2.371708245126284*rdxF[0]*f[1]*nu[32]; 
  out[14] += 9.486832980505138*rdxF[0]*f[26]*nu[50]+10.86853255964208*rdxF[0]*f[6]*nu[50]+7.115124735378852*rdxF[0]*f[14]*nu[40]+5.303300858899106*rdxF[0]*f[3]*nu[40]+4.743416490252569*rdxF[0]*f[6]*nu[34]+2.371708245126284*rdxF[0]*f[3]*nu[32]; 
  out[16] += 1.620185174601965*rdxF[0]*f[9]*nu[50]+2.371708245126284*rdxF[0]*f[16]*nu[40]+1.060660171779821*rdxF[0]*f[9]*nu[34]; 
  out[18] += 26.08879069638912*rdxF[0]*f[8]*nu[50]+14.8492424049175*f[0]*rdxF[0]*nu[50]+14.23024947075771*rdxF[0]*f[18]*nu[40]+18.1142209327368*rdxF[0]*f[2]*nu[40]+10.86853255964208*rdxF[0]*f[8]*nu[34]+9.721111047611789*f[0]*rdxF[0]*nu[34]+8.100925873009823*rdxF[0]*f[2]*nu[32]; 
  out[20] += 1.620185174601965*rdxF[0]*f[13]*nu[50]+2.371708245126284*rdxF[0]*f[20]*nu[40]+1.060660171779821*rdxF[0]*f[13]*nu[34]; 
  out[21] += 9.48683298050514*rdxF[0]*f[30]*nu[50]+10.86853255964208*rdxF[0]*f[10]*nu[50]+7.115124735378852*rdxF[0]*f[21]*nu[40]+5.303300858899105*rdxF[0]*f[5]*nu[40]+4.743416490252569*rdxF[0]*f[10]*nu[34]+2.371708245126284*rdxF[0]*f[5]*nu[32]; 
  out[22] += 1.620185174601965*rdxF[0]*f[15]*nu[50]+2.371708245126284*rdxF[0]*f[22]*nu[40]+1.060660171779821*rdxF[0]*f[15]*nu[34]; 
  out[23] += 1.620185174601965*rdxF[0]*f[17]*nu[50]+2.371708245126284*rdxF[0]*f[23]*nu[40]+1.060660171779821*rdxF[0]*f[17]*nu[34]; 
  out[24] += 26.08879069638912*rdxF[0]*f[12]*nu[50]+14.8492424049175*rdxF[0]*f[1]*nu[50]+14.23024947075771*rdxF[0]*f[24]*nu[40]+18.1142209327368*rdxF[0]*f[4]*nu[40]+10.86853255964208*rdxF[0]*f[12]*nu[34]+9.721111047611787*rdxF[0]*f[1]*nu[34]+8.100925873009823*rdxF[0]*f[4]*nu[32]; 
  out[26] += 26.08879069638912*rdxF[0]*f[14]*nu[50]+14.8492424049175*rdxF[0]*f[3]*nu[50]+14.23024947075771*rdxF[0]*f[26]*nu[40]+18.1142209327368*rdxF[0]*f[6]*nu[40]+10.86853255964208*rdxF[0]*f[14]*nu[34]+9.721111047611787*rdxF[0]*f[3]*nu[34]+8.100925873009823*rdxF[0]*f[6]*nu[32]; 
  out[28] += 1.620185174601965*rdxF[0]*f[19]*nu[50]+2.371708245126284*rdxF[0]*f[28]*nu[40]+1.060660171779821*rdxF[0]*f[19]*nu[34]; 
  out[29] += 1.620185174601965*rdxF[0]*f[25]*nu[50]+2.371708245126284*rdxF[0]*f[29]*nu[40]+1.060660171779821*rdxF[0]*f[25]*nu[34]; 
  out[30] += 26.08879069638912*rdxF[0]*f[21]*nu[50]+14.8492424049175*rdxF[0]*f[5]*nu[50]+14.23024947075771*rdxF[0]*f[30]*nu[40]+18.1142209327368*rdxF[0]*f[10]*nu[40]+10.86853255964208*rdxF[0]*f[21]*nu[34]+9.721111047611789*rdxF[0]*f[5]*nu[34]+8.100925873009823*rdxF[0]*f[10]*nu[32]; 
  out[31] += 1.620185174601965*rdxF[0]*f[27]*nu[50]+2.371708245126284*rdxF[0]*f[31]*nu[40]+1.060660171779821*rdxF[0]*f[27]*nu[34]; 

  return (rdxF[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol3xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[96]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[2] += 3.622844186547359*rdxF[0]*f[8]*nu[50]+1.620185174601965*f[0]*rdxF[0]*nu[50]+2.371708245126284*rdxF[0]*f[2]*nu[40]+1.060660171779821*f[0]*rdxF[0]*nu[34]; 
  out[3] += 3.622844186547359*rdxF[1]*f[9]*nu[83]+1.620185174601965*f[0]*rdxF[1]*nu[83]+2.371708245126284*rdxF[1]*f[3]*nu[73]+1.060660171779821*f[0]*rdxF[1]*nu[67]; 
  out[4] += 3.622844186547359*rdxF[0]*f[12]*nu[50]+1.620185174601965*rdxF[0]*f[1]*nu[50]+2.371708245126284*rdxF[0]*f[4]*nu[40]+1.060660171779821*rdxF[0]*f[1]*nu[34]; 
  out[5] += 3.622844186547359*rdxF[1]*f[15]*nu[83]+1.620185174601965*f[1]*rdxF[1]*nu[83]+2.371708245126284*rdxF[1]*f[5]*nu[73]+1.060660171779821*f[1]*rdxF[1]*nu[67]; 
  out[6] += 3.622844186547359*rdxF[1]*f[16]*nu[83]+1.620185174601965*rdxF[1]*f[2]*nu[83]+2.371708245126284*rdxF[1]*f[6]*nu[73]+1.060660171779821*rdxF[1]*f[2]*nu[67]+3.622844186547359*rdxF[0]*f[14]*nu[50]+1.620185174601965*rdxF[0]*f[3]*nu[50]+2.371708245126284*rdxF[0]*f[6]*nu[40]+1.060660171779821*rdxF[0]*f[3]*nu[34]; 
  out[8] += 9.48683298050514*rdxF[0]*f[18]*nu[50]+10.86853255964208*rdxF[0]*f[2]*nu[50]+7.115124735378852*rdxF[0]*f[8]*nu[40]+5.303300858899105*f[0]*rdxF[0]*nu[40]+4.743416490252569*rdxF[0]*f[2]*nu[34]+2.371708245126284*f[0]*rdxF[0]*nu[32]; 
  out[9] += 9.48683298050514*rdxF[1]*f[19]*nu[83]+10.86853255964208*rdxF[1]*f[3]*nu[83]+7.115124735378852*rdxF[1]*f[9]*nu[73]+5.303300858899105*f[0]*rdxF[1]*nu[73]+4.743416490252569*rdxF[1]*f[3]*nu[67]+2.371708245126284*f[0]*rdxF[1]*nu[64]; 
  out[10] += 3.622844186547359*rdxF[1]*f[22]*nu[83]+1.620185174601965*rdxF[1]*f[4]*nu[83]+2.371708245126284*rdxF[1]*f[10]*nu[73]+1.060660171779821*rdxF[1]*f[4]*nu[67]+3.622844186547359*rdxF[0]*f[21]*nu[50]+1.620185174601965*rdxF[0]*f[5]*nu[50]+2.371708245126284*rdxF[0]*f[10]*nu[40]+1.060660171779821*rdxF[0]*f[5]*nu[34]; 
  out[11] += 1.620185174601965*rdxF[0]*f[7]*nu[50]+2.371708245126284*rdxF[0]*f[11]*nu[40]+1.060660171779821*rdxF[0]*f[7]*nu[34]; 
  out[12] += 9.486832980505138*rdxF[0]*f[24]*nu[50]+10.86853255964208*rdxF[0]*f[4]*nu[50]+7.115124735378852*rdxF[0]*f[12]*nu[40]+5.303300858899106*rdxF[0]*f[1]*nu[40]+4.743416490252569*rdxF[0]*f[4]*nu[34]+2.371708245126284*rdxF[0]*f[1]*nu[32]; 
  out[13] += 1.620185174601965*rdxF[1]*f[7]*nu[83]+2.371708245126284*rdxF[1]*f[13]*nu[73]+1.060660171779821*rdxF[1]*f[7]*nu[67]; 
  out[14] += 1.620185174601965*rdxF[1]*f[8]*nu[83]+2.371708245126284*rdxF[1]*f[14]*nu[73]+1.060660171779821*rdxF[1]*f[8]*nu[67]+9.486832980505138*rdxF[0]*f[26]*nu[50]+10.86853255964208*rdxF[0]*f[6]*nu[50]+7.115124735378852*rdxF[0]*f[14]*nu[40]+5.303300858899106*rdxF[0]*f[3]*nu[40]+4.743416490252569*rdxF[0]*f[6]*nu[34]+2.371708245126284*rdxF[0]*f[3]*nu[32]; 
  out[15] += 9.486832980505138*rdxF[1]*f[27]*nu[83]+10.86853255964208*rdxF[1]*f[5]*nu[83]+7.115124735378852*rdxF[1]*f[15]*nu[73]+5.303300858899106*f[1]*rdxF[1]*nu[73]+4.743416490252569*rdxF[1]*f[5]*nu[67]+2.371708245126284*f[1]*rdxF[1]*nu[64]; 
  out[16] += 9.486832980505138*rdxF[1]*f[28]*nu[83]+10.86853255964208*rdxF[1]*f[6]*nu[83]+7.115124735378852*rdxF[1]*f[16]*nu[73]+5.303300858899106*rdxF[1]*f[2]*nu[73]+4.743416490252569*rdxF[1]*f[6]*nu[67]+2.371708245126284*rdxF[1]*f[2]*nu[64]+1.620185174601965*rdxF[0]*f[9]*nu[50]+2.371708245126284*rdxF[0]*f[16]*nu[40]+1.060660171779821*rdxF[0]*f[9]*nu[34]; 
  out[18] += 26.08879069638912*rdxF[0]*f[8]*nu[50]+14.8492424049175*f[0]*rdxF[0]*nu[50]+14.23024947075771*rdxF[0]*f[18]*nu[40]+18.1142209327368*rdxF[0]*f[2]*nu[40]+10.86853255964208*rdxF[0]*f[8]*nu[34]+9.721111047611789*f[0]*rdxF[0]*nu[34]+8.100925873009823*rdxF[0]*f[2]*nu[32]; 
  out[19] += 26.08879069638912*rdxF[1]*f[9]*nu[83]+14.8492424049175*f[0]*rdxF[1]*nu[83]+14.23024947075771*rdxF[1]*f[19]*nu[73]+18.1142209327368*rdxF[1]*f[3]*nu[73]+10.86853255964208*rdxF[1]*f[9]*nu[67]+9.721111047611789*f[0]*rdxF[1]*nu[67]+8.100925873009823*rdxF[1]*f[3]*nu[64]; 
  out[20] += 1.620185174601965*rdxF[1]*f[11]*nu[83]+2.371708245126284*rdxF[1]*f[20]*nu[73]+1.060660171779821*rdxF[1]*f[11]*nu[67]+1.620185174601965*rdxF[0]*f[13]*nu[50]+2.371708245126284*rdxF[0]*f[20]*nu[40]+1.060660171779821*rdxF[0]*f[13]*nu[34]; 
  out[21] += 1.620185174601965*rdxF[1]*f[12]*nu[83]+2.371708245126284*rdxF[1]*f[21]*nu[73]+1.060660171779821*rdxF[1]*f[12]*nu[67]+9.48683298050514*rdxF[0]*f[30]*nu[50]+10.86853255964208*rdxF[0]*f[10]*nu[50]+7.115124735378852*rdxF[0]*f[21]*nu[40]+5.303300858899105*rdxF[0]*f[5]*nu[40]+4.743416490252569*rdxF[0]*f[10]*nu[34]+2.371708245126284*rdxF[0]*f[5]*nu[32]; 
  out[22] += 9.48683298050514*rdxF[1]*f[31]*nu[83]+10.86853255964208*rdxF[1]*f[10]*nu[83]+7.115124735378852*rdxF[1]*f[22]*nu[73]+5.303300858899105*rdxF[1]*f[4]*nu[73]+4.743416490252569*rdxF[1]*f[10]*nu[67]+2.371708245126284*rdxF[1]*f[4]*nu[64]+1.620185174601965*rdxF[0]*f[15]*nu[50]+2.371708245126284*rdxF[0]*f[22]*nu[40]+1.060660171779821*rdxF[0]*f[15]*nu[34]; 
  out[23] += 1.620185174601965*rdxF[0]*f[17]*nu[50]+2.371708245126284*rdxF[0]*f[23]*nu[40]+1.060660171779821*rdxF[0]*f[17]*nu[34]; 
  out[24] += 26.08879069638912*rdxF[0]*f[12]*nu[50]+14.8492424049175*rdxF[0]*f[1]*nu[50]+14.23024947075771*rdxF[0]*f[24]*nu[40]+18.1142209327368*rdxF[0]*f[4]*nu[40]+10.86853255964208*rdxF[0]*f[12]*nu[34]+9.721111047611787*rdxF[0]*f[1]*nu[34]+8.100925873009823*rdxF[0]*f[4]*nu[32]; 
  out[25] += 1.620185174601965*rdxF[1]*f[17]*nu[83]+2.371708245126284*rdxF[1]*f[25]*nu[73]+1.060660171779821*rdxF[1]*f[17]*nu[67]; 
  out[26] += 1.620185174601965*rdxF[1]*f[18]*nu[83]+2.371708245126284*rdxF[1]*f[26]*nu[73]+1.060660171779821*rdxF[1]*f[18]*nu[67]+26.08879069638912*rdxF[0]*f[14]*nu[50]+14.8492424049175*rdxF[0]*f[3]*nu[50]+14.23024947075771*rdxF[0]*f[26]*nu[40]+18.1142209327368*rdxF[0]*f[6]*nu[40]+10.86853255964208*rdxF[0]*f[14]*nu[34]+9.721111047611787*rdxF[0]*f[3]*nu[34]+8.100925873009823*rdxF[0]*f[6]*nu[32]; 
  out[27] += 26.08879069638912*rdxF[1]*f[15]*nu[83]+14.8492424049175*f[1]*rdxF[1]*nu[83]+14.23024947075771*rdxF[1]*f[27]*nu[73]+18.1142209327368*rdxF[1]*f[5]*nu[73]+10.86853255964208*rdxF[1]*f[15]*nu[67]+9.721111047611787*f[1]*rdxF[1]*nu[67]+8.100925873009823*rdxF[1]*f[5]*nu[64]; 
  out[28] += 26.08879069638912*rdxF[1]*f[16]*nu[83]+14.8492424049175*rdxF[1]*f[2]*nu[83]+14.23024947075771*rdxF[1]*f[28]*nu[73]+18.1142209327368*rdxF[1]*f[6]*nu[73]+10.86853255964208*rdxF[1]*f[16]*nu[67]+9.721111047611787*rdxF[1]*f[2]*nu[67]+8.100925873009823*rdxF[1]*f[6]*nu[64]+1.620185174601965*rdxF[0]*f[19]*nu[50]+2.371708245126284*rdxF[0]*f[28]*nu[40]+1.060660171779821*rdxF[0]*f[19]*nu[34]; 
  out[29] += 1.620185174601965*rdxF[1]*f[23]*nu[83]+2.371708245126284*rdxF[1]*f[29]*nu[73]+1.060660171779821*rdxF[1]*f[23]*nu[67]+1.620185174601965*rdxF[0]*f[25]*nu[50]+2.371708245126284*rdxF[0]*f[29]*nu[40]+1.060660171779821*rdxF[0]*f[25]*nu[34]; 
  out[30] += 1.620185174601965*rdxF[1]*f[24]*nu[83]+2.371708245126284*rdxF[1]*f[30]*nu[73]+1.060660171779821*rdxF[1]*f[24]*nu[67]+26.08879069638912*rdxF[0]*f[21]*nu[50]+14.8492424049175*rdxF[0]*f[5]*nu[50]+14.23024947075771*rdxF[0]*f[30]*nu[40]+18.1142209327368*rdxF[0]*f[10]*nu[40]+10.86853255964208*rdxF[0]*f[21]*nu[34]+9.721111047611789*rdxF[0]*f[5]*nu[34]+8.100925873009823*rdxF[0]*f[10]*nu[32]; 
  out[31] += 26.08879069638912*rdxF[1]*f[22]*nu[83]+14.8492424049175*rdxF[1]*f[4]*nu[83]+14.23024947075771*rdxF[1]*f[31]*nu[73]+18.1142209327368*rdxF[1]*f[10]*nu[73]+10.86853255964208*rdxF[1]*f[22]*nu[67]+9.721111047611789*rdxF[1]*f[4]*nu[67]+8.100925873009823*rdxF[1]*f[10]*nu[64]+1.620185174601965*rdxF[0]*f[27]*nu[50]+2.371708245126284*rdxF[0]*f[31]*nu[40]+1.060660171779821*rdxF[0]*f[27]*nu[34]; 

  return (rdxF[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+rdxF[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol3xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[96]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 

  out[3] += 3.622844186547359*rdxF[0]*f[9]*nu[83]+1.620185174601965*f[0]*rdxF[0]*nu[83]+2.371708245126284*rdxF[0]*f[3]*nu[73]+1.060660171779821*f[0]*rdxF[0]*nu[67]; 
  out[5] += 3.622844186547359*rdxF[0]*f[15]*nu[83]+1.620185174601965*rdxF[0]*f[1]*nu[83]+2.371708245126284*rdxF[0]*f[5]*nu[73]+1.060660171779821*rdxF[0]*f[1]*nu[67]; 
  out[6] += 3.622844186547359*rdxF[0]*f[16]*nu[83]+1.620185174601965*rdxF[0]*f[2]*nu[83]+2.371708245126284*rdxF[0]*f[6]*nu[73]+1.060660171779821*rdxF[0]*f[2]*nu[67]; 
  out[9] += 9.48683298050514*rdxF[0]*f[19]*nu[83]+10.86853255964208*rdxF[0]*f[3]*nu[83]+7.115124735378852*rdxF[0]*f[9]*nu[73]+5.303300858899105*f[0]*rdxF[0]*nu[73]+4.743416490252569*rdxF[0]*f[3]*nu[67]+2.371708245126284*f[0]*rdxF[0]*nu[64]; 
  out[10] += 3.622844186547359*rdxF[0]*f[22]*nu[83]+1.620185174601965*rdxF[0]*f[4]*nu[83]+2.371708245126284*rdxF[0]*f[10]*nu[73]+1.060660171779821*rdxF[0]*f[4]*nu[67]; 
  out[13] += 1.620185174601965*rdxF[0]*f[7]*nu[83]+2.371708245126284*rdxF[0]*f[13]*nu[73]+1.060660171779821*rdxF[0]*f[7]*nu[67]; 
  out[14] += 1.620185174601965*rdxF[0]*f[8]*nu[83]+2.371708245126284*rdxF[0]*f[14]*nu[73]+1.060660171779821*rdxF[0]*f[8]*nu[67]; 
  out[15] += 9.486832980505138*rdxF[0]*f[27]*nu[83]+10.86853255964208*rdxF[0]*f[5]*nu[83]+7.115124735378852*rdxF[0]*f[15]*nu[73]+5.303300858899106*rdxF[0]*f[1]*nu[73]+4.743416490252569*rdxF[0]*f[5]*nu[67]+2.371708245126284*rdxF[0]*f[1]*nu[64]; 
  out[16] += 9.486832980505138*rdxF[0]*f[28]*nu[83]+10.86853255964208*rdxF[0]*f[6]*nu[83]+7.115124735378852*rdxF[0]*f[16]*nu[73]+5.303300858899106*rdxF[0]*f[2]*nu[73]+4.743416490252569*rdxF[0]*f[6]*nu[67]+2.371708245126284*rdxF[0]*f[2]*nu[64]; 
  out[19] += 26.08879069638912*rdxF[0]*f[9]*nu[83]+14.8492424049175*f[0]*rdxF[0]*nu[83]+14.23024947075771*rdxF[0]*f[19]*nu[73]+18.1142209327368*rdxF[0]*f[3]*nu[73]+10.86853255964208*rdxF[0]*f[9]*nu[67]+9.721111047611789*f[0]*rdxF[0]*nu[67]+8.100925873009823*rdxF[0]*f[3]*nu[64]; 
  out[20] += 1.620185174601965*rdxF[0]*f[11]*nu[83]+2.371708245126284*rdxF[0]*f[20]*nu[73]+1.060660171779821*rdxF[0]*f[11]*nu[67]; 
  out[21] += 1.620185174601965*rdxF[0]*f[12]*nu[83]+2.371708245126284*rdxF[0]*f[21]*nu[73]+1.060660171779821*rdxF[0]*f[12]*nu[67]; 
  out[22] += 9.48683298050514*rdxF[0]*f[31]*nu[83]+10.86853255964208*rdxF[0]*f[10]*nu[83]+7.115124735378852*rdxF[0]*f[22]*nu[73]+5.303300858899105*rdxF[0]*f[4]*nu[73]+4.743416490252569*rdxF[0]*f[10]*nu[67]+2.371708245126284*rdxF[0]*f[4]*nu[64]; 
  out[25] += 1.620185174601965*rdxF[0]*f[17]*nu[83]+2.371708245126284*rdxF[0]*f[25]*nu[73]+1.060660171779821*rdxF[0]*f[17]*nu[67]; 
  out[26] += 1.620185174601965*rdxF[0]*f[18]*nu[83]+2.371708245126284*rdxF[0]*f[26]*nu[73]+1.060660171779821*rdxF[0]*f[18]*nu[67]; 
  out[27] += 26.08879069638912*rdxF[0]*f[15]*nu[83]+14.8492424049175*rdxF[0]*f[1]*nu[83]+14.23024947075771*rdxF[0]*f[27]*nu[73]+18.1142209327368*rdxF[0]*f[5]*nu[73]+10.86853255964208*rdxF[0]*f[15]*nu[67]+9.721111047611787*rdxF[0]*f[1]*nu[67]+8.100925873009823*rdxF[0]*f[5]*nu[64]; 
  out[28] += 26.08879069638912*rdxF[0]*f[16]*nu[83]+14.8492424049175*rdxF[0]*f[2]*nu[83]+14.23024947075771*rdxF[0]*f[28]*nu[73]+18.1142209327368*rdxF[0]*f[6]*nu[73]+10.86853255964208*rdxF[0]*f[16]*nu[67]+9.721111047611787*rdxF[0]*f[2]*nu[67]+8.100925873009823*rdxF[0]*f[6]*nu[64]; 
  out[29] += 1.620185174601965*rdxF[0]*f[23]*nu[83]+2.371708245126284*rdxF[0]*f[29]*nu[73]+1.060660171779821*rdxF[0]*f[23]*nu[67]; 
  out[30] += 1.620185174601965*rdxF[0]*f[24]*nu[83]+2.371708245126284*rdxF[0]*f[30]*nu[73]+1.060660171779821*rdxF[0]*f[24]*nu[67]; 
  out[31] += 26.08879069638912*rdxF[0]*f[22]*nu[83]+14.8492424049175*rdxF[0]*f[4]*nu[83]+14.23024947075771*rdxF[0]*f[31]*nu[73]+18.1142209327368*rdxF[0]*f[10]*nu[73]+10.86853255964208*rdxF[0]*f[22]*nu[67]+9.721111047611789*rdxF[0]*f[4]*nu[67]+8.100925873009823*rdxF[0]*f[10]*nu[64]; 

  return (rdxF[0]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*1.142857142857143;

} 
void ConstDiffusionCFLfreqMin3xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin3xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin3xSerP3_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+kxSq[2]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin3xSerP3_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin3xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin3xSerP3_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+kxSq[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin3xSerP3_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP3_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+kxSq[2]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP3_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP3_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+kxSq[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP3_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP3_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+kxSq[2]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP3_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP3_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[41])-0.3952847075210473*nu[40]-0.3952847075210473*nu[39]+0.3535533905932737*nu[32])+kxSq[1]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP3_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[96]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[73])-0.3952847075210473*nu[72]-0.3952847075210473*nu[71]+0.3535533905932737*nu[64]))*2.285714285714286);

} 
