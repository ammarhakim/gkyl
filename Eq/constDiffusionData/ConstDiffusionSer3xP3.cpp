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
