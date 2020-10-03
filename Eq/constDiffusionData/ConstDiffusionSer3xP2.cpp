#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[17] += 6.708203932499369*rdxSq4nu[0]*f[6]; 

  return (rdxSq4nu[0])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[17] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[18] += 6.708203932499369*rdxSq4nu[1]*f[5]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[17] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[18] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[19] += 6.708203932499369*rdxSq4nu[2]*f[4]; 

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[2]*dx[2]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[16] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[17] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[19] += 6.708203932499369*rdxSq4nu[1]*f[4]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 

  out[8] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[12] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[18] += 6.708203932499369*rdxSq4nu[0]*f[5]; 

  return (rdxSq4nu[0])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[2]*dx[2]); 

  out[8] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[12] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[16] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[18] += 6.708203932499369*rdxSq4nu[0]*f[5]; 
  out[19] += 6.708203932499369*rdxSq4nu[1]*f[4]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[2]*dx[2]); 

  out[9] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[15] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[16] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[4]; 

  return (rdxSq4nu[0])*0.9;

} 
