#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol2xMaxP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[4] += 6.708203932499369*f[0]*rdxSq4nu[0]; 

  return (rdxSq4nu[0])*0.9;

} 
double ConstDiffusionVol2xMaxP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[4] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[5] += 6.708203932499369*f[0]*rdxSq4nu[1]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol2xMaxP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 

  out[5] += 6.708203932499369*f[0]*rdxSq4nu[0]; 

  return (rdxSq4nu[0])*0.9;

} 
