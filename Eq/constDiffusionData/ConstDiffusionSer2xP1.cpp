#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol2xSerP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 


  return (rdxSq4nu[0])*0.6666666666666666;

} 
double ConstDiffusionVol2xSerP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxSq4nu[0]+rdxSq4nu[1])*0.6666666666666666;

} 
double ConstDiffusionVol2xSerP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 


  return (rdxSq4nu[0])*0.6666666666666666;

} 
