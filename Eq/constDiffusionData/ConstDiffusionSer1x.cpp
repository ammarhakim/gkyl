#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dxv[1]: Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 


  return rdxSq4nu[0]*0.5*1.333333333333333; 

} 
double ConstDiffusionVol1xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dxv[1]: Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[2] += 6.708203932499369*f[0]*rdxSq4nu[0]; 

  return rdxSq4nu[0]*0.5*1.8; 

} 
double ConstDiffusionVol1xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dxv[1]: Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[2] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[3] += 22.9128784747792*rdxSq4nu[0]*f[1]; 

  return rdxSq4nu[0]*0.5*2.285714285714286; 

} 
