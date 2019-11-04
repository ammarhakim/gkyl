#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dxv[3]: Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 


  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5*1.333333333333333; 

} 
double ConstDiffusionVol3xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dxv[3]: Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

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

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5*1.8; 

} 
double ConstDiffusionVol3xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dxv[3]: Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[17] += 22.9128784747792*rdxSq4nu[0]*f[1]; 
  out[18] += 22.9128784747792*rdxSq4nu[1]*f[2]; 
  out[19] += 22.9128784747792*rdxSq4nu[2]*f[3]; 
  out[20] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[21] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[22] += 6.708203932499369*rdxSq4nu[2]*f[4]; 
  out[23] += 22.9128784747792*rdxSq4nu[0]*f[4]; 
  out[24] += 22.9128784747792*rdxSq4nu[1]*f[4]; 
  out[25] += 22.9128784747792*rdxSq4nu[0]*f[5]; 
  out[26] += 22.9128784747792*rdxSq4nu[1]*f[6]; 
  out[27] += 22.9128784747792*rdxSq4nu[2]*f[5]; 
  out[28] += 22.9128784747792*rdxSq4nu[2]*f[6]; 
  out[29] += 22.9128784747792*rdxSq4nu[0]*f[10]; 
  out[30] += 22.9128784747792*rdxSq4nu[1]*f[10]; 
  out[31] += 22.9128784747792*rdxSq4nu[2]*f[10]; 

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5*2.285714285714286; 

} 
