#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 


return rdxSq4nu[0]*0.5; 

} 
double ConstDiffusionVol1xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[2] += 6.708203932499369*f[0]*rdxSq4nu[0]; 

return rdxSq4nu[0]*0.5; 

} 
double ConstDiffusionVol1xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[2] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[3] += 22.9128784747792*rdxSq4nu[0]*f[1]; 

return rdxSq4nu[0]*0.5; 

} 
double ConstDiffusionVol1xMaxP4(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[2] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[3] += 22.9128784747792*rdxSq4nu[0]*f[1]; 
  out[4] += 46.95742752749558*rdxSq4nu[0]*f[2]+30.0*f[0]*rdxSq4nu[0]; 

return rdxSq4nu[0]*0.5; 

} 
