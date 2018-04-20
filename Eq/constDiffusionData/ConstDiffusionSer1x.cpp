#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[1] += -3.0*rdxSq4nu[0]*f[1]; 

return rdxSq4nu[0]*0.5; 

} 
double ConstDiffusionVol1xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[1] += -3.0*rdxSq4nu[0]*f[1]; 
  out[2] += -15.0*rdxSq4nu[0]*f[2]; 

return rdxSq4nu[0]*0.5; 

} 
double ConstDiffusionVol1xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[1] += (-4.58257569495584*rdxSq4nu[0]*f[3])-3.0*rdxSq4nu[0]*f[1]; 
  out[2] += -15.0*rdxSq4nu[0]*f[2]; 
  out[3] += (-42.0*rdxSq4nu[0]*f[3])-4.58257569495584*rdxSq4nu[0]*f[1]; 

return rdxSq4nu[0]*0.5; 

} 
double ConstDiffusionVol1xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 

  out[1] += (-4.58257569495584*rdxSq4nu[0]*f[3])-3.0*rdxSq4nu[0]*f[1]; 
  out[2] += (-20.12461179749811*rdxSq4nu[0]*f[4])-15.0*rdxSq4nu[0]*f[2]; 
  out[3] += (-42.0*rdxSq4nu[0]*f[3])-4.58257569495584*rdxSq4nu[0]*f[1]; 
  out[4] += (-90.0*rdxSq4nu[0]*f[4])-20.12461179749811*rdxSq4nu[0]*f[2]; 

return rdxSq4nu[0]*0.5; 

} 
