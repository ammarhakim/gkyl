#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol2xMaxP1(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 

  out[1] += -3.0*rdxSq4nu[0]*f[1]; 
  out[2] += -3.0*rdxSq4nu[1]*f[2]; 

return (rdxSq4nu[0]+rdxSq4nu[1])*0.5; 

} 
double ConstDiffusionVol2xMaxP2(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 

  out[1] += -3.0*rdxSq4nu[0]*f[1]; 
  out[2] += -3.0*rdxSq4nu[1]*f[2]; 
  out[3] += (-3.0*rdxSq4nu[1]*f[3])-3.0*rdxSq4nu[0]*f[3]; 
  out[4] += -15.0*rdxSq4nu[0]*f[4]; 
  out[5] += -15.0*rdxSq4nu[1]*f[5]; 

return (rdxSq4nu[0]+rdxSq4nu[1])*0.5; 

} 
double ConstDiffusionVol2xMaxP3(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 

  out[1] += (-4.58257569495584*rdxSq4nu[0]*f[8])-3.0*rdxSq4nu[0]*f[1]; 
  out[2] += (-4.58257569495584*rdxSq4nu[1]*f[9])-3.0*rdxSq4nu[1]*f[2]; 
  out[3] += (-3.0*rdxSq4nu[1]*f[3])-3.0*rdxSq4nu[0]*f[3]; 
  out[4] += -15.0*rdxSq4nu[0]*f[4]; 
  out[5] += -15.0*rdxSq4nu[1]*f[5]; 
  out[6] += (-3.0*rdxSq4nu[1]*f[6])-15.0*rdxSq4nu[0]*f[6]; 
  out[7] += (-15.0*rdxSq4nu[1]*f[7])-3.0*rdxSq4nu[0]*f[7]; 
  out[8] += (-42.0*rdxSq4nu[0]*f[8])-4.58257569495584*rdxSq4nu[0]*f[1]; 
  out[9] += (-42.0*rdxSq4nu[1]*f[9])-4.58257569495584*rdxSq4nu[1]*f[2]; 

return (rdxSq4nu[0]+rdxSq4nu[1])*0.5; 

} 
double ConstDiffusionVol2xMaxP4(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 

  out[1] += (-4.58257569495584*rdxSq4nu[0]*f[8])-3.0*rdxSq4nu[0]*f[1]; 
  out[2] += (-4.58257569495584*rdxSq4nu[1]*f[9])-3.0*rdxSq4nu[1]*f[2]; 
  out[3] += (-4.58257569495584*rdxSq4nu[1]*f[12])-4.58257569495584*rdxSq4nu[0]*f[11]-3.0*rdxSq4nu[1]*f[3]-3.0*rdxSq4nu[0]*f[3]; 
  out[4] += (-20.12461179749811*rdxSq4nu[0]*f[13])-15.0*rdxSq4nu[0]*f[4]; 
  out[5] += (-20.12461179749811*rdxSq4nu[1]*f[14])-15.0*rdxSq4nu[1]*f[5]; 
  out[6] += (-3.0*rdxSq4nu[1]*f[6])-15.0*rdxSq4nu[0]*f[6]; 
  out[7] += (-15.0*rdxSq4nu[1]*f[7])-3.0*rdxSq4nu[0]*f[7]; 
  out[8] += (-42.0*rdxSq4nu[0]*f[8])-4.58257569495584*rdxSq4nu[0]*f[1]; 
  out[9] += (-42.0*rdxSq4nu[1]*f[9])-4.58257569495584*rdxSq4nu[1]*f[2]; 
  out[10] += (-15.0*rdxSq4nu[1]*f[10])-15.0*rdxSq4nu[0]*f[10]; 
  out[11] += (-3.0*rdxSq4nu[1]*f[11])-42.0*rdxSq4nu[0]*f[11]-4.58257569495584*rdxSq4nu[0]*f[3]; 
  out[12] += (-42.0*rdxSq4nu[1]*f[12])-3.0*rdxSq4nu[0]*f[12]-4.58257569495584*rdxSq4nu[1]*f[3]; 
  out[13] += (-90.0*rdxSq4nu[0]*f[13])-20.12461179749811*rdxSq4nu[0]*f[4]; 
  out[14] += (-90.0*rdxSq4nu[1]*f[14])-20.12461179749811*rdxSq4nu[1]*f[5]; 

return (rdxSq4nu[0]+rdxSq4nu[1])*0.5; 

} 
