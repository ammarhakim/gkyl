#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xSerP1(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[1]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 

  out[1] += 3.0*rdx[0]*f[1]*nu; 

return std::nu*rdx[0]*0.5; 

} 
double ConstDiffusionVol1xSerP2(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[1]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 

  out[1] += 3.0*rdx[0]*f[1]*nu; 
  out[2] += 15.0*rdx[0]*f[2]*nu; 

return std::nu*rdx[0]*0.5; 

} 
double ConstDiffusionVol1xSerP3(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[1]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 

  out[1] += 4.58257569495584*rdx[0]*f[3]*nu+3.0*rdx[0]*f[1]*nu; 
  out[2] += 15.0*rdx[0]*f[2]*nu; 
  out[3] += 42.0*rdx[0]*f[3]*nu+4.58257569495584*rdx[0]*f[1]*nu; 

return std::nu*rdx[0]*0.5; 

} 
double ConstDiffusionVol1xSerP4(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[1]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 

  out[1] += 4.58257569495584*rdx[0]*f[3]*nu+3.0*rdx[0]*f[1]*nu; 
  out[2] += 20.12461179749811*rdx[0]*f[4]*nu+15.0*rdx[0]*f[2]*nu; 
  out[3] += 42.0*rdx[0]*f[3]*nu+4.58257569495584*rdx[0]*f[1]*nu; 
  out[4] += 90.0*rdx[0]*f[4]*nu+20.12461179749811*rdx[0]*f[2]*nu; 

return std::nu*rdx[0]*0.5; 

} 
