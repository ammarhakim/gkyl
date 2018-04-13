#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol2xSerP1(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[2]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 

  out[1] += 3.0*rdx[0]*f[1]*nu; 
  out[2] += 3.0*rdx[1]*f[2]*nu; 
  out[3] += 3.0*rdx[1]*f[3]*nu+3.0*rdx[0]*f[3]*nu; 

return std::nu*(rdx[0]+rdx[1])*0.5; 

} 
double ConstDiffusionVol2xSerP2(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[2]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 

  out[1] += 3.0*rdx[0]*f[1]*nu; 
  out[2] += 3.0*rdx[1]*f[2]*nu; 
  out[3] += 3.0*rdx[1]*f[3]*nu+3.0*rdx[0]*f[3]*nu; 
  out[4] += 15.0*rdx[0]*f[4]*nu; 
  out[5] += 15.0*rdx[1]*f[5]*nu; 
  out[6] += 3.0*rdx[1]*f[6]*nu+15.0*rdx[0]*f[6]*nu; 
  out[7] += 15.0*rdx[1]*f[7]*nu+3.0*rdx[0]*f[7]*nu; 

return std::nu*(rdx[0]+rdx[1])*0.5; 

} 
double ConstDiffusionVol2xSerP3(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[2]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 

  out[1] += 4.58257569495584*rdx[0]*f[8]*nu+3.0*rdx[0]*f[1]*nu; 
  out[2] += 4.58257569495584*rdx[1]*f[9]*nu+3.0*rdx[1]*f[2]*nu; 
  out[3] += 4.58257569495584*rdx[1]*f[11]*nu+4.58257569495584*rdx[0]*f[10]*nu+3.0*rdx[1]*f[3]*nu+3.0*rdx[0]*f[3]*nu; 
  out[4] += 15.0*rdx[0]*f[4]*nu; 
  out[5] += 15.0*rdx[1]*f[5]*nu; 
  out[6] += 3.0*rdx[1]*f[6]*nu+15.0*rdx[0]*f[6]*nu; 
  out[7] += 15.0*rdx[1]*f[7]*nu+3.0*rdx[0]*f[7]*nu; 
  out[8] += 42.0*rdx[0]*f[8]*nu+4.58257569495584*rdx[0]*f[1]*nu; 
  out[9] += 42.0*rdx[1]*f[9]*nu+4.58257569495584*rdx[1]*f[2]*nu; 
  out[10] += 3.0*rdx[1]*f[10]*nu+42.0*rdx[0]*f[10]*nu+4.58257569495584*rdx[0]*f[3]*nu; 
  out[11] += 42.0*rdx[1]*f[11]*nu+3.0*rdx[0]*f[11]*nu+4.58257569495584*rdx[1]*f[3]*nu; 

return std::nu*(rdx[0]+rdx[1])*0.5; 

} 
double ConstDiffusionVol2xSerP4(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[2]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 

  out[1] += 4.58257569495584*rdx[0]*f[8]*nu+3.0*rdx[0]*f[1]*nu; 
  out[2] += 4.58257569495584*rdx[1]*f[9]*nu+3.0*rdx[1]*f[2]*nu; 
  out[3] += 4.58257569495584*rdx[1]*f[12]*nu+4.58257569495584*rdx[0]*f[11]*nu+3.0*rdx[1]*f[3]*nu+3.0*rdx[0]*f[3]*nu; 
  out[4] += 20.12461179749811*rdx[0]*f[13]*nu+15.0*rdx[0]*f[4]*nu; 
  out[5] += 20.12461179749811*rdx[1]*f[14]*nu+15.0*rdx[1]*f[5]*nu; 
  out[6] += 20.1246117974981*rdx[0]*f[15]*nu+3.0*rdx[1]*f[6]*nu+15.0*rdx[0]*f[6]*nu; 
  out[7] += 20.1246117974981*rdx[1]*f[16]*nu+15.0*rdx[1]*f[7]*nu+3.0*rdx[0]*f[7]*nu; 
  out[8] += 42.0*rdx[0]*f[8]*nu+4.58257569495584*rdx[0]*f[1]*nu; 
  out[9] += 42.0*rdx[1]*f[9]*nu+4.58257569495584*rdx[1]*f[2]*nu; 
  out[10] += 15.0*rdx[1]*f[10]*nu+15.0*rdx[0]*f[10]*nu; 
  out[11] += 3.0*rdx[1]*f[11]*nu+42.0*rdx[0]*f[11]*nu+4.58257569495584*rdx[0]*f[3]*nu; 
  out[12] += 42.0*rdx[1]*f[12]*nu+3.0*rdx[0]*f[12]*nu+4.58257569495584*rdx[1]*f[3]*nu; 
  out[13] += 90.0*rdx[0]*f[13]*nu+20.12461179749811*rdx[0]*f[4]*nu; 
  out[14] += 90.0*rdx[1]*f[14]*nu+20.12461179749811*rdx[1]*f[5]*nu; 
  out[15] += 3.0*rdx[1]*f[15]*nu+90.0*rdx[0]*f[15]*nu+20.1246117974981*rdx[0]*f[6]*nu; 
  out[16] += 90.0*rdx[1]*f[16]*nu+3.0*rdx[0]*f[16]*nu+20.1246117974981*rdx[1]*f[7]*nu; 

return std::nu*(rdx[0]+rdx[1])*0.5; 

} 
