#include <VlasovModDecl.h> 
double VlasovVolStream1x3vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 

  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[5] += 3.464101615137754*f[2]*w0dx0+f[0]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+f[7]*dv0dx0; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[9]*dv0dx0; 
  out[11] += 3.464101615137754*f[7]*w0dx0+f[3]*dv0dx0; 
  out[12] += 3.464101615137754*f[9]*w0dx0+f[4]*dv0dx0; 
  out[13] += 3.464101615137754*f[10]*w0dx0+f[14]*dv0dx0; 
  out[15] += 3.464101615137754*f[14]*w0dx0+f[10]*dv0dx0; 
  return std::abs(w0dx0)+dv0dx0/2; 
} 
