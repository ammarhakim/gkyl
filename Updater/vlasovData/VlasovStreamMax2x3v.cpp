#include <VlasovModDecl.h> 
void VlasovVolStream2x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
} 
void VlasovVolStream2x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[6] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[9]*dv1dx1+f[8]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w0dx0+0.8944271909999159*f[18]*dv0dx0+f[0]*dv0dx0; 
  out[8] += 3.464101615137754*f[3]*w1dx1+f[11]*dv1dx1; 
  out[9] += 3.464101615137754*f[4]*w0dx0+f[11]*dv0dx0; 
  out[10] += 3.464101615137754*f[4]*w1dx1+0.8944271909999159*f[19]*dv1dx1+f[0]*dv1dx1; 
  out[12] += 3.464101615137754*f[5]*w0dx0+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[5]*w1dx1+f[15]*dv1dx1; 
  out[16] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[7]*dv0dx0; 
  out[17] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[10]*dv1dx1; 
} 
