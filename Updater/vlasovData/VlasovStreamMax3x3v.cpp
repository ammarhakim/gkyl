#include <VlasovModDecl.h> 
void VlasovVolStream3x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[4]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[5]*dv1dx1; 
  out[3] += 3.464101615137754*f[0]*w2dx2+f[6]*dv2dx2; 
} 
void VlasovVolStream3x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[4]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[5]*dv1dx1; 
  out[3] += 3.464101615137754*f[0]*w2dx2+f[6]*dv2dx2; 
  out[7] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[13]*dv1dx1+f[11]*dv0dx0; 
  out[8] += 3.464101615137754*f[1]*w2dx2+3.464101615137754*f[3]*w0dx0+f[17]*dv2dx2+f[12]*dv0dx0; 
  out[9] += 3.464101615137754*f[2]*w2dx2+3.464101615137754*f[3]*w1dx1+f[18]*dv2dx2+f[15]*dv1dx1; 
  out[10] += 3.464101615137754*f[4]*w0dx0+0.8944271909999159*f[25]*dv0dx0+f[0]*dv0dx0; 
  out[11] += 3.464101615137754*f[4]*w1dx1+f[16]*dv1dx1; 
  out[12] += 3.464101615137754*f[4]*w2dx2+f[20]*dv2dx2; 
  out[13] += 3.464101615137754*f[5]*w0dx0+f[16]*dv0dx0; 
  out[14] += 3.464101615137754*f[5]*w1dx1+0.8944271909999159*f[26]*dv1dx1+f[0]*dv1dx1; 
  out[15] += 3.464101615137754*f[5]*w2dx2+f[21]*dv2dx2; 
  out[17] += 3.464101615137754*f[6]*w0dx0+f[20]*dv0dx0; 
  out[18] += 3.464101615137754*f[6]*w1dx1+f[21]*dv1dx1; 
  out[19] += 3.464101615137754*f[6]*w2dx2+0.8944271909999159*f[27]*dv2dx2+f[0]*dv2dx2; 
  out[22] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[23] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[14]*dv1dx1; 
  out[24] += 7.745966692414834*f[3]*w2dx2+2.23606797749979*f[19]*dv2dx2; 
} 
