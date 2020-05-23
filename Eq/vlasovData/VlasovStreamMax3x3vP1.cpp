#include <VlasovModDecl.h> 
double VlasovVolStream3x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
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
  return std::abs(w0dx0)+std::abs(w1dx1)+std::abs(w2dx2)+dv0dx0/2+dv1dx1/2+dv2dx2/2; 
} 
