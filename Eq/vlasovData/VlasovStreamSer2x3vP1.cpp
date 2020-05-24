#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVolStream2x3vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 

  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[6] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[9]*dv1dx1+f[8]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w0dx0+f[0]*dv0dx0; 
  out[8] += 3.464101615137754*f[3]*w1dx1+f[11]*dv1dx1; 
  out[9] += 3.464101615137754*f[4]*w0dx0+f[11]*dv0dx0; 
  out[10] += 3.464101615137754*f[4]*w1dx1+f[0]*dv1dx1; 
  out[12] += 3.464101615137754*f[5]*w0dx0+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[5]*w1dx1+f[15]*dv1dx1; 
  out[16] += 3.464101615137754*f[7]*w1dx1+3.464101615137754*f[8]*w0dx0+f[18]*dv1dx1+f[2]*dv0dx0; 
  out[17] += 3.464101615137754*f[9]*w1dx1+3.464101615137754*f[10]*w0dx0+f[1]*dv1dx1+f[19]*dv0dx0; 
  out[18] += 3.464101615137754*f[11]*w0dx0+f[4]*dv0dx0; 
  out[19] += 3.464101615137754*f[11]*w1dx1+f[3]*dv1dx1; 
  out[20] += 3.464101615137754*f[12]*w1dx1+3.464101615137754*f[13]*w0dx0+f[23]*dv1dx1+f[22]*dv0dx0; 
  out[21] += 3.464101615137754*f[14]*w0dx0+f[5]*dv0dx0; 
  out[22] += 3.464101615137754*f[14]*w1dx1+f[25]*dv1dx1; 
  out[23] += 3.464101615137754*f[15]*w0dx0+f[25]*dv0dx0; 
  out[24] += 3.464101615137754*f[15]*w1dx1+f[5]*dv1dx1; 
  out[26] += 3.464101615137754*f[18]*w1dx1+3.464101615137754*f[19]*w0dx0+f[7]*dv1dx1+f[10]*dv0dx0; 
  out[27] += 3.464101615137754*f[21]*w1dx1+3.464101615137754*f[22]*w0dx0+f[29]*dv1dx1+f[13]*dv0dx0; 
  out[28] += 3.464101615137754*f[23]*w1dx1+3.464101615137754*f[24]*w0dx0+f[12]*dv1dx1+f[30]*dv0dx0; 
  out[29] += 3.464101615137754*f[25]*w0dx0+f[15]*dv0dx0; 
  out[30] += 3.464101615137754*f[25]*w1dx1+f[14]*dv1dx1; 
  out[31] += 3.464101615137754*f[29]*w1dx1+3.464101615137754*f[30]*w0dx0+f[21]*dv1dx1+f[24]*dv0dx0; 

  return std::abs(w0dx0)+std::abs(w1dx1)+dv0dx0/2+dv1dx1/2; 
} 
