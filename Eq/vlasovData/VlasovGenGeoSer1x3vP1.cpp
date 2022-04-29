#include <VlasovModDecl.h> 
__host__ __device__ double VlasovGenGeoVol1x3vSerP1(const double *w, const double *dxv, const double *alphaGeo, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alphaGeo:  Components for tangent basis vectors.
  // f:         Input distribution function.
  // out:       Incremented output.

  double w1Ddx0  = w[1]/dxv[0]; 
  double dv1Ddx0 = dxv[1]/dxv[0]; 
  const double *alpha0 = &alphaGeo[0]; 
  const double *alpha1 = &alphaGeo[16]; 
  const double *alpha2 = &alphaGeo[32]; 

  double alpha_mid = 0.0; 

  alpha_mid += 2*std::abs((0.25*alpha0[0])/dxv[0]); 

  out[1] += (0.8660254037844386*(alpha0[15]*f[15]+alpha0[14]*f[14]+alpha0[13]*f[13]+alpha0[12]*f[12]+alpha0[11]*f[11]+alpha0[10]*f[10]+alpha0[9]*f[9]+alpha0[8]*f[8]+alpha0[7]*f[7]+alpha0[6]*f[6]+alpha0[5]*f[5]+alpha0[4]*f[4]+alpha0[3]*f[3]+alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]))/dxv[0]; 
  out[5] += (0.8660254037844386*(alpha0[13]*f[15]+f[13]*alpha0[15]+alpha0[10]*f[14]+f[10]*alpha0[14]+alpha0[8]*f[12]+f[8]*alpha0[12]+alpha0[6]*f[11]+f[6]*alpha0[11]+alpha0[4]*f[9]+f[4]*alpha0[9]+alpha0[3]*f[7]+f[3]*alpha0[7]+alpha0[1]*f[5]+f[1]*alpha0[5]+alpha0[0]*f[2]+f[0]*alpha0[2]))/dxv[0]; 
  out[6] += (0.8660254037844386*(alpha0[12]*f[15]+f[12]*alpha0[15]+alpha0[9]*f[14]+f[9]*alpha0[14]+alpha0[8]*f[13]+f[8]*alpha0[13]+alpha0[5]*f[11]+f[5]*alpha0[11]+alpha0[4]*f[10]+f[4]*alpha0[10]+alpha0[2]*f[7]+f[2]*alpha0[7]+alpha0[1]*f[6]+f[1]*alpha0[6]+alpha0[0]*f[3]+f[0]*alpha0[3]))/dxv[0]; 
  out[8] += (0.8660254037844386*(alpha0[11]*f[15]+f[11]*alpha0[15]+alpha0[7]*f[14]+f[7]*alpha0[14]+alpha0[6]*f[13]+f[6]*alpha0[13]+alpha0[5]*f[12]+f[5]*alpha0[12]+alpha0[3]*f[10]+f[3]*alpha0[10]+alpha0[2]*f[9]+f[2]*alpha0[9]+alpha0[1]*f[8]+f[1]*alpha0[8]+alpha0[0]*f[4]+f[0]*alpha0[4]))/dxv[0]; 
  out[11] += (0.8660254037844386*(alpha0[8]*f[15]+f[8]*alpha0[15]+alpha0[4]*f[14]+f[4]*alpha0[14]+alpha0[12]*f[13]+f[12]*alpha0[13]+alpha0[1]*f[11]+f[1]*alpha0[11]+alpha0[9]*f[10]+f[9]*alpha0[10]+alpha0[0]*f[7]+f[0]*alpha0[7]+alpha0[5]*f[6]+f[5]*alpha0[6]+alpha0[2]*f[3]+f[2]*alpha0[3]))/dxv[0]; 
  out[12] += (0.8660254037844386*(alpha0[6]*f[15]+f[6]*alpha0[15]+alpha0[3]*f[14]+f[3]*alpha0[14]+alpha0[11]*f[13]+f[11]*alpha0[13]+alpha0[1]*f[12]+f[1]*alpha0[12]+alpha0[7]*f[10]+f[7]*alpha0[10]+alpha0[0]*f[9]+f[0]*alpha0[9]+alpha0[5]*f[8]+f[5]*alpha0[8]+alpha0[2]*f[4]+f[2]*alpha0[4]))/dxv[0]; 
  out[13] += (0.8660254037844386*(alpha0[5]*f[15]+f[5]*alpha0[15]+alpha0[2]*f[14]+f[2]*alpha0[14]+alpha0[1]*f[13]+f[1]*alpha0[13]+alpha0[11]*f[12]+f[11]*alpha0[12]+alpha0[0]*f[10]+f[0]*alpha0[10]+alpha0[7]*f[9]+f[7]*alpha0[9]+alpha0[6]*f[8]+f[6]*alpha0[8]+alpha0[3]*f[4]+f[3]*alpha0[4]))/dxv[0]; 
  out[15] += (0.8660254037844386*(alpha0[1]*f[15]+f[1]*alpha0[15]+alpha0[0]*f[14]+f[0]*alpha0[14]+alpha0[5]*f[13]+f[5]*alpha0[13]+alpha0[6]*f[12]+f[6]*alpha0[12]+alpha0[8]*f[11]+f[8]*alpha0[11]+alpha0[2]*f[10]+f[2]*alpha0[10]+alpha0[3]*f[9]+f[3]*alpha0[9]+alpha0[4]*f[7]+f[4]*alpha0[7]))/dxv[0]; 

  return alpha_mid; 
} 
