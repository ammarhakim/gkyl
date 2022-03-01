#include <VlasovModDecl.h> 
__host__ __device__ double VlasovGenGeoVol1x3vSerP1(const double *w, const double *dxv, const double *alphaGeo, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alphaGeo:  Components for tangent basis vectors.
  // f:         Input distribution function.
  // out:       Incremented output.

  const double *alpha0 = &alphaGeo[0]; 

  double alpha_mid = 0.0; 

  alpha_mid += std::abs(0.25*alpha0[1]); 

  out[1] += 0.4330127018922193*(f[15]*alpha0[16]+f[14]*alpha0[15]+f[13]*alpha0[14]+f[12]*alpha0[13]+f[11]*alpha0[12]+f[10]*alpha0[11]+f[9]*alpha0[10]+f[8]*alpha0[9]+f[7]*alpha0[8]+f[6]*alpha0[7]+f[5]*alpha0[6]+f[4]*alpha0[5]+f[3]*alpha0[4]+f[2]*alpha0[3]+f[1]*alpha0[2]+f[0]*alpha0[1]); 
  out[5] += 0.4330127018922193*(f[13]*alpha0[16]+alpha0[14]*f[15]+f[10]*alpha0[15]+alpha0[11]*f[14]+f[8]*alpha0[13]+alpha0[9]*f[12]+f[6]*alpha0[12]+alpha0[7]*f[11]+f[4]*alpha0[10]+alpha0[5]*f[9]+f[3]*alpha0[8]+alpha0[4]*f[7]+f[1]*alpha0[6]+alpha0[2]*f[5]+f[0]*alpha0[3]+alpha0[1]*f[2]); 
  out[6] += 0.4330127018922193*(f[12]*alpha0[16]+alpha0[13]*f[15]+f[9]*alpha0[15]+alpha0[10]*f[14]+f[8]*alpha0[14]+alpha0[9]*f[13]+f[5]*alpha0[12]+alpha0[6]*f[11]+f[4]*alpha0[11]+alpha0[5]*f[10]+f[2]*alpha0[8]+alpha0[3]*f[7]+f[1]*alpha0[7]+alpha0[2]*f[6]+f[0]*alpha0[4]+alpha0[1]*f[3]); 
  out[8] += 0.4330127018922193*(f[11]*alpha0[16]+alpha0[12]*f[15]+f[7]*alpha0[15]+alpha0[8]*f[14]+f[6]*alpha0[14]+alpha0[7]*f[13]+f[5]*alpha0[13]+alpha0[6]*f[12]+f[3]*alpha0[11]+alpha0[4]*f[10]+f[2]*alpha0[10]+alpha0[3]*f[9]+f[1]*alpha0[9]+alpha0[2]*f[8]+f[0]*alpha0[5]+alpha0[1]*f[4]); 
  out[11] += 0.4330127018922193*(f[8]*alpha0[16]+alpha0[9]*f[15]+f[4]*alpha0[15]+alpha0[5]*f[14]+f[12]*alpha0[14]+alpha0[13]*f[13]+f[1]*alpha0[12]+alpha0[2]*f[11]+f[9]*alpha0[11]+alpha0[10]*f[10]+f[0]*alpha0[8]+alpha0[1]*f[7]+f[5]*alpha0[7]+alpha0[6]*f[6]+f[2]*alpha0[4]+alpha0[3]*f[3]); 
  out[12] += 0.4330127018922193*(f[6]*alpha0[16]+alpha0[7]*f[15]+f[3]*alpha0[15]+alpha0[4]*f[14]+f[11]*alpha0[14]+alpha0[12]*f[13]+f[1]*alpha0[13]+alpha0[2]*f[12]+f[7]*alpha0[11]+alpha0[8]*f[10]+f[0]*alpha0[10]+alpha0[1]*f[9]+f[5]*alpha0[9]+alpha0[6]*f[8]+f[2]*alpha0[5]+alpha0[3]*f[4]); 
  out[13] += 0.4330127018922193*(f[5]*alpha0[16]+alpha0[6]*f[15]+f[2]*alpha0[15]+alpha0[3]*f[14]+f[1]*alpha0[14]+alpha0[2]*f[13]+f[11]*alpha0[13]+alpha0[12]*f[12]+f[0]*alpha0[11]+alpha0[1]*f[10]+f[7]*alpha0[10]+alpha0[8]*f[9]+f[6]*alpha0[9]+alpha0[7]*f[8]+f[3]*alpha0[5]+alpha0[4]*f[4]); 
  out[15] += 0.4330127018922193*(f[1]*alpha0[16]+alpha0[2]*f[15]+f[0]*alpha0[15]+alpha0[1]*f[14]+f[5]*alpha0[14]+alpha0[6]*f[13]+f[6]*alpha0[13]+alpha0[7]*f[12]+f[8]*alpha0[12]+alpha0[9]*f[11]+f[2]*alpha0[11]+alpha0[3]*f[10]+f[3]*alpha0[10]+alpha0[4]*f[9]+f[4]*alpha0[8]+alpha0[5]*f[7]); 

  return alpha_mid; 
} 
