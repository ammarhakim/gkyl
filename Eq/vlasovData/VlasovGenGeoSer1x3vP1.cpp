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
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 

  out[1] += (0.8660254037844386*(alpha0[15]*f[15]+alpha0[14]*f[14]+alpha0[13]*f[13]+alpha0[12]*f[12]+alpha0[11]*f[11]+alpha0[10]*f[10]+alpha0[9]*f[9]+alpha0[8]*f[8]+alpha0[7]*f[7]+alpha0[6]*f[6]+alpha0[5]*f[5]+alpha0[4]*f[4]+alpha0[3]*f[3]+alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]))/dxv[0]; 
  out[5] += (0.8660254037844386*(alpha0[13]*f[15]+f[13]*alpha0[15]+alpha0[10]*f[14]+f[10]*alpha0[14]+alpha0[8]*f[12]+f[8]*alpha0[12]+alpha0[6]*f[11]+f[6]*alpha0[11]+alpha0[4]*f[9]+f[4]*alpha0[9]+alpha0[3]*f[7]+f[3]*alpha0[7]+alpha0[1]*f[5]+f[1]*alpha0[5]+alpha0[0]*f[2]+f[0]*alpha0[2]))/dxv[0]; 
  out[6] += (0.8660254037844386*(alpha0[12]*f[15]+f[12]*alpha0[15]+alpha0[9]*f[14]+f[9]*alpha0[14]+alpha0[8]*f[13]+f[8]*alpha0[13]+alpha0[5]*f[11]+f[5]*alpha0[11]+alpha0[4]*f[10]+f[4]*alpha0[10]+alpha0[2]*f[7]+f[2]*alpha0[7]+alpha0[1]*f[6]+f[1]*alpha0[6]+alpha0[0]*f[3]+f[0]*alpha0[3]))/dxv[0]; 
  out[8] += (0.8660254037844386*(alpha0[11]*f[15]+f[11]*alpha0[15]+alpha0[7]*f[14]+f[7]*alpha0[14]+alpha0[6]*f[13]+f[6]*alpha0[13]+alpha0[5]*f[12]+f[5]*alpha0[12]+alpha0[3]*f[10]+f[3]*alpha0[10]+alpha0[2]*f[9]+f[2]*alpha0[9]+alpha0[1]*f[8]+f[1]*alpha0[8]+alpha0[0]*f[4]+f[0]*alpha0[4]))/dxv[0]; 
  out[11] += (0.8660254037844386*(alpha0[8]*f[15]+f[8]*alpha0[15]+alpha0[4]*f[14]+f[4]*alpha0[14]+alpha0[12]*f[13]+f[12]*alpha0[13]+alpha0[1]*f[11]+f[1]*alpha0[11]+alpha0[9]*f[10]+f[9]*alpha0[10]+alpha0[0]*f[7]+f[0]*alpha0[7]+alpha0[5]*f[6]+f[5]*alpha0[6]+alpha0[2]*f[3]+f[2]*alpha0[3]))/dxv[0]; 
  out[12] += (0.8660254037844386*(alpha0[6]*f[15]+f[6]*alpha0[15]+alpha0[3]*f[14]+f[3]*alpha0[14]+alpha0[11]*f[13]+f[11]*alpha0[13]+alpha0[1]*f[12]+f[1]*alpha0[12]+alpha0[7]*f[10]+f[7]*alpha0[10]+alpha0[0]*f[9]+f[0]*alpha0[9]+alpha0[5]*f[8]+f[5]*alpha0[8]+alpha0[2]*f[4]+f[2]*alpha0[4]))/dxv[0]; 
  out[13] += (0.8660254037844386*(alpha0[5]*f[15]+f[5]*alpha0[15]+alpha0[2]*f[14]+f[2]*alpha0[14]+alpha0[1]*f[13]+f[1]*alpha0[13]+alpha0[11]*f[12]+f[11]*alpha0[12]+alpha0[0]*f[10]+f[0]*alpha0[10]+alpha0[7]*f[9]+f[7]*alpha0[9]+alpha0[6]*f[8]+f[6]*alpha0[8]+alpha0[3]*f[4]+f[3]*alpha0[4]))/dxv[0]; 
  out[15] += (0.8660254037844386*(alpha0[1]*f[15]+f[1]*alpha0[15]+alpha0[0]*f[14]+f[0]*alpha0[14]+alpha0[5]*f[13]+f[5]*alpha0[13]+alpha0[6]*f[12]+f[6]*alpha0[12]+alpha0[8]*f[11]+f[8]*alpha0[11]+alpha0[2]*f[10]+f[2]*alpha0[10]+alpha0[3]*f[9]+f[3]*alpha0[9]+alpha0[4]*f[7]+f[4]*alpha0[7]))/dxv[0]; 

  double alphax[16]; 
  alphax[0] = (2.0*alpha0[0])/dxv[0]; 
  alphax[1] = (2.0*alpha0[1])/dxv[0]; 
  alphax[2] = (2.0*alpha0[2])/dxv[0]; 
  alphax[3] = (2.0*alpha0[3])/dxv[0]; 
  alphax[4] = (2.0*alpha0[4])/dxv[0]; 
  alphax[5] = (2.0*alpha0[5])/dxv[0]; 
  alphax[6] = (2.0*alpha0[6])/dxv[0]; 
  alphax[7] = (2.0*alpha0[7])/dxv[0]; 
  alphax[8] = (2.0*alpha0[8])/dxv[0]; 
  alphax[9] = (2.0*alpha0[9])/dxv[0]; 
  alphax[10] = (2.0*alpha0[10])/dxv[0]; 
  alphax[11] = (2.0*alpha0[11])/dxv[0]; 
  alphax[12] = (2.0*alpha0[12])/dxv[0]; 
  alphax[13] = (2.0*alpha0[13])/dxv[0]; 
  alphax[14] = (2.0*alpha0[14])/dxv[0]; 
  alphax[15] = (2.0*alpha0[15])/dxv[0]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.4330127018922193*alphax[15]-0.25*alphax[14]-0.4330127018922193*(alphax[13]+alphax[12]+alphax[11])+0.25*(alphax[10]+alphax[9])+0.4330127018922193*alphax[8]+0.25*alphax[7]+0.4330127018922193*(alphax[6]+alphax[5])-0.25*(alphax[4]+alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[15])+0.25*alphax[14]-0.4330127018922193*alphax[13]+0.4330127018922193*(alphax[12]+alphax[11])+0.25*alphax[10]-0.25*alphax[9]+0.4330127018922193*alphax[8]-0.25*alphax[7]+0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]-0.25*(alphax[4]+alphax[3])+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[15])+0.25*alphax[14]+0.4330127018922193*alphax[13]-0.4330127018922193*alphax[12]+0.4330127018922193*alphax[11]-0.25*alphax[10]+0.25*alphax[9]+0.4330127018922193*alphax[8]-0.25*alphax[7]-0.4330127018922193*alphax[6]+0.4330127018922193*alphax[5]-0.25*alphax[4]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[15]-0.25*alphax[14]+0.4330127018922193*(alphax[13]+alphax[12])-0.4330127018922193*alphax[11]-0.25*(alphax[10]+alphax[9])+0.4330127018922193*alphax[8]+0.25*alphax[7]-0.4330127018922193*(alphax[6]+alphax[5])-0.25*alphax[4]+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[15])+0.25*alphax[14]+0.4330127018922193*(alphax[13]+alphax[12])-0.4330127018922193*alphax[11]-0.25*(alphax[10]+alphax[9])-0.4330127018922193*alphax[8]+0.25*alphax[7]+0.4330127018922193*(alphax[6]+alphax[5])+0.25*alphax[4]-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[15]-0.25*alphax[14]+0.4330127018922193*alphax[13]-0.4330127018922193*alphax[12]+0.4330127018922193*alphax[11]-0.25*alphax[10]+0.25*alphax[9]-0.4330127018922193*alphax[8]-0.25*alphax[7]+0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]+0.25*alphax[4]-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[15]-0.25*alphax[14]-0.4330127018922193*alphax[13]+0.4330127018922193*(alphax[12]+alphax[11])+0.25*alphax[10]-0.25*alphax[9]-0.4330127018922193*alphax[8]-0.25*alphax[7]-0.4330127018922193*alphax[6]+0.4330127018922193*alphax[5]+0.25*(alphax[4]+alphax[3])-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[15])+0.25*alphax[14]-0.4330127018922193*(alphax[13]+alphax[12]+alphax[11])+0.25*(alphax[10]+alphax[9])-0.4330127018922193*alphax[8]+0.25*alphax[7]-0.4330127018922193*(alphax[6]+alphax[5])+0.25*(alphax[4]+alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.4330127018922193*alphax[15])-0.25*alphax[14]+0.4330127018922193*(alphax[13]+alphax[12]+alphax[11])+0.25*(alphax[10]+alphax[9])-0.4330127018922193*alphax[8]+0.25*alphax[7]-0.4330127018922193*(alphax[6]+alphax[5])-0.25*(alphax[4]+alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[15]+0.25*alphax[14]+0.4330127018922193*alphax[13]-0.4330127018922193*(alphax[12]+alphax[11])+0.25*alphax[10]-0.25*alphax[9]-0.4330127018922193*alphax[8]-0.25*alphax[7]-0.4330127018922193*alphax[6]+0.4330127018922193*alphax[5]-0.25*(alphax[4]+alphax[3])+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[15]+0.25*alphax[14]-0.4330127018922193*alphax[13]+0.4330127018922193*alphax[12]-0.4330127018922193*alphax[11]-0.25*alphax[10]+0.25*alphax[9]-0.4330127018922193*alphax[8]-0.25*alphax[7]+0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]-0.25*alphax[4]+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[15])-0.25*alphax[14]-0.4330127018922193*(alphax[13]+alphax[12])+0.4330127018922193*alphax[11]-0.25*(alphax[10]+alphax[9])-0.4330127018922193*alphax[8]+0.25*alphax[7]+0.4330127018922193*(alphax[6]+alphax[5])-0.25*alphax[4]+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[15]+0.25*alphax[14]-0.4330127018922193*(alphax[13]+alphax[12])+0.4330127018922193*alphax[11]-0.25*(alphax[10]+alphax[9])+0.4330127018922193*alphax[8]+0.25*alphax[7]-0.4330127018922193*(alphax[6]+alphax[5])+0.25*alphax[4]-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[15])-0.25*alphax[14]-0.4330127018922193*alphax[13]+0.4330127018922193*alphax[12]-0.4330127018922193*alphax[11]-0.25*alphax[10]+0.25*alphax[9]+0.4330127018922193*alphax[8]-0.25*alphax[7]-0.4330127018922193*alphax[6]+0.4330127018922193*alphax[5]+0.25*alphax[4]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[15])-0.25*alphax[14]+0.4330127018922193*alphax[13]-0.4330127018922193*(alphax[12]+alphax[11])+0.25*alphax[10]-0.25*alphax[9]+0.4330127018922193*alphax[8]-0.25*alphax[7]+0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]+0.25*(alphax[4]+alphax[3])-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[15]+0.25*alphax[14]+0.4330127018922193*(alphax[13]+alphax[12]+alphax[11])+0.25*(alphax[10]+alphax[9])+0.4330127018922193*alphax[8]+0.25*alphax[7]+0.4330127018922193*(alphax[6]+alphax[5])+0.25*(alphax[4]+alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  return cflFreq;
} 
