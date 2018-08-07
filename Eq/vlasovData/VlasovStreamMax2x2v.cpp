#include <VlasovModDecl.h> 
double VlasovVolStream2x2vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  double alpha0[5]; 

  double alpha1[5]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[3] = 2.309401076758503*dv0dx0; 

  alpha1[0] = 8.0*w1dx1; 
  alpha1[4] = 2.309401076758503*dv1dx1; 

  out[1] += 0.4330127018922193*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1); 
} 
double VlasovVolStream2x2vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  double alpha0[15]; 

  double alpha1[15]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[3] = 2.309401076758503*dv0dx0; 

  alpha1[0] = 8.0*w1dx1; 
  alpha1[4] = 2.309401076758503*dv1dx1; 

  out[1] += 0.4330127018922193*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha1[4]*f[8]+alpha0[3]*f[7]+alpha0[0]*f[2]+alpha1[0]*f[1]); 
  out[6] += 0.3872983346207416*alpha0[3]*f[13]+0.4330127018922193*(alpha0[0]*f[3]+f[0]*alpha0[3]); 
  out[7] += 0.4330127018922193*(alpha1[4]*f[10]+alpha1[0]*f[3]); 
  out[8] += 0.4330127018922193*(alpha0[3]*f[10]+alpha0[0]*f[4]); 
  out[9] += 0.3872983346207416*alpha1[4]*f[14]+0.4330127018922193*(alpha1[0]*f[4]+f[0]*alpha1[4]); 
  out[11] += 0.9682458365518543*(alpha0[3]*f[6]+alpha0[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha1[4]*f[9]+alpha1[0]*f[2]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1); 
} 
double VlasovVolStream2x2vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  double alpha0[35]; 

  double alpha1[35]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[3] = 2.309401076758503*dv0dx0; 

  alpha1[0] = 8.0*w1dx1; 
  alpha1[4] = 2.309401076758503*dv1dx1; 

  out[1] += 0.4330127018922193*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha1[4]*f[8]+alpha0[3]*f[7]+alpha0[0]*f[2]+alpha1[0]*f[1]); 
  out[6] += 0.3872983346207416*alpha0[3]*f[13]+0.4330127018922193*(alpha0[0]*f[3]+f[0]*alpha0[3]); 
  out[7] += 0.4330127018922193*(alpha1[4]*f[10]+alpha1[0]*f[3]); 
  out[8] += 0.4330127018922193*(alpha0[3]*f[10]+alpha0[0]*f[4]); 
  out[9] += 0.3872983346207416*alpha1[4]*f[14]+0.4330127018922193*(alpha1[0]*f[4]+f[0]*alpha1[4]); 
  out[11] += 0.9682458365518543*(alpha0[3]*f[6]+alpha0[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha1[4]*f[9]+alpha1[0]*f[2]); 
  out[15] += 0.3872983346207416*alpha0[3]*f[24]+0.4330127018922193*(alpha1[4]*f[17]+alpha0[0]*f[7]+alpha1[0]*f[6]+f[2]*alpha0[3]); 
  out[16] += 0.3872983346207416*alpha1[4]*f[28]+0.4330127018922193*(alpha0[3]*f[18]+alpha0[0]*f[9]+alpha1[0]*f[8]+f[1]*alpha1[4]); 
  out[17] += 0.3872983346207416*alpha0[3]*f[27]+0.4330127018922193*(alpha0[0]*f[10]+alpha0[3]*f[4]); 
  out[18] += 0.3872983346207416*alpha1[4]*f[30]+0.4330127018922193*(alpha1[0]*f[10]+f[3]*alpha1[4]); 
  out[19] += 0.4330127018922193*alpha1[4]*f[25]+0.9682458365518543*alpha0[3]*f[15]+0.4330127018922193*alpha1[0]*f[11]+0.9682458365518543*alpha0[0]*f[5]; 
  out[20] += 0.4330127018922193*alpha0[3]*f[22]+0.9682458365518543*alpha1[4]*f[16]+0.4330127018922193*alpha0[0]*f[12]+0.9682458365518543*alpha1[0]*f[5]; 
  out[21] += 0.8660254037844386*alpha0[3]*f[23]+0.9682458365518543*(alpha0[0]*f[6]+f[1]*alpha0[3]); 
  out[22] += 0.9682458365518543*(alpha1[4]*f[18]+alpha1[0]*f[7]); 
  out[23] += 0.3803194146278324*alpha0[3]*f[33]+0.4330127018922193*alpha0[0]*f[13]+0.3872983346207416*alpha0[3]*f[3]; 
  out[24] += 0.4330127018922193*(alpha1[4]*f[27]+alpha1[0]*f[13]); 
  out[25] += 0.9682458365518543*(alpha0[3]*f[17]+alpha0[0]*f[8]); 
  out[26] += 0.8660254037844386*alpha1[4]*f[29]+0.9682458365518543*(alpha1[0]*f[9]+f[2]*alpha1[4]); 
  out[28] += 0.4330127018922193*(alpha0[3]*f[30]+alpha0[0]*f[14]); 
  out[29] += 0.3803194146278324*alpha1[4]*f[34]+0.4330127018922193*alpha1[0]*f[14]+0.3872983346207416*alpha1[4]*f[4]; 
  out[31] += 1.479019945774904*(alpha0[3]*f[21]+alpha0[0]*f[11])+0.6614378277661477*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[32] += 1.479019945774904*(alpha1[4]*f[26]+alpha1[0]*f[12])+0.6614378277661477*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1); 
} 
