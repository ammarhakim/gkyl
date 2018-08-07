#include <VlasovModDecl.h> 
double VlasovVolStream1x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[5]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[2] = 2.309401076758503*dv0dx0; 

  out[1] += 0.4330127018922193*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[15]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[2] = 2.309401076758503*dv0dx0; 

  out[1] += 0.4330127018922193*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[5] += 0.3872983346207416*alpha0[2]*f[12]+0.4330127018922193*(alpha0[0]*f[2]+f[0]*alpha0[2]); 
  out[6] += 0.4330127018922193*(alpha0[2]*f[7]+alpha0[0]*f[3]); 
  out[8] += 0.4330127018922193*(alpha0[2]*f[9]+alpha0[0]*f[4]); 
  out[11] += 0.9682458365518543*(alpha0[2]*f[5]+alpha0[0]*f[1]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x3vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[35]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[2] = 2.309401076758503*dv0dx0; 

  out[1] += 0.4330127018922193*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[5] += 0.3872983346207416*alpha0[2]*f[12]+0.4330127018922193*(alpha0[0]*f[2]+f[0]*alpha0[2]); 
  out[6] += 0.4330127018922193*(alpha0[2]*f[7]+alpha0[0]*f[3]); 
  out[8] += 0.4330127018922193*(alpha0[2]*f[9]+alpha0[0]*f[4]); 
  out[11] += 0.9682458365518543*(alpha0[2]*f[5]+alpha0[0]*f[1]); 
  out[15] += 0.3872983346207416*alpha0[2]*f[22]+0.4330127018922193*(alpha0[0]*f[7]+alpha0[2]*f[3]); 
  out[16] += 0.3872983346207416*alpha0[2]*f[26]+0.4330127018922193*(alpha0[0]*f[9]+alpha0[2]*f[4]); 
  out[17] += 0.4330127018922193*(alpha0[2]*f[18]+alpha0[0]*f[10]); 
  out[19] += 0.8660254037844386*alpha0[2]*f[20]+0.9682458365518543*(alpha0[0]*f[5]+f[1]*alpha0[2]); 
  out[20] += 0.3803194146278324*alpha0[2]*f[32]+0.4330127018922193*alpha0[0]*f[12]+0.3872983346207416*alpha0[2]*f[2]; 
  out[21] += 0.9682458365518543*(alpha0[2]*f[15]+alpha0[0]*f[6]); 
  out[23] += 0.4330127018922193*(alpha0[2]*f[24]+alpha0[0]*f[13]); 
  out[25] += 0.9682458365518543*(alpha0[2]*f[16]+alpha0[0]*f[8]); 
  out[28] += 0.4330127018922193*(alpha0[2]*f[29]+alpha0[0]*f[14]); 
  out[31] += 1.479019945774904*(alpha0[2]*f[19]+alpha0[0]*f[11])+0.6614378277661477*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
