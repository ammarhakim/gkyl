#include <VlasovModDecl.h> 
__host__ __device__ double VlasovNeutralVol1x3vSerP1(const double *w, const double *dxv, const double *boA, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // boA:       Input body acceleration.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *Fo0 = &boA[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *Fo1 = &boA[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double *Fo2 = &boA[4]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[16]; 
  double alpha_vdim[48]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.828427124746191*Fo0[0]*dv10; 
  alpha_vdim[1] = 2.828427124746191*Fo0[1]*dv10; 
  alpha_mid += std::abs(0.125*alpha_vdim[0]); 

  alpha_vdim[16] = 2.828427124746191*Fo1[0]*dv11; 
  alpha_vdim[17] = 2.828427124746191*Fo1[1]*dv11; 
  alpha_mid += std::abs(0.125*alpha_vdim[16]); 

  alpha_vdim[32] = 2.828427124746191*Fo2[0]*dv12; 
  alpha_vdim[33] = 2.828427124746191*Fo2[1]*dv12; 
  alpha_mid += std::abs(0.125*alpha_vdim[32]); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[1]*alpha_vdim[17]+f[0]*alpha_vdim[16]); 
  out[4] += 0.4330127018922193*(f[1]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[5] += 0.4330127018922193*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.4330127018922193*(f[0]*alpha_vdim[17]+f[1]*alpha_vdim[16]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(f[5]*alpha_vdim[17]+f[2]*alpha_vdim[16]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[8] += 0.4330127018922193*(f[0]*alpha_vdim[33]+f[1]*alpha_vdim[32]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[5]*alpha_vdim[33]+f[2]*alpha_vdim[32]+alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[6]*alpha_vdim[33]+f[3]*alpha_vdim[32]+f[8]*alpha_vdim[17]+f[4]*alpha_vdim[16]); 
  out[11] += 0.4330127018922193*(f[2]*alpha_vdim[17]+f[5]*alpha_vdim[16]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[12] += 0.4330127018922193*(f[2]*alpha_vdim[33]+f[5]*alpha_vdim[32]+alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]); 
  out[13] += 0.4330127018922193*(f[3]*alpha_vdim[33]+f[6]*alpha_vdim[32]+f[4]*alpha_vdim[17]+f[8]*alpha_vdim[16]+alpha_cdim[2]*f[14]+alpha_cdim[0]*f[10]); 
  out[14] += 0.4330127018922193*(f[11]*alpha_vdim[33]+f[7]*alpha_vdim[32]+f[12]*alpha_vdim[17]+f[9]*alpha_vdim[16]+alpha_vdim[1]*f[13]+alpha_vdim[0]*f[10]); 
  out[15] += 0.4330127018922193*(f[7]*alpha_vdim[33]+f[11]*alpha_vdim[32]+f[9]*alpha_vdim[17]+f[12]*alpha_vdim[16]+alpha_cdim[0]*f[14]+alpha_vdim[0]*f[13]+(alpha_cdim[2]+alpha_vdim[1])*f[10]); 

  return alpha_mid; 
} 
