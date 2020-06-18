#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVol3x3vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // EM:        Input EM-field.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  const double dv10 = 2/dxv[3]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv11 = 2/dxv[4]; 
  const double *E1 = &EM[4]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv12 = 2/dxv[5]; 
  const double *E2 = &EM[8]; 
  const double dv3 = dxv[5], wv3 = w[5]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[21]; 
  double alpha_vdim[21]; 

  alpha_cdim[0] = 16.0*w0dx0; 
  alpha_cdim[4] = 4.618802153517007*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_cdim[7] = 16.0*w1dx1; 
  alpha_cdim[12] = 4.618802153517007*dv1dx1; 
  alpha_mid += std::abs(w1dx1)+0.5*dv1dx1; 

  alpha_cdim[14] = 16.0*w2dx2; 
  alpha_cdim[20] = 4.618802153517007*dv2dx2; 
  alpha_mid += std::abs(w2dx2)+0.5*dv2dx2; 

  alpha_vdim[0] = dv10*(2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3); 
  alpha_vdim[1] = dv10*(2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3); 
  alpha_vdim[2] = dv10*(2.828427124746191*(B2[2]*wv2+E0[2])-2.828427124746191*B1[2]*wv3); 
  alpha_vdim[3] = dv10*(2.828427124746191*(B2[3]*wv2+E0[3])-2.828427124746191*B1[3]*wv3); 
  alpha_vdim[5] = 0.8164965809277261*B2[0]*dv10*dv2; 
  alpha_vdim[6] = -0.8164965809277261*B1[0]*dv10*dv3; 
  alpha_mid += std::abs(0.0625*alpha_vdim[0]); 

  alpha_vdim[7] = dv11*(2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]); 
  alpha_vdim[8] = dv11*(2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]); 
  alpha_vdim[9] = dv11*(2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]); 
  alpha_vdim[10] = dv11*(2.828427124746191*B0[3]*wv3-2.828427124746191*B2[3]*wv1+2.828427124746191*E1[3]); 
  alpha_vdim[11] = -0.8164965809277261*B2[0]*dv1*dv11; 
  alpha_vdim[13] = 0.8164965809277261*B0[0]*dv11*dv3; 
  alpha_mid += std::abs(0.0625*alpha_vdim[7]); 

  alpha_vdim[14] = dv12*(2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2); 
  alpha_vdim[15] = dv12*(2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2); 
  alpha_vdim[16] = dv12*(2.828427124746191*(B1[2]*wv1+E2[2])-2.828427124746191*B0[2]*wv2); 
  alpha_vdim[17] = dv12*(2.828427124746191*(B1[3]*wv1+E2[3])-2.828427124746191*B0[3]*wv2); 
  alpha_vdim[18] = 0.8164965809277261*B1[0]*dv1*dv12; 
  alpha_vdim[19] = -0.8164965809277261*B0[0]*dv12*dv2; 
  alpha_mid += std::abs(0.0625*alpha_vdim[14]); 

  out[1] += 0.2165063509461096*(alpha_cdim[4]*f[4]+alpha_cdim[0]*f[0]); 
  out[2] += 0.2165063509461096*(f[5]*alpha_cdim[12]+f[0]*alpha_cdim[7]); 
  out[3] += 0.2165063509461096*(f[6]*alpha_cdim[20]+f[0]*alpha_cdim[14]); 
  out[4] += 0.2165063509461096*(alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[6]*alpha_vdim[13]+f[4]*alpha_vdim[11]+f[3]*alpha_vdim[10]+f[2]*alpha_vdim[9]+f[1]*alpha_vdim[8]+f[0]*alpha_vdim[7]); 
  out[6] += 0.2165063509461096*(f[5]*alpha_vdim[19]+f[4]*alpha_vdim[18]+f[3]*alpha_vdim[17]+f[2]*alpha_vdim[16]+f[1]*alpha_vdim[15]+f[0]*alpha_vdim[14]); 

  return alpha_mid; 
} 
