#include <VlasovModDecl.h> 
double VlasovVol2x3vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[3]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv12 = 2/dxv[4]; 
  const double *E2 = &EM[6]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[12]; 
  double alpha_vdim[18]; 

  alpha_cdim[0] = 11.31370849898477*w0dx0; 
  alpha_cdim[3] = 3.265986323710906*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_cdim[6] = 11.31370849898477*w1dx1; 
  alpha_cdim[10] = 3.265986323710906*dv1dx1; 
  alpha_mid += std::abs(w1dx1)+0.5*dv1dx1; 
  alpha_vdim[0] = dv10*(2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3); 
  alpha_vdim[1] = dv10*(2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3); 
  alpha_vdim[2] = dv10*(2.828427124746191*(B2[2]*wv2+E0[2])-2.828427124746191*B1[2]*wv3); 
  alpha_vdim[4] = 0.8164965809277261*B2[0]*dv10*dv2; 
  alpha_vdim[5] = -0.8164965809277261*B1[0]*dv10*dv3; 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[0]); 
  alpha_vdim[6] = dv11*(2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]); 
  alpha_vdim[7] = dv11*(2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]); 
  alpha_vdim[8] = dv11*(2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]); 
  alpha_vdim[9] = -0.8164965809277261*B2[0]*dv1*dv11; 
  alpha_vdim[11] = 0.8164965809277261*B0[0]*dv11*dv3; 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[6]); 
  alpha_vdim[12] = dv12*(2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2); 
  alpha_vdim[13] = dv12*(2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2); 
  alpha_vdim[14] = dv12*(2.828427124746191*(B1[2]*wv1+E2[2])-2.828427124746191*B0[2]*wv2); 
  alpha_vdim[15] = 0.8164965809277261*B1[0]*dv1*dv12; 
  alpha_vdim[16] = -0.8164965809277261*B0[0]*dv12*dv2; 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[12]); 
  out[1] += 0.3061862178478971*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.3061862178478971*(f[4]*alpha_cdim[10]+f[0]*alpha_cdim[6]); 
  out[3] += 0.3061862178478971*(alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[5]*alpha_vdim[11]+f[3]*alpha_vdim[9]+f[2]*alpha_vdim[8]+f[1]*alpha_vdim[7]+f[0]*alpha_vdim[6]); 
  out[5] += 0.3061862178478971*(f[4]*alpha_vdim[16]+f[3]*alpha_vdim[15]+f[2]*alpha_vdim[14]+f[1]*alpha_vdim[13]+f[0]*alpha_vdim[12]); 

  return alpha_mid; 
} 
