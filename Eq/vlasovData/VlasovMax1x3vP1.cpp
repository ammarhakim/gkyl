#include <VlasovModDecl.h> 
double VlasovVol1x3vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double *E2 = &EM[4]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[5]; 
  double alpha_vdim[15]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_vdim[0] = dv10*(2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3); 
  alpha_vdim[1] = dv10*(2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3); 
  alpha_vdim[3] = 0.8164965809277261*B2[0]*dv10*dv2; 
  alpha_vdim[4] = -0.8164965809277261*B1[0]*dv10*dv3; 
  alpha_mid += std::abs(0.125*alpha_vdim[0]); 
  alpha_vdim[5] = dv11*(2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]); 
  alpha_vdim[6] = dv11*(2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]); 
  alpha_vdim[7] = -0.8164965809277261*B2[0]*dv1*dv11; 
  alpha_vdim[9] = 0.8164965809277261*B0[0]*dv11*dv3; 
  alpha_mid += std::abs(0.125*alpha_vdim[5]); 
  alpha_vdim[10] = dv12*(2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2); 
  alpha_vdim[11] = dv12*(2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2); 
  alpha_vdim[12] = 0.8164965809277261*B1[0]*dv1*dv12; 
  alpha_vdim[13] = -0.8164965809277261*B0[0]*dv12*dv2; 
  alpha_mid += std::abs(0.125*alpha_vdim[10]); 
  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[4]*alpha_vdim[9]+f[2]*alpha_vdim[7]+f[1]*alpha_vdim[6]+f[0]*alpha_vdim[5]); 
  out[4] += 0.4330127018922193*(f[3]*alpha_vdim[13]+f[2]*alpha_vdim[12]+f[1]*alpha_vdim[11]+f[0]*alpha_vdim[10]); 

  return alpha_mid; 
} 
