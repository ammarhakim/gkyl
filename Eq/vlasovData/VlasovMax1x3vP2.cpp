#include <VlasovModDecl.h> 
double VlasovVol1x3vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[3]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double *E2 = &EM[6]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[15]; 
  double alpha_vdim[45]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_vdim[0] = dv10*(2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3); 
  alpha_vdim[1] = dv10*(2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3); 
  alpha_vdim[3] = 0.8164965809277261*B2[0]*dv10*dv2; 
  alpha_vdim[4] = -0.8164965809277261*B1[0]*dv10*dv3; 
  alpha_vdim[6] = 0.8164965809277261*B2[1]*dv10*dv2; 
  alpha_vdim[8] = -0.8164965809277261*B1[1]*dv10*dv3; 
  alpha_vdim[11] = dv10*(2.828427124746191*(B2[2]*wv2+E0[2])-2.828427124746191*B1[2]*wv3); 
  alpha_mid += std::abs(0.125*alpha_vdim[0]-0.1397542485937369*alpha_vdim[11]); 
  alpha_vdim[15] = dv11*(2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]); 
  alpha_vdim[16] = dv11*(2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]); 
  alpha_vdim[17] = -0.8164965809277261*B2[0]*dv1*dv11; 
  alpha_vdim[19] = 0.8164965809277261*B0[0]*dv11*dv3; 
  alpha_vdim[20] = -0.8164965809277261*B2[1]*dv1*dv11; 
  alpha_vdim[23] = 0.8164965809277261*B0[1]*dv11*dv3; 
  alpha_vdim[26] = dv11*(2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]); 
  alpha_mid += std::abs(0.125*alpha_vdim[15]-0.1397542485937369*alpha_vdim[26]); 
  alpha_vdim[30] = dv12*(2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2); 
  alpha_vdim[31] = dv12*(2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2); 
  alpha_vdim[32] = 0.8164965809277261*B1[0]*dv1*dv12; 
  alpha_vdim[33] = -0.8164965809277261*B0[0]*dv12*dv2; 
  alpha_vdim[35] = 0.8164965809277261*B1[1]*dv1*dv12; 
  alpha_vdim[36] = -0.8164965809277261*B0[1]*dv12*dv2; 
  alpha_vdim[41] = dv12*(2.828427124746191*(B1[2]*wv1+E2[2])-2.828427124746191*B0[2]*wv2); 
  alpha_mid += std::abs(0.125*alpha_vdim[30]-0.1397542485937369*alpha_vdim[41]); 
  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[11]*f[11]+alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[11]*alpha_vdim[26]+f[8]*alpha_vdim[23]+f[5]*alpha_vdim[20]+f[4]*alpha_vdim[19]+f[2]*alpha_vdim[17]+f[1]*alpha_vdim[16]+f[0]*alpha_vdim[15]); 
  out[4] += 0.4330127018922193*(f[11]*alpha_vdim[41]+f[6]*alpha_vdim[36]+f[5]*alpha_vdim[35]+f[3]*alpha_vdim[33]+f[2]*alpha_vdim[32]+f[1]*alpha_vdim[31]+f[0]*alpha_vdim[30]); 
  out[5] += 0.3872983346207416*(alpha_cdim[2]*f[12]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.4330127018922193*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.3872983346207416*f[1]*alpha_vdim[26]+0.4330127018922193*(f[4]*alpha_vdim[23]+f[2]*alpha_vdim[20]+f[8]*alpha_vdim[19]+f[5]*alpha_vdim[17])+0.3872983346207416*f[11]*alpha_vdim[16]+0.4330127018922193*(f[0]*alpha_vdim[16]+f[1]*alpha_vdim[15]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(f[1]*alpha_vdim[20]+f[9]*alpha_vdim[19])+0.3872983346207416*f[12]*alpha_vdim[17]+0.4330127018922193*(f[0]*alpha_vdim[17]+f[5]*alpha_vdim[16]+f[2]*alpha_vdim[15])+0.3872983346207416*alpha_vdim[3]*f[13]+0.4330127018922193*(alpha_vdim[4]*f[10]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[8] += 0.3872983346207416*f[1]*alpha_vdim[41]+0.4330127018922193*(f[3]*alpha_vdim[36]+f[2]*alpha_vdim[35]+f[6]*alpha_vdim[33]+f[5]*alpha_vdim[32])+0.3872983346207416*f[11]*alpha_vdim[31]+0.4330127018922193*(f[0]*alpha_vdim[31]+f[1]*alpha_vdim[30]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[1]*alpha_vdim[35]+f[7]*alpha_vdim[33])+0.3872983346207416*f[12]*alpha_vdim[32]+0.4330127018922193*(f[0]*alpha_vdim[32]+f[5]*alpha_vdim[31]+f[2]*alpha_vdim[30])+0.3872983346207416*alpha_vdim[4]*f[14]+0.4330127018922193*(alpha_vdim[3]*f[10]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[10] += 0.4330127018922193*f[1]*alpha_vdim[36]+0.3872983346207416*f[13]*alpha_vdim[33]+0.4330127018922193*(f[0]*alpha_vdim[33]+f[7]*alpha_vdim[32]+f[6]*alpha_vdim[31]+f[3]*alpha_vdim[30]+f[1]*alpha_vdim[23])+0.3872983346207416*f[14]*alpha_vdim[19]+0.4330127018922193*(f[0]*alpha_vdim[19]+f[9]*alpha_vdim[17]+f[8]*alpha_vdim[16]+f[4]*alpha_vdim[15]); 
  out[11] += 0.9682458365518543*(alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_vdim[4]*f[9]+alpha_vdim[3]*f[7]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[13] += 0.9682458365518543*(f[10]*alpha_vdim[19]+f[7]*alpha_vdim[17]+f[6]*alpha_vdim[16]+f[3]*alpha_vdim[15]); 
  out[14] += 0.9682458365518543*(f[10]*alpha_vdim[33]+f[9]*alpha_vdim[32]+f[8]*alpha_vdim[31]+f[4]*alpha_vdim[30]); 

  return alpha_mid; 
} 
