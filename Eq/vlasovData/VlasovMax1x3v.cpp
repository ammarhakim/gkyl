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
double VlasovVol1x3vMaxP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[4]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double *E2 = &EM[8]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[35]; 
  double alpha_vdim[105]; 

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
  alpha_vdim[21] = 0.8164965809277261*B2[2]*dv10*dv2; 
  alpha_vdim[25] = -0.8164965809277261*B1[2]*dv10*dv3; 
  alpha_vdim[31] = dv10*(2.828427124746191*(B2[3]*wv2+E0[3])-2.828427124746191*B1[3]*wv3); 
  alpha_mid += std::abs(0.125*alpha_vdim[0]-0.1397542485937369*alpha_vdim[11]); 
  alpha_vdim[35] = dv11*(2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]); 
  alpha_vdim[36] = dv11*(2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]); 
  alpha_vdim[37] = -0.8164965809277261*B2[0]*dv1*dv11; 
  alpha_vdim[39] = 0.8164965809277261*B0[0]*dv11*dv3; 
  alpha_vdim[40] = -0.8164965809277261*B2[1]*dv1*dv11; 
  alpha_vdim[43] = 0.8164965809277261*B0[1]*dv11*dv3; 
  alpha_vdim[46] = dv11*(2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]); 
  alpha_vdim[54] = -0.8164965809277261*B2[2]*dv1*dv11; 
  alpha_vdim[60] = 0.8164965809277261*B0[2]*dv11*dv3; 
  alpha_vdim[66] = dv11*(2.828427124746191*B0[3]*wv3-2.828427124746191*B2[3]*wv1+2.828427124746191*E1[3]); 
  alpha_mid += std::abs(0.125*alpha_vdim[35]-0.1397542485937369*alpha_vdim[46]); 
  alpha_vdim[70] = dv12*(2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2); 
  alpha_vdim[71] = dv12*(2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2); 
  alpha_vdim[72] = 0.8164965809277261*B1[0]*dv1*dv12; 
  alpha_vdim[73] = -0.8164965809277261*B0[0]*dv12*dv2; 
  alpha_vdim[75] = 0.8164965809277261*B1[1]*dv1*dv12; 
  alpha_vdim[76] = -0.8164965809277261*B0[1]*dv12*dv2; 
  alpha_vdim[81] = dv12*(2.828427124746191*(B1[2]*wv1+E2[2])-2.828427124746191*B0[2]*wv2); 
  alpha_vdim[89] = 0.8164965809277261*B1[2]*dv1*dv12; 
  alpha_vdim[91] = -0.8164965809277261*B0[2]*dv12*dv2; 
  alpha_vdim[101] = dv12*(2.828427124746191*(B1[3]*wv1+E2[3])-2.828427124746191*B0[3]*wv2); 
  alpha_mid += std::abs(0.125*alpha_vdim[70]-0.1397542485937369*alpha_vdim[81]); 
  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[31]*f[31]+alpha_vdim[25]*f[25]+alpha_vdim[21]*f[21]+alpha_vdim[11]*f[11]+alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[31]*alpha_vdim[66]+f[25]*alpha_vdim[60]+f[19]*alpha_vdim[54]+f[11]*alpha_vdim[46]+f[8]*alpha_vdim[43]+f[5]*alpha_vdim[40]+f[4]*alpha_vdim[39]+f[2]*alpha_vdim[37]+f[1]*alpha_vdim[36]+f[0]*alpha_vdim[35]); 
  out[4] += 0.4330127018922193*(f[31]*alpha_vdim[101]+f[21]*alpha_vdim[91]+f[19]*alpha_vdim[89]+f[11]*alpha_vdim[81]+f[6]*alpha_vdim[76]+f[5]*alpha_vdim[75]+f[3]*alpha_vdim[73]+f[2]*alpha_vdim[72]+f[1]*alpha_vdim[71]+f[0]*alpha_vdim[70]); 
  out[5] += 0.3803194146278324*(alpha_vdim[11]*f[31]+f[11]*alpha_vdim[31])+0.3872983346207416*(alpha_vdim[8]*f[25]+f[8]*alpha_vdim[25]+alpha_vdim[6]*f[21]+f[6]*alpha_vdim[21]+alpha_cdim[2]*f[12]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.4330127018922193*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[3]*f[6]+f[3]*alpha_vdim[6]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.3803194146278324*f[11]*alpha_vdim[66]+0.3872983346207416*(f[8]*alpha_vdim[60]+f[5]*alpha_vdim[54])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[46]+(0.3872983346207416*f[25]+0.4330127018922193*f[4])*alpha_vdim[43]+0.3872983346207416*f[19]*alpha_vdim[40]+0.4330127018922193*(f[2]*alpha_vdim[40]+f[8]*alpha_vdim[39]+f[5]*alpha_vdim[37])+0.3872983346207416*f[11]*alpha_vdim[36]+0.4330127018922193*(f[0]*alpha_vdim[36]+f[1]*alpha_vdim[35]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(f[11]*alpha_vdim[54]+f[19]*alpha_vdim[46]+f[16]*alpha_vdim[43])+0.3872983346207416*f[20]*alpha_vdim[40]+0.4330127018922193*(f[1]*alpha_vdim[40]+f[9]*alpha_vdim[39])+0.3872983346207416*f[12]*alpha_vdim[37]+0.4330127018922193*(f[0]*alpha_vdim[37]+f[5]*alpha_vdim[36]+f[2]*alpha_vdim[35])+0.3872983346207416*alpha_vdim[6]*f[23]+0.4330127018922193*(alpha_vdim[11]*f[21]+f[11]*alpha_vdim[21]+alpha_vdim[8]*f[17])+0.3872983346207416*alpha_vdim[3]*f[13]+0.4330127018922193*(alpha_vdim[4]*f[10]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[8] += 0.3803194146278324*f[11]*alpha_vdim[101]+0.3872983346207416*(f[6]*alpha_vdim[91]+f[5]*alpha_vdim[89])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[81]+(0.3872983346207416*f[21]+0.4330127018922193*f[3])*alpha_vdim[76]+0.3872983346207416*f[19]*alpha_vdim[75]+0.4330127018922193*(f[2]*alpha_vdim[75]+f[6]*alpha_vdim[73]+f[5]*alpha_vdim[72])+0.3872983346207416*f[11]*alpha_vdim[71]+0.4330127018922193*(f[0]*alpha_vdim[71]+f[1]*alpha_vdim[70]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[11]*alpha_vdim[89]+f[19]*alpha_vdim[81]+f[15]*alpha_vdim[76])+0.3872983346207416*f[20]*alpha_vdim[75]+0.4330127018922193*(f[1]*alpha_vdim[75]+f[7]*alpha_vdim[73])+0.3872983346207416*f[12]*alpha_vdim[72]+0.4330127018922193*(f[0]*alpha_vdim[72]+f[5]*alpha_vdim[71]+f[2]*alpha_vdim[70])+0.3872983346207416*alpha_vdim[8]*f[28]+0.4330127018922193*(alpha_vdim[11]*f[25]+f[11]*alpha_vdim[25]+alpha_vdim[6]*f[17])+0.3872983346207416*alpha_vdim[4]*f[14]+0.4330127018922193*(alpha_vdim[3]*f[10]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[10] += 0.4330127018922193*(f[11]*alpha_vdim[91]+f[21]*alpha_vdim[81])+0.3872983346207416*f[23]*alpha_vdim[76]+0.4330127018922193*(f[1]*alpha_vdim[76]+f[15]*alpha_vdim[75])+0.3872983346207416*f[13]*alpha_vdim[73]+0.4330127018922193*(f[0]*alpha_vdim[73]+f[7]*alpha_vdim[72]+f[6]*alpha_vdim[71]+f[3]*alpha_vdim[70]+f[11]*alpha_vdim[60]+f[25]*alpha_vdim[46])+0.3872983346207416*f[28]*alpha_vdim[43]+0.4330127018922193*(f[1]*alpha_vdim[43]+f[16]*alpha_vdim[40])+0.3872983346207416*f[14]*alpha_vdim[39]+0.4330127018922193*(f[0]*alpha_vdim[39]+f[9]*alpha_vdim[37]+f[8]*alpha_vdim[36]+f[4]*alpha_vdim[35]); 
  out[11] += 0.9682458365518543*(alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_vdim[11]*f[19]+alpha_vdim[8]*f[16]+alpha_vdim[6]*f[15]+alpha_vdim[4]*f[9]+alpha_vdim[3]*f[7]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[13] += 0.9682458365518543*(f[21]*alpha_vdim[46]+f[17]*alpha_vdim[43]+f[15]*alpha_vdim[40]+f[10]*alpha_vdim[39]+f[7]*alpha_vdim[37]+f[6]*alpha_vdim[36]+f[3]*alpha_vdim[35]); 
  out[14] += 0.9682458365518543*(f[25]*alpha_vdim[81]+f[17]*alpha_vdim[76]+f[16]*alpha_vdim[75]+f[10]*alpha_vdim[73]+f[9]*alpha_vdim[72]+f[8]*alpha_vdim[71]+f[4]*alpha_vdim[70]); 
  out[15] += 0.3803194146278324*f[19]*alpha_vdim[66]+0.3872983346207416*f[16]*alpha_vdim[60]+(0.3803194146278324*f[31]+0.3464101615137755*f[20])*alpha_vdim[54]+0.3872983346207416*(f[1]*alpha_vdim[54]+f[5]*alpha_vdim[46])+0.4330127018922193*f[9]*alpha_vdim[43]+0.3872983346207416*(f[12]+f[11])*alpha_vdim[40]+0.4330127018922193*(f[0]*alpha_vdim[40]+f[16]*alpha_vdim[39])+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[37]+0.3872983346207416*f[19]*alpha_vdim[36]+0.4330127018922193*(f[2]*alpha_vdim[36]+f[5]*alpha_vdim[35])+0.3803194146278324*(alpha_vdim[21]*f[31]+f[21]*alpha_vdim[31])+0.3872983346207416*f[17]*alpha_vdim[25]+0.3464101615137755*alpha_vdim[21]*f[23]+0.3872983346207416*(alpha_vdim[3]*f[23]+alpha_cdim[2]*f[22]+alpha_vdim[1]*f[21]+f[1]*alpha_vdim[21])+0.4330127018922193*alpha_vdim[4]*f[17]+0.3872983346207416*(alpha_vdim[6]*(f[13]+f[11])+f[6]*alpha_vdim[11])+0.4330127018922193*(alpha_vdim[8]*f[10]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+f[0]*alpha_vdim[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[16] += 0.3803194146278324*f[19]*alpha_vdim[101]+0.3872983346207416*f[15]*alpha_vdim[91]+(0.3803194146278324*f[31]+0.3464101615137755*f[20])*alpha_vdim[89]+0.3872983346207416*(f[1]*alpha_vdim[89]+f[5]*alpha_vdim[81])+0.4330127018922193*f[7]*alpha_vdim[76]+0.3872983346207416*(f[12]+f[11])*alpha_vdim[75]+0.4330127018922193*(f[0]*alpha_vdim[75]+f[15]*alpha_vdim[73])+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[72]+0.3872983346207416*f[19]*alpha_vdim[71]+0.4330127018922193*(f[2]*alpha_vdim[71]+f[5]*alpha_vdim[70])+0.3803194146278324*(alpha_vdim[25]*f[31]+f[25]*alpha_vdim[31])+0.3464101615137755*alpha_vdim[25]*f[28]+0.3872983346207416*(alpha_vdim[4]*f[28]+alpha_cdim[2]*f[26]+alpha_vdim[1]*f[25]+f[1]*alpha_vdim[25])+f[17]*(0.3872983346207416*alpha_vdim[21]+0.4330127018922193*alpha_vdim[3])+0.3872983346207416*(alpha_vdim[8]*(f[14]+f[11])+f[8]*alpha_vdim[11])+0.4330127018922193*(alpha_vdim[6]*f[10]+alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]); 
  out[17] += 0.3803194146278324*f[21]*alpha_vdim[101]+(0.3803194146278324*f[31]+0.3464101615137755*f[23])*alpha_vdim[91]+0.3872983346207416*(f[1]*alpha_vdim[91]+f[15]*alpha_vdim[89]+f[6]*alpha_vdim[81]+(f[13]+f[11])*alpha_vdim[76])+0.4330127018922193*(f[0]*alpha_vdim[76]+f[7]*alpha_vdim[75])+0.3872983346207416*f[23]*alpha_vdim[73]+0.4330127018922193*(f[1]*alpha_vdim[73]+f[15]*alpha_vdim[72])+0.3872983346207416*f[21]*alpha_vdim[71]+0.4330127018922193*(f[3]*alpha_vdim[71]+f[6]*alpha_vdim[70])+0.3803194146278324*f[25]*alpha_vdim[66]+(0.3803194146278324*f[31]+0.3464101615137755*f[28])*alpha_vdim[60]+0.3872983346207416*(f[1]*alpha_vdim[60]+f[16]*alpha_vdim[54]+f[8]*alpha_vdim[46]+(f[14]+f[11])*alpha_vdim[43])+0.4330127018922193*(f[0]*alpha_vdim[43]+f[9]*alpha_vdim[40])+0.3872983346207416*f[28]*alpha_vdim[39]+0.4330127018922193*(f[1]*alpha_vdim[39]+f[16]*alpha_vdim[37])+0.3872983346207416*f[25]*alpha_vdim[36]+0.4330127018922193*(f[4]*alpha_vdim[36]+f[8]*alpha_vdim[35]+alpha_cdim[2]*f[18]+alpha_cdim[0]*f[10]); 
  out[18] += 0.4330127018922193*(f[19]*alpha_vdim[91]+f[21]*alpha_vdim[89]+f[5]*alpha_vdim[76]+f[6]*alpha_vdim[75])+(0.3872983346207416*f[24]+0.4330127018922193*f[2])*alpha_vdim[73]+0.3872983346207416*f[22]*alpha_vdim[72]+0.4330127018922193*(f[3]*alpha_vdim[72]+f[15]*alpha_vdim[71]+f[7]*alpha_vdim[70]+f[19]*alpha_vdim[60]+f[25]*alpha_vdim[54]+f[5]*alpha_vdim[43]+f[8]*alpha_vdim[40])+(0.3872983346207416*f[29]+0.4330127018922193*f[2])*alpha_vdim[39]+0.3872983346207416*f[26]*alpha_vdim[37]+0.4330127018922193*(f[4]*alpha_vdim[37]+f[16]*alpha_vdim[36]+f[9]*alpha_vdim[35])+0.3872983346207416*(alpha_vdim[4]*f[30]+alpha_vdim[3]*f[27])+0.4330127018922193*(alpha_vdim[21]*f[25]+f[21]*alpha_vdim[25]+alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10]+alpha_vdim[6]*f[8]+f[6]*alpha_vdim[8]+alpha_vdim[3]*f[4]+f[3]*alpha_vdim[4]); 
  out[19] += 0.2581988897471612*alpha_vdim[31]*f[31]+0.3803194146278324*(alpha_vdim[1]*f[31]+f[1]*alpha_vdim[31])+0.276641667586244*alpha_vdim[25]*f[25]+0.4330127018922193*(alpha_vdim[4]*f[25]+f[4]*alpha_vdim[25])+0.276641667586244*alpha_vdim[21]*f[21]+0.4330127018922193*(alpha_vdim[3]*f[21]+f[3]*alpha_vdim[21])+0.8660254037844386*alpha_cdim[2]*f[20]+0.276641667586244*alpha_vdim[11]*f[11]+0.4330127018922193*(alpha_vdim[0]*f[11]+f[0]*alpha_vdim[11])+0.3872983346207416*(alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6])+0.9682458365518543*alpha_cdim[0]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1]); 
  out[20] += 0.3803194146278324*alpha_cdim[2]*f[32]+0.8504200642707612*f[19]*alpha_vdim[31]+0.8660254037844386*(f[16]*alpha_vdim[25]+f[15]*alpha_vdim[21]+alpha_vdim[1]*f[19])+0.9682458365518543*(alpha_vdim[4]*f[16]+alpha_vdim[3]*f[15])+0.4330127018922193*alpha_cdim[0]*f[12]+0.8660254037844386*f[5]*alpha_vdim[11]+0.9682458365518543*(alpha_vdim[8]*f[9]+alpha_vdim[6]*f[7]+alpha_vdim[0]*f[5])+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[2]; 
  out[21] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_vdim[66]+(0.276641667586244*f[25]+0.4330127018922193*f[4])*alpha_vdim[60]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_vdim[54]+(0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_vdim[46]+0.3872983346207416*(f[8]*alpha_vdim[43]+f[5]*alpha_vdim[40])+0.4330127018922193*(f[25]*alpha_vdim[39]+f[19]*alpha_vdim[37])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[36]+0.4330127018922193*f[11]*alpha_vdim[35]+0.9682458365518543*(alpha_cdim[2]*f[15]+alpha_cdim[0]*f[6]); 
  out[22] += 0.3872983346207416*(f[19]*alpha_vdim[54]+f[5]*alpha_vdim[40])+0.4330127018922193*f[26]*alpha_vdim[39]+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alpha_vdim[37]+0.4330127018922193*(f[20]*alpha_vdim[36]+f[12]*alpha_vdim[35])+0.8660254037844386*alpha_vdim[3]*f[24]+0.9682458365518543*(f[19]*alpha_vdim[21]+alpha_vdim[4]*f[18]+alpha_vdim[1]*f[15]+alpha_vdim[0]*f[7]+f[5]*alpha_vdim[6]+f[2]*alpha_vdim[3]); 
  out[23] += 0.8504200642707612*f[21]*alpha_vdim[66]+0.8660254037844386*(f[17]*alpha_vdim[60]+f[15]*alpha_vdim[54]+f[6]*alpha_vdim[46])+0.9682458365518543*(f[10]*alpha_vdim[43]+f[7]*alpha_vdim[40]+f[17]*alpha_vdim[39]+f[15]*alpha_vdim[37])+0.8660254037844386*f[21]*alpha_vdim[36]+0.9682458365518543*(f[3]*alpha_vdim[36]+f[6]*alpha_vdim[35])+0.4330127018922193*(alpha_cdim[2]*f[24]+alpha_cdim[0]*f[13]); 
  out[24] += 0.9682458365518543*(f[21]*alpha_vdim[54]+f[6]*alpha_vdim[40]+f[18]*alpha_vdim[39])+0.8660254037844386*f[22]*alpha_vdim[37]+0.9682458365518543*(f[3]*alpha_vdim[37]+f[15]*alpha_vdim[36]+f[7]*alpha_vdim[35])+0.3803194146278324*alpha_vdim[3]*f[33]+0.4330127018922193*(alpha_vdim[4]*f[27]+alpha_vdim[1]*f[23])+0.3872983346207416*alpha_vdim[21]*f[21]+0.4330127018922193*alpha_vdim[0]*f[13]+0.3872983346207416*(alpha_vdim[6]*f[6]+alpha_vdim[3]*f[3]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_vdim[101]+(0.276641667586244*f[21]+0.4330127018922193*f[3])*alpha_vdim[91]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_vdim[89]+(0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_vdim[81]+0.3872983346207416*(f[6]*alpha_vdim[76]+f[5]*alpha_vdim[75])+0.4330127018922193*(f[21]*alpha_vdim[73]+f[19]*alpha_vdim[72])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[71]+0.4330127018922193*f[11]*alpha_vdim[70]+0.9682458365518543*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[8]); 
  out[26] += 0.3872983346207416*(f[19]*alpha_vdim[89]+f[5]*alpha_vdim[75])+0.4330127018922193*f[22]*alpha_vdim[73]+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alpha_vdim[72]+0.4330127018922193*(f[20]*alpha_vdim[71]+f[12]*alpha_vdim[70])+0.8660254037844386*alpha_vdim[4]*f[29]+0.9682458365518543*(f[19]*alpha_vdim[25]+alpha_vdim[3]*f[18]+alpha_vdim[1]*f[16]+alpha_vdim[0]*f[9]+f[5]*alpha_vdim[8]+f[2]*alpha_vdim[4]); 
  out[27] += 0.3872983346207416*(f[21]*alpha_vdim[91]+f[6]*alpha_vdim[76])+(0.3803194146278324*f[33]+0.3872983346207416*f[3])*alpha_vdim[73]+0.4330127018922193*(f[24]*alpha_vdim[72]+f[23]*alpha_vdim[71]+f[13]*alpha_vdim[70])+0.9682458365518543*(f[21]*alpha_vdim[60]+f[6]*alpha_vdim[43])+0.8660254037844386*f[30]*alpha_vdim[39]+0.9682458365518543*(f[3]*alpha_vdim[39]+f[18]*alpha_vdim[37]+f[17]*alpha_vdim[36]+f[10]*alpha_vdim[35]); 
  out[28] += 0.8504200642707612*f[25]*alpha_vdim[101]+0.8660254037844386*(f[17]*alpha_vdim[91]+f[16]*alpha_vdim[89]+f[8]*alpha_vdim[81])+0.9682458365518543*(f[10]*alpha_vdim[76]+f[9]*alpha_vdim[75]+f[17]*alpha_vdim[73]+f[16]*alpha_vdim[72])+0.8660254037844386*f[25]*alpha_vdim[71]+0.9682458365518543*(f[4]*alpha_vdim[71]+f[8]*alpha_vdim[70])+0.4330127018922193*(alpha_cdim[2]*f[29]+alpha_cdim[0]*f[14]); 
  out[29] += 0.9682458365518543*(f[25]*alpha_vdim[89]+f[8]*alpha_vdim[75]+f[18]*alpha_vdim[73])+0.8660254037844386*f[26]*alpha_vdim[72]+0.9682458365518543*(f[4]*alpha_vdim[72]+f[16]*alpha_vdim[71]+f[9]*alpha_vdim[70])+0.3803194146278324*alpha_vdim[4]*f[34]+0.4330127018922193*(alpha_vdim[3]*f[30]+alpha_vdim[1]*f[28])+0.3872983346207416*alpha_vdim[25]*f[25]+0.4330127018922193*alpha_vdim[0]*f[14]+0.3872983346207416*(alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]); 
  out[30] += 0.9682458365518543*(f[25]*alpha_vdim[91]+f[8]*alpha_vdim[76])+0.8660254037844386*f[27]*alpha_vdim[73]+0.9682458365518543*(f[4]*alpha_vdim[73]+f[18]*alpha_vdim[72]+f[17]*alpha_vdim[71]+f[10]*alpha_vdim[70])+0.3872983346207416*(f[25]*alpha_vdim[60]+f[8]*alpha_vdim[43])+(0.3803194146278324*f[34]+0.3872983346207416*f[4])*alpha_vdim[39]+0.4330127018922193*(f[29]*alpha_vdim[37]+f[28]*alpha_vdim[36]+f[14]*alpha_vdim[35]); 
  out[31] += 1.479019945774904*(alpha_cdim[2]*f[19]+alpha_cdim[0]*f[11])+0.6614378277661477*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[32] += 0.6614378277661477*alpha_vdim[31]*f[31]+1.479019945774904*alpha_vdim[4]*f[26]+0.6614378277661477*alpha_vdim[25]*f[25]+1.479019945774904*alpha_vdim[3]*f[22]+0.6614378277661477*alpha_vdim[21]*f[21]+1.479019945774904*(alpha_vdim[1]*f[20]+alpha_vdim[0]*f[12])+0.6614378277661477*(alpha_vdim[11]*f[11]+alpha_vdim[8]*f[8]+alpha_vdim[6]*f[6]+alpha_vdim[4]*f[4]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[33] += 0.6614378277661477*(f[31]*alpha_vdim[66]+f[25]*alpha_vdim[60]+f[19]*alpha_vdim[54]+f[11]*alpha_vdim[46]+f[8]*alpha_vdim[43]+f[5]*alpha_vdim[40])+(1.479019945774904*f[27]+0.6614378277661477*f[4])*alpha_vdim[39]+(1.479019945774904*f[24]+0.6614378277661477*f[2])*alpha_vdim[37]+(1.479019945774904*f[23]+0.6614378277661477*f[1])*alpha_vdim[36]+(1.479019945774904*f[13]+0.6614378277661477*f[0])*alpha_vdim[35]; 
  out[34] += 0.6614378277661477*(f[31]*alpha_vdim[101]+f[21]*alpha_vdim[91]+f[19]*alpha_vdim[89]+f[11]*alpha_vdim[81]+f[6]*alpha_vdim[76]+f[5]*alpha_vdim[75])+(1.479019945774904*f[30]+0.6614378277661477*f[3])*alpha_vdim[73]+(1.479019945774904*f[29]+0.6614378277661477*f[2])*alpha_vdim[72]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_vdim[71]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alpha_vdim[70]; 
return alpha_mid; 

} 
