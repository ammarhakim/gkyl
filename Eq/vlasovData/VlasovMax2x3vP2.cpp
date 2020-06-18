#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVol2x3vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // EM:        Input EM-field.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[6]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv12 = 2/dxv[4]; 
  const double *E2 = &EM[12]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[42]; 
  double alpha_vdim[63]; 

  alpha_cdim[0] = 11.31370849898477*w0dx0; 
  alpha_cdim[3] = 3.265986323710906*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_cdim[21] = 11.31370849898477*w1dx1; 
  alpha_cdim[25] = 3.265986323710906*dv1dx1; 
  alpha_mid += std::abs(w1dx1)+0.5*dv1dx1; 

  alpha_vdim[0] = dv10*(2.828427124746191*(B2[0]*wv2+E0[0])-2.828427124746191*B1[0]*wv3); 
  alpha_vdim[1] = dv10*(2.828427124746191*(B2[1]*wv2+E0[1])-2.828427124746191*B1[1]*wv3); 
  alpha_vdim[2] = dv10*(2.828427124746191*(B2[2]*wv2+E0[2])-2.828427124746191*B1[2]*wv3); 
  alpha_vdim[4] = 0.8164965809277261*B2[0]*dv10*dv2; 
  alpha_vdim[5] = -0.8164965809277261*B1[0]*dv10*dv3; 
  alpha_vdim[6] = dv10*(2.828427124746191*(B2[3]*wv2+E0[3])-2.828427124746191*B1[3]*wv3); 
  alpha_vdim[9] = 0.8164965809277261*B2[1]*dv10*dv2; 
  alpha_vdim[10] = 0.8164965809277261*B2[2]*dv10*dv2; 
  alpha_vdim[12] = -0.8164965809277261*B1[1]*dv10*dv3; 
  alpha_vdim[13] = -0.8164965809277261*B1[2]*dv10*dv3; 
  alpha_vdim[16] = dv10*(2.828427124746191*(B2[4]*wv2+E0[4])-2.828427124746191*B1[4]*wv3); 
  alpha_vdim[17] = dv10*(2.828427124746191*(B2[5]*wv2+E0[5])-2.828427124746191*B1[5]*wv3); 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[0]-0.09882117688026182*(alpha_vdim[17]+alpha_vdim[16])); 

  alpha_vdim[21] = dv11*(2.828427124746191*B0[0]*wv3-2.828427124746191*B2[0]*wv1+2.828427124746191*E1[0]); 
  alpha_vdim[22] = dv11*(2.828427124746191*B0[1]*wv3-2.828427124746191*B2[1]*wv1+2.828427124746191*E1[1]); 
  alpha_vdim[23] = dv11*(2.828427124746191*B0[2]*wv3-2.828427124746191*B2[2]*wv1+2.828427124746191*E1[2]); 
  alpha_vdim[24] = -0.8164965809277261*B2[0]*dv1*dv11; 
  alpha_vdim[26] = 0.8164965809277261*B0[0]*dv11*dv3; 
  alpha_vdim[27] = dv11*(2.828427124746191*B0[3]*wv3-2.828427124746191*B2[3]*wv1+2.828427124746191*E1[3]); 
  alpha_vdim[28] = -0.8164965809277261*B2[1]*dv1*dv11; 
  alpha_vdim[29] = -0.8164965809277261*B2[2]*dv1*dv11; 
  alpha_vdim[33] = 0.8164965809277261*B0[1]*dv11*dv3; 
  alpha_vdim[34] = 0.8164965809277261*B0[2]*dv11*dv3; 
  alpha_vdim[37] = dv11*(2.828427124746191*B0[4]*wv3-2.828427124746191*B2[4]*wv1+2.828427124746191*E1[4]); 
  alpha_vdim[38] = dv11*(2.828427124746191*B0[5]*wv3-2.828427124746191*B2[5]*wv1+2.828427124746191*E1[5]); 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[21]-0.09882117688026182*(alpha_vdim[38]+alpha_vdim[37])); 

  alpha_vdim[42] = dv12*(2.828427124746191*(B1[0]*wv1+E2[0])-2.828427124746191*B0[0]*wv2); 
  alpha_vdim[43] = dv12*(2.828427124746191*(B1[1]*wv1+E2[1])-2.828427124746191*B0[1]*wv2); 
  alpha_vdim[44] = dv12*(2.828427124746191*(B1[2]*wv1+E2[2])-2.828427124746191*B0[2]*wv2); 
  alpha_vdim[45] = 0.8164965809277261*B1[0]*dv1*dv12; 
  alpha_vdim[46] = -0.8164965809277261*B0[0]*dv12*dv2; 
  alpha_vdim[48] = dv12*(2.828427124746191*(B1[3]*wv1+E2[3])-2.828427124746191*B0[3]*wv2); 
  alpha_vdim[49] = 0.8164965809277261*B1[1]*dv1*dv12; 
  alpha_vdim[50] = 0.8164965809277261*B1[2]*dv1*dv12; 
  alpha_vdim[51] = -0.8164965809277261*B0[1]*dv12*dv2; 
  alpha_vdim[52] = -0.8164965809277261*B0[2]*dv12*dv2; 
  alpha_vdim[58] = dv12*(2.828427124746191*(B1[4]*wv1+E2[4])-2.828427124746191*B0[4]*wv2); 
  alpha_vdim[59] = dv12*(2.828427124746191*(B1[5]*wv1+E2[5])-2.828427124746191*B0[5]*wv2); 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[42]-0.09882117688026182*(alpha_vdim[59]+alpha_vdim[58])); 

  out[1] += 0.3061862178478971*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.3061862178478971*(f[4]*alpha_cdim[25]+f[0]*alpha_cdim[21]); 
  out[3] += 0.3061862178478971*(alpha_vdim[17]*f[17]+alpha_vdim[16]*f[16]+alpha_vdim[13]*f[13]+alpha_vdim[12]*f[12]+alpha_vdim[10]*f[10]+alpha_vdim[9]*f[9]+alpha_vdim[6]*f[6]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[17]*alpha_vdim[38]+f[16]*alpha_vdim[37]+f[13]*alpha_vdim[34]+f[12]*alpha_vdim[33]+f[8]*alpha_vdim[29]+f[7]*alpha_vdim[28]+f[6]*alpha_vdim[27]+f[5]*alpha_vdim[26]+f[3]*alpha_vdim[24]+f[2]*alpha_vdim[23]+f[1]*alpha_vdim[22]+f[0]*alpha_vdim[21]); 
  out[5] += 0.3061862178478971*(f[17]*alpha_vdim[59]+f[16]*alpha_vdim[58]+f[10]*alpha_vdim[52]+f[9]*alpha_vdim[51]+f[8]*alpha_vdim[50]+f[7]*alpha_vdim[49]+f[6]*alpha_vdim[48]+f[4]*alpha_vdim[46]+f[3]*alpha_vdim[45]+f[2]*alpha_vdim[44]+f[1]*alpha_vdim[43]+f[0]*alpha_vdim[42]); 
  out[6] += 0.3061862178478971*(f[9]*alpha_cdim[25]+f[1]*alpha_cdim[21]+alpha_cdim[3]*f[8]+alpha_cdim[0]*f[2]); 
  out[7] += 0.273861278752583*(alpha_cdim[3]*f[18]+alpha_vdim[1]*f[16]+f[1]*alpha_vdim[16])+0.3061862178478971*(alpha_vdim[5]*f[12]+f[5]*alpha_vdim[12]+alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[2]*f[6]+f[2]*alpha_vdim[6]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[8] += 0.3061862178478971*(f[11]*alpha_cdim[25]+f[3]*alpha_cdim[21])+0.273861278752583*(alpha_vdim[2]*f[17]+f[2]*alpha_vdim[17])+0.3061862178478971*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13]+alpha_vdim[4]*f[10]+f[4]*alpha_vdim[10]+alpha_vdim[1]*f[6]+f[1]*alpha_vdim[6]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[9] += 0.273861278752583*f[1]*alpha_vdim[37]+0.3061862178478971*(f[5]*alpha_vdim[33]+f[3]*alpha_vdim[28]+f[2]*alpha_vdim[27]+f[12]*alpha_vdim[26]+f[7]*alpha_vdim[24]+f[6]*alpha_vdim[23])+0.273861278752583*f[16]*alpha_vdim[22]+0.3061862178478971*(f[0]*alpha_vdim[22]+f[1]*alpha_vdim[21]+alpha_cdim[3]*f[11]+alpha_cdim[0]*f[4]); 
  out[10] += 0.273861278752583*f[2]*alpha_vdim[38]+0.3061862178478971*(f[5]*alpha_vdim[34]+f[3]*alpha_vdim[29]+f[1]*alpha_vdim[27]+f[13]*alpha_vdim[26])+0.273861278752583*f[19]*alpha_cdim[25]+0.3061862178478971*(f[0]*alpha_cdim[25]+f[8]*alpha_vdim[24])+0.273861278752583*f[17]*alpha_vdim[23]+0.3061862178478971*(f[0]*alpha_vdim[23]+f[6]*alpha_vdim[22]+f[2]*alpha_vdim[21]+f[4]*alpha_cdim[21]); 
  out[11] += 0.3061862178478971*(f[2]*alpha_vdim[29]+f[1]*alpha_vdim[28]+f[14]*alpha_vdim[26])+0.273861278752583*f[18]*alpha_vdim[24]+0.3061862178478971*(f[0]*alpha_vdim[24]+f[8]*alpha_vdim[23]+f[7]*alpha_vdim[22]+f[3]*alpha_vdim[21])+0.273861278752583*alpha_vdim[4]*f[19]+0.3061862178478971*(alpha_vdim[5]*f[15]+alpha_vdim[2]*f[10]+f[2]*alpha_vdim[10]+alpha_vdim[1]*f[9]+f[1]*alpha_vdim[9]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[12] += 0.273861278752583*f[1]*alpha_vdim[58]+0.3061862178478971*(f[4]*alpha_vdim[51]+f[3]*alpha_vdim[49]+f[2]*alpha_vdim[48]+f[9]*alpha_vdim[46]+f[7]*alpha_vdim[45]+f[6]*alpha_vdim[44])+0.273861278752583*f[16]*alpha_vdim[43]+0.3061862178478971*(f[0]*alpha_vdim[43]+f[1]*alpha_vdim[42]+alpha_cdim[3]*f[14]+alpha_cdim[0]*f[5]); 
  out[13] += 0.273861278752583*f[2]*alpha_vdim[59]+0.3061862178478971*(f[4]*alpha_vdim[52]+f[3]*alpha_vdim[50]+f[1]*alpha_vdim[48]+f[10]*alpha_vdim[46]+f[8]*alpha_vdim[45])+0.273861278752583*f[17]*alpha_vdim[44]+0.3061862178478971*(f[0]*alpha_vdim[44]+f[6]*alpha_vdim[43]+f[2]*alpha_vdim[42]+f[15]*alpha_cdim[25]+f[5]*alpha_cdim[21]); 
  out[14] += 0.3061862178478971*(f[2]*alpha_vdim[50]+f[1]*alpha_vdim[49]+f[11]*alpha_vdim[46])+0.273861278752583*f[18]*alpha_vdim[45]+0.3061862178478971*(f[0]*alpha_vdim[45]+f[8]*alpha_vdim[44]+f[7]*alpha_vdim[43]+f[3]*alpha_vdim[42])+0.273861278752583*alpha_vdim[5]*f[20]+0.3061862178478971*(alpha_vdim[4]*f[15]+alpha_vdim[2]*f[13]+f[2]*alpha_vdim[13]+alpha_vdim[1]*f[12]+f[1]*alpha_vdim[12]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]); 
  out[15] += 0.3061862178478971*(f[2]*alpha_vdim[52]+f[1]*alpha_vdim[51])+0.273861278752583*f[19]*alpha_vdim[46]+0.3061862178478971*(f[0]*alpha_vdim[46]+f[11]*alpha_vdim[45]+f[10]*alpha_vdim[44]+f[9]*alpha_vdim[43]+f[4]*alpha_vdim[42]+f[2]*alpha_vdim[34]+f[1]*alpha_vdim[33])+0.273861278752583*f[20]*alpha_vdim[26]+0.3061862178478971*(f[0]*alpha_vdim[26]+f[14]*alpha_vdim[24]+f[13]*alpha_vdim[23]+f[12]*alpha_vdim[22]+f[5]*alpha_vdim[21]); 
  out[16] += 0.6846531968814573*(alpha_cdim[3]*f[7]+alpha_cdim[0]*f[1]); 
  out[17] += 0.6846531968814573*(f[10]*alpha_cdim[25]+f[2]*alpha_cdim[21]); 
  out[18] += 0.6846531968814573*(alpha_vdim[5]*f[14]+alpha_vdim[4]*f[11]+alpha_vdim[2]*f[8]+alpha_vdim[1]*f[7]+alpha_vdim[0]*f[3]); 
  out[19] += 0.6846531968814573*(f[15]*alpha_vdim[26]+f[11]*alpha_vdim[24]+f[10]*alpha_vdim[23]+f[9]*alpha_vdim[22]+f[4]*alpha_vdim[21]); 
  out[20] += 0.6846531968814573*(f[15]*alpha_vdim[46]+f[14]*alpha_vdim[45]+f[13]*alpha_vdim[44]+f[12]*alpha_vdim[43]+f[5]*alpha_vdim[42]); 

  return alpha_mid; 
} 
