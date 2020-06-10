//#include <VlasovModDecl.h> 
//__host__ __device__ double VlasovVol2x3vSerP1_v1(const int nC, const int nC_EM, const double *w, const double *dxv, const double *EM, const double *f, double *out) 
__host__ __device__ double VlasovVol2x3vSerP1_v1(const int nC, const int nC_EM,  const double* __restrict__ w,
                                                 const double* __restrict__ dxv,
                                                 const double* __restrict__ EM, const double* __restrict__ f, double *out)
    
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[nC_EM*0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[nC_EM*4]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv12 = 2/dxv[4]; 
  const double *E2 = &EM[nC_EM*8]; 
  const double dv3 = dxv[4], wv3 = w[4]; 

  const double *B0 = &EM[nC_EM*12]; 
  const double *B1 = &EM[nC_EM*16]; 
  const double *B2 = &EM[nC_EM*20]; 
  double alpha_mid = 0.0; 
  double alpha_cdim[64]; 
  double alpha_vdim[96]; 

  alpha_cdim[0] = 11.31370849898477*w0dx0; 
  alpha_cdim[3] = 3.265986323710906*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_cdim[32] = 11.31370849898477*w1dx1; 
  alpha_cdim[36] = 3.265986323710906*dv1dx1; 
  alpha_mid += std::abs(w1dx1)+0.5*dv1dx1; 

  alpha_vdim[0] = dv10*(2.828427124746191*(B2[nC_EM*0]*wv2+E0[nC_EM*0])-2.828427124746191*B1[nC_EM*0]*wv3); 
  alpha_vdim[1] = dv10*(2.828427124746191*(B2[nC_EM*1]*wv2+E0[nC_EM*1])-2.828427124746191*B1[nC_EM*1]*wv3); 
  alpha_vdim[2] = dv10*(2.828427124746191*(B2[nC_EM*2]*wv2+E0[nC_EM*2])-2.828427124746191*B1[nC_EM*2]*wv3); 
  alpha_vdim[4] = 0.8164965809277261*B2[nC_EM*0]*dv10*dv2; 
  alpha_vdim[5] = -0.8164965809277261*B1[nC_EM*0]*dv10*dv3; 
  alpha_vdim[6] = dv10*(2.828427124746191*(B2[nC_EM*3]*wv2+E0[nC_EM*3])-2.828427124746191*B1[nC_EM*3]*wv3); 
  alpha_vdim[9] = 0.8164965809277261*B2[nC_EM*1]*dv10*dv2; 
  alpha_vdim[10] = 0.8164965809277261*B2[nC_EM*2]*dv10*dv2; 
  alpha_vdim[12] = -0.8164965809277261*B1[nC_EM*1]*dv10*dv3; 
  alpha_vdim[13] = -0.8164965809277261*B1[nC_EM*2]*dv10*dv3; 
  alpha_vdim[17] = 0.8164965809277261*B2[nC_EM*3]*dv10*dv2; 
  alpha_vdim[20] = -0.8164965809277261*B1[nC_EM*3]*dv10*dv3; 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[0]); 

  alpha_vdim[32] = dv11*(2.828427124746191*B0[nC_EM*0]*wv3-2.828427124746191*B2[nC_EM*0]*wv1+2.828427124746191*E1[nC_EM*0]); 
  alpha_vdim[33] = dv11*(2.828427124746191*B0[nC_EM*1]*wv3-2.828427124746191*B2[nC_EM*1]*wv1+2.828427124746191*E1[nC_EM*1]); 
  alpha_vdim[34] = dv11*(2.828427124746191*B0[nC_EM*2]*wv3-2.828427124746191*B2[nC_EM*2]*wv1+2.828427124746191*E1[nC_EM*2]); 
  alpha_vdim[35] = -0.8164965809277261*B2[nC_EM*0]*dv1*dv11; 
  alpha_vdim[37] = 0.8164965809277261*B0[nC_EM*0]*dv11*dv3; 
  alpha_vdim[38] = dv11*(2.828427124746191*B0[nC_EM*3]*wv3-2.828427124746191*B2[nC_EM*3]*wv1+2.828427124746191*E1[nC_EM*3]); 
  alpha_vdim[39] = -0.8164965809277261*B2[nC_EM*1]*dv1*dv11; 
  alpha_vdim[40] = -0.8164965809277261*B2[nC_EM*2]*dv1*dv11; 
  alpha_vdim[44] = 0.8164965809277261*B0[nC_EM*1]*dv11*dv3; 
  alpha_vdim[45] = 0.8164965809277261*B0[nC_EM*2]*dv11*dv3; 
  alpha_vdim[48] = -0.8164965809277261*B2[nC_EM*3]*dv1*dv11; 
  alpha_vdim[52] = 0.8164965809277261*B0[nC_EM*3]*dv11*dv3; 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[32]); 

  alpha_vdim[64] = dv12*(2.828427124746191*(B1[nC_EM*0]*wv1+E2[nC_EM*0])-2.828427124746191*B0[nC_EM*0]*wv2); 
  alpha_vdim[65] = dv12*(2.828427124746191*(B1[nC_EM*1]*wv1+E2[nC_EM*1])-2.828427124746191*B0[nC_EM*1]*wv2); 
  alpha_vdim[66] = dv12*(2.828427124746191*(B1[nC_EM*2]*wv1+E2[nC_EM*2])-2.828427124746191*B0[nC_EM*2]*wv2); 
  alpha_vdim[67] = 0.8164965809277261*B1[nC_EM*0]*dv1*dv12; 
  alpha_vdim[68] = -0.8164965809277261*B0[nC_EM*0]*dv12*dv2; 
  alpha_vdim[70] = dv12*(2.828427124746191*(B1[nC_EM*3]*wv1+E2[nC_EM*3])-2.828427124746191*B0[nC_EM*3]*wv2); 
  alpha_vdim[71] = 0.8164965809277261*B1[nC_EM*1]*dv1*dv12; 
  alpha_vdim[72] = 0.8164965809277261*B1[nC_EM*2]*dv1*dv12; 
  alpha_vdim[73] = -0.8164965809277261*B0[nC_EM*1]*dv12*dv2; 
  alpha_vdim[74] = -0.8164965809277261*B0[nC_EM*2]*dv12*dv2; 
  alpha_vdim[80] = 0.8164965809277261*B1[nC_EM*3]*dv1*dv12; 
  alpha_vdim[81] = -0.8164965809277261*B0[nC_EM*3]*dv12*dv2; 
  alpha_mid += std::abs(0.0883883476483184*alpha_vdim[64]); 

  // nC is number of Cells, i.e., the distance of buffer between out[0] and out[1], or f[0] and f[1]
  
  out[1*nC] += 0.3061862178478971*(alpha_cdim[3]*f[3*nC]+alpha_cdim[0]*f[0]); 
  out[2*nC] += 0.3061862178478971*(f[4*nC]*alpha_cdim[36]+f[0]*alpha_cdim[32]); 
  out[3*nC] += 0.3061862178478971*(alpha_vdim[20]*f[20*nC]+alpha_vdim[17]*f[17*nC]+alpha_vdim[13]*f[13*nC]+alpha_vdim[12]*f[12*nC]+alpha_vdim[10]*f[10*nC]+alpha_vdim[9]*f[9*nC]+alpha_vdim[6]*f[6*nC]+alpha_vdim[5]*f[5*nC]+alpha_vdim[4]*f[4*nC]+alpha_vdim[2]*f[2*nC]+alpha_vdim[1]*f[1*nC]+alpha_vdim[0]*f[0]); 
  out[4*nC] += 0.3061862178478971*(f[20*nC]*alpha_vdim[52]+f[16*nC]*alpha_vdim[48]+f[13*nC]*alpha_vdim[45]+f[12*nC]*alpha_vdim[44]+f[8*nC]*alpha_vdim[40]+f[7*nC]*alpha_vdim[39]+f[6*nC]*alpha_vdim[38]+f[5*nC]*alpha_vdim[37]+f[3*nC]*alpha_vdim[35]+f[2*nC]*alpha_vdim[34]+f[1*nC]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[5*nC] += 0.3061862178478971*(f[17*nC]*alpha_vdim[81]+f[16*nC]*alpha_vdim[80]+f[10*nC]*alpha_vdim[74]+f[9*nC]*alpha_vdim[73]+f[8*nC]*alpha_vdim[72]+f[7*nC]*alpha_vdim[71]+f[6*nC]*alpha_vdim[70]+f[4*nC]*alpha_vdim[68]+f[3*nC]*alpha_vdim[67]+f[2*nC]*alpha_vdim[66]+f[1*nC]*alpha_vdim[65]+f[0]*alpha_vdim[64]); 
  out[6*nC] += 0.3061862178478971*(f[9*nC]*alpha_cdim[36]+f[1*nC]*alpha_cdim[32]+alpha_cdim[3]*f[8*nC]+alpha_cdim[0]*f[2*nC]); 
  out[7*nC] += 0.3061862178478971*(alpha_vdim[13]*f[20*nC]+f[13*nC]*alpha_vdim[20]+alpha_vdim[10]*f[17*nC]+f[10*nC]*alpha_vdim[17]+alpha_vdim[5]*f[12*nC]+f[5*nC]*alpha_vdim[12]+alpha_vdim[4]*f[9*nC]+f[4*nC]*alpha_vdim[9]+alpha_vdim[2]*f[6*nC]+f[2*nC]*alpha_vdim[6]+alpha_cdim[0]*f[3*nC]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1*nC]+f[0]*alpha_vdim[1]); 
  out[8*nC] += 0.3061862178478971*(f[11*nC]*alpha_cdim[36]+f[3*nC]*alpha_cdim[32]+alpha_vdim[12]*f[20*nC]+f[12*nC]*alpha_vdim[20]+alpha_vdim[9]*f[17*nC]+f[9*nC]*alpha_vdim[17]+alpha_vdim[5]*f[13*nC]+f[5*nC]*alpha_vdim[13]+alpha_vdim[4]*f[10*nC]+f[4*nC]*alpha_vdim[10]+alpha_vdim[1]*f[6*nC]+f[1*nC]*alpha_vdim[6]+alpha_vdim[0]*f[2*nC]+f[0]*alpha_vdim[2]); 
  out[9*nC] += 0.3061862178478971*(f[13*nC]*alpha_vdim[52]+f[8*nC]*alpha_vdim[48]+f[20*nC]*alpha_vdim[45]+f[5*nC]*alpha_vdim[44]+f[16*nC]*alpha_vdim[40]+f[3*nC]*alpha_vdim[39]+f[2*nC]*alpha_vdim[38]+f[12*nC]*alpha_vdim[37]+f[7*nC]*alpha_vdim[35]+f[6*nC]*alpha_vdim[34]+f[0]*alpha_vdim[33]+f[1*nC]*alpha_vdim[32]+alpha_cdim[3]*f[11*nC]+alpha_cdim[0]*f[4*nC]); 
  out[10*nC] += 0.3061862178478971*(f[12*nC]*alpha_vdim[52]+f[7*nC]*alpha_vdim[48]+f[5*nC]*alpha_vdim[45]+f[20*nC]*alpha_vdim[44]+f[3*nC]*alpha_vdim[40]+f[16*nC]*alpha_vdim[39]+f[1*nC]*alpha_vdim[38]+f[13*nC]*alpha_vdim[37]+f[0]*alpha_cdim[36]+f[8*nC]*alpha_vdim[35]+f[0]*alpha_vdim[34]+f[6*nC]*alpha_vdim[33]+f[2*nC]*alpha_vdim[32]+f[4*nC]*alpha_cdim[32]); 
  out[11*nC] += 0.3061862178478971*(f[27*nC]*alpha_vdim[52]+f[6*nC]*alpha_vdim[48]+f[22*nC]*alpha_vdim[45]+f[21*nC]*alpha_vdim[44]+f[2*nC]*alpha_vdim[40]+f[1*nC]*alpha_vdim[39]+f[16*nC]*alpha_vdim[38]+f[14*nC]*alpha_vdim[37]+f[0]*alpha_vdim[35]+f[8*nC]*alpha_vdim[34]+f[7*nC]*alpha_vdim[33]+f[3*nC]*alpha_vdim[32]+alpha_vdim[20]*f[28*nC]+alpha_vdim[13]*f[24*nC]+alpha_vdim[12]*f[23*nC]+alpha_vdim[6]*f[17*nC]+f[6*nC]*alpha_vdim[17]+alpha_vdim[5]*f[15*nC]+alpha_vdim[2]*f[10*nC]+f[2*nC]*alpha_vdim[10]+alpha_vdim[1]*f[9*nC]+f[1*nC]*alpha_vdim[9]+alpha_vdim[0]*f[4*nC]+f[0]*alpha_vdim[4]); 
  out[12*nC] += 0.3061862178478971*(f[10*nC]*alpha_vdim[81]+f[8*nC]*alpha_vdim[80]+f[17*nC]*alpha_vdim[74]+f[4*nC]*alpha_vdim[73]+f[16*nC]*alpha_vdim[72]+f[3*nC]*alpha_vdim[71]+f[2*nC]*alpha_vdim[70]+f[9*nC]*alpha_vdim[68]+f[7*nC]*alpha_vdim[67]+f[6*nC]*alpha_vdim[66]+f[0]*alpha_vdim[65]+f[1*nC]*alpha_vdim[64]+alpha_cdim[3]*f[14*nC]+alpha_cdim[0]*f[5*nC]); 
  out[13*nC] += 0.3061862178478971*(f[9*nC]*alpha_vdim[81]+f[7*nC]*alpha_vdim[80]+f[4*nC]*alpha_vdim[74]+f[17*nC]*alpha_vdim[73]+f[3*nC]*alpha_vdim[72]+f[16*nC]*alpha_vdim[71]+f[1*nC]*alpha_vdim[70]+f[10*nC]*alpha_vdim[68]+f[8*nC]*alpha_vdim[67]+f[0]*alpha_vdim[66]+f[6*nC]*alpha_vdim[65]+f[2*nC]*alpha_vdim[64]+f[15*nC]*alpha_cdim[36]+f[5*nC]*alpha_cdim[32]); 
  out[14*nC] += 0.3061862178478971*(f[26*nC]*alpha_vdim[81]+f[6*nC]*alpha_vdim[80]+f[19*nC]*alpha_vdim[74]+f[18*nC]*alpha_vdim[73]+f[2*nC]*alpha_vdim[72]+f[1*nC]*alpha_vdim[71]+f[16*nC]*alpha_vdim[70]+f[11*nC]*alpha_vdim[68]+f[0]*alpha_vdim[67]+f[8*nC]*alpha_vdim[66]+f[7*nC]*alpha_vdim[65]+f[3*nC]*alpha_vdim[64]+alpha_vdim[17]*f[28*nC]+alpha_vdim[10]*f[24*nC]+alpha_vdim[9]*f[23*nC]+alpha_vdim[6]*f[20*nC]+f[6*nC]*alpha_vdim[20]+alpha_vdim[4]*f[15*nC]+alpha_vdim[2]*f[13*nC]+f[2*nC]*alpha_vdim[13]+alpha_vdim[1]*f[12*nC]+f[1*nC]*alpha_vdim[12]+alpha_vdim[0]*f[5*nC]+f[0]*alpha_vdim[5]); 
  out[15*nC] += 0.3061862178478971*(f[6*nC]*alpha_vdim[81]+f[26*nC]*alpha_vdim[80]+f[2*nC]*alpha_vdim[74]+f[1*nC]*alpha_vdim[73]+f[19*nC]*alpha_vdim[72]+f[18*nC]*alpha_vdim[71]+f[17*nC]*alpha_vdim[70]+f[0]*alpha_vdim[68]+f[11*nC]*alpha_vdim[67]+f[10*nC]*alpha_vdim[66]+f[9*nC]*alpha_vdim[65]+f[4*nC]*alpha_vdim[64]+f[6*nC]*alpha_vdim[52]+f[27*nC]*alpha_vdim[48]+f[2*nC]*alpha_vdim[45]+f[1*nC]*alpha_vdim[44]+f[22*nC]*alpha_vdim[40]+f[21*nC]*alpha_vdim[39]+f[20*nC]*alpha_vdim[38]+f[0]*alpha_vdim[37]+f[14*nC]*alpha_vdim[35]+f[13*nC]*alpha_vdim[34]+f[12*nC]*alpha_vdim[33]+f[5*nC]*alpha_vdim[32]); 
  out[16*nC] += 0.3061862178478971*(f[18*nC]*alpha_cdim[36]+f[7*nC]*alpha_cdim[32]+alpha_vdim[5]*f[20*nC]+f[5*nC]*alpha_vdim[20]+alpha_vdim[4]*f[17*nC]+f[4*nC]*alpha_vdim[17]+alpha_vdim[12]*f[13*nC]+f[12*nC]*alpha_vdim[13]+alpha_vdim[9]*f[10*nC]+f[9*nC]*alpha_vdim[10]+alpha_cdim[0]*f[8*nC]+alpha_vdim[0]*f[6*nC]+f[0]*alpha_vdim[6]+f[2*nC]*(alpha_cdim[3]+alpha_vdim[1])+f[1*nC]*alpha_vdim[2]); 
  out[17*nC] += 0.3061862178478971*(f[5*nC]*alpha_vdim[52]+f[3*nC]*alpha_vdim[48]+f[12*nC]*alpha_vdim[45]+f[13*nC]*alpha_vdim[44]+f[7*nC]*alpha_vdim[40]+f[8*nC]*alpha_vdim[39]+f[0]*alpha_vdim[38]+f[20*nC]*alpha_vdim[37]+f[1*nC]*alpha_cdim[36]+f[16*nC]*alpha_vdim[35]+f[1*nC]*alpha_vdim[34]+f[2*nC]*alpha_vdim[33]+f[6*nC]*alpha_vdim[32]+f[9*nC]*alpha_cdim[32]+alpha_cdim[3]*f[19*nC]+alpha_cdim[0]*f[10*nC]); 
  out[18*nC] += 0.3061862178478971*(f[22*nC]*alpha_vdim[52]+f[2*nC]*alpha_vdim[48]+f[27*nC]*alpha_vdim[45]+f[14*nC]*alpha_vdim[44]+f[6*nC]*alpha_vdim[40]+f[0]*alpha_vdim[39]+f[8*nC]*alpha_vdim[38]+f[21*nC]*alpha_vdim[37]+f[1*nC]*alpha_vdim[35]+f[16*nC]*alpha_vdim[34]+f[3*nC]*alpha_vdim[33]+f[7*nC]*alpha_vdim[32]+alpha_vdim[13]*f[28*nC]+alpha_vdim[20]*f[24*nC]+alpha_vdim[5]*f[23*nC]+alpha_vdim[2]*f[17*nC]+f[2*nC]*alpha_vdim[17]+alpha_vdim[12]*f[15*nC]+alpha_cdim[0]*f[11*nC]+alpha_vdim[6]*f[10*nC]+f[6*nC]*alpha_vdim[10]+alpha_vdim[0]*f[9*nC]+f[0]*alpha_vdim[9]+(alpha_cdim[3]+alpha_vdim[1])*f[4*nC]+f[1*nC]*alpha_vdim[4]); 
  out[19*nC] += 0.3061862178478971*(f[21*nC]*alpha_vdim[52]+f[1*nC]*alpha_vdim[48]+f[14*nC]*alpha_vdim[45]+f[27*nC]*alpha_vdim[44]+f[0]*alpha_vdim[40]+f[6*nC]*alpha_vdim[39]+f[7*nC]*alpha_vdim[38]+f[22*nC]*alpha_vdim[37]+f[3*nC]*alpha_cdim[36]+f[2*nC]*alpha_vdim[35]+f[3*nC]*alpha_vdim[34]+f[16*nC]*alpha_vdim[33]+f[8*nC]*alpha_vdim[32]+f[11*nC]*alpha_cdim[32]+alpha_vdim[12]*f[28*nC]+alpha_vdim[5]*f[24*nC]+alpha_vdim[20]*f[23*nC]+alpha_vdim[1]*f[17*nC]+f[1*nC]*alpha_vdim[17]+alpha_vdim[13]*f[15*nC]+alpha_vdim[0]*f[10*nC]+f[0]*alpha_vdim[10]+alpha_vdim[6]*f[9*nC]+f[6*nC]*alpha_vdim[9]+alpha_vdim[2]*f[4*nC]+f[2*nC]*alpha_vdim[4]); 
  out[20*nC] += 0.3061862178478971*(f[4*nC]*alpha_vdim[81]+f[3*nC]*alpha_vdim[80]+f[9*nC]*alpha_vdim[74]+f[10*nC]*alpha_vdim[73]+f[7*nC]*alpha_vdim[72]+f[8*nC]*alpha_vdim[71]+f[0]*alpha_vdim[70]+f[17*nC]*alpha_vdim[68]+f[16*nC]*alpha_vdim[67]+f[1*nC]*alpha_vdim[66]+f[2*nC]*alpha_vdim[65]+f[6*nC]*alpha_vdim[64]+f[23*nC]*alpha_cdim[36]+f[12*nC]*alpha_cdim[32]+alpha_cdim[3]*f[22*nC]+alpha_cdim[0]*f[13*nC]); 
  out[21*nC] += 0.3061862178478971*(f[19*nC]*alpha_vdim[81]+f[2*nC]*alpha_vdim[80]+f[26*nC]*alpha_vdim[74]+f[11*nC]*alpha_vdim[73]+f[6*nC]*alpha_vdim[72]+f[0]*alpha_vdim[71]+f[8*nC]*alpha_vdim[70]+f[18*nC]*alpha_vdim[68]+f[1*nC]*alpha_vdim[67]+f[16*nC]*alpha_vdim[66]+f[3*nC]*alpha_vdim[65]+f[7*nC]*alpha_vdim[64]+alpha_vdim[10]*f[28*nC]+alpha_vdim[17]*f[24*nC]+alpha_vdim[4]*f[23*nC]+alpha_vdim[2]*f[20*nC]+f[2*nC]*alpha_vdim[20]+alpha_vdim[9]*f[15*nC]+alpha_cdim[0]*f[14*nC]+alpha_vdim[6]*f[13*nC]+f[6*nC]*alpha_vdim[13]+alpha_vdim[0]*f[12*nC]+f[0]*alpha_vdim[12]+(alpha_cdim[3]+alpha_vdim[1])*f[5*nC]+f[1*nC]*alpha_vdim[5]); 
  out[22*nC] += 0.3061862178478971*(f[18*nC]*alpha_vdim[81]+f[1*nC]*alpha_vdim[80]+f[11*nC]*alpha_vdim[74]+f[26*nC]*alpha_vdim[73]+f[0]*alpha_vdim[72]+f[6*nC]*alpha_vdim[71]+f[7*nC]*alpha_vdim[70]+f[19*nC]*alpha_vdim[68]+f[2*nC]*alpha_vdim[67]+f[3*nC]*alpha_vdim[66]+f[16*nC]*alpha_vdim[65]+f[8*nC]*alpha_vdim[64]+f[25*nC]*alpha_cdim[36]+f[14*nC]*alpha_cdim[32]+alpha_vdim[9]*f[28*nC]+alpha_vdim[4]*f[24*nC]+alpha_vdim[17]*f[23*nC]+alpha_vdim[1]*f[20*nC]+f[1*nC]*alpha_vdim[20]+alpha_vdim[10]*f[15*nC]+alpha_vdim[0]*f[13*nC]+f[0]*alpha_vdim[13]+alpha_vdim[6]*f[12*nC]+f[6*nC]*alpha_vdim[12]+alpha_vdim[2]*f[5*nC]+f[2*nC]*alpha_vdim[5]); 
  out[23*nC] += 0.3061862178478971*(f[2*nC]*alpha_vdim[81]+f[19*nC]*alpha_vdim[80]+f[6*nC]*alpha_vdim[74]+f[0]*alpha_vdim[73]+f[26*nC]*alpha_vdim[72]+f[11*nC]*alpha_vdim[71]+f[10*nC]*alpha_vdim[70]+f[1*nC]*alpha_vdim[68]+f[18*nC]*alpha_vdim[67]+f[17*nC]*alpha_vdim[66]+f[4*nC]*alpha_vdim[65]+f[9*nC]*alpha_vdim[64]+f[2*nC]*alpha_vdim[52]+f[22*nC]*alpha_vdim[48]+f[6*nC]*alpha_vdim[45]+f[0]*alpha_vdim[44]+f[27*nC]*alpha_vdim[40]+f[14*nC]*alpha_vdim[39]+f[13*nC]*alpha_vdim[38]+f[1*nC]*alpha_vdim[37]+f[21*nC]*alpha_vdim[35]+f[20*nC]*alpha_vdim[34]+f[5*nC]*alpha_vdim[33]+f[12*nC]*alpha_vdim[32]+alpha_cdim[3]*f[25*nC]+alpha_cdim[0]*f[15*nC]); 
  out[24*nC] += 0.3061862178478971*(f[1*nC]*alpha_vdim[81]+f[18*nC]*alpha_vdim[80]+f[0]*alpha_vdim[74]+f[6*nC]*alpha_vdim[73]+f[11*nC]*alpha_vdim[72]+f[26*nC]*alpha_vdim[71]+f[9*nC]*alpha_vdim[70]+f[2*nC]*alpha_vdim[68]+f[19*nC]*alpha_vdim[67]+f[4*nC]*alpha_vdim[66]+f[17*nC]*alpha_vdim[65]+f[10*nC]*alpha_vdim[64]+f[1*nC]*alpha_vdim[52]+f[21*nC]*alpha_vdim[48]+f[0]*alpha_vdim[45]+f[6*nC]*alpha_vdim[44]+f[14*nC]*alpha_vdim[40]+f[27*nC]*alpha_vdim[39]+f[12*nC]*alpha_vdim[38]+f[2*nC]*alpha_vdim[37]+f[5*nC]*alpha_cdim[36]+f[22*nC]*alpha_vdim[35]+f[5*nC]*alpha_vdim[34]+f[20*nC]*alpha_vdim[33]+f[13*nC]*alpha_vdim[32]+f[15*nC]*alpha_cdim[32]); 
  out[25*nC] += 0.3061862178478971*(f[16*nC]*alpha_vdim[81]+f[17*nC]*alpha_vdim[80]+f[8*nC]*alpha_vdim[74]+f[7*nC]*alpha_vdim[73]+f[10*nC]*alpha_vdim[72]+f[9*nC]*alpha_vdim[71]+f[26*nC]*alpha_vdim[70]+f[3*nC]*alpha_vdim[68]+f[4*nC]*alpha_vdim[67]+f[19*nC]*alpha_vdim[66]+f[18*nC]*alpha_vdim[65]+f[11*nC]*alpha_vdim[64]+f[16*nC]*alpha_vdim[52]+f[20*nC]*alpha_vdim[48]+f[8*nC]*alpha_vdim[45]+f[7*nC]*alpha_vdim[44]+f[13*nC]*alpha_vdim[40]+f[12*nC]*alpha_vdim[39]+f[27*nC]*alpha_vdim[38]+f[3*nC]*alpha_vdim[37]+f[5*nC]*alpha_vdim[35]+f[22*nC]*alpha_vdim[34]+f[21*nC]*alpha_vdim[33]+f[14*nC]*alpha_vdim[32]+alpha_vdim[6]*f[28*nC]+alpha_vdim[2]*f[24*nC]+alpha_vdim[1]*f[23*nC]+alpha_vdim[17]*f[20*nC]+f[17*nC]*alpha_vdim[20]+alpha_vdim[0]*f[15*nC]+alpha_vdim[10]*f[13*nC]+f[10*nC]*alpha_vdim[13]+alpha_vdim[9]*f[12*nC]+f[9*nC]*alpha_vdim[12]+alpha_vdim[4]*f[5*nC]+f[4*nC]*alpha_vdim[5]); 
  out[26*nC] += 0.3061862178478971*(f[14*nC]*alpha_vdim[52]+f[0]*alpha_vdim[48]+f[21*nC]*alpha_vdim[45]+f[22*nC]*alpha_vdim[44]+f[1*nC]*alpha_vdim[40]+f[2*nC]*alpha_vdim[39]+f[3*nC]*alpha_vdim[38]+f[27*nC]*alpha_vdim[37]+f[7*nC]*alpha_cdim[36]+f[6*nC]*alpha_vdim[35]+f[7*nC]*alpha_vdim[34]+f[8*nC]*alpha_vdim[33]+f[16*nC]*alpha_vdim[32]+f[18*nC]*alpha_cdim[32]+alpha_vdim[5]*f[28*nC]+alpha_vdim[12]*f[24*nC]+alpha_vdim[13]*f[23*nC]+f[15*nC]*alpha_vdim[20]+alpha_cdim[0]*f[19*nC]+alpha_vdim[0]*f[17*nC]+f[0]*alpha_vdim[17]+(alpha_cdim[3]+alpha_vdim[1])*f[10*nC]+f[1*nC]*alpha_vdim[10]+alpha_vdim[2]*f[9*nC]+f[2*nC]*alpha_vdim[9]+alpha_vdim[4]*f[6*nC]+f[4*nC]*alpha_vdim[6]); 
  out[27*nC] += 0.3061862178478971*(f[11*nC]*alpha_vdim[81]+f[0]*alpha_vdim[80]+f[18*nC]*alpha_vdim[74]+f[19*nC]*alpha_vdim[73]+f[1*nC]*alpha_vdim[72]+f[2*nC]*alpha_vdim[71]+f[3*nC]*alpha_vdim[70]+f[26*nC]*alpha_vdim[68]+f[6*nC]*alpha_vdim[67]+f[7*nC]*alpha_vdim[66]+f[8*nC]*alpha_vdim[65]+f[16*nC]*alpha_vdim[64]+f[29*nC]*alpha_cdim[36]+f[21*nC]*alpha_cdim[32]+alpha_vdim[4]*f[28*nC]+alpha_vdim[9]*f[24*nC]+alpha_vdim[10]*f[23*nC]+alpha_cdim[0]*f[22*nC]+alpha_vdim[0]*f[20*nC]+f[0]*alpha_vdim[20]+f[15*nC]*alpha_vdim[17]+(alpha_cdim[3]+alpha_vdim[1])*f[13*nC]+f[1*nC]*alpha_vdim[13]+alpha_vdim[2]*f[12*nC]+f[2*nC]*alpha_vdim[12]+alpha_vdim[5]*f[6*nC]+f[5*nC]*alpha_vdim[6]); 
  out[28*nC] += 0.3061862178478971*(f[0]*alpha_vdim[81]+f[11*nC]*alpha_vdim[80]+f[1*nC]*alpha_vdim[74]+f[2*nC]*alpha_vdim[73]+f[18*nC]*alpha_vdim[72]+f[19*nC]*alpha_vdim[71]+f[4*nC]*alpha_vdim[70]+f[6*nC]*alpha_vdim[68]+f[26*nC]*alpha_vdim[67]+f[9*nC]*alpha_vdim[66]+f[10*nC]*alpha_vdim[65]+f[17*nC]*alpha_vdim[64]+f[0]*alpha_vdim[52]+f[14*nC]*alpha_vdim[48]+f[1*nC]*alpha_vdim[45]+f[2*nC]*alpha_vdim[44]+f[21*nC]*alpha_vdim[40]+f[22*nC]*alpha_vdim[39]+f[5*nC]*alpha_vdim[38]+f[6*nC]*alpha_vdim[37]+f[12*nC]*alpha_cdim[36]+f[27*nC]*alpha_vdim[35]+f[12*nC]*alpha_vdim[34]+f[13*nC]*alpha_vdim[33]+f[20*nC]*alpha_vdim[32]+f[23*nC]*alpha_cdim[32]+alpha_cdim[3]*f[30*nC]+alpha_cdim[0]*f[24*nC]); 
  out[29*nC] += 0.3061862178478971*(f[8*nC]*alpha_vdim[81]+f[10*nC]*alpha_vdim[80]+f[16*nC]*alpha_vdim[74]+f[3*nC]*alpha_vdim[73]+f[17*nC]*alpha_vdim[72]+f[4*nC]*alpha_vdim[71]+f[19*nC]*alpha_vdim[70]+f[7*nC]*alpha_vdim[68]+f[9*nC]*alpha_vdim[67]+f[26*nC]*alpha_vdim[66]+f[11*nC]*alpha_vdim[65]+f[18*nC]*alpha_vdim[64]+f[8*nC]*alpha_vdim[52]+f[13*nC]*alpha_vdim[48]+f[16*nC]*alpha_vdim[45]+f[3*nC]*alpha_vdim[44]+f[20*nC]*alpha_vdim[40]+f[5*nC]*alpha_vdim[39]+f[22*nC]*alpha_vdim[38]+f[7*nC]*alpha_vdim[37]+f[12*nC]*alpha_vdim[35]+f[27*nC]*alpha_vdim[34]+f[14*nC]*alpha_vdim[33]+f[21*nC]*alpha_vdim[32]+alpha_vdim[2]*f[28*nC]+alpha_cdim[0]*f[25*nC]+alpha_vdim[6]*f[24*nC]+alpha_vdim[0]*f[23*nC]+alpha_vdim[10]*f[20*nC]+f[10*nC]*alpha_vdim[20]+alpha_vdim[13]*f[17*nC]+f[13*nC]*alpha_vdim[17]+(alpha_cdim[3]+alpha_vdim[1])*f[15*nC]+alpha_vdim[4]*f[12*nC]+f[4*nC]*alpha_vdim[12]+alpha_vdim[5]*f[9*nC]+f[5*nC]*alpha_vdim[9]); 
  out[30*nC] += 0.3061862178478971*(f[7*nC]*alpha_vdim[81]+f[9*nC]*alpha_vdim[80]+f[3*nC]*alpha_vdim[74]+f[16*nC]*alpha_vdim[73]+f[4*nC]*alpha_vdim[72]+f[17*nC]*alpha_vdim[71]+f[18*nC]*alpha_vdim[70]+f[8*nC]*alpha_vdim[68]+f[10*nC]*alpha_vdim[67]+f[11*nC]*alpha_vdim[66]+f[26*nC]*alpha_vdim[65]+f[19*nC]*alpha_vdim[64]+f[7*nC]*alpha_vdim[52]+f[12*nC]*alpha_vdim[48]+f[3*nC]*alpha_vdim[45]+f[16*nC]*alpha_vdim[44]+f[5*nC]*alpha_vdim[40]+f[20*nC]*alpha_vdim[39]+f[21*nC]*alpha_vdim[38]+f[8*nC]*alpha_vdim[37]+f[14*nC]*alpha_cdim[36]+f[13*nC]*alpha_vdim[35]+f[14*nC]*alpha_vdim[34]+f[27*nC]*alpha_vdim[33]+f[22*nC]*alpha_vdim[32]+f[25*nC]*alpha_cdim[32]+alpha_vdim[1]*f[28*nC]+alpha_vdim[0]*f[24*nC]+alpha_vdim[6]*f[23*nC]+alpha_vdim[9]*f[20*nC]+f[9*nC]*alpha_vdim[20]+alpha_vdim[12]*f[17*nC]+f[12*nC]*alpha_vdim[17]+alpha_vdim[2]*f[15*nC]+alpha_vdim[4]*f[13*nC]+f[4*nC]*alpha_vdim[13]+alpha_vdim[5]*f[10*nC]+f[5*nC]*alpha_vdim[10]); 
  out[31*nC] += 0.3061862178478971*(f[3*nC]*alpha_vdim[81]+f[4*nC]*alpha_vdim[80]+f[7*nC]*alpha_vdim[74]+f[8*nC]*alpha_vdim[73]+f[9*nC]*alpha_vdim[72]+f[10*nC]*alpha_vdim[71]+f[11*nC]*alpha_vdim[70]+f[16*nC]*alpha_vdim[68]+f[17*nC]*alpha_vdim[67]+f[18*nC]*alpha_vdim[66]+f[19*nC]*alpha_vdim[65]+f[26*nC]*alpha_vdim[64]+f[3*nC]*alpha_vdim[52]+f[5*nC]*alpha_vdim[48]+f[7*nC]*alpha_vdim[45]+f[8*nC]*alpha_vdim[44]+f[12*nC]*alpha_vdim[40]+f[13*nC]*alpha_vdim[39]+f[14*nC]*alpha_vdim[38]+f[16*nC]*alpha_vdim[37]+f[21*nC]*alpha_cdim[36]+f[20*nC]*alpha_vdim[35]+f[21*nC]*alpha_vdim[34]+f[22*nC]*alpha_vdim[33]+f[27*nC]*alpha_vdim[32]+f[29*nC]*alpha_cdim[32]+alpha_cdim[0]*f[30*nC]+alpha_vdim[0]*f[28*nC]+(alpha_cdim[3]+alpha_vdim[1])*f[24*nC]+alpha_vdim[2]*f[23*nC]+alpha_vdim[4]*f[20*nC]+f[4*nC]*alpha_vdim[20]+alpha_vdim[5]*f[17*nC]+f[5*nC]*alpha_vdim[17]+alpha_vdim[6]*f[15*nC]+alpha_vdim[9]*f[13*nC]+f[9*nC]*alpha_vdim[13]+alpha_vdim[10]*f[12*nC]+f[10*nC]*alpha_vdim[12]); 

  return alpha_mid; 
} 
