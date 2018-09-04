#include <VlasovModDecl.h> 
double VlasovVol2x2vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double *B2 = &EM[15]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[10]; 
  double alpha_vdim[10]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[3] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_cdim[5] = 8.0*w1dx1; 
  alpha_cdim[9] = 2.309401076758503*dv1dx1; 
  alpha_mid += std::abs(w1dx1)+0.5*dv1dx1; 
  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[2] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_vdim[4] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_mid += std::abs(0.125*alpha_vdim[0]); 
  alpha_vdim[5] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[6] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[7] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_vdim[8] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_mid += std::abs(0.125*alpha_vdim[5]); 
  out[1] += 0.4330127018922193*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[4]*alpha_cdim[9]+f[0]*alpha_cdim[5]); 
  out[3] += 0.4330127018922193*(alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[3]*alpha_vdim[8]+f[2]*alpha_vdim[7]+f[1]*alpha_vdim[6]+f[0]*alpha_vdim[5]); 
return alpha_mid; 

} 
double VlasovVol2x2vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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
  const double *E1 = &EM[6]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B2 = &EM[30]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[30]; 
  double alpha_vdim[30]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[3] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_cdim[15] = 8.0*w1dx1; 
  alpha_cdim[19] = 2.309401076758503*dv1dx1; 
  alpha_mid += std::abs(w1dx1)+0.5*dv1dx1; 
  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[2] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_vdim[4] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 2.0*dv10*(B2[3]*wv2+E0[3]); 
  alpha_vdim[8] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[9] = 0.5773502691896258*B2[2]*dv10*dv2; 
  alpha_vdim[11] = 2.0*dv10*(B2[4]*wv2+E0[4]); 
  alpha_vdim[12] = 2.0*dv10*(B2[5]*wv2+E0[5]); 
  alpha_mid += std::abs(0.125*alpha_vdim[0]-0.1397542485937369*(alpha_vdim[12]+alpha_vdim[11])); 
  alpha_vdim[15] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[16] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[17] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_vdim[18] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[20] = dv11*(2.0*E1[3]-2.0*B2[3]*wv1); 
  alpha_vdim[21] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[22] = -0.5773502691896258*B2[2]*dv1*dv11; 
  alpha_vdim[26] = dv11*(2.0*E1[4]-2.0*B2[4]*wv1); 
  alpha_vdim[27] = dv11*(2.0*E1[5]-2.0*B2[5]*wv1); 
  alpha_mid += std::abs(0.125*alpha_vdim[15]-0.1397542485937369*(alpha_vdim[27]+alpha_vdim[26])); 
  out[1] += 0.4330127018922193*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[4]*alpha_cdim[19]+f[0]*alpha_cdim[15]); 
  out[3] += 0.4330127018922193*(alpha_vdim[12]*f[12]+alpha_vdim[11]*f[11]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[12]*alpha_vdim[27]+f[11]*alpha_vdim[26]+f[7]*alpha_vdim[22]+f[6]*alpha_vdim[21]+f[5]*alpha_vdim[20]+f[3]*alpha_vdim[18]+f[2]*alpha_vdim[17]+f[1]*alpha_vdim[16]+f[0]*alpha_vdim[15]); 
  out[5] += 0.4330127018922193*(f[8]*alpha_cdim[19]+f[1]*alpha_cdim[15]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.3872983346207416*(alpha_cdim[3]*f[13]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.4330127018922193*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[7] += 0.4330127018922193*(f[10]*alpha_cdim[19]+f[3]*alpha_cdim[15])+0.3872983346207416*(alpha_vdim[2]*f[12]+f[2]*alpha_vdim[12])+0.4330127018922193*(alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[8] += 0.3872983346207416*f[1]*alpha_vdim[26]+0.4330127018922193*(f[3]*alpha_vdim[21]+f[2]*alpha_vdim[20]+f[6]*alpha_vdim[18]+f[5]*alpha_vdim[17])+0.3872983346207416*f[11]*alpha_vdim[16]+0.4330127018922193*(f[0]*alpha_vdim[16]+f[1]*alpha_vdim[15]+alpha_cdim[3]*f[10]+alpha_cdim[0]*f[4]); 
  out[9] += 0.3872983346207416*f[2]*alpha_vdim[27]+0.4330127018922193*(f[3]*alpha_vdim[22]+f[1]*alpha_vdim[20])+0.3872983346207416*f[14]*alpha_cdim[19]+0.4330127018922193*(f[0]*alpha_cdim[19]+f[7]*alpha_vdim[18])+0.3872983346207416*f[12]*alpha_vdim[17]+0.4330127018922193*(f[0]*alpha_vdim[17]+f[5]*alpha_vdim[16]+f[2]*alpha_vdim[15]+f[4]*alpha_cdim[15]); 
  out[10] += 0.4330127018922193*(f[2]*alpha_vdim[22]+f[1]*alpha_vdim[21])+0.3872983346207416*f[13]*alpha_vdim[18]+0.4330127018922193*(f[0]*alpha_vdim[18]+f[7]*alpha_vdim[17]+f[6]*alpha_vdim[16]+f[3]*alpha_vdim[15])+0.3872983346207416*alpha_vdim[4]*f[14]+0.4330127018922193*(alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[11] += 0.9682458365518543*(alpha_cdim[3]*f[6]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(f[9]*alpha_cdim[19]+f[2]*alpha_cdim[15]); 
  out[13] += 0.9682458365518543*(alpha_vdim[4]*f[10]+alpha_vdim[2]*f[7]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[14] += 0.9682458365518543*(f[10]*alpha_vdim[18]+f[9]*alpha_vdim[17]+f[8]*alpha_vdim[16]+f[4]*alpha_vdim[15]); 
return alpha_mid; 

} 
double VlasovVol2x2vMaxP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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
  const double *E1 = &EM[10]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B2 = &EM[50]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[70]; 
  double alpha_vdim[70]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[3] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_cdim[35] = 8.0*w1dx1; 
  alpha_cdim[39] = 2.309401076758503*dv1dx1; 
  alpha_mid += std::abs(w1dx1)+0.5*dv1dx1; 
  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[2] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_vdim[4] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 2.0*dv10*(B2[3]*wv2+E0[3]); 
  alpha_vdim[8] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[9] = 0.5773502691896258*B2[2]*dv10*dv2; 
  alpha_vdim[11] = 2.0*dv10*(B2[4]*wv2+E0[4]); 
  alpha_vdim[12] = 2.0*dv10*(B2[5]*wv2+E0[5]); 
  alpha_vdim[16] = 0.5773502691896258*B2[3]*dv10*dv2; 
  alpha_vdim[19] = 2.0*dv10*(B2[6]*wv2+E0[6]); 
  alpha_vdim[20] = 2.0*dv10*(B2[7]*wv2+E0[7]); 
  alpha_vdim[25] = 0.5773502691896258*B2[4]*dv10*dv2; 
  alpha_vdim[26] = 0.5773502691896258*B2[5]*dv10*dv2; 
  alpha_vdim[31] = 2.0*dv10*(B2[8]*wv2+E0[8]); 
  alpha_vdim[32] = 2.0*dv10*(B2[9]*wv2+E0[9]); 
  alpha_mid += std::abs(0.125*alpha_vdim[0]-0.1397542485937369*(alpha_vdim[12]+alpha_vdim[11])); 
  alpha_vdim[35] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[36] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[37] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_vdim[38] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[40] = dv11*(2.0*E1[3]-2.0*B2[3]*wv1); 
  alpha_vdim[41] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[42] = -0.5773502691896258*B2[2]*dv1*dv11; 
  alpha_vdim[46] = dv11*(2.0*E1[4]-2.0*B2[4]*wv1); 
  alpha_vdim[47] = dv11*(2.0*E1[5]-2.0*B2[5]*wv1); 
  alpha_vdim[50] = -0.5773502691896258*B2[3]*dv1*dv11; 
  alpha_vdim[54] = dv11*(2.0*E1[6]-2.0*B2[6]*wv1); 
  alpha_vdim[55] = dv11*(2.0*E1[7]-2.0*B2[7]*wv1); 
  alpha_vdim[56] = -0.5773502691896258*B2[4]*dv1*dv11; 
  alpha_vdim[57] = -0.5773502691896258*B2[5]*dv1*dv11; 
  alpha_vdim[66] = dv11*(2.0*E1[8]-2.0*B2[8]*wv1); 
  alpha_vdim[67] = dv11*(2.0*E1[9]-2.0*B2[9]*wv1); 
  alpha_mid += std::abs(0.125*alpha_vdim[35]-0.1397542485937369*(alpha_vdim[47]+alpha_vdim[46])); 
  out[1] += 0.4330127018922193*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(f[4]*alpha_cdim[39]+f[0]*alpha_cdim[35]); 
  out[3] += 0.4330127018922193*(alpha_vdim[32]*f[32]+alpha_vdim[31]*f[31]+alpha_vdim[26]*f[26]+alpha_vdim[25]*f[25]+alpha_vdim[20]*f[20]+alpha_vdim[19]*f[19]+alpha_vdim[16]*f[16]+alpha_vdim[12]*f[12]+alpha_vdim[11]*f[11]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[32]*alpha_vdim[67]+f[31]*alpha_vdim[66]+f[22]*alpha_vdim[57]+f[21]*alpha_vdim[56]+f[20]*alpha_vdim[55]+f[19]*alpha_vdim[54]+f[15]*alpha_vdim[50]+f[12]*alpha_vdim[47]+f[11]*alpha_vdim[46]+f[7]*alpha_vdim[42]+f[6]*alpha_vdim[41]+f[5]*alpha_vdim[40]+f[3]*alpha_vdim[38]+f[2]*alpha_vdim[37]+f[1]*alpha_vdim[36]+f[0]*alpha_vdim[35]); 
  out[5] += 0.4330127018922193*(f[8]*alpha_cdim[39]+f[1]*alpha_cdim[35]+alpha_cdim[3]*f[7]+alpha_cdim[0]*f[2]); 
  out[6] += 0.3803194146278324*(alpha_vdim[11]*f[31]+f[11]*alpha_vdim[31])+0.3872983346207416*(alpha_vdim[8]*f[25]+f[8]*alpha_vdim[25])+0.4330127018922193*(alpha_vdim[12]*f[20]+f[12]*alpha_vdim[20])+0.3872983346207416*(alpha_vdim[5]*f[19]+f[5]*alpha_vdim[19])+0.4330127018922193*(alpha_vdim[9]*f[16]+f[9]*alpha_vdim[16])+0.3872983346207416*(alpha_cdim[3]*f[13]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.4330127018922193*(alpha_vdim[4]*f[8]+f[4]*alpha_vdim[8]+alpha_vdim[2]*f[5]+f[2]*alpha_vdim[5]+alpha_cdim[0]*f[3]+f[0]*alpha_cdim[3]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[7] += 0.4330127018922193*(f[10]*alpha_cdim[39]+f[3]*alpha_cdim[35])+0.3803194146278324*(alpha_vdim[12]*f[32]+f[12]*alpha_vdim[32])+0.3872983346207416*(alpha_vdim[9]*f[26]+f[9]*alpha_vdim[26]+alpha_vdim[5]*f[20]+f[5]*alpha_vdim[20])+0.4330127018922193*(alpha_vdim[11]*f[19]+f[11]*alpha_vdim[19]+alpha_vdim[8]*f[16]+f[8]*alpha_vdim[16])+0.3872983346207416*(alpha_vdim[2]*f[12]+f[2]*alpha_vdim[12])+0.4330127018922193*(alpha_vdim[4]*f[9]+f[4]*alpha_vdim[9]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[2]+f[0]*alpha_vdim[2]); 
  out[8] += 0.3803194146278324*f[11]*alpha_vdim[66]+0.3872983346207416*f[6]*alpha_vdim[56]+0.4330127018922193*f[12]*alpha_vdim[55]+0.3872983346207416*f[5]*alpha_vdim[54]+0.4330127018922193*(f[7]*alpha_vdim[50]+f[20]*alpha_vdim[47])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[46]+0.4330127018922193*f[15]*alpha_vdim[42]+(0.3872983346207416*f[21]+0.4330127018922193*f[3])*alpha_vdim[41]+0.3872983346207416*f[19]*alpha_vdim[40]+0.4330127018922193*(f[2]*alpha_vdim[40]+f[6]*alpha_vdim[38]+f[5]*alpha_vdim[37])+0.3872983346207416*f[11]*alpha_vdim[36]+0.4330127018922193*(f[0]*alpha_vdim[36]+f[1]*alpha_vdim[35]+alpha_cdim[3]*f[10]+alpha_cdim[0]*f[4]); 
  out[9] += 0.3803194146278324*f[12]*alpha_vdim[67]+0.3872983346207416*(f[7]*alpha_vdim[57]+f[5]*alpha_vdim[55])+0.4330127018922193*(f[11]*alpha_vdim[54]+f[6]*alpha_vdim[50])+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alpha_vdim[47]+0.4330127018922193*f[19]*alpha_vdim[46]+0.3872983346207416*f[22]*alpha_vdim[42]+0.4330127018922193*(f[3]*alpha_vdim[42]+f[15]*alpha_vdim[41])+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[40]+0.3872983346207416*f[14]*alpha_cdim[39]+0.4330127018922193*(f[0]*alpha_cdim[39]+f[7]*alpha_vdim[38])+0.3872983346207416*f[12]*alpha_vdim[37]+0.4330127018922193*(f[0]*alpha_vdim[37]+f[5]*alpha_vdim[36]+f[2]*alpha_vdim[35]+f[4]*alpha_cdim[35]); 
  out[10] += 0.4330127018922193*(f[12]*alpha_vdim[57]+f[11]*alpha_vdim[56]+f[5]*alpha_vdim[50]+f[22]*alpha_vdim[47]+f[21]*alpha_vdim[46])+(0.3872983346207416*f[24]+0.4330127018922193*f[2])*alpha_vdim[42]+0.3872983346207416*f[23]*alpha_vdim[41]+0.4330127018922193*(f[1]*alpha_vdim[41]+f[15]*alpha_vdim[40])+0.3872983346207416*f[13]*alpha_vdim[38]+0.4330127018922193*(f[0]*alpha_vdim[38]+f[7]*alpha_vdim[37]+f[6]*alpha_vdim[36]+f[3]*alpha_vdim[35])+0.3872983346207416*(alpha_vdim[9]*f[29]+alpha_vdim[8]*f[28])+0.4330127018922193*(alpha_vdim[12]*f[26]+f[12]*alpha_vdim[26]+alpha_vdim[11]*f[25]+f[11]*alpha_vdim[25]+alpha_vdim[5]*f[16]+f[5]*alpha_vdim[16])+0.3872983346207416*alpha_vdim[4]*f[14]+0.4330127018922193*(alpha_vdim[2]*f[9]+f[2]*alpha_vdim[9]+alpha_vdim[1]*f[8]+f[1]*alpha_vdim[8]+alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4]); 
  out[11] += 0.9682458365518543*(alpha_cdim[3]*f[6]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(f[9]*alpha_cdim[39]+f[2]*alpha_cdim[35]); 
  out[13] += 0.9682458365518543*(alpha_vdim[12]*f[22]+alpha_vdim[11]*f[21]+alpha_vdim[9]*f[18]+alpha_vdim[8]*f[17]+alpha_vdim[5]*f[15]+alpha_vdim[4]*f[10]+alpha_vdim[2]*f[7]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[14] += 0.9682458365518543*(f[26]*alpha_vdim[47]+f[25]*alpha_vdim[46]+f[18]*alpha_vdim[42]+f[17]*alpha_vdim[41]+f[16]*alpha_vdim[40]+f[10]*alpha_vdim[38]+f[9]*alpha_vdim[37]+f[8]*alpha_vdim[36]+f[4]*alpha_vdim[35]); 
  out[15] += 0.4330127018922193*(f[17]*alpha_cdim[39]+f[6]*alpha_cdim[35])+0.3803194146278324*(alpha_vdim[20]*f[32]+f[20]*alpha_vdim[32]+alpha_vdim[19]*f[31]+f[19]*alpha_vdim[31])+0.3872983346207416*(alpha_vdim[16]*f[26]+f[16]*alpha_vdim[26]+alpha_vdim[16]*f[25]+f[16]*alpha_vdim[25]+alpha_cdim[3]*f[24])+(0.3464101615137755*alpha_vdim[19]+0.3872983346207416*alpha_vdim[2])*f[20]+0.3464101615137755*f[19]*alpha_vdim[20]+0.3872983346207416*(f[2]*alpha_vdim[20]+alpha_vdim[1]*f[19]+f[1]*alpha_vdim[19])+0.4330127018922193*(alpha_vdim[4]*f[16]+f[4]*alpha_vdim[16])+0.3872983346207416*(alpha_vdim[5]*f[12]+f[5]*alpha_vdim[12]+alpha_vdim[5]*f[11]+f[5]*alpha_vdim[11])+0.4330127018922193*(alpha_vdim[8]*f[9]+f[8]*alpha_vdim[9]+alpha_cdim[0]*f[7]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+f[2]*(alpha_cdim[3]+alpha_vdim[1])+f[1]*alpha_vdim[2]); 
  out[16] += 0.3803194146278324*(f[20]*alpha_vdim[67]+f[19]*alpha_vdim[66])+0.3872983346207416*f[15]*(alpha_vdim[57]+alpha_vdim[56])+(0.3803194146278324*f[32]+0.3464101615137755*f[19]+0.3872983346207416*f[2])*alpha_vdim[55]+(0.3803194146278324*f[31]+0.3464101615137755*f[20]+0.3872983346207416*f[1])*alpha_vdim[54]+(0.3872983346207416*(f[22]+f[21])+0.4330127018922193*f[3])*alpha_vdim[50]+0.3872983346207416*f[5]*(alpha_vdim[47]+alpha_vdim[46])+0.4330127018922193*(f[6]*alpha_vdim[42]+f[7]*alpha_vdim[41])+(0.3872983346207416*(f[12]+f[11])+0.4330127018922193*f[0])*alpha_vdim[40]+0.3872983346207416*f[28]*alpha_cdim[39]+0.4330127018922193*(f[1]*alpha_cdim[39]+f[15]*alpha_vdim[38])+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alpha_vdim[37]+0.3872983346207416*f[19]*alpha_vdim[36]+0.4330127018922193*(f[2]*alpha_vdim[36]+f[5]*alpha_vdim[35]+f[8]*alpha_cdim[35]+alpha_cdim[3]*f[18]+alpha_cdim[0]*f[9]); 
  out[17] += 0.3803194146278324*f[21]*alpha_vdim[66]+0.4330127018922193*f[20]*alpha_vdim[57]+(0.3803194146278324*f[31]+0.3464101615137755*f[23]+0.3872983346207416*f[1])*alpha_vdim[56]+0.4330127018922193*f[22]*alpha_vdim[55]+0.3872983346207416*f[15]*alpha_vdim[54]+(0.3872983346207416*(f[24]+f[19])+0.4330127018922193*f[2])*alpha_vdim[50]+0.3872983346207416*f[6]*alpha_vdim[46]+0.4330127018922193*f[5]*alpha_vdim[42]+0.3872983346207416*(f[13]+f[11])*alpha_vdim[41]+0.4330127018922193*(f[0]*alpha_vdim[41]+f[7]*alpha_vdim[40])+0.3872983346207416*f[23]*alpha_vdim[38]+0.4330127018922193*(f[1]*alpha_vdim[38]+f[15]*alpha_vdim[37])+0.3872983346207416*f[21]*alpha_vdim[36]+0.4330127018922193*(f[3]*alpha_vdim[36]+f[6]*alpha_vdim[35])+0.3803194146278324*(alpha_vdim[25]*f[31]+f[25]*alpha_vdim[31])+0.3872983346207416*alpha_vdim[16]*f[29]+0.3464101615137755*alpha_vdim[25]*f[28]+0.3872983346207416*(alpha_vdim[4]*f[28]+alpha_cdim[3]*f[27])+0.4330127018922193*(alpha_vdim[20]*f[26]+f[20]*alpha_vdim[26])+0.3872983346207416*(alpha_vdim[1]*f[25]+f[1]*alpha_vdim[25]+alpha_vdim[16]*f[19]+f[16]*alpha_vdim[19])+0.4330127018922193*(alpha_vdim[2]*f[16]+f[2]*alpha_vdim[16])+0.3872983346207416*(alpha_vdim[8]*(f[14]+f[11])+f[8]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*f[10]+alpha_vdim[5]*f[9]+f[5]*alpha_vdim[9]+alpha_vdim[0]*f[8]+f[0]*alpha_vdim[8]+(alpha_cdim[3]+alpha_vdim[1])*f[4]+f[1]*alpha_vdim[4]); 
  out[18] += 0.3803194146278324*f[22]*alpha_vdim[67]+(0.3803194146278324*f[32]+0.3464101615137755*f[24]+0.3872983346207416*f[2])*alpha_vdim[57]+0.4330127018922193*f[19]*alpha_vdim[56]+0.3872983346207416*f[15]*alpha_vdim[55]+0.4330127018922193*f[21]*alpha_vdim[54]+(0.3872983346207416*(f[23]+f[20])+0.4330127018922193*f[1])*alpha_vdim[50]+0.3872983346207416*(f[7]*alpha_vdim[47]+(f[13]+f[12])*alpha_vdim[42])+0.4330127018922193*(f[0]*alpha_vdim[42]+f[5]*alpha_vdim[41]+f[6]*alpha_vdim[40])+(0.3872983346207416*f[30]+0.4330127018922193*f[3])*alpha_cdim[39]+(0.3872983346207416*f[24]+0.4330127018922193*f[2])*alpha_vdim[38]+0.3872983346207416*f[22]*alpha_vdim[37]+0.4330127018922193*(f[3]*alpha_vdim[37]+f[15]*alpha_vdim[36]+f[7]*alpha_vdim[35]+f[10]*alpha_cdim[35])+0.3803194146278324*(alpha_vdim[26]*f[32]+f[26]*alpha_vdim[32])+0.3464101615137755*alpha_vdim[26]*f[29]+0.3872983346207416*(alpha_vdim[4]*f[29]+alpha_vdim[16]*f[28]+alpha_vdim[2]*f[26]+f[2]*alpha_vdim[26])+0.4330127018922193*(alpha_vdim[19]*f[25]+f[19]*alpha_vdim[25])+0.3872983346207416*(alpha_vdim[16]*f[20]+f[16]*alpha_vdim[20])+0.4330127018922193*(alpha_vdim[1]*f[16]+f[1]*alpha_vdim[16])+0.3872983346207416*(alpha_vdim[9]*(f[14]+f[12])+f[9]*alpha_vdim[12])+0.4330127018922193*(alpha_vdim[0]*f[9]+f[0]*alpha_vdim[9]+alpha_vdim[5]*f[8]+f[5]*alpha_vdim[8]+alpha_vdim[2]*f[4]+f[2]*alpha_vdim[4]); 
  out[19] += 0.4330127018922193*(f[25]*alpha_cdim[39]+f[11]*alpha_cdim[35])+0.9682458365518543*(alpha_cdim[3]*f[15]+alpha_cdim[0]*f[5]); 
  out[20] += 0.9682458365518543*(f[16]*alpha_cdim[39]+f[5]*alpha_cdim[35])+0.4330127018922193*(alpha_cdim[3]*f[22]+alpha_cdim[0]*f[12]); 
  out[21] += 0.2581988897471612*alpha_vdim[31]*f[31]+0.3803194146278324*(alpha_vdim[1]*f[31]+f[1]*alpha_vdim[31])+0.276641667586244*alpha_vdim[25]*f[25]+0.4330127018922193*(alpha_vdim[4]*f[25]+f[4]*alpha_vdim[25])+0.8660254037844386*alpha_cdim[3]*f[23]+0.3872983346207416*alpha_vdim[20]*f[20]+0.276641667586244*alpha_vdim[19]*f[19]+0.4330127018922193*(alpha_vdim[2]*f[19]+f[2]*alpha_vdim[19])+0.3872983346207416*alpha_vdim[16]*f[16]+0.276641667586244*alpha_vdim[11]*f[11]+0.4330127018922193*(alpha_vdim[0]*f[11]+f[0]*alpha_vdim[11])+0.3872983346207416*alpha_vdim[8]*f[8]+0.9682458365518543*alpha_cdim[0]*f[6]+0.3872983346207416*alpha_vdim[5]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[3]+0.3872983346207416*alpha_vdim[1]); 
  out[22] += 0.9682458365518543*(f[18]*alpha_cdim[39]+f[7]*alpha_cdim[35])+0.2581988897471612*alpha_vdim[32]*f[32]+0.3803194146278324*(alpha_vdim[2]*f[32]+f[2]*alpha_vdim[32])+0.276641667586244*alpha_vdim[26]*f[26]+0.4330127018922193*(alpha_vdim[4]*f[26]+f[4]*alpha_vdim[26])+0.276641667586244*alpha_vdim[20]*f[20]+0.4330127018922193*(alpha_vdim[1]*f[20]+f[1]*alpha_vdim[20])+0.3872983346207416*(alpha_vdim[19]*f[19]+alpha_vdim[16]*f[16])+0.276641667586244*alpha_vdim[12]*f[12]+0.4330127018922193*(alpha_vdim[0]*f[12]+f[0]*alpha_vdim[12])+0.3872983346207416*(alpha_vdim[9]*f[9]+alpha_vdim[5]*f[5]+alpha_vdim[2]*f[2]); 
  out[23] += 0.3803194146278324*alpha_cdim[3]*f[33]+0.8504200642707612*f[21]*alpha_vdim[31]+0.8660254037844386*f[17]*alpha_vdim[25]+0.9682458365518543*alpha_vdim[20]*f[22]+0.8660254037844386*(alpha_vdim[1]*f[21]+f[15]*alpha_vdim[19])+0.9682458365518543*(alpha_vdim[16]*f[18]+alpha_vdim[4]*f[17]+alpha_vdim[2]*f[15])+0.4330127018922193*alpha_cdim[0]*f[13]+0.8660254037844386*f[6]*alpha_vdim[11]+0.9682458365518543*(alpha_vdim[8]*f[10]+alpha_vdim[5]*f[7]+alpha_vdim[0]*f[6])+(0.3872983346207416*alpha_cdim[3]+0.9682458365518543*alpha_vdim[1])*f[3]; 
  out[24] += 0.4330127018922193*(f[27]*alpha_cdim[39]+f[13]*alpha_cdim[35])+0.8504200642707612*f[22]*alpha_vdim[32]+0.8660254037844386*(f[18]*alpha_vdim[26]+alpha_vdim[2]*f[22])+0.9682458365518543*alpha_vdim[19]*f[21]+0.8660254037844386*f[15]*alpha_vdim[20]+0.9682458365518543*(alpha_vdim[4]*f[18]+alpha_vdim[16]*f[17]+alpha_vdim[1]*f[15])+0.8660254037844386*f[7]*alpha_vdim[12]+0.9682458365518543*(alpha_vdim[9]*f[10]+alpha_vdim[0]*f[7]+alpha_vdim[5]*f[6]+alpha_vdim[2]*f[3]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_vdim[66]+(0.276641667586244*f[21]+0.4330127018922193*f[3])*alpha_vdim[56]+0.3872983346207416*f[20]*alpha_vdim[55]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_vdim[54]+0.3872983346207416*f[15]*alpha_vdim[50]+(0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_vdim[46]+0.3872983346207416*(f[6]*alpha_vdim[41]+f[5]*alpha_vdim[40])+0.4330127018922193*(f[21]*alpha_vdim[38]+f[19]*alpha_vdim[37])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[36]+0.4330127018922193*f[11]*alpha_vdim[35]+0.9682458365518543*(alpha_cdim[3]*f[17]+alpha_cdim[0]*f[8]); 
  out[26] += (0.2581988897471612*f[32]+0.3803194146278324*f[2])*alpha_vdim[67]+(0.276641667586244*f[22]+0.4330127018922193*f[3])*alpha_vdim[57]+(0.276641667586244*f[20]+0.4330127018922193*f[1])*alpha_vdim[55]+0.3872983346207416*(f[19]*alpha_vdim[54]+f[15]*alpha_vdim[50])+(0.276641667586244*f[12]+0.4330127018922193*f[0])*alpha_vdim[47]+0.3872983346207416*(f[7]*alpha_vdim[42]+f[5]*alpha_vdim[40])+(0.8660254037844386*f[29]+0.9682458365518543*f[2])*alpha_cdim[39]+0.4330127018922193*f[22]*alpha_vdim[38]+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alpha_vdim[37]+0.4330127018922193*(f[20]*alpha_vdim[36]+f[12]*alpha_vdim[35])+0.9682458365518543*f[9]*alpha_cdim[35]; 
  out[27] += 0.3872983346207416*(f[22]*alpha_vdim[57]+f[21]*alpha_vdim[56]+f[15]*alpha_vdim[50]+f[7]*alpha_vdim[42]+f[6]*alpha_vdim[41])+(0.3803194146278324*f[33]+0.3872983346207416*f[3])*alpha_vdim[38]+0.4330127018922193*(f[24]*alpha_vdim[37]+f[23]*alpha_vdim[36]+f[13]*alpha_vdim[35])+0.8660254037844386*alpha_vdim[4]*f[30]+0.9682458365518543*(f[22]*alpha_vdim[26]+f[21]*alpha_vdim[25]+alpha_vdim[2]*f[18]+alpha_vdim[1]*f[17]+f[15]*alpha_vdim[16]+alpha_vdim[0]*f[10]+f[7]*alpha_vdim[9]+f[6]*alpha_vdim[8]+f[3]*alpha_vdim[4]); 
  out[28] += 0.8504200642707612*f[25]*alpha_vdim[66]+0.8660254037844386*f[17]*alpha_vdim[56]+0.9682458365518543*f[26]*alpha_vdim[55]+0.8660254037844386*f[16]*alpha_vdim[54]+0.9682458365518543*f[18]*alpha_vdim[50]+0.8660254037844386*f[8]*alpha_vdim[46]+0.9682458365518543*(f[10]*alpha_vdim[41]+f[9]*alpha_vdim[40]+f[17]*alpha_vdim[38]+f[16]*alpha_vdim[37])+0.8660254037844386*f[25]*alpha_vdim[36]+0.9682458365518543*(f[4]*alpha_vdim[36]+f[8]*alpha_vdim[35])+0.4330127018922193*(alpha_cdim[3]*f[30]+alpha_cdim[0]*f[14]); 
  out[29] += 0.8504200642707612*f[26]*alpha_vdim[67]+0.8660254037844386*(f[18]*alpha_vdim[57]+f[16]*alpha_vdim[55])+0.9682458365518543*(f[25]*alpha_vdim[54]+f[17]*alpha_vdim[50])+0.8660254037844386*f[9]*alpha_vdim[47]+0.9682458365518543*(f[10]*alpha_vdim[42]+f[8]*alpha_vdim[40])+(0.3803194146278324*f[34]+0.3872983346207416*f[4])*alpha_cdim[39]+0.9682458365518543*f[18]*alpha_vdim[38]+0.8660254037844386*f[26]*alpha_vdim[37]+0.9682458365518543*(f[4]*alpha_vdim[37]+f[16]*alpha_vdim[36]+f[9]*alpha_vdim[35])+0.4330127018922193*f[14]*alpha_cdim[35]; 
  out[30] += 0.9682458365518543*(f[26]*alpha_vdim[57]+f[25]*alpha_vdim[56]+f[16]*alpha_vdim[50]+f[9]*alpha_vdim[42]+f[8]*alpha_vdim[41])+0.8660254037844386*f[27]*alpha_vdim[38]+0.9682458365518543*(f[4]*alpha_vdim[38]+f[18]*alpha_vdim[37]+f[17]*alpha_vdim[36]+f[10]*alpha_vdim[35])+0.3803194146278324*alpha_vdim[4]*f[34]+0.4330127018922193*(alpha_vdim[2]*f[29]+alpha_vdim[1]*f[28])+0.3872983346207416*(alpha_vdim[26]*f[26]+alpha_vdim[25]*f[25]+alpha_vdim[16]*f[16])+0.4330127018922193*alpha_vdim[0]*f[14]+0.3872983346207416*(alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[4]*f[4]); 
  out[31] += 1.479019945774904*(alpha_cdim[3]*f[21]+alpha_cdim[0]*f[11])+0.6614378277661477*(alpha_cdim[3]*f[3]+alpha_cdim[0]*f[0]); 
  out[32] += (1.479019945774904*f[26]+0.6614378277661477*f[4])*alpha_cdim[39]+(1.479019945774904*f[12]+0.6614378277661477*f[0])*alpha_cdim[35]; 
  out[33] += 0.6614378277661477*(alpha_vdim[32]*f[32]+alpha_vdim[31]*f[31])+1.479019945774904*alpha_vdim[4]*f[27]+0.6614378277661477*(alpha_vdim[26]*f[26]+alpha_vdim[25]*f[25])+1.479019945774904*(alpha_vdim[2]*f[24]+alpha_vdim[1]*f[23])+0.6614378277661477*(alpha_vdim[20]*f[20]+alpha_vdim[19]*f[19]+alpha_vdim[16]*f[16])+1.479019945774904*alpha_vdim[0]*f[13]+0.6614378277661477*(alpha_vdim[12]*f[12]+alpha_vdim[11]*f[11]+alpha_vdim[9]*f[9]+alpha_vdim[8]*f[8]+alpha_vdim[5]*f[5]+alpha_vdim[4]*f[4]+alpha_vdim[2]*f[2]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[34] += 0.6614378277661477*(f[32]*alpha_vdim[67]+f[31]*alpha_vdim[66]+f[22]*alpha_vdim[57]+f[21]*alpha_vdim[56]+f[20]*alpha_vdim[55]+f[19]*alpha_vdim[54]+f[15]*alpha_vdim[50]+f[12]*alpha_vdim[47]+f[11]*alpha_vdim[46]+f[7]*alpha_vdim[42]+f[6]*alpha_vdim[41]+f[5]*alpha_vdim[40])+(1.479019945774904*f[30]+0.6614378277661477*f[3])*alpha_vdim[38]+(1.479019945774904*f[29]+0.6614378277661477*f[2])*alpha_vdim[37]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_vdim[36]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alpha_vdim[35]; 
return alpha_mid; 

} 
