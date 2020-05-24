#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVol2x2vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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
