#include <VlasovModDecl.h> 
double VlasovVol1x2vSerP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double *B2 = &EM[10]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[8]; 
  double alpha_vdim[16]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]); 
  alpha_vdim[8] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[9] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[10] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[12] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[8]); 
  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[4]*alpha_vdim[12]+f[2]*alpha_vdim[10]+f[1]*alpha_vdim[9]+f[0]*alpha_vdim[8]); 
  out[4] += 0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[2]*alpha_vdim[12]+f[4]*alpha_vdim[10]+f[0]*alpha_vdim[9]+f[1]*alpha_vdim[8]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[1]*alpha_vdim[12]+f[0]*alpha_vdim[10]+f[4]*alpha_vdim[9]+f[2]*alpha_vdim[8]+alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 0.6123724356957944*(f[0]*alpha_vdim[12]+f[1]*alpha_vdim[10]+f[2]*alpha_vdim[9]+f[4]*alpha_vdim[8]+alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
return alpha_mid; 

} 
double VlasovVol1x2vSerP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double *B2 = &EM[15]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[20]; 
  double alpha_vdim[40]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[7] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_vdim[13] = 0.5773502691896258*B2[2]*dv10*dv2; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[7]); 
  alpha_vdim[20] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[21] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[22] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[24] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[27] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_vdim[31] = -0.5773502691896258*B2[2]*dv1*dv11; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[20]-0.1976423537605236*alpha_vdim[27]); 
  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[13]*f[13]+alpha_vdim[7]*f[7]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[11]*alpha_vdim[31]+f[7]*alpha_vdim[27]+f[4]*alpha_vdim[24]+f[2]*alpha_vdim[22]+f[1]*alpha_vdim[21]+f[0]*alpha_vdim[20]); 
  out[4] += 0.5477225575051661*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13]+alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.5477225575051661*(f[4]*alpha_vdim[31]+f[1]*alpha_vdim[27]+f[11]*alpha_vdim[24])+0.6123724356957944*(f[2]*alpha_vdim[24]+f[4]*alpha_vdim[22])+0.5477225575051661*f[7]*alpha_vdim[21]+0.6123724356957944*(f[0]*alpha_vdim[21]+f[1]*alpha_vdim[20]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[7]*alpha_vdim[31]+f[11]*alpha_vdim[27])+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[24]+0.5477225575051661*f[8]*alpha_vdim[22]+0.6123724356957944*(f[0]*alpha_vdim[22]+f[4]*alpha_vdim[21]+f[2]*alpha_vdim[20])+0.5477225575051661*alpha_vdim[5]*f[15]+0.6123724356957944*(alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13])+0.5477225575051661*alpha_vdim[3]*f[9]+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[13]*f[17]+alpha_vdim[7]*f[11]+alpha_vdim[5]*f[10]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[17]*alpha_vdim[31]+f[13]*alpha_vdim[27]+f[10]*alpha_vdim[24]+f[6]*alpha_vdim[22]+f[5]*alpha_vdim[21]+f[3]*alpha_vdim[20]); 
  out[10] += 0.4898979485566357*f[12]*alpha_vdim[31]+0.5477225575051661*(f[1]*alpha_vdim[31]+f[4]*alpha_vdim[27])+(0.5477225575051661*(f[8]+f[7])+0.6123724356957944*f[0])*alpha_vdim[24]+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[22]+0.5477225575051661*f[11]*alpha_vdim[21]+0.6123724356957944*(f[2]*alpha_vdim[21]+f[4]*alpha_vdim[20])+0.4898979485566357*alpha_vdim[13]*f[15]+0.5477225575051661*(alpha_vdim[3]*f[15]+alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[5]*(f[9]+f[7])+f[5]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[11] += 0.3912303982179757*alpha_vdim[13]*f[13]+0.6123724356957944*(alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13])+1.224744871391589*alpha_cdim[2]*f[12]+0.3912303982179757*alpha_vdim[7]*f[7]+0.6123724356957944*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+0.5477225575051661*alpha_vdim[5]*f[5]+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.224744871391589*(alpha_vdim[5]*f[17]+f[10]*alpha_vdim[13]+alpha_vdim[1]*f[11])+1.369306393762915*alpha_vdim[3]*f[10]+0.6123724356957944*alpha_cdim[0]*f[8]+1.224744871391589*f[4]*alpha_vdim[7]+1.369306393762915*(alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += (0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[31]+(0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[27]+0.5477225575051661*f[4]*alpha_vdim[24]+0.6123724356957944*f[11]*alpha_vdim[22]+0.5477225575051661*f[1]*alpha_vdim[21]+0.6123724356957944*f[7]*alpha_vdim[20]+1.369306393762915*(alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += 0.5477225575051661*(f[11]*alpha_vdim[31]+f[4]*alpha_vdim[24]+f[2]*alpha_vdim[22])+0.6123724356957944*(f[12]*alpha_vdim[21]+f[8]*alpha_vdim[20])+1.224744871391589*alpha_vdim[5]*f[19]+1.369306393762915*alpha_vdim[7]*f[17]+1.224744871391589*alpha_vdim[3]*f[16]+1.369306393762915*(f[11]*alpha_vdim[13]+alpha_vdim[1]*f[10]+alpha_vdim[0]*f[6]+f[4]*alpha_vdim[5]+f[2]*alpha_vdim[3]); 
  out[15] += 1.224744871391589*(f[10]*alpha_vdim[31]+f[5]*alpha_vdim[27]+f[17]*alpha_vdim[24])+1.369306393762915*(f[6]*alpha_vdim[24]+f[10]*alpha_vdim[22])+1.224744871391589*f[13]*alpha_vdim[21]+1.369306393762915*(f[3]*alpha_vdim[21]+f[5]*alpha_vdim[20])+0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += 1.369306393762915*(f[13]*alpha_vdim[31]+f[17]*alpha_vdim[27])+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[24]+1.224744871391589*f[14]*alpha_vdim[22]+1.369306393762915*(f[3]*alpha_vdim[22]+f[10]*alpha_vdim[21]+f[6]*alpha_vdim[20])+0.6123724356957944*alpha_vdim[1]*f[15]+0.5477225575051661*alpha_vdim[13]*f[13]+0.6123724356957944*alpha_vdim[0]*f[9]+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[17] += (0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[31]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[27]+(0.4898979485566357*f[12]+0.5477225575051661*f[1])*alpha_vdim[24]+0.6123724356957944*f[7]*alpha_vdim[22]+0.5477225575051661*f[4]*alpha_vdim[21]+0.6123724356957944*f[11]*alpha_vdim[20]+1.224744871391589*alpha_cdim[2]*f[18]+0.4898979485566357*alpha_vdim[5]*f[15]+(0.3912303982179757*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[13]+(0.5477225575051661*f[9]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[13]+1.369306393762915*alpha_cdim[0]*f[10]+0.6123724356957944*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+1.369306393762915*alpha_cdim[2]*f[5]+0.5477225575051661*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]); 
  out[18] += 0.4898979485566357*f[4]*alpha_vdim[31]+0.5477225575051661*f[12]*alpha_vdim[27]+0.4898979485566357*f[11]*alpha_vdim[24]+0.5477225575051661*(f[2]*alpha_vdim[24]+f[4]*alpha_vdim[22])+0.6123724356957944*(f[8]*alpha_vdim[21]+f[12]*alpha_vdim[20])+1.095445115010332*alpha_vdim[13]*f[19]+1.224744871391589*(alpha_vdim[3]*f[19]+alpha_vdim[1]*f[17]+alpha_vdim[5]*f[16])+0.6123724356957944*alpha_cdim[0]*f[14]+1.224744871391589*(f[4]*alpha_vdim[13]+alpha_vdim[5]*f[11])+(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])*f[10]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]); 
  out[19] += 1.095445115010332*f[18]*alpha_vdim[31]+1.224744871391589*(f[5]*alpha_vdim[31]+f[10]*alpha_vdim[27])+(1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*alpha_vdim[24]+(1.224744871391589*f[18]+1.369306393762915*f[5])*alpha_vdim[22]+1.224744871391589*f[17]*alpha_vdim[21]+1.369306393762915*(f[6]*alpha_vdim[21]+f[10]*alpha_vdim[20])+0.6123724356957944*alpha_cdim[0]*f[16]+(0.5477225575051661*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[15]+0.4898979485566357*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+0.6123724356957944*(alpha_cdim[2]+alpha_vdim[1])*f[9]+0.5477225575051661*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 
return alpha_mid; 

} 
double VlasovVol1x2vSerP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double *B2 = &EM[20]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[32]; 
  double alpha_vdim[64]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 
  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[7] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_vdim[13] = 0.5773502691896258*B2[2]*dv10*dv2; 
  alpha_vdim[17] = 2.0*dv10*(B2[3]*wv2+E0[3]); 
  alpha_vdim[25] = 0.5773502691896258*B2[3]*dv10*dv2; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[7]); 
  alpha_vdim[32] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[33] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[34] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[36] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[39] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_vdim[43] = -0.5773502691896258*B2[2]*dv1*dv11; 
  alpha_vdim[49] = dv11*(2.0*E1[3]-2.0*B2[3]*wv1); 
  alpha_vdim[55] = -0.5773502691896258*B2[3]*dv1*dv11; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[32]-0.1976423537605236*alpha_vdim[39]); 
  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[25]*f[25]+alpha_vdim[17]*f[17]+alpha_vdim[13]*f[13]+alpha_vdim[7]*f[7]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[23]*alpha_vdim[55]+f[17]*alpha_vdim[49]+f[11]*alpha_vdim[43]+f[7]*alpha_vdim[39]+f[4]*alpha_vdim[36]+f[2]*alpha_vdim[34]+f[1]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[4] += 0.537852874200477*(alpha_vdim[13]*f[25]+f[13]*alpha_vdim[25]+alpha_vdim[7]*f[17]+f[7]*alpha_vdim[17])+0.5477225575051661*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13]+alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.537852874200477*(f[11]*alpha_vdim[55]+f[7]*alpha_vdim[49])+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alpha_vdim[43]+0.537852874200477*f[17]*alpha_vdim[39]+0.5477225575051661*(f[1]*alpha_vdim[39]+f[11]*alpha_vdim[36])+0.6123724356957944*(f[2]*alpha_vdim[36]+f[4]*alpha_vdim[34])+0.5477225575051661*f[7]*alpha_vdim[33]+0.6123724356957944*(f[0]*alpha_vdim[33]+f[1]*alpha_vdim[32]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[17]*alpha_vdim[55]+f[23]*alpha_vdim[49]+f[7]*alpha_vdim[43]+f[11]*alpha_vdim[39])+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[36]+0.5477225575051661*f[8]*alpha_vdim[34]+0.6123724356957944*(f[0]*alpha_vdim[34]+f[4]*alpha_vdim[33]+f[2]*alpha_vdim[32]+alpha_vdim[17]*f[25]+f[17]*alpha_vdim[25])+0.5477225575051661*alpha_vdim[5]*f[15]+0.6123724356957944*(alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13])+0.5477225575051661*alpha_vdim[3]*f[9]+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[25]*f[29]+alpha_vdim[17]*f[23]+alpha_vdim[13]*f[20]+alpha_vdim[7]*f[11]+alpha_vdim[5]*f[10]+alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[29]*alpha_vdim[55]+f[25]*alpha_vdim[49]+f[20]*alpha_vdim[43]+f[13]*alpha_vdim[39]+f[10]*alpha_vdim[36]+f[6]*alpha_vdim[34]+f[5]*alpha_vdim[33]+f[3]*alpha_vdim[32]); 
  out[10] += 0.537852874200477*(f[7]*alpha_vdim[55]+f[11]*alpha_vdim[49])+(0.537852874200477*f[17]+0.4898979485566357*f[12]+0.5477225575051661*f[1])*alpha_vdim[43]+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alpha_vdim[39]+(0.5477225575051661*(f[8]+f[7])+0.6123724356957944*f[0])*alpha_vdim[36]+(0.5477225575051661*f[12]+0.6123724356957944*f[1])*alpha_vdim[34]+0.5477225575051661*f[11]*alpha_vdim[33]+0.6123724356957944*(f[2]*alpha_vdim[33]+f[4]*alpha_vdim[32])+0.537852874200477*(alpha_vdim[7]*f[25]+f[7]*alpha_vdim[25]+alpha_vdim[13]*f[17]+f[13]*alpha_vdim[17])+0.4898979485566357*alpha_vdim[13]*f[15]+0.5477225575051661*(alpha_vdim[3]*f[15]+alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[5]*(f[9]+f[7])+f[5]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]+f[1]*alpha_vdim[3]); 
  out[11] += 0.3651483716701108*alpha_vdim[25]*f[25]+0.537852874200477*(alpha_vdim[5]*f[25]+f[5]*alpha_vdim[25])+0.3651483716701108*alpha_vdim[17]*f[17]+0.537852874200477*(alpha_vdim[1]*f[17]+f[1]*alpha_vdim[17])+0.3912303982179757*alpha_vdim[13]*f[13]+0.6123724356957944*(alpha_vdim[3]*f[13]+f[3]*alpha_vdim[13])+1.224744871391589*alpha_cdim[2]*f[12]+0.3912303982179757*alpha_vdim[7]*f[7]+0.6123724356957944*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+0.5477225575051661*alpha_vdim[5]*f[5]+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.202675588605909*(alpha_vdim[13]*f[29]+f[20]*alpha_vdim[25]+alpha_vdim[7]*f[23])+1.224744871391589*alpha_vdim[5]*f[20]+0.537852874200477*alpha_cdim[2]*f[18]+1.202675588605909*f[11]*alpha_vdim[17]+1.224744871391589*(f[10]*alpha_vdim[13]+alpha_vdim[1]*f[11])+1.369306393762915*alpha_vdim[3]*f[10]+0.6123724356957944*alpha_cdim[0]*f[8]+1.224744871391589*f[4]*alpha_vdim[7]+1.369306393762915*(alpha_vdim[5]*f[6]+alpha_vdim[0]*f[4])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += (0.3651483716701108*f[23]+0.537852874200477*f[4])*alpha_vdim[55]+(0.3651483716701108*f[17]+0.537852874200477*f[1])*alpha_vdim[49]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[43]+(0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[39]+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alpha_vdim[36]+0.6123724356957944*f[11]*alpha_vdim[34]+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alpha_vdim[33]+0.6123724356957944*f[7]*alpha_vdim[32]+1.369306393762915*(alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += 0.5477225575051661*(f[23]*alpha_vdim[55]+f[11]*alpha_vdim[43])+(0.537852874200477*f[24]+0.5477225575051661*f[4])*alpha_vdim[36]+(0.537852874200477*f[18]+0.5477225575051661*f[2])*alpha_vdim[34]+0.6123724356957944*(f[12]*alpha_vdim[33]+f[8]*alpha_vdim[32])+1.369306393762915*(alpha_vdim[17]*f[29]+f[23]*alpha_vdim[25])+1.224744871391589*alpha_vdim[5]*f[22]+1.369306393762915*alpha_vdim[7]*f[20]+1.224744871391589*alpha_vdim[3]*f[16]+1.369306393762915*(f[11]*alpha_vdim[13]+alpha_vdim[1]*f[10]+alpha_vdim[0]*f[6]+f[4]*alpha_vdim[5]+f[2]*alpha_vdim[3]); 
  out[15] += 1.202675588605909*(f[20]*alpha_vdim[55]+f[13]*alpha_vdim[49])+(1.202675588605909*f[29]+1.224744871391589*f[10])*alpha_vdim[43]+1.202675588605909*f[25]*alpha_vdim[39]+1.224744871391589*(f[5]*alpha_vdim[39]+f[20]*alpha_vdim[36])+1.369306393762915*(f[6]*alpha_vdim[36]+f[10]*alpha_vdim[34])+1.224744871391589*f[13]*alpha_vdim[33]+1.369306393762915*(f[3]*alpha_vdim[33]+f[5]*alpha_vdim[32])+0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += 1.369306393762915*(f[25]*alpha_vdim[55]+f[29]*alpha_vdim[49]+f[13]*alpha_vdim[43]+f[20]*alpha_vdim[39])+(1.224744871391589*f[21]+1.369306393762915*f[5])*alpha_vdim[36]+1.224744871391589*f[14]*alpha_vdim[34]+1.369306393762915*(f[3]*alpha_vdim[34]+f[10]*alpha_vdim[33]+f[6]*alpha_vdim[32])+0.537852874200477*alpha_vdim[5]*f[27]+0.5477225575051661*alpha_vdim[25]*f[25]+0.537852874200477*alpha_vdim[3]*f[19]+0.6123724356957944*alpha_vdim[1]*f[15]+0.5477225575051661*alpha_vdim[13]*f[13]+0.6123724356957944*alpha_vdim[0]*f[9]+0.5477225575051661*(alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]); 
  out[17] += 2.091650066335188*(alpha_cdim[2]*f[11]+alpha_cdim[0]*f[7])+0.9354143466934851*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[18] += 0.9354143466934851*alpha_vdim[25]*f[25]+2.091650066335188*alpha_vdim[5]*f[21]+0.9354143466934851*alpha_vdim[17]*f[17]+2.091650066335188*alpha_vdim[3]*f[14]+0.9354143466934851*alpha_vdim[13]*f[13]+2.091650066335188*(alpha_vdim[1]*f[12]+alpha_vdim[0]*f[8])+0.9354143466934851*(alpha_vdim[7]*f[7]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[19] += 0.9354143466934851*(f[23]*alpha_vdim[55]+f[17]*alpha_vdim[49]+f[11]*alpha_vdim[43]+f[7]*alpha_vdim[39])+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_vdim[36]+(2.091650066335188*f[16]+0.9354143466934851*f[2])*alpha_vdim[34]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alpha_vdim[33]+(2.091650066335188*f[9]+0.9354143466934851*f[0])*alpha_vdim[32]; 
  out[20] += (0.3651483716701108*f[17]+0.4810702354423638*f[12]+0.537852874200477*f[1])*alpha_vdim[55]+(0.3651483716701108*f[23]+0.537852874200477*f[4])*alpha_vdim[49]+(0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[43]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[39]+(0.537852874200477*f[17]+0.4898979485566357*f[12]+0.5477225575051661*f[1])*alpha_vdim[36]+0.6123724356957944*f[7]*alpha_vdim[34]+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alpha_vdim[33]+0.6123724356957944*f[11]*alpha_vdim[32]+(0.3651483716701108*alpha_vdim[17]+0.537852874200477*alpha_vdim[1])*f[25]+(0.3651483716701108*f[17]+0.4810702354423638*f[15]+0.537852874200477*f[1])*alpha_vdim[25]+1.224744871391589*alpha_cdim[2]*f[21]+0.537852874200477*(alpha_vdim[5]*f[17]+f[5]*alpha_vdim[17])+0.4898979485566357*alpha_vdim[5]*f[15]+(0.3912303982179757*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[13]+(0.5477225575051661*f[9]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[13]+1.369306393762915*alpha_cdim[0]*f[10]+0.6123724356957944*(alpha_vdim[3]*f[7]+f[3]*alpha_vdim[7])+1.369306393762915*alpha_cdim[2]*f[5]+0.5477225575051661*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]); 
  out[21] += 0.4810702354423638*f[11]*alpha_vdim[55]+(0.4810702354423638*(f[24]+f[23])+0.4898979485566357*f[4])*alpha_vdim[43]+0.5477225575051661*f[12]*alpha_vdim[39]+(0.537852874200477*f[18]+0.4898979485566357*f[11]+0.5477225575051661*f[2])*alpha_vdim[36]+(0.537852874200477*f[24]+0.5477225575051661*f[4])*alpha_vdim[34]+0.6123724356957944*(f[8]*alpha_vdim[33]+f[12]*alpha_vdim[32])+1.202675588605909*alpha_vdim[7]*f[29]+0.537852874200477*alpha_cdim[2]*f[26]+1.202675588605909*(f[11]*alpha_vdim[25]+alpha_vdim[13]*f[23])+(1.095445115010332*alpha_vdim[13]+1.224744871391589*alpha_vdim[3])*f[22]+1.202675588605909*alpha_vdim[17]*f[20]+1.224744871391589*(alpha_vdim[1]*f[20]+alpha_vdim[5]*f[16])+0.6123724356957944*alpha_cdim[0]*f[14]+1.224744871391589*(f[4]*alpha_vdim[13]+alpha_vdim[5]*f[11])+(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])*f[10]+0.5477225575051661*alpha_cdim[2]*f[6]+1.369306393762915*(alpha_vdim[1]*f[6]+f[2]*alpha_vdim[5]+alpha_vdim[3]*f[4]); 
  out[22] += 1.202675588605909*(f[13]*alpha_vdim[55]+f[20]*alpha_vdim[49])+(1.202675588605909*f[25]+1.095445115010332*f[21]+1.224744871391589*f[5])*alpha_vdim[43]+(1.202675588605909*f[29]+1.224744871391589*f[10])*alpha_vdim[39]+(1.224744871391589*(f[14]+f[13])+1.369306393762915*f[3])*alpha_vdim[36]+(1.224744871391589*f[21]+1.369306393762915*f[5])*alpha_vdim[34]+1.224744871391589*f[20]*alpha_vdim[33]+1.369306393762915*(f[6]*alpha_vdim[33]+f[10]*alpha_vdim[32])+(0.4810702354423638*alpha_vdim[13]+0.537852874200477*alpha_vdim[3])*f[27]+0.4810702354423638*(alpha_vdim[13]*f[25]+f[13]*alpha_vdim[25])+0.537852874200477*alpha_vdim[5]*f[19]+0.6123724356957944*alpha_cdim[0]*f[16]+(0.5477225575051661*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[15]+0.4898979485566357*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+0.6123724356957944*(alpha_cdim[2]+alpha_vdim[1])*f[9]+0.5477225575051661*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]); 
  out[23] += (0.3651483716701108*alpha_vdim[13]+0.6123724356957944*alpha_vdim[3])*f[25]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alpha_vdim[25]+(0.3651483716701108*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_vdim[17]+0.537852874200477*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+2.091650066335188*alpha_cdim[0]*f[11]+alpha_cdim[2]*(0.8366600265340755*f[8]+2.091650066335188*f[7])+0.537852874200477*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.9354143466934851*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[24] += 0.8215838362577489*(alpha_vdim[13]*f[25]+f[13]*alpha_vdim[25])+(1.870828693386971*alpha_vdim[13]+2.091650066335188*alpha_vdim[3])*f[21]+0.6123724356957944*alpha_cdim[0]*f[18]+0.8215838362577489*(alpha_vdim[7]*f[17]+f[7]*alpha_vdim[17])+2.091650066335188*alpha_vdim[5]*f[14]+0.8366600265340755*(alpha_vdim[5]*f[13]+f[5]*alpha_vdim[13])+(1.870828693386971*alpha_vdim[7]+2.091650066335188*alpha_vdim[0])*f[12]+(0.537852874200477*alpha_cdim[2]+2.091650066335188*alpha_vdim[1])*f[8]+0.8366600265340755*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.9354143466934851*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[25] += (0.3651483716701108*f[11]+0.6123724356957944*f[2])*alpha_vdim[55]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_vdim[49]+(0.3651483716701108*f[23]+0.537852874200477*f[4])*alpha_vdim[43]+0.3651483716701108*f[17]*alpha_vdim[39]+0.537852874200477*(f[1]*alpha_vdim[39]+f[11]*alpha_vdim[36])+0.6123724356957944*f[23]*alpha_vdim[34]+0.537852874200477*f[7]*alpha_vdim[33]+0.6123724356957944*f[17]*alpha_vdim[32]+2.091650066335188*(alpha_cdim[2]*f[20]+alpha_cdim[0]*f[13])+0.9354143466934851*(alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[26] += 0.537852874200477*(f[12]*alpha_vdim[36]+f[8]*alpha_vdim[34])+0.6123724356957944*(f[24]*alpha_vdim[33]+f[18]*alpha_vdim[32])+0.9354143466934851*(alpha_vdim[17]*f[25]+f[17]*alpha_vdim[25])+2.091650066335188*alpha_vdim[1]*f[21]+0.8366600265340755*alpha_vdim[5]*f[15]+2.091650066335188*alpha_vdim[0]*f[14]+0.9354143466934851*(alpha_vdim[7]*f[13]+f[7]*alpha_vdim[13])+2.091650066335188*alpha_vdim[5]*f[12]+alpha_vdim[3]*(0.8366600265340755*f[9]+2.091650066335188*f[8])+0.9354143466934851*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[27] += 0.8215838362577489*(f[11]*alpha_vdim[55]+f[7]*alpha_vdim[49])+(0.8215838362577489*f[23]+1.870828693386971*f[22]+0.8366600265340755*f[4])*alpha_vdim[43]+(0.8215838362577489*f[17]+1.870828693386971*f[15]+0.8366600265340755*f[1])*alpha_vdim[39]+(2.091650066335188*f[16]+0.8366600265340755*f[11]+0.9354143466934851*f[2])*alpha_vdim[36]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_vdim[34]+(2.091650066335188*f[9]+0.8366600265340755*f[7]+0.9354143466934851*f[0])*alpha_vdim[33]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alpha_vdim[32]+0.6123724356957944*(alpha_cdim[2]*f[28]+alpha_cdim[0]*f[19]); 
  out[28] += 0.9354143466934851*(f[17]*alpha_vdim[55]+f[23]*alpha_vdim[49]+f[7]*alpha_vdim[43]+f[11]*alpha_vdim[39])+(2.091650066335188*f[15]+0.8366600265340755*f[12]+0.9354143466934851*f[1])*alpha_vdim[36]+(2.091650066335188*f[9]+0.8366600265340755*f[8]+0.9354143466934851*f[0])*alpha_vdim[34]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_vdim[33]+(2.091650066335188*f[16]+0.9354143466934851*f[2])*alpha_vdim[32]+0.6123724356957944*(alpha_vdim[1]*f[27]+alpha_vdim[0]*f[19])+0.537852874200477*(alpha_vdim[5]*f[15]+alpha_vdim[3]*f[9]); 
  out[29] += (0.5477225575051661*f[8]+0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_vdim[55]+(0.3651483716701108*f[11]+0.6123724356957944*f[2])*alpha_vdim[49]+(0.3651483716701108*f[17]+0.4810702354423638*f[12]+0.537852874200477*f[1])*alpha_vdim[43]+0.3651483716701108*f[23]*alpha_vdim[39]+0.537852874200477*(f[4]*alpha_vdim[39]+f[7]*alpha_vdim[36])+0.6123724356957944*f[17]*alpha_vdim[34]+0.537852874200477*f[11]*alpha_vdim[33]+0.6123724356957944*f[23]*alpha_vdim[32]+(0.3651483716701108*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[25]+(0.5477225575051661*f[9]+0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_vdim[25]+2.091650066335188*alpha_cdim[0]*f[20]+(0.3651483716701108*alpha_vdim[13]+0.6123724356957944*alpha_vdim[3])*f[17]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alpha_vdim[17]+0.4810702354423638*alpha_vdim[13]*f[15]+alpha_cdim[2]*(0.8366600265340755*f[14]+2.091650066335188*f[13])+0.537852874200477*(alpha_vdim[1]*f[13]+f[1]*alpha_vdim[13]+alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7])+0.9354143466934851*(alpha_cdim[0]*f[6]+alpha_cdim[2]*f[3]); 
  out[30] += 0.4810702354423638*f[12]*alpha_vdim[43]+0.5477225575051661*f[24]*alpha_vdim[39]+0.537852874200477*(f[8]*alpha_vdim[36]+f[12]*alpha_vdim[34])+0.6123724356957944*(f[18]*alpha_vdim[33]+f[24]*alpha_vdim[32]+alpha_cdim[0]*f[26])+0.8215838362577489*(alpha_vdim[7]*f[25]+f[7]*alpha_vdim[25])+(1.870828693386971*alpha_vdim[7]+2.091650066335188*alpha_vdim[0])*f[21]+0.8215838362577489*(alpha_vdim[13]*f[17]+f[13]*alpha_vdim[17])+(0.7483314773547884*alpha_vdim[13]+0.8366600265340755*alpha_vdim[3])*f[15]+0.537852874200477*alpha_cdim[2]*f[14]+alpha_vdim[1]*(2.091650066335188*f[14]+0.8366600265340755*f[13])+(1.870828693386971*f[12]+0.8366600265340755*f[1])*alpha_vdim[13]+2.091650066335188*alpha_vdim[3]*f[12]+alpha_vdim[5]*(0.8366600265340755*f[9]+2.091650066335188*f[8])+0.8366600265340755*(alpha_vdim[5]*f[7]+f[5]*alpha_vdim[7])+0.9354143466934851*(alpha_vdim[0]*f[5]+f[0]*alpha_vdim[5]+alpha_vdim[1]*f[3]+f[1]*alpha_vdim[3]); 
  out[31] += 0.8215838362577489*(f[7]*alpha_vdim[55]+f[11]*alpha_vdim[49])+(0.8215838362577489*f[17]+1.870828693386971*f[15]+0.7483314773547884*f[12]+0.8366600265340755*f[1])*alpha_vdim[43]+(0.8215838362577489*f[23]+1.870828693386971*f[22]+0.8366600265340755*f[4])*alpha_vdim[39]+(2.091650066335188*f[9]+0.8366600265340755*(f[8]+f[7])+0.9354143466934851*f[0])*alpha_vdim[36]+(2.091650066335188*f[15]+0.8366600265340755*f[12]+0.9354143466934851*f[1])*alpha_vdim[34]+(2.091650066335188*f[16]+0.8366600265340755*f[11]+0.9354143466934851*f[2])*alpha_vdim[33]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_vdim[32]+0.6123724356957944*alpha_cdim[0]*f[28]+0.5477225575051661*alpha_vdim[7]*f[27]+0.6123724356957944*(alpha_vdim[0]*f[27]+(alpha_cdim[2]+alpha_vdim[1])*f[19])+0.4810702354423638*alpha_vdim[13]*f[15]+0.537852874200477*(alpha_vdim[3]*f[15]+alpha_vdim[5]*f[9]); 
return alpha_mid; 

} 
