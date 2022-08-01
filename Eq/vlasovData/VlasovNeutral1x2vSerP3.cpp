#include <VlasovModDecl.h> 
__host__ __device__ double VlasovNeutralVol1x2vSerP3(const double *w, const double *dxv, const double *boA, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // boA:       Input body acceleration.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *Fo0 = &boA[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *Fo1 = &boA[4]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[32]; 
  double alpha_vdim[64]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.0*Fo0[0]*dv10; 
  alpha_vdim[1] = 2.0*Fo0[1]*dv10; 
  alpha_vdim[7] = 2.0*Fo0[2]*dv10; 
  alpha_vdim[17] = 2.0*Fo0[3]*dv10; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[7]); 

  alpha_vdim[32] = 2.0*Fo1[0]*dv11; 
  alpha_vdim[33] = 2.0*Fo1[1]*dv11; 
  alpha_vdim[39] = 2.0*Fo1[2]*dv11; 
  alpha_vdim[49] = 2.0*Fo1[3]*dv11; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[32]-0.1976423537605236*alpha_vdim[39]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[17]*f[17]+alpha_vdim[7]*f[7]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[17]*alpha_vdim[49]+f[7]*alpha_vdim[39]+f[1]*alpha_vdim[33]+f[0]*alpha_vdim[32]); 
  out[4] += 0.537852874200477*(alpha_vdim[7]*f[17]+f[7]*alpha_vdim[17])+0.5477225575051661*(alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.537852874200477*(f[7]*alpha_vdim[49]+f[17]*alpha_vdim[39])+0.5477225575051661*(f[1]*alpha_vdim[39]+f[7]*alpha_vdim[33])+0.6123724356957944*(f[0]*alpha_vdim[33]+f[1]*alpha_vdim[32]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[23]*alpha_vdim[49]+f[11]*alpha_vdim[39]+f[4]*alpha_vdim[33]+f[2]*alpha_vdim[32]+alpha_vdim[17]*f[25]+alpha_vdim[7]*f[13]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[17]*f[23]+alpha_vdim[7]*f[11]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[25]*alpha_vdim[49]+f[13]*alpha_vdim[39]+f[5]*alpha_vdim[33]+f[3]*alpha_vdim[32]); 
  out[10] += 0.537852874200477*(f[11]*alpha_vdim[49]+f[23]*alpha_vdim[39])+0.5477225575051661*(f[4]*alpha_vdim[39]+f[11]*alpha_vdim[33])+0.6123724356957944*(f[2]*alpha_vdim[33]+f[4]*alpha_vdim[32])+0.537852874200477*(alpha_vdim[7]*f[25]+f[13]*alpha_vdim[17])+0.5477225575051661*(alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+f[5]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[11] += 0.3651483716701108*alpha_vdim[17]*f[17]+0.537852874200477*(alpha_vdim[1]*f[17]+f[1]*alpha_vdim[17])+1.224744871391589*alpha_cdim[2]*f[12]+0.3912303982179757*alpha_vdim[7]*f[7]+0.6123724356957944*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.202675588605909*alpha_vdim[7]*f[23]+0.537852874200477*alpha_cdim[2]*f[18]+f[11]*(1.202675588605909*alpha_vdim[17]+1.224744871391589*alpha_vdim[1])+0.6123724356957944*alpha_cdim[0]*f[8]+f[4]*(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += (0.3651483716701108*f[17]+0.537852874200477*f[1])*alpha_vdim[49]+(0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[39]+(0.537852874200477*f[17]+0.5477225575051661*f[1])*alpha_vdim[33]+0.6123724356957944*f[7]*alpha_vdim[32]+1.369306393762915*(alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += 0.6123724356957944*(f[12]*alpha_vdim[33]+f[8]*alpha_vdim[32])+1.369306393762915*(alpha_vdim[17]*f[29]+alpha_vdim[7]*f[20]+alpha_vdim[1]*f[10]+alpha_vdim[0]*f[6]); 
  out[15] += 1.202675588605909*(f[13]*alpha_vdim[49]+f[25]*alpha_vdim[39])+1.224744871391589*(f[5]*alpha_vdim[39]+f[13]*alpha_vdim[33])+1.369306393762915*(f[3]*alpha_vdim[33]+f[5]*alpha_vdim[32])+0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += 1.369306393762915*(f[29]*alpha_vdim[49]+f[20]*alpha_vdim[39]+f[10]*alpha_vdim[33]+f[6]*alpha_vdim[32])+0.6123724356957944*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[9]); 
  out[17] += 2.091650066335188*(alpha_cdim[2]*f[11]+alpha_cdim[0]*f[7])+0.9354143466934851*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[18] += 0.9354143466934851*alpha_vdim[17]*f[17]+2.091650066335188*(alpha_vdim[1]*f[12]+alpha_vdim[0]*f[8])+0.9354143466934851*(alpha_vdim[7]*f[7]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[19] += 0.9354143466934851*(f[17]*alpha_vdim[49]+f[7]*alpha_vdim[39])+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alpha_vdim[33]+(2.091650066335188*f[9]+0.9354143466934851*f[0])*alpha_vdim[32]; 
  out[20] += (0.3651483716701108*f[23]+0.537852874200477*f[4])*alpha_vdim[49]+(0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[39]+(0.537852874200477*f[23]+0.5477225575051661*f[4])*alpha_vdim[33]+0.6123724356957944*f[11]*alpha_vdim[32]+(0.3651483716701108*alpha_vdim[17]+0.537852874200477*alpha_vdim[1])*f[25]+1.224744871391589*alpha_cdim[2]*f[21]+0.537852874200477*f[5]*alpha_vdim[17]+(0.3912303982179757*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[13]+1.369306393762915*alpha_cdim[0]*f[10]+0.6123724356957944*f[3]*alpha_vdim[7]+(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1])*f[5]; 
  out[21] += 0.5477225575051661*f[12]*alpha_vdim[39]+0.6123724356957944*(f[8]*alpha_vdim[33]+f[12]*alpha_vdim[32])+1.202675588605909*alpha_vdim[7]*f[29]+0.537852874200477*alpha_cdim[2]*f[26]+(1.202675588605909*alpha_vdim[17]+1.224744871391589*alpha_vdim[1])*f[20]+0.6123724356957944*alpha_cdim[0]*f[14]+(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])*f[10]+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[6]; 
  out[22] += 1.202675588605909*(f[20]*alpha_vdim[49]+f[29]*alpha_vdim[39])+1.224744871391589*(f[10]*alpha_vdim[39]+f[20]*alpha_vdim[33])+1.369306393762915*(f[6]*alpha_vdim[33]+f[10]*alpha_vdim[32])+0.6123724356957944*alpha_cdim[0]*f[16]+0.5477225575051661*alpha_vdim[7]*f[15]+0.6123724356957944*(alpha_vdim[0]*f[15]+(alpha_cdim[2]+alpha_vdim[1])*f[9]); 
  out[23] += (0.3651483716701108*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_vdim[17]+2.091650066335188*alpha_cdim[0]*f[11]+alpha_cdim[2]*(0.8366600265340755*f[8]+2.091650066335188*f[7])+0.537852874200477*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.9354143466934851*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[24] += 0.6123724356957944*alpha_cdim[0]*f[18]+0.8215838362577489*(alpha_vdim[7]*f[17]+f[7]*alpha_vdim[17])+(1.870828693386971*alpha_vdim[7]+2.091650066335188*alpha_vdim[0])*f[12]+(0.537852874200477*alpha_cdim[2]+2.091650066335188*alpha_vdim[1])*f[8]+0.8366600265340755*(alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.9354143466934851*(alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[25] += (0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha_vdim[49]+0.3651483716701108*f[17]*alpha_vdim[39]+0.537852874200477*(f[1]*alpha_vdim[39]+f[7]*alpha_vdim[33])+0.6123724356957944*f[17]*alpha_vdim[32]+2.091650066335188*(alpha_cdim[2]*f[20]+alpha_cdim[0]*f[13])+0.9354143466934851*(alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[26] += 0.6123724356957944*(f[24]*alpha_vdim[33]+f[18]*alpha_vdim[32])+0.9354143466934851*alpha_vdim[17]*f[25]+2.091650066335188*(alpha_vdim[1]*f[21]+alpha_vdim[0]*f[14])+0.9354143466934851*(alpha_vdim[7]*f[13]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[3]); 
  out[27] += 0.8215838362577489*f[7]*alpha_vdim[49]+(0.8215838362577489*f[17]+1.870828693386971*f[15]+0.8366600265340755*f[1])*alpha_vdim[39]+(2.091650066335188*f[9]+0.8366600265340755*f[7]+0.9354143466934851*f[0])*alpha_vdim[33]+(2.091650066335188*f[15]+0.9354143466934851*f[1])*alpha_vdim[32]+0.6123724356957944*(alpha_cdim[2]*f[28]+alpha_cdim[0]*f[19]); 
  out[28] += 0.9354143466934851*(f[23]*alpha_vdim[49]+f[11]*alpha_vdim[39])+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_vdim[33]+(2.091650066335188*f[16]+0.9354143466934851*f[2])*alpha_vdim[32]+0.6123724356957944*(alpha_vdim[1]*f[27]+alpha_vdim[0]*f[19]); 
  out[29] += (0.3651483716701108*f[11]+0.6123724356957944*f[2])*alpha_vdim[49]+0.3651483716701108*f[23]*alpha_vdim[39]+0.537852874200477*(f[4]*alpha_vdim[39]+f[11]*alpha_vdim[33])+0.6123724356957944*f[23]*alpha_vdim[32]+(0.3651483716701108*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[25]+2.091650066335188*alpha_cdim[0]*f[20]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alpha_vdim[17]+alpha_cdim[2]*(0.8366600265340755*f[14]+2.091650066335188*f[13])+0.537852874200477*(alpha_vdim[1]*f[13]+f[5]*alpha_vdim[7])+0.9354143466934851*(alpha_cdim[0]*f[6]+alpha_cdim[2]*f[3]); 
  out[30] += 0.5477225575051661*f[24]*alpha_vdim[39]+0.6123724356957944*(f[18]*alpha_vdim[33]+f[24]*alpha_vdim[32]+alpha_cdim[0]*f[26])+0.8215838362577489*alpha_vdim[7]*f[25]+(1.870828693386971*alpha_vdim[7]+2.091650066335188*alpha_vdim[0])*f[21]+0.8215838362577489*f[13]*alpha_vdim[17]+(0.537852874200477*alpha_cdim[2]+2.091650066335188*alpha_vdim[1])*f[14]+0.8366600265340755*(alpha_vdim[1]*f[13]+f[5]*alpha_vdim[7])+0.9354143466934851*(alpha_vdim[0]*f[5]+alpha_vdim[1]*f[3]); 
  out[31] += 0.8215838362577489*f[11]*alpha_vdim[49]+(0.8215838362577489*f[23]+1.870828693386971*f[22]+0.8366600265340755*f[4])*alpha_vdim[39]+(2.091650066335188*f[16]+0.8366600265340755*f[11]+0.9354143466934851*f[2])*alpha_vdim[33]+(2.091650066335188*f[22]+0.9354143466934851*f[4])*alpha_vdim[32]+0.6123724356957944*alpha_cdim[0]*f[28]+0.5477225575051661*alpha_vdim[7]*f[27]+0.6123724356957944*(alpha_vdim[0]*f[27]+(alpha_cdim[2]+alpha_vdim[1])*f[19]); 

  return alpha_mid; 
} 
