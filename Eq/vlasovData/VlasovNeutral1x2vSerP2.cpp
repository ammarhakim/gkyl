#include <VlasovModDecl.h> 
__host__ __device__ double VlasovNeutralVol1x2vSerP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // EM:        Input EM-field.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *Fo0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *Fo1 = &EM[3]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[20]; 
  double alpha_vdim[40]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.0*Fo0[0]*dv10; 
  alpha_vdim[1] = 2.0*Fo0[1]*dv10; 
  alpha_vdim[7] = 2.0*Fo0[2]*dv10; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[7]); 

  alpha_vdim[20] = 2.0*Fo1[0]*dv11; 
  alpha_vdim[21] = 2.0*Fo1[1]*dv11; 
  alpha_vdim[27] = 2.0*Fo1[2]*dv11; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[20]-0.1976423537605236*alpha_vdim[27]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[7]*f[7]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alpha_vdim[27]+f[1]*alpha_vdim[21]+f[0]*alpha_vdim[20]); 
  out[4] += 0.5477225575051661*(alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.5477225575051661*(f[1]*alpha_vdim[27]+f[7]*alpha_vdim[21])+0.6123724356957944*(f[0]*alpha_vdim[21]+f[1]*alpha_vdim[20]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[11]*alpha_vdim[27]+f[4]*alpha_vdim[21]+f[2]*alpha_vdim[20]+alpha_vdim[7]*f[13]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[7]*f[11]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[13]*alpha_vdim[27]+f[5]*alpha_vdim[21]+f[3]*alpha_vdim[20]); 
  out[10] += 0.5477225575051661*(f[4]*alpha_vdim[27]+f[11]*alpha_vdim[21])+0.6123724356957944*(f[2]*alpha_vdim[21]+f[4]*alpha_vdim[20])+0.5477225575051661*(alpha_cdim[2]*f[14]+alpha_vdim[1]*f[13]+f[5]*alpha_vdim[7])+0.6123724356957944*(alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[11] += 1.224744871391589*alpha_cdim[2]*f[12]+0.3912303982179757*alpha_vdim[7]*f[7]+0.6123724356957944*(alpha_vdim[0]*f[7]+f[0]*alpha_vdim[7])+1.369306393762915*alpha_cdim[0]*f[4]+f[1]*(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1]); 
  out[12] += 1.224744871391589*alpha_vdim[1]*f[11]+0.6123724356957944*alpha_cdim[0]*f[8]+f[4]*(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[2]; 
  out[13] += (0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha_vdim[27]+0.5477225575051661*f[1]*alpha_vdim[21]+0.6123724356957944*f[7]*alpha_vdim[20]+1.369306393762915*(alpha_cdim[2]*f[10]+alpha_cdim[0]*f[5]); 
  out[14] += 0.6123724356957944*(f[12]*alpha_vdim[21]+f[8]*alpha_vdim[20])+1.369306393762915*(alpha_vdim[7]*f[17]+alpha_vdim[1]*f[10]+alpha_vdim[0]*f[6]); 
  out[15] += 1.224744871391589*(f[5]*alpha_vdim[27]+f[13]*alpha_vdim[21])+1.369306393762915*(f[3]*alpha_vdim[21]+f[5]*alpha_vdim[20])+0.6123724356957944*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[9]); 
  out[16] += 1.369306393762915*(f[17]*alpha_vdim[27]+f[10]*alpha_vdim[21]+f[6]*alpha_vdim[20])+0.6123724356957944*(alpha_vdim[1]*f[15]+alpha_vdim[0]*f[9]); 
  out[17] += (0.3912303982179757*f[11]+0.6123724356957944*f[2])*alpha_vdim[27]+0.5477225575051661*f[4]*alpha_vdim[21]+0.6123724356957944*f[11]*alpha_vdim[20]+1.224744871391589*alpha_cdim[2]*f[18]+(0.3912303982179757*alpha_vdim[7]+0.6123724356957944*alpha_vdim[0])*f[13]+1.369306393762915*alpha_cdim[0]*f[10]+0.6123724356957944*f[3]*alpha_vdim[7]+(1.369306393762915*alpha_cdim[2]+0.5477225575051661*alpha_vdim[1])*f[5]; 
  out[18] += 0.5477225575051661*f[12]*alpha_vdim[27]+0.6123724356957944*(f[8]*alpha_vdim[21]+f[12]*alpha_vdim[20])+1.224744871391589*alpha_vdim[1]*f[17]+0.6123724356957944*alpha_cdim[0]*f[14]+(1.224744871391589*alpha_vdim[7]+1.369306393762915*alpha_vdim[0])*f[10]+(0.5477225575051661*alpha_cdim[2]+1.369306393762915*alpha_vdim[1])*f[6]; 
  out[19] += 1.224744871391589*(f[10]*alpha_vdim[27]+f[17]*alpha_vdim[21])+1.369306393762915*(f[6]*alpha_vdim[21]+f[10]*alpha_vdim[20])+0.6123724356957944*alpha_cdim[0]*f[16]+0.5477225575051661*alpha_vdim[7]*f[15]+0.6123724356957944*(alpha_vdim[0]*f[15]+(alpha_cdim[2]+alpha_vdim[1])*f[9]); 

  return alpha_mid; 
} 
