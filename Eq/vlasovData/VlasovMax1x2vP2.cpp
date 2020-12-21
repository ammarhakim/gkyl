#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVol1x2vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // EM:        Input EM-field.
  // f:         Input distribution function.
  // out:       Incremented output.
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
  double alpha_cdim[10]; 
  double alpha_vdim[20]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_vdim[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha_vdim[7] = 2.0*dv10*(B2[2]*wv2+E0[2]); 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]-0.1976423537605236*alpha_vdim[7]); 

  alpha_vdim[10] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[11] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[12] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_vdim[14] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha_vdim[17] = dv11*(2.0*E1[2]-2.0*B2[2]*wv1); 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[10]-0.1976423537605236*alpha_vdim[17]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[7]*f[7]+alpha_vdim[5]*f[5]+alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alpha_vdim[17]+f[4]*alpha_vdim[14]+f[2]*alpha_vdim[12]+f[1]*alpha_vdim[11]+f[0]*alpha_vdim[10]); 
  out[4] += 0.5477225575051661*(alpha_cdim[2]*f[8]+alpha_vdim[1]*f[7]+f[1]*alpha_vdim[7])+0.6123724356957944*(alpha_vdim[3]*f[5]+f[3]*alpha_vdim[5]+alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.5477225575051661*f[1]*alpha_vdim[17]+0.6123724356957944*(f[2]*alpha_vdim[14]+f[4]*alpha_vdim[12])+0.5477225575051661*f[7]*alpha_vdim[11]+0.6123724356957944*(f[0]*alpha_vdim[11]+f[1]*alpha_vdim[10]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*f[1]*alpha_vdim[14]+0.5477225575051661*f[8]*alpha_vdim[12]+0.6123724356957944*(f[0]*alpha_vdim[12]+f[4]*alpha_vdim[11]+f[2]*alpha_vdim[10])+0.5477225575051661*alpha_vdim[3]*f[9]+0.6123724356957944*(alpha_vdim[1]*f[5]+f[1]*alpha_vdim[5]+alpha_vdim[0]*f[3]+f[0]*alpha_vdim[3]); 
  out[7] += 1.369306393762915*(alpha_cdim[2]*f[4]+alpha_cdim[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha_vdim[3]*f[6]+alpha_vdim[1]*f[4]+alpha_vdim[0]*f[2]); 
  out[9] += 1.369306393762915*(f[6]*alpha_vdim[12]+f[5]*alpha_vdim[11]+f[3]*alpha_vdim[10]); 

  return alpha_mid; 
} 
