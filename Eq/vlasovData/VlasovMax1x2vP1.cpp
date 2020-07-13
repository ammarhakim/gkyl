#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVol1x2vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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
  double alpha_cdim[4]; 
  double alpha_vdim[8]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.0*dv10*(B2[0]*wv2+E0[0]); 
  alpha_vdim[1] = 2.0*dv10*(B2[1]*wv2+E0[1]); 
  alpha_vdim[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[4] = dv11*(2.0*E1[0]-2.0*B2[0]*wv1); 
  alpha_vdim[5] = dv11*(2.0*E1[1]-2.0*B2[1]*wv1); 
  alpha_vdim[6] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[4]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[3]*f[3]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[2]*alpha_vdim[6]+f[1]*alpha_vdim[5]+f[0]*alpha_vdim[4]); 

  return alpha_mid; 
} 
