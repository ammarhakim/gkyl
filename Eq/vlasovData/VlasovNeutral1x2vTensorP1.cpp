#include <VlasovModDecl.h> 
__host__ __device__ double VlasovNeutralVol1x2vTensorP1(const double *w, const double *dxv, const double *boA, const double *f, double *out) 
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
  const double *Fo1 = &boA[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[8]; 
  double alpha_vdim[16]; 

  alpha_cdim[0] = 5.656854249492382*w0dx0; 
  alpha_cdim[2] = 1.632993161855453*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.0*Fo0[0]*dv10; 
  alpha_vdim[1] = 2.0*Fo0[1]*dv10; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[0]); 

  alpha_vdim[8] = 2.0*Fo1[0]*dv11; 
  alpha_vdim[9] = 2.0*Fo1[1]*dv11; 
  alpha_mid += std::abs(0.1767766952966368*alpha_vdim[8]); 

  out[1] += 0.6123724356957944*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[1]*alpha_vdim[9]+f[0]*alpha_vdim[8]); 
  out[4] += 0.6123724356957944*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[5] += 0.6123724356957944*(f[0]*alpha_vdim[9]+f[1]*alpha_vdim[8]+alpha_cdim[2]*f[6]+alpha_cdim[0]*f[3]); 
  out[6] += 0.6123724356957944*(f[4]*alpha_vdim[9]+f[2]*alpha_vdim[8]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[2]*alpha_vdim[9]+f[4]*alpha_vdim[8]+alpha_cdim[0]*f[6]+alpha_vdim[0]*f[5]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 

  return alpha_mid; 
} 
