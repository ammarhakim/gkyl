#include <VlasovModDecl.h> 

__host__ __device__ double VlasovPhiVol1x1vSerP1(const double *w, const double *dxv, const double *phi, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // phi:       Input phi-field.
  // EM:        Input (external) EM field vectors.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dv0dx0 = dxv[1]/dxv[0]; 
  const double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2./dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *E0 = &EM[0]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[4]; 
  double alpha_vdim[4]; 

  alpha_cdim[0] = 4.0*w0dx0; 
  alpha_cdim[2] = 1.154700538379252*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = (1.414213562373095*E0[0]-2.449489742783178*phi[1])*dv10; 
  alpha_vdim[1] = 1.414213562373095*E0[1]*dv10; 
  alpha_mid += std::abs(0.25*alpha_vdim[0]); 

  out[1] += 0.8660254037844386*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 

  return alpha_mid; 
} 

__host__ __device__ double VlasovPhiBextVol1x1vSerP1(const double *w, const double *dxv, const double *phi, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // phi:       Input phi-field.
  // EM:        Input (external) EM field vectors.
  // f:         Input distribution function.
  // out:       Incremented output.
  const double dv0dx0 = dxv[1]/dxv[0]; 
  const double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2./dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *E0 = &EM[0]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[4]; 
  double alpha_vdim[4]; 

  alpha_cdim[0] = 4.0*w0dx0; 
  alpha_cdim[2] = 1.154700538379252*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = (1.414213562373095*E0[0]-2.449489742783178*phi[1])*dv10; 
  alpha_vdim[1] = 1.414213562373095*E0[1]*dv10; 
  alpha_mid += std::abs(0.25*alpha_vdim[0]); 

  out[1] += 0.8660254037844386*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 

  return alpha_mid; 
} 

