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
