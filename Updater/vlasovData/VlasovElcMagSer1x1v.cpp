#include <VlasovModDecl.h> 
void VlasovVolElcMag1x1vSerP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  out[2] += 1.224744871391589*f[1]*E0[1]*dv10+1.224744871391589*f[0]*E0[0]*dv10; 
  out[3] += 1.224744871391589*f[0]*E0[1]*dv10+1.224744871391589*E0[0]*f[1]*dv10; 
} 
void VlasovVolElcMag1x1vSerP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  out[2] += 1.224744871391589*E0[2]*f[4]*dv10+1.224744871391589*f[1]*E0[1]*dv10+1.224744871391589*f[0]*E0[0]*dv10; 
  out[3] += 1.095445115010332*E0[1]*f[4]*dv10+1.095445115010332*f[1]*E0[2]*dv10+1.224744871391589*f[0]*E0[1]*dv10+1.224744871391589*E0[0]*f[1]*dv10; 
  out[5] += 2.738612787525831*E0[2]*f[6]*dv10+2.738612787525831*E0[1]*f[3]*dv10+2.738612787525831*E0[0]*f[2]*dv10; 
  out[6] += 0.7824607964359517*E0[2]*f[4]*dv10+1.224744871391589*E0[0]*f[4]*dv10+1.224744871391589*f[0]*E0[2]*dv10+1.095445115010332*f[1]*E0[1]*dv10; 
  out[7] += 2.449489742783178*E0[1]*f[6]*dv10+2.449489742783178*E0[2]*f[3]*dv10+2.738612787525831*E0[0]*f[3]*dv10+2.738612787525831*E0[1]*f[2]*dv10; 
} 
