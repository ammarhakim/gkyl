#include <VlasovModDecl.h> 
void VlasovVolElc1x3vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &E[0]; 
  out[2] += 1.224744871391589*E0[1]*f[1]*dv10+1.224744871391589*E0[0]*f[0]*dv10; 
} 
void VlasovVolElc1x3vMaxP2(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &E[0]; 
  out[2] += 1.224744871391589*E0[2]*f[11]*dv10+1.224744871391589*E0[1]*f[1]*dv10+1.224744871391589*E0[0]*f[0]*dv10; 
  out[5] += 1.095445115010332*E0[1]*f[11]*dv10+1.095445115010332*f[1]*E0[2]*dv10+1.224744871391589*E0[0]*f[1]*dv10+1.224744871391589*f[0]*E0[1]*dv10; 
  out[7] += 1.224744871391589*E0[1]*f[6]*dv10+1.224744871391589*E0[0]*f[3]*dv10; 
  out[9] += 1.224744871391589*E0[1]*f[8]*dv10+1.224744871391589*E0[0]*f[4]*dv10; 
  out[12] += 2.738612787525831*E0[1]*f[5]*dv10+2.738612787525831*E0[0]*f[2]*dv10; 
} 
