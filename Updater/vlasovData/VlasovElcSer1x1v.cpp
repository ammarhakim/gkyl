#include <VlasovModDecl.h> 
void VlasovVolElc1x1vSerP1(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[1]; 
  out[2] += 1.224744871391589*E[1]*f[1]*dv10+1.224744871391589*E[0]*f[0]*dv10; 
  out[3] += 1.224744871391589*E[0]*f[1]*dv10+1.224744871391589*f[0]*E[1]*dv10; 
} 
void VlasovVolElc1x1vSerP2(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input electric-field/distribution function. out: Incremented output 
  double dv10 = 2/dxv[1]; 
  out[2] += 1.224744871391589*E[2]*f[4]*dv10+1.224744871391589*E[1]*f[1]*dv10+1.224744871391589*E[0]*f[0]*dv10; 
  out[3] += 1.095445115010332*E[1]*f[4]*dv10+1.095445115010332*f[1]*E[2]*dv10+1.224744871391589*E[0]*f[1]*dv10+1.224744871391589*f[0]*E[1]*dv10; 
  out[5] += 2.738612787525831*E[2]*f[6]*dv10+2.738612787525831*E[1]*f[3]*dv10+2.738612787525831*E[0]*f[2]*dv10; 
  out[6] += 0.7824607964359517*E[2]*f[4]*dv10+1.224744871391589*E[0]*f[4]*dv10+1.224744871391589*f[0]*E[2]*dv10+1.095445115010332*E[1]*f[1]*dv10; 
  out[7] += 2.449489742783178*E[1]*f[6]*dv10+2.449489742783178*E[2]*f[3]*dv10+2.738612787525831*E[0]*f[3]*dv10+2.738612787525831*E[1]*f[2]*dv10; 
} 
