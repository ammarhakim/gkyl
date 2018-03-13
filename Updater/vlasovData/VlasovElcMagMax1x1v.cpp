#include <VlasovModDecl.h> 
void VlasovVolElcMag1x1vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 

  double incr0[3]; 

  for(unsigned int i=0; i<3; ++i){ 

    incr0[i]=0.0; 

  }; 

  incr0[2] = 1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
} 
void VlasovVolElcMag1x1vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 

  double incr0[6]; 

  for(unsigned int i=0; i<6; ++i){ 

    incr0[i]=0.0; 

  }; 

  incr0[2] = 1.224744871391589*abar0[2]*f[4]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[3] = 1.095445115010332*abar0[1]*f[4]+1.095445115010332*f[1]*abar0[2]+1.224744871391589*abar0[0]*f[1]+1.224744871391589*f[0]*abar0[1]; 
  incr0[5] = 2.738612787525831*abar0[1]*f[3]+2.738612787525831*abar0[0]*f[2]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
  out[3] += incr0[3]*dv10; 
  out[4] += incr0[4]*dv10; 
  out[5] += incr0[5]*dv10; 
} 
