#include <VlasovModDecl.h> 
double VlasovVolElcMag1x1vSerP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  double incr0[4]; 

  for(unsigned int i=0; i<4; ++i){ 

    incr0[i]=0.0; 

  }; 

  const double amid1 = 0.7071067811865475*abar0[0]; 
  incr0[2] = 1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[3] = 1.224744871391589*abar0[0]*f[1]+1.224744871391589*f[0]*abar0[1]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
  out[3] += incr0[3]*dv10; 
return std::abs(amid1)/dv1; 
} 
double VlasovVolElcMag1x1vSerP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  double incr0[8]; 

  for(unsigned int i=0; i<8; ++i){ 

    incr0[i]=0.0; 

  }; 

  const double amid1 = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr0[2] = 1.224744871391589*abar0[2]*f[4]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[3] = abar0[1]*(1.095445115010332*f[4]+1.224744871391589*f[0])+1.095445115010332*f[1]*abar0[2]+1.224744871391589*abar0[0]*f[1]; 
  incr0[5] = 2.738612787525831*abar0[2]*f[6]+2.738612787525831*abar0[1]*f[3]+2.738612787525831*abar0[0]*f[2]; 
  incr0[6] = 1.224744871391589*abar0[0]*f[4]+abar0[2]*(0.7824607964359517*f[4]+1.224744871391589*f[0])+1.095445115010332*abar0[1]*f[1]; 
  incr0[7] = abar0[1]*(2.449489742783178*f[6]+2.738612787525831*f[2])+2.449489742783178*abar0[2]*f[3]+2.738612787525831*abar0[0]*f[3]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
  out[3] += incr0[3]*dv10; 
  out[4] += incr0[4]*dv10; 
  out[5] += incr0[5]*dv10; 
  out[6] += incr0[6]*dv10; 
  out[7] += incr0[7]*dv10; 
return std::abs(amid1)/dv1; 
} 
