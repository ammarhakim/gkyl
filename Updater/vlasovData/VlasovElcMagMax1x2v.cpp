#include <VlasovModDecl.h> 
void VlasovVolElcMag1x2vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double abar0[2]; 

  double abar1[2]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar1[0] = E1[0]-wv1*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 
  abar1[1] = E1[1]-wv1*B2[1]; 

  double incr0[4]; 

  double incr1[4]; 

  for(unsigned int i=0; i<4; ++i){ 

    incr0[i]=0.0; 

    incr1[i]=0.0; 

  }; 

  incr0[2] = 0.3535533905932737*B2[0]*f[3]*dv2+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr1[3] = (-0.3535533905932737*B2[0]*f[2]*dv1)+1.224744871391589*abar1[1]*f[1]+1.224744871391589*abar1[0]*f[0]; 

  out[0] += incr1[0]*dv11+incr0[0]*dv10; 
  out[1] += incr1[1]*dv11+incr0[1]*dv10; 
  out[2] += incr1[2]*dv11+incr0[2]*dv10; 
  out[3] += incr1[3]*dv11+incr0[3]*dv10; 
} 
void VlasovVolElcMag1x2vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[3]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double abar0[3]; 

  double abar1[3]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar1[0] = E1[0]-wv1*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 
  abar1[1] = E1[1]-wv1*B2[1]; 
  abar0[2] = E0[2]+wv2*B2[2]; 
  abar1[2] = E1[2]-wv1*B2[2]; 

  double incr0[10]; 

  double incr1[10]; 

  for(unsigned int i=0; i<10; ++i){ 

    incr0[i]=0.0; 

    incr1[i]=0.0; 

  }; 

  incr0[2] = 0.3535533905932737*B2[1]*f[5]*dv2+0.3535533905932737*B2[0]*f[3]*dv2+1.224744871391589*abar0[2]*f[7]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[4] = 0.3162277660168379*B2[2]*f[5]*dv2+0.3535533905932737*B2[0]*f[5]*dv2+0.3535533905932737*B2[1]*f[3]*dv2+1.095445115010332*abar0[1]*f[7]+1.095445115010332*f[1]*abar0[2]+1.224744871391589*abar0[0]*f[1]+1.224744871391589*f[0]*abar0[1]; 
  incr0[6] = 0.3162277660168379*B2[0]*f[9]*dv2+0.3535533905932737*B2[2]*f[7]*dv2+0.3535533905932737*f[1]*B2[1]*dv2+0.3535533905932737*f[0]*B2[0]*dv2+1.224744871391589*abar0[1]*f[5]+1.224744871391589*abar0[0]*f[3]; 
  incr0[8] = 0.7905694150420947*B2[0]*f[6]*dv2+2.738612787525831*abar0[1]*f[4]+2.738612787525831*abar0[0]*f[2]; 
  incr1[3] = (-0.3535533905932737*B2[1]*f[4]*dv1)-0.3535533905932737*B2[0]*f[2]*dv1+1.224744871391589*abar1[2]*f[7]+1.224744871391589*abar1[1]*f[1]+1.224744871391589*abar1[0]*f[0]; 
  incr1[5] = (-0.3162277660168379*B2[2]*f[4]*dv1)-0.3535533905932737*B2[0]*f[4]*dv1-0.3535533905932737*B2[1]*f[2]*dv1+1.095445115010332*abar1[1]*f[7]+1.095445115010332*f[1]*abar1[2]+1.224744871391589*abar1[0]*f[1]+1.224744871391589*f[0]*abar1[1]; 
  incr1[6] = (-0.3162277660168379*B2[0]*f[8]*dv1)-0.3535533905932737*B2[2]*f[7]*dv1-0.3535533905932737*f[1]*B2[1]*dv1-0.3535533905932737*f[0]*B2[0]*dv1+1.224744871391589*abar1[1]*f[4]+1.224744871391589*abar1[0]*f[2]; 
  incr1[9] = (-0.7905694150420947*B2[0]*f[6]*dv1)+2.738612787525831*abar1[1]*f[5]+2.738612787525831*abar1[0]*f[3]; 

  out[0] += incr1[0]*dv11+incr0[0]*dv10; 
  out[1] += incr1[1]*dv11+incr0[1]*dv10; 
  out[2] += incr1[2]*dv11+incr0[2]*dv10; 
  out[3] += incr1[3]*dv11+incr0[3]*dv10; 
  out[4] += incr1[4]*dv11+incr0[4]*dv10; 
  out[5] += incr1[5]*dv11+incr0[5]*dv10; 
  out[6] += incr1[6]*dv11+incr0[6]*dv10; 
  out[7] += incr1[7]*dv11+incr0[7]*dv10; 
  out[8] += incr1[8]*dv11+incr0[8]*dv10; 
  out[9] += incr1[9]*dv11+incr0[9]*dv10; 
} 
