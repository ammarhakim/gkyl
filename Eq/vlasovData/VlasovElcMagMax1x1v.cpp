#include <VlasovModDecl.h> 
double VlasovVolElcMag1x1vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double amid1 = 0.7071067811865475*abar0[0]; 
  incr0[2] = 1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
return std::abs(amid1)/dv1; 
} 
double VlasovVolElcMag1x1vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double amid1 = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr0[2] = 1.224744871391589*abar0[2]*f[4]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[3] = 1.095445115010332*abar0[1]*f[4]+1.095445115010332*f[1]*abar0[2]+1.224744871391589*abar0[0]*f[1]+1.224744871391589*f[0]*abar0[1]; 
  incr0[5] = 2.738612787525831*abar0[1]*f[3]+2.738612787525831*abar0[0]*f[2]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
  out[3] += incr0[3]*dv10; 
  out[4] += incr0[4]*dv10; 
  out[5] += incr0[5]*dv10; 
return std::abs(amid1)/dv1; 
} 
double VlasovVolElcMag1x1vMaxP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double abar0[4]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 
  abar0[3] = E0[3]; 

  double incr0[10]; 

  for(unsigned int i=0; i<10; ++i){ 

    incr0[i]=0.0; 

  }; 

  const double amid1 = 0.7071067811865475*abar0[0]-0.7905694150420947*abar0[2]; 
  incr0[2] = 1.224744871391589*abar0[3]*f[8]+1.224744871391589*abar0[2]*f[4]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[3] = 1.075705748400954*abar0[2]*f[8]+1.075705748400954*abar0[3]*f[4]+1.095445115010332*abar0[1]*f[4]+1.095445115010332*f[1]*abar0[2]+1.224744871391589*abar0[0]*f[1]+1.224744871391589*f[0]*abar0[1]; 
  incr0[5] = 2.738612787525831*abar0[2]*f[6]+2.738612787525831*abar0[1]*f[3]+2.738612787525831*abar0[0]*f[2]; 
  incr0[6] = 0.7302967433402215*abar0[3]*f[8]+1.075705748400954*abar0[1]*f[8]+0.7824607964359517*abar0[2]*f[4]+1.224744871391589*abar0[0]*f[4]+1.075705748400954*f[1]*abar0[3]+1.224744871391589*f[0]*abar0[2]+1.095445115010332*abar0[1]*f[1]; 
  incr0[7] = 2.405351177211819*abar0[3]*f[6]+2.449489742783178*abar0[1]*f[6]+2.449489742783178*abar0[2]*f[3]+2.738612787525831*abar0[0]*f[3]+2.738612787525831*abar0[1]*f[2]; 
  incr0[9] = 1.870828693386971*abar0[3]*f[8]+4.183300132670379*abar0[1]*f[7]+4.183300132670378*abar0[0]*f[5]+1.870828693386971*abar0[2]*f[4]+1.870828693386971*abar0[1]*f[1]+1.870828693386971*abar0[0]*f[0]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
  out[3] += incr0[3]*dv10; 
  out[4] += incr0[4]*dv10; 
  out[5] += incr0[5]*dv10; 
  out[6] += incr0[6]*dv10; 
  out[7] += incr0[7]*dv10; 
  out[8] += incr0[8]*dv10; 
  out[9] += incr0[9]*dv10; 
return std::abs(amid1)/dv1; 
} 
double VlasovVolElcMag1x1vMaxP4(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  const double *B0 = &EM[15]; 
  const double *B1 = &EM[20]; 
  const double *B2 = &EM[25]; 

  double abar0[5]; 


  abar0[0] = E0[0]; 
  abar0[1] = E0[1]; 
  abar0[2] = E0[2]; 
  abar0[3] = E0[3]; 
  abar0[4] = E0[4]; 

  double incr0[15]; 

  for(unsigned int i=0; i<15; ++i){ 

    incr0[i]=0.0; 

  }; 

  const double amid1 = 0.7954951288348656*abar0[4]-0.7905694150420947*abar0[2]+0.7071067811865475*abar0[0]; 
  incr0[2] = 1.224744871391589*abar0[4]*f[13]+1.224744871391589*abar0[3]*f[8]+1.224744871391589*abar0[2]*f[4]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[3] = 1.069044967649698*abar0[3]*f[13]+1.069044967649698*abar0[4]*f[8]+1.075705748400954*abar0[2]*f[8]+1.075705748400954*abar0[3]*f[4]+1.095445115010332*abar0[1]*f[4]+1.095445115010332*f[1]*abar0[2]+1.224744871391589*abar0[0]*f[1]+1.224744871391589*f[0]*abar0[1]; 
  incr0[5] = 2.73861278752583*abar0[3]*f[11]+2.738612787525831*abar0[2]*f[6]+2.738612787525831*abar0[1]*f[3]+2.738612787525831*abar0[0]*f[2]; 
  incr0[6] = 0.7113279967599562*abar0[4]*f[13]+1.049781318335648*abar0[2]*f[13]+0.7302967433402215*abar0[3]*f[8]+1.075705748400954*abar0[1]*f[8]+1.049781318335648*abar0[4]*f[4]+0.7824607964359517*abar0[2]*f[4]+1.224744871391589*abar0[0]*f[4]+1.075705748400954*f[1]*abar0[3]+1.224744871391589*f[0]*abar0[2]+1.095445115010332*abar0[1]*f[1]; 
  incr0[7] = 2.390457218668788*abar0[4]*f[11]+2.405351177211819*abar0[2]*f[11]+2.405351177211819*abar0[3]*f[6]+2.449489742783178*abar0[1]*f[6]+2.449489742783178*abar0[2]*f[3]+2.738612787525831*abar0[0]*f[3]+2.738612787525831*abar0[1]*f[2]; 
  incr0[9] = 1.870828693386971*abar0[4]*f[13]+4.183300132670378*abar0[2]*f[10]+1.870828693386971*abar0[3]*f[8]+4.183300132670379*abar0[1]*f[7]+4.183300132670378*abar0[0]*f[5]+1.870828693386971*abar0[2]*f[4]+1.870828693386971*abar0[1]*f[1]+1.870828693386971*abar0[0]*f[0]; 
  incr0[10] = 1.632993161855452*abar0[3]*f[11]+2.405351177211819*abar0[1]*f[11]+2.347382389307855*abar0[4]*f[6]+1.749635530559413*abar0[2]*f[6]+2.738612787525831*abar0[0]*f[6]+2.405351177211819*abar0[3]*f[3]+2.449489742783178*abar0[1]*f[3]+2.738612787525831*abar0[2]*f[2]; 
  incr0[11] = 0.6680426571226847*abar0[3]*f[13]+1.069044967649698*abar0[1]*f[13]+0.6680426571226847*abar0[4]*f[8]+0.7302967433402215*abar0[2]*f[8]+1.224744871391589*abar0[0]*f[8]+0.7302967433402215*abar0[3]*f[4]+1.075705748400954*abar0[1]*f[4]+1.069044967649698*f[1]*abar0[4]+1.224744871391589*f[0]*abar0[3]+1.075705748400954*f[1]*abar0[2]; 
  incr0[12] = 1.632993161855452*abar0[3]*f[13]+3.674234614174766*abar0[3]*f[10]+3.741657386773942*abar0[1]*f[10]+1.632993161855452*abar0[4]*f[8]+1.643167672515498*abar0[2]*f[8]+3.741657386773941*abar0[2]*f[7]+4.183300132670377*abar0[0]*f[7]+4.183300132670378*abar0[1]*f[5]+1.643167672515498*abar0[3]*f[4]+1.673320053068151*abar0[1]*f[4]+1.673320053068151*f[1]*abar0[2]+1.870828693386971*abar0[0]*f[1]+1.870828693386971*f[0]*abar0[1]; 
  incr0[14] = 5.612486080160911*abar0[1]*f[12]+3.674234614174766*abar0[3]*f[11]+5.612486080160912*abar0[0]*f[9]+3.674234614174767*abar0[2]*f[6]+3.674234614174766*abar0[1]*f[3]+3.674234614174766*abar0[0]*f[2]; 

  out[0] += incr0[0]*dv10; 
  out[1] += incr0[1]*dv10; 
  out[2] += incr0[2]*dv10; 
  out[3] += incr0[3]*dv10; 
  out[4] += incr0[4]*dv10; 
  out[5] += incr0[5]*dv10; 
  out[6] += incr0[6]*dv10; 
  out[7] += incr0[7]*dv10; 
  out[8] += incr0[8]*dv10; 
  out[9] += incr0[9]*dv10; 
  out[10] += incr0[10]*dv10; 
  out[11] += incr0[11]*dv10; 
  out[12] += incr0[12]*dv10; 
  out[13] += incr0[13]*dv10; 
  out[14] += incr0[14]*dv10; 
return std::abs(amid1)/dv1; 
} 
