#include <VlasovModDecl.h> 
double VlasovVolElcMag1x2vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double amid1 = 0.3535533905932737*B2[0]*dv2+0.7071067811865475*abar0[0]; 
  incr0[2] = 0.3535533905932737*B2[0]*f[3]*dv2+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  const double amid2 = 0.7071067811865475*abar1[0]-0.3535533905932737*B2[0]*dv1; 
  incr1[3] = (-0.3535533905932737*B2[0]*f[2]*dv1)+1.224744871391589*abar1[1]*f[1]+1.224744871391589*abar1[0]*f[0]; 

  out[0] += incr1[0]*dv11+incr0[0]*dv10; 
  out[1] += incr1[1]*dv11+incr0[1]*dv10; 
  out[2] += incr1[2]*dv11+incr0[2]*dv10; 
  out[3] += incr1[3]*dv11+incr0[3]*dv10; 
return std::abs(amid1)/dv1+std::abs(amid2)/dv2; 
} 
double VlasovVolElcMag1x2vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
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

  const double amid1 = (-0.3952847075210473*B2[2]*dv2)+0.3535533905932737*B2[0]*dv2-0.7905694150420947*abar0[2]+0.7071067811865475*abar0[0]; 
  incr0[2] = (0.3535533905932737*B2[1]*f[5]+0.3535533905932737*B2[0]*f[3])*dv2+1.224744871391589*abar0[2]*f[7]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[4] = (0.3162277660168379*B2[2]*f[5]+0.3535533905932737*B2[0]*f[5]+0.3535533905932737*B2[1]*f[3])*dv2+abar0[1]*(1.095445115010332*f[7]+1.224744871391589*f[0])+1.095445115010332*f[1]*abar0[2]+1.224744871391589*abar0[0]*f[1]; 
  incr0[6] = (0.3162277660168379*B2[0]*f[9]+0.3535533905932737*B2[2]*f[7]+0.3535533905932737*f[1]*B2[1]+0.3535533905932737*f[0]*B2[0])*dv2+1.224744871391589*abar0[1]*f[5]+1.224744871391589*abar0[0]*f[3]; 
  incr0[8] = 0.7905694150420947*B2[0]*f[6]*dv2+2.738612787525831*abar0[1]*f[4]+2.738612787525831*abar0[0]*f[2]; 
  const double amid2 = 0.3952847075210473*B2[2]*dv1-0.3535533905932737*B2[0]*dv1-0.7905694150420947*abar1[2]+0.7071067811865475*abar1[0]; 
  incr1[3] = ((-0.3535533905932737*B2[1]*f[4])-0.3535533905932737*B2[0]*f[2])*dv1+1.224744871391589*abar1[2]*f[7]+1.224744871391589*abar1[1]*f[1]+1.224744871391589*abar1[0]*f[0]; 
  incr1[5] = ((-0.3162277660168379*B2[2]*f[4])-0.3535533905932737*B2[0]*f[4]-0.3535533905932737*B2[1]*f[2])*dv1+abar1[1]*(1.095445115010332*f[7]+1.224744871391589*f[0])+1.095445115010332*f[1]*abar1[2]+1.224744871391589*abar1[0]*f[1]; 
  incr1[6] = ((-0.3162277660168379*B2[0]*f[8])-0.3535533905932737*B2[2]*f[7]-0.3535533905932737*f[1]*B2[1]-0.3535533905932737*f[0]*B2[0])*dv1+1.224744871391589*abar1[1]*f[4]+1.224744871391589*abar1[0]*f[2]; 
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
return std::abs(amid1)/dv1+std::abs(amid2)/dv2; 
} 
double VlasovVolElcMag1x2vMaxP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[4]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double abar0[4]; 

  double abar1[4]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar1[0] = E1[0]-wv1*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 
  abar1[1] = E1[1]-wv1*B2[1]; 
  abar0[2] = E0[2]+wv2*B2[2]; 
  abar1[2] = E1[2]-wv1*B2[2]; 
  abar0[3] = E0[3]+wv2*B2[3]; 
  abar1[3] = E1[3]-wv1*B2[3]; 

  double incr0[20]; 

  double incr1[20]; 

  for(unsigned int i=0; i<20; ++i){ 

    incr0[i]=0.0; 

    incr1[i]=0.0; 

  }; 

  const double amid1 = (-0.3952847075210473*B2[2]*dv2)+0.3535533905932737*B2[0]*dv2-0.7905694150420947*abar0[2]+0.7071067811865475*abar0[0]; 
  incr0[2] = (0.3535533905932737*B2[2]*f[13]+0.3535533905932737*B2[1]*f[5]+0.3535533905932737*B2[0]*f[3])*dv2+1.224744871391589*abar0[3]*f[17]+1.224744871391589*abar0[2]*f[7]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[4] = (0.3105295017040593*B2[3]*f[13]+0.3162277660168379*B2[1]*f[13]+0.3162277660168379*B2[2]*f[5]+0.3535533905932737*B2[0]*f[5]+0.3535533905932737*B2[1]*f[3])*dv2+abar0[2]*(1.075705748400954*f[17]+1.095445115010332*f[1])+abar0[1]*(1.095445115010332*f[7]+1.224744871391589*f[0])+1.075705748400954*abar0[3]*f[7]+1.224744871391589*abar0[0]*f[1]; 
  incr0[6] = (0.3535533905932737*B2[3]*f[17]+0.3162277660168379*B2[1]*f[15]+0.3162277660168379*B2[0]*f[9]+0.3535533905932737*B2[2]*f[7]+0.3535533905932737*f[1]*B2[1]+0.3535533905932737*f[0]*B2[0])*dv2+1.224744871391589*abar0[2]*f[13]+1.224744871391589*abar0[1]*f[5]+1.224744871391589*abar0[0]*f[3]; 
  incr0[8] = (0.7905694150420947*B2[1]*f[10]+0.7905694150420947*B2[0]*f[6])*dv2+2.738612787525831*abar0[2]*f[11]+2.738612787525831*abar0[1]*f[4]+2.738612787525831*abar0[0]*f[2]; 
  incr0[10] = (0.3105295017040592*B2[2]*f[17]+0.282842712474619*B2[2]*f[15]+0.3162277660168379*B2[0]*f[15]+0.3162277660168379*B2[1]*f[9]+0.3105295017040592*B2[3]*f[7]+0.3162277660168379*B2[1]*f[7]+0.3162277660168379*f[1]*B2[2]+0.3535533905932737*f[0]*B2[1]+0.3535533905932737*B2[0]*f[1])*dv2+abar0[1]*(1.095445115010332*f[13]+1.224744871391589*f[3])+1.075705748400954*abar0[3]*f[13]+1.095445115010332*abar0[2]*f[5]+1.224744871391589*abar0[0]*f[5]; 
  incr0[11] = (0.2258769757263128*B2[2]*f[13]+0.3535533905932737*B2[0]*f[13]+0.3105295017040593*B2[3]*f[5]+0.3162277660168379*B2[1]*f[5]+0.3535533905932737*B2[2]*f[3])*dv2+abar0[1]*(1.075705748400954*f[17]+1.095445115010332*f[1])+abar0[3]*(0.7302967433402215*f[17]+1.075705748400954*f[1])+1.224744871391589*abar0[0]*f[7]+abar0[2]*(0.7824607964359517*f[7]+1.224744871391589*f[0]); 
  incr0[12] = (0.7071067811865475*B2[2]*f[10]+0.7905694150420948*B2[0]*f[10]+0.7905694150420948*B2[1]*f[6])*dv2+abar0[1]*(2.449489742783178*f[11]+2.738612787525831*f[2])+2.405351177211819*abar0[3]*f[11]+2.449489742783178*abar0[2]*f[4]+2.738612787525831*abar0[0]*f[4]; 
  incr0[14] = (0.7071067811865475*B2[0]*f[16]+0.7905694150420947*B2[2]*f[11]+0.7905694150420948*B2[1]*f[4]+0.7905694150420948*B2[0]*f[2])*dv2+2.738612787525831*abar0[1]*f[10]+2.738612787525831*abar0[0]*f[6]; 
  incr0[16] = (0.3105295017040593*B2[0]*f[19]+0.3162277660168379*B2[2]*f[13]+0.3162277660168379*B2[1]*f[5]+0.3162277660168379*B2[0]*f[3])*dv2+1.224744871391589*abar0[1]*f[15]+1.224744871391589*abar0[0]*f[9]; 
  incr0[18] = (1.20761472884912*B2[0]*f[14]+0.5400617248673215*B2[2]*f[13]+0.5400617248673216*B2[1]*f[5]+0.5400617248673216*B2[0]*f[3])*dv2+1.870828693386971*abar0[3]*f[17]+abar0[1]*(4.183300132670379*f[12]+1.870828693386971*f[1])+abar0[0]*(4.183300132670378*f[8]+1.870828693386971*f[0])+1.870828693386971*abar0[2]*f[7]; 
  const double amid2 = 0.3952847075210473*B2[2]*dv1-0.3535533905932737*B2[0]*dv1-0.7905694150420947*abar1[2]+0.7071067811865475*abar1[0]; 
  incr1[3] = ((-0.3535533905932737*B2[2]*f[11])-0.3535533905932737*B2[1]*f[4]-0.3535533905932737*B2[0]*f[2])*dv1+1.224744871391589*abar1[3]*f[17]+1.224744871391589*abar1[2]*f[7]+1.224744871391589*abar1[1]*f[1]+1.224744871391589*abar1[0]*f[0]; 
  incr1[5] = ((-0.3105295017040593*B2[3]*f[11])-0.3162277660168379*B2[1]*f[11]-0.3162277660168379*B2[2]*f[4]-0.3535533905932737*B2[0]*f[4]-0.3535533905932737*B2[1]*f[2])*dv1+abar1[2]*(1.075705748400954*f[17]+1.095445115010332*f[1])+abar1[1]*(1.095445115010332*f[7]+1.224744871391589*f[0])+1.075705748400954*abar1[3]*f[7]+1.224744871391589*abar1[0]*f[1]; 
  incr1[6] = ((-0.3535533905932737*B2[3]*f[17])-0.3162277660168379*B2[1]*f[12]-0.3162277660168379*B2[0]*f[8]-0.3535533905932737*B2[2]*f[7]-0.3535533905932737*f[1]*B2[1]-0.3535533905932737*f[0]*B2[0])*dv1+1.224744871391589*abar1[2]*f[11]+1.224744871391589*abar1[1]*f[4]+1.224744871391589*abar1[0]*f[2]; 
  incr1[9] = ((-0.7905694150420947*B2[1]*f[10])-0.7905694150420947*B2[0]*f[6])*dv1+2.738612787525831*abar1[2]*f[13]+2.738612787525831*abar1[1]*f[5]+2.738612787525831*abar1[0]*f[3]; 
  incr1[10] = ((-0.3105295017040592*B2[2]*f[17])-0.282842712474619*B2[2]*f[12]-0.3162277660168379*B2[0]*f[12]-0.3162277660168379*B2[1]*f[8]-0.3105295017040592*B2[3]*f[7]-0.3162277660168379*B2[1]*f[7]-0.3162277660168379*f[1]*B2[2]-0.3535533905932737*f[0]*B2[1]-0.3535533905932737*B2[0]*f[1])*dv1+abar1[1]*(1.095445115010332*f[11]+1.224744871391589*f[2])+1.075705748400954*abar1[3]*f[11]+1.095445115010332*abar1[2]*f[4]+1.224744871391589*abar1[0]*f[4]; 
  incr1[13] = ((-0.2258769757263128*B2[2]*f[11])-0.3535533905932737*B2[0]*f[11]-0.3105295017040593*B2[3]*f[4]-0.3162277660168379*B2[1]*f[4]-0.3535533905932737*f[2]*B2[2])*dv1+abar1[1]*(1.075705748400954*f[17]+1.095445115010332*f[1])+abar1[3]*(0.7302967433402215*f[17]+1.075705748400954*f[1])+1.224744871391589*abar1[0]*f[7]+abar1[2]*(0.7824607964359517*f[7]+1.224744871391589*f[0]); 
  incr1[14] = ((-0.3105295017040593*B2[0]*f[18])-0.3162277660168379*B2[2]*f[11]-0.3162277660168379*B2[1]*f[4]-0.3162277660168379*B2[0]*f[2])*dv1+1.224744871391589*abar1[1]*f[12]+1.224744871391589*abar1[0]*f[8]; 
  incr1[15] = ((-0.7071067811865475*B2[2]*f[10])-0.7905694150420948*B2[0]*f[10]-0.7905694150420948*B2[1]*f[6])*dv1+abar1[1]*(2.449489742783178*f[13]+2.738612787525831*f[3])+2.405351177211819*abar1[3]*f[13]+2.449489742783178*abar1[2]*f[5]+2.738612787525831*abar1[0]*f[5]; 
  incr1[16] = ((-0.7071067811865475*B2[0]*f[14])-0.7905694150420947*B2[2]*f[13]-0.7905694150420948*B2[1]*f[5]-0.7905694150420948*B2[0]*f[3])*dv1+2.738612787525831*abar1[1]*f[10]+2.738612787525831*abar1[0]*f[6]; 
  incr1[19] = ((-1.20761472884912*B2[0]*f[16])-0.5400617248673215*B2[2]*f[11]-0.5400617248673216*B2[1]*f[4]-0.5400617248673216*B2[0]*f[2])*dv1+1.870828693386971*abar1[3]*f[17]+abar1[1]*(4.183300132670379*f[15]+1.870828693386971*f[1])+abar1[0]*(4.183300132670378*f[9]+1.870828693386971*f[0])+1.870828693386971*abar1[2]*f[7]; 

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
  out[10] += incr1[10]*dv11+incr0[10]*dv10; 
  out[11] += incr1[11]*dv11+incr0[11]*dv10; 
  out[12] += incr1[12]*dv11+incr0[12]*dv10; 
  out[13] += incr1[13]*dv11+incr0[13]*dv10; 
  out[14] += incr1[14]*dv11+incr0[14]*dv10; 
  out[15] += incr1[15]*dv11+incr0[15]*dv10; 
  out[16] += incr1[16]*dv11+incr0[16]*dv10; 
  out[17] += incr1[17]*dv11+incr0[17]*dv10; 
  out[18] += incr1[18]*dv11+incr0[18]*dv10; 
  out[19] += incr1[19]*dv11+incr0[19]*dv10; 
return std::abs(amid1)/dv1+std::abs(amid2)/dv2; 
} 
double VlasovVolElcMag1x2vMaxP4(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[5]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[15]; 
  const double *B1 = &EM[20]; 
  const double *B2 = &EM[25]; 

  double abar0[5]; 

  double abar1[5]; 


  abar0[0] = E0[0]+wv2*B2[0]; 
  abar1[0] = E1[0]-wv1*B2[0]; 
  abar0[1] = E0[1]+wv2*B2[1]; 
  abar1[1] = E1[1]-wv1*B2[1]; 
  abar0[2] = E0[2]+wv2*B2[2]; 
  abar1[2] = E1[2]-wv1*B2[2]; 
  abar0[3] = E0[3]+wv2*B2[3]; 
  abar1[3] = E1[3]-wv1*B2[3]; 
  abar0[4] = E0[4]+wv2*B2[4]; 
  abar1[4] = E1[4]-wv1*B2[4]; 

  double incr0[35]; 

  double incr1[35]; 

  for(unsigned int i=0; i<35; ++i){ 

    incr0[i]=0.0; 

    incr1[i]=0.0; 

  }; 

  const double amid1 = 0.3977475644174328*B2[4]*dv2-0.3952847075210473*B2[2]*dv2+0.3535533905932737*B2[0]*dv2+0.7954951288348656*abar0[4]-0.7905694150420947*abar0[2]+0.7071067811865475*abar0[0]; 
  incr0[2] = (0.3535533905932737*B2[3]*f[28]+0.3535533905932737*B2[2]*f[13]+0.3535533905932737*B2[1]*f[5]+0.3535533905932737*B2[0]*f[3])*dv2+1.224744871391589*abar0[4]*f[32]+1.224744871391589*abar0[3]*f[17]+1.224744871391589*abar0[2]*f[7]+1.224744871391589*abar0[1]*f[1]+1.224744871391589*abar0[0]*f[0]; 
  incr0[4] = (0.3086066999241839*B2[4]*f[28]+0.3105295017040593*B2[2]*f[28]+0.3105295017040593*B2[3]*f[13]+0.3162277660168379*B2[1]*f[13]+0.3162277660168379*B2[2]*f[5]+0.3535533905932737*B2[0]*f[5]+0.3535533905932737*B2[1]*f[3])*dv2+abar0[3]*(1.069044967649698*f[32]+1.075705748400954*f[7])+abar0[2]*(1.075705748400954*f[17]+1.095445115010332*f[1])+1.069044967649698*abar0[4]*f[17]+abar0[1]*(1.095445115010332*f[7]+1.224744871391589*f[0])+1.224744871391589*abar0[0]*f[1]; 
  incr0[6] = (0.3535533905932737*B2[4]*f[32]+0.3162277660168379*B2[2]*f[24]+0.3535533905932737*B2[3]*f[17]+0.3162277660168379*B2[1]*f[15]+0.3162277660168379*B2[0]*f[9]+0.3535533905932737*B2[2]*f[7]+0.3535533905932737*f[1]*B2[1]+0.3535533905932737*f[0]*B2[0])*dv2+1.224744871391589*abar0[3]*f[28]+1.224744871391589*abar0[2]*f[13]+1.224744871391589*abar0[1]*f[5]+1.224744871391589*abar0[0]*f[3]; 
  incr0[8] = (0.7905694150420947*B2[2]*f[20]+0.7905694150420947*B2[1]*f[10]+0.7905694150420947*B2[0]*f[6])*dv2+2.73861278752583*abar0[3]*f[26]+2.738612787525831*abar0[2]*f[11]+2.738612787525831*abar0[1]*f[4]+2.738612787525831*abar0[0]*f[2]; 
  incr0[10] = (0.3086066999241838*B2[3]*f[32]+0.2777460299317654*B2[3]*f[24]+0.2828427124746191*B2[1]*f[24]+0.3086066999241838*B2[4]*f[17]+0.3105295017040592*B2[2]*f[17]+0.282842712474619*B2[2]*f[15]+0.3162277660168379*B2[0]*f[15]+0.3162277660168379*B2[1]*f[9]+0.3105295017040592*B2[3]*f[7]+0.3162277660168379*B2[1]*f[7]+0.3162277660168379*f[1]*B2[2]+0.3535533905932737*f[0]*B2[1]+0.3535533905932737*B2[0]*f[1])*dv2+abar0[2]*(1.075705748400954*f[28]+1.095445115010332*f[5])+1.069044967649698*abar0[4]*f[28]+abar0[1]*(1.095445115010332*f[13]+1.224744871391589*f[3])+1.075705748400954*abar0[3]*f[13]+1.224744871391589*abar0[0]*f[5]; 
  incr0[11] = (0.2108185106778919*B2[3]*f[28]+0.3105295017040593*B2[1]*f[28]+0.3030457633656632*B2[4]*f[13]+0.2258769757263128*B2[2]*f[13]+0.3535533905932737*B2[0]*f[13]+0.3105295017040593*B2[3]*f[5]+0.3162277660168379*B2[1]*f[5]+0.3535533905932737*B2[2]*f[3])*dv2+abar0[2]*(1.049781318335648*f[32]+0.7824607964359517*f[7]+1.224744871391589*f[0])+abar0[4]*(0.7113279967599562*f[32]+1.049781318335648*f[7])+abar0[1]*(1.075705748400954*f[17]+1.095445115010332*f[1])+abar0[3]*(0.7302967433402215*f[17]+1.075705748400954*f[1])+1.224744871391589*abar0[0]*f[7]; 
  incr0[12] = (0.6943650748294133*B2[3]*f[20]+0.7071067811865475*B2[1]*f[20]+0.7071067811865475*B2[2]*f[10]+0.7905694150420948*B2[0]*f[10]+0.7905694150420948*B2[1]*f[6])*dv2+abar0[2]*(2.405351177211819*f[26]+2.449489742783178*f[4])+2.390457218668788*abar0[4]*f[26]+abar0[1]*(2.449489742783178*f[11]+2.738612787525831*f[2])+2.405351177211819*abar0[3]*f[11]+2.738612787525831*abar0[0]*f[4]; 
  incr0[14] = (0.7905694150420945*B2[3]*f[26]+0.7071067811865475*B2[1]*f[22]+0.7071067811865475*B2[0]*f[16]+0.7905694150420947*B2[2]*f[11]+0.7905694150420948*B2[1]*f[4]+0.7905694150420948*B2[0]*f[2])*dv2+2.738612787525831*abar0[2]*f[20]+2.738612787525831*abar0[1]*f[10]+2.738612787525831*abar0[0]*f[6]; 
  incr0[16] = (0.3105295017040593*B2[1]*f[30]+0.3162277660168379*B2[3]*f[28]+0.3105295017040593*B2[0]*f[19]+0.3162277660168379*B2[2]*f[13]+0.3162277660168379*B2[1]*f[5]+0.3162277660168379*B2[0]*f[3])*dv2+1.224744871391589*abar0[2]*f[24]+1.224744871391589*abar0[1]*f[15]+1.224744871391589*abar0[0]*f[9]; 
  incr0[18] = (0.5400617248673215*B2[3]*f[28]+1.20761472884912*B2[1]*f[21]+1.20761472884912*B2[0]*f[14]+0.5400617248673215*B2[2]*f[13]+0.5400617248673216*B2[1]*f[5]+0.5400617248673216*B2[0]*f[3])*dv2+1.870828693386971*abar0[4]*f[32]+abar0[2]*(4.183300132670378*f[23]+1.870828693386971*f[7])+1.870828693386971*abar0[3]*f[17]+abar0[1]*(4.183300132670379*f[12]+1.870828693386971*f[1])+abar0[0]*(4.183300132670378*f[8]+1.870828693386971*f[0]); 
  incr0[20] = (0.205342705205739*B2[4]*f[32]+0.3030457633656632*B2[2]*f[32]+0.2710523708715754*B2[4]*f[24]+0.2020305089104422*B2[2]*f[24]+0.3162277660168379*B2[0]*f[24]+0.210818510677892*B2[3]*f[17]+0.3105295017040592*B2[1]*f[17]+0.2777460299317654*B2[3]*f[15]+0.282842712474619*B2[1]*f[15]+0.3162277660168379*B2[2]*f[9]+0.3030457633656632*B2[4]*f[7]+0.2258769757263128*B2[2]*f[7]+0.3535533905932737*B2[0]*f[7]+0.3105295017040592*f[1]*B2[3]+0.3535533905932737*f[0]*B2[2]+0.3162277660168379*f[1]*B2[1])*dv2+abar0[1]*(1.075705748400954*f[28]+1.095445115010332*f[5])+abar0[3]*(0.7302967433402215*f[28]+1.075705748400954*f[5])+1.049781318335648*abar0[4]*f[13]+1.224744871391589*abar0[0]*f[13]+abar0[2]*(0.7824607964359517*f[13]+1.224744871391589*f[3]); 
  incr0[21] = (0.6900655593423543*B2[4]*f[26]+0.6943650748294133*B2[2]*f[26]+0.6324555320336759*B2[2]*f[22]+0.7071067811865475*B2[0]*f[22]+0.7071067811865475*B2[1]*f[16]+0.6943650748294133*B2[3]*f[11]+0.7071067811865475*B2[1]*f[11]+0.7071067811865475*B2[2]*f[4]+0.7905694150420947*B2[0]*f[4]+0.7905694150420947*B2[1]*f[2])*dv2+abar0[1]*(2.449489742783178*f[20]+2.738612787525831*f[6])+2.405351177211819*abar0[3]*f[20]+2.449489742783178*abar0[2]*f[10]+2.738612787525831*abar0[0]*f[10]; 
  incr0[22] = (0.2777460299317654*B2[2]*f[30]+0.3105295017040593*B2[0]*f[30]+0.2760262237369417*B2[4]*f[28]+0.2777460299317654*B2[2]*f[28]+0.3105295017040592*B2[1]*f[19]+0.2777460299317654*B2[3]*f[13]+0.282842712474619*B2[1]*f[13]+0.2828427124746191*B2[2]*f[5]+0.3162277660168379*B2[0]*f[5]+0.3162277660168379*B2[1]*f[3])*dv2+abar0[1]*(1.095445115010332*f[24]+1.224744871391589*f[9])+1.075705748400954*abar0[3]*f[24]+1.095445115010332*abar0[2]*f[15]+1.224744871391589*abar0[0]*f[15]; 
  incr0[23] = (0.6776309271789384*B2[4]*f[20]+0.5050762722761053*B2[2]*f[20]+0.7905694150420947*B2[0]*f[20]+0.6943650748294133*B2[3]*f[10]+0.7071067811865475*B2[1]*f[10]+0.7905694150420947*B2[2]*f[6])*dv2+abar0[1]*(2.405351177211819*f[26]+2.449489742783178*f[4])+abar0[3]*(1.632993161855452*f[26]+2.405351177211819*f[4])+abar0[2]*(1.749635530559413*f[11]+2.738612787525831*f[2])+2.347382389307855*abar0[4]*f[11]+2.738612787525831*abar0[0]*f[11]; 
  incr0[25] = (0.6943650748294133*B2[0]*f[31]+0.7071067811865475*B2[2]*f[20]+0.7071067811865475*B2[1]*f[10]+0.7071067811865475*B2[0]*f[6])*dv2+2.738612787525831*abar0[1]*f[22]+2.738612787525831*abar0[0]*f[16]; 
  incr0[26] = (0.1928473039599675*B2[4]*f[28]+0.210818510677892*B2[2]*f[28]+0.3535533905932737*B2[0]*f[28]+0.2108185106778919*B2[3]*f[13]+0.3105295017040593*B2[1]*f[13]+0.3086066999241839*B2[4]*f[5]+0.3105295017040593*B2[2]*f[5]+0.3535533905932737*f[3]*B2[3])*dv2+abar0[1]*(1.069044967649698*f[32]+1.075705748400954*f[7])+abar0[3]*(0.6680426571226847*f[32]+0.7302967433402215*f[7]+1.224744871391589*f[0])+1.224744871391589*abar0[0]*f[17]+abar0[2]*(0.7302967433402215*f[17]+1.075705748400954*f[1])+abar0[4]*(0.6680426571226847*f[17]+1.069044967649698*f[1]); 
  incr0[27] = (0.4714045207910317*B2[4]*f[28]+0.4743416490252568*B2[2]*f[28]+1.080123449734643*B2[2]*f[21]+1.20761472884912*B2[0]*f[21]+1.20761472884912*B2[1]*f[14]+0.4743416490252568*B2[3]*f[13]+0.4830458915396481*B2[1]*f[13]+0.4830458915396479*B2[2]*f[5]+0.5400617248673215*B2[0]*f[5]+0.5400617248673215*B2[1]*f[3])*dv2+abar0[3]*(1.632993161855452*f[32]+3.674234614174766*f[23]+1.643167672515498*f[7])+abar0[1]*(3.741657386773942*f[23]+4.183300132670378*f[8]+1.673320053068151*f[7]+1.870828693386971*f[0])+abar0[2]*(1.643167672515498*f[17]+3.741657386773941*f[12]+1.673320053068151*f[1])+1.632993161855452*abar0[4]*f[17]+abar0[0]*(4.183300132670377*f[12]+1.870828693386971*f[1]); 
  incr0[29] = (0.5400617248673215*B2[4]*f[32]+1.080123449734643*B2[0]*f[25]+0.4830458915396479*B2[2]*f[24]+1.20761472884912*B2[2]*f[23]+0.5400617248673215*B2[3]*f[17]+0.4830458915396481*B2[1]*f[15]+1.20761472884912*B2[1]*f[12]+0.4830458915396479*B2[0]*f[9]+1.20761472884912*B2[0]*f[8]+0.5400617248673215*B2[2]*f[7]+0.5400617248673215*f[1]*B2[1]+0.5400617248673215*f[0]*B2[0])*dv2+1.870828693386971*abar0[3]*f[28]+abar0[1]*(4.183300132670378*f[21]+1.870828693386971*f[5])+abar0[0]*(4.183300132670377*f[14]+1.870828693386971*f[3])+1.87082869338697*abar0[2]*f[13]; 
  incr0[31] = (0.3086066999241839*B2[0]*f[34]+0.3105295017040593*B2[2]*f[24]+0.3105295017040593*B2[1]*f[15]+0.3105295017040593*B2[0]*f[9])*dv2+1.224744871391589*abar0[1]*f[30]+1.224744871391589*abar0[0]*f[19]; 
  incr0[33] = (1.620185174601965*B2[0]*f[29]+1.060660171779821*B2[2]*f[20]+1.060660171779821*B2[1]*f[10]+1.060660171779821*B2[0]*f[6])*dv2+abar0[1]*(5.612486080160911*f[27]+3.674234614174766*f[4])+3.674234614174766*abar0[3]*f[26]+abar0[0]*(5.612486080160912*f[18]+3.674234614174766*f[2])+3.674234614174767*abar0[2]*f[11]; 
  const double amid2 = (-0.3977475644174328*B2[4]*dv1)+0.3952847075210473*B2[2]*dv1-0.3535533905932737*B2[0]*dv1+0.7954951288348656*abar1[4]-0.7905694150420947*abar1[2]+0.7071067811865475*abar1[0]; 
  incr1[3] = ((-0.3535533905932737*B2[3]*f[26])-0.3535533905932737*B2[2]*f[11]-0.3535533905932737*B2[1]*f[4]-0.3535533905932737*B2[0]*f[2])*dv1+1.224744871391589*abar1[4]*f[32]+1.224744871391589*abar1[3]*f[17]+1.224744871391589*abar1[2]*f[7]+1.224744871391589*abar1[1]*f[1]+1.224744871391589*abar1[0]*f[0]; 
  incr1[5] = ((-0.3086066999241839*B2[4]*f[26])-0.3105295017040593*B2[2]*f[26]-0.3105295017040593*B2[3]*f[11]-0.3162277660168379*B2[1]*f[11]-0.3162277660168379*B2[2]*f[4]-0.3535533905932737*B2[0]*f[4]-0.3535533905932737*B2[1]*f[2])*dv1+abar1[3]*(1.069044967649698*f[32]+1.075705748400954*f[7])+abar1[2]*(1.075705748400954*f[17]+1.095445115010332*f[1])+1.069044967649698*abar1[4]*f[17]+abar1[1]*(1.095445115010332*f[7]+1.224744871391589*f[0])+1.224744871391589*abar1[0]*f[1]; 
  incr1[6] = ((-0.3535533905932737*B2[4]*f[32])-0.3162277660168379*B2[2]*f[23]-0.3535533905932737*B2[3]*f[17]-0.3162277660168379*B2[1]*f[12]-0.3162277660168379*B2[0]*f[8]-0.3535533905932737*B2[2]*f[7]-0.3535533905932737*f[1]*B2[1]-0.3535533905932737*f[0]*B2[0])*dv1+1.224744871391589*abar1[3]*f[26]+1.224744871391589*abar1[2]*f[11]+1.224744871391589*abar1[1]*f[4]+1.224744871391589*abar1[0]*f[2]; 
  incr1[9] = ((-0.7905694150420947*B2[2]*f[20])-0.7905694150420947*B2[1]*f[10]-0.7905694150420947*B2[0]*f[6])*dv1+2.73861278752583*abar1[3]*f[28]+2.738612787525831*abar1[2]*f[13]+2.738612787525831*abar1[1]*f[5]+2.738612787525831*abar1[0]*f[3]; 
  incr1[10] = ((-0.3086066999241838*B2[3]*f[32])-0.2777460299317654*B2[3]*f[23]-0.2828427124746191*B2[1]*f[23]-0.3086066999241838*B2[4]*f[17]-0.3105295017040592*B2[2]*f[17]-0.282842712474619*B2[2]*f[12]-0.3162277660168379*B2[0]*f[12]-0.3162277660168379*B2[1]*f[8]-0.3105295017040592*B2[3]*f[7]-0.3162277660168379*B2[1]*f[7]-0.3162277660168379*f[1]*B2[2]-0.3535533905932737*f[0]*B2[1]-0.3535533905932737*B2[0]*f[1])*dv1+abar1[2]*(1.075705748400954*f[26]+1.095445115010332*f[4])+1.069044967649698*abar1[4]*f[26]+abar1[1]*(1.095445115010332*f[11]+1.224744871391589*f[2])+1.075705748400954*abar1[3]*f[11]+1.224744871391589*abar1[0]*f[4]; 
  incr1[13] = ((-0.2108185106778919*B2[3]*f[26])-0.3105295017040593*B2[1]*f[26]-0.3030457633656632*B2[4]*f[11]-0.2258769757263128*B2[2]*f[11]-0.3535533905932737*B2[0]*f[11]-0.3105295017040593*B2[3]*f[4]-0.3162277660168379*B2[1]*f[4]-0.3535533905932737*f[2]*B2[2])*dv1+abar1[2]*(1.049781318335648*f[32]+0.7824607964359517*f[7]+1.224744871391589*f[0])+abar1[4]*(0.7113279967599562*f[32]+1.049781318335648*f[7])+abar1[1]*(1.075705748400954*f[17]+1.095445115010332*f[1])+abar1[3]*(0.7302967433402215*f[17]+1.075705748400954*f[1])+1.224744871391589*abar1[0]*f[7]; 
  incr1[14] = ((-0.3105295017040593*B2[1]*f[27])-0.3162277660168379*B2[3]*f[26]-0.3105295017040593*B2[0]*f[18]-0.3162277660168379*B2[2]*f[11]-0.3162277660168379*B2[1]*f[4]-0.3162277660168379*B2[0]*f[2])*dv1+1.224744871391589*abar1[2]*f[23]+1.224744871391589*abar1[1]*f[12]+1.224744871391589*abar1[0]*f[8]; 
  incr1[15] = ((-0.6943650748294133*B2[3]*f[20])-0.7071067811865475*B2[1]*f[20]-0.7071067811865475*B2[2]*f[10]-0.7905694150420948*B2[0]*f[10]-0.7905694150420948*B2[1]*f[6])*dv1+abar1[2]*(2.405351177211819*f[28]+2.449489742783178*f[5])+2.390457218668788*abar1[4]*f[28]+abar1[1]*(2.449489742783178*f[13]+2.738612787525831*f[3])+2.405351177211819*abar1[3]*f[13]+2.738612787525831*abar1[0]*f[5]; 
  incr1[16] = ((-0.7905694150420945*B2[3]*f[28])-0.7071067811865475*B2[1]*f[21]-0.7071067811865475*B2[0]*f[14]-0.7905694150420947*B2[2]*f[13]-0.7905694150420948*B2[1]*f[5]-0.7905694150420948*B2[0]*f[3])*dv1+2.738612787525831*abar1[2]*f[20]+2.738612787525831*abar1[1]*f[10]+2.738612787525831*abar1[0]*f[6]; 
  incr1[19] = ((-0.5400617248673215*B2[3]*f[26])-1.20761472884912*B2[1]*f[22]-1.20761472884912*B2[0]*f[16]-0.5400617248673215*B2[2]*f[11]-0.5400617248673216*B2[1]*f[4]-0.5400617248673216*B2[0]*f[2])*dv1+1.870828693386971*abar1[4]*f[32]+abar1[2]*(4.183300132670378*f[24]+1.870828693386971*f[7])+1.870828693386971*abar1[3]*f[17]+abar1[1]*(4.183300132670379*f[15]+1.870828693386971*f[1])+abar1[0]*(4.183300132670378*f[9]+1.870828693386971*f[0]); 
  incr1[20] = ((-0.205342705205739*B2[4]*f[32])-0.3030457633656632*B2[2]*f[32]-0.2710523708715754*B2[4]*f[23]-0.2020305089104422*B2[2]*f[23]-0.3162277660168379*B2[0]*f[23]-0.210818510677892*B2[3]*f[17]-0.3105295017040592*B2[1]*f[17]-0.2777460299317654*B2[3]*f[12]-0.282842712474619*B2[1]*f[12]-0.3162277660168379*B2[2]*f[8]-0.3030457633656632*B2[4]*f[7]-0.2258769757263128*B2[2]*f[7]-0.3535533905932737*B2[0]*f[7]-0.3105295017040592*f[1]*B2[3]-0.3535533905932737*f[0]*B2[2]-0.3162277660168379*f[1]*B2[1])*dv1+abar1[1]*(1.075705748400954*f[26]+1.095445115010332*f[4])+abar1[3]*(0.7302967433402215*f[26]+1.075705748400954*f[4])+1.049781318335648*abar1[4]*f[11]+1.224744871391589*abar1[0]*f[11]+abar1[2]*(0.7824607964359517*f[11]+1.224744871391589*f[2]); 
  incr1[21] = ((-0.2777460299317654*B2[2]*f[27])-0.3105295017040593*B2[0]*f[27]-0.2760262237369417*B2[4]*f[26]-0.2777460299317654*B2[2]*f[26]-0.3105295017040592*B2[1]*f[18]-0.2777460299317654*B2[3]*f[11]-0.282842712474619*B2[1]*f[11]-0.2828427124746191*B2[2]*f[4]-0.3162277660168379*B2[0]*f[4]-0.3162277660168379*B2[1]*f[2])*dv1+abar1[1]*(1.095445115010332*f[23]+1.224744871391589*f[8])+1.075705748400954*abar1[3]*f[23]+1.095445115010332*abar1[2]*f[12]+1.224744871391589*abar1[0]*f[12]; 
  incr1[22] = ((-0.6900655593423543*B2[4]*f[28])-0.6943650748294133*B2[2]*f[28]-0.6324555320336759*B2[2]*f[21]-0.7071067811865475*B2[0]*f[21]-0.7071067811865475*B2[1]*f[14]-0.6943650748294133*B2[3]*f[13]-0.7071067811865475*B2[1]*f[13]-0.7071067811865475*B2[2]*f[5]-0.7905694150420947*B2[0]*f[5]-0.7905694150420947*B2[1]*f[3])*dv1+abar1[1]*(2.449489742783178*f[20]+2.738612787525831*f[6])+2.405351177211819*abar1[3]*f[20]+2.449489742783178*abar1[2]*f[10]+2.738612787525831*abar1[0]*f[10]; 
  incr1[24] = ((-0.6776309271789384*B2[4]*f[20])-0.5050762722761053*B2[2]*f[20]-0.7905694150420947*B2[0]*f[20]-0.6943650748294133*B2[3]*f[10]-0.7071067811865475*B2[1]*f[10]-0.7905694150420947*B2[2]*f[6])*dv1+abar1[1]*(2.405351177211819*f[28]+2.449489742783178*f[5])+abar1[3]*(1.632993161855452*f[28]+2.405351177211819*f[5])+abar1[2]*(1.749635530559413*f[13]+2.738612787525831*f[3])+2.347382389307855*abar1[4]*f[13]+2.738612787525831*abar1[0]*f[13]; 
  incr1[25] = ((-0.6943650748294133*B2[0]*f[29])-0.7071067811865475*B2[2]*f[20]-0.7071067811865475*B2[1]*f[10]-0.7071067811865475*B2[0]*f[6])*dv1+2.738612787525831*abar1[1]*f[21]+2.738612787525831*abar1[0]*f[14]; 
  incr1[28] = ((-0.1928473039599675*B2[4]*f[26])-0.210818510677892*B2[2]*f[26]-0.3535533905932737*B2[0]*f[26]-0.2108185106778919*B2[3]*f[11]-0.3105295017040593*B2[1]*f[11]-0.3086066999241839*f[4]*B2[4]-0.3105295017040593*B2[2]*f[4]-0.3535533905932737*f[2]*B2[3])*dv1+abar1[1]*(1.069044967649698*f[32]+1.075705748400954*f[7])+abar1[3]*(0.6680426571226847*f[32]+0.7302967433402215*f[7]+1.224744871391589*f[0])+1.224744871391589*abar1[0]*f[17]+abar1[2]*(0.7302967433402215*f[17]+1.075705748400954*f[1])+abar1[4]*(0.6680426571226847*f[17]+1.069044967649698*f[1]); 
  incr1[29] = ((-0.3086066999241839*B2[0]*f[33])-0.3105295017040593*B2[2]*f[23]-0.3105295017040593*B2[1]*f[12]-0.3105295017040593*B2[0]*f[8])*dv1+1.224744871391589*abar1[1]*f[27]+1.224744871391589*abar1[0]*f[18]; 
  incr1[30] = ((-0.4714045207910317*B2[4]*f[26])-0.4743416490252568*B2[2]*f[26]-1.080123449734643*B2[2]*f[22]-1.20761472884912*B2[0]*f[22]-1.20761472884912*B2[1]*f[16]-0.4743416490252568*B2[3]*f[11]-0.4830458915396481*B2[1]*f[11]-0.4830458915396479*B2[2]*f[4]-0.5400617248673215*B2[0]*f[4]-0.5400617248673215*B2[1]*f[2])*dv1+abar1[3]*(1.632993161855452*f[32]+3.674234614174766*f[24]+1.643167672515498*f[7])+abar1[1]*(3.741657386773942*f[24]+4.183300132670378*f[9]+1.673320053068151*f[7]+1.870828693386971*f[0])+abar1[2]*(1.643167672515498*f[17]+3.741657386773941*f[15]+1.673320053068151*f[1])+1.632993161855452*abar1[4]*f[17]+abar1[0]*(4.183300132670377*f[15]+1.870828693386971*f[1]); 
  incr1[31] = ((-0.5400617248673215*B2[4]*f[32])-1.080123449734643*B2[0]*f[25]-1.20761472884912*B2[2]*f[24]-0.4830458915396479*B2[2]*f[23]-0.5400617248673215*B2[3]*f[17]-1.20761472884912*B2[1]*f[15]-0.4830458915396481*B2[1]*f[12]-1.20761472884912*B2[0]*f[9]-0.4830458915396479*B2[0]*f[8]-0.5400617248673215*B2[2]*f[7]-0.5400617248673215*f[1]*B2[1]-0.5400617248673215*f[0]*B2[0])*dv1+1.870828693386971*abar1[3]*f[26]+abar1[1]*(4.183300132670378*f[22]+1.870828693386971*f[4])+abar1[0]*(4.183300132670377*f[16]+1.870828693386971*f[2])+1.87082869338697*abar1[2]*f[11]; 
  incr1[34] = ((-1.620185174601965*B2[0]*f[31])-1.060660171779821*B2[2]*f[20]-1.060660171779821*B2[1]*f[10]-1.060660171779821*B2[0]*f[6])*dv1+abar1[1]*(5.612486080160911*f[30]+3.674234614174766*f[5])+3.674234614174766*abar1[3]*f[28]+abar1[0]*(5.612486080160912*f[19]+3.674234614174766*f[3])+3.674234614174767*abar1[2]*f[13]; 

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
  out[10] += incr1[10]*dv11+incr0[10]*dv10; 
  out[11] += incr1[11]*dv11+incr0[11]*dv10; 
  out[12] += incr1[12]*dv11+incr0[12]*dv10; 
  out[13] += incr1[13]*dv11+incr0[13]*dv10; 
  out[14] += incr1[14]*dv11+incr0[14]*dv10; 
  out[15] += incr1[15]*dv11+incr0[15]*dv10; 
  out[16] += incr1[16]*dv11+incr0[16]*dv10; 
  out[17] += incr1[17]*dv11+incr0[17]*dv10; 
  out[18] += incr1[18]*dv11+incr0[18]*dv10; 
  out[19] += incr1[19]*dv11+incr0[19]*dv10; 
  out[20] += incr1[20]*dv11+incr0[20]*dv10; 
  out[21] += incr1[21]*dv11+incr0[21]*dv10; 
  out[22] += incr1[22]*dv11+incr0[22]*dv10; 
  out[23] += incr1[23]*dv11+incr0[23]*dv10; 
  out[24] += incr1[24]*dv11+incr0[24]*dv10; 
  out[25] += incr1[25]*dv11+incr0[25]*dv10; 
  out[26] += incr1[26]*dv11+incr0[26]*dv10; 
  out[27] += incr1[27]*dv11+incr0[27]*dv10; 
  out[28] += incr1[28]*dv11+incr0[28]*dv10; 
  out[29] += incr1[29]*dv11+incr0[29]*dv10; 
  out[30] += incr1[30]*dv11+incr0[30]*dv10; 
  out[31] += incr1[31]*dv11+incr0[31]*dv10; 
  out[32] += incr1[32]*dv11+incr0[32]*dv10; 
  out[33] += incr1[33]*dv11+incr0[33]*dv10; 
  out[34] += incr1[34]*dv11+incr0[34]*dv10; 
return std::abs(amid1)/dv1+std::abs(amid2)/dv2; 
} 
