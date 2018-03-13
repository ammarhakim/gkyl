#include <VlasovModDecl.h> 
void VlasovVolElcMag2x2vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[3]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

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

  double incr0[5]; 

  double incr1[5]; 

  for(unsigned int i=0; i<5; ++i){ 

    incr0[i]=0.0; 

    incr1[i]=0.0; 

  }; 

  incr0[3] = 0.25*B2[0]*f[4]*dv2+0.8660254037844386*abar0[2]*f[2]+0.8660254037844386*abar0[1]*f[1]+0.8660254037844386*abar0[0]*f[0]; 
  incr1[4] = (-0.25*B2[0]*f[3]*dv1)+0.8660254037844386*abar1[2]*f[2]+0.8660254037844386*abar1[1]*f[1]+0.8660254037844386*abar1[0]*f[0]; 

  out[0] += incr1[0]*dv11+incr0[0]*dv10; 
  out[1] += incr1[1]*dv11+incr0[1]*dv10; 
  out[2] += incr1[2]*dv11+incr0[2]*dv10; 
  out[3] += incr1[3]*dv11+incr0[3]*dv10; 
  out[4] += incr1[4]*dv11+incr0[4]*dv10; 
} 
void VlasovVolElcMag2x2vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[6]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double abar0[6]; 

  double abar1[6]; 


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
  abar0[5] = E0[5]+wv2*B2[5]; 
  abar1[5] = E1[5]-wv1*B2[5]; 

  double incr0[15]; 

  double incr1[15]; 

  for(unsigned int i=0; i<15; ++i){ 

    incr0[i]=0.0; 

    incr1[i]=0.0; 

  }; 

  incr0[3] = 0.25*B2[2]*f[9]*dv2+0.25*B2[1]*f[8]*dv2+0.25*B2[0]*f[4]*dv2+0.8660254037844386*abar0[5]*f[12]+0.8660254037844386*abar0[4]*f[11]+0.8660254037844386*abar0[3]*f[5]+0.8660254037844386*abar0[2]*f[2]+0.8660254037844386*abar0[1]*f[1]+0.8660254037844386*abar0[0]*f[0]; 
  incr0[6] = 0.25*B2[3]*f[9]*dv2+0.223606797749979*B2[4]*f[8]*dv2+0.25*B2[0]*f[8]*dv2+0.25*B2[1]*f[4]*dv2+0.7745966692414833*abar0[1]*f[11]+0.8660254037844386*abar0[2]*f[5]+0.7745966692414833*f[1]*abar0[4]+0.8660254037844386*f[2]*abar0[3]+0.8660254037844386*abar0[0]*f[1]+0.8660254037844386*f[0]*abar0[1]; 
  incr0[7] = 0.223606797749979*B2[5]*f[9]*dv2+0.25*B2[0]*f[9]*dv2+0.25*B2[3]*f[8]*dv2+0.25*B2[2]*f[4]*dv2+0.7745966692414833*abar0[2]*f[12]+0.8660254037844386*abar0[1]*f[5]+0.7745966692414833*f[2]*abar0[5]+0.8660254037844386*f[1]*abar0[3]+0.8660254037844386*abar0[0]*f[2]+0.8660254037844386*f[0]*abar0[2]; 
  incr0[10] = 0.223606797749979*B2[0]*f[14]*dv2+0.25*B2[5]*f[12]*dv2+0.25*B2[4]*f[11]*dv2+0.25*B2[3]*f[5]*dv2+0.25*f[2]*B2[2]*dv2+0.25*f[1]*B2[1]*dv2+0.25*f[0]*B2[0]*dv2+0.8660254037844386*abar0[2]*f[9]+0.8660254037844386*abar0[1]*f[8]+0.8660254037844386*abar0[0]*f[4]; 
  incr0[13] = 0.5590169943749475*B2[0]*f[10]*dv2+1.936491673103709*abar0[2]*f[7]+1.936491673103709*abar0[1]*f[6]+1.936491673103709*abar0[0]*f[3]; 
  incr1[4] = (-0.25*B2[2]*f[7]*dv1)-0.25*B2[1]*f[6]*dv1-0.25*B2[0]*f[3]*dv1+0.8660254037844386*abar1[5]*f[12]+0.8660254037844386*abar1[4]*f[11]+0.8660254037844386*abar1[3]*f[5]+0.8660254037844386*abar1[2]*f[2]+0.8660254037844386*abar1[1]*f[1]+0.8660254037844386*abar1[0]*f[0]; 
  incr1[8] = (-0.25*B2[3]*f[7]*dv1)-0.223606797749979*B2[4]*f[6]*dv1-0.25*B2[0]*f[6]*dv1-0.25*B2[1]*f[3]*dv1+0.7745966692414833*abar1[1]*f[11]+0.8660254037844386*abar1[2]*f[5]+0.7745966692414833*f[1]*abar1[4]+0.8660254037844386*f[2]*abar1[3]+0.8660254037844386*abar1[0]*f[1]+0.8660254037844386*f[0]*abar1[1]; 
  incr1[9] = (-0.223606797749979*B2[5]*f[7]*dv1)-0.25*B2[0]*f[7]*dv1-0.25*B2[3]*f[6]*dv1-0.25*B2[2]*f[3]*dv1+0.7745966692414833*abar1[2]*f[12]+0.8660254037844386*abar1[1]*f[5]+0.7745966692414833*f[2]*abar1[5]+0.8660254037844386*f[1]*abar1[3]+0.8660254037844386*abar1[0]*f[2]+0.8660254037844386*f[0]*abar1[2]; 
  incr1[10] = (-0.223606797749979*B2[0]*f[13]*dv1)-0.25*B2[5]*f[12]*dv1-0.25*B2[4]*f[11]*dv1-0.25*B2[3]*f[5]*dv1-0.25*f[2]*B2[2]*dv1-0.25*f[1]*B2[1]*dv1-0.25*f[0]*B2[0]*dv1+0.8660254037844386*abar1[2]*f[7]+0.8660254037844386*abar1[1]*f[6]+0.8660254037844386*abar1[0]*f[3]; 
  incr1[14] = (-0.5590169943749475*B2[0]*f[10]*dv1)+1.936491673103709*abar1[2]*f[9]+1.936491673103709*abar1[1]*f[8]+1.936491673103709*abar1[0]*f[4]; 

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
} 
