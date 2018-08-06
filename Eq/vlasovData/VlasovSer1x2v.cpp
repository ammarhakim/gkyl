#include <VlasovModDecl.h> 
double VlasovVol1x2vSerP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[2]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[6]; 
  const double *B1 = &EM[8]; 
  const double *B2 = &EM[10]; 

  double alpha0[8]; 

  double alpha1[8]; 

  double alpha2[8]; 

  alpha0[0] = 5.656854249492382*w0dx0; 
  alpha0[2] = 1.632993161855453*dv0dx0; 

  alpha1[0] = dv10*(2.0*B2[0]*wv2+2.0*E0[0]); 
  alpha1[1] = dv10*(2.0*B2[1]*wv2+2.0*E0[1]); 
  alpha1[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha1[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha2[0] = 2.0*E1[0]*dv11-2.0*B2[0]*dv11*wv1; 
  alpha2[1] = 2.0*E1[1]*dv11-2.0*B2[1]*dv11*wv1; 
  alpha2[2] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha2[4] = -0.5773502691896258*B2[1]*dv1*dv11; 
  const double amid1 = 0.3535533905932737*alpha1[0]; 
  const double amid2 = 0.3535533905932737*alpha2[0]; 
  out[1] += 0.6123724356957944*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha1[5]*f[5]+alpha1[3]*f[3]+alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[3] += 0.6123724356957944*(alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[4] += 0.6123724356957944*(alpha1[3]*f[5]+f[3]*alpha1[5]+alpha0[0]*f[2]+f[0]*alpha0[2]+alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[5] += 0.6123724356957944*(alpha0[2]*f[6]+alpha2[2]*f[4]+f[2]*alpha2[4]+alpha0[0]*f[3]+alpha2[0]*f[1]+f[0]*alpha2[1]); 
  out[6] += 0.6123724356957944*(alpha1[1]*f[5]+f[1]*alpha1[5]+alpha2[1]*f[4]+f[1]*alpha2[4]+alpha1[0]*f[3]+f[0]*alpha1[3]+alpha2[0]*f[2]+f[0]*alpha2[2]); 
  out[7] += 0.6123724356957944*(alpha0[0]*f[6]+alpha1[0]*f[5]+f[0]*alpha1[5]+alpha2[0]*f[4]+f[0]*alpha2[4]+(alpha0[2]+alpha1[1])*f[3]+f[1]*alpha1[3]+alpha2[1]*f[2]+f[1]*alpha2[2]); 
return std::abs(w0dx0)+0.5*(dv0dx0+std::abs(amid1)+std::abs(amid2)); 
} 
double VlasovVol1x2vSerP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[3]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double alpha0[20]; 

  double alpha1[20]; 

  double alpha2[20]; 

  alpha0[0] = 5.656854249492382*w0dx0; 
  alpha0[2] = 1.632993161855453*dv0dx0; 

  alpha1[0] = dv10*(2.0*B2[0]*wv2+2.0*E0[0]); 
  alpha1[1] = dv10*(2.0*B2[1]*wv2+2.0*E0[1]); 
  alpha1[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha1[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha1[7] = dv10*(2.0*B2[2]*wv2+2.0*E0[2]); 
  alpha1[13] = 0.5773502691896257*B2[2]*dv10*dv2; 
  alpha2[0] = 2.0*E1[0]*dv11-2.0*B2[0]*dv11*wv1; 
  alpha2[1] = 2.0*E1[1]*dv11-2.0*B2[1]*dv11*wv1; 
  alpha2[2] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha2[4] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha2[7] = 2.0*E1[2]*dv11-2.0*B2[2]*dv11*wv1; 
  alpha2[11] = -0.5773502691896257*B2[2]*dv1*dv11; 
  const double amid1 = 0.3535533905932737*alpha1[0]-0.3952847075210473*alpha1[7]; 
  const double amid2 = 0.3535533905932737*alpha2[0]-0.3952847075210473*alpha2[7]; 
  out[1] += 0.6123724356957944*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha1[13]*f[13]+alpha1[7]*f[7]+alpha1[5]*f[5]+alpha1[3]*f[3]+alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[3] += 0.6123724356957944*(alpha2[11]*f[11]+alpha2[7]*f[7]+alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[4] += 0.5477225575051661*(alpha1[5]*f[13]+f[5]*alpha1[13]+alpha0[2]*f[8]+alpha1[1]*f[7]+f[1]*alpha1[7])+0.6123724356957944*(alpha1[3]*f[5]+f[3]*alpha1[5]+alpha0[0]*f[2]+f[0]*alpha0[2]+alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[5] += 0.5477225575051661*(alpha2[4]*f[11]+f[4]*alpha2[11]+alpha2[1]*f[7]+f[1]*alpha2[7])+0.6123724356957944*(alpha0[2]*f[6]+alpha2[2]*f[4]+f[2]*alpha2[4]+alpha0[0]*f[3]+alpha2[0]*f[1]+f[0]*alpha2[1]); 
  out[6] += 0.5477225575051661*alpha1[5]*f[15]+0.6123724356957944*(alpha1[7]*f[13]+f[7]*alpha1[13])+0.5477225575051661*alpha2[4]*f[12]+0.6123724356957944*(alpha2[7]*f[11]+f[7]*alpha2[11])+0.5477225575051661*(alpha1[3]*f[9]+alpha2[2]*f[8])+0.6123724356957944*(alpha1[1]*f[5]+f[1]*alpha1[5]+alpha2[1]*f[4]+f[1]*alpha2[4]+alpha1[0]*f[3]+f[0]*alpha1[3]+alpha2[0]*f[2]+f[0]*alpha2[2]); 
  out[7] += 1.369306393762915*(alpha0[2]*f[4]+alpha0[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha1[13]*f[17]+alpha1[7]*f[11]+alpha1[5]*f[10]+alpha1[3]*f[6]+alpha1[1]*f[4]+alpha1[0]*f[2]); 
  out[9] += 1.369306393762915*(alpha2[11]*f[17]+alpha2[7]*f[13]+alpha2[4]*f[10]+alpha2[2]*f[6]+alpha2[1]*f[5]+alpha2[0]*f[3]); 
  out[10] += 0.4898979485566357*alpha1[13]*f[15]+0.5477225575051661*(alpha1[3]*f[15]+alpha0[2]*f[14]+alpha1[1]*f[13]+f[1]*alpha1[13])+0.4898979485566357*alpha2[11]*f[12]+0.5477225575051661*(alpha2[2]*f[12]+alpha2[1]*f[11]+f[1]*alpha2[11]+alpha1[5]*f[9]+alpha2[4]*f[8]+(alpha1[5]+alpha2[4])*f[7]+f[4]*alpha2[7]+f[5]*alpha1[7])+0.6123724356957944*(alpha0[0]*f[6]+alpha1[0]*f[5]+f[0]*alpha1[5]+alpha2[0]*f[4]+f[0]*alpha2[4]+(alpha0[2]+alpha1[1])*f[3]+f[1]*alpha1[3]+alpha2[1]*f[2]+f[1]*alpha2[2]); 
  out[11] += 0.3912303982179757*alpha1[13]*f[13]+0.6123724356957944*(alpha1[3]*f[13]+f[3]*alpha1[13])+1.224744871391589*alpha0[2]*f[12]+0.3912303982179757*alpha1[7]*f[7]+0.6123724356957944*(alpha1[0]*f[7]+f[0]*alpha1[7])+0.5477225575051661*alpha1[5]*f[5]+1.369306393762915*alpha0[0]*f[4]+f[1]*(1.369306393762915*alpha0[2]+0.5477225575051661*alpha1[1]); 
  out[12] += 1.224744871391589*(alpha1[5]*f[17]+f[10]*alpha1[13]+alpha1[1]*f[11])+1.369306393762915*alpha1[3]*f[10]+0.6123724356957944*alpha0[0]*f[8]+1.224744871391589*f[4]*alpha1[7]+1.369306393762915*(alpha1[5]*f[6]+alpha1[0]*f[4])+(0.5477225575051661*alpha0[2]+1.369306393762915*alpha1[1])*f[2]; 
  out[13] += 0.3912303982179757*alpha2[11]*f[11]+0.6123724356957944*(alpha2[2]*f[11]+f[2]*alpha2[11])+1.369306393762915*alpha0[2]*f[10]+0.3912303982179757*alpha2[7]*f[7]+0.6123724356957944*(alpha2[0]*f[7]+f[0]*alpha2[7])+1.369306393762915*alpha0[0]*f[5]+0.5477225575051661*(alpha2[4]*f[4]+alpha2[1]*f[1]); 
  out[14] += 1.224744871391589*alpha1[5]*f[19]+1.369306393762915*alpha1[7]*f[17]+1.224744871391589*alpha1[3]*f[16]+1.369306393762915*f[11]*alpha1[13]+0.6123724356957944*alpha2[1]*f[12]+0.5477225575051661*alpha2[11]*f[11]+1.369306393762915*alpha1[1]*f[10]+0.6123724356957944*alpha2[0]*f[8]+1.369306393762915*alpha1[0]*f[6]+f[4]*(1.369306393762915*alpha1[5]+0.5477225575051661*alpha2[4])+f[2]*(1.369306393762915*alpha1[3]+0.5477225575051661*alpha2[2]); 
  out[15] += 1.224744871391589*alpha2[4]*f[17]+0.6123724356957944*alpha0[2]*f[16]+1.224744871391589*alpha2[1]*f[13]+f[10]*(1.224744871391589*alpha2[11]+1.369306393762915*alpha2[2])+0.6123724356957944*alpha0[0]*f[9]+1.224744871391589*f[5]*alpha2[7]+1.369306393762915*(alpha2[4]*f[6]+alpha2[0]*f[5]+alpha2[1]*f[3]); 
  out[16] += 1.224744871391589*alpha2[4]*f[18]+1.369306393762915*alpha2[7]*f[17]+0.6123724356957944*alpha1[1]*f[15]+1.224744871391589*alpha2[2]*f[14]+0.5477225575051661*alpha1[13]*f[13]+1.369306393762915*(alpha2[11]*f[13]+alpha2[1]*f[10])+0.6123724356957944*alpha1[0]*f[9]+1.369306393762915*alpha2[0]*f[6]+(0.5477225575051661*alpha1[5]+1.369306393762915*alpha2[4])*f[5]+(0.5477225575051661*alpha1[3]+1.369306393762915*alpha2[2])*f[3]; 
  out[17] += 1.224744871391589*alpha0[2]*f[18]+0.4898979485566357*alpha1[5]*f[15]+(0.3912303982179757*alpha1[7]+0.6123724356957944*alpha1[0])*f[13]+(0.5477225575051661*f[9]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha1[13]+0.4898979485566357*alpha2[4]*f[12]+(0.3912303982179757*alpha2[7]+0.6123724356957944*alpha2[0])*f[11]+(0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha2[11]+1.369306393762915*alpha0[0]*f[10]+0.6123724356957944*((alpha1[3]+alpha2[2])*f[7]+f[2]*alpha2[7]+f[3]*alpha1[7])+1.369306393762915*alpha0[2]*f[5]+0.5477225575051661*(alpha1[1]*f[5]+f[1]*alpha1[5]+alpha2[1]*f[4]+f[1]*alpha2[4]); 
  out[18] += 1.095445115010332*alpha1[13]*f[19]+1.224744871391589*(alpha1[3]*f[19]+alpha1[1]*f[17]+alpha1[5]*f[16])+0.6123724356957944*alpha0[0]*f[14]+1.224744871391589*f[4]*alpha1[13]+(0.5477225575051661*alpha2[7]+0.6123724356957944*alpha2[0])*f[12]+1.224744871391589*alpha1[5]*f[11]+0.4898979485566357*(alpha2[4]*f[11]+f[4]*alpha2[11])+(1.224744871391589*alpha1[7]+1.369306393762915*alpha1[0])*f[10]+0.6123724356957944*alpha2[1]*f[8]+0.5477225575051661*alpha0[2]*f[6]+1.369306393762915*(alpha1[1]*f[6]+f[2]*alpha1[5]+alpha1[3]*f[4])+0.5477225575051661*(alpha2[2]*f[4]+f[2]*alpha2[4]); 
  out[19] += 1.095445115010332*alpha2[11]*f[18]+1.224744871391589*(alpha2[2]*f[18]+alpha2[1]*f[17])+0.6123724356957944*alpha0[0]*f[16]+(0.5477225575051661*alpha1[7]+0.6123724356957944*alpha1[0])*f[15]+1.224744871391589*alpha2[4]*f[14]+(0.4898979485566357*alpha1[5]+1.224744871391589*alpha2[4])*f[13]+f[5]*(0.4898979485566357*alpha1[13]+1.224744871391589*alpha2[11])+(1.224744871391589*alpha2[7]+1.369306393762915*alpha2[0])*f[10]+0.6123724356957944*(alpha0[2]+alpha1[1])*f[9]+1.369306393762915*alpha2[1]*f[6]+(0.5477225575051661*alpha1[3]+1.369306393762915*alpha2[2])*f[5]+f[3]*(0.5477225575051661*alpha1[5]+1.369306393762915*alpha2[4]); 
return std::abs(w0dx0)+0.5*(dv0dx0+std::abs(amid1)+std::abs(amid2)); 
} 
double VlasovVol1x2vSerP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *E1 = &EM[4]; 
  const double dv2 = dxv[2], wv2 = w[2]; 

  const double *B0 = &EM[12]; 
  const double *B1 = &EM[16]; 
  const double *B2 = &EM[20]; 

  double alpha0[32]; 

  double alpha1[32]; 

  double alpha2[32]; 

  alpha0[0] = 5.656854249492382*w0dx0; 
  alpha0[2] = 1.632993161855453*dv0dx0; 

  alpha1[0] = dv10*(2.0*B2[0]*wv2+2.0*E0[0]); 
  alpha1[1] = dv10*(2.0*B2[1]*wv2+2.0*E0[1]); 
  alpha1[3] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha1[5] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha1[7] = dv10*(2.0*B2[2]*wv2+2.0*E0[2]); 
  alpha1[13] = 0.5773502691896257*B2[2]*dv10*dv2; 
  alpha1[17] = dv10*(2.0*B2[3]*wv2+2.0*E0[3]); 
  alpha1[25] = 0.5773502691896256*B2[3]*dv10*dv2; 
  alpha2[0] = 2.0*E1[0]*dv11-2.0*B2[0]*dv11*wv1; 
  alpha2[1] = 2.0*E1[1]*dv11-2.0*B2[1]*dv11*wv1; 
  alpha2[2] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha2[4] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha2[7] = 2.0*E1[2]*dv11-2.0*B2[2]*dv11*wv1; 
  alpha2[11] = -0.5773502691896257*B2[2]*dv1*dv11; 
  alpha2[17] = 2.0*E1[3]*dv11-2.0*B2[3]*dv11*wv1; 
  alpha2[23] = -0.5773502691896256*B2[3]*dv1*dv11; 
  const double amid1 = 0.3535533905932737*alpha1[0]-0.3952847075210473*alpha1[7]; 
  const double amid2 = 0.3535533905932737*alpha2[0]-0.3952847075210473*alpha2[7]; 
  out[1] += 0.6123724356957944*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[2] += 0.6123724356957944*(alpha1[25]*f[25]+alpha1[17]*f[17]+alpha1[13]*f[13]+alpha1[7]*f[7]+alpha1[5]*f[5]+alpha1[3]*f[3]+alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[3] += 0.6123724356957944*(alpha2[23]*f[23]+alpha2[17]*f[17]+alpha2[11]*f[11]+alpha2[7]*f[7]+alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[4] += 0.537852874200477*(alpha1[13]*f[25]+f[13]*alpha1[25]+alpha1[7]*f[17]+f[7]*alpha1[17])+0.5477225575051661*(alpha1[5]*f[13]+f[5]*alpha1[13]+alpha0[2]*f[8]+alpha1[1]*f[7]+f[1]*alpha1[7])+0.6123724356957944*(alpha1[3]*f[5]+f[3]*alpha1[5]+alpha0[0]*f[2]+f[0]*alpha0[2]+alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[5] += 0.537852874200477*(alpha2[11]*f[23]+f[11]*alpha2[23]+alpha2[7]*f[17]+f[7]*alpha2[17])+0.5477225575051661*(alpha2[4]*f[11]+f[4]*alpha2[11]+alpha2[1]*f[7]+f[1]*alpha2[7])+0.6123724356957944*(alpha0[2]*f[6]+alpha2[2]*f[4]+f[2]*alpha2[4]+alpha0[0]*f[3]+alpha2[0]*f[1]+f[0]*alpha2[1]); 
  out[6] += 0.6123724356957944*(alpha1[17]*f[25]+f[17]*alpha1[25]+alpha2[17]*f[23]+f[17]*alpha2[23])+0.5477225575051661*alpha1[5]*f[15]+0.6123724356957944*(alpha1[7]*f[13]+f[7]*alpha1[13])+0.5477225575051661*alpha2[4]*f[12]+0.6123724356957944*(alpha2[7]*f[11]+f[7]*alpha2[11])+0.5477225575051661*(alpha1[3]*f[9]+alpha2[2]*f[8])+0.6123724356957944*(alpha1[1]*f[5]+f[1]*alpha1[5]+alpha2[1]*f[4]+f[1]*alpha2[4]+alpha1[0]*f[3]+f[0]*alpha1[3]+alpha2[0]*f[2]+f[0]*alpha2[2]); 
  out[7] += 1.369306393762915*(alpha0[2]*f[4]+alpha0[0]*f[1]); 
  out[8] += 1.369306393762915*(alpha1[25]*f[29]+alpha1[17]*f[23]+alpha1[13]*f[20]+alpha1[7]*f[11]+alpha1[5]*f[10]+alpha1[3]*f[6]+alpha1[1]*f[4]+alpha1[0]*f[2]); 
  out[9] += 1.369306393762915*(alpha2[23]*f[29]+alpha2[17]*f[25]+alpha2[11]*f[20]+alpha2[7]*f[13]+alpha2[4]*f[10]+alpha2[2]*f[6]+alpha2[1]*f[5]+alpha2[0]*f[3]); 
  out[10] += 0.537852874200477*(alpha1[7]*f[25]+f[7]*alpha1[25]+alpha2[7]*f[23]+f[7]*alpha2[23]+(alpha1[13]+alpha2[11])*f[17]+f[11]*alpha2[17]+f[13]*alpha1[17])+0.4898979485566357*alpha1[13]*f[15]+0.5477225575051661*(alpha1[3]*f[15]+alpha0[2]*f[14]+alpha1[1]*f[13]+f[1]*alpha1[13])+0.4898979485566357*alpha2[11]*f[12]+0.5477225575051661*(alpha2[2]*f[12]+alpha2[1]*f[11]+f[1]*alpha2[11]+alpha1[5]*f[9]+alpha2[4]*f[8]+(alpha1[5]+alpha2[4])*f[7]+f[4]*alpha2[7]+f[5]*alpha1[7])+0.6123724356957944*(alpha0[0]*f[6]+alpha1[0]*f[5]+f[0]*alpha1[5]+alpha2[0]*f[4]+f[0]*alpha2[4]+(alpha0[2]+alpha1[1])*f[3]+f[1]*alpha1[3]+alpha2[1]*f[2]+f[1]*alpha2[2]); 
  out[11] += 0.3651483716701108*alpha1[25]*f[25]+0.537852874200477*(alpha1[5]*f[25]+f[5]*alpha1[25])+0.3651483716701108*alpha1[17]*f[17]+0.537852874200477*(alpha1[1]*f[17]+f[1]*alpha1[17])+0.3912303982179757*alpha1[13]*f[13]+0.6123724356957944*(alpha1[3]*f[13]+f[3]*alpha1[13])+1.224744871391589*alpha0[2]*f[12]+0.3912303982179757*alpha1[7]*f[7]+0.6123724356957944*(alpha1[0]*f[7]+f[0]*alpha1[7])+0.5477225575051661*alpha1[5]*f[5]+1.369306393762915*alpha0[0]*f[4]+f[1]*(1.369306393762915*alpha0[2]+0.5477225575051661*alpha1[1]); 
  out[12] += 1.202675588605909*(alpha1[13]*f[29]+f[20]*alpha1[25]+alpha1[7]*f[23])+1.224744871391589*alpha1[5]*f[20]+0.537852874200477*alpha0[2]*f[18]+1.202675588605909*f[11]*alpha1[17]+1.224744871391589*(f[10]*alpha1[13]+alpha1[1]*f[11])+1.369306393762915*alpha1[3]*f[10]+0.6123724356957944*alpha0[0]*f[8]+1.224744871391589*f[4]*alpha1[7]+1.369306393762915*(alpha1[5]*f[6]+alpha1[0]*f[4])+(0.5477225575051661*alpha0[2]+1.369306393762915*alpha1[1])*f[2]; 
  out[13] += 0.3651483716701108*alpha2[23]*f[23]+0.537852874200477*(alpha2[4]*f[23]+f[4]*alpha2[23])+0.3651483716701108*alpha2[17]*f[17]+0.537852874200477*(alpha2[1]*f[17]+f[1]*alpha2[17])+0.3912303982179757*alpha2[11]*f[11]+0.6123724356957944*(alpha2[2]*f[11]+f[2]*alpha2[11])+1.369306393762915*alpha0[2]*f[10]+0.3912303982179757*alpha2[7]*f[7]+0.6123724356957944*(alpha2[0]*f[7]+f[0]*alpha2[7])+1.369306393762915*alpha0[0]*f[5]+0.5477225575051661*(alpha2[4]*f[4]+alpha2[1]*f[1]); 
  out[14] += 1.369306393762915*(alpha1[17]*f[29]+f[23]*alpha1[25])+0.537852874200477*alpha2[4]*f[24]+0.5477225575051661*alpha2[23]*f[23]+1.224744871391589*alpha1[5]*f[22]+1.369306393762915*alpha1[7]*f[20]+0.537852874200477*alpha2[2]*f[18]+1.224744871391589*alpha1[3]*f[16]+1.369306393762915*f[11]*alpha1[13]+0.6123724356957944*alpha2[1]*f[12]+0.5477225575051661*alpha2[11]*f[11]+1.369306393762915*alpha1[1]*f[10]+0.6123724356957944*alpha2[0]*f[8]+1.369306393762915*alpha1[0]*f[6]+f[4]*(1.369306393762915*alpha1[5]+0.5477225575051661*alpha2[4])+f[2]*(1.369306393762915*alpha1[3]+0.5477225575051661*alpha2[2]); 
  out[15] += 1.202675588605909*(alpha2[11]*f[29]+alpha2[7]*f[25])+f[20]*(1.202675588605909*alpha2[23]+1.224744871391589*alpha2[4])+1.202675588605909*f[13]*alpha2[17]+0.6123724356957944*alpha0[2]*f[16]+1.224744871391589*alpha2[1]*f[13]+f[10]*(1.224744871391589*alpha2[11]+1.369306393762915*alpha2[2])+0.6123724356957944*alpha0[0]*f[9]+1.224744871391589*f[5]*alpha2[7]+1.369306393762915*(alpha2[4]*f[6]+alpha2[0]*f[5]+alpha2[1]*f[3]); 
  out[16] += 1.369306393762915*alpha2[17]*f[29]+0.537852874200477*alpha1[5]*f[27]+(0.5477225575051661*alpha1[25]+1.369306393762915*alpha2[23])*f[25]+1.224744871391589*alpha2[4]*f[21]+1.369306393762915*alpha2[7]*f[20]+0.537852874200477*alpha1[3]*f[19]+0.6123724356957944*alpha1[1]*f[15]+1.224744871391589*alpha2[2]*f[14]+0.5477225575051661*alpha1[13]*f[13]+1.369306393762915*(alpha2[11]*f[13]+alpha2[1]*f[10])+0.6123724356957944*alpha1[0]*f[9]+1.369306393762915*alpha2[0]*f[6]+(0.5477225575051661*alpha1[5]+1.369306393762915*alpha2[4])*f[5]+(0.5477225575051661*alpha1[3]+1.369306393762915*alpha2[2])*f[3]; 
  out[17] += 2.091650066335188*(alpha0[2]*f[11]+alpha0[0]*f[7])+0.9354143466934851*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[18] += 0.9354143466934851*alpha1[25]*f[25]+2.091650066335188*alpha1[5]*f[21]+0.9354143466934851*alpha1[17]*f[17]+2.091650066335188*alpha1[3]*f[14]+0.9354143466934851*alpha1[13]*f[13]+2.091650066335188*(alpha1[1]*f[12]+alpha1[0]*f[8])+0.9354143466934851*(alpha1[7]*f[7]+alpha1[5]*f[5]+alpha1[3]*f[3]+alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[19] += 0.9354143466934851*alpha2[23]*f[23]+2.091650066335188*alpha2[4]*f[22]+0.9354143466934851*alpha2[17]*f[17]+2.091650066335188*(alpha2[2]*f[16]+alpha2[1]*f[15])+0.9354143466934851*alpha2[11]*f[11]+2.091650066335188*alpha2[0]*f[9]+0.9354143466934851*(alpha2[7]*f[7]+alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[20] += (0.3651483716701108*alpha1[17]+0.537852874200477*alpha1[1])*f[25]+(0.3651483716701108*f[17]+0.4810702354423638*f[15]+0.537852874200477*f[1])*alpha1[25]+(0.3651483716701108*alpha2[17]+0.537852874200477*alpha2[1])*f[23]+(0.3651483716701108*f[17]+0.4810702354423638*f[12]+0.537852874200477*f[1])*alpha2[23]+1.224744871391589*alpha0[2]*f[21]+0.537852874200477*((alpha1[5]+alpha2[4])*f[17]+f[4]*alpha2[17]+f[5]*alpha1[17])+0.4898979485566357*alpha1[5]*f[15]+(0.3912303982179757*alpha1[7]+0.6123724356957944*alpha1[0])*f[13]+(0.5477225575051661*f[9]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha1[13]+0.4898979485566357*alpha2[4]*f[12]+(0.3912303982179757*alpha2[7]+0.6123724356957944*alpha2[0])*f[11]+(0.5477225575051661*f[8]+0.3912303982179757*f[7]+0.6123724356957944*f[0])*alpha2[11]+1.369306393762915*alpha0[0]*f[10]+0.6123724356957944*((alpha1[3]+alpha2[2])*f[7]+f[2]*alpha2[7]+f[3]*alpha1[7])+1.369306393762915*alpha0[2]*f[5]+0.5477225575051661*(alpha1[1]*f[5]+f[1]*alpha1[5]+alpha2[1]*f[4]+f[1]*alpha2[4]); 
  out[21] += 1.202675588605909*alpha1[7]*f[29]+0.537852874200477*alpha0[2]*f[26]+1.202675588605909*f[11]*alpha1[25]+(0.4810702354423638*alpha2[11]+0.537852874200477*alpha2[2])*f[24]+1.202675588605909*alpha1[13]*f[23]+0.4810702354423638*(alpha2[11]*f[23]+f[11]*alpha2[23])+(1.095445115010332*alpha1[13]+1.224744871391589*alpha1[3])*f[22]+(1.202675588605909*alpha1[17]+1.224744871391589*alpha1[1])*f[20]+0.537852874200477*alpha2[4]*f[18]+1.224744871391589*alpha1[5]*f[16]+0.6123724356957944*alpha0[0]*f[14]+1.224744871391589*f[4]*alpha1[13]+(0.5477225575051661*alpha2[7]+0.6123724356957944*alpha2[0])*f[12]+1.224744871391589*alpha1[5]*f[11]+0.4898979485566357*(alpha2[4]*f[11]+f[4]*alpha2[11])+(1.224744871391589*alpha1[7]+1.369306393762915*alpha1[0])*f[10]+0.6123724356957944*alpha2[1]*f[8]+0.5477225575051661*alpha0[2]*f[6]+1.369306393762915*(alpha1[1]*f[6]+f[2]*alpha1[5]+alpha1[3]*f[4])+0.5477225575051661*(alpha2[2]*f[4]+f[2]*alpha2[4]); 
  out[22] += 1.202675588605909*alpha2[7]*f[29]+(0.4810702354423638*alpha1[13]+0.537852874200477*alpha1[3])*f[27]+(0.4810702354423638*alpha1[13]+1.202675588605909*alpha2[11])*f[25]+f[13]*(0.4810702354423638*alpha1[25]+1.202675588605909*alpha2[23])+(1.095445115010332*alpha2[11]+1.224744871391589*alpha2[2])*f[21]+(1.202675588605909*alpha2[17]+1.224744871391589*alpha2[1])*f[20]+0.537852874200477*alpha1[5]*f[19]+0.6123724356957944*alpha0[0]*f[16]+(0.5477225575051661*alpha1[7]+0.6123724356957944*alpha1[0])*f[15]+1.224744871391589*alpha2[4]*f[14]+(0.4898979485566357*alpha1[5]+1.224744871391589*alpha2[4])*f[13]+f[5]*(0.4898979485566357*alpha1[13]+1.224744871391589*alpha2[11])+(1.224744871391589*alpha2[7]+1.369306393762915*alpha2[0])*f[10]+0.6123724356957944*(alpha0[2]+alpha1[1])*f[9]+1.369306393762915*alpha2[1]*f[6]+(0.5477225575051661*alpha1[3]+1.369306393762915*alpha2[2])*f[5]+f[3]*(0.5477225575051661*alpha1[5]+1.369306393762915*alpha2[4]); 
  out[23] += (0.3651483716701108*alpha1[13]+0.6123724356957944*alpha1[3])*f[25]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alpha1[25]+(0.3651483716701108*alpha1[7]+0.6123724356957944*alpha1[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha1[17]+0.537852874200477*(alpha1[5]*f[13]+f[5]*alpha1[13])+2.091650066335188*alpha0[0]*f[11]+alpha0[2]*(0.8366600265340755*f[8]+2.091650066335188*f[7])+0.537852874200477*(alpha1[1]*f[7]+f[1]*alpha1[7])+0.9354143466934851*(alpha0[0]*f[2]+f[0]*alpha0[2]); 
  out[24] += 0.8215838362577489*(alpha1[13]*f[25]+f[13]*alpha1[25])+(1.870828693386971*alpha1[13]+2.091650066335188*alpha1[3])*f[21]+0.6123724356957944*alpha0[0]*f[18]+0.8215838362577489*(alpha1[7]*f[17]+f[7]*alpha1[17])+2.091650066335188*alpha1[5]*f[14]+0.8366600265340755*(alpha1[5]*f[13]+f[5]*alpha1[13])+(1.870828693386971*alpha1[7]+2.091650066335188*alpha1[0])*f[12]+(0.537852874200477*alpha0[2]+2.091650066335188*alpha1[1])*f[8]+0.8366600265340755*(alpha1[1]*f[7]+f[1]*alpha1[7])+0.9354143466934851*(alpha1[3]*f[5]+f[3]*alpha1[5]+alpha1[0]*f[1]+f[0]*alpha1[1]); 
  out[25] += (0.3651483716701108*alpha2[11]+0.6123724356957944*alpha2[2])*f[23]+(0.3651483716701108*f[11]+0.6123724356957944*f[2])*alpha2[23]+2.091650066335188*alpha0[2]*f[20]+(0.3651483716701108*alpha2[7]+0.6123724356957944*alpha2[0])*f[17]+(0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha2[17]+2.091650066335188*alpha0[0]*f[13]+0.537852874200477*(alpha2[4]*f[11]+f[4]*alpha2[11]+alpha2[1]*f[7]+f[1]*alpha2[7])+0.9354143466934851*(alpha0[2]*f[6]+alpha0[0]*f[3]); 
  out[26] += 0.9354143466934851*(alpha1[17]*f[25]+f[17]*alpha1[25])+0.6123724356957944*alpha2[1]*f[24]+2.091650066335188*alpha1[1]*f[21]+0.6123724356957944*alpha2[0]*f[18]+0.8366600265340755*alpha1[5]*f[15]+2.091650066335188*alpha1[0]*f[14]+0.9354143466934851*(alpha1[7]*f[13]+f[7]*alpha1[13])+(2.091650066335188*alpha1[5]+0.537852874200477*alpha2[4])*f[12]+0.8366600265340755*alpha1[3]*f[9]+(2.091650066335188*alpha1[3]+0.537852874200477*alpha2[2])*f[8]+0.9354143466934851*(alpha1[1]*f[5]+f[1]*alpha1[5]+alpha1[0]*f[3]+f[0]*alpha1[3]); 
  out[27] += 0.6123724356957944*alpha0[2]*f[28]+0.8215838362577489*(alpha2[11]*f[23]+f[11]*alpha2[23])+(1.870828693386971*alpha2[11]+2.091650066335188*alpha2[2])*f[22]+0.6123724356957944*alpha0[0]*f[19]+0.8215838362577489*(alpha2[7]*f[17]+f[7]*alpha2[17])+2.091650066335188*alpha2[4]*f[16]+(1.870828693386971*alpha2[7]+2.091650066335188*alpha2[0])*f[15]+0.8366600265340755*(alpha2[4]*f[11]+f[4]*alpha2[11])+2.091650066335188*alpha2[1]*f[9]+0.8366600265340755*(alpha2[1]*f[7]+f[1]*alpha2[7])+0.9354143466934851*(alpha2[2]*f[4]+f[2]*alpha2[4]+alpha2[0]*f[1]+f[0]*alpha2[1]); 
  out[28] += 0.6123724356957944*alpha1[1]*f[27]+0.9354143466934851*(alpha2[17]*f[23]+f[17]*alpha2[23])+2.091650066335188*alpha2[1]*f[22]+0.6123724356957944*alpha1[0]*f[19]+2.091650066335188*alpha2[0]*f[16]+0.537852874200477*alpha1[5]*f[15]+alpha2[4]*(2.091650066335188*f[15]+0.8366600265340755*f[12])+0.9354143466934851*(alpha2[7]*f[11]+f[7]*alpha2[11])+0.537852874200477*alpha1[3]*f[9]+alpha2[2]*(2.091650066335188*f[9]+0.8366600265340755*f[8])+0.9354143466934851*(alpha2[1]*f[4]+f[1]*alpha2[4]+alpha2[0]*f[2]+f[0]*alpha2[2]); 
  out[29] += (0.3651483716701108*alpha1[7]+0.6123724356957944*alpha1[0])*f[25]+(0.5477225575051661*f[9]+0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha1[25]+(0.3651483716701108*alpha2[7]+0.6123724356957944*alpha2[0])*f[23]+(0.5477225575051661*f[8]+0.3651483716701108*f[7]+0.6123724356957944*f[0])*alpha2[23]+2.091650066335188*alpha0[0]*f[20]+(0.3651483716701108*(alpha1[13]+alpha2[11])+0.6123724356957944*(alpha1[3]+alpha2[2]))*f[17]+(0.3651483716701108*f[11]+0.6123724356957944*f[2])*alpha2[17]+(0.3651483716701108*f[13]+0.6123724356957944*f[3])*alpha1[17]+0.4810702354423638*alpha1[13]*f[15]+alpha0[2]*(0.8366600265340755*f[14]+2.091650066335188*f[13])+0.537852874200477*(alpha1[1]*f[13]+f[1]*alpha1[13])+0.4810702354423638*alpha2[11]*f[12]+0.537852874200477*(alpha2[1]*f[11]+f[1]*alpha2[11]+(alpha1[5]+alpha2[4])*f[7]+f[4]*alpha2[7]+f[5]*alpha1[7])+0.9354143466934851*(alpha0[0]*f[6]+alpha0[2]*f[3]); 
  out[30] += 0.6123724356957944*alpha0[0]*f[26]+0.8215838362577489*(alpha1[7]*f[25]+f[7]*alpha1[25])+(0.5477225575051661*alpha2[7]+0.6123724356957944*alpha2[0])*f[24]+(1.870828693386971*alpha1[7]+2.091650066335188*alpha1[0])*f[21]+0.6123724356957944*alpha2[1]*f[18]+0.8215838362577489*(alpha1[13]*f[17]+f[13]*alpha1[17])+(0.7483314773547884*alpha1[13]+0.8366600265340755*alpha1[3])*f[15]+0.537852874200477*alpha0[2]*f[14]+alpha1[1]*(2.091650066335188*f[14]+0.8366600265340755*f[13])+(1.870828693386971*f[12]+0.8366600265340755*f[1])*alpha1[13]+(0.4810702354423638*alpha2[11]+2.091650066335188*alpha1[3]+0.537852874200477*alpha2[2])*f[12]+0.8366600265340755*alpha1[5]*f[9]+(2.091650066335188*alpha1[5]+0.537852874200477*alpha2[4])*f[8]+0.8366600265340755*(alpha1[5]*f[7]+f[5]*alpha1[7])+0.9354143466934851*(alpha1[0]*f[5]+f[0]*alpha1[5]+alpha1[1]*f[3]+f[1]*alpha1[3]); 
  out[31] += 0.6123724356957944*alpha0[0]*f[28]+(0.5477225575051661*alpha1[7]+0.6123724356957944*alpha1[0])*f[27]+0.8215838362577489*(alpha2[7]*f[23]+f[7]*alpha2[23])+(1.870828693386971*alpha2[7]+2.091650066335188*alpha2[0])*f[22]+0.6123724356957944*(alpha0[2]+alpha1[1])*f[19]+0.8215838362577489*(alpha2[11]*f[17]+f[11]*alpha2[17])+2.091650066335188*alpha2[1]*f[16]+(0.4810702354423638*alpha1[13]+1.870828693386971*alpha2[11]+0.537852874200477*alpha1[3]+2.091650066335188*alpha2[2])*f[15]+0.7483314773547884*alpha2[11]*f[12]+0.8366600265340755*(alpha2[2]*f[12]+alpha2[1]*f[11]+f[1]*alpha2[11])+(0.537852874200477*alpha1[5]+2.091650066335188*alpha2[4])*f[9]+0.8366600265340755*(alpha2[4]*(f[8]+f[7])+f[4]*alpha2[7])+0.9354143466934851*(alpha2[0]*f[4]+f[0]*alpha2[4]+alpha2[1]*f[2]+f[1]*alpha2[2]); 
return std::abs(w0dx0)+0.5*(dv0dx0+std::abs(amid1)+std::abs(amid2)); 
} 
