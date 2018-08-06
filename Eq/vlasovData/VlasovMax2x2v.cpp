#include <VlasovModDecl.h> 
double VlasovVol2x2vMaxP1(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[3]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[9]; 
  const double *B1 = &EM[12]; 
  const double *B2 = &EM[15]; 

  double alpha0[5]; 

  double alpha1[5]; 

  double alpha2[5]; 

  double alpha3[5]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[3] = 2.309401076758503*dv0dx0; 

  alpha1[0] = 8.0*w1dx1; 
  alpha1[4] = 2.309401076758503*dv1dx1; 

  alpha2[0] = dv10*(2.0*B2[0]*wv2+2.0*E0[0]); 
  alpha2[1] = dv10*(2.0*B2[1]*wv2+2.0*E0[1]); 
  alpha2[2] = dv10*(2.0*B2[2]*wv2+2.0*E0[2]); 
  alpha2[4] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha3[0] = 2.0*E1[0]*dv11-2.0*B2[0]*dv11*wv1; 
  alpha3[1] = 2.0*E1[1]*dv11-2.0*B2[1]*dv11*wv1; 
  alpha3[2] = 2.0*E1[2]*dv11-2.0*B2[2]*dv11*wv1; 
  alpha3[3] = -0.5773502691896258*B2[0]*dv1*dv11; 
  const double amid1 = 0.25*alpha2[0]; 
  const double amid2 = 0.25*alpha3[0]; 
  out[1] += 0.4330127018922193*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[3] += 0.4330127018922193*(alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[4] += 0.4330127018922193*(alpha3[3]*f[3]+alpha3[2]*f[2]+alpha3[1]*f[1]+alpha3[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1+std::abs(amid1)+std::abs(amid2)); 
} 
double VlasovVol2x2vMaxP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[6]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[18]; 
  const double *B1 = &EM[24]; 
  const double *B2 = &EM[30]; 

  double alpha0[15]; 

  double alpha1[15]; 

  double alpha2[15]; 

  double alpha3[15]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[3] = 2.309401076758503*dv0dx0; 

  alpha1[0] = 8.0*w1dx1; 
  alpha1[4] = 2.309401076758503*dv1dx1; 

  alpha2[0] = dv10*(2.0*B2[0]*wv2+2.0*E0[0]); 
  alpha2[1] = dv10*(2.0*B2[1]*wv2+2.0*E0[1]); 
  alpha2[2] = dv10*(2.0*B2[2]*wv2+2.0*E0[2]); 
  alpha2[4] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha2[5] = dv10*(2.0*B2[3]*wv2+2.0*E0[3]); 
  alpha2[8] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha2[9] = 0.5773502691896258*B2[2]*dv10*dv2; 
  alpha2[11] = dv10*(2.0*B2[4]*wv2+2.0*E0[4]); 
  alpha2[12] = dv10*(2.0*B2[5]*wv2+2.0*E0[5]); 
  alpha3[0] = 2.0*E1[0]*dv11-2.0*B2[0]*dv11*wv1; 
  alpha3[1] = 2.0*E1[1]*dv11-2.0*B2[1]*dv11*wv1; 
  alpha3[2] = 2.0*E1[2]*dv11-2.0*B2[2]*dv11*wv1; 
  alpha3[3] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha3[5] = 2.0*E1[3]*dv11-2.0*B2[3]*dv11*wv1; 
  alpha3[6] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha3[7] = -0.5773502691896258*B2[2]*dv1*dv11; 
  alpha3[11] = 2.0*E1[4]*dv11-2.0*B2[4]*dv11*wv1; 
  alpha3[12] = 2.0*E1[5]*dv11-2.0*B2[5]*dv11*wv1; 
  const double amid1 = (-0.2795084971874737*alpha2[12])-0.2795084971874737*alpha2[11]+0.25*alpha2[0]; 
  const double amid2 = (-0.2795084971874737*alpha3[12])-0.2795084971874737*alpha3[11]+0.25*alpha3[0]; 
  out[1] += 0.4330127018922193*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[3] += 0.4330127018922193*(alpha2[12]*f[12]+alpha2[11]*f[11]+alpha2[9]*f[9]+alpha2[8]*f[8]+alpha2[5]*f[5]+alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[4] += 0.4330127018922193*(alpha3[12]*f[12]+alpha3[11]*f[11]+alpha3[7]*f[7]+alpha3[6]*f[6]+alpha3[5]*f[5]+alpha3[3]*f[3]+alpha3[2]*f[2]+alpha3[1]*f[1]+alpha3[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha1[4]*f[8]+alpha0[3]*f[7]+alpha0[0]*f[2]+alpha1[0]*f[1]); 
  out[6] += 0.3872983346207416*(alpha0[3]*f[13]+alpha2[1]*f[11]+f[1]*alpha2[11])+0.4330127018922193*(alpha2[4]*f[8]+f[4]*alpha2[8]+alpha2[2]*f[5]+f[2]*alpha2[5]+alpha0[0]*f[3]+f[0]*alpha0[3]+alpha2[0]*f[1]+f[0]*alpha2[1]); 
  out[7] += 0.3872983346207416*(alpha2[2]*f[12]+f[2]*alpha2[12])+0.4330127018922193*(alpha1[4]*f[10]+alpha2[4]*f[9]+f[4]*alpha2[9]+alpha2[1]*f[5]+f[1]*alpha2[5]+alpha1[0]*f[3]+alpha2[0]*f[2]+f[0]*alpha2[2]); 
  out[8] += 0.3872983346207416*(alpha3[1]*f[11]+f[1]*alpha3[11])+0.4330127018922193*(alpha0[3]*f[10]+alpha3[3]*f[6]+f[3]*alpha3[6]+alpha3[2]*f[5]+f[2]*alpha3[5]+alpha0[0]*f[4]+alpha3[0]*f[1]+f[0]*alpha3[1]); 
  out[9] += 0.3872983346207416*(alpha1[4]*f[14]+alpha3[2]*f[12]+f[2]*alpha3[12])+0.4330127018922193*(alpha3[3]*f[7]+f[3]*alpha3[7]+alpha3[1]*f[5]+f[1]*alpha3[5]+alpha1[0]*f[4]+f[0]*alpha1[4]+alpha3[0]*f[2]+f[0]*alpha3[2]); 
  out[10] += 0.3872983346207416*(alpha2[4]*f[14]+alpha3[3]*f[13])+0.4330127018922193*(alpha2[2]*f[9]+f[2]*alpha2[9]+alpha2[1]*f[8]+f[1]*alpha2[8]+alpha3[2]*f[7]+f[2]*alpha3[7]+alpha3[1]*f[6]+f[1]*alpha3[6]+alpha2[0]*f[4]+f[0]*alpha2[4]+alpha3[0]*f[3]+f[0]*alpha3[3]); 
  out[11] += 0.9682458365518543*(alpha0[3]*f[6]+alpha0[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha1[4]*f[9]+alpha1[0]*f[2]); 
  out[13] += 0.9682458365518543*(alpha2[4]*f[10]+alpha2[2]*f[7]+alpha2[1]*f[6]+alpha2[0]*f[3]); 
  out[14] += 0.9682458365518543*(alpha3[3]*f[10]+alpha3[2]*f[9]+alpha3[1]*f[8]+alpha3[0]*f[4]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1+std::abs(amid1)+std::abs(amid2)); 
} 
double VlasovVol2x2vMaxP3(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. EM/f: Input EM-field/distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  const double dv10 = 2/dxv[2]; 
  const double *E0 = &EM[0]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv11 = 2/dxv[3]; 
  const double *E1 = &EM[10]; 
  const double dv2 = dxv[3], wv2 = w[3]; 

  const double *B0 = &EM[30]; 
  const double *B1 = &EM[40]; 
  const double *B2 = &EM[50]; 

  double alpha0[35]; 

  double alpha1[35]; 

  double alpha2[35]; 

  double alpha3[35]; 

  alpha0[0] = 8.0*w0dx0; 
  alpha0[3] = 2.309401076758503*dv0dx0; 

  alpha1[0] = 8.0*w1dx1; 
  alpha1[4] = 2.309401076758503*dv1dx1; 

  alpha2[0] = dv10*(2.0*B2[0]*wv2+2.0*E0[0]); 
  alpha2[1] = dv10*(2.0*B2[1]*wv2+2.0*E0[1]); 
  alpha2[2] = dv10*(2.0*B2[2]*wv2+2.0*E0[2]); 
  alpha2[4] = 0.5773502691896258*B2[0]*dv10*dv2; 
  alpha2[5] = dv10*(2.0*B2[3]*wv2+2.0*E0[3]); 
  alpha2[8] = 0.5773502691896258*B2[1]*dv10*dv2; 
  alpha2[9] = 0.5773502691896258*B2[2]*dv10*dv2; 
  alpha2[11] = dv10*(2.0*B2[4]*wv2+2.0*E0[4]); 
  alpha2[12] = dv10*(2.0*B2[5]*wv2+2.0*E0[5]); 
  alpha2[16] = 0.5773502691896258*B2[3]*dv10*dv2; 
  alpha2[19] = dv10*(2.0*B2[6]*wv2+2.0*E0[6]); 
  alpha2[20] = dv10*(2.0*B2[7]*wv2+2.0*E0[7]); 
  alpha2[25] = 0.5773502691896257*B2[4]*dv10*dv2; 
  alpha2[26] = 0.5773502691896257*B2[5]*dv10*dv2; 
  alpha2[31] = dv10*(2.0*B2[8]*wv2+2.0*E0[8]); 
  alpha2[32] = dv10*(2.0*B2[9]*wv2+2.0*E0[9]); 
  alpha3[0] = 2.0*E1[0]*dv11-2.0*B2[0]*dv11*wv1; 
  alpha3[1] = 2.0*E1[1]*dv11-2.0*B2[1]*dv11*wv1; 
  alpha3[2] = 2.0*E1[2]*dv11-2.0*B2[2]*dv11*wv1; 
  alpha3[3] = -0.5773502691896258*B2[0]*dv1*dv11; 
  alpha3[5] = 2.0*E1[3]*dv11-2.0*B2[3]*dv11*wv1; 
  alpha3[6] = -0.5773502691896258*B2[1]*dv1*dv11; 
  alpha3[7] = -0.5773502691896258*B2[2]*dv1*dv11; 
  alpha3[11] = 2.0*E1[4]*dv11-2.0*B2[4]*dv11*wv1; 
  alpha3[12] = 2.0*E1[5]*dv11-2.0*B2[5]*dv11*wv1; 
  alpha3[15] = -0.5773502691896258*B2[3]*dv1*dv11; 
  alpha3[19] = 2.0*E1[6]*dv11-2.0*B2[6]*dv11*wv1; 
  alpha3[20] = 2.0*E1[7]*dv11-2.0*B2[7]*dv11*wv1; 
  alpha3[21] = -0.5773502691896257*B2[4]*dv1*dv11; 
  alpha3[22] = -0.5773502691896257*B2[5]*dv1*dv11; 
  alpha3[31] = 2.0*E1[8]*dv11-2.0*B2[8]*dv11*wv1; 
  alpha3[32] = 2.0*E1[9]*dv11-2.0*B2[9]*dv11*wv1; 
  const double amid1 = (-0.2795084971874737*alpha2[12])-0.2795084971874737*alpha2[11]+0.25*alpha2[0]; 
  const double amid2 = (-0.2795084971874737*alpha3[12])-0.2795084971874737*alpha3[11]+0.25*alpha3[0]; 
  out[1] += 0.4330127018922193*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[3] += 0.4330127018922193*(alpha2[32]*f[32]+alpha2[31]*f[31]+alpha2[26]*f[26]+alpha2[25]*f[25]+alpha2[20]*f[20]+alpha2[19]*f[19]+alpha2[16]*f[16]+alpha2[12]*f[12]+alpha2[11]*f[11]+alpha2[9]*f[9]+alpha2[8]*f[8]+alpha2[5]*f[5]+alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[4] += 0.4330127018922193*(alpha3[32]*f[32]+alpha3[31]*f[31]+alpha3[22]*f[22]+alpha3[21]*f[21]+alpha3[20]*f[20]+alpha3[19]*f[19]+alpha3[15]*f[15]+alpha3[12]*f[12]+alpha3[11]*f[11]+alpha3[7]*f[7]+alpha3[6]*f[6]+alpha3[5]*f[5]+alpha3[3]*f[3]+alpha3[2]*f[2]+alpha3[1]*f[1]+alpha3[0]*f[0]); 
  out[5] += 0.4330127018922193*(alpha1[4]*f[8]+alpha0[3]*f[7]+alpha0[0]*f[2]+alpha1[0]*f[1]); 
  out[6] += 0.3803194146278324*(alpha2[11]*f[31]+f[11]*alpha2[31])+0.3872983346207416*(alpha2[8]*f[25]+f[8]*alpha2[25])+0.4330127018922193*(alpha2[12]*f[20]+f[12]*alpha2[20])+0.3872983346207416*(alpha2[5]*f[19]+f[5]*alpha2[19])+0.4330127018922193*(alpha2[9]*f[16]+f[9]*alpha2[16])+0.3872983346207416*(alpha0[3]*f[13]+alpha2[1]*f[11]+f[1]*alpha2[11])+0.4330127018922193*(alpha2[4]*f[8]+f[4]*alpha2[8]+alpha2[2]*f[5]+f[2]*alpha2[5]+alpha0[0]*f[3]+f[0]*alpha0[3]+alpha2[0]*f[1]+f[0]*alpha2[1]); 
  out[7] += 0.3803194146278324*(alpha2[12]*f[32]+f[12]*alpha2[32])+0.3872983346207416*(alpha2[9]*f[26]+f[9]*alpha2[26]+alpha2[5]*f[20]+f[5]*alpha2[20])+0.4330127018922193*(alpha2[11]*f[19]+f[11]*alpha2[19]+alpha2[8]*f[16]+f[8]*alpha2[16])+0.3872983346207416*(alpha2[2]*f[12]+f[2]*alpha2[12])+0.4330127018922193*(alpha1[4]*f[10]+alpha2[4]*f[9]+f[4]*alpha2[9]+alpha2[1]*f[5]+f[1]*alpha2[5]+alpha1[0]*f[3]+alpha2[0]*f[2]+f[0]*alpha2[2]); 
  out[8] += 0.3803194146278324*(alpha3[11]*f[31]+f[11]*alpha3[31])+0.3872983346207416*(alpha3[6]*f[21]+f[6]*alpha3[21])+0.4330127018922193*(alpha3[12]*f[20]+f[12]*alpha3[20])+0.3872983346207416*(alpha3[5]*f[19]+f[5]*alpha3[19])+0.4330127018922193*(alpha3[7]*f[15]+f[7]*alpha3[15])+0.3872983346207416*(alpha3[1]*f[11]+f[1]*alpha3[11])+0.4330127018922193*(alpha0[3]*f[10]+alpha3[3]*f[6]+f[3]*alpha3[6]+alpha3[2]*f[5]+f[2]*alpha3[5]+alpha0[0]*f[4]+alpha3[0]*f[1]+f[0]*alpha3[1]); 
  out[9] += 0.3803194146278324*(alpha3[12]*f[32]+f[12]*alpha3[32])+0.3872983346207416*(alpha3[7]*f[22]+f[7]*alpha3[22]+alpha3[5]*f[20]+f[5]*alpha3[20])+0.4330127018922193*(alpha3[11]*f[19]+f[11]*alpha3[19]+alpha3[6]*f[15]+f[6]*alpha3[15])+0.3872983346207416*(alpha1[4]*f[14]+alpha3[2]*f[12]+f[2]*alpha3[12])+0.4330127018922193*(alpha3[3]*f[7]+f[3]*alpha3[7]+alpha3[1]*f[5]+f[1]*alpha3[5]+alpha1[0]*f[4]+f[0]*alpha1[4]+alpha3[0]*f[2]+f[0]*alpha3[2]); 
  out[10] += 0.3872983346207416*(alpha2[9]*f[29]+alpha2[8]*f[28])+0.4330127018922193*(alpha2[12]*f[26]+f[12]*alpha2[26]+alpha2[11]*f[25]+f[11]*alpha2[25])+0.3872983346207416*(alpha3[7]*f[24]+alpha3[6]*f[23])+0.4330127018922193*(alpha3[12]*f[22]+f[12]*alpha3[22]+alpha3[11]*f[21]+f[11]*alpha3[21]+alpha2[5]*f[16]+f[5]*alpha2[16]+alpha3[5]*f[15]+f[5]*alpha3[15])+0.3872983346207416*(alpha2[4]*f[14]+alpha3[3]*f[13])+0.4330127018922193*(alpha2[2]*f[9]+f[2]*alpha2[9]+alpha2[1]*f[8]+f[1]*alpha2[8]+alpha3[2]*f[7]+f[2]*alpha3[7]+alpha3[1]*f[6]+f[1]*alpha3[6]+alpha2[0]*f[4]+f[0]*alpha2[4]+alpha3[0]*f[3]+f[0]*alpha3[3]); 
  out[11] += 0.9682458365518543*(alpha0[3]*f[6]+alpha0[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha1[4]*f[9]+alpha1[0]*f[2]); 
  out[13] += 0.9682458365518543*(alpha2[12]*f[22]+alpha2[11]*f[21]+alpha2[9]*f[18]+alpha2[8]*f[17]+alpha2[5]*f[15]+alpha2[4]*f[10]+alpha2[2]*f[7]+alpha2[1]*f[6]+alpha2[0]*f[3]); 
  out[14] += 0.9682458365518543*(alpha3[12]*f[26]+alpha3[11]*f[25]+alpha3[7]*f[18]+alpha3[6]*f[17]+alpha3[5]*f[16]+alpha3[3]*f[10]+alpha3[2]*f[9]+alpha3[1]*f[8]+alpha3[0]*f[4]); 
  out[15] += 0.3803194146278324*(alpha2[20]*f[32]+f[20]*alpha2[32]+alpha2[19]*f[31]+f[19]*alpha2[31])+0.3872983346207416*(alpha2[16]*f[26]+f[16]*alpha2[26]+alpha2[16]*f[25]+f[16]*alpha2[25]+alpha0[3]*f[24])+(0.3464101615137755*alpha2[19]+0.3872983346207416*alpha2[2])*f[20]+0.3464101615137755*f[19]*alpha2[20]+0.3872983346207416*(f[2]*alpha2[20]+alpha2[1]*f[19]+f[1]*alpha2[19])+0.4330127018922193*(alpha1[4]*f[17]+alpha2[4]*f[16]+f[4]*alpha2[16])+0.3872983346207416*(alpha2[5]*f[12]+f[5]*alpha2[12]+alpha2[5]*f[11]+f[5]*alpha2[11])+0.4330127018922193*(alpha2[8]*f[9]+f[8]*alpha2[9]+alpha0[0]*f[7]+alpha1[0]*f[6]+alpha2[0]*f[5]+f[0]*alpha2[5]+f[2]*(alpha0[3]+alpha2[1])+f[1]*alpha2[2]); 
  out[16] += 0.3803194146278324*(alpha3[20]*f[32]+f[20]*alpha3[32]+alpha3[19]*f[31]+f[19]*alpha3[31])+0.3872983346207416*(alpha1[4]*f[28]+alpha3[15]*f[22]+f[15]*alpha3[22]+alpha3[15]*f[21]+f[15]*alpha3[21])+(0.3464101615137755*alpha3[19]+0.3872983346207416*alpha3[2])*f[20]+0.3464101615137755*f[19]*alpha3[20]+0.3872983346207416*(f[2]*alpha3[20]+alpha3[1]*f[19]+f[1]*alpha3[19])+0.4330127018922193*(alpha0[3]*f[18]+alpha3[3]*f[15]+f[3]*alpha3[15])+0.3872983346207416*(alpha3[5]*f[12]+f[5]*alpha3[12]+alpha3[5]*f[11]+f[5]*alpha3[11])+0.4330127018922193*(alpha0[0]*f[9]+alpha1[0]*f[8]+alpha3[6]*f[7]+f[6]*alpha3[7]+alpha3[0]*f[5]+f[0]*alpha3[5]+f[1]*alpha1[4]+alpha3[1]*f[2]+f[1]*alpha3[2]); 
  out[17] += 0.3803194146278324*((alpha2[25]+alpha3[21])*f[31]+f[21]*alpha3[31]+f[25]*alpha2[31])+0.3872983346207416*alpha2[16]*f[29]+0.3464101615137755*alpha2[25]*f[28]+0.3872983346207416*(alpha2[4]*f[28]+alpha0[3]*f[27])+0.4330127018922193*(alpha2[20]*f[26]+f[20]*alpha2[26])+0.3872983346207416*(alpha2[1]*f[25]+f[1]*alpha2[25]+alpha3[15]*f[24])+(0.3464101615137755*alpha3[21]+0.3872983346207416*alpha3[3])*f[23]+0.4330127018922193*(alpha3[20]*f[22]+f[20]*alpha3[22])+0.3872983346207416*(alpha3[1]*f[21]+f[1]*alpha3[21]+(alpha2[16]+alpha3[15])*f[19]+f[15]*alpha3[19]+f[16]*alpha2[19])+0.4330127018922193*(alpha2[2]*f[16]+f[2]*alpha2[16]+alpha3[2]*f[15]+f[2]*alpha3[15])+0.3872983346207416*(alpha2[8]*f[14]+alpha3[6]*f[13]+(alpha2[8]+alpha3[6])*f[11]+f[6]*alpha3[11]+f[8]*alpha2[11])+0.4330127018922193*(alpha0[0]*f[10]+alpha2[5]*f[9]+f[5]*alpha2[9]+alpha2[0]*f[8]+f[0]*alpha2[8]+alpha3[5]*f[7]+f[5]*alpha3[7]+alpha3[0]*f[6]+f[0]*alpha3[6]+(alpha0[3]+alpha2[1])*f[4]+f[1]*alpha2[4]+alpha3[1]*f[3]+f[1]*alpha3[3]); 
  out[18] += 0.3803194146278324*((alpha2[26]+alpha3[22])*f[32]+f[22]*alpha3[32]+f[26]*alpha2[32])+0.3872983346207416*alpha1[4]*f[30]+0.3464101615137755*alpha2[26]*f[29]+0.3872983346207416*(alpha2[4]*f[29]+alpha2[16]*f[28]+alpha2[2]*f[26]+f[2]*alpha2[26])+0.4330127018922193*(alpha2[19]*f[25]+f[19]*alpha2[25])+0.3464101615137755*alpha3[22]*f[24]+0.3872983346207416*(alpha3[3]*f[24]+alpha3[15]*f[23]+alpha3[2]*f[22]+f[2]*alpha3[22])+0.4330127018922193*(alpha3[19]*f[21]+f[19]*alpha3[21])+0.3872983346207416*((alpha2[16]+alpha3[15])*f[20]+f[15]*alpha3[20]+f[16]*alpha2[20])+0.4330127018922193*(alpha2[1]*f[16]+f[1]*alpha2[16]+alpha3[1]*f[15]+f[1]*alpha3[15])+0.3872983346207416*(alpha2[9]*f[14]+alpha3[7]*f[13]+(alpha2[9]+alpha3[7])*f[12]+f[7]*alpha3[12]+f[9]*alpha2[12])+0.4330127018922193*(alpha1[0]*f[10]+alpha2[0]*f[9]+f[0]*alpha2[9]+alpha2[5]*f[8]+f[5]*alpha2[8]+alpha3[0]*f[7]+f[0]*alpha3[7]+alpha3[5]*f[6]+f[5]*alpha3[6]+alpha2[2]*f[4]+f[2]*alpha2[4]+f[3]*(alpha1[4]+alpha3[2])+f[2]*alpha3[3]); 
  out[19] += 0.4330127018922193*alpha1[4]*f[25]+0.9682458365518543*alpha0[3]*f[15]+0.4330127018922193*alpha1[0]*f[11]+0.9682458365518543*alpha0[0]*f[5]; 
  out[20] += 0.4330127018922193*alpha0[3]*f[22]+0.9682458365518543*alpha1[4]*f[16]+0.4330127018922193*alpha0[0]*f[12]+0.9682458365518543*alpha1[0]*f[5]; 
  out[21] += 0.2581988897471612*alpha2[31]*f[31]+0.3803194146278324*(alpha2[1]*f[31]+f[1]*alpha2[31])+0.276641667586244*alpha2[25]*f[25]+0.4330127018922193*(alpha2[4]*f[25]+f[4]*alpha2[25])+0.8660254037844386*alpha0[3]*f[23]+0.3872983346207416*alpha2[20]*f[20]+0.276641667586244*alpha2[19]*f[19]+0.4330127018922193*(alpha2[2]*f[19]+f[2]*alpha2[19])+0.3872983346207416*alpha2[16]*f[16]+0.276641667586244*alpha2[11]*f[11]+0.4330127018922193*(alpha2[0]*f[11]+f[0]*alpha2[11])+0.3872983346207416*alpha2[8]*f[8]+0.9682458365518543*alpha0[0]*f[6]+0.3872983346207416*alpha2[5]*f[5]+f[1]*(0.9682458365518543*alpha0[3]+0.3872983346207416*alpha2[1]); 
  out[22] += 0.2581988897471612*alpha2[32]*f[32]+0.3803194146278324*(alpha2[2]*f[32]+f[2]*alpha2[32])+0.276641667586244*alpha2[26]*f[26]+0.4330127018922193*(alpha2[4]*f[26]+f[4]*alpha2[26])+0.276641667586244*alpha2[20]*f[20]+0.4330127018922193*(alpha2[1]*f[20]+f[1]*alpha2[20])+0.3872983346207416*alpha2[19]*f[19]+0.9682458365518543*alpha1[4]*f[18]+0.3872983346207416*alpha2[16]*f[16]+0.276641667586244*alpha2[12]*f[12]+0.4330127018922193*(alpha2[0]*f[12]+f[0]*alpha2[12])+0.3872983346207416*alpha2[9]*f[9]+0.9682458365518543*alpha1[0]*f[7]+0.3872983346207416*(alpha2[5]*f[5]+alpha2[2]*f[2]); 
  out[23] += 0.3803194146278324*alpha0[3]*f[33]+0.8504200642707612*f[21]*alpha2[31]+0.8660254037844386*f[17]*alpha2[25]+0.9682458365518543*alpha2[20]*f[22]+0.8660254037844386*(alpha2[1]*f[21]+f[15]*alpha2[19])+0.9682458365518543*(alpha2[16]*f[18]+alpha2[4]*f[17]+alpha2[2]*f[15])+0.4330127018922193*alpha0[0]*f[13]+0.8660254037844386*f[6]*alpha2[11]+0.9682458365518543*(alpha2[8]*f[10]+alpha2[5]*f[7]+alpha2[0]*f[6])+(0.3872983346207416*alpha0[3]+0.9682458365518543*alpha2[1])*f[3]; 
  out[24] += 0.8504200642707612*f[22]*alpha2[32]+0.4330127018922193*alpha1[4]*f[27]+0.8660254037844386*(f[18]*alpha2[26]+alpha2[2]*f[22])+0.9682458365518543*alpha2[19]*f[21]+0.8660254037844386*f[15]*alpha2[20]+0.9682458365518543*(alpha2[4]*f[18]+alpha2[16]*f[17]+alpha2[1]*f[15])+0.4330127018922193*alpha1[0]*f[13]+0.8660254037844386*f[7]*alpha2[12]+0.9682458365518543*(alpha2[9]*f[10]+alpha2[0]*f[7]+alpha2[5]*f[6]+alpha2[2]*f[3]); 
  out[25] += 0.2581988897471612*alpha3[31]*f[31]+0.3803194146278324*(alpha3[1]*f[31]+f[1]*alpha3[31])+0.276641667586244*alpha3[21]*f[21]+0.4330127018922193*(alpha3[3]*f[21]+f[3]*alpha3[21])+0.3872983346207416*alpha3[20]*f[20]+0.276641667586244*alpha3[19]*f[19]+0.4330127018922193*(alpha3[2]*f[19]+f[2]*alpha3[19])+0.9682458365518543*alpha0[3]*f[17]+0.3872983346207416*alpha3[15]*f[15]+0.276641667586244*alpha3[11]*f[11]+0.4330127018922193*(alpha3[0]*f[11]+f[0]*alpha3[11])+0.9682458365518543*alpha0[0]*f[8]+0.3872983346207416*(alpha3[6]*f[6]+alpha3[5]*f[5]+alpha3[1]*f[1]); 
  out[26] += 0.2581988897471612*alpha3[32]*f[32]+0.3803194146278324*(alpha3[2]*f[32]+f[2]*alpha3[32])+0.8660254037844386*alpha1[4]*f[29]+0.276641667586244*alpha3[22]*f[22]+0.4330127018922193*(alpha3[3]*f[22]+f[3]*alpha3[22])+0.276641667586244*alpha3[20]*f[20]+0.4330127018922193*(alpha3[1]*f[20]+f[1]*alpha3[20])+0.3872983346207416*(alpha3[19]*f[19]+alpha3[15]*f[15])+0.276641667586244*alpha3[12]*f[12]+0.4330127018922193*(alpha3[0]*f[12]+f[0]*alpha3[12])+0.9682458365518543*alpha1[0]*f[9]+0.3872983346207416*(alpha3[7]*f[7]+alpha3[5]*f[5])+f[2]*(0.9682458365518543*alpha1[4]+0.3872983346207416*alpha3[2]); 
  out[27] += 0.3803194146278324*alpha3[3]*f[33]+0.8660254037844386*alpha2[4]*f[30]+0.9682458365518543*(f[22]*alpha2[26]+f[21]*alpha2[25])+0.4330127018922193*(alpha3[2]*f[24]+alpha3[1]*f[23])+0.3872983346207416*(alpha3[22]*f[22]+alpha3[21]*f[21])+0.9682458365518543*(alpha2[2]*f[18]+alpha2[1]*f[17])+f[15]*(0.9682458365518543*alpha2[16]+0.3872983346207416*alpha3[15])+0.4330127018922193*alpha3[0]*f[13]+0.9682458365518543*(alpha2[0]*f[10]+f[7]*alpha2[9]+f[6]*alpha2[8])+0.3872983346207416*(alpha3[7]*f[7]+alpha3[6]*f[6])+f[3]*(0.9682458365518543*alpha2[4]+0.3872983346207416*alpha3[3]); 
  out[28] += 0.8504200642707612*f[25]*alpha3[31]+0.4330127018922193*alpha0[3]*f[30]+0.9682458365518543*alpha3[20]*f[26]+0.8660254037844386*(alpha3[1]*f[25]+f[17]*alpha3[21]+f[16]*alpha3[19])+0.9682458365518543*(alpha3[15]*f[18]+alpha3[3]*f[17]+alpha3[2]*f[16])+0.4330127018922193*alpha0[0]*f[14]+0.8660254037844386*f[8]*alpha3[11]+0.9682458365518543*(alpha3[6]*f[10]+alpha3[5]*f[9]+alpha3[0]*f[8]+alpha3[1]*f[4]); 
  out[29] += 0.3803194146278324*alpha1[4]*f[34]+f[26]*(0.8504200642707612*alpha3[32]+0.8660254037844386*alpha3[2])+0.9682458365518543*alpha3[19]*f[25]+0.8660254037844386*(f[18]*alpha3[22]+f[16]*alpha3[20])+0.9682458365518543*(alpha3[3]*f[18]+alpha3[15]*f[17]+alpha3[1]*f[16])+0.4330127018922193*alpha1[0]*f[14]+0.8660254037844386*f[9]*alpha3[12]+0.9682458365518543*(alpha3[7]*f[10]+alpha3[0]*f[9]+alpha3[5]*f[8])+(0.3872983346207416*alpha1[4]+0.9682458365518543*alpha3[2])*f[4]; 
  out[30] += 0.3803194146278324*alpha2[4]*f[34]+0.4330127018922193*(alpha2[2]*f[29]+alpha2[1]*f[28])+0.8660254037844386*alpha3[3]*f[27]+(0.3872983346207416*alpha2[26]+0.9682458365518543*alpha3[22])*f[26]+0.3872983346207416*alpha2[25]*f[25]+0.9682458365518543*(alpha3[21]*f[25]+alpha3[2]*f[18]+alpha3[1]*f[17])+(0.3872983346207416*alpha2[16]+0.9682458365518543*alpha3[15])*f[16]+0.4330127018922193*alpha2[0]*f[14]+0.9682458365518543*alpha3[0]*f[10]+(0.3872983346207416*alpha2[9]+0.9682458365518543*alpha3[7])*f[9]+(0.3872983346207416*alpha2[8]+0.9682458365518543*alpha3[6])*f[8]+(0.3872983346207416*alpha2[4]+0.9682458365518543*alpha3[3])*f[4]; 
  out[31] += 1.479019945774904*(alpha0[3]*f[21]+alpha0[0]*f[11])+0.6614378277661477*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[32] += 1.479019945774904*(alpha1[4]*f[26]+alpha1[0]*f[12])+0.6614378277661477*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[33] += 0.6614378277661477*(alpha2[32]*f[32]+alpha2[31]*f[31])+1.479019945774904*alpha2[4]*f[27]+0.6614378277661477*(alpha2[26]*f[26]+alpha2[25]*f[25])+1.479019945774904*(alpha2[2]*f[24]+alpha2[1]*f[23])+0.6614378277661477*(alpha2[20]*f[20]+alpha2[19]*f[19]+alpha2[16]*f[16])+1.479019945774904*alpha2[0]*f[13]+0.6614378277661477*(alpha2[12]*f[12]+alpha2[11]*f[11]+alpha2[9]*f[9]+alpha2[8]*f[8]+alpha2[5]*f[5]+alpha2[4]*f[4]+alpha2[2]*f[2]+alpha2[1]*f[1]+alpha2[0]*f[0]); 
  out[34] += 0.6614378277661477*(alpha3[32]*f[32]+alpha3[31]*f[31])+1.479019945774904*(alpha3[3]*f[30]+alpha3[2]*f[29]+alpha3[1]*f[28])+0.6614378277661477*(alpha3[22]*f[22]+alpha3[21]*f[21]+alpha3[20]*f[20]+alpha3[19]*f[19]+alpha3[15]*f[15])+1.479019945774904*alpha3[0]*f[14]+0.6614378277661477*(alpha3[12]*f[12]+alpha3[11]*f[11]+alpha3[7]*f[7]+alpha3[6]*f[6]+alpha3[5]*f[5]+alpha3[3]*f[3]+alpha3[2]*f[2]+alpha3[1]*f[1]+alpha3[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1+std::abs(amid1)+std::abs(amid2)); 
} 
