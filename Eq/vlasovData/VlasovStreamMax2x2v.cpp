#include <VlasovModDecl.h> 
double VlasovVolStream2x2vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
return std::abs(w0dx0)+std::abs(w1dx1)+dv0dx0/2+dv1dx1/2; 
} 
double VlasovVolStream2x2vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[5] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[8]*dv1dx1+f[7]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+(0.8944271909999159*f[13]+f[0])*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w1dx1+f[10]*dv1dx1; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[10]*dv0dx0; 
  out[9] += 3.464101615137754*f[4]*w1dx1+(0.8944271909999159*f[14]+f[0])*dv1dx1; 
  out[11] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[6]*dv0dx0; 
  out[12] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[9]*dv1dx1; 
return std::abs(w0dx0)+std::abs(w1dx1)+dv0dx0/2+dv1dx1/2; 
} 
double VlasovVolStream2x2vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[5] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[8]*dv1dx1+f[7]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+(0.8944271909999159*f[13]+f[0])*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w1dx1+f[10]*dv1dx1; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[10]*dv0dx0; 
  out[9] += 3.464101615137754*f[4]*w1dx1+(0.8944271909999159*f[14]+f[0])*dv1dx1; 
  out[11] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[6]*dv0dx0; 
  out[12] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[9]*dv1dx1; 
  out[15] += 3.464101615137754*f[6]*w1dx1+3.464101615137754*f[7]*w0dx0+f[17]*dv1dx1+(0.8944271909999161*f[24]+f[2])*dv0dx0; 
  out[16] += 3.464101615137754*f[8]*w1dx1+3.464101615137754*f[9]*w0dx0+(0.8944271909999161*f[28]+f[1])*dv1dx1+f[18]*dv0dx0; 
  out[17] += 3.464101615137754*f[10]*w0dx0+(0.8944271909999161*f[27]+f[4])*dv0dx0; 
  out[18] += 3.464101615137754*f[10]*w1dx1+(0.8944271909999161*f[30]+f[3])*dv1dx1; 
  out[19] += 3.464101615137755*f[11]*w1dx1+7.745966692414834*f[5]*w0dx0+f[25]*dv1dx1+2.23606797749979*f[15]*dv0dx0; 
  out[20] += 7.745966692414834*f[5]*w1dx1+3.464101615137755*f[12]*w0dx0+2.23606797749979*f[16]*dv1dx1+f[22]*dv0dx0; 
  out[21] += 7.745966692414834*f[6]*w0dx0+(2.0*f[23]+2.23606797749979*f[1])*dv0dx0; 
  out[22] += 7.745966692414834*f[7]*w1dx1+2.23606797749979*f[18]*dv1dx1; 
  out[23] += 3.464101615137755*f[13]*w0dx0+(0.8783100656536798*f[33]+0.8944271909999161*f[3])*dv0dx0; 
  out[24] += 3.464101615137755*f[13]*w1dx1+f[27]*dv1dx1; 
  out[25] += 7.745966692414834*f[8]*w0dx0+2.23606797749979*f[17]*dv0dx0; 
  out[26] += 7.745966692414834*f[9]*w1dx1+(2.0*f[29]+2.23606797749979*f[2])*dv1dx1; 
  out[28] += 3.464101615137755*f[14]*w0dx0+f[30]*dv0dx0; 
  out[29] += 3.464101615137755*f[14]*w1dx1+(0.8783100656536798*f[34]+0.8944271909999161*f[4])*dv1dx1; 
  out[31] += (11.83215956619923*f[11]+5.291502622129181*f[0])*w0dx0+(3.415650255319866*f[21]+1.527525231651947*f[3])*dv0dx0; 
  out[32] += (11.83215956619923*f[12]+5.291502622129181*f[0])*w1dx1+(3.415650255319866*f[26]+1.527525231651947*f[4])*dv1dx1; 
return std::abs(w0dx0)+std::abs(w1dx1)+dv0dx0/2+dv1dx1/2; 
} 
