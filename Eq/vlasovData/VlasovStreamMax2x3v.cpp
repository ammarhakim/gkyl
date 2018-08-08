#include <VlasovModDecl.h> 
double VlasovVolStream2x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
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
double VlasovVolStream2x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[6] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[9]*dv1dx1+f[8]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w0dx0+(0.8944271909999159*f[18]+f[0])*dv0dx0; 
  out[8] += 3.464101615137754*f[3]*w1dx1+f[11]*dv1dx1; 
  out[9] += 3.464101615137754*f[4]*w0dx0+f[11]*dv0dx0; 
  out[10] += 3.464101615137754*f[4]*w1dx1+(0.8944271909999159*f[19]+f[0])*dv1dx1; 
  out[12] += 3.464101615137754*f[5]*w0dx0+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[5]*w1dx1+f[15]*dv1dx1; 
  out[16] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[7]*dv0dx0; 
  out[17] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[10]*dv1dx1; 
return std::abs(w0dx0)+std::abs(w1dx1)+dv0dx0/2+dv1dx1/2; 
} 
double VlasovVolStream2x3vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[6] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[9]*dv1dx1+f[8]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w0dx0+(0.8944271909999159*f[18]+f[0])*dv0dx0; 
  out[8] += 3.464101615137754*f[3]*w1dx1+f[11]*dv1dx1; 
  out[9] += 3.464101615137754*f[4]*w0dx0+f[11]*dv0dx0; 
  out[10] += 3.464101615137754*f[4]*w1dx1+(0.8944271909999159*f[19]+f[0])*dv1dx1; 
  out[12] += 3.464101615137754*f[5]*w0dx0+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[5]*w1dx1+f[15]*dv1dx1; 
  out[16] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[7]*dv0dx0; 
  out[17] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[10]*dv1dx1; 
  out[21] += 3.464101615137754*f[7]*w1dx1+3.464101615137754*f[8]*w0dx0+f[23]*dv1dx1+(0.8944271909999161*f[36]+f[2])*dv0dx0; 
  out[22] += 3.464101615137754*f[9]*w1dx1+3.464101615137754*f[10]*w0dx0+(0.8944271909999161*f[40]+f[1])*dv1dx1+f[24]*dv0dx0; 
  out[23] += 3.464101615137754*f[11]*w0dx0+(0.8944271909999161*f[39]+f[4])*dv0dx0; 
  out[24] += 3.464101615137754*f[11]*w1dx1+(0.8944271909999161*f[42]+f[3])*dv1dx1; 
  out[25] += 3.464101615137754*f[12]*w1dx1+3.464101615137754*f[13]*w0dx0+f[28]*dv1dx1+f[27]*dv0dx0; 
  out[26] += 3.464101615137754*f[14]*w0dx0+(0.8944271909999161*f[45]+f[5])*dv0dx0; 
  out[27] += 3.464101615137754*f[14]*w1dx1+f[30]*dv1dx1; 
  out[28] += 3.464101615137754*f[15]*w0dx0+f[30]*dv0dx0; 
  out[29] += 3.464101615137754*f[15]*w1dx1+(0.8944271909999161*f[46]+f[5])*dv1dx1; 
  out[31] += 3.464101615137755*f[16]*w1dx1+7.745966692414834*f[6]*w0dx0+f[37]*dv1dx1+2.23606797749979*f[21]*dv0dx0; 
  out[32] += 7.745966692414834*f[6]*w1dx1+3.464101615137755*f[17]*w0dx0+2.23606797749979*f[22]*dv1dx1+f[34]*dv0dx0; 
  out[33] += 7.745966692414834*f[7]*w0dx0+(2.0*f[35]+2.23606797749979*f[1])*dv0dx0; 
  out[34] += 7.745966692414834*f[8]*w1dx1+2.23606797749979*f[24]*dv1dx1; 
  out[35] += 3.464101615137755*f[18]*w0dx0+(0.8783100656536798*f[53]+0.8944271909999161*f[3])*dv0dx0; 
  out[36] += 3.464101615137755*f[18]*w1dx1+f[39]*dv1dx1; 
  out[37] += 7.745966692414834*f[9]*w0dx0+2.23606797749979*f[23]*dv0dx0; 
  out[38] += 7.745966692414834*f[10]*w1dx1+(2.0*f[41]+2.23606797749979*f[2])*dv1dx1; 
  out[40] += 3.464101615137755*f[19]*w0dx0+f[42]*dv0dx0; 
  out[41] += 3.464101615137755*f[19]*w1dx1+(0.8783100656536798*f[54]+0.8944271909999161*f[4])*dv1dx1; 
  out[43] += 7.745966692414834*f[12]*w0dx0+2.23606797749979*f[26]*dv0dx0; 
  out[44] += 7.745966692414834*f[13]*w1dx1+2.23606797749979*f[29]*dv1dx1; 
  out[47] += 3.464101615137755*f[20]*w0dx0+f[49]*dv0dx0; 
  out[48] += 3.464101615137755*f[20]*w1dx1+f[50]*dv1dx1; 
  out[51] += (11.83215956619923*f[16]+5.291502622129181*f[0])*w0dx0+(3.415650255319866*f[33]+1.527525231651947*f[3])*dv0dx0; 
  out[52] += (11.83215956619923*f[17]+5.291502622129181*f[0])*w1dx1+(3.415650255319866*f[38]+1.527525231651947*f[4])*dv1dx1; 
return std::abs(w0dx0)+std::abs(w1dx1)+dv0dx0/2+dv1dx1/2; 
} 
