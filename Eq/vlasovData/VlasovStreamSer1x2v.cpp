#include <VlasovModDecl.h> 
double VlasovVolStream1x2vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+f[0]*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 3.464101615137754*f[6]*w0dx0+f[3]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x2vSerP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[8]+f[0])*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[4]*dv0dx0; 
  out[10] += 3.464101615137754*f[6]*w0dx0+(0.8944271909999161*f[14]+f[3])*dv0dx0; 
  out[11] += 7.745966692414834*f[4]*w0dx0+(2.0*f[12]+2.23606797749979*f[1])*dv0dx0; 
  out[12] += 3.464101615137755*f[8]*w0dx0+0.8944271909999161*f[2]*dv0dx0; 
  out[13] += 7.745966692414834*f[5]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[15] += 3.464101615137755*f[9]*w0dx0+f[16]*dv0dx0; 
  out[17] += 7.745966692414834*f[10]*w0dx0+(2.0*f[18]+2.23606797749979*f[5])*dv0dx0; 
  out[18] += 3.464101615137755*f[14]*w0dx0+0.8944271909999159*f[6]*dv0dx0; 
  out[19] += 3.464101615137755*f[16]*w0dx0+f[9]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x2vSerP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[8]+f[0])*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[4]*dv0dx0; 
  out[10] += 3.464101615137754*f[6]*w0dx0+(0.8944271909999161*f[14]+f[3])*dv0dx0; 
  out[11] += 7.745966692414834*f[4]*w0dx0+(2.0*f[12]+2.23606797749979*f[1])*dv0dx0; 
  out[12] += 3.464101615137755*f[8]*w0dx0+(0.8783100656536798*f[18]+0.8944271909999161*f[2])*dv0dx0; 
  out[13] += 7.745966692414834*f[5]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[15] += 3.464101615137755*f[9]*w0dx0+f[16]*dv0dx0; 
  out[17] += (11.83215956619923*f[7]+5.291502622129181*f[0])*w0dx0+(3.415650255319866*f[11]+1.527525231651947*f[2])*dv0dx0; 
  out[20] += 7.745966692414834*f[10]*w0dx0+(2.0*f[21]+2.23606797749979*f[5])*dv0dx0; 
  out[21] += 3.464101615137755*f[14]*w0dx0+(0.8783100656536798*f[26]+0.8944271909999159*f[6])*dv0dx0; 
  out[22] += 3.464101615137755*f[16]*w0dx0+f[9]*dv0dx0; 
  out[23] += (11.83215956619923*f[11]+5.291502622129181*f[2])*w0dx0+(1.366260102127946*f[8]+3.415650255319866*f[7]+1.527525231651947*f[0])*dv0dx0; 
  out[24] += 3.464101615137754*f[18]*w0dx0+0.8783100656536798*f[8]*dv0dx0; 
  out[25] += (11.83215956619923*f[13]+5.291502622129181*f[3])*w0dx0+(3.415650255319866*f[20]+1.527525231651947*f[6])*dv0dx0; 
  out[27] += 3.464101615137754*f[19]*w0dx0+f[28]*dv0dx0; 
  out[29] += (11.83215956619923*f[20]+5.291502622129181*f[6])*w0dx0+(1.366260102127946*f[14]+3.415650255319866*f[13]+1.527525231651947*f[3])*dv0dx0; 
  out[30] += 3.464101615137754*f[26]*w0dx0+0.8783100656536798*f[14]*dv0dx0; 
  out[31] += 3.464101615137754*f[28]*w0dx0+f[19]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
