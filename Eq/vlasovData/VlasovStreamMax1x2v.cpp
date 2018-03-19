#include <VlasovModDecl.h> 
double VlasovVolStream1x2vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
return std:abs(w0dx0); 
} 
double VlasovVolStream1x2vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[8]*dv0dx0+f[0]*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[4]*dv0dx0; 
return std:abs(w0dx0); 
} 
double VlasovVolStream1x2vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[8]*dv0dx0+f[0]*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[4]*dv0dx0; 
  out[10] += 3.464101615137754*f[6]*w0dx0+0.8944271909999161*f[14]*dv0dx0+f[3]*dv0dx0; 
  out[11] += 7.745966692414834*f[4]*w0dx0+2.0*f[12]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[12] += 3.464101615137755*f[8]*w0dx0+0.8783100656536798*f[18]*dv0dx0+0.8944271909999161*f[2]*dv0dx0; 
  out[13] += 7.745966692414834*f[5]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[15] += 3.464101615137755*f[9]*w0dx0+f[16]*dv0dx0; 
  out[17] += 11.83215956619923*f[7]*w0dx0+5.291502622129181*f[0]*w0dx0+3.415650255319866*f[11]*dv0dx0+1.527525231651947*f[2]*dv0dx0; 
return std:abs(w0dx0); 
} 
double VlasovVolStream1x2vMaxP4(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[8]*dv0dx0+f[0]*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[4]*dv0dx0; 
  out[10] += 3.464101615137754*f[6]*w0dx0+0.8944271909999161*f[14]*dv0dx0+f[3]*dv0dx0; 
  out[11] += 7.745966692414834*f[4]*w0dx0+2.0*f[12]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[12] += 3.464101615137755*f[8]*w0dx0+0.8783100656536798*f[18]*dv0dx0+0.8944271909999161*f[2]*dv0dx0; 
  out[13] += 7.745966692414834*f[5]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[15] += 3.464101615137755*f[9]*w0dx0+f[16]*dv0dx0; 
  out[17] += 11.83215956619923*f[7]*w0dx0+5.291502622129181*f[0]*w0dx0+3.415650255319866*f[11]*dv0dx0+1.527525231651947*f[2]*dv0dx0; 
  out[20] += 7.745966692414834*f[10]*w0dx0+2.0*f[21]*dv0dx0+2.23606797749979*f[5]*dv0dx0; 
  out[21] += 3.464101615137755*f[14]*w0dx0+0.8783100656536798*f[29]*dv0dx0+0.8944271909999159*f[6]*dv0dx0; 
  out[22] += 3.464101615137755*f[16]*w0dx0+0.8944271909999159*f[25]*dv0dx0+f[9]*dv0dx0; 
  out[23] += 7.745966692414834*f[12]*w0dx0+1.963961012123931*f[27]*dv0dx0+2.0*f[4]*dv0dx0; 
  out[24] += 7.745966692414834*f[15]*w0dx0+2.23606797749979*f[22]*dv0dx0; 
  out[26] += 11.83215956619923*f[11]*w0dx0+5.291502622129181*f[2]*w0dx0+3.055050463303893*f[23]*dv0dx0+1.366260102127946*f[8]*dv0dx0+3.415650255319866*f[7]*dv0dx0+1.527525231651947*f[0]*dv0dx0; 
  out[27] += 3.464101615137754*f[18]*w0dx0+0.8728715609439696*f[33]*dv0dx0+0.8783100656536798*f[8]*dv0dx0; 
  out[28] += 11.83215956619923*f[13]*w0dx0+5.291502622129181*f[3]*w0dx0+3.415650255319866*f[20]*dv0dx0+1.527525231651947*f[6]*dv0dx0; 
  out[30] += 3.464101615137754*f[19]*w0dx0+f[31]*dv0dx0; 
  out[32] += 15.87450786638754*f[17]*w0dx0+10.39230484541326*f[1]*w0dx0+4.58257569495584*f[26]*dv0dx0+3.0*f[4]*dv0dx0; 
return std:abs(w0dx0); 
} 
