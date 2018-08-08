#include <VlasovModDecl.h> 
double VlasovVolStream1x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[5] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[12]+f[0])*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+f[7]*dv0dx0; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[9]*dv0dx0; 
  out[11] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[5]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x3vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[5] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[12]+f[0])*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+f[7]*dv0dx0; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[9]*dv0dx0; 
  out[11] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[5]*dv0dx0; 
  out[15] += 3.464101615137754*f[7]*w0dx0+(0.8944271909999161*f[22]+f[3])*dv0dx0; 
  out[16] += 3.464101615137754*f[9]*w0dx0+(0.8944271909999161*f[26]+f[4])*dv0dx0; 
  out[17] += 3.464101615137754*f[10]*w0dx0+f[18]*dv0dx0; 
  out[19] += 7.745966692414834*f[5]*w0dx0+(2.0*f[20]+2.23606797749979*f[1])*dv0dx0; 
  out[20] += 3.464101615137755*f[12]*w0dx0+(0.8783100656536798*f[32]+0.8944271909999161*f[2])*dv0dx0; 
  out[21] += 7.745966692414834*f[6]*w0dx0+2.23606797749979*f[15]*dv0dx0; 
  out[23] += 3.464101615137755*f[13]*w0dx0+f[24]*dv0dx0; 
  out[25] += 7.745966692414834*f[8]*w0dx0+2.23606797749979*f[16]*dv0dx0; 
  out[28] += 3.464101615137755*f[14]*w0dx0+f[29]*dv0dx0; 
  out[31] += (11.83215956619923*f[11]+5.291502622129181*f[0])*w0dx0+(3.415650255319866*f[19]+1.527525231651947*f[2])*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
