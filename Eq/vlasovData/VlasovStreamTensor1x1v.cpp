#include <VlasovModDecl.h> 
double VlasovVolStream1x1vTensorP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[3] += 3.464101615137754*f[2]*w0dx0+f[0]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x1vTensorP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[3] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[5]+f[0])*dv0dx0; 
  out[4] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[3]*dv0dx0; 
  out[6] += 7.745966692414834*f[3]*w0dx0+(2.0*f[7]+2.23606797749979*f[1])*dv0dx0; 
  out[7] += 3.464101615137755*f[5]*w0dx0+0.8944271909999161*f[2]*dv0dx0; 
  out[8] += 7.745966692414834*f[7]*w0dx0+2.0*f[3]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x1vTensorP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[3] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[5]+f[0])*dv0dx0; 
  out[4] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[3]*dv0dx0; 
  out[6] += 7.745966692414834*f[3]*w0dx0+(2.0*f[7]+2.23606797749979*f[1])*dv0dx0; 
  out[7] += 3.464101615137755*f[5]*w0dx0+(0.8783100656536798*f[9]+0.8944271909999161*f[2])*dv0dx0; 
  out[8] += (11.83215956619923*f[4]+5.291502622129181*f[0])*w0dx0+(3.415650255319866*f[6]+1.527525231651947*f[2])*dv0dx0; 
  out[10] += 7.745966692414834*f[7]*w0dx0+(1.963961012123931*f[12]+2.0*f[3])*dv0dx0; 
  out[11] += (11.83215956619923*f[6]+5.291502622129181*f[2])*w0dx0+(3.055050463303893*f[10]+1.366260102127946*f[5]+3.415650255319866*f[4]+1.527525231651947*f[0])*dv0dx0; 
  out[12] += 3.464101615137754*f[9]*w0dx0+0.8783100656536798*f[5]*dv0dx0; 
  out[13] += (11.83215956619923*f[10]+5.291502622129181*f[5])*w0dx0+(3.0*f[14]+1.341640786499874*f[9]+3.055050463303893*f[6]+1.366260102127947*f[2])*dv0dx0; 
  out[14] += 7.745966692414833*f[12]*w0dx0+1.963961012123931*f[7]*dv0dx0; 
  out[15] += (11.83215956619923*f[14]+5.291502622129181*f[9])*w0dx0+(3.0*f[10]+1.341640786499874*f[5])*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
