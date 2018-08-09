#include <VlasovModDecl.h> 
double VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x1vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[3] += 3.464101615137754*f[2]*w0dx0+(0.8944271909999159*f[5]+f[0])*dv0dx0; 
  out[4] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[3]*dv0dx0; 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x1vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
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
return std::abs(w0dx0)+dv0dx0/2; 
} 
