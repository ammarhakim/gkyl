#include <VlasovModDecl.h> 
double VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[3]; 

  alpha0[0] = 4.0*w0dx0; 
  alpha0[2] = 1.154700538379252*dv0dx0; 

  out[1] += 0.8660254037844386*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x1vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[6]; 

  alpha0[0] = 4.0*w0dx0; 
  alpha0[2] = 1.154700538379252*dv0dx0; 

  out[1] += 0.8660254037844386*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[3] += 0.7745966692414833*alpha0[2]*f[5]+0.8660254037844386*(alpha0[0]*f[2]+f[0]*alpha0[2]); 
  out[4] += 1.936491673103709*(alpha0[2]*f[3]+alpha0[0]*f[1]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x1vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[10]; 

  alpha0[0] = 4.0*w0dx0; 
  alpha0[2] = 1.154700538379252*dv0dx0; 

  out[1] += 0.8660254037844386*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[3] += 0.7745966692414833*alpha0[2]*f[5]+0.8660254037844386*(alpha0[0]*f[2]+f[0]*alpha0[2]); 
  out[4] += 1.936491673103709*(alpha0[2]*f[3]+alpha0[0]*f[1]); 
  out[6] += 1.732050807568877*alpha0[2]*f[7]+1.936491673103709*(alpha0[0]*f[3]+f[1]*alpha0[2]); 
  out[7] += 0.7606388292556648*alpha0[2]*f[9]+0.8660254037844386*alpha0[0]*f[5]+0.7745966692414833*alpha0[2]*f[2]; 
  out[8] += 2.958039891549809*(alpha0[2]*f[6]+alpha0[0]*f[4])+1.322875655532295*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
