#include <VlasovModDecl.h> 
double VlasovVolStream1x2vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[4]; 

  alpha0[0] = 5.656854249492382*w0dx0; 
  alpha0[2] = 1.632993161855453*dv0dx0; 

  out[1] += 0.6123724356957944*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x2vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[10]; 

  alpha0[0] = 5.656854249492382*w0dx0; 
  alpha0[2] = 1.632993161855453*dv0dx0; 

  out[1] += 0.6123724356957944*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[4] += 0.5477225575051661*alpha0[2]*f[8]+0.6123724356957944*(alpha0[0]*f[2]+f[0]*alpha0[2]); 
  out[5] += 0.6123724356957944*(alpha0[2]*f[6]+alpha0[0]*f[3]); 
  out[7] += 1.369306393762915*(alpha0[2]*f[4]+alpha0[0]*f[1]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
double VlasovVolStream1x2vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  double alpha0[20]; 

  alpha0[0] = 5.656854249492382*w0dx0; 
  alpha0[2] = 1.632993161855453*dv0dx0; 

  out[1] += 0.6123724356957944*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
  out[4] += 0.5477225575051661*alpha0[2]*f[8]+0.6123724356957944*(alpha0[0]*f[2]+f[0]*alpha0[2]); 
  out[5] += 0.6123724356957944*(alpha0[2]*f[6]+alpha0[0]*f[3]); 
  out[7] += 1.369306393762915*(alpha0[2]*f[4]+alpha0[0]*f[1]); 
  out[10] += 0.5477225575051661*alpha0[2]*f[14]+0.6123724356957944*(alpha0[0]*f[6]+alpha0[2]*f[3]); 
  out[11] += 1.224744871391589*alpha0[2]*f[12]+1.369306393762915*(alpha0[0]*f[4]+f[1]*alpha0[2]); 
  out[12] += 0.537852874200477*alpha0[2]*f[18]+0.6123724356957944*alpha0[0]*f[8]+0.5477225575051661*alpha0[2]*f[2]; 
  out[13] += 1.369306393762915*(alpha0[2]*f[10]+alpha0[0]*f[5]); 
  out[15] += 0.6123724356957944*(alpha0[2]*f[16]+alpha0[0]*f[9]); 
  out[17] += 2.091650066335188*(alpha0[2]*f[11]+alpha0[0]*f[7])+0.9354143466934851*(alpha0[2]*f[2]+alpha0[0]*f[0]); 
return std::abs(w0dx0)+dv0dx0/2; 
} 
