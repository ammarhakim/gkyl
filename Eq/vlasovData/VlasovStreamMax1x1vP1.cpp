#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVolStream1x1vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 

  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 

  return std::abs(w0dx0)+dv0dx0/2; 
} 
