#include <VlasovModDecl.h> 
__host__ __device__ double VlasovVolForce1x1vMaxP1(const double *w, const double *dxv, const double *E, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. E/f: Input force/distribution function. out: Incremented output 
  double dv10 = 2/dxv[1]; 
  const double *E0 = &E[0]; 
  out[2] += (1.224744871391589*E0[1]*f[1]+1.224744871391589*E0[0]*f[0])*dv10; 

  return std::abs(E0[0])*dv10*0.5; 
} 
