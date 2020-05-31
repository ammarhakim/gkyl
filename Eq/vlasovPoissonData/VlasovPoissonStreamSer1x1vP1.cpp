#include <VlasovPoissonModDecl.h> 
__host__ __device__ double VlasovPoissonVolStream1x1vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f: Input distribution function.
  // out: Incremented output.
  double w1Ddx0  = w[1]/dxv[0]; 
  double dv1Ddx0 = dxv[1]/dxv[0]; 

  out[1] += 3.464101615137754*f[0]*w1Ddx0+f[2]*dv1Ddx0; 
  out[3] += 3.464101615137754*f[2]*w1Ddx0+f[0]*dv1Ddx0; 

  return std::abs(w0Ddx0)+0.5*dv0Ddx0;
} 
