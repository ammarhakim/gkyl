#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xMaxP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
