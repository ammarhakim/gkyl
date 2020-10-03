#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xSerP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 


  return (rdxSq4nu[0])*0.6666666666666666;

} 
