#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
double ConstHyperDiffusion4Vol1xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol1xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[2]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 2.121320343559642*f[0]*rdxF[0]*nu[1]; 

  return (rdxF[0]*(0.7071067811865475*nu[0]))*0.6666666666666666;

} 
void ConstDiffusionCFLfreqMin1xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[2]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.7071067811865475*nu[0]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin1xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[2]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.7071067811865475*nu[0]))*1.333333333333333);

} 
