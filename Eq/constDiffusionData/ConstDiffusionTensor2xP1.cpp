#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol2xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol2xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol2xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol2xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol2xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol2xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol2xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[8]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 1.5*f[0]*rdxF[0]*nu[1]; 
  out[3] += 1.5*rdxF[0]*nu[1]*f[2]; 

  return (rdxF[0]*(0.5*nu[0]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol2xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[8]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 1.5*f[0]*rdxF[0]*nu[1]; 
  out[2] += 1.5*f[0]*rdxF[1]*nu[6]; 
  out[3] += 1.5*f[1]*rdxF[1]*nu[6]+1.5*rdxF[0]*nu[1]*f[2]; 

  return (rdxF[0]*(0.5*nu[0])+rdxF[1]*(0.5*nu[4]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol2xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[8]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 1.5*f[0]*rdxF[0]*nu[6]; 
  out[3] += 1.5*rdxF[0]*f[1]*nu[6]; 

  return (rdxF[0]*(0.5*nu[4]))*0.6666666666666666;

} 
void ConstDiffusionCFLfreqMin2xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[8]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.5*nu[0]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin2xTensorP1_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[8]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.5*nu[0])+kxSq[1]*(0.5*nu[4]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin2xTensorP1_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[8]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.5*nu[4]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin2xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[8]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.5*nu[0]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin2xTensorP1_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[8]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.5*nu[0])+kxSq[1]*(0.5*nu[4]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin2xTensorP1_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[8]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.5*nu[4]))*1.333333333333333);

} 
