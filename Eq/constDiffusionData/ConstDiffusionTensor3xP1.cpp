#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol3xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol3xTensorP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol3xTensorP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol3xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol3xTensorP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol3xTensorP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol3xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol3xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol3xTensorP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol3xTensorP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol3xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol3xTensorP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol3xTensorP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol3xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[4] += 1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[7] += 1.060660171779821*rdxF[0]*nu[1]*f[6]; 

  return (rdxF[0]*(0.3535533905932737*nu[0]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol3xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[2] += 1.060660171779821*f[0]*rdxF[1]*nu[10]; 
  out[4] += 1.060660171779821*f[1]*rdxF[1]*nu[10]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 1.060660171779821*rdxF[1]*f[3]*nu[10]; 
  out[7] += 1.060660171779821*rdxF[1]*f[5]*nu[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 

  return (rdxF[0]*(0.3535533905932737*nu[0])+rdxF[1]*(0.3535533905932737*nu[8]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol3xTensorP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 

  out[1] += 1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[2] += 1.060660171779821*f[0]*rdxF[1]*nu[10]; 
  out[3] += 1.060660171779821*f[0]*rdxF[2]*nu[19]; 
  out[4] += 1.060660171779821*f[1]*rdxF[1]*nu[10]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 1.060660171779821*f[1]*rdxF[2]*nu[19]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 1.060660171779821*f[2]*rdxF[2]*nu[19]+1.060660171779821*rdxF[1]*f[3]*nu[10]; 
  out[7] += 1.060660171779821*rdxF[2]*f[4]*nu[19]+1.060660171779821*rdxF[1]*f[5]*nu[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 

  return (rdxF[0]*(0.3535533905932737*nu[0])+rdxF[1]*(0.3535533905932737*nu[8])+rdxF[2]*(0.3535533905932737*nu[16]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol3xTensorP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[1] += 1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[3] += 1.060660171779821*f[0]*rdxF[1]*nu[19]; 
  out[4] += 1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 1.060660171779821*f[1]*rdxF[1]*nu[19]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 1.060660171779821*rdxF[1]*f[2]*nu[19]; 
  out[7] += 1.060660171779821*rdxF[1]*f[4]*nu[19]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 

  return (rdxF[0]*(0.3535533905932737*nu[0])+rdxF[1]*(0.3535533905932737*nu[16]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol3xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 1.060660171779821*f[0]*rdxF[0]*nu[10]; 
  out[4] += 1.060660171779821*rdxF[0]*f[1]*nu[10]; 
  out[6] += 1.060660171779821*rdxF[0]*f[3]*nu[10]; 
  out[7] += 1.060660171779821*rdxF[0]*f[5]*nu[10]; 

  return (rdxF[0]*(0.3535533905932737*nu[8]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol3xTensorP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[2] += 1.060660171779821*f[0]*rdxF[0]*nu[10]; 
  out[3] += 1.060660171779821*f[0]*rdxF[1]*nu[19]; 
  out[4] += 1.060660171779821*rdxF[0]*f[1]*nu[10]; 
  out[5] += 1.060660171779821*f[1]*rdxF[1]*nu[19]; 
  out[6] += 1.060660171779821*rdxF[1]*f[2]*nu[19]+1.060660171779821*rdxF[0]*f[3]*nu[10]; 
  out[7] += 1.060660171779821*rdxF[1]*f[4]*nu[19]+1.060660171779821*rdxF[0]*f[5]*nu[10]; 

  return (rdxF[0]*(0.3535533905932737*nu[8])+rdxF[1]*(0.3535533905932737*nu[16]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol3xTensorP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 

  out[3] += 1.060660171779821*f[0]*rdxF[0]*nu[19]; 
  out[5] += 1.060660171779821*rdxF[0]*f[1]*nu[19]; 
  out[6] += 1.060660171779821*rdxF[0]*f[2]*nu[19]; 
  out[7] += 1.060660171779821*rdxF[0]*f[4]*nu[19]; 

  return (rdxF[0]*(0.3535533905932737*nu[16]))*0.6666666666666666;

} 
void ConstDiffusionCFLfreqMin3xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin3xTensorP1_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0])+kxSq[1]*(0.3535533905932737*nu[8]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin3xTensorP1_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0])+kxSq[1]*(0.3535533905932737*nu[8])+kxSq[2]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin3xTensorP1_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0])+kxSq[1]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin3xTensorP1_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[8]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin3xTensorP1_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[8])+kxSq[1]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin3xTensorP1_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin3xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin3xTensorP1_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0])+kxSq[1]*(0.3535533905932737*nu[8]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin3xTensorP1_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0])+kxSq[1]*(0.3535533905932737*nu[8])+kxSq[2]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin3xTensorP1_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[0])+kxSq[1]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin3xTensorP1_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[8]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin3xTensorP1_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[8])+kxSq[1]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin3xTensorP1_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.3535533905932737*nu[16]))*1.333333333333333);

} 
