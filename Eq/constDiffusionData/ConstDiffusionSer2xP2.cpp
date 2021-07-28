#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol2xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[4] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[6] += 6.708203932499369*rdxFnu[0]*f[2]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol2xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[4] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[5] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[6] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[7] += 6.708203932499369*f[1]*rdxFnu[1]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol2xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[5] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[7] += 6.708203932499369*rdxFnu[0]*f[1]; 

  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol2xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol2xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol2xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol2xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol2xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol2xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVarCoeffVol2xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[16]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 3.354101966249685*rdxF[0]*f[1]*nu[4]+1.5*f[0]*rdxF[0]*nu[1]; 
  out[3] += 3.354101966249685*rdxF[0]*f[3]*nu[4]+1.5*rdxF[0]*nu[1]*f[2]; 
  out[4] += 10.06230589874905*rdxF[0]*f[4]*nu[4]+7.5*f[0]*rdxF[0]*nu[4]+6.708203932499369*rdxF[0]*f[1]*nu[1]+3.354101966249685*f[0]*nu[0]*rdxF[0]; 
  out[6] += 10.06230589874905*rdxF[0]*nu[4]*f[6]+7.500000000000001*rdxF[0]*f[2]*nu[4]+6.708203932499369*rdxF[0]*nu[1]*f[3]+3.354101966249684*nu[0]*rdxF[0]*f[2]; 
  out[7] += 3.354101966249685*rdxF[0]*nu[4]*f[7]+1.5*rdxF[0]*nu[1]*f[5]; 

  return (rdxF[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*0.9;

} 
double ConstDiffusionVarCoeffVol2xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[16]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 3.354101966249685*rdxF[0]*f[1]*nu[4]+1.5*f[0]*rdxF[0]*nu[1]; 
  out[2] += 3.354101966249685*rdxF[1]*f[2]*nu[13]+1.5*f[0]*rdxF[1]*nu[10]; 
  out[3] += 3.354101966249685*rdxF[1]*f[3]*nu[13]+1.5*f[1]*rdxF[1]*nu[10]+3.354101966249685*rdxF[0]*f[3]*nu[4]+1.5*rdxF[0]*nu[1]*f[2]; 
  out[4] += 10.06230589874905*rdxF[0]*f[4]*nu[4]+7.5*f[0]*rdxF[0]*nu[4]+6.708203932499369*rdxF[0]*f[1]*nu[1]+3.354101966249685*f[0]*nu[0]*rdxF[0]; 
  out[5] += 10.06230589874905*rdxF[1]*f[5]*nu[13]+7.5*f[0]*rdxF[1]*nu[13]+6.708203932499369*rdxF[1]*f[2]*nu[10]+3.354101966249685*f[0]*rdxF[1]*nu[8]; 
  out[6] += 3.354101966249685*rdxF[1]*f[6]*nu[13]+1.5*rdxF[1]*f[4]*nu[10]+10.06230589874905*rdxF[0]*nu[4]*f[6]+7.500000000000001*rdxF[0]*f[2]*nu[4]+6.708203932499369*rdxF[0]*nu[1]*f[3]+3.354101966249684*nu[0]*rdxF[0]*f[2]; 
  out[7] += 10.06230589874905*rdxF[1]*f[7]*nu[13]+7.500000000000001*f[1]*rdxF[1]*nu[13]+6.708203932499369*rdxF[1]*f[3]*nu[10]+3.354101966249684*f[1]*rdxF[1]*nu[8]+3.354101966249685*rdxF[0]*nu[4]*f[7]+1.5*rdxF[0]*nu[1]*f[5]; 

  return (rdxF[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+rdxF[1]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*0.9;

} 
double ConstDiffusionVarCoeffVol2xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[16]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 3.354101966249685*rdxF[0]*f[2]*nu[13]+1.5*f[0]*rdxF[0]*nu[10]; 
  out[3] += 3.354101966249685*rdxF[0]*f[3]*nu[13]+1.5*rdxF[0]*f[1]*nu[10]; 
  out[5] += 10.06230589874905*rdxF[0]*f[5]*nu[13]+7.5*f[0]*rdxF[0]*nu[13]+6.708203932499369*rdxF[0]*f[2]*nu[10]+3.354101966249685*f[0]*rdxF[0]*nu[8]; 
  out[6] += 3.354101966249685*rdxF[0]*f[6]*nu[13]+1.5*rdxF[0]*f[4]*nu[10]; 
  out[7] += 10.06230589874905*rdxF[0]*f[7]*nu[13]+7.500000000000001*rdxF[0]*f[1]*nu[13]+6.708203932499369*rdxF[0]*f[3]*nu[10]+3.354101966249684*rdxF[0]*f[1]*nu[8]; 

  return (rdxF[0]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*0.9;

} 
void ConstDiffusionCFLfreqMin2xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*1.8);

} 
void ConstDiffusionCFLfreqMin2xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+kxSq[1]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*1.8);

} 
void ConstDiffusionCFLfreqMin2xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin2xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin2xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+kxSq[1]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin2xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin2xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin2xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+kxSq[1]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin2xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[16]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[13])-0.5590169943749475*nu[12]+0.5*nu[8]))*1.8);

} 
