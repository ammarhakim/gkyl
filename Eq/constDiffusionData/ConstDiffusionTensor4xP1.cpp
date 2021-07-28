#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol4xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[3] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[2] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstDiffusionVol4xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[3] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.6666666666666666;

} 
double ConstHyperDiffusion4Vol4xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 0.75*f[0]*rdxF[0]*nu[18]; 
  out[5] += 0.75*rdxF[0]*f[1]*nu[18]; 
  out[7] += 0.75*rdxF[0]*f[3]*nu[18]; 
  out[9] += 0.75*rdxF[0]*f[4]*nu[18]; 
  out[11] += 0.75*rdxF[0]*f[6]*nu[18]; 
  out[12] += 0.75*rdxF[0]*f[8]*nu[18]; 
  out[14] += 0.75*rdxF[0]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[0]*f[13]*nu[18]; 

  return (rdxF[0]*(0.25*nu[16]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[2] += 0.75*f[0]*rdxF[0]*nu[18]; 
  out[3] += 0.75*f[0]*rdxF[1]*nu[35]; 
  out[5] += 0.75*rdxF[0]*f[1]*nu[18]; 
  out[6] += 0.75*f[1]*rdxF[1]*nu[35]; 
  out[7] += 0.75*rdxF[1]*f[2]*nu[35]+0.75*rdxF[0]*f[3]*nu[18]; 
  out[9] += 0.75*rdxF[0]*f[4]*nu[18]; 
  out[10] += 0.75*rdxF[1]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[1]*f[5]*nu[35]+0.75*rdxF[0]*f[6]*nu[18]; 
  out[12] += 0.75*rdxF[0]*f[8]*nu[18]; 
  out[13] += 0.75*rdxF[1]*f[8]*nu[35]; 
  out[14] += 0.75*rdxF[1]*f[9]*nu[35]+0.75*rdxF[0]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[1]*f[12]*nu[35]+0.75*rdxF[0]*f[13]*nu[18]; 

  return (rdxF[0]*(0.25*nu[16])+rdxF[1]*(0.25*nu[32]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[2] += 0.75*f[0]*rdxF[0]*nu[18]; 
  out[3] += 0.75*f[0]*rdxF[1]*nu[35]; 
  out[4] += 0.75*f[0]*rdxF[2]*nu[52]; 
  out[5] += 0.75*rdxF[0]*f[1]*nu[18]; 
  out[6] += 0.75*f[1]*rdxF[1]*nu[35]; 
  out[7] += 0.75*rdxF[1]*f[2]*nu[35]+0.75*rdxF[0]*f[3]*nu[18]; 
  out[8] += 0.75*f[1]*rdxF[2]*nu[52]; 
  out[9] += 0.75*f[2]*rdxF[2]*nu[52]+0.75*rdxF[0]*f[4]*nu[18]; 
  out[10] += 0.75*rdxF[2]*f[3]*nu[52]+0.75*rdxF[1]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[1]*f[5]*nu[35]+0.75*rdxF[0]*f[6]*nu[18]; 
  out[12] += 0.75*rdxF[2]*f[5]*nu[52]+0.75*rdxF[0]*f[8]*nu[18]; 
  out[13] += 0.75*rdxF[2]*f[6]*nu[52]+0.75*rdxF[1]*f[8]*nu[35]; 
  out[14] += 0.75*rdxF[2]*f[7]*nu[52]+0.75*rdxF[1]*f[9]*nu[35]+0.75*rdxF[0]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[2]*f[11]*nu[52]+0.75*rdxF[1]*f[12]*nu[35]+0.75*rdxF[0]*f[13]*nu[18]; 

  return (rdxF[0]*(0.25*nu[16])+rdxF[1]*(0.25*nu[32])+rdxF[2]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 
  rdxF[3] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 0.75*f[0]*rdxF[1]*nu[18]; 
  out[3] += 0.75*f[0]*rdxF[2]*nu[35]; 
  out[4] += 0.75*f[0]*rdxF[3]*nu[52]; 
  out[5] += 0.75*f[1]*rdxF[1]*nu[18]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*f[1]*rdxF[2]*nu[35]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 0.75*f[2]*rdxF[2]*nu[35]+0.75*rdxF[1]*f[3]*nu[18]; 
  out[8] += 0.75*f[1]*rdxF[3]*nu[52]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 0.75*f[2]*rdxF[3]*nu[52]+0.75*rdxF[1]*f[4]*nu[18]; 
  out[10] += 0.75*f[3]*rdxF[3]*nu[52]+0.75*rdxF[2]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[2]*f[5]*nu[35]+0.75*rdxF[1]*f[6]*nu[18]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[3]*f[5]*nu[52]+0.75*rdxF[1]*f[8]*nu[18]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[3]*f[6]*nu[52]+0.75*rdxF[2]*f[8]*nu[35]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[14] += 0.75*rdxF[3]*f[7]*nu[52]+0.75*rdxF[2]*f[9]*nu[35]+0.75*rdxF[1]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[3]*f[11]*nu[52]+0.75*rdxF[2]*f[12]*nu[35]+0.75*rdxF[1]*f[13]*nu[18]+0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0])+rdxF[1]*(0.25*nu[16])+rdxF[2]*(0.25*nu[32])+rdxF[3]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 0.75*f[0]*rdxF[1]*nu[18]; 
  out[3] += 0.75*f[0]*rdxF[2]*nu[35]; 
  out[5] += 0.75*f[1]*rdxF[1]*nu[18]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*f[1]*rdxF[2]*nu[35]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 0.75*f[2]*rdxF[2]*nu[35]+0.75*rdxF[1]*f[3]*nu[18]; 
  out[8] += 0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 0.75*rdxF[1]*f[4]*nu[18]; 
  out[10] += 0.75*rdxF[2]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[2]*f[5]*nu[35]+0.75*rdxF[1]*f[6]*nu[18]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[1]*f[8]*nu[18]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[2]*f[8]*nu[35]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[14] += 0.75*rdxF[2]*f[9]*nu[35]+0.75*rdxF[1]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[2]*f[12]*nu[35]+0.75*rdxF[1]*f[13]*nu[18]+0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0])+rdxF[1]*(0.25*nu[16])+rdxF[2]*(0.25*nu[32]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[2] += 0.75*f[0]*rdxF[0]*nu[18]; 
  out[4] += 0.75*f[0]*rdxF[1]*nu[52]; 
  out[5] += 0.75*rdxF[0]*f[1]*nu[18]; 
  out[7] += 0.75*rdxF[0]*f[3]*nu[18]; 
  out[8] += 0.75*f[1]*rdxF[1]*nu[52]; 
  out[9] += 0.75*rdxF[1]*f[2]*nu[52]+0.75*rdxF[0]*f[4]*nu[18]; 
  out[10] += 0.75*rdxF[1]*f[3]*nu[52]; 
  out[11] += 0.75*rdxF[0]*f[6]*nu[18]; 
  out[12] += 0.75*rdxF[1]*f[5]*nu[52]+0.75*rdxF[0]*f[8]*nu[18]; 
  out[13] += 0.75*rdxF[1]*f[6]*nu[52]; 
  out[14] += 0.75*rdxF[1]*f[7]*nu[52]+0.75*rdxF[0]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[1]*f[11]*nu[52]+0.75*rdxF[0]*f[13]*nu[18]; 

  return (rdxF[0]*(0.25*nu[16])+rdxF[1]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 0.75*f[0]*rdxF[1]*nu[18]; 
  out[4] += 0.75*f[0]*rdxF[2]*nu[52]; 
  out[5] += 0.75*f[1]*rdxF[1]*nu[18]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 0.75*rdxF[1]*f[3]*nu[18]; 
  out[8] += 0.75*f[1]*rdxF[2]*nu[52]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 0.75*f[2]*rdxF[2]*nu[52]+0.75*rdxF[1]*f[4]*nu[18]; 
  out[10] += 0.75*rdxF[2]*f[3]*nu[52]; 
  out[11] += 0.75*rdxF[1]*f[6]*nu[18]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[2]*f[5]*nu[52]+0.75*rdxF[1]*f[8]*nu[18]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[2]*f[6]*nu[52]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[14] += 0.75*rdxF[2]*f[7]*nu[52]+0.75*rdxF[1]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[2]*f[11]*nu[52]+0.75*rdxF[1]*f[13]*nu[18]+0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0])+rdxF[1]*(0.25*nu[16])+rdxF[2]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 0.75*f[0]*rdxF[1]*nu[18]; 
  out[5] += 0.75*f[1]*rdxF[1]*nu[18]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 0.75*rdxF[1]*f[3]*nu[18]; 
  out[8] += 0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 0.75*rdxF[1]*f[4]*nu[18]; 
  out[11] += 0.75*rdxF[1]*f[6]*nu[18]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[1]*f[8]*nu[18]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[0]*nu[1]*f[10]; 
  out[14] += 0.75*rdxF[1]*f[10]*nu[18]; 
  out[15] += 0.75*rdxF[1]*f[13]*nu[18]+0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0])+rdxF[1]*(0.25*nu[16]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 

  out[3] += 0.75*f[0]*rdxF[0]*nu[35]; 
  out[6] += 0.75*rdxF[0]*f[1]*nu[35]; 
  out[7] += 0.75*rdxF[0]*f[2]*nu[35]; 
  out[10] += 0.75*rdxF[0]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[0]*f[5]*nu[35]; 
  out[13] += 0.75*rdxF[0]*f[8]*nu[35]; 
  out[14] += 0.75*rdxF[0]*f[9]*nu[35]; 
  out[15] += 0.75*rdxF[0]*f[12]*nu[35]; 

  return (rdxF[0]*(0.25*nu[32]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[3] += 0.75*f[0]*rdxF[0]*nu[35]; 
  out[4] += 0.75*f[0]*rdxF[1]*nu[52]; 
  out[6] += 0.75*rdxF[0]*f[1]*nu[35]; 
  out[7] += 0.75*rdxF[0]*f[2]*nu[35]; 
  out[8] += 0.75*f[1]*rdxF[1]*nu[52]; 
  out[9] += 0.75*rdxF[1]*f[2]*nu[52]; 
  out[10] += 0.75*rdxF[1]*f[3]*nu[52]+0.75*rdxF[0]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[0]*f[5]*nu[35]; 
  out[12] += 0.75*rdxF[1]*f[5]*nu[52]; 
  out[13] += 0.75*rdxF[1]*f[6]*nu[52]+0.75*rdxF[0]*f[8]*nu[35]; 
  out[14] += 0.75*rdxF[1]*f[7]*nu[52]+0.75*rdxF[0]*f[9]*nu[35]; 
  out[15] += 0.75*rdxF[1]*f[11]*nu[52]+0.75*rdxF[0]*f[12]*nu[35]; 

  return (rdxF[0]*(0.25*nu[32])+rdxF[1]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.75*f[0]*rdxF[1]*nu[35]; 
  out[4] += 0.75*f[0]*rdxF[2]*nu[52]; 
  out[5] += 0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*f[1]*rdxF[1]*nu[35]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 0.75*rdxF[1]*f[2]*nu[35]; 
  out[8] += 0.75*f[1]*rdxF[2]*nu[52]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 0.75*f[2]*rdxF[2]*nu[52]; 
  out[10] += 0.75*rdxF[2]*f[3]*nu[52]+0.75*rdxF[1]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[1]*f[5]*nu[35]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[2]*f[5]*nu[52]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[2]*f[6]*nu[52]+0.75*rdxF[1]*f[8]*nu[35]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[14] += 0.75*rdxF[2]*f[7]*nu[52]+0.75*rdxF[1]*f[9]*nu[35]; 
  out[15] += 0.75*rdxF[2]*f[11]*nu[52]+0.75*rdxF[1]*f[12]*nu[35]+0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0])+rdxF[1]*(0.25*nu[32])+rdxF[2]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[3] += 0.75*f[0]*rdxF[1]*nu[35]; 
  out[5] += 0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*f[1]*rdxF[1]*nu[35]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 0.75*rdxF[1]*f[2]*nu[35]; 
  out[8] += 0.75*rdxF[0]*nu[1]*f[4]; 
  out[10] += 0.75*rdxF[1]*f[4]*nu[35]; 
  out[11] += 0.75*rdxF[1]*f[5]*nu[35]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[1]*f[8]*nu[35]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[14] += 0.75*rdxF[1]*f[9]*nu[35]; 
  out[15] += 0.75*rdxF[1]*f[12]*nu[35]+0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0])+rdxF[1]*(0.25*nu[32]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[3]*dx[3]); 

  out[4] += 0.75*f[0]*rdxF[0]*nu[52]; 
  out[8] += 0.75*rdxF[0]*f[1]*nu[52]; 
  out[9] += 0.75*rdxF[0]*f[2]*nu[52]; 
  out[10] += 0.75*rdxF[0]*f[3]*nu[52]; 
  out[12] += 0.75*rdxF[0]*f[5]*nu[52]; 
  out[13] += 0.75*rdxF[0]*f[6]*nu[52]; 
  out[14] += 0.75*rdxF[0]*f[7]*nu[52]; 
  out[15] += 0.75*rdxF[0]*f[11]*nu[52]; 

  return (rdxF[0]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[4] += 0.75*f[0]*rdxF[1]*nu[52]; 
  out[5] += 0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.75*f[1]*rdxF[1]*nu[52]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 0.75*rdxF[1]*f[2]*nu[52]; 
  out[10] += 0.75*rdxF[1]*f[3]*nu[52]; 
  out[11] += 0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[1]*f[5]*nu[52]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[1]*f[6]*nu[52]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[14] += 0.75*rdxF[1]*f[7]*nu[52]; 
  out[15] += 0.75*rdxF[1]*f[11]*nu[52]+0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0])+rdxF[1]*(0.25*nu[48]))*0.6666666666666666;

} 
double ConstDiffusionVarCoeffVol4xTensorP1_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[64]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 0.75*f[0]*rdxF[0]*nu[1]; 
  out[5] += 0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 0.75*rdxF[0]*nu[1]*f[3]; 
  out[8] += 0.75*rdxF[0]*nu[1]*f[4]; 
  out[11] += 0.75*rdxF[0]*nu[1]*f[7]; 
  out[12] += 0.75*rdxF[0]*nu[1]*f[9]; 
  out[13] += 0.75*rdxF[0]*nu[1]*f[10]; 
  out[15] += 0.75*rdxF[0]*nu[1]*f[14]; 

  return (rdxF[0]*(0.25*nu[0]))*0.6666666666666666;

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16])+kxSq[1]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16])+kxSq[1]*(0.25*nu[32])+kxSq[2]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16])+kxSq[2]*(0.25*nu[32])+kxSq[3]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16])+kxSq[2]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16])+kxSq[1]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16])+kxSq[2]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[32])+kxSq[1]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[32])+kxSq[2]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstDiffusionCFLfreqMin4xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16])+kxSq[1]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16])+kxSq[1]*(0.25*nu[32])+kxSq[2]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16])+kxSq[2]*(0.25*nu[32])+kxSq[3]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16])+kxSq[2]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[16])+kxSq[1]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16])+kxSq[2]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[16]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[32])+kxSq[1]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[32])+kxSq[2]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[32]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0])+kxSq[1]*(0.25*nu[48]))*1.333333333333333);

} 
void ConstHyperDiffusion4CFLfreqMin4xTensorP1_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[64]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.25*nu[0]))*1.333333333333333);

} 
