#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol1xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[2] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[3] += 22.9128784747792*rdxFnu[0]*f[1]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol1xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol1xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[1]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol1xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[1]:   Cell-center coordinates.
  // dx[1]:  Cell spacing.
  // nu[4]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 7.24568837309472*rdxF[0]*f[2]*nu[3]+3.24037034920393*f[0]*rdxF[0]*nu[3]+4.743416490252569*rdxF[0]*f[1]*nu[2]+2.121320343559642*f[0]*rdxF[0]*nu[1]; 
  out[2] += 18.97366596101028*rdxF[0]*f[3]*nu[3]+21.73706511928415*rdxF[0]*f[1]*nu[3]+14.23024947075771*rdxF[0]*f[2]*nu[2]+10.60660171779821*f[0]*rdxF[0]*nu[2]+9.48683298050514*rdxF[0]*f[1]*nu[1]+4.743416490252569*f[0]*nu[0]*rdxF[0]; 
  out[3] += 52.17758139277825*rdxF[0]*f[2]*nu[3]+29.698484809835*f[0]*rdxF[0]*nu[3]+28.46049894151542*rdxF[0]*nu[2]*f[3]+36.22844186547361*rdxF[0]*f[1]*nu[2]+21.73706511928415*rdxF[0]*nu[1]*f[2]+19.44222209522358*f[0]*rdxF[0]*nu[1]+16.20185174601965*nu[0]*rdxF[0]*f[1]; 

  return (rdxF[0]*(0.7071067811865475*nu[0]-0.7905694150420947*nu[2]))*1.142857142857143;

} 
void ConstDiffusionCFLfreqMin1xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[4]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.7071067811865475*nu[0]-0.7905694150420947*nu[2]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin1xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[4]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.7071067811865475*nu[0]-0.7905694150420947*nu[2]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin1xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[4]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*(0.7071067811865475*nu[0]-0.7905694150420947*nu[2]))*2.285714285714286);

} 
