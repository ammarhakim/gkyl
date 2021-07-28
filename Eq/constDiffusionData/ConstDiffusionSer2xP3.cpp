#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol2xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[8] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[10] += 22.9128784747792*rdxFnu[0]*f[3]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVol2xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[8] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[9] += 22.9128784747792*rdxFnu[1]*f[2]; 
  out[10] += 22.9128784747792*rdxFnu[0]*f[3]; 
  out[11] += 22.9128784747792*rdxFnu[1]*f[3]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol2xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[9] += 22.9128784747792*rdxFnu[0]*f[2]; 
  out[11] += 22.9128784747792*rdxFnu[0]*f[3]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol2xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol2xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol2xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol2xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol2xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol2xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[2]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol2xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 5.1234753829798*rdxF[0]*f[4]*nu[8]+2.29128784747792*f[0]*rdxF[0]*nu[8]+3.354101966249685*rdxF[0]*f[1]*nu[4]+1.5*f[0]*rdxF[0]*nu[1]; 
  out[3] += 5.1234753829798*rdxF[0]*f[6]*nu[8]+2.29128784747792*rdxF[0]*f[2]*nu[8]+3.354101966249685*rdxF[0]*f[3]*nu[4]+1.5*rdxF[0]*nu[1]*f[2]; 
  out[4] += 13.41640786499874*rdxF[0]*f[8]*nu[8]+15.3704261489394*rdxF[0]*f[1]*nu[8]+10.06230589874905*rdxF[0]*f[4]*nu[4]+7.5*f[0]*rdxF[0]*nu[4]+6.708203932499369*rdxF[0]*f[1]*nu[1]+3.354101966249685*f[0]*nu[0]*rdxF[0]; 
  out[6] += 13.41640786499874*rdxF[0]*nu[8]*f[10]+15.3704261489394*rdxF[0]*f[3]*nu[8]+10.06230589874905*rdxF[0]*nu[4]*f[6]+7.500000000000001*rdxF[0]*f[2]*nu[4]+6.708203932499369*rdxF[0]*nu[1]*f[3]+3.354101966249684*nu[0]*rdxF[0]*f[2]; 
  out[7] += 2.29128784747792*rdxF[0]*f[5]*nu[8]+3.354101966249685*rdxF[0]*nu[4]*f[7]+1.5*rdxF[0]*nu[1]*f[5]; 
  out[8] += 36.89512162874653*rdxF[0]*f[4]*nu[8]+21.0*f[0]*rdxF[0]*nu[8]+20.12461179749811*rdxF[0]*nu[4]*f[8]+25.617376914899*rdxF[0]*f[1]*nu[4]+15.3704261489394*rdxF[0]*nu[1]*f[4]+13.74772708486752*f[0]*rdxF[0]*nu[1]+11.4564392373896*nu[0]*rdxF[0]*f[1]; 
  out[10] += 20.12461179749811*rdxF[0]*nu[4]*f[10]+36.89512162874653*rdxF[0]*f[6]*nu[8]+21.0*rdxF[0]*f[2]*nu[8]+15.3704261489394*rdxF[0]*nu[1]*f[6]+25.617376914899*rdxF[0]*f[3]*nu[4]+11.4564392373896*nu[0]*rdxF[0]*f[3]+13.74772708486752*rdxF[0]*nu[1]*f[2]; 
  out[11] += 3.354101966249685*rdxF[0]*nu[4]*f[11]+2.29128784747792*rdxF[0]*nu[8]*f[9]+1.5*rdxF[0]*nu[1]*f[9]; 

  return (rdxF[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol2xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 5.1234753829798*rdxF[0]*f[4]*nu[8]+2.29128784747792*f[0]*rdxF[0]*nu[8]+3.354101966249685*rdxF[0]*f[1]*nu[4]+1.5*f[0]*rdxF[0]*nu[1]; 
  out[2] += 5.1234753829798*rdxF[1]*f[5]*nu[21]+2.29128784747792*f[0]*rdxF[1]*nu[21]+3.354101966249685*rdxF[1]*f[2]*nu[17]+1.5*f[0]*rdxF[1]*nu[14]; 
  out[3] += 5.1234753829798*rdxF[1]*f[7]*nu[21]+2.29128784747792*f[1]*rdxF[1]*nu[21]+3.354101966249685*rdxF[1]*f[3]*nu[17]+1.5*f[1]*rdxF[1]*nu[14]+5.1234753829798*rdxF[0]*f[6]*nu[8]+2.29128784747792*rdxF[0]*f[2]*nu[8]+3.354101966249685*rdxF[0]*f[3]*nu[4]+1.5*rdxF[0]*nu[1]*f[2]; 
  out[4] += 13.41640786499874*rdxF[0]*f[8]*nu[8]+15.3704261489394*rdxF[0]*f[1]*nu[8]+10.06230589874905*rdxF[0]*f[4]*nu[4]+7.5*f[0]*rdxF[0]*nu[4]+6.708203932499369*rdxF[0]*f[1]*nu[1]+3.354101966249685*f[0]*nu[0]*rdxF[0]; 
  out[5] += 13.41640786499874*rdxF[1]*f[9]*nu[21]+15.3704261489394*rdxF[1]*f[2]*nu[21]+10.06230589874905*rdxF[1]*f[5]*nu[17]+7.5*f[0]*rdxF[1]*nu[17]+6.708203932499369*rdxF[1]*f[2]*nu[14]+3.354101966249685*f[0]*rdxF[1]*nu[12]; 
  out[6] += 2.29128784747792*rdxF[1]*f[4]*nu[21]+3.354101966249685*rdxF[1]*f[6]*nu[17]+1.5*rdxF[1]*f[4]*nu[14]+13.41640786499874*rdxF[0]*nu[8]*f[10]+15.3704261489394*rdxF[0]*f[3]*nu[8]+10.06230589874905*rdxF[0]*nu[4]*f[6]+7.500000000000001*rdxF[0]*f[2]*nu[4]+6.708203932499369*rdxF[0]*nu[1]*f[3]+3.354101966249684*nu[0]*rdxF[0]*f[2]; 
  out[7] += 13.41640786499874*rdxF[1]*f[11]*nu[21]+15.3704261489394*rdxF[1]*f[3]*nu[21]+10.06230589874905*rdxF[1]*f[7]*nu[17]+7.500000000000001*f[1]*rdxF[1]*nu[17]+6.708203932499369*rdxF[1]*f[3]*nu[14]+3.354101966249684*f[1]*rdxF[1]*nu[12]+2.29128784747792*rdxF[0]*f[5]*nu[8]+3.354101966249685*rdxF[0]*nu[4]*f[7]+1.5*rdxF[0]*nu[1]*f[5]; 
  out[8] += 36.89512162874653*rdxF[0]*f[4]*nu[8]+21.0*f[0]*rdxF[0]*nu[8]+20.12461179749811*rdxF[0]*nu[4]*f[8]+25.617376914899*rdxF[0]*f[1]*nu[4]+15.3704261489394*rdxF[0]*nu[1]*f[4]+13.74772708486752*f[0]*rdxF[0]*nu[1]+11.4564392373896*nu[0]*rdxF[0]*f[1]; 
  out[9] += 36.89512162874653*rdxF[1]*f[5]*nu[21]+21.0*f[0]*rdxF[1]*nu[21]+20.12461179749811*rdxF[1]*f[9]*nu[17]+25.617376914899*rdxF[1]*f[2]*nu[17]+15.3704261489394*rdxF[1]*f[5]*nu[14]+13.74772708486752*f[0]*rdxF[1]*nu[14]+11.4564392373896*rdxF[1]*f[2]*nu[12]; 
  out[10] += 2.29128784747792*rdxF[1]*f[8]*nu[21]+3.354101966249685*rdxF[1]*f[10]*nu[17]+1.5*rdxF[1]*f[8]*nu[14]+20.12461179749811*rdxF[0]*nu[4]*f[10]+36.89512162874653*rdxF[0]*f[6]*nu[8]+21.0*rdxF[0]*f[2]*nu[8]+15.3704261489394*rdxF[0]*nu[1]*f[6]+25.617376914899*rdxF[0]*f[3]*nu[4]+11.4564392373896*nu[0]*rdxF[0]*f[3]+13.74772708486752*rdxF[0]*nu[1]*f[2]; 
  out[11] += 36.89512162874653*rdxF[1]*f[7]*nu[21]+21.0*f[1]*rdxF[1]*nu[21]+20.12461179749811*rdxF[1]*f[11]*nu[17]+25.617376914899*rdxF[1]*f[3]*nu[17]+15.3704261489394*rdxF[1]*f[7]*nu[14]+13.74772708486752*f[1]*rdxF[1]*nu[14]+11.4564392373896*rdxF[1]*f[3]*nu[12]+3.354101966249685*rdxF[0]*nu[4]*f[11]+2.29128784747792*rdxF[0]*nu[8]*f[9]+1.5*rdxF[0]*nu[1]*f[9]; 

  return (rdxF[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+rdxF[1]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol2xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[2]:   Cell-center coordinates.
  // dx[2]:  Cell spacing.
  // nu[24]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 5.1234753829798*rdxF[0]*f[5]*nu[21]+2.29128784747792*f[0]*rdxF[0]*nu[21]+3.354101966249685*rdxF[0]*f[2]*nu[17]+1.5*f[0]*rdxF[0]*nu[14]; 
  out[3] += 5.1234753829798*rdxF[0]*f[7]*nu[21]+2.29128784747792*rdxF[0]*f[1]*nu[21]+3.354101966249685*rdxF[0]*f[3]*nu[17]+1.5*rdxF[0]*f[1]*nu[14]; 
  out[5] += 13.41640786499874*rdxF[0]*f[9]*nu[21]+15.3704261489394*rdxF[0]*f[2]*nu[21]+10.06230589874905*rdxF[0]*f[5]*nu[17]+7.5*f[0]*rdxF[0]*nu[17]+6.708203932499369*rdxF[0]*f[2]*nu[14]+3.354101966249685*f[0]*rdxF[0]*nu[12]; 
  out[6] += 2.29128784747792*rdxF[0]*f[4]*nu[21]+3.354101966249685*rdxF[0]*f[6]*nu[17]+1.5*rdxF[0]*f[4]*nu[14]; 
  out[7] += 13.41640786499874*rdxF[0]*f[11]*nu[21]+15.3704261489394*rdxF[0]*f[3]*nu[21]+10.06230589874905*rdxF[0]*f[7]*nu[17]+7.500000000000001*rdxF[0]*f[1]*nu[17]+6.708203932499369*rdxF[0]*f[3]*nu[14]+3.354101966249684*rdxF[0]*f[1]*nu[12]; 
  out[9] += 36.89512162874653*rdxF[0]*f[5]*nu[21]+21.0*f[0]*rdxF[0]*nu[21]+20.12461179749811*rdxF[0]*f[9]*nu[17]+25.617376914899*rdxF[0]*f[2]*nu[17]+15.3704261489394*rdxF[0]*f[5]*nu[14]+13.74772708486752*f[0]*rdxF[0]*nu[14]+11.4564392373896*rdxF[0]*f[2]*nu[12]; 
  out[10] += 2.29128784747792*rdxF[0]*f[8]*nu[21]+3.354101966249685*rdxF[0]*f[10]*nu[17]+1.5*rdxF[0]*f[8]*nu[14]; 
  out[11] += 36.89512162874653*rdxF[0]*f[7]*nu[21]+21.0*rdxF[0]*f[1]*nu[21]+20.12461179749811*rdxF[0]*f[11]*nu[17]+25.617376914899*rdxF[0]*f[3]*nu[17]+15.3704261489394*rdxF[0]*f[7]*nu[14]+13.74772708486752*rdxF[0]*f[1]*nu[14]+11.4564392373896*rdxF[0]*f[3]*nu[12]; 

  return (rdxF[0]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*1.142857142857143;

} 
void ConstDiffusionCFLfreqMin2xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin2xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+kxSq[1]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin2xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin2xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin2xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+kxSq[1]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin2xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin2xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin2xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[5])-0.5590169943749475*nu[4]+0.5*nu[0])+kxSq[1]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin2xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[24]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.5590169943749475*nu[17])-0.5590169943749475*nu[16]+0.5*nu[12]))*2.285714285714286);

} 
