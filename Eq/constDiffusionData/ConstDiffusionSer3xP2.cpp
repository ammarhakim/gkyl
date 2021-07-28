#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[18] += 6.708203932499369*rdxFnu[1]*f[5]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[18] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[19] += 6.708203932499369*rdxFnu[2]*f[4]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[7] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[11] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[13] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[16] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[17] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[19] += 6.708203932499369*rdxFnu[1]*f[4]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[8] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[12] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[18] += 6.708203932499369*rdxFnu[0]*f[5]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[8] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[9] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[12] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[14] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[16] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[18] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[19] += 6.708203932499369*rdxFnu[1]*f[4]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol3xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[9] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[15] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[16] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[4]; 

  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol3xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol3xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol3xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion4Vol3xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol3xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol3xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol3xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol3xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol3xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol3xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion6Vol3xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol3xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol3xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol3xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[3]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVarCoeffVol3xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[60]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[4] += 2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[7] += 7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[10] += 2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[15] += 2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[17] += 7.115124735378852*rdxF[0]*nu[7]*f[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[18] += 2.371708245126284*rdxF[0]*nu[7]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[19] += 2.371708245126284*rdxF[0]*nu[7]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*0.9;

} 
double ConstDiffusionVarCoeffVol3xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[60]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[2] += 2.371708245126284*rdxF[1]*f[2]*nu[28]+1.060660171779821*f[0]*rdxF[1]*nu[22]; 
  out[4] += 2.371708245126284*rdxF[1]*f[4]*nu[28]+1.060660171779821*f[1]*rdxF[1]*nu[22]+2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 2.371708245126284*rdxF[1]*f[6]*nu[28]+1.060660171779821*rdxF[1]*f[3]*nu[22]; 
  out[7] += 7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[8] += 7.115124735378852*rdxF[1]*f[8]*nu[28]+5.303300858899105*f[0]*rdxF[1]*nu[28]+4.743416490252569*rdxF[1]*f[2]*nu[22]+2.371708245126284*f[0]*rdxF[1]*nu[20]; 
  out[10] += 2.371708245126284*rdxF[1]*f[10]*nu[28]+1.060660171779821*rdxF[1]*f[5]*nu[22]+2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 2.371708245126284*rdxF[1]*f[11]*nu[28]+1.060660171779821*rdxF[1]*f[7]*nu[22]+7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 7.115124735378852*rdxF[1]*f[12]*nu[28]+5.303300858899106*f[1]*rdxF[1]*nu[28]+4.743416490252569*rdxF[1]*f[4]*nu[22]+2.371708245126284*f[1]*rdxF[1]*nu[20]+2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[14] += 7.115124735378852*rdxF[1]*f[14]*nu[28]+5.303300858899106*rdxF[1]*f[3]*nu[28]+4.743416490252569*rdxF[1]*f[6]*nu[22]+2.371708245126284*rdxF[1]*f[3]*nu[20]; 
  out[15] += 2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[16] += 2.371708245126284*rdxF[1]*f[16]*nu[28]+1.060660171779821*rdxF[1]*f[9]*nu[22]; 
  out[17] += 2.371708245126284*rdxF[1]*f[17]*nu[28]+1.060660171779821*rdxF[1]*f[13]*nu[22]+7.115124735378852*rdxF[0]*nu[7]*f[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[18] += 7.115124735378852*rdxF[1]*f[18]*nu[28]+5.303300858899105*rdxF[1]*f[5]*nu[28]+4.743416490252569*rdxF[1]*f[10]*nu[22]+2.371708245126284*rdxF[1]*f[5]*nu[20]+2.371708245126284*rdxF[0]*nu[7]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[19] += 2.371708245126284*rdxF[1]*f[19]*nu[28]+1.060660171779821*rdxF[1]*f[15]*nu[22]+2.371708245126284*rdxF[0]*nu[7]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+rdxF[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*0.9;

} 
double ConstDiffusionVarCoeffVol3xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[60]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 

  out[1] += 2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[2] += 2.371708245126284*rdxF[1]*f[2]*nu[28]+1.060660171779821*f[0]*rdxF[1]*nu[22]; 
  out[3] += 2.371708245126284*rdxF[2]*f[3]*nu[49]+1.060660171779821*f[0]*rdxF[2]*nu[43]; 
  out[4] += 2.371708245126284*rdxF[1]*f[4]*nu[28]+1.060660171779821*f[1]*rdxF[1]*nu[22]+2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 2.371708245126284*rdxF[2]*f[5]*nu[49]+1.060660171779821*f[1]*rdxF[2]*nu[43]+2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 2.371708245126284*rdxF[2]*f[6]*nu[49]+1.060660171779821*f[2]*rdxF[2]*nu[43]+2.371708245126284*rdxF[1]*f[6]*nu[28]+1.060660171779821*rdxF[1]*f[3]*nu[22]; 
  out[7] += 7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[8] += 7.115124735378852*rdxF[1]*f[8]*nu[28]+5.303300858899105*f[0]*rdxF[1]*nu[28]+4.743416490252569*rdxF[1]*f[2]*nu[22]+2.371708245126284*f[0]*rdxF[1]*nu[20]; 
  out[9] += 7.115124735378852*rdxF[2]*f[9]*nu[49]+5.303300858899105*f[0]*rdxF[2]*nu[49]+4.743416490252569*rdxF[2]*f[3]*nu[43]+2.371708245126284*f[0]*rdxF[2]*nu[40]; 
  out[10] += 2.371708245126284*rdxF[2]*f[10]*nu[49]+1.060660171779821*rdxF[2]*f[4]*nu[43]+2.371708245126284*rdxF[1]*f[10]*nu[28]+1.060660171779821*rdxF[1]*f[5]*nu[22]+2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 2.371708245126284*rdxF[1]*f[11]*nu[28]+1.060660171779821*rdxF[1]*f[7]*nu[22]+7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 7.115124735378852*rdxF[1]*f[12]*nu[28]+5.303300858899106*f[1]*rdxF[1]*nu[28]+4.743416490252569*rdxF[1]*f[4]*nu[22]+2.371708245126284*f[1]*rdxF[1]*nu[20]+2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 2.371708245126284*rdxF[2]*f[13]*nu[49]+1.060660171779821*rdxF[2]*f[7]*nu[43]+7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[14] += 2.371708245126284*rdxF[2]*f[14]*nu[49]+1.060660171779821*rdxF[2]*f[8]*nu[43]+7.115124735378852*rdxF[1]*f[14]*nu[28]+5.303300858899106*rdxF[1]*f[3]*nu[28]+4.743416490252569*rdxF[1]*f[6]*nu[22]+2.371708245126284*rdxF[1]*f[3]*nu[20]; 
  out[15] += 7.115124735378852*rdxF[2]*f[15]*nu[49]+5.303300858899106*f[1]*rdxF[2]*nu[49]+4.743416490252569*rdxF[2]*f[5]*nu[43]+2.371708245126284*f[1]*rdxF[2]*nu[40]+2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[16] += 7.115124735378852*rdxF[2]*f[16]*nu[49]+5.303300858899106*f[2]*rdxF[2]*nu[49]+4.743416490252569*rdxF[2]*f[6]*nu[43]+2.371708245126284*f[2]*rdxF[2]*nu[40]+2.371708245126284*rdxF[1]*f[16]*nu[28]+1.060660171779821*rdxF[1]*f[9]*nu[22]; 
  out[17] += 2.371708245126284*rdxF[2]*f[17]*nu[49]+1.060660171779821*rdxF[2]*f[11]*nu[43]+2.371708245126284*rdxF[1]*f[17]*nu[28]+1.060660171779821*rdxF[1]*f[13]*nu[22]+7.115124735378852*rdxF[0]*nu[7]*f[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[18] += 2.371708245126284*rdxF[2]*f[18]*nu[49]+1.060660171779821*rdxF[2]*f[12]*nu[43]+7.115124735378852*rdxF[1]*f[18]*nu[28]+5.303300858899105*rdxF[1]*f[5]*nu[28]+4.743416490252569*rdxF[1]*f[10]*nu[22]+2.371708245126284*rdxF[1]*f[5]*nu[20]+2.371708245126284*rdxF[0]*nu[7]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[19] += 7.115124735378852*rdxF[2]*f[19]*nu[49]+5.303300858899105*rdxF[2]*f[4]*nu[49]+4.743416490252569*rdxF[2]*f[10]*nu[43]+2.371708245126284*rdxF[2]*f[4]*nu[40]+2.371708245126284*rdxF[1]*f[19]*nu[28]+1.060660171779821*rdxF[1]*f[15]*nu[22]+2.371708245126284*rdxF[0]*nu[7]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+rdxF[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+rdxF[2]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*0.9;

} 
double ConstDiffusionVarCoeffVol3xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[60]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[1] += 2.371708245126284*rdxF[0]*f[1]*nu[7]+1.060660171779821*f[0]*rdxF[0]*nu[1]; 
  out[3] += 2.371708245126284*rdxF[1]*f[3]*nu[49]+1.060660171779821*f[0]*rdxF[1]*nu[43]; 
  out[4] += 2.371708245126284*rdxF[0]*f[4]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[2]; 
  out[5] += 2.371708245126284*rdxF[1]*f[5]*nu[49]+1.060660171779821*f[1]*rdxF[1]*nu[43]+2.371708245126284*rdxF[0]*f[5]*nu[7]+1.060660171779821*rdxF[0]*nu[1]*f[3]; 
  out[6] += 2.371708245126284*rdxF[1]*f[6]*nu[49]+1.060660171779821*rdxF[1]*f[2]*nu[43]; 
  out[7] += 7.115124735378852*rdxF[0]*f[7]*nu[7]+5.303300858899105*f[0]*rdxF[0]*nu[7]+4.743416490252569*rdxF[0]*f[1]*nu[1]+2.371708245126284*f[0]*nu[0]*rdxF[0]; 
  out[9] += 7.115124735378852*rdxF[1]*f[9]*nu[49]+5.303300858899105*f[0]*rdxF[1]*nu[49]+4.743416490252569*rdxF[1]*f[3]*nu[43]+2.371708245126284*f[0]*rdxF[1]*nu[40]; 
  out[10] += 2.371708245126284*rdxF[1]*f[10]*nu[49]+1.060660171779821*rdxF[1]*f[4]*nu[43]+2.371708245126284*rdxF[0]*nu[7]*f[10]+1.060660171779821*rdxF[0]*nu[1]*f[6]; 
  out[11] += 7.115124735378852*rdxF[0]*nu[7]*f[11]+5.303300858899106*rdxF[0]*f[2]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[4]+2.371708245126284*nu[0]*rdxF[0]*f[2]; 
  out[12] += 2.371708245126284*rdxF[0]*nu[7]*f[12]+1.060660171779821*rdxF[0]*nu[1]*f[8]; 
  out[13] += 2.371708245126284*rdxF[1]*f[13]*nu[49]+1.060660171779821*rdxF[1]*f[7]*nu[43]+7.115124735378852*rdxF[0]*nu[7]*f[13]+5.303300858899106*rdxF[0]*f[3]*nu[7]+4.743416490252569*rdxF[0]*nu[1]*f[5]+2.371708245126284*nu[0]*rdxF[0]*f[3]; 
  out[14] += 2.371708245126284*rdxF[1]*f[14]*nu[49]+1.060660171779821*rdxF[1]*f[8]*nu[43]; 
  out[15] += 7.115124735378852*rdxF[1]*f[15]*nu[49]+5.303300858899106*f[1]*rdxF[1]*nu[49]+4.743416490252569*rdxF[1]*f[5]*nu[43]+2.371708245126284*f[1]*rdxF[1]*nu[40]+2.371708245126284*rdxF[0]*nu[7]*f[15]+1.060660171779821*rdxF[0]*nu[1]*f[9]; 
  out[16] += 7.115124735378852*rdxF[1]*f[16]*nu[49]+5.303300858899106*rdxF[1]*f[2]*nu[49]+4.743416490252569*rdxF[1]*f[6]*nu[43]+2.371708245126284*rdxF[1]*f[2]*nu[40]; 
  out[17] += 2.371708245126284*rdxF[1]*f[17]*nu[49]+1.060660171779821*rdxF[1]*f[11]*nu[43]+7.115124735378852*rdxF[0]*nu[7]*f[17]+4.743416490252569*rdxF[0]*nu[1]*f[10]+5.303300858899105*rdxF[0]*f[6]*nu[7]+2.371708245126284*nu[0]*rdxF[0]*f[6]; 
  out[18] += 2.371708245126284*rdxF[1]*f[18]*nu[49]+1.060660171779821*rdxF[1]*f[12]*nu[43]+2.371708245126284*rdxF[0]*nu[7]*f[18]+1.060660171779821*rdxF[0]*nu[1]*f[14]; 
  out[19] += 7.115124735378852*rdxF[1]*f[19]*nu[49]+5.303300858899105*rdxF[1]*f[4]*nu[49]+4.743416490252569*rdxF[1]*f[10]*nu[43]+2.371708245126284*rdxF[1]*f[4]*nu[40]+2.371708245126284*rdxF[0]*nu[7]*f[19]+1.060660171779821*rdxF[0]*nu[1]*f[16]; 

  return (rdxF[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+rdxF[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*0.9;

} 
double ConstDiffusionVarCoeffVol3xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[60]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 2.371708245126284*rdxF[0]*f[2]*nu[28]+1.060660171779821*f[0]*rdxF[0]*nu[22]; 
  out[4] += 2.371708245126284*rdxF[0]*f[4]*nu[28]+1.060660171779821*rdxF[0]*f[1]*nu[22]; 
  out[6] += 2.371708245126284*rdxF[0]*f[6]*nu[28]+1.060660171779821*rdxF[0]*f[3]*nu[22]; 
  out[8] += 7.115124735378852*rdxF[0]*f[8]*nu[28]+5.303300858899105*f[0]*rdxF[0]*nu[28]+4.743416490252569*rdxF[0]*f[2]*nu[22]+2.371708245126284*f[0]*rdxF[0]*nu[20]; 
  out[10] += 2.371708245126284*rdxF[0]*f[10]*nu[28]+1.060660171779821*rdxF[0]*f[5]*nu[22]; 
  out[11] += 2.371708245126284*rdxF[0]*f[11]*nu[28]+1.060660171779821*rdxF[0]*f[7]*nu[22]; 
  out[12] += 7.115124735378852*rdxF[0]*f[12]*nu[28]+5.303300858899106*rdxF[0]*f[1]*nu[28]+4.743416490252569*rdxF[0]*f[4]*nu[22]+2.371708245126284*rdxF[0]*f[1]*nu[20]; 
  out[14] += 7.115124735378852*rdxF[0]*f[14]*nu[28]+5.303300858899106*rdxF[0]*f[3]*nu[28]+4.743416490252569*rdxF[0]*f[6]*nu[22]+2.371708245126284*rdxF[0]*f[3]*nu[20]; 
  out[16] += 2.371708245126284*rdxF[0]*f[16]*nu[28]+1.060660171779821*rdxF[0]*f[9]*nu[22]; 
  out[17] += 2.371708245126284*rdxF[0]*f[17]*nu[28]+1.060660171779821*rdxF[0]*f[13]*nu[22]; 
  out[18] += 7.115124735378852*rdxF[0]*f[18]*nu[28]+5.303300858899105*rdxF[0]*f[5]*nu[28]+4.743416490252569*rdxF[0]*f[10]*nu[22]+2.371708245126284*rdxF[0]*f[5]*nu[20]; 
  out[19] += 2.371708245126284*rdxF[0]*f[19]*nu[28]+1.060660171779821*rdxF[0]*f[15]*nu[22]; 

  return (rdxF[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*0.9;

} 
double ConstDiffusionVarCoeffVol3xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[60]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[2] += 2.371708245126284*rdxF[0]*f[2]*nu[28]+1.060660171779821*f[0]*rdxF[0]*nu[22]; 
  out[3] += 2.371708245126284*rdxF[1]*f[3]*nu[49]+1.060660171779821*f[0]*rdxF[1]*nu[43]; 
  out[4] += 2.371708245126284*rdxF[0]*f[4]*nu[28]+1.060660171779821*rdxF[0]*f[1]*nu[22]; 
  out[5] += 2.371708245126284*rdxF[1]*f[5]*nu[49]+1.060660171779821*f[1]*rdxF[1]*nu[43]; 
  out[6] += 2.371708245126284*rdxF[1]*f[6]*nu[49]+1.060660171779821*rdxF[1]*f[2]*nu[43]+2.371708245126284*rdxF[0]*f[6]*nu[28]+1.060660171779821*rdxF[0]*f[3]*nu[22]; 
  out[8] += 7.115124735378852*rdxF[0]*f[8]*nu[28]+5.303300858899105*f[0]*rdxF[0]*nu[28]+4.743416490252569*rdxF[0]*f[2]*nu[22]+2.371708245126284*f[0]*rdxF[0]*nu[20]; 
  out[9] += 7.115124735378852*rdxF[1]*f[9]*nu[49]+5.303300858899105*f[0]*rdxF[1]*nu[49]+4.743416490252569*rdxF[1]*f[3]*nu[43]+2.371708245126284*f[0]*rdxF[1]*nu[40]; 
  out[10] += 2.371708245126284*rdxF[1]*f[10]*nu[49]+1.060660171779821*rdxF[1]*f[4]*nu[43]+2.371708245126284*rdxF[0]*f[10]*nu[28]+1.060660171779821*rdxF[0]*f[5]*nu[22]; 
  out[11] += 2.371708245126284*rdxF[0]*f[11]*nu[28]+1.060660171779821*rdxF[0]*f[7]*nu[22]; 
  out[12] += 7.115124735378852*rdxF[0]*f[12]*nu[28]+5.303300858899106*rdxF[0]*f[1]*nu[28]+4.743416490252569*rdxF[0]*f[4]*nu[22]+2.371708245126284*rdxF[0]*f[1]*nu[20]; 
  out[13] += 2.371708245126284*rdxF[1]*f[13]*nu[49]+1.060660171779821*rdxF[1]*f[7]*nu[43]; 
  out[14] += 2.371708245126284*rdxF[1]*f[14]*nu[49]+1.060660171779821*rdxF[1]*f[8]*nu[43]+7.115124735378852*rdxF[0]*f[14]*nu[28]+5.303300858899106*rdxF[0]*f[3]*nu[28]+4.743416490252569*rdxF[0]*f[6]*nu[22]+2.371708245126284*rdxF[0]*f[3]*nu[20]; 
  out[15] += 7.115124735378852*rdxF[1]*f[15]*nu[49]+5.303300858899106*f[1]*rdxF[1]*nu[49]+4.743416490252569*rdxF[1]*f[5]*nu[43]+2.371708245126284*f[1]*rdxF[1]*nu[40]; 
  out[16] += 7.115124735378852*rdxF[1]*f[16]*nu[49]+5.303300858899106*rdxF[1]*f[2]*nu[49]+4.743416490252569*rdxF[1]*f[6]*nu[43]+2.371708245126284*rdxF[1]*f[2]*nu[40]+2.371708245126284*rdxF[0]*f[16]*nu[28]+1.060660171779821*rdxF[0]*f[9]*nu[22]; 
  out[17] += 2.371708245126284*rdxF[1]*f[17]*nu[49]+1.060660171779821*rdxF[1]*f[11]*nu[43]+2.371708245126284*rdxF[0]*f[17]*nu[28]+1.060660171779821*rdxF[0]*f[13]*nu[22]; 
  out[18] += 2.371708245126284*rdxF[1]*f[18]*nu[49]+1.060660171779821*rdxF[1]*f[12]*nu[43]+7.115124735378852*rdxF[0]*f[18]*nu[28]+5.303300858899105*rdxF[0]*f[5]*nu[28]+4.743416490252569*rdxF[0]*f[10]*nu[22]+2.371708245126284*rdxF[0]*f[5]*nu[20]; 
  out[19] += 7.115124735378852*rdxF[1]*f[19]*nu[49]+5.303300858899105*rdxF[1]*f[4]*nu[49]+4.743416490252569*rdxF[1]*f[10]*nu[43]+2.371708245126284*rdxF[1]*f[4]*nu[40]+2.371708245126284*rdxF[0]*f[19]*nu[28]+1.060660171779821*rdxF[0]*f[15]*nu[22]; 

  return (rdxF[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+rdxF[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*0.9;

} 
double ConstDiffusionVarCoeffVol3xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[3]:   Cell-center coordinates.
  // dx[3]:  Cell spacing.
  // nu[60]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 

  out[3] += 2.371708245126284*rdxF[0]*f[3]*nu[49]+1.060660171779821*f[0]*rdxF[0]*nu[43]; 
  out[5] += 2.371708245126284*rdxF[0]*f[5]*nu[49]+1.060660171779821*rdxF[0]*f[1]*nu[43]; 
  out[6] += 2.371708245126284*rdxF[0]*f[6]*nu[49]+1.060660171779821*rdxF[0]*f[2]*nu[43]; 
  out[9] += 7.115124735378852*rdxF[0]*f[9]*nu[49]+5.303300858899105*f[0]*rdxF[0]*nu[49]+4.743416490252569*rdxF[0]*f[3]*nu[43]+2.371708245126284*f[0]*rdxF[0]*nu[40]; 
  out[10] += 2.371708245126284*rdxF[0]*f[10]*nu[49]+1.060660171779821*rdxF[0]*f[4]*nu[43]; 
  out[13] += 2.371708245126284*rdxF[0]*f[13]*nu[49]+1.060660171779821*rdxF[0]*f[7]*nu[43]; 
  out[14] += 2.371708245126284*rdxF[0]*f[14]*nu[49]+1.060660171779821*rdxF[0]*f[8]*nu[43]; 
  out[15] += 7.115124735378852*rdxF[0]*f[15]*nu[49]+5.303300858899106*rdxF[0]*f[1]*nu[49]+4.743416490252569*rdxF[0]*f[5]*nu[43]+2.371708245126284*rdxF[0]*f[1]*nu[40]; 
  out[16] += 7.115124735378852*rdxF[0]*f[16]*nu[49]+5.303300858899106*rdxF[0]*f[2]*nu[49]+4.743416490252569*rdxF[0]*f[6]*nu[43]+2.371708245126284*rdxF[0]*f[2]*nu[40]; 
  out[17] += 2.371708245126284*rdxF[0]*f[17]*nu[49]+1.060660171779821*rdxF[0]*f[11]*nu[43]; 
  out[18] += 2.371708245126284*rdxF[0]*f[18]*nu[49]+1.060660171779821*rdxF[0]*f[12]*nu[43]; 
  out[19] += 7.115124735378852*rdxF[0]*f[19]*nu[49]+5.303300858899105*rdxF[0]*f[4]*nu[49]+4.743416490252569*rdxF[0]*f[10]*nu[43]+2.371708245126284*rdxF[0]*f[4]*nu[40]; 

  return (rdxF[0]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*0.9;

} 
void ConstDiffusionCFLfreqMin3xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*1.8);

} 
void ConstDiffusionCFLfreqMin3xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*1.8);

} 
void ConstDiffusionCFLfreqMin3xSerP2_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+kxSq[2]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstDiffusionCFLfreqMin3xSerP2_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstDiffusionCFLfreqMin3xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*1.8);

} 
void ConstDiffusionCFLfreqMin3xSerP2_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+kxSq[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstDiffusionCFLfreqMin3xSerP2_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP2_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+kxSq[2]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP2_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP2_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+kxSq[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin3xSerP2_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP2_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+kxSq[2]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP2_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[9])-0.3952847075210473*nu[8]-0.3952847075210473*nu[7]+0.3535533905932737*nu[0])+kxSq[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP2_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[29])-0.3952847075210473*nu[28]-0.3952847075210473*nu[27]+0.3535533905932737*nu[20])+kxSq[1]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin3xSerP2_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[60]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.3952847075210473*nu[49])-0.3952847075210473*nu[48]-0.3952847075210473*nu[47]+0.3535533905932737*nu[40]))*1.8);

} 
