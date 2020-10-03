#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol4xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 

  out[12] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[20] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[26] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[33] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[36] += 6.708203932499369*rdxSq4nu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[45] += 6.708203932499369*rdxSq4nu[0]*f[17]; 

  return (rdxSq4nu[0])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[2]*dx[2]); 

  out[12] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[20] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[24] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[26] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[33] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[34] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[36] += 6.708203932499369*rdxSq4nu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[1]*f[9]; 
  out[45] += 6.708203932499369*rdxSq4nu[0]*f[17]; 
  out[46] += 6.708203932499369*rdxSq4nu[1]*f[16]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[2]*dx[2]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dx[3]*dx[3]); 

  out[12] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[14] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[20] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[24] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[26] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[29] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[30] += 6.708203932499369*rdxSq4nu[2]*f[3]; 
  out[33] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[34] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[36] += 6.708203932499369*rdxSq4nu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[1]*f[9]; 
  out[41] += 6.708203932499369*rdxSq4nu[2]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[2]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[2]*f[7]; 
  out[45] += 6.708203932499369*rdxSq4nu[0]*f[17]; 
  out[46] += 6.708203932499369*rdxSq4nu[1]*f[16]; 
  out[47] += 6.708203932499369*rdxSq4nu[2]*f[15]; 

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[4]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxSq4nu[3] = 4.0*nu[3]/(dx[3]*dx[3]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[14] += 6.708203932499369*f[0]*rdxSq4nu[3]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[24] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[27] += 6.708203932499369*rdxSq4nu[2]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxSq4nu[3]; 
  out[29] += 6.708203932499369*f[2]*rdxSq4nu[3]; 
  out[30] += 6.708203932499369*f[3]*rdxSq4nu[3]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxSq4nu[1]*f[6]; 
  out[34] += 6.708203932499369*rdxSq4nu[2]*f[5]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxSq4nu[1]*f[10]; 
  out[39] += 6.708203932499369*rdxSq4nu[2]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[2]*f[9]; 
  out[41] += 6.708203932499369*rdxSq4nu[3]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[3]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[3]*f[7]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxSq4nu[1]*f[17]; 
  out[46] += 6.708203932499369*rdxSq4nu[2]*f[16]; 
  out[47] += 6.708203932499369*rdxSq4nu[3]*f[15]; 

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2]+rdxSq4nu[3])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[24] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[27] += 6.708203932499369*rdxSq4nu[2]*f[4]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxSq4nu[1]*f[6]; 
  out[34] += 6.708203932499369*rdxSq4nu[2]*f[5]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxSq4nu[1]*f[10]; 
  out[39] += 6.708203932499369*rdxSq4nu[2]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[2]*f[9]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxSq4nu[1]*f[17]; 
  out[46] += 6.708203932499369*rdxSq4nu[2]*f[16]; 

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[1]*dx[1]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[3]*dx[3]); 

  out[12] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[14] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[20] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[26] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[29] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[30] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[33] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[36] += 6.708203932499369*rdxSq4nu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[41] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[1]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[1]*f[7]; 
  out[45] += 6.708203932499369*rdxSq4nu[0]*f[17]; 
  out[47] += 6.708203932499369*rdxSq4nu[1]*f[15]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dx[3]*dx[3]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[14] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[29] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[30] += 6.708203932499369*rdxSq4nu[2]*f[3]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxSq4nu[1]*f[6]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxSq4nu[1]*f[10]; 
  out[41] += 6.708203932499369*rdxSq4nu[2]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[2]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[2]*f[7]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxSq4nu[1]*f[17]; 
  out[47] += 6.708203932499369*rdxSq4nu[2]*f[15]; 

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxSq4nu[1]*f[6]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxSq4nu[1]*f[10]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxSq4nu[1]*f[17]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[2]*dx[2]); 

  out[13] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[23] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[24] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[27] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[34] += 6.708203932499369*rdxSq4nu[0]*f[5]; 
  out[39] += 6.708203932499369*rdxSq4nu[0]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[46] += 6.708203932499369*rdxSq4nu[0]*f[16]; 

  return (rdxSq4nu[0])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[2]*dx[2]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[3]*dx[3]); 

  out[13] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[14] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[23] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[24] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[27] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[29] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[30] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[34] += 6.708203932499369*rdxSq4nu[0]*f[5]; 
  out[39] += 6.708203932499369*rdxSq4nu[0]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[41] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[1]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[1]*f[7]; 
  out[46] += 6.708203932499369*rdxSq4nu[0]*f[16]; 
  out[47] += 6.708203932499369*rdxSq4nu[1]*f[15]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[2]*dx[2]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dx[3]*dx[3]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[14] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[24] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[29] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[30] += 6.708203932499369*rdxSq4nu[2]*f[3]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[34] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[1]*f[9]; 
  out[41] += 6.708203932499369*rdxSq4nu[2]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[2]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[2]*f[7]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 
  out[46] += 6.708203932499369*rdxSq4nu[1]*f[16]; 
  out[47] += 6.708203932499369*rdxSq4nu[2]*f[15]; 

  return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[2]*dx[2]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[24] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxSq4nu[1]*f[4]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[34] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxSq4nu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxSq4nu[1]*f[9]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 
  out[46] += 6.708203932499369*rdxSq4nu[1]*f[16]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[3]*dx[3]); 

  out[14] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[28] += 6.708203932499369*rdxSq4nu[0]*f[1]; 
  out[29] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[30] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[41] += 6.708203932499369*rdxSq4nu[0]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[47] += 6.708203932499369*rdxSq4nu[0]*f[15]; 

  return (rdxSq4nu[0])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[2]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dx[3]*dx[3]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[14] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[29] += 6.708203932499369*rdxSq4nu[1]*f[2]; 
  out[30] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[41] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[42] += 6.708203932499369*rdxSq4nu[1]*f[6]; 
  out[43] += 6.708203932499369*rdxSq4nu[1]*f[7]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 
  out[47] += 6.708203932499369*rdxSq4nu[1]*f[15]; 

  return (rdxSq4nu[0]+rdxSq4nu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxSq4nu[1]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[11] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[19] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[25] += 6.708203932499369*rdxSq4nu[0]*f[4]; 
  out[32] += 6.708203932499369*rdxSq4nu[0]*f[7]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxSq4nu[0]*f[10]; 
  out[44] += 6.708203932499369*rdxSq4nu[0]*f[18]; 

  return (rdxSq4nu[0])*0.9;

} 
