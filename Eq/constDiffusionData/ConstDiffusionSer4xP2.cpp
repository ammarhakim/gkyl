#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol4xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[12] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[33] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[45] += 6.708203932499369*rdxFnu[0]*f[17]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[12] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[24] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[33] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[34] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[45] += 6.708203932499369*rdxFnu[0]*f[17]; 
  out[46] += 6.708203932499369*rdxFnu[1]*f[16]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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

  out[12] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[14] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[24] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[29] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[30] += 6.708203932499369*rdxFnu[2]*f[3]; 
  out[33] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[34] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[2]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[2]*f[7]; 
  out[45] += 6.708203932499369*rdxFnu[0]*f[17]; 
  out[46] += 6.708203932499369*rdxFnu[1]*f[16]; 
  out[47] += 6.708203932499369*rdxFnu[2]*f[15]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[13] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[14] += 6.708203932499369*f[0]*rdxFnu[3]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[24] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[27] += 6.708203932499369*rdxFnu[2]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxFnu[3]; 
  out[29] += 6.708203932499369*f[2]*rdxFnu[3]; 
  out[30] += 6.708203932499369*f[3]*rdxFnu[3]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[34] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[39] += 6.708203932499369*rdxFnu[2]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[2]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[3]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[3]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[3]*f[7]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxFnu[1]*f[17]; 
  out[46] += 6.708203932499369*rdxFnu[2]*f[16]; 
  out[47] += 6.708203932499369*rdxFnu[3]*f[15]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[13] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[24] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[27] += 6.708203932499369*rdxFnu[2]*f[4]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[34] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[39] += 6.708203932499369*rdxFnu[2]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[2]*f[9]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxFnu[1]*f[17]; 
  out[46] += 6.708203932499369*rdxFnu[2]*f[16]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[1]/(dx[1]*dx[1]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 

  out[12] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[14] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[20] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[22] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[26] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[29] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[30] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[33] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[38] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[41] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[1]*f[7]; 
  out[45] += 6.708203932499369*rdxFnu[0]*f[17]; 
  out[47] += 6.708203932499369*rdxFnu[1]*f[15]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[14] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[29] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[30] += 6.708203932499369*rdxFnu[2]*f[3]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[41] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[2]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[2]*f[7]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxFnu[1]*f[17]; 
  out[47] += 6.708203932499369*rdxFnu[2]*f[15]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[1]/(dx[1]*dx[1]); 

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[12] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[20] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[22] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[26] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[33] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[36] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[45] += 6.708203932499369*rdxFnu[1]*f[17]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[13] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[24] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[27] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[34] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[46] += 6.708203932499369*rdxFnu[0]*f[16]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[2]/(dx[2]*dx[2]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 

  out[13] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[14] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[23] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[24] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[27] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[29] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[30] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[34] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[1]*f[7]; 
  out[46] += 6.708203932499369*rdxFnu[0]*f[16]; 
  out[47] += 6.708203932499369*rdxFnu[1]*f[15]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[14] += 6.708203932499369*f[0]*rdxFnu[2]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[24] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxFnu[2]; 
  out[29] += 6.708203932499369*f[2]*rdxFnu[2]; 
  out[30] += 6.708203932499369*rdxFnu[2]*f[3]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[34] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[2]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[2]*f[7]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[46] += 6.708203932499369*rdxFnu[1]*f[16]; 
  out[47] += 6.708203932499369*rdxFnu[2]*f[15]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[2]/(dx[2]*dx[2]); 

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[13] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[23] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[24] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[27] += 6.708203932499369*rdxFnu[1]*f[4]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[34] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[39] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[46] += 6.708203932499369*rdxFnu[1]*f[16]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[3]/(dx[3]*dx[3]); 

  out[14] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[28] += 6.708203932499369*rdxFnu[0]*f[1]; 
  out[29] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[30] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[47] += 6.708203932499369*rdxFnu[0]*f[15]; 

  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 
  rdxFnu[1] = 4.0*nu[3]/(dx[3]*dx[3]); 

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[14] += 6.708203932499369*f[0]*rdxFnu[1]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[28] += 6.708203932499369*f[1]*rdxFnu[1]; 
  out[29] += 6.708203932499369*rdxFnu[1]*f[2]; 
  out[30] += 6.708203932499369*rdxFnu[1]*f[3]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[41] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[42] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[43] += 6.708203932499369*rdxFnu[1]*f[7]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[47] += 6.708203932499369*rdxFnu[1]*f[15]; 

  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstDiffusionVol4xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 4.0*nu[0]/(dx[0]*dx[0]); 

  out[11] += 6.708203932499369*f[0]*rdxFnu[0]; 
  out[19] += 6.708203932499369*rdxFnu[0]*f[2]; 
  out[21] += 6.708203932499369*rdxFnu[0]*f[3]; 
  out[25] += 6.708203932499369*rdxFnu[0]*f[4]; 
  out[32] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[35] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[18]; 

  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion4Vol4xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[4]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[3] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[2] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[3]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[2] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*0.9;

} 
double ConstHyperDiffusion6Vol4xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 1.677050983124842*rdxF[0]*f[2]*nu[60]+0.75*f[0]*rdxF[0]*nu[50]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[60]+0.75*rdxF[0]*f[1]*nu[50]; 
  out[7] += 1.677050983124842*rdxF[0]*f[7]*nu[60]+0.75*rdxF[0]*f[3]*nu[50]; 
  out[9] += 1.677050983124842*rdxF[0]*f[9]*nu[60]+0.75*rdxF[0]*f[4]*nu[50]; 
  out[12] += 5.031152949374527*rdxF[0]*f[12]*nu[60]+3.75*f[0]*rdxF[0]*nu[60]+3.354101966249685*rdxF[0]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[0]*nu[48]; 
  out[15] += 1.677050983124842*rdxF[0]*f[15]*nu[60]+0.75*rdxF[0]*f[6]*nu[50]; 
  out[16] += 1.677050983124842*rdxF[0]*f[16]*nu[60]+0.75*rdxF[0]*f[8]*nu[50]; 
  out[18] += 1.677050983124842*rdxF[0]*f[18]*nu[60]+0.75*rdxF[0]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[0]*f[19]*nu[60]+0.75*rdxF[0]*f[11]*nu[50]; 
  out[20] += 5.031152949374527*rdxF[0]*f[20]*nu[60]+3.75*rdxF[0]*f[1]*nu[60]+3.354101966249684*rdxF[0]*f[5]*nu[50]+1.677050983124842*rdxF[0]*f[1]*nu[48]; 
  out[22] += 5.031152949374527*rdxF[0]*f[22]*nu[60]+3.75*rdxF[0]*f[3]*nu[60]+3.354101966249684*rdxF[0]*f[7]*nu[50]+1.677050983124842*rdxF[0]*f[3]*nu[48]; 
  out[24] += 1.677050983124842*rdxF[0]*f[24]*nu[60]+0.75*rdxF[0]*f[13]*nu[50]; 
  out[26] += 5.031152949374527*rdxF[0]*f[26]*nu[60]+3.75*rdxF[0]*f[4]*nu[60]+3.354101966249684*rdxF[0]*f[9]*nu[50]+1.677050983124842*rdxF[0]*f[4]*nu[48]; 
  out[29] += 1.677050983124842*rdxF[0]*f[29]*nu[60]+0.75*rdxF[0]*f[14]*nu[50]; 
  out[31] += 1.677050983124842*rdxF[0]*f[31]*nu[60]+0.75*rdxF[0]*f[17]*nu[50]; 
  out[32] += 1.677050983124842*rdxF[0]*f[32]*nu[60]+0.75*rdxF[0]*f[21]*nu[50]; 
  out[33] += 5.031152949374527*rdxF[0]*f[33]*nu[60]+3.75*rdxF[0]*f[6]*nu[60]+3.354101966249685*rdxF[0]*f[15]*nu[50]+1.677050983124842*rdxF[0]*f[6]*nu[48]; 
  out[34] += 1.677050983124842*rdxF[0]*f[34]*nu[60]+0.75*rdxF[0]*f[23]*nu[50]; 
  out[35] += 1.677050983124842*rdxF[0]*f[35]*nu[60]+0.75*rdxF[0]*f[25]*nu[50]; 
  out[36] += 5.031152949374527*rdxF[0]*f[36]*nu[60]+3.75*rdxF[0]*f[8]*nu[60]+3.354101966249685*rdxF[0]*f[16]*nu[50]+1.677050983124842*rdxF[0]*f[8]*nu[48]; 
  out[38] += 5.031152949374527*rdxF[0]*f[38]*nu[60]+3.75*rdxF[0]*f[10]*nu[60]+3.354101966249685*rdxF[0]*f[18]*nu[50]+1.677050983124842*rdxF[0]*f[10]*nu[48]; 
  out[40] += 1.677050983124842*rdxF[0]*f[40]*nu[60]+0.75*rdxF[0]*f[27]*nu[50]; 
  out[41] += 1.677050983124842*rdxF[0]*f[41]*nu[60]+0.75*rdxF[0]*f[28]*nu[50]; 
  out[43] += 1.677050983124842*rdxF[0]*f[43]*nu[60]+0.75*rdxF[0]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[0]*f[44]*nu[60]+0.75*rdxF[0]*f[37]*nu[50]; 
  out[45] += 5.031152949374527*rdxF[0]*f[45]*nu[60]+3.75*rdxF[0]*f[17]*nu[60]+3.354101966249684*rdxF[0]*f[31]*nu[50]+1.677050983124842*rdxF[0]*f[17]*nu[48]; 
  out[46] += 1.677050983124842*rdxF[0]*f[46]*nu[60]+0.75*rdxF[0]*f[39]*nu[50]; 
  out[47] += 1.677050983124842*rdxF[0]*f[47]*nu[60]+0.75*rdxF[0]*f[42]*nu[50]; 

  return (rdxF[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[2] += 1.677050983124842*rdxF[0]*f[2]*nu[60]+0.75*f[0]*rdxF[0]*nu[50]; 
  out[3] += 1.677050983124842*rdxF[1]*f[3]*nu[109]+0.75*f[0]*rdxF[1]*nu[99]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[60]+0.75*rdxF[0]*f[1]*nu[50]; 
  out[6] += 1.677050983124842*rdxF[1]*f[6]*nu[109]+0.75*f[1]*rdxF[1]*nu[99]; 
  out[7] += 1.677050983124842*rdxF[1]*f[7]*nu[109]+0.75*rdxF[1]*f[2]*nu[99]+1.677050983124842*rdxF[0]*f[7]*nu[60]+0.75*rdxF[0]*f[3]*nu[50]; 
  out[9] += 1.677050983124842*rdxF[0]*f[9]*nu[60]+0.75*rdxF[0]*f[4]*nu[50]; 
  out[10] += 1.677050983124842*rdxF[1]*f[10]*nu[109]+0.75*rdxF[1]*f[4]*nu[99]; 
  out[12] += 5.031152949374527*rdxF[0]*f[12]*nu[60]+3.75*f[0]*rdxF[0]*nu[60]+3.354101966249685*rdxF[0]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[0]*nu[48]; 
  out[13] += 5.031152949374527*rdxF[1]*f[13]*nu[109]+3.75*f[0]*rdxF[1]*nu[109]+3.354101966249685*rdxF[1]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[1]*nu[96]; 
  out[15] += 1.677050983124842*rdxF[1]*f[15]*nu[109]+0.75*rdxF[1]*f[5]*nu[99]+1.677050983124842*rdxF[0]*f[15]*nu[60]+0.75*rdxF[0]*f[6]*nu[50]; 
  out[16] += 1.677050983124842*rdxF[0]*f[16]*nu[60]+0.75*rdxF[0]*f[8]*nu[50]; 
  out[17] += 1.677050983124842*rdxF[1]*f[17]*nu[109]+0.75*rdxF[1]*f[8]*nu[99]; 
  out[18] += 1.677050983124842*rdxF[1]*f[18]*nu[109]+0.75*rdxF[1]*f[9]*nu[99]+1.677050983124842*rdxF[0]*f[18]*nu[60]+0.75*rdxF[0]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[0]*f[19]*nu[60]+0.75*rdxF[0]*f[11]*nu[50]; 
  out[20] += 5.031152949374527*rdxF[0]*f[20]*nu[60]+3.75*rdxF[0]*f[1]*nu[60]+3.354101966249684*rdxF[0]*f[5]*nu[50]+1.677050983124842*rdxF[0]*f[1]*nu[48]; 
  out[21] += 1.677050983124842*rdxF[1]*f[21]*nu[109]+0.75*rdxF[1]*f[11]*nu[99]; 
  out[22] += 1.677050983124842*rdxF[1]*f[22]*nu[109]+0.75*rdxF[1]*f[12]*nu[99]+5.031152949374527*rdxF[0]*f[22]*nu[60]+3.75*rdxF[0]*f[3]*nu[60]+3.354101966249684*rdxF[0]*f[7]*nu[50]+1.677050983124842*rdxF[0]*f[3]*nu[48]; 
  out[23] += 5.031152949374527*rdxF[1]*f[23]*nu[109]+3.75*f[1]*rdxF[1]*nu[109]+3.354101966249684*rdxF[1]*f[6]*nu[99]+1.677050983124842*f[1]*rdxF[1]*nu[96]; 
  out[24] += 5.031152949374527*rdxF[1]*f[24]*nu[109]+3.75*rdxF[1]*f[2]*nu[109]+3.354101966249684*rdxF[1]*f[7]*nu[99]+1.677050983124842*rdxF[1]*f[2]*nu[96]+1.677050983124842*rdxF[0]*f[24]*nu[60]+0.75*rdxF[0]*f[13]*nu[50]; 
  out[26] += 5.031152949374527*rdxF[0]*f[26]*nu[60]+3.75*rdxF[0]*f[4]*nu[60]+3.354101966249684*rdxF[0]*f[9]*nu[50]+1.677050983124842*rdxF[0]*f[4]*nu[48]; 
  out[27] += 5.031152949374527*rdxF[1]*f[27]*nu[109]+3.75*rdxF[1]*f[4]*nu[109]+3.354101966249684*rdxF[1]*f[10]*nu[99]+1.677050983124842*rdxF[1]*f[4]*nu[96]; 
  out[29] += 1.677050983124842*rdxF[0]*f[29]*nu[60]+0.75*rdxF[0]*f[14]*nu[50]; 
  out[30] += 1.677050983124842*rdxF[1]*f[30]*nu[109]+0.75*rdxF[1]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[1]*f[31]*nu[109]+0.75*rdxF[1]*f[16]*nu[99]+1.677050983124842*rdxF[0]*f[31]*nu[60]+0.75*rdxF[0]*f[17]*nu[50]; 
  out[32] += 1.677050983124842*rdxF[1]*f[32]*nu[109]+0.75*rdxF[1]*f[19]*nu[99]+1.677050983124842*rdxF[0]*f[32]*nu[60]+0.75*rdxF[0]*f[21]*nu[50]; 
  out[33] += 1.677050983124842*rdxF[1]*f[33]*nu[109]+0.75*rdxF[1]*f[20]*nu[99]+5.031152949374527*rdxF[0]*f[33]*nu[60]+3.75*rdxF[0]*f[6]*nu[60]+3.354101966249685*rdxF[0]*f[15]*nu[50]+1.677050983124842*rdxF[0]*f[6]*nu[48]; 
  out[34] += 5.031152949374527*rdxF[1]*f[34]*nu[109]+3.75*rdxF[1]*f[5]*nu[109]+3.354101966249685*rdxF[1]*f[15]*nu[99]+1.677050983124842*rdxF[1]*f[5]*nu[96]+1.677050983124842*rdxF[0]*f[34]*nu[60]+0.75*rdxF[0]*f[23]*nu[50]; 
  out[35] += 1.677050983124842*rdxF[0]*f[35]*nu[60]+0.75*rdxF[0]*f[25]*nu[50]; 
  out[36] += 5.031152949374527*rdxF[0]*f[36]*nu[60]+3.75*rdxF[0]*f[8]*nu[60]+3.354101966249685*rdxF[0]*f[16]*nu[50]+1.677050983124842*rdxF[0]*f[8]*nu[48]; 
  out[37] += 1.677050983124842*rdxF[1]*f[37]*nu[109]+0.75*rdxF[1]*f[25]*nu[99]; 
  out[38] += 1.677050983124842*rdxF[1]*f[38]*nu[109]+0.75*rdxF[1]*f[26]*nu[99]+5.031152949374527*rdxF[0]*f[38]*nu[60]+3.75*rdxF[0]*f[10]*nu[60]+3.354101966249685*rdxF[0]*f[18]*nu[50]+1.677050983124842*rdxF[0]*f[10]*nu[48]; 
  out[39] += 5.031152949374527*rdxF[1]*f[39]*nu[109]+3.75*rdxF[1]*f[8]*nu[109]+3.354101966249685*rdxF[1]*f[17]*nu[99]+1.677050983124842*rdxF[1]*f[8]*nu[96]; 
  out[40] += 5.031152949374527*rdxF[1]*f[40]*nu[109]+3.75*rdxF[1]*f[9]*nu[109]+3.354101966249685*rdxF[1]*f[18]*nu[99]+1.677050983124842*rdxF[1]*f[9]*nu[96]+1.677050983124842*rdxF[0]*f[40]*nu[60]+0.75*rdxF[0]*f[27]*nu[50]; 
  out[41] += 1.677050983124842*rdxF[0]*f[41]*nu[60]+0.75*rdxF[0]*f[28]*nu[50]; 
  out[42] += 1.677050983124842*rdxF[1]*f[42]*nu[109]+0.75*rdxF[1]*f[28]*nu[99]; 
  out[43] += 1.677050983124842*rdxF[1]*f[43]*nu[109]+0.75*rdxF[1]*f[29]*nu[99]+1.677050983124842*rdxF[0]*f[43]*nu[60]+0.75*rdxF[0]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[1]*f[44]*nu[109]+0.75*rdxF[1]*f[35]*nu[99]+1.677050983124842*rdxF[0]*f[44]*nu[60]+0.75*rdxF[0]*f[37]*nu[50]; 
  out[45] += 1.677050983124842*rdxF[1]*f[45]*nu[109]+0.75*rdxF[1]*f[36]*nu[99]+5.031152949374527*rdxF[0]*f[45]*nu[60]+3.75*rdxF[0]*f[17]*nu[60]+3.354101966249684*rdxF[0]*f[31]*nu[50]+1.677050983124842*rdxF[0]*f[17]*nu[48]; 
  out[46] += 5.031152949374527*rdxF[1]*f[46]*nu[109]+3.75*rdxF[1]*f[16]*nu[109]+3.354101966249684*rdxF[1]*f[31]*nu[99]+1.677050983124842*rdxF[1]*f[16]*nu[96]+1.677050983124842*rdxF[0]*f[46]*nu[60]+0.75*rdxF[0]*f[39]*nu[50]; 
  out[47] += 1.677050983124842*rdxF[1]*f[47]*nu[109]+0.75*rdxF[1]*f[41]*nu[99]+1.677050983124842*rdxF[0]*f[47]*nu[60]+0.75*rdxF[0]*f[42]*nu[50]; 

  return (rdxF[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+rdxF[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[2] += 1.677050983124842*rdxF[0]*f[2]*nu[60]+0.75*f[0]*rdxF[0]*nu[50]; 
  out[3] += 1.677050983124842*rdxF[1]*f[3]*nu[109]+0.75*f[0]*rdxF[1]*nu[99]; 
  out[4] += 1.677050983124842*rdxF[2]*f[4]*nu[158]+0.75*f[0]*rdxF[2]*nu[148]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[60]+0.75*rdxF[0]*f[1]*nu[50]; 
  out[6] += 1.677050983124842*rdxF[1]*f[6]*nu[109]+0.75*f[1]*rdxF[1]*nu[99]; 
  out[7] += 1.677050983124842*rdxF[1]*f[7]*nu[109]+0.75*rdxF[1]*f[2]*nu[99]+1.677050983124842*rdxF[0]*f[7]*nu[60]+0.75*rdxF[0]*f[3]*nu[50]; 
  out[8] += 1.677050983124842*rdxF[2]*f[8]*nu[158]+0.75*f[1]*rdxF[2]*nu[148]; 
  out[9] += 1.677050983124842*rdxF[2]*f[9]*nu[158]+0.75*f[2]*rdxF[2]*nu[148]+1.677050983124842*rdxF[0]*f[9]*nu[60]+0.75*rdxF[0]*f[4]*nu[50]; 
  out[10] += 1.677050983124842*rdxF[2]*f[10]*nu[158]+0.75*rdxF[2]*f[3]*nu[148]+1.677050983124842*rdxF[1]*f[10]*nu[109]+0.75*rdxF[1]*f[4]*nu[99]; 
  out[12] += 5.031152949374527*rdxF[0]*f[12]*nu[60]+3.75*f[0]*rdxF[0]*nu[60]+3.354101966249685*rdxF[0]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[0]*nu[48]; 
  out[13] += 5.031152949374527*rdxF[1]*f[13]*nu[109]+3.75*f[0]*rdxF[1]*nu[109]+3.354101966249685*rdxF[1]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[1]*nu[96]; 
  out[14] += 5.031152949374527*rdxF[2]*f[14]*nu[158]+3.75*f[0]*rdxF[2]*nu[158]+3.354101966249685*rdxF[2]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[2]*nu[144]; 
  out[15] += 1.677050983124842*rdxF[1]*f[15]*nu[109]+0.75*rdxF[1]*f[5]*nu[99]+1.677050983124842*rdxF[0]*f[15]*nu[60]+0.75*rdxF[0]*f[6]*nu[50]; 
  out[16] += 1.677050983124842*rdxF[2]*f[16]*nu[158]+0.75*rdxF[2]*f[5]*nu[148]+1.677050983124842*rdxF[0]*f[16]*nu[60]+0.75*rdxF[0]*f[8]*nu[50]; 
  out[17] += 1.677050983124842*rdxF[2]*f[17]*nu[158]+0.75*rdxF[2]*f[6]*nu[148]+1.677050983124842*rdxF[1]*f[17]*nu[109]+0.75*rdxF[1]*f[8]*nu[99]; 
  out[18] += 1.677050983124842*rdxF[2]*f[18]*nu[158]+0.75*rdxF[2]*f[7]*nu[148]+1.677050983124842*rdxF[1]*f[18]*nu[109]+0.75*rdxF[1]*f[9]*nu[99]+1.677050983124842*rdxF[0]*f[18]*nu[60]+0.75*rdxF[0]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[0]*f[19]*nu[60]+0.75*rdxF[0]*f[11]*nu[50]; 
  out[20] += 5.031152949374527*rdxF[0]*f[20]*nu[60]+3.75*rdxF[0]*f[1]*nu[60]+3.354101966249684*rdxF[0]*f[5]*nu[50]+1.677050983124842*rdxF[0]*f[1]*nu[48]; 
  out[21] += 1.677050983124842*rdxF[1]*f[21]*nu[109]+0.75*rdxF[1]*f[11]*nu[99]; 
  out[22] += 1.677050983124842*rdxF[1]*f[22]*nu[109]+0.75*rdxF[1]*f[12]*nu[99]+5.031152949374527*rdxF[0]*f[22]*nu[60]+3.75*rdxF[0]*f[3]*nu[60]+3.354101966249684*rdxF[0]*f[7]*nu[50]+1.677050983124842*rdxF[0]*f[3]*nu[48]; 
  out[23] += 5.031152949374527*rdxF[1]*f[23]*nu[109]+3.75*f[1]*rdxF[1]*nu[109]+3.354101966249684*rdxF[1]*f[6]*nu[99]+1.677050983124842*f[1]*rdxF[1]*nu[96]; 
  out[24] += 5.031152949374527*rdxF[1]*f[24]*nu[109]+3.75*rdxF[1]*f[2]*nu[109]+3.354101966249684*rdxF[1]*f[7]*nu[99]+1.677050983124842*rdxF[1]*f[2]*nu[96]+1.677050983124842*rdxF[0]*f[24]*nu[60]+0.75*rdxF[0]*f[13]*nu[50]; 
  out[25] += 1.677050983124842*rdxF[2]*f[25]*nu[158]+0.75*rdxF[2]*f[11]*nu[148]; 
  out[26] += 1.677050983124842*rdxF[2]*f[26]*nu[158]+0.75*rdxF[2]*f[12]*nu[148]+5.031152949374527*rdxF[0]*f[26]*nu[60]+3.75*rdxF[0]*f[4]*nu[60]+3.354101966249684*rdxF[0]*f[9]*nu[50]+1.677050983124842*rdxF[0]*f[4]*nu[48]; 
  out[27] += 1.677050983124842*rdxF[2]*f[27]*nu[158]+0.75*rdxF[2]*f[13]*nu[148]+5.031152949374527*rdxF[1]*f[27]*nu[109]+3.75*rdxF[1]*f[4]*nu[109]+3.354101966249684*rdxF[1]*f[10]*nu[99]+1.677050983124842*rdxF[1]*f[4]*nu[96]; 
  out[28] += 5.031152949374527*rdxF[2]*f[28]*nu[158]+3.75*f[1]*rdxF[2]*nu[158]+3.354101966249684*rdxF[2]*f[8]*nu[148]+1.677050983124842*f[1]*rdxF[2]*nu[144]; 
  out[29] += 5.031152949374527*rdxF[2]*f[29]*nu[158]+3.75*f[2]*rdxF[2]*nu[158]+3.354101966249684*rdxF[2]*f[9]*nu[148]+1.677050983124842*f[2]*rdxF[2]*nu[144]+1.677050983124842*rdxF[0]*f[29]*nu[60]+0.75*rdxF[0]*f[14]*nu[50]; 
  out[30] += 5.031152949374527*rdxF[2]*f[30]*nu[158]+3.75*rdxF[2]*f[3]*nu[158]+3.354101966249684*rdxF[2]*f[10]*nu[148]+1.677050983124842*rdxF[2]*f[3]*nu[144]+1.677050983124842*rdxF[1]*f[30]*nu[109]+0.75*rdxF[1]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[2]*f[31]*nu[158]+0.75*rdxF[2]*f[15]*nu[148]+1.677050983124842*rdxF[1]*f[31]*nu[109]+0.75*rdxF[1]*f[16]*nu[99]+1.677050983124842*rdxF[0]*f[31]*nu[60]+0.75*rdxF[0]*f[17]*nu[50]; 
  out[32] += 1.677050983124842*rdxF[1]*f[32]*nu[109]+0.75*rdxF[1]*f[19]*nu[99]+1.677050983124842*rdxF[0]*f[32]*nu[60]+0.75*rdxF[0]*f[21]*nu[50]; 
  out[33] += 1.677050983124842*rdxF[1]*f[33]*nu[109]+0.75*rdxF[1]*f[20]*nu[99]+5.031152949374527*rdxF[0]*f[33]*nu[60]+3.75*rdxF[0]*f[6]*nu[60]+3.354101966249685*rdxF[0]*f[15]*nu[50]+1.677050983124842*rdxF[0]*f[6]*nu[48]; 
  out[34] += 5.031152949374527*rdxF[1]*f[34]*nu[109]+3.75*rdxF[1]*f[5]*nu[109]+3.354101966249685*rdxF[1]*f[15]*nu[99]+1.677050983124842*rdxF[1]*f[5]*nu[96]+1.677050983124842*rdxF[0]*f[34]*nu[60]+0.75*rdxF[0]*f[23]*nu[50]; 
  out[35] += 1.677050983124842*rdxF[2]*f[35]*nu[158]+0.75*rdxF[2]*f[19]*nu[148]+1.677050983124842*rdxF[0]*f[35]*nu[60]+0.75*rdxF[0]*f[25]*nu[50]; 
  out[36] += 1.677050983124842*rdxF[2]*f[36]*nu[158]+0.75*rdxF[2]*f[20]*nu[148]+5.031152949374527*rdxF[0]*f[36]*nu[60]+3.75*rdxF[0]*f[8]*nu[60]+3.354101966249685*rdxF[0]*f[16]*nu[50]+1.677050983124842*rdxF[0]*f[8]*nu[48]; 
  out[37] += 1.677050983124842*rdxF[2]*f[37]*nu[158]+0.75*rdxF[2]*f[21]*nu[148]+1.677050983124842*rdxF[1]*f[37]*nu[109]+0.75*rdxF[1]*f[25]*nu[99]; 
  out[38] += 1.677050983124842*rdxF[2]*f[38]*nu[158]+0.75*rdxF[2]*f[22]*nu[148]+1.677050983124842*rdxF[1]*f[38]*nu[109]+0.75*rdxF[1]*f[26]*nu[99]+5.031152949374527*rdxF[0]*f[38]*nu[60]+3.75*rdxF[0]*f[10]*nu[60]+3.354101966249685*rdxF[0]*f[18]*nu[50]+1.677050983124842*rdxF[0]*f[10]*nu[48]; 
  out[39] += 1.677050983124842*rdxF[2]*f[39]*nu[158]+0.75*rdxF[2]*f[23]*nu[148]+5.031152949374527*rdxF[1]*f[39]*nu[109]+3.75*rdxF[1]*f[8]*nu[109]+3.354101966249685*rdxF[1]*f[17]*nu[99]+1.677050983124842*rdxF[1]*f[8]*nu[96]; 
  out[40] += 1.677050983124842*rdxF[2]*f[40]*nu[158]+0.75*rdxF[2]*f[24]*nu[148]+5.031152949374527*rdxF[1]*f[40]*nu[109]+3.75*rdxF[1]*f[9]*nu[109]+3.354101966249685*rdxF[1]*f[18]*nu[99]+1.677050983124842*rdxF[1]*f[9]*nu[96]+1.677050983124842*rdxF[0]*f[40]*nu[60]+0.75*rdxF[0]*f[27]*nu[50]; 
  out[41] += 5.031152949374527*rdxF[2]*f[41]*nu[158]+3.75*rdxF[2]*f[5]*nu[158]+3.354101966249685*rdxF[2]*f[16]*nu[148]+1.677050983124842*rdxF[2]*f[5]*nu[144]+1.677050983124842*rdxF[0]*f[41]*nu[60]+0.75*rdxF[0]*f[28]*nu[50]; 
  out[42] += 5.031152949374527*rdxF[2]*f[42]*nu[158]+3.75*rdxF[2]*f[6]*nu[158]+3.354101966249685*rdxF[2]*f[17]*nu[148]+1.677050983124842*rdxF[2]*f[6]*nu[144]+1.677050983124842*rdxF[1]*f[42]*nu[109]+0.75*rdxF[1]*f[28]*nu[99]; 
  out[43] += 5.031152949374527*rdxF[2]*f[43]*nu[158]+3.75*rdxF[2]*f[7]*nu[158]+3.354101966249685*rdxF[2]*f[18]*nu[148]+1.677050983124842*rdxF[2]*f[7]*nu[144]+1.677050983124842*rdxF[1]*f[43]*nu[109]+0.75*rdxF[1]*f[29]*nu[99]+1.677050983124842*rdxF[0]*f[43]*nu[60]+0.75*rdxF[0]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[2]*f[44]*nu[158]+0.75*rdxF[2]*f[32]*nu[148]+1.677050983124842*rdxF[1]*f[44]*nu[109]+0.75*rdxF[1]*f[35]*nu[99]+1.677050983124842*rdxF[0]*f[44]*nu[60]+0.75*rdxF[0]*f[37]*nu[50]; 
  out[45] += 1.677050983124842*rdxF[2]*f[45]*nu[158]+0.75*rdxF[2]*f[33]*nu[148]+1.677050983124842*rdxF[1]*f[45]*nu[109]+0.75*rdxF[1]*f[36]*nu[99]+5.031152949374527*rdxF[0]*f[45]*nu[60]+3.75*rdxF[0]*f[17]*nu[60]+3.354101966249684*rdxF[0]*f[31]*nu[50]+1.677050983124842*rdxF[0]*f[17]*nu[48]; 
  out[46] += 1.677050983124842*rdxF[2]*f[46]*nu[158]+0.75*rdxF[2]*f[34]*nu[148]+5.031152949374527*rdxF[1]*f[46]*nu[109]+3.75*rdxF[1]*f[16]*nu[109]+3.354101966249684*rdxF[1]*f[31]*nu[99]+1.677050983124842*rdxF[1]*f[16]*nu[96]+1.677050983124842*rdxF[0]*f[46]*nu[60]+0.75*rdxF[0]*f[39]*nu[50]; 
  out[47] += 5.031152949374527*rdxF[2]*f[47]*nu[158]+3.75*rdxF[2]*f[15]*nu[158]+3.354101966249684*rdxF[2]*f[31]*nu[148]+1.677050983124842*rdxF[2]*f[15]*nu[144]+1.677050983124842*rdxF[1]*f[47]*nu[109]+0.75*rdxF[1]*f[41]*nu[99]+1.677050983124842*rdxF[0]*f[47]*nu[60]+0.75*rdxF[0]*f[42]*nu[50]; 

  return (rdxF[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+rdxF[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+rdxF[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 
  rdxF[3] = 4.0/(dx[3]*dx[3]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 1.677050983124842*rdxF[1]*f[2]*nu[60]+0.75*f[0]*rdxF[1]*nu[50]; 
  out[3] += 1.677050983124842*rdxF[2]*f[3]*nu[109]+0.75*f[0]*rdxF[2]*nu[99]; 
  out[4] += 1.677050983124842*rdxF[3]*f[4]*nu[158]+0.75*f[0]*rdxF[3]*nu[148]; 
  out[5] += 1.677050983124842*rdxF[1]*f[5]*nu[60]+0.75*f[1]*rdxF[1]*nu[50]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[2]*f[6]*nu[109]+0.75*f[1]*rdxF[2]*nu[99]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 1.677050983124842*rdxF[2]*f[7]*nu[109]+0.75*f[2]*rdxF[2]*nu[99]+1.677050983124842*rdxF[1]*f[7]*nu[60]+0.75*rdxF[1]*f[3]*nu[50]; 
  out[8] += 1.677050983124842*rdxF[3]*f[8]*nu[158]+0.75*f[1]*rdxF[3]*nu[148]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 1.677050983124842*rdxF[3]*f[9]*nu[158]+0.75*f[2]*rdxF[3]*nu[148]+1.677050983124842*rdxF[1]*f[9]*nu[60]+0.75*rdxF[1]*f[4]*nu[50]; 
  out[10] += 1.677050983124842*rdxF[3]*f[10]*nu[158]+0.75*f[3]*rdxF[3]*nu[148]+1.677050983124842*rdxF[2]*f[10]*nu[109]+0.75*rdxF[2]*f[4]*nu[99]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 5.031152949374527*rdxF[1]*f[12]*nu[60]+3.75*f[0]*rdxF[1]*nu[60]+3.354101966249685*rdxF[1]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[1]*nu[48]; 
  out[13] += 5.031152949374527*rdxF[2]*f[13]*nu[109]+3.75*f[0]*rdxF[2]*nu[109]+3.354101966249685*rdxF[2]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[2]*nu[96]; 
  out[14] += 5.031152949374527*rdxF[3]*f[14]*nu[158]+3.75*f[0]*rdxF[3]*nu[158]+3.354101966249685*rdxF[3]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[3]*nu[144]; 
  out[15] += 1.677050983124842*rdxF[2]*f[15]*nu[109]+0.75*rdxF[2]*f[5]*nu[99]+1.677050983124842*rdxF[1]*f[15]*nu[60]+0.75*rdxF[1]*f[6]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[3]*f[16]*nu[158]+0.75*rdxF[3]*f[5]*nu[148]+1.677050983124842*rdxF[1]*f[16]*nu[60]+0.75*rdxF[1]*f[8]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[3]*f[17]*nu[158]+0.75*rdxF[3]*f[6]*nu[148]+1.677050983124842*rdxF[2]*f[17]*nu[109]+0.75*rdxF[2]*f[8]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 1.677050983124842*rdxF[3]*f[18]*nu[158]+0.75*rdxF[3]*f[7]*nu[148]+1.677050983124842*rdxF[2]*f[18]*nu[109]+0.75*rdxF[2]*f[9]*nu[99]+1.677050983124842*rdxF[1]*f[18]*nu[60]+0.75*rdxF[1]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[1]*f[19]*nu[60]+0.75*rdxF[1]*f[11]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 5.031152949374527*rdxF[1]*f[20]*nu[60]+3.75*f[1]*rdxF[1]*nu[60]+3.354101966249684*rdxF[1]*f[5]*nu[50]+1.677050983124842*f[1]*rdxF[1]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.677050983124842*rdxF[2]*f[21]*nu[109]+0.75*rdxF[2]*f[11]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.677050983124842*rdxF[2]*f[22]*nu[109]+0.75*rdxF[2]*f[12]*nu[99]+5.031152949374527*rdxF[1]*f[22]*nu[60]+3.75*rdxF[1]*f[3]*nu[60]+3.354101966249684*rdxF[1]*f[7]*nu[50]+1.677050983124842*rdxF[1]*f[3]*nu[48]; 
  out[23] += 5.031152949374527*rdxF[2]*f[23]*nu[109]+3.75*f[1]*rdxF[2]*nu[109]+3.354101966249684*rdxF[2]*f[6]*nu[99]+1.677050983124842*f[1]*rdxF[2]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 5.031152949374527*rdxF[2]*f[24]*nu[109]+3.75*f[2]*rdxF[2]*nu[109]+3.354101966249684*rdxF[2]*f[7]*nu[99]+1.677050983124842*f[2]*rdxF[2]*nu[96]+1.677050983124842*rdxF[1]*f[24]*nu[60]+0.75*rdxF[1]*f[13]*nu[50]; 
  out[25] += 1.677050983124842*rdxF[3]*f[25]*nu[158]+0.75*rdxF[3]*f[11]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.677050983124842*rdxF[3]*f[26]*nu[158]+0.75*rdxF[3]*f[12]*nu[148]+5.031152949374527*rdxF[1]*f[26]*nu[60]+3.75*rdxF[1]*f[4]*nu[60]+3.354101966249684*rdxF[1]*f[9]*nu[50]+1.677050983124842*rdxF[1]*f[4]*nu[48]; 
  out[27] += 1.677050983124842*rdxF[3]*f[27]*nu[158]+0.75*rdxF[3]*f[13]*nu[148]+5.031152949374527*rdxF[2]*f[27]*nu[109]+3.75*rdxF[2]*f[4]*nu[109]+3.354101966249684*rdxF[2]*f[10]*nu[99]+1.677050983124842*rdxF[2]*f[4]*nu[96]; 
  out[28] += 5.031152949374527*rdxF[3]*f[28]*nu[158]+3.75*f[1]*rdxF[3]*nu[158]+3.354101966249684*rdxF[3]*f[8]*nu[148]+1.677050983124842*f[1]*rdxF[3]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 5.031152949374527*rdxF[3]*f[29]*nu[158]+3.75*f[2]*rdxF[3]*nu[158]+3.354101966249684*rdxF[3]*f[9]*nu[148]+1.677050983124842*f[2]*rdxF[3]*nu[144]+1.677050983124842*rdxF[1]*f[29]*nu[60]+0.75*rdxF[1]*f[14]*nu[50]; 
  out[30] += 5.031152949374527*rdxF[3]*f[30]*nu[158]+3.75*f[3]*rdxF[3]*nu[158]+3.354101966249684*rdxF[3]*f[10]*nu[148]+1.677050983124842*f[3]*rdxF[3]*nu[144]+1.677050983124842*rdxF[2]*f[30]*nu[109]+0.75*rdxF[2]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[3]*f[31]*nu[158]+0.75*rdxF[3]*f[15]*nu[148]+1.677050983124842*rdxF[2]*f[31]*nu[109]+0.75*rdxF[2]*f[16]*nu[99]+1.677050983124842*rdxF[1]*f[31]*nu[60]+0.75*rdxF[1]*f[17]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 1.677050983124842*rdxF[2]*f[32]*nu[109]+0.75*rdxF[2]*f[19]*nu[99]+1.677050983124842*rdxF[1]*f[32]*nu[60]+0.75*rdxF[1]*f[21]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 1.677050983124842*rdxF[2]*f[33]*nu[109]+0.75*rdxF[2]*f[20]*nu[99]+5.031152949374527*rdxF[1]*f[33]*nu[60]+3.75*rdxF[1]*f[6]*nu[60]+3.354101966249685*rdxF[1]*f[15]*nu[50]+1.677050983124842*rdxF[1]*f[6]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 5.031152949374527*rdxF[2]*f[34]*nu[109]+3.75*rdxF[2]*f[5]*nu[109]+3.354101966249685*rdxF[2]*f[15]*nu[99]+1.677050983124842*rdxF[2]*f[5]*nu[96]+1.677050983124842*rdxF[1]*f[34]*nu[60]+0.75*rdxF[1]*f[23]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 1.677050983124842*rdxF[3]*f[35]*nu[158]+0.75*rdxF[3]*f[19]*nu[148]+1.677050983124842*rdxF[1]*f[35]*nu[60]+0.75*rdxF[1]*f[25]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 1.677050983124842*rdxF[3]*f[36]*nu[158]+0.75*rdxF[3]*f[20]*nu[148]+5.031152949374527*rdxF[1]*f[36]*nu[60]+3.75*rdxF[1]*f[8]*nu[60]+3.354101966249685*rdxF[1]*f[16]*nu[50]+1.677050983124842*rdxF[1]*f[8]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 1.677050983124842*rdxF[3]*f[37]*nu[158]+0.75*rdxF[3]*f[21]*nu[148]+1.677050983124842*rdxF[2]*f[37]*nu[109]+0.75*rdxF[2]*f[25]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[38] += 1.677050983124842*rdxF[3]*f[38]*nu[158]+0.75*rdxF[3]*f[22]*nu[148]+1.677050983124842*rdxF[2]*f[38]*nu[109]+0.75*rdxF[2]*f[26]*nu[99]+5.031152949374527*rdxF[1]*f[38]*nu[60]+3.75*rdxF[1]*f[10]*nu[60]+3.354101966249685*rdxF[1]*f[18]*nu[50]+1.677050983124842*rdxF[1]*f[10]*nu[48]; 
  out[39] += 1.677050983124842*rdxF[3]*f[39]*nu[158]+0.75*rdxF[3]*f[23]*nu[148]+5.031152949374527*rdxF[2]*f[39]*nu[109]+3.75*rdxF[2]*f[8]*nu[109]+3.354101966249685*rdxF[2]*f[17]*nu[99]+1.677050983124842*rdxF[2]*f[8]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[40] += 1.677050983124842*rdxF[3]*f[40]*nu[158]+0.75*rdxF[3]*f[24]*nu[148]+5.031152949374527*rdxF[2]*f[40]*nu[109]+3.75*rdxF[2]*f[9]*nu[109]+3.354101966249685*rdxF[2]*f[18]*nu[99]+1.677050983124842*rdxF[2]*f[9]*nu[96]+1.677050983124842*rdxF[1]*f[40]*nu[60]+0.75*rdxF[1]*f[27]*nu[50]; 
  out[41] += 5.031152949374527*rdxF[3]*f[41]*nu[158]+3.75*rdxF[3]*f[5]*nu[158]+3.354101966249685*rdxF[3]*f[16]*nu[148]+1.677050983124842*rdxF[3]*f[5]*nu[144]+1.677050983124842*rdxF[1]*f[41]*nu[60]+0.75*rdxF[1]*f[28]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 5.031152949374527*rdxF[3]*f[42]*nu[158]+3.75*rdxF[3]*f[6]*nu[158]+3.354101966249685*rdxF[3]*f[17]*nu[148]+1.677050983124842*rdxF[3]*f[6]*nu[144]+1.677050983124842*rdxF[2]*f[42]*nu[109]+0.75*rdxF[2]*f[28]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[43] += 5.031152949374527*rdxF[3]*f[43]*nu[158]+3.75*rdxF[3]*f[7]*nu[158]+3.354101966249685*rdxF[3]*f[18]*nu[148]+1.677050983124842*rdxF[3]*f[7]*nu[144]+1.677050983124842*rdxF[2]*f[43]*nu[109]+0.75*rdxF[2]*f[29]*nu[99]+1.677050983124842*rdxF[1]*f[43]*nu[60]+0.75*rdxF[1]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[3]*f[44]*nu[158]+0.75*rdxF[3]*f[32]*nu[148]+1.677050983124842*rdxF[2]*f[44]*nu[109]+0.75*rdxF[2]*f[35]*nu[99]+1.677050983124842*rdxF[1]*f[44]*nu[60]+0.75*rdxF[1]*f[37]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 1.677050983124842*rdxF[3]*f[45]*nu[158]+0.75*rdxF[3]*f[33]*nu[148]+1.677050983124842*rdxF[2]*f[45]*nu[109]+0.75*rdxF[2]*f[36]*nu[99]+5.031152949374527*rdxF[1]*f[45]*nu[60]+3.75*rdxF[1]*f[17]*nu[60]+3.354101966249684*rdxF[1]*f[31]*nu[50]+1.677050983124842*rdxF[1]*f[17]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 1.677050983124842*rdxF[3]*f[46]*nu[158]+0.75*rdxF[3]*f[34]*nu[148]+5.031152949374527*rdxF[2]*f[46]*nu[109]+3.75*rdxF[2]*f[16]*nu[109]+3.354101966249684*rdxF[2]*f[31]*nu[99]+1.677050983124842*rdxF[2]*f[16]*nu[96]+1.677050983124842*rdxF[1]*f[46]*nu[60]+0.75*rdxF[1]*f[39]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 5.031152949374527*rdxF[3]*f[47]*nu[158]+3.75*rdxF[3]*f[15]*nu[158]+3.354101966249684*rdxF[3]*f[31]*nu[148]+1.677050983124842*rdxF[3]*f[15]*nu[144]+1.677050983124842*rdxF[2]*f[47]*nu[109]+0.75*rdxF[2]*f[41]*nu[99]+1.677050983124842*rdxF[1]*f[47]*nu[60]+0.75*rdxF[1]*f[42]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+rdxF[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+rdxF[3]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 1.677050983124842*rdxF[1]*f[2]*nu[60]+0.75*f[0]*rdxF[1]*nu[50]; 
  out[3] += 1.677050983124842*rdxF[2]*f[3]*nu[109]+0.75*f[0]*rdxF[2]*nu[99]; 
  out[5] += 1.677050983124842*rdxF[1]*f[5]*nu[60]+0.75*f[1]*rdxF[1]*nu[50]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[2]*f[6]*nu[109]+0.75*f[1]*rdxF[2]*nu[99]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 1.677050983124842*rdxF[2]*f[7]*nu[109]+0.75*f[2]*rdxF[2]*nu[99]+1.677050983124842*rdxF[1]*f[7]*nu[60]+0.75*rdxF[1]*f[3]*nu[50]; 
  out[8] += 1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 1.677050983124842*rdxF[1]*f[9]*nu[60]+0.75*rdxF[1]*f[4]*nu[50]; 
  out[10] += 1.677050983124842*rdxF[2]*f[10]*nu[109]+0.75*rdxF[2]*f[4]*nu[99]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 5.031152949374527*rdxF[1]*f[12]*nu[60]+3.75*f[0]*rdxF[1]*nu[60]+3.354101966249685*rdxF[1]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[1]*nu[48]; 
  out[13] += 5.031152949374527*rdxF[2]*f[13]*nu[109]+3.75*f[0]*rdxF[2]*nu[109]+3.354101966249685*rdxF[2]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[2]*nu[96]; 
  out[15] += 1.677050983124842*rdxF[2]*f[15]*nu[109]+0.75*rdxF[2]*f[5]*nu[99]+1.677050983124842*rdxF[1]*f[15]*nu[60]+0.75*rdxF[1]*f[6]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[1]*f[16]*nu[60]+0.75*rdxF[1]*f[8]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[2]*f[17]*nu[109]+0.75*rdxF[2]*f[8]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 1.677050983124842*rdxF[2]*f[18]*nu[109]+0.75*rdxF[2]*f[9]*nu[99]+1.677050983124842*rdxF[1]*f[18]*nu[60]+0.75*rdxF[1]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[1]*f[19]*nu[60]+0.75*rdxF[1]*f[11]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 5.031152949374527*rdxF[1]*f[20]*nu[60]+3.75*f[1]*rdxF[1]*nu[60]+3.354101966249684*rdxF[1]*f[5]*nu[50]+1.677050983124842*f[1]*rdxF[1]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.677050983124842*rdxF[2]*f[21]*nu[109]+0.75*rdxF[2]*f[11]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.677050983124842*rdxF[2]*f[22]*nu[109]+0.75*rdxF[2]*f[12]*nu[99]+5.031152949374527*rdxF[1]*f[22]*nu[60]+3.75*rdxF[1]*f[3]*nu[60]+3.354101966249684*rdxF[1]*f[7]*nu[50]+1.677050983124842*rdxF[1]*f[3]*nu[48]; 
  out[23] += 5.031152949374527*rdxF[2]*f[23]*nu[109]+3.75*f[1]*rdxF[2]*nu[109]+3.354101966249684*rdxF[2]*f[6]*nu[99]+1.677050983124842*f[1]*rdxF[2]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 5.031152949374527*rdxF[2]*f[24]*nu[109]+3.75*f[2]*rdxF[2]*nu[109]+3.354101966249684*rdxF[2]*f[7]*nu[99]+1.677050983124842*f[2]*rdxF[2]*nu[96]+1.677050983124842*rdxF[1]*f[24]*nu[60]+0.75*rdxF[1]*f[13]*nu[50]; 
  out[25] += 5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 5.031152949374527*rdxF[1]*f[26]*nu[60]+3.75*rdxF[1]*f[4]*nu[60]+3.354101966249684*rdxF[1]*f[9]*nu[50]+1.677050983124842*rdxF[1]*f[4]*nu[48]; 
  out[27] += 5.031152949374527*rdxF[2]*f[27]*nu[109]+3.75*rdxF[2]*f[4]*nu[109]+3.354101966249684*rdxF[2]*f[10]*nu[99]+1.677050983124842*rdxF[2]*f[4]*nu[96]; 
  out[28] += 1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 1.677050983124842*rdxF[1]*f[29]*nu[60]+0.75*rdxF[1]*f[14]*nu[50]; 
  out[30] += 1.677050983124842*rdxF[2]*f[30]*nu[109]+0.75*rdxF[2]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[2]*f[31]*nu[109]+0.75*rdxF[2]*f[16]*nu[99]+1.677050983124842*rdxF[1]*f[31]*nu[60]+0.75*rdxF[1]*f[17]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 1.677050983124842*rdxF[2]*f[32]*nu[109]+0.75*rdxF[2]*f[19]*nu[99]+1.677050983124842*rdxF[1]*f[32]*nu[60]+0.75*rdxF[1]*f[21]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 1.677050983124842*rdxF[2]*f[33]*nu[109]+0.75*rdxF[2]*f[20]*nu[99]+5.031152949374527*rdxF[1]*f[33]*nu[60]+3.75*rdxF[1]*f[6]*nu[60]+3.354101966249685*rdxF[1]*f[15]*nu[50]+1.677050983124842*rdxF[1]*f[6]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 5.031152949374527*rdxF[2]*f[34]*nu[109]+3.75*rdxF[2]*f[5]*nu[109]+3.354101966249685*rdxF[2]*f[15]*nu[99]+1.677050983124842*rdxF[2]*f[5]*nu[96]+1.677050983124842*rdxF[1]*f[34]*nu[60]+0.75*rdxF[1]*f[23]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 1.677050983124842*rdxF[1]*f[35]*nu[60]+0.75*rdxF[1]*f[25]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 5.031152949374527*rdxF[1]*f[36]*nu[60]+3.75*rdxF[1]*f[8]*nu[60]+3.354101966249685*rdxF[1]*f[16]*nu[50]+1.677050983124842*rdxF[1]*f[8]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 1.677050983124842*rdxF[2]*f[37]*nu[109]+0.75*rdxF[2]*f[25]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[38] += 1.677050983124842*rdxF[2]*f[38]*nu[109]+0.75*rdxF[2]*f[26]*nu[99]+5.031152949374527*rdxF[1]*f[38]*nu[60]+3.75*rdxF[1]*f[10]*nu[60]+3.354101966249685*rdxF[1]*f[18]*nu[50]+1.677050983124842*rdxF[1]*f[10]*nu[48]; 
  out[39] += 5.031152949374527*rdxF[2]*f[39]*nu[109]+3.75*rdxF[2]*f[8]*nu[109]+3.354101966249685*rdxF[2]*f[17]*nu[99]+1.677050983124842*rdxF[2]*f[8]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[40] += 5.031152949374527*rdxF[2]*f[40]*nu[109]+3.75*rdxF[2]*f[9]*nu[109]+3.354101966249685*rdxF[2]*f[18]*nu[99]+1.677050983124842*rdxF[2]*f[9]*nu[96]+1.677050983124842*rdxF[1]*f[40]*nu[60]+0.75*rdxF[1]*f[27]*nu[50]; 
  out[41] += 1.677050983124842*rdxF[1]*f[41]*nu[60]+0.75*rdxF[1]*f[28]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 1.677050983124842*rdxF[2]*f[42]*nu[109]+0.75*rdxF[2]*f[28]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[43] += 1.677050983124842*rdxF[2]*f[43]*nu[109]+0.75*rdxF[2]*f[29]*nu[99]+1.677050983124842*rdxF[1]*f[43]*nu[60]+0.75*rdxF[1]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[2]*f[44]*nu[109]+0.75*rdxF[2]*f[35]*nu[99]+1.677050983124842*rdxF[1]*f[44]*nu[60]+0.75*rdxF[1]*f[37]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 1.677050983124842*rdxF[2]*f[45]*nu[109]+0.75*rdxF[2]*f[36]*nu[99]+5.031152949374527*rdxF[1]*f[45]*nu[60]+3.75*rdxF[1]*f[17]*nu[60]+3.354101966249684*rdxF[1]*f[31]*nu[50]+1.677050983124842*rdxF[1]*f[17]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 5.031152949374527*rdxF[2]*f[46]*nu[109]+3.75*rdxF[2]*f[16]*nu[109]+3.354101966249684*rdxF[2]*f[31]*nu[99]+1.677050983124842*rdxF[2]*f[16]*nu[96]+1.677050983124842*rdxF[1]*f[46]*nu[60]+0.75*rdxF[1]*f[39]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 1.677050983124842*rdxF[2]*f[47]*nu[109]+0.75*rdxF[2]*f[41]*nu[99]+1.677050983124842*rdxF[1]*f[47]*nu[60]+0.75*rdxF[1]*f[42]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+rdxF[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[2] += 1.677050983124842*rdxF[0]*f[2]*nu[60]+0.75*f[0]*rdxF[0]*nu[50]; 
  out[4] += 1.677050983124842*rdxF[1]*f[4]*nu[158]+0.75*f[0]*rdxF[1]*nu[148]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[60]+0.75*rdxF[0]*f[1]*nu[50]; 
  out[7] += 1.677050983124842*rdxF[0]*f[7]*nu[60]+0.75*rdxF[0]*f[3]*nu[50]; 
  out[8] += 1.677050983124842*rdxF[1]*f[8]*nu[158]+0.75*f[1]*rdxF[1]*nu[148]; 
  out[9] += 1.677050983124842*rdxF[1]*f[9]*nu[158]+0.75*rdxF[1]*f[2]*nu[148]+1.677050983124842*rdxF[0]*f[9]*nu[60]+0.75*rdxF[0]*f[4]*nu[50]; 
  out[10] += 1.677050983124842*rdxF[1]*f[10]*nu[158]+0.75*rdxF[1]*f[3]*nu[148]; 
  out[12] += 5.031152949374527*rdxF[0]*f[12]*nu[60]+3.75*f[0]*rdxF[0]*nu[60]+3.354101966249685*rdxF[0]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[0]*nu[48]; 
  out[14] += 5.031152949374527*rdxF[1]*f[14]*nu[158]+3.75*f[0]*rdxF[1]*nu[158]+3.354101966249685*rdxF[1]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[1]*nu[144]; 
  out[15] += 1.677050983124842*rdxF[0]*f[15]*nu[60]+0.75*rdxF[0]*f[6]*nu[50]; 
  out[16] += 1.677050983124842*rdxF[1]*f[16]*nu[158]+0.75*rdxF[1]*f[5]*nu[148]+1.677050983124842*rdxF[0]*f[16]*nu[60]+0.75*rdxF[0]*f[8]*nu[50]; 
  out[17] += 1.677050983124842*rdxF[1]*f[17]*nu[158]+0.75*rdxF[1]*f[6]*nu[148]; 
  out[18] += 1.677050983124842*rdxF[1]*f[18]*nu[158]+0.75*rdxF[1]*f[7]*nu[148]+1.677050983124842*rdxF[0]*f[18]*nu[60]+0.75*rdxF[0]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[0]*f[19]*nu[60]+0.75*rdxF[0]*f[11]*nu[50]; 
  out[20] += 5.031152949374527*rdxF[0]*f[20]*nu[60]+3.75*rdxF[0]*f[1]*nu[60]+3.354101966249684*rdxF[0]*f[5]*nu[50]+1.677050983124842*rdxF[0]*f[1]*nu[48]; 
  out[22] += 5.031152949374527*rdxF[0]*f[22]*nu[60]+3.75*rdxF[0]*f[3]*nu[60]+3.354101966249684*rdxF[0]*f[7]*nu[50]+1.677050983124842*rdxF[0]*f[3]*nu[48]; 
  out[24] += 1.677050983124842*rdxF[0]*f[24]*nu[60]+0.75*rdxF[0]*f[13]*nu[50]; 
  out[25] += 1.677050983124842*rdxF[1]*f[25]*nu[158]+0.75*rdxF[1]*f[11]*nu[148]; 
  out[26] += 1.677050983124842*rdxF[1]*f[26]*nu[158]+0.75*rdxF[1]*f[12]*nu[148]+5.031152949374527*rdxF[0]*f[26]*nu[60]+3.75*rdxF[0]*f[4]*nu[60]+3.354101966249684*rdxF[0]*f[9]*nu[50]+1.677050983124842*rdxF[0]*f[4]*nu[48]; 
  out[27] += 1.677050983124842*rdxF[1]*f[27]*nu[158]+0.75*rdxF[1]*f[13]*nu[148]; 
  out[28] += 5.031152949374527*rdxF[1]*f[28]*nu[158]+3.75*f[1]*rdxF[1]*nu[158]+3.354101966249684*rdxF[1]*f[8]*nu[148]+1.677050983124842*f[1]*rdxF[1]*nu[144]; 
  out[29] += 5.031152949374527*rdxF[1]*f[29]*nu[158]+3.75*rdxF[1]*f[2]*nu[158]+3.354101966249684*rdxF[1]*f[9]*nu[148]+1.677050983124842*rdxF[1]*f[2]*nu[144]+1.677050983124842*rdxF[0]*f[29]*nu[60]+0.75*rdxF[0]*f[14]*nu[50]; 
  out[30] += 5.031152949374527*rdxF[1]*f[30]*nu[158]+3.75*rdxF[1]*f[3]*nu[158]+3.354101966249684*rdxF[1]*f[10]*nu[148]+1.677050983124842*rdxF[1]*f[3]*nu[144]; 
  out[31] += 1.677050983124842*rdxF[1]*f[31]*nu[158]+0.75*rdxF[1]*f[15]*nu[148]+1.677050983124842*rdxF[0]*f[31]*nu[60]+0.75*rdxF[0]*f[17]*nu[50]; 
  out[32] += 1.677050983124842*rdxF[0]*f[32]*nu[60]+0.75*rdxF[0]*f[21]*nu[50]; 
  out[33] += 5.031152949374527*rdxF[0]*f[33]*nu[60]+3.75*rdxF[0]*f[6]*nu[60]+3.354101966249685*rdxF[0]*f[15]*nu[50]+1.677050983124842*rdxF[0]*f[6]*nu[48]; 
  out[34] += 1.677050983124842*rdxF[0]*f[34]*nu[60]+0.75*rdxF[0]*f[23]*nu[50]; 
  out[35] += 1.677050983124842*rdxF[1]*f[35]*nu[158]+0.75*rdxF[1]*f[19]*nu[148]+1.677050983124842*rdxF[0]*f[35]*nu[60]+0.75*rdxF[0]*f[25]*nu[50]; 
  out[36] += 1.677050983124842*rdxF[1]*f[36]*nu[158]+0.75*rdxF[1]*f[20]*nu[148]+5.031152949374527*rdxF[0]*f[36]*nu[60]+3.75*rdxF[0]*f[8]*nu[60]+3.354101966249685*rdxF[0]*f[16]*nu[50]+1.677050983124842*rdxF[0]*f[8]*nu[48]; 
  out[37] += 1.677050983124842*rdxF[1]*f[37]*nu[158]+0.75*rdxF[1]*f[21]*nu[148]; 
  out[38] += 1.677050983124842*rdxF[1]*f[38]*nu[158]+0.75*rdxF[1]*f[22]*nu[148]+5.031152949374527*rdxF[0]*f[38]*nu[60]+3.75*rdxF[0]*f[10]*nu[60]+3.354101966249685*rdxF[0]*f[18]*nu[50]+1.677050983124842*rdxF[0]*f[10]*nu[48]; 
  out[39] += 1.677050983124842*rdxF[1]*f[39]*nu[158]+0.75*rdxF[1]*f[23]*nu[148]; 
  out[40] += 1.677050983124842*rdxF[1]*f[40]*nu[158]+0.75*rdxF[1]*f[24]*nu[148]+1.677050983124842*rdxF[0]*f[40]*nu[60]+0.75*rdxF[0]*f[27]*nu[50]; 
  out[41] += 5.031152949374527*rdxF[1]*f[41]*nu[158]+3.75*rdxF[1]*f[5]*nu[158]+3.354101966249685*rdxF[1]*f[16]*nu[148]+1.677050983124842*rdxF[1]*f[5]*nu[144]+1.677050983124842*rdxF[0]*f[41]*nu[60]+0.75*rdxF[0]*f[28]*nu[50]; 
  out[42] += 5.031152949374527*rdxF[1]*f[42]*nu[158]+3.75*rdxF[1]*f[6]*nu[158]+3.354101966249685*rdxF[1]*f[17]*nu[148]+1.677050983124842*rdxF[1]*f[6]*nu[144]; 
  out[43] += 5.031152949374527*rdxF[1]*f[43]*nu[158]+3.75*rdxF[1]*f[7]*nu[158]+3.354101966249685*rdxF[1]*f[18]*nu[148]+1.677050983124842*rdxF[1]*f[7]*nu[144]+1.677050983124842*rdxF[0]*f[43]*nu[60]+0.75*rdxF[0]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[1]*f[44]*nu[158]+0.75*rdxF[1]*f[32]*nu[148]+1.677050983124842*rdxF[0]*f[44]*nu[60]+0.75*rdxF[0]*f[37]*nu[50]; 
  out[45] += 1.677050983124842*rdxF[1]*f[45]*nu[158]+0.75*rdxF[1]*f[33]*nu[148]+5.031152949374527*rdxF[0]*f[45]*nu[60]+3.75*rdxF[0]*f[17]*nu[60]+3.354101966249684*rdxF[0]*f[31]*nu[50]+1.677050983124842*rdxF[0]*f[17]*nu[48]; 
  out[46] += 1.677050983124842*rdxF[1]*f[46]*nu[158]+0.75*rdxF[1]*f[34]*nu[148]+1.677050983124842*rdxF[0]*f[46]*nu[60]+0.75*rdxF[0]*f[39]*nu[50]; 
  out[47] += 5.031152949374527*rdxF[1]*f[47]*nu[158]+3.75*rdxF[1]*f[15]*nu[158]+3.354101966249684*rdxF[1]*f[31]*nu[148]+1.677050983124842*rdxF[1]*f[15]*nu[144]+1.677050983124842*rdxF[0]*f[47]*nu[60]+0.75*rdxF[0]*f[42]*nu[50]; 

  return (rdxF[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+rdxF[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 1.677050983124842*rdxF[1]*f[2]*nu[60]+0.75*f[0]*rdxF[1]*nu[50]; 
  out[4] += 1.677050983124842*rdxF[2]*f[4]*nu[158]+0.75*f[0]*rdxF[2]*nu[148]; 
  out[5] += 1.677050983124842*rdxF[1]*f[5]*nu[60]+0.75*f[1]*rdxF[1]*nu[50]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 1.677050983124842*rdxF[1]*f[7]*nu[60]+0.75*rdxF[1]*f[3]*nu[50]; 
  out[8] += 1.677050983124842*rdxF[2]*f[8]*nu[158]+0.75*f[1]*rdxF[2]*nu[148]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 1.677050983124842*rdxF[2]*f[9]*nu[158]+0.75*f[2]*rdxF[2]*nu[148]+1.677050983124842*rdxF[1]*f[9]*nu[60]+0.75*rdxF[1]*f[4]*nu[50]; 
  out[10] += 1.677050983124842*rdxF[2]*f[10]*nu[158]+0.75*rdxF[2]*f[3]*nu[148]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 5.031152949374527*rdxF[1]*f[12]*nu[60]+3.75*f[0]*rdxF[1]*nu[60]+3.354101966249685*rdxF[1]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[1]*nu[48]; 
  out[14] += 5.031152949374527*rdxF[2]*f[14]*nu[158]+3.75*f[0]*rdxF[2]*nu[158]+3.354101966249685*rdxF[2]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[2]*nu[144]; 
  out[15] += 1.677050983124842*rdxF[1]*f[15]*nu[60]+0.75*rdxF[1]*f[6]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[2]*f[16]*nu[158]+0.75*rdxF[2]*f[5]*nu[148]+1.677050983124842*rdxF[1]*f[16]*nu[60]+0.75*rdxF[1]*f[8]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[2]*f[17]*nu[158]+0.75*rdxF[2]*f[6]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 1.677050983124842*rdxF[2]*f[18]*nu[158]+0.75*rdxF[2]*f[7]*nu[148]+1.677050983124842*rdxF[1]*f[18]*nu[60]+0.75*rdxF[1]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[1]*f[19]*nu[60]+0.75*rdxF[1]*f[11]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 5.031152949374527*rdxF[1]*f[20]*nu[60]+3.75*f[1]*rdxF[1]*nu[60]+3.354101966249684*rdxF[1]*f[5]*nu[50]+1.677050983124842*f[1]*rdxF[1]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 5.031152949374527*rdxF[1]*f[22]*nu[60]+3.75*rdxF[1]*f[3]*nu[60]+3.354101966249684*rdxF[1]*f[7]*nu[50]+1.677050983124842*rdxF[1]*f[3]*nu[48]; 
  out[23] += 1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 1.677050983124842*rdxF[1]*f[24]*nu[60]+0.75*rdxF[1]*f[13]*nu[50]; 
  out[25] += 1.677050983124842*rdxF[2]*f[25]*nu[158]+0.75*rdxF[2]*f[11]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.677050983124842*rdxF[2]*f[26]*nu[158]+0.75*rdxF[2]*f[12]*nu[148]+5.031152949374527*rdxF[1]*f[26]*nu[60]+3.75*rdxF[1]*f[4]*nu[60]+3.354101966249684*rdxF[1]*f[9]*nu[50]+1.677050983124842*rdxF[1]*f[4]*nu[48]; 
  out[27] += 1.677050983124842*rdxF[2]*f[27]*nu[158]+0.75*rdxF[2]*f[13]*nu[148]; 
  out[28] += 5.031152949374527*rdxF[2]*f[28]*nu[158]+3.75*f[1]*rdxF[2]*nu[158]+3.354101966249684*rdxF[2]*f[8]*nu[148]+1.677050983124842*f[1]*rdxF[2]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 5.031152949374527*rdxF[2]*f[29]*nu[158]+3.75*f[2]*rdxF[2]*nu[158]+3.354101966249684*rdxF[2]*f[9]*nu[148]+1.677050983124842*f[2]*rdxF[2]*nu[144]+1.677050983124842*rdxF[1]*f[29]*nu[60]+0.75*rdxF[1]*f[14]*nu[50]; 
  out[30] += 5.031152949374527*rdxF[2]*f[30]*nu[158]+3.75*rdxF[2]*f[3]*nu[158]+3.354101966249684*rdxF[2]*f[10]*nu[148]+1.677050983124842*rdxF[2]*f[3]*nu[144]; 
  out[31] += 1.677050983124842*rdxF[2]*f[31]*nu[158]+0.75*rdxF[2]*f[15]*nu[148]+1.677050983124842*rdxF[1]*f[31]*nu[60]+0.75*rdxF[1]*f[17]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 1.677050983124842*rdxF[1]*f[32]*nu[60]+0.75*rdxF[1]*f[21]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 5.031152949374527*rdxF[1]*f[33]*nu[60]+3.75*rdxF[1]*f[6]*nu[60]+3.354101966249685*rdxF[1]*f[15]*nu[50]+1.677050983124842*rdxF[1]*f[6]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 1.677050983124842*rdxF[1]*f[34]*nu[60]+0.75*rdxF[1]*f[23]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 1.677050983124842*rdxF[2]*f[35]*nu[158]+0.75*rdxF[2]*f[19]*nu[148]+1.677050983124842*rdxF[1]*f[35]*nu[60]+0.75*rdxF[1]*f[25]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 1.677050983124842*rdxF[2]*f[36]*nu[158]+0.75*rdxF[2]*f[20]*nu[148]+5.031152949374527*rdxF[1]*f[36]*nu[60]+3.75*rdxF[1]*f[8]*nu[60]+3.354101966249685*rdxF[1]*f[16]*nu[50]+1.677050983124842*rdxF[1]*f[8]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 1.677050983124842*rdxF[2]*f[37]*nu[158]+0.75*rdxF[2]*f[21]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[38] += 1.677050983124842*rdxF[2]*f[38]*nu[158]+0.75*rdxF[2]*f[22]*nu[148]+5.031152949374527*rdxF[1]*f[38]*nu[60]+3.75*rdxF[1]*f[10]*nu[60]+3.354101966249685*rdxF[1]*f[18]*nu[50]+1.677050983124842*rdxF[1]*f[10]*nu[48]; 
  out[39] += 1.677050983124842*rdxF[2]*f[39]*nu[158]+0.75*rdxF[2]*f[23]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[40] += 1.677050983124842*rdxF[2]*f[40]*nu[158]+0.75*rdxF[2]*f[24]*nu[148]+1.677050983124842*rdxF[1]*f[40]*nu[60]+0.75*rdxF[1]*f[27]*nu[50]; 
  out[41] += 5.031152949374527*rdxF[2]*f[41]*nu[158]+3.75*rdxF[2]*f[5]*nu[158]+3.354101966249685*rdxF[2]*f[16]*nu[148]+1.677050983124842*rdxF[2]*f[5]*nu[144]+1.677050983124842*rdxF[1]*f[41]*nu[60]+0.75*rdxF[1]*f[28]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 5.031152949374527*rdxF[2]*f[42]*nu[158]+3.75*rdxF[2]*f[6]*nu[158]+3.354101966249685*rdxF[2]*f[17]*nu[148]+1.677050983124842*rdxF[2]*f[6]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[43] += 5.031152949374527*rdxF[2]*f[43]*nu[158]+3.75*rdxF[2]*f[7]*nu[158]+3.354101966249685*rdxF[2]*f[18]*nu[148]+1.677050983124842*rdxF[2]*f[7]*nu[144]+1.677050983124842*rdxF[1]*f[43]*nu[60]+0.75*rdxF[1]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[2]*f[44]*nu[158]+0.75*rdxF[2]*f[32]*nu[148]+1.677050983124842*rdxF[1]*f[44]*nu[60]+0.75*rdxF[1]*f[37]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 1.677050983124842*rdxF[2]*f[45]*nu[158]+0.75*rdxF[2]*f[33]*nu[148]+5.031152949374527*rdxF[1]*f[45]*nu[60]+3.75*rdxF[1]*f[17]*nu[60]+3.354101966249684*rdxF[1]*f[31]*nu[50]+1.677050983124842*rdxF[1]*f[17]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 1.677050983124842*rdxF[2]*f[46]*nu[158]+0.75*rdxF[2]*f[34]*nu[148]+1.677050983124842*rdxF[1]*f[46]*nu[60]+0.75*rdxF[1]*f[39]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 5.031152949374527*rdxF[2]*f[47]*nu[158]+3.75*rdxF[2]*f[15]*nu[158]+3.354101966249684*rdxF[2]*f[31]*nu[148]+1.677050983124842*rdxF[2]*f[15]*nu[144]+1.677050983124842*rdxF[1]*f[47]*nu[60]+0.75*rdxF[1]*f[42]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+rdxF[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 1.677050983124842*rdxF[1]*f[2]*nu[60]+0.75*f[0]*rdxF[1]*nu[50]; 
  out[5] += 1.677050983124842*rdxF[1]*f[5]*nu[60]+0.75*f[1]*rdxF[1]*nu[50]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 1.677050983124842*rdxF[1]*f[7]*nu[60]+0.75*rdxF[1]*f[3]*nu[50]; 
  out[8] += 1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 1.677050983124842*rdxF[1]*f[9]*nu[60]+0.75*rdxF[1]*f[4]*nu[50]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 5.031152949374527*rdxF[1]*f[12]*nu[60]+3.75*f[0]*rdxF[1]*nu[60]+3.354101966249685*rdxF[1]*f[2]*nu[50]+1.677050983124842*f[0]*rdxF[1]*nu[48]; 
  out[15] += 1.677050983124842*rdxF[1]*f[15]*nu[60]+0.75*rdxF[1]*f[6]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[1]*f[16]*nu[60]+0.75*rdxF[1]*f[8]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 1.677050983124842*rdxF[1]*f[18]*nu[60]+0.75*rdxF[1]*f[10]*nu[50]; 
  out[19] += 1.677050983124842*rdxF[1]*f[19]*nu[60]+0.75*rdxF[1]*f[11]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 5.031152949374527*rdxF[1]*f[20]*nu[60]+3.75*f[1]*rdxF[1]*nu[60]+3.354101966249684*rdxF[1]*f[5]*nu[50]+1.677050983124842*f[1]*rdxF[1]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 5.031152949374527*rdxF[1]*f[22]*nu[60]+3.75*rdxF[1]*f[3]*nu[60]+3.354101966249684*rdxF[1]*f[7]*nu[50]+1.677050983124842*rdxF[1]*f[3]*nu[48]; 
  out[23] += 1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 1.677050983124842*rdxF[1]*f[24]*nu[60]+0.75*rdxF[1]*f[13]*nu[50]; 
  out[25] += 5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 5.031152949374527*rdxF[1]*f[26]*nu[60]+3.75*rdxF[1]*f[4]*nu[60]+3.354101966249684*rdxF[1]*f[9]*nu[50]+1.677050983124842*rdxF[1]*f[4]*nu[48]; 
  out[28] += 1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 1.677050983124842*rdxF[1]*f[29]*nu[60]+0.75*rdxF[1]*f[14]*nu[50]; 
  out[31] += 1.677050983124842*rdxF[1]*f[31]*nu[60]+0.75*rdxF[1]*f[17]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 1.677050983124842*rdxF[1]*f[32]*nu[60]+0.75*rdxF[1]*f[21]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 5.031152949374527*rdxF[1]*f[33]*nu[60]+3.75*rdxF[1]*f[6]*nu[60]+3.354101966249685*rdxF[1]*f[15]*nu[50]+1.677050983124842*rdxF[1]*f[6]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 1.677050983124842*rdxF[1]*f[34]*nu[60]+0.75*rdxF[1]*f[23]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 1.677050983124842*rdxF[1]*f[35]*nu[60]+0.75*rdxF[1]*f[25]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 5.031152949374527*rdxF[1]*f[36]*nu[60]+3.75*rdxF[1]*f[8]*nu[60]+3.354101966249685*rdxF[1]*f[16]*nu[50]+1.677050983124842*rdxF[1]*f[8]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[38] += 5.031152949374527*rdxF[1]*f[38]*nu[60]+3.75*rdxF[1]*f[10]*nu[60]+3.354101966249685*rdxF[1]*f[18]*nu[50]+1.677050983124842*rdxF[1]*f[10]*nu[48]; 
  out[39] += 1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[40] += 1.677050983124842*rdxF[1]*f[40]*nu[60]+0.75*rdxF[1]*f[27]*nu[50]; 
  out[41] += 1.677050983124842*rdxF[1]*f[41]*nu[60]+0.75*rdxF[1]*f[28]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[43] += 1.677050983124842*rdxF[1]*f[43]*nu[60]+0.75*rdxF[1]*f[30]*nu[50]; 
  out[44] += 1.677050983124842*rdxF[1]*f[44]*nu[60]+0.75*rdxF[1]*f[37]*nu[50]+5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 5.031152949374527*rdxF[1]*f[45]*nu[60]+3.75*rdxF[1]*f[17]*nu[60]+3.354101966249684*rdxF[1]*f[31]*nu[50]+1.677050983124842*rdxF[1]*f[17]*nu[48]+1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 1.677050983124842*rdxF[1]*f[46]*nu[60]+0.75*rdxF[1]*f[39]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 1.677050983124842*rdxF[1]*f[47]*nu[60]+0.75*rdxF[1]*f[42]*nu[50]+1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 

  out[3] += 1.677050983124842*rdxF[0]*f[3]*nu[109]+0.75*f[0]*rdxF[0]*nu[99]; 
  out[6] += 1.677050983124842*rdxF[0]*f[6]*nu[109]+0.75*rdxF[0]*f[1]*nu[99]; 
  out[7] += 1.677050983124842*rdxF[0]*f[7]*nu[109]+0.75*rdxF[0]*f[2]*nu[99]; 
  out[10] += 1.677050983124842*rdxF[0]*f[10]*nu[109]+0.75*rdxF[0]*f[4]*nu[99]; 
  out[13] += 5.031152949374527*rdxF[0]*f[13]*nu[109]+3.75*f[0]*rdxF[0]*nu[109]+3.354101966249685*rdxF[0]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[0]*nu[96]; 
  out[15] += 1.677050983124842*rdxF[0]*f[15]*nu[109]+0.75*rdxF[0]*f[5]*nu[99]; 
  out[17] += 1.677050983124842*rdxF[0]*f[17]*nu[109]+0.75*rdxF[0]*f[8]*nu[99]; 
  out[18] += 1.677050983124842*rdxF[0]*f[18]*nu[109]+0.75*rdxF[0]*f[9]*nu[99]; 
  out[21] += 1.677050983124842*rdxF[0]*f[21]*nu[109]+0.75*rdxF[0]*f[11]*nu[99]; 
  out[22] += 1.677050983124842*rdxF[0]*f[22]*nu[109]+0.75*rdxF[0]*f[12]*nu[99]; 
  out[23] += 5.031152949374527*rdxF[0]*f[23]*nu[109]+3.75*rdxF[0]*f[1]*nu[109]+3.354101966249684*rdxF[0]*f[6]*nu[99]+1.677050983124842*rdxF[0]*f[1]*nu[96]; 
  out[24] += 5.031152949374527*rdxF[0]*f[24]*nu[109]+3.75*rdxF[0]*f[2]*nu[109]+3.354101966249684*rdxF[0]*f[7]*nu[99]+1.677050983124842*rdxF[0]*f[2]*nu[96]; 
  out[27] += 5.031152949374527*rdxF[0]*f[27]*nu[109]+3.75*rdxF[0]*f[4]*nu[109]+3.354101966249684*rdxF[0]*f[10]*nu[99]+1.677050983124842*rdxF[0]*f[4]*nu[96]; 
  out[30] += 1.677050983124842*rdxF[0]*f[30]*nu[109]+0.75*rdxF[0]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[0]*f[31]*nu[109]+0.75*rdxF[0]*f[16]*nu[99]; 
  out[32] += 1.677050983124842*rdxF[0]*f[32]*nu[109]+0.75*rdxF[0]*f[19]*nu[99]; 
  out[33] += 1.677050983124842*rdxF[0]*f[33]*nu[109]+0.75*rdxF[0]*f[20]*nu[99]; 
  out[34] += 5.031152949374527*rdxF[0]*f[34]*nu[109]+3.75*rdxF[0]*f[5]*nu[109]+3.354101966249685*rdxF[0]*f[15]*nu[99]+1.677050983124842*rdxF[0]*f[5]*nu[96]; 
  out[37] += 1.677050983124842*rdxF[0]*f[37]*nu[109]+0.75*rdxF[0]*f[25]*nu[99]; 
  out[38] += 1.677050983124842*rdxF[0]*f[38]*nu[109]+0.75*rdxF[0]*f[26]*nu[99]; 
  out[39] += 5.031152949374527*rdxF[0]*f[39]*nu[109]+3.75*rdxF[0]*f[8]*nu[109]+3.354101966249685*rdxF[0]*f[17]*nu[99]+1.677050983124842*rdxF[0]*f[8]*nu[96]; 
  out[40] += 5.031152949374527*rdxF[0]*f[40]*nu[109]+3.75*rdxF[0]*f[9]*nu[109]+3.354101966249685*rdxF[0]*f[18]*nu[99]+1.677050983124842*rdxF[0]*f[9]*nu[96]; 
  out[42] += 1.677050983124842*rdxF[0]*f[42]*nu[109]+0.75*rdxF[0]*f[28]*nu[99]; 
  out[43] += 1.677050983124842*rdxF[0]*f[43]*nu[109]+0.75*rdxF[0]*f[29]*nu[99]; 
  out[44] += 1.677050983124842*rdxF[0]*f[44]*nu[109]+0.75*rdxF[0]*f[35]*nu[99]; 
  out[45] += 1.677050983124842*rdxF[0]*f[45]*nu[109]+0.75*rdxF[0]*f[36]*nu[99]; 
  out[46] += 5.031152949374527*rdxF[0]*f[46]*nu[109]+3.75*rdxF[0]*f[16]*nu[109]+3.354101966249684*rdxF[0]*f[31]*nu[99]+1.677050983124842*rdxF[0]*f[16]*nu[96]; 
  out[47] += 1.677050983124842*rdxF[0]*f[47]*nu[109]+0.75*rdxF[0]*f[41]*nu[99]; 

  return (rdxF[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[3] += 1.677050983124842*rdxF[0]*f[3]*nu[109]+0.75*f[0]*rdxF[0]*nu[99]; 
  out[4] += 1.677050983124842*rdxF[1]*f[4]*nu[158]+0.75*f[0]*rdxF[1]*nu[148]; 
  out[6] += 1.677050983124842*rdxF[0]*f[6]*nu[109]+0.75*rdxF[0]*f[1]*nu[99]; 
  out[7] += 1.677050983124842*rdxF[0]*f[7]*nu[109]+0.75*rdxF[0]*f[2]*nu[99]; 
  out[8] += 1.677050983124842*rdxF[1]*f[8]*nu[158]+0.75*f[1]*rdxF[1]*nu[148]; 
  out[9] += 1.677050983124842*rdxF[1]*f[9]*nu[158]+0.75*rdxF[1]*f[2]*nu[148]; 
  out[10] += 1.677050983124842*rdxF[1]*f[10]*nu[158]+0.75*rdxF[1]*f[3]*nu[148]+1.677050983124842*rdxF[0]*f[10]*nu[109]+0.75*rdxF[0]*f[4]*nu[99]; 
  out[13] += 5.031152949374527*rdxF[0]*f[13]*nu[109]+3.75*f[0]*rdxF[0]*nu[109]+3.354101966249685*rdxF[0]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[0]*nu[96]; 
  out[14] += 5.031152949374527*rdxF[1]*f[14]*nu[158]+3.75*f[0]*rdxF[1]*nu[158]+3.354101966249685*rdxF[1]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[1]*nu[144]; 
  out[15] += 1.677050983124842*rdxF[0]*f[15]*nu[109]+0.75*rdxF[0]*f[5]*nu[99]; 
  out[16] += 1.677050983124842*rdxF[1]*f[16]*nu[158]+0.75*rdxF[1]*f[5]*nu[148]; 
  out[17] += 1.677050983124842*rdxF[1]*f[17]*nu[158]+0.75*rdxF[1]*f[6]*nu[148]+1.677050983124842*rdxF[0]*f[17]*nu[109]+0.75*rdxF[0]*f[8]*nu[99]; 
  out[18] += 1.677050983124842*rdxF[1]*f[18]*nu[158]+0.75*rdxF[1]*f[7]*nu[148]+1.677050983124842*rdxF[0]*f[18]*nu[109]+0.75*rdxF[0]*f[9]*nu[99]; 
  out[21] += 1.677050983124842*rdxF[0]*f[21]*nu[109]+0.75*rdxF[0]*f[11]*nu[99]; 
  out[22] += 1.677050983124842*rdxF[0]*f[22]*nu[109]+0.75*rdxF[0]*f[12]*nu[99]; 
  out[23] += 5.031152949374527*rdxF[0]*f[23]*nu[109]+3.75*rdxF[0]*f[1]*nu[109]+3.354101966249684*rdxF[0]*f[6]*nu[99]+1.677050983124842*rdxF[0]*f[1]*nu[96]; 
  out[24] += 5.031152949374527*rdxF[0]*f[24]*nu[109]+3.75*rdxF[0]*f[2]*nu[109]+3.354101966249684*rdxF[0]*f[7]*nu[99]+1.677050983124842*rdxF[0]*f[2]*nu[96]; 
  out[25] += 1.677050983124842*rdxF[1]*f[25]*nu[158]+0.75*rdxF[1]*f[11]*nu[148]; 
  out[26] += 1.677050983124842*rdxF[1]*f[26]*nu[158]+0.75*rdxF[1]*f[12]*nu[148]; 
  out[27] += 1.677050983124842*rdxF[1]*f[27]*nu[158]+0.75*rdxF[1]*f[13]*nu[148]+5.031152949374527*rdxF[0]*f[27]*nu[109]+3.75*rdxF[0]*f[4]*nu[109]+3.354101966249684*rdxF[0]*f[10]*nu[99]+1.677050983124842*rdxF[0]*f[4]*nu[96]; 
  out[28] += 5.031152949374527*rdxF[1]*f[28]*nu[158]+3.75*f[1]*rdxF[1]*nu[158]+3.354101966249684*rdxF[1]*f[8]*nu[148]+1.677050983124842*f[1]*rdxF[1]*nu[144]; 
  out[29] += 5.031152949374527*rdxF[1]*f[29]*nu[158]+3.75*rdxF[1]*f[2]*nu[158]+3.354101966249684*rdxF[1]*f[9]*nu[148]+1.677050983124842*rdxF[1]*f[2]*nu[144]; 
  out[30] += 5.031152949374527*rdxF[1]*f[30]*nu[158]+3.75*rdxF[1]*f[3]*nu[158]+3.354101966249684*rdxF[1]*f[10]*nu[148]+1.677050983124842*rdxF[1]*f[3]*nu[144]+1.677050983124842*rdxF[0]*f[30]*nu[109]+0.75*rdxF[0]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[1]*f[31]*nu[158]+0.75*rdxF[1]*f[15]*nu[148]+1.677050983124842*rdxF[0]*f[31]*nu[109]+0.75*rdxF[0]*f[16]*nu[99]; 
  out[32] += 1.677050983124842*rdxF[0]*f[32]*nu[109]+0.75*rdxF[0]*f[19]*nu[99]; 
  out[33] += 1.677050983124842*rdxF[0]*f[33]*nu[109]+0.75*rdxF[0]*f[20]*nu[99]; 
  out[34] += 5.031152949374527*rdxF[0]*f[34]*nu[109]+3.75*rdxF[0]*f[5]*nu[109]+3.354101966249685*rdxF[0]*f[15]*nu[99]+1.677050983124842*rdxF[0]*f[5]*nu[96]; 
  out[35] += 1.677050983124842*rdxF[1]*f[35]*nu[158]+0.75*rdxF[1]*f[19]*nu[148]; 
  out[36] += 1.677050983124842*rdxF[1]*f[36]*nu[158]+0.75*rdxF[1]*f[20]*nu[148]; 
  out[37] += 1.677050983124842*rdxF[1]*f[37]*nu[158]+0.75*rdxF[1]*f[21]*nu[148]+1.677050983124842*rdxF[0]*f[37]*nu[109]+0.75*rdxF[0]*f[25]*nu[99]; 
  out[38] += 1.677050983124842*rdxF[1]*f[38]*nu[158]+0.75*rdxF[1]*f[22]*nu[148]+1.677050983124842*rdxF[0]*f[38]*nu[109]+0.75*rdxF[0]*f[26]*nu[99]; 
  out[39] += 1.677050983124842*rdxF[1]*f[39]*nu[158]+0.75*rdxF[1]*f[23]*nu[148]+5.031152949374527*rdxF[0]*f[39]*nu[109]+3.75*rdxF[0]*f[8]*nu[109]+3.354101966249685*rdxF[0]*f[17]*nu[99]+1.677050983124842*rdxF[0]*f[8]*nu[96]; 
  out[40] += 1.677050983124842*rdxF[1]*f[40]*nu[158]+0.75*rdxF[1]*f[24]*nu[148]+5.031152949374527*rdxF[0]*f[40]*nu[109]+3.75*rdxF[0]*f[9]*nu[109]+3.354101966249685*rdxF[0]*f[18]*nu[99]+1.677050983124842*rdxF[0]*f[9]*nu[96]; 
  out[41] += 5.031152949374527*rdxF[1]*f[41]*nu[158]+3.75*rdxF[1]*f[5]*nu[158]+3.354101966249685*rdxF[1]*f[16]*nu[148]+1.677050983124842*rdxF[1]*f[5]*nu[144]; 
  out[42] += 5.031152949374527*rdxF[1]*f[42]*nu[158]+3.75*rdxF[1]*f[6]*nu[158]+3.354101966249685*rdxF[1]*f[17]*nu[148]+1.677050983124842*rdxF[1]*f[6]*nu[144]+1.677050983124842*rdxF[0]*f[42]*nu[109]+0.75*rdxF[0]*f[28]*nu[99]; 
  out[43] += 5.031152949374527*rdxF[1]*f[43]*nu[158]+3.75*rdxF[1]*f[7]*nu[158]+3.354101966249685*rdxF[1]*f[18]*nu[148]+1.677050983124842*rdxF[1]*f[7]*nu[144]+1.677050983124842*rdxF[0]*f[43]*nu[109]+0.75*rdxF[0]*f[29]*nu[99]; 
  out[44] += 1.677050983124842*rdxF[1]*f[44]*nu[158]+0.75*rdxF[1]*f[32]*nu[148]+1.677050983124842*rdxF[0]*f[44]*nu[109]+0.75*rdxF[0]*f[35]*nu[99]; 
  out[45] += 1.677050983124842*rdxF[1]*f[45]*nu[158]+0.75*rdxF[1]*f[33]*nu[148]+1.677050983124842*rdxF[0]*f[45]*nu[109]+0.75*rdxF[0]*f[36]*nu[99]; 
  out[46] += 1.677050983124842*rdxF[1]*f[46]*nu[158]+0.75*rdxF[1]*f[34]*nu[148]+5.031152949374527*rdxF[0]*f[46]*nu[109]+3.75*rdxF[0]*f[16]*nu[109]+3.354101966249684*rdxF[0]*f[31]*nu[99]+1.677050983124842*rdxF[0]*f[16]*nu[96]; 
  out[47] += 5.031152949374527*rdxF[1]*f[47]*nu[158]+3.75*rdxF[1]*f[15]*nu[158]+3.354101966249684*rdxF[1]*f[31]*nu[148]+1.677050983124842*rdxF[1]*f[15]*nu[144]+1.677050983124842*rdxF[0]*f[47]*nu[109]+0.75*rdxF[0]*f[41]*nu[99]; 

  return (rdxF[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+rdxF[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[3] += 1.677050983124842*rdxF[1]*f[3]*nu[109]+0.75*f[0]*rdxF[1]*nu[99]; 
  out[4] += 1.677050983124842*rdxF[2]*f[4]*nu[158]+0.75*f[0]*rdxF[2]*nu[148]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[1]*f[6]*nu[109]+0.75*f[1]*rdxF[1]*nu[99]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 1.677050983124842*rdxF[1]*f[7]*nu[109]+0.75*rdxF[1]*f[2]*nu[99]; 
  out[8] += 1.677050983124842*rdxF[2]*f[8]*nu[158]+0.75*f[1]*rdxF[2]*nu[148]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 1.677050983124842*rdxF[2]*f[9]*nu[158]+0.75*f[2]*rdxF[2]*nu[148]; 
  out[10] += 1.677050983124842*rdxF[2]*f[10]*nu[158]+0.75*rdxF[2]*f[3]*nu[148]+1.677050983124842*rdxF[1]*f[10]*nu[109]+0.75*rdxF[1]*f[4]*nu[99]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[13] += 5.031152949374527*rdxF[1]*f[13]*nu[109]+3.75*f[0]*rdxF[1]*nu[109]+3.354101966249685*rdxF[1]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[1]*nu[96]; 
  out[14] += 5.031152949374527*rdxF[2]*f[14]*nu[158]+3.75*f[0]*rdxF[2]*nu[158]+3.354101966249685*rdxF[2]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[2]*nu[144]; 
  out[15] += 1.677050983124842*rdxF[1]*f[15]*nu[109]+0.75*rdxF[1]*f[5]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[2]*f[16]*nu[158]+0.75*rdxF[2]*f[5]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[2]*f[17]*nu[158]+0.75*rdxF[2]*f[6]*nu[148]+1.677050983124842*rdxF[1]*f[17]*nu[109]+0.75*rdxF[1]*f[8]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 1.677050983124842*rdxF[2]*f[18]*nu[158]+0.75*rdxF[2]*f[7]*nu[148]+1.677050983124842*rdxF[1]*f[18]*nu[109]+0.75*rdxF[1]*f[9]*nu[99]; 
  out[19] += 5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.677050983124842*rdxF[1]*f[21]*nu[109]+0.75*rdxF[1]*f[11]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.677050983124842*rdxF[1]*f[22]*nu[109]+0.75*rdxF[1]*f[12]*nu[99]; 
  out[23] += 5.031152949374527*rdxF[1]*f[23]*nu[109]+3.75*f[1]*rdxF[1]*nu[109]+3.354101966249684*rdxF[1]*f[6]*nu[99]+1.677050983124842*f[1]*rdxF[1]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 5.031152949374527*rdxF[1]*f[24]*nu[109]+3.75*rdxF[1]*f[2]*nu[109]+3.354101966249684*rdxF[1]*f[7]*nu[99]+1.677050983124842*rdxF[1]*f[2]*nu[96]; 
  out[25] += 1.677050983124842*rdxF[2]*f[25]*nu[158]+0.75*rdxF[2]*f[11]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.677050983124842*rdxF[2]*f[26]*nu[158]+0.75*rdxF[2]*f[12]*nu[148]; 
  out[27] += 1.677050983124842*rdxF[2]*f[27]*nu[158]+0.75*rdxF[2]*f[13]*nu[148]+5.031152949374527*rdxF[1]*f[27]*nu[109]+3.75*rdxF[1]*f[4]*nu[109]+3.354101966249684*rdxF[1]*f[10]*nu[99]+1.677050983124842*rdxF[1]*f[4]*nu[96]; 
  out[28] += 5.031152949374527*rdxF[2]*f[28]*nu[158]+3.75*f[1]*rdxF[2]*nu[158]+3.354101966249684*rdxF[2]*f[8]*nu[148]+1.677050983124842*f[1]*rdxF[2]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 5.031152949374527*rdxF[2]*f[29]*nu[158]+3.75*f[2]*rdxF[2]*nu[158]+3.354101966249684*rdxF[2]*f[9]*nu[148]+1.677050983124842*f[2]*rdxF[2]*nu[144]; 
  out[30] += 5.031152949374527*rdxF[2]*f[30]*nu[158]+3.75*rdxF[2]*f[3]*nu[158]+3.354101966249684*rdxF[2]*f[10]*nu[148]+1.677050983124842*rdxF[2]*f[3]*nu[144]+1.677050983124842*rdxF[1]*f[30]*nu[109]+0.75*rdxF[1]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[2]*f[31]*nu[158]+0.75*rdxF[2]*f[15]*nu[148]+1.677050983124842*rdxF[1]*f[31]*nu[109]+0.75*rdxF[1]*f[16]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 1.677050983124842*rdxF[1]*f[32]*nu[109]+0.75*rdxF[1]*f[19]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 1.677050983124842*rdxF[1]*f[33]*nu[109]+0.75*rdxF[1]*f[20]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 5.031152949374527*rdxF[1]*f[34]*nu[109]+3.75*rdxF[1]*f[5]*nu[109]+3.354101966249685*rdxF[1]*f[15]*nu[99]+1.677050983124842*rdxF[1]*f[5]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 1.677050983124842*rdxF[2]*f[35]*nu[158]+0.75*rdxF[2]*f[19]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 1.677050983124842*rdxF[2]*f[36]*nu[158]+0.75*rdxF[2]*f[20]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 1.677050983124842*rdxF[2]*f[37]*nu[158]+0.75*rdxF[2]*f[21]*nu[148]+1.677050983124842*rdxF[1]*f[37]*nu[109]+0.75*rdxF[1]*f[25]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[38] += 1.677050983124842*rdxF[2]*f[38]*nu[158]+0.75*rdxF[2]*f[22]*nu[148]+1.677050983124842*rdxF[1]*f[38]*nu[109]+0.75*rdxF[1]*f[26]*nu[99]; 
  out[39] += 1.677050983124842*rdxF[2]*f[39]*nu[158]+0.75*rdxF[2]*f[23]*nu[148]+5.031152949374527*rdxF[1]*f[39]*nu[109]+3.75*rdxF[1]*f[8]*nu[109]+3.354101966249685*rdxF[1]*f[17]*nu[99]+1.677050983124842*rdxF[1]*f[8]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[40] += 1.677050983124842*rdxF[2]*f[40]*nu[158]+0.75*rdxF[2]*f[24]*nu[148]+5.031152949374527*rdxF[1]*f[40]*nu[109]+3.75*rdxF[1]*f[9]*nu[109]+3.354101966249685*rdxF[1]*f[18]*nu[99]+1.677050983124842*rdxF[1]*f[9]*nu[96]; 
  out[41] += 5.031152949374527*rdxF[2]*f[41]*nu[158]+3.75*rdxF[2]*f[5]*nu[158]+3.354101966249685*rdxF[2]*f[16]*nu[148]+1.677050983124842*rdxF[2]*f[5]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 5.031152949374527*rdxF[2]*f[42]*nu[158]+3.75*rdxF[2]*f[6]*nu[158]+3.354101966249685*rdxF[2]*f[17]*nu[148]+1.677050983124842*rdxF[2]*f[6]*nu[144]+1.677050983124842*rdxF[1]*f[42]*nu[109]+0.75*rdxF[1]*f[28]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[43] += 5.031152949374527*rdxF[2]*f[43]*nu[158]+3.75*rdxF[2]*f[7]*nu[158]+3.354101966249685*rdxF[2]*f[18]*nu[148]+1.677050983124842*rdxF[2]*f[7]*nu[144]+1.677050983124842*rdxF[1]*f[43]*nu[109]+0.75*rdxF[1]*f[29]*nu[99]; 
  out[44] += 1.677050983124842*rdxF[2]*f[44]*nu[158]+0.75*rdxF[2]*f[32]*nu[148]+1.677050983124842*rdxF[1]*f[44]*nu[109]+0.75*rdxF[1]*f[35]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 1.677050983124842*rdxF[2]*f[45]*nu[158]+0.75*rdxF[2]*f[33]*nu[148]+1.677050983124842*rdxF[1]*f[45]*nu[109]+0.75*rdxF[1]*f[36]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 1.677050983124842*rdxF[2]*f[46]*nu[158]+0.75*rdxF[2]*f[34]*nu[148]+5.031152949374527*rdxF[1]*f[46]*nu[109]+3.75*rdxF[1]*f[16]*nu[109]+3.354101966249684*rdxF[1]*f[31]*nu[99]+1.677050983124842*rdxF[1]*f[16]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 5.031152949374527*rdxF[2]*f[47]*nu[158]+3.75*rdxF[2]*f[15]*nu[158]+3.354101966249684*rdxF[2]*f[31]*nu[148]+1.677050983124842*rdxF[2]*f[15]*nu[144]+1.677050983124842*rdxF[1]*f[47]*nu[109]+0.75*rdxF[1]*f[41]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+rdxF[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[3] += 1.677050983124842*rdxF[1]*f[3]*nu[109]+0.75*f[0]*rdxF[1]*nu[99]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[1]*f[6]*nu[109]+0.75*f[1]*rdxF[1]*nu[99]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 1.677050983124842*rdxF[1]*f[7]*nu[109]+0.75*rdxF[1]*f[2]*nu[99]; 
  out[8] += 1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[10] += 1.677050983124842*rdxF[1]*f[10]*nu[109]+0.75*rdxF[1]*f[4]*nu[99]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[13] += 5.031152949374527*rdxF[1]*f[13]*nu[109]+3.75*f[0]*rdxF[1]*nu[109]+3.354101966249685*rdxF[1]*f[3]*nu[99]+1.677050983124842*f[0]*rdxF[1]*nu[96]; 
  out[15] += 1.677050983124842*rdxF[1]*f[15]*nu[109]+0.75*rdxF[1]*f[5]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[1]*f[17]*nu[109]+0.75*rdxF[1]*f[8]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 1.677050983124842*rdxF[1]*f[18]*nu[109]+0.75*rdxF[1]*f[9]*nu[99]; 
  out[19] += 5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.677050983124842*rdxF[1]*f[21]*nu[109]+0.75*rdxF[1]*f[11]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.677050983124842*rdxF[1]*f[22]*nu[109]+0.75*rdxF[1]*f[12]*nu[99]; 
  out[23] += 5.031152949374527*rdxF[1]*f[23]*nu[109]+3.75*f[1]*rdxF[1]*nu[109]+3.354101966249684*rdxF[1]*f[6]*nu[99]+1.677050983124842*f[1]*rdxF[1]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 5.031152949374527*rdxF[1]*f[24]*nu[109]+3.75*rdxF[1]*f[2]*nu[109]+3.354101966249684*rdxF[1]*f[7]*nu[99]+1.677050983124842*rdxF[1]*f[2]*nu[96]; 
  out[25] += 5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[27] += 5.031152949374527*rdxF[1]*f[27]*nu[109]+3.75*rdxF[1]*f[4]*nu[109]+3.354101966249684*rdxF[1]*f[10]*nu[99]+1.677050983124842*rdxF[1]*f[4]*nu[96]; 
  out[28] += 1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[30] += 1.677050983124842*rdxF[1]*f[30]*nu[109]+0.75*rdxF[1]*f[14]*nu[99]; 
  out[31] += 1.677050983124842*rdxF[1]*f[31]*nu[109]+0.75*rdxF[1]*f[16]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 1.677050983124842*rdxF[1]*f[32]*nu[109]+0.75*rdxF[1]*f[19]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 1.677050983124842*rdxF[1]*f[33]*nu[109]+0.75*rdxF[1]*f[20]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 5.031152949374527*rdxF[1]*f[34]*nu[109]+3.75*rdxF[1]*f[5]*nu[109]+3.354101966249685*rdxF[1]*f[15]*nu[99]+1.677050983124842*rdxF[1]*f[5]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 1.677050983124842*rdxF[1]*f[37]*nu[109]+0.75*rdxF[1]*f[25]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[38] += 1.677050983124842*rdxF[1]*f[38]*nu[109]+0.75*rdxF[1]*f[26]*nu[99]; 
  out[39] += 5.031152949374527*rdxF[1]*f[39]*nu[109]+3.75*rdxF[1]*f[8]*nu[109]+3.354101966249685*rdxF[1]*f[17]*nu[99]+1.677050983124842*rdxF[1]*f[8]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[40] += 5.031152949374527*rdxF[1]*f[40]*nu[109]+3.75*rdxF[1]*f[9]*nu[109]+3.354101966249685*rdxF[1]*f[18]*nu[99]+1.677050983124842*rdxF[1]*f[9]*nu[96]; 
  out[41] += 1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 1.677050983124842*rdxF[1]*f[42]*nu[109]+0.75*rdxF[1]*f[28]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[43] += 1.677050983124842*rdxF[1]*f[43]*nu[109]+0.75*rdxF[1]*f[29]*nu[99]; 
  out[44] += 1.677050983124842*rdxF[1]*f[44]*nu[109]+0.75*rdxF[1]*f[35]*nu[99]+5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 1.677050983124842*rdxF[1]*f[45]*nu[109]+0.75*rdxF[1]*f[36]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 5.031152949374527*rdxF[1]*f[46]*nu[109]+3.75*rdxF[1]*f[16]*nu[109]+3.354101966249684*rdxF[1]*f[31]*nu[99]+1.677050983124842*rdxF[1]*f[16]*nu[96]+1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 1.677050983124842*rdxF[1]*f[47]*nu[109]+0.75*rdxF[1]*f[41]*nu[99]+1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[3]*dx[3]); 

  out[4] += 1.677050983124842*rdxF[0]*f[4]*nu[158]+0.75*f[0]*rdxF[0]*nu[148]; 
  out[8] += 1.677050983124842*rdxF[0]*f[8]*nu[158]+0.75*rdxF[0]*f[1]*nu[148]; 
  out[9] += 1.677050983124842*rdxF[0]*f[9]*nu[158]+0.75*rdxF[0]*f[2]*nu[148]; 
  out[10] += 1.677050983124842*rdxF[0]*f[10]*nu[158]+0.75*rdxF[0]*f[3]*nu[148]; 
  out[14] += 5.031152949374527*rdxF[0]*f[14]*nu[158]+3.75*f[0]*rdxF[0]*nu[158]+3.354101966249685*rdxF[0]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[0]*nu[144]; 
  out[16] += 1.677050983124842*rdxF[0]*f[16]*nu[158]+0.75*rdxF[0]*f[5]*nu[148]; 
  out[17] += 1.677050983124842*rdxF[0]*f[17]*nu[158]+0.75*rdxF[0]*f[6]*nu[148]; 
  out[18] += 1.677050983124842*rdxF[0]*f[18]*nu[158]+0.75*rdxF[0]*f[7]*nu[148]; 
  out[25] += 1.677050983124842*rdxF[0]*f[25]*nu[158]+0.75*rdxF[0]*f[11]*nu[148]; 
  out[26] += 1.677050983124842*rdxF[0]*f[26]*nu[158]+0.75*rdxF[0]*f[12]*nu[148]; 
  out[27] += 1.677050983124842*rdxF[0]*f[27]*nu[158]+0.75*rdxF[0]*f[13]*nu[148]; 
  out[28] += 5.031152949374527*rdxF[0]*f[28]*nu[158]+3.75*rdxF[0]*f[1]*nu[158]+3.354101966249684*rdxF[0]*f[8]*nu[148]+1.677050983124842*rdxF[0]*f[1]*nu[144]; 
  out[29] += 5.031152949374527*rdxF[0]*f[29]*nu[158]+3.75*rdxF[0]*f[2]*nu[158]+3.354101966249684*rdxF[0]*f[9]*nu[148]+1.677050983124842*rdxF[0]*f[2]*nu[144]; 
  out[30] += 5.031152949374527*rdxF[0]*f[30]*nu[158]+3.75*rdxF[0]*f[3]*nu[158]+3.354101966249684*rdxF[0]*f[10]*nu[148]+1.677050983124842*rdxF[0]*f[3]*nu[144]; 
  out[31] += 1.677050983124842*rdxF[0]*f[31]*nu[158]+0.75*rdxF[0]*f[15]*nu[148]; 
  out[35] += 1.677050983124842*rdxF[0]*f[35]*nu[158]+0.75*rdxF[0]*f[19]*nu[148]; 
  out[36] += 1.677050983124842*rdxF[0]*f[36]*nu[158]+0.75*rdxF[0]*f[20]*nu[148]; 
  out[37] += 1.677050983124842*rdxF[0]*f[37]*nu[158]+0.75*rdxF[0]*f[21]*nu[148]; 
  out[38] += 1.677050983124842*rdxF[0]*f[38]*nu[158]+0.75*rdxF[0]*f[22]*nu[148]; 
  out[39] += 1.677050983124842*rdxF[0]*f[39]*nu[158]+0.75*rdxF[0]*f[23]*nu[148]; 
  out[40] += 1.677050983124842*rdxF[0]*f[40]*nu[158]+0.75*rdxF[0]*f[24]*nu[148]; 
  out[41] += 5.031152949374527*rdxF[0]*f[41]*nu[158]+3.75*rdxF[0]*f[5]*nu[158]+3.354101966249685*rdxF[0]*f[16]*nu[148]+1.677050983124842*rdxF[0]*f[5]*nu[144]; 
  out[42] += 5.031152949374527*rdxF[0]*f[42]*nu[158]+3.75*rdxF[0]*f[6]*nu[158]+3.354101966249685*rdxF[0]*f[17]*nu[148]+1.677050983124842*rdxF[0]*f[6]*nu[144]; 
  out[43] += 5.031152949374527*rdxF[0]*f[43]*nu[158]+3.75*rdxF[0]*f[7]*nu[158]+3.354101966249685*rdxF[0]*f[18]*nu[148]+1.677050983124842*rdxF[0]*f[7]*nu[144]; 
  out[44] += 1.677050983124842*rdxF[0]*f[44]*nu[158]+0.75*rdxF[0]*f[32]*nu[148]; 
  out[45] += 1.677050983124842*rdxF[0]*f[45]*nu[158]+0.75*rdxF[0]*f[33]*nu[148]; 
  out[46] += 1.677050983124842*rdxF[0]*f[46]*nu[158]+0.75*rdxF[0]*f[34]*nu[148]; 
  out[47] += 5.031152949374527*rdxF[0]*f[47]*nu[158]+3.75*rdxF[0]*f[15]*nu[158]+3.354101966249684*rdxF[0]*f[31]*nu[148]+1.677050983124842*rdxF[0]*f[15]*nu[144]; 

  return (rdxF[0]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[4] += 1.677050983124842*rdxF[1]*f[4]*nu[158]+0.75*f[0]*rdxF[1]*nu[148]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[8] += 1.677050983124842*rdxF[1]*f[8]*nu[158]+0.75*f[1]*rdxF[1]*nu[148]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 1.677050983124842*rdxF[1]*f[9]*nu[158]+0.75*rdxF[1]*f[2]*nu[148]; 
  out[10] += 1.677050983124842*rdxF[1]*f[10]*nu[158]+0.75*rdxF[1]*f[3]*nu[148]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[14] += 5.031152949374527*rdxF[1]*f[14]*nu[158]+3.75*f[0]*rdxF[1]*nu[158]+3.354101966249685*rdxF[1]*f[4]*nu[148]+1.677050983124842*f[0]*rdxF[1]*nu[144]; 
  out[15] += 1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[1]*f[16]*nu[158]+0.75*rdxF[1]*f[5]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[1]*f[17]*nu[158]+0.75*rdxF[1]*f[6]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 1.677050983124842*rdxF[1]*f[18]*nu[158]+0.75*rdxF[1]*f[7]*nu[148]; 
  out[19] += 5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[23] += 1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[25] += 1.677050983124842*rdxF[1]*f[25]*nu[158]+0.75*rdxF[1]*f[11]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.677050983124842*rdxF[1]*f[26]*nu[158]+0.75*rdxF[1]*f[12]*nu[148]; 
  out[27] += 1.677050983124842*rdxF[1]*f[27]*nu[158]+0.75*rdxF[1]*f[13]*nu[148]; 
  out[28] += 5.031152949374527*rdxF[1]*f[28]*nu[158]+3.75*f[1]*rdxF[1]*nu[158]+3.354101966249684*rdxF[1]*f[8]*nu[148]+1.677050983124842*f[1]*rdxF[1]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 5.031152949374527*rdxF[1]*f[29]*nu[158]+3.75*rdxF[1]*f[2]*nu[158]+3.354101966249684*rdxF[1]*f[9]*nu[148]+1.677050983124842*rdxF[1]*f[2]*nu[144]; 
  out[30] += 5.031152949374527*rdxF[1]*f[30]*nu[158]+3.75*rdxF[1]*f[3]*nu[158]+3.354101966249684*rdxF[1]*f[10]*nu[148]+1.677050983124842*rdxF[1]*f[3]*nu[144]; 
  out[31] += 1.677050983124842*rdxF[1]*f[31]*nu[158]+0.75*rdxF[1]*f[15]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 1.677050983124842*rdxF[1]*f[35]*nu[158]+0.75*rdxF[1]*f[19]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 1.677050983124842*rdxF[1]*f[36]*nu[158]+0.75*rdxF[1]*f[20]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 1.677050983124842*rdxF[1]*f[37]*nu[158]+0.75*rdxF[1]*f[21]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[38] += 1.677050983124842*rdxF[1]*f[38]*nu[158]+0.75*rdxF[1]*f[22]*nu[148]; 
  out[39] += 1.677050983124842*rdxF[1]*f[39]*nu[158]+0.75*rdxF[1]*f[23]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[40] += 1.677050983124842*rdxF[1]*f[40]*nu[158]+0.75*rdxF[1]*f[24]*nu[148]; 
  out[41] += 5.031152949374527*rdxF[1]*f[41]*nu[158]+3.75*rdxF[1]*f[5]*nu[158]+3.354101966249685*rdxF[1]*f[16]*nu[148]+1.677050983124842*rdxF[1]*f[5]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 5.031152949374527*rdxF[1]*f[42]*nu[158]+3.75*rdxF[1]*f[6]*nu[158]+3.354101966249685*rdxF[1]*f[17]*nu[148]+1.677050983124842*rdxF[1]*f[6]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[43] += 5.031152949374527*rdxF[1]*f[43]*nu[158]+3.75*rdxF[1]*f[7]*nu[158]+3.354101966249685*rdxF[1]*f[18]*nu[148]+1.677050983124842*rdxF[1]*f[7]*nu[144]; 
  out[44] += 1.677050983124842*rdxF[1]*f[44]*nu[158]+0.75*rdxF[1]*f[32]*nu[148]+5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 1.677050983124842*rdxF[1]*f[45]*nu[158]+0.75*rdxF[1]*f[33]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 1.677050983124842*rdxF[1]*f[46]*nu[158]+0.75*rdxF[1]*f[34]*nu[148]+1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 5.031152949374527*rdxF[1]*f[47]*nu[158]+3.75*rdxF[1]*f[15]*nu[158]+3.354101966249684*rdxF[1]*f[31]*nu[148]+1.677050983124842*rdxF[1]*f[15]*nu[144]+1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*0.9;

} 
double ConstDiffusionVarCoeffVol4xSerP2_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[192]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[5] += 1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[8] += 1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[11] += 5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[15] += 1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[19] += 5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[23] += 1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[25] += 5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[28] += 1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[31] += 1.677050983124842*rdxF[0]*nu[11]*f[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[32] += 5.031152949374527*rdxF[0]*nu[11]*f[32]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[33] += 1.677050983124842*rdxF[0]*nu[11]*f[33]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[34] += 1.677050983124842*rdxF[0]*nu[11]*f[34]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[35] += 5.031152949374527*rdxF[0]*nu[11]*f[35]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[36] += 1.677050983124842*rdxF[0]*nu[11]*f[36]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[37] += 5.031152949374527*rdxF[0]*nu[11]*f[37]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[39] += 1.677050983124842*rdxF[0]*nu[11]*f[39]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[41] += 1.677050983124842*rdxF[0]*nu[11]*f[41]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[42] += 1.677050983124842*rdxF[0]*nu[11]*f[42]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[44] += 5.031152949374527*rdxF[0]*nu[11]*f[44]+3.354101966249684*rdxF[0]*nu[1]*f[31]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[45] += 1.677050983124842*rdxF[0]*nu[11]*f[45]+0.75*rdxF[0]*nu[1]*f[38]; 
  out[46] += 1.677050983124842*rdxF[0]*nu[11]*f[46]+0.75*rdxF[0]*nu[1]*f[40]; 
  out[47] += 1.677050983124842*rdxF[0]*nu[11]*f[47]+0.75*rdxF[0]*nu[1]*f[43]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*0.9;

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[3]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstDiffusionCFLfreqMin4xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[3]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[3]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[62])-0.2795084971874737*nu[61]-0.2795084971874737*nu[60]-0.2795084971874737*nu[59]+0.25*nu[48]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96])+kxSq[2]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[110])-0.2795084971874737*nu[109]-0.2795084971874737*nu[108]-0.2795084971874737*nu[107]+0.25*nu[96]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[158])-0.2795084971874737*nu[157]-0.2795084971874737*nu[156]-0.2795084971874737*nu[155]+0.25*nu[144]))*1.8);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP2_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[192]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*1.8);

} 
