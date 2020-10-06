#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol4xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[32] += 22.9128784747792*rdxFnu[0]*f[2]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[40] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[42] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[49] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[51] += 22.9128784747792*rdxFnu[0]*f[7]; 
  out[55] += 22.9128784747792*rdxFnu[0]*f[9]; 
  out[61] += 6.708203932499369*rdxFnu[0]*f[17]; 
  out[65] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[68] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[70] += 22.9128784747792*rdxFnu[0]*f[18]; 
  out[77] += 22.9128784747792*rdxFnu[0]*f[35]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[32] += 22.9128784747792*rdxFnu[0]*f[2]; 
  out[33] += 22.9128784747792*rdxFnu[1]*f[3]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[40] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[42] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[43] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[49] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[51] += 22.9128784747792*rdxFnu[0]*f[7]; 
  out[52] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[55] += 22.9128784747792*rdxFnu[0]*f[9]; 
  out[56] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[61] += 6.708203932499369*rdxFnu[0]*f[17]; 
  out[62] += 6.708203932499369*rdxFnu[1]*f[16]; 
  out[65] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[66] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[68] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[70] += 22.9128784747792*rdxFnu[0]*f[18]; 
  out[71] += 22.9128784747792*rdxFnu[1]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[77] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[78] += 22.9128784747792*rdxFnu[1]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[32] += 22.9128784747792*rdxFnu[0]*f[2]; 
  out[33] += 22.9128784747792*rdxFnu[1]*f[3]; 
  out[34] += 22.9128784747792*rdxFnu[2]*f[4]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[40] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[42] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[43] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[45] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[2]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[2]*f[7]; 
  out[49] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[51] += 22.9128784747792*rdxFnu[0]*f[7]; 
  out[52] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[55] += 22.9128784747792*rdxFnu[0]*f[9]; 
  out[56] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[57] += 22.9128784747792*rdxFnu[2]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[2]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[2]*f[10]; 
  out[61] += 6.708203932499369*rdxFnu[0]*f[17]; 
  out[62] += 6.708203932499369*rdxFnu[1]*f[16]; 
  out[63] += 6.708203932499369*rdxFnu[2]*f[15]; 
  out[65] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[66] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[68] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[70] += 22.9128784747792*rdxFnu[0]*f[18]; 
  out[71] += 22.9128784747792*rdxFnu[1]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[73] += 22.9128784747792*rdxFnu[2]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[2]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[2]*f[18]; 
  out[77] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[78] += 22.9128784747792*rdxFnu[1]*f[35]; 
  out[79] += 22.9128784747792*rdxFnu[2]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[32] += 22.9128784747792*rdxFnu[1]*f[2]; 
  out[33] += 22.9128784747792*rdxFnu[2]*f[3]; 
  out[34] += 22.9128784747792*rdxFnu[3]*f[4]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[37] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[38] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[42] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[43] += 6.708203932499369*rdxFnu[2]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[2]*f[9]; 
  out[45] += 6.708203932499369*rdxFnu[3]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[3]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[3]*f[7]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[49] += 22.9128784747792*rdxFnu[1]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[51] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[52] += 22.9128784747792*rdxFnu[2]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[2]*f[7]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[55] += 22.9128784747792*rdxFnu[1]*f[9]; 
  out[56] += 22.9128784747792*rdxFnu[2]*f[10]; 
  out[57] += 22.9128784747792*rdxFnu[3]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[3]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[3]*f[10]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[61] += 6.708203932499369*rdxFnu[1]*f[17]; 
  out[62] += 6.708203932499369*rdxFnu[2]*f[16]; 
  out[63] += 6.708203932499369*rdxFnu[3]*f[15]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[65] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[66] += 22.9128784747792*rdxFnu[2]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[68] += 22.9128784747792*rdxFnu[1]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[70] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[71] += 22.9128784747792*rdxFnu[2]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[2]*f[18]; 
  out[73] += 22.9128784747792*rdxFnu[3]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[3]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[3]*f[18]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[77] += 22.9128784747792*rdxFnu[1]*f[35]; 
  out[78] += 22.9128784747792*rdxFnu[2]*f[35]; 
  out[79] += 22.9128784747792*rdxFnu[3]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[32] += 22.9128784747792*rdxFnu[1]*f[2]; 
  out[33] += 22.9128784747792*rdxFnu[2]*f[3]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[37] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[38] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[42] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[43] += 6.708203932499369*rdxFnu[2]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[2]*f[9]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[49] += 22.9128784747792*rdxFnu[1]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[51] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[52] += 22.9128784747792*rdxFnu[2]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[2]*f[7]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[55] += 22.9128784747792*rdxFnu[1]*f[9]; 
  out[56] += 22.9128784747792*rdxFnu[2]*f[10]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[61] += 6.708203932499369*rdxFnu[1]*f[17]; 
  out[62] += 6.708203932499369*rdxFnu[2]*f[16]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[65] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[66] += 22.9128784747792*rdxFnu[2]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[68] += 22.9128784747792*rdxFnu[1]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[70] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[71] += 22.9128784747792*rdxFnu[2]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[2]*f[18]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[77] += 22.9128784747792*rdxFnu[1]*f[35]; 
  out[78] += 22.9128784747792*rdxFnu[2]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[32] += 22.9128784747792*rdxFnu[0]*f[2]; 
  out[34] += 22.9128784747792*rdxFnu[1]*f[4]; 
  out[37] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[40] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[42] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[45] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[1]*f[7]; 
  out[49] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[51] += 22.9128784747792*rdxFnu[0]*f[7]; 
  out[55] += 22.9128784747792*rdxFnu[0]*f[9]; 
  out[57] += 22.9128784747792*rdxFnu[1]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[1]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[61] += 6.708203932499369*rdxFnu[0]*f[17]; 
  out[63] += 6.708203932499369*rdxFnu[1]*f[15]; 
  out[65] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[68] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[70] += 22.9128784747792*rdxFnu[0]*f[18]; 
  out[73] += 22.9128784747792*rdxFnu[1]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[1]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[77] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[79] += 22.9128784747792*rdxFnu[1]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[32] += 22.9128784747792*rdxFnu[1]*f[2]; 
  out[34] += 22.9128784747792*rdxFnu[2]*f[4]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[37] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[42] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[45] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[2]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[2]*f[7]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[49] += 22.9128784747792*rdxFnu[1]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[51] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[55] += 22.9128784747792*rdxFnu[1]*f[9]; 
  out[57] += 22.9128784747792*rdxFnu[2]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[2]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[2]*f[10]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[61] += 6.708203932499369*rdxFnu[1]*f[17]; 
  out[63] += 6.708203932499369*rdxFnu[2]*f[15]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[65] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[68] += 22.9128784747792*rdxFnu[1]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[70] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[73] += 22.9128784747792*rdxFnu[2]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[2]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[2]*f[18]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[77] += 22.9128784747792*rdxFnu[1]*f[35]; 
  out[79] += 22.9128784747792*rdxFnu[2]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[32] += 22.9128784747792*rdxFnu[1]*f[2]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[37] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[40] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[42] += 6.708203932499369*rdxFnu[1]*f[10]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[49] += 22.9128784747792*rdxFnu[1]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[51] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[55] += 22.9128784747792*rdxFnu[1]*f[9]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[61] += 6.708203932499369*rdxFnu[1]*f[17]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[65] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[68] += 22.9128784747792*rdxFnu[1]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[70] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[77] += 22.9128784747792*rdxFnu[1]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[33] += 22.9128784747792*rdxFnu[0]*f[3]; 
  out[38] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[43] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[52] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[0]*f[7]; 
  out[56] += 22.9128784747792*rdxFnu[0]*f[10]; 
  out[62] += 6.708203932499369*rdxFnu[0]*f[16]; 
  out[66] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[71] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[0]*f[18]; 
  out[78] += 22.9128784747792*rdxFnu[0]*f[35]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[33] += 22.9128784747792*rdxFnu[0]*f[3]; 
  out[34] += 22.9128784747792*rdxFnu[1]*f[4]; 
  out[38] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[43] += 6.708203932499369*rdxFnu[0]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[45] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[1]*f[7]; 
  out[52] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[0]*f[7]; 
  out[56] += 22.9128784747792*rdxFnu[0]*f[10]; 
  out[57] += 22.9128784747792*rdxFnu[1]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[1]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[62] += 6.708203932499369*rdxFnu[0]*f[16]; 
  out[63] += 6.708203932499369*rdxFnu[1]*f[15]; 
  out[66] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[71] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[0]*f[18]; 
  out[73] += 22.9128784747792*rdxFnu[1]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[1]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[78] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[79] += 22.9128784747792*rdxFnu[1]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[33] += 22.9128784747792*rdxFnu[1]*f[3]; 
  out[34] += 22.9128784747792*rdxFnu[2]*f[4]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[43] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[45] += 6.708203932499369*rdxFnu[2]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[2]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[2]*f[7]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[52] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[56] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[57] += 22.9128784747792*rdxFnu[2]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[2]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[2]*f[10]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[62] += 6.708203932499369*rdxFnu[1]*f[16]; 
  out[63] += 6.708203932499369*rdxFnu[2]*f[15]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[66] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[71] += 22.9128784747792*rdxFnu[1]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[73] += 22.9128784747792*rdxFnu[2]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[2]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[2]*f[18]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[78] += 22.9128784747792*rdxFnu[1]*f[35]; 
  out[79] += 22.9128784747792*rdxFnu[2]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[33] += 22.9128784747792*rdxFnu[1]*f[3]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[38] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[43] += 6.708203932499369*rdxFnu[1]*f[8]; 
  out[44] += 6.708203932499369*rdxFnu[1]*f[9]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[52] += 22.9128784747792*rdxFnu[1]*f[6]; 
  out[53] += 22.9128784747792*rdxFnu[1]*f[7]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[56] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[62] += 6.708203932499369*rdxFnu[1]*f[16]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[66] += 22.9128784747792*rdxFnu[1]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[71] += 22.9128784747792*rdxFnu[1]*f[17]; 
  out[72] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[78] += 22.9128784747792*rdxFnu[1]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[34] += 22.9128784747792*rdxFnu[0]*f[4]; 
  out[45] += 6.708203932499369*rdxFnu[0]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[0]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[57] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[0]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[0]*f[10]; 
  out[63] += 6.708203932499369*rdxFnu[0]*f[15]; 
  out[73] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[0]*f[18]; 
  out[79] += 22.9128784747792*rdxFnu[0]*f[35]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[34] += 22.9128784747792*rdxFnu[1]*f[4]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[45] += 6.708203932499369*rdxFnu[1]*f[5]; 
  out[46] += 6.708203932499369*rdxFnu[1]*f[6]; 
  out[47] += 6.708203932499369*rdxFnu[1]*f[7]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[57] += 22.9128784747792*rdxFnu[1]*f[8]; 
  out[58] += 22.9128784747792*rdxFnu[1]*f[9]; 
  out[59] += 22.9128784747792*rdxFnu[1]*f[10]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[63] += 6.708203932499369*rdxFnu[1]*f[15]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[73] += 22.9128784747792*rdxFnu[1]*f[16]; 
  out[74] += 22.9128784747792*rdxFnu[1]*f[17]; 
  out[75] += 22.9128784747792*rdxFnu[1]*f[18]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 
  out[79] += 22.9128784747792*rdxFnu[1]*f[35]; 

  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstDiffusionVol4xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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
  out[31] += 22.9128784747792*rdxFnu[0]*f[1]; 
  out[36] += 6.708203932499369*rdxFnu[0]*f[7]; 
  out[39] += 6.708203932499369*rdxFnu[0]*f[9]; 
  out[41] += 6.708203932499369*rdxFnu[0]*f[10]; 
  out[48] += 22.9128784747792*rdxFnu[0]*f[5]; 
  out[50] += 22.9128784747792*rdxFnu[0]*f[6]; 
  out[54] += 22.9128784747792*rdxFnu[0]*f[8]; 
  out[60] += 6.708203932499369*rdxFnu[0]*f[18]; 
  out[64] += 22.9128784747792*rdxFnu[0]*f[15]; 
  out[67] += 22.9128784747792*rdxFnu[0]*f[16]; 
  out[69] += 22.9128784747792*rdxFnu[0]*f[17]; 
  out[76] += 22.9128784747792*rdxFnu[0]*f[35]; 

  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 16.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion4Vol4xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 16.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2]+rdxFnu[3])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 
  rdxFnu[1] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[1]/(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 
  rdxFnu[1] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
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


  return (rdxFnu[0]+rdxFnu[1]+rdxFnu[2])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[2]/(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[2]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 
  rdxFnu[1] = 64.0*nu[3]/(dx[3]*dx[3]*dx[3]*dx[3]*dx[3]*dx[3]); 


  return (rdxFnu[0]+rdxFnu[1])*1.142857142857143;

} 
double ConstHyperDiffusion6Vol4xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[4]:  diffusion coefficient (collisionality).
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxFnu[1]; 
  rdxFnu[0] = 64.0*nu[0]/(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0]); 


  return (rdxFnu[0])*1.142857142857143;

} 
