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
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs2(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 

  out[2] += 2.5617376914899*rdxF[0]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[0]*nu[112]+1.677050983124842*rdxF[0]*f[2]*nu[92]+0.75*f[0]*rdxF[0]*nu[82]; 
  out[5] += 2.5617376914899*rdxF[0]*f[20]*nu[112]+1.14564392373896*rdxF[0]*f[1]*nu[112]+1.677050983124842*rdxF[0]*f[5]*nu[92]+0.75*rdxF[0]*f[1]*nu[82]; 
  out[7] += 2.5617376914899*rdxF[0]*f[22]*nu[112]+1.14564392373896*rdxF[0]*f[3]*nu[112]+1.677050983124842*rdxF[0]*f[7]*nu[92]+0.75*rdxF[0]*f[3]*nu[82]; 
  out[9] += 2.5617376914899*rdxF[0]*f[26]*nu[112]+1.14564392373896*rdxF[0]*f[4]*nu[112]+1.677050983124842*rdxF[0]*f[9]*nu[92]+0.75*rdxF[0]*f[4]*nu[82]; 
  out[12] += 6.708203932499369*rdxF[0]*f[32]*nu[112]+7.685213074469699*rdxF[0]*f[2]*nu[112]+5.031152949374527*rdxF[0]*f[12]*nu[92]+3.75*f[0]*rdxF[0]*nu[92]+3.354101966249685*rdxF[0]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[0]*nu[80]; 
  out[15] += 2.5617376914899*rdxF[0]*f[37]*nu[112]+1.14564392373896*rdxF[0]*f[6]*nu[112]+1.677050983124842*rdxF[0]*f[15]*nu[92]+0.75*rdxF[0]*f[6]*nu[82]; 
  out[16] += 2.5617376914899*rdxF[0]*f[40]*nu[112]+1.14564392373896*rdxF[0]*f[8]*nu[112]+1.677050983124842*rdxF[0]*f[16]*nu[92]+0.75*rdxF[0]*f[8]*nu[82]; 
  out[18] += 2.5617376914899*rdxF[0]*f[42]*nu[112]+1.14564392373896*rdxF[0]*f[10]*nu[112]+1.677050983124842*rdxF[0]*f[18]*nu[92]+0.75*rdxF[0]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[0]*f[11]*nu[112]+1.677050983124842*rdxF[0]*f[19]*nu[92]+0.75*rdxF[0]*f[11]*nu[82]; 
  out[20] += 6.708203932499369*rdxF[0]*f[49]*nu[112]+7.6852130744697*rdxF[0]*f[5]*nu[112]+5.031152949374527*rdxF[0]*f[20]*nu[92]+3.75*rdxF[0]*f[1]*nu[92]+3.354101966249684*rdxF[0]*f[5]*nu[82]+1.677050983124842*rdxF[0]*f[1]*nu[80]; 
  out[22] += 6.708203932499369*rdxF[0]*f[51]*nu[112]+7.6852130744697*rdxF[0]*f[7]*nu[112]+5.031152949374527*rdxF[0]*f[22]*nu[92]+3.75*rdxF[0]*f[3]*nu[92]+3.354101966249684*rdxF[0]*f[7]*nu[82]+1.677050983124842*rdxF[0]*f[3]*nu[80]; 
  out[24] += 1.14564392373896*rdxF[0]*f[13]*nu[112]+1.677050983124842*rdxF[0]*f[24]*nu[92]+0.75*rdxF[0]*f[13]*nu[82]; 
  out[26] += 6.708203932499369*rdxF[0]*f[55]*nu[112]+7.6852130744697*rdxF[0]*f[9]*nu[112]+5.031152949374527*rdxF[0]*f[26]*nu[92]+3.75*rdxF[0]*f[4]*nu[92]+3.354101966249684*rdxF[0]*f[9]*nu[82]+1.677050983124842*rdxF[0]*f[4]*nu[80]; 
  out[29] += 1.14564392373896*rdxF[0]*f[14]*nu[112]+1.677050983124842*rdxF[0]*f[29]*nu[92]+0.75*rdxF[0]*f[14]*nu[82]; 
  out[32] += 18.44756081437327*rdxF[0]*f[12]*nu[112]+10.5*f[0]*rdxF[0]*nu[112]+10.06230589874905*rdxF[0]*f[32]*nu[92]+12.8086884574495*rdxF[0]*f[2]*nu[92]+7.685213074469699*rdxF[0]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[0]*nu[82]+5.7282196186948*rdxF[0]*f[2]*nu[80]; 
  out[35] += 2.5617376914899*rdxF[0]*f[61]*nu[112]+1.14564392373896*rdxF[0]*f[17]*nu[112]+1.677050983124842*rdxF[0]*f[35]*nu[92]+0.75*rdxF[0]*f[17]*nu[82]; 
  out[36] += 1.14564392373896*rdxF[0]*f[21]*nu[112]+1.677050983124842*rdxF[0]*f[36]*nu[92]+0.75*rdxF[0]*f[21]*nu[82]; 
  out[37] += 6.708203932499369*rdxF[0]*f[65]*nu[112]+7.685213074469699*rdxF[0]*f[15]*nu[112]+5.031152949374527*rdxF[0]*f[37]*nu[92]+3.75*rdxF[0]*f[6]*nu[92]+3.354101966249685*rdxF[0]*f[15]*nu[82]+1.677050983124842*rdxF[0]*f[6]*nu[80]; 
  out[38] += 1.14564392373896*rdxF[0]*f[23]*nu[112]+1.677050983124842*rdxF[0]*f[38]*nu[92]+0.75*rdxF[0]*f[23]*nu[82]; 
  out[39] += 1.14564392373896*rdxF[0]*f[25]*nu[112]+1.677050983124842*rdxF[0]*f[39]*nu[92]+0.75*rdxF[0]*f[25]*nu[82]; 
  out[40] += 6.708203932499369*rdxF[0]*f[68]*nu[112]+7.685213074469699*rdxF[0]*f[16]*nu[112]+5.031152949374527*rdxF[0]*f[40]*nu[92]+3.75*rdxF[0]*f[8]*nu[92]+3.354101966249685*rdxF[0]*f[16]*nu[82]+1.677050983124842*rdxF[0]*f[8]*nu[80]; 
  out[42] += 6.708203932499369*rdxF[0]*f[70]*nu[112]+7.685213074469699*rdxF[0]*f[18]*nu[112]+5.031152949374527*rdxF[0]*f[42]*nu[92]+3.75*rdxF[0]*f[10]*nu[92]+3.354101966249685*rdxF[0]*f[18]*nu[82]+1.677050983124842*rdxF[0]*f[10]*nu[80]; 
  out[44] += 1.14564392373896*rdxF[0]*f[27]*nu[112]+1.677050983124842*rdxF[0]*f[44]*nu[92]+0.75*rdxF[0]*f[27]*nu[82]; 
  out[45] += 1.14564392373896*rdxF[0]*f[28]*nu[112]+1.677050983124842*rdxF[0]*f[45]*nu[92]+0.75*rdxF[0]*f[28]*nu[82]; 
  out[47] += 1.14564392373896*rdxF[0]*f[30]*nu[112]+1.677050983124842*rdxF[0]*f[47]*nu[92]+0.75*rdxF[0]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[0]*f[31]*nu[112]+1.677050983124842*rdxF[0]*f[48]*nu[92]+0.7499999999999999*rdxF[0]*f[31]*nu[82]; 
  out[49] += 18.44756081437327*rdxF[0]*f[20]*nu[112]+10.5*rdxF[0]*f[1]*nu[112]+10.06230589874905*rdxF[0]*f[49]*nu[92]+12.8086884574495*rdxF[0]*f[5]*nu[92]+7.685213074469698*rdxF[0]*f[20]*nu[82]+6.873863542433758*rdxF[0]*f[1]*nu[82]+5.7282196186948*rdxF[0]*f[5]*nu[80]; 
  out[51] += 18.44756081437327*rdxF[0]*f[22]*nu[112]+10.5*rdxF[0]*f[3]*nu[112]+10.06230589874905*rdxF[0]*f[51]*nu[92]+12.8086884574495*rdxF[0]*f[7]*nu[92]+7.685213074469698*rdxF[0]*f[22]*nu[82]+6.873863542433758*rdxF[0]*f[3]*nu[82]+5.7282196186948*rdxF[0]*f[7]*nu[80]; 
  out[53] += 1.14564392373896*rdxF[0]*f[33]*nu[112]+1.677050983124842*rdxF[0]*f[53]*nu[92]+0.7499999999999999*rdxF[0]*f[33]*nu[82]; 
  out[55] += 18.44756081437327*rdxF[0]*f[26]*nu[112]+10.5*rdxF[0]*f[4]*nu[112]+10.06230589874905*rdxF[0]*f[55]*nu[92]+12.8086884574495*rdxF[0]*f[9]*nu[92]+7.685213074469698*rdxF[0]*f[26]*nu[82]+6.873863542433758*rdxF[0]*f[4]*nu[82]+5.7282196186948*rdxF[0]*f[9]*nu[80]; 
  out[58] += 1.14564392373896*rdxF[0]*f[34]*nu[112]+1.677050983124842*rdxF[0]*f[58]*nu[92]+0.7499999999999999*rdxF[0]*f[34]*nu[82]; 
  out[60] += 1.14564392373896*rdxF[0]*f[41]*nu[112]+1.677050983124842*rdxF[0]*f[60]*nu[92]+0.75*rdxF[0]*f[41]*nu[82]; 
  out[61] += 6.708203932499369*rdxF[0]*f[77]*nu[112]+7.6852130744697*rdxF[0]*f[35]*nu[112]+5.031152949374527*rdxF[0]*f[61]*nu[92]+3.75*rdxF[0]*f[17]*nu[92]+3.354101966249684*rdxF[0]*f[35]*nu[82]+1.677050983124842*rdxF[0]*f[17]*nu[80]; 
  out[62] += 1.14564392373896*rdxF[0]*f[43]*nu[112]+1.677050983124842*rdxF[0]*f[62]*nu[92]+0.75*rdxF[0]*f[43]*nu[82]; 
  out[63] += 1.14564392373896*rdxF[0]*f[46]*nu[112]+1.677050983124842*rdxF[0]*f[63]*nu[92]+0.75*rdxF[0]*f[46]*nu[82]; 
  out[64] += 1.14564392373896*rdxF[0]*f[50]*nu[112]+1.677050983124842*rdxF[0]*f[64]*nu[92]+0.7499999999999999*rdxF[0]*f[50]*nu[82]; 
  out[65] += 18.44756081437327*rdxF[0]*f[37]*nu[112]+10.5*rdxF[0]*f[6]*nu[112]+10.06230589874905*rdxF[0]*f[65]*nu[92]+12.8086884574495*rdxF[0]*f[15]*nu[92]+7.685213074469699*rdxF[0]*f[37]*nu[82]+6.873863542433759*rdxF[0]*f[6]*nu[82]+5.7282196186948*rdxF[0]*f[15]*nu[80]; 
  out[66] += 1.14564392373896*rdxF[0]*f[52]*nu[112]+1.677050983124842*rdxF[0]*f[66]*nu[92]+0.7499999999999999*rdxF[0]*f[52]*nu[82]; 
  out[67] += 1.14564392373896*rdxF[0]*f[54]*nu[112]+1.677050983124842*rdxF[0]*f[67]*nu[92]+0.7499999999999999*rdxF[0]*f[54]*nu[82]; 
  out[68] += 18.44756081437327*rdxF[0]*f[40]*nu[112]+10.5*rdxF[0]*f[8]*nu[112]+10.06230589874905*rdxF[0]*f[68]*nu[92]+12.8086884574495*rdxF[0]*f[16]*nu[92]+7.685213074469699*rdxF[0]*f[40]*nu[82]+6.873863542433759*rdxF[0]*f[8]*nu[82]+5.7282196186948*rdxF[0]*f[16]*nu[80]; 
  out[70] += 18.44756081437327*rdxF[0]*f[42]*nu[112]+10.5*rdxF[0]*f[10]*nu[112]+10.06230589874905*rdxF[0]*f[70]*nu[92]+12.8086884574495*rdxF[0]*f[18]*nu[92]+7.685213074469699*rdxF[0]*f[42]*nu[82]+6.873863542433759*rdxF[0]*f[10]*nu[82]+5.7282196186948*rdxF[0]*f[18]*nu[80]; 
  out[72] += 1.14564392373896*rdxF[0]*f[56]*nu[112]+1.677050983124842*rdxF[0]*f[72]*nu[92]+0.7499999999999999*rdxF[0]*f[56]*nu[82]; 
  out[73] += 1.14564392373896*rdxF[0]*f[57]*nu[112]+1.677050983124842*rdxF[0]*f[73]*nu[92]+0.7499999999999999*rdxF[0]*f[57]*nu[82]; 
  out[75] += 1.14564392373896*rdxF[0]*f[59]*nu[112]+1.677050983124842*rdxF[0]*f[75]*nu[92]+0.7499999999999999*rdxF[0]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[0]*f[69]*nu[112]+1.677050983124842*rdxF[0]*f[76]*nu[92]+0.7499999999999999*rdxF[0]*f[69]*nu[82]; 
  out[77] += 18.44756081437327*rdxF[0]*f[61]*nu[112]+10.5*rdxF[0]*f[17]*nu[112]+10.06230589874905*rdxF[0]*f[77]*nu[92]+12.8086884574495*rdxF[0]*f[35]*nu[92]+7.685213074469698*rdxF[0]*f[61]*nu[82]+6.873863542433758*rdxF[0]*f[17]*nu[82]+5.7282196186948*rdxF[0]*f[35]*nu[80]; 
  out[78] += 1.14564392373896*rdxF[0]*f[71]*nu[112]+1.677050983124842*rdxF[0]*f[78]*nu[92]+0.7499999999999999*rdxF[0]*f[71]*nu[82]; 
  out[79] += 1.14564392373896*rdxF[0]*f[74]*nu[112]+1.677050983124842*rdxF[0]*f[79]*nu[92]+0.7499999999999999*rdxF[0]*f[74]*nu[82]; 

  return (rdxF[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs23(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[2] += 2.5617376914899*rdxF[0]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[0]*nu[112]+1.677050983124842*rdxF[0]*f[2]*nu[92]+0.75*f[0]*rdxF[0]*nu[82]; 
  out[3] += 2.5617376914899*rdxF[1]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[3]*nu[173]+0.75*f[0]*rdxF[1]*nu[163]; 
  out[5] += 2.5617376914899*rdxF[0]*f[20]*nu[112]+1.14564392373896*rdxF[0]*f[1]*nu[112]+1.677050983124842*rdxF[0]*f[5]*nu[92]+0.75*rdxF[0]*f[1]*nu[82]; 
  out[6] += 2.5617376914899*rdxF[1]*f[23]*nu[193]+1.14564392373896*f[1]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[6]*nu[173]+0.75*f[1]*rdxF[1]*nu[163]; 
  out[7] += 2.5617376914899*rdxF[1]*f[24]*nu[193]+1.14564392373896*rdxF[1]*f[2]*nu[193]+1.677050983124842*rdxF[1]*f[7]*nu[173]+0.75*rdxF[1]*f[2]*nu[163]+2.5617376914899*rdxF[0]*f[22]*nu[112]+1.14564392373896*rdxF[0]*f[3]*nu[112]+1.677050983124842*rdxF[0]*f[7]*nu[92]+0.75*rdxF[0]*f[3]*nu[82]; 
  out[9] += 2.5617376914899*rdxF[0]*f[26]*nu[112]+1.14564392373896*rdxF[0]*f[4]*nu[112]+1.677050983124842*rdxF[0]*f[9]*nu[92]+0.75*rdxF[0]*f[4]*nu[82]; 
  out[10] += 2.5617376914899*rdxF[1]*f[27]*nu[193]+1.14564392373896*rdxF[1]*f[4]*nu[193]+1.677050983124842*rdxF[1]*f[10]*nu[173]+0.75*rdxF[1]*f[4]*nu[163]; 
  out[12] += 6.708203932499369*rdxF[0]*f[32]*nu[112]+7.685213074469699*rdxF[0]*f[2]*nu[112]+5.031152949374527*rdxF[0]*f[12]*nu[92]+3.75*f[0]*rdxF[0]*nu[92]+3.354101966249685*rdxF[0]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[0]*nu[80]; 
  out[13] += 6.708203932499369*rdxF[1]*f[33]*nu[193]+7.685213074469699*rdxF[1]*f[3]*nu[193]+5.031152949374527*rdxF[1]*f[13]*nu[173]+3.75*f[0]*rdxF[1]*nu[173]+3.354101966249685*rdxF[1]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[1]*nu[160]; 
  out[15] += 2.5617376914899*rdxF[1]*f[38]*nu[193]+1.14564392373896*rdxF[1]*f[5]*nu[193]+1.677050983124842*rdxF[1]*f[15]*nu[173]+0.75*rdxF[1]*f[5]*nu[163]+2.5617376914899*rdxF[0]*f[37]*nu[112]+1.14564392373896*rdxF[0]*f[6]*nu[112]+1.677050983124842*rdxF[0]*f[15]*nu[92]+0.75*rdxF[0]*f[6]*nu[82]; 
  out[16] += 2.5617376914899*rdxF[0]*f[40]*nu[112]+1.14564392373896*rdxF[0]*f[8]*nu[112]+1.677050983124842*rdxF[0]*f[16]*nu[92]+0.75*rdxF[0]*f[8]*nu[82]; 
  out[17] += 2.5617376914899*rdxF[1]*f[43]*nu[193]+1.14564392373896*rdxF[1]*f[8]*nu[193]+1.677050983124842*rdxF[1]*f[17]*nu[173]+0.75*rdxF[1]*f[8]*nu[163]; 
  out[18] += 2.5617376914899*rdxF[1]*f[44]*nu[193]+1.14564392373896*rdxF[1]*f[9]*nu[193]+1.677050983124842*rdxF[1]*f[18]*nu[173]+0.75*rdxF[1]*f[9]*nu[163]+2.5617376914899*rdxF[0]*f[42]*nu[112]+1.14564392373896*rdxF[0]*f[10]*nu[112]+1.677050983124842*rdxF[0]*f[18]*nu[92]+0.75*rdxF[0]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[0]*f[11]*nu[112]+1.677050983124842*rdxF[0]*f[19]*nu[92]+0.75*rdxF[0]*f[11]*nu[82]; 
  out[20] += 6.708203932499369*rdxF[0]*f[49]*nu[112]+7.6852130744697*rdxF[0]*f[5]*nu[112]+5.031152949374527*rdxF[0]*f[20]*nu[92]+3.75*rdxF[0]*f[1]*nu[92]+3.354101966249684*rdxF[0]*f[5]*nu[82]+1.677050983124842*rdxF[0]*f[1]*nu[80]; 
  out[21] += 1.14564392373896*rdxF[1]*f[11]*nu[193]+1.677050983124842*rdxF[1]*f[21]*nu[173]+0.75*rdxF[1]*f[11]*nu[163]; 
  out[22] += 1.14564392373896*rdxF[1]*f[12]*nu[193]+1.677050983124842*rdxF[1]*f[22]*nu[173]+0.75*rdxF[1]*f[12]*nu[163]+6.708203932499369*rdxF[0]*f[51]*nu[112]+7.6852130744697*rdxF[0]*f[7]*nu[112]+5.031152949374527*rdxF[0]*f[22]*nu[92]+3.75*rdxF[0]*f[3]*nu[92]+3.354101966249684*rdxF[0]*f[7]*nu[82]+1.677050983124842*rdxF[0]*f[3]*nu[80]; 
  out[23] += 6.708203932499369*rdxF[1]*f[52]*nu[193]+7.6852130744697*rdxF[1]*f[6]*nu[193]+5.031152949374527*rdxF[1]*f[23]*nu[173]+3.75*f[1]*rdxF[1]*nu[173]+3.354101966249684*rdxF[1]*f[6]*nu[163]+1.677050983124842*f[1]*rdxF[1]*nu[160]; 
  out[24] += 6.708203932499369*rdxF[1]*f[53]*nu[193]+7.6852130744697*rdxF[1]*f[7]*nu[193]+5.031152949374527*rdxF[1]*f[24]*nu[173]+3.75*rdxF[1]*f[2]*nu[173]+3.354101966249684*rdxF[1]*f[7]*nu[163]+1.677050983124842*rdxF[1]*f[2]*nu[160]+1.14564392373896*rdxF[0]*f[13]*nu[112]+1.677050983124842*rdxF[0]*f[24]*nu[92]+0.75*rdxF[0]*f[13]*nu[82]; 
  out[26] += 6.708203932499369*rdxF[0]*f[55]*nu[112]+7.6852130744697*rdxF[0]*f[9]*nu[112]+5.031152949374527*rdxF[0]*f[26]*nu[92]+3.75*rdxF[0]*f[4]*nu[92]+3.354101966249684*rdxF[0]*f[9]*nu[82]+1.677050983124842*rdxF[0]*f[4]*nu[80]; 
  out[27] += 6.708203932499369*rdxF[1]*f[56]*nu[193]+7.6852130744697*rdxF[1]*f[10]*nu[193]+5.031152949374527*rdxF[1]*f[27]*nu[173]+3.75*rdxF[1]*f[4]*nu[173]+3.354101966249684*rdxF[1]*f[10]*nu[163]+1.677050983124842*rdxF[1]*f[4]*nu[160]; 
  out[29] += 1.14564392373896*rdxF[0]*f[14]*nu[112]+1.677050983124842*rdxF[0]*f[29]*nu[92]+0.75*rdxF[0]*f[14]*nu[82]; 
  out[30] += 1.14564392373896*rdxF[1]*f[14]*nu[193]+1.677050983124842*rdxF[1]*f[30]*nu[173]+0.75*rdxF[1]*f[14]*nu[163]; 
  out[32] += 18.44756081437327*rdxF[0]*f[12]*nu[112]+10.5*f[0]*rdxF[0]*nu[112]+10.06230589874905*rdxF[0]*f[32]*nu[92]+12.8086884574495*rdxF[0]*f[2]*nu[92]+7.685213074469699*rdxF[0]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[0]*nu[82]+5.7282196186948*rdxF[0]*f[2]*nu[80]; 
  out[33] += 18.44756081437327*rdxF[1]*f[13]*nu[193]+10.5*f[0]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[33]*nu[173]+12.8086884574495*rdxF[1]*f[3]*nu[173]+7.685213074469699*rdxF[1]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[3]*nu[160]; 
  out[35] += 2.5617376914899*rdxF[1]*f[62]*nu[193]+1.14564392373896*rdxF[1]*f[16]*nu[193]+1.677050983124842*rdxF[1]*f[35]*nu[173]+0.75*rdxF[1]*f[16]*nu[163]+2.5617376914899*rdxF[0]*f[61]*nu[112]+1.14564392373896*rdxF[0]*f[17]*nu[112]+1.677050983124842*rdxF[0]*f[35]*nu[92]+0.75*rdxF[0]*f[17]*nu[82]; 
  out[36] += 1.14564392373896*rdxF[1]*f[19]*nu[193]+1.677050983124842*rdxF[1]*f[36]*nu[173]+0.75*rdxF[1]*f[19]*nu[163]+1.14564392373896*rdxF[0]*f[21]*nu[112]+1.677050983124842*rdxF[0]*f[36]*nu[92]+0.75*rdxF[0]*f[21]*nu[82]; 
  out[37] += 1.14564392373896*rdxF[1]*f[20]*nu[193]+1.677050983124842*rdxF[1]*f[37]*nu[173]+0.75*rdxF[1]*f[20]*nu[163]+6.708203932499369*rdxF[0]*f[65]*nu[112]+7.685213074469699*rdxF[0]*f[15]*nu[112]+5.031152949374527*rdxF[0]*f[37]*nu[92]+3.75*rdxF[0]*f[6]*nu[92]+3.354101966249685*rdxF[0]*f[15]*nu[82]+1.677050983124842*rdxF[0]*f[6]*nu[80]; 
  out[38] += 6.708203932499369*rdxF[1]*f[66]*nu[193]+7.685213074469699*rdxF[1]*f[15]*nu[193]+5.031152949374527*rdxF[1]*f[38]*nu[173]+3.75*rdxF[1]*f[5]*nu[173]+3.354101966249685*rdxF[1]*f[15]*nu[163]+1.677050983124842*rdxF[1]*f[5]*nu[160]+1.14564392373896*rdxF[0]*f[23]*nu[112]+1.677050983124842*rdxF[0]*f[38]*nu[92]+0.75*rdxF[0]*f[23]*nu[82]; 
  out[39] += 1.14564392373896*rdxF[0]*f[25]*nu[112]+1.677050983124842*rdxF[0]*f[39]*nu[92]+0.75*rdxF[0]*f[25]*nu[82]; 
  out[40] += 6.708203932499369*rdxF[0]*f[68]*nu[112]+7.685213074469699*rdxF[0]*f[16]*nu[112]+5.031152949374527*rdxF[0]*f[40]*nu[92]+3.75*rdxF[0]*f[8]*nu[92]+3.354101966249685*rdxF[0]*f[16]*nu[82]+1.677050983124842*rdxF[0]*f[8]*nu[80]; 
  out[41] += 1.14564392373896*rdxF[1]*f[25]*nu[193]+1.677050983124842*rdxF[1]*f[41]*nu[173]+0.75*rdxF[1]*f[25]*nu[163]; 
  out[42] += 1.14564392373896*rdxF[1]*f[26]*nu[193]+1.677050983124842*rdxF[1]*f[42]*nu[173]+0.75*rdxF[1]*f[26]*nu[163]+6.708203932499369*rdxF[0]*f[70]*nu[112]+7.685213074469699*rdxF[0]*f[18]*nu[112]+5.031152949374527*rdxF[0]*f[42]*nu[92]+3.75*rdxF[0]*f[10]*nu[92]+3.354101966249685*rdxF[0]*f[18]*nu[82]+1.677050983124842*rdxF[0]*f[10]*nu[80]; 
  out[43] += 6.708203932499369*rdxF[1]*f[71]*nu[193]+7.685213074469699*rdxF[1]*f[17]*nu[193]+5.031152949374527*rdxF[1]*f[43]*nu[173]+3.75*rdxF[1]*f[8]*nu[173]+3.354101966249685*rdxF[1]*f[17]*nu[163]+1.677050983124842*rdxF[1]*f[8]*nu[160]; 
  out[44] += 6.708203932499369*rdxF[1]*f[72]*nu[193]+7.685213074469699*rdxF[1]*f[18]*nu[193]+5.031152949374527*rdxF[1]*f[44]*nu[173]+3.75*rdxF[1]*f[9]*nu[173]+3.354101966249685*rdxF[1]*f[18]*nu[163]+1.677050983124842*rdxF[1]*f[9]*nu[160]+1.14564392373896*rdxF[0]*f[27]*nu[112]+1.677050983124842*rdxF[0]*f[44]*nu[92]+0.75*rdxF[0]*f[27]*nu[82]; 
  out[45] += 1.14564392373896*rdxF[0]*f[28]*nu[112]+1.677050983124842*rdxF[0]*f[45]*nu[92]+0.75*rdxF[0]*f[28]*nu[82]; 
  out[46] += 1.14564392373896*rdxF[1]*f[28]*nu[193]+1.677050983124842*rdxF[1]*f[46]*nu[173]+0.75*rdxF[1]*f[28]*nu[163]; 
  out[47] += 1.14564392373896*rdxF[1]*f[29]*nu[193]+1.677050983124842*rdxF[1]*f[47]*nu[173]+0.75*rdxF[1]*f[29]*nu[163]+1.14564392373896*rdxF[0]*f[30]*nu[112]+1.677050983124842*rdxF[0]*f[47]*nu[92]+0.75*rdxF[0]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[0]*f[31]*nu[112]+1.677050983124842*rdxF[0]*f[48]*nu[92]+0.7499999999999999*rdxF[0]*f[31]*nu[82]; 
  out[49] += 18.44756081437327*rdxF[0]*f[20]*nu[112]+10.5*rdxF[0]*f[1]*nu[112]+10.06230589874905*rdxF[0]*f[49]*nu[92]+12.8086884574495*rdxF[0]*f[5]*nu[92]+7.685213074469698*rdxF[0]*f[20]*nu[82]+6.873863542433758*rdxF[0]*f[1]*nu[82]+5.7282196186948*rdxF[0]*f[5]*nu[80]; 
  out[50] += 1.14564392373896*rdxF[1]*f[31]*nu[193]+1.677050983124842*rdxF[1]*f[50]*nu[173]+0.7499999999999999*rdxF[1]*f[31]*nu[163]; 
  out[51] += 1.14564392373896*rdxF[1]*f[32]*nu[193]+1.677050983124842*rdxF[1]*f[51]*nu[173]+0.7499999999999999*rdxF[1]*f[32]*nu[163]+18.44756081437327*rdxF[0]*f[22]*nu[112]+10.5*rdxF[0]*f[3]*nu[112]+10.06230589874905*rdxF[0]*f[51]*nu[92]+12.8086884574495*rdxF[0]*f[7]*nu[92]+7.685213074469698*rdxF[0]*f[22]*nu[82]+6.873863542433758*rdxF[0]*f[3]*nu[82]+5.7282196186948*rdxF[0]*f[7]*nu[80]; 
  out[52] += 18.44756081437327*rdxF[1]*f[23]*nu[193]+10.5*f[1]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[52]*nu[173]+12.8086884574495*rdxF[1]*f[6]*nu[173]+7.685213074469698*rdxF[1]*f[23]*nu[163]+6.873863542433758*f[1]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[6]*nu[160]; 
  out[53] += 18.44756081437327*rdxF[1]*f[24]*nu[193]+10.5*rdxF[1]*f[2]*nu[193]+10.06230589874905*rdxF[1]*f[53]*nu[173]+12.8086884574495*rdxF[1]*f[7]*nu[173]+7.685213074469698*rdxF[1]*f[24]*nu[163]+6.873863542433758*rdxF[1]*f[2]*nu[163]+5.7282196186948*rdxF[1]*f[7]*nu[160]+1.14564392373896*rdxF[0]*f[33]*nu[112]+1.677050983124842*rdxF[0]*f[53]*nu[92]+0.7499999999999999*rdxF[0]*f[33]*nu[82]; 
  out[55] += 18.44756081437327*rdxF[0]*f[26]*nu[112]+10.5*rdxF[0]*f[4]*nu[112]+10.06230589874905*rdxF[0]*f[55]*nu[92]+12.8086884574495*rdxF[0]*f[9]*nu[92]+7.685213074469698*rdxF[0]*f[26]*nu[82]+6.873863542433758*rdxF[0]*f[4]*nu[82]+5.7282196186948*rdxF[0]*f[9]*nu[80]; 
  out[56] += 18.44756081437327*rdxF[1]*f[27]*nu[193]+10.5*rdxF[1]*f[4]*nu[193]+10.06230589874905*rdxF[1]*f[56]*nu[173]+12.8086884574495*rdxF[1]*f[10]*nu[173]+7.685213074469698*rdxF[1]*f[27]*nu[163]+6.873863542433758*rdxF[1]*f[4]*nu[163]+5.7282196186948*rdxF[1]*f[10]*nu[160]; 
  out[58] += 1.14564392373896*rdxF[0]*f[34]*nu[112]+1.677050983124842*rdxF[0]*f[58]*nu[92]+0.7499999999999999*rdxF[0]*f[34]*nu[82]; 
  out[59] += 1.14564392373896*rdxF[1]*f[34]*nu[193]+1.677050983124842*rdxF[1]*f[59]*nu[173]+0.7499999999999999*rdxF[1]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[1]*f[39]*nu[193]+1.677050983124842*rdxF[1]*f[60]*nu[173]+0.75*rdxF[1]*f[39]*nu[163]+1.14564392373896*rdxF[0]*f[41]*nu[112]+1.677050983124842*rdxF[0]*f[60]*nu[92]+0.75*rdxF[0]*f[41]*nu[82]; 
  out[61] += 1.14564392373896*rdxF[1]*f[40]*nu[193]+1.677050983124842*rdxF[1]*f[61]*nu[173]+0.75*rdxF[1]*f[40]*nu[163]+6.708203932499369*rdxF[0]*f[77]*nu[112]+7.6852130744697*rdxF[0]*f[35]*nu[112]+5.031152949374527*rdxF[0]*f[61]*nu[92]+3.75*rdxF[0]*f[17]*nu[92]+3.354101966249684*rdxF[0]*f[35]*nu[82]+1.677050983124842*rdxF[0]*f[17]*nu[80]; 
  out[62] += 6.708203932499369*rdxF[1]*f[78]*nu[193]+7.6852130744697*rdxF[1]*f[35]*nu[193]+5.031152949374527*rdxF[1]*f[62]*nu[173]+3.75*rdxF[1]*f[16]*nu[173]+3.354101966249684*rdxF[1]*f[35]*nu[163]+1.677050983124842*rdxF[1]*f[16]*nu[160]+1.14564392373896*rdxF[0]*f[43]*nu[112]+1.677050983124842*rdxF[0]*f[62]*nu[92]+0.75*rdxF[0]*f[43]*nu[82]; 
  out[63] += 1.14564392373896*rdxF[1]*f[45]*nu[193]+1.677050983124842*rdxF[1]*f[63]*nu[173]+0.75*rdxF[1]*f[45]*nu[163]+1.14564392373896*rdxF[0]*f[46]*nu[112]+1.677050983124842*rdxF[0]*f[63]*nu[92]+0.75*rdxF[0]*f[46]*nu[82]; 
  out[64] += 1.14564392373896*rdxF[1]*f[48]*nu[193]+1.677050983124842*rdxF[1]*f[64]*nu[173]+0.7499999999999999*rdxF[1]*f[48]*nu[163]+1.14564392373896*rdxF[0]*f[50]*nu[112]+1.677050983124842*rdxF[0]*f[64]*nu[92]+0.7499999999999999*rdxF[0]*f[50]*nu[82]; 
  out[65] += 1.14564392373896*rdxF[1]*f[49]*nu[193]+1.677050983124842*rdxF[1]*f[65]*nu[173]+0.7499999999999999*rdxF[1]*f[49]*nu[163]+18.44756081437327*rdxF[0]*f[37]*nu[112]+10.5*rdxF[0]*f[6]*nu[112]+10.06230589874905*rdxF[0]*f[65]*nu[92]+12.8086884574495*rdxF[0]*f[15]*nu[92]+7.685213074469699*rdxF[0]*f[37]*nu[82]+6.873863542433759*rdxF[0]*f[6]*nu[82]+5.7282196186948*rdxF[0]*f[15]*nu[80]; 
  out[66] += 18.44756081437327*rdxF[1]*f[38]*nu[193]+10.5*rdxF[1]*f[5]*nu[193]+10.06230589874905*rdxF[1]*f[66]*nu[173]+12.8086884574495*rdxF[1]*f[15]*nu[173]+7.685213074469699*rdxF[1]*f[38]*nu[163]+6.873863542433759*rdxF[1]*f[5]*nu[163]+5.7282196186948*rdxF[1]*f[15]*nu[160]+1.14564392373896*rdxF[0]*f[52]*nu[112]+1.677050983124842*rdxF[0]*f[66]*nu[92]+0.7499999999999999*rdxF[0]*f[52]*nu[82]; 
  out[67] += 1.14564392373896*rdxF[0]*f[54]*nu[112]+1.677050983124842*rdxF[0]*f[67]*nu[92]+0.7499999999999999*rdxF[0]*f[54]*nu[82]; 
  out[68] += 18.44756081437327*rdxF[0]*f[40]*nu[112]+10.5*rdxF[0]*f[8]*nu[112]+10.06230589874905*rdxF[0]*f[68]*nu[92]+12.8086884574495*rdxF[0]*f[16]*nu[92]+7.685213074469699*rdxF[0]*f[40]*nu[82]+6.873863542433759*rdxF[0]*f[8]*nu[82]+5.7282196186948*rdxF[0]*f[16]*nu[80]; 
  out[69] += 1.14564392373896*rdxF[1]*f[54]*nu[193]+1.677050983124842*rdxF[1]*f[69]*nu[173]+0.7499999999999999*rdxF[1]*f[54]*nu[163]; 
  out[70] += 1.14564392373896*rdxF[1]*f[55]*nu[193]+1.677050983124842*rdxF[1]*f[70]*nu[173]+0.7499999999999999*rdxF[1]*f[55]*nu[163]+18.44756081437327*rdxF[0]*f[42]*nu[112]+10.5*rdxF[0]*f[10]*nu[112]+10.06230589874905*rdxF[0]*f[70]*nu[92]+12.8086884574495*rdxF[0]*f[18]*nu[92]+7.685213074469699*rdxF[0]*f[42]*nu[82]+6.873863542433759*rdxF[0]*f[10]*nu[82]+5.7282196186948*rdxF[0]*f[18]*nu[80]; 
  out[71] += 18.44756081437327*rdxF[1]*f[43]*nu[193]+10.5*rdxF[1]*f[8]*nu[193]+10.06230589874905*rdxF[1]*f[71]*nu[173]+12.8086884574495*rdxF[1]*f[17]*nu[173]+7.685213074469699*rdxF[1]*f[43]*nu[163]+6.873863542433759*rdxF[1]*f[8]*nu[163]+5.7282196186948*rdxF[1]*f[17]*nu[160]; 
  out[72] += 18.44756081437327*rdxF[1]*f[44]*nu[193]+10.5*rdxF[1]*f[9]*nu[193]+10.06230589874905*rdxF[1]*f[72]*nu[173]+12.8086884574495*rdxF[1]*f[18]*nu[173]+7.685213074469699*rdxF[1]*f[44]*nu[163]+6.873863542433759*rdxF[1]*f[9]*nu[163]+5.7282196186948*rdxF[1]*f[18]*nu[160]+1.14564392373896*rdxF[0]*f[56]*nu[112]+1.677050983124842*rdxF[0]*f[72]*nu[92]+0.7499999999999999*rdxF[0]*f[56]*nu[82]; 
  out[73] += 1.14564392373896*rdxF[0]*f[57]*nu[112]+1.677050983124842*rdxF[0]*f[73]*nu[92]+0.7499999999999999*rdxF[0]*f[57]*nu[82]; 
  out[74] += 1.14564392373896*rdxF[1]*f[57]*nu[193]+1.677050983124842*rdxF[1]*f[74]*nu[173]+0.7499999999999999*rdxF[1]*f[57]*nu[163]; 
  out[75] += 1.14564392373896*rdxF[1]*f[58]*nu[193]+1.677050983124842*rdxF[1]*f[75]*nu[173]+0.7499999999999999*rdxF[1]*f[58]*nu[163]+1.14564392373896*rdxF[0]*f[59]*nu[112]+1.677050983124842*rdxF[0]*f[75]*nu[92]+0.7499999999999999*rdxF[0]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[1]*f[67]*nu[193]+1.677050983124842*rdxF[1]*f[76]*nu[173]+0.7499999999999999*rdxF[1]*f[67]*nu[163]+1.14564392373896*rdxF[0]*f[69]*nu[112]+1.677050983124842*rdxF[0]*f[76]*nu[92]+0.7499999999999999*rdxF[0]*f[69]*nu[82]; 
  out[77] += 1.14564392373896*rdxF[1]*f[68]*nu[193]+1.677050983124842*rdxF[1]*f[77]*nu[173]+0.7499999999999999*rdxF[1]*f[68]*nu[163]+18.44756081437327*rdxF[0]*f[61]*nu[112]+10.5*rdxF[0]*f[17]*nu[112]+10.06230589874905*rdxF[0]*f[77]*nu[92]+12.8086884574495*rdxF[0]*f[35]*nu[92]+7.685213074469698*rdxF[0]*f[61]*nu[82]+6.873863542433758*rdxF[0]*f[17]*nu[82]+5.7282196186948*rdxF[0]*f[35]*nu[80]; 
  out[78] += 18.44756081437327*rdxF[1]*f[62]*nu[193]+10.5*rdxF[1]*f[16]*nu[193]+10.06230589874905*rdxF[1]*f[78]*nu[173]+12.8086884574495*rdxF[1]*f[35]*nu[173]+7.685213074469698*rdxF[1]*f[62]*nu[163]+6.873863542433758*rdxF[1]*f[16]*nu[163]+5.7282196186948*rdxF[1]*f[35]*nu[160]+1.14564392373896*rdxF[0]*f[71]*nu[112]+1.677050983124842*rdxF[0]*f[78]*nu[92]+0.7499999999999999*rdxF[0]*f[71]*nu[82]; 
  out[79] += 1.14564392373896*rdxF[1]*f[73]*nu[193]+1.677050983124842*rdxF[1]*f[79]*nu[173]+0.7499999999999999*rdxF[1]*f[73]*nu[163]+1.14564392373896*rdxF[0]*f[74]*nu[112]+1.677050983124842*rdxF[0]*f[79]*nu[92]+0.7499999999999999*rdxF[0]*f[74]*nu[82]; 

  return (rdxF[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+rdxF[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[2] += 2.5617376914899*rdxF[0]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[0]*nu[112]+1.677050983124842*rdxF[0]*f[2]*nu[92]+0.75*f[0]*rdxF[0]*nu[82]; 
  out[3] += 2.5617376914899*rdxF[1]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[3]*nu[173]+0.75*f[0]*rdxF[1]*nu[163]; 
  out[4] += 2.5617376914899*rdxF[2]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[4]*nu[254]+0.75*f[0]*rdxF[2]*nu[244]; 
  out[5] += 2.5617376914899*rdxF[0]*f[20]*nu[112]+1.14564392373896*rdxF[0]*f[1]*nu[112]+1.677050983124842*rdxF[0]*f[5]*nu[92]+0.75*rdxF[0]*f[1]*nu[82]; 
  out[6] += 2.5617376914899*rdxF[1]*f[23]*nu[193]+1.14564392373896*f[1]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[6]*nu[173]+0.75*f[1]*rdxF[1]*nu[163]; 
  out[7] += 2.5617376914899*rdxF[1]*f[24]*nu[193]+1.14564392373896*rdxF[1]*f[2]*nu[193]+1.677050983124842*rdxF[1]*f[7]*nu[173]+0.75*rdxF[1]*f[2]*nu[163]+2.5617376914899*rdxF[0]*f[22]*nu[112]+1.14564392373896*rdxF[0]*f[3]*nu[112]+1.677050983124842*rdxF[0]*f[7]*nu[92]+0.75*rdxF[0]*f[3]*nu[82]; 
  out[8] += 2.5617376914899*rdxF[2]*f[28]*nu[274]+1.14564392373896*f[1]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[8]*nu[254]+0.75*f[1]*rdxF[2]*nu[244]; 
  out[9] += 2.5617376914899*rdxF[2]*f[29]*nu[274]+1.14564392373896*f[2]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[9]*nu[254]+0.75*f[2]*rdxF[2]*nu[244]+2.5617376914899*rdxF[0]*f[26]*nu[112]+1.14564392373896*rdxF[0]*f[4]*nu[112]+1.677050983124842*rdxF[0]*f[9]*nu[92]+0.75*rdxF[0]*f[4]*nu[82]; 
  out[10] += 2.5617376914899*rdxF[2]*f[30]*nu[274]+1.14564392373896*rdxF[2]*f[3]*nu[274]+1.677050983124842*rdxF[2]*f[10]*nu[254]+0.75*rdxF[2]*f[3]*nu[244]+2.5617376914899*rdxF[1]*f[27]*nu[193]+1.14564392373896*rdxF[1]*f[4]*nu[193]+1.677050983124842*rdxF[1]*f[10]*nu[173]+0.75*rdxF[1]*f[4]*nu[163]; 
  out[12] += 6.708203932499369*rdxF[0]*f[32]*nu[112]+7.685213074469699*rdxF[0]*f[2]*nu[112]+5.031152949374527*rdxF[0]*f[12]*nu[92]+3.75*f[0]*rdxF[0]*nu[92]+3.354101966249685*rdxF[0]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[0]*nu[80]; 
  out[13] += 6.708203932499369*rdxF[1]*f[33]*nu[193]+7.685213074469699*rdxF[1]*f[3]*nu[193]+5.031152949374527*rdxF[1]*f[13]*nu[173]+3.75*f[0]*rdxF[1]*nu[173]+3.354101966249685*rdxF[1]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[1]*nu[160]; 
  out[14] += 6.708203932499369*rdxF[2]*f[34]*nu[274]+7.685213074469699*rdxF[2]*f[4]*nu[274]+5.031152949374527*rdxF[2]*f[14]*nu[254]+3.75*f[0]*rdxF[2]*nu[254]+3.354101966249685*rdxF[2]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[2]*nu[240]; 
  out[15] += 2.5617376914899*rdxF[1]*f[38]*nu[193]+1.14564392373896*rdxF[1]*f[5]*nu[193]+1.677050983124842*rdxF[1]*f[15]*nu[173]+0.75*rdxF[1]*f[5]*nu[163]+2.5617376914899*rdxF[0]*f[37]*nu[112]+1.14564392373896*rdxF[0]*f[6]*nu[112]+1.677050983124842*rdxF[0]*f[15]*nu[92]+0.75*rdxF[0]*f[6]*nu[82]; 
  out[16] += 2.5617376914899*rdxF[2]*f[45]*nu[274]+1.14564392373896*rdxF[2]*f[5]*nu[274]+1.677050983124842*rdxF[2]*f[16]*nu[254]+0.75*rdxF[2]*f[5]*nu[244]+2.5617376914899*rdxF[0]*f[40]*nu[112]+1.14564392373896*rdxF[0]*f[8]*nu[112]+1.677050983124842*rdxF[0]*f[16]*nu[92]+0.75*rdxF[0]*f[8]*nu[82]; 
  out[17] += 2.5617376914899*rdxF[2]*f[46]*nu[274]+1.14564392373896*rdxF[2]*f[6]*nu[274]+1.677050983124842*rdxF[2]*f[17]*nu[254]+0.75*rdxF[2]*f[6]*nu[244]+2.5617376914899*rdxF[1]*f[43]*nu[193]+1.14564392373896*rdxF[1]*f[8]*nu[193]+1.677050983124842*rdxF[1]*f[17]*nu[173]+0.75*rdxF[1]*f[8]*nu[163]; 
  out[18] += 2.5617376914899*rdxF[2]*f[47]*nu[274]+1.14564392373896*rdxF[2]*f[7]*nu[274]+1.677050983124842*rdxF[2]*f[18]*nu[254]+0.75*rdxF[2]*f[7]*nu[244]+2.5617376914899*rdxF[1]*f[44]*nu[193]+1.14564392373896*rdxF[1]*f[9]*nu[193]+1.677050983124842*rdxF[1]*f[18]*nu[173]+0.75*rdxF[1]*f[9]*nu[163]+2.5617376914899*rdxF[0]*f[42]*nu[112]+1.14564392373896*rdxF[0]*f[10]*nu[112]+1.677050983124842*rdxF[0]*f[18]*nu[92]+0.75*rdxF[0]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[0]*f[11]*nu[112]+1.677050983124842*rdxF[0]*f[19]*nu[92]+0.75*rdxF[0]*f[11]*nu[82]; 
  out[20] += 6.708203932499369*rdxF[0]*f[49]*nu[112]+7.6852130744697*rdxF[0]*f[5]*nu[112]+5.031152949374527*rdxF[0]*f[20]*nu[92]+3.75*rdxF[0]*f[1]*nu[92]+3.354101966249684*rdxF[0]*f[5]*nu[82]+1.677050983124842*rdxF[0]*f[1]*nu[80]; 
  out[21] += 1.14564392373896*rdxF[1]*f[11]*nu[193]+1.677050983124842*rdxF[1]*f[21]*nu[173]+0.75*rdxF[1]*f[11]*nu[163]; 
  out[22] += 1.14564392373896*rdxF[1]*f[12]*nu[193]+1.677050983124842*rdxF[1]*f[22]*nu[173]+0.75*rdxF[1]*f[12]*nu[163]+6.708203932499369*rdxF[0]*f[51]*nu[112]+7.6852130744697*rdxF[0]*f[7]*nu[112]+5.031152949374527*rdxF[0]*f[22]*nu[92]+3.75*rdxF[0]*f[3]*nu[92]+3.354101966249684*rdxF[0]*f[7]*nu[82]+1.677050983124842*rdxF[0]*f[3]*nu[80]; 
  out[23] += 6.708203932499369*rdxF[1]*f[52]*nu[193]+7.6852130744697*rdxF[1]*f[6]*nu[193]+5.031152949374527*rdxF[1]*f[23]*nu[173]+3.75*f[1]*rdxF[1]*nu[173]+3.354101966249684*rdxF[1]*f[6]*nu[163]+1.677050983124842*f[1]*rdxF[1]*nu[160]; 
  out[24] += 6.708203932499369*rdxF[1]*f[53]*nu[193]+7.6852130744697*rdxF[1]*f[7]*nu[193]+5.031152949374527*rdxF[1]*f[24]*nu[173]+3.75*rdxF[1]*f[2]*nu[173]+3.354101966249684*rdxF[1]*f[7]*nu[163]+1.677050983124842*rdxF[1]*f[2]*nu[160]+1.14564392373896*rdxF[0]*f[13]*nu[112]+1.677050983124842*rdxF[0]*f[24]*nu[92]+0.75*rdxF[0]*f[13]*nu[82]; 
  out[25] += 1.14564392373896*rdxF[2]*f[11]*nu[274]+1.677050983124842*rdxF[2]*f[25]*nu[254]+0.75*rdxF[2]*f[11]*nu[244]; 
  out[26] += 1.14564392373896*rdxF[2]*f[12]*nu[274]+1.677050983124842*rdxF[2]*f[26]*nu[254]+0.75*rdxF[2]*f[12]*nu[244]+6.708203932499369*rdxF[0]*f[55]*nu[112]+7.6852130744697*rdxF[0]*f[9]*nu[112]+5.031152949374527*rdxF[0]*f[26]*nu[92]+3.75*rdxF[0]*f[4]*nu[92]+3.354101966249684*rdxF[0]*f[9]*nu[82]+1.677050983124842*rdxF[0]*f[4]*nu[80]; 
  out[27] += 1.14564392373896*rdxF[2]*f[13]*nu[274]+1.677050983124842*rdxF[2]*f[27]*nu[254]+0.75*rdxF[2]*f[13]*nu[244]+6.708203932499369*rdxF[1]*f[56]*nu[193]+7.6852130744697*rdxF[1]*f[10]*nu[193]+5.031152949374527*rdxF[1]*f[27]*nu[173]+3.75*rdxF[1]*f[4]*nu[173]+3.354101966249684*rdxF[1]*f[10]*nu[163]+1.677050983124842*rdxF[1]*f[4]*nu[160]; 
  out[28] += 6.708203932499369*rdxF[2]*f[57]*nu[274]+7.6852130744697*rdxF[2]*f[8]*nu[274]+5.031152949374527*rdxF[2]*f[28]*nu[254]+3.75*f[1]*rdxF[2]*nu[254]+3.354101966249684*rdxF[2]*f[8]*nu[244]+1.677050983124842*f[1]*rdxF[2]*nu[240]; 
  out[29] += 6.708203932499369*rdxF[2]*f[58]*nu[274]+7.6852130744697*rdxF[2]*f[9]*nu[274]+5.031152949374527*rdxF[2]*f[29]*nu[254]+3.75*f[2]*rdxF[2]*nu[254]+3.354101966249684*rdxF[2]*f[9]*nu[244]+1.677050983124842*f[2]*rdxF[2]*nu[240]+1.14564392373896*rdxF[0]*f[14]*nu[112]+1.677050983124842*rdxF[0]*f[29]*nu[92]+0.75*rdxF[0]*f[14]*nu[82]; 
  out[30] += 6.708203932499369*rdxF[2]*f[59]*nu[274]+7.6852130744697*rdxF[2]*f[10]*nu[274]+5.031152949374527*rdxF[2]*f[30]*nu[254]+3.75*rdxF[2]*f[3]*nu[254]+3.354101966249684*rdxF[2]*f[10]*nu[244]+1.677050983124842*rdxF[2]*f[3]*nu[240]+1.14564392373896*rdxF[1]*f[14]*nu[193]+1.677050983124842*rdxF[1]*f[30]*nu[173]+0.75*rdxF[1]*f[14]*nu[163]; 
  out[32] += 18.44756081437327*rdxF[0]*f[12]*nu[112]+10.5*f[0]*rdxF[0]*nu[112]+10.06230589874905*rdxF[0]*f[32]*nu[92]+12.8086884574495*rdxF[0]*f[2]*nu[92]+7.685213074469699*rdxF[0]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[0]*nu[82]+5.7282196186948*rdxF[0]*f[2]*nu[80]; 
  out[33] += 18.44756081437327*rdxF[1]*f[13]*nu[193]+10.5*f[0]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[33]*nu[173]+12.8086884574495*rdxF[1]*f[3]*nu[173]+7.685213074469699*rdxF[1]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[3]*nu[160]; 
  out[34] += 18.44756081437327*rdxF[2]*f[14]*nu[274]+10.5*f[0]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[34]*nu[254]+12.8086884574495*rdxF[2]*f[4]*nu[254]+7.685213074469699*rdxF[2]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[2]*f[63]*nu[274]+1.14564392373896*rdxF[2]*f[15]*nu[274]+1.677050983124842*rdxF[2]*f[35]*nu[254]+0.75*rdxF[2]*f[15]*nu[244]+2.5617376914899*rdxF[1]*f[62]*nu[193]+1.14564392373896*rdxF[1]*f[16]*nu[193]+1.677050983124842*rdxF[1]*f[35]*nu[173]+0.75*rdxF[1]*f[16]*nu[163]+2.5617376914899*rdxF[0]*f[61]*nu[112]+1.14564392373896*rdxF[0]*f[17]*nu[112]+1.677050983124842*rdxF[0]*f[35]*nu[92]+0.75*rdxF[0]*f[17]*nu[82]; 
  out[36] += 1.14564392373896*rdxF[1]*f[19]*nu[193]+1.677050983124842*rdxF[1]*f[36]*nu[173]+0.75*rdxF[1]*f[19]*nu[163]+1.14564392373896*rdxF[0]*f[21]*nu[112]+1.677050983124842*rdxF[0]*f[36]*nu[92]+0.75*rdxF[0]*f[21]*nu[82]; 
  out[37] += 1.14564392373896*rdxF[1]*f[20]*nu[193]+1.677050983124842*rdxF[1]*f[37]*nu[173]+0.75*rdxF[1]*f[20]*nu[163]+6.708203932499369*rdxF[0]*f[65]*nu[112]+7.685213074469699*rdxF[0]*f[15]*nu[112]+5.031152949374527*rdxF[0]*f[37]*nu[92]+3.75*rdxF[0]*f[6]*nu[92]+3.354101966249685*rdxF[0]*f[15]*nu[82]+1.677050983124842*rdxF[0]*f[6]*nu[80]; 
  out[38] += 6.708203932499369*rdxF[1]*f[66]*nu[193]+7.685213074469699*rdxF[1]*f[15]*nu[193]+5.031152949374527*rdxF[1]*f[38]*nu[173]+3.75*rdxF[1]*f[5]*nu[173]+3.354101966249685*rdxF[1]*f[15]*nu[163]+1.677050983124842*rdxF[1]*f[5]*nu[160]+1.14564392373896*rdxF[0]*f[23]*nu[112]+1.677050983124842*rdxF[0]*f[38]*nu[92]+0.75*rdxF[0]*f[23]*nu[82]; 
  out[39] += 1.14564392373896*rdxF[2]*f[19]*nu[274]+1.677050983124842*rdxF[2]*f[39]*nu[254]+0.75*rdxF[2]*f[19]*nu[244]+1.14564392373896*rdxF[0]*f[25]*nu[112]+1.677050983124842*rdxF[0]*f[39]*nu[92]+0.75*rdxF[0]*f[25]*nu[82]; 
  out[40] += 1.14564392373896*rdxF[2]*f[20]*nu[274]+1.677050983124842*rdxF[2]*f[40]*nu[254]+0.75*rdxF[2]*f[20]*nu[244]+6.708203932499369*rdxF[0]*f[68]*nu[112]+7.685213074469699*rdxF[0]*f[16]*nu[112]+5.031152949374527*rdxF[0]*f[40]*nu[92]+3.75*rdxF[0]*f[8]*nu[92]+3.354101966249685*rdxF[0]*f[16]*nu[82]+1.677050983124842*rdxF[0]*f[8]*nu[80]; 
  out[41] += 1.14564392373896*rdxF[2]*f[21]*nu[274]+1.677050983124842*rdxF[2]*f[41]*nu[254]+0.75*rdxF[2]*f[21]*nu[244]+1.14564392373896*rdxF[1]*f[25]*nu[193]+1.677050983124842*rdxF[1]*f[41]*nu[173]+0.75*rdxF[1]*f[25]*nu[163]; 
  out[42] += 1.14564392373896*rdxF[2]*f[22]*nu[274]+1.677050983124842*rdxF[2]*f[42]*nu[254]+0.75*rdxF[2]*f[22]*nu[244]+1.14564392373896*rdxF[1]*f[26]*nu[193]+1.677050983124842*rdxF[1]*f[42]*nu[173]+0.75*rdxF[1]*f[26]*nu[163]+6.708203932499369*rdxF[0]*f[70]*nu[112]+7.685213074469699*rdxF[0]*f[18]*nu[112]+5.031152949374527*rdxF[0]*f[42]*nu[92]+3.75*rdxF[0]*f[10]*nu[92]+3.354101966249685*rdxF[0]*f[18]*nu[82]+1.677050983124842*rdxF[0]*f[10]*nu[80]; 
  out[43] += 1.14564392373896*rdxF[2]*f[23]*nu[274]+1.677050983124842*rdxF[2]*f[43]*nu[254]+0.75*rdxF[2]*f[23]*nu[244]+6.708203932499369*rdxF[1]*f[71]*nu[193]+7.685213074469699*rdxF[1]*f[17]*nu[193]+5.031152949374527*rdxF[1]*f[43]*nu[173]+3.75*rdxF[1]*f[8]*nu[173]+3.354101966249685*rdxF[1]*f[17]*nu[163]+1.677050983124842*rdxF[1]*f[8]*nu[160]; 
  out[44] += 1.14564392373896*rdxF[2]*f[24]*nu[274]+1.677050983124842*rdxF[2]*f[44]*nu[254]+0.75*rdxF[2]*f[24]*nu[244]+6.708203932499369*rdxF[1]*f[72]*nu[193]+7.685213074469699*rdxF[1]*f[18]*nu[193]+5.031152949374527*rdxF[1]*f[44]*nu[173]+3.75*rdxF[1]*f[9]*nu[173]+3.354101966249685*rdxF[1]*f[18]*nu[163]+1.677050983124842*rdxF[1]*f[9]*nu[160]+1.14564392373896*rdxF[0]*f[27]*nu[112]+1.677050983124842*rdxF[0]*f[44]*nu[92]+0.75*rdxF[0]*f[27]*nu[82]; 
  out[45] += 6.708203932499369*rdxF[2]*f[73]*nu[274]+7.685213074469699*rdxF[2]*f[16]*nu[274]+5.031152949374527*rdxF[2]*f[45]*nu[254]+3.75*rdxF[2]*f[5]*nu[254]+3.354101966249685*rdxF[2]*f[16]*nu[244]+1.677050983124842*rdxF[2]*f[5]*nu[240]+1.14564392373896*rdxF[0]*f[28]*nu[112]+1.677050983124842*rdxF[0]*f[45]*nu[92]+0.75*rdxF[0]*f[28]*nu[82]; 
  out[46] += 6.708203932499369*rdxF[2]*f[74]*nu[274]+7.685213074469699*rdxF[2]*f[17]*nu[274]+5.031152949374527*rdxF[2]*f[46]*nu[254]+3.75*rdxF[2]*f[6]*nu[254]+3.354101966249685*rdxF[2]*f[17]*nu[244]+1.677050983124842*rdxF[2]*f[6]*nu[240]+1.14564392373896*rdxF[1]*f[28]*nu[193]+1.677050983124842*rdxF[1]*f[46]*nu[173]+0.75*rdxF[1]*f[28]*nu[163]; 
  out[47] += 6.708203932499369*rdxF[2]*f[75]*nu[274]+7.685213074469699*rdxF[2]*f[18]*nu[274]+5.031152949374527*rdxF[2]*f[47]*nu[254]+3.75*rdxF[2]*f[7]*nu[254]+3.354101966249685*rdxF[2]*f[18]*nu[244]+1.677050983124842*rdxF[2]*f[7]*nu[240]+1.14564392373896*rdxF[1]*f[29]*nu[193]+1.677050983124842*rdxF[1]*f[47]*nu[173]+0.75*rdxF[1]*f[29]*nu[163]+1.14564392373896*rdxF[0]*f[30]*nu[112]+1.677050983124842*rdxF[0]*f[47]*nu[92]+0.75*rdxF[0]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[0]*f[31]*nu[112]+1.677050983124842*rdxF[0]*f[48]*nu[92]+0.7499999999999999*rdxF[0]*f[31]*nu[82]; 
  out[49] += 18.44756081437327*rdxF[0]*f[20]*nu[112]+10.5*rdxF[0]*f[1]*nu[112]+10.06230589874905*rdxF[0]*f[49]*nu[92]+12.8086884574495*rdxF[0]*f[5]*nu[92]+7.685213074469698*rdxF[0]*f[20]*nu[82]+6.873863542433758*rdxF[0]*f[1]*nu[82]+5.7282196186948*rdxF[0]*f[5]*nu[80]; 
  out[50] += 1.14564392373896*rdxF[1]*f[31]*nu[193]+1.677050983124842*rdxF[1]*f[50]*nu[173]+0.7499999999999999*rdxF[1]*f[31]*nu[163]; 
  out[51] += 1.14564392373896*rdxF[1]*f[32]*nu[193]+1.677050983124842*rdxF[1]*f[51]*nu[173]+0.7499999999999999*rdxF[1]*f[32]*nu[163]+18.44756081437327*rdxF[0]*f[22]*nu[112]+10.5*rdxF[0]*f[3]*nu[112]+10.06230589874905*rdxF[0]*f[51]*nu[92]+12.8086884574495*rdxF[0]*f[7]*nu[92]+7.685213074469698*rdxF[0]*f[22]*nu[82]+6.873863542433758*rdxF[0]*f[3]*nu[82]+5.7282196186948*rdxF[0]*f[7]*nu[80]; 
  out[52] += 18.44756081437327*rdxF[1]*f[23]*nu[193]+10.5*f[1]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[52]*nu[173]+12.8086884574495*rdxF[1]*f[6]*nu[173]+7.685213074469698*rdxF[1]*f[23]*nu[163]+6.873863542433758*f[1]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[6]*nu[160]; 
  out[53] += 18.44756081437327*rdxF[1]*f[24]*nu[193]+10.5*rdxF[1]*f[2]*nu[193]+10.06230589874905*rdxF[1]*f[53]*nu[173]+12.8086884574495*rdxF[1]*f[7]*nu[173]+7.685213074469698*rdxF[1]*f[24]*nu[163]+6.873863542433758*rdxF[1]*f[2]*nu[163]+5.7282196186948*rdxF[1]*f[7]*nu[160]+1.14564392373896*rdxF[0]*f[33]*nu[112]+1.677050983124842*rdxF[0]*f[53]*nu[92]+0.7499999999999999*rdxF[0]*f[33]*nu[82]; 
  out[54] += 1.14564392373896*rdxF[2]*f[31]*nu[274]+1.677050983124842*rdxF[2]*f[54]*nu[254]+0.7499999999999999*rdxF[2]*f[31]*nu[244]; 
  out[55] += 1.14564392373896*rdxF[2]*f[32]*nu[274]+1.677050983124842*rdxF[2]*f[55]*nu[254]+0.7499999999999999*rdxF[2]*f[32]*nu[244]+18.44756081437327*rdxF[0]*f[26]*nu[112]+10.5*rdxF[0]*f[4]*nu[112]+10.06230589874905*rdxF[0]*f[55]*nu[92]+12.8086884574495*rdxF[0]*f[9]*nu[92]+7.685213074469698*rdxF[0]*f[26]*nu[82]+6.873863542433758*rdxF[0]*f[4]*nu[82]+5.7282196186948*rdxF[0]*f[9]*nu[80]; 
  out[56] += 1.14564392373896*rdxF[2]*f[33]*nu[274]+1.677050983124842*rdxF[2]*f[56]*nu[254]+0.7499999999999999*rdxF[2]*f[33]*nu[244]+18.44756081437327*rdxF[1]*f[27]*nu[193]+10.5*rdxF[1]*f[4]*nu[193]+10.06230589874905*rdxF[1]*f[56]*nu[173]+12.8086884574495*rdxF[1]*f[10]*nu[173]+7.685213074469698*rdxF[1]*f[27]*nu[163]+6.873863542433758*rdxF[1]*f[4]*nu[163]+5.7282196186948*rdxF[1]*f[10]*nu[160]; 
  out[57] += 18.44756081437327*rdxF[2]*f[28]*nu[274]+10.5*f[1]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[57]*nu[254]+12.8086884574495*rdxF[2]*f[8]*nu[254]+7.685213074469698*rdxF[2]*f[28]*nu[244]+6.873863542433758*f[1]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[8]*nu[240]; 
  out[58] += 18.44756081437327*rdxF[2]*f[29]*nu[274]+10.5*f[2]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[58]*nu[254]+12.8086884574495*rdxF[2]*f[9]*nu[254]+7.685213074469698*rdxF[2]*f[29]*nu[244]+6.873863542433758*f[2]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[9]*nu[240]+1.14564392373896*rdxF[0]*f[34]*nu[112]+1.677050983124842*rdxF[0]*f[58]*nu[92]+0.7499999999999999*rdxF[0]*f[34]*nu[82]; 
  out[59] += 18.44756081437327*rdxF[2]*f[30]*nu[274]+10.5*rdxF[2]*f[3]*nu[274]+10.06230589874905*rdxF[2]*f[59]*nu[254]+12.8086884574495*rdxF[2]*f[10]*nu[254]+7.685213074469698*rdxF[2]*f[30]*nu[244]+6.873863542433758*rdxF[2]*f[3]*nu[244]+5.7282196186948*rdxF[2]*f[10]*nu[240]+1.14564392373896*rdxF[1]*f[34]*nu[193]+1.677050983124842*rdxF[1]*f[59]*nu[173]+0.7499999999999999*rdxF[1]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[2]*f[36]*nu[274]+1.677050983124842*rdxF[2]*f[60]*nu[254]+0.75*rdxF[2]*f[36]*nu[244]+1.14564392373896*rdxF[1]*f[39]*nu[193]+1.677050983124842*rdxF[1]*f[60]*nu[173]+0.75*rdxF[1]*f[39]*nu[163]+1.14564392373896*rdxF[0]*f[41]*nu[112]+1.677050983124842*rdxF[0]*f[60]*nu[92]+0.75*rdxF[0]*f[41]*nu[82]; 
  out[61] += 1.14564392373896*rdxF[2]*f[37]*nu[274]+1.677050983124842*rdxF[2]*f[61]*nu[254]+0.75*rdxF[2]*f[37]*nu[244]+1.14564392373896*rdxF[1]*f[40]*nu[193]+1.677050983124842*rdxF[1]*f[61]*nu[173]+0.75*rdxF[1]*f[40]*nu[163]+6.708203932499369*rdxF[0]*f[77]*nu[112]+7.6852130744697*rdxF[0]*f[35]*nu[112]+5.031152949374527*rdxF[0]*f[61]*nu[92]+3.75*rdxF[0]*f[17]*nu[92]+3.354101966249684*rdxF[0]*f[35]*nu[82]+1.677050983124842*rdxF[0]*f[17]*nu[80]; 
  out[62] += 1.14564392373896*rdxF[2]*f[38]*nu[274]+1.677050983124842*rdxF[2]*f[62]*nu[254]+0.75*rdxF[2]*f[38]*nu[244]+6.708203932499369*rdxF[1]*f[78]*nu[193]+7.6852130744697*rdxF[1]*f[35]*nu[193]+5.031152949374527*rdxF[1]*f[62]*nu[173]+3.75*rdxF[1]*f[16]*nu[173]+3.354101966249684*rdxF[1]*f[35]*nu[163]+1.677050983124842*rdxF[1]*f[16]*nu[160]+1.14564392373896*rdxF[0]*f[43]*nu[112]+1.677050983124842*rdxF[0]*f[62]*nu[92]+0.75*rdxF[0]*f[43]*nu[82]; 
  out[63] += 6.708203932499369*rdxF[2]*f[79]*nu[274]+7.6852130744697*rdxF[2]*f[35]*nu[274]+5.031152949374527*rdxF[2]*f[63]*nu[254]+3.75*rdxF[2]*f[15]*nu[254]+3.354101966249684*rdxF[2]*f[35]*nu[244]+1.677050983124842*rdxF[2]*f[15]*nu[240]+1.14564392373896*rdxF[1]*f[45]*nu[193]+1.677050983124842*rdxF[1]*f[63]*nu[173]+0.75*rdxF[1]*f[45]*nu[163]+1.14564392373896*rdxF[0]*f[46]*nu[112]+1.677050983124842*rdxF[0]*f[63]*nu[92]+0.75*rdxF[0]*f[46]*nu[82]; 
  out[64] += 1.14564392373896*rdxF[1]*f[48]*nu[193]+1.677050983124842*rdxF[1]*f[64]*nu[173]+0.7499999999999999*rdxF[1]*f[48]*nu[163]+1.14564392373896*rdxF[0]*f[50]*nu[112]+1.677050983124842*rdxF[0]*f[64]*nu[92]+0.7499999999999999*rdxF[0]*f[50]*nu[82]; 
  out[65] += 1.14564392373896*rdxF[1]*f[49]*nu[193]+1.677050983124842*rdxF[1]*f[65]*nu[173]+0.7499999999999999*rdxF[1]*f[49]*nu[163]+18.44756081437327*rdxF[0]*f[37]*nu[112]+10.5*rdxF[0]*f[6]*nu[112]+10.06230589874905*rdxF[0]*f[65]*nu[92]+12.8086884574495*rdxF[0]*f[15]*nu[92]+7.685213074469699*rdxF[0]*f[37]*nu[82]+6.873863542433759*rdxF[0]*f[6]*nu[82]+5.7282196186948*rdxF[0]*f[15]*nu[80]; 
  out[66] += 18.44756081437327*rdxF[1]*f[38]*nu[193]+10.5*rdxF[1]*f[5]*nu[193]+10.06230589874905*rdxF[1]*f[66]*nu[173]+12.8086884574495*rdxF[1]*f[15]*nu[173]+7.685213074469699*rdxF[1]*f[38]*nu[163]+6.873863542433759*rdxF[1]*f[5]*nu[163]+5.7282196186948*rdxF[1]*f[15]*nu[160]+1.14564392373896*rdxF[0]*f[52]*nu[112]+1.677050983124842*rdxF[0]*f[66]*nu[92]+0.7499999999999999*rdxF[0]*f[52]*nu[82]; 
  out[67] += 1.14564392373896*rdxF[2]*f[48]*nu[274]+1.677050983124842*rdxF[2]*f[67]*nu[254]+0.7499999999999999*rdxF[2]*f[48]*nu[244]+1.14564392373896*rdxF[0]*f[54]*nu[112]+1.677050983124842*rdxF[0]*f[67]*nu[92]+0.7499999999999999*rdxF[0]*f[54]*nu[82]; 
  out[68] += 1.14564392373896*rdxF[2]*f[49]*nu[274]+1.677050983124842*rdxF[2]*f[68]*nu[254]+0.7499999999999999*rdxF[2]*f[49]*nu[244]+18.44756081437327*rdxF[0]*f[40]*nu[112]+10.5*rdxF[0]*f[8]*nu[112]+10.06230589874905*rdxF[0]*f[68]*nu[92]+12.8086884574495*rdxF[0]*f[16]*nu[92]+7.685213074469699*rdxF[0]*f[40]*nu[82]+6.873863542433759*rdxF[0]*f[8]*nu[82]+5.7282196186948*rdxF[0]*f[16]*nu[80]; 
  out[69] += 1.14564392373896*rdxF[2]*f[50]*nu[274]+1.677050983124842*rdxF[2]*f[69]*nu[254]+0.7499999999999999*rdxF[2]*f[50]*nu[244]+1.14564392373896*rdxF[1]*f[54]*nu[193]+1.677050983124842*rdxF[1]*f[69]*nu[173]+0.7499999999999999*rdxF[1]*f[54]*nu[163]; 
  out[70] += 1.14564392373896*rdxF[2]*f[51]*nu[274]+1.677050983124842*rdxF[2]*f[70]*nu[254]+0.7499999999999999*rdxF[2]*f[51]*nu[244]+1.14564392373896*rdxF[1]*f[55]*nu[193]+1.677050983124842*rdxF[1]*f[70]*nu[173]+0.7499999999999999*rdxF[1]*f[55]*nu[163]+18.44756081437327*rdxF[0]*f[42]*nu[112]+10.5*rdxF[0]*f[10]*nu[112]+10.06230589874905*rdxF[0]*f[70]*nu[92]+12.8086884574495*rdxF[0]*f[18]*nu[92]+7.685213074469699*rdxF[0]*f[42]*nu[82]+6.873863542433759*rdxF[0]*f[10]*nu[82]+5.7282196186948*rdxF[0]*f[18]*nu[80]; 
  out[71] += 1.14564392373896*rdxF[2]*f[52]*nu[274]+1.677050983124842*rdxF[2]*f[71]*nu[254]+0.7499999999999999*rdxF[2]*f[52]*nu[244]+18.44756081437327*rdxF[1]*f[43]*nu[193]+10.5*rdxF[1]*f[8]*nu[193]+10.06230589874905*rdxF[1]*f[71]*nu[173]+12.8086884574495*rdxF[1]*f[17]*nu[173]+7.685213074469699*rdxF[1]*f[43]*nu[163]+6.873863542433759*rdxF[1]*f[8]*nu[163]+5.7282196186948*rdxF[1]*f[17]*nu[160]; 
  out[72] += 1.14564392373896*rdxF[2]*f[53]*nu[274]+1.677050983124842*rdxF[2]*f[72]*nu[254]+0.7499999999999999*rdxF[2]*f[53]*nu[244]+18.44756081437327*rdxF[1]*f[44]*nu[193]+10.5*rdxF[1]*f[9]*nu[193]+10.06230589874905*rdxF[1]*f[72]*nu[173]+12.8086884574495*rdxF[1]*f[18]*nu[173]+7.685213074469699*rdxF[1]*f[44]*nu[163]+6.873863542433759*rdxF[1]*f[9]*nu[163]+5.7282196186948*rdxF[1]*f[18]*nu[160]+1.14564392373896*rdxF[0]*f[56]*nu[112]+1.677050983124842*rdxF[0]*f[72]*nu[92]+0.7499999999999999*rdxF[0]*f[56]*nu[82]; 
  out[73] += 18.44756081437327*rdxF[2]*f[45]*nu[274]+10.5*rdxF[2]*f[5]*nu[274]+10.06230589874905*rdxF[2]*f[73]*nu[254]+12.8086884574495*rdxF[2]*f[16]*nu[254]+7.685213074469699*rdxF[2]*f[45]*nu[244]+6.873863542433759*rdxF[2]*f[5]*nu[244]+5.7282196186948*rdxF[2]*f[16]*nu[240]+1.14564392373896*rdxF[0]*f[57]*nu[112]+1.677050983124842*rdxF[0]*f[73]*nu[92]+0.7499999999999999*rdxF[0]*f[57]*nu[82]; 
  out[74] += 18.44756081437327*rdxF[2]*f[46]*nu[274]+10.5*rdxF[2]*f[6]*nu[274]+10.06230589874905*rdxF[2]*f[74]*nu[254]+12.8086884574495*rdxF[2]*f[17]*nu[254]+7.685213074469699*rdxF[2]*f[46]*nu[244]+6.873863542433759*rdxF[2]*f[6]*nu[244]+5.7282196186948*rdxF[2]*f[17]*nu[240]+1.14564392373896*rdxF[1]*f[57]*nu[193]+1.677050983124842*rdxF[1]*f[74]*nu[173]+0.7499999999999999*rdxF[1]*f[57]*nu[163]; 
  out[75] += 18.44756081437327*rdxF[2]*f[47]*nu[274]+10.5*rdxF[2]*f[7]*nu[274]+10.06230589874905*rdxF[2]*f[75]*nu[254]+12.8086884574495*rdxF[2]*f[18]*nu[254]+7.685213074469699*rdxF[2]*f[47]*nu[244]+6.873863542433759*rdxF[2]*f[7]*nu[244]+5.7282196186948*rdxF[2]*f[18]*nu[240]+1.14564392373896*rdxF[1]*f[58]*nu[193]+1.677050983124842*rdxF[1]*f[75]*nu[173]+0.7499999999999999*rdxF[1]*f[58]*nu[163]+1.14564392373896*rdxF[0]*f[59]*nu[112]+1.677050983124842*rdxF[0]*f[75]*nu[92]+0.7499999999999999*rdxF[0]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[2]*f[64]*nu[274]+1.677050983124842*rdxF[2]*f[76]*nu[254]+0.7499999999999999*rdxF[2]*f[64]*nu[244]+1.14564392373896*rdxF[1]*f[67]*nu[193]+1.677050983124842*rdxF[1]*f[76]*nu[173]+0.7499999999999999*rdxF[1]*f[67]*nu[163]+1.14564392373896*rdxF[0]*f[69]*nu[112]+1.677050983124842*rdxF[0]*f[76]*nu[92]+0.7499999999999999*rdxF[0]*f[69]*nu[82]; 
  out[77] += 1.14564392373896*rdxF[2]*f[65]*nu[274]+1.677050983124842*rdxF[2]*f[77]*nu[254]+0.7499999999999999*rdxF[2]*f[65]*nu[244]+1.14564392373896*rdxF[1]*f[68]*nu[193]+1.677050983124842*rdxF[1]*f[77]*nu[173]+0.7499999999999999*rdxF[1]*f[68]*nu[163]+18.44756081437327*rdxF[0]*f[61]*nu[112]+10.5*rdxF[0]*f[17]*nu[112]+10.06230589874905*rdxF[0]*f[77]*nu[92]+12.8086884574495*rdxF[0]*f[35]*nu[92]+7.685213074469698*rdxF[0]*f[61]*nu[82]+6.873863542433758*rdxF[0]*f[17]*nu[82]+5.7282196186948*rdxF[0]*f[35]*nu[80]; 
  out[78] += 1.14564392373896*rdxF[2]*f[66]*nu[274]+1.677050983124842*rdxF[2]*f[78]*nu[254]+0.7499999999999999*rdxF[2]*f[66]*nu[244]+18.44756081437327*rdxF[1]*f[62]*nu[193]+10.5*rdxF[1]*f[16]*nu[193]+10.06230589874905*rdxF[1]*f[78]*nu[173]+12.8086884574495*rdxF[1]*f[35]*nu[173]+7.685213074469698*rdxF[1]*f[62]*nu[163]+6.873863542433758*rdxF[1]*f[16]*nu[163]+5.7282196186948*rdxF[1]*f[35]*nu[160]+1.14564392373896*rdxF[0]*f[71]*nu[112]+1.677050983124842*rdxF[0]*f[78]*nu[92]+0.7499999999999999*rdxF[0]*f[71]*nu[82]; 
  out[79] += 18.44756081437327*rdxF[2]*f[63]*nu[274]+10.5*rdxF[2]*f[15]*nu[274]+10.06230589874905*rdxF[2]*f[79]*nu[254]+12.8086884574495*rdxF[2]*f[35]*nu[254]+7.685213074469698*rdxF[2]*f[63]*nu[244]+6.873863542433758*rdxF[2]*f[15]*nu[244]+5.7282196186948*rdxF[2]*f[35]*nu[240]+1.14564392373896*rdxF[1]*f[73]*nu[193]+1.677050983124842*rdxF[1]*f[79]*nu[173]+0.7499999999999999*rdxF[1]*f[73]*nu[163]+1.14564392373896*rdxF[0]*f[74]*nu[112]+1.677050983124842*rdxF[0]*f[79]*nu[92]+0.7499999999999999*rdxF[0]*f[74]*nu[82]; 

  return (rdxF[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+rdxF[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+rdxF[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs1234(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[4]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 
  rdxF[3] = 4.0/(dx[3]*dx[3]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 2.5617376914899*rdxF[1]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[2]*nu[92]+0.75*f[0]*rdxF[1]*nu[82]; 
  out[3] += 2.5617376914899*rdxF[2]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[2]*nu[193]+1.677050983124842*rdxF[2]*f[3]*nu[173]+0.75*f[0]*rdxF[2]*nu[163]; 
  out[4] += 2.5617376914899*rdxF[3]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[3]*nu[274]+1.677050983124842*rdxF[3]*f[4]*nu[254]+0.75*f[0]*rdxF[3]*nu[244]; 
  out[5] += 2.5617376914899*rdxF[1]*f[20]*nu[112]+1.14564392373896*f[1]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[5]*nu[92]+0.75*f[1]*rdxF[1]*nu[82]+2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[2]*f[23]*nu[193]+1.14564392373896*f[1]*rdxF[2]*nu[193]+1.677050983124842*rdxF[2]*f[6]*nu[173]+0.75*f[1]*rdxF[2]*nu[163]+2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 2.5617376914899*rdxF[2]*f[24]*nu[193]+1.14564392373896*f[2]*rdxF[2]*nu[193]+1.677050983124842*rdxF[2]*f[7]*nu[173]+0.75*f[2]*rdxF[2]*nu[163]+2.5617376914899*rdxF[1]*f[22]*nu[112]+1.14564392373896*rdxF[1]*f[3]*nu[112]+1.677050983124842*rdxF[1]*f[7]*nu[92]+0.75*rdxF[1]*f[3]*nu[82]; 
  out[8] += 2.5617376914899*rdxF[3]*f[28]*nu[274]+1.14564392373896*f[1]*rdxF[3]*nu[274]+1.677050983124842*rdxF[3]*f[8]*nu[254]+0.75*f[1]*rdxF[3]*nu[244]+2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 2.5617376914899*rdxF[3]*f[29]*nu[274]+1.14564392373896*f[2]*rdxF[3]*nu[274]+1.677050983124842*rdxF[3]*f[9]*nu[254]+0.75*f[2]*rdxF[3]*nu[244]+2.5617376914899*rdxF[1]*f[26]*nu[112]+1.14564392373896*rdxF[1]*f[4]*nu[112]+1.677050983124842*rdxF[1]*f[9]*nu[92]+0.75*rdxF[1]*f[4]*nu[82]; 
  out[10] += 2.5617376914899*rdxF[3]*f[30]*nu[274]+1.14564392373896*f[3]*rdxF[3]*nu[274]+1.677050983124842*rdxF[3]*f[10]*nu[254]+0.75*f[3]*rdxF[3]*nu[244]+2.5617376914899*rdxF[2]*f[27]*nu[193]+1.14564392373896*rdxF[2]*f[4]*nu[193]+1.677050983124842*rdxF[2]*f[10]*nu[173]+0.75*rdxF[2]*f[4]*nu[163]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 6.708203932499369*rdxF[1]*f[32]*nu[112]+7.685213074469699*rdxF[1]*f[2]*nu[112]+5.031152949374527*rdxF[1]*f[12]*nu[92]+3.75*f[0]*rdxF[1]*nu[92]+3.354101966249685*rdxF[1]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[1]*nu[80]; 
  out[13] += 6.708203932499369*rdxF[2]*f[33]*nu[193]+7.685213074469699*rdxF[2]*f[3]*nu[193]+5.031152949374527*rdxF[2]*f[13]*nu[173]+3.75*f[0]*rdxF[2]*nu[173]+3.354101966249685*rdxF[2]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[2]*nu[160]; 
  out[14] += 6.708203932499369*rdxF[3]*f[34]*nu[274]+7.685213074469699*rdxF[3]*f[4]*nu[274]+5.031152949374527*rdxF[3]*f[14]*nu[254]+3.75*f[0]*rdxF[3]*nu[254]+3.354101966249685*rdxF[3]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[3]*nu[240]; 
  out[15] += 2.5617376914899*rdxF[2]*f[38]*nu[193]+1.14564392373896*rdxF[2]*f[5]*nu[193]+1.677050983124842*rdxF[2]*f[15]*nu[173]+0.75*rdxF[2]*f[5]*nu[163]+2.5617376914899*rdxF[1]*f[37]*nu[112]+1.14564392373896*rdxF[1]*f[6]*nu[112]+1.677050983124842*rdxF[1]*f[15]*nu[92]+0.75*rdxF[1]*f[6]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[3]*f[45]*nu[274]+1.14564392373896*rdxF[3]*f[5]*nu[274]+1.677050983124842*rdxF[3]*f[16]*nu[254]+0.75*rdxF[3]*f[5]*nu[244]+2.5617376914899*rdxF[1]*f[40]*nu[112]+1.14564392373896*rdxF[1]*f[8]*nu[112]+1.677050983124842*rdxF[1]*f[16]*nu[92]+0.75*rdxF[1]*f[8]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[3]*f[46]*nu[274]+1.14564392373896*rdxF[3]*f[6]*nu[274]+1.677050983124842*rdxF[3]*f[17]*nu[254]+0.75*rdxF[3]*f[6]*nu[244]+2.5617376914899*rdxF[2]*f[43]*nu[193]+1.14564392373896*rdxF[2]*f[8]*nu[193]+1.677050983124842*rdxF[2]*f[17]*nu[173]+0.75*rdxF[2]*f[8]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 2.5617376914899*rdxF[3]*f[47]*nu[274]+1.14564392373896*rdxF[3]*f[7]*nu[274]+1.677050983124842*rdxF[3]*f[18]*nu[254]+0.75*rdxF[3]*f[7]*nu[244]+2.5617376914899*rdxF[2]*f[44]*nu[193]+1.14564392373896*rdxF[2]*f[9]*nu[193]+1.677050983124842*rdxF[2]*f[18]*nu[173]+0.75*rdxF[2]*f[9]*nu[163]+2.5617376914899*rdxF[1]*f[42]*nu[112]+1.14564392373896*rdxF[1]*f[10]*nu[112]+1.677050983124842*rdxF[1]*f[18]*nu[92]+0.75*rdxF[1]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[1]*f[11]*nu[112]+1.677050983124842*rdxF[1]*f[19]*nu[92]+0.75*rdxF[1]*f[11]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 6.708203932499369*rdxF[1]*f[49]*nu[112]+7.6852130744697*rdxF[1]*f[5]*nu[112]+5.031152949374527*rdxF[1]*f[20]*nu[92]+3.75*f[1]*rdxF[1]*nu[92]+3.354101966249684*rdxF[1]*f[5]*nu[82]+1.677050983124842*f[1]*rdxF[1]*nu[80]+1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.14564392373896*rdxF[2]*f[11]*nu[193]+1.677050983124842*rdxF[2]*f[21]*nu[173]+0.75*rdxF[2]*f[11]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.14564392373896*rdxF[2]*f[12]*nu[193]+1.677050983124842*rdxF[2]*f[22]*nu[173]+0.75*rdxF[2]*f[12]*nu[163]+6.708203932499369*rdxF[1]*f[51]*nu[112]+7.6852130744697*rdxF[1]*f[7]*nu[112]+5.031152949374527*rdxF[1]*f[22]*nu[92]+3.75*rdxF[1]*f[3]*nu[92]+3.354101966249684*rdxF[1]*f[7]*nu[82]+1.677050983124842*rdxF[1]*f[3]*nu[80]; 
  out[23] += 6.708203932499369*rdxF[2]*f[52]*nu[193]+7.6852130744697*rdxF[2]*f[6]*nu[193]+5.031152949374527*rdxF[2]*f[23]*nu[173]+3.75*f[1]*rdxF[2]*nu[173]+3.354101966249684*rdxF[2]*f[6]*nu[163]+1.677050983124842*f[1]*rdxF[2]*nu[160]+1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 6.708203932499369*rdxF[2]*f[53]*nu[193]+7.6852130744697*rdxF[2]*f[7]*nu[193]+5.031152949374527*rdxF[2]*f[24]*nu[173]+3.75*f[2]*rdxF[2]*nu[173]+3.354101966249684*rdxF[2]*f[7]*nu[163]+1.677050983124842*f[2]*rdxF[2]*nu[160]+1.14564392373896*rdxF[1]*f[13]*nu[112]+1.677050983124842*rdxF[1]*f[24]*nu[92]+0.75*rdxF[1]*f[13]*nu[82]; 
  out[25] += 1.14564392373896*rdxF[3]*f[11]*nu[274]+1.677050983124842*rdxF[3]*f[25]*nu[254]+0.75*rdxF[3]*f[11]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.14564392373896*rdxF[3]*f[12]*nu[274]+1.677050983124842*rdxF[3]*f[26]*nu[254]+0.75*rdxF[3]*f[12]*nu[244]+6.708203932499369*rdxF[1]*f[55]*nu[112]+7.6852130744697*rdxF[1]*f[9]*nu[112]+5.031152949374527*rdxF[1]*f[26]*nu[92]+3.75*rdxF[1]*f[4]*nu[92]+3.354101966249684*rdxF[1]*f[9]*nu[82]+1.677050983124842*rdxF[1]*f[4]*nu[80]; 
  out[27] += 1.14564392373896*rdxF[3]*f[13]*nu[274]+1.677050983124842*rdxF[3]*f[27]*nu[254]+0.75*rdxF[3]*f[13]*nu[244]+6.708203932499369*rdxF[2]*f[56]*nu[193]+7.6852130744697*rdxF[2]*f[10]*nu[193]+5.031152949374527*rdxF[2]*f[27]*nu[173]+3.75*rdxF[2]*f[4]*nu[173]+3.354101966249684*rdxF[2]*f[10]*nu[163]+1.677050983124842*rdxF[2]*f[4]*nu[160]; 
  out[28] += 6.708203932499369*rdxF[3]*f[57]*nu[274]+7.6852130744697*rdxF[3]*f[8]*nu[274]+5.031152949374527*rdxF[3]*f[28]*nu[254]+3.75*f[1]*rdxF[3]*nu[254]+3.354101966249684*rdxF[3]*f[8]*nu[244]+1.677050983124842*f[1]*rdxF[3]*nu[240]+1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 6.708203932499369*rdxF[3]*f[58]*nu[274]+7.6852130744697*rdxF[3]*f[9]*nu[274]+5.031152949374527*rdxF[3]*f[29]*nu[254]+3.75*f[2]*rdxF[3]*nu[254]+3.354101966249684*rdxF[3]*f[9]*nu[244]+1.677050983124842*f[2]*rdxF[3]*nu[240]+1.14564392373896*rdxF[1]*f[14]*nu[112]+1.677050983124842*rdxF[1]*f[29]*nu[92]+0.75*rdxF[1]*f[14]*nu[82]; 
  out[30] += 6.708203932499369*rdxF[3]*f[59]*nu[274]+7.6852130744697*rdxF[3]*f[10]*nu[274]+5.031152949374527*rdxF[3]*f[30]*nu[254]+3.75*f[3]*rdxF[3]*nu[254]+3.354101966249684*rdxF[3]*f[10]*nu[244]+1.677050983124842*f[3]*rdxF[3]*nu[240]+1.14564392373896*rdxF[2]*f[14]*nu[193]+1.677050983124842*rdxF[2]*f[30]*nu[173]+0.75*rdxF[2]*f[14]*nu[163]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[32] += 18.44756081437327*rdxF[1]*f[12]*nu[112]+10.5*f[0]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[32]*nu[92]+12.8086884574495*rdxF[1]*f[2]*nu[92]+7.685213074469699*rdxF[1]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[2]*nu[80]; 
  out[33] += 18.44756081437327*rdxF[2]*f[13]*nu[193]+10.5*f[0]*rdxF[2]*nu[193]+10.06230589874905*rdxF[2]*f[33]*nu[173]+12.8086884574495*rdxF[2]*f[3]*nu[173]+7.685213074469699*rdxF[2]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[2]*nu[163]+5.7282196186948*rdxF[2]*f[3]*nu[160]; 
  out[34] += 18.44756081437327*rdxF[3]*f[14]*nu[274]+10.5*f[0]*rdxF[3]*nu[274]+10.06230589874905*rdxF[3]*f[34]*nu[254]+12.8086884574495*rdxF[3]*f[4]*nu[254]+7.685213074469699*rdxF[3]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[3]*nu[244]+5.7282196186948*rdxF[3]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[3]*f[63]*nu[274]+1.14564392373896*rdxF[3]*f[15]*nu[274]+1.677050983124842*rdxF[3]*f[35]*nu[254]+0.75*rdxF[3]*f[15]*nu[244]+2.5617376914899*rdxF[2]*f[62]*nu[193]+1.14564392373896*rdxF[2]*f[16]*nu[193]+1.677050983124842*rdxF[2]*f[35]*nu[173]+0.75*rdxF[2]*f[16]*nu[163]+2.5617376914899*rdxF[1]*f[61]*nu[112]+1.14564392373896*rdxF[1]*f[17]*nu[112]+1.677050983124842*rdxF[1]*f[35]*nu[92]+0.75*rdxF[1]*f[17]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 1.14564392373896*rdxF[2]*f[19]*nu[193]+1.677050983124842*rdxF[2]*f[36]*nu[173]+0.75*rdxF[2]*f[19]*nu[163]+1.14564392373896*rdxF[1]*f[21]*nu[112]+1.677050983124842*rdxF[1]*f[36]*nu[92]+0.75*rdxF[1]*f[21]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 1.14564392373896*rdxF[2]*f[20]*nu[193]+1.677050983124842*rdxF[2]*f[37]*nu[173]+0.75*rdxF[2]*f[20]*nu[163]+6.708203932499369*rdxF[1]*f[65]*nu[112]+7.685213074469699*rdxF[1]*f[15]*nu[112]+5.031152949374527*rdxF[1]*f[37]*nu[92]+3.75*rdxF[1]*f[6]*nu[92]+3.354101966249685*rdxF[1]*f[15]*nu[82]+1.677050983124842*rdxF[1]*f[6]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 6.708203932499369*rdxF[2]*f[66]*nu[193]+7.685213074469699*rdxF[2]*f[15]*nu[193]+5.031152949374527*rdxF[2]*f[38]*nu[173]+3.75*rdxF[2]*f[5]*nu[173]+3.354101966249685*rdxF[2]*f[15]*nu[163]+1.677050983124842*rdxF[2]*f[5]*nu[160]+1.14564392373896*rdxF[1]*f[23]*nu[112]+1.677050983124842*rdxF[1]*f[38]*nu[92]+0.75*rdxF[1]*f[23]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 1.14564392373896*rdxF[3]*f[19]*nu[274]+1.677050983124842*rdxF[3]*f[39]*nu[254]+0.75*rdxF[3]*f[19]*nu[244]+1.14564392373896*rdxF[1]*f[25]*nu[112]+1.677050983124842*rdxF[1]*f[39]*nu[92]+0.75*rdxF[1]*f[25]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 1.14564392373896*rdxF[3]*f[20]*nu[274]+1.677050983124842*rdxF[3]*f[40]*nu[254]+0.75*rdxF[3]*f[20]*nu[244]+6.708203932499369*rdxF[1]*f[68]*nu[112]+7.685213074469699*rdxF[1]*f[16]*nu[112]+5.031152949374527*rdxF[1]*f[40]*nu[92]+3.75*rdxF[1]*f[8]*nu[92]+3.354101966249685*rdxF[1]*f[16]*nu[82]+1.677050983124842*rdxF[1]*f[8]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 1.14564392373896*rdxF[3]*f[21]*nu[274]+1.677050983124842*rdxF[3]*f[41]*nu[254]+0.75*rdxF[3]*f[21]*nu[244]+1.14564392373896*rdxF[2]*f[25]*nu[193]+1.677050983124842*rdxF[2]*f[41]*nu[173]+0.75*rdxF[2]*f[25]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[42] += 1.14564392373896*rdxF[3]*f[22]*nu[274]+1.677050983124842*rdxF[3]*f[42]*nu[254]+0.75*rdxF[3]*f[22]*nu[244]+1.14564392373896*rdxF[2]*f[26]*nu[193]+1.677050983124842*rdxF[2]*f[42]*nu[173]+0.75*rdxF[2]*f[26]*nu[163]+6.708203932499369*rdxF[1]*f[70]*nu[112]+7.685213074469699*rdxF[1]*f[18]*nu[112]+5.031152949374527*rdxF[1]*f[42]*nu[92]+3.75*rdxF[1]*f[10]*nu[92]+3.354101966249685*rdxF[1]*f[18]*nu[82]+1.677050983124842*rdxF[1]*f[10]*nu[80]; 
  out[43] += 1.14564392373896*rdxF[3]*f[23]*nu[274]+1.677050983124842*rdxF[3]*f[43]*nu[254]+0.75*rdxF[3]*f[23]*nu[244]+6.708203932499369*rdxF[2]*f[71]*nu[193]+7.685213074469699*rdxF[2]*f[17]*nu[193]+5.031152949374527*rdxF[2]*f[43]*nu[173]+3.75*rdxF[2]*f[8]*nu[173]+3.354101966249685*rdxF[2]*f[17]*nu[163]+1.677050983124842*rdxF[2]*f[8]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[44] += 1.14564392373896*rdxF[3]*f[24]*nu[274]+1.677050983124842*rdxF[3]*f[44]*nu[254]+0.75*rdxF[3]*f[24]*nu[244]+6.708203932499369*rdxF[2]*f[72]*nu[193]+7.685213074469699*rdxF[2]*f[18]*nu[193]+5.031152949374527*rdxF[2]*f[44]*nu[173]+3.75*rdxF[2]*f[9]*nu[173]+3.354101966249685*rdxF[2]*f[18]*nu[163]+1.677050983124842*rdxF[2]*f[9]*nu[160]+1.14564392373896*rdxF[1]*f[27]*nu[112]+1.677050983124842*rdxF[1]*f[44]*nu[92]+0.75*rdxF[1]*f[27]*nu[82]; 
  out[45] += 6.708203932499369*rdxF[3]*f[73]*nu[274]+7.685213074469699*rdxF[3]*f[16]*nu[274]+5.031152949374527*rdxF[3]*f[45]*nu[254]+3.75*rdxF[3]*f[5]*nu[254]+3.354101966249685*rdxF[3]*f[16]*nu[244]+1.677050983124842*rdxF[3]*f[5]*nu[240]+1.14564392373896*rdxF[1]*f[28]*nu[112]+1.677050983124842*rdxF[1]*f[45]*nu[92]+0.75*rdxF[1]*f[28]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 6.708203932499369*rdxF[3]*f[74]*nu[274]+7.685213074469699*rdxF[3]*f[17]*nu[274]+5.031152949374527*rdxF[3]*f[46]*nu[254]+3.75*rdxF[3]*f[6]*nu[254]+3.354101966249685*rdxF[3]*f[17]*nu[244]+1.677050983124842*rdxF[3]*f[6]*nu[240]+1.14564392373896*rdxF[2]*f[28]*nu[193]+1.677050983124842*rdxF[2]*f[46]*nu[173]+0.75*rdxF[2]*f[28]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[47] += 6.708203932499369*rdxF[3]*f[75]*nu[274]+7.685213074469699*rdxF[3]*f[18]*nu[274]+5.031152949374527*rdxF[3]*f[47]*nu[254]+3.75*rdxF[3]*f[7]*nu[254]+3.354101966249685*rdxF[3]*f[18]*nu[244]+1.677050983124842*rdxF[3]*f[7]*nu[240]+1.14564392373896*rdxF[2]*f[29]*nu[193]+1.677050983124842*rdxF[2]*f[47]*nu[173]+0.75*rdxF[2]*f[29]*nu[163]+1.14564392373896*rdxF[1]*f[30]*nu[112]+1.677050983124842*rdxF[1]*f[47]*nu[92]+0.75*rdxF[1]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[1]*f[31]*nu[112]+1.677050983124842*rdxF[1]*f[48]*nu[92]+0.7499999999999999*rdxF[1]*f[31]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 18.44756081437327*rdxF[1]*f[20]*nu[112]+10.5*f[1]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[49]*nu[92]+12.8086884574495*rdxF[1]*f[5]*nu[92]+7.685213074469698*rdxF[1]*f[20]*nu[82]+6.873863542433758*f[1]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[5]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 1.14564392373896*rdxF[2]*f[31]*nu[193]+1.677050983124842*rdxF[2]*f[50]*nu[173]+0.7499999999999999*rdxF[2]*f[31]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[51] += 1.14564392373896*rdxF[2]*f[32]*nu[193]+1.677050983124842*rdxF[2]*f[51]*nu[173]+0.7499999999999999*rdxF[2]*f[32]*nu[163]+18.44756081437327*rdxF[1]*f[22]*nu[112]+10.5*rdxF[1]*f[3]*nu[112]+10.06230589874905*rdxF[1]*f[51]*nu[92]+12.8086884574495*rdxF[1]*f[7]*nu[92]+7.685213074469698*rdxF[1]*f[22]*nu[82]+6.873863542433758*rdxF[1]*f[3]*nu[82]+5.7282196186948*rdxF[1]*f[7]*nu[80]; 
  out[52] += 18.44756081437327*rdxF[2]*f[23]*nu[193]+10.5*f[1]*rdxF[2]*nu[193]+10.06230589874905*rdxF[2]*f[52]*nu[173]+12.8086884574495*rdxF[2]*f[6]*nu[173]+7.685213074469698*rdxF[2]*f[23]*nu[163]+6.873863542433758*f[1]*rdxF[2]*nu[163]+5.7282196186948*rdxF[2]*f[6]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[53] += 18.44756081437327*rdxF[2]*f[24]*nu[193]+10.5*f[2]*rdxF[2]*nu[193]+10.06230589874905*rdxF[2]*f[53]*nu[173]+12.8086884574495*rdxF[2]*f[7]*nu[173]+7.685213074469698*rdxF[2]*f[24]*nu[163]+6.873863542433758*f[2]*rdxF[2]*nu[163]+5.7282196186948*rdxF[2]*f[7]*nu[160]+1.14564392373896*rdxF[1]*f[33]*nu[112]+1.677050983124842*rdxF[1]*f[53]*nu[92]+0.7499999999999999*rdxF[1]*f[33]*nu[82]; 
  out[54] += 1.14564392373896*rdxF[3]*f[31]*nu[274]+1.677050983124842*rdxF[3]*f[54]*nu[254]+0.7499999999999999*rdxF[3]*f[31]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[55] += 1.14564392373896*rdxF[3]*f[32]*nu[274]+1.677050983124842*rdxF[3]*f[55]*nu[254]+0.7499999999999999*rdxF[3]*f[32]*nu[244]+18.44756081437327*rdxF[1]*f[26]*nu[112]+10.5*rdxF[1]*f[4]*nu[112]+10.06230589874905*rdxF[1]*f[55]*nu[92]+12.8086884574495*rdxF[1]*f[9]*nu[92]+7.685213074469698*rdxF[1]*f[26]*nu[82]+6.873863542433758*rdxF[1]*f[4]*nu[82]+5.7282196186948*rdxF[1]*f[9]*nu[80]; 
  out[56] += 1.14564392373896*rdxF[3]*f[33]*nu[274]+1.677050983124842*rdxF[3]*f[56]*nu[254]+0.7499999999999999*rdxF[3]*f[33]*nu[244]+18.44756081437327*rdxF[2]*f[27]*nu[193]+10.5*rdxF[2]*f[4]*nu[193]+10.06230589874905*rdxF[2]*f[56]*nu[173]+12.8086884574495*rdxF[2]*f[10]*nu[173]+7.685213074469698*rdxF[2]*f[27]*nu[163]+6.873863542433758*rdxF[2]*f[4]*nu[163]+5.7282196186948*rdxF[2]*f[10]*nu[160]; 
  out[57] += 18.44756081437327*rdxF[3]*f[28]*nu[274]+10.5*f[1]*rdxF[3]*nu[274]+10.06230589874905*rdxF[3]*f[57]*nu[254]+12.8086884574495*rdxF[3]*f[8]*nu[254]+7.685213074469698*rdxF[3]*f[28]*nu[244]+6.873863542433758*f[1]*rdxF[3]*nu[244]+5.7282196186948*rdxF[3]*f[8]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[58] += 18.44756081437327*rdxF[3]*f[29]*nu[274]+10.5*f[2]*rdxF[3]*nu[274]+10.06230589874905*rdxF[3]*f[58]*nu[254]+12.8086884574495*rdxF[3]*f[9]*nu[254]+7.685213074469698*rdxF[3]*f[29]*nu[244]+6.873863542433758*f[2]*rdxF[3]*nu[244]+5.7282196186948*rdxF[3]*f[9]*nu[240]+1.14564392373896*rdxF[1]*f[34]*nu[112]+1.677050983124842*rdxF[1]*f[58]*nu[92]+0.7499999999999999*rdxF[1]*f[34]*nu[82]; 
  out[59] += 18.44756081437327*rdxF[3]*f[30]*nu[274]+10.5*f[3]*rdxF[3]*nu[274]+10.06230589874905*rdxF[3]*f[59]*nu[254]+12.8086884574495*rdxF[3]*f[10]*nu[254]+7.685213074469698*rdxF[3]*f[30]*nu[244]+6.873863542433758*f[3]*rdxF[3]*nu[244]+5.7282196186948*rdxF[3]*f[10]*nu[240]+1.14564392373896*rdxF[2]*f[34]*nu[193]+1.677050983124842*rdxF[2]*f[59]*nu[173]+0.7499999999999999*rdxF[2]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[3]*f[36]*nu[274]+1.677050983124842*rdxF[3]*f[60]*nu[254]+0.75*rdxF[3]*f[36]*nu[244]+1.14564392373896*rdxF[2]*f[39]*nu[193]+1.677050983124842*rdxF[2]*f[60]*nu[173]+0.75*rdxF[2]*f[39]*nu[163]+1.14564392373896*rdxF[1]*f[41]*nu[112]+1.677050983124842*rdxF[1]*f[60]*nu[92]+0.75*rdxF[1]*f[41]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 1.14564392373896*rdxF[3]*f[37]*nu[274]+1.677050983124842*rdxF[3]*f[61]*nu[254]+0.75*rdxF[3]*f[37]*nu[244]+1.14564392373896*rdxF[2]*f[40]*nu[193]+1.677050983124842*rdxF[2]*f[61]*nu[173]+0.75*rdxF[2]*f[40]*nu[163]+6.708203932499369*rdxF[1]*f[77]*nu[112]+7.6852130744697*rdxF[1]*f[35]*nu[112]+5.031152949374527*rdxF[1]*f[61]*nu[92]+3.75*rdxF[1]*f[17]*nu[92]+3.354101966249684*rdxF[1]*f[35]*nu[82]+1.677050983124842*rdxF[1]*f[17]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 1.14564392373896*rdxF[3]*f[38]*nu[274]+1.677050983124842*rdxF[3]*f[62]*nu[254]+0.75*rdxF[3]*f[38]*nu[244]+6.708203932499369*rdxF[2]*f[78]*nu[193]+7.6852130744697*rdxF[2]*f[35]*nu[193]+5.031152949374527*rdxF[2]*f[62]*nu[173]+3.75*rdxF[2]*f[16]*nu[173]+3.354101966249684*rdxF[2]*f[35]*nu[163]+1.677050983124842*rdxF[2]*f[16]*nu[160]+1.14564392373896*rdxF[1]*f[43]*nu[112]+1.677050983124842*rdxF[1]*f[62]*nu[92]+0.75*rdxF[1]*f[43]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 6.708203932499369*rdxF[3]*f[79]*nu[274]+7.6852130744697*rdxF[3]*f[35]*nu[274]+5.031152949374527*rdxF[3]*f[63]*nu[254]+3.75*rdxF[3]*f[15]*nu[254]+3.354101966249684*rdxF[3]*f[35]*nu[244]+1.677050983124842*rdxF[3]*f[15]*nu[240]+1.14564392373896*rdxF[2]*f[45]*nu[193]+1.677050983124842*rdxF[2]*f[63]*nu[173]+0.75*rdxF[2]*f[45]*nu[163]+1.14564392373896*rdxF[1]*f[46]*nu[112]+1.677050983124842*rdxF[1]*f[63]*nu[92]+0.75*rdxF[1]*f[46]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 1.14564392373896*rdxF[2]*f[48]*nu[193]+1.677050983124842*rdxF[2]*f[64]*nu[173]+0.7499999999999999*rdxF[2]*f[48]*nu[163]+1.14564392373896*rdxF[1]*f[50]*nu[112]+1.677050983124842*rdxF[1]*f[64]*nu[92]+0.7499999999999999*rdxF[1]*f[50]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 1.14564392373896*rdxF[2]*f[49]*nu[193]+1.677050983124842*rdxF[2]*f[65]*nu[173]+0.7499999999999999*rdxF[2]*f[49]*nu[163]+18.44756081437327*rdxF[1]*f[37]*nu[112]+10.5*rdxF[1]*f[6]*nu[112]+10.06230589874905*rdxF[1]*f[65]*nu[92]+12.8086884574495*rdxF[1]*f[15]*nu[92]+7.685213074469699*rdxF[1]*f[37]*nu[82]+6.873863542433759*rdxF[1]*f[6]*nu[82]+5.7282196186948*rdxF[1]*f[15]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 18.44756081437327*rdxF[2]*f[38]*nu[193]+10.5*rdxF[2]*f[5]*nu[193]+10.06230589874905*rdxF[2]*f[66]*nu[173]+12.8086884574495*rdxF[2]*f[15]*nu[173]+7.685213074469699*rdxF[2]*f[38]*nu[163]+6.873863542433759*rdxF[2]*f[5]*nu[163]+5.7282196186948*rdxF[2]*f[15]*nu[160]+1.14564392373896*rdxF[1]*f[52]*nu[112]+1.677050983124842*rdxF[1]*f[66]*nu[92]+0.7499999999999999*rdxF[1]*f[52]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 1.14564392373896*rdxF[3]*f[48]*nu[274]+1.677050983124842*rdxF[3]*f[67]*nu[254]+0.7499999999999999*rdxF[3]*f[48]*nu[244]+1.14564392373896*rdxF[1]*f[54]*nu[112]+1.677050983124842*rdxF[1]*f[67]*nu[92]+0.7499999999999999*rdxF[1]*f[54]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 1.14564392373896*rdxF[3]*f[49]*nu[274]+1.677050983124842*rdxF[3]*f[68]*nu[254]+0.7499999999999999*rdxF[3]*f[49]*nu[244]+18.44756081437327*rdxF[1]*f[40]*nu[112]+10.5*rdxF[1]*f[8]*nu[112]+10.06230589874905*rdxF[1]*f[68]*nu[92]+12.8086884574495*rdxF[1]*f[16]*nu[92]+7.685213074469699*rdxF[1]*f[40]*nu[82]+6.873863542433759*rdxF[1]*f[8]*nu[82]+5.7282196186948*rdxF[1]*f[16]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 1.14564392373896*rdxF[3]*f[50]*nu[274]+1.677050983124842*rdxF[3]*f[69]*nu[254]+0.7499999999999999*rdxF[3]*f[50]*nu[244]+1.14564392373896*rdxF[2]*f[54]*nu[193]+1.677050983124842*rdxF[2]*f[69]*nu[173]+0.7499999999999999*rdxF[2]*f[54]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[70] += 1.14564392373896*rdxF[3]*f[51]*nu[274]+1.677050983124842*rdxF[3]*f[70]*nu[254]+0.7499999999999999*rdxF[3]*f[51]*nu[244]+1.14564392373896*rdxF[2]*f[55]*nu[193]+1.677050983124842*rdxF[2]*f[70]*nu[173]+0.7499999999999999*rdxF[2]*f[55]*nu[163]+18.44756081437327*rdxF[1]*f[42]*nu[112]+10.5*rdxF[1]*f[10]*nu[112]+10.06230589874905*rdxF[1]*f[70]*nu[92]+12.8086884574495*rdxF[1]*f[18]*nu[92]+7.685213074469699*rdxF[1]*f[42]*nu[82]+6.873863542433759*rdxF[1]*f[10]*nu[82]+5.7282196186948*rdxF[1]*f[18]*nu[80]; 
  out[71] += 1.14564392373896*rdxF[3]*f[52]*nu[274]+1.677050983124842*rdxF[3]*f[71]*nu[254]+0.7499999999999999*rdxF[3]*f[52]*nu[244]+18.44756081437327*rdxF[2]*f[43]*nu[193]+10.5*rdxF[2]*f[8]*nu[193]+10.06230589874905*rdxF[2]*f[71]*nu[173]+12.8086884574495*rdxF[2]*f[17]*nu[173]+7.685213074469699*rdxF[2]*f[43]*nu[163]+6.873863542433759*rdxF[2]*f[8]*nu[163]+5.7282196186948*rdxF[2]*f[17]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[72] += 1.14564392373896*rdxF[3]*f[53]*nu[274]+1.677050983124842*rdxF[3]*f[72]*nu[254]+0.7499999999999999*rdxF[3]*f[53]*nu[244]+18.44756081437327*rdxF[2]*f[44]*nu[193]+10.5*rdxF[2]*f[9]*nu[193]+10.06230589874905*rdxF[2]*f[72]*nu[173]+12.8086884574495*rdxF[2]*f[18]*nu[173]+7.685213074469699*rdxF[2]*f[44]*nu[163]+6.873863542433759*rdxF[2]*f[9]*nu[163]+5.7282196186948*rdxF[2]*f[18]*nu[160]+1.14564392373896*rdxF[1]*f[56]*nu[112]+1.677050983124842*rdxF[1]*f[72]*nu[92]+0.7499999999999999*rdxF[1]*f[56]*nu[82]; 
  out[73] += 18.44756081437327*rdxF[3]*f[45]*nu[274]+10.5*rdxF[3]*f[5]*nu[274]+10.06230589874905*rdxF[3]*f[73]*nu[254]+12.8086884574495*rdxF[3]*f[16]*nu[254]+7.685213074469699*rdxF[3]*f[45]*nu[244]+6.873863542433759*rdxF[3]*f[5]*nu[244]+5.7282196186948*rdxF[3]*f[16]*nu[240]+1.14564392373896*rdxF[1]*f[57]*nu[112]+1.677050983124842*rdxF[1]*f[73]*nu[92]+0.7499999999999999*rdxF[1]*f[57]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 18.44756081437327*rdxF[3]*f[46]*nu[274]+10.5*rdxF[3]*f[6]*nu[274]+10.06230589874905*rdxF[3]*f[74]*nu[254]+12.8086884574495*rdxF[3]*f[17]*nu[254]+7.685213074469699*rdxF[3]*f[46]*nu[244]+6.873863542433759*rdxF[3]*f[6]*nu[244]+5.7282196186948*rdxF[3]*f[17]*nu[240]+1.14564392373896*rdxF[2]*f[57]*nu[193]+1.677050983124842*rdxF[2]*f[74]*nu[173]+0.7499999999999999*rdxF[2]*f[57]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[75] += 18.44756081437327*rdxF[3]*f[47]*nu[274]+10.5*rdxF[3]*f[7]*nu[274]+10.06230589874905*rdxF[3]*f[75]*nu[254]+12.8086884574495*rdxF[3]*f[18]*nu[254]+7.685213074469699*rdxF[3]*f[47]*nu[244]+6.873863542433759*rdxF[3]*f[7]*nu[244]+5.7282196186948*rdxF[3]*f[18]*nu[240]+1.14564392373896*rdxF[2]*f[58]*nu[193]+1.677050983124842*rdxF[2]*f[75]*nu[173]+0.7499999999999999*rdxF[2]*f[58]*nu[163]+1.14564392373896*rdxF[1]*f[59]*nu[112]+1.677050983124842*rdxF[1]*f[75]*nu[92]+0.7499999999999999*rdxF[1]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[3]*f[64]*nu[274]+1.677050983124842*rdxF[3]*f[76]*nu[254]+0.7499999999999999*rdxF[3]*f[64]*nu[244]+1.14564392373896*rdxF[2]*f[67]*nu[193]+1.677050983124842*rdxF[2]*f[76]*nu[173]+0.7499999999999999*rdxF[2]*f[67]*nu[163]+1.14564392373896*rdxF[1]*f[69]*nu[112]+1.677050983124842*rdxF[1]*f[76]*nu[92]+0.7499999999999999*rdxF[1]*f[69]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 1.14564392373896*rdxF[3]*f[65]*nu[274]+1.677050983124842*rdxF[3]*f[77]*nu[254]+0.7499999999999999*rdxF[3]*f[65]*nu[244]+1.14564392373896*rdxF[2]*f[68]*nu[193]+1.677050983124842*rdxF[2]*f[77]*nu[173]+0.7499999999999999*rdxF[2]*f[68]*nu[163]+18.44756081437327*rdxF[1]*f[61]*nu[112]+10.5*rdxF[1]*f[17]*nu[112]+10.06230589874905*rdxF[1]*f[77]*nu[92]+12.8086884574495*rdxF[1]*f[35]*nu[92]+7.685213074469698*rdxF[1]*f[61]*nu[82]+6.873863542433758*rdxF[1]*f[17]*nu[82]+5.7282196186948*rdxF[1]*f[35]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 1.14564392373896*rdxF[3]*f[66]*nu[274]+1.677050983124842*rdxF[3]*f[78]*nu[254]+0.7499999999999999*rdxF[3]*f[66]*nu[244]+18.44756081437327*rdxF[2]*f[62]*nu[193]+10.5*rdxF[2]*f[16]*nu[193]+10.06230589874905*rdxF[2]*f[78]*nu[173]+12.8086884574495*rdxF[2]*f[35]*nu[173]+7.685213074469698*rdxF[2]*f[62]*nu[163]+6.873863542433758*rdxF[2]*f[16]*nu[163]+5.7282196186948*rdxF[2]*f[35]*nu[160]+1.14564392373896*rdxF[1]*f[71]*nu[112]+1.677050983124842*rdxF[1]*f[78]*nu[92]+0.7499999999999999*rdxF[1]*f[71]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 18.44756081437327*rdxF[3]*f[63]*nu[274]+10.5*rdxF[3]*f[15]*nu[274]+10.06230589874905*rdxF[3]*f[79]*nu[254]+12.8086884574495*rdxF[3]*f[35]*nu[254]+7.685213074469698*rdxF[3]*f[63]*nu[244]+6.873863542433758*rdxF[3]*f[15]*nu[244]+5.7282196186948*rdxF[3]*f[35]*nu[240]+1.14564392373896*rdxF[2]*f[73]*nu[193]+1.677050983124842*rdxF[2]*f[79]*nu[173]+0.7499999999999999*rdxF[2]*f[73]*nu[163]+1.14564392373896*rdxF[1]*f[74]*nu[112]+1.677050983124842*rdxF[1]*f[79]*nu[92]+0.7499999999999999*rdxF[1]*f[74]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+rdxF[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+rdxF[3]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs123(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[2]*dx[2]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 2.5617376914899*rdxF[1]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[2]*nu[92]+0.75*f[0]*rdxF[1]*nu[82]; 
  out[3] += 2.5617376914899*rdxF[2]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[2]*nu[193]+1.677050983124842*rdxF[2]*f[3]*nu[173]+0.75*f[0]*rdxF[2]*nu[163]; 
  out[5] += 2.5617376914899*rdxF[1]*f[20]*nu[112]+1.14564392373896*f[1]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[5]*nu[92]+0.75*f[1]*rdxF[1]*nu[82]+2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[2]*f[23]*nu[193]+1.14564392373896*f[1]*rdxF[2]*nu[193]+1.677050983124842*rdxF[2]*f[6]*nu[173]+0.75*f[1]*rdxF[2]*nu[163]+2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 2.5617376914899*rdxF[2]*f[24]*nu[193]+1.14564392373896*f[2]*rdxF[2]*nu[193]+1.677050983124842*rdxF[2]*f[7]*nu[173]+0.75*f[2]*rdxF[2]*nu[163]+2.5617376914899*rdxF[1]*f[22]*nu[112]+1.14564392373896*rdxF[1]*f[3]*nu[112]+1.677050983124842*rdxF[1]*f[7]*nu[92]+0.75*rdxF[1]*f[3]*nu[82]; 
  out[8] += 2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 2.5617376914899*rdxF[1]*f[26]*nu[112]+1.14564392373896*rdxF[1]*f[4]*nu[112]+1.677050983124842*rdxF[1]*f[9]*nu[92]+0.75*rdxF[1]*f[4]*nu[82]; 
  out[10] += 2.5617376914899*rdxF[2]*f[27]*nu[193]+1.14564392373896*rdxF[2]*f[4]*nu[193]+1.677050983124842*rdxF[2]*f[10]*nu[173]+0.75*rdxF[2]*f[4]*nu[163]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 6.708203932499369*rdxF[1]*f[32]*nu[112]+7.685213074469699*rdxF[1]*f[2]*nu[112]+5.031152949374527*rdxF[1]*f[12]*nu[92]+3.75*f[0]*rdxF[1]*nu[92]+3.354101966249685*rdxF[1]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[1]*nu[80]; 
  out[13] += 6.708203932499369*rdxF[2]*f[33]*nu[193]+7.685213074469699*rdxF[2]*f[3]*nu[193]+5.031152949374527*rdxF[2]*f[13]*nu[173]+3.75*f[0]*rdxF[2]*nu[173]+3.354101966249685*rdxF[2]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[2]*nu[160]; 
  out[15] += 2.5617376914899*rdxF[2]*f[38]*nu[193]+1.14564392373896*rdxF[2]*f[5]*nu[193]+1.677050983124842*rdxF[2]*f[15]*nu[173]+0.75*rdxF[2]*f[5]*nu[163]+2.5617376914899*rdxF[1]*f[37]*nu[112]+1.14564392373896*rdxF[1]*f[6]*nu[112]+1.677050983124842*rdxF[1]*f[15]*nu[92]+0.75*rdxF[1]*f[6]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[1]*f[40]*nu[112]+1.14564392373896*rdxF[1]*f[8]*nu[112]+1.677050983124842*rdxF[1]*f[16]*nu[92]+0.75*rdxF[1]*f[8]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[2]*f[43]*nu[193]+1.14564392373896*rdxF[2]*f[8]*nu[193]+1.677050983124842*rdxF[2]*f[17]*nu[173]+0.75*rdxF[2]*f[8]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 2.5617376914899*rdxF[2]*f[44]*nu[193]+1.14564392373896*rdxF[2]*f[9]*nu[193]+1.677050983124842*rdxF[2]*f[18]*nu[173]+0.75*rdxF[2]*f[9]*nu[163]+2.5617376914899*rdxF[1]*f[42]*nu[112]+1.14564392373896*rdxF[1]*f[10]*nu[112]+1.677050983124842*rdxF[1]*f[18]*nu[92]+0.75*rdxF[1]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[1]*f[11]*nu[112]+1.677050983124842*rdxF[1]*f[19]*nu[92]+0.75*rdxF[1]*f[11]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 6.708203932499369*rdxF[1]*f[49]*nu[112]+7.6852130744697*rdxF[1]*f[5]*nu[112]+5.031152949374527*rdxF[1]*f[20]*nu[92]+3.75*f[1]*rdxF[1]*nu[92]+3.354101966249684*rdxF[1]*f[5]*nu[82]+1.677050983124842*f[1]*rdxF[1]*nu[80]+1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.14564392373896*rdxF[2]*f[11]*nu[193]+1.677050983124842*rdxF[2]*f[21]*nu[173]+0.75*rdxF[2]*f[11]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.14564392373896*rdxF[2]*f[12]*nu[193]+1.677050983124842*rdxF[2]*f[22]*nu[173]+0.75*rdxF[2]*f[12]*nu[163]+6.708203932499369*rdxF[1]*f[51]*nu[112]+7.6852130744697*rdxF[1]*f[7]*nu[112]+5.031152949374527*rdxF[1]*f[22]*nu[92]+3.75*rdxF[1]*f[3]*nu[92]+3.354101966249684*rdxF[1]*f[7]*nu[82]+1.677050983124842*rdxF[1]*f[3]*nu[80]; 
  out[23] += 6.708203932499369*rdxF[2]*f[52]*nu[193]+7.6852130744697*rdxF[2]*f[6]*nu[193]+5.031152949374527*rdxF[2]*f[23]*nu[173]+3.75*f[1]*rdxF[2]*nu[173]+3.354101966249684*rdxF[2]*f[6]*nu[163]+1.677050983124842*f[1]*rdxF[2]*nu[160]+1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 6.708203932499369*rdxF[2]*f[53]*nu[193]+7.6852130744697*rdxF[2]*f[7]*nu[193]+5.031152949374527*rdxF[2]*f[24]*nu[173]+3.75*f[2]*rdxF[2]*nu[173]+3.354101966249684*rdxF[2]*f[7]*nu[163]+1.677050983124842*f[2]*rdxF[2]*nu[160]+1.14564392373896*rdxF[1]*f[13]*nu[112]+1.677050983124842*rdxF[1]*f[24]*nu[92]+0.75*rdxF[1]*f[13]*nu[82]; 
  out[25] += 6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 6.708203932499369*rdxF[1]*f[55]*nu[112]+7.6852130744697*rdxF[1]*f[9]*nu[112]+5.031152949374527*rdxF[1]*f[26]*nu[92]+3.75*rdxF[1]*f[4]*nu[92]+3.354101966249684*rdxF[1]*f[9]*nu[82]+1.677050983124842*rdxF[1]*f[4]*nu[80]; 
  out[27] += 6.708203932499369*rdxF[2]*f[56]*nu[193]+7.6852130744697*rdxF[2]*f[10]*nu[193]+5.031152949374527*rdxF[2]*f[27]*nu[173]+3.75*rdxF[2]*f[4]*nu[173]+3.354101966249684*rdxF[2]*f[10]*nu[163]+1.677050983124842*rdxF[2]*f[4]*nu[160]; 
  out[28] += 1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 1.14564392373896*rdxF[1]*f[14]*nu[112]+1.677050983124842*rdxF[1]*f[29]*nu[92]+0.75*rdxF[1]*f[14]*nu[82]; 
  out[30] += 1.14564392373896*rdxF[2]*f[14]*nu[193]+1.677050983124842*rdxF[2]*f[30]*nu[173]+0.75*rdxF[2]*f[14]*nu[163]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[32] += 18.44756081437327*rdxF[1]*f[12]*nu[112]+10.5*f[0]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[32]*nu[92]+12.8086884574495*rdxF[1]*f[2]*nu[92]+7.685213074469699*rdxF[1]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[2]*nu[80]; 
  out[33] += 18.44756081437327*rdxF[2]*f[13]*nu[193]+10.5*f[0]*rdxF[2]*nu[193]+10.06230589874905*rdxF[2]*f[33]*nu[173]+12.8086884574495*rdxF[2]*f[3]*nu[173]+7.685213074469699*rdxF[2]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[2]*nu[163]+5.7282196186948*rdxF[2]*f[3]*nu[160]; 
  out[35] += 2.5617376914899*rdxF[2]*f[62]*nu[193]+1.14564392373896*rdxF[2]*f[16]*nu[193]+1.677050983124842*rdxF[2]*f[35]*nu[173]+0.75*rdxF[2]*f[16]*nu[163]+2.5617376914899*rdxF[1]*f[61]*nu[112]+1.14564392373896*rdxF[1]*f[17]*nu[112]+1.677050983124842*rdxF[1]*f[35]*nu[92]+0.75*rdxF[1]*f[17]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 1.14564392373896*rdxF[2]*f[19]*nu[193]+1.677050983124842*rdxF[2]*f[36]*nu[173]+0.75*rdxF[2]*f[19]*nu[163]+1.14564392373896*rdxF[1]*f[21]*nu[112]+1.677050983124842*rdxF[1]*f[36]*nu[92]+0.75*rdxF[1]*f[21]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 1.14564392373896*rdxF[2]*f[20]*nu[193]+1.677050983124842*rdxF[2]*f[37]*nu[173]+0.75*rdxF[2]*f[20]*nu[163]+6.708203932499369*rdxF[1]*f[65]*nu[112]+7.685213074469699*rdxF[1]*f[15]*nu[112]+5.031152949374527*rdxF[1]*f[37]*nu[92]+3.75*rdxF[1]*f[6]*nu[92]+3.354101966249685*rdxF[1]*f[15]*nu[82]+1.677050983124842*rdxF[1]*f[6]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 6.708203932499369*rdxF[2]*f[66]*nu[193]+7.685213074469699*rdxF[2]*f[15]*nu[193]+5.031152949374527*rdxF[2]*f[38]*nu[173]+3.75*rdxF[2]*f[5]*nu[173]+3.354101966249685*rdxF[2]*f[15]*nu[163]+1.677050983124842*rdxF[2]*f[5]*nu[160]+1.14564392373896*rdxF[1]*f[23]*nu[112]+1.677050983124842*rdxF[1]*f[38]*nu[92]+0.75*rdxF[1]*f[23]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 1.14564392373896*rdxF[1]*f[25]*nu[112]+1.677050983124842*rdxF[1]*f[39]*nu[92]+0.75*rdxF[1]*f[25]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 6.708203932499369*rdxF[1]*f[68]*nu[112]+7.685213074469699*rdxF[1]*f[16]*nu[112]+5.031152949374527*rdxF[1]*f[40]*nu[92]+3.75*rdxF[1]*f[8]*nu[92]+3.354101966249685*rdxF[1]*f[16]*nu[82]+1.677050983124842*rdxF[1]*f[8]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 1.14564392373896*rdxF[2]*f[25]*nu[193]+1.677050983124842*rdxF[2]*f[41]*nu[173]+0.75*rdxF[2]*f[25]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[42] += 1.14564392373896*rdxF[2]*f[26]*nu[193]+1.677050983124842*rdxF[2]*f[42]*nu[173]+0.75*rdxF[2]*f[26]*nu[163]+6.708203932499369*rdxF[1]*f[70]*nu[112]+7.685213074469699*rdxF[1]*f[18]*nu[112]+5.031152949374527*rdxF[1]*f[42]*nu[92]+3.75*rdxF[1]*f[10]*nu[92]+3.354101966249685*rdxF[1]*f[18]*nu[82]+1.677050983124842*rdxF[1]*f[10]*nu[80]; 
  out[43] += 6.708203932499369*rdxF[2]*f[71]*nu[193]+7.685213074469699*rdxF[2]*f[17]*nu[193]+5.031152949374527*rdxF[2]*f[43]*nu[173]+3.75*rdxF[2]*f[8]*nu[173]+3.354101966249685*rdxF[2]*f[17]*nu[163]+1.677050983124842*rdxF[2]*f[8]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[44] += 6.708203932499369*rdxF[2]*f[72]*nu[193]+7.685213074469699*rdxF[2]*f[18]*nu[193]+5.031152949374527*rdxF[2]*f[44]*nu[173]+3.75*rdxF[2]*f[9]*nu[173]+3.354101966249685*rdxF[2]*f[18]*nu[163]+1.677050983124842*rdxF[2]*f[9]*nu[160]+1.14564392373896*rdxF[1]*f[27]*nu[112]+1.677050983124842*rdxF[1]*f[44]*nu[92]+0.75*rdxF[1]*f[27]*nu[82]; 
  out[45] += 1.14564392373896*rdxF[1]*f[28]*nu[112]+1.677050983124842*rdxF[1]*f[45]*nu[92]+0.75*rdxF[1]*f[28]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 1.14564392373896*rdxF[2]*f[28]*nu[193]+1.677050983124842*rdxF[2]*f[46]*nu[173]+0.75*rdxF[2]*f[28]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[47] += 1.14564392373896*rdxF[2]*f[29]*nu[193]+1.677050983124842*rdxF[2]*f[47]*nu[173]+0.75*rdxF[2]*f[29]*nu[163]+1.14564392373896*rdxF[1]*f[30]*nu[112]+1.677050983124842*rdxF[1]*f[47]*nu[92]+0.75*rdxF[1]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[1]*f[31]*nu[112]+1.677050983124842*rdxF[1]*f[48]*nu[92]+0.7499999999999999*rdxF[1]*f[31]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 18.44756081437327*rdxF[1]*f[20]*nu[112]+10.5*f[1]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[49]*nu[92]+12.8086884574495*rdxF[1]*f[5]*nu[92]+7.685213074469698*rdxF[1]*f[20]*nu[82]+6.873863542433758*f[1]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[5]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 1.14564392373896*rdxF[2]*f[31]*nu[193]+1.677050983124842*rdxF[2]*f[50]*nu[173]+0.7499999999999999*rdxF[2]*f[31]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[51] += 1.14564392373896*rdxF[2]*f[32]*nu[193]+1.677050983124842*rdxF[2]*f[51]*nu[173]+0.7499999999999999*rdxF[2]*f[32]*nu[163]+18.44756081437327*rdxF[1]*f[22]*nu[112]+10.5*rdxF[1]*f[3]*nu[112]+10.06230589874905*rdxF[1]*f[51]*nu[92]+12.8086884574495*rdxF[1]*f[7]*nu[92]+7.685213074469698*rdxF[1]*f[22]*nu[82]+6.873863542433758*rdxF[1]*f[3]*nu[82]+5.7282196186948*rdxF[1]*f[7]*nu[80]; 
  out[52] += 18.44756081437327*rdxF[2]*f[23]*nu[193]+10.5*f[1]*rdxF[2]*nu[193]+10.06230589874905*rdxF[2]*f[52]*nu[173]+12.8086884574495*rdxF[2]*f[6]*nu[173]+7.685213074469698*rdxF[2]*f[23]*nu[163]+6.873863542433758*f[1]*rdxF[2]*nu[163]+5.7282196186948*rdxF[2]*f[6]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[53] += 18.44756081437327*rdxF[2]*f[24]*nu[193]+10.5*f[2]*rdxF[2]*nu[193]+10.06230589874905*rdxF[2]*f[53]*nu[173]+12.8086884574495*rdxF[2]*f[7]*nu[173]+7.685213074469698*rdxF[2]*f[24]*nu[163]+6.873863542433758*f[2]*rdxF[2]*nu[163]+5.7282196186948*rdxF[2]*f[7]*nu[160]+1.14564392373896*rdxF[1]*f[33]*nu[112]+1.677050983124842*rdxF[1]*f[53]*nu[92]+0.7499999999999999*rdxF[1]*f[33]*nu[82]; 
  out[54] += 10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[55] += 18.44756081437327*rdxF[1]*f[26]*nu[112]+10.5*rdxF[1]*f[4]*nu[112]+10.06230589874905*rdxF[1]*f[55]*nu[92]+12.8086884574495*rdxF[1]*f[9]*nu[92]+7.685213074469698*rdxF[1]*f[26]*nu[82]+6.873863542433758*rdxF[1]*f[4]*nu[82]+5.7282196186948*rdxF[1]*f[9]*nu[80]; 
  out[56] += 18.44756081437327*rdxF[2]*f[27]*nu[193]+10.5*rdxF[2]*f[4]*nu[193]+10.06230589874905*rdxF[2]*f[56]*nu[173]+12.8086884574495*rdxF[2]*f[10]*nu[173]+7.685213074469698*rdxF[2]*f[27]*nu[163]+6.873863542433758*rdxF[2]*f[4]*nu[163]+5.7282196186948*rdxF[2]*f[10]*nu[160]; 
  out[57] += 1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[58] += 1.14564392373896*rdxF[1]*f[34]*nu[112]+1.677050983124842*rdxF[1]*f[58]*nu[92]+0.7499999999999999*rdxF[1]*f[34]*nu[82]; 
  out[59] += 1.14564392373896*rdxF[2]*f[34]*nu[193]+1.677050983124842*rdxF[2]*f[59]*nu[173]+0.7499999999999999*rdxF[2]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[2]*f[39]*nu[193]+1.677050983124842*rdxF[2]*f[60]*nu[173]+0.75*rdxF[2]*f[39]*nu[163]+1.14564392373896*rdxF[1]*f[41]*nu[112]+1.677050983124842*rdxF[1]*f[60]*nu[92]+0.75*rdxF[1]*f[41]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 1.14564392373896*rdxF[2]*f[40]*nu[193]+1.677050983124842*rdxF[2]*f[61]*nu[173]+0.75*rdxF[2]*f[40]*nu[163]+6.708203932499369*rdxF[1]*f[77]*nu[112]+7.6852130744697*rdxF[1]*f[35]*nu[112]+5.031152949374527*rdxF[1]*f[61]*nu[92]+3.75*rdxF[1]*f[17]*nu[92]+3.354101966249684*rdxF[1]*f[35]*nu[82]+1.677050983124842*rdxF[1]*f[17]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 6.708203932499369*rdxF[2]*f[78]*nu[193]+7.6852130744697*rdxF[2]*f[35]*nu[193]+5.031152949374527*rdxF[2]*f[62]*nu[173]+3.75*rdxF[2]*f[16]*nu[173]+3.354101966249684*rdxF[2]*f[35]*nu[163]+1.677050983124842*rdxF[2]*f[16]*nu[160]+1.14564392373896*rdxF[1]*f[43]*nu[112]+1.677050983124842*rdxF[1]*f[62]*nu[92]+0.75*rdxF[1]*f[43]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 1.14564392373896*rdxF[2]*f[45]*nu[193]+1.677050983124842*rdxF[2]*f[63]*nu[173]+0.75*rdxF[2]*f[45]*nu[163]+1.14564392373896*rdxF[1]*f[46]*nu[112]+1.677050983124842*rdxF[1]*f[63]*nu[92]+0.75*rdxF[1]*f[46]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 1.14564392373896*rdxF[2]*f[48]*nu[193]+1.677050983124842*rdxF[2]*f[64]*nu[173]+0.7499999999999999*rdxF[2]*f[48]*nu[163]+1.14564392373896*rdxF[1]*f[50]*nu[112]+1.677050983124842*rdxF[1]*f[64]*nu[92]+0.7499999999999999*rdxF[1]*f[50]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 1.14564392373896*rdxF[2]*f[49]*nu[193]+1.677050983124842*rdxF[2]*f[65]*nu[173]+0.7499999999999999*rdxF[2]*f[49]*nu[163]+18.44756081437327*rdxF[1]*f[37]*nu[112]+10.5*rdxF[1]*f[6]*nu[112]+10.06230589874905*rdxF[1]*f[65]*nu[92]+12.8086884574495*rdxF[1]*f[15]*nu[92]+7.685213074469699*rdxF[1]*f[37]*nu[82]+6.873863542433759*rdxF[1]*f[6]*nu[82]+5.7282196186948*rdxF[1]*f[15]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 18.44756081437327*rdxF[2]*f[38]*nu[193]+10.5*rdxF[2]*f[5]*nu[193]+10.06230589874905*rdxF[2]*f[66]*nu[173]+12.8086884574495*rdxF[2]*f[15]*nu[173]+7.685213074469699*rdxF[2]*f[38]*nu[163]+6.873863542433759*rdxF[2]*f[5]*nu[163]+5.7282196186948*rdxF[2]*f[15]*nu[160]+1.14564392373896*rdxF[1]*f[52]*nu[112]+1.677050983124842*rdxF[1]*f[66]*nu[92]+0.7499999999999999*rdxF[1]*f[52]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 1.14564392373896*rdxF[1]*f[54]*nu[112]+1.677050983124842*rdxF[1]*f[67]*nu[92]+0.7499999999999999*rdxF[1]*f[54]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 18.44756081437327*rdxF[1]*f[40]*nu[112]+10.5*rdxF[1]*f[8]*nu[112]+10.06230589874905*rdxF[1]*f[68]*nu[92]+12.8086884574495*rdxF[1]*f[16]*nu[92]+7.685213074469699*rdxF[1]*f[40]*nu[82]+6.873863542433759*rdxF[1]*f[8]*nu[82]+5.7282196186948*rdxF[1]*f[16]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 1.14564392373896*rdxF[2]*f[54]*nu[193]+1.677050983124842*rdxF[2]*f[69]*nu[173]+0.7499999999999999*rdxF[2]*f[54]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[70] += 1.14564392373896*rdxF[2]*f[55]*nu[193]+1.677050983124842*rdxF[2]*f[70]*nu[173]+0.7499999999999999*rdxF[2]*f[55]*nu[163]+18.44756081437327*rdxF[1]*f[42]*nu[112]+10.5*rdxF[1]*f[10]*nu[112]+10.06230589874905*rdxF[1]*f[70]*nu[92]+12.8086884574495*rdxF[1]*f[18]*nu[92]+7.685213074469699*rdxF[1]*f[42]*nu[82]+6.873863542433759*rdxF[1]*f[10]*nu[82]+5.7282196186948*rdxF[1]*f[18]*nu[80]; 
  out[71] += 18.44756081437327*rdxF[2]*f[43]*nu[193]+10.5*rdxF[2]*f[8]*nu[193]+10.06230589874905*rdxF[2]*f[71]*nu[173]+12.8086884574495*rdxF[2]*f[17]*nu[173]+7.685213074469699*rdxF[2]*f[43]*nu[163]+6.873863542433759*rdxF[2]*f[8]*nu[163]+5.7282196186948*rdxF[2]*f[17]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[72] += 18.44756081437327*rdxF[2]*f[44]*nu[193]+10.5*rdxF[2]*f[9]*nu[193]+10.06230589874905*rdxF[2]*f[72]*nu[173]+12.8086884574495*rdxF[2]*f[18]*nu[173]+7.685213074469699*rdxF[2]*f[44]*nu[163]+6.873863542433759*rdxF[2]*f[9]*nu[163]+5.7282196186948*rdxF[2]*f[18]*nu[160]+1.14564392373896*rdxF[1]*f[56]*nu[112]+1.677050983124842*rdxF[1]*f[72]*nu[92]+0.7499999999999999*rdxF[1]*f[56]*nu[82]; 
  out[73] += 1.14564392373896*rdxF[1]*f[57]*nu[112]+1.677050983124842*rdxF[1]*f[73]*nu[92]+0.7499999999999999*rdxF[1]*f[57]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 1.14564392373896*rdxF[2]*f[57]*nu[193]+1.677050983124842*rdxF[2]*f[74]*nu[173]+0.7499999999999999*rdxF[2]*f[57]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[75] += 1.14564392373896*rdxF[2]*f[58]*nu[193]+1.677050983124842*rdxF[2]*f[75]*nu[173]+0.7499999999999999*rdxF[2]*f[58]*nu[163]+1.14564392373896*rdxF[1]*f[59]*nu[112]+1.677050983124842*rdxF[1]*f[75]*nu[92]+0.7499999999999999*rdxF[1]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[2]*f[67]*nu[193]+1.677050983124842*rdxF[2]*f[76]*nu[173]+0.7499999999999999*rdxF[2]*f[67]*nu[163]+1.14564392373896*rdxF[1]*f[69]*nu[112]+1.677050983124842*rdxF[1]*f[76]*nu[92]+0.7499999999999999*rdxF[1]*f[69]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 1.14564392373896*rdxF[2]*f[68]*nu[193]+1.677050983124842*rdxF[2]*f[77]*nu[173]+0.7499999999999999*rdxF[2]*f[68]*nu[163]+18.44756081437327*rdxF[1]*f[61]*nu[112]+10.5*rdxF[1]*f[17]*nu[112]+10.06230589874905*rdxF[1]*f[77]*nu[92]+12.8086884574495*rdxF[1]*f[35]*nu[92]+7.685213074469698*rdxF[1]*f[61]*nu[82]+6.873863542433758*rdxF[1]*f[17]*nu[82]+5.7282196186948*rdxF[1]*f[35]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 18.44756081437327*rdxF[2]*f[62]*nu[193]+10.5*rdxF[2]*f[16]*nu[193]+10.06230589874905*rdxF[2]*f[78]*nu[173]+12.8086884574495*rdxF[2]*f[35]*nu[173]+7.685213074469698*rdxF[2]*f[62]*nu[163]+6.873863542433758*rdxF[2]*f[16]*nu[163]+5.7282196186948*rdxF[2]*f[35]*nu[160]+1.14564392373896*rdxF[1]*f[71]*nu[112]+1.677050983124842*rdxF[1]*f[78]*nu[92]+0.7499999999999999*rdxF[1]*f[71]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 1.14564392373896*rdxF[2]*f[73]*nu[193]+1.677050983124842*rdxF[2]*f[79]*nu[173]+0.7499999999999999*rdxF[2]*f[73]*nu[163]+1.14564392373896*rdxF[1]*f[74]*nu[112]+1.677050983124842*rdxF[1]*f[79]*nu[92]+0.7499999999999999*rdxF[1]*f[74]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+rdxF[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs24(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[1]*dx[1]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[2] += 2.5617376914899*rdxF[0]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[0]*nu[112]+1.677050983124842*rdxF[0]*f[2]*nu[92]+0.75*f[0]*rdxF[0]*nu[82]; 
  out[4] += 2.5617376914899*rdxF[1]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[1]*nu[274]+1.677050983124842*rdxF[1]*f[4]*nu[254]+0.75*f[0]*rdxF[1]*nu[244]; 
  out[5] += 2.5617376914899*rdxF[0]*f[20]*nu[112]+1.14564392373896*rdxF[0]*f[1]*nu[112]+1.677050983124842*rdxF[0]*f[5]*nu[92]+0.75*rdxF[0]*f[1]*nu[82]; 
  out[7] += 2.5617376914899*rdxF[0]*f[22]*nu[112]+1.14564392373896*rdxF[0]*f[3]*nu[112]+1.677050983124842*rdxF[0]*f[7]*nu[92]+0.75*rdxF[0]*f[3]*nu[82]; 
  out[8] += 2.5617376914899*rdxF[1]*f[28]*nu[274]+1.14564392373896*f[1]*rdxF[1]*nu[274]+1.677050983124842*rdxF[1]*f[8]*nu[254]+0.75*f[1]*rdxF[1]*nu[244]; 
  out[9] += 2.5617376914899*rdxF[1]*f[29]*nu[274]+1.14564392373896*rdxF[1]*f[2]*nu[274]+1.677050983124842*rdxF[1]*f[9]*nu[254]+0.75*rdxF[1]*f[2]*nu[244]+2.5617376914899*rdxF[0]*f[26]*nu[112]+1.14564392373896*rdxF[0]*f[4]*nu[112]+1.677050983124842*rdxF[0]*f[9]*nu[92]+0.75*rdxF[0]*f[4]*nu[82]; 
  out[10] += 2.5617376914899*rdxF[1]*f[30]*nu[274]+1.14564392373896*rdxF[1]*f[3]*nu[274]+1.677050983124842*rdxF[1]*f[10]*nu[254]+0.75*rdxF[1]*f[3]*nu[244]; 
  out[12] += 6.708203932499369*rdxF[0]*f[32]*nu[112]+7.685213074469699*rdxF[0]*f[2]*nu[112]+5.031152949374527*rdxF[0]*f[12]*nu[92]+3.75*f[0]*rdxF[0]*nu[92]+3.354101966249685*rdxF[0]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[0]*nu[80]; 
  out[14] += 6.708203932499369*rdxF[1]*f[34]*nu[274]+7.685213074469699*rdxF[1]*f[4]*nu[274]+5.031152949374527*rdxF[1]*f[14]*nu[254]+3.75*f[0]*rdxF[1]*nu[254]+3.354101966249685*rdxF[1]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[1]*nu[240]; 
  out[15] += 2.5617376914899*rdxF[0]*f[37]*nu[112]+1.14564392373896*rdxF[0]*f[6]*nu[112]+1.677050983124842*rdxF[0]*f[15]*nu[92]+0.75*rdxF[0]*f[6]*nu[82]; 
  out[16] += 2.5617376914899*rdxF[1]*f[45]*nu[274]+1.14564392373896*rdxF[1]*f[5]*nu[274]+1.677050983124842*rdxF[1]*f[16]*nu[254]+0.75*rdxF[1]*f[5]*nu[244]+2.5617376914899*rdxF[0]*f[40]*nu[112]+1.14564392373896*rdxF[0]*f[8]*nu[112]+1.677050983124842*rdxF[0]*f[16]*nu[92]+0.75*rdxF[0]*f[8]*nu[82]; 
  out[17] += 2.5617376914899*rdxF[1]*f[46]*nu[274]+1.14564392373896*rdxF[1]*f[6]*nu[274]+1.677050983124842*rdxF[1]*f[17]*nu[254]+0.75*rdxF[1]*f[6]*nu[244]; 
  out[18] += 2.5617376914899*rdxF[1]*f[47]*nu[274]+1.14564392373896*rdxF[1]*f[7]*nu[274]+1.677050983124842*rdxF[1]*f[18]*nu[254]+0.75*rdxF[1]*f[7]*nu[244]+2.5617376914899*rdxF[0]*f[42]*nu[112]+1.14564392373896*rdxF[0]*f[10]*nu[112]+1.677050983124842*rdxF[0]*f[18]*nu[92]+0.75*rdxF[0]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[0]*f[11]*nu[112]+1.677050983124842*rdxF[0]*f[19]*nu[92]+0.75*rdxF[0]*f[11]*nu[82]; 
  out[20] += 6.708203932499369*rdxF[0]*f[49]*nu[112]+7.6852130744697*rdxF[0]*f[5]*nu[112]+5.031152949374527*rdxF[0]*f[20]*nu[92]+3.75*rdxF[0]*f[1]*nu[92]+3.354101966249684*rdxF[0]*f[5]*nu[82]+1.677050983124842*rdxF[0]*f[1]*nu[80]; 
  out[22] += 6.708203932499369*rdxF[0]*f[51]*nu[112]+7.6852130744697*rdxF[0]*f[7]*nu[112]+5.031152949374527*rdxF[0]*f[22]*nu[92]+3.75*rdxF[0]*f[3]*nu[92]+3.354101966249684*rdxF[0]*f[7]*nu[82]+1.677050983124842*rdxF[0]*f[3]*nu[80]; 
  out[24] += 1.14564392373896*rdxF[0]*f[13]*nu[112]+1.677050983124842*rdxF[0]*f[24]*nu[92]+0.75*rdxF[0]*f[13]*nu[82]; 
  out[25] += 1.14564392373896*rdxF[1]*f[11]*nu[274]+1.677050983124842*rdxF[1]*f[25]*nu[254]+0.75*rdxF[1]*f[11]*nu[244]; 
  out[26] += 1.14564392373896*rdxF[1]*f[12]*nu[274]+1.677050983124842*rdxF[1]*f[26]*nu[254]+0.75*rdxF[1]*f[12]*nu[244]+6.708203932499369*rdxF[0]*f[55]*nu[112]+7.6852130744697*rdxF[0]*f[9]*nu[112]+5.031152949374527*rdxF[0]*f[26]*nu[92]+3.75*rdxF[0]*f[4]*nu[92]+3.354101966249684*rdxF[0]*f[9]*nu[82]+1.677050983124842*rdxF[0]*f[4]*nu[80]; 
  out[27] += 1.14564392373896*rdxF[1]*f[13]*nu[274]+1.677050983124842*rdxF[1]*f[27]*nu[254]+0.75*rdxF[1]*f[13]*nu[244]; 
  out[28] += 6.708203932499369*rdxF[1]*f[57]*nu[274]+7.6852130744697*rdxF[1]*f[8]*nu[274]+5.031152949374527*rdxF[1]*f[28]*nu[254]+3.75*f[1]*rdxF[1]*nu[254]+3.354101966249684*rdxF[1]*f[8]*nu[244]+1.677050983124842*f[1]*rdxF[1]*nu[240]; 
  out[29] += 6.708203932499369*rdxF[1]*f[58]*nu[274]+7.6852130744697*rdxF[1]*f[9]*nu[274]+5.031152949374527*rdxF[1]*f[29]*nu[254]+3.75*rdxF[1]*f[2]*nu[254]+3.354101966249684*rdxF[1]*f[9]*nu[244]+1.677050983124842*rdxF[1]*f[2]*nu[240]+1.14564392373896*rdxF[0]*f[14]*nu[112]+1.677050983124842*rdxF[0]*f[29]*nu[92]+0.75*rdxF[0]*f[14]*nu[82]; 
  out[30] += 6.708203932499369*rdxF[1]*f[59]*nu[274]+7.6852130744697*rdxF[1]*f[10]*nu[274]+5.031152949374527*rdxF[1]*f[30]*nu[254]+3.75*rdxF[1]*f[3]*nu[254]+3.354101966249684*rdxF[1]*f[10]*nu[244]+1.677050983124842*rdxF[1]*f[3]*nu[240]; 
  out[32] += 18.44756081437327*rdxF[0]*f[12]*nu[112]+10.5*f[0]*rdxF[0]*nu[112]+10.06230589874905*rdxF[0]*f[32]*nu[92]+12.8086884574495*rdxF[0]*f[2]*nu[92]+7.685213074469699*rdxF[0]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[0]*nu[82]+5.7282196186948*rdxF[0]*f[2]*nu[80]; 
  out[34] += 18.44756081437327*rdxF[1]*f[14]*nu[274]+10.5*f[0]*rdxF[1]*nu[274]+10.06230589874905*rdxF[1]*f[34]*nu[254]+12.8086884574495*rdxF[1]*f[4]*nu[254]+7.685213074469699*rdxF[1]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[1]*nu[244]+5.7282196186948*rdxF[1]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[1]*f[63]*nu[274]+1.14564392373896*rdxF[1]*f[15]*nu[274]+1.677050983124842*rdxF[1]*f[35]*nu[254]+0.75*rdxF[1]*f[15]*nu[244]+2.5617376914899*rdxF[0]*f[61]*nu[112]+1.14564392373896*rdxF[0]*f[17]*nu[112]+1.677050983124842*rdxF[0]*f[35]*nu[92]+0.75*rdxF[0]*f[17]*nu[82]; 
  out[36] += 1.14564392373896*rdxF[0]*f[21]*nu[112]+1.677050983124842*rdxF[0]*f[36]*nu[92]+0.75*rdxF[0]*f[21]*nu[82]; 
  out[37] += 6.708203932499369*rdxF[0]*f[65]*nu[112]+7.685213074469699*rdxF[0]*f[15]*nu[112]+5.031152949374527*rdxF[0]*f[37]*nu[92]+3.75*rdxF[0]*f[6]*nu[92]+3.354101966249685*rdxF[0]*f[15]*nu[82]+1.677050983124842*rdxF[0]*f[6]*nu[80]; 
  out[38] += 1.14564392373896*rdxF[0]*f[23]*nu[112]+1.677050983124842*rdxF[0]*f[38]*nu[92]+0.75*rdxF[0]*f[23]*nu[82]; 
  out[39] += 1.14564392373896*rdxF[1]*f[19]*nu[274]+1.677050983124842*rdxF[1]*f[39]*nu[254]+0.75*rdxF[1]*f[19]*nu[244]+1.14564392373896*rdxF[0]*f[25]*nu[112]+1.677050983124842*rdxF[0]*f[39]*nu[92]+0.75*rdxF[0]*f[25]*nu[82]; 
  out[40] += 1.14564392373896*rdxF[1]*f[20]*nu[274]+1.677050983124842*rdxF[1]*f[40]*nu[254]+0.75*rdxF[1]*f[20]*nu[244]+6.708203932499369*rdxF[0]*f[68]*nu[112]+7.685213074469699*rdxF[0]*f[16]*nu[112]+5.031152949374527*rdxF[0]*f[40]*nu[92]+3.75*rdxF[0]*f[8]*nu[92]+3.354101966249685*rdxF[0]*f[16]*nu[82]+1.677050983124842*rdxF[0]*f[8]*nu[80]; 
  out[41] += 1.14564392373896*rdxF[1]*f[21]*nu[274]+1.677050983124842*rdxF[1]*f[41]*nu[254]+0.75*rdxF[1]*f[21]*nu[244]; 
  out[42] += 1.14564392373896*rdxF[1]*f[22]*nu[274]+1.677050983124842*rdxF[1]*f[42]*nu[254]+0.75*rdxF[1]*f[22]*nu[244]+6.708203932499369*rdxF[0]*f[70]*nu[112]+7.685213074469699*rdxF[0]*f[18]*nu[112]+5.031152949374527*rdxF[0]*f[42]*nu[92]+3.75*rdxF[0]*f[10]*nu[92]+3.354101966249685*rdxF[0]*f[18]*nu[82]+1.677050983124842*rdxF[0]*f[10]*nu[80]; 
  out[43] += 1.14564392373896*rdxF[1]*f[23]*nu[274]+1.677050983124842*rdxF[1]*f[43]*nu[254]+0.75*rdxF[1]*f[23]*nu[244]; 
  out[44] += 1.14564392373896*rdxF[1]*f[24]*nu[274]+1.677050983124842*rdxF[1]*f[44]*nu[254]+0.75*rdxF[1]*f[24]*nu[244]+1.14564392373896*rdxF[0]*f[27]*nu[112]+1.677050983124842*rdxF[0]*f[44]*nu[92]+0.75*rdxF[0]*f[27]*nu[82]; 
  out[45] += 6.708203932499369*rdxF[1]*f[73]*nu[274]+7.685213074469699*rdxF[1]*f[16]*nu[274]+5.031152949374527*rdxF[1]*f[45]*nu[254]+3.75*rdxF[1]*f[5]*nu[254]+3.354101966249685*rdxF[1]*f[16]*nu[244]+1.677050983124842*rdxF[1]*f[5]*nu[240]+1.14564392373896*rdxF[0]*f[28]*nu[112]+1.677050983124842*rdxF[0]*f[45]*nu[92]+0.75*rdxF[0]*f[28]*nu[82]; 
  out[46] += 6.708203932499369*rdxF[1]*f[74]*nu[274]+7.685213074469699*rdxF[1]*f[17]*nu[274]+5.031152949374527*rdxF[1]*f[46]*nu[254]+3.75*rdxF[1]*f[6]*nu[254]+3.354101966249685*rdxF[1]*f[17]*nu[244]+1.677050983124842*rdxF[1]*f[6]*nu[240]; 
  out[47] += 6.708203932499369*rdxF[1]*f[75]*nu[274]+7.685213074469699*rdxF[1]*f[18]*nu[274]+5.031152949374527*rdxF[1]*f[47]*nu[254]+3.75*rdxF[1]*f[7]*nu[254]+3.354101966249685*rdxF[1]*f[18]*nu[244]+1.677050983124842*rdxF[1]*f[7]*nu[240]+1.14564392373896*rdxF[0]*f[30]*nu[112]+1.677050983124842*rdxF[0]*f[47]*nu[92]+0.75*rdxF[0]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[0]*f[31]*nu[112]+1.677050983124842*rdxF[0]*f[48]*nu[92]+0.7499999999999999*rdxF[0]*f[31]*nu[82]; 
  out[49] += 18.44756081437327*rdxF[0]*f[20]*nu[112]+10.5*rdxF[0]*f[1]*nu[112]+10.06230589874905*rdxF[0]*f[49]*nu[92]+12.8086884574495*rdxF[0]*f[5]*nu[92]+7.685213074469698*rdxF[0]*f[20]*nu[82]+6.873863542433758*rdxF[0]*f[1]*nu[82]+5.7282196186948*rdxF[0]*f[5]*nu[80]; 
  out[51] += 18.44756081437327*rdxF[0]*f[22]*nu[112]+10.5*rdxF[0]*f[3]*nu[112]+10.06230589874905*rdxF[0]*f[51]*nu[92]+12.8086884574495*rdxF[0]*f[7]*nu[92]+7.685213074469698*rdxF[0]*f[22]*nu[82]+6.873863542433758*rdxF[0]*f[3]*nu[82]+5.7282196186948*rdxF[0]*f[7]*nu[80]; 
  out[53] += 1.14564392373896*rdxF[0]*f[33]*nu[112]+1.677050983124842*rdxF[0]*f[53]*nu[92]+0.7499999999999999*rdxF[0]*f[33]*nu[82]; 
  out[54] += 1.14564392373896*rdxF[1]*f[31]*nu[274]+1.677050983124842*rdxF[1]*f[54]*nu[254]+0.7499999999999999*rdxF[1]*f[31]*nu[244]; 
  out[55] += 1.14564392373896*rdxF[1]*f[32]*nu[274]+1.677050983124842*rdxF[1]*f[55]*nu[254]+0.7499999999999999*rdxF[1]*f[32]*nu[244]+18.44756081437327*rdxF[0]*f[26]*nu[112]+10.5*rdxF[0]*f[4]*nu[112]+10.06230589874905*rdxF[0]*f[55]*nu[92]+12.8086884574495*rdxF[0]*f[9]*nu[92]+7.685213074469698*rdxF[0]*f[26]*nu[82]+6.873863542433758*rdxF[0]*f[4]*nu[82]+5.7282196186948*rdxF[0]*f[9]*nu[80]; 
  out[56] += 1.14564392373896*rdxF[1]*f[33]*nu[274]+1.677050983124842*rdxF[1]*f[56]*nu[254]+0.7499999999999999*rdxF[1]*f[33]*nu[244]; 
  out[57] += 18.44756081437327*rdxF[1]*f[28]*nu[274]+10.5*f[1]*rdxF[1]*nu[274]+10.06230589874905*rdxF[1]*f[57]*nu[254]+12.8086884574495*rdxF[1]*f[8]*nu[254]+7.685213074469698*rdxF[1]*f[28]*nu[244]+6.873863542433758*f[1]*rdxF[1]*nu[244]+5.7282196186948*rdxF[1]*f[8]*nu[240]; 
  out[58] += 18.44756081437327*rdxF[1]*f[29]*nu[274]+10.5*rdxF[1]*f[2]*nu[274]+10.06230589874905*rdxF[1]*f[58]*nu[254]+12.8086884574495*rdxF[1]*f[9]*nu[254]+7.685213074469698*rdxF[1]*f[29]*nu[244]+6.873863542433758*rdxF[1]*f[2]*nu[244]+5.7282196186948*rdxF[1]*f[9]*nu[240]+1.14564392373896*rdxF[0]*f[34]*nu[112]+1.677050983124842*rdxF[0]*f[58]*nu[92]+0.7499999999999999*rdxF[0]*f[34]*nu[82]; 
  out[59] += 18.44756081437327*rdxF[1]*f[30]*nu[274]+10.5*rdxF[1]*f[3]*nu[274]+10.06230589874905*rdxF[1]*f[59]*nu[254]+12.8086884574495*rdxF[1]*f[10]*nu[254]+7.685213074469698*rdxF[1]*f[30]*nu[244]+6.873863542433758*rdxF[1]*f[3]*nu[244]+5.7282196186948*rdxF[1]*f[10]*nu[240]; 
  out[60] += 1.14564392373896*rdxF[1]*f[36]*nu[274]+1.677050983124842*rdxF[1]*f[60]*nu[254]+0.75*rdxF[1]*f[36]*nu[244]+1.14564392373896*rdxF[0]*f[41]*nu[112]+1.677050983124842*rdxF[0]*f[60]*nu[92]+0.75*rdxF[0]*f[41]*nu[82]; 
  out[61] += 1.14564392373896*rdxF[1]*f[37]*nu[274]+1.677050983124842*rdxF[1]*f[61]*nu[254]+0.75*rdxF[1]*f[37]*nu[244]+6.708203932499369*rdxF[0]*f[77]*nu[112]+7.6852130744697*rdxF[0]*f[35]*nu[112]+5.031152949374527*rdxF[0]*f[61]*nu[92]+3.75*rdxF[0]*f[17]*nu[92]+3.354101966249684*rdxF[0]*f[35]*nu[82]+1.677050983124842*rdxF[0]*f[17]*nu[80]; 
  out[62] += 1.14564392373896*rdxF[1]*f[38]*nu[274]+1.677050983124842*rdxF[1]*f[62]*nu[254]+0.75*rdxF[1]*f[38]*nu[244]+1.14564392373896*rdxF[0]*f[43]*nu[112]+1.677050983124842*rdxF[0]*f[62]*nu[92]+0.75*rdxF[0]*f[43]*nu[82]; 
  out[63] += 6.708203932499369*rdxF[1]*f[79]*nu[274]+7.6852130744697*rdxF[1]*f[35]*nu[274]+5.031152949374527*rdxF[1]*f[63]*nu[254]+3.75*rdxF[1]*f[15]*nu[254]+3.354101966249684*rdxF[1]*f[35]*nu[244]+1.677050983124842*rdxF[1]*f[15]*nu[240]+1.14564392373896*rdxF[0]*f[46]*nu[112]+1.677050983124842*rdxF[0]*f[63]*nu[92]+0.75*rdxF[0]*f[46]*nu[82]; 
  out[64] += 1.14564392373896*rdxF[0]*f[50]*nu[112]+1.677050983124842*rdxF[0]*f[64]*nu[92]+0.7499999999999999*rdxF[0]*f[50]*nu[82]; 
  out[65] += 18.44756081437327*rdxF[0]*f[37]*nu[112]+10.5*rdxF[0]*f[6]*nu[112]+10.06230589874905*rdxF[0]*f[65]*nu[92]+12.8086884574495*rdxF[0]*f[15]*nu[92]+7.685213074469699*rdxF[0]*f[37]*nu[82]+6.873863542433759*rdxF[0]*f[6]*nu[82]+5.7282196186948*rdxF[0]*f[15]*nu[80]; 
  out[66] += 1.14564392373896*rdxF[0]*f[52]*nu[112]+1.677050983124842*rdxF[0]*f[66]*nu[92]+0.7499999999999999*rdxF[0]*f[52]*nu[82]; 
  out[67] += 1.14564392373896*rdxF[1]*f[48]*nu[274]+1.677050983124842*rdxF[1]*f[67]*nu[254]+0.7499999999999999*rdxF[1]*f[48]*nu[244]+1.14564392373896*rdxF[0]*f[54]*nu[112]+1.677050983124842*rdxF[0]*f[67]*nu[92]+0.7499999999999999*rdxF[0]*f[54]*nu[82]; 
  out[68] += 1.14564392373896*rdxF[1]*f[49]*nu[274]+1.677050983124842*rdxF[1]*f[68]*nu[254]+0.7499999999999999*rdxF[1]*f[49]*nu[244]+18.44756081437327*rdxF[0]*f[40]*nu[112]+10.5*rdxF[0]*f[8]*nu[112]+10.06230589874905*rdxF[0]*f[68]*nu[92]+12.8086884574495*rdxF[0]*f[16]*nu[92]+7.685213074469699*rdxF[0]*f[40]*nu[82]+6.873863542433759*rdxF[0]*f[8]*nu[82]+5.7282196186948*rdxF[0]*f[16]*nu[80]; 
  out[69] += 1.14564392373896*rdxF[1]*f[50]*nu[274]+1.677050983124842*rdxF[1]*f[69]*nu[254]+0.7499999999999999*rdxF[1]*f[50]*nu[244]; 
  out[70] += 1.14564392373896*rdxF[1]*f[51]*nu[274]+1.677050983124842*rdxF[1]*f[70]*nu[254]+0.7499999999999999*rdxF[1]*f[51]*nu[244]+18.44756081437327*rdxF[0]*f[42]*nu[112]+10.5*rdxF[0]*f[10]*nu[112]+10.06230589874905*rdxF[0]*f[70]*nu[92]+12.8086884574495*rdxF[0]*f[18]*nu[92]+7.685213074469699*rdxF[0]*f[42]*nu[82]+6.873863542433759*rdxF[0]*f[10]*nu[82]+5.7282196186948*rdxF[0]*f[18]*nu[80]; 
  out[71] += 1.14564392373896*rdxF[1]*f[52]*nu[274]+1.677050983124842*rdxF[1]*f[71]*nu[254]+0.7499999999999999*rdxF[1]*f[52]*nu[244]; 
  out[72] += 1.14564392373896*rdxF[1]*f[53]*nu[274]+1.677050983124842*rdxF[1]*f[72]*nu[254]+0.7499999999999999*rdxF[1]*f[53]*nu[244]+1.14564392373896*rdxF[0]*f[56]*nu[112]+1.677050983124842*rdxF[0]*f[72]*nu[92]+0.7499999999999999*rdxF[0]*f[56]*nu[82]; 
  out[73] += 18.44756081437327*rdxF[1]*f[45]*nu[274]+10.5*rdxF[1]*f[5]*nu[274]+10.06230589874905*rdxF[1]*f[73]*nu[254]+12.8086884574495*rdxF[1]*f[16]*nu[254]+7.685213074469699*rdxF[1]*f[45]*nu[244]+6.873863542433759*rdxF[1]*f[5]*nu[244]+5.7282196186948*rdxF[1]*f[16]*nu[240]+1.14564392373896*rdxF[0]*f[57]*nu[112]+1.677050983124842*rdxF[0]*f[73]*nu[92]+0.7499999999999999*rdxF[0]*f[57]*nu[82]; 
  out[74] += 18.44756081437327*rdxF[1]*f[46]*nu[274]+10.5*rdxF[1]*f[6]*nu[274]+10.06230589874905*rdxF[1]*f[74]*nu[254]+12.8086884574495*rdxF[1]*f[17]*nu[254]+7.685213074469699*rdxF[1]*f[46]*nu[244]+6.873863542433759*rdxF[1]*f[6]*nu[244]+5.7282196186948*rdxF[1]*f[17]*nu[240]; 
  out[75] += 18.44756081437327*rdxF[1]*f[47]*nu[274]+10.5*rdxF[1]*f[7]*nu[274]+10.06230589874905*rdxF[1]*f[75]*nu[254]+12.8086884574495*rdxF[1]*f[18]*nu[254]+7.685213074469699*rdxF[1]*f[47]*nu[244]+6.873863542433759*rdxF[1]*f[7]*nu[244]+5.7282196186948*rdxF[1]*f[18]*nu[240]+1.14564392373896*rdxF[0]*f[59]*nu[112]+1.677050983124842*rdxF[0]*f[75]*nu[92]+0.7499999999999999*rdxF[0]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[1]*f[64]*nu[274]+1.677050983124842*rdxF[1]*f[76]*nu[254]+0.7499999999999999*rdxF[1]*f[64]*nu[244]+1.14564392373896*rdxF[0]*f[69]*nu[112]+1.677050983124842*rdxF[0]*f[76]*nu[92]+0.7499999999999999*rdxF[0]*f[69]*nu[82]; 
  out[77] += 1.14564392373896*rdxF[1]*f[65]*nu[274]+1.677050983124842*rdxF[1]*f[77]*nu[254]+0.7499999999999999*rdxF[1]*f[65]*nu[244]+18.44756081437327*rdxF[0]*f[61]*nu[112]+10.5*rdxF[0]*f[17]*nu[112]+10.06230589874905*rdxF[0]*f[77]*nu[92]+12.8086884574495*rdxF[0]*f[35]*nu[92]+7.685213074469698*rdxF[0]*f[61]*nu[82]+6.873863542433758*rdxF[0]*f[17]*nu[82]+5.7282196186948*rdxF[0]*f[35]*nu[80]; 
  out[78] += 1.14564392373896*rdxF[1]*f[66]*nu[274]+1.677050983124842*rdxF[1]*f[78]*nu[254]+0.7499999999999999*rdxF[1]*f[66]*nu[244]+1.14564392373896*rdxF[0]*f[71]*nu[112]+1.677050983124842*rdxF[0]*f[78]*nu[92]+0.7499999999999999*rdxF[0]*f[71]*nu[82]; 
  out[79] += 18.44756081437327*rdxF[1]*f[63]*nu[274]+10.5*rdxF[1]*f[15]*nu[274]+10.06230589874905*rdxF[1]*f[79]*nu[254]+12.8086884574495*rdxF[1]*f[35]*nu[254]+7.685213074469698*rdxF[1]*f[63]*nu[244]+6.873863542433758*rdxF[1]*f[15]*nu[244]+5.7282196186948*rdxF[1]*f[35]*nu[240]+1.14564392373896*rdxF[0]*f[74]*nu[112]+1.677050983124842*rdxF[0]*f[79]*nu[92]+0.7499999999999999*rdxF[0]*f[74]*nu[82]; 

  return (rdxF[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+rdxF[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs124(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 2.5617376914899*rdxF[1]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[2]*nu[92]+0.75*f[0]*rdxF[1]*nu[82]; 
  out[4] += 2.5617376914899*rdxF[2]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[4]*nu[254]+0.75*f[0]*rdxF[2]*nu[244]; 
  out[5] += 2.5617376914899*rdxF[1]*f[20]*nu[112]+1.14564392373896*f[1]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[5]*nu[92]+0.75*f[1]*rdxF[1]*nu[82]+2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 2.5617376914899*rdxF[1]*f[22]*nu[112]+1.14564392373896*rdxF[1]*f[3]*nu[112]+1.677050983124842*rdxF[1]*f[7]*nu[92]+0.75*rdxF[1]*f[3]*nu[82]; 
  out[8] += 2.5617376914899*rdxF[2]*f[28]*nu[274]+1.14564392373896*f[1]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[8]*nu[254]+0.75*f[1]*rdxF[2]*nu[244]+2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 2.5617376914899*rdxF[2]*f[29]*nu[274]+1.14564392373896*f[2]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[9]*nu[254]+0.75*f[2]*rdxF[2]*nu[244]+2.5617376914899*rdxF[1]*f[26]*nu[112]+1.14564392373896*rdxF[1]*f[4]*nu[112]+1.677050983124842*rdxF[1]*f[9]*nu[92]+0.75*rdxF[1]*f[4]*nu[82]; 
  out[10] += 2.5617376914899*rdxF[2]*f[30]*nu[274]+1.14564392373896*rdxF[2]*f[3]*nu[274]+1.677050983124842*rdxF[2]*f[10]*nu[254]+0.75*rdxF[2]*f[3]*nu[244]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 6.708203932499369*rdxF[1]*f[32]*nu[112]+7.685213074469699*rdxF[1]*f[2]*nu[112]+5.031152949374527*rdxF[1]*f[12]*nu[92]+3.75*f[0]*rdxF[1]*nu[92]+3.354101966249685*rdxF[1]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[1]*nu[80]; 
  out[14] += 6.708203932499369*rdxF[2]*f[34]*nu[274]+7.685213074469699*rdxF[2]*f[4]*nu[274]+5.031152949374527*rdxF[2]*f[14]*nu[254]+3.75*f[0]*rdxF[2]*nu[254]+3.354101966249685*rdxF[2]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[2]*nu[240]; 
  out[15] += 2.5617376914899*rdxF[1]*f[37]*nu[112]+1.14564392373896*rdxF[1]*f[6]*nu[112]+1.677050983124842*rdxF[1]*f[15]*nu[92]+0.75*rdxF[1]*f[6]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[2]*f[45]*nu[274]+1.14564392373896*rdxF[2]*f[5]*nu[274]+1.677050983124842*rdxF[2]*f[16]*nu[254]+0.75*rdxF[2]*f[5]*nu[244]+2.5617376914899*rdxF[1]*f[40]*nu[112]+1.14564392373896*rdxF[1]*f[8]*nu[112]+1.677050983124842*rdxF[1]*f[16]*nu[92]+0.75*rdxF[1]*f[8]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[2]*f[46]*nu[274]+1.14564392373896*rdxF[2]*f[6]*nu[274]+1.677050983124842*rdxF[2]*f[17]*nu[254]+0.75*rdxF[2]*f[6]*nu[244]+2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 2.5617376914899*rdxF[2]*f[47]*nu[274]+1.14564392373896*rdxF[2]*f[7]*nu[274]+1.677050983124842*rdxF[2]*f[18]*nu[254]+0.75*rdxF[2]*f[7]*nu[244]+2.5617376914899*rdxF[1]*f[42]*nu[112]+1.14564392373896*rdxF[1]*f[10]*nu[112]+1.677050983124842*rdxF[1]*f[18]*nu[92]+0.75*rdxF[1]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[1]*f[11]*nu[112]+1.677050983124842*rdxF[1]*f[19]*nu[92]+0.75*rdxF[1]*f[11]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 6.708203932499369*rdxF[1]*f[49]*nu[112]+7.6852130744697*rdxF[1]*f[5]*nu[112]+5.031152949374527*rdxF[1]*f[20]*nu[92]+3.75*f[1]*rdxF[1]*nu[92]+3.354101966249684*rdxF[1]*f[5]*nu[82]+1.677050983124842*f[1]*rdxF[1]*nu[80]+1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 6.708203932499369*rdxF[1]*f[51]*nu[112]+7.6852130744697*rdxF[1]*f[7]*nu[112]+5.031152949374527*rdxF[1]*f[22]*nu[92]+3.75*rdxF[1]*f[3]*nu[92]+3.354101966249684*rdxF[1]*f[7]*nu[82]+1.677050983124842*rdxF[1]*f[3]*nu[80]; 
  out[23] += 1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 1.14564392373896*rdxF[1]*f[13]*nu[112]+1.677050983124842*rdxF[1]*f[24]*nu[92]+0.75*rdxF[1]*f[13]*nu[82]; 
  out[25] += 1.14564392373896*rdxF[2]*f[11]*nu[274]+1.677050983124842*rdxF[2]*f[25]*nu[254]+0.75*rdxF[2]*f[11]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.14564392373896*rdxF[2]*f[12]*nu[274]+1.677050983124842*rdxF[2]*f[26]*nu[254]+0.75*rdxF[2]*f[12]*nu[244]+6.708203932499369*rdxF[1]*f[55]*nu[112]+7.6852130744697*rdxF[1]*f[9]*nu[112]+5.031152949374527*rdxF[1]*f[26]*nu[92]+3.75*rdxF[1]*f[4]*nu[92]+3.354101966249684*rdxF[1]*f[9]*nu[82]+1.677050983124842*rdxF[1]*f[4]*nu[80]; 
  out[27] += 1.14564392373896*rdxF[2]*f[13]*nu[274]+1.677050983124842*rdxF[2]*f[27]*nu[254]+0.75*rdxF[2]*f[13]*nu[244]; 
  out[28] += 6.708203932499369*rdxF[2]*f[57]*nu[274]+7.6852130744697*rdxF[2]*f[8]*nu[274]+5.031152949374527*rdxF[2]*f[28]*nu[254]+3.75*f[1]*rdxF[2]*nu[254]+3.354101966249684*rdxF[2]*f[8]*nu[244]+1.677050983124842*f[1]*rdxF[2]*nu[240]+1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 6.708203932499369*rdxF[2]*f[58]*nu[274]+7.6852130744697*rdxF[2]*f[9]*nu[274]+5.031152949374527*rdxF[2]*f[29]*nu[254]+3.75*f[2]*rdxF[2]*nu[254]+3.354101966249684*rdxF[2]*f[9]*nu[244]+1.677050983124842*f[2]*rdxF[2]*nu[240]+1.14564392373896*rdxF[1]*f[14]*nu[112]+1.677050983124842*rdxF[1]*f[29]*nu[92]+0.75*rdxF[1]*f[14]*nu[82]; 
  out[30] += 6.708203932499369*rdxF[2]*f[59]*nu[274]+7.6852130744697*rdxF[2]*f[10]*nu[274]+5.031152949374527*rdxF[2]*f[30]*nu[254]+3.75*rdxF[2]*f[3]*nu[254]+3.354101966249684*rdxF[2]*f[10]*nu[244]+1.677050983124842*rdxF[2]*f[3]*nu[240]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[32] += 18.44756081437327*rdxF[1]*f[12]*nu[112]+10.5*f[0]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[32]*nu[92]+12.8086884574495*rdxF[1]*f[2]*nu[92]+7.685213074469699*rdxF[1]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[2]*nu[80]; 
  out[34] += 18.44756081437327*rdxF[2]*f[14]*nu[274]+10.5*f[0]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[34]*nu[254]+12.8086884574495*rdxF[2]*f[4]*nu[254]+7.685213074469699*rdxF[2]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[2]*f[63]*nu[274]+1.14564392373896*rdxF[2]*f[15]*nu[274]+1.677050983124842*rdxF[2]*f[35]*nu[254]+0.75*rdxF[2]*f[15]*nu[244]+2.5617376914899*rdxF[1]*f[61]*nu[112]+1.14564392373896*rdxF[1]*f[17]*nu[112]+1.677050983124842*rdxF[1]*f[35]*nu[92]+0.75*rdxF[1]*f[17]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 1.14564392373896*rdxF[1]*f[21]*nu[112]+1.677050983124842*rdxF[1]*f[36]*nu[92]+0.75*rdxF[1]*f[21]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 6.708203932499369*rdxF[1]*f[65]*nu[112]+7.685213074469699*rdxF[1]*f[15]*nu[112]+5.031152949374527*rdxF[1]*f[37]*nu[92]+3.75*rdxF[1]*f[6]*nu[92]+3.354101966249685*rdxF[1]*f[15]*nu[82]+1.677050983124842*rdxF[1]*f[6]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 1.14564392373896*rdxF[1]*f[23]*nu[112]+1.677050983124842*rdxF[1]*f[38]*nu[92]+0.75*rdxF[1]*f[23]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 1.14564392373896*rdxF[2]*f[19]*nu[274]+1.677050983124842*rdxF[2]*f[39]*nu[254]+0.75*rdxF[2]*f[19]*nu[244]+1.14564392373896*rdxF[1]*f[25]*nu[112]+1.677050983124842*rdxF[1]*f[39]*nu[92]+0.75*rdxF[1]*f[25]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 1.14564392373896*rdxF[2]*f[20]*nu[274]+1.677050983124842*rdxF[2]*f[40]*nu[254]+0.75*rdxF[2]*f[20]*nu[244]+6.708203932499369*rdxF[1]*f[68]*nu[112]+7.685213074469699*rdxF[1]*f[16]*nu[112]+5.031152949374527*rdxF[1]*f[40]*nu[92]+3.75*rdxF[1]*f[8]*nu[92]+3.354101966249685*rdxF[1]*f[16]*nu[82]+1.677050983124842*rdxF[1]*f[8]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 1.14564392373896*rdxF[2]*f[21]*nu[274]+1.677050983124842*rdxF[2]*f[41]*nu[254]+0.75*rdxF[2]*f[21]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[42] += 1.14564392373896*rdxF[2]*f[22]*nu[274]+1.677050983124842*rdxF[2]*f[42]*nu[254]+0.75*rdxF[2]*f[22]*nu[244]+6.708203932499369*rdxF[1]*f[70]*nu[112]+7.685213074469699*rdxF[1]*f[18]*nu[112]+5.031152949374527*rdxF[1]*f[42]*nu[92]+3.75*rdxF[1]*f[10]*nu[92]+3.354101966249685*rdxF[1]*f[18]*nu[82]+1.677050983124842*rdxF[1]*f[10]*nu[80]; 
  out[43] += 1.14564392373896*rdxF[2]*f[23]*nu[274]+1.677050983124842*rdxF[2]*f[43]*nu[254]+0.75*rdxF[2]*f[23]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[44] += 1.14564392373896*rdxF[2]*f[24]*nu[274]+1.677050983124842*rdxF[2]*f[44]*nu[254]+0.75*rdxF[2]*f[24]*nu[244]+1.14564392373896*rdxF[1]*f[27]*nu[112]+1.677050983124842*rdxF[1]*f[44]*nu[92]+0.75*rdxF[1]*f[27]*nu[82]; 
  out[45] += 6.708203932499369*rdxF[2]*f[73]*nu[274]+7.685213074469699*rdxF[2]*f[16]*nu[274]+5.031152949374527*rdxF[2]*f[45]*nu[254]+3.75*rdxF[2]*f[5]*nu[254]+3.354101966249685*rdxF[2]*f[16]*nu[244]+1.677050983124842*rdxF[2]*f[5]*nu[240]+1.14564392373896*rdxF[1]*f[28]*nu[112]+1.677050983124842*rdxF[1]*f[45]*nu[92]+0.75*rdxF[1]*f[28]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 6.708203932499369*rdxF[2]*f[74]*nu[274]+7.685213074469699*rdxF[2]*f[17]*nu[274]+5.031152949374527*rdxF[2]*f[46]*nu[254]+3.75*rdxF[2]*f[6]*nu[254]+3.354101966249685*rdxF[2]*f[17]*nu[244]+1.677050983124842*rdxF[2]*f[6]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[47] += 6.708203932499369*rdxF[2]*f[75]*nu[274]+7.685213074469699*rdxF[2]*f[18]*nu[274]+5.031152949374527*rdxF[2]*f[47]*nu[254]+3.75*rdxF[2]*f[7]*nu[254]+3.354101966249685*rdxF[2]*f[18]*nu[244]+1.677050983124842*rdxF[2]*f[7]*nu[240]+1.14564392373896*rdxF[1]*f[30]*nu[112]+1.677050983124842*rdxF[1]*f[47]*nu[92]+0.75*rdxF[1]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[1]*f[31]*nu[112]+1.677050983124842*rdxF[1]*f[48]*nu[92]+0.7499999999999999*rdxF[1]*f[31]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 18.44756081437327*rdxF[1]*f[20]*nu[112]+10.5*f[1]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[49]*nu[92]+12.8086884574495*rdxF[1]*f[5]*nu[92]+7.685213074469698*rdxF[1]*f[20]*nu[82]+6.873863542433758*f[1]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[5]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[51] += 18.44756081437327*rdxF[1]*f[22]*nu[112]+10.5*rdxF[1]*f[3]*nu[112]+10.06230589874905*rdxF[1]*f[51]*nu[92]+12.8086884574495*rdxF[1]*f[7]*nu[92]+7.685213074469698*rdxF[1]*f[22]*nu[82]+6.873863542433758*rdxF[1]*f[3]*nu[82]+5.7282196186948*rdxF[1]*f[7]*nu[80]; 
  out[52] += 1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[53] += 1.14564392373896*rdxF[1]*f[33]*nu[112]+1.677050983124842*rdxF[1]*f[53]*nu[92]+0.7499999999999999*rdxF[1]*f[33]*nu[82]; 
  out[54] += 1.14564392373896*rdxF[2]*f[31]*nu[274]+1.677050983124842*rdxF[2]*f[54]*nu[254]+0.7499999999999999*rdxF[2]*f[31]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[55] += 1.14564392373896*rdxF[2]*f[32]*nu[274]+1.677050983124842*rdxF[2]*f[55]*nu[254]+0.7499999999999999*rdxF[2]*f[32]*nu[244]+18.44756081437327*rdxF[1]*f[26]*nu[112]+10.5*rdxF[1]*f[4]*nu[112]+10.06230589874905*rdxF[1]*f[55]*nu[92]+12.8086884574495*rdxF[1]*f[9]*nu[92]+7.685213074469698*rdxF[1]*f[26]*nu[82]+6.873863542433758*rdxF[1]*f[4]*nu[82]+5.7282196186948*rdxF[1]*f[9]*nu[80]; 
  out[56] += 1.14564392373896*rdxF[2]*f[33]*nu[274]+1.677050983124842*rdxF[2]*f[56]*nu[254]+0.7499999999999999*rdxF[2]*f[33]*nu[244]; 
  out[57] += 18.44756081437327*rdxF[2]*f[28]*nu[274]+10.5*f[1]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[57]*nu[254]+12.8086884574495*rdxF[2]*f[8]*nu[254]+7.685213074469698*rdxF[2]*f[28]*nu[244]+6.873863542433758*f[1]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[8]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[58] += 18.44756081437327*rdxF[2]*f[29]*nu[274]+10.5*f[2]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[58]*nu[254]+12.8086884574495*rdxF[2]*f[9]*nu[254]+7.685213074469698*rdxF[2]*f[29]*nu[244]+6.873863542433758*f[2]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[9]*nu[240]+1.14564392373896*rdxF[1]*f[34]*nu[112]+1.677050983124842*rdxF[1]*f[58]*nu[92]+0.7499999999999999*rdxF[1]*f[34]*nu[82]; 
  out[59] += 18.44756081437327*rdxF[2]*f[30]*nu[274]+10.5*rdxF[2]*f[3]*nu[274]+10.06230589874905*rdxF[2]*f[59]*nu[254]+12.8086884574495*rdxF[2]*f[10]*nu[254]+7.685213074469698*rdxF[2]*f[30]*nu[244]+6.873863542433758*rdxF[2]*f[3]*nu[244]+5.7282196186948*rdxF[2]*f[10]*nu[240]; 
  out[60] += 1.14564392373896*rdxF[2]*f[36]*nu[274]+1.677050983124842*rdxF[2]*f[60]*nu[254]+0.75*rdxF[2]*f[36]*nu[244]+1.14564392373896*rdxF[1]*f[41]*nu[112]+1.677050983124842*rdxF[1]*f[60]*nu[92]+0.75*rdxF[1]*f[41]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 1.14564392373896*rdxF[2]*f[37]*nu[274]+1.677050983124842*rdxF[2]*f[61]*nu[254]+0.75*rdxF[2]*f[37]*nu[244]+6.708203932499369*rdxF[1]*f[77]*nu[112]+7.6852130744697*rdxF[1]*f[35]*nu[112]+5.031152949374527*rdxF[1]*f[61]*nu[92]+3.75*rdxF[1]*f[17]*nu[92]+3.354101966249684*rdxF[1]*f[35]*nu[82]+1.677050983124842*rdxF[1]*f[17]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 1.14564392373896*rdxF[2]*f[38]*nu[274]+1.677050983124842*rdxF[2]*f[62]*nu[254]+0.75*rdxF[2]*f[38]*nu[244]+1.14564392373896*rdxF[1]*f[43]*nu[112]+1.677050983124842*rdxF[1]*f[62]*nu[92]+0.75*rdxF[1]*f[43]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 6.708203932499369*rdxF[2]*f[79]*nu[274]+7.6852130744697*rdxF[2]*f[35]*nu[274]+5.031152949374527*rdxF[2]*f[63]*nu[254]+3.75*rdxF[2]*f[15]*nu[254]+3.354101966249684*rdxF[2]*f[35]*nu[244]+1.677050983124842*rdxF[2]*f[15]*nu[240]+1.14564392373896*rdxF[1]*f[46]*nu[112]+1.677050983124842*rdxF[1]*f[63]*nu[92]+0.75*rdxF[1]*f[46]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 1.14564392373896*rdxF[1]*f[50]*nu[112]+1.677050983124842*rdxF[1]*f[64]*nu[92]+0.7499999999999999*rdxF[1]*f[50]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 18.44756081437327*rdxF[1]*f[37]*nu[112]+10.5*rdxF[1]*f[6]*nu[112]+10.06230589874905*rdxF[1]*f[65]*nu[92]+12.8086884574495*rdxF[1]*f[15]*nu[92]+7.685213074469699*rdxF[1]*f[37]*nu[82]+6.873863542433759*rdxF[1]*f[6]*nu[82]+5.7282196186948*rdxF[1]*f[15]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 1.14564392373896*rdxF[1]*f[52]*nu[112]+1.677050983124842*rdxF[1]*f[66]*nu[92]+0.7499999999999999*rdxF[1]*f[52]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 1.14564392373896*rdxF[2]*f[48]*nu[274]+1.677050983124842*rdxF[2]*f[67]*nu[254]+0.7499999999999999*rdxF[2]*f[48]*nu[244]+1.14564392373896*rdxF[1]*f[54]*nu[112]+1.677050983124842*rdxF[1]*f[67]*nu[92]+0.7499999999999999*rdxF[1]*f[54]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 1.14564392373896*rdxF[2]*f[49]*nu[274]+1.677050983124842*rdxF[2]*f[68]*nu[254]+0.7499999999999999*rdxF[2]*f[49]*nu[244]+18.44756081437327*rdxF[1]*f[40]*nu[112]+10.5*rdxF[1]*f[8]*nu[112]+10.06230589874905*rdxF[1]*f[68]*nu[92]+12.8086884574495*rdxF[1]*f[16]*nu[92]+7.685213074469699*rdxF[1]*f[40]*nu[82]+6.873863542433759*rdxF[1]*f[8]*nu[82]+5.7282196186948*rdxF[1]*f[16]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 1.14564392373896*rdxF[2]*f[50]*nu[274]+1.677050983124842*rdxF[2]*f[69]*nu[254]+0.7499999999999999*rdxF[2]*f[50]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[70] += 1.14564392373896*rdxF[2]*f[51]*nu[274]+1.677050983124842*rdxF[2]*f[70]*nu[254]+0.7499999999999999*rdxF[2]*f[51]*nu[244]+18.44756081437327*rdxF[1]*f[42]*nu[112]+10.5*rdxF[1]*f[10]*nu[112]+10.06230589874905*rdxF[1]*f[70]*nu[92]+12.8086884574495*rdxF[1]*f[18]*nu[92]+7.685213074469699*rdxF[1]*f[42]*nu[82]+6.873863542433759*rdxF[1]*f[10]*nu[82]+5.7282196186948*rdxF[1]*f[18]*nu[80]; 
  out[71] += 1.14564392373896*rdxF[2]*f[52]*nu[274]+1.677050983124842*rdxF[2]*f[71]*nu[254]+0.7499999999999999*rdxF[2]*f[52]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[72] += 1.14564392373896*rdxF[2]*f[53]*nu[274]+1.677050983124842*rdxF[2]*f[72]*nu[254]+0.7499999999999999*rdxF[2]*f[53]*nu[244]+1.14564392373896*rdxF[1]*f[56]*nu[112]+1.677050983124842*rdxF[1]*f[72]*nu[92]+0.7499999999999999*rdxF[1]*f[56]*nu[82]; 
  out[73] += 18.44756081437327*rdxF[2]*f[45]*nu[274]+10.5*rdxF[2]*f[5]*nu[274]+10.06230589874905*rdxF[2]*f[73]*nu[254]+12.8086884574495*rdxF[2]*f[16]*nu[254]+7.685213074469699*rdxF[2]*f[45]*nu[244]+6.873863542433759*rdxF[2]*f[5]*nu[244]+5.7282196186948*rdxF[2]*f[16]*nu[240]+1.14564392373896*rdxF[1]*f[57]*nu[112]+1.677050983124842*rdxF[1]*f[73]*nu[92]+0.7499999999999999*rdxF[1]*f[57]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 18.44756081437327*rdxF[2]*f[46]*nu[274]+10.5*rdxF[2]*f[6]*nu[274]+10.06230589874905*rdxF[2]*f[74]*nu[254]+12.8086884574495*rdxF[2]*f[17]*nu[254]+7.685213074469699*rdxF[2]*f[46]*nu[244]+6.873863542433759*rdxF[2]*f[6]*nu[244]+5.7282196186948*rdxF[2]*f[17]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[75] += 18.44756081437327*rdxF[2]*f[47]*nu[274]+10.5*rdxF[2]*f[7]*nu[274]+10.06230589874905*rdxF[2]*f[75]*nu[254]+12.8086884574495*rdxF[2]*f[18]*nu[254]+7.685213074469699*rdxF[2]*f[47]*nu[244]+6.873863542433759*rdxF[2]*f[7]*nu[244]+5.7282196186948*rdxF[2]*f[18]*nu[240]+1.14564392373896*rdxF[1]*f[59]*nu[112]+1.677050983124842*rdxF[1]*f[75]*nu[92]+0.7499999999999999*rdxF[1]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[2]*f[64]*nu[274]+1.677050983124842*rdxF[2]*f[76]*nu[254]+0.7499999999999999*rdxF[2]*f[64]*nu[244]+1.14564392373896*rdxF[1]*f[69]*nu[112]+1.677050983124842*rdxF[1]*f[76]*nu[92]+0.7499999999999999*rdxF[1]*f[69]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 1.14564392373896*rdxF[2]*f[65]*nu[274]+1.677050983124842*rdxF[2]*f[77]*nu[254]+0.7499999999999999*rdxF[2]*f[65]*nu[244]+18.44756081437327*rdxF[1]*f[61]*nu[112]+10.5*rdxF[1]*f[17]*nu[112]+10.06230589874905*rdxF[1]*f[77]*nu[92]+12.8086884574495*rdxF[1]*f[35]*nu[92]+7.685213074469698*rdxF[1]*f[61]*nu[82]+6.873863542433758*rdxF[1]*f[17]*nu[82]+5.7282196186948*rdxF[1]*f[35]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 1.14564392373896*rdxF[2]*f[66]*nu[274]+1.677050983124842*rdxF[2]*f[78]*nu[254]+0.7499999999999999*rdxF[2]*f[66]*nu[244]+1.14564392373896*rdxF[1]*f[71]*nu[112]+1.677050983124842*rdxF[1]*f[78]*nu[92]+0.7499999999999999*rdxF[1]*f[71]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 18.44756081437327*rdxF[2]*f[63]*nu[274]+10.5*rdxF[2]*f[15]*nu[274]+10.06230589874905*rdxF[2]*f[79]*nu[254]+12.8086884574495*rdxF[2]*f[35]*nu[254]+7.685213074469698*rdxF[2]*f[63]*nu[244]+6.873863542433758*rdxF[2]*f[15]*nu[244]+5.7282196186948*rdxF[2]*f[35]*nu[240]+1.14564392373896*rdxF[1]*f[74]*nu[112]+1.677050983124842*rdxF[1]*f[79]*nu[92]+0.7499999999999999*rdxF[1]*f[74]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+rdxF[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs12(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[1]*dx[1]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[2] += 2.5617376914899*rdxF[1]*f[12]*nu[112]+1.14564392373896*f[0]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[2]*nu[92]+0.75*f[0]*rdxF[1]*nu[82]; 
  out[5] += 2.5617376914899*rdxF[1]*f[20]*nu[112]+1.14564392373896*f[1]*rdxF[1]*nu[112]+1.677050983124842*rdxF[1]*f[5]*nu[92]+0.75*f[1]*rdxF[1]*nu[82]+2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 2.5617376914899*rdxF[1]*f[22]*nu[112]+1.14564392373896*rdxF[1]*f[3]*nu[112]+1.677050983124842*rdxF[1]*f[7]*nu[92]+0.75*rdxF[1]*f[3]*nu[82]; 
  out[8] += 2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 2.5617376914899*rdxF[1]*f[26]*nu[112]+1.14564392373896*rdxF[1]*f[4]*nu[112]+1.677050983124842*rdxF[1]*f[9]*nu[92]+0.75*rdxF[1]*f[4]*nu[82]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[12] += 6.708203932499369*rdxF[1]*f[32]*nu[112]+7.685213074469699*rdxF[1]*f[2]*nu[112]+5.031152949374527*rdxF[1]*f[12]*nu[92]+3.75*f[0]*rdxF[1]*nu[92]+3.354101966249685*rdxF[1]*f[2]*nu[82]+1.677050983124842*f[0]*rdxF[1]*nu[80]; 
  out[15] += 2.5617376914899*rdxF[1]*f[37]*nu[112]+1.14564392373896*rdxF[1]*f[6]*nu[112]+1.677050983124842*rdxF[1]*f[15]*nu[92]+0.75*rdxF[1]*f[6]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[1]*f[40]*nu[112]+1.14564392373896*rdxF[1]*f[8]*nu[112]+1.677050983124842*rdxF[1]*f[16]*nu[92]+0.75*rdxF[1]*f[8]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 2.5617376914899*rdxF[1]*f[42]*nu[112]+1.14564392373896*rdxF[1]*f[10]*nu[112]+1.677050983124842*rdxF[1]*f[18]*nu[92]+0.75*rdxF[1]*f[10]*nu[82]; 
  out[19] += 1.14564392373896*rdxF[1]*f[11]*nu[112]+1.677050983124842*rdxF[1]*f[19]*nu[92]+0.75*rdxF[1]*f[11]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 6.708203932499369*rdxF[1]*f[49]*nu[112]+7.6852130744697*rdxF[1]*f[5]*nu[112]+5.031152949374527*rdxF[1]*f[20]*nu[92]+3.75*f[1]*rdxF[1]*nu[92]+3.354101966249684*rdxF[1]*f[5]*nu[82]+1.677050983124842*f[1]*rdxF[1]*nu[80]+1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 6.708203932499369*rdxF[1]*f[51]*nu[112]+7.6852130744697*rdxF[1]*f[7]*nu[112]+5.031152949374527*rdxF[1]*f[22]*nu[92]+3.75*rdxF[1]*f[3]*nu[92]+3.354101966249684*rdxF[1]*f[7]*nu[82]+1.677050983124842*rdxF[1]*f[3]*nu[80]; 
  out[23] += 1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 1.14564392373896*rdxF[1]*f[13]*nu[112]+1.677050983124842*rdxF[1]*f[24]*nu[92]+0.75*rdxF[1]*f[13]*nu[82]; 
  out[25] += 6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 6.708203932499369*rdxF[1]*f[55]*nu[112]+7.6852130744697*rdxF[1]*f[9]*nu[112]+5.031152949374527*rdxF[1]*f[26]*nu[92]+3.75*rdxF[1]*f[4]*nu[92]+3.354101966249684*rdxF[1]*f[9]*nu[82]+1.677050983124842*rdxF[1]*f[4]*nu[80]; 
  out[28] += 1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 1.14564392373896*rdxF[1]*f[14]*nu[112]+1.677050983124842*rdxF[1]*f[29]*nu[92]+0.75*rdxF[1]*f[14]*nu[82]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[32] += 18.44756081437327*rdxF[1]*f[12]*nu[112]+10.5*f[0]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[32]*nu[92]+12.8086884574495*rdxF[1]*f[2]*nu[92]+7.685213074469699*rdxF[1]*f[12]*nu[82]+6.873863542433759*f[0]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[2]*nu[80]; 
  out[35] += 2.5617376914899*rdxF[1]*f[61]*nu[112]+1.14564392373896*rdxF[1]*f[17]*nu[112]+1.677050983124842*rdxF[1]*f[35]*nu[92]+0.75*rdxF[1]*f[17]*nu[82]+2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 1.14564392373896*rdxF[1]*f[21]*nu[112]+1.677050983124842*rdxF[1]*f[36]*nu[92]+0.75*rdxF[1]*f[21]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 6.708203932499369*rdxF[1]*f[65]*nu[112]+7.685213074469699*rdxF[1]*f[15]*nu[112]+5.031152949374527*rdxF[1]*f[37]*nu[92]+3.75*rdxF[1]*f[6]*nu[92]+3.354101966249685*rdxF[1]*f[15]*nu[82]+1.677050983124842*rdxF[1]*f[6]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 1.14564392373896*rdxF[1]*f[23]*nu[112]+1.677050983124842*rdxF[1]*f[38]*nu[92]+0.75*rdxF[1]*f[23]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 1.14564392373896*rdxF[1]*f[25]*nu[112]+1.677050983124842*rdxF[1]*f[39]*nu[92]+0.75*rdxF[1]*f[25]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 6.708203932499369*rdxF[1]*f[68]*nu[112]+7.685213074469699*rdxF[1]*f[16]*nu[112]+5.031152949374527*rdxF[1]*f[40]*nu[92]+3.75*rdxF[1]*f[8]*nu[92]+3.354101966249685*rdxF[1]*f[16]*nu[82]+1.677050983124842*rdxF[1]*f[8]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[42] += 6.708203932499369*rdxF[1]*f[70]*nu[112]+7.685213074469699*rdxF[1]*f[18]*nu[112]+5.031152949374527*rdxF[1]*f[42]*nu[92]+3.75*rdxF[1]*f[10]*nu[92]+3.354101966249685*rdxF[1]*f[18]*nu[82]+1.677050983124842*rdxF[1]*f[10]*nu[80]; 
  out[43] += 1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[44] += 1.14564392373896*rdxF[1]*f[27]*nu[112]+1.677050983124842*rdxF[1]*f[44]*nu[92]+0.75*rdxF[1]*f[27]*nu[82]; 
  out[45] += 1.14564392373896*rdxF[1]*f[28]*nu[112]+1.677050983124842*rdxF[1]*f[45]*nu[92]+0.75*rdxF[1]*f[28]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[47] += 1.14564392373896*rdxF[1]*f[30]*nu[112]+1.677050983124842*rdxF[1]*f[47]*nu[92]+0.75*rdxF[1]*f[30]*nu[82]; 
  out[48] += 1.14564392373896*rdxF[1]*f[31]*nu[112]+1.677050983124842*rdxF[1]*f[48]*nu[92]+0.7499999999999999*rdxF[1]*f[31]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 18.44756081437327*rdxF[1]*f[20]*nu[112]+10.5*f[1]*rdxF[1]*nu[112]+10.06230589874905*rdxF[1]*f[49]*nu[92]+12.8086884574495*rdxF[1]*f[5]*nu[92]+7.685213074469698*rdxF[1]*f[20]*nu[82]+6.873863542433758*f[1]*rdxF[1]*nu[82]+5.7282196186948*rdxF[1]*f[5]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[51] += 18.44756081437327*rdxF[1]*f[22]*nu[112]+10.5*rdxF[1]*f[3]*nu[112]+10.06230589874905*rdxF[1]*f[51]*nu[92]+12.8086884574495*rdxF[1]*f[7]*nu[92]+7.685213074469698*rdxF[1]*f[22]*nu[82]+6.873863542433758*rdxF[1]*f[3]*nu[82]+5.7282196186948*rdxF[1]*f[7]*nu[80]; 
  out[52] += 1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[53] += 1.14564392373896*rdxF[1]*f[33]*nu[112]+1.677050983124842*rdxF[1]*f[53]*nu[92]+0.7499999999999999*rdxF[1]*f[33]*nu[82]; 
  out[54] += 10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[55] += 18.44756081437327*rdxF[1]*f[26]*nu[112]+10.5*rdxF[1]*f[4]*nu[112]+10.06230589874905*rdxF[1]*f[55]*nu[92]+12.8086884574495*rdxF[1]*f[9]*nu[92]+7.685213074469698*rdxF[1]*f[26]*nu[82]+6.873863542433758*rdxF[1]*f[4]*nu[82]+5.7282196186948*rdxF[1]*f[9]*nu[80]; 
  out[57] += 1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[58] += 1.14564392373896*rdxF[1]*f[34]*nu[112]+1.677050983124842*rdxF[1]*f[58]*nu[92]+0.7499999999999999*rdxF[1]*f[34]*nu[82]; 
  out[60] += 1.14564392373896*rdxF[1]*f[41]*nu[112]+1.677050983124842*rdxF[1]*f[60]*nu[92]+0.75*rdxF[1]*f[41]*nu[82]+6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 6.708203932499369*rdxF[1]*f[77]*nu[112]+7.6852130744697*rdxF[1]*f[35]*nu[112]+5.031152949374527*rdxF[1]*f[61]*nu[92]+3.75*rdxF[1]*f[17]*nu[92]+3.354101966249684*rdxF[1]*f[35]*nu[82]+1.677050983124842*rdxF[1]*f[17]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 1.14564392373896*rdxF[1]*f[43]*nu[112]+1.677050983124842*rdxF[1]*f[62]*nu[92]+0.75*rdxF[1]*f[43]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 1.14564392373896*rdxF[1]*f[46]*nu[112]+1.677050983124842*rdxF[1]*f[63]*nu[92]+0.75*rdxF[1]*f[46]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 1.14564392373896*rdxF[1]*f[50]*nu[112]+1.677050983124842*rdxF[1]*f[64]*nu[92]+0.7499999999999999*rdxF[1]*f[50]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 18.44756081437327*rdxF[1]*f[37]*nu[112]+10.5*rdxF[1]*f[6]*nu[112]+10.06230589874905*rdxF[1]*f[65]*nu[92]+12.8086884574495*rdxF[1]*f[15]*nu[92]+7.685213074469699*rdxF[1]*f[37]*nu[82]+6.873863542433759*rdxF[1]*f[6]*nu[82]+5.7282196186948*rdxF[1]*f[15]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 1.14564392373896*rdxF[1]*f[52]*nu[112]+1.677050983124842*rdxF[1]*f[66]*nu[92]+0.7499999999999999*rdxF[1]*f[52]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 1.14564392373896*rdxF[1]*f[54]*nu[112]+1.677050983124842*rdxF[1]*f[67]*nu[92]+0.7499999999999999*rdxF[1]*f[54]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 18.44756081437327*rdxF[1]*f[40]*nu[112]+10.5*rdxF[1]*f[8]*nu[112]+10.06230589874905*rdxF[1]*f[68]*nu[92]+12.8086884574495*rdxF[1]*f[16]*nu[92]+7.685213074469699*rdxF[1]*f[40]*nu[82]+6.873863542433759*rdxF[1]*f[8]*nu[82]+5.7282196186948*rdxF[1]*f[16]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[70] += 18.44756081437327*rdxF[1]*f[42]*nu[112]+10.5*rdxF[1]*f[10]*nu[112]+10.06230589874905*rdxF[1]*f[70]*nu[92]+12.8086884574495*rdxF[1]*f[18]*nu[92]+7.685213074469699*rdxF[1]*f[42]*nu[82]+6.873863542433759*rdxF[1]*f[10]*nu[82]+5.7282196186948*rdxF[1]*f[18]*nu[80]; 
  out[71] += 1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[72] += 1.14564392373896*rdxF[1]*f[56]*nu[112]+1.677050983124842*rdxF[1]*f[72]*nu[92]+0.7499999999999999*rdxF[1]*f[56]*nu[82]; 
  out[73] += 1.14564392373896*rdxF[1]*f[57]*nu[112]+1.677050983124842*rdxF[1]*f[73]*nu[92]+0.7499999999999999*rdxF[1]*f[57]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[75] += 1.14564392373896*rdxF[1]*f[59]*nu[112]+1.677050983124842*rdxF[1]*f[75]*nu[92]+0.7499999999999999*rdxF[1]*f[59]*nu[82]; 
  out[76] += 1.14564392373896*rdxF[1]*f[69]*nu[112]+1.677050983124842*rdxF[1]*f[76]*nu[92]+0.7499999999999999*rdxF[1]*f[69]*nu[82]+10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 18.44756081437327*rdxF[1]*f[61]*nu[112]+10.5*rdxF[1]*f[17]*nu[112]+10.06230589874905*rdxF[1]*f[77]*nu[92]+12.8086884574495*rdxF[1]*f[35]*nu[92]+7.685213074469698*rdxF[1]*f[61]*nu[82]+6.873863542433758*rdxF[1]*f[17]*nu[82]+5.7282196186948*rdxF[1]*f[35]*nu[80]+1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 1.14564392373896*rdxF[1]*f[71]*nu[112]+1.677050983124842*rdxF[1]*f[78]*nu[92]+0.7499999999999999*rdxF[1]*f[71]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 1.14564392373896*rdxF[1]*f[74]*nu[112]+1.677050983124842*rdxF[1]*f[79]*nu[92]+0.7499999999999999*rdxF[1]*f[74]*nu[82]+1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs3(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 

  out[3] += 2.5617376914899*rdxF[0]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[0]*nu[193]+1.677050983124842*rdxF[0]*f[3]*nu[173]+0.75*f[0]*rdxF[0]*nu[163]; 
  out[6] += 2.5617376914899*rdxF[0]*f[23]*nu[193]+1.14564392373896*rdxF[0]*f[1]*nu[193]+1.677050983124842*rdxF[0]*f[6]*nu[173]+0.75*rdxF[0]*f[1]*nu[163]; 
  out[7] += 2.5617376914899*rdxF[0]*f[24]*nu[193]+1.14564392373896*rdxF[0]*f[2]*nu[193]+1.677050983124842*rdxF[0]*f[7]*nu[173]+0.75*rdxF[0]*f[2]*nu[163]; 
  out[10] += 2.5617376914899*rdxF[0]*f[27]*nu[193]+1.14564392373896*rdxF[0]*f[4]*nu[193]+1.677050983124842*rdxF[0]*f[10]*nu[173]+0.75*rdxF[0]*f[4]*nu[163]; 
  out[13] += 6.708203932499369*rdxF[0]*f[33]*nu[193]+7.685213074469699*rdxF[0]*f[3]*nu[193]+5.031152949374527*rdxF[0]*f[13]*nu[173]+3.75*f[0]*rdxF[0]*nu[173]+3.354101966249685*rdxF[0]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[0]*nu[160]; 
  out[15] += 2.5617376914899*rdxF[0]*f[38]*nu[193]+1.14564392373896*rdxF[0]*f[5]*nu[193]+1.677050983124842*rdxF[0]*f[15]*nu[173]+0.75*rdxF[0]*f[5]*nu[163]; 
  out[17] += 2.5617376914899*rdxF[0]*f[43]*nu[193]+1.14564392373896*rdxF[0]*f[8]*nu[193]+1.677050983124842*rdxF[0]*f[17]*nu[173]+0.75*rdxF[0]*f[8]*nu[163]; 
  out[18] += 2.5617376914899*rdxF[0]*f[44]*nu[193]+1.14564392373896*rdxF[0]*f[9]*nu[193]+1.677050983124842*rdxF[0]*f[18]*nu[173]+0.75*rdxF[0]*f[9]*nu[163]; 
  out[21] += 1.14564392373896*rdxF[0]*f[11]*nu[193]+1.677050983124842*rdxF[0]*f[21]*nu[173]+0.75*rdxF[0]*f[11]*nu[163]; 
  out[22] += 1.14564392373896*rdxF[0]*f[12]*nu[193]+1.677050983124842*rdxF[0]*f[22]*nu[173]+0.75*rdxF[0]*f[12]*nu[163]; 
  out[23] += 6.708203932499369*rdxF[0]*f[52]*nu[193]+7.6852130744697*rdxF[0]*f[6]*nu[193]+5.031152949374527*rdxF[0]*f[23]*nu[173]+3.75*rdxF[0]*f[1]*nu[173]+3.354101966249684*rdxF[0]*f[6]*nu[163]+1.677050983124842*rdxF[0]*f[1]*nu[160]; 
  out[24] += 6.708203932499369*rdxF[0]*f[53]*nu[193]+7.6852130744697*rdxF[0]*f[7]*nu[193]+5.031152949374527*rdxF[0]*f[24]*nu[173]+3.75*rdxF[0]*f[2]*nu[173]+3.354101966249684*rdxF[0]*f[7]*nu[163]+1.677050983124842*rdxF[0]*f[2]*nu[160]; 
  out[27] += 6.708203932499369*rdxF[0]*f[56]*nu[193]+7.6852130744697*rdxF[0]*f[10]*nu[193]+5.031152949374527*rdxF[0]*f[27]*nu[173]+3.75*rdxF[0]*f[4]*nu[173]+3.354101966249684*rdxF[0]*f[10]*nu[163]+1.677050983124842*rdxF[0]*f[4]*nu[160]; 
  out[30] += 1.14564392373896*rdxF[0]*f[14]*nu[193]+1.677050983124842*rdxF[0]*f[30]*nu[173]+0.75*rdxF[0]*f[14]*nu[163]; 
  out[33] += 18.44756081437327*rdxF[0]*f[13]*nu[193]+10.5*f[0]*rdxF[0]*nu[193]+10.06230589874905*rdxF[0]*f[33]*nu[173]+12.8086884574495*rdxF[0]*f[3]*nu[173]+7.685213074469699*rdxF[0]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[0]*nu[163]+5.7282196186948*rdxF[0]*f[3]*nu[160]; 
  out[35] += 2.5617376914899*rdxF[0]*f[62]*nu[193]+1.14564392373896*rdxF[0]*f[16]*nu[193]+1.677050983124842*rdxF[0]*f[35]*nu[173]+0.75*rdxF[0]*f[16]*nu[163]; 
  out[36] += 1.14564392373896*rdxF[0]*f[19]*nu[193]+1.677050983124842*rdxF[0]*f[36]*nu[173]+0.75*rdxF[0]*f[19]*nu[163]; 
  out[37] += 1.14564392373896*rdxF[0]*f[20]*nu[193]+1.677050983124842*rdxF[0]*f[37]*nu[173]+0.75*rdxF[0]*f[20]*nu[163]; 
  out[38] += 6.708203932499369*rdxF[0]*f[66]*nu[193]+7.685213074469699*rdxF[0]*f[15]*nu[193]+5.031152949374527*rdxF[0]*f[38]*nu[173]+3.75*rdxF[0]*f[5]*nu[173]+3.354101966249685*rdxF[0]*f[15]*nu[163]+1.677050983124842*rdxF[0]*f[5]*nu[160]; 
  out[41] += 1.14564392373896*rdxF[0]*f[25]*nu[193]+1.677050983124842*rdxF[0]*f[41]*nu[173]+0.75*rdxF[0]*f[25]*nu[163]; 
  out[42] += 1.14564392373896*rdxF[0]*f[26]*nu[193]+1.677050983124842*rdxF[0]*f[42]*nu[173]+0.75*rdxF[0]*f[26]*nu[163]; 
  out[43] += 6.708203932499369*rdxF[0]*f[71]*nu[193]+7.685213074469699*rdxF[0]*f[17]*nu[193]+5.031152949374527*rdxF[0]*f[43]*nu[173]+3.75*rdxF[0]*f[8]*nu[173]+3.354101966249685*rdxF[0]*f[17]*nu[163]+1.677050983124842*rdxF[0]*f[8]*nu[160]; 
  out[44] += 6.708203932499369*rdxF[0]*f[72]*nu[193]+7.685213074469699*rdxF[0]*f[18]*nu[193]+5.031152949374527*rdxF[0]*f[44]*nu[173]+3.75*rdxF[0]*f[9]*nu[173]+3.354101966249685*rdxF[0]*f[18]*nu[163]+1.677050983124842*rdxF[0]*f[9]*nu[160]; 
  out[46] += 1.14564392373896*rdxF[0]*f[28]*nu[193]+1.677050983124842*rdxF[0]*f[46]*nu[173]+0.75*rdxF[0]*f[28]*nu[163]; 
  out[47] += 1.14564392373896*rdxF[0]*f[29]*nu[193]+1.677050983124842*rdxF[0]*f[47]*nu[173]+0.75*rdxF[0]*f[29]*nu[163]; 
  out[50] += 1.14564392373896*rdxF[0]*f[31]*nu[193]+1.677050983124842*rdxF[0]*f[50]*nu[173]+0.7499999999999999*rdxF[0]*f[31]*nu[163]; 
  out[51] += 1.14564392373896*rdxF[0]*f[32]*nu[193]+1.677050983124842*rdxF[0]*f[51]*nu[173]+0.7499999999999999*rdxF[0]*f[32]*nu[163]; 
  out[52] += 18.44756081437327*rdxF[0]*f[23]*nu[193]+10.5*rdxF[0]*f[1]*nu[193]+10.06230589874905*rdxF[0]*f[52]*nu[173]+12.8086884574495*rdxF[0]*f[6]*nu[173]+7.685213074469698*rdxF[0]*f[23]*nu[163]+6.873863542433758*rdxF[0]*f[1]*nu[163]+5.7282196186948*rdxF[0]*f[6]*nu[160]; 
  out[53] += 18.44756081437327*rdxF[0]*f[24]*nu[193]+10.5*rdxF[0]*f[2]*nu[193]+10.06230589874905*rdxF[0]*f[53]*nu[173]+12.8086884574495*rdxF[0]*f[7]*nu[173]+7.685213074469698*rdxF[0]*f[24]*nu[163]+6.873863542433758*rdxF[0]*f[2]*nu[163]+5.7282196186948*rdxF[0]*f[7]*nu[160]; 
  out[56] += 18.44756081437327*rdxF[0]*f[27]*nu[193]+10.5*rdxF[0]*f[4]*nu[193]+10.06230589874905*rdxF[0]*f[56]*nu[173]+12.8086884574495*rdxF[0]*f[10]*nu[173]+7.685213074469698*rdxF[0]*f[27]*nu[163]+6.873863542433758*rdxF[0]*f[4]*nu[163]+5.7282196186948*rdxF[0]*f[10]*nu[160]; 
  out[59] += 1.14564392373896*rdxF[0]*f[34]*nu[193]+1.677050983124842*rdxF[0]*f[59]*nu[173]+0.7499999999999999*rdxF[0]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[0]*f[39]*nu[193]+1.677050983124842*rdxF[0]*f[60]*nu[173]+0.75*rdxF[0]*f[39]*nu[163]; 
  out[61] += 1.14564392373896*rdxF[0]*f[40]*nu[193]+1.677050983124842*rdxF[0]*f[61]*nu[173]+0.75*rdxF[0]*f[40]*nu[163]; 
  out[62] += 6.708203932499369*rdxF[0]*f[78]*nu[193]+7.6852130744697*rdxF[0]*f[35]*nu[193]+5.031152949374527*rdxF[0]*f[62]*nu[173]+3.75*rdxF[0]*f[16]*nu[173]+3.354101966249684*rdxF[0]*f[35]*nu[163]+1.677050983124842*rdxF[0]*f[16]*nu[160]; 
  out[63] += 1.14564392373896*rdxF[0]*f[45]*nu[193]+1.677050983124842*rdxF[0]*f[63]*nu[173]+0.75*rdxF[0]*f[45]*nu[163]; 
  out[64] += 1.14564392373896*rdxF[0]*f[48]*nu[193]+1.677050983124842*rdxF[0]*f[64]*nu[173]+0.7499999999999999*rdxF[0]*f[48]*nu[163]; 
  out[65] += 1.14564392373896*rdxF[0]*f[49]*nu[193]+1.677050983124842*rdxF[0]*f[65]*nu[173]+0.7499999999999999*rdxF[0]*f[49]*nu[163]; 
  out[66] += 18.44756081437327*rdxF[0]*f[38]*nu[193]+10.5*rdxF[0]*f[5]*nu[193]+10.06230589874905*rdxF[0]*f[66]*nu[173]+12.8086884574495*rdxF[0]*f[15]*nu[173]+7.685213074469699*rdxF[0]*f[38]*nu[163]+6.873863542433759*rdxF[0]*f[5]*nu[163]+5.7282196186948*rdxF[0]*f[15]*nu[160]; 
  out[69] += 1.14564392373896*rdxF[0]*f[54]*nu[193]+1.677050983124842*rdxF[0]*f[69]*nu[173]+0.7499999999999999*rdxF[0]*f[54]*nu[163]; 
  out[70] += 1.14564392373896*rdxF[0]*f[55]*nu[193]+1.677050983124842*rdxF[0]*f[70]*nu[173]+0.7499999999999999*rdxF[0]*f[55]*nu[163]; 
  out[71] += 18.44756081437327*rdxF[0]*f[43]*nu[193]+10.5*rdxF[0]*f[8]*nu[193]+10.06230589874905*rdxF[0]*f[71]*nu[173]+12.8086884574495*rdxF[0]*f[17]*nu[173]+7.685213074469699*rdxF[0]*f[43]*nu[163]+6.873863542433759*rdxF[0]*f[8]*nu[163]+5.7282196186948*rdxF[0]*f[17]*nu[160]; 
  out[72] += 18.44756081437327*rdxF[0]*f[44]*nu[193]+10.5*rdxF[0]*f[9]*nu[193]+10.06230589874905*rdxF[0]*f[72]*nu[173]+12.8086884574495*rdxF[0]*f[18]*nu[173]+7.685213074469699*rdxF[0]*f[44]*nu[163]+6.873863542433759*rdxF[0]*f[9]*nu[163]+5.7282196186948*rdxF[0]*f[18]*nu[160]; 
  out[74] += 1.14564392373896*rdxF[0]*f[57]*nu[193]+1.677050983124842*rdxF[0]*f[74]*nu[173]+0.7499999999999999*rdxF[0]*f[57]*nu[163]; 
  out[75] += 1.14564392373896*rdxF[0]*f[58]*nu[193]+1.677050983124842*rdxF[0]*f[75]*nu[173]+0.7499999999999999*rdxF[0]*f[58]*nu[163]; 
  out[76] += 1.14564392373896*rdxF[0]*f[67]*nu[193]+1.677050983124842*rdxF[0]*f[76]*nu[173]+0.7499999999999999*rdxF[0]*f[67]*nu[163]; 
  out[77] += 1.14564392373896*rdxF[0]*f[68]*nu[193]+1.677050983124842*rdxF[0]*f[77]*nu[173]+0.7499999999999999*rdxF[0]*f[68]*nu[163]; 
  out[78] += 18.44756081437327*rdxF[0]*f[62]*nu[193]+10.5*rdxF[0]*f[16]*nu[193]+10.06230589874905*rdxF[0]*f[78]*nu[173]+12.8086884574495*rdxF[0]*f[35]*nu[173]+7.685213074469698*rdxF[0]*f[62]*nu[163]+6.873863542433758*rdxF[0]*f[16]*nu[163]+5.7282196186948*rdxF[0]*f[35]*nu[160]; 
  out[79] += 1.14564392373896*rdxF[0]*f[73]*nu[193]+1.677050983124842*rdxF[0]*f[79]*nu[173]+0.7499999999999999*rdxF[0]*f[73]*nu[163]; 

  return (rdxF[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs34(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[2]*dx[2]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[3] += 2.5617376914899*rdxF[0]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[0]*nu[193]+1.677050983124842*rdxF[0]*f[3]*nu[173]+0.75*f[0]*rdxF[0]*nu[163]; 
  out[4] += 2.5617376914899*rdxF[1]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[1]*nu[274]+1.677050983124842*rdxF[1]*f[4]*nu[254]+0.75*f[0]*rdxF[1]*nu[244]; 
  out[6] += 2.5617376914899*rdxF[0]*f[23]*nu[193]+1.14564392373896*rdxF[0]*f[1]*nu[193]+1.677050983124842*rdxF[0]*f[6]*nu[173]+0.75*rdxF[0]*f[1]*nu[163]; 
  out[7] += 2.5617376914899*rdxF[0]*f[24]*nu[193]+1.14564392373896*rdxF[0]*f[2]*nu[193]+1.677050983124842*rdxF[0]*f[7]*nu[173]+0.75*rdxF[0]*f[2]*nu[163]; 
  out[8] += 2.5617376914899*rdxF[1]*f[28]*nu[274]+1.14564392373896*f[1]*rdxF[1]*nu[274]+1.677050983124842*rdxF[1]*f[8]*nu[254]+0.75*f[1]*rdxF[1]*nu[244]; 
  out[9] += 2.5617376914899*rdxF[1]*f[29]*nu[274]+1.14564392373896*rdxF[1]*f[2]*nu[274]+1.677050983124842*rdxF[1]*f[9]*nu[254]+0.75*rdxF[1]*f[2]*nu[244]; 
  out[10] += 2.5617376914899*rdxF[1]*f[30]*nu[274]+1.14564392373896*rdxF[1]*f[3]*nu[274]+1.677050983124842*rdxF[1]*f[10]*nu[254]+0.75*rdxF[1]*f[3]*nu[244]+2.5617376914899*rdxF[0]*f[27]*nu[193]+1.14564392373896*rdxF[0]*f[4]*nu[193]+1.677050983124842*rdxF[0]*f[10]*nu[173]+0.75*rdxF[0]*f[4]*nu[163]; 
  out[13] += 6.708203932499369*rdxF[0]*f[33]*nu[193]+7.685213074469699*rdxF[0]*f[3]*nu[193]+5.031152949374527*rdxF[0]*f[13]*nu[173]+3.75*f[0]*rdxF[0]*nu[173]+3.354101966249685*rdxF[0]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[0]*nu[160]; 
  out[14] += 6.708203932499369*rdxF[1]*f[34]*nu[274]+7.685213074469699*rdxF[1]*f[4]*nu[274]+5.031152949374527*rdxF[1]*f[14]*nu[254]+3.75*f[0]*rdxF[1]*nu[254]+3.354101966249685*rdxF[1]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[1]*nu[240]; 
  out[15] += 2.5617376914899*rdxF[0]*f[38]*nu[193]+1.14564392373896*rdxF[0]*f[5]*nu[193]+1.677050983124842*rdxF[0]*f[15]*nu[173]+0.75*rdxF[0]*f[5]*nu[163]; 
  out[16] += 2.5617376914899*rdxF[1]*f[45]*nu[274]+1.14564392373896*rdxF[1]*f[5]*nu[274]+1.677050983124842*rdxF[1]*f[16]*nu[254]+0.75*rdxF[1]*f[5]*nu[244]; 
  out[17] += 2.5617376914899*rdxF[1]*f[46]*nu[274]+1.14564392373896*rdxF[1]*f[6]*nu[274]+1.677050983124842*rdxF[1]*f[17]*nu[254]+0.75*rdxF[1]*f[6]*nu[244]+2.5617376914899*rdxF[0]*f[43]*nu[193]+1.14564392373896*rdxF[0]*f[8]*nu[193]+1.677050983124842*rdxF[0]*f[17]*nu[173]+0.75*rdxF[0]*f[8]*nu[163]; 
  out[18] += 2.5617376914899*rdxF[1]*f[47]*nu[274]+1.14564392373896*rdxF[1]*f[7]*nu[274]+1.677050983124842*rdxF[1]*f[18]*nu[254]+0.75*rdxF[1]*f[7]*nu[244]+2.5617376914899*rdxF[0]*f[44]*nu[193]+1.14564392373896*rdxF[0]*f[9]*nu[193]+1.677050983124842*rdxF[0]*f[18]*nu[173]+0.75*rdxF[0]*f[9]*nu[163]; 
  out[21] += 1.14564392373896*rdxF[0]*f[11]*nu[193]+1.677050983124842*rdxF[0]*f[21]*nu[173]+0.75*rdxF[0]*f[11]*nu[163]; 
  out[22] += 1.14564392373896*rdxF[0]*f[12]*nu[193]+1.677050983124842*rdxF[0]*f[22]*nu[173]+0.75*rdxF[0]*f[12]*nu[163]; 
  out[23] += 6.708203932499369*rdxF[0]*f[52]*nu[193]+7.6852130744697*rdxF[0]*f[6]*nu[193]+5.031152949374527*rdxF[0]*f[23]*nu[173]+3.75*rdxF[0]*f[1]*nu[173]+3.354101966249684*rdxF[0]*f[6]*nu[163]+1.677050983124842*rdxF[0]*f[1]*nu[160]; 
  out[24] += 6.708203932499369*rdxF[0]*f[53]*nu[193]+7.6852130744697*rdxF[0]*f[7]*nu[193]+5.031152949374527*rdxF[0]*f[24]*nu[173]+3.75*rdxF[0]*f[2]*nu[173]+3.354101966249684*rdxF[0]*f[7]*nu[163]+1.677050983124842*rdxF[0]*f[2]*nu[160]; 
  out[25] += 1.14564392373896*rdxF[1]*f[11]*nu[274]+1.677050983124842*rdxF[1]*f[25]*nu[254]+0.75*rdxF[1]*f[11]*nu[244]; 
  out[26] += 1.14564392373896*rdxF[1]*f[12]*nu[274]+1.677050983124842*rdxF[1]*f[26]*nu[254]+0.75*rdxF[1]*f[12]*nu[244]; 
  out[27] += 1.14564392373896*rdxF[1]*f[13]*nu[274]+1.677050983124842*rdxF[1]*f[27]*nu[254]+0.75*rdxF[1]*f[13]*nu[244]+6.708203932499369*rdxF[0]*f[56]*nu[193]+7.6852130744697*rdxF[0]*f[10]*nu[193]+5.031152949374527*rdxF[0]*f[27]*nu[173]+3.75*rdxF[0]*f[4]*nu[173]+3.354101966249684*rdxF[0]*f[10]*nu[163]+1.677050983124842*rdxF[0]*f[4]*nu[160]; 
  out[28] += 6.708203932499369*rdxF[1]*f[57]*nu[274]+7.6852130744697*rdxF[1]*f[8]*nu[274]+5.031152949374527*rdxF[1]*f[28]*nu[254]+3.75*f[1]*rdxF[1]*nu[254]+3.354101966249684*rdxF[1]*f[8]*nu[244]+1.677050983124842*f[1]*rdxF[1]*nu[240]; 
  out[29] += 6.708203932499369*rdxF[1]*f[58]*nu[274]+7.6852130744697*rdxF[1]*f[9]*nu[274]+5.031152949374527*rdxF[1]*f[29]*nu[254]+3.75*rdxF[1]*f[2]*nu[254]+3.354101966249684*rdxF[1]*f[9]*nu[244]+1.677050983124842*rdxF[1]*f[2]*nu[240]; 
  out[30] += 6.708203932499369*rdxF[1]*f[59]*nu[274]+7.6852130744697*rdxF[1]*f[10]*nu[274]+5.031152949374527*rdxF[1]*f[30]*nu[254]+3.75*rdxF[1]*f[3]*nu[254]+3.354101966249684*rdxF[1]*f[10]*nu[244]+1.677050983124842*rdxF[1]*f[3]*nu[240]+1.14564392373896*rdxF[0]*f[14]*nu[193]+1.677050983124842*rdxF[0]*f[30]*nu[173]+0.75*rdxF[0]*f[14]*nu[163]; 
  out[33] += 18.44756081437327*rdxF[0]*f[13]*nu[193]+10.5*f[0]*rdxF[0]*nu[193]+10.06230589874905*rdxF[0]*f[33]*nu[173]+12.8086884574495*rdxF[0]*f[3]*nu[173]+7.685213074469699*rdxF[0]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[0]*nu[163]+5.7282196186948*rdxF[0]*f[3]*nu[160]; 
  out[34] += 18.44756081437327*rdxF[1]*f[14]*nu[274]+10.5*f[0]*rdxF[1]*nu[274]+10.06230589874905*rdxF[1]*f[34]*nu[254]+12.8086884574495*rdxF[1]*f[4]*nu[254]+7.685213074469699*rdxF[1]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[1]*nu[244]+5.7282196186948*rdxF[1]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[1]*f[63]*nu[274]+1.14564392373896*rdxF[1]*f[15]*nu[274]+1.677050983124842*rdxF[1]*f[35]*nu[254]+0.75*rdxF[1]*f[15]*nu[244]+2.5617376914899*rdxF[0]*f[62]*nu[193]+1.14564392373896*rdxF[0]*f[16]*nu[193]+1.677050983124842*rdxF[0]*f[35]*nu[173]+0.75*rdxF[0]*f[16]*nu[163]; 
  out[36] += 1.14564392373896*rdxF[0]*f[19]*nu[193]+1.677050983124842*rdxF[0]*f[36]*nu[173]+0.75*rdxF[0]*f[19]*nu[163]; 
  out[37] += 1.14564392373896*rdxF[0]*f[20]*nu[193]+1.677050983124842*rdxF[0]*f[37]*nu[173]+0.75*rdxF[0]*f[20]*nu[163]; 
  out[38] += 6.708203932499369*rdxF[0]*f[66]*nu[193]+7.685213074469699*rdxF[0]*f[15]*nu[193]+5.031152949374527*rdxF[0]*f[38]*nu[173]+3.75*rdxF[0]*f[5]*nu[173]+3.354101966249685*rdxF[0]*f[15]*nu[163]+1.677050983124842*rdxF[0]*f[5]*nu[160]; 
  out[39] += 1.14564392373896*rdxF[1]*f[19]*nu[274]+1.677050983124842*rdxF[1]*f[39]*nu[254]+0.75*rdxF[1]*f[19]*nu[244]; 
  out[40] += 1.14564392373896*rdxF[1]*f[20]*nu[274]+1.677050983124842*rdxF[1]*f[40]*nu[254]+0.75*rdxF[1]*f[20]*nu[244]; 
  out[41] += 1.14564392373896*rdxF[1]*f[21]*nu[274]+1.677050983124842*rdxF[1]*f[41]*nu[254]+0.75*rdxF[1]*f[21]*nu[244]+1.14564392373896*rdxF[0]*f[25]*nu[193]+1.677050983124842*rdxF[0]*f[41]*nu[173]+0.75*rdxF[0]*f[25]*nu[163]; 
  out[42] += 1.14564392373896*rdxF[1]*f[22]*nu[274]+1.677050983124842*rdxF[1]*f[42]*nu[254]+0.75*rdxF[1]*f[22]*nu[244]+1.14564392373896*rdxF[0]*f[26]*nu[193]+1.677050983124842*rdxF[0]*f[42]*nu[173]+0.75*rdxF[0]*f[26]*nu[163]; 
  out[43] += 1.14564392373896*rdxF[1]*f[23]*nu[274]+1.677050983124842*rdxF[1]*f[43]*nu[254]+0.75*rdxF[1]*f[23]*nu[244]+6.708203932499369*rdxF[0]*f[71]*nu[193]+7.685213074469699*rdxF[0]*f[17]*nu[193]+5.031152949374527*rdxF[0]*f[43]*nu[173]+3.75*rdxF[0]*f[8]*nu[173]+3.354101966249685*rdxF[0]*f[17]*nu[163]+1.677050983124842*rdxF[0]*f[8]*nu[160]; 
  out[44] += 1.14564392373896*rdxF[1]*f[24]*nu[274]+1.677050983124842*rdxF[1]*f[44]*nu[254]+0.75*rdxF[1]*f[24]*nu[244]+6.708203932499369*rdxF[0]*f[72]*nu[193]+7.685213074469699*rdxF[0]*f[18]*nu[193]+5.031152949374527*rdxF[0]*f[44]*nu[173]+3.75*rdxF[0]*f[9]*nu[173]+3.354101966249685*rdxF[0]*f[18]*nu[163]+1.677050983124842*rdxF[0]*f[9]*nu[160]; 
  out[45] += 6.708203932499369*rdxF[1]*f[73]*nu[274]+7.685213074469699*rdxF[1]*f[16]*nu[274]+5.031152949374527*rdxF[1]*f[45]*nu[254]+3.75*rdxF[1]*f[5]*nu[254]+3.354101966249685*rdxF[1]*f[16]*nu[244]+1.677050983124842*rdxF[1]*f[5]*nu[240]; 
  out[46] += 6.708203932499369*rdxF[1]*f[74]*nu[274]+7.685213074469699*rdxF[1]*f[17]*nu[274]+5.031152949374527*rdxF[1]*f[46]*nu[254]+3.75*rdxF[1]*f[6]*nu[254]+3.354101966249685*rdxF[1]*f[17]*nu[244]+1.677050983124842*rdxF[1]*f[6]*nu[240]+1.14564392373896*rdxF[0]*f[28]*nu[193]+1.677050983124842*rdxF[0]*f[46]*nu[173]+0.75*rdxF[0]*f[28]*nu[163]; 
  out[47] += 6.708203932499369*rdxF[1]*f[75]*nu[274]+7.685213074469699*rdxF[1]*f[18]*nu[274]+5.031152949374527*rdxF[1]*f[47]*nu[254]+3.75*rdxF[1]*f[7]*nu[254]+3.354101966249685*rdxF[1]*f[18]*nu[244]+1.677050983124842*rdxF[1]*f[7]*nu[240]+1.14564392373896*rdxF[0]*f[29]*nu[193]+1.677050983124842*rdxF[0]*f[47]*nu[173]+0.75*rdxF[0]*f[29]*nu[163]; 
  out[50] += 1.14564392373896*rdxF[0]*f[31]*nu[193]+1.677050983124842*rdxF[0]*f[50]*nu[173]+0.7499999999999999*rdxF[0]*f[31]*nu[163]; 
  out[51] += 1.14564392373896*rdxF[0]*f[32]*nu[193]+1.677050983124842*rdxF[0]*f[51]*nu[173]+0.7499999999999999*rdxF[0]*f[32]*nu[163]; 
  out[52] += 18.44756081437327*rdxF[0]*f[23]*nu[193]+10.5*rdxF[0]*f[1]*nu[193]+10.06230589874905*rdxF[0]*f[52]*nu[173]+12.8086884574495*rdxF[0]*f[6]*nu[173]+7.685213074469698*rdxF[0]*f[23]*nu[163]+6.873863542433758*rdxF[0]*f[1]*nu[163]+5.7282196186948*rdxF[0]*f[6]*nu[160]; 
  out[53] += 18.44756081437327*rdxF[0]*f[24]*nu[193]+10.5*rdxF[0]*f[2]*nu[193]+10.06230589874905*rdxF[0]*f[53]*nu[173]+12.8086884574495*rdxF[0]*f[7]*nu[173]+7.685213074469698*rdxF[0]*f[24]*nu[163]+6.873863542433758*rdxF[0]*f[2]*nu[163]+5.7282196186948*rdxF[0]*f[7]*nu[160]; 
  out[54] += 1.14564392373896*rdxF[1]*f[31]*nu[274]+1.677050983124842*rdxF[1]*f[54]*nu[254]+0.7499999999999999*rdxF[1]*f[31]*nu[244]; 
  out[55] += 1.14564392373896*rdxF[1]*f[32]*nu[274]+1.677050983124842*rdxF[1]*f[55]*nu[254]+0.7499999999999999*rdxF[1]*f[32]*nu[244]; 
  out[56] += 1.14564392373896*rdxF[1]*f[33]*nu[274]+1.677050983124842*rdxF[1]*f[56]*nu[254]+0.7499999999999999*rdxF[1]*f[33]*nu[244]+18.44756081437327*rdxF[0]*f[27]*nu[193]+10.5*rdxF[0]*f[4]*nu[193]+10.06230589874905*rdxF[0]*f[56]*nu[173]+12.8086884574495*rdxF[0]*f[10]*nu[173]+7.685213074469698*rdxF[0]*f[27]*nu[163]+6.873863542433758*rdxF[0]*f[4]*nu[163]+5.7282196186948*rdxF[0]*f[10]*nu[160]; 
  out[57] += 18.44756081437327*rdxF[1]*f[28]*nu[274]+10.5*f[1]*rdxF[1]*nu[274]+10.06230589874905*rdxF[1]*f[57]*nu[254]+12.8086884574495*rdxF[1]*f[8]*nu[254]+7.685213074469698*rdxF[1]*f[28]*nu[244]+6.873863542433758*f[1]*rdxF[1]*nu[244]+5.7282196186948*rdxF[1]*f[8]*nu[240]; 
  out[58] += 18.44756081437327*rdxF[1]*f[29]*nu[274]+10.5*rdxF[1]*f[2]*nu[274]+10.06230589874905*rdxF[1]*f[58]*nu[254]+12.8086884574495*rdxF[1]*f[9]*nu[254]+7.685213074469698*rdxF[1]*f[29]*nu[244]+6.873863542433758*rdxF[1]*f[2]*nu[244]+5.7282196186948*rdxF[1]*f[9]*nu[240]; 
  out[59] += 18.44756081437327*rdxF[1]*f[30]*nu[274]+10.5*rdxF[1]*f[3]*nu[274]+10.06230589874905*rdxF[1]*f[59]*nu[254]+12.8086884574495*rdxF[1]*f[10]*nu[254]+7.685213074469698*rdxF[1]*f[30]*nu[244]+6.873863542433758*rdxF[1]*f[3]*nu[244]+5.7282196186948*rdxF[1]*f[10]*nu[240]+1.14564392373896*rdxF[0]*f[34]*nu[193]+1.677050983124842*rdxF[0]*f[59]*nu[173]+0.7499999999999999*rdxF[0]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[1]*f[36]*nu[274]+1.677050983124842*rdxF[1]*f[60]*nu[254]+0.75*rdxF[1]*f[36]*nu[244]+1.14564392373896*rdxF[0]*f[39]*nu[193]+1.677050983124842*rdxF[0]*f[60]*nu[173]+0.75*rdxF[0]*f[39]*nu[163]; 
  out[61] += 1.14564392373896*rdxF[1]*f[37]*nu[274]+1.677050983124842*rdxF[1]*f[61]*nu[254]+0.75*rdxF[1]*f[37]*nu[244]+1.14564392373896*rdxF[0]*f[40]*nu[193]+1.677050983124842*rdxF[0]*f[61]*nu[173]+0.75*rdxF[0]*f[40]*nu[163]; 
  out[62] += 1.14564392373896*rdxF[1]*f[38]*nu[274]+1.677050983124842*rdxF[1]*f[62]*nu[254]+0.75*rdxF[1]*f[38]*nu[244]+6.708203932499369*rdxF[0]*f[78]*nu[193]+7.6852130744697*rdxF[0]*f[35]*nu[193]+5.031152949374527*rdxF[0]*f[62]*nu[173]+3.75*rdxF[0]*f[16]*nu[173]+3.354101966249684*rdxF[0]*f[35]*nu[163]+1.677050983124842*rdxF[0]*f[16]*nu[160]; 
  out[63] += 6.708203932499369*rdxF[1]*f[79]*nu[274]+7.6852130744697*rdxF[1]*f[35]*nu[274]+5.031152949374527*rdxF[1]*f[63]*nu[254]+3.75*rdxF[1]*f[15]*nu[254]+3.354101966249684*rdxF[1]*f[35]*nu[244]+1.677050983124842*rdxF[1]*f[15]*nu[240]+1.14564392373896*rdxF[0]*f[45]*nu[193]+1.677050983124842*rdxF[0]*f[63]*nu[173]+0.75*rdxF[0]*f[45]*nu[163]; 
  out[64] += 1.14564392373896*rdxF[0]*f[48]*nu[193]+1.677050983124842*rdxF[0]*f[64]*nu[173]+0.7499999999999999*rdxF[0]*f[48]*nu[163]; 
  out[65] += 1.14564392373896*rdxF[0]*f[49]*nu[193]+1.677050983124842*rdxF[0]*f[65]*nu[173]+0.7499999999999999*rdxF[0]*f[49]*nu[163]; 
  out[66] += 18.44756081437327*rdxF[0]*f[38]*nu[193]+10.5*rdxF[0]*f[5]*nu[193]+10.06230589874905*rdxF[0]*f[66]*nu[173]+12.8086884574495*rdxF[0]*f[15]*nu[173]+7.685213074469699*rdxF[0]*f[38]*nu[163]+6.873863542433759*rdxF[0]*f[5]*nu[163]+5.7282196186948*rdxF[0]*f[15]*nu[160]; 
  out[67] += 1.14564392373896*rdxF[1]*f[48]*nu[274]+1.677050983124842*rdxF[1]*f[67]*nu[254]+0.7499999999999999*rdxF[1]*f[48]*nu[244]; 
  out[68] += 1.14564392373896*rdxF[1]*f[49]*nu[274]+1.677050983124842*rdxF[1]*f[68]*nu[254]+0.7499999999999999*rdxF[1]*f[49]*nu[244]; 
  out[69] += 1.14564392373896*rdxF[1]*f[50]*nu[274]+1.677050983124842*rdxF[1]*f[69]*nu[254]+0.7499999999999999*rdxF[1]*f[50]*nu[244]+1.14564392373896*rdxF[0]*f[54]*nu[193]+1.677050983124842*rdxF[0]*f[69]*nu[173]+0.7499999999999999*rdxF[0]*f[54]*nu[163]; 
  out[70] += 1.14564392373896*rdxF[1]*f[51]*nu[274]+1.677050983124842*rdxF[1]*f[70]*nu[254]+0.7499999999999999*rdxF[1]*f[51]*nu[244]+1.14564392373896*rdxF[0]*f[55]*nu[193]+1.677050983124842*rdxF[0]*f[70]*nu[173]+0.7499999999999999*rdxF[0]*f[55]*nu[163]; 
  out[71] += 1.14564392373896*rdxF[1]*f[52]*nu[274]+1.677050983124842*rdxF[1]*f[71]*nu[254]+0.7499999999999999*rdxF[1]*f[52]*nu[244]+18.44756081437327*rdxF[0]*f[43]*nu[193]+10.5*rdxF[0]*f[8]*nu[193]+10.06230589874905*rdxF[0]*f[71]*nu[173]+12.8086884574495*rdxF[0]*f[17]*nu[173]+7.685213074469699*rdxF[0]*f[43]*nu[163]+6.873863542433759*rdxF[0]*f[8]*nu[163]+5.7282196186948*rdxF[0]*f[17]*nu[160]; 
  out[72] += 1.14564392373896*rdxF[1]*f[53]*nu[274]+1.677050983124842*rdxF[1]*f[72]*nu[254]+0.7499999999999999*rdxF[1]*f[53]*nu[244]+18.44756081437327*rdxF[0]*f[44]*nu[193]+10.5*rdxF[0]*f[9]*nu[193]+10.06230589874905*rdxF[0]*f[72]*nu[173]+12.8086884574495*rdxF[0]*f[18]*nu[173]+7.685213074469699*rdxF[0]*f[44]*nu[163]+6.873863542433759*rdxF[0]*f[9]*nu[163]+5.7282196186948*rdxF[0]*f[18]*nu[160]; 
  out[73] += 18.44756081437327*rdxF[1]*f[45]*nu[274]+10.5*rdxF[1]*f[5]*nu[274]+10.06230589874905*rdxF[1]*f[73]*nu[254]+12.8086884574495*rdxF[1]*f[16]*nu[254]+7.685213074469699*rdxF[1]*f[45]*nu[244]+6.873863542433759*rdxF[1]*f[5]*nu[244]+5.7282196186948*rdxF[1]*f[16]*nu[240]; 
  out[74] += 18.44756081437327*rdxF[1]*f[46]*nu[274]+10.5*rdxF[1]*f[6]*nu[274]+10.06230589874905*rdxF[1]*f[74]*nu[254]+12.8086884574495*rdxF[1]*f[17]*nu[254]+7.685213074469699*rdxF[1]*f[46]*nu[244]+6.873863542433759*rdxF[1]*f[6]*nu[244]+5.7282196186948*rdxF[1]*f[17]*nu[240]+1.14564392373896*rdxF[0]*f[57]*nu[193]+1.677050983124842*rdxF[0]*f[74]*nu[173]+0.7499999999999999*rdxF[0]*f[57]*nu[163]; 
  out[75] += 18.44756081437327*rdxF[1]*f[47]*nu[274]+10.5*rdxF[1]*f[7]*nu[274]+10.06230589874905*rdxF[1]*f[75]*nu[254]+12.8086884574495*rdxF[1]*f[18]*nu[254]+7.685213074469699*rdxF[1]*f[47]*nu[244]+6.873863542433759*rdxF[1]*f[7]*nu[244]+5.7282196186948*rdxF[1]*f[18]*nu[240]+1.14564392373896*rdxF[0]*f[58]*nu[193]+1.677050983124842*rdxF[0]*f[75]*nu[173]+0.7499999999999999*rdxF[0]*f[58]*nu[163]; 
  out[76] += 1.14564392373896*rdxF[1]*f[64]*nu[274]+1.677050983124842*rdxF[1]*f[76]*nu[254]+0.7499999999999999*rdxF[1]*f[64]*nu[244]+1.14564392373896*rdxF[0]*f[67]*nu[193]+1.677050983124842*rdxF[0]*f[76]*nu[173]+0.7499999999999999*rdxF[0]*f[67]*nu[163]; 
  out[77] += 1.14564392373896*rdxF[1]*f[65]*nu[274]+1.677050983124842*rdxF[1]*f[77]*nu[254]+0.7499999999999999*rdxF[1]*f[65]*nu[244]+1.14564392373896*rdxF[0]*f[68]*nu[193]+1.677050983124842*rdxF[0]*f[77]*nu[173]+0.7499999999999999*rdxF[0]*f[68]*nu[163]; 
  out[78] += 1.14564392373896*rdxF[1]*f[66]*nu[274]+1.677050983124842*rdxF[1]*f[78]*nu[254]+0.7499999999999999*rdxF[1]*f[66]*nu[244]+18.44756081437327*rdxF[0]*f[62]*nu[193]+10.5*rdxF[0]*f[16]*nu[193]+10.06230589874905*rdxF[0]*f[78]*nu[173]+12.8086884574495*rdxF[0]*f[35]*nu[173]+7.685213074469698*rdxF[0]*f[62]*nu[163]+6.873863542433758*rdxF[0]*f[16]*nu[163]+5.7282196186948*rdxF[0]*f[35]*nu[160]; 
  out[79] += 18.44756081437327*rdxF[1]*f[63]*nu[274]+10.5*rdxF[1]*f[15]*nu[274]+10.06230589874905*rdxF[1]*f[79]*nu[254]+12.8086884574495*rdxF[1]*f[35]*nu[254]+7.685213074469698*rdxF[1]*f[63]*nu[244]+6.873863542433758*rdxF[1]*f[15]*nu[244]+5.7282196186948*rdxF[1]*f[35]*nu[240]+1.14564392373896*rdxF[0]*f[73]*nu[193]+1.677050983124842*rdxF[0]*f[79]*nu[173]+0.7499999999999999*rdxF[0]*f[73]*nu[163]; 

  return (rdxF[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+rdxF[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs134(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[3]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 
  rdxF[2] = 4.0/(dx[3]*dx[3]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[3] += 2.5617376914899*rdxF[1]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[3]*nu[173]+0.75*f[0]*rdxF[1]*nu[163]; 
  out[4] += 2.5617376914899*rdxF[2]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[4]*nu[254]+0.75*f[0]*rdxF[2]*nu[244]; 
  out[5] += 2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[1]*f[23]*nu[193]+1.14564392373896*f[1]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[6]*nu[173]+0.75*f[1]*rdxF[1]*nu[163]+2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 2.5617376914899*rdxF[1]*f[24]*nu[193]+1.14564392373896*rdxF[1]*f[2]*nu[193]+1.677050983124842*rdxF[1]*f[7]*nu[173]+0.75*rdxF[1]*f[2]*nu[163]; 
  out[8] += 2.5617376914899*rdxF[2]*f[28]*nu[274]+1.14564392373896*f[1]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[8]*nu[254]+0.75*f[1]*rdxF[2]*nu[244]+2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 2.5617376914899*rdxF[2]*f[29]*nu[274]+1.14564392373896*f[2]*rdxF[2]*nu[274]+1.677050983124842*rdxF[2]*f[9]*nu[254]+0.75*f[2]*rdxF[2]*nu[244]; 
  out[10] += 2.5617376914899*rdxF[2]*f[30]*nu[274]+1.14564392373896*rdxF[2]*f[3]*nu[274]+1.677050983124842*rdxF[2]*f[10]*nu[254]+0.75*rdxF[2]*f[3]*nu[244]+2.5617376914899*rdxF[1]*f[27]*nu[193]+1.14564392373896*rdxF[1]*f[4]*nu[193]+1.677050983124842*rdxF[1]*f[10]*nu[173]+0.75*rdxF[1]*f[4]*nu[163]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[13] += 6.708203932499369*rdxF[1]*f[33]*nu[193]+7.685213074469699*rdxF[1]*f[3]*nu[193]+5.031152949374527*rdxF[1]*f[13]*nu[173]+3.75*f[0]*rdxF[1]*nu[173]+3.354101966249685*rdxF[1]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[1]*nu[160]; 
  out[14] += 6.708203932499369*rdxF[2]*f[34]*nu[274]+7.685213074469699*rdxF[2]*f[4]*nu[274]+5.031152949374527*rdxF[2]*f[14]*nu[254]+3.75*f[0]*rdxF[2]*nu[254]+3.354101966249685*rdxF[2]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[2]*nu[240]; 
  out[15] += 2.5617376914899*rdxF[1]*f[38]*nu[193]+1.14564392373896*rdxF[1]*f[5]*nu[193]+1.677050983124842*rdxF[1]*f[15]*nu[173]+0.75*rdxF[1]*f[5]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[2]*f[45]*nu[274]+1.14564392373896*rdxF[2]*f[5]*nu[274]+1.677050983124842*rdxF[2]*f[16]*nu[254]+0.75*rdxF[2]*f[5]*nu[244]+2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[2]*f[46]*nu[274]+1.14564392373896*rdxF[2]*f[6]*nu[274]+1.677050983124842*rdxF[2]*f[17]*nu[254]+0.75*rdxF[2]*f[6]*nu[244]+2.5617376914899*rdxF[1]*f[43]*nu[193]+1.14564392373896*rdxF[1]*f[8]*nu[193]+1.677050983124842*rdxF[1]*f[17]*nu[173]+0.75*rdxF[1]*f[8]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 2.5617376914899*rdxF[2]*f[47]*nu[274]+1.14564392373896*rdxF[2]*f[7]*nu[274]+1.677050983124842*rdxF[2]*f[18]*nu[254]+0.75*rdxF[2]*f[7]*nu[244]+2.5617376914899*rdxF[1]*f[44]*nu[193]+1.14564392373896*rdxF[1]*f[9]*nu[193]+1.677050983124842*rdxF[1]*f[18]*nu[173]+0.75*rdxF[1]*f[9]*nu[163]; 
  out[19] += 6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.14564392373896*rdxF[1]*f[11]*nu[193]+1.677050983124842*rdxF[1]*f[21]*nu[173]+0.75*rdxF[1]*f[11]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.14564392373896*rdxF[1]*f[12]*nu[193]+1.677050983124842*rdxF[1]*f[22]*nu[173]+0.75*rdxF[1]*f[12]*nu[163]; 
  out[23] += 6.708203932499369*rdxF[1]*f[52]*nu[193]+7.6852130744697*rdxF[1]*f[6]*nu[193]+5.031152949374527*rdxF[1]*f[23]*nu[173]+3.75*f[1]*rdxF[1]*nu[173]+3.354101966249684*rdxF[1]*f[6]*nu[163]+1.677050983124842*f[1]*rdxF[1]*nu[160]+1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 6.708203932499369*rdxF[1]*f[53]*nu[193]+7.6852130744697*rdxF[1]*f[7]*nu[193]+5.031152949374527*rdxF[1]*f[24]*nu[173]+3.75*rdxF[1]*f[2]*nu[173]+3.354101966249684*rdxF[1]*f[7]*nu[163]+1.677050983124842*rdxF[1]*f[2]*nu[160]; 
  out[25] += 1.14564392373896*rdxF[2]*f[11]*nu[274]+1.677050983124842*rdxF[2]*f[25]*nu[254]+0.75*rdxF[2]*f[11]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.14564392373896*rdxF[2]*f[12]*nu[274]+1.677050983124842*rdxF[2]*f[26]*nu[254]+0.75*rdxF[2]*f[12]*nu[244]; 
  out[27] += 1.14564392373896*rdxF[2]*f[13]*nu[274]+1.677050983124842*rdxF[2]*f[27]*nu[254]+0.75*rdxF[2]*f[13]*nu[244]+6.708203932499369*rdxF[1]*f[56]*nu[193]+7.6852130744697*rdxF[1]*f[10]*nu[193]+5.031152949374527*rdxF[1]*f[27]*nu[173]+3.75*rdxF[1]*f[4]*nu[173]+3.354101966249684*rdxF[1]*f[10]*nu[163]+1.677050983124842*rdxF[1]*f[4]*nu[160]; 
  out[28] += 6.708203932499369*rdxF[2]*f[57]*nu[274]+7.6852130744697*rdxF[2]*f[8]*nu[274]+5.031152949374527*rdxF[2]*f[28]*nu[254]+3.75*f[1]*rdxF[2]*nu[254]+3.354101966249684*rdxF[2]*f[8]*nu[244]+1.677050983124842*f[1]*rdxF[2]*nu[240]+1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 6.708203932499369*rdxF[2]*f[58]*nu[274]+7.6852130744697*rdxF[2]*f[9]*nu[274]+5.031152949374527*rdxF[2]*f[29]*nu[254]+3.75*f[2]*rdxF[2]*nu[254]+3.354101966249684*rdxF[2]*f[9]*nu[244]+1.677050983124842*f[2]*rdxF[2]*nu[240]; 
  out[30] += 6.708203932499369*rdxF[2]*f[59]*nu[274]+7.6852130744697*rdxF[2]*f[10]*nu[274]+5.031152949374527*rdxF[2]*f[30]*nu[254]+3.75*rdxF[2]*f[3]*nu[254]+3.354101966249684*rdxF[2]*f[10]*nu[244]+1.677050983124842*rdxF[2]*f[3]*nu[240]+1.14564392373896*rdxF[1]*f[14]*nu[193]+1.677050983124842*rdxF[1]*f[30]*nu[173]+0.75*rdxF[1]*f[14]*nu[163]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[33] += 18.44756081437327*rdxF[1]*f[13]*nu[193]+10.5*f[0]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[33]*nu[173]+12.8086884574495*rdxF[1]*f[3]*nu[173]+7.685213074469699*rdxF[1]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[3]*nu[160]; 
  out[34] += 18.44756081437327*rdxF[2]*f[14]*nu[274]+10.5*f[0]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[34]*nu[254]+12.8086884574495*rdxF[2]*f[4]*nu[254]+7.685213074469699*rdxF[2]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[2]*f[63]*nu[274]+1.14564392373896*rdxF[2]*f[15]*nu[274]+1.677050983124842*rdxF[2]*f[35]*nu[254]+0.75*rdxF[2]*f[15]*nu[244]+2.5617376914899*rdxF[1]*f[62]*nu[193]+1.14564392373896*rdxF[1]*f[16]*nu[193]+1.677050983124842*rdxF[1]*f[35]*nu[173]+0.75*rdxF[1]*f[16]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 1.14564392373896*rdxF[1]*f[19]*nu[193]+1.677050983124842*rdxF[1]*f[36]*nu[173]+0.75*rdxF[1]*f[19]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 1.14564392373896*rdxF[1]*f[20]*nu[193]+1.677050983124842*rdxF[1]*f[37]*nu[173]+0.75*rdxF[1]*f[20]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 6.708203932499369*rdxF[1]*f[66]*nu[193]+7.685213074469699*rdxF[1]*f[15]*nu[193]+5.031152949374527*rdxF[1]*f[38]*nu[173]+3.75*rdxF[1]*f[5]*nu[173]+3.354101966249685*rdxF[1]*f[15]*nu[163]+1.677050983124842*rdxF[1]*f[5]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 1.14564392373896*rdxF[2]*f[19]*nu[274]+1.677050983124842*rdxF[2]*f[39]*nu[254]+0.75*rdxF[2]*f[19]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 1.14564392373896*rdxF[2]*f[20]*nu[274]+1.677050983124842*rdxF[2]*f[40]*nu[254]+0.75*rdxF[2]*f[20]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 1.14564392373896*rdxF[2]*f[21]*nu[274]+1.677050983124842*rdxF[2]*f[41]*nu[254]+0.75*rdxF[2]*f[21]*nu[244]+1.14564392373896*rdxF[1]*f[25]*nu[193]+1.677050983124842*rdxF[1]*f[41]*nu[173]+0.75*rdxF[1]*f[25]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[42] += 1.14564392373896*rdxF[2]*f[22]*nu[274]+1.677050983124842*rdxF[2]*f[42]*nu[254]+0.75*rdxF[2]*f[22]*nu[244]+1.14564392373896*rdxF[1]*f[26]*nu[193]+1.677050983124842*rdxF[1]*f[42]*nu[173]+0.75*rdxF[1]*f[26]*nu[163]; 
  out[43] += 1.14564392373896*rdxF[2]*f[23]*nu[274]+1.677050983124842*rdxF[2]*f[43]*nu[254]+0.75*rdxF[2]*f[23]*nu[244]+6.708203932499369*rdxF[1]*f[71]*nu[193]+7.685213074469699*rdxF[1]*f[17]*nu[193]+5.031152949374527*rdxF[1]*f[43]*nu[173]+3.75*rdxF[1]*f[8]*nu[173]+3.354101966249685*rdxF[1]*f[17]*nu[163]+1.677050983124842*rdxF[1]*f[8]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[44] += 1.14564392373896*rdxF[2]*f[24]*nu[274]+1.677050983124842*rdxF[2]*f[44]*nu[254]+0.75*rdxF[2]*f[24]*nu[244]+6.708203932499369*rdxF[1]*f[72]*nu[193]+7.685213074469699*rdxF[1]*f[18]*nu[193]+5.031152949374527*rdxF[1]*f[44]*nu[173]+3.75*rdxF[1]*f[9]*nu[173]+3.354101966249685*rdxF[1]*f[18]*nu[163]+1.677050983124842*rdxF[1]*f[9]*nu[160]; 
  out[45] += 6.708203932499369*rdxF[2]*f[73]*nu[274]+7.685213074469699*rdxF[2]*f[16]*nu[274]+5.031152949374527*rdxF[2]*f[45]*nu[254]+3.75*rdxF[2]*f[5]*nu[254]+3.354101966249685*rdxF[2]*f[16]*nu[244]+1.677050983124842*rdxF[2]*f[5]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 6.708203932499369*rdxF[2]*f[74]*nu[274]+7.685213074469699*rdxF[2]*f[17]*nu[274]+5.031152949374527*rdxF[2]*f[46]*nu[254]+3.75*rdxF[2]*f[6]*nu[254]+3.354101966249685*rdxF[2]*f[17]*nu[244]+1.677050983124842*rdxF[2]*f[6]*nu[240]+1.14564392373896*rdxF[1]*f[28]*nu[193]+1.677050983124842*rdxF[1]*f[46]*nu[173]+0.75*rdxF[1]*f[28]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[47] += 6.708203932499369*rdxF[2]*f[75]*nu[274]+7.685213074469699*rdxF[2]*f[18]*nu[274]+5.031152949374527*rdxF[2]*f[47]*nu[254]+3.75*rdxF[2]*f[7]*nu[254]+3.354101966249685*rdxF[2]*f[18]*nu[244]+1.677050983124842*rdxF[2]*f[7]*nu[240]+1.14564392373896*rdxF[1]*f[29]*nu[193]+1.677050983124842*rdxF[1]*f[47]*nu[173]+0.75*rdxF[1]*f[29]*nu[163]; 
  out[48] += 10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 1.14564392373896*rdxF[1]*f[31]*nu[193]+1.677050983124842*rdxF[1]*f[50]*nu[173]+0.7499999999999999*rdxF[1]*f[31]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[51] += 1.14564392373896*rdxF[1]*f[32]*nu[193]+1.677050983124842*rdxF[1]*f[51]*nu[173]+0.7499999999999999*rdxF[1]*f[32]*nu[163]; 
  out[52] += 18.44756081437327*rdxF[1]*f[23]*nu[193]+10.5*f[1]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[52]*nu[173]+12.8086884574495*rdxF[1]*f[6]*nu[173]+7.685213074469698*rdxF[1]*f[23]*nu[163]+6.873863542433758*f[1]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[6]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[53] += 18.44756081437327*rdxF[1]*f[24]*nu[193]+10.5*rdxF[1]*f[2]*nu[193]+10.06230589874905*rdxF[1]*f[53]*nu[173]+12.8086884574495*rdxF[1]*f[7]*nu[173]+7.685213074469698*rdxF[1]*f[24]*nu[163]+6.873863542433758*rdxF[1]*f[2]*nu[163]+5.7282196186948*rdxF[1]*f[7]*nu[160]; 
  out[54] += 1.14564392373896*rdxF[2]*f[31]*nu[274]+1.677050983124842*rdxF[2]*f[54]*nu[254]+0.7499999999999999*rdxF[2]*f[31]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[55] += 1.14564392373896*rdxF[2]*f[32]*nu[274]+1.677050983124842*rdxF[2]*f[55]*nu[254]+0.7499999999999999*rdxF[2]*f[32]*nu[244]; 
  out[56] += 1.14564392373896*rdxF[2]*f[33]*nu[274]+1.677050983124842*rdxF[2]*f[56]*nu[254]+0.7499999999999999*rdxF[2]*f[33]*nu[244]+18.44756081437327*rdxF[1]*f[27]*nu[193]+10.5*rdxF[1]*f[4]*nu[193]+10.06230589874905*rdxF[1]*f[56]*nu[173]+12.8086884574495*rdxF[1]*f[10]*nu[173]+7.685213074469698*rdxF[1]*f[27]*nu[163]+6.873863542433758*rdxF[1]*f[4]*nu[163]+5.7282196186948*rdxF[1]*f[10]*nu[160]; 
  out[57] += 18.44756081437327*rdxF[2]*f[28]*nu[274]+10.5*f[1]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[57]*nu[254]+12.8086884574495*rdxF[2]*f[8]*nu[254]+7.685213074469698*rdxF[2]*f[28]*nu[244]+6.873863542433758*f[1]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[8]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[58] += 18.44756081437327*rdxF[2]*f[29]*nu[274]+10.5*f[2]*rdxF[2]*nu[274]+10.06230589874905*rdxF[2]*f[58]*nu[254]+12.8086884574495*rdxF[2]*f[9]*nu[254]+7.685213074469698*rdxF[2]*f[29]*nu[244]+6.873863542433758*f[2]*rdxF[2]*nu[244]+5.7282196186948*rdxF[2]*f[9]*nu[240]; 
  out[59] += 18.44756081437327*rdxF[2]*f[30]*nu[274]+10.5*rdxF[2]*f[3]*nu[274]+10.06230589874905*rdxF[2]*f[59]*nu[254]+12.8086884574495*rdxF[2]*f[10]*nu[254]+7.685213074469698*rdxF[2]*f[30]*nu[244]+6.873863542433758*rdxF[2]*f[3]*nu[244]+5.7282196186948*rdxF[2]*f[10]*nu[240]+1.14564392373896*rdxF[1]*f[34]*nu[193]+1.677050983124842*rdxF[1]*f[59]*nu[173]+0.7499999999999999*rdxF[1]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[2]*f[36]*nu[274]+1.677050983124842*rdxF[2]*f[60]*nu[254]+0.75*rdxF[2]*f[36]*nu[244]+1.14564392373896*rdxF[1]*f[39]*nu[193]+1.677050983124842*rdxF[1]*f[60]*nu[173]+0.75*rdxF[1]*f[39]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 1.14564392373896*rdxF[2]*f[37]*nu[274]+1.677050983124842*rdxF[2]*f[61]*nu[254]+0.75*rdxF[2]*f[37]*nu[244]+1.14564392373896*rdxF[1]*f[40]*nu[193]+1.677050983124842*rdxF[1]*f[61]*nu[173]+0.75*rdxF[1]*f[40]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 1.14564392373896*rdxF[2]*f[38]*nu[274]+1.677050983124842*rdxF[2]*f[62]*nu[254]+0.75*rdxF[2]*f[38]*nu[244]+6.708203932499369*rdxF[1]*f[78]*nu[193]+7.6852130744697*rdxF[1]*f[35]*nu[193]+5.031152949374527*rdxF[1]*f[62]*nu[173]+3.75*rdxF[1]*f[16]*nu[173]+3.354101966249684*rdxF[1]*f[35]*nu[163]+1.677050983124842*rdxF[1]*f[16]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 6.708203932499369*rdxF[2]*f[79]*nu[274]+7.6852130744697*rdxF[2]*f[35]*nu[274]+5.031152949374527*rdxF[2]*f[63]*nu[254]+3.75*rdxF[2]*f[15]*nu[254]+3.354101966249684*rdxF[2]*f[35]*nu[244]+1.677050983124842*rdxF[2]*f[15]*nu[240]+1.14564392373896*rdxF[1]*f[45]*nu[193]+1.677050983124842*rdxF[1]*f[63]*nu[173]+0.75*rdxF[1]*f[45]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 1.14564392373896*rdxF[1]*f[48]*nu[193]+1.677050983124842*rdxF[1]*f[64]*nu[173]+0.7499999999999999*rdxF[1]*f[48]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 1.14564392373896*rdxF[1]*f[49]*nu[193]+1.677050983124842*rdxF[1]*f[65]*nu[173]+0.7499999999999999*rdxF[1]*f[49]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 18.44756081437327*rdxF[1]*f[38]*nu[193]+10.5*rdxF[1]*f[5]*nu[193]+10.06230589874905*rdxF[1]*f[66]*nu[173]+12.8086884574495*rdxF[1]*f[15]*nu[173]+7.685213074469699*rdxF[1]*f[38]*nu[163]+6.873863542433759*rdxF[1]*f[5]*nu[163]+5.7282196186948*rdxF[1]*f[15]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 1.14564392373896*rdxF[2]*f[48]*nu[274]+1.677050983124842*rdxF[2]*f[67]*nu[254]+0.7499999999999999*rdxF[2]*f[48]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 1.14564392373896*rdxF[2]*f[49]*nu[274]+1.677050983124842*rdxF[2]*f[68]*nu[254]+0.7499999999999999*rdxF[2]*f[49]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 1.14564392373896*rdxF[2]*f[50]*nu[274]+1.677050983124842*rdxF[2]*f[69]*nu[254]+0.7499999999999999*rdxF[2]*f[50]*nu[244]+1.14564392373896*rdxF[1]*f[54]*nu[193]+1.677050983124842*rdxF[1]*f[69]*nu[173]+0.7499999999999999*rdxF[1]*f[54]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[70] += 1.14564392373896*rdxF[2]*f[51]*nu[274]+1.677050983124842*rdxF[2]*f[70]*nu[254]+0.7499999999999999*rdxF[2]*f[51]*nu[244]+1.14564392373896*rdxF[1]*f[55]*nu[193]+1.677050983124842*rdxF[1]*f[70]*nu[173]+0.7499999999999999*rdxF[1]*f[55]*nu[163]; 
  out[71] += 1.14564392373896*rdxF[2]*f[52]*nu[274]+1.677050983124842*rdxF[2]*f[71]*nu[254]+0.7499999999999999*rdxF[2]*f[52]*nu[244]+18.44756081437327*rdxF[1]*f[43]*nu[193]+10.5*rdxF[1]*f[8]*nu[193]+10.06230589874905*rdxF[1]*f[71]*nu[173]+12.8086884574495*rdxF[1]*f[17]*nu[173]+7.685213074469699*rdxF[1]*f[43]*nu[163]+6.873863542433759*rdxF[1]*f[8]*nu[163]+5.7282196186948*rdxF[1]*f[17]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[72] += 1.14564392373896*rdxF[2]*f[53]*nu[274]+1.677050983124842*rdxF[2]*f[72]*nu[254]+0.7499999999999999*rdxF[2]*f[53]*nu[244]+18.44756081437327*rdxF[1]*f[44]*nu[193]+10.5*rdxF[1]*f[9]*nu[193]+10.06230589874905*rdxF[1]*f[72]*nu[173]+12.8086884574495*rdxF[1]*f[18]*nu[173]+7.685213074469699*rdxF[1]*f[44]*nu[163]+6.873863542433759*rdxF[1]*f[9]*nu[163]+5.7282196186948*rdxF[1]*f[18]*nu[160]; 
  out[73] += 18.44756081437327*rdxF[2]*f[45]*nu[274]+10.5*rdxF[2]*f[5]*nu[274]+10.06230589874905*rdxF[2]*f[73]*nu[254]+12.8086884574495*rdxF[2]*f[16]*nu[254]+7.685213074469699*rdxF[2]*f[45]*nu[244]+6.873863542433759*rdxF[2]*f[5]*nu[244]+5.7282196186948*rdxF[2]*f[16]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 18.44756081437327*rdxF[2]*f[46]*nu[274]+10.5*rdxF[2]*f[6]*nu[274]+10.06230589874905*rdxF[2]*f[74]*nu[254]+12.8086884574495*rdxF[2]*f[17]*nu[254]+7.685213074469699*rdxF[2]*f[46]*nu[244]+6.873863542433759*rdxF[2]*f[6]*nu[244]+5.7282196186948*rdxF[2]*f[17]*nu[240]+1.14564392373896*rdxF[1]*f[57]*nu[193]+1.677050983124842*rdxF[1]*f[74]*nu[173]+0.7499999999999999*rdxF[1]*f[57]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[75] += 18.44756081437327*rdxF[2]*f[47]*nu[274]+10.5*rdxF[2]*f[7]*nu[274]+10.06230589874905*rdxF[2]*f[75]*nu[254]+12.8086884574495*rdxF[2]*f[18]*nu[254]+7.685213074469699*rdxF[2]*f[47]*nu[244]+6.873863542433759*rdxF[2]*f[7]*nu[244]+5.7282196186948*rdxF[2]*f[18]*nu[240]+1.14564392373896*rdxF[1]*f[58]*nu[193]+1.677050983124842*rdxF[1]*f[75]*nu[173]+0.7499999999999999*rdxF[1]*f[58]*nu[163]; 
  out[76] += 1.14564392373896*rdxF[2]*f[64]*nu[274]+1.677050983124842*rdxF[2]*f[76]*nu[254]+0.7499999999999999*rdxF[2]*f[64]*nu[244]+1.14564392373896*rdxF[1]*f[67]*nu[193]+1.677050983124842*rdxF[1]*f[76]*nu[173]+0.7499999999999999*rdxF[1]*f[67]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 1.14564392373896*rdxF[2]*f[65]*nu[274]+1.677050983124842*rdxF[2]*f[77]*nu[254]+0.7499999999999999*rdxF[2]*f[65]*nu[244]+1.14564392373896*rdxF[1]*f[68]*nu[193]+1.677050983124842*rdxF[1]*f[77]*nu[173]+0.7499999999999999*rdxF[1]*f[68]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 1.14564392373896*rdxF[2]*f[66]*nu[274]+1.677050983124842*rdxF[2]*f[78]*nu[254]+0.7499999999999999*rdxF[2]*f[66]*nu[244]+18.44756081437327*rdxF[1]*f[62]*nu[193]+10.5*rdxF[1]*f[16]*nu[193]+10.06230589874905*rdxF[1]*f[78]*nu[173]+12.8086884574495*rdxF[1]*f[35]*nu[173]+7.685213074469698*rdxF[1]*f[62]*nu[163]+6.873863542433758*rdxF[1]*f[16]*nu[163]+5.7282196186948*rdxF[1]*f[35]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 18.44756081437327*rdxF[2]*f[63]*nu[274]+10.5*rdxF[2]*f[15]*nu[274]+10.06230589874905*rdxF[2]*f[79]*nu[254]+12.8086884574495*rdxF[2]*f[35]*nu[254]+7.685213074469698*rdxF[2]*f[63]*nu[244]+6.873863542433758*rdxF[2]*f[15]*nu[244]+5.7282196186948*rdxF[2]*f[35]*nu[240]+1.14564392373896*rdxF[1]*f[73]*nu[193]+1.677050983124842*rdxF[1]*f[79]*nu[173]+0.7499999999999999*rdxF[1]*f[73]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+rdxF[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs13(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[2]*dx[2]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[3] += 2.5617376914899*rdxF[1]*f[13]*nu[193]+1.14564392373896*f[0]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[3]*nu[173]+0.75*f[0]*rdxF[1]*nu[163]; 
  out[5] += 2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[1]*f[23]*nu[193]+1.14564392373896*f[1]*rdxF[1]*nu[193]+1.677050983124842*rdxF[1]*f[6]*nu[173]+0.75*f[1]*rdxF[1]*nu[163]+2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[7] += 2.5617376914899*rdxF[1]*f[24]*nu[193]+1.14564392373896*rdxF[1]*f[2]*nu[193]+1.677050983124842*rdxF[1]*f[7]*nu[173]+0.75*rdxF[1]*f[2]*nu[163]; 
  out[8] += 2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[10] += 2.5617376914899*rdxF[1]*f[27]*nu[193]+1.14564392373896*rdxF[1]*f[4]*nu[193]+1.677050983124842*rdxF[1]*f[10]*nu[173]+0.75*rdxF[1]*f[4]*nu[163]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[13] += 6.708203932499369*rdxF[1]*f[33]*nu[193]+7.685213074469699*rdxF[1]*f[3]*nu[193]+5.031152949374527*rdxF[1]*f[13]*nu[173]+3.75*f[0]*rdxF[1]*nu[173]+3.354101966249685*rdxF[1]*f[3]*nu[163]+1.677050983124842*f[0]*rdxF[1]*nu[160]; 
  out[15] += 2.5617376914899*rdxF[1]*f[38]*nu[193]+1.14564392373896*rdxF[1]*f[5]*nu[193]+1.677050983124842*rdxF[1]*f[15]*nu[173]+0.75*rdxF[1]*f[5]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[1]*f[43]*nu[193]+1.14564392373896*rdxF[1]*f[8]*nu[193]+1.677050983124842*rdxF[1]*f[17]*nu[173]+0.75*rdxF[1]*f[8]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 2.5617376914899*rdxF[1]*f[44]*nu[193]+1.14564392373896*rdxF[1]*f[9]*nu[193]+1.677050983124842*rdxF[1]*f[18]*nu[173]+0.75*rdxF[1]*f[9]*nu[163]; 
  out[19] += 6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 1.14564392373896*rdxF[1]*f[11]*nu[193]+1.677050983124842*rdxF[1]*f[21]*nu[173]+0.75*rdxF[1]*f[11]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[22] += 1.14564392373896*rdxF[1]*f[12]*nu[193]+1.677050983124842*rdxF[1]*f[22]*nu[173]+0.75*rdxF[1]*f[12]*nu[163]; 
  out[23] += 6.708203932499369*rdxF[1]*f[52]*nu[193]+7.6852130744697*rdxF[1]*f[6]*nu[193]+5.031152949374527*rdxF[1]*f[23]*nu[173]+3.75*f[1]*rdxF[1]*nu[173]+3.354101966249684*rdxF[1]*f[6]*nu[163]+1.677050983124842*f[1]*rdxF[1]*nu[160]+1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[24] += 6.708203932499369*rdxF[1]*f[53]*nu[193]+7.6852130744697*rdxF[1]*f[7]*nu[193]+5.031152949374527*rdxF[1]*f[24]*nu[173]+3.75*rdxF[1]*f[2]*nu[173]+3.354101966249684*rdxF[1]*f[7]*nu[163]+1.677050983124842*rdxF[1]*f[2]*nu[160]; 
  out[25] += 6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[27] += 6.708203932499369*rdxF[1]*f[56]*nu[193]+7.6852130744697*rdxF[1]*f[10]*nu[193]+5.031152949374527*rdxF[1]*f[27]*nu[173]+3.75*rdxF[1]*f[4]*nu[173]+3.354101966249684*rdxF[1]*f[10]*nu[163]+1.677050983124842*rdxF[1]*f[4]*nu[160]; 
  out[28] += 1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[30] += 1.14564392373896*rdxF[1]*f[14]*nu[193]+1.677050983124842*rdxF[1]*f[30]*nu[173]+0.75*rdxF[1]*f[14]*nu[163]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[33] += 18.44756081437327*rdxF[1]*f[13]*nu[193]+10.5*f[0]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[33]*nu[173]+12.8086884574495*rdxF[1]*f[3]*nu[173]+7.685213074469699*rdxF[1]*f[13]*nu[163]+6.873863542433759*f[0]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[3]*nu[160]; 
  out[35] += 2.5617376914899*rdxF[1]*f[62]*nu[193]+1.14564392373896*rdxF[1]*f[16]*nu[193]+1.677050983124842*rdxF[1]*f[35]*nu[173]+0.75*rdxF[1]*f[16]*nu[163]+2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 1.14564392373896*rdxF[1]*f[19]*nu[193]+1.677050983124842*rdxF[1]*f[36]*nu[173]+0.75*rdxF[1]*f[19]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 1.14564392373896*rdxF[1]*f[20]*nu[193]+1.677050983124842*rdxF[1]*f[37]*nu[173]+0.75*rdxF[1]*f[20]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 6.708203932499369*rdxF[1]*f[66]*nu[193]+7.685213074469699*rdxF[1]*f[15]*nu[193]+5.031152949374527*rdxF[1]*f[38]*nu[173]+3.75*rdxF[1]*f[5]*nu[173]+3.354101966249685*rdxF[1]*f[15]*nu[163]+1.677050983124842*rdxF[1]*f[5]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 1.14564392373896*rdxF[1]*f[25]*nu[193]+1.677050983124842*rdxF[1]*f[41]*nu[173]+0.75*rdxF[1]*f[25]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[42] += 1.14564392373896*rdxF[1]*f[26]*nu[193]+1.677050983124842*rdxF[1]*f[42]*nu[173]+0.75*rdxF[1]*f[26]*nu[163]; 
  out[43] += 6.708203932499369*rdxF[1]*f[71]*nu[193]+7.685213074469699*rdxF[1]*f[17]*nu[193]+5.031152949374527*rdxF[1]*f[43]*nu[173]+3.75*rdxF[1]*f[8]*nu[173]+3.354101966249685*rdxF[1]*f[17]*nu[163]+1.677050983124842*rdxF[1]*f[8]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[44] += 6.708203932499369*rdxF[1]*f[72]*nu[193]+7.685213074469699*rdxF[1]*f[18]*nu[193]+5.031152949374527*rdxF[1]*f[44]*nu[173]+3.75*rdxF[1]*f[9]*nu[173]+3.354101966249685*rdxF[1]*f[18]*nu[163]+1.677050983124842*rdxF[1]*f[9]*nu[160]; 
  out[45] += 1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 1.14564392373896*rdxF[1]*f[28]*nu[193]+1.677050983124842*rdxF[1]*f[46]*nu[173]+0.75*rdxF[1]*f[28]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[47] += 1.14564392373896*rdxF[1]*f[29]*nu[193]+1.677050983124842*rdxF[1]*f[47]*nu[173]+0.75*rdxF[1]*f[29]*nu[163]; 
  out[48] += 10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 1.14564392373896*rdxF[1]*f[31]*nu[193]+1.677050983124842*rdxF[1]*f[50]*nu[173]+0.7499999999999999*rdxF[1]*f[31]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[51] += 1.14564392373896*rdxF[1]*f[32]*nu[193]+1.677050983124842*rdxF[1]*f[51]*nu[173]+0.7499999999999999*rdxF[1]*f[32]*nu[163]; 
  out[52] += 18.44756081437327*rdxF[1]*f[23]*nu[193]+10.5*f[1]*rdxF[1]*nu[193]+10.06230589874905*rdxF[1]*f[52]*nu[173]+12.8086884574495*rdxF[1]*f[6]*nu[173]+7.685213074469698*rdxF[1]*f[23]*nu[163]+6.873863542433758*f[1]*rdxF[1]*nu[163]+5.7282196186948*rdxF[1]*f[6]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[53] += 18.44756081437327*rdxF[1]*f[24]*nu[193]+10.5*rdxF[1]*f[2]*nu[193]+10.06230589874905*rdxF[1]*f[53]*nu[173]+12.8086884574495*rdxF[1]*f[7]*nu[173]+7.685213074469698*rdxF[1]*f[24]*nu[163]+6.873863542433758*rdxF[1]*f[2]*nu[163]+5.7282196186948*rdxF[1]*f[7]*nu[160]; 
  out[54] += 10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[56] += 18.44756081437327*rdxF[1]*f[27]*nu[193]+10.5*rdxF[1]*f[4]*nu[193]+10.06230589874905*rdxF[1]*f[56]*nu[173]+12.8086884574495*rdxF[1]*f[10]*nu[173]+7.685213074469698*rdxF[1]*f[27]*nu[163]+6.873863542433758*rdxF[1]*f[4]*nu[163]+5.7282196186948*rdxF[1]*f[10]*nu[160]; 
  out[57] += 1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[59] += 1.14564392373896*rdxF[1]*f[34]*nu[193]+1.677050983124842*rdxF[1]*f[59]*nu[173]+0.7499999999999999*rdxF[1]*f[34]*nu[163]; 
  out[60] += 1.14564392373896*rdxF[1]*f[39]*nu[193]+1.677050983124842*rdxF[1]*f[60]*nu[173]+0.75*rdxF[1]*f[39]*nu[163]+6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 1.14564392373896*rdxF[1]*f[40]*nu[193]+1.677050983124842*rdxF[1]*f[61]*nu[173]+0.75*rdxF[1]*f[40]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 6.708203932499369*rdxF[1]*f[78]*nu[193]+7.6852130744697*rdxF[1]*f[35]*nu[193]+5.031152949374527*rdxF[1]*f[62]*nu[173]+3.75*rdxF[1]*f[16]*nu[173]+3.354101966249684*rdxF[1]*f[35]*nu[163]+1.677050983124842*rdxF[1]*f[16]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 1.14564392373896*rdxF[1]*f[45]*nu[193]+1.677050983124842*rdxF[1]*f[63]*nu[173]+0.75*rdxF[1]*f[45]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 1.14564392373896*rdxF[1]*f[48]*nu[193]+1.677050983124842*rdxF[1]*f[64]*nu[173]+0.7499999999999999*rdxF[1]*f[48]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 1.14564392373896*rdxF[1]*f[49]*nu[193]+1.677050983124842*rdxF[1]*f[65]*nu[173]+0.7499999999999999*rdxF[1]*f[49]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 18.44756081437327*rdxF[1]*f[38]*nu[193]+10.5*rdxF[1]*f[5]*nu[193]+10.06230589874905*rdxF[1]*f[66]*nu[173]+12.8086884574495*rdxF[1]*f[15]*nu[173]+7.685213074469699*rdxF[1]*f[38]*nu[163]+6.873863542433759*rdxF[1]*f[5]*nu[163]+5.7282196186948*rdxF[1]*f[15]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 1.14564392373896*rdxF[1]*f[54]*nu[193]+1.677050983124842*rdxF[1]*f[69]*nu[173]+0.7499999999999999*rdxF[1]*f[54]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[70] += 1.14564392373896*rdxF[1]*f[55]*nu[193]+1.677050983124842*rdxF[1]*f[70]*nu[173]+0.7499999999999999*rdxF[1]*f[55]*nu[163]; 
  out[71] += 18.44756081437327*rdxF[1]*f[43]*nu[193]+10.5*rdxF[1]*f[8]*nu[193]+10.06230589874905*rdxF[1]*f[71]*nu[173]+12.8086884574495*rdxF[1]*f[17]*nu[173]+7.685213074469699*rdxF[1]*f[43]*nu[163]+6.873863542433759*rdxF[1]*f[8]*nu[163]+5.7282196186948*rdxF[1]*f[17]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[72] += 18.44756081437327*rdxF[1]*f[44]*nu[193]+10.5*rdxF[1]*f[9]*nu[193]+10.06230589874905*rdxF[1]*f[72]*nu[173]+12.8086884574495*rdxF[1]*f[18]*nu[173]+7.685213074469699*rdxF[1]*f[44]*nu[163]+6.873863542433759*rdxF[1]*f[9]*nu[163]+5.7282196186948*rdxF[1]*f[18]*nu[160]; 
  out[73] += 1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 1.14564392373896*rdxF[1]*f[57]*nu[193]+1.677050983124842*rdxF[1]*f[74]*nu[173]+0.7499999999999999*rdxF[1]*f[57]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[75] += 1.14564392373896*rdxF[1]*f[58]*nu[193]+1.677050983124842*rdxF[1]*f[75]*nu[173]+0.7499999999999999*rdxF[1]*f[58]*nu[163]; 
  out[76] += 1.14564392373896*rdxF[1]*f[67]*nu[193]+1.677050983124842*rdxF[1]*f[76]*nu[173]+0.7499999999999999*rdxF[1]*f[67]*nu[163]+10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 1.14564392373896*rdxF[1]*f[68]*nu[193]+1.677050983124842*rdxF[1]*f[77]*nu[173]+0.7499999999999999*rdxF[1]*f[68]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 18.44756081437327*rdxF[1]*f[62]*nu[193]+10.5*rdxF[1]*f[16]*nu[193]+10.06230589874905*rdxF[1]*f[78]*nu[173]+12.8086884574495*rdxF[1]*f[35]*nu[173]+7.685213074469698*rdxF[1]*f[62]*nu[163]+6.873863542433758*rdxF[1]*f[16]*nu[163]+5.7282196186948*rdxF[1]*f[35]*nu[160]+1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 1.14564392373896*rdxF[1]*f[73]*nu[193]+1.677050983124842*rdxF[1]*f[79]*nu[173]+0.7499999999999999*rdxF[1]*f[73]*nu[163]+1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs4(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[3]*dx[3]); 

  out[4] += 2.5617376914899*rdxF[0]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[0]*nu[274]+1.677050983124842*rdxF[0]*f[4]*nu[254]+0.75*f[0]*rdxF[0]*nu[244]; 
  out[8] += 2.5617376914899*rdxF[0]*f[28]*nu[274]+1.14564392373896*rdxF[0]*f[1]*nu[274]+1.677050983124842*rdxF[0]*f[8]*nu[254]+0.75*rdxF[0]*f[1]*nu[244]; 
  out[9] += 2.5617376914899*rdxF[0]*f[29]*nu[274]+1.14564392373896*rdxF[0]*f[2]*nu[274]+1.677050983124842*rdxF[0]*f[9]*nu[254]+0.75*rdxF[0]*f[2]*nu[244]; 
  out[10] += 2.5617376914899*rdxF[0]*f[30]*nu[274]+1.14564392373896*rdxF[0]*f[3]*nu[274]+1.677050983124842*rdxF[0]*f[10]*nu[254]+0.75*rdxF[0]*f[3]*nu[244]; 
  out[14] += 6.708203932499369*rdxF[0]*f[34]*nu[274]+7.685213074469699*rdxF[0]*f[4]*nu[274]+5.031152949374527*rdxF[0]*f[14]*nu[254]+3.75*f[0]*rdxF[0]*nu[254]+3.354101966249685*rdxF[0]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[0]*nu[240]; 
  out[16] += 2.5617376914899*rdxF[0]*f[45]*nu[274]+1.14564392373896*rdxF[0]*f[5]*nu[274]+1.677050983124842*rdxF[0]*f[16]*nu[254]+0.75*rdxF[0]*f[5]*nu[244]; 
  out[17] += 2.5617376914899*rdxF[0]*f[46]*nu[274]+1.14564392373896*rdxF[0]*f[6]*nu[274]+1.677050983124842*rdxF[0]*f[17]*nu[254]+0.75*rdxF[0]*f[6]*nu[244]; 
  out[18] += 2.5617376914899*rdxF[0]*f[47]*nu[274]+1.14564392373896*rdxF[0]*f[7]*nu[274]+1.677050983124842*rdxF[0]*f[18]*nu[254]+0.75*rdxF[0]*f[7]*nu[244]; 
  out[25] += 1.14564392373896*rdxF[0]*f[11]*nu[274]+1.677050983124842*rdxF[0]*f[25]*nu[254]+0.75*rdxF[0]*f[11]*nu[244]; 
  out[26] += 1.14564392373896*rdxF[0]*f[12]*nu[274]+1.677050983124842*rdxF[0]*f[26]*nu[254]+0.75*rdxF[0]*f[12]*nu[244]; 
  out[27] += 1.14564392373896*rdxF[0]*f[13]*nu[274]+1.677050983124842*rdxF[0]*f[27]*nu[254]+0.75*rdxF[0]*f[13]*nu[244]; 
  out[28] += 6.708203932499369*rdxF[0]*f[57]*nu[274]+7.6852130744697*rdxF[0]*f[8]*nu[274]+5.031152949374527*rdxF[0]*f[28]*nu[254]+3.75*rdxF[0]*f[1]*nu[254]+3.354101966249684*rdxF[0]*f[8]*nu[244]+1.677050983124842*rdxF[0]*f[1]*nu[240]; 
  out[29] += 6.708203932499369*rdxF[0]*f[58]*nu[274]+7.6852130744697*rdxF[0]*f[9]*nu[274]+5.031152949374527*rdxF[0]*f[29]*nu[254]+3.75*rdxF[0]*f[2]*nu[254]+3.354101966249684*rdxF[0]*f[9]*nu[244]+1.677050983124842*rdxF[0]*f[2]*nu[240]; 
  out[30] += 6.708203932499369*rdxF[0]*f[59]*nu[274]+7.6852130744697*rdxF[0]*f[10]*nu[274]+5.031152949374527*rdxF[0]*f[30]*nu[254]+3.75*rdxF[0]*f[3]*nu[254]+3.354101966249684*rdxF[0]*f[10]*nu[244]+1.677050983124842*rdxF[0]*f[3]*nu[240]; 
  out[34] += 18.44756081437327*rdxF[0]*f[14]*nu[274]+10.5*f[0]*rdxF[0]*nu[274]+10.06230589874905*rdxF[0]*f[34]*nu[254]+12.8086884574495*rdxF[0]*f[4]*nu[254]+7.685213074469699*rdxF[0]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[0]*nu[244]+5.7282196186948*rdxF[0]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[0]*f[63]*nu[274]+1.14564392373896*rdxF[0]*f[15]*nu[274]+1.677050983124842*rdxF[0]*f[35]*nu[254]+0.75*rdxF[0]*f[15]*nu[244]; 
  out[39] += 1.14564392373896*rdxF[0]*f[19]*nu[274]+1.677050983124842*rdxF[0]*f[39]*nu[254]+0.75*rdxF[0]*f[19]*nu[244]; 
  out[40] += 1.14564392373896*rdxF[0]*f[20]*nu[274]+1.677050983124842*rdxF[0]*f[40]*nu[254]+0.75*rdxF[0]*f[20]*nu[244]; 
  out[41] += 1.14564392373896*rdxF[0]*f[21]*nu[274]+1.677050983124842*rdxF[0]*f[41]*nu[254]+0.75*rdxF[0]*f[21]*nu[244]; 
  out[42] += 1.14564392373896*rdxF[0]*f[22]*nu[274]+1.677050983124842*rdxF[0]*f[42]*nu[254]+0.75*rdxF[0]*f[22]*nu[244]; 
  out[43] += 1.14564392373896*rdxF[0]*f[23]*nu[274]+1.677050983124842*rdxF[0]*f[43]*nu[254]+0.75*rdxF[0]*f[23]*nu[244]; 
  out[44] += 1.14564392373896*rdxF[0]*f[24]*nu[274]+1.677050983124842*rdxF[0]*f[44]*nu[254]+0.75*rdxF[0]*f[24]*nu[244]; 
  out[45] += 6.708203932499369*rdxF[0]*f[73]*nu[274]+7.685213074469699*rdxF[0]*f[16]*nu[274]+5.031152949374527*rdxF[0]*f[45]*nu[254]+3.75*rdxF[0]*f[5]*nu[254]+3.354101966249685*rdxF[0]*f[16]*nu[244]+1.677050983124842*rdxF[0]*f[5]*nu[240]; 
  out[46] += 6.708203932499369*rdxF[0]*f[74]*nu[274]+7.685213074469699*rdxF[0]*f[17]*nu[274]+5.031152949374527*rdxF[0]*f[46]*nu[254]+3.75*rdxF[0]*f[6]*nu[254]+3.354101966249685*rdxF[0]*f[17]*nu[244]+1.677050983124842*rdxF[0]*f[6]*nu[240]; 
  out[47] += 6.708203932499369*rdxF[0]*f[75]*nu[274]+7.685213074469699*rdxF[0]*f[18]*nu[274]+5.031152949374527*rdxF[0]*f[47]*nu[254]+3.75*rdxF[0]*f[7]*nu[254]+3.354101966249685*rdxF[0]*f[18]*nu[244]+1.677050983124842*rdxF[0]*f[7]*nu[240]; 
  out[54] += 1.14564392373896*rdxF[0]*f[31]*nu[274]+1.677050983124842*rdxF[0]*f[54]*nu[254]+0.7499999999999999*rdxF[0]*f[31]*nu[244]; 
  out[55] += 1.14564392373896*rdxF[0]*f[32]*nu[274]+1.677050983124842*rdxF[0]*f[55]*nu[254]+0.7499999999999999*rdxF[0]*f[32]*nu[244]; 
  out[56] += 1.14564392373896*rdxF[0]*f[33]*nu[274]+1.677050983124842*rdxF[0]*f[56]*nu[254]+0.7499999999999999*rdxF[0]*f[33]*nu[244]; 
  out[57] += 18.44756081437327*rdxF[0]*f[28]*nu[274]+10.5*rdxF[0]*f[1]*nu[274]+10.06230589874905*rdxF[0]*f[57]*nu[254]+12.8086884574495*rdxF[0]*f[8]*nu[254]+7.685213074469698*rdxF[0]*f[28]*nu[244]+6.873863542433758*rdxF[0]*f[1]*nu[244]+5.7282196186948*rdxF[0]*f[8]*nu[240]; 
  out[58] += 18.44756081437327*rdxF[0]*f[29]*nu[274]+10.5*rdxF[0]*f[2]*nu[274]+10.06230589874905*rdxF[0]*f[58]*nu[254]+12.8086884574495*rdxF[0]*f[9]*nu[254]+7.685213074469698*rdxF[0]*f[29]*nu[244]+6.873863542433758*rdxF[0]*f[2]*nu[244]+5.7282196186948*rdxF[0]*f[9]*nu[240]; 
  out[59] += 18.44756081437327*rdxF[0]*f[30]*nu[274]+10.5*rdxF[0]*f[3]*nu[274]+10.06230589874905*rdxF[0]*f[59]*nu[254]+12.8086884574495*rdxF[0]*f[10]*nu[254]+7.685213074469698*rdxF[0]*f[30]*nu[244]+6.873863542433758*rdxF[0]*f[3]*nu[244]+5.7282196186948*rdxF[0]*f[10]*nu[240]; 
  out[60] += 1.14564392373896*rdxF[0]*f[36]*nu[274]+1.677050983124842*rdxF[0]*f[60]*nu[254]+0.75*rdxF[0]*f[36]*nu[244]; 
  out[61] += 1.14564392373896*rdxF[0]*f[37]*nu[274]+1.677050983124842*rdxF[0]*f[61]*nu[254]+0.75*rdxF[0]*f[37]*nu[244]; 
  out[62] += 1.14564392373896*rdxF[0]*f[38]*nu[274]+1.677050983124842*rdxF[0]*f[62]*nu[254]+0.75*rdxF[0]*f[38]*nu[244]; 
  out[63] += 6.708203932499369*rdxF[0]*f[79]*nu[274]+7.6852130744697*rdxF[0]*f[35]*nu[274]+5.031152949374527*rdxF[0]*f[63]*nu[254]+3.75*rdxF[0]*f[15]*nu[254]+3.354101966249684*rdxF[0]*f[35]*nu[244]+1.677050983124842*rdxF[0]*f[15]*nu[240]; 
  out[67] += 1.14564392373896*rdxF[0]*f[48]*nu[274]+1.677050983124842*rdxF[0]*f[67]*nu[254]+0.7499999999999999*rdxF[0]*f[48]*nu[244]; 
  out[68] += 1.14564392373896*rdxF[0]*f[49]*nu[274]+1.677050983124842*rdxF[0]*f[68]*nu[254]+0.7499999999999999*rdxF[0]*f[49]*nu[244]; 
  out[69] += 1.14564392373896*rdxF[0]*f[50]*nu[274]+1.677050983124842*rdxF[0]*f[69]*nu[254]+0.7499999999999999*rdxF[0]*f[50]*nu[244]; 
  out[70] += 1.14564392373896*rdxF[0]*f[51]*nu[274]+1.677050983124842*rdxF[0]*f[70]*nu[254]+0.7499999999999999*rdxF[0]*f[51]*nu[244]; 
  out[71] += 1.14564392373896*rdxF[0]*f[52]*nu[274]+1.677050983124842*rdxF[0]*f[71]*nu[254]+0.7499999999999999*rdxF[0]*f[52]*nu[244]; 
  out[72] += 1.14564392373896*rdxF[0]*f[53]*nu[274]+1.677050983124842*rdxF[0]*f[72]*nu[254]+0.7499999999999999*rdxF[0]*f[53]*nu[244]; 
  out[73] += 18.44756081437327*rdxF[0]*f[45]*nu[274]+10.5*rdxF[0]*f[5]*nu[274]+10.06230589874905*rdxF[0]*f[73]*nu[254]+12.8086884574495*rdxF[0]*f[16]*nu[254]+7.685213074469699*rdxF[0]*f[45]*nu[244]+6.873863542433759*rdxF[0]*f[5]*nu[244]+5.7282196186948*rdxF[0]*f[16]*nu[240]; 
  out[74] += 18.44756081437327*rdxF[0]*f[46]*nu[274]+10.5*rdxF[0]*f[6]*nu[274]+10.06230589874905*rdxF[0]*f[74]*nu[254]+12.8086884574495*rdxF[0]*f[17]*nu[254]+7.685213074469699*rdxF[0]*f[46]*nu[244]+6.873863542433759*rdxF[0]*f[6]*nu[244]+5.7282196186948*rdxF[0]*f[17]*nu[240]; 
  out[75] += 18.44756081437327*rdxF[0]*f[47]*nu[274]+10.5*rdxF[0]*f[7]*nu[274]+10.06230589874905*rdxF[0]*f[75]*nu[254]+12.8086884574495*rdxF[0]*f[18]*nu[254]+7.685213074469699*rdxF[0]*f[47]*nu[244]+6.873863542433759*rdxF[0]*f[7]*nu[244]+5.7282196186948*rdxF[0]*f[18]*nu[240]; 
  out[76] += 1.14564392373896*rdxF[0]*f[64]*nu[274]+1.677050983124842*rdxF[0]*f[76]*nu[254]+0.7499999999999999*rdxF[0]*f[64]*nu[244]; 
  out[77] += 1.14564392373896*rdxF[0]*f[65]*nu[274]+1.677050983124842*rdxF[0]*f[77]*nu[254]+0.7499999999999999*rdxF[0]*f[65]*nu[244]; 
  out[78] += 1.14564392373896*rdxF[0]*f[66]*nu[274]+1.677050983124842*rdxF[0]*f[78]*nu[254]+0.7499999999999999*rdxF[0]*f[66]*nu[244]; 
  out[79] += 18.44756081437327*rdxF[0]*f[63]*nu[274]+10.5*rdxF[0]*f[15]*nu[274]+10.06230589874905*rdxF[0]*f[79]*nu[254]+12.8086884574495*rdxF[0]*f[35]*nu[254]+7.685213074469698*rdxF[0]*f[63]*nu[244]+6.873863542433758*rdxF[0]*f[15]*nu[244]+5.7282196186948*rdxF[0]*f[35]*nu[240]; 

  return (rdxF[0]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs14(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[2]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 
  rdxF[1] = 4.0/(dx[3]*dx[3]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[4] += 2.5617376914899*rdxF[1]*f[14]*nu[274]+1.14564392373896*f[0]*rdxF[1]*nu[274]+1.677050983124842*rdxF[1]*f[4]*nu[254]+0.75*f[0]*rdxF[1]*nu[244]; 
  out[5] += 2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[8] += 2.5617376914899*rdxF[1]*f[28]*nu[274]+1.14564392373896*f[1]*rdxF[1]*nu[274]+1.677050983124842*rdxF[1]*f[8]*nu[254]+0.75*f[1]*rdxF[1]*nu[244]+2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[9] += 2.5617376914899*rdxF[1]*f[29]*nu[274]+1.14564392373896*rdxF[1]*f[2]*nu[274]+1.677050983124842*rdxF[1]*f[9]*nu[254]+0.75*rdxF[1]*f[2]*nu[244]; 
  out[10] += 2.5617376914899*rdxF[1]*f[30]*nu[274]+1.14564392373896*rdxF[1]*f[3]*nu[274]+1.677050983124842*rdxF[1]*f[10]*nu[254]+0.75*rdxF[1]*f[3]*nu[244]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[14] += 6.708203932499369*rdxF[1]*f[34]*nu[274]+7.685213074469699*rdxF[1]*f[4]*nu[274]+5.031152949374527*rdxF[1]*f[14]*nu[254]+3.75*f[0]*rdxF[1]*nu[254]+3.354101966249685*rdxF[1]*f[4]*nu[244]+1.677050983124842*f[0]*rdxF[1]*nu[240]; 
  out[15] += 2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[1]*f[45]*nu[274]+1.14564392373896*rdxF[1]*f[5]*nu[274]+1.677050983124842*rdxF[1]*f[16]*nu[254]+0.75*rdxF[1]*f[5]*nu[244]+2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[1]*f[46]*nu[274]+1.14564392373896*rdxF[1]*f[6]*nu[274]+1.677050983124842*rdxF[1]*f[17]*nu[254]+0.75*rdxF[1]*f[6]*nu[244]+2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[18] += 2.5617376914899*rdxF[1]*f[47]*nu[274]+1.14564392373896*rdxF[1]*f[7]*nu[274]+1.677050983124842*rdxF[1]*f[18]*nu[254]+0.75*rdxF[1]*f[7]*nu[244]; 
  out[19] += 6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[23] += 1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[25] += 1.14564392373896*rdxF[1]*f[11]*nu[274]+1.677050983124842*rdxF[1]*f[25]*nu[254]+0.75*rdxF[1]*f[11]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[26] += 1.14564392373896*rdxF[1]*f[12]*nu[274]+1.677050983124842*rdxF[1]*f[26]*nu[254]+0.75*rdxF[1]*f[12]*nu[244]; 
  out[27] += 1.14564392373896*rdxF[1]*f[13]*nu[274]+1.677050983124842*rdxF[1]*f[27]*nu[254]+0.75*rdxF[1]*f[13]*nu[244]; 
  out[28] += 6.708203932499369*rdxF[1]*f[57]*nu[274]+7.6852130744697*rdxF[1]*f[8]*nu[274]+5.031152949374527*rdxF[1]*f[28]*nu[254]+3.75*f[1]*rdxF[1]*nu[254]+3.354101966249684*rdxF[1]*f[8]*nu[244]+1.677050983124842*f[1]*rdxF[1]*nu[240]+1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[29] += 6.708203932499369*rdxF[1]*f[58]*nu[274]+7.6852130744697*rdxF[1]*f[9]*nu[274]+5.031152949374527*rdxF[1]*f[29]*nu[254]+3.75*rdxF[1]*f[2]*nu[254]+3.354101966249684*rdxF[1]*f[9]*nu[244]+1.677050983124842*rdxF[1]*f[2]*nu[240]; 
  out[30] += 6.708203932499369*rdxF[1]*f[59]*nu[274]+7.6852130744697*rdxF[1]*f[10]*nu[274]+5.031152949374527*rdxF[1]*f[30]*nu[254]+3.75*rdxF[1]*f[3]*nu[254]+3.354101966249684*rdxF[1]*f[10]*nu[244]+1.677050983124842*rdxF[1]*f[3]*nu[240]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[34] += 18.44756081437327*rdxF[1]*f[14]*nu[274]+10.5*f[0]*rdxF[1]*nu[274]+10.06230589874905*rdxF[1]*f[34]*nu[254]+12.8086884574495*rdxF[1]*f[4]*nu[254]+7.685213074469699*rdxF[1]*f[14]*nu[244]+6.873863542433759*f[0]*rdxF[1]*nu[244]+5.7282196186948*rdxF[1]*f[4]*nu[240]; 
  out[35] += 2.5617376914899*rdxF[1]*f[63]*nu[274]+1.14564392373896*rdxF[1]*f[15]*nu[274]+1.677050983124842*rdxF[1]*f[35]*nu[254]+0.75*rdxF[1]*f[15]*nu[244]+2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 1.14564392373896*rdxF[1]*f[19]*nu[274]+1.677050983124842*rdxF[1]*f[39]*nu[254]+0.75*rdxF[1]*f[19]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 1.14564392373896*rdxF[1]*f[20]*nu[274]+1.677050983124842*rdxF[1]*f[40]*nu[254]+0.75*rdxF[1]*f[20]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 1.14564392373896*rdxF[1]*f[21]*nu[274]+1.677050983124842*rdxF[1]*f[41]*nu[254]+0.75*rdxF[1]*f[21]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[42] += 1.14564392373896*rdxF[1]*f[22]*nu[274]+1.677050983124842*rdxF[1]*f[42]*nu[254]+0.75*rdxF[1]*f[22]*nu[244]; 
  out[43] += 1.14564392373896*rdxF[1]*f[23]*nu[274]+1.677050983124842*rdxF[1]*f[43]*nu[254]+0.75*rdxF[1]*f[23]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[44] += 1.14564392373896*rdxF[1]*f[24]*nu[274]+1.677050983124842*rdxF[1]*f[44]*nu[254]+0.75*rdxF[1]*f[24]*nu[244]; 
  out[45] += 6.708203932499369*rdxF[1]*f[73]*nu[274]+7.685213074469699*rdxF[1]*f[16]*nu[274]+5.031152949374527*rdxF[1]*f[45]*nu[254]+3.75*rdxF[1]*f[5]*nu[254]+3.354101966249685*rdxF[1]*f[16]*nu[244]+1.677050983124842*rdxF[1]*f[5]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 6.708203932499369*rdxF[1]*f[74]*nu[274]+7.685213074469699*rdxF[1]*f[17]*nu[274]+5.031152949374527*rdxF[1]*f[46]*nu[254]+3.75*rdxF[1]*f[6]*nu[254]+3.354101966249685*rdxF[1]*f[17]*nu[244]+1.677050983124842*rdxF[1]*f[6]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[47] += 6.708203932499369*rdxF[1]*f[75]*nu[274]+7.685213074469699*rdxF[1]*f[18]*nu[274]+5.031152949374527*rdxF[1]*f[47]*nu[254]+3.75*rdxF[1]*f[7]*nu[254]+3.354101966249685*rdxF[1]*f[18]*nu[244]+1.677050983124842*rdxF[1]*f[7]*nu[240]; 
  out[48] += 10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[52] += 1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[54] += 1.14564392373896*rdxF[1]*f[31]*nu[274]+1.677050983124842*rdxF[1]*f[54]*nu[254]+0.7499999999999999*rdxF[1]*f[31]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[55] += 1.14564392373896*rdxF[1]*f[32]*nu[274]+1.677050983124842*rdxF[1]*f[55]*nu[254]+0.7499999999999999*rdxF[1]*f[32]*nu[244]; 
  out[56] += 1.14564392373896*rdxF[1]*f[33]*nu[274]+1.677050983124842*rdxF[1]*f[56]*nu[254]+0.7499999999999999*rdxF[1]*f[33]*nu[244]; 
  out[57] += 18.44756081437327*rdxF[1]*f[28]*nu[274]+10.5*f[1]*rdxF[1]*nu[274]+10.06230589874905*rdxF[1]*f[57]*nu[254]+12.8086884574495*rdxF[1]*f[8]*nu[254]+7.685213074469698*rdxF[1]*f[28]*nu[244]+6.873863542433758*f[1]*rdxF[1]*nu[244]+5.7282196186948*rdxF[1]*f[8]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[58] += 18.44756081437327*rdxF[1]*f[29]*nu[274]+10.5*rdxF[1]*f[2]*nu[274]+10.06230589874905*rdxF[1]*f[58]*nu[254]+12.8086884574495*rdxF[1]*f[9]*nu[254]+7.685213074469698*rdxF[1]*f[29]*nu[244]+6.873863542433758*rdxF[1]*f[2]*nu[244]+5.7282196186948*rdxF[1]*f[9]*nu[240]; 
  out[59] += 18.44756081437327*rdxF[1]*f[30]*nu[274]+10.5*rdxF[1]*f[3]*nu[274]+10.06230589874905*rdxF[1]*f[59]*nu[254]+12.8086884574495*rdxF[1]*f[10]*nu[254]+7.685213074469698*rdxF[1]*f[30]*nu[244]+6.873863542433758*rdxF[1]*f[3]*nu[244]+5.7282196186948*rdxF[1]*f[10]*nu[240]; 
  out[60] += 1.14564392373896*rdxF[1]*f[36]*nu[274]+1.677050983124842*rdxF[1]*f[60]*nu[254]+0.75*rdxF[1]*f[36]*nu[244]+6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 1.14564392373896*rdxF[1]*f[37]*nu[274]+1.677050983124842*rdxF[1]*f[61]*nu[254]+0.75*rdxF[1]*f[37]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 1.14564392373896*rdxF[1]*f[38]*nu[274]+1.677050983124842*rdxF[1]*f[62]*nu[254]+0.75*rdxF[1]*f[38]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 6.708203932499369*rdxF[1]*f[79]*nu[274]+7.6852130744697*rdxF[1]*f[35]*nu[274]+5.031152949374527*rdxF[1]*f[63]*nu[254]+3.75*rdxF[1]*f[15]*nu[254]+3.354101966249684*rdxF[1]*f[35]*nu[244]+1.677050983124842*rdxF[1]*f[15]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 1.14564392373896*rdxF[1]*f[48]*nu[274]+1.677050983124842*rdxF[1]*f[67]*nu[254]+0.7499999999999999*rdxF[1]*f[48]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 1.14564392373896*rdxF[1]*f[49]*nu[274]+1.677050983124842*rdxF[1]*f[68]*nu[254]+0.7499999999999999*rdxF[1]*f[49]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 1.14564392373896*rdxF[1]*f[50]*nu[274]+1.677050983124842*rdxF[1]*f[69]*nu[254]+0.7499999999999999*rdxF[1]*f[50]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[70] += 1.14564392373896*rdxF[1]*f[51]*nu[274]+1.677050983124842*rdxF[1]*f[70]*nu[254]+0.7499999999999999*rdxF[1]*f[51]*nu[244]; 
  out[71] += 1.14564392373896*rdxF[1]*f[52]*nu[274]+1.677050983124842*rdxF[1]*f[71]*nu[254]+0.7499999999999999*rdxF[1]*f[52]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[72] += 1.14564392373896*rdxF[1]*f[53]*nu[274]+1.677050983124842*rdxF[1]*f[72]*nu[254]+0.7499999999999999*rdxF[1]*f[53]*nu[244]; 
  out[73] += 18.44756081437327*rdxF[1]*f[45]*nu[274]+10.5*rdxF[1]*f[5]*nu[274]+10.06230589874905*rdxF[1]*f[73]*nu[254]+12.8086884574495*rdxF[1]*f[16]*nu[254]+7.685213074469699*rdxF[1]*f[45]*nu[244]+6.873863542433759*rdxF[1]*f[5]*nu[244]+5.7282196186948*rdxF[1]*f[16]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 18.44756081437327*rdxF[1]*f[46]*nu[274]+10.5*rdxF[1]*f[6]*nu[274]+10.06230589874905*rdxF[1]*f[74]*nu[254]+12.8086884574495*rdxF[1]*f[17]*nu[254]+7.685213074469699*rdxF[1]*f[46]*nu[244]+6.873863542433759*rdxF[1]*f[6]*nu[244]+5.7282196186948*rdxF[1]*f[17]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[75] += 18.44756081437327*rdxF[1]*f[47]*nu[274]+10.5*rdxF[1]*f[7]*nu[274]+10.06230589874905*rdxF[1]*f[75]*nu[254]+12.8086884574495*rdxF[1]*f[18]*nu[254]+7.685213074469699*rdxF[1]*f[47]*nu[244]+6.873863542433759*rdxF[1]*f[7]*nu[244]+5.7282196186948*rdxF[1]*f[18]*nu[240]; 
  out[76] += 1.14564392373896*rdxF[1]*f[64]*nu[274]+1.677050983124842*rdxF[1]*f[76]*nu[254]+0.7499999999999999*rdxF[1]*f[64]*nu[244]+10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 1.14564392373896*rdxF[1]*f[65]*nu[274]+1.677050983124842*rdxF[1]*f[77]*nu[254]+0.7499999999999999*rdxF[1]*f[65]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 1.14564392373896*rdxF[1]*f[66]*nu[274]+1.677050983124842*rdxF[1]*f[78]*nu[254]+0.7499999999999999*rdxF[1]*f[66]*nu[244]+1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 18.44756081437327*rdxF[1]*f[63]*nu[274]+10.5*rdxF[1]*f[15]*nu[274]+10.06230589874905*rdxF[1]*f[79]*nu[254]+12.8086884574495*rdxF[1]*f[35]*nu[254]+7.685213074469698*rdxF[1]*f[63]*nu[244]+6.873863542433758*rdxF[1]*f[15]*nu[244]+5.7282196186948*rdxF[1]*f[35]*nu[240]+1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+rdxF[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*1.142857142857143;

} 
double ConstDiffusionVarCoeffVol4xSerP3_diffDirs1(const double *w, const double *dx, const double *nu, const double *f, double *out) 
{ 
  // w[4]:   Cell-center coordinates.
  // dx[4]:  Cell spacing.
  // nu[320]:  diffusion coefficient.
  // f:      Input distribution function.
  // out:    Incremented output 
  double rdxF[1]; 
  rdxF[0] = 4.0/(dx[0]*dx[0]); 

  out[1] += 2.5617376914899*rdxF[0]*f[11]*nu[31]+1.14564392373896*f[0]*rdxF[0]*nu[31]+1.677050983124842*rdxF[0]*f[1]*nu[11]+0.75*f[0]*rdxF[0]*nu[1]; 
  out[5] += 2.5617376914899*rdxF[0]*f[19]*nu[31]+1.14564392373896*rdxF[0]*f[2]*nu[31]+1.677050983124842*rdxF[0]*f[5]*nu[11]+0.75*rdxF[0]*nu[1]*f[2]; 
  out[6] += 2.5617376914899*rdxF[0]*f[21]*nu[31]+1.14564392373896*rdxF[0]*f[3]*nu[31]+1.677050983124842*rdxF[0]*f[6]*nu[11]+0.75*rdxF[0]*nu[1]*f[3]; 
  out[8] += 2.5617376914899*rdxF[0]*f[25]*nu[31]+1.14564392373896*rdxF[0]*f[4]*nu[31]+1.677050983124842*rdxF[0]*f[8]*nu[11]+0.75*rdxF[0]*nu[1]*f[4]; 
  out[11] += 6.708203932499369*rdxF[0]*f[31]*nu[31]+7.685213074469699*rdxF[0]*f[1]*nu[31]+5.031152949374527*rdxF[0]*f[11]*nu[11]+3.75*f[0]*rdxF[0]*nu[11]+3.354101966249685*rdxF[0]*f[1]*nu[1]+1.677050983124842*f[0]*nu[0]*rdxF[0]; 
  out[15] += 2.5617376914899*rdxF[0]*nu[31]*f[36]+1.14564392373896*rdxF[0]*f[7]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[15]+0.75*rdxF[0]*nu[1]*f[7]; 
  out[16] += 2.5617376914899*rdxF[0]*nu[31]*f[39]+1.14564392373896*rdxF[0]*f[9]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[16]+0.75*rdxF[0]*nu[1]*f[9]; 
  out[17] += 2.5617376914899*rdxF[0]*nu[31]*f[41]+1.14564392373896*rdxF[0]*f[10]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[17]+0.75*rdxF[0]*nu[1]*f[10]; 
  out[19] += 6.708203932499369*rdxF[0]*nu[31]*f[48]+7.6852130744697*rdxF[0]*f[5]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[19]+3.75*rdxF[0]*f[2]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[5]+1.677050983124842*nu[0]*rdxF[0]*f[2]; 
  out[20] += 1.14564392373896*rdxF[0]*f[12]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[20]+0.75*rdxF[0]*nu[1]*f[12]; 
  out[21] += 6.708203932499369*rdxF[0]*nu[31]*f[50]+7.6852130744697*rdxF[0]*f[6]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[21]+3.75*rdxF[0]*f[3]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[6]+1.677050983124842*nu[0]*rdxF[0]*f[3]; 
  out[23] += 1.14564392373896*rdxF[0]*f[13]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[23]+0.75*rdxF[0]*nu[1]*f[13]; 
  out[25] += 6.708203932499369*rdxF[0]*nu[31]*f[54]+7.6852130744697*rdxF[0]*f[8]*nu[31]+5.031152949374527*rdxF[0]*nu[11]*f[25]+3.75*rdxF[0]*f[4]*nu[11]+3.354101966249684*rdxF[0]*nu[1]*f[8]+1.677050983124842*nu[0]*rdxF[0]*f[4]; 
  out[28] += 1.14564392373896*rdxF[0]*f[14]*nu[31]+1.677050983124842*rdxF[0]*nu[11]*f[28]+0.75*rdxF[0]*nu[1]*f[14]; 
  out[31] += 18.44756081437327*rdxF[0]*f[11]*nu[31]+10.5*f[0]*rdxF[0]*nu[31]+10.06230589874905*rdxF[0]*nu[11]*f[31]+12.8086884574495*rdxF[0]*f[1]*nu[11]+7.685213074469699*rdxF[0]*nu[1]*f[11]+6.873863542433759*f[0]*rdxF[0]*nu[1]+5.7282196186948*nu[0]*rdxF[0]*f[1]; 
  out[35] += 2.5617376914899*rdxF[0]*nu[31]*f[60]+1.677050983124842*rdxF[0]*nu[11]*f[35]+1.14564392373896*rdxF[0]*f[18]*nu[31]+0.75*rdxF[0]*nu[1]*f[18]; 
  out[36] += 6.708203932499369*rdxF[0]*nu[31]*f[64]+5.031152949374527*rdxF[0]*nu[11]*f[36]+7.685213074469699*rdxF[0]*f[15]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[15]+3.75*rdxF[0]*f[7]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[7]; 
  out[37] += 1.677050983124842*rdxF[0]*nu[11]*f[37]+1.14564392373896*rdxF[0]*f[22]*nu[31]+0.75*rdxF[0]*nu[1]*f[22]; 
  out[38] += 1.677050983124842*rdxF[0]*nu[11]*f[38]+1.14564392373896*rdxF[0]*f[24]*nu[31]+0.75*rdxF[0]*nu[1]*f[24]; 
  out[39] += 6.708203932499369*rdxF[0]*nu[31]*f[67]+5.031152949374527*rdxF[0]*nu[11]*f[39]+7.685213074469699*rdxF[0]*f[16]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[16]+3.75*rdxF[0]*f[9]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[9]; 
  out[40] += 1.677050983124842*rdxF[0]*nu[11]*f[40]+1.14564392373896*rdxF[0]*f[26]*nu[31]+0.75*rdxF[0]*nu[1]*f[26]; 
  out[41] += 6.708203932499369*rdxF[0]*nu[31]*f[69]+5.031152949374527*rdxF[0]*nu[11]*f[41]+7.685213074469699*rdxF[0]*f[17]*nu[31]+3.354101966249685*rdxF[0]*nu[1]*f[17]+3.75*rdxF[0]*f[10]*nu[11]+1.677050983124842*nu[0]*rdxF[0]*f[10]; 
  out[43] += 1.677050983124842*rdxF[0]*nu[11]*f[43]+1.14564392373896*rdxF[0]*f[27]*nu[31]+0.75*rdxF[0]*nu[1]*f[27]; 
  out[45] += 1.677050983124842*rdxF[0]*nu[11]*f[45]+1.14564392373896*rdxF[0]*f[29]*nu[31]+0.75*rdxF[0]*nu[1]*f[29]; 
  out[46] += 1.677050983124842*rdxF[0]*nu[11]*f[46]+1.14564392373896*rdxF[0]*f[30]*nu[31]+0.75*rdxF[0]*nu[1]*f[30]; 
  out[48] += 10.06230589874905*rdxF[0]*nu[11]*f[48]+18.44756081437327*rdxF[0]*f[19]*nu[31]+10.5*rdxF[0]*f[2]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[19]+12.8086884574495*rdxF[0]*f[5]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[5]+6.873863542433758*rdxF[0]*nu[1]*f[2]; 
  out[49] += 1.677050983124842*rdxF[0]*nu[11]*f[49]+1.14564392373896*rdxF[0]*nu[31]*f[32]+0.7499999999999999*rdxF[0]*nu[1]*f[32]; 
  out[50] += 10.06230589874905*rdxF[0]*nu[11]*f[50]+18.44756081437327*rdxF[0]*f[21]*nu[31]+10.5*rdxF[0]*f[3]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[21]+12.8086884574495*rdxF[0]*f[6]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[6]+6.873863542433758*rdxF[0]*nu[1]*f[3]; 
  out[52] += 1.677050983124842*rdxF[0]*nu[11]*f[52]+1.14564392373896*rdxF[0]*nu[31]*f[33]+0.7499999999999999*rdxF[0]*nu[1]*f[33]; 
  out[54] += 10.06230589874905*rdxF[0]*nu[11]*f[54]+18.44756081437327*rdxF[0]*f[25]*nu[31]+10.5*rdxF[0]*f[4]*nu[31]+7.685213074469698*rdxF[0]*nu[1]*f[25]+12.8086884574495*rdxF[0]*f[8]*nu[11]+5.7282196186948*nu[0]*rdxF[0]*f[8]+6.873863542433758*rdxF[0]*nu[1]*f[4]; 
  out[57] += 1.677050983124842*rdxF[0]*nu[11]*f[57]+1.14564392373896*rdxF[0]*nu[31]*f[34]+0.7499999999999999*rdxF[0]*nu[1]*f[34]; 
  out[60] += 6.708203932499369*rdxF[0]*nu[31]*f[76]+5.031152949374527*rdxF[0]*nu[11]*f[60]+7.6852130744697*rdxF[0]*nu[31]*f[35]+3.354101966249684*rdxF[0]*nu[1]*f[35]+3.75*rdxF[0]*nu[11]*f[18]+1.677050983124842*nu[0]*rdxF[0]*f[18]; 
  out[61] += 1.677050983124842*rdxF[0]*nu[11]*f[61]+1.14564392373896*rdxF[0]*nu[31]*f[42]+0.75*rdxF[0]*nu[1]*f[42]; 
  out[62] += 1.677050983124842*rdxF[0]*nu[11]*f[62]+1.14564392373896*rdxF[0]*nu[31]*f[44]+0.75*rdxF[0]*nu[1]*f[44]; 
  out[63] += 1.677050983124842*rdxF[0]*nu[11]*f[63]+1.14564392373896*rdxF[0]*nu[31]*f[47]+0.75*rdxF[0]*nu[1]*f[47]; 
  out[64] += 10.06230589874905*rdxF[0]*nu[11]*f[64]+18.44756081437327*rdxF[0]*nu[31]*f[36]+7.685213074469699*rdxF[0]*nu[1]*f[36]+10.5*rdxF[0]*f[7]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[15]+5.7282196186948*nu[0]*rdxF[0]*f[15]+6.873863542433759*rdxF[0]*nu[1]*f[7]; 
  out[65] += 1.677050983124842*rdxF[0]*nu[11]*f[65]+1.14564392373896*rdxF[0]*nu[31]*f[51]+0.7499999999999999*rdxF[0]*nu[1]*f[51]; 
  out[66] += 1.677050983124842*rdxF[0]*nu[11]*f[66]+1.14564392373896*rdxF[0]*nu[31]*f[53]+0.7499999999999999*rdxF[0]*nu[1]*f[53]; 
  out[67] += 10.06230589874905*rdxF[0]*nu[11]*f[67]+18.44756081437327*rdxF[0]*nu[31]*f[39]+7.685213074469699*rdxF[0]*nu[1]*f[39]+10.5*rdxF[0]*f[9]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[16]+5.7282196186948*nu[0]*rdxF[0]*f[16]+6.873863542433759*rdxF[0]*nu[1]*f[9]; 
  out[68] += 1.677050983124842*rdxF[0]*nu[11]*f[68]+1.14564392373896*rdxF[0]*nu[31]*f[55]+0.7499999999999999*rdxF[0]*nu[1]*f[55]; 
  out[69] += 10.06230589874905*rdxF[0]*nu[11]*f[69]+18.44756081437327*rdxF[0]*nu[31]*f[41]+7.685213074469699*rdxF[0]*nu[1]*f[41]+10.5*rdxF[0]*f[10]*nu[31]+12.8086884574495*rdxF[0]*nu[11]*f[17]+5.7282196186948*nu[0]*rdxF[0]*f[17]+6.873863542433759*rdxF[0]*nu[1]*f[10]; 
  out[71] += 1.677050983124842*rdxF[0]*nu[11]*f[71]+1.14564392373896*rdxF[0]*nu[31]*f[56]+0.7499999999999999*rdxF[0]*nu[1]*f[56]; 
  out[73] += 1.677050983124842*rdxF[0]*nu[11]*f[73]+1.14564392373896*rdxF[0]*nu[31]*f[58]+0.7499999999999999*rdxF[0]*nu[1]*f[58]; 
  out[74] += 1.677050983124842*rdxF[0]*nu[11]*f[74]+1.14564392373896*rdxF[0]*nu[31]*f[59]+0.7499999999999999*rdxF[0]*nu[1]*f[59]; 
  out[76] += 10.06230589874905*rdxF[0]*nu[11]*f[76]+18.44756081437327*rdxF[0]*nu[31]*f[60]+7.685213074469698*rdxF[0]*nu[1]*f[60]+12.8086884574495*rdxF[0]*nu[11]*f[35]+5.7282196186948*nu[0]*rdxF[0]*f[35]+10.5*rdxF[0]*f[18]*nu[31]+6.873863542433758*rdxF[0]*nu[1]*f[18]; 
  out[77] += 1.677050983124842*rdxF[0]*nu[11]*f[77]+1.14564392373896*rdxF[0]*nu[31]*f[70]+0.7499999999999999*rdxF[0]*nu[1]*f[70]; 
  out[78] += 1.677050983124842*rdxF[0]*nu[11]*f[78]+1.14564392373896*rdxF[0]*nu[31]*f[72]+0.7499999999999999*rdxF[0]*nu[1]*f[72]; 
  out[79] += 1.677050983124842*rdxF[0]*nu[11]*f[79]+1.14564392373896*rdxF[0]*nu[31]*f[75]+0.7499999999999999*rdxF[0]*nu[1]*f[75]; 

  return (rdxF[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*1.142857142857143;

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[3]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstDiffusionCFLfreqMin4xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[3]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion4CFLfreqMin4xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs2(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs23(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs1234(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[4]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[3] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[3]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs123(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs24(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs124(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs12(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]*Lx[1]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[94])-0.2795084971874737*nu[93]-0.2795084971874737*nu[92]-0.2795084971874737*nu[91]+0.25*nu[80]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs3(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs34(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs134(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[3]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 
  kxSq[2] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160])+kxSq[2]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs13(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]*Lx[2]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[174])-0.2795084971874737*nu[173]-0.2795084971874737*nu[172]-0.2795084971874737*nu[171]+0.25*nu[160]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs4(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs14(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[2]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 
  kxSq[1] = 39.47841760435743/(Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]*Lx[3]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0])+kxSq[1]*((-0.2795084971874737*nu[254])-0.2795084971874737*nu[253]-0.2795084971874737*nu[252]-0.2795084971874737*nu[251]+0.25*nu[240]))*2.285714285714286);

} 
void ConstHyperDiffusion6CFLfreqMin4xSerP3_diffDirs1(const double *Lx, const double *nu, double *cflFreq) 
{ 
  // Lx[vdim]:  domain length.
  // nu[320]:  diffusion coefficient (collisionality).
  double kxSq[1]; 
  kxSq[0] = 39.47841760435743/(Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]*Lx[0]); 

  cflFreq[0] = fmax(cflFreq[0],(kxSq[0]*((-0.2795084971874737*nu[14])-0.2795084971874737*nu[13]-0.2795084971874737*nu[12]-0.2795084971874737*nu[11]+0.25*nu[0]))*2.285714285714286);

} 
