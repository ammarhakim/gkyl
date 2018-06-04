#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 


return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
double ConstDiffusionVol3xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[17] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[18] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[19] += 6.708203932499369*rdxSq4nu[2]*f[4]; 

return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
double ConstDiffusionVol3xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[17] += 22.9128784747792*rdxSq4nu[0]*f[1]; 
  out[18] += 22.9128784747792*rdxSq4nu[1]*f[2]; 
  out[19] += 22.9128784747792*rdxSq4nu[2]*f[3]; 
  out[20] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[21] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[22] += 6.708203932499369*rdxSq4nu[2]*f[4]; 
  out[23] += 22.9128784747792*rdxSq4nu[0]*f[4]; 
  out[24] += 22.9128784747792*rdxSq4nu[1]*f[4]; 
  out[25] += 22.9128784747792*rdxSq4nu[0]*f[5]; 
  out[26] += 22.9128784747792*rdxSq4nu[1]*f[6]; 
  out[27] += 22.9128784747792*rdxSq4nu[2]*f[5]; 
  out[28] += 22.9128784747792*rdxSq4nu[2]*f[6]; 
  out[29] += 22.9128784747792*rdxSq4nu[0]*f[10]; 
  out[30] += 22.9128784747792*rdxSq4nu[1]*f[10]; 
  out[31] += 22.9128784747792*rdxSq4nu[2]*f[10]; 

return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
double ConstDiffusionVol3xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[7] += 6.708203932499369*f[0]*rdxSq4nu[0]; 
  out[8] += 6.708203932499369*f[0]*rdxSq4nu[1]; 
  out[9] += 6.708203932499369*f[0]*rdxSq4nu[2]; 
  out[11] += 6.708203932499369*rdxSq4nu[0]*f[2]; 
  out[12] += 6.708203932499369*f[1]*rdxSq4nu[1]; 
  out[13] += 6.708203932499369*rdxSq4nu[0]*f[3]; 
  out[14] += 6.708203932499369*rdxSq4nu[1]*f[3]; 
  out[15] += 6.708203932499369*f[1]*rdxSq4nu[2]; 
  out[16] += 6.708203932499369*f[2]*rdxSq4nu[2]; 
  out[17] += 22.9128784747792*rdxSq4nu[0]*f[1]; 
  out[18] += 22.9128784747792*rdxSq4nu[1]*f[2]; 
  out[19] += 22.9128784747792*rdxSq4nu[2]*f[3]; 
  out[20] += 6.708203932499369*rdxSq4nu[0]*f[6]; 
  out[21] += 6.708203932499369*rdxSq4nu[1]*f[5]; 
  out[22] += 6.708203932499369*rdxSq4nu[2]*f[4]; 
  out[23] += 6.708203932499369*rdxSq4nu[0]*f[8]+6.708203932499369*rdxSq4nu[1]*f[7]; 
  out[24] += 6.708203932499369*rdxSq4nu[0]*f[9]+6.708203932499369*rdxSq4nu[2]*f[7]; 
  out[25] += 6.708203932499369*rdxSq4nu[1]*f[9]+6.708203932499369*rdxSq4nu[2]*f[8]; 
  out[26] += 22.9128784747792*rdxSq4nu[0]*f[4]; 
  out[27] += 22.9128784747792*rdxSq4nu[1]*f[4]; 
  out[28] += 22.9128784747792*rdxSq4nu[0]*f[5]; 
  out[29] += 22.9128784747792*rdxSq4nu[1]*f[6]; 
  out[30] += 22.9128784747792*rdxSq4nu[2]*f[5]; 
  out[31] += 22.9128784747792*rdxSq4nu[2]*f[6]; 
  out[32] += 46.95742752749558*rdxSq4nu[0]*f[7]+30.0*f[0]*rdxSq4nu[0]; 
  out[33] += 46.95742752749558*rdxSq4nu[1]*f[8]+30.0*f[0]*rdxSq4nu[1]; 
  out[34] += 46.95742752749558*rdxSq4nu[2]*f[9]+30.0*f[0]*rdxSq4nu[2]; 
  out[35] += 6.708203932499369*rdxSq4nu[0]*f[14]+6.708203932499369*rdxSq4nu[1]*f[13]; 
  out[36] += 6.708203932499369*rdxSq4nu[0]*f[16]+6.708203932499369*rdxSq4nu[2]*f[11]; 
  out[37] += 6.708203932499369*rdxSq4nu[1]*f[15]+6.708203932499369*rdxSq4nu[2]*f[12]; 
  out[38] += 22.9128784747792*rdxSq4nu[0]*f[10]; 
  out[39] += 22.9128784747792*rdxSq4nu[1]*f[10]; 
  out[40] += 22.9128784747792*rdxSq4nu[2]*f[10]; 
  out[41] += 46.95742752749558*rdxSq4nu[0]*f[11]+30.0*rdxSq4nu[0]*f[2]; 
  out[42] += 46.95742752749558*rdxSq4nu[1]*f[12]+30.0*f[1]*rdxSq4nu[1]; 
  out[43] += 46.95742752749558*rdxSq4nu[0]*f[13]+30.0*rdxSq4nu[0]*f[3]; 
  out[44] += 46.95742752749558*rdxSq4nu[1]*f[14]+30.0*rdxSq4nu[1]*f[3]; 
  out[45] += 46.95742752749558*rdxSq4nu[2]*f[15]+30.0*f[1]*rdxSq4nu[2]; 
  out[46] += 46.95742752749558*rdxSq4nu[2]*f[16]+30.0*f[2]*rdxSq4nu[2]; 
  out[47] += 46.95742752749558*rdxSq4nu[0]*f[20]+30.0*rdxSq4nu[0]*f[6]; 
  out[48] += 46.95742752749558*rdxSq4nu[1]*f[21]+30.0*rdxSq4nu[1]*f[5]; 
  out[49] += 46.95742752749558*rdxSq4nu[2]*f[22]+30.0*rdxSq4nu[2]*f[4]; 

return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
