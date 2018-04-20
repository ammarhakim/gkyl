#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xSerP1(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[1] += -3.0*rdxSq4nu[0]*f[1]; 
  out[2] += -3.0*rdxSq4nu[1]*f[2]; 
  out[3] += -3.0*rdxSq4nu[2]*f[3]; 
  out[4] += (-3.0*rdxSq4nu[1]*f[4])-3.0*rdxSq4nu[0]*f[4]; 
  out[5] += (-3.0*rdxSq4nu[2]*f[5])-3.0*rdxSq4nu[0]*f[5]; 
  out[6] += (-3.0*rdxSq4nu[2]*f[6])-3.0*rdxSq4nu[1]*f[6]; 
  out[7] += (-3.0*rdxSq4nu[2]*f[7])-3.0*rdxSq4nu[1]*f[7]-3.0*rdxSq4nu[0]*f[7]; 

return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
double ConstDiffusionVol3xSerP2(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[1] += -3.0*rdxSq4nu[0]*f[1]; 
  out[2] += -3.0*rdxSq4nu[1]*f[2]; 
  out[3] += -3.0*rdxSq4nu[2]*f[3]; 
  out[4] += (-3.0*rdxSq4nu[1]*f[4])-3.0*rdxSq4nu[0]*f[4]; 
  out[5] += (-3.0*rdxSq4nu[2]*f[5])-3.0*rdxSq4nu[0]*f[5]; 
  out[6] += (-3.0*rdxSq4nu[2]*f[6])-3.0*rdxSq4nu[1]*f[6]; 
  out[7] += -15.0*rdxSq4nu[0]*f[7]; 
  out[8] += -15.0*rdxSq4nu[1]*f[8]; 
  out[9] += -15.0*rdxSq4nu[2]*f[9]; 
  out[10] += (-3.0*rdxSq4nu[2]*f[10])-3.0*rdxSq4nu[1]*f[10]-3.0*rdxSq4nu[0]*f[10]; 
  out[11] += (-3.0*rdxSq4nu[1]*f[11])-15.0*rdxSq4nu[0]*f[11]; 
  out[12] += (-15.0*rdxSq4nu[1]*f[12])-3.0*rdxSq4nu[0]*f[12]; 
  out[13] += (-3.0*rdxSq4nu[2]*f[13])-15.0*rdxSq4nu[0]*f[13]; 
  out[14] += (-3.0*rdxSq4nu[2]*f[14])-15.0*rdxSq4nu[1]*f[14]; 
  out[15] += (-15.0*rdxSq4nu[2]*f[15])-3.0*rdxSq4nu[0]*f[15]; 
  out[16] += (-15.0*rdxSq4nu[2]*f[16])-3.0*rdxSq4nu[1]*f[16]; 
  out[17] += (-3.0*rdxSq4nu[2]*f[17])-3.0*rdxSq4nu[1]*f[17]-15.0*rdxSq4nu[0]*f[17]; 
  out[18] += (-3.0*rdxSq4nu[2]*f[18])-15.0*rdxSq4nu[1]*f[18]-3.0*rdxSq4nu[0]*f[18]; 
  out[19] += (-15.0*rdxSq4nu[2]*f[19])-3.0*rdxSq4nu[1]*f[19]-3.0*rdxSq4nu[0]*f[19]; 

return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
double ConstDiffusionVol3xSerP3(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[1] += (-4.58257569495584*rdxSq4nu[0]*f[17])-3.0*rdxSq4nu[0]*f[1]; 
  out[2] += (-4.58257569495584*rdxSq4nu[1]*f[18])-3.0*rdxSq4nu[1]*f[2]; 
  out[3] += (-4.58257569495584*rdxSq4nu[2]*f[19])-3.0*rdxSq4nu[2]*f[3]; 
  out[4] += (-4.58257569495584*rdxSq4nu[1]*f[24])-4.58257569495584*rdxSq4nu[0]*f[23]-3.0*rdxSq4nu[1]*f[4]-3.0*rdxSq4nu[0]*f[4]; 
  out[5] += (-4.58257569495584*rdxSq4nu[2]*f[27])-4.58257569495584*rdxSq4nu[0]*f[25]-3.0*rdxSq4nu[2]*f[5]-3.0*rdxSq4nu[0]*f[5]; 
  out[6] += (-4.58257569495584*rdxSq4nu[2]*f[28])-4.58257569495584*rdxSq4nu[1]*f[26]-3.0*rdxSq4nu[2]*f[6]-3.0*rdxSq4nu[1]*f[6]; 
  out[7] += -15.0*rdxSq4nu[0]*f[7]; 
  out[8] += -15.0*rdxSq4nu[1]*f[8]; 
  out[9] += -15.0*rdxSq4nu[2]*f[9]; 
  out[10] += (-4.58257569495584*rdxSq4nu[2]*f[31])-4.58257569495584*rdxSq4nu[1]*f[30]-4.58257569495584*rdxSq4nu[0]*f[29]-3.0*rdxSq4nu[2]*f[10]-3.0*rdxSq4nu[1]*f[10]-3.0*rdxSq4nu[0]*f[10]; 
  out[11] += (-3.0*rdxSq4nu[1]*f[11])-15.0*rdxSq4nu[0]*f[11]; 
  out[12] += (-15.0*rdxSq4nu[1]*f[12])-3.0*rdxSq4nu[0]*f[12]; 
  out[13] += (-3.0*rdxSq4nu[2]*f[13])-15.0*rdxSq4nu[0]*f[13]; 
  out[14] += (-3.0*rdxSq4nu[2]*f[14])-15.0*rdxSq4nu[1]*f[14]; 
  out[15] += (-15.0*rdxSq4nu[2]*f[15])-3.0*rdxSq4nu[0]*f[15]; 
  out[16] += (-15.0*rdxSq4nu[2]*f[16])-3.0*rdxSq4nu[1]*f[16]; 
  out[17] += (-42.0*rdxSq4nu[0]*f[17])-4.58257569495584*rdxSq4nu[0]*f[1]; 
  out[18] += (-42.0*rdxSq4nu[1]*f[18])-4.58257569495584*rdxSq4nu[1]*f[2]; 
  out[19] += (-42.0*rdxSq4nu[2]*f[19])-4.58257569495584*rdxSq4nu[2]*f[3]; 
  out[20] += (-3.0*rdxSq4nu[2]*f[20])-3.0*rdxSq4nu[1]*f[20]-15.0*rdxSq4nu[0]*f[20]; 
  out[21] += (-3.0*rdxSq4nu[2]*f[21])-15.0*rdxSq4nu[1]*f[21]-3.0*rdxSq4nu[0]*f[21]; 
  out[22] += (-15.0*rdxSq4nu[2]*f[22])-3.0*rdxSq4nu[1]*f[22]-3.0*rdxSq4nu[0]*f[22]; 
  out[23] += (-3.0*rdxSq4nu[1]*f[23])-42.0*rdxSq4nu[0]*f[23]-4.58257569495584*rdxSq4nu[0]*f[4]; 
  out[24] += (-42.0*rdxSq4nu[1]*f[24])-3.0*rdxSq4nu[0]*f[24]-4.58257569495584*rdxSq4nu[1]*f[4]; 
  out[25] += (-3.0*rdxSq4nu[2]*f[25])-42.0*rdxSq4nu[0]*f[25]-4.58257569495584*rdxSq4nu[0]*f[5]; 
  out[26] += (-3.0*rdxSq4nu[2]*f[26])-42.0*rdxSq4nu[1]*f[26]-4.58257569495584*rdxSq4nu[1]*f[6]; 
  out[27] += (-42.0*rdxSq4nu[2]*f[27])-3.0*rdxSq4nu[0]*f[27]-4.58257569495584*rdxSq4nu[2]*f[5]; 
  out[28] += (-42.0*rdxSq4nu[2]*f[28])-3.0*rdxSq4nu[1]*f[28]-4.58257569495584*rdxSq4nu[2]*f[6]; 
  out[29] += (-3.0*rdxSq4nu[2]*f[29])-3.0*rdxSq4nu[1]*f[29]-42.0*rdxSq4nu[0]*f[29]-4.58257569495584*rdxSq4nu[0]*f[10]; 
  out[30] += (-3.0*rdxSq4nu[2]*f[30])-42.0*rdxSq4nu[1]*f[30]-3.0*rdxSq4nu[0]*f[30]-4.58257569495584*rdxSq4nu[1]*f[10]; 
  out[31] += (-42.0*rdxSq4nu[2]*f[31])-3.0*rdxSq4nu[1]*f[31]-3.0*rdxSq4nu[0]*f[31]-4.58257569495584*rdxSq4nu[2]*f[10]; 

return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
double ConstDiffusionVol3xSerP4(const double *w, const double *dxv, const double *nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). f: Input distribution function. out: Incremented output 
  double rdxSq4nu[3]; 
  rdxSq4nu[0] = 4.0*nu[0]/(dxv[0]*dxv[0]); 
  rdxSq4nu[1] = 4.0*nu[1]/(dxv[1]*dxv[1]); 
  rdxSq4nu[2] = 4.0*nu[2]/(dxv[2]*dxv[2]); 

  out[1] += (-4.58257569495584*rdxSq4nu[0]*f[17])-3.0*rdxSq4nu[0]*f[1]; 
  out[2] += (-4.58257569495584*rdxSq4nu[1]*f[18])-3.0*rdxSq4nu[1]*f[2]; 
  out[3] += (-4.58257569495584*rdxSq4nu[2]*f[19])-3.0*rdxSq4nu[2]*f[3]; 
  out[4] += (-4.58257569495584*rdxSq4nu[1]*f[27])-4.58257569495584*rdxSq4nu[0]*f[26]-3.0*rdxSq4nu[1]*f[4]-3.0*rdxSq4nu[0]*f[4]; 
  out[5] += (-4.58257569495584*rdxSq4nu[2]*f[30])-4.58257569495584*rdxSq4nu[0]*f[28]-3.0*rdxSq4nu[2]*f[5]-3.0*rdxSq4nu[0]*f[5]; 
  out[6] += (-4.58257569495584*rdxSq4nu[2]*f[31])-4.58257569495584*rdxSq4nu[1]*f[29]-3.0*rdxSq4nu[2]*f[6]-3.0*rdxSq4nu[1]*f[6]; 
  out[7] += (-20.12461179749811*rdxSq4nu[0]*f[32])-15.0*rdxSq4nu[0]*f[7]; 
  out[8] += (-20.12461179749811*rdxSq4nu[1]*f[33])-15.0*rdxSq4nu[1]*f[8]; 
  out[9] += (-20.12461179749811*rdxSq4nu[2]*f[34])-15.0*rdxSq4nu[2]*f[9]; 
  out[10] += (-4.58257569495584*rdxSq4nu[2]*f[40])-4.58257569495584*rdxSq4nu[1]*f[39]-4.58257569495584*rdxSq4nu[0]*f[38]-3.0*rdxSq4nu[2]*f[10]-3.0*rdxSq4nu[1]*f[10]-3.0*rdxSq4nu[0]*f[10]; 
  out[11] += (-20.1246117974981*rdxSq4nu[0]*f[41])-3.0*rdxSq4nu[1]*f[11]-15.0*rdxSq4nu[0]*f[11]; 
  out[12] += (-20.1246117974981*rdxSq4nu[1]*f[42])-15.0*rdxSq4nu[1]*f[12]-3.0*rdxSq4nu[0]*f[12]; 
  out[13] += (-20.1246117974981*rdxSq4nu[0]*f[43])-3.0*rdxSq4nu[2]*f[13]-15.0*rdxSq4nu[0]*f[13]; 
  out[14] += (-20.1246117974981*rdxSq4nu[1]*f[44])-3.0*rdxSq4nu[2]*f[14]-15.0*rdxSq4nu[1]*f[14]; 
  out[15] += (-20.1246117974981*rdxSq4nu[2]*f[45])-15.0*rdxSq4nu[2]*f[15]-3.0*rdxSq4nu[0]*f[15]; 
  out[16] += (-20.1246117974981*rdxSq4nu[2]*f[46])-15.0*rdxSq4nu[2]*f[16]-3.0*rdxSq4nu[1]*f[16]; 
  out[17] += (-42.0*rdxSq4nu[0]*f[17])-4.58257569495584*rdxSq4nu[0]*f[1]; 
  out[18] += (-42.0*rdxSq4nu[1]*f[18])-4.58257569495584*rdxSq4nu[1]*f[2]; 
  out[19] += (-42.0*rdxSq4nu[2]*f[19])-4.58257569495584*rdxSq4nu[2]*f[3]; 
  out[20] += (-20.12461179749811*rdxSq4nu[0]*f[47])-3.0*rdxSq4nu[2]*f[20]-3.0*rdxSq4nu[1]*f[20]-15.0*rdxSq4nu[0]*f[20]; 
  out[21] += (-20.12461179749811*rdxSq4nu[1]*f[48])-3.0*rdxSq4nu[2]*f[21]-15.0*rdxSq4nu[1]*f[21]-3.0*rdxSq4nu[0]*f[21]; 
  out[22] += (-20.12461179749811*rdxSq4nu[2]*f[49])-15.0*rdxSq4nu[2]*f[22]-3.0*rdxSq4nu[1]*f[22]-3.0*rdxSq4nu[0]*f[22]; 
  out[23] += (-15.0*rdxSq4nu[1]*f[23])-15.0*rdxSq4nu[0]*f[23]; 
  out[24] += (-15.0*rdxSq4nu[2]*f[24])-15.0*rdxSq4nu[0]*f[24]; 
  out[25] += (-15.0*rdxSq4nu[2]*f[25])-15.0*rdxSq4nu[1]*f[25]; 
  out[26] += (-3.0*rdxSq4nu[1]*f[26])-42.0*rdxSq4nu[0]*f[26]-4.58257569495584*rdxSq4nu[0]*f[4]; 
  out[27] += (-42.0*rdxSq4nu[1]*f[27])-3.0*rdxSq4nu[0]*f[27]-4.58257569495584*rdxSq4nu[1]*f[4]; 
  out[28] += (-3.0*rdxSq4nu[2]*f[28])-42.0*rdxSq4nu[0]*f[28]-4.58257569495584*rdxSq4nu[0]*f[5]; 
  out[29] += (-3.0*rdxSq4nu[2]*f[29])-42.0*rdxSq4nu[1]*f[29]-4.58257569495584*rdxSq4nu[1]*f[6]; 
  out[30] += (-42.0*rdxSq4nu[2]*f[30])-3.0*rdxSq4nu[0]*f[30]-4.58257569495584*rdxSq4nu[2]*f[5]; 
  out[31] += (-42.0*rdxSq4nu[2]*f[31])-3.0*rdxSq4nu[1]*f[31]-4.58257569495584*rdxSq4nu[2]*f[6]; 
  out[32] += (-90.0*rdxSq4nu[0]*f[32])-20.12461179749811*rdxSq4nu[0]*f[7]; 
  out[33] += (-90.0*rdxSq4nu[1]*f[33])-20.12461179749811*rdxSq4nu[1]*f[8]; 
  out[34] += (-90.0*rdxSq4nu[2]*f[34])-20.12461179749811*rdxSq4nu[2]*f[9]; 
  out[35] += (-3.0*rdxSq4nu[2]*f[35])-15.0*rdxSq4nu[1]*f[35]-15.0*rdxSq4nu[0]*f[35]; 
  out[36] += (-15.0*rdxSq4nu[2]*f[36])-3.0*rdxSq4nu[1]*f[36]-15.0*rdxSq4nu[0]*f[36]; 
  out[37] += (-15.0*rdxSq4nu[2]*f[37])-15.0*rdxSq4nu[1]*f[37]-3.0*rdxSq4nu[0]*f[37]; 
  out[38] += (-3.0*rdxSq4nu[2]*f[38])-3.0*rdxSq4nu[1]*f[38]-42.0*rdxSq4nu[0]*f[38]-4.58257569495584*rdxSq4nu[0]*f[10]; 
  out[39] += (-3.0*rdxSq4nu[2]*f[39])-42.0*rdxSq4nu[1]*f[39]-3.0*rdxSq4nu[0]*f[39]-4.58257569495584*rdxSq4nu[1]*f[10]; 
  out[40] += (-42.0*rdxSq4nu[2]*f[40])-3.0*rdxSq4nu[1]*f[40]-3.0*rdxSq4nu[0]*f[40]-4.58257569495584*rdxSq4nu[2]*f[10]; 
  out[41] += (-3.0*rdxSq4nu[1]*f[41])-90.0*rdxSq4nu[0]*f[41]-20.1246117974981*rdxSq4nu[0]*f[11]; 
  out[42] += (-90.0*rdxSq4nu[1]*f[42])-3.0*rdxSq4nu[0]*f[42]-20.1246117974981*rdxSq4nu[1]*f[12]; 
  out[43] += (-3.0*rdxSq4nu[2]*f[43])-90.0*rdxSq4nu[0]*f[43]-20.1246117974981*rdxSq4nu[0]*f[13]; 
  out[44] += (-3.0*rdxSq4nu[2]*f[44])-90.0*rdxSq4nu[1]*f[44]-20.1246117974981*rdxSq4nu[1]*f[14]; 
  out[45] += (-90.0*rdxSq4nu[2]*f[45])-3.0*rdxSq4nu[0]*f[45]-20.1246117974981*rdxSq4nu[2]*f[15]; 
  out[46] += (-90.0*rdxSq4nu[2]*f[46])-3.0*rdxSq4nu[1]*f[46]-20.1246117974981*rdxSq4nu[2]*f[16]; 
  out[47] += (-3.0*rdxSq4nu[2]*f[47])-3.0*rdxSq4nu[1]*f[47]-90.0*rdxSq4nu[0]*f[47]-20.12461179749811*rdxSq4nu[0]*f[20]; 
  out[48] += (-3.0*rdxSq4nu[2]*f[48])-90.0*rdxSq4nu[1]*f[48]-3.0*rdxSq4nu[0]*f[48]-20.12461179749811*rdxSq4nu[1]*f[21]; 
  out[49] += (-90.0*rdxSq4nu[2]*f[49])-3.0*rdxSq4nu[1]*f[49]-3.0*rdxSq4nu[0]*f[49]-20.12461179749811*rdxSq4nu[2]*f[22]; 

return (rdxSq4nu[0]+rdxSq4nu[1]+rdxSq4nu[2])*0.5; 

} 
