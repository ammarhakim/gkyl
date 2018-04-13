#include <ConstDiffusionModDecl.h> 
double ConstDiffusionVol3xSerP1(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[3]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 
  rdx[2] = 4.0/(dxv[2]*dxv[2]); 

  out[1] += 3.0*rdx[0]*f[1]*nu; 
  out[2] += 3.0*rdx[1]*f[2]*nu; 
  out[3] += 3.0*rdx[2]*f[3]*nu; 
  out[4] += 3.0*rdx[1]*f[4]*nu+3.0*rdx[0]*f[4]*nu; 
  out[5] += 3.0*rdx[2]*f[5]*nu+3.0*rdx[0]*f[5]*nu; 
  out[6] += 3.0*rdx[2]*f[6]*nu+3.0*rdx[1]*f[6]*nu; 
  out[7] += 3.0*rdx[2]*f[7]*nu+3.0*rdx[1]*f[7]*nu+3.0*rdx[0]*f[7]*nu; 

return nu*(rdx[0]+rdx[1]+rdx[2])*0.5; 

} 
double ConstDiffusionVol3xSerP2(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[3]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 
  rdx[2] = 4.0/(dxv[2]*dxv[2]); 

  out[1] += 3.0*rdx[0]*f[1]*nu; 
  out[2] += 3.0*rdx[1]*f[2]*nu; 
  out[3] += 3.0*rdx[2]*f[3]*nu; 
  out[4] += 3.0*rdx[1]*f[4]*nu+3.0*rdx[0]*f[4]*nu; 
  out[5] += 3.0*rdx[2]*f[5]*nu+3.0*rdx[0]*f[5]*nu; 
  out[6] += 3.0*rdx[2]*f[6]*nu+3.0*rdx[1]*f[6]*nu; 
  out[7] += 15.0*rdx[0]*f[7]*nu; 
  out[8] += 15.0*rdx[1]*f[8]*nu; 
  out[9] += 15.0*rdx[2]*f[9]*nu; 
  out[10] += 3.0*rdx[2]*f[10]*nu+3.0*rdx[1]*f[10]*nu+3.0*rdx[0]*f[10]*nu; 
  out[11] += 3.0*rdx[1]*f[11]*nu+15.0*rdx[0]*f[11]*nu; 
  out[12] += 15.0*rdx[1]*f[12]*nu+3.0*rdx[0]*f[12]*nu; 
  out[13] += 3.0*rdx[2]*f[13]*nu+15.0*rdx[0]*f[13]*nu; 
  out[14] += 3.0*rdx[2]*f[14]*nu+15.0*rdx[1]*f[14]*nu; 
  out[15] += 15.0*rdx[2]*f[15]*nu+3.0*rdx[0]*f[15]*nu; 
  out[16] += 15.0*rdx[2]*f[16]*nu+3.0*rdx[1]*f[16]*nu; 
  out[17] += 3.0*rdx[2]*f[17]*nu+3.0*rdx[1]*f[17]*nu+15.0*rdx[0]*f[17]*nu; 
  out[18] += 3.0*rdx[2]*f[18]*nu+15.0*rdx[1]*f[18]*nu+3.0*rdx[0]*f[18]*nu; 
  out[19] += 15.0*rdx[2]*f[19]*nu+3.0*rdx[1]*f[19]*nu+3.0*rdx[0]*f[19]*nu; 

return nu*(rdx[0]+rdx[1]+rdx[2])*0.5; 

} 
double ConstDiffusionVol3xSerP3(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[3]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 
  rdx[2] = 4.0/(dxv[2]*dxv[2]); 

  out[1] += 4.58257569495584*rdx[0]*f[17]*nu+3.0*rdx[0]*f[1]*nu; 
  out[2] += 4.58257569495584*rdx[1]*f[18]*nu+3.0*rdx[1]*f[2]*nu; 
  out[3] += 4.58257569495584*rdx[2]*f[19]*nu+3.0*rdx[2]*f[3]*nu; 
  out[4] += 4.58257569495584*rdx[1]*f[24]*nu+4.58257569495584*rdx[0]*f[23]*nu+3.0*rdx[1]*f[4]*nu+3.0*rdx[0]*f[4]*nu; 
  out[5] += 4.58257569495584*rdx[2]*f[27]*nu+4.58257569495584*rdx[0]*f[25]*nu+3.0*rdx[2]*f[5]*nu+3.0*rdx[0]*f[5]*nu; 
  out[6] += 4.58257569495584*rdx[2]*f[28]*nu+4.58257569495584*rdx[1]*f[26]*nu+3.0*rdx[2]*f[6]*nu+3.0*rdx[1]*f[6]*nu; 
  out[7] += 15.0*rdx[0]*f[7]*nu; 
  out[8] += 15.0*rdx[1]*f[8]*nu; 
  out[9] += 15.0*rdx[2]*f[9]*nu; 
  out[10] += 4.58257569495584*rdx[2]*f[31]*nu+4.58257569495584*rdx[1]*f[30]*nu+4.58257569495584*rdx[0]*f[29]*nu+3.0*rdx[2]*f[10]*nu+3.0*rdx[1]*f[10]*nu+3.0*rdx[0]*f[10]*nu; 
  out[11] += 3.0*rdx[1]*f[11]*nu+15.0*rdx[0]*f[11]*nu; 
  out[12] += 15.0*rdx[1]*f[12]*nu+3.0*rdx[0]*f[12]*nu; 
  out[13] += 3.0*rdx[2]*f[13]*nu+15.0*rdx[0]*f[13]*nu; 
  out[14] += 3.0*rdx[2]*f[14]*nu+15.0*rdx[1]*f[14]*nu; 
  out[15] += 15.0*rdx[2]*f[15]*nu+3.0*rdx[0]*f[15]*nu; 
  out[16] += 15.0*rdx[2]*f[16]*nu+3.0*rdx[1]*f[16]*nu; 
  out[17] += 42.0*rdx[0]*f[17]*nu+4.58257569495584*rdx[0]*f[1]*nu; 
  out[18] += 42.0*rdx[1]*f[18]*nu+4.58257569495584*rdx[1]*f[2]*nu; 
  out[19] += 42.0*rdx[2]*f[19]*nu+4.58257569495584*rdx[2]*f[3]*nu; 
  out[20] += 3.0*rdx[2]*f[20]*nu+3.0*rdx[1]*f[20]*nu+15.0*rdx[0]*f[20]*nu; 
  out[21] += 3.0*rdx[2]*f[21]*nu+15.0*rdx[1]*f[21]*nu+3.0*rdx[0]*f[21]*nu; 
  out[22] += 15.0*rdx[2]*f[22]*nu+3.0*rdx[1]*f[22]*nu+3.0*rdx[0]*f[22]*nu; 
  out[23] += 3.0*rdx[1]*f[23]*nu+42.0*rdx[0]*f[23]*nu+4.58257569495584*rdx[0]*f[4]*nu; 
  out[24] += 42.0*rdx[1]*f[24]*nu+3.0*rdx[0]*f[24]*nu+4.58257569495584*rdx[1]*f[4]*nu; 
  out[25] += 3.0*rdx[2]*f[25]*nu+42.0*rdx[0]*f[25]*nu+4.58257569495584*rdx[0]*f[5]*nu; 
  out[26] += 3.0*rdx[2]*f[26]*nu+42.0*rdx[1]*f[26]*nu+4.58257569495584*rdx[1]*f[6]*nu; 
  out[27] += 42.0*rdx[2]*f[27]*nu+3.0*rdx[0]*f[27]*nu+4.58257569495584*rdx[2]*f[5]*nu; 
  out[28] += 42.0*rdx[2]*f[28]*nu+3.0*rdx[1]*f[28]*nu+4.58257569495584*rdx[2]*f[6]*nu; 
  out[29] += 3.0*rdx[2]*f[29]*nu+3.0*rdx[1]*f[29]*nu+42.0*rdx[0]*f[29]*nu+4.58257569495584*rdx[0]*f[10]*nu; 
  out[30] += 3.0*rdx[2]*f[30]*nu+42.0*rdx[1]*f[30]*nu+3.0*rdx[0]*f[30]*nu+4.58257569495584*rdx[1]*f[10]*nu; 
  out[31] += 42.0*rdx[2]*f[31]*nu+3.0*rdx[1]*f[31]*nu+3.0*rdx[0]*f[31]*nu+4.58257569495584*rdx[2]*f[10]*nu; 

return nu*(rdx[0]+rdx[1]+rdx[2])*0.5; 

} 
double ConstDiffusionVol3xSerP4(const double *w, const double *dxv, const double nu, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double rdx[3]; 
  rdx[0] = 4.0/(dxv[0]*dxv[0]); 
  rdx[1] = 4.0/(dxv[1]*dxv[1]); 
  rdx[2] = 4.0/(dxv[2]*dxv[2]); 

  out[1] += 4.58257569495584*rdx[0]*f[17]*nu+3.0*rdx[0]*f[1]*nu; 
  out[2] += 4.58257569495584*rdx[1]*f[18]*nu+3.0*rdx[1]*f[2]*nu; 
  out[3] += 4.58257569495584*rdx[2]*f[19]*nu+3.0*rdx[2]*f[3]*nu; 
  out[4] += 4.58257569495584*rdx[1]*f[27]*nu+4.58257569495584*rdx[0]*f[26]*nu+3.0*rdx[1]*f[4]*nu+3.0*rdx[0]*f[4]*nu; 
  out[5] += 4.58257569495584*rdx[2]*f[30]*nu+4.58257569495584*rdx[0]*f[28]*nu+3.0*rdx[2]*f[5]*nu+3.0*rdx[0]*f[5]*nu; 
  out[6] += 4.58257569495584*rdx[2]*f[31]*nu+4.58257569495584*rdx[1]*f[29]*nu+3.0*rdx[2]*f[6]*nu+3.0*rdx[1]*f[6]*nu; 
  out[7] += 20.12461179749811*rdx[0]*f[32]*nu+15.0*rdx[0]*f[7]*nu; 
  out[8] += 20.12461179749811*rdx[1]*f[33]*nu+15.0*rdx[1]*f[8]*nu; 
  out[9] += 20.12461179749811*rdx[2]*f[34]*nu+15.0*rdx[2]*f[9]*nu; 
  out[10] += 4.58257569495584*rdx[2]*f[40]*nu+4.58257569495584*rdx[1]*f[39]*nu+4.58257569495584*rdx[0]*f[38]*nu+3.0*rdx[2]*f[10]*nu+3.0*rdx[1]*f[10]*nu+3.0*rdx[0]*f[10]*nu; 
  out[11] += 20.1246117974981*rdx[0]*f[41]*nu+3.0*rdx[1]*f[11]*nu+15.0*rdx[0]*f[11]*nu; 
  out[12] += 20.1246117974981*rdx[1]*f[42]*nu+15.0*rdx[1]*f[12]*nu+3.0*rdx[0]*f[12]*nu; 
  out[13] += 20.1246117974981*rdx[0]*f[43]*nu+3.0*rdx[2]*f[13]*nu+15.0*rdx[0]*f[13]*nu; 
  out[14] += 20.1246117974981*rdx[1]*f[44]*nu+3.0*rdx[2]*f[14]*nu+15.0*rdx[1]*f[14]*nu; 
  out[15] += 20.1246117974981*rdx[2]*f[45]*nu+15.0*rdx[2]*f[15]*nu+3.0*rdx[0]*f[15]*nu; 
  out[16] += 20.1246117974981*rdx[2]*f[46]*nu+15.0*rdx[2]*f[16]*nu+3.0*rdx[1]*f[16]*nu; 
  out[17] += 42.0*rdx[0]*f[17]*nu+4.58257569495584*rdx[0]*f[1]*nu; 
  out[18] += 42.0*rdx[1]*f[18]*nu+4.58257569495584*rdx[1]*f[2]*nu; 
  out[19] += 42.0*rdx[2]*f[19]*nu+4.58257569495584*rdx[2]*f[3]*nu; 
  out[20] += 20.12461179749811*rdx[0]*f[47]*nu+3.0*rdx[2]*f[20]*nu+3.0*rdx[1]*f[20]*nu+15.0*rdx[0]*f[20]*nu; 
  out[21] += 20.12461179749811*rdx[1]*f[48]*nu+3.0*rdx[2]*f[21]*nu+15.0*rdx[1]*f[21]*nu+3.0*rdx[0]*f[21]*nu; 
  out[22] += 20.12461179749811*rdx[2]*f[49]*nu+15.0*rdx[2]*f[22]*nu+3.0*rdx[1]*f[22]*nu+3.0*rdx[0]*f[22]*nu; 
  out[23] += 15.0*rdx[1]*f[23]*nu+15.0*rdx[0]*f[23]*nu; 
  out[24] += 15.0*rdx[2]*f[24]*nu+15.0*rdx[0]*f[24]*nu; 
  out[25] += 15.0*rdx[2]*f[25]*nu+15.0*rdx[1]*f[25]*nu; 
  out[26] += 3.0*rdx[1]*f[26]*nu+42.0*rdx[0]*f[26]*nu+4.58257569495584*rdx[0]*f[4]*nu; 
  out[27] += 42.0*rdx[1]*f[27]*nu+3.0*rdx[0]*f[27]*nu+4.58257569495584*rdx[1]*f[4]*nu; 
  out[28] += 3.0*rdx[2]*f[28]*nu+42.0*rdx[0]*f[28]*nu+4.58257569495584*rdx[0]*f[5]*nu; 
  out[29] += 3.0*rdx[2]*f[29]*nu+42.0*rdx[1]*f[29]*nu+4.58257569495584*rdx[1]*f[6]*nu; 
  out[30] += 42.0*rdx[2]*f[30]*nu+3.0*rdx[0]*f[30]*nu+4.58257569495584*rdx[2]*f[5]*nu; 
  out[31] += 42.0*rdx[2]*f[31]*nu+3.0*rdx[1]*f[31]*nu+4.58257569495584*rdx[2]*f[6]*nu; 
  out[32] += 90.0*rdx[0]*f[32]*nu+20.12461179749811*rdx[0]*f[7]*nu; 
  out[33] += 90.0*rdx[1]*f[33]*nu+20.12461179749811*rdx[1]*f[8]*nu; 
  out[34] += 90.0*rdx[2]*f[34]*nu+20.12461179749811*rdx[2]*f[9]*nu; 
  out[35] += 3.0*rdx[2]*f[35]*nu+15.0*rdx[1]*f[35]*nu+15.0*rdx[0]*f[35]*nu; 
  out[36] += 15.0*rdx[2]*f[36]*nu+3.0*rdx[1]*f[36]*nu+15.0*rdx[0]*f[36]*nu; 
  out[37] += 15.0*rdx[2]*f[37]*nu+15.0*rdx[1]*f[37]*nu+3.0*rdx[0]*f[37]*nu; 
  out[38] += 3.0*rdx[2]*f[38]*nu+3.0*rdx[1]*f[38]*nu+42.0*rdx[0]*f[38]*nu+4.58257569495584*rdx[0]*f[10]*nu; 
  out[39] += 3.0*rdx[2]*f[39]*nu+42.0*rdx[1]*f[39]*nu+3.0*rdx[0]*f[39]*nu+4.58257569495584*rdx[1]*f[10]*nu; 
  out[40] += 42.0*rdx[2]*f[40]*nu+3.0*rdx[1]*f[40]*nu+3.0*rdx[0]*f[40]*nu+4.58257569495584*rdx[2]*f[10]*nu; 
  out[41] += 3.0*rdx[1]*f[41]*nu+90.0*rdx[0]*f[41]*nu+20.1246117974981*rdx[0]*f[11]*nu; 
  out[42] += 90.0*rdx[1]*f[42]*nu+3.0*rdx[0]*f[42]*nu+20.1246117974981*rdx[1]*f[12]*nu; 
  out[43] += 3.0*rdx[2]*f[43]*nu+90.0*rdx[0]*f[43]*nu+20.1246117974981*rdx[0]*f[13]*nu; 
  out[44] += 3.0*rdx[2]*f[44]*nu+90.0*rdx[1]*f[44]*nu+20.1246117974981*rdx[1]*f[14]*nu; 
  out[45] += 90.0*rdx[2]*f[45]*nu+3.0*rdx[0]*f[45]*nu+20.1246117974981*rdx[2]*f[15]*nu; 
  out[46] += 90.0*rdx[2]*f[46]*nu+3.0*rdx[1]*f[46]*nu+20.1246117974981*rdx[2]*f[16]*nu; 
  out[47] += 3.0*rdx[2]*f[47]*nu+3.0*rdx[1]*f[47]*nu+90.0*rdx[0]*f[47]*nu+20.12461179749811*rdx[0]*f[20]*nu; 
  out[48] += 3.0*rdx[2]*f[48]*nu+90.0*rdx[1]*f[48]*nu+3.0*rdx[0]*f[48]*nu+20.12461179749811*rdx[1]*f[21]*nu; 
  out[49] += 90.0*rdx[2]*f[49]*nu+3.0*rdx[1]*f[49]*nu+3.0*rdx[0]*f[49]*nu+20.12461179749811*rdx[2]*f[22]*nu; 

return nu*(rdx[0]+rdx[1]+rdx[2])*0.5; 

} 
