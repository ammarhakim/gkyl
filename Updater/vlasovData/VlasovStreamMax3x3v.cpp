#include <VlasovModDecl.h> 
void VlasovVolStream3x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[4]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[5]*dv1dx1; 
  out[3] += 3.464101615137754*f[0]*w2dx2+f[6]*dv2dx2; 
} 
void VlasovVolStream3x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[4]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[5]*dv1dx1; 
  out[3] += 3.464101615137754*f[0]*w2dx2+f[6]*dv2dx2; 
  out[7] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[13]*dv1dx1+f[11]*dv0dx0; 
  out[8] += 3.464101615137754*f[1]*w2dx2+3.464101615137754*f[3]*w0dx0+f[17]*dv2dx2+f[12]*dv0dx0; 
  out[9] += 3.464101615137754*f[2]*w2dx2+3.464101615137754*f[3]*w1dx1+f[18]*dv2dx2+f[15]*dv1dx1; 
  out[10] += 3.464101615137754*f[4]*w0dx0+0.8944271909999159*f[25]*dv0dx0+f[0]*dv0dx0; 
  out[11] += 3.464101615137754*f[4]*w1dx1+f[16]*dv1dx1; 
  out[12] += 3.464101615137754*f[4]*w2dx2+f[20]*dv2dx2; 
  out[13] += 3.464101615137754*f[5]*w0dx0+f[16]*dv0dx0; 
  out[14] += 3.464101615137754*f[5]*w1dx1+0.8944271909999159*f[26]*dv1dx1+f[0]*dv1dx1; 
  out[15] += 3.464101615137754*f[5]*w2dx2+f[21]*dv2dx2; 
  out[17] += 3.464101615137754*f[6]*w0dx0+f[20]*dv0dx0; 
  out[18] += 3.464101615137754*f[6]*w1dx1+f[21]*dv1dx1; 
  out[19] += 3.464101615137754*f[6]*w2dx2+0.8944271909999159*f[27]*dv2dx2+f[0]*dv2dx2; 
  out[22] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[23] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[14]*dv1dx1; 
  out[24] += 7.745966692414834*f[3]*w2dx2+2.23606797749979*f[19]*dv2dx2; 
} 
void VlasovVolStream3x3vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[4]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[5]*dv1dx1; 
  out[3] += 3.464101615137754*f[0]*w2dx2+f[6]*dv2dx2; 
  out[7] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[13]*dv1dx1+f[11]*dv0dx0; 
  out[8] += 3.464101615137754*f[1]*w2dx2+3.464101615137754*f[3]*w0dx0+f[17]*dv2dx2+f[12]*dv0dx0; 
  out[9] += 3.464101615137754*f[2]*w2dx2+3.464101615137754*f[3]*w1dx1+f[18]*dv2dx2+f[15]*dv1dx1; 
  out[10] += 3.464101615137754*f[4]*w0dx0+0.8944271909999159*f[25]*dv0dx0+f[0]*dv0dx0; 
  out[11] += 3.464101615137754*f[4]*w1dx1+f[16]*dv1dx1; 
  out[12] += 3.464101615137754*f[4]*w2dx2+f[20]*dv2dx2; 
  out[13] += 3.464101615137754*f[5]*w0dx0+f[16]*dv0dx0; 
  out[14] += 3.464101615137754*f[5]*w1dx1+0.8944271909999159*f[26]*dv1dx1+f[0]*dv1dx1; 
  out[15] += 3.464101615137754*f[5]*w2dx2+f[21]*dv2dx2; 
  out[17] += 3.464101615137754*f[6]*w0dx0+f[20]*dv0dx0; 
  out[18] += 3.464101615137754*f[6]*w1dx1+f[21]*dv1dx1; 
  out[19] += 3.464101615137754*f[6]*w2dx2+0.8944271909999159*f[27]*dv2dx2+f[0]*dv2dx2; 
  out[22] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[23] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[14]*dv1dx1; 
  out[24] += 7.745966692414834*f[3]*w2dx2+2.23606797749979*f[19]*dv2dx2; 
  out[28] += 3.464101615137754*f[7]*w2dx2+3.464101615137754*f[8]*w1dx1+3.464101615137754*f[9]*w0dx0+f[38]*dv2dx2+f[33]*dv1dx1+f[31]*dv0dx0; 
  out[29] += 3.464101615137754*f[10]*w1dx1+3.464101615137754*f[11]*w0dx0+f[35]*dv1dx1+0.8944271909999161*f[58]*dv0dx0+f[2]*dv0dx0; 
  out[30] += 3.464101615137754*f[10]*w2dx2+3.464101615137754*f[12]*w0dx0+f[41]*dv2dx2+0.8944271909999161*f[59]*dv0dx0+f[3]*dv0dx0; 
  out[31] += 3.464101615137754*f[11]*w2dx2+3.464101615137754*f[12]*w1dx1+f[42]*dv2dx2+f[37]*dv1dx1; 
  out[32] += 3.464101615137754*f[13]*w1dx1+3.464101615137754*f[14]*w0dx0+0.8944271909999161*f[64]*dv1dx1+f[1]*dv1dx1+f[36]*dv0dx0; 
  out[33] += 3.464101615137754*f[13]*w2dx2+3.464101615137754*f[15]*w0dx0+f[44]*dv2dx2+f[37]*dv0dx0; 
  out[34] += 3.464101615137754*f[14]*w2dx2+3.464101615137754*f[15]*w1dx1+f[45]*dv2dx2+0.8944271909999161*f[66]*dv1dx1+f[3]*dv1dx1; 
  out[35] += 3.464101615137754*f[16]*w0dx0+0.8944271909999161*f[63]*dv0dx0+f[5]*dv0dx0; 
  out[36] += 3.464101615137754*f[16]*w1dx1+0.8944271909999161*f[67]*dv1dx1+f[4]*dv1dx1; 
  out[37] += 3.464101615137754*f[16]*w2dx2+f[47]*dv2dx2; 
  out[38] += 3.464101615137754*f[17]*w1dx1+3.464101615137754*f[18]*w0dx0+f[44]*dv1dx1+f[42]*dv0dx0; 
  out[39] += 3.464101615137754*f[17]*w2dx2+3.464101615137754*f[19]*w0dx0+0.8944271909999161*f[73]*dv2dx2+f[1]*dv2dx2+f[43]*dv0dx0; 
  out[40] += 3.464101615137754*f[18]*w2dx2+3.464101615137754*f[19]*w1dx1+0.8944271909999161*f[74]*dv2dx2+f[2]*dv2dx2+f[46]*dv1dx1; 
  out[41] += 3.464101615137754*f[20]*w0dx0+0.8944271909999161*f[71]*dv0dx0+f[6]*dv0dx0; 
  out[42] += 3.464101615137754*f[20]*w1dx1+f[47]*dv1dx1; 
  out[43] += 3.464101615137754*f[20]*w2dx2+0.8944271909999161*f[76]*dv2dx2+f[4]*dv2dx2; 
  out[44] += 3.464101615137754*f[21]*w0dx0+f[47]*dv0dx0; 
  out[45] += 3.464101615137754*f[21]*w1dx1+0.8944271909999161*f[72]*dv1dx1+f[6]*dv1dx1; 
  out[46] += 3.464101615137754*f[21]*w2dx2+0.8944271909999161*f[77]*dv2dx2+f[5]*dv2dx2; 
  out[48] += 3.464101615137755*f[22]*w1dx1+7.745966692414834*f[7]*w0dx0+f[60]*dv1dx1+2.23606797749979*f[29]*dv0dx0; 
  out[49] += 7.745966692414834*f[7]*w1dx1+3.464101615137755*f[23]*w0dx0+2.23606797749979*f[32]*dv1dx1+f[55]*dv0dx0; 
  out[50] += 3.464101615137755*f[22]*w2dx2+7.745966692414834*f[8]*w0dx0+f[68]*dv2dx2+2.23606797749979*f[30]*dv0dx0; 
  out[51] += 3.464101615137755*f[23]*w2dx2+7.745966692414834*f[9]*w1dx1+f[69]*dv2dx2+2.23606797749979*f[34]*dv1dx1; 
  out[52] += 7.745966692414834*f[8]*w2dx2+3.464101615137755*f[24]*w0dx0+2.23606797749979*f[39]*dv2dx2+f[56]*dv0dx0; 
  out[53] += 7.745966692414834*f[9]*w2dx2+3.464101615137755*f[24]*w1dx1+2.23606797749979*f[40]*dv2dx2+f[62]*dv1dx1; 
  out[54] += 7.745966692414834*f[10]*w0dx0+2.0*f[57]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[55] += 7.745966692414834*f[11]*w1dx1+2.23606797749979*f[36]*dv1dx1; 
  out[56] += 7.745966692414834*f[12]*w2dx2+2.23606797749979*f[43]*dv2dx2; 
  out[57] += 3.464101615137755*f[25]*w0dx0+0.8783100656536798*f[81]*dv0dx0+0.8944271909999161*f[4]*dv0dx0; 
  out[58] += 3.464101615137755*f[25]*w1dx1+f[63]*dv1dx1; 
  out[59] += 3.464101615137755*f[25]*w2dx2+f[71]*dv2dx2; 
  out[60] += 7.745966692414834*f[13]*w0dx0+2.23606797749979*f[35]*dv0dx0; 
  out[61] += 7.745966692414834*f[14]*w1dx1+2.0*f[65]*dv1dx1+2.23606797749979*f[2]*dv1dx1; 
  out[62] += 7.745966692414834*f[15]*w2dx2+2.23606797749979*f[46]*dv2dx2; 
  out[64] += 3.464101615137755*f[26]*w0dx0+f[67]*dv0dx0; 
  out[65] += 3.464101615137755*f[26]*w1dx1+0.8783100656536798*f[82]*dv1dx1+0.8944271909999161*f[5]*dv1dx1; 
  out[66] += 3.464101615137755*f[26]*w2dx2+f[72]*dv2dx2; 
  out[68] += 7.745966692414834*f[17]*w0dx0+2.23606797749979*f[41]*dv0dx0; 
  out[69] += 7.745966692414834*f[18]*w1dx1+2.23606797749979*f[45]*dv1dx1; 
  out[70] += 7.745966692414834*f[19]*w2dx2+2.0*f[75]*dv2dx2+2.23606797749979*f[3]*dv2dx2; 
  out[73] += 3.464101615137755*f[27]*w0dx0+f[76]*dv0dx0; 
  out[74] += 3.464101615137755*f[27]*w1dx1+f[77]*dv1dx1; 
  out[75] += 3.464101615137755*f[27]*w2dx2+0.8783100656536798*f[83]*dv2dx2+0.8944271909999161*f[6]*dv2dx2; 
  out[78] += 11.83215956619923*f[22]*w0dx0+5.291502622129181*f[0]*w0dx0+3.415650255319866*f[54]*dv0dx0+1.527525231651947*f[4]*dv0dx0; 
  out[79] += 11.83215956619923*f[23]*w1dx1+5.291502622129181*f[0]*w1dx1+3.415650255319866*f[61]*dv1dx1+1.527525231651947*f[5]*dv1dx1; 
  out[80] += 11.83215956619923*f[24]*w2dx2+5.291502622129181*f[0]*w2dx2+3.415650255319866*f[70]*dv2dx2+1.527525231651947*f[6]*dv2dx2; 
} 
void VlasovVolStream3x3vMaxP4(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[4]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[5]*dv1dx1; 
  out[3] += 3.464101615137754*f[0]*w2dx2+f[6]*dv2dx2; 
  out[7] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[13]*dv1dx1+f[11]*dv0dx0; 
  out[8] += 3.464101615137754*f[1]*w2dx2+3.464101615137754*f[3]*w0dx0+f[17]*dv2dx2+f[12]*dv0dx0; 
  out[9] += 3.464101615137754*f[2]*w2dx2+3.464101615137754*f[3]*w1dx1+f[18]*dv2dx2+f[15]*dv1dx1; 
  out[10] += 3.464101615137754*f[4]*w0dx0+0.8944271909999159*f[25]*dv0dx0+f[0]*dv0dx0; 
  out[11] += 3.464101615137754*f[4]*w1dx1+f[16]*dv1dx1; 
  out[12] += 3.464101615137754*f[4]*w2dx2+f[20]*dv2dx2; 
  out[13] += 3.464101615137754*f[5]*w0dx0+f[16]*dv0dx0; 
  out[14] += 3.464101615137754*f[5]*w1dx1+0.8944271909999159*f[26]*dv1dx1+f[0]*dv1dx1; 
  out[15] += 3.464101615137754*f[5]*w2dx2+f[21]*dv2dx2; 
  out[17] += 3.464101615137754*f[6]*w0dx0+f[20]*dv0dx0; 
  out[18] += 3.464101615137754*f[6]*w1dx1+f[21]*dv1dx1; 
  out[19] += 3.464101615137754*f[6]*w2dx2+0.8944271909999159*f[27]*dv2dx2+f[0]*dv2dx2; 
  out[22] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[23] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[14]*dv1dx1; 
  out[24] += 7.745966692414834*f[3]*w2dx2+2.23606797749979*f[19]*dv2dx2; 
  out[28] += 3.464101615137754*f[7]*w2dx2+3.464101615137754*f[8]*w1dx1+3.464101615137754*f[9]*w0dx0+f[38]*dv2dx2+f[33]*dv1dx1+f[31]*dv0dx0; 
  out[29] += 3.464101615137754*f[10]*w1dx1+3.464101615137754*f[11]*w0dx0+f[35]*dv1dx1+0.8944271909999161*f[58]*dv0dx0+f[2]*dv0dx0; 
  out[30] += 3.464101615137754*f[10]*w2dx2+3.464101615137754*f[12]*w0dx0+f[41]*dv2dx2+0.8944271909999161*f[59]*dv0dx0+f[3]*dv0dx0; 
  out[31] += 3.464101615137754*f[11]*w2dx2+3.464101615137754*f[12]*w1dx1+f[42]*dv2dx2+f[37]*dv1dx1; 
  out[32] += 3.464101615137754*f[13]*w1dx1+3.464101615137754*f[14]*w0dx0+0.8944271909999161*f[64]*dv1dx1+f[1]*dv1dx1+f[36]*dv0dx0; 
  out[33] += 3.464101615137754*f[13]*w2dx2+3.464101615137754*f[15]*w0dx0+f[44]*dv2dx2+f[37]*dv0dx0; 
  out[34] += 3.464101615137754*f[14]*w2dx2+3.464101615137754*f[15]*w1dx1+f[45]*dv2dx2+0.8944271909999161*f[66]*dv1dx1+f[3]*dv1dx1; 
  out[35] += 3.464101615137754*f[16]*w0dx0+0.8944271909999161*f[63]*dv0dx0+f[5]*dv0dx0; 
  out[36] += 3.464101615137754*f[16]*w1dx1+0.8944271909999161*f[67]*dv1dx1+f[4]*dv1dx1; 
  out[37] += 3.464101615137754*f[16]*w2dx2+f[47]*dv2dx2; 
  out[38] += 3.464101615137754*f[17]*w1dx1+3.464101615137754*f[18]*w0dx0+f[44]*dv1dx1+f[42]*dv0dx0; 
  out[39] += 3.464101615137754*f[17]*w2dx2+3.464101615137754*f[19]*w0dx0+0.8944271909999161*f[73]*dv2dx2+f[1]*dv2dx2+f[43]*dv0dx0; 
  out[40] += 3.464101615137754*f[18]*w2dx2+3.464101615137754*f[19]*w1dx1+0.8944271909999161*f[74]*dv2dx2+f[2]*dv2dx2+f[46]*dv1dx1; 
  out[41] += 3.464101615137754*f[20]*w0dx0+0.8944271909999161*f[71]*dv0dx0+f[6]*dv0dx0; 
  out[42] += 3.464101615137754*f[20]*w1dx1+f[47]*dv1dx1; 
  out[43] += 3.464101615137754*f[20]*w2dx2+0.8944271909999161*f[76]*dv2dx2+f[4]*dv2dx2; 
  out[44] += 3.464101615137754*f[21]*w0dx0+f[47]*dv0dx0; 
  out[45] += 3.464101615137754*f[21]*w1dx1+0.8944271909999161*f[72]*dv1dx1+f[6]*dv1dx1; 
  out[46] += 3.464101615137754*f[21]*w2dx2+0.8944271909999161*f[77]*dv2dx2+f[5]*dv2dx2; 
  out[48] += 3.464101615137755*f[22]*w1dx1+7.745966692414834*f[7]*w0dx0+f[60]*dv1dx1+2.23606797749979*f[29]*dv0dx0; 
  out[49] += 7.745966692414834*f[7]*w1dx1+3.464101615137755*f[23]*w0dx0+2.23606797749979*f[32]*dv1dx1+f[55]*dv0dx0; 
  out[50] += 3.464101615137755*f[22]*w2dx2+7.745966692414834*f[8]*w0dx0+f[68]*dv2dx2+2.23606797749979*f[30]*dv0dx0; 
  out[51] += 3.464101615137755*f[23]*w2dx2+7.745966692414834*f[9]*w1dx1+f[69]*dv2dx2+2.23606797749979*f[34]*dv1dx1; 
  out[52] += 7.745966692414834*f[8]*w2dx2+3.464101615137755*f[24]*w0dx0+2.23606797749979*f[39]*dv2dx2+f[56]*dv0dx0; 
  out[53] += 7.745966692414834*f[9]*w2dx2+3.464101615137755*f[24]*w1dx1+2.23606797749979*f[40]*dv2dx2+f[62]*dv1dx1; 
  out[54] += 7.745966692414834*f[10]*w0dx0+2.0*f[57]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[55] += 7.745966692414834*f[11]*w1dx1+2.23606797749979*f[36]*dv1dx1; 
  out[56] += 7.745966692414834*f[12]*w2dx2+2.23606797749979*f[43]*dv2dx2; 
  out[57] += 3.464101615137755*f[25]*w0dx0+0.8783100656536798*f[81]*dv0dx0+0.8944271909999161*f[4]*dv0dx0; 
  out[58] += 3.464101615137755*f[25]*w1dx1+f[63]*dv1dx1; 
  out[59] += 3.464101615137755*f[25]*w2dx2+f[71]*dv2dx2; 
  out[60] += 7.745966692414834*f[13]*w0dx0+2.23606797749979*f[35]*dv0dx0; 
  out[61] += 7.745966692414834*f[14]*w1dx1+2.0*f[65]*dv1dx1+2.23606797749979*f[2]*dv1dx1; 
  out[62] += 7.745966692414834*f[15]*w2dx2+2.23606797749979*f[46]*dv2dx2; 
  out[64] += 3.464101615137755*f[26]*w0dx0+f[67]*dv0dx0; 
  out[65] += 3.464101615137755*f[26]*w1dx1+0.8783100656536798*f[82]*dv1dx1+0.8944271909999161*f[5]*dv1dx1; 
  out[66] += 3.464101615137755*f[26]*w2dx2+f[72]*dv2dx2; 
  out[68] += 7.745966692414834*f[17]*w0dx0+2.23606797749979*f[41]*dv0dx0; 
  out[69] += 7.745966692414834*f[18]*w1dx1+2.23606797749979*f[45]*dv1dx1; 
  out[70] += 7.745966692414834*f[19]*w2dx2+2.0*f[75]*dv2dx2+2.23606797749979*f[3]*dv2dx2; 
  out[73] += 3.464101615137755*f[27]*w0dx0+f[76]*dv0dx0; 
  out[74] += 3.464101615137755*f[27]*w1dx1+f[77]*dv1dx1; 
  out[75] += 3.464101615137755*f[27]*w2dx2+0.8783100656536798*f[83]*dv2dx2+0.8944271909999161*f[6]*dv2dx2; 
  out[78] += 11.83215956619923*f[22]*w0dx0+5.291502622129181*f[0]*w0dx0+3.415650255319866*f[54]*dv0dx0+1.527525231651947*f[4]*dv0dx0; 
  out[79] += 11.83215956619923*f[23]*w1dx1+5.291502622129181*f[0]*w1dx1+3.415650255319866*f[61]*dv1dx1+1.527525231651947*f[5]*dv1dx1; 
  out[80] += 11.83215956619923*f[24]*w2dx2+5.291502622129181*f[0]*w2dx2+3.415650255319866*f[70]*dv2dx2+1.527525231651947*f[6]*dv2dx2; 
  out[84] += 3.464101615137754*f[29]*w2dx2+3.464101615137754*f[30]*w1dx1+3.464101615137754*f[31]*w0dx0+f[90]*dv2dx2+f[87]*dv1dx1+0.8944271909999159*f[110]*dv0dx0+f[9]*dv0dx0; 
  out[85] += 3.464101615137754*f[32]*w2dx2+3.464101615137754*f[33]*w1dx1+3.464101615137754*f[34]*w0dx0+f[93]*dv2dx2+0.8944271909999159*f[124]*dv1dx1+f[8]*dv1dx1+f[88]*dv0dx0; 
  out[86] += 3.464101615137754*f[35]*w1dx1+3.464101615137754*f[36]*w0dx0+0.8944271909999159*f[126]*dv1dx1+f[10]*dv1dx1+0.8944271909999159*f[121]*dv0dx0+f[14]*dv0dx0; 
  out[87] += 3.464101615137754*f[35]*w2dx2+3.464101615137754*f[37]*w0dx0+f[96]*dv2dx2+0.8944271909999159*f[122]*dv0dx0+f[15]*dv0dx0; 
  out[88] += 3.464101615137754*f[36]*w2dx2+3.464101615137754*f[37]*w1dx1+f[97]*dv2dx2+0.8944271909999159*f[128]*dv1dx1+f[12]*dv1dx1; 
  out[89] += 3.464101615137754*f[38]*w2dx2+3.464101615137754*f[39]*w1dx1+3.464101615137754*f[40]*w0dx0+0.8944271909999159*f[149]*dv2dx2+f[7]*dv2dx2+f[94]*dv1dx1+f[92]*dv0dx0; 
  out[90] += 3.464101615137754*f[41]*w1dx1+3.464101615137754*f[42]*w0dx0+f[96]*dv1dx1+0.8944271909999159*f[139]*dv0dx0+f[18]*dv0dx0; 
  out[91] += 3.464101615137754*f[41]*w2dx2+3.464101615137754*f[43]*w0dx0+0.8944271909999159*f[152]*dv2dx2+f[10]*dv2dx2+0.8944271909999159*f[140]*dv0dx0+f[19]*dv0dx0; 
  out[92] += 3.464101615137754*f[42]*w2dx2+3.464101615137754*f[43]*w1dx1+0.8944271909999159*f[153]*dv2dx2+f[11]*dv2dx2+f[98]*dv1dx1; 
  out[93] += 3.464101615137754*f[44]*w1dx1+3.464101615137754*f[45]*w0dx0+0.8944271909999159*f[145]*dv1dx1+f[17]*dv1dx1+f[97]*dv0dx0; 
  out[94] += 3.464101615137754*f[44]*w2dx2+3.464101615137754*f[46]*w0dx0+0.8944271909999159*f[155]*dv2dx2+f[13]*dv2dx2+f[98]*dv0dx0; 
  out[95] += 3.464101615137754*f[45]*w2dx2+3.464101615137754*f[46]*w1dx1+0.8944271909999159*f[156]*dv2dx2+f[14]*dv2dx2+0.8944271909999159*f[147]*dv1dx1+f[19]*dv1dx1; 
  out[96] += 3.464101615137754*f[47]*w0dx0+0.8944271909999159*f[144]*dv0dx0+f[21]*dv0dx0; 
  out[97] += 3.464101615137754*f[47]*w1dx1+0.8944271909999159*f[148]*dv1dx1+f[20]*dv1dx1; 
  out[98] += 3.464101615137754*f[47]*w2dx2+0.8944271909999159*f[158]*dv2dx2+f[16]*dv2dx2; 
  out[99] += 3.464101615137755*f[48]*w2dx2+3.464101615137755*f[50]*w1dx1+7.745966692414834*f[28]*w0dx0+f[129]*dv2dx2+f[113]*dv1dx1+2.23606797749979*f[84]*dv0dx0; 
  out[100] += 3.464101615137755*f[49]*w2dx2+7.745966692414834*f[28]*w1dx1+3.464101615137755*f[51]*w0dx0+f[130]*dv2dx2+2.23606797749979*f[85]*dv1dx1+f[105]*dv0dx0; 
  out[101] += 7.745966692414834*f[28]*w2dx2+3.464101615137755*f[52]*w1dx1+3.464101615137755*f[53]*w0dx0+2.23606797749979*f[89]*dv2dx2+f[115]*dv1dx1+f[107]*dv0dx0; 
  out[102] += 3.464101615137755*f[54]*w1dx1+7.745966692414834*f[29]*w0dx0+f[117]*dv1dx1+2.0*f[108]*dv0dx0+2.23606797749979*f[7]*dv0dx0; 
  out[103] += 7.745966692414834*f[29]*w1dx1+3.464101615137755*f[55]*w0dx0+2.23606797749979*f[86]*dv1dx1+0.8944271909999159*f[163]*dv0dx0+f[23]*dv0dx0; 
  out[104] += 3.464101615137755*f[54]*w2dx2+7.745966692414834*f[30]*w0dx0+f[135]*dv2dx2+2.0*f[109]*dv0dx0+2.23606797749979*f[8]*dv0dx0; 
  out[105] += 3.464101615137755*f[55]*w2dx2+7.745966692414834*f[31]*w1dx1+f[136]*dv2dx2+2.23606797749979*f[88]*dv1dx1; 
  out[106] += 7.745966692414834*f[30]*w2dx2+3.464101615137755*f[56]*w0dx0+2.23606797749979*f[91]*dv2dx2+0.8944271909999159*f[164]*dv0dx0+f[24]*dv0dx0; 
  out[107] += 7.745966692414834*f[31]*w2dx2+3.464101615137755*f[56]*w1dx1+2.23606797749979*f[92]*dv2dx2+f[119]*dv1dx1; 
  out[108] += 3.464101615137755*f[57]*w1dx1+3.464101615137755*f[58]*w0dx0+f[120]*dv1dx1+0.8783100656536798*f[184]*dv0dx0+0.8944271909999159*f[11]*dv0dx0; 
  out[109] += 3.464101615137755*f[57]*w2dx2+3.464101615137755*f[59]*w0dx0+f[138]*dv2dx2+0.8783100656536798*f[185]*dv0dx0+0.8944271909999159*f[12]*dv0dx0; 
  out[110] += 3.464101615137755*f[58]*w2dx2+3.464101615137755*f[59]*w1dx1+f[139]*dv2dx2+f[122]*dv1dx1; 
  out[111] += 3.464101615137755*f[60]*w1dx1+7.745966692414834*f[32]*w0dx0+0.8944271909999159*f[165]*dv1dx1+f[22]*dv1dx1+2.23606797749979*f[86]*dv0dx0; 
  out[112] += 7.745966692414834*f[32]*w1dx1+3.464101615137755*f[61]*w0dx0+2.0*f[123]*dv1dx1+2.23606797749979*f[7]*dv1dx1+f[118]*dv0dx0; 
  out[113] += 3.464101615137755*f[60]*w2dx2+7.745966692414834*f[33]*w0dx0+f[141]*dv2dx2+2.23606797749979*f[87]*dv0dx0; 
  out[114] += 3.464101615137755*f[61]*w2dx2+7.745966692414834*f[34]*w1dx1+f[142]*dv2dx2+2.0*f[125]*dv1dx1+2.23606797749979*f[9]*dv1dx1; 
  out[115] += 7.745966692414834*f[33]*w2dx2+3.464101615137755*f[62]*w0dx0+2.23606797749979*f[94]*dv2dx2+f[119]*dv0dx0; 
  out[116] += 7.745966692414834*f[34]*w2dx2+3.464101615137755*f[62]*w1dx1+2.23606797749979*f[95]*dv2dx2+0.8944271909999159*f[167]*dv1dx1+f[24]*dv1dx1; 
  out[117] += 7.745966692414834*f[35]*w0dx0+2.0*f[120]*dv0dx0+2.23606797749979*f[13]*dv0dx0; 
  out[118] += 7.745966692414834*f[36]*w1dx1+2.0*f[127]*dv1dx1+2.23606797749979*f[11]*dv1dx1; 
  out[119] += 7.745966692414834*f[37]*w2dx2+2.23606797749979*f[98]*dv2dx2; 
  out[120] += 3.464101615137755*f[63]*w0dx0+0.8783100656536798*f[189]*dv0dx0+0.8944271909999159*f[16]*dv0dx0; 
  out[121] += 3.464101615137755*f[63]*w1dx1+0.8944271909999159*f[168]*dv1dx1+f[25]*dv1dx1; 
  out[122] += 3.464101615137755*f[63]*w2dx2+f[144]*dv2dx2; 
  out[123] += 3.464101615137755*f[64]*w1dx1+3.464101615137755*f[65]*w0dx0+0.8783100656536798*f[190]*dv1dx1+0.8944271909999159*f[13]*dv1dx1+f[127]*dv0dx0; 
  out[124] += 3.464101615137755*f[64]*w2dx2+3.464101615137755*f[66]*w0dx0+f[145]*dv2dx2+f[128]*dv0dx0; 
  out[125] += 3.464101615137755*f[65]*w2dx2+3.464101615137755*f[66]*w1dx1+f[146]*dv2dx2+0.8783100656536798*f[192]*dv1dx1+0.8944271909999159*f[15]*dv1dx1; 
  out[126] += 3.464101615137755*f[67]*w0dx0+0.8944271909999159*f[168]*dv0dx0+f[26]*dv0dx0; 
  out[127] += 3.464101615137755*f[67]*w1dx1+0.8783100656536798*f[193]*dv1dx1+0.8944271909999159*f[16]*dv1dx1; 
  out[128] += 3.464101615137755*f[67]*w2dx2+f[148]*dv2dx2; 
  out[129] += 3.464101615137755*f[68]*w1dx1+7.745966692414834*f[38]*w0dx0+f[141]*dv1dx1+2.23606797749979*f[90]*dv0dx0; 
  out[130] += 7.745966692414834*f[38]*w1dx1+3.464101615137755*f[69]*w0dx0+2.23606797749979*f[93]*dv1dx1+f[136]*dv0dx0; 
  out[131] += 3.464101615137755*f[68]*w2dx2+7.745966692414834*f[39]*w0dx0+0.8944271909999159*f[169]*dv2dx2+f[22]*dv2dx2+2.23606797749979*f[91]*dv0dx0; 
  out[132] += 3.464101615137755*f[69]*w2dx2+7.745966692414834*f[40]*w1dx1+0.8944271909999159*f[170]*dv2dx2+f[23]*dv2dx2+2.23606797749979*f[95]*dv1dx1; 
  out[133] += 7.745966692414834*f[39]*w2dx2+3.464101615137755*f[70]*w0dx0+2.0*f[150]*dv2dx2+2.23606797749979*f[8]*dv2dx2+f[137]*dv0dx0; 
  out[134] += 7.745966692414834*f[40]*w2dx2+3.464101615137755*f[70]*w1dx1+2.0*f[151]*dv2dx2+2.23606797749979*f[9]*dv2dx2+f[143]*dv1dx1; 
  out[135] += 7.745966692414834*f[41]*w0dx0+2.0*f[138]*dv0dx0+2.23606797749979*f[17]*dv0dx0; 
  out[136] += 7.745966692414834*f[42]*w1dx1+2.23606797749979*f[97]*dv1dx1; 
  out[137] += 7.745966692414834*f[43]*w2dx2+2.0*f[154]*dv2dx2+2.23606797749979*f[12]*dv2dx2; 
  out[138] += 3.464101615137755*f[71]*w0dx0+0.8783100656536798*f[197]*dv0dx0+0.8944271909999159*f[20]*dv0dx0; 
  out[139] += 3.464101615137755*f[71]*w1dx1+f[144]*dv1dx1; 
  out[140] += 3.464101615137755*f[71]*w2dx2+0.8944271909999159*f[172]*dv2dx2+f[25]*dv2dx2; 
  out[141] += 7.745966692414834*f[44]*w0dx0+2.23606797749979*f[96]*dv0dx0; 
  out[142] += 7.745966692414834*f[45]*w1dx1+2.0*f[146]*dv1dx1+2.23606797749979*f[18]*dv1dx1; 
  out[143] += 7.745966692414834*f[46]*w2dx2+2.0*f[157]*dv2dx2+2.23606797749979*f[15]*dv2dx2; 
  out[145] += 3.464101615137755*f[72]*w0dx0+f[148]*dv0dx0; 
  out[146] += 3.464101615137755*f[72]*w1dx1+0.8783100656536798*f[198]*dv1dx1+0.8944271909999159*f[21]*dv1dx1; 
  out[147] += 3.464101615137755*f[72]*w2dx2+0.8944271909999159*f[173]*dv2dx2+f[26]*dv2dx2; 
  out[149] += 3.464101615137755*f[73]*w1dx1+3.464101615137755*f[74]*w0dx0+f[155]*dv1dx1+f[153]*dv0dx0; 
  out[150] += 3.464101615137755*f[73]*w2dx2+3.464101615137755*f[75]*w0dx0+0.8783100656536798*f[199]*dv2dx2+0.8944271909999159*f[17]*dv2dx2+f[154]*dv0dx0; 
  out[151] += 3.464101615137755*f[74]*w2dx2+3.464101615137755*f[75]*w1dx1+0.8783100656536798*f[200]*dv2dx2+0.8944271909999159*f[18]*dv2dx2+f[157]*dv1dx1; 
  out[152] += 3.464101615137755*f[76]*w0dx0+0.8944271909999159*f[172]*dv0dx0+f[27]*dv0dx0; 
  out[153] += 3.464101615137755*f[76]*w1dx1+f[158]*dv1dx1; 
  out[154] += 3.464101615137755*f[76]*w2dx2+0.8783100656536798*f[202]*dv2dx2+0.8944271909999159*f[20]*dv2dx2; 
  out[155] += 3.464101615137755*f[77]*w0dx0+f[158]*dv0dx0; 
  out[156] += 3.464101615137755*f[77]*w1dx1+0.8944271909999159*f[173]*dv1dx1+f[27]*dv1dx1; 
  out[157] += 3.464101615137755*f[77]*w2dx2+0.8783100656536798*f[203]*dv2dx2+0.8944271909999159*f[21]*dv2dx2; 
  out[159] += 7.745966692414834*f[48]*w1dx1+7.745966692414834*f[49]*w0dx0+2.23606797749979*f[111]*dv1dx1+2.23606797749979*f[103]*dv0dx0; 
  out[160] += 7.745966692414834*f[50]*w2dx2+7.745966692414834*f[52]*w0dx0+2.23606797749979*f[131]*dv2dx2+2.23606797749979*f[106]*dv0dx0; 
  out[161] += 7.745966692414834*f[51]*w2dx2+7.745966692414834*f[53]*w1dx1+2.23606797749979*f[132]*dv2dx2+2.23606797749979*f[116]*dv1dx1; 
  out[162] += 7.745966692414834*f[57]*w0dx0+1.963961012123931*f[183]*dv0dx0+2.0*f[10]*dv0dx0; 
  out[163] += 7.745966692414834*f[58]*w1dx1+2.23606797749979*f[121]*dv1dx1; 
  out[164] += 7.745966692414834*f[59]*w2dx2+2.23606797749979*f[140]*dv2dx2; 
  out[165] += 7.745966692414834*f[64]*w0dx0+2.23606797749979*f[126]*dv0dx0; 
  out[166] += 7.745966692414834*f[65]*w1dx1+1.963961012123931*f[191]*dv1dx1+2.0*f[14]*dv1dx1; 
  out[167] += 7.745966692414834*f[66]*w2dx2+2.23606797749979*f[147]*dv2dx2; 
  out[169] += 7.745966692414834*f[73]*w0dx0+2.23606797749979*f[152]*dv0dx0; 
  out[170] += 7.745966692414834*f[74]*w1dx1+2.23606797749979*f[156]*dv1dx1; 
  out[171] += 7.745966692414834*f[75]*w2dx2+1.963961012123931*f[201]*dv2dx2+2.0*f[19]*dv2dx2; 
  out[174] += 3.464101615137754*f[78]*w1dx1+11.83215956619923*f[48]*w0dx0+5.291502622129181*f[2]*w0dx0+f[186]*dv1dx1+3.415650255319866*f[102]*dv0dx0+1.527525231651947*f[11]*dv0dx0; 
  out[175] += 11.83215956619923*f[49]*w1dx1+5.291502622129181*f[1]*w1dx1+3.464101615137754*f[79]*w0dx0+3.415650255319866*f[112]*dv1dx1+1.527525231651947*f[13]*dv1dx1+f[181]*dv0dx0; 
  out[176] += 3.464101615137754*f[78]*w2dx2+11.83215956619923*f[50]*w0dx0+5.291502622129181*f[3]*w0dx0+f[194]*dv2dx2+3.415650255319866*f[104]*dv0dx0+1.527525231651947*f[12]*dv0dx0; 
  out[177] += 3.464101615137754*f[79]*w2dx2+11.83215956619923*f[51]*w1dx1+5.291502622129181*f[3]*w1dx1+f[195]*dv2dx2+3.415650255319866*f[114]*dv1dx1+1.527525231651947*f[15]*dv1dx1; 
  out[178] += 11.83215956619923*f[52]*w2dx2+5.291502622129181*f[1]*w2dx2+3.464101615137754*f[80]*w0dx0+3.415650255319866*f[133]*dv2dx2+1.527525231651947*f[17]*dv2dx2+f[182]*dv0dx0; 
  out[179] += 11.83215956619923*f[53]*w2dx2+5.291502622129181*f[2]*w2dx2+3.464101615137754*f[80]*w1dx1+3.415650255319866*f[134]*dv2dx2+1.527525231651947*f[18]*dv2dx2+f[188]*dv1dx1; 
  out[180] += 11.83215956619923*f[54]*w0dx0+5.291502622129181*f[4]*w0dx0+3.055050463303893*f[162]*dv0dx0+1.366260102127946*f[25]*dv0dx0+3.415650255319866*f[22]*dv0dx0+1.527525231651947*f[0]*dv0dx0; 
  out[181] += 11.83215956619923*f[55]*w1dx1+5.291502622129181*f[4]*w1dx1+3.415650255319866*f[118]*dv1dx1+1.527525231651947*f[16]*dv1dx1; 
  out[182] += 11.83215956619923*f[56]*w2dx2+5.291502622129181*f[4]*w2dx2+3.415650255319866*f[137]*dv2dx2+1.527525231651947*f[20]*dv2dx2; 
  out[183] += 3.464101615137754*f[81]*w0dx0+0.8728715609439696*f[207]*dv0dx0+0.8783100656536798*f[25]*dv0dx0; 
  out[184] += 3.464101615137754*f[81]*w1dx1+f[189]*dv1dx1; 
  out[185] += 3.464101615137754*f[81]*w2dx2+f[197]*dv2dx2; 
  out[186] += 11.83215956619923*f[60]*w0dx0+5.291502622129181*f[5]*w0dx0+3.415650255319866*f[117]*dv0dx0+1.527525231651947*f[16]*dv0dx0; 
  out[187] += 11.83215956619923*f[61]*w1dx1+5.291502622129181*f[5]*w1dx1+3.055050463303893*f[166]*dv1dx1+1.366260102127946*f[26]*dv1dx1+3.415650255319866*f[23]*dv1dx1+1.527525231651947*f[0]*dv1dx1; 
  out[188] += 11.83215956619923*f[62]*w2dx2+5.291502622129181*f[5]*w2dx2+3.415650255319866*f[143]*dv2dx2+1.527525231651947*f[21]*dv2dx2; 
  out[190] += 3.464101615137754*f[82]*w0dx0+f[193]*dv0dx0; 
  out[191] += 3.464101615137754*f[82]*w1dx1+0.8728715609439696*f[208]*dv1dx1+0.8783100656536798*f[26]*dv1dx1; 
  out[192] += 3.464101615137754*f[82]*w2dx2+f[198]*dv2dx2; 
  out[194] += 11.83215956619923*f[68]*w0dx0+5.291502622129181*f[6]*w0dx0+3.415650255319866*f[135]*dv0dx0+1.527525231651947*f[20]*dv0dx0; 
  out[195] += 11.83215956619923*f[69]*w1dx1+5.291502622129181*f[6]*w1dx1+3.415650255319866*f[142]*dv1dx1+1.527525231651947*f[21]*dv1dx1; 
  out[196] += 11.83215956619923*f[70]*w2dx2+5.291502622129181*f[6]*w2dx2+3.055050463303893*f[171]*dv2dx2+1.366260102127946*f[27]*dv2dx2+3.415650255319866*f[24]*dv2dx2+1.527525231651947*f[0]*dv2dx2; 
  out[199] += 3.464101615137754*f[83]*w0dx0+f[202]*dv0dx0; 
  out[200] += 3.464101615137754*f[83]*w1dx1+f[203]*dv1dx1; 
  out[201] += 3.464101615137754*f[83]*w2dx2+0.8728715609439696*f[209]*dv2dx2+0.8783100656536798*f[27]*dv2dx2; 
  out[204] += 15.87450786638754*f[78]*w0dx0+10.39230484541326*f[1]*w0dx0+4.58257569495584*f[180]*dv0dx0+3.0*f[10]*dv0dx0; 
  out[205] += 15.87450786638754*f[79]*w1dx1+10.39230484541326*f[2]*w1dx1+4.58257569495584*f[187]*dv1dx1+3.0*f[14]*dv1dx1; 
  out[206] += 15.87450786638754*f[80]*w2dx2+10.39230484541326*f[3]*w2dx2+4.58257569495584*f[196]*dv2dx2+3.0*f[19]*dv2dx2; 
} 
