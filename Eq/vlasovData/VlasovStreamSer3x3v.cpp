#include <VlasovModDecl.h> 
double VlasovVolStream3x3vSerP1(const double *w, const double *dxv, const double *f, double *out) 
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
  out[10] += 3.464101615137754*f[4]*w0dx0+f[0]*dv0dx0; 
  out[11] += 3.464101615137754*f[4]*w1dx1+f[16]*dv1dx1; 
  out[12] += 3.464101615137754*f[4]*w2dx2+f[20]*dv2dx2; 
  out[13] += 3.464101615137754*f[5]*w0dx0+f[16]*dv0dx0; 
  out[14] += 3.464101615137754*f[5]*w1dx1+f[0]*dv1dx1; 
  out[15] += 3.464101615137754*f[5]*w2dx2+f[21]*dv2dx2; 
  out[17] += 3.464101615137754*f[6]*w0dx0+f[20]*dv0dx0; 
  out[18] += 3.464101615137754*f[6]*w1dx1+f[21]*dv1dx1; 
  out[19] += 3.464101615137754*f[6]*w2dx2+f[0]*dv2dx2; 
  out[22] += 3.464101615137754*f[7]*w2dx2+3.464101615137754*f[8]*w1dx1+3.464101615137754*f[9]*w0dx0+f[32]*dv2dx2+f[27]*dv1dx1+f[25]*dv0dx0; 
  out[23] += 3.464101615137754*f[10]*w1dx1+3.464101615137754*f[11]*w0dx0+f[29]*dv1dx1+f[2]*dv0dx0; 
  out[24] += 3.464101615137754*f[10]*w2dx2+3.464101615137754*f[12]*w0dx0+f[35]*dv2dx2+f[3]*dv0dx0; 
  out[25] += 3.464101615137754*f[11]*w2dx2+3.464101615137754*f[12]*w1dx1+f[36]*dv2dx2+f[31]*dv1dx1; 
  out[26] += 3.464101615137754*f[13]*w1dx1+3.464101615137754*f[14]*w0dx0+f[1]*dv1dx1+f[30]*dv0dx0; 
  out[27] += 3.464101615137754*f[13]*w2dx2+3.464101615137754*f[15]*w0dx0+f[38]*dv2dx2+f[31]*dv0dx0; 
  out[28] += 3.464101615137754*f[14]*w2dx2+3.464101615137754*f[15]*w1dx1+f[39]*dv2dx2+f[3]*dv1dx1; 
  out[29] += 3.464101615137754*f[16]*w0dx0+f[5]*dv0dx0; 
  out[30] += 3.464101615137754*f[16]*w1dx1+f[4]*dv1dx1; 
  out[31] += 3.464101615137754*f[16]*w2dx2+f[41]*dv2dx2; 
  out[32] += 3.464101615137754*f[17]*w1dx1+3.464101615137754*f[18]*w0dx0+f[38]*dv1dx1+f[36]*dv0dx0; 
  out[33] += 3.464101615137754*f[17]*w2dx2+3.464101615137754*f[19]*w0dx0+f[1]*dv2dx2+f[37]*dv0dx0; 
  out[34] += 3.464101615137754*f[18]*w2dx2+3.464101615137754*f[19]*w1dx1+f[2]*dv2dx2+f[40]*dv1dx1; 
  out[35] += 3.464101615137754*f[20]*w0dx0+f[6]*dv0dx0; 
  out[36] += 3.464101615137754*f[20]*w1dx1+f[41]*dv1dx1; 
  out[37] += 3.464101615137754*f[20]*w2dx2+f[4]*dv2dx2; 
  out[38] += 3.464101615137754*f[21]*w0dx0+f[41]*dv0dx0; 
  out[39] += 3.464101615137754*f[21]*w1dx1+f[6]*dv1dx1; 
  out[40] += 3.464101615137754*f[21]*w2dx2+f[5]*dv2dx2; 
  out[42] += 3.464101615137754*f[23]*w2dx2+3.464101615137754*f[24]*w1dx1+3.464101615137754*f[25]*w0dx0+f[48]*dv2dx2+f[45]*dv1dx1+f[9]*dv0dx0; 
  out[43] += 3.464101615137754*f[26]*w2dx2+3.464101615137754*f[27]*w1dx1+3.464101615137754*f[28]*w0dx0+f[51]*dv2dx2+f[8]*dv1dx1+f[46]*dv0dx0; 
  out[44] += 3.464101615137754*f[29]*w1dx1+3.464101615137754*f[30]*w0dx0+f[10]*dv1dx1+f[14]*dv0dx0; 
  out[45] += 3.464101615137754*f[29]*w2dx2+3.464101615137754*f[31]*w0dx0+f[54]*dv2dx2+f[15]*dv0dx0; 
  out[46] += 3.464101615137754*f[30]*w2dx2+3.464101615137754*f[31]*w1dx1+f[55]*dv2dx2+f[12]*dv1dx1; 
  out[47] += 3.464101615137754*f[32]*w2dx2+3.464101615137754*f[33]*w1dx1+3.464101615137754*f[34]*w0dx0+f[7]*dv2dx2+f[52]*dv1dx1+f[50]*dv0dx0; 
  out[48] += 3.464101615137754*f[35]*w1dx1+3.464101615137754*f[36]*w0dx0+f[54]*dv1dx1+f[18]*dv0dx0; 
  out[49] += 3.464101615137754*f[35]*w2dx2+3.464101615137754*f[37]*w0dx0+f[10]*dv2dx2+f[19]*dv0dx0; 
  out[50] += 3.464101615137754*f[36]*w2dx2+3.464101615137754*f[37]*w1dx1+f[11]*dv2dx2+f[56]*dv1dx1; 
  out[51] += 3.464101615137754*f[38]*w1dx1+3.464101615137754*f[39]*w0dx0+f[17]*dv1dx1+f[55]*dv0dx0; 
  out[52] += 3.464101615137754*f[38]*w2dx2+3.464101615137754*f[40]*w0dx0+f[13]*dv2dx2+f[56]*dv0dx0; 
  out[53] += 3.464101615137754*f[39]*w2dx2+3.464101615137754*f[40]*w1dx1+f[14]*dv2dx2+f[19]*dv1dx1; 
  out[54] += 3.464101615137754*f[41]*w0dx0+f[21]*dv0dx0; 
  out[55] += 3.464101615137754*f[41]*w1dx1+f[20]*dv1dx1; 
  out[56] += 3.464101615137754*f[41]*w2dx2+f[16]*dv2dx2; 
  out[57] += 3.464101615137754*f[44]*w2dx2+3.464101615137754*f[45]*w1dx1+3.464101615137754*f[46]*w0dx0+f[60]*dv2dx2+f[24]*dv1dx1+f[28]*dv0dx0; 
  out[58] += 3.464101615137754*f[48]*w2dx2+3.464101615137754*f[49]*w1dx1+3.464101615137754*f[50]*w0dx0+f[23]*dv2dx2+f[61]*dv1dx1+f[34]*dv0dx0; 
  out[59] += 3.464101615137754*f[51]*w2dx2+3.464101615137754*f[52]*w1dx1+3.464101615137754*f[53]*w0dx0+f[26]*dv2dx2+f[33]*dv1dx1+f[62]*dv0dx0; 
  out[60] += 3.464101615137754*f[54]*w1dx1+3.464101615137754*f[55]*w0dx0+f[35]*dv1dx1+f[39]*dv0dx0; 
  out[61] += 3.464101615137754*f[54]*w2dx2+3.464101615137754*f[56]*w0dx0+f[29]*dv2dx2+f[40]*dv0dx0; 
  out[62] += 3.464101615137754*f[55]*w2dx2+3.464101615137754*f[56]*w1dx1+f[30]*dv2dx2+f[37]*dv1dx1; 
  out[63] += 3.464101615137754*f[60]*w2dx2+3.464101615137754*f[61]*w1dx1+3.464101615137754*f[62]*w0dx0+f[44]*dv2dx2+f[49]*dv1dx1+f[53]*dv0dx0; 
return std::abs(w0dx0)+std::abs(w1dx1)+std::abs(w2dx2); 
} 
double VlasovVolStream3x3vSerP2(const double *w, const double *dxv, const double *f, double *out) 
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
  out[57] += 3.464101615137755*f[25]*w0dx0+0.8944271909999161*f[4]*dv0dx0; 
  out[58] += 3.464101615137755*f[25]*w1dx1+f[63]*dv1dx1; 
  out[59] += 3.464101615137755*f[25]*w2dx2+f[71]*dv2dx2; 
  out[60] += 7.745966692414834*f[13]*w0dx0+2.23606797749979*f[35]*dv0dx0; 
  out[61] += 7.745966692414834*f[14]*w1dx1+2.0*f[65]*dv1dx1+2.23606797749979*f[2]*dv1dx1; 
  out[62] += 7.745966692414834*f[15]*w2dx2+2.23606797749979*f[46]*dv2dx2; 
  out[64] += 3.464101615137755*f[26]*w0dx0+f[67]*dv0dx0; 
  out[65] += 3.464101615137755*f[26]*w1dx1+0.8944271909999161*f[5]*dv1dx1; 
  out[66] += 3.464101615137755*f[26]*w2dx2+f[72]*dv2dx2; 
  out[68] += 7.745966692414834*f[17]*w0dx0+2.23606797749979*f[41]*dv0dx0; 
  out[69] += 7.745966692414834*f[18]*w1dx1+2.23606797749979*f[45]*dv1dx1; 
  out[70] += 7.745966692414834*f[19]*w2dx2+2.0*f[75]*dv2dx2+2.23606797749979*f[3]*dv2dx2; 
  out[73] += 3.464101615137755*f[27]*w0dx0+f[76]*dv0dx0; 
  out[74] += 3.464101615137755*f[27]*w1dx1+f[77]*dv1dx1; 
  out[75] += 3.464101615137755*f[27]*w2dx2+0.8944271909999161*f[6]*dv2dx2; 
  out[78] += 3.464101615137754*f[29]*w2dx2+3.464101615137754*f[30]*w1dx1+3.464101615137754*f[31]*w0dx0+f[84]*dv2dx2+f[81]*dv1dx1+0.8944271909999159*f[104]*dv0dx0+f[9]*dv0dx0; 
  out[79] += 3.464101615137754*f[32]*w2dx2+3.464101615137754*f[33]*w1dx1+3.464101615137754*f[34]*w0dx0+f[87]*dv2dx2+0.8944271909999159*f[118]*dv1dx1+f[8]*dv1dx1+f[82]*dv0dx0; 
  out[80] += 3.464101615137754*f[35]*w1dx1+3.464101615137754*f[36]*w0dx0+0.8944271909999159*f[120]*dv1dx1+f[10]*dv1dx1+0.8944271909999159*f[115]*dv0dx0+f[14]*dv0dx0; 
  out[81] += 3.464101615137754*f[35]*w2dx2+3.464101615137754*f[37]*w0dx0+f[90]*dv2dx2+0.8944271909999159*f[116]*dv0dx0+f[15]*dv0dx0; 
  out[82] += 3.464101615137754*f[36]*w2dx2+3.464101615137754*f[37]*w1dx1+f[91]*dv2dx2+0.8944271909999159*f[122]*dv1dx1+f[12]*dv1dx1; 
  out[83] += 3.464101615137754*f[38]*w2dx2+3.464101615137754*f[39]*w1dx1+3.464101615137754*f[40]*w0dx0+0.8944271909999159*f[143]*dv2dx2+f[7]*dv2dx2+f[88]*dv1dx1+f[86]*dv0dx0; 
  out[84] += 3.464101615137754*f[41]*w1dx1+3.464101615137754*f[42]*w0dx0+f[90]*dv1dx1+0.8944271909999159*f[133]*dv0dx0+f[18]*dv0dx0; 
  out[85] += 3.464101615137754*f[41]*w2dx2+3.464101615137754*f[43]*w0dx0+0.8944271909999159*f[146]*dv2dx2+f[10]*dv2dx2+0.8944271909999159*f[134]*dv0dx0+f[19]*dv0dx0; 
  out[86] += 3.464101615137754*f[42]*w2dx2+3.464101615137754*f[43]*w1dx1+0.8944271909999159*f[147]*dv2dx2+f[11]*dv2dx2+f[92]*dv1dx1; 
  out[87] += 3.464101615137754*f[44]*w1dx1+3.464101615137754*f[45]*w0dx0+0.8944271909999159*f[139]*dv1dx1+f[17]*dv1dx1+f[91]*dv0dx0; 
  out[88] += 3.464101615137754*f[44]*w2dx2+3.464101615137754*f[46]*w0dx0+0.8944271909999159*f[149]*dv2dx2+f[13]*dv2dx2+f[92]*dv0dx0; 
  out[89] += 3.464101615137754*f[45]*w2dx2+3.464101615137754*f[46]*w1dx1+0.8944271909999159*f[150]*dv2dx2+f[14]*dv2dx2+0.8944271909999159*f[141]*dv1dx1+f[19]*dv1dx1; 
  out[90] += 3.464101615137754*f[47]*w0dx0+0.8944271909999159*f[138]*dv0dx0+f[21]*dv0dx0; 
  out[91] += 3.464101615137754*f[47]*w1dx1+0.8944271909999159*f[142]*dv1dx1+f[20]*dv1dx1; 
  out[92] += 3.464101615137754*f[47]*w2dx2+0.8944271909999159*f[152]*dv2dx2+f[16]*dv2dx2; 
  out[93] += 3.464101615137755*f[48]*w2dx2+3.464101615137755*f[50]*w1dx1+7.745966692414834*f[28]*w0dx0+f[123]*dv2dx2+f[107]*dv1dx1+2.23606797749979*f[78]*dv0dx0; 
  out[94] += 3.464101615137755*f[49]*w2dx2+7.745966692414834*f[28]*w1dx1+3.464101615137755*f[51]*w0dx0+f[124]*dv2dx2+2.23606797749979*f[79]*dv1dx1+f[99]*dv0dx0; 
  out[95] += 7.745966692414834*f[28]*w2dx2+3.464101615137755*f[52]*w1dx1+3.464101615137755*f[53]*w0dx0+2.23606797749979*f[83]*dv2dx2+f[109]*dv1dx1+f[101]*dv0dx0; 
  out[96] += 3.464101615137755*f[54]*w1dx1+7.745966692414834*f[29]*w0dx0+f[111]*dv1dx1+2.0*f[102]*dv0dx0+2.23606797749979*f[7]*dv0dx0; 
  out[97] += 7.745966692414834*f[29]*w1dx1+3.464101615137755*f[55]*w0dx0+2.23606797749979*f[80]*dv1dx1+f[23]*dv0dx0; 
  out[98] += 3.464101615137755*f[54]*w2dx2+7.745966692414834*f[30]*w0dx0+f[129]*dv2dx2+2.0*f[103]*dv0dx0+2.23606797749979*f[8]*dv0dx0; 
  out[99] += 3.464101615137755*f[55]*w2dx2+7.745966692414834*f[31]*w1dx1+f[130]*dv2dx2+2.23606797749979*f[82]*dv1dx1; 
  out[100] += 7.745966692414834*f[30]*w2dx2+3.464101615137755*f[56]*w0dx0+2.23606797749979*f[85]*dv2dx2+f[24]*dv0dx0; 
  out[101] += 7.745966692414834*f[31]*w2dx2+3.464101615137755*f[56]*w1dx1+2.23606797749979*f[86]*dv2dx2+f[113]*dv1dx1; 
  out[102] += 3.464101615137755*f[57]*w1dx1+3.464101615137755*f[58]*w0dx0+f[114]*dv1dx1+0.8944271909999159*f[11]*dv0dx0; 
  out[103] += 3.464101615137755*f[57]*w2dx2+3.464101615137755*f[59]*w0dx0+f[132]*dv2dx2+0.8944271909999159*f[12]*dv0dx0; 
  out[104] += 3.464101615137755*f[58]*w2dx2+3.464101615137755*f[59]*w1dx1+f[133]*dv2dx2+f[116]*dv1dx1; 
  out[105] += 3.464101615137755*f[60]*w1dx1+7.745966692414834*f[32]*w0dx0+f[22]*dv1dx1+2.23606797749979*f[80]*dv0dx0; 
  out[106] += 7.745966692414834*f[32]*w1dx1+3.464101615137755*f[61]*w0dx0+2.0*f[117]*dv1dx1+2.23606797749979*f[7]*dv1dx1+f[112]*dv0dx0; 
  out[107] += 3.464101615137755*f[60]*w2dx2+7.745966692414834*f[33]*w0dx0+f[135]*dv2dx2+2.23606797749979*f[81]*dv0dx0; 
  out[108] += 3.464101615137755*f[61]*w2dx2+7.745966692414834*f[34]*w1dx1+f[136]*dv2dx2+2.0*f[119]*dv1dx1+2.23606797749979*f[9]*dv1dx1; 
  out[109] += 7.745966692414834*f[33]*w2dx2+3.464101615137755*f[62]*w0dx0+2.23606797749979*f[88]*dv2dx2+f[113]*dv0dx0; 
  out[110] += 7.745966692414834*f[34]*w2dx2+3.464101615137755*f[62]*w1dx1+2.23606797749979*f[89]*dv2dx2+f[24]*dv1dx1; 
  out[111] += 7.745966692414834*f[35]*w0dx0+2.0*f[114]*dv0dx0+2.23606797749979*f[13]*dv0dx0; 
  out[112] += 7.745966692414834*f[36]*w1dx1+2.0*f[121]*dv1dx1+2.23606797749979*f[11]*dv1dx1; 
  out[113] += 7.745966692414834*f[37]*w2dx2+2.23606797749979*f[92]*dv2dx2; 
  out[114] += 3.464101615137755*f[63]*w0dx0+0.8944271909999159*f[16]*dv0dx0; 
  out[115] += 3.464101615137755*f[63]*w1dx1+f[25]*dv1dx1; 
  out[116] += 3.464101615137755*f[63]*w2dx2+f[138]*dv2dx2; 
  out[117] += 3.464101615137755*f[64]*w1dx1+3.464101615137755*f[65]*w0dx0+0.8944271909999159*f[13]*dv1dx1+f[121]*dv0dx0; 
  out[118] += 3.464101615137755*f[64]*w2dx2+3.464101615137755*f[66]*w0dx0+f[139]*dv2dx2+f[122]*dv0dx0; 
  out[119] += 3.464101615137755*f[65]*w2dx2+3.464101615137755*f[66]*w1dx1+f[140]*dv2dx2+0.8944271909999159*f[15]*dv1dx1; 
  out[120] += 3.464101615137755*f[67]*w0dx0+f[26]*dv0dx0; 
  out[121] += 3.464101615137755*f[67]*w1dx1+0.8944271909999159*f[16]*dv1dx1; 
  out[122] += 3.464101615137755*f[67]*w2dx2+f[142]*dv2dx2; 
  out[123] += 3.464101615137755*f[68]*w1dx1+7.745966692414834*f[38]*w0dx0+f[135]*dv1dx1+2.23606797749979*f[84]*dv0dx0; 
  out[124] += 7.745966692414834*f[38]*w1dx1+3.464101615137755*f[69]*w0dx0+2.23606797749979*f[87]*dv1dx1+f[130]*dv0dx0; 
  out[125] += 3.464101615137755*f[68]*w2dx2+7.745966692414834*f[39]*w0dx0+f[22]*dv2dx2+2.23606797749979*f[85]*dv0dx0; 
  out[126] += 3.464101615137755*f[69]*w2dx2+7.745966692414834*f[40]*w1dx1+f[23]*dv2dx2+2.23606797749979*f[89]*dv1dx1; 
  out[127] += 7.745966692414834*f[39]*w2dx2+3.464101615137755*f[70]*w0dx0+2.0*f[144]*dv2dx2+2.23606797749979*f[8]*dv2dx2+f[131]*dv0dx0; 
  out[128] += 7.745966692414834*f[40]*w2dx2+3.464101615137755*f[70]*w1dx1+2.0*f[145]*dv2dx2+2.23606797749979*f[9]*dv2dx2+f[137]*dv1dx1; 
  out[129] += 7.745966692414834*f[41]*w0dx0+2.0*f[132]*dv0dx0+2.23606797749979*f[17]*dv0dx0; 
  out[130] += 7.745966692414834*f[42]*w1dx1+2.23606797749979*f[91]*dv1dx1; 
  out[131] += 7.745966692414834*f[43]*w2dx2+2.0*f[148]*dv2dx2+2.23606797749979*f[12]*dv2dx2; 
  out[132] += 3.464101615137755*f[71]*w0dx0+0.8944271909999159*f[20]*dv0dx0; 
  out[133] += 3.464101615137755*f[71]*w1dx1+f[138]*dv1dx1; 
  out[134] += 3.464101615137755*f[71]*w2dx2+f[25]*dv2dx2; 
  out[135] += 7.745966692414834*f[44]*w0dx0+2.23606797749979*f[90]*dv0dx0; 
  out[136] += 7.745966692414834*f[45]*w1dx1+2.0*f[140]*dv1dx1+2.23606797749979*f[18]*dv1dx1; 
  out[137] += 7.745966692414834*f[46]*w2dx2+2.0*f[151]*dv2dx2+2.23606797749979*f[15]*dv2dx2; 
  out[139] += 3.464101615137755*f[72]*w0dx0+f[142]*dv0dx0; 
  out[140] += 3.464101615137755*f[72]*w1dx1+0.8944271909999159*f[21]*dv1dx1; 
  out[141] += 3.464101615137755*f[72]*w2dx2+f[26]*dv2dx2; 
  out[143] += 3.464101615137755*f[73]*w1dx1+3.464101615137755*f[74]*w0dx0+f[149]*dv1dx1+f[147]*dv0dx0; 
  out[144] += 3.464101615137755*f[73]*w2dx2+3.464101615137755*f[75]*w0dx0+0.8944271909999159*f[17]*dv2dx2+f[148]*dv0dx0; 
  out[145] += 3.464101615137755*f[74]*w2dx2+3.464101615137755*f[75]*w1dx1+0.8944271909999159*f[18]*dv2dx2+f[151]*dv1dx1; 
  out[146] += 3.464101615137755*f[76]*w0dx0+f[27]*dv0dx0; 
  out[147] += 3.464101615137755*f[76]*w1dx1+f[152]*dv1dx1; 
  out[148] += 3.464101615137755*f[76]*w2dx2+0.8944271909999159*f[20]*dv2dx2; 
  out[149] += 3.464101615137755*f[77]*w0dx0+f[152]*dv0dx0; 
  out[150] += 3.464101615137755*f[77]*w1dx1+f[27]*dv1dx1; 
  out[151] += 3.464101615137755*f[77]*w2dx2+0.8944271909999159*f[21]*dv2dx2; 
  out[153] += 3.464101615137754*f[80]*w2dx2+3.464101615137754*f[81]*w1dx1+3.464101615137754*f[82]*w0dx0+f[156]*dv2dx2+0.8944271909999161*f[177]*dv1dx1+f[30]*dv1dx1+0.8944271909999161*f[174]*dv0dx0+f[34]*dv0dx0; 
  out[154] += 3.464101615137754*f[84]*w2dx2+3.464101615137754*f[85]*w1dx1+3.464101615137754*f[86]*w0dx0+0.8944271909999161*f[210]*dv2dx2+f[29]*dv2dx2+f[157]*dv1dx1+0.8944271909999161*f[190]*dv0dx0+f[40]*dv0dx0; 
  out[155] += 3.464101615137754*f[87]*w2dx2+3.464101615137754*f[88]*w1dx1+3.464101615137754*f[89]*w0dx0+0.8944271909999161*f[213]*dv2dx2+f[32]*dv2dx2+0.8944271909999161*f[204]*dv1dx1+f[39]*dv1dx1+f[158]*dv0dx0; 
  out[156] += 3.464101615137754*f[90]*w1dx1+3.464101615137754*f[91]*w0dx0+0.8944271909999161*f[206]*dv1dx1+f[41]*dv1dx1+0.8944271909999161*f[201]*dv0dx0+f[45]*dv0dx0; 
  out[157] += 3.464101615137754*f[90]*w2dx2+3.464101615137754*f[92]*w0dx0+0.8944271909999161*f[216]*dv2dx2+f[35]*dv2dx2+0.8944271909999161*f[202]*dv0dx0+f[46]*dv0dx0; 
  out[158] += 3.464101615137754*f[91]*w2dx2+3.464101615137754*f[92]*w1dx1+0.8944271909999161*f[217]*dv2dx2+f[36]*dv2dx2+0.8944271909999161*f[208]*dv1dx1+f[43]*dv1dx1; 
  out[159] += 3.464101615137755*f[96]*w2dx2+3.464101615137755*f[98]*w1dx1+7.745966692414834*f[78]*w0dx0+f[182]*dv2dx2+f[168]*dv1dx1+2.0*f[162]*dv0dx0+2.23606797749979*f[28]*dv0dx0; 
  out[160] += 3.464101615137755*f[97]*w2dx2+7.745966692414834*f[78]*w1dx1+3.464101615137755*f[99]*w0dx0+f[183]*dv2dx2+2.23606797749979*f[153]*dv1dx1+f[51]*dv0dx0; 
  out[161] += 7.745966692414834*f[78]*w2dx2+3.464101615137755*f[100]*w1dx1+3.464101615137755*f[101]*w0dx0+2.23606797749979*f[154]*dv2dx2+f[170]*dv1dx1+f[53]*dv0dx0; 
  out[162] += 3.464101615137755*f[102]*w2dx2+3.464101615137755*f[103]*w1dx1+3.464101615137755*f[104]*w0dx0+f[188]*dv2dx2+f[173]*dv1dx1+0.8944271909999161*f[31]*dv0dx0; 
  out[163] += 3.464101615137755*f[105]*w2dx2+3.464101615137755*f[107]*w1dx1+7.745966692414834*f[79]*w0dx0+f[191]*dv2dx2+f[50]*dv1dx1+2.23606797749979*f[153]*dv0dx0; 
  out[164] += 3.464101615137755*f[106]*w2dx2+7.745966692414834*f[79]*w1dx1+3.464101615137755*f[108]*w0dx0+f[192]*dv2dx2+2.0*f[175]*dv1dx1+2.23606797749979*f[28]*dv1dx1+f[169]*dv0dx0; 
  out[165] += 7.745966692414834*f[79]*w2dx2+3.464101615137755*f[109]*w1dx1+3.464101615137755*f[110]*w0dx0+2.23606797749979*f[155]*dv2dx2+f[52]*dv1dx1+f[171]*dv0dx0; 
  out[166] += 3.464101615137755*f[111]*w1dx1+7.745966692414834*f[80]*w0dx0+f[54]*dv1dx1+2.0*f[172]*dv0dx0+2.23606797749979*f[32]*dv0dx0; 
  out[167] += 7.745966692414834*f[80]*w1dx1+3.464101615137755*f[112]*w0dx0+2.0*f[176]*dv1dx1+2.23606797749979*f[29]*dv1dx1+f[61]*dv0dx0; 
  out[168] += 3.464101615137755*f[111]*w2dx2+7.745966692414834*f[81]*w0dx0+f[197]*dv2dx2+2.0*f[173]*dv0dx0+2.23606797749979*f[33]*dv0dx0; 
  out[169] += 3.464101615137755*f[112]*w2dx2+7.745966692414834*f[82]*w1dx1+f[198]*dv2dx2+2.0*f[178]*dv1dx1+2.23606797749979*f[31]*dv1dx1; 
  out[170] += 7.745966692414834*f[81]*w2dx2+3.464101615137755*f[113]*w0dx0+2.23606797749979*f[157]*dv2dx2+f[62]*dv0dx0; 
  out[171] += 7.745966692414834*f[82]*w2dx2+3.464101615137755*f[113]*w1dx1+2.23606797749979*f[158]*dv2dx2+f[56]*dv1dx1; 
  out[172] += 3.464101615137755*f[114]*w1dx1+3.464101615137755*f[115]*w0dx0+f[57]*dv1dx1+0.8944271909999161*f[36]*dv0dx0; 
  out[173] += 3.464101615137755*f[114]*w2dx2+3.464101615137755*f[116]*w0dx0+f[200]*dv2dx2+0.8944271909999161*f[37]*dv0dx0; 
  out[174] += 3.464101615137755*f[115]*w2dx2+3.464101615137755*f[116]*w1dx1+f[201]*dv2dx2+f[59]*dv1dx1; 
  out[175] += 3.464101615137755*f[117]*w2dx2+3.464101615137755*f[118]*w1dx1+3.464101615137755*f[119]*w0dx0+f[203]*dv2dx2+0.8944271909999161*f[33]*dv1dx1+f[178]*dv0dx0; 
  out[176] += 3.464101615137755*f[120]*w1dx1+3.464101615137755*f[121]*w0dx0+0.8944271909999161*f[35]*dv1dx1+f[65]*dv0dx0; 
  out[177] += 3.464101615137755*f[120]*w2dx2+3.464101615137755*f[122]*w0dx0+f[206]*dv2dx2+f[66]*dv0dx0; 
  out[178] += 3.464101615137755*f[121]*w2dx2+3.464101615137755*f[122]*w1dx1+f[207]*dv2dx2+0.8944271909999161*f[37]*dv1dx1; 
  out[179] += 3.464101615137755*f[123]*w2dx2+3.464101615137755*f[125]*w1dx1+7.745966692414834*f[83]*w0dx0+f[48]*dv2dx2+f[193]*dv1dx1+2.23606797749979*f[154]*dv0dx0; 
  out[180] += 3.464101615137755*f[124]*w2dx2+7.745966692414834*f[83]*w1dx1+3.464101615137755*f[126]*w0dx0+f[49]*dv2dx2+2.23606797749979*f[155]*dv1dx1+f[185]*dv0dx0; 
  out[181] += 7.745966692414834*f[83]*w2dx2+3.464101615137755*f[127]*w1dx1+3.464101615137755*f[128]*w0dx0+2.0*f[209]*dv2dx2+2.23606797749979*f[28]*dv2dx2+f[195]*dv1dx1+f[187]*dv0dx0; 
  out[182] += 3.464101615137755*f[129]*w1dx1+7.745966692414834*f[84]*w0dx0+f[197]*dv1dx1+2.0*f[188]*dv0dx0+2.23606797749979*f[38]*dv0dx0; 
  out[183] += 7.745966692414834*f[84]*w1dx1+3.464101615137755*f[130]*w0dx0+2.23606797749979*f[156]*dv1dx1+f[69]*dv0dx0; 
  out[184] += 3.464101615137755*f[129]*w2dx2+7.745966692414834*f[85]*w0dx0+f[54]*dv2dx2+2.0*f[189]*dv0dx0+2.23606797749979*f[39]*dv0dx0; 
  out[185] += 3.464101615137755*f[130]*w2dx2+7.745966692414834*f[86]*w1dx1+f[55]*dv2dx2+2.23606797749979*f[158]*dv1dx1; 
  out[186] += 7.745966692414834*f[85]*w2dx2+3.464101615137755*f[131]*w0dx0+2.0*f[211]*dv2dx2+2.23606797749979*f[30]*dv2dx2+f[70]*dv0dx0; 
  out[187] += 7.745966692414834*f[86]*w2dx2+3.464101615137755*f[131]*w1dx1+2.0*f[212]*dv2dx2+2.23606797749979*f[31]*dv2dx2+f[199]*dv1dx1; 
  out[188] += 3.464101615137755*f[132]*w1dx1+3.464101615137755*f[133]*w0dx0+f[200]*dv1dx1+0.8944271909999161*f[42]*dv0dx0; 
  out[189] += 3.464101615137755*f[132]*w2dx2+3.464101615137755*f[134]*w0dx0+f[57]*dv2dx2+0.8944271909999161*f[43]*dv0dx0; 
  out[190] += 3.464101615137755*f[133]*w2dx2+3.464101615137755*f[134]*w1dx1+f[58]*dv2dx2+f[202]*dv1dx1; 
  out[191] += 3.464101615137755*f[135]*w1dx1+7.745966692414834*f[87]*w0dx0+f[68]*dv1dx1+2.23606797749979*f[156]*dv0dx0; 
  out[192] += 7.745966692414834*f[87]*w1dx1+3.464101615137755*f[136]*w0dx0+2.0*f[203]*dv1dx1+2.23606797749979*f[38]*dv1dx1+f[198]*dv0dx0; 
  out[193] += 3.464101615137755*f[135]*w2dx2+7.745966692414834*f[88]*w0dx0+f[60]*dv2dx2+2.23606797749979*f[157]*dv0dx0; 
  out[194] += 3.464101615137755*f[136]*w2dx2+7.745966692414834*f[89]*w1dx1+f[61]*dv2dx2+2.0*f[205]*dv1dx1+2.23606797749979*f[40]*dv1dx1; 
  out[195] += 7.745966692414834*f[88]*w2dx2+3.464101615137755*f[137]*w0dx0+2.0*f[214]*dv2dx2+2.23606797749979*f[33]*dv2dx2+f[199]*dv0dx0; 
  out[196] += 7.745966692414834*f[89]*w2dx2+3.464101615137755*f[137]*w1dx1+2.0*f[215]*dv2dx2+2.23606797749979*f[34]*dv2dx2+f[70]*dv1dx1; 
  out[197] += 7.745966692414834*f[90]*w0dx0+2.0*f[200]*dv0dx0+2.23606797749979*f[44]*dv0dx0; 
  out[198] += 7.745966692414834*f[91]*w1dx1+2.0*f[207]*dv1dx1+2.23606797749979*f[42]*dv1dx1; 
  out[199] += 7.745966692414834*f[92]*w2dx2+2.0*f[218]*dv2dx2+2.23606797749979*f[37]*dv2dx2; 
  out[200] += 3.464101615137755*f[138]*w0dx0+0.8944271909999161*f[47]*dv0dx0; 
  out[201] += 3.464101615137755*f[138]*w1dx1+f[71]*dv1dx1; 
  out[202] += 3.464101615137755*f[138]*w2dx2+f[63]*dv2dx2; 
  out[203] += 3.464101615137755*f[139]*w1dx1+3.464101615137755*f[140]*w0dx0+0.8944271909999161*f[44]*dv1dx1+f[207]*dv0dx0; 
  out[204] += 3.464101615137755*f[139]*w2dx2+3.464101615137755*f[141]*w0dx0+f[64]*dv2dx2+f[208]*dv0dx0; 
  out[205] += 3.464101615137755*f[140]*w2dx2+3.464101615137755*f[141]*w1dx1+f[65]*dv2dx2+0.8944271909999161*f[46]*dv1dx1; 
  out[206] += 3.464101615137755*f[142]*w0dx0+f[72]*dv0dx0; 
  out[207] += 3.464101615137755*f[142]*w1dx1+0.8944271909999161*f[47]*dv1dx1; 
  out[208] += 3.464101615137755*f[142]*w2dx2+f[67]*dv2dx2; 
  out[209] += 3.464101615137755*f[143]*w2dx2+3.464101615137755*f[144]*w1dx1+3.464101615137755*f[145]*w0dx0+0.8944271909999161*f[38]*dv2dx2+f[214]*dv1dx1+f[212]*dv0dx0; 
  out[210] += 3.464101615137755*f[146]*w1dx1+3.464101615137755*f[147]*w0dx0+f[216]*dv1dx1+f[74]*dv0dx0; 
  out[211] += 3.464101615137755*f[146]*w2dx2+3.464101615137755*f[148]*w0dx0+0.8944271909999161*f[41]*dv2dx2+f[75]*dv0dx0; 
  out[212] += 3.464101615137755*f[147]*w2dx2+3.464101615137755*f[148]*w1dx1+0.8944271909999161*f[42]*dv2dx2+f[218]*dv1dx1; 
  out[213] += 3.464101615137755*f[149]*w1dx1+3.464101615137755*f[150]*w0dx0+f[73]*dv1dx1+f[217]*dv0dx0; 
  out[214] += 3.464101615137755*f[149]*w2dx2+3.464101615137755*f[151]*w0dx0+0.8944271909999161*f[44]*dv2dx2+f[218]*dv0dx0; 
  out[215] += 3.464101615137755*f[150]*w2dx2+3.464101615137755*f[151]*w1dx1+0.8944271909999161*f[45]*dv2dx2+f[75]*dv1dx1; 
  out[216] += 3.464101615137755*f[152]*w0dx0+f[77]*dv0dx0; 
  out[217] += 3.464101615137755*f[152]*w1dx1+f[76]*dv1dx1; 
  out[218] += 3.464101615137755*f[152]*w2dx2+0.8944271909999161*f[47]*dv2dx2; 
  out[219] += 3.464101615137754*f[156]*w2dx2+3.464101615137754*f[157]*w1dx1+3.464101615137754*f[158]*w0dx0+0.8944271909999159*f[247]*dv2dx2+f[80]*dv2dx2+0.8944271909999159*f[243]*dv1dx1+f[85]*dv1dx1+0.8944271909999159*f[240]*dv0dx0+f[89]*dv0dx0; 
  out[220] += 3.464101615137755*f[166]*w2dx2+3.464101615137755*f[168]*w1dx1+7.745966692414834*f[153]*w0dx0+f[232]*dv2dx2+f[98]*dv1dx1+2.0*f[223]*dv0dx0+2.23606797749979*f[79]*dv0dx0; 
  out[221] += 3.464101615137755*f[167]*w2dx2+7.745966692414834*f[153]*w1dx1+3.464101615137755*f[169]*w0dx0+f[233]*dv2dx2+2.0*f[224]*dv1dx1+2.23606797749979*f[78]*dv1dx1+f[108]*dv0dx0; 
  out[222] += 7.745966692414834*f[153]*w2dx2+3.464101615137755*f[170]*w1dx1+3.464101615137755*f[171]*w0dx0+2.23606797749979*f[219]*dv2dx2+f[100]*dv1dx1+f[110]*dv0dx0; 
  out[223] += 3.464101615137755*f[172]*w2dx2+3.464101615137755*f[173]*w1dx1+3.464101615137755*f[174]*w0dx0+f[238]*dv2dx2+f[103]*dv1dx1+0.8944271909999159*f[82]*dv0dx0; 
  out[224] += 3.464101615137755*f[176]*w2dx2+3.464101615137755*f[177]*w1dx1+3.464101615137755*f[178]*w0dx0+f[242]*dv2dx2+0.8944271909999159*f[81]*dv1dx1+f[119]*dv0dx0; 
  out[225] += 3.464101615137755*f[182]*w2dx2+3.464101615137755*f[184]*w1dx1+7.745966692414834*f[154]*w0dx0+f[96]*dv2dx2+f[234]*dv1dx1+2.0*f[228]*dv0dx0+2.23606797749979*f[83]*dv0dx0; 
  out[226] += 3.464101615137755*f[183]*w2dx2+7.745966692414834*f[154]*w1dx1+3.464101615137755*f[185]*w0dx0+f[97]*dv2dx2+2.23606797749979*f[219]*dv1dx1+f[126]*dv0dx0; 
  out[227] += 7.745966692414834*f[154]*w2dx2+3.464101615137755*f[186]*w1dx1+3.464101615137755*f[187]*w0dx0+2.0*f[245]*dv2dx2+2.23606797749979*f[78]*dv2dx2+f[236]*dv1dx1+f[128]*dv0dx0; 
  out[228] += 3.464101615137755*f[188]*w2dx2+3.464101615137755*f[189]*w1dx1+3.464101615137755*f[190]*w0dx0+f[102]*dv2dx2+f[239]*dv1dx1+0.8944271909999159*f[86]*dv0dx0; 
  out[229] += 3.464101615137755*f[191]*w2dx2+3.464101615137755*f[193]*w1dx1+7.745966692414834*f[155]*w0dx0+f[105]*dv2dx2+f[125]*dv1dx1+2.23606797749979*f[219]*dv0dx0; 
  out[230] += 3.464101615137755*f[192]*w2dx2+7.745966692414834*f[155]*w1dx1+3.464101615137755*f[194]*w0dx0+f[106]*dv2dx2+2.0*f[241]*dv1dx1+2.23606797749979*f[83]*dv1dx1+f[235]*dv0dx0; 
  out[231] += 7.745966692414834*f[155]*w2dx2+3.464101615137755*f[195]*w1dx1+3.464101615137755*f[196]*w0dx0+2.0*f[246]*dv2dx2+2.23606797749979*f[79]*dv2dx2+f[127]*dv1dx1+f[237]*dv0dx0; 
  out[232] += 3.464101615137755*f[197]*w1dx1+7.745966692414834*f[156]*w0dx0+f[129]*dv1dx1+2.0*f[238]*dv0dx0+2.23606797749979*f[87]*dv0dx0; 
  out[233] += 7.745966692414834*f[156]*w1dx1+3.464101615137755*f[198]*w0dx0+2.0*f[242]*dv1dx1+2.23606797749979*f[84]*dv1dx1+f[136]*dv0dx0; 
  out[234] += 3.464101615137755*f[197]*w2dx2+7.745966692414834*f[157]*w0dx0+f[111]*dv2dx2+2.0*f[239]*dv0dx0+2.23606797749979*f[88]*dv0dx0; 
  out[235] += 3.464101615137755*f[198]*w2dx2+7.745966692414834*f[158]*w1dx1+f[112]*dv2dx2+2.0*f[244]*dv1dx1+2.23606797749979*f[86]*dv1dx1; 
  out[236] += 7.745966692414834*f[157]*w2dx2+3.464101615137755*f[199]*w0dx0+2.0*f[248]*dv2dx2+2.23606797749979*f[81]*dv2dx2+f[137]*dv0dx0; 
  out[237] += 7.745966692414834*f[158]*w2dx2+3.464101615137755*f[199]*w1dx1+2.0*f[249]*dv2dx2+2.23606797749979*f[82]*dv2dx2+f[131]*dv1dx1; 
  out[238] += 3.464101615137755*f[200]*w1dx1+3.464101615137755*f[201]*w0dx0+f[132]*dv1dx1+0.8944271909999159*f[91]*dv0dx0; 
  out[239] += 3.464101615137755*f[200]*w2dx2+3.464101615137755*f[202]*w0dx0+f[114]*dv2dx2+0.8944271909999159*f[92]*dv0dx0; 
  out[240] += 3.464101615137755*f[201]*w2dx2+3.464101615137755*f[202]*w1dx1+f[115]*dv2dx2+f[134]*dv1dx1; 
  out[241] += 3.464101615137755*f[203]*w2dx2+3.464101615137755*f[204]*w1dx1+3.464101615137755*f[205]*w0dx0+f[117]*dv2dx2+0.8944271909999159*f[88]*dv1dx1+f[244]*dv0dx0; 
  out[242] += 3.464101615137755*f[206]*w1dx1+3.464101615137755*f[207]*w0dx0+0.8944271909999159*f[90]*dv1dx1+f[140]*dv0dx0; 
  out[243] += 3.464101615137755*f[206]*w2dx2+3.464101615137755*f[208]*w0dx0+f[120]*dv2dx2+f[141]*dv0dx0; 
  out[244] += 3.464101615137755*f[207]*w2dx2+3.464101615137755*f[208]*w1dx1+f[121]*dv2dx2+0.8944271909999159*f[92]*dv1dx1; 
  out[245] += 3.464101615137755*f[210]*w2dx2+3.464101615137755*f[211]*w1dx1+3.464101615137755*f[212]*w0dx0+0.8944271909999159*f[84]*dv2dx2+f[248]*dv1dx1+f[145]*dv0dx0; 
  out[246] += 3.464101615137755*f[213]*w2dx2+3.464101615137755*f[214]*w1dx1+3.464101615137755*f[215]*w0dx0+0.8944271909999159*f[87]*dv2dx2+f[144]*dv1dx1+f[249]*dv0dx0; 
  out[247] += 3.464101615137755*f[216]*w1dx1+3.464101615137755*f[217]*w0dx0+f[146]*dv1dx1+f[150]*dv0dx0; 
  out[248] += 3.464101615137755*f[216]*w2dx2+3.464101615137755*f[218]*w0dx0+0.8944271909999159*f[90]*dv2dx2+f[151]*dv0dx0; 
  out[249] += 3.464101615137755*f[217]*w2dx2+3.464101615137755*f[218]*w1dx1+0.8944271909999159*f[91]*dv2dx2+f[148]*dv1dx1; 
  out[250] += 3.464101615137755*f[232]*w2dx2+3.464101615137755*f[234]*w1dx1+7.745966692414834*f[219]*w0dx0+f[166]*dv2dx2+f[184]*dv1dx1+2.0*f[253]*dv0dx0+2.23606797749979*f[155]*dv0dx0; 
  out[251] += 3.464101615137755*f[233]*w2dx2+7.745966692414834*f[219]*w1dx1+3.464101615137755*f[235]*w0dx0+f[167]*dv2dx2+2.0*f[254]*dv1dx1+2.23606797749979*f[154]*dv1dx1+f[194]*dv0dx0; 
  out[252] += 7.745966692414834*f[219]*w2dx2+3.464101615137755*f[236]*w1dx1+3.464101615137755*f[237]*w0dx0+2.0*f[255]*dv2dx2+2.23606797749979*f[153]*dv2dx2+f[186]*dv1dx1+f[196]*dv0dx0; 
  out[253] += 3.464101615137755*f[238]*w2dx2+3.464101615137755*f[239]*w1dx1+3.464101615137755*f[240]*w0dx0+f[172]*dv2dx2+f[189]*dv1dx1+0.8944271909999161*f[158]*dv0dx0; 
  out[254] += 3.464101615137755*f[242]*w2dx2+3.464101615137755*f[243]*w1dx1+3.464101615137755*f[244]*w0dx0+f[176]*dv2dx2+0.8944271909999161*f[157]*dv1dx1+f[205]*dv0dx0; 
  out[255] += 3.464101615137755*f[247]*w2dx2+3.464101615137755*f[248]*w1dx1+3.464101615137755*f[249]*w0dx0+0.8944271909999161*f[156]*dv2dx2+f[211]*dv1dx1+f[215]*dv0dx0; 
return std::abs(w0dx0)+std::abs(w1dx1)+std::abs(w2dx2); 
} 
