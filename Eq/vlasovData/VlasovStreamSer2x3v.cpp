#include <VlasovModDecl.h> 
double VlasovVolStream2x3vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[6] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[9]*dv1dx1+f[8]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w0dx0+f[0]*dv0dx0; 
  out[8] += 3.464101615137754*f[3]*w1dx1+f[11]*dv1dx1; 
  out[9] += 3.464101615137754*f[4]*w0dx0+f[11]*dv0dx0; 
  out[10] += 3.464101615137754*f[4]*w1dx1+f[0]*dv1dx1; 
  out[12] += 3.464101615137754*f[5]*w0dx0+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[5]*w1dx1+f[15]*dv1dx1; 
  out[16] += 3.464101615137754*f[7]*w1dx1+3.464101615137754*f[8]*w0dx0+f[18]*dv1dx1+f[2]*dv0dx0; 
  out[17] += 3.464101615137754*f[9]*w1dx1+3.464101615137754*f[10]*w0dx0+f[1]*dv1dx1+f[19]*dv0dx0; 
  out[18] += 3.464101615137754*f[11]*w0dx0+f[4]*dv0dx0; 
  out[19] += 3.464101615137754*f[11]*w1dx1+f[3]*dv1dx1; 
  out[20] += 3.464101615137754*f[12]*w1dx1+3.464101615137754*f[13]*w0dx0+f[23]*dv1dx1+f[22]*dv0dx0; 
  out[21] += 3.464101615137754*f[14]*w0dx0+f[5]*dv0dx0; 
  out[22] += 3.464101615137754*f[14]*w1dx1+f[25]*dv1dx1; 
  out[23] += 3.464101615137754*f[15]*w0dx0+f[25]*dv0dx0; 
  out[24] += 3.464101615137754*f[15]*w1dx1+f[5]*dv1dx1; 
  out[26] += 3.464101615137754*f[18]*w1dx1+3.464101615137754*f[19]*w0dx0+f[7]*dv1dx1+f[10]*dv0dx0; 
  out[27] += 3.464101615137754*f[21]*w1dx1+3.464101615137754*f[22]*w0dx0+f[29]*dv1dx1+f[13]*dv0dx0; 
  out[28] += 3.464101615137754*f[23]*w1dx1+3.464101615137754*f[24]*w0dx0+f[12]*dv1dx1+f[30]*dv0dx0; 
  out[29] += 3.464101615137754*f[25]*w0dx0+f[15]*dv0dx0; 
  out[30] += 3.464101615137754*f[25]*w1dx1+f[14]*dv1dx1; 
  out[31] += 3.464101615137754*f[29]*w1dx1+3.464101615137754*f[30]*w0dx0+f[21]*dv1dx1+f[24]*dv0dx0; 
return std::abs(w0dx0)+std::abs(w1dx1); 
} 
double VlasovVolStream2x3vSerP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[6] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[9]*dv1dx1+f[8]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w0dx0+0.8944271909999159*f[18]*dv0dx0+f[0]*dv0dx0; 
  out[8] += 3.464101615137754*f[3]*w1dx1+f[11]*dv1dx1; 
  out[9] += 3.464101615137754*f[4]*w0dx0+f[11]*dv0dx0; 
  out[10] += 3.464101615137754*f[4]*w1dx1+0.8944271909999159*f[19]*dv1dx1+f[0]*dv1dx1; 
  out[12] += 3.464101615137754*f[5]*w0dx0+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[5]*w1dx1+f[15]*dv1dx1; 
  out[16] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[7]*dv0dx0; 
  out[17] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[10]*dv1dx1; 
  out[21] += 3.464101615137754*f[7]*w1dx1+3.464101615137754*f[8]*w0dx0+f[23]*dv1dx1+0.8944271909999161*f[36]*dv0dx0+f[2]*dv0dx0; 
  out[22] += 3.464101615137754*f[9]*w1dx1+3.464101615137754*f[10]*w0dx0+0.8944271909999161*f[40]*dv1dx1+f[1]*dv1dx1+f[24]*dv0dx0; 
  out[23] += 3.464101615137754*f[11]*w0dx0+0.8944271909999161*f[39]*dv0dx0+f[4]*dv0dx0; 
  out[24] += 3.464101615137754*f[11]*w1dx1+0.8944271909999161*f[42]*dv1dx1+f[3]*dv1dx1; 
  out[25] += 3.464101615137754*f[12]*w1dx1+3.464101615137754*f[13]*w0dx0+f[28]*dv1dx1+f[27]*dv0dx0; 
  out[26] += 3.464101615137754*f[14]*w0dx0+0.8944271909999161*f[45]*dv0dx0+f[5]*dv0dx0; 
  out[27] += 3.464101615137754*f[14]*w1dx1+f[30]*dv1dx1; 
  out[28] += 3.464101615137754*f[15]*w0dx0+f[30]*dv0dx0; 
  out[29] += 3.464101615137754*f[15]*w1dx1+0.8944271909999161*f[46]*dv1dx1+f[5]*dv1dx1; 
  out[31] += 3.464101615137755*f[16]*w1dx1+7.745966692414834*f[6]*w0dx0+f[37]*dv1dx1+2.23606797749979*f[21]*dv0dx0; 
  out[32] += 7.745966692414834*f[6]*w1dx1+3.464101615137755*f[17]*w0dx0+2.23606797749979*f[22]*dv1dx1+f[34]*dv0dx0; 
  out[33] += 7.745966692414834*f[7]*w0dx0+2.0*f[35]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[34] += 7.745966692414834*f[8]*w1dx1+2.23606797749979*f[24]*dv1dx1; 
  out[35] += 3.464101615137755*f[18]*w0dx0+0.8944271909999161*f[3]*dv0dx0; 
  out[36] += 3.464101615137755*f[18]*w1dx1+f[39]*dv1dx1; 
  out[37] += 7.745966692414834*f[9]*w0dx0+2.23606797749979*f[23]*dv0dx0; 
  out[38] += 7.745966692414834*f[10]*w1dx1+2.0*f[41]*dv1dx1+2.23606797749979*f[2]*dv1dx1; 
  out[40] += 3.464101615137755*f[19]*w0dx0+f[42]*dv0dx0; 
  out[41] += 3.464101615137755*f[19]*w1dx1+0.8944271909999161*f[4]*dv1dx1; 
  out[43] += 7.745966692414834*f[12]*w0dx0+2.23606797749979*f[26]*dv0dx0; 
  out[44] += 7.745966692414834*f[13]*w1dx1+2.23606797749979*f[29]*dv1dx1; 
  out[47] += 3.464101615137755*f[20]*w0dx0+f[49]*dv0dx0; 
  out[48] += 3.464101615137755*f[20]*w1dx1+f[50]*dv1dx1; 
  out[51] += 3.464101615137754*f[23]*w1dx1+3.464101615137754*f[24]*w0dx0+0.8944271909999159*f[66]*dv1dx1+f[7]*dv1dx1+0.8944271909999159*f[64]*dv0dx0+f[10]*dv0dx0; 
  out[52] += 3.464101615137754*f[26]*w1dx1+3.464101615137754*f[27]*w0dx0+f[54]*dv1dx1+0.8944271909999159*f[73]*dv0dx0+f[13]*dv0dx0; 
  out[53] += 3.464101615137754*f[28]*w1dx1+3.464101615137754*f[29]*w0dx0+0.8944271909999159*f[77]*dv1dx1+f[12]*dv1dx1+f[55]*dv0dx0; 
  out[54] += 3.464101615137754*f[30]*w0dx0+0.8944271909999159*f[76]*dv0dx0+f[15]*dv0dx0; 
  out[55] += 3.464101615137754*f[30]*w1dx1+0.8944271909999159*f[79]*dv1dx1+f[14]*dv1dx1; 
  out[56] += 3.464101615137755*f[33]*w1dx1+7.745966692414834*f[21]*w0dx0+f[61]*dv1dx1+2.0*f[58]*dv0dx0+2.23606797749979*f[6]*dv0dx0; 
  out[57] += 7.745966692414834*f[21]*w1dx1+3.464101615137755*f[34]*w0dx0+2.23606797749979*f[51]*dv1dx1+f[17]*dv0dx0; 
  out[58] += 3.464101615137755*f[35]*w1dx1+3.464101615137755*f[36]*w0dx0+f[63]*dv1dx1+0.8944271909999159*f[8]*dv0dx0; 
  out[59] += 3.464101615137755*f[37]*w1dx1+7.745966692414834*f[22]*w0dx0+f[16]*dv1dx1+2.23606797749979*f[51]*dv0dx0; 
  out[60] += 7.745966692414834*f[22]*w1dx1+3.464101615137755*f[38]*w0dx0+2.0*f[65]*dv1dx1+2.23606797749979*f[6]*dv1dx1+f[62]*dv0dx0; 
  out[61] += 7.745966692414834*f[23]*w0dx0+2.0*f[63]*dv0dx0+2.23606797749979*f[9]*dv0dx0; 
  out[62] += 7.745966692414834*f[24]*w1dx1+2.0*f[67]*dv1dx1+2.23606797749979*f[8]*dv1dx1; 
  out[63] += 3.464101615137755*f[39]*w0dx0+0.8944271909999159*f[11]*dv0dx0; 
  out[64] += 3.464101615137755*f[39]*w1dx1+f[18]*dv1dx1; 
  out[65] += 3.464101615137755*f[40]*w1dx1+3.464101615137755*f[41]*w0dx0+0.8944271909999159*f[9]*dv1dx1+f[67]*dv0dx0; 
  out[66] += 3.464101615137755*f[42]*w0dx0+f[19]*dv0dx0; 
  out[67] += 3.464101615137755*f[42]*w1dx1+0.8944271909999159*f[11]*dv1dx1; 
  out[68] += 3.464101615137755*f[43]*w1dx1+7.745966692414834*f[25]*w0dx0+f[74]*dv1dx1+2.23606797749979*f[52]*dv0dx0; 
  out[69] += 7.745966692414834*f[25]*w1dx1+3.464101615137755*f[44]*w0dx0+2.23606797749979*f[53]*dv1dx1+f[71]*dv0dx0; 
  out[70] += 7.745966692414834*f[26]*w0dx0+2.0*f[72]*dv0dx0+2.23606797749979*f[12]*dv0dx0; 
  out[71] += 7.745966692414834*f[27]*w1dx1+2.23606797749979*f[55]*dv1dx1; 
  out[72] += 3.464101615137755*f[45]*w0dx0+0.8944271909999159*f[14]*dv0dx0; 
  out[73] += 3.464101615137755*f[45]*w1dx1+f[76]*dv1dx1; 
  out[74] += 7.745966692414834*f[28]*w0dx0+2.23606797749979*f[54]*dv0dx0; 
  out[75] += 7.745966692414834*f[29]*w1dx1+2.0*f[78]*dv1dx1+2.23606797749979*f[13]*dv1dx1; 
  out[77] += 3.464101615137755*f[46]*w0dx0+f[79]*dv0dx0; 
  out[78] += 3.464101615137755*f[46]*w1dx1+0.8944271909999159*f[15]*dv1dx1; 
  out[80] += 3.464101615137755*f[47]*w1dx1+3.464101615137755*f[48]*w0dx0+f[83]*dv1dx1+f[82]*dv0dx0; 
  out[81] += 3.464101615137755*f[49]*w0dx0+f[20]*dv0dx0; 
  out[82] += 3.464101615137755*f[49]*w1dx1+f[85]*dv1dx1; 
  out[83] += 3.464101615137755*f[50]*w0dx0+f[85]*dv0dx0; 
  out[84] += 3.464101615137755*f[50]*w1dx1+f[20]*dv1dx1; 
  out[86] += 3.464101615137754*f[54]*w1dx1+3.464101615137754*f[55]*w0dx0+0.8944271909999161*f[101]*dv1dx1+f[26]*dv1dx1+0.8944271909999161*f[99]*dv0dx0+f[29]*dv0dx0; 
  out[87] += 3.464101615137755*f[61]*w1dx1+7.745966692414834*f[51]*w0dx0+f[33]*dv1dx1+2.0*f[89]*dv0dx0+2.23606797749979*f[22]*dv0dx0; 
  out[88] += 7.745966692414834*f[51]*w1dx1+3.464101615137755*f[62]*w0dx0+2.0*f[90]*dv1dx1+2.23606797749979*f[21]*dv1dx1+f[38]*dv0dx0; 
  out[89] += 3.464101615137755*f[63]*w1dx1+3.464101615137755*f[64]*w0dx0+f[35]*dv1dx1+0.8944271909999161*f[24]*dv0dx0; 
  out[90] += 3.464101615137755*f[66]*w1dx1+3.464101615137755*f[67]*w0dx0+0.8944271909999161*f[23]*dv1dx1+f[41]*dv0dx0; 
  out[91] += 3.464101615137755*f[70]*w1dx1+7.745966692414834*f[52]*w0dx0+f[96]*dv1dx1+2.0*f[93]*dv0dx0+2.23606797749979*f[25]*dv0dx0; 
  out[92] += 7.745966692414834*f[52]*w1dx1+3.464101615137755*f[71]*w0dx0+2.23606797749979*f[86]*dv1dx1+f[44]*dv0dx0; 
  out[93] += 3.464101615137755*f[72]*w1dx1+3.464101615137755*f[73]*w0dx0+f[98]*dv1dx1+0.8944271909999161*f[27]*dv0dx0; 
  out[94] += 3.464101615137755*f[74]*w1dx1+7.745966692414834*f[53]*w0dx0+f[43]*dv1dx1+2.23606797749979*f[86]*dv0dx0; 
  out[95] += 7.745966692414834*f[53]*w1dx1+3.464101615137755*f[75]*w0dx0+2.0*f[100]*dv1dx1+2.23606797749979*f[25]*dv1dx1+f[97]*dv0dx0; 
  out[96] += 7.745966692414834*f[54]*w0dx0+2.0*f[98]*dv0dx0+2.23606797749979*f[28]*dv0dx0; 
  out[97] += 7.745966692414834*f[55]*w1dx1+2.0*f[102]*dv1dx1+2.23606797749979*f[27]*dv1dx1; 
  out[98] += 3.464101615137755*f[76]*w0dx0+0.8944271909999161*f[30]*dv0dx0; 
  out[99] += 3.464101615137755*f[76]*w1dx1+f[45]*dv1dx1; 
  out[100] += 3.464101615137755*f[77]*w1dx1+3.464101615137755*f[78]*w0dx0+0.8944271909999161*f[28]*dv1dx1+f[102]*dv0dx0; 
  out[101] += 3.464101615137755*f[79]*w0dx0+f[46]*dv0dx0; 
  out[102] += 3.464101615137755*f[79]*w1dx1+0.8944271909999161*f[30]*dv1dx1; 
  out[103] += 3.464101615137755*f[81]*w1dx1+3.464101615137755*f[82]*w0dx0+f[105]*dv1dx1+f[48]*dv0dx0; 
  out[104] += 3.464101615137755*f[83]*w1dx1+3.464101615137755*f[84]*w0dx0+f[47]*dv1dx1+f[106]*dv0dx0; 
  out[105] += 3.464101615137755*f[85]*w0dx0+f[50]*dv0dx0; 
  out[106] += 3.464101615137755*f[85]*w1dx1+f[49]*dv1dx1; 
  out[107] += 3.464101615137755*f[96]*w1dx1+7.745966692414834*f[86]*w0dx0+f[70]*dv1dx1+2.0*f[109]*dv0dx0+2.23606797749979*f[53]*dv0dx0; 
  out[108] += 7.745966692414834*f[86]*w1dx1+3.464101615137755*f[97]*w0dx0+2.0*f[110]*dv1dx1+2.23606797749979*f[52]*dv1dx1+f[75]*dv0dx0; 
  out[109] += 3.464101615137755*f[98]*w1dx1+3.464101615137755*f[99]*w0dx0+f[72]*dv1dx1+0.8944271909999159*f[55]*dv0dx0; 
  out[110] += 3.464101615137755*f[101]*w1dx1+3.464101615137755*f[102]*w0dx0+0.8944271909999159*f[54]*dv1dx1+f[78]*dv0dx0; 
  out[111] += 3.464101615137755*f[105]*w1dx1+3.464101615137755*f[106]*w0dx0+f[81]*dv1dx1+f[84]*dv0dx0; 
return std::abs(w0dx0)+std::abs(w1dx1); 
} 
