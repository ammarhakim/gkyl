#include <VlasovModDecl.h> 
void VlasovVolStream2x2vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[5] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[8]*dv1dx1+f[7]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+f[0]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w1dx1+f[10]*dv1dx1; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[10]*dv0dx0; 
  out[9] += 3.464101615137754*f[4]*w1dx1+f[0]*dv1dx1; 
  out[11] += 3.464101615137754*f[6]*w1dx1+3.464101615137754*f[7]*w0dx0+f[13]*dv1dx1+f[2]*dv0dx0; 
  out[12] += 3.464101615137754*f[8]*w1dx1+3.464101615137754*f[9]*w0dx0+f[1]*dv1dx1+f[14]*dv0dx0; 
  out[13] += 3.464101615137754*f[10]*w0dx0+f[4]*dv0dx0; 
  out[14] += 3.464101615137754*f[10]*w1dx1+f[3]*dv1dx1; 
  out[15] += 3.464101615137754*f[13]*w1dx1+3.464101615137754*f[14]*w0dx0+f[6]*dv1dx1+f[9]*dv0dx0; 
} 
void VlasovVolStream2x2vSerP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[3]*dv0dx0; 
  out[2] += 3.464101615137754*f[0]*w1dx1+f[4]*dv1dx1; 
  out[5] += 3.464101615137754*f[1]*w1dx1+3.464101615137754*f[2]*w0dx0+f[8]*dv1dx1+f[7]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+0.8944271909999159*f[13]*dv0dx0+f[0]*dv0dx0; 
  out[7] += 3.464101615137754*f[3]*w1dx1+f[10]*dv1dx1; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[10]*dv0dx0; 
  out[9] += 3.464101615137754*f[4]*w1dx1+0.8944271909999159*f[14]*dv1dx1+f[0]*dv1dx1; 
  out[11] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[6]*dv0dx0; 
  out[12] += 7.745966692414834*f[2]*w1dx1+2.23606797749979*f[9]*dv1dx1; 
  out[15] += 3.464101615137754*f[6]*w1dx1+3.464101615137754*f[7]*w0dx0+f[17]*dv1dx1+0.8944271909999159*f[24]*dv0dx0+f[2]*dv0dx0; 
  out[16] += 3.464101615137754*f[8]*w1dx1+3.464101615137754*f[9]*w0dx0+0.8944271909999159*f[28]*dv1dx1+f[1]*dv1dx1+f[18]*dv0dx0; 
  out[17] += 3.464101615137754*f[10]*w0dx0+0.8944271909999159*f[27]*dv0dx0+f[4]*dv0dx0; 
  out[18] += 3.464101615137754*f[10]*w1dx1+0.8944271909999159*f[30]*dv1dx1+f[3]*dv1dx1; 
  out[19] += 3.464101615137754*f[11]*w1dx1+7.745966692414834*f[5]*w0dx0+f[25]*dv1dx1+2.23606797749979*f[15]*dv0dx0; 
  out[20] += 7.745966692414834*f[5]*w1dx1+3.464101615137754*f[12]*w0dx0+2.23606797749979*f[16]*dv1dx1+f[22]*dv0dx0; 
  out[21] += 7.745966692414834*f[6]*w0dx0+2.0*f[23]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[22] += 7.745966692414834*f[7]*w1dx1+2.23606797749979*f[18]*dv1dx1; 
  out[23] += 3.464101615137754*f[13]*w0dx0+0.8944271909999159*f[3]*dv0dx0; 
  out[24] += 3.464101615137754*f[13]*w1dx1+f[27]*dv1dx1; 
  out[25] += 7.745966692414834*f[8]*w0dx0+2.23606797749979*f[17]*dv0dx0; 
  out[26] += 7.745966692414834*f[9]*w1dx1+2.0*f[29]*dv1dx1+2.23606797749979*f[2]*dv1dx1; 
  out[28] += 3.464101615137754*f[14]*w0dx0+f[30]*dv0dx0; 
  out[29] += 3.464101615137754*f[14]*w1dx1+0.8944271909999159*f[4]*dv1dx1; 
  out[31] += 3.464101615137754*f[17]*w1dx1+3.464101615137754*f[18]*w0dx0+0.8944271909999159*f[42]*dv1dx1+f[6]*dv1dx1+0.8944271909999159*f[40]*dv0dx0+f[9]*dv0dx0; 
  out[32] += 3.464101615137754*f[21]*w1dx1+7.745966692414834*f[15]*w0dx0+f[37]*dv1dx1+2.0*f[34]*dv0dx0+2.23606797749979*f[5]*dv0dx0; 
  out[33] += 7.745966692414834*f[15]*w1dx1+3.464101615137754*f[22]*w0dx0+2.23606797749979*f[31]*dv1dx1+f[12]*dv0dx0; 
  out[34] += 3.464101615137754*f[23]*w1dx1+3.464101615137754*f[24]*w0dx0+f[39]*dv1dx1+0.8944271909999159*f[7]*dv0dx0; 
  out[35] += 3.464101615137754*f[25]*w1dx1+7.745966692414834*f[16]*w0dx0+f[11]*dv1dx1+2.23606797749979*f[31]*dv0dx0; 
  out[36] += 7.745966692414834*f[16]*w1dx1+3.464101615137754*f[26]*w0dx0+2.0*f[41]*dv1dx1+2.23606797749979*f[5]*dv1dx1+f[38]*dv0dx0; 
  out[37] += 7.745966692414834*f[17]*w0dx0+2.0*f[39]*dv0dx0+2.23606797749979*f[8]*dv0dx0; 
  out[38] += 7.745966692414834*f[18]*w1dx1+2.0*f[43]*dv1dx1+2.23606797749979*f[7]*dv1dx1; 
  out[39] += 3.464101615137754*f[27]*w0dx0+0.8944271909999159*f[10]*dv0dx0; 
  out[40] += 3.464101615137754*f[27]*w1dx1+f[13]*dv1dx1; 
  out[41] += 3.464101615137754*f[28]*w1dx1+3.464101615137754*f[29]*w0dx0+0.8944271909999159*f[8]*dv1dx1+f[43]*dv0dx0; 
  out[42] += 3.464101615137754*f[30]*w0dx0+f[14]*dv0dx0; 
  out[43] += 3.464101615137754*f[30]*w1dx1+0.8944271909999159*f[10]*dv1dx1; 
  out[44] += 3.464101615137754*f[37]*w1dx1+7.745966692414834*f[31]*w0dx0+f[21]*dv1dx1+2.0*f[46]*dv0dx0+2.23606797749979*f[16]*dv0dx0; 
  out[45] += 7.745966692414834*f[31]*w1dx1+3.464101615137754*f[38]*w0dx0+2.0*f[47]*dv1dx1+2.23606797749979*f[15]*dv1dx1+f[26]*dv0dx0; 
  out[46] += 3.464101615137754*f[39]*w1dx1+3.464101615137754*f[40]*w0dx0+f[23]*dv1dx1+0.8944271909999159*f[18]*dv0dx0; 
  out[47] += 3.464101615137754*f[42]*w1dx1+3.464101615137754*f[43]*w0dx0+0.8944271909999159*f[17]*dv1dx1+f[29]*dv0dx0; 
} 
