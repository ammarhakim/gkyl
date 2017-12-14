#include <VlasovModDecl.h> 
void VlasovVolStream1x3vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[5] += 3.464101615137754*f[2]*w0dx0+f[0]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+f[7]*dv0dx0; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[9]*dv0dx0; 
  out[11] += 3.464101615137754*f[7]*w0dx0+f[3]*dv0dx0; 
  out[12] += 3.464101615137754*f[9]*w0dx0+f[4]*dv0dx0; 
  out[13] += 3.464101615137754*f[10]*w0dx0+f[14]*dv0dx0; 
  out[15] += 3.464101615137754*f[14]*w0dx0+f[10]*dv0dx0; 
} 
void VlasovVolStream1x3vSerP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[5] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[12]*dv0dx0+f[0]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+f[7]*dv0dx0; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[9]*dv0dx0; 
  out[11] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[5]*dv0dx0; 
  out[15] += 3.464101615137754*f[7]*w0dx0+0.8944271909999159*f[22]*dv0dx0+f[3]*dv0dx0; 
  out[16] += 3.464101615137754*f[9]*w0dx0+0.8944271909999159*f[26]*dv0dx0+f[4]*dv0dx0; 
  out[17] += 3.464101615137754*f[10]*w0dx0+f[18]*dv0dx0; 
  out[19] += 7.745966692414834*f[5]*w0dx0+2.0*f[20]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[20] += 3.464101615137754*f[12]*w0dx0+0.8944271909999159*f[2]*dv0dx0; 
  out[21] += 7.745966692414834*f[6]*w0dx0+2.23606797749979*f[15]*dv0dx0; 
  out[23] += 3.464101615137754*f[13]*w0dx0+f[24]*dv0dx0; 
  out[25] += 7.745966692414834*f[8]*w0dx0+2.23606797749979*f[16]*dv0dx0; 
  out[28] += 3.464101615137754*f[14]*w0dx0+f[29]*dv0dx0; 
  out[31] += 3.464101615137754*f[18]*w0dx0+0.8944271909999159*f[38]*dv0dx0+f[10]*dv0dx0; 
  out[32] += 7.745966692414834*f[15]*w0dx0+2.0*f[33]*dv0dx0+2.23606797749979*f[6]*dv0dx0; 
  out[33] += 3.464101615137754*f[22]*w0dx0+0.8944271909999159*f[7]*dv0dx0; 
  out[34] += 3.464101615137754*f[24]*w0dx0+f[13]*dv0dx0; 
  out[35] += 7.745966692414834*f[16]*w0dx0+2.0*f[36]*dv0dx0+2.23606797749979*f[8]*dv0dx0; 
  out[36] += 3.464101615137754*f[26]*w0dx0+0.8944271909999159*f[9]*dv0dx0; 
  out[37] += 7.745966692414834*f[17]*w0dx0+2.23606797749979*f[31]*dv0dx0; 
  out[39] += 3.464101615137754*f[27]*w0dx0+f[40]*dv0dx0; 
  out[41] += 3.464101615137754*f[29]*w0dx0+f[14]*dv0dx0; 
  out[42] += 3.464101615137754*f[30]*w0dx0+f[43]*dv0dx0; 
  out[44] += 7.745966692414834*f[31]*w0dx0+2.0*f[45]*dv0dx0+2.23606797749979*f[17]*dv0dx0; 
  out[45] += 3.464101615137754*f[38]*w0dx0+0.8944271909999159*f[18]*dv0dx0; 
  out[46] += 3.464101615137754*f[40]*w0dx0+f[27]*dv0dx0; 
  out[47] += 3.464101615137754*f[43]*w0dx0+f[30]*dv0dx0; 
} 
