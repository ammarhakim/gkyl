void VlasovVolStream1x2vSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+f[0]*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 3.464101615137754*f[6]*w0dx0+f[3]*dv0dx0; 
} 
void VlasovVolStream1x2vSerP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[4] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[8]*dv0dx0+f[0]*dv0dx0; 
  out[5] += 3.464101615137754*f[3]*w0dx0+f[6]*dv0dx0; 
  out[7] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[4]*dv0dx0; 
  out[10] += 3.464101615137754*f[6]*w0dx0+0.8944271909999159*f[14]*dv0dx0+f[3]*dv0dx0; 
  out[11] += 7.745966692414834*f[4]*w0dx0+2.0*f[12]*dv0dx0+2.23606797749979*f[1]*dv0dx0; 
  out[12] += 3.464101615137754*f[8]*w0dx0+0.8944271909999159*f[2]*dv0dx0; 
  out[13] += 7.745966692414834*f[5]*w0dx0+2.23606797749979*f[10]*dv0dx0; 
  out[15] += 3.464101615137754*f[9]*w0dx0+f[16]*dv0dx0; 
  out[17] += 7.745966692414834*f[10]*w0dx0+2.0*f[18]*dv0dx0+2.23606797749979*f[5]*dv0dx0; 
  out[18] += 3.464101615137754*f[14]*w0dx0+0.8944271909999159*f[6]*dv0dx0; 
  out[19] += 3.464101615137754*f[16]*w0dx0+f[9]*dv0dx0; 
} 
