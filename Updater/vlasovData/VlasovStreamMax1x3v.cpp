void VlasovVolStream1x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
} 
void VlasovVolStream1x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  out[1] += 3.464101615137754*f[0]*w0dx0+f[2]*dv0dx0; 
  out[5] += 3.464101615137754*f[2]*w0dx0+0.8944271909999159*f[12]*dv0dx0+f[0]*dv0dx0; 
  out[6] += 3.464101615137754*f[3]*w0dx0+f[7]*dv0dx0; 
  out[8] += 3.464101615137754*f[4]*w0dx0+f[9]*dv0dx0; 
  out[11] += 7.745966692414834*f[1]*w0dx0+2.23606797749979*f[5]*dv0dx0; 
} 
