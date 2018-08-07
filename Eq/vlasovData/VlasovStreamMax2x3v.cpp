#include <VlasovModDecl.h> 
double VlasovVolStream2x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  double alpha0[6]; 

  double alpha1[6]; 

  alpha0[0] = 11.31370849898477*w0dx0; 
  alpha0[3] = 3.265986323710906*dv0dx0; 

  alpha1[0] = 11.31370849898477*w1dx1; 
  alpha1[4] = 3.265986323710906*dv1dx1; 

  out[1] += 0.3061862178478971*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.3061862178478971*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1); 
} 
double VlasovVolStream2x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  double alpha0[21]; 

  double alpha1[21]; 

  alpha0[0] = 11.31370849898477*w0dx0; 
  alpha0[3] = 3.265986323710906*dv0dx0; 

  alpha1[0] = 11.31370849898477*w1dx1; 
  alpha1[4] = 3.265986323710906*dv1dx1; 

  out[1] += 0.3061862178478971*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.3061862178478971*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[6] += 0.3061862178478971*(alpha1[4]*f[9]+alpha0[3]*f[8]+alpha0[0]*f[2]+alpha1[0]*f[1]); 
  out[7] += 0.273861278752583*alpha0[3]*f[18]+0.3061862178478971*(alpha0[0]*f[3]+f[0]*alpha0[3]); 
  out[8] += 0.3061862178478971*(alpha1[4]*f[11]+alpha1[0]*f[3]); 
  out[9] += 0.3061862178478971*(alpha0[3]*f[11]+alpha0[0]*f[4]); 
  out[10] += 0.273861278752583*alpha1[4]*f[19]+0.3061862178478971*(alpha1[0]*f[4]+f[0]*alpha1[4]); 
  out[12] += 0.3061862178478971*(alpha0[3]*f[14]+alpha0[0]*f[5]); 
  out[13] += 0.3061862178478971*(alpha1[4]*f[15]+alpha1[0]*f[5]); 
  out[16] += 0.6846531968814573*(alpha0[3]*f[7]+alpha0[0]*f[1]); 
  out[17] += 0.6846531968814573*(alpha1[4]*f[10]+alpha1[0]*f[2]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1); 
} 
double VlasovVolStream2x3vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[2]/dxv[0]; 
  double w0dx0 = w[2]/dxv[0]; 
  double dv1dx1 = dxv[3]/dxv[1]; 
  double w1dx1 = w[3]/dxv[1]; 
  double alpha0[56]; 

  double alpha1[56]; 

  alpha0[0] = 11.31370849898477*w0dx0; 
  alpha0[3] = 3.265986323710906*dv0dx0; 

  alpha1[0] = 11.31370849898477*w1dx1; 
  alpha1[4] = 3.265986323710906*dv1dx1; 

  out[1] += 0.3061862178478971*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[2] += 0.3061862178478971*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
  out[6] += 0.3061862178478971*(alpha1[4]*f[9]+alpha0[3]*f[8]+alpha0[0]*f[2]+alpha1[0]*f[1]); 
  out[7] += 0.273861278752583*alpha0[3]*f[18]+0.3061862178478971*(alpha0[0]*f[3]+f[0]*alpha0[3]); 
  out[8] += 0.3061862178478971*(alpha1[4]*f[11]+alpha1[0]*f[3]); 
  out[9] += 0.3061862178478971*(alpha0[3]*f[11]+alpha0[0]*f[4]); 
  out[10] += 0.273861278752583*alpha1[4]*f[19]+0.3061862178478971*(alpha1[0]*f[4]+f[0]*alpha1[4]); 
  out[12] += 0.3061862178478971*(alpha0[3]*f[14]+alpha0[0]*f[5]); 
  out[13] += 0.3061862178478971*(alpha1[4]*f[15]+alpha1[0]*f[5]); 
  out[16] += 0.6846531968814573*(alpha0[3]*f[7]+alpha0[0]*f[1]); 
  out[17] += 0.6846531968814573*(alpha1[4]*f[10]+alpha1[0]*f[2]); 
  out[21] += 0.273861278752583*alpha0[3]*f[36]+0.3061862178478971*(alpha1[4]*f[23]+alpha0[0]*f[8]+alpha1[0]*f[7]+f[2]*alpha0[3]); 
  out[22] += 0.273861278752583*alpha1[4]*f[40]+0.3061862178478971*(alpha0[3]*f[24]+alpha0[0]*f[10]+alpha1[0]*f[9]+f[1]*alpha1[4]); 
  out[23] += 0.273861278752583*alpha0[3]*f[39]+0.3061862178478971*(alpha0[0]*f[11]+alpha0[3]*f[4]); 
  out[24] += 0.273861278752583*alpha1[4]*f[42]+0.3061862178478971*(alpha1[0]*f[11]+f[3]*alpha1[4]); 
  out[25] += 0.3061862178478971*(alpha1[4]*f[28]+alpha0[3]*f[27]+alpha0[0]*f[13]+alpha1[0]*f[12]); 
  out[26] += 0.273861278752583*alpha0[3]*f[45]+0.3061862178478971*(alpha0[0]*f[14]+alpha0[3]*f[5]); 
  out[27] += 0.3061862178478971*(alpha1[4]*f[30]+alpha1[0]*f[14]); 
  out[28] += 0.3061862178478971*(alpha0[3]*f[30]+alpha0[0]*f[15]); 
  out[29] += 0.273861278752583*alpha1[4]*f[46]+0.3061862178478971*(alpha1[0]*f[15]+alpha1[4]*f[5]); 
  out[31] += 0.3061862178478971*alpha1[4]*f[37]+0.6846531968814573*alpha0[3]*f[21]+0.3061862178478971*alpha1[0]*f[16]+0.6846531968814573*alpha0[0]*f[6]; 
  out[32] += 0.3061862178478971*alpha0[3]*f[34]+0.6846531968814573*alpha1[4]*f[22]+0.3061862178478971*alpha0[0]*f[17]+0.6846531968814573*alpha1[0]*f[6]; 
  out[33] += 0.6123724356957944*alpha0[3]*f[35]+0.6846531968814573*(alpha0[0]*f[7]+f[1]*alpha0[3]); 
  out[34] += 0.6846531968814573*(alpha1[4]*f[24]+alpha1[0]*f[8]); 
  out[35] += 0.2689264371002384*alpha0[3]*f[53]+0.3061862178478971*alpha0[0]*f[18]+0.273861278752583*alpha0[3]*f[3]; 
  out[36] += 0.3061862178478971*(alpha1[4]*f[39]+alpha1[0]*f[18]); 
  out[37] += 0.6846531968814573*(alpha0[3]*f[23]+alpha0[0]*f[9]); 
  out[38] += 0.6123724356957944*alpha1[4]*f[41]+0.6846531968814573*(alpha1[0]*f[10]+f[2]*alpha1[4]); 
  out[40] += 0.3061862178478971*(alpha0[3]*f[42]+alpha0[0]*f[19]); 
  out[41] += 0.2689264371002384*alpha1[4]*f[54]+0.3061862178478971*alpha1[0]*f[19]+0.273861278752583*alpha1[4]*f[4]; 
  out[43] += 0.6846531968814573*(alpha0[3]*f[26]+alpha0[0]*f[12]); 
  out[44] += 0.6846531968814573*(alpha1[4]*f[29]+alpha1[0]*f[13]); 
  out[47] += 0.3061862178478971*(alpha0[3]*f[49]+alpha0[0]*f[20]); 
  out[48] += 0.3061862178478971*(alpha1[4]*f[50]+alpha1[0]*f[20]); 
  out[51] += 1.045825033167594*(alpha0[3]*f[33]+alpha0[0]*f[16])+0.4677071733467425*(alpha0[3]*f[3]+alpha0[0]*f[0]); 
  out[52] += 1.045825033167594*(alpha1[4]*f[38]+alpha1[0]*f[17])+0.4677071733467425*(alpha1[4]*f[4]+alpha1[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+0.5*(dv0dx0+dv1dx1); 
} 
