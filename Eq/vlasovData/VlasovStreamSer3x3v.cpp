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
  double alpha0[64]; 

  double alpha1[64]; 

  double alpha2[64]; 

  alpha0[0] = 16.0*w0dx0; 
  alpha0[4] = 4.618802153517007*dv0dx0; 

  alpha1[0] = 16.0*w1dx1; 
  alpha1[5] = 4.618802153517007*dv1dx1; 

  alpha2[0] = 16.0*w2dx2; 
  alpha2[6] = 4.618802153517007*dv2dx2; 

  out[1] += 0.2165063509461096*(alpha0[4]*f[4]+alpha0[0]*f[0]); 
  out[2] += 0.2165063509461096*(alpha1[5]*f[5]+alpha1[0]*f[0]); 
  out[3] += 0.2165063509461096*(alpha2[6]*f[6]+alpha2[0]*f[0]); 
  out[7] += 0.2165063509461096*(alpha1[5]*f[13]+alpha0[4]*f[11]+alpha0[0]*f[2]+alpha1[0]*f[1]); 
  out[8] += 0.2165063509461096*(alpha2[6]*f[17]+alpha0[4]*f[12]+alpha0[0]*f[3]+alpha2[0]*f[1]); 
  out[9] += 0.2165063509461096*(alpha2[6]*f[18]+alpha1[5]*f[15]+alpha1[0]*f[3]+alpha2[0]*f[2]); 
  out[10] += 0.2165063509461096*(alpha0[0]*f[4]+f[0]*alpha0[4]); 
  out[11] += 0.2165063509461096*(alpha1[5]*f[16]+alpha1[0]*f[4]); 
  out[12] += 0.2165063509461096*(alpha2[6]*f[20]+alpha2[0]*f[4]); 
  out[13] += 0.2165063509461096*(alpha0[4]*f[16]+alpha0[0]*f[5]); 
  out[14] += 0.2165063509461096*(alpha1[0]*f[5]+f[0]*alpha1[5]); 
  out[15] += 0.2165063509461096*(alpha2[6]*f[21]+alpha2[0]*f[5]); 
  out[17] += 0.2165063509461096*(alpha0[4]*f[20]+alpha0[0]*f[6]); 
  out[18] += 0.2165063509461096*(alpha1[5]*f[21]+alpha1[0]*f[6]); 
  out[19] += 0.2165063509461096*(alpha2[0]*f[6]+f[0]*alpha2[6]); 
  out[22] += 0.2165063509461096*(alpha2[6]*f[32]+alpha1[5]*f[27]+alpha0[4]*f[25]+alpha0[0]*f[9]+alpha1[0]*f[8]+alpha2[0]*f[7]); 
  out[23] += 0.2165063509461096*(alpha1[5]*f[29]+alpha0[0]*f[11]+alpha1[0]*f[10]+f[2]*alpha0[4]); 
  out[24] += 0.2165063509461096*(alpha2[6]*f[35]+alpha0[0]*f[12]+alpha2[0]*f[10]+f[3]*alpha0[4]); 
  out[25] += 0.2165063509461096*(alpha2[6]*f[36]+alpha1[5]*f[31]+alpha1[0]*f[12]+alpha2[0]*f[11]); 
  out[26] += 0.2165063509461096*(alpha0[4]*f[30]+alpha0[0]*f[14]+alpha1[0]*f[13]+f[1]*alpha1[5]); 
  out[27] += 0.2165063509461096*(alpha2[6]*f[38]+alpha0[4]*f[31]+alpha0[0]*f[15]+alpha2[0]*f[13]); 
  out[28] += 0.2165063509461096*(alpha2[6]*f[39]+alpha1[0]*f[15]+alpha2[0]*f[14]+f[3]*alpha1[5]); 
  out[29] += 0.2165063509461096*(alpha0[0]*f[16]+alpha0[4]*f[5]); 
  out[30] += 0.2165063509461096*(alpha1[0]*f[16]+f[4]*alpha1[5]); 
  out[31] += 0.2165063509461096*(alpha2[6]*f[41]+alpha2[0]*f[16]); 
  out[32] += 0.2165063509461096*(alpha1[5]*f[38]+alpha0[4]*f[36]+alpha0[0]*f[18]+alpha1[0]*f[17]); 
  out[33] += 0.2165063509461096*(alpha0[4]*f[37]+alpha0[0]*f[19]+alpha2[0]*f[17]+f[1]*alpha2[6]); 
  out[34] += 0.2165063509461096*(alpha1[5]*f[40]+alpha1[0]*f[19]+alpha2[0]*f[18]+f[2]*alpha2[6]); 
  out[35] += 0.2165063509461096*(alpha0[0]*f[20]+alpha0[4]*f[6]); 
  out[36] += 0.2165063509461096*(alpha1[5]*f[41]+alpha1[0]*f[20]); 
  out[37] += 0.2165063509461096*(alpha2[0]*f[20]+f[4]*alpha2[6]); 
  out[38] += 0.2165063509461096*(alpha0[4]*f[41]+alpha0[0]*f[21]); 
  out[39] += 0.2165063509461096*(alpha1[0]*f[21]+alpha1[5]*f[6]); 
  out[40] += 0.2165063509461096*(alpha2[0]*f[21]+f[5]*alpha2[6]); 
  out[42] += 0.2165063509461096*(alpha2[6]*f[48]+alpha1[5]*f[45]+alpha0[0]*f[25]+alpha1[0]*f[24]+alpha2[0]*f[23]+alpha0[4]*f[9]); 
  out[43] += 0.2165063509461096*(alpha2[6]*f[51]+alpha0[4]*f[46]+alpha0[0]*f[28]+alpha1[0]*f[27]+alpha2[0]*f[26]+alpha1[5]*f[8]); 
  out[44] += 0.2165063509461096*(alpha0[0]*f[30]+alpha1[0]*f[29]+alpha0[4]*f[14]+alpha1[5]*f[10]); 
  out[45] += 0.2165063509461096*(alpha2[6]*f[54]+alpha0[0]*f[31]+alpha2[0]*f[29]+alpha0[4]*f[15]); 
  out[46] += 0.2165063509461096*(alpha2[6]*f[55]+alpha1[0]*f[31]+alpha2[0]*f[30]+alpha1[5]*f[12]); 
  out[47] += 0.2165063509461096*(alpha1[5]*f[52]+alpha0[4]*f[50]+alpha0[0]*f[34]+alpha1[0]*f[33]+alpha2[0]*f[32]+alpha2[6]*f[7]); 
  out[48] += 0.2165063509461096*(alpha1[5]*f[54]+alpha0[0]*f[36]+alpha1[0]*f[35]+alpha0[4]*f[18]); 
  out[49] += 0.2165063509461096*(alpha0[0]*f[37]+alpha2[0]*f[35]+alpha0[4]*f[19]+alpha2[6]*f[10]); 
  out[50] += 0.2165063509461096*(alpha1[5]*f[56]+alpha1[0]*f[37]+alpha2[0]*f[36]+alpha2[6]*f[11]); 
  out[51] += 0.2165063509461096*(alpha0[4]*f[55]+alpha0[0]*f[39]+alpha1[0]*f[38]+alpha1[5]*f[17]); 
  out[52] += 0.2165063509461096*(alpha0[4]*f[56]+alpha0[0]*f[40]+alpha2[0]*f[38]+alpha2[6]*f[13]); 
  out[53] += 0.2165063509461096*(alpha1[0]*f[40]+alpha2[0]*f[39]+alpha1[5]*f[19]+alpha2[6]*f[14]); 
  out[54] += 0.2165063509461096*(alpha0[0]*f[41]+alpha0[4]*f[21]); 
  out[55] += 0.2165063509461096*(alpha1[0]*f[41]+alpha1[5]*f[20]); 
  out[56] += 0.2165063509461096*(alpha2[0]*f[41]+alpha2[6]*f[16]); 
  out[57] += 0.2165063509461096*(alpha2[6]*f[60]+alpha0[0]*f[46]+alpha1[0]*f[45]+alpha2[0]*f[44]+alpha0[4]*f[28]+alpha1[5]*f[24]); 
  out[58] += 0.2165063509461096*(alpha1[5]*f[61]+alpha0[0]*f[50]+alpha1[0]*f[49]+alpha2[0]*f[48]+alpha0[4]*f[34]+alpha2[6]*f[23]); 
  out[59] += 0.2165063509461096*(alpha0[4]*f[62]+alpha0[0]*f[53]+alpha1[0]*f[52]+alpha2[0]*f[51]+alpha1[5]*f[33]+alpha2[6]*f[26]); 
  out[60] += 0.2165063509461096*(alpha0[0]*f[55]+alpha1[0]*f[54]+alpha0[4]*f[39]+alpha1[5]*f[35]); 
  out[61] += 0.2165063509461096*(alpha0[0]*f[56]+alpha2[0]*f[54]+alpha0[4]*f[40]+alpha2[6]*f[29]); 
  out[62] += 0.2165063509461096*(alpha1[0]*f[56]+alpha2[0]*f[55]+alpha1[5]*f[37]+alpha2[6]*f[30]); 
  out[63] += 0.2165063509461096*(alpha0[0]*f[62]+alpha1[0]*f[61]+alpha2[0]*f[60]+alpha0[4]*f[53]+alpha1[5]*f[49]+alpha2[6]*f[44]); 
return std::abs(w0dx0)+std::abs(w1dx1)+std::abs(w2dx2)+0.5*(dv0dx0+dv1dx1+dv2dx2); 
} 
