#include <VlasovModDecl.h> 
double VlasovVolStream3x3vMaxP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  double alpha0[7]; 

  double alpha1[7]; 

  double alpha2[7]; 

  alpha0[0] = 16.0*w0dx0; 
  alpha0[4] = 4.618802153517007*dv0dx0; 

  alpha1[0] = 16.0*w1dx1; 
  alpha1[5] = 4.618802153517007*dv1dx1; 

  alpha2[0] = 16.0*w2dx2; 
  alpha2[6] = 4.618802153517007*dv2dx2; 

  out[1] += 0.2165063509461096*(alpha0[4]*f[4]+alpha0[0]*f[0]); 
  out[2] += 0.2165063509461096*(alpha1[5]*f[5]+alpha1[0]*f[0]); 
  out[3] += 0.2165063509461096*(alpha2[6]*f[6]+alpha2[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+std::abs(w2dx2)+0.5*(dv0dx0+dv1dx1+dv2dx2); 
} 
double VlasovVolStream3x3vMaxP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  double alpha0[28]; 

  double alpha1[28]; 

  double alpha2[28]; 

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
  out[10] += 0.1936491673103708*alpha0[4]*f[25]+0.2165063509461096*(alpha0[0]*f[4]+f[0]*alpha0[4]); 
  out[11] += 0.2165063509461096*(alpha1[5]*f[16]+alpha1[0]*f[4]); 
  out[12] += 0.2165063509461096*(alpha2[6]*f[20]+alpha2[0]*f[4]); 
  out[13] += 0.2165063509461096*(alpha0[4]*f[16]+alpha0[0]*f[5]); 
  out[14] += 0.1936491673103708*alpha1[5]*f[26]+0.2165063509461096*(alpha1[0]*f[5]+f[0]*alpha1[5]); 
  out[15] += 0.2165063509461096*(alpha2[6]*f[21]+alpha2[0]*f[5]); 
  out[17] += 0.2165063509461096*(alpha0[4]*f[20]+alpha0[0]*f[6]); 
  out[18] += 0.2165063509461096*(alpha1[5]*f[21]+alpha1[0]*f[6]); 
  out[19] += 0.1936491673103708*alpha2[6]*f[27]+0.2165063509461096*(alpha2[0]*f[6]+f[0]*alpha2[6]); 
  out[22] += 0.4841229182759271*(alpha0[4]*f[10]+alpha0[0]*f[1]); 
  out[23] += 0.4841229182759271*(alpha1[5]*f[14]+alpha1[0]*f[2]); 
  out[24] += 0.4841229182759271*(alpha2[6]*f[19]+alpha2[0]*f[3]); 
return std::abs(w0dx0)+std::abs(w1dx1)+std::abs(w2dx2)+0.5*(dv0dx0+dv1dx1+dv2dx2); 
} 
double VlasovVolStream3x3vMaxP3(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. f: Input distribution function. out: Incremented output 
  double dv0dx0 = dxv[3]/dxv[0]; 
  double w0dx0 = w[3]/dxv[0]; 
  double dv1dx1 = dxv[4]/dxv[1]; 
  double w1dx1 = w[4]/dxv[1]; 
  double dv2dx2 = dxv[5]/dxv[2]; 
  double w2dx2 = w[5]/dxv[2]; 
  double alpha0[84]; 

  double alpha1[84]; 

  double alpha2[84]; 

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
  out[10] += 0.1936491673103708*alpha0[4]*f[25]+0.2165063509461096*(alpha0[0]*f[4]+f[0]*alpha0[4]); 
  out[11] += 0.2165063509461096*(alpha1[5]*f[16]+alpha1[0]*f[4]); 
  out[12] += 0.2165063509461096*(alpha2[6]*f[20]+alpha2[0]*f[4]); 
  out[13] += 0.2165063509461096*(alpha0[4]*f[16]+alpha0[0]*f[5]); 
  out[14] += 0.1936491673103708*alpha1[5]*f[26]+0.2165063509461096*(alpha1[0]*f[5]+f[0]*alpha1[5]); 
  out[15] += 0.2165063509461096*(alpha2[6]*f[21]+alpha2[0]*f[5]); 
  out[17] += 0.2165063509461096*(alpha0[4]*f[20]+alpha0[0]*f[6]); 
  out[18] += 0.2165063509461096*(alpha1[5]*f[21]+alpha1[0]*f[6]); 
  out[19] += 0.1936491673103708*alpha2[6]*f[27]+0.2165063509461096*(alpha2[0]*f[6]+f[0]*alpha2[6]); 
  out[22] += 0.4841229182759271*(alpha0[4]*f[10]+alpha0[0]*f[1]); 
  out[23] += 0.4841229182759271*(alpha1[5]*f[14]+alpha1[0]*f[2]); 
  out[24] += 0.4841229182759271*(alpha2[6]*f[19]+alpha2[0]*f[3]); 
  out[28] += 0.2165063509461096*(alpha2[6]*f[38]+alpha1[5]*f[33]+alpha0[4]*f[31]+alpha0[0]*f[9]+alpha1[0]*f[8]+alpha2[0]*f[7]); 
  out[29] += 0.1936491673103708*alpha0[4]*f[58]+0.2165063509461096*(alpha1[5]*f[35]+alpha0[0]*f[11]+alpha1[0]*f[10]+f[2]*alpha0[4]); 
  out[30] += 0.1936491673103708*alpha0[4]*f[59]+0.2165063509461096*(alpha2[6]*f[41]+alpha0[0]*f[12]+alpha2[0]*f[10]+f[3]*alpha0[4]); 
  out[31] += 0.2165063509461096*(alpha2[6]*f[42]+alpha1[5]*f[37]+alpha1[0]*f[12]+alpha2[0]*f[11]); 
  out[32] += 0.1936491673103708*alpha1[5]*f[64]+0.2165063509461096*(alpha0[4]*f[36]+alpha0[0]*f[14]+alpha1[0]*f[13]+f[1]*alpha1[5]); 
  out[33] += 0.2165063509461096*(alpha2[6]*f[44]+alpha0[4]*f[37]+alpha0[0]*f[15]+alpha2[0]*f[13]); 
  out[34] += 0.1936491673103708*alpha1[5]*f[66]+0.2165063509461096*(alpha2[6]*f[45]+alpha1[0]*f[15]+alpha2[0]*f[14]+f[3]*alpha1[5]); 
  out[35] += 0.1936491673103708*alpha0[4]*f[63]+0.2165063509461096*(alpha0[0]*f[16]+alpha0[4]*f[5]); 
  out[36] += 0.1936491673103708*alpha1[5]*f[67]+0.2165063509461096*(alpha1[0]*f[16]+f[4]*alpha1[5]); 
  out[37] += 0.2165063509461096*(alpha2[6]*f[47]+alpha2[0]*f[16]); 
  out[38] += 0.2165063509461096*(alpha1[5]*f[44]+alpha0[4]*f[42]+alpha0[0]*f[18]+alpha1[0]*f[17]); 
  out[39] += 0.1936491673103708*alpha2[6]*f[73]+0.2165063509461096*(alpha0[4]*f[43]+alpha0[0]*f[19]+alpha2[0]*f[17]+f[1]*alpha2[6]); 
  out[40] += 0.1936491673103708*alpha2[6]*f[74]+0.2165063509461096*(alpha1[5]*f[46]+alpha1[0]*f[19]+alpha2[0]*f[18]+f[2]*alpha2[6]); 
  out[41] += 0.1936491673103708*alpha0[4]*f[71]+0.2165063509461096*(alpha0[0]*f[20]+alpha0[4]*f[6]); 
  out[42] += 0.2165063509461096*(alpha1[5]*f[47]+alpha1[0]*f[20]); 
  out[43] += 0.1936491673103708*alpha2[6]*f[76]+0.2165063509461096*(alpha2[0]*f[20]+f[4]*alpha2[6]); 
  out[44] += 0.2165063509461096*(alpha0[4]*f[47]+alpha0[0]*f[21]); 
  out[45] += 0.1936491673103708*alpha1[5]*f[72]+0.2165063509461096*(alpha1[0]*f[21]+alpha1[5]*f[6]); 
  out[46] += 0.1936491673103708*alpha2[6]*f[77]+0.2165063509461096*(alpha2[0]*f[21]+f[5]*alpha2[6]); 
  out[48] += 0.2165063509461096*alpha1[5]*f[60]+0.4841229182759271*alpha0[4]*f[29]+0.2165063509461096*alpha1[0]*f[22]+0.4841229182759271*alpha0[0]*f[7]; 
  out[49] += 0.2165063509461096*alpha0[4]*f[55]+0.4841229182759271*alpha1[5]*f[32]+0.2165063509461096*alpha0[0]*f[23]+0.4841229182759271*alpha1[0]*f[7]; 
  out[50] += 0.2165063509461096*alpha2[6]*f[68]+0.4841229182759271*alpha0[4]*f[30]+0.2165063509461096*alpha2[0]*f[22]+0.4841229182759271*alpha0[0]*f[8]; 
  out[51] += 0.2165063509461096*alpha2[6]*f[69]+0.4841229182759271*alpha1[5]*f[34]+0.2165063509461096*alpha2[0]*f[23]+0.4841229182759271*alpha1[0]*f[9]; 
  out[52] += 0.2165063509461096*alpha0[4]*f[56]+0.4841229182759271*alpha2[6]*f[39]+0.2165063509461096*alpha0[0]*f[24]+0.4841229182759271*alpha2[0]*f[8]; 
  out[53] += 0.2165063509461096*alpha1[5]*f[62]+0.4841229182759271*alpha2[6]*f[40]+0.2165063509461096*alpha1[0]*f[24]+0.4841229182759271*alpha2[0]*f[9]; 
  out[54] += 0.4330127018922193*alpha0[4]*f[57]+0.4841229182759271*(alpha0[0]*f[10]+f[1]*alpha0[4]); 
  out[55] += 0.4841229182759271*(alpha1[5]*f[36]+alpha1[0]*f[11]); 
  out[56] += 0.4841229182759271*(alpha2[6]*f[43]+alpha2[0]*f[12]); 
  out[57] += 0.1901597073139162*alpha0[4]*f[81]+0.2165063509461096*alpha0[0]*f[25]+0.1936491673103708*alpha0[4]*f[4]; 
  out[58] += 0.2165063509461096*(alpha1[5]*f[63]+alpha1[0]*f[25]); 
  out[59] += 0.2165063509461096*(alpha2[6]*f[71]+alpha2[0]*f[25]); 
  out[60] += 0.4841229182759271*(alpha0[4]*f[35]+alpha0[0]*f[13]); 
  out[61] += 0.4330127018922193*alpha1[5]*f[65]+0.4841229182759271*(alpha1[0]*f[14]+f[2]*alpha1[5]); 
  out[62] += 0.4841229182759271*(alpha2[6]*f[46]+alpha2[0]*f[15]); 
  out[64] += 0.2165063509461096*(alpha0[4]*f[67]+alpha0[0]*f[26]); 
  out[65] += 0.1901597073139162*alpha1[5]*f[82]+0.2165063509461096*alpha1[0]*f[26]+0.1936491673103708*alpha1[5]*f[5]; 
  out[66] += 0.2165063509461096*(alpha2[6]*f[72]+alpha2[0]*f[26]); 
  out[68] += 0.4841229182759271*(alpha0[4]*f[41]+alpha0[0]*f[17]); 
  out[69] += 0.4841229182759271*(alpha1[5]*f[45]+alpha1[0]*f[18]); 
  out[70] += 0.4330127018922193*alpha2[6]*f[75]+0.4841229182759271*(alpha2[0]*f[19]+f[3]*alpha2[6]); 
  out[73] += 0.2165063509461096*(alpha0[4]*f[76]+alpha0[0]*f[27]); 
  out[74] += 0.2165063509461096*(alpha1[5]*f[77]+alpha1[0]*f[27]); 
  out[75] += 0.1901597073139162*alpha2[6]*f[83]+0.2165063509461096*alpha2[0]*f[27]+0.1936491673103708*alpha2[6]*f[6]; 
  out[78] += 0.7395099728874521*(alpha0[4]*f[54]+alpha0[0]*f[22])+0.3307189138830738*(alpha0[4]*f[4]+alpha0[0]*f[0]); 
  out[79] += 0.7395099728874521*(alpha1[5]*f[61]+alpha1[0]*f[23])+0.3307189138830738*(alpha1[5]*f[5]+alpha1[0]*f[0]); 
  out[80] += 0.7395099728874521*(alpha2[6]*f[70]+alpha2[0]*f[24])+0.3307189138830738*(alpha2[6]*f[6]+alpha2[0]*f[0]); 
return std::abs(w0dx0)+std::abs(w1dx1)+std::abs(w2dx2)+0.5*(dv0dx0+dv1dx1+dv2dx2); 
} 
