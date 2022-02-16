#include <VlasovModDecl.h> 
__host__ __device__ double VlasovNeutralVol1x3vSerP2(const double *w, const double *dxv, const double *EM, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // EM:        Input EM-field.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *Fo0 = &EM[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *Fo1 = &EM[3]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double *Fo2 = &EM[6]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[48]; 
  double alpha_vdim[144]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.828427124746191*Fo0[0]*dv10; 
  alpha_vdim[1] = 2.828427124746191*Fo0[1]*dv10; 
  alpha_vdim[11] = 2.828427124746191*Fo0[2]*dv10; 
  alpha_mid += std::abs(0.125*alpha_vdim[0]-0.1397542485937369*alpha_vdim[11]); 

  alpha_vdim[48] = 2.828427124746191*Fo1[0]*dv11; 
  alpha_vdim[49] = 2.828427124746191*Fo1[1]*dv11; 
  alpha_vdim[59] = 2.828427124746191*Fo1[2]*dv11; 
  alpha_mid += std::abs(0.125*alpha_vdim[48]-0.1397542485937369*alpha_vdim[59]); 

  alpha_vdim[96] = 2.828427124746191*Fo2[0]*dv12; 
  alpha_vdim[97] = 2.828427124746191*Fo2[1]*dv12; 
  alpha_vdim[107] = 2.828427124746191*Fo2[2]*dv12; 
  alpha_mid += std::abs(0.125*alpha_vdim[96]-0.1397542485937369*alpha_vdim[107]); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[11]*f[11]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[11]*alpha_vdim[59]+f[1]*alpha_vdim[49]+f[0]*alpha_vdim[48]); 
  out[4] += 0.4330127018922193*(f[11]*alpha_vdim[107]+f[1]*alpha_vdim[97]+f[0]*alpha_vdim[96]); 
  out[5] += 0.3872983346207416*(alpha_cdim[2]*f[12]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.3872983346207416*(f[1]*alpha_vdim[59]+f[11]*alpha_vdim[49])+0.4330127018922193*(f[0]*alpha_vdim[49]+f[1]*alpha_vdim[48]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(f[19]*alpha_vdim[59]+f[5]*alpha_vdim[49]+f[2]*alpha_vdim[48]+alpha_vdim[11]*f[21]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[8] += 0.3872983346207416*(f[1]*alpha_vdim[107]+f[11]*alpha_vdim[97])+0.4330127018922193*(f[0]*alpha_vdim[97]+f[1]*alpha_vdim[96]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[19]*alpha_vdim[107]+f[5]*alpha_vdim[97]+f[2]*alpha_vdim[96]+alpha_vdim[11]*f[25]+alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[21]*alpha_vdim[107]+f[6]*alpha_vdim[97]+f[3]*alpha_vdim[96]+f[25]*alpha_vdim[59]+f[8]*alpha_vdim[49]+f[4]*alpha_vdim[48]); 
  out[11] += 0.9682458365518543*(alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_vdim[11]*f[19]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[13] += 0.9682458365518543*(f[21]*alpha_vdim[59]+f[6]*alpha_vdim[49]+f[3]*alpha_vdim[48]); 
  out[14] += 0.9682458365518543*(f[25]*alpha_vdim[107]+f[8]*alpha_vdim[97]+f[4]*alpha_vdim[96]); 
  out[15] += 0.3872983346207416*(f[5]*alpha_vdim[59]+f[19]*alpha_vdim[49])+0.4330127018922193*(f[2]*alpha_vdim[49]+f[5]*alpha_vdim[48])+0.3872983346207416*(alpha_cdim[2]*f[22]+alpha_vdim[1]*f[21]+f[6]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[16] += 0.3872983346207416*(f[5]*alpha_vdim[107]+f[19]*alpha_vdim[97])+0.4330127018922193*(f[2]*alpha_vdim[97]+f[5]*alpha_vdim[96])+0.3872983346207416*(alpha_cdim[2]*f[26]+alpha_vdim[1]*f[25]+f[8]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]); 
  out[17] += 0.3872983346207416*(f[6]*alpha_vdim[107]+f[21]*alpha_vdim[97])+0.4330127018922193*(f[3]*alpha_vdim[97]+f[6]*alpha_vdim[96])+0.3872983346207416*(f[8]*alpha_vdim[59]+f[25]*alpha_vdim[49])+0.4330127018922193*(f[4]*alpha_vdim[49]+f[8]*alpha_vdim[48]+alpha_cdim[2]*f[18]+alpha_cdim[0]*f[10]); 
  out[18] += 0.4330127018922193*(f[32]*alpha_vdim[107]+f[15]*alpha_vdim[97]+f[7]*alpha_vdim[96]+f[35]*alpha_vdim[59]+f[16]*alpha_vdim[49]+f[9]*alpha_vdim[48]+alpha_vdim[11]*f[37]+alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10]); 
  out[19] += 0.8660254037844386*alpha_cdim[2]*f[20]+0.276641667586244*alpha_vdim[11]*f[11]+0.4330127018922193*(alpha_vdim[0]*f[11]+f[0]*alpha_vdim[11])+0.9682458365518543*alpha_cdim[0]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1]); 
  out[20] += 0.8660254037844386*alpha_vdim[1]*f[19]+0.4330127018922193*alpha_cdim[0]*f[12]+f[5]*(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[2]; 
  out[21] += (0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_vdim[59]+0.3872983346207416*f[1]*alpha_vdim[49]+0.4330127018922193*f[11]*alpha_vdim[48]+0.9682458365518543*(alpha_cdim[2]*f[15]+alpha_cdim[0]*f[6]); 
  out[22] += 0.4330127018922193*(f[20]*alpha_vdim[49]+f[12]*alpha_vdim[48])+0.9682458365518543*(alpha_vdim[11]*f[32]+alpha_vdim[1]*f[15]+alpha_vdim[0]*f[7]); 
  out[23] += 0.8660254037844386*(f[6]*alpha_vdim[59]+f[21]*alpha_vdim[49])+0.9682458365518543*(f[3]*alpha_vdim[49]+f[6]*alpha_vdim[48])+0.4330127018922193*(alpha_cdim[2]*f[24]+alpha_cdim[0]*f[13]); 
  out[24] += 0.9682458365518543*(f[32]*alpha_vdim[59]+f[15]*alpha_vdim[49]+f[7]*alpha_vdim[48])+0.4330127018922193*(alpha_vdim[1]*f[23]+alpha_vdim[0]*f[13]); 
  out[25] += (0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_vdim[107]+0.3872983346207416*f[1]*alpha_vdim[97]+0.4330127018922193*f[11]*alpha_vdim[96]+0.9682458365518543*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[8]); 
  out[26] += 0.4330127018922193*(f[20]*alpha_vdim[97]+f[12]*alpha_vdim[96])+0.9682458365518543*(alpha_vdim[11]*f[35]+alpha_vdim[1]*f[16]+alpha_vdim[0]*f[9]); 
  out[27] += 0.4330127018922193*(f[23]*alpha_vdim[97]+f[13]*alpha_vdim[96])+0.9682458365518543*(f[37]*alpha_vdim[59]+f[17]*alpha_vdim[49]+f[10]*alpha_vdim[48]); 
  out[28] += 0.8660254037844386*(f[8]*alpha_vdim[107]+f[25]*alpha_vdim[97])+0.9682458365518543*(f[4]*alpha_vdim[97]+f[8]*alpha_vdim[96])+0.4330127018922193*(alpha_cdim[2]*f[29]+alpha_cdim[0]*f[14]); 
  out[29] += 0.9682458365518543*(f[35]*alpha_vdim[107]+f[16]*alpha_vdim[97]+f[9]*alpha_vdim[96])+0.4330127018922193*(alpha_vdim[1]*f[28]+alpha_vdim[0]*f[14]); 
  out[30] += 0.9682458365518543*(f[37]*alpha_vdim[107]+f[17]*alpha_vdim[97]+f[10]*alpha_vdim[96])+0.4330127018922193*(f[28]*alpha_vdim[49]+f[14]*alpha_vdim[48]); 
  out[31] += 0.3872983346207416*(f[15]*alpha_vdim[107]+f[32]*alpha_vdim[97])+0.4330127018922193*(f[7]*alpha_vdim[97]+f[15]*alpha_vdim[96])+0.3872983346207416*(f[16]*alpha_vdim[59]+f[35]*alpha_vdim[49])+0.4330127018922193*(f[9]*alpha_vdim[49]+f[16]*alpha_vdim[48])+0.3872983346207416*(alpha_cdim[2]*f[38]+alpha_vdim[1]*f[37])+0.4330127018922193*alpha_cdim[0]*f[18]+0.3872983346207416*alpha_vdim[11]*f[17]+0.4330127018922193*(alpha_vdim[0]*f[17]+(alpha_cdim[2]+alpha_vdim[1])*f[10]); 
  out[32] += (0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_vdim[59]+0.3872983346207416*f[5]*alpha_vdim[49]+0.4330127018922193*f[19]*alpha_vdim[48]+0.8660254037844386*alpha_cdim[2]*f[33]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[21]+0.9682458365518543*alpha_cdim[0]*f[15]+0.4330127018922193*f[3]*alpha_vdim[11]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[6]; 
  out[33] += 0.3872983346207416*f[20]*alpha_vdim[59]+0.4330127018922193*(f[12]*alpha_vdim[49]+f[20]*alpha_vdim[48])+0.8660254037844386*alpha_vdim[1]*f[32]+0.4330127018922193*alpha_cdim[0]*f[22]+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[15]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[7]; 
  out[34] += 0.8660254037844386*(f[15]*alpha_vdim[59]+f[32]*alpha_vdim[49])+0.9682458365518543*(f[7]*alpha_vdim[49]+f[15]*alpha_vdim[48])+0.4330127018922193*alpha_cdim[0]*f[24]+0.3872983346207416*alpha_vdim[11]*f[23]+0.4330127018922193*(alpha_vdim[0]*f[23]+(alpha_cdim[2]+alpha_vdim[1])*f[13]); 
  out[35] += (0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_vdim[107]+0.3872983346207416*f[5]*alpha_vdim[97]+0.4330127018922193*f[19]*alpha_vdim[96]+0.8660254037844386*alpha_cdim[2]*f[36]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[25]+0.9682458365518543*alpha_cdim[0]*f[16]+0.4330127018922193*f[4]*alpha_vdim[11]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[8]; 
  out[36] += 0.3872983346207416*f[20]*alpha_vdim[107]+0.4330127018922193*(f[12]*alpha_vdim[97]+f[20]*alpha_vdim[96])+0.8660254037844386*alpha_vdim[1]*f[35]+0.4330127018922193*alpha_cdim[0]*f[26]+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[16]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[9]; 
  out[37] += (0.276641667586244*f[21]+0.4330127018922193*f[3])*alpha_vdim[107]+0.3872983346207416*f[6]*alpha_vdim[97]+0.4330127018922193*f[21]*alpha_vdim[96]+(0.276641667586244*f[25]+0.4330127018922193*f[4])*alpha_vdim[59]+0.3872983346207416*f[8]*alpha_vdim[49]+0.4330127018922193*f[25]*alpha_vdim[48]+0.9682458365518543*(alpha_cdim[2]*f[31]+alpha_cdim[0]*f[17]); 
  out[38] += 0.4330127018922193*(f[33]*alpha_vdim[97]+f[22]*alpha_vdim[96]+f[36]*alpha_vdim[49]+f[26]*alpha_vdim[48])+0.9682458365518543*(alpha_vdim[11]*f[44]+alpha_vdim[1]*f[31]+alpha_vdim[0]*f[18]); 
  out[39] += 0.3872983346207416*f[23]*alpha_vdim[107]+0.4330127018922193*(f[13]*alpha_vdim[97]+f[23]*alpha_vdim[96])+0.8660254037844386*(f[17]*alpha_vdim[59]+f[37]*alpha_vdim[49])+0.9682458365518543*(f[10]*alpha_vdim[49]+f[17]*alpha_vdim[48])+0.4330127018922193*(alpha_cdim[2]*f[40]+alpha_cdim[0]*f[27]); 
  out[40] += 0.4330127018922193*(f[34]*alpha_vdim[97]+f[24]*alpha_vdim[96])+0.9682458365518543*(f[44]*alpha_vdim[59]+f[31]*alpha_vdim[49]+f[18]*alpha_vdim[48])+0.4330127018922193*(alpha_vdim[1]*f[39]+alpha_vdim[0]*f[27]); 
  out[41] += 0.8660254037844386*(f[16]*alpha_vdim[107]+f[35]*alpha_vdim[97])+0.9682458365518543*(f[9]*alpha_vdim[97]+f[16]*alpha_vdim[96])+0.4330127018922193*alpha_cdim[0]*f[29]+0.3872983346207416*alpha_vdim[11]*f[28]+0.4330127018922193*(alpha_vdim[0]*f[28]+(alpha_cdim[2]+alpha_vdim[1])*f[14]); 
  out[42] += 0.8660254037844386*(f[17]*alpha_vdim[107]+f[37]*alpha_vdim[97])+0.9682458365518543*(f[10]*alpha_vdim[97]+f[17]*alpha_vdim[96])+0.3872983346207416*f[28]*alpha_vdim[59]+0.4330127018922193*(f[14]*alpha_vdim[49]+f[28]*alpha_vdim[48]+alpha_cdim[2]*f[43]+alpha_cdim[0]*f[30]); 
  out[43] += 0.9682458365518543*(f[44]*alpha_vdim[107]+f[31]*alpha_vdim[97]+f[18]*alpha_vdim[96])+0.4330127018922193*(f[41]*alpha_vdim[49]+f[29]*alpha_vdim[48]+alpha_vdim[1]*f[42]+alpha_vdim[0]*f[30]); 
  out[44] += (0.276641667586244*f[32]+0.4330127018922193*f[7])*alpha_vdim[107]+0.3872983346207416*f[15]*alpha_vdim[97]+0.4330127018922193*f[32]*alpha_vdim[96]+(0.276641667586244*f[35]+0.4330127018922193*f[9])*alpha_vdim[59]+0.3872983346207416*f[16]*alpha_vdim[49]+0.4330127018922193*f[35]*alpha_vdim[48]+0.8660254037844386*alpha_cdim[2]*f[45]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[37]+0.9682458365518543*alpha_cdim[0]*f[31]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[17]+0.4330127018922193*f[10]*alpha_vdim[11]; 
  out[45] += 0.3872983346207416*f[33]*alpha_vdim[107]+0.4330127018922193*(f[22]*alpha_vdim[97]+f[33]*alpha_vdim[96])+0.3872983346207416*f[36]*alpha_vdim[59]+0.4330127018922193*(f[26]*alpha_vdim[49]+f[36]*alpha_vdim[48])+0.8660254037844386*alpha_vdim[1]*f[44]+0.4330127018922193*alpha_cdim[0]*f[38]+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[31]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[18]; 
  out[46] += 0.3872983346207416*f[34]*alpha_vdim[107]+0.4330127018922193*(f[24]*alpha_vdim[97]+f[34]*alpha_vdim[96])+0.8660254037844386*(f[31]*alpha_vdim[59]+f[44]*alpha_vdim[49])+0.9682458365518543*(f[18]*alpha_vdim[49]+f[31]*alpha_vdim[48])+0.4330127018922193*alpha_cdim[0]*f[40]+0.3872983346207416*alpha_vdim[11]*f[39]+0.4330127018922193*(alpha_vdim[0]*f[39]+(alpha_cdim[2]+alpha_vdim[1])*f[27]); 
  out[47] += 0.8660254037844386*(f[31]*alpha_vdim[107]+f[44]*alpha_vdim[97])+0.9682458365518543*(f[18]*alpha_vdim[97]+f[31]*alpha_vdim[96])+0.3872983346207416*f[41]*alpha_vdim[59]+0.4330127018922193*(f[29]*alpha_vdim[49]+f[41]*alpha_vdim[48]+alpha_cdim[0]*f[43])+0.3872983346207416*alpha_vdim[11]*f[42]+0.4330127018922193*(alpha_vdim[0]*f[42]+(alpha_cdim[2]+alpha_vdim[1])*f[30]); 

  return alpha_mid; 
} 
