#include <VlasovModDecl.h> 
__host__ __device__ double VlasovNeutralVol1x3vSerP3(const double *w, const double *dxv, const double *boA, const double *f, double *out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // boA:       Input body acceleration.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *Fo0 = &boA[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv11 = 2/dxv[2]; 
  const double *Fo1 = &boA[4]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv12 = 2/dxv[3]; 
  const double *Fo2 = &boA[8]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  double alpha_mid = 0.0; 
  double alpha_cdim[80]; 
  double alpha_vdim[240]; 

  alpha_cdim[0] = 8.0*w0dx0; 
  alpha_cdim[2] = 2.309401076758503*dv0dx0; 
  alpha_mid += std::abs(w0dx0)+0.5*dv0dx0; 

  alpha_vdim[0] = 2.828427124746191*Fo0[0]*dv10; 
  alpha_vdim[1] = 2.828427124746191*Fo0[1]*dv10; 
  alpha_vdim[11] = 2.828427124746191*Fo0[2]*dv10; 
  alpha_vdim[31] = 2.828427124746191*Fo0[3]*dv10; 
  alpha_mid += std::abs(0.125*alpha_vdim[0]-0.1397542485937369*alpha_vdim[11]); 

  alpha_vdim[80] = 2.828427124746191*Fo1[0]*dv11; 
  alpha_vdim[81] = 2.828427124746191*Fo1[1]*dv11; 
  alpha_vdim[91] = 2.828427124746191*Fo1[2]*dv11; 
  alpha_vdim[111] = 2.828427124746191*Fo1[3]*dv11; 
  alpha_mid += std::abs(0.125*alpha_vdim[80]-0.1397542485937369*alpha_vdim[91]); 

  alpha_vdim[160] = 2.828427124746191*Fo2[0]*dv12; 
  alpha_vdim[161] = 2.828427124746191*Fo2[1]*dv12; 
  alpha_vdim[171] = 2.828427124746191*Fo2[2]*dv12; 
  alpha_vdim[191] = 2.828427124746191*Fo2[3]*dv12; 
  alpha_mid += std::abs(0.125*alpha_vdim[160]-0.1397542485937369*alpha_vdim[171]); 

  out[1] += 0.4330127018922193*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.4330127018922193*(alpha_vdim[31]*f[31]+alpha_vdim[11]*f[11]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[31]*alpha_vdim[111]+f[11]*alpha_vdim[91]+f[1]*alpha_vdim[81]+f[0]*alpha_vdim[80]); 
  out[4] += 0.4330127018922193*(f[31]*alpha_vdim[191]+f[11]*alpha_vdim[171]+f[1]*alpha_vdim[161]+f[0]*alpha_vdim[160]); 
  out[5] += 0.3803194146278324*(alpha_vdim[11]*f[31]+f[11]*alpha_vdim[31])+0.3872983346207416*(alpha_cdim[2]*f[12]+alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[6] += 0.3803194146278324*(f[11]*alpha_vdim[111]+f[31]*alpha_vdim[91])+0.3872983346207416*(f[1]*alpha_vdim[91]+f[11]*alpha_vdim[81])+0.4330127018922193*(f[0]*alpha_vdim[81]+f[1]*alpha_vdim[80]+alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[7] += 0.4330127018922193*(f[48]*alpha_vdim[111]+f[19]*alpha_vdim[91]+f[5]*alpha_vdim[81]+f[2]*alpha_vdim[80]+alpha_vdim[31]*f[50]+alpha_vdim[11]*f[21]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[8] += 0.3803194146278324*(f[11]*alpha_vdim[191]+f[31]*alpha_vdim[171])+0.3872983346207416*(f[1]*alpha_vdim[171]+f[11]*alpha_vdim[161])+0.4330127018922193*(f[0]*alpha_vdim[161]+f[1]*alpha_vdim[160]+alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[9] += 0.4330127018922193*(f[48]*alpha_vdim[191]+f[19]*alpha_vdim[171]+f[5]*alpha_vdim[161]+f[2]*alpha_vdim[160]+alpha_vdim[31]*f[54]+alpha_vdim[11]*f[25]+alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[50]*alpha_vdim[191]+f[21]*alpha_vdim[171]+f[6]*alpha_vdim[161]+f[3]*alpha_vdim[160]+f[54]*alpha_vdim[111]+f[25]*alpha_vdim[91]+f[8]*alpha_vdim[81]+f[4]*alpha_vdim[80]); 
  out[11] += 0.9682458365518543*(alpha_cdim[2]*f[5]+alpha_cdim[0]*f[1]); 
  out[12] += 0.9682458365518543*(alpha_vdim[31]*f[48]+alpha_vdim[11]*f[19]+alpha_vdim[1]*f[5]+alpha_vdim[0]*f[2]); 
  out[13] += 0.9682458365518543*(f[50]*alpha_vdim[111]+f[21]*alpha_vdim[91]+f[6]*alpha_vdim[81]+f[3]*alpha_vdim[80]); 
  out[14] += 0.9682458365518543*(f[54]*alpha_vdim[191]+f[25]*alpha_vdim[171]+f[8]*alpha_vdim[161]+f[4]*alpha_vdim[160]); 
  out[15] += 0.3803194146278324*(f[19]*alpha_vdim[111]+f[48]*alpha_vdim[91])+0.3872983346207416*(f[5]*alpha_vdim[91]+f[19]*alpha_vdim[81])+0.4330127018922193*(f[2]*alpha_vdim[81]+f[5]*alpha_vdim[80])+0.3803194146278324*(alpha_vdim[11]*f[50]+f[21]*alpha_vdim[31])+0.3872983346207416*(alpha_cdim[2]*f[22]+alpha_vdim[1]*f[21]+f[6]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*f[7]+alpha_vdim[0]*f[6]+(alpha_cdim[2]+alpha_vdim[1])*f[3]); 
  out[16] += 0.3803194146278324*(f[19]*alpha_vdim[191]+f[48]*alpha_vdim[171])+0.3872983346207416*(f[5]*alpha_vdim[171]+f[19]*alpha_vdim[161])+0.4330127018922193*(f[2]*alpha_vdim[161]+f[5]*alpha_vdim[160])+0.3803194146278324*(alpha_vdim[11]*f[54]+f[25]*alpha_vdim[31])+0.3872983346207416*(alpha_cdim[2]*f[26]+alpha_vdim[1]*f[25]+f[8]*alpha_vdim[11])+0.4330127018922193*(alpha_cdim[0]*f[9]+alpha_vdim[0]*f[8]+(alpha_cdim[2]+alpha_vdim[1])*f[4]); 
  out[17] += 0.3803194146278324*(f[21]*alpha_vdim[191]+f[50]*alpha_vdim[171])+0.3872983346207416*(f[6]*alpha_vdim[171]+f[21]*alpha_vdim[161])+0.4330127018922193*(f[3]*alpha_vdim[161]+f[6]*alpha_vdim[160])+0.3803194146278324*(f[25]*alpha_vdim[111]+f[54]*alpha_vdim[91])+0.3872983346207416*(f[8]*alpha_vdim[91]+f[25]*alpha_vdim[81])+0.4330127018922193*(f[4]*alpha_vdim[81]+f[8]*alpha_vdim[80]+alpha_cdim[2]*f[18]+alpha_cdim[0]*f[10]); 
  out[18] += 0.4330127018922193*(f[64]*alpha_vdim[191]+f[36]*alpha_vdim[171]+f[15]*alpha_vdim[161]+f[7]*alpha_vdim[160]+f[67]*alpha_vdim[111]+f[39]*alpha_vdim[91]+f[16]*alpha_vdim[81]+f[9]*alpha_vdim[80]+alpha_vdim[31]*f[69]+alpha_vdim[11]*f[41]+alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10]); 
  out[19] += 0.2581988897471612*alpha_vdim[31]*f[31]+0.3803194146278324*(alpha_vdim[1]*f[31]+f[1]*alpha_vdim[31])+0.8660254037844386*alpha_cdim[2]*f[20]+0.276641667586244*alpha_vdim[11]*f[11]+0.4330127018922193*(alpha_vdim[0]*f[11]+f[0]*alpha_vdim[11])+0.9682458365518543*alpha_cdim[0]*f[5]+f[1]*(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1]); 
  out[20] += 0.8504200642707612*alpha_vdim[11]*f[48]+0.3803194146278324*alpha_cdim[2]*f[32]+f[19]*(0.8504200642707612*alpha_vdim[31]+0.8660254037844386*alpha_vdim[1])+0.4330127018922193*alpha_cdim[0]*f[12]+f[5]*(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[2]; 
  out[21] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_vdim[111]+(0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_vdim[91]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[81]+0.4330127018922193*f[11]*alpha_vdim[80]+0.9682458365518543*(alpha_cdim[2]*f[15]+alpha_cdim[0]*f[6]); 
  out[22] += 0.4330127018922193*(f[20]*alpha_vdim[81]+f[12]*alpha_vdim[80])+0.9682458365518543*(alpha_vdim[31]*f[64]+alpha_vdim[11]*f[36]+alpha_vdim[1]*f[15]+alpha_vdim[0]*f[7]); 
  out[23] += 0.8504200642707612*(f[21]*alpha_vdim[111]+f[50]*alpha_vdim[91])+0.8660254037844386*(f[6]*alpha_vdim[91]+f[21]*alpha_vdim[81])+0.9682458365518543*(f[3]*alpha_vdim[81]+f[6]*alpha_vdim[80])+0.4330127018922193*(alpha_cdim[2]*f[24]+alpha_cdim[0]*f[13]); 
  out[24] += 0.9682458365518543*(f[64]*alpha_vdim[111]+f[36]*alpha_vdim[91]+f[15]*alpha_vdim[81]+f[7]*alpha_vdim[80])+0.4330127018922193*(alpha_vdim[1]*f[23]+alpha_vdim[0]*f[13]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alpha_vdim[191]+(0.276641667586244*f[11]+0.4330127018922193*f[0])*alpha_vdim[171]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alpha_vdim[161]+0.4330127018922193*f[11]*alpha_vdim[160]+0.9682458365518543*(alpha_cdim[2]*f[16]+alpha_cdim[0]*f[8]); 
  out[26] += 0.4330127018922193*(f[20]*alpha_vdim[161]+f[12]*alpha_vdim[160])+0.9682458365518543*(alpha_vdim[31]*f[67]+alpha_vdim[11]*f[39]+alpha_vdim[1]*f[16]+alpha_vdim[0]*f[9]); 
  out[27] += 0.4330127018922193*(f[23]*alpha_vdim[161]+f[13]*alpha_vdim[160])+0.9682458365518543*(f[69]*alpha_vdim[111]+f[41]*alpha_vdim[91]+f[17]*alpha_vdim[81]+f[10]*alpha_vdim[80]); 
  out[28] += 0.8504200642707612*(f[25]*alpha_vdim[191]+f[54]*alpha_vdim[171])+0.8660254037844386*(f[8]*alpha_vdim[171]+f[25]*alpha_vdim[161])+0.9682458365518543*(f[4]*alpha_vdim[161]+f[8]*alpha_vdim[160])+0.4330127018922193*(alpha_cdim[2]*f[29]+alpha_cdim[0]*f[14]); 
  out[29] += 0.9682458365518543*(f[67]*alpha_vdim[191]+f[39]*alpha_vdim[171]+f[16]*alpha_vdim[161]+f[9]*alpha_vdim[160])+0.4330127018922193*(alpha_vdim[1]*f[28]+alpha_vdim[0]*f[14]); 
  out[30] += 0.9682458365518543*(f[69]*alpha_vdim[191]+f[41]*alpha_vdim[171]+f[17]*alpha_vdim[161]+f[10]*alpha_vdim[160])+0.4330127018922193*(f[28]*alpha_vdim[81]+f[14]*alpha_vdim[80]); 
  out[31] += 1.479019945774904*(alpha_cdim[2]*f[19]+alpha_cdim[0]*f[11])+0.6614378277661477*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[32] += 0.6614378277661477*alpha_vdim[31]*f[31]+1.479019945774904*(alpha_vdim[1]*f[20]+alpha_vdim[0]*f[12])+0.6614378277661477*(alpha_vdim[11]*f[11]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[33] += 0.6614378277661477*(f[31]*alpha_vdim[111]+f[11]*alpha_vdim[91])+(1.479019945774904*f[23]+0.6614378277661477*f[1])*alpha_vdim[81]+(1.479019945774904*f[13]+0.6614378277661477*f[0])*alpha_vdim[80]; 
  out[34] += 0.6614378277661477*(f[31]*alpha_vdim[191]+f[11]*alpha_vdim[171])+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_vdim[161]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alpha_vdim[160]; 
  out[35] += 0.3803194146278324*(f[36]*alpha_vdim[191]+f[64]*alpha_vdim[171])+0.3872983346207416*(f[15]*alpha_vdim[171]+f[36]*alpha_vdim[161])+0.4330127018922193*(f[7]*alpha_vdim[161]+f[15]*alpha_vdim[160])+0.3803194146278324*(f[39]*alpha_vdim[111]+f[67]*alpha_vdim[91])+0.3872983346207416*(f[16]*alpha_vdim[91]+f[39]*alpha_vdim[81])+0.4330127018922193*(f[9]*alpha_vdim[81]+f[16]*alpha_vdim[80])+0.3803194146278324*alpha_vdim[11]*f[69]+0.3872983346207416*alpha_cdim[2]*f[42]+(0.3803194146278324*alpha_vdim[31]+0.3872983346207416*alpha_vdim[1])*f[41]+0.4330127018922193*alpha_cdim[0]*f[18]+0.3872983346207416*alpha_vdim[11]*f[17]+0.4330127018922193*(alpha_vdim[0]*f[17]+(alpha_cdim[2]+alpha_vdim[1])*f[10]); 
  out[36] += (0.2581988897471612*f[48]+0.3803194146278324*f[5])*alpha_vdim[111]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_vdim[91]+(0.3803194146278324*f[48]+0.3872983346207416*f[5])*alpha_vdim[81]+0.4330127018922193*f[19]*alpha_vdim[80]+(0.2581988897471612*alpha_vdim[31]+0.3803194146278324*alpha_vdim[1])*f[50]+0.8660254037844386*alpha_cdim[2]*f[37]+0.3803194146278324*f[6]*alpha_vdim[31]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[21]+0.9682458365518543*alpha_cdim[0]*f[15]+0.4330127018922193*f[3]*alpha_vdim[11]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[6]; 
  out[37] += 0.3872983346207416*f[20]*alpha_vdim[91]+0.4330127018922193*(f[12]*alpha_vdim[81]+f[20]*alpha_vdim[80])+0.8504200642707612*alpha_vdim[11]*f[64]+0.3803194146278324*alpha_cdim[2]*f[51]+(0.8504200642707612*alpha_vdim[31]+0.8660254037844386*alpha_vdim[1])*f[36]+0.4330127018922193*alpha_cdim[0]*f[22]+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[15]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[7]; 
  out[38] += 0.8504200642707612*(f[36]*alpha_vdim[111]+f[64]*alpha_vdim[91])+0.8660254037844386*(f[15]*alpha_vdim[91]+f[36]*alpha_vdim[81])+0.9682458365518543*(f[7]*alpha_vdim[81]+f[15]*alpha_vdim[80])+0.4330127018922193*alpha_cdim[0]*f[24]+0.3872983346207416*alpha_vdim[11]*f[23]+0.4330127018922193*(alpha_vdim[0]*f[23]+(alpha_cdim[2]+alpha_vdim[1])*f[13]); 
  out[39] += (0.2581988897471612*f[48]+0.3803194146278324*f[5])*alpha_vdim[191]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alpha_vdim[171]+(0.3803194146278324*f[48]+0.3872983346207416*f[5])*alpha_vdim[161]+0.4330127018922193*f[19]*alpha_vdim[160]+(0.2581988897471612*alpha_vdim[31]+0.3803194146278324*alpha_vdim[1])*f[54]+0.8660254037844386*alpha_cdim[2]*f[40]+0.3803194146278324*f[8]*alpha_vdim[31]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[25]+0.9682458365518543*alpha_cdim[0]*f[16]+0.4330127018922193*f[4]*alpha_vdim[11]+(0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])*f[8]; 
  out[40] += 0.3872983346207416*f[20]*alpha_vdim[171]+0.4330127018922193*(f[12]*alpha_vdim[161]+f[20]*alpha_vdim[160])+0.8504200642707612*alpha_vdim[11]*f[67]+0.3803194146278324*alpha_cdim[2]*f[55]+(0.8504200642707612*alpha_vdim[31]+0.8660254037844386*alpha_vdim[1])*f[39]+0.4330127018922193*alpha_cdim[0]*f[26]+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[16]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[9]; 
  out[41] += (0.2581988897471612*f[50]+0.3803194146278324*f[6])*alpha_vdim[191]+(0.276641667586244*f[21]+0.4330127018922193*f[3])*alpha_vdim[171]+(0.3803194146278324*f[50]+0.3872983346207416*f[6])*alpha_vdim[161]+0.4330127018922193*f[21]*alpha_vdim[160]+(0.2581988897471612*f[54]+0.3803194146278324*f[8])*alpha_vdim[111]+(0.276641667586244*f[25]+0.4330127018922193*f[4])*alpha_vdim[91]+(0.3803194146278324*f[54]+0.3872983346207416*f[8])*alpha_vdim[81]+0.4330127018922193*f[25]*alpha_vdim[80]+0.9682458365518543*(alpha_cdim[2]*f[35]+alpha_cdim[0]*f[17]); 
  out[42] += 0.4330127018922193*(f[37]*alpha_vdim[161]+f[22]*alpha_vdim[160]+f[40]*alpha_vdim[81]+f[26]*alpha_vdim[80])+0.9682458365518543*(alpha_vdim[31]*f[76]+alpha_vdim[11]*f[60]+alpha_vdim[1]*f[35]+alpha_vdim[0]*f[18]); 
  out[43] += 0.3872983346207416*f[23]*alpha_vdim[171]+0.4330127018922193*(f[13]*alpha_vdim[161]+f[23]*alpha_vdim[160])+0.8504200642707612*(f[41]*alpha_vdim[111]+f[69]*alpha_vdim[91])+0.8660254037844386*(f[17]*alpha_vdim[91]+f[41]*alpha_vdim[81])+0.9682458365518543*(f[10]*alpha_vdim[81]+f[17]*alpha_vdim[80])+0.4330127018922193*(alpha_cdim[2]*f[44]+alpha_cdim[0]*f[27]); 
  out[44] += 0.4330127018922193*(f[38]*alpha_vdim[161]+f[24]*alpha_vdim[160])+0.9682458365518543*(f[76]*alpha_vdim[111]+f[60]*alpha_vdim[91]+f[35]*alpha_vdim[81]+f[18]*alpha_vdim[80])+0.4330127018922193*(alpha_vdim[1]*f[43]+alpha_vdim[0]*f[27]); 
  out[45] += 0.8504200642707612*(f[39]*alpha_vdim[191]+f[67]*alpha_vdim[171])+0.8660254037844386*(f[16]*alpha_vdim[171]+f[39]*alpha_vdim[161])+0.9682458365518543*(f[9]*alpha_vdim[161]+f[16]*alpha_vdim[160])+0.4330127018922193*alpha_cdim[0]*f[29]+0.3872983346207416*alpha_vdim[11]*f[28]+0.4330127018922193*(alpha_vdim[0]*f[28]+(alpha_cdim[2]+alpha_vdim[1])*f[14]); 
  out[46] += 0.8504200642707612*(f[41]*alpha_vdim[191]+f[69]*alpha_vdim[171])+0.8660254037844386*(f[17]*alpha_vdim[171]+f[41]*alpha_vdim[161])+0.9682458365518543*(f[10]*alpha_vdim[161]+f[17]*alpha_vdim[160])+0.3872983346207416*f[28]*alpha_vdim[91]+0.4330127018922193*(f[14]*alpha_vdim[81]+f[28]*alpha_vdim[80]+alpha_cdim[2]*f[47]+alpha_cdim[0]*f[30]); 
  out[47] += 0.9682458365518543*(f[76]*alpha_vdim[191]+f[60]*alpha_vdim[171]+f[35]*alpha_vdim[161]+f[18]*alpha_vdim[160])+0.4330127018922193*(f[45]*alpha_vdim[81]+f[29]*alpha_vdim[80]+alpha_vdim[1]*f[46]+alpha_vdim[0]*f[30]); 
  out[48] += (0.2581988897471612*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[31]+(0.2581988897471612*f[11]+0.4330127018922193*f[0])*alpha_vdim[31]+1.479019945774904*alpha_cdim[0]*f[19]+alpha_cdim[2]*(0.5916079783099616*f[12]+1.479019945774904*f[11])+0.3803194146278324*(alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.6614378277661477*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]); 
  out[49] += 0.4330127018922193*alpha_cdim[0]*f[32]+0.5809475019311124*(alpha_vdim[11]*f[31]+f[11]*alpha_vdim[31])+(1.322875655532295*alpha_vdim[11]+1.479019945774904*alpha_vdim[0])*f[20]+(0.3803194146278324*alpha_cdim[2]+1.479019945774904*alpha_vdim[1])*f[12]+0.5916079783099616*(alpha_vdim[1]*f[11]+f[1]*alpha_vdim[11])+0.6614378277661477*(alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[50] += (0.2581988897471612*f[11]+0.4330127018922193*f[0])*alpha_vdim[111]+0.2581988897471612*f[31]*alpha_vdim[91]+0.3803194146278324*(f[1]*alpha_vdim[91]+f[11]*alpha_vdim[81])+0.4330127018922193*f[31]*alpha_vdim[80]+1.479019945774904*(alpha_cdim[2]*f[36]+alpha_cdim[0]*f[21])+0.6614378277661477*(alpha_cdim[2]*f[7]+alpha_cdim[0]*f[3]); 
  out[51] += 0.4330127018922193*(f[49]*alpha_vdim[81]+f[32]*alpha_vdim[80])+0.6614378277661477*alpha_vdim[31]*f[50]+1.479019945774904*(alpha_vdim[1]*f[37]+alpha_vdim[0]*f[22])+0.6614378277661477*(alpha_vdim[11]*f[21]+alpha_vdim[1]*f[6]+alpha_vdim[0]*f[3]); 
  out[52] += 0.5809475019311124*f[11]*alpha_vdim[111]+(0.5809475019311124*f[31]+1.322875655532295*f[23]+0.5916079783099616*f[1])*alpha_vdim[91]+(1.479019945774904*f[13]+0.5916079783099616*f[11]+0.6614378277661477*f[0])*alpha_vdim[81]+(1.479019945774904*f[23]+0.6614378277661477*f[1])*alpha_vdim[80]+0.4330127018922193*(alpha_cdim[2]*f[53]+alpha_cdim[0]*f[33]); 
  out[53] += 0.6614378277661477*(f[48]*alpha_vdim[111]+f[19]*alpha_vdim[91])+(1.479019945774904*f[38]+0.6614378277661477*f[5])*alpha_vdim[81]+(1.479019945774904*f[24]+0.6614378277661477*f[2])*alpha_vdim[80]+0.4330127018922193*(alpha_vdim[1]*f[52]+alpha_vdim[0]*f[33]); 
  out[54] += (0.2581988897471612*f[11]+0.4330127018922193*f[0])*alpha_vdim[191]+0.2581988897471612*f[31]*alpha_vdim[171]+0.3803194146278324*(f[1]*alpha_vdim[171]+f[11]*alpha_vdim[161])+0.4330127018922193*f[31]*alpha_vdim[160]+1.479019945774904*(alpha_cdim[2]*f[39]+alpha_cdim[0]*f[25])+0.6614378277661477*(alpha_cdim[2]*f[9]+alpha_cdim[0]*f[4]); 
  out[55] += 0.4330127018922193*(f[49]*alpha_vdim[161]+f[32]*alpha_vdim[160])+0.6614378277661477*alpha_vdim[31]*f[54]+1.479019945774904*(alpha_vdim[1]*f[40]+alpha_vdim[0]*f[26])+0.6614378277661477*(alpha_vdim[11]*f[25]+alpha_vdim[1]*f[8]+alpha_vdim[0]*f[4]); 
  out[56] += 0.4330127018922193*(f[52]*alpha_vdim[161]+f[33]*alpha_vdim[160])+0.6614378277661477*(f[54]*alpha_vdim[111]+f[25]*alpha_vdim[91])+(1.479019945774904*f[43]+0.6614378277661477*f[8])*alpha_vdim[81]+(1.479019945774904*f[27]+0.6614378277661477*f[4])*alpha_vdim[80]; 
  out[57] += 0.5809475019311124*f[11]*alpha_vdim[191]+(0.5809475019311124*f[31]+1.322875655532295*f[28]+0.5916079783099616*f[1])*alpha_vdim[171]+(1.479019945774904*f[14]+0.5916079783099616*f[11]+0.6614378277661477*f[0])*alpha_vdim[161]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alpha_vdim[160]+0.4330127018922193*(alpha_cdim[2]*f[58]+alpha_cdim[0]*f[34]); 
  out[58] += 0.6614378277661477*(f[48]*alpha_vdim[191]+f[19]*alpha_vdim[171])+(1.479019945774904*f[45]+0.6614378277661477*f[5])*alpha_vdim[161]+(1.479019945774904*f[29]+0.6614378277661477*f[2])*alpha_vdim[160]+0.4330127018922193*(alpha_vdim[1]*f[57]+alpha_vdim[0]*f[34]); 
  out[59] += 0.6614378277661477*(f[50]*alpha_vdim[191]+f[21]*alpha_vdim[171])+(1.479019945774904*f[46]+0.6614378277661477*f[6])*alpha_vdim[161]+(1.479019945774904*f[30]+0.6614378277661477*f[3])*alpha_vdim[160]+0.4330127018922193*(f[57]*alpha_vdim[81]+f[34]*alpha_vdim[80]); 
  out[60] += (0.2581988897471612*f[64]+0.3803194146278324*f[15])*alpha_vdim[191]+(0.276641667586244*f[36]+0.4330127018922193*f[7])*alpha_vdim[171]+(0.3803194146278324*f[64]+0.3872983346207416*f[15])*alpha_vdim[161]+0.4330127018922193*f[36]*alpha_vdim[160]+(0.2581988897471612*f[67]+0.3803194146278324*f[16])*alpha_vdim[111]+(0.276641667586244*f[39]+0.4330127018922193*f[9])*alpha_vdim[91]+(0.3803194146278324*f[67]+0.3872983346207416*f[16])*alpha_vdim[81]+0.4330127018922193*f[39]*alpha_vdim[80]+(0.2581988897471612*alpha_vdim[31]+0.3803194146278324*alpha_vdim[1])*f[69]+0.8660254037844386*alpha_cdim[2]*f[61]+(0.276641667586244*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[41]+0.9682458365518543*alpha_cdim[0]*f[35]+f[17]*(0.3803194146278324*alpha_vdim[31]+0.9682458365518543*alpha_cdim[2]+0.3872983346207416*alpha_vdim[1])+0.4330127018922193*f[10]*alpha_vdim[11]; 
  out[61] += 0.3872983346207416*f[37]*alpha_vdim[171]+0.4330127018922193*(f[22]*alpha_vdim[161]+f[37]*alpha_vdim[160])+0.3872983346207416*f[40]*alpha_vdim[91]+0.4330127018922193*(f[26]*alpha_vdim[81]+f[40]*alpha_vdim[80])+0.8504200642707612*alpha_vdim[11]*f[76]+0.3803194146278324*alpha_cdim[2]*f[70]+(0.8504200642707612*alpha_vdim[31]+0.8660254037844386*alpha_vdim[1])*f[60]+0.4330127018922193*alpha_cdim[0]*f[42]+(0.8660254037844386*alpha_vdim[11]+0.9682458365518543*alpha_vdim[0])*f[35]+(0.3872983346207416*alpha_cdim[2]+0.9682458365518543*alpha_vdim[1])*f[18]; 
  out[62] += 0.3872983346207416*f[38]*alpha_vdim[171]+0.4330127018922193*(f[24]*alpha_vdim[161]+f[38]*alpha_vdim[160])+0.8504200642707612*(f[60]*alpha_vdim[111]+f[76]*alpha_vdim[91])+0.8660254037844386*(f[35]*alpha_vdim[91]+f[60]*alpha_vdim[81])+0.9682458365518543*(f[18]*alpha_vdim[81]+f[35]*alpha_vdim[80])+0.4330127018922193*alpha_cdim[0]*f[44]+0.3872983346207416*alpha_vdim[11]*f[43]+0.4330127018922193*(alpha_vdim[0]*f[43]+(alpha_cdim[2]+alpha_vdim[1])*f[27]); 
  out[63] += 0.8504200642707612*(f[60]*alpha_vdim[191]+f[76]*alpha_vdim[171])+0.8660254037844386*(f[35]*alpha_vdim[171]+f[60]*alpha_vdim[161])+0.9682458365518543*(f[18]*alpha_vdim[161]+f[35]*alpha_vdim[160])+0.3872983346207416*f[45]*alpha_vdim[91]+0.4330127018922193*(f[29]*alpha_vdim[81]+f[45]*alpha_vdim[80]+alpha_cdim[0]*f[47])+0.3872983346207416*alpha_vdim[11]*f[46]+0.4330127018922193*(alpha_vdim[0]*f[46]+(alpha_cdim[2]+alpha_vdim[1])*f[30]); 
  out[64] += (0.2581988897471612*f[19]+0.4330127018922193*f[2])*alpha_vdim[111]+0.2581988897471612*f[48]*alpha_vdim[91]+0.3803194146278324*(f[5]*alpha_vdim[91]+f[19]*alpha_vdim[81])+0.4330127018922193*f[48]*alpha_vdim[80]+(0.2581988897471612*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[50]+1.479019945774904*alpha_cdim[0]*f[36]+(0.2581988897471612*f[21]+0.4330127018922193*f[3])*alpha_vdim[31]+alpha_cdim[2]*(0.5916079783099616*f[22]+1.479019945774904*f[21])+0.3803194146278324*(alpha_vdim[1]*f[21]+f[6]*alpha_vdim[11])+0.6614378277661477*(alpha_cdim[0]*f[7]+alpha_cdim[2]*f[3]); 
  out[65] += 0.3872983346207416*f[49]*alpha_vdim[91]+0.4330127018922193*(f[32]*alpha_vdim[81]+f[49]*alpha_vdim[80]+alpha_cdim[0]*f[51])+0.5809475019311124*alpha_vdim[11]*f[50]+(1.322875655532295*alpha_vdim[11]+1.479019945774904*alpha_vdim[0])*f[37]+0.5809475019311124*f[21]*alpha_vdim[31]+(0.3803194146278324*alpha_cdim[2]+1.479019945774904*alpha_vdim[1])*f[22]+0.5916079783099616*(alpha_vdim[1]*f[21]+f[6]*alpha_vdim[11])+0.6614378277661477*(alpha_vdim[0]*f[6]+alpha_vdim[1]*f[3]); 
  out[66] += 0.5809475019311124*f[19]*alpha_vdim[111]+(0.5809475019311124*f[48]+1.322875655532295*f[38]+0.5916079783099616*f[5])*alpha_vdim[91]+(1.479019945774904*f[24]+0.5916079783099616*f[19]+0.6614378277661477*f[2])*alpha_vdim[81]+(1.479019945774904*f[38]+0.6614378277661477*f[5])*alpha_vdim[80]+0.4330127018922193*alpha_cdim[0]*f[53]+0.3872983346207416*alpha_vdim[11]*f[52]+0.4330127018922193*(alpha_vdim[0]*f[52]+(alpha_cdim[2]+alpha_vdim[1])*f[33]); 
  out[67] += (0.2581988897471612*f[19]+0.4330127018922193*f[2])*alpha_vdim[191]+0.2581988897471612*f[48]*alpha_vdim[171]+0.3803194146278324*(f[5]*alpha_vdim[171]+f[19]*alpha_vdim[161])+0.4330127018922193*f[48]*alpha_vdim[160]+(0.2581988897471612*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[54]+1.479019945774904*alpha_cdim[0]*f[39]+(0.2581988897471612*f[25]+0.4330127018922193*f[4])*alpha_vdim[31]+alpha_cdim[2]*(0.5916079783099616*f[26]+1.479019945774904*f[25])+0.3803194146278324*(alpha_vdim[1]*f[25]+f[8]*alpha_vdim[11])+0.6614378277661477*(alpha_cdim[0]*f[9]+alpha_cdim[2]*f[4]); 
  out[68] += 0.3872983346207416*f[49]*alpha_vdim[171]+0.4330127018922193*(f[32]*alpha_vdim[161]+f[49]*alpha_vdim[160]+alpha_cdim[0]*f[55])+0.5809475019311124*alpha_vdim[11]*f[54]+(1.322875655532295*alpha_vdim[11]+1.479019945774904*alpha_vdim[0])*f[40]+0.5809475019311124*f[25]*alpha_vdim[31]+(0.3803194146278324*alpha_cdim[2]+1.479019945774904*alpha_vdim[1])*f[26]+0.5916079783099616*(alpha_vdim[1]*f[25]+f[8]*alpha_vdim[11])+0.6614378277661477*(alpha_vdim[0]*f[8]+alpha_vdim[1]*f[4]); 
  out[69] += (0.2581988897471612*f[21]+0.4330127018922193*f[3])*alpha_vdim[191]+0.2581988897471612*f[50]*alpha_vdim[171]+0.3803194146278324*(f[6]*alpha_vdim[171]+f[21]*alpha_vdim[161])+0.4330127018922193*f[50]*alpha_vdim[160]+(0.2581988897471612*f[25]+0.4330127018922193*f[4])*alpha_vdim[111]+0.2581988897471612*f[54]*alpha_vdim[91]+0.3803194146278324*(f[8]*alpha_vdim[91]+f[25]*alpha_vdim[81])+0.4330127018922193*f[54]*alpha_vdim[80]+1.479019945774904*(alpha_cdim[2]*f[60]+alpha_cdim[0]*f[41])+0.6614378277661477*(alpha_cdim[2]*f[18]+alpha_cdim[0]*f[10]); 
  out[70] += 0.4330127018922193*(f[65]*alpha_vdim[161]+f[51]*alpha_vdim[160]+f[68]*alpha_vdim[81]+f[55]*alpha_vdim[80])+0.6614378277661477*alpha_vdim[31]*f[69]+1.479019945774904*(alpha_vdim[1]*f[61]+alpha_vdim[0]*f[42])+0.6614378277661477*(alpha_vdim[11]*f[41]+alpha_vdim[1]*f[17]+alpha_vdim[0]*f[10]); 
  out[71] += 0.3872983346207416*f[52]*alpha_vdim[171]+0.4330127018922193*(f[33]*alpha_vdim[161]+f[52]*alpha_vdim[160])+0.5809475019311124*f[25]*alpha_vdim[111]+(0.5809475019311124*f[54]+1.322875655532295*f[43]+0.5916079783099616*f[8])*alpha_vdim[91]+(1.479019945774904*f[27]+0.5916079783099616*f[25]+0.6614378277661477*f[4])*alpha_vdim[81]+(1.479019945774904*f[43]+0.6614378277661477*f[8])*alpha_vdim[80]+0.4330127018922193*(alpha_cdim[2]*f[72]+alpha_cdim[0]*f[56]); 
  out[72] += 0.4330127018922193*(f[66]*alpha_vdim[161]+f[53]*alpha_vdim[160])+0.6614378277661477*(f[67]*alpha_vdim[111]+f[39]*alpha_vdim[91])+(1.479019945774904*f[62]+0.6614378277661477*f[16])*alpha_vdim[81]+(1.479019945774904*f[44]+0.6614378277661477*f[9])*alpha_vdim[80]+0.4330127018922193*(alpha_vdim[1]*f[71]+alpha_vdim[0]*f[56]); 
  out[73] += 0.5809475019311124*f[19]*alpha_vdim[191]+(0.5809475019311124*f[48]+1.322875655532295*f[45]+0.5916079783099616*f[5])*alpha_vdim[171]+(1.479019945774904*f[29]+0.5916079783099616*f[19]+0.6614378277661477*f[2])*alpha_vdim[161]+(1.479019945774904*f[45]+0.6614378277661477*f[5])*alpha_vdim[160]+0.4330127018922193*alpha_cdim[0]*f[58]+0.3872983346207416*alpha_vdim[11]*f[57]+0.4330127018922193*(alpha_vdim[0]*f[57]+(alpha_cdim[2]+alpha_vdim[1])*f[34]); 
  out[74] += 0.5809475019311124*f[21]*alpha_vdim[191]+(0.5809475019311124*f[50]+1.322875655532295*f[46]+0.5916079783099616*f[6])*alpha_vdim[171]+(1.479019945774904*f[30]+0.5916079783099616*f[21]+0.6614378277661477*f[3])*alpha_vdim[161]+(1.479019945774904*f[46]+0.6614378277661477*f[6])*alpha_vdim[160]+0.3872983346207416*f[57]*alpha_vdim[91]+0.4330127018922193*(f[34]*alpha_vdim[81]+f[57]*alpha_vdim[80]+alpha_cdim[2]*f[75]+alpha_cdim[0]*f[59]); 
  out[75] += 0.6614378277661477*(f[64]*alpha_vdim[191]+f[36]*alpha_vdim[171])+(1.479019945774904*f[63]+0.6614378277661477*f[15])*alpha_vdim[161]+(1.479019945774904*f[47]+0.6614378277661477*f[7])*alpha_vdim[160]+0.4330127018922193*(f[73]*alpha_vdim[81]+f[58]*alpha_vdim[80]+alpha_vdim[1]*f[74]+alpha_vdim[0]*f[59]); 
  out[76] += (0.2581988897471612*f[36]+0.4330127018922193*f[7])*alpha_vdim[191]+0.2581988897471612*f[64]*alpha_vdim[171]+0.3803194146278324*(f[15]*alpha_vdim[171]+f[36]*alpha_vdim[161])+0.4330127018922193*f[64]*alpha_vdim[160]+(0.2581988897471612*f[39]+0.4330127018922193*f[9])*alpha_vdim[111]+0.2581988897471612*f[67]*alpha_vdim[91]+0.3803194146278324*(f[16]*alpha_vdim[91]+f[39]*alpha_vdim[81])+0.4330127018922193*f[67]*alpha_vdim[80]+(0.2581988897471612*alpha_vdim[11]+0.4330127018922193*alpha_vdim[0])*f[69]+1.479019945774904*alpha_cdim[0]*f[60]+0.5916079783099616*alpha_cdim[2]*f[42]+(0.2581988897471612*alpha_vdim[31]+1.479019945774904*alpha_cdim[2]+0.3803194146278324*alpha_vdim[1])*f[41]+0.4330127018922193*f[10]*alpha_vdim[31]+0.6614378277661477*alpha_cdim[0]*f[18]+0.3803194146278324*alpha_vdim[11]*f[17]+0.6614378277661477*alpha_cdim[2]*f[10]; 
  out[77] += 0.3872983346207416*f[65]*alpha_vdim[171]+0.4330127018922193*(f[51]*alpha_vdim[161]+f[65]*alpha_vdim[160])+0.3872983346207416*f[68]*alpha_vdim[91]+0.4330127018922193*(f[55]*alpha_vdim[81]+f[68]*alpha_vdim[80]+alpha_cdim[0]*f[70])+0.5809475019311124*alpha_vdim[11]*f[69]+(1.322875655532295*alpha_vdim[11]+1.479019945774904*alpha_vdim[0])*f[61]+(0.3803194146278324*alpha_cdim[2]+1.479019945774904*alpha_vdim[1])*f[42]+0.5809475019311124*alpha_vdim[31]*f[41]+0.5916079783099616*(alpha_vdim[1]*f[41]+alpha_vdim[11]*f[17])+0.6614378277661477*(alpha_vdim[0]*f[17]+alpha_vdim[1]*f[10]); 
  out[78] += 0.3872983346207416*f[66]*alpha_vdim[171]+0.4330127018922193*(f[53]*alpha_vdim[161]+f[66]*alpha_vdim[160])+0.5809475019311124*f[39]*alpha_vdim[111]+(0.5809475019311124*f[67]+1.322875655532295*f[62]+0.5916079783099616*f[16])*alpha_vdim[91]+(1.479019945774904*f[44]+0.5916079783099616*f[39]+0.6614378277661477*f[9])*alpha_vdim[81]+(1.479019945774904*f[62]+0.6614378277661477*f[16])*alpha_vdim[80]+0.4330127018922193*alpha_cdim[0]*f[72]+0.3872983346207416*alpha_vdim[11]*f[71]+0.4330127018922193*(alpha_vdim[0]*f[71]+(alpha_cdim[2]+alpha_vdim[1])*f[56]); 
  out[79] += 0.5809475019311124*f[36]*alpha_vdim[191]+(0.5809475019311124*f[64]+1.322875655532295*f[63]+0.5916079783099616*f[15])*alpha_vdim[171]+(1.479019945774904*f[47]+0.5916079783099616*f[36]+0.6614378277661477*f[7])*alpha_vdim[161]+(1.479019945774904*f[63]+0.6614378277661477*f[15])*alpha_vdim[160]+0.3872983346207416*f[73]*alpha_vdim[91]+0.4330127018922193*(f[58]*alpha_vdim[81]+f[73]*alpha_vdim[80]+alpha_cdim[0]*f[75])+0.3872983346207416*alpha_vdim[11]*f[74]+0.4330127018922193*(alpha_vdim[0]*f[74]+(alpha_cdim[2]+alpha_vdim[1])*f[59]); 

  return alpha_mid; 
} 