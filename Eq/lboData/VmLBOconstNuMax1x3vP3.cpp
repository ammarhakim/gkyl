#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x3vMaxP3(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvz2 = 2.0/dxv[3]; 
  const double rdvzSq4 = 4.0/(dxv[3]*dxv[3]); 

  double alphaDrag[105]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-4.0*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -1.154700538379252*dxv[1]*nuSum*rdvx2; 
  alphaDrag[11] = 2.828427124746191*nuUSum[2]*rdvx2; 
  alphaDrag[31] = 2.828427124746191*nuUSum[3]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[35] = (2.828427124746191*nuUSum[4]-4.0*w[2]*nuSum)*rdvy2; 
  alphaDrag[36] = 2.828427124746191*nuUSum[5]*rdvy2; 
  alphaDrag[38] = -1.154700538379252*dxv[2]*nuSum*rdvy2; 
  alphaDrag[46] = 2.828427124746191*nuUSum[6]*rdvy2; 
  alphaDrag[66] = 2.828427124746191*nuUSum[7]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[70] = (2.828427124746191*nuUSum[8]-4.0*w[3]*nuSum)*rdvz2; 
  alphaDrag[71] = 2.828427124746191*nuUSum[9]*rdvz2; 
  alphaDrag[74] = -1.154700538379252*dxv[3]*nuSum*rdvz2; 
  alphaDrag[81] = 2.828427124746191*nuUSum[10]*rdvz2; 
  alphaDrag[101] = 2.828427124746191*nuUSum[11]*rdvz2; 

  double facDiff[4]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 
  facDiff[3] = nuVtSqSum[3]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alphaDrag[31]*f[31]+alphaDrag[11]*f[11]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[31]*alphaDrag[66]+f[11]*alphaDrag[46]+f[3]*alphaDrag[38]+f[1]*alphaDrag[36]+f[0]*alphaDrag[35]); 
  out[4] += 0.4330127018922193*(f[31]*alphaDrag[101]+f[11]*alphaDrag[81]+f[4]*alphaDrag[74]+f[1]*alphaDrag[71]+f[0]*alphaDrag[70]); 
  out[5] += 0.3803194146278324*(alphaDrag[11]*f[31]+f[11]*alphaDrag[31])+0.3872983346207416*(alphaDrag[1]*f[11]+f[1]*alphaDrag[11])+0.4330127018922193*(alphaDrag[2]*f[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.3803194146278324*f[11]*alphaDrag[66]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alphaDrag[46]+0.4330127018922193*f[6]*alphaDrag[38]+0.3872983346207416*f[11]*alphaDrag[36]+0.4330127018922193*(f[0]*alphaDrag[36]+f[1]*alphaDrag[35]); 
  out[7] += 0.4330127018922193*(f[19]*alphaDrag[46]+f[7]*alphaDrag[38]+f[5]*alphaDrag[36]+f[2]*alphaDrag[35]+alphaDrag[11]*f[21]+alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]); 
  out[8] += 0.3803194146278324*f[11]*alphaDrag[101]+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alphaDrag[81]+0.4330127018922193*f[8]*alphaDrag[74]+0.3872983346207416*f[11]*alphaDrag[71]+0.4330127018922193*(f[0]*alphaDrag[71]+f[1]*alphaDrag[70]); 
  out[9] += 0.4330127018922193*(f[19]*alphaDrag[81]+f[9]*alphaDrag[74]+f[5]*alphaDrag[71]+f[2]*alphaDrag[70]+alphaDrag[11]*f[25]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[21]*alphaDrag[81]+f[10]*alphaDrag[74]+f[6]*alphaDrag[71]+f[3]*alphaDrag[70]+f[25]*alphaDrag[46]+f[10]*alphaDrag[38]+f[8]*alphaDrag[36]+f[4]*alphaDrag[35]); 
  out[12] += 4.743416490252569*(facDiff[3]*f[31]+facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+0.9682458365518543*alphaDrag[11]*f[19]+0.8660254037844386*alphaDrag[2]*f[12]+0.9682458365518543*(alphaDrag[1]*f[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[13] += 4.743416490252569*(facDiff[3]*f[31]+facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+0.9682458365518543*f[21]*alphaDrag[46]+0.8660254037844386*f[13]*alphaDrag[38]+0.9682458365518543*(f[0]*alphaDrag[38]+f[6]*alphaDrag[36]+f[3]*alphaDrag[35]); 
  out[14] += 4.743416490252569*(facDiff[3]*f[31]+facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvzSq4+0.9682458365518543*f[25]*alphaDrag[81]+0.8660254037844386*f[14]*alphaDrag[74]+0.9682458365518543*(f[0]*alphaDrag[74]+f[8]*alphaDrag[71]+f[4]*alphaDrag[70]); 
  out[15] += 0.3803194146278324*f[19]*alphaDrag[66]+0.3872983346207416*f[5]*alphaDrag[46]+0.4330127018922193*f[15]*alphaDrag[38]+0.3872983346207416*f[19]*alphaDrag[36]+0.4330127018922193*(f[2]*alphaDrag[36]+f[5]*alphaDrag[35])+f[21]*(0.3803194146278324*alphaDrag[31]+0.3872983346207416*alphaDrag[1])+0.4330127018922193*alphaDrag[2]*f[15]+0.3872983346207416*f[6]*alphaDrag[11]+0.4330127018922193*(alphaDrag[0]*f[6]+alphaDrag[1]*f[3]); 
  out[16] += 0.3803194146278324*f[19]*alphaDrag[101]+0.3872983346207416*f[5]*alphaDrag[81]+0.4330127018922193*f[16]*alphaDrag[74]+0.3872983346207416*f[19]*alphaDrag[71]+0.4330127018922193*(f[2]*alphaDrag[71]+f[5]*alphaDrag[70])+f[25]*(0.3803194146278324*alphaDrag[31]+0.3872983346207416*alphaDrag[1])+0.4330127018922193*alphaDrag[2]*f[16]+0.3872983346207416*f[8]*alphaDrag[11]+0.4330127018922193*(alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[17] += 0.3803194146278324*f[21]*alphaDrag[101]+0.3872983346207416*f[6]*alphaDrag[81]+0.4330127018922193*f[17]*alphaDrag[74]+0.3872983346207416*f[21]*alphaDrag[71]+0.4330127018922193*(f[3]*alphaDrag[71]+f[6]*alphaDrag[70])+0.3803194146278324*f[25]*alphaDrag[66]+0.3872983346207416*f[8]*alphaDrag[46]+0.4330127018922193*f[17]*alphaDrag[38]+0.3872983346207416*f[25]*alphaDrag[36]+0.4330127018922193*(f[4]*alphaDrag[36]+f[8]*alphaDrag[35]); 
  out[18] += 0.4330127018922193*(f[18]*alphaDrag[74]+f[15]*alphaDrag[71]+f[7]*alphaDrag[70]+f[18]*alphaDrag[38]+f[16]*alphaDrag[36]+f[9]*alphaDrag[35]+alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[0]*f[10]); 
  out[19] += 0.2581988897471612*alphaDrag[31]*f[31]+0.3803194146278324*(alphaDrag[1]*f[31]+f[1]*alphaDrag[31])+0.4330127018922193*alphaDrag[2]*f[19]+0.276641667586244*alphaDrag[11]*f[11]+0.4330127018922193*(alphaDrag[0]*f[11]+f[0]*alphaDrag[11])+0.3872983346207416*alphaDrag[1]*f[1]; 
  out[20] += (4.166190448976479*facDiff[2]*f[31]+4.242640687119286*(facDiff[1]*f[11]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[11]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+0.8504200642707612*f[19]*alphaDrag[31]+0.8660254037844386*(alphaDrag[2]*f[20]+alphaDrag[1]*f[19]+f[5]*alphaDrag[11])+0.9682458365518543*(alphaDrag[0]*f[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[21] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alphaDrag[66]+0.276641667586244*f[11]*alphaDrag[46]+0.4330127018922193*(f[0]*alphaDrag[46]+f[21]*alphaDrag[38])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alphaDrag[36]+0.4330127018922193*f[11]*alphaDrag[35]; 
  out[22] += 4.743416490252569*(facDiff[2]*f[21]+facDiff[1]*f[6]+facDiff[0]*f[3])*rdvxSq4+0.4330127018922193*(f[22]*alphaDrag[38]+f[20]*alphaDrag[36]+f[12]*alphaDrag[35])+0.8660254037844386*alphaDrag[2]*f[22]+0.9682458365518543*(alphaDrag[1]*f[15]+alphaDrag[0]*f[7]+alphaDrag[2]*f[3]); 
  out[23] += (4.166190448976479*facDiff[2]*f[31]+4.242640687119286*(facDiff[1]*f[11]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[11]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4+0.8504200642707612*f[21]*alphaDrag[66]+0.8660254037844386*f[6]*alphaDrag[46]+(0.8660254037844386*f[23]+0.9682458365518543*f[1])*alphaDrag[38]+0.8660254037844386*f[21]*alphaDrag[36]+0.9682458365518543*(f[3]*alphaDrag[36]+f[6]*alphaDrag[35]); 
  out[24] += 4.743416490252569*(facDiff[2]*f[19]+facDiff[1]*f[5]+facDiff[0]*f[2])*rdvySq4+0.8660254037844386*f[24]*alphaDrag[38]+0.9682458365518543*(f[2]*alphaDrag[38]+f[15]*alphaDrag[36]+f[7]*alphaDrag[35])+0.4330127018922193*(alphaDrag[2]*f[24]+alphaDrag[1]*f[23]+alphaDrag[0]*f[13]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alphaDrag[101]+0.276641667586244*f[11]*alphaDrag[81]+0.4330127018922193*(f[0]*alphaDrag[81]+f[25]*alphaDrag[74])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alphaDrag[71]+0.4330127018922193*f[11]*alphaDrag[70]; 
  out[26] += 4.743416490252569*(facDiff[2]*f[25]+facDiff[1]*f[8]+facDiff[0]*f[4])*rdvxSq4+0.4330127018922193*(f[26]*alphaDrag[74]+f[20]*alphaDrag[71]+f[12]*alphaDrag[70])+0.8660254037844386*alphaDrag[2]*f[26]+0.9682458365518543*(alphaDrag[1]*f[16]+alphaDrag[0]*f[9]+alphaDrag[2]*f[4]); 
  out[27] += 4.743416490252569*(facDiff[2]*f[25]+facDiff[1]*f[8]+facDiff[0]*f[4])*rdvySq4+0.4330127018922193*(f[27]*alphaDrag[74]+f[23]*alphaDrag[71]+f[13]*alphaDrag[70])+0.8660254037844386*f[27]*alphaDrag[38]+0.9682458365518543*(f[4]*alphaDrag[38]+f[17]*alphaDrag[36]+f[10]*alphaDrag[35]); 
  out[28] += (4.166190448976479*facDiff[2]*f[31]+4.242640687119286*(facDiff[1]*f[11]+f[1]*facDiff[2])+4.166190448976479*facDiff[3]*f[11]+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvzSq4+0.8504200642707612*f[25]*alphaDrag[101]+0.8660254037844386*f[8]*alphaDrag[81]+(0.8660254037844386*f[28]+0.9682458365518543*f[1])*alphaDrag[74]+0.8660254037844386*f[25]*alphaDrag[71]+0.9682458365518543*(f[4]*alphaDrag[71]+f[8]*alphaDrag[70]); 
  out[29] += 4.743416490252569*(facDiff[2]*f[19]+facDiff[1]*f[5]+facDiff[0]*f[2])*rdvzSq4+0.8660254037844386*f[29]*alphaDrag[74]+0.9682458365518543*(f[2]*alphaDrag[74]+f[16]*alphaDrag[71]+f[9]*alphaDrag[70])+0.4330127018922193*(alphaDrag[2]*f[29]+alphaDrag[1]*f[28]+alphaDrag[0]*f[14]); 
  out[30] += 4.743416490252569*(facDiff[2]*f[21]+facDiff[1]*f[6]+facDiff[0]*f[3])*rdvzSq4+0.8660254037844386*f[30]*alphaDrag[74]+0.9682458365518543*(f[3]*alphaDrag[74]+f[17]*alphaDrag[71]+f[10]*alphaDrag[70])+0.4330127018922193*(f[30]*alphaDrag[38]+f[28]*alphaDrag[36]+f[14]*alphaDrag[35]); 
  out[32] += 16.20185174601965*(facDiff[2]*f[19]+facDiff[1]*f[5]+facDiff[0]*f[2])*rdvxSq4+1.299038105676658*alphaDrag[2]*f[32]+0.6614378277661477*alphaDrag[31]*f[31]+1.479019945774904*(alphaDrag[1]*f[20]+alphaDrag[0]*f[12])+0.6614378277661477*alphaDrag[11]*f[11]+1.984313483298443*alphaDrag[2]*f[2]+0.6614378277661477*(alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[33] += 16.20185174601965*(facDiff[2]*f[21]+facDiff[1]*f[6]+facDiff[0]*f[3])*rdvySq4+0.6614378277661477*(f[31]*alphaDrag[66]+f[11]*alphaDrag[46])+(1.299038105676658*f[33]+1.984313483298443*f[3])*alphaDrag[38]+(1.479019945774904*f[23]+0.6614378277661477*f[1])*alphaDrag[36]+(1.479019945774904*f[13]+0.6614378277661477*f[0])*alphaDrag[35]; 
  out[34] += 16.20185174601965*(facDiff[2]*f[25]+facDiff[1]*f[8]+facDiff[0]*f[4])*rdvzSq4+0.6614378277661477*(f[31]*alphaDrag[101]+f[11]*alphaDrag[81])+(1.299038105676658*f[34]+1.984313483298443*f[4])*alphaDrag[74]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alphaDrag[71]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alphaDrag[70]; 

  return std::abs(0.125*alphaDrag[0]-0.1397542485937369*alphaDrag[11])+std::abs(0.125*alphaDrag[35]-0.1397542485937369*alphaDrag[46])+std::abs(0.125*alphaDrag[70]-0.1397542485937369*alphaDrag[81])+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvxSq4)+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvySq4)+std::abs((1.616244071283538*facDiff[0]-1.807015805810503*facDiff[2])*rdvzSq4); 

} 
