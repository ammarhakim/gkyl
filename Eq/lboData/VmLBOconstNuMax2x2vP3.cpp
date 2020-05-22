#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x2vMaxP3(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[2]; 
  const double rdvxSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvy2 = 2.0/dxv[3]; 
  const double rdvySq4 = 4.0/(dxv[3]*dxv[3]); 

  double alphaDrag[70]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-4.0*w[2]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.0*nuUSum[2]*rdvx2; 
  alphaDrag[3] = -1.154700538379252*dxv[2]*nuSum*rdvx2; 
  alphaDrag[5] = 2.0*nuUSum[3]*rdvx2; 
  alphaDrag[11] = 2.0*nuUSum[4]*rdvx2; 
  alphaDrag[12] = 2.0*nuUSum[5]*rdvx2; 
  alphaDrag[19] = 2.0*nuUSum[6]*rdvx2; 
  alphaDrag[20] = 2.0*nuUSum[7]*rdvx2; 
  alphaDrag[31] = 2.0*nuUSum[8]*rdvx2; 
  alphaDrag[32] = 2.0*nuUSum[9]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[35] = (2.0*nuUSum[10]-4.0*w[3]*nuSum)*rdvy2; 
  alphaDrag[36] = 2.0*nuUSum[11]*rdvy2; 
  alphaDrag[37] = 2.0*nuUSum[12]*rdvy2; 
  alphaDrag[39] = -1.154700538379252*dxv[3]*nuSum*rdvy2; 
  alphaDrag[40] = 2.0*nuUSum[13]*rdvy2; 
  alphaDrag[46] = 2.0*nuUSum[14]*rdvy2; 
  alphaDrag[47] = 2.0*nuUSum[15]*rdvy2; 
  alphaDrag[54] = 2.0*nuUSum[16]*rdvy2; 
  alphaDrag[55] = 2.0*nuUSum[17]*rdvy2; 
  alphaDrag[66] = 2.0*nuUSum[18]*rdvy2; 
  alphaDrag[67] = 2.0*nuUSum[19]*rdvy2; 

  double facDiff[10]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 
  facDiff[3] = nuVtSqSum[3]; 
  facDiff[4] = nuVtSqSum[4]; 
  facDiff[5] = nuVtSqSum[5]; 
  facDiff[6] = nuVtSqSum[6]; 
  facDiff[7] = nuVtSqSum[7]; 
  facDiff[8] = nuVtSqSum[8]; 
  facDiff[9] = nuVtSqSum[9]; 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.4330127018922193*(alphaDrag[32]*f[32]+alphaDrag[31]*f[31]+alphaDrag[20]*f[20]+alphaDrag[19]*f[19]+alphaDrag[12]*f[12]+alphaDrag[11]*f[11]+alphaDrag[5]*f[5]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[32]*alphaDrag[67]+f[31]*alphaDrag[66]+f[20]*alphaDrag[55]+f[19]*alphaDrag[54]+f[12]*alphaDrag[47]+f[11]*alphaDrag[46]+f[5]*alphaDrag[40]+f[4]*alphaDrag[39]+f[2]*alphaDrag[37]+f[1]*alphaDrag[36]+f[0]*alphaDrag[35]); 
  out[6] += 0.3803194146278324*(alphaDrag[11]*f[31]+f[11]*alphaDrag[31])+0.4330127018922193*(alphaDrag[12]*f[20]+f[12]*alphaDrag[20])+0.3872983346207416*(alphaDrag[5]*f[19]+f[5]*alphaDrag[19]+alphaDrag[1]*f[11]+f[1]*alphaDrag[11])+0.4330127018922193*(alphaDrag[3]*f[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[7] += 0.3803194146278324*(alphaDrag[12]*f[32]+f[12]*alphaDrag[32])+0.3872983346207416*(alphaDrag[5]*f[20]+f[5]*alphaDrag[20])+0.4330127018922193*(alphaDrag[11]*f[19]+f[11]*alphaDrag[19])+0.3872983346207416*(alphaDrag[2]*f[12]+f[2]*alphaDrag[12])+0.4330127018922193*(alphaDrag[3]*f[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[8] += 0.3803194146278324*f[11]*alphaDrag[66]+0.4330127018922193*f[12]*alphaDrag[55]+0.3872983346207416*f[5]*alphaDrag[54]+0.4330127018922193*f[20]*alphaDrag[47]+0.3803194146278324*f[31]*alphaDrag[46]+0.3872983346207416*(f[1]*alphaDrag[46]+f[19]*alphaDrag[40])+0.4330127018922193*(f[2]*alphaDrag[40]+f[8]*alphaDrag[39]+f[5]*alphaDrag[37])+0.3872983346207416*f[11]*alphaDrag[36]+0.4330127018922193*(f[0]*alphaDrag[36]+f[1]*alphaDrag[35]); 
  out[9] += 0.3803194146278324*f[12]*alphaDrag[67]+0.3872983346207416*f[5]*alphaDrag[55]+0.4330127018922193*f[11]*alphaDrag[54]+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alphaDrag[47]+0.4330127018922193*f[19]*alphaDrag[46]+0.3872983346207416*f[20]*alphaDrag[40]+0.4330127018922193*(f[1]*alphaDrag[40]+f[9]*alphaDrag[39])+0.3872983346207416*f[12]*alphaDrag[37]+0.4330127018922193*(f[0]*alphaDrag[37]+f[5]*alphaDrag[36]+f[2]*alphaDrag[35]); 
  out[10] += 0.4330127018922193*(f[22]*alphaDrag[47]+f[21]*alphaDrag[46]+f[15]*alphaDrag[40]+f[10]*alphaDrag[39]+f[7]*alphaDrag[37]+f[6]*alphaDrag[36]+f[3]*alphaDrag[35]+alphaDrag[12]*f[26]+alphaDrag[11]*f[25]+alphaDrag[5]*f[16]+alphaDrag[3]*f[10]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[13] += 3.354101966249685*(facDiff[9]*f[32]+facDiff[8]*f[31]+facDiff[7]*f[20]+facDiff[6]*f[19]+facDiff[5]*f[12]+facDiff[4]*f[11]+facDiff[3]*f[5]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+0.9682458365518543*(alphaDrag[12]*f[22]+alphaDrag[11]*f[21]+alphaDrag[5]*f[15])+0.8660254037844386*alphaDrag[3]*f[13]+0.9682458365518543*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[14] += 3.354101966249685*(facDiff[9]*f[32]+facDiff[8]*f[31]+facDiff[7]*f[20]+facDiff[6]*f[19]+facDiff[5]*f[12]+facDiff[4]*f[11]+facDiff[3]*f[5]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+0.9682458365518543*(f[26]*alphaDrag[47]+f[25]*alphaDrag[46]+f[16]*alphaDrag[40])+0.8660254037844386*f[14]*alphaDrag[39]+0.9682458365518543*(f[0]*alphaDrag[39]+f[9]*alphaDrag[37]+f[8]*alphaDrag[36]+f[4]*alphaDrag[35]); 
  out[15] += 0.3803194146278324*(alphaDrag[20]*f[32]+f[20]*alphaDrag[32]+alphaDrag[19]*f[31]+f[19]*alphaDrag[31])+(0.3464101615137755*alphaDrag[19]+0.3872983346207416*alphaDrag[2])*f[20]+0.3464101615137755*f[19]*alphaDrag[20]+0.3872983346207416*(f[2]*alphaDrag[20]+alphaDrag[1]*f[19]+f[1]*alphaDrag[19])+0.4330127018922193*alphaDrag[3]*f[15]+0.3872983346207416*(alphaDrag[5]*f[12]+f[5]*alphaDrag[12]+alphaDrag[5]*f[11]+f[5]*alphaDrag[11])+0.4330127018922193*(alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[16] += 0.3803194146278324*(f[20]*alphaDrag[67]+f[19]*alphaDrag[66])+(0.3803194146278324*f[32]+0.3464101615137755*f[19]+0.3872983346207416*f[2])*alphaDrag[55]+(0.3803194146278324*f[31]+0.3464101615137755*f[20])*alphaDrag[54]+0.3872983346207416*(f[1]*alphaDrag[54]+f[5]*(alphaDrag[47]+alphaDrag[46])+(f[12]+f[11])*alphaDrag[40])+0.4330127018922193*(f[0]*alphaDrag[40]+f[16]*alphaDrag[39])+(0.3872983346207416*f[20]+0.4330127018922193*f[1])*alphaDrag[37]+0.3872983346207416*f[19]*alphaDrag[36]+0.4330127018922193*(f[2]*alphaDrag[36]+f[5]*alphaDrag[35]); 
  out[17] += 0.3803194146278324*f[21]*alphaDrag[66]+0.4330127018922193*f[22]*alphaDrag[55]+0.3872983346207416*(f[15]*alphaDrag[54]+f[6]*alphaDrag[46])+0.4330127018922193*(f[7]*alphaDrag[40]+f[17]*alphaDrag[39]+f[15]*alphaDrag[37])+0.3872983346207416*f[21]*alphaDrag[36]+0.4330127018922193*(f[3]*alphaDrag[36]+f[6]*alphaDrag[35])+0.3803194146278324*f[25]*alphaDrag[31]+0.4330127018922193*alphaDrag[20]*f[26]+0.3872983346207416*(alphaDrag[1]*f[25]+f[16]*alphaDrag[19])+0.4330127018922193*(alphaDrag[3]*f[17]+alphaDrag[2]*f[16])+0.3872983346207416*f[8]*alphaDrag[11]+0.4330127018922193*(alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[18] += 0.3803194146278324*f[22]*alphaDrag[67]+0.3872983346207416*f[15]*alphaDrag[55]+0.4330127018922193*f[21]*alphaDrag[54]+0.3872983346207416*f[7]*alphaDrag[47]+0.4330127018922193*(f[6]*alphaDrag[40]+f[18]*alphaDrag[39])+0.3872983346207416*f[22]*alphaDrag[37]+0.4330127018922193*(f[3]*alphaDrag[37]+f[15]*alphaDrag[36]+f[7]*alphaDrag[35])+f[26]*(0.3803194146278324*alphaDrag[32]+0.3872983346207416*alphaDrag[2])+0.4330127018922193*alphaDrag[19]*f[25]+0.3872983346207416*f[16]*alphaDrag[20]+0.4330127018922193*(alphaDrag[3]*f[18]+alphaDrag[1]*f[16])+0.3872983346207416*f[9]*alphaDrag[12]+0.4330127018922193*(alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[21] += 0.2581988897471612*alphaDrag[31]*f[31]+0.3803194146278324*(alphaDrag[1]*f[31]+f[1]*alphaDrag[31])+0.4330127018922193*alphaDrag[3]*f[21]+0.3872983346207416*alphaDrag[20]*f[20]+0.276641667586244*alphaDrag[19]*f[19]+0.4330127018922193*(alphaDrag[2]*f[19]+f[2]*alphaDrag[19])+0.276641667586244*alphaDrag[11]*f[11]+0.4330127018922193*(alphaDrag[0]*f[11]+f[0]*alphaDrag[11])+0.3872983346207416*(alphaDrag[5]*f[5]+alphaDrag[1]*f[1]); 
  out[22] += 0.2581988897471612*alphaDrag[32]*f[32]+0.3803194146278324*(alphaDrag[2]*f[32]+f[2]*alphaDrag[32])+0.4330127018922193*alphaDrag[3]*f[22]+0.276641667586244*alphaDrag[20]*f[20]+0.4330127018922193*(alphaDrag[1]*f[20]+f[1]*alphaDrag[20])+0.3872983346207416*alphaDrag[19]*f[19]+0.276641667586244*alphaDrag[12]*f[12]+0.4330127018922193*(alphaDrag[0]*f[12]+f[0]*alphaDrag[12])+0.3872983346207416*(alphaDrag[5]*f[5]+alphaDrag[2]*f[2]); 
  out[23] += (2.945941518185896*facDiff[4]*f[31]+3.354101966249685*facDiff[5]*f[20]+3.0*facDiff[3]*f[19]+3.354101966249685*facDiff[7]*f[12]+(2.945941518185896*facDiff[8]+3.0*facDiff[1])*f[11]+f[5]*(3.0*facDiff[6]+3.354101966249685*facDiff[2])+3.0*f[1]*facDiff[4]+3.354101966249685*(f[2]*facDiff[3]+f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+0.8504200642707612*f[21]*alphaDrag[31]+0.8660254037844386*alphaDrag[3]*f[23]+0.9682458365518543*alphaDrag[20]*f[22]+0.8660254037844386*alphaDrag[1]*f[21]+f[15]*(0.8660254037844386*alphaDrag[19]+0.9682458365518543*alphaDrag[2])+0.8660254037844386*f[6]*alphaDrag[11]+0.9682458365518543*(alphaDrag[5]*f[7]+alphaDrag[0]*f[6]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[24] += (2.945941518185896*facDiff[5]*f[32]+3.0*facDiff[3]*f[20]+3.354101966249685*facDiff[4]*f[19]+(2.945941518185896*facDiff[9]+3.0*facDiff[2])*f[12]+3.354101966249685*facDiff[6]*f[11]+3.0*(f[5]*facDiff[7]+f[2]*facDiff[5])+3.354101966249685*(facDiff[1]*f[5]+f[1]*facDiff[3]+f[0]*facDiff[2]+facDiff[0]*f[2]))*rdvxSq4+0.8504200642707612*f[22]*alphaDrag[32]+0.8660254037844386*(alphaDrag[3]*f[24]+alphaDrag[2]*f[22])+0.9682458365518543*alphaDrag[19]*f[21]+f[15]*(0.8660254037844386*alphaDrag[20]+0.9682458365518543*alphaDrag[1])+0.8660254037844386*f[7]*alphaDrag[12]+0.9682458365518543*(alphaDrag[0]*f[7]+alphaDrag[5]*f[6]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[25] += (0.2581988897471612*f[31]+0.3803194146278324*f[1])*alphaDrag[66]+0.3872983346207416*f[20]*alphaDrag[55]+(0.276641667586244*f[19]+0.4330127018922193*f[2])*alphaDrag[54]+(0.276641667586244*f[11]+0.4330127018922193*f[0])*alphaDrag[46]+0.3872983346207416*f[5]*alphaDrag[40]+0.4330127018922193*(f[25]*alphaDrag[39]+f[19]*alphaDrag[37])+(0.3803194146278324*f[31]+0.3872983346207416*f[1])*alphaDrag[36]+0.4330127018922193*f[11]*alphaDrag[35]; 
  out[26] += (0.2581988897471612*f[32]+0.3803194146278324*f[2])*alphaDrag[67]+(0.276641667586244*f[20]+0.4330127018922193*f[1])*alphaDrag[55]+0.3872983346207416*f[19]*alphaDrag[54]+(0.276641667586244*f[12]+0.4330127018922193*f[0])*alphaDrag[47]+0.3872983346207416*f[5]*alphaDrag[40]+0.4330127018922193*f[26]*alphaDrag[39]+(0.3803194146278324*f[32]+0.3872983346207416*f[2])*alphaDrag[37]+0.4330127018922193*(f[20]*alphaDrag[36]+f[12]*alphaDrag[35]); 
  out[27] += 3.354101966249685*(facDiff[5]*f[26]+facDiff[4]*f[25]+facDiff[3]*f[16]+facDiff[2]*f[9]+facDiff[1]*f[8]+facDiff[0]*f[4])*rdvxSq4+0.4330127018922193*(f[27]*alphaDrag[39]+f[24]*alphaDrag[37]+f[23]*alphaDrag[36]+f[13]*alphaDrag[35])+0.8660254037844386*alphaDrag[3]*f[27]+0.9682458365518543*(alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[0]*f[10]+alphaDrag[3]*f[4]); 
  out[28] += (2.945941518185896*facDiff[4]*f[31]+3.354101966249685*facDiff[5]*f[20]+3.0*facDiff[3]*f[19]+3.354101966249685*facDiff[7]*f[12]+(2.945941518185896*facDiff[8]+3.0*facDiff[1])*f[11]+f[5]*(3.0*facDiff[6]+3.354101966249685*facDiff[2])+3.0*f[1]*facDiff[4]+3.354101966249685*(f[2]*facDiff[3]+f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4+0.8504200642707612*f[25]*alphaDrag[66]+0.9682458365518543*f[26]*alphaDrag[55]+0.8660254037844386*(f[16]*alphaDrag[54]+f[8]*alphaDrag[46])+0.9682458365518543*f[9]*alphaDrag[40]+0.8660254037844386*f[28]*alphaDrag[39]+0.9682458365518543*(f[1]*alphaDrag[39]+f[16]*alphaDrag[37])+0.8660254037844386*f[25]*alphaDrag[36]+0.9682458365518543*(f[4]*alphaDrag[36]+f[8]*alphaDrag[35]); 
  out[29] += (2.945941518185896*facDiff[5]*f[32]+3.0*facDiff[3]*f[20]+3.354101966249685*facDiff[4]*f[19]+(2.945941518185896*facDiff[9]+3.0*facDiff[2])*f[12]+3.354101966249685*facDiff[6]*f[11]+3.0*(f[5]*facDiff[7]+f[2]*facDiff[5])+3.354101966249685*(facDiff[1]*f[5]+f[1]*facDiff[3]+f[0]*facDiff[2]+facDiff[0]*f[2]))*rdvySq4+0.8504200642707612*f[26]*alphaDrag[67]+0.8660254037844386*f[16]*alphaDrag[55]+0.9682458365518543*f[25]*alphaDrag[54]+0.8660254037844386*f[9]*alphaDrag[47]+0.9682458365518543*f[8]*alphaDrag[40]+(0.8660254037844386*f[29]+0.9682458365518543*f[2])*alphaDrag[39]+0.8660254037844386*f[26]*alphaDrag[37]+0.9682458365518543*(f[4]*alphaDrag[37]+f[16]*alphaDrag[36]+f[9]*alphaDrag[35]); 
  out[30] += 3.354101966249685*(facDiff[5]*f[22]+facDiff[4]*f[21]+facDiff[3]*f[15]+facDiff[2]*f[7]+facDiff[1]*f[6]+facDiff[0]*f[3])*rdvySq4+0.8660254037844386*f[30]*alphaDrag[39]+0.9682458365518543*(f[3]*alphaDrag[39]+f[18]*alphaDrag[37]+f[17]*alphaDrag[36]+f[10]*alphaDrag[35])+0.4330127018922193*(alphaDrag[3]*f[30]+alphaDrag[2]*f[29]+alphaDrag[1]*f[28]+alphaDrag[0]*f[14]); 
  out[33] += 11.4564392373896*(facDiff[5]*f[22]+facDiff[4]*f[21]+facDiff[3]*f[15]+facDiff[2]*f[7]+facDiff[1]*f[6]+facDiff[0]*f[3])*rdvxSq4+1.299038105676658*alphaDrag[3]*f[33]+0.6614378277661477*(alphaDrag[32]*f[32]+alphaDrag[31]*f[31])+1.479019945774904*(alphaDrag[2]*f[24]+alphaDrag[1]*f[23])+0.6614378277661477*(alphaDrag[20]*f[20]+alphaDrag[19]*f[19])+1.479019945774904*alphaDrag[0]*f[13]+0.6614378277661477*(alphaDrag[12]*f[12]+alphaDrag[11]*f[11]+alphaDrag[5]*f[5])+1.984313483298443*alphaDrag[3]*f[3]+0.6614378277661477*(alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[34] += 11.4564392373896*(facDiff[5]*f[26]+facDiff[4]*f[25]+facDiff[3]*f[16]+facDiff[2]*f[9]+facDiff[1]*f[8]+facDiff[0]*f[4])*rdvySq4+0.6614378277661477*(f[32]*alphaDrag[67]+f[31]*alphaDrag[66]+f[20]*alphaDrag[55]+f[19]*alphaDrag[54]+f[12]*alphaDrag[47]+f[11]*alphaDrag[46]+f[5]*alphaDrag[40])+(1.299038105676658*f[34]+1.984313483298443*f[4])*alphaDrag[39]+(1.479019945774904*f[29]+0.6614378277661477*f[2])*alphaDrag[37]+(1.479019945774904*f[28]+0.6614378277661477*f[1])*alphaDrag[36]+(1.479019945774904*f[14]+0.6614378277661477*f[0])*alphaDrag[35]; 

  return std::abs(0.125*alphaDrag[0]-0.1397542485937369*(alphaDrag[12]+alphaDrag[11]))+std::abs(0.125*alphaDrag[35]-0.1397542485937369*(alphaDrag[47]+alphaDrag[46]))+std::abs((1.142857142857143*facDiff[0]-1.27775312999988*(facDiff[5]+facDiff[4]))*rdvxSq4)+std::abs((1.142857142857143*facDiff[0]-1.27775312999988*(facDiff[5]+facDiff[4]))*rdvySq4); 

} 
