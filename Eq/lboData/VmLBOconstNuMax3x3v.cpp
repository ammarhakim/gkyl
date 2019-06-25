#include <VmLBOModDecl.h> 
double VmLBOconstNuVol3x3vMaxP1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[6]:      Cell-center coordinates. 
  // dxv[6]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[3]; 
  const double rdvxSq4 = 4.0/(dxv[3]*dxv[3]); 
  const double rdvy2 = 2.0/dxv[4]; 
  const double rdvySq4 = 4.0/(dxv[4]*dxv[4]); 
  const double rdvz2 = 2.0/dxv[5]; 
  const double rdvzSq4 = 4.0/(dxv[5]*dxv[5]); 

  double alphaDrag[21]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = rdvx2*(2.828427124746191*nuUSum[0]-8.0*w[3]*nuSum); 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.828427124746191*nuUSum[2]*rdvx2; 
  alphaDrag[3] = 2.828427124746191*nuUSum[3]*rdvx2; 
  alphaDrag[4] = -2.309401076758503*dxv[3]*rdvx2*nuSum; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[7] = rdvy2*(2.828427124746191*nuUSum[4]-8.0*w[4]*nuSum); 
  alphaDrag[8] = 2.828427124746191*nuUSum[5]*rdvy2; 
  alphaDrag[9] = 2.828427124746191*nuUSum[6]*rdvy2; 
  alphaDrag[10] = 2.828427124746191*nuUSum[7]*rdvy2; 
  alphaDrag[12] = -2.309401076758503*dxv[4]*rdvy2*nuSum; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[14] = rdvz2*(2.828427124746191*nuUSum[8]-8.0*w[5]*nuSum); 
  alphaDrag[15] = 2.828427124746191*nuUSum[9]*rdvz2; 
  alphaDrag[16] = 2.828427124746191*nuUSum[10]*rdvz2; 
  alphaDrag[17] = 2.828427124746191*nuUSum[11]*rdvz2; 
  alphaDrag[20] = -2.309401076758503*dxv[5]*rdvz2*nuSum; 

  // Put together updates due to drag and diffusion terms.
  out[4] += 0.2165063509461096*(f[4]*alphaDrag[4]+f[3]*alphaDrag[3]+f[2]*alphaDrag[2]+f[1]*alphaDrag[1]+f[0]*alphaDrag[0]); 
  out[5] += 0.2165063509461096*(f[5]*alphaDrag[12]+f[3]*alphaDrag[10]+f[2]*alphaDrag[9]+f[1]*alphaDrag[8]+f[0]*alphaDrag[7]); 
  out[6] += 0.2165063509461096*(f[6]*alphaDrag[20]+f[3]*alphaDrag[17]+f[2]*alphaDrag[16]+f[1]*alphaDrag[15]+f[0]*alphaDrag[14]); 

  return std::abs(0.0625*alphaDrag[0])+std::abs(0.0625*alphaDrag[7])+std::abs(0.0625*alphaDrag[14])+std::abs(0.4714045207910317*nuVtSqSum[0]*rdvxSq4)+std::abs(0.4714045207910317*nuVtSqSum[0]*rdvySq4)+std::abs(0.4714045207910317*nuVtSqSum[0]*rdvzSq4); 

} 
double VmLBOconstNuVol3x3vMaxP2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[6]:      Cell-center coordinates. 
  // dxv[6]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[3]; 
  const double rdvxSq4 = 4.0/(dxv[3]*dxv[3]); 
  const double rdvy2 = 2.0/dxv[4]; 
  const double rdvySq4 = 4.0/(dxv[4]*dxv[4]); 
  const double rdvz2 = 2.0/dxv[5]; 
  const double rdvzSq4 = 4.0/(dxv[5]*dxv[5]); 

  double alphaDrag[84]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = rdvx2*(2.828427124746191*nuUSum[0]-8.0*w[3]*nuSum); 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.828427124746191*nuUSum[2]*rdvx2; 
  alphaDrag[3] = 2.828427124746191*nuUSum[3]*rdvx2; 
  alphaDrag[4] = -2.309401076758503*dxv[3]*rdvx2*nuSum; 
  alphaDrag[7] = 2.828427124746191*nuUSum[4]*rdvx2; 
  alphaDrag[8] = 2.828427124746191*nuUSum[5]*rdvx2; 
  alphaDrag[9] = 2.828427124746191*nuUSum[6]*rdvx2; 
  alphaDrag[22] = 2.828427124746191*nuUSum[7]*rdvx2; 
  alphaDrag[23] = 2.828427124746191*nuUSum[8]*rdvx2; 
  alphaDrag[24] = 2.828427124746191*nuUSum[9]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[28] = rdvy2*(2.828427124746191*nuUSum[10]-8.0*w[4]*nuSum); 
  alphaDrag[29] = 2.828427124746191*nuUSum[11]*rdvy2; 
  alphaDrag[30] = 2.828427124746191*nuUSum[12]*rdvy2; 
  alphaDrag[31] = 2.828427124746191*nuUSum[13]*rdvy2; 
  alphaDrag[33] = -2.309401076758503*dxv[4]*rdvy2*nuSum; 
  alphaDrag[35] = 2.828427124746191*nuUSum[14]*rdvy2; 
  alphaDrag[36] = 2.828427124746191*nuUSum[15]*rdvy2; 
  alphaDrag[37] = 2.828427124746191*nuUSum[16]*rdvy2; 
  alphaDrag[50] = 2.828427124746191*nuUSum[17]*rdvy2; 
  alphaDrag[51] = 2.828427124746191*nuUSum[18]*rdvy2; 
  alphaDrag[52] = 2.828427124746191*nuUSum[19]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[56] = rdvz2*(2.828427124746191*nuUSum[20]-8.0*w[5]*nuSum); 
  alphaDrag[57] = 2.828427124746191*nuUSum[21]*rdvz2; 
  alphaDrag[58] = 2.828427124746191*nuUSum[22]*rdvz2; 
  alphaDrag[59] = 2.828427124746191*nuUSum[23]*rdvz2; 
  alphaDrag[62] = -2.309401076758503*dxv[5]*rdvz2*nuSum; 
  alphaDrag[63] = 2.828427124746191*nuUSum[24]*rdvz2; 
  alphaDrag[64] = 2.828427124746191*nuUSum[25]*rdvz2; 
  alphaDrag[65] = 2.828427124746191*nuUSum[26]*rdvz2; 
  alphaDrag[78] = 2.828427124746191*nuUSum[27]*rdvz2; 
  alphaDrag[79] = 2.828427124746191*nuUSum[28]*rdvz2; 
  alphaDrag[80] = 2.828427124746191*nuUSum[29]*rdvz2; 

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
  out[4] += 0.2165063509461096*(f[24]*alphaDrag[24]+f[23]*alphaDrag[23]+f[22]*alphaDrag[22]+f[9]*alphaDrag[9]+f[8]*alphaDrag[8]+f[7]*alphaDrag[7]+f[4]*alphaDrag[4]+f[3]*alphaDrag[3]+f[2]*alphaDrag[2]+f[1]*alphaDrag[1]+f[0]*alphaDrag[0]); 
  out[5] += 0.2165063509461096*(f[24]*alphaDrag[52]+f[23]*alphaDrag[51]+f[22]*alphaDrag[50]+f[9]*alphaDrag[37]+f[8]*alphaDrag[36]+f[7]*alphaDrag[35]+f[5]*alphaDrag[33]+f[3]*alphaDrag[31]+f[2]*alphaDrag[30]+f[1]*alphaDrag[29]+f[0]*alphaDrag[28]); 
  out[6] += 0.2165063509461096*(f[24]*alphaDrag[80]+f[23]*alphaDrag[79]+f[22]*alphaDrag[78]+f[9]*alphaDrag[65]+f[8]*alphaDrag[64]+f[7]*alphaDrag[63]+f[6]*alphaDrag[62]+f[3]*alphaDrag[59]+f[2]*alphaDrag[58]+f[1]*alphaDrag[57]+f[0]*alphaDrag[56]); 
  out[10] += 0.1936491673103708*(f[1]*alphaDrag[22]+alphaDrag[1]*f[22])+0.2165063509461096*(alphaDrag[4]*f[10]+f[3]*alphaDrag[8]+alphaDrag[3]*f[8]+f[2]*alphaDrag[7]+alphaDrag[2]*f[7]+f[0]*alphaDrag[1]+alphaDrag[0]*f[1]); 
  out[11] += 0.1936491673103708*(f[2]*alphaDrag[23]+alphaDrag[2]*f[23])+0.2165063509461096*(alphaDrag[4]*f[11]+f[3]*alphaDrag[9]+alphaDrag[3]*f[9]+f[1]*alphaDrag[7]+alphaDrag[1]*f[7]+f[0]*alphaDrag[2]+alphaDrag[0]*f[2]); 
  out[12] += 0.1936491673103708*(f[3]*alphaDrag[24]+alphaDrag[3]*f[24])+0.2165063509461096*(alphaDrag[4]*f[12]+f[2]*alphaDrag[9]+alphaDrag[2]*f[9]+f[1]*alphaDrag[8]+alphaDrag[1]*f[8]+f[0]*alphaDrag[3]+alphaDrag[0]*f[3]); 
  out[13] += 0.1936491673103708*f[1]*alphaDrag[50]+0.2165063509461096*(f[3]*alphaDrag[36]+f[2]*alphaDrag[35]+f[13]*alphaDrag[33]+f[8]*alphaDrag[31]+f[7]*alphaDrag[30])+0.1936491673103708*f[22]*alphaDrag[29]+0.2165063509461096*(f[0]*alphaDrag[29]+f[1]*alphaDrag[28]); 
  out[14] += 0.1936491673103708*f[2]*alphaDrag[51]+0.2165063509461096*(f[3]*alphaDrag[37]+f[1]*alphaDrag[35]+f[14]*alphaDrag[33]+f[9]*alphaDrag[31])+0.1936491673103708*f[23]*alphaDrag[30]+0.2165063509461096*(f[0]*alphaDrag[30]+f[7]*alphaDrag[29]+f[2]*alphaDrag[28]); 
  out[15] += 0.1936491673103708*f[3]*alphaDrag[52]+0.2165063509461096*(f[2]*alphaDrag[37]+f[1]*alphaDrag[36]+f[15]*alphaDrag[33])+0.1936491673103708*f[24]*alphaDrag[31]+0.2165063509461096*(f[0]*alphaDrag[31]+f[9]*alphaDrag[30]+f[8]*alphaDrag[29]+f[3]*alphaDrag[28]); 
  out[16] += 0.2165063509461096*(f[16]*alphaDrag[33]+f[12]*alphaDrag[31]+f[11]*alphaDrag[30]+f[10]*alphaDrag[29]+f[4]*alphaDrag[28]+alphaDrag[4]*f[16]+alphaDrag[3]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[5]); 
  out[17] += 0.1936491673103708*f[1]*alphaDrag[78]+0.2165063509461096*(f[3]*alphaDrag[64]+f[2]*alphaDrag[63]+f[17]*alphaDrag[62]+f[8]*alphaDrag[59]+f[7]*alphaDrag[58])+0.1936491673103708*f[22]*alphaDrag[57]+0.2165063509461096*(f[0]*alphaDrag[57]+f[1]*alphaDrag[56]); 
  out[18] += 0.1936491673103708*f[2]*alphaDrag[79]+0.2165063509461096*(f[3]*alphaDrag[65]+f[1]*alphaDrag[63]+f[18]*alphaDrag[62]+f[9]*alphaDrag[59])+0.1936491673103708*f[23]*alphaDrag[58]+0.2165063509461096*(f[0]*alphaDrag[58]+f[7]*alphaDrag[57]+f[2]*alphaDrag[56]); 
  out[19] += 0.1936491673103708*f[3]*alphaDrag[80]+0.2165063509461096*(f[2]*alphaDrag[65]+f[1]*alphaDrag[64]+f[19]*alphaDrag[62])+0.1936491673103708*f[24]*alphaDrag[59]+0.2165063509461096*(f[0]*alphaDrag[59]+f[9]*alphaDrag[58]+f[8]*alphaDrag[57]+f[3]*alphaDrag[56]); 
  out[20] += 0.2165063509461096*(f[20]*alphaDrag[62]+f[12]*alphaDrag[59]+f[11]*alphaDrag[58]+f[10]*alphaDrag[57]+f[4]*alphaDrag[56]+alphaDrag[4]*f[20]+alphaDrag[3]*f[19]+alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[0]*f[6]); 
  out[21] += 0.2165063509461096*(f[21]*alphaDrag[62]+f[15]*alphaDrag[59]+f[14]*alphaDrag[58]+f[13]*alphaDrag[57]+f[5]*alphaDrag[56]+f[21]*alphaDrag[33]+f[19]*alphaDrag[31]+f[18]*alphaDrag[30]+f[17]*alphaDrag[29]+f[6]*alphaDrag[28]); 
  out[25] += 2.371708245126284*(facDiff[9]*f[24]+facDiff[8]*f[23]+facDiff[7]*f[22]+facDiff[6]*f[9]+facDiff[5]*f[8]+facDiff[4]*f[7]+f[3]*facDiff[3]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+0.4330127018922193*alphaDrag[4]*f[25]+0.4841229182759271*(alphaDrag[3]*f[12]+alphaDrag[2]*f[11]+alphaDrag[1]*f[10]+f[0]*alphaDrag[4]+alphaDrag[0]*f[4]); 
  out[26] += 2.371708245126284*(facDiff[9]*f[24]+facDiff[8]*f[23]+facDiff[7]*f[22]+facDiff[6]*f[9]+facDiff[5]*f[8]+facDiff[4]*f[7]+f[3]*facDiff[3]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+0.4330127018922193*f[26]*alphaDrag[33]+0.4841229182759271*(f[0]*alphaDrag[33]+f[15]*alphaDrag[31]+f[14]*alphaDrag[30]+f[13]*alphaDrag[29]+f[5]*alphaDrag[28]); 
  out[27] += 2.371708245126284*(facDiff[9]*f[24]+facDiff[8]*f[23]+facDiff[7]*f[22]+facDiff[6]*f[9]+facDiff[5]*f[8]+facDiff[4]*f[7]+f[3]*facDiff[3]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvzSq4+0.4330127018922193*f[27]*alphaDrag[62]+0.4841229182759271*(f[0]*alphaDrag[62]+f[19]*alphaDrag[59]+f[18]*alphaDrag[58]+f[17]*alphaDrag[57]+f[6]*alphaDrag[56]); 

  return std::abs(0.0625*alphaDrag[0]-0.06987712429686843*(alphaDrag[24]+alphaDrag[23]+alphaDrag[22]))+std::abs(0.0625*alphaDrag[28]-0.06987712429686843*(alphaDrag[52]+alphaDrag[51]+alphaDrag[50]))+std::abs(0.0625*alphaDrag[56]-0.06987712429686843*(alphaDrag[80]+alphaDrag[79]+alphaDrag[78]))+std::abs((0.6363961030678926*facDiff[0]-0.711512473537885*(facDiff[9]+facDiff[8]+facDiff[7]))*rdvxSq4)+std::abs((0.6363961030678926*facDiff[0]-0.711512473537885*(facDiff[9]+facDiff[8]+facDiff[7]))*rdvySq4)+std::abs((0.6363961030678926*facDiff[0]-0.711512473537885*(facDiff[9]+facDiff[8]+facDiff[7]))*rdvzSq4); 

} 
