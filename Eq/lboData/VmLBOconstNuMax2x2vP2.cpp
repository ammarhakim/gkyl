#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x2vMaxP2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
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

  double alphaDrag[30]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-4.0*w[2]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.0*nuUSum[2]*rdvx2; 
  alphaDrag[3] = -1.154700538379252*dxv[2]*nuSum*rdvx2; 
  alphaDrag[5] = 2.0*nuUSum[3]*rdvx2; 
  alphaDrag[11] = 2.0*nuUSum[4]*rdvx2; 
  alphaDrag[12] = 2.0*nuUSum[5]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[15] = (2.0*nuUSum[6]-4.0*w[3]*nuSum)*rdvy2; 
  alphaDrag[16] = 2.0*nuUSum[7]*rdvy2; 
  alphaDrag[17] = 2.0*nuUSum[8]*rdvy2; 
  alphaDrag[19] = -1.154700538379252*dxv[3]*nuSum*rdvy2; 
  alphaDrag[20] = 2.0*nuUSum[9]*rdvy2; 
  alphaDrag[26] = 2.0*nuUSum[10]*rdvy2; 
  alphaDrag[27] = 2.0*nuUSum[11]*rdvy2; 

  double facDiff[6]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 
  facDiff[3] = nuVtSqSum[3]; 
  facDiff[4] = nuVtSqSum[4]; 
  facDiff[5] = nuVtSqSum[5]; 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.4330127018922193*(alphaDrag[12]*f[12]+alphaDrag[11]*f[11]+alphaDrag[5]*f[5]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.4330127018922193*(f[12]*alphaDrag[27]+f[11]*alphaDrag[26]+f[5]*alphaDrag[20]+f[4]*alphaDrag[19]+f[2]*alphaDrag[17]+f[1]*alphaDrag[16]+f[0]*alphaDrag[15]); 
  out[6] += 0.3872983346207416*(alphaDrag[1]*f[11]+f[1]*alphaDrag[11])+0.4330127018922193*(alphaDrag[3]*f[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[7] += 0.3872983346207416*(alphaDrag[2]*f[12]+f[2]*alphaDrag[12])+0.4330127018922193*(alphaDrag[3]*f[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[8] += 0.3872983346207416*f[1]*alphaDrag[26]+0.4330127018922193*(f[2]*alphaDrag[20]+f[8]*alphaDrag[19]+f[5]*alphaDrag[17])+0.3872983346207416*f[11]*alphaDrag[16]+0.4330127018922193*(f[0]*alphaDrag[16]+f[1]*alphaDrag[15]); 
  out[9] += 0.3872983346207416*f[2]*alphaDrag[27]+0.4330127018922193*(f[1]*alphaDrag[20]+f[9]*alphaDrag[19])+0.3872983346207416*f[12]*alphaDrag[17]+0.4330127018922193*(f[0]*alphaDrag[17]+f[5]*alphaDrag[16]+f[2]*alphaDrag[15]); 
  out[10] += 0.4330127018922193*(f[10]*alphaDrag[19]+f[7]*alphaDrag[17]+f[6]*alphaDrag[16]+f[3]*alphaDrag[15]+alphaDrag[3]*f[10]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[13] += 3.354101966249685*(facDiff[5]*f[12]+facDiff[4]*f[11]+facDiff[3]*f[5]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+0.8660254037844386*alphaDrag[3]*f[13]+0.9682458365518543*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[14] += 3.354101966249685*(facDiff[5]*f[12]+facDiff[4]*f[11]+facDiff[3]*f[5]+f[2]*facDiff[2]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+0.8660254037844386*f[14]*alphaDrag[19]+0.9682458365518543*(f[0]*alphaDrag[19]+f[9]*alphaDrag[17]+f[8]*alphaDrag[16]+f[4]*alphaDrag[15]); 

  return std::abs(0.125*alphaDrag[0]-0.1397542485937369*(alphaDrag[12]+alphaDrag[11]))+std::abs(0.125*alphaDrag[15]-0.1397542485937369*(alphaDrag[27]+alphaDrag[26]))+std::abs((0.9*facDiff[0]-1.006230589874905*(facDiff[5]+facDiff[4]))*rdvxSq4)+std::abs((0.9*facDiff[0]-1.006230589874905*(facDiff[5]+facDiff[4]))*rdvySq4); 

} 
