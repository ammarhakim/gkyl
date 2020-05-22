#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x3vMaxP2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
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

  double alphaDrag[45]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-4.0*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -1.154700538379252*dxv[1]*nuSum*rdvx2; 
  alphaDrag[11] = 2.828427124746191*nuUSum[2]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[15] = (2.828427124746191*nuUSum[3]-4.0*w[2]*nuSum)*rdvy2; 
  alphaDrag[16] = 2.828427124746191*nuUSum[4]*rdvy2; 
  alphaDrag[18] = -1.154700538379252*dxv[2]*nuSum*rdvy2; 
  alphaDrag[26] = 2.828427124746191*nuUSum[5]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[30] = (2.828427124746191*nuUSum[6]-4.0*w[3]*nuSum)*rdvz2; 
  alphaDrag[31] = 2.828427124746191*nuUSum[7]*rdvz2; 
  alphaDrag[34] = -1.154700538379252*dxv[3]*nuSum*rdvz2; 
  alphaDrag[41] = 2.828427124746191*nuUSum[8]*rdvz2; 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alphaDrag[11]*f[11]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[11]*alphaDrag[26]+f[3]*alphaDrag[18]+f[1]*alphaDrag[16]+f[0]*alphaDrag[15]); 
  out[4] += 0.4330127018922193*(f[11]*alphaDrag[41]+f[4]*alphaDrag[34]+f[1]*alphaDrag[31]+f[0]*alphaDrag[30]); 
  out[5] += 0.3872983346207416*(alphaDrag[1]*f[11]+f[1]*alphaDrag[11])+0.4330127018922193*(alphaDrag[2]*f[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.3872983346207416*f[1]*alphaDrag[26]+0.4330127018922193*f[6]*alphaDrag[18]+0.3872983346207416*f[11]*alphaDrag[16]+0.4330127018922193*(f[0]*alphaDrag[16]+f[1]*alphaDrag[15]); 
  out[7] += 0.4330127018922193*(f[7]*alphaDrag[18]+f[5]*alphaDrag[16]+f[2]*alphaDrag[15]+alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]); 
  out[8] += 0.3872983346207416*f[1]*alphaDrag[41]+0.4330127018922193*f[8]*alphaDrag[34]+0.3872983346207416*f[11]*alphaDrag[31]+0.4330127018922193*(f[0]*alphaDrag[31]+f[1]*alphaDrag[30]); 
  out[9] += 0.4330127018922193*(f[9]*alphaDrag[34]+f[5]*alphaDrag[31]+f[2]*alphaDrag[30]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[10]*alphaDrag[34]+f[6]*alphaDrag[31]+f[3]*alphaDrag[30]+f[10]*alphaDrag[18]+f[8]*alphaDrag[16]+f[4]*alphaDrag[15]); 
  out[12] += 4.743416490252569*(facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+0.8660254037844386*alphaDrag[2]*f[12]+0.9682458365518543*(alphaDrag[1]*f[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[13] += 4.743416490252569*(facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+0.8660254037844386*f[13]*alphaDrag[18]+0.9682458365518543*(f[0]*alphaDrag[18]+f[6]*alphaDrag[16]+f[3]*alphaDrag[15]); 
  out[14] += 4.743416490252569*(facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvzSq4+0.8660254037844386*f[14]*alphaDrag[34]+0.9682458365518543*(f[0]*alphaDrag[34]+f[8]*alphaDrag[31]+f[4]*alphaDrag[30]); 

  return std::abs(0.125*alphaDrag[0]-0.1397542485937369*alphaDrag[11])+std::abs(0.125*alphaDrag[15]-0.1397542485937369*alphaDrag[26])+std::abs(0.125*alphaDrag[30]-0.1397542485937369*alphaDrag[41])+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4)+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4)+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvzSq4); 

} 
