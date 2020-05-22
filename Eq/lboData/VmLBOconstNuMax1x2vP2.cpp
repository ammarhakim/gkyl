#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vMaxP2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 

  double alphaDrag[20]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*nuSum*rdvx2; 
  alphaDrag[7] = 2.0*nuUSum[2]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[10] = (2.0*nuUSum[3]-2.828427124746191*w[2]*nuSum)*rdvy2; 
  alphaDrag[11] = 2.0*nuUSum[4]*rdvy2; 
  alphaDrag[13] = -0.8164965809277261*dxv[2]*nuSum*rdvy2; 
  alphaDrag[17] = 2.0*nuUSum[5]*rdvy2; 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[7]*f[7]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alphaDrag[17]+f[3]*alphaDrag[13]+f[1]*alphaDrag[11]+f[0]*alphaDrag[10]); 
  out[4] += 0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.5477225575051661*f[1]*alphaDrag[17]+0.6123724356957944*f[5]*alphaDrag[13]+0.5477225575051661*f[7]*alphaDrag[11]+0.6123724356957944*(f[0]*alphaDrag[11]+f[1]*alphaDrag[10]); 
  out[6] += 0.6123724356957944*(f[6]*alphaDrag[13]+f[4]*alphaDrag[11]+f[2]*alphaDrag[10]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+1.224744871391589*f[9]*alphaDrag[13]+1.369306393762915*(f[0]*alphaDrag[13]+f[5]*alphaDrag[11]+f[3]*alphaDrag[10]); 

  return std::abs(0.1767766952966368*alphaDrag[0]-0.1976423537605236*alphaDrag[7])+std::abs(0.1767766952966368*alphaDrag[10]-0.1976423537605236*alphaDrag[17])+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4)+std::abs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4); 

} 
