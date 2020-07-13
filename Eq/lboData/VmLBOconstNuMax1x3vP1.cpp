#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x3vMaxP1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
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

  double alphaDrag[15]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-4.0*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -1.154700538379252*dxv[1]*nuSum*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[5] = (2.828427124746191*nuUSum[2]-4.0*w[2]*nuSum)*rdvy2; 
  alphaDrag[6] = 2.828427124746191*nuUSum[3]*rdvy2; 
  alphaDrag[8] = -1.154700538379252*dxv[2]*nuSum*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[10] = (2.828427124746191*nuUSum[4]-4.0*w[3]*nuSum)*rdvz2; 
  alphaDrag[11] = 2.828427124746191*nuUSum[5]*rdvz2; 
  alphaDrag[14] = -1.154700538379252*dxv[3]*nuSum*rdvz2; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[3]*alphaDrag[8]+f[1]*alphaDrag[6]+f[0]*alphaDrag[5]); 
  out[4] += 0.4330127018922193*(f[4]*alphaDrag[14]+f[1]*alphaDrag[11]+f[0]*alphaDrag[10]); 

  return std::abs(0.125*alphaDrag[0])+std::abs(0.125*alphaDrag[5])+std::abs(0.125*alphaDrag[10])+std::abs(0.9428090415820636*nuVtSqSum[0]*rdvxSq4)+std::abs(0.9428090415820636*nuVtSqSum[0]*rdvySq4)+std::abs(0.9428090415820636*nuVtSqSum[0]*rdvzSq4); 

} 
