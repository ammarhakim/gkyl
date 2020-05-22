#include <VmLBOModDecl.h> 
double VmLBOconstNuVol1x2vMaxP1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
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

  double alphaDrag[8]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*dxv[1]*nuSum*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[4] = (2.0*nuUSum[2]-2.828427124746191*w[2]*nuSum)*rdvy2; 
  alphaDrag[5] = 2.0*nuUSum[3]*rdvy2; 
  alphaDrag[7] = -0.8164965809277261*dxv[2]*nuSum*rdvy2; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[3]*alphaDrag[7]+f[1]*alphaDrag[5]+f[0]*alphaDrag[4]); 

  return std::abs(0.1767766952966368*alphaDrag[0])+std::abs(0.1767766952966368*alphaDrag[4])+std::abs(0.9428090415820636*nuVtSqSum[0]*rdvxSq4)+std::abs(0.9428090415820636*nuVtSqSum[0]*rdvySq4); 

} 
