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
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-8.0*w[3]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.828427124746191*nuUSum[2]*rdvx2; 
  alphaDrag[3] = 2.828427124746191*nuUSum[3]*rdvx2; 
  alphaDrag[4] = -2.309401076758503*dxv[3]*nuSum*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[7] = (2.828427124746191*nuUSum[4]-8.0*w[4]*nuSum)*rdvy2; 
  alphaDrag[8] = 2.828427124746191*nuUSum[5]*rdvy2; 
  alphaDrag[9] = 2.828427124746191*nuUSum[6]*rdvy2; 
  alphaDrag[10] = 2.828427124746191*nuUSum[7]*rdvy2; 
  alphaDrag[12] = -2.309401076758503*dxv[4]*nuSum*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[14] = (2.828427124746191*nuUSum[8]-8.0*w[5]*nuSum)*rdvz2; 
  alphaDrag[15] = 2.828427124746191*nuUSum[9]*rdvz2; 
  alphaDrag[16] = 2.828427124746191*nuUSum[10]*rdvz2; 
  alphaDrag[17] = 2.828427124746191*nuUSum[11]*rdvz2; 
  alphaDrag[20] = -2.309401076758503*dxv[5]*nuSum*rdvz2; 

  // Put together updates due to drag and diffusion terms.
  out[4] += 0.2165063509461096*(alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[5]*alphaDrag[12]+f[3]*alphaDrag[10]+f[2]*alphaDrag[9]+f[1]*alphaDrag[8]+f[0]*alphaDrag[7]); 
  out[6] += 0.2165063509461096*(f[6]*alphaDrag[20]+f[3]*alphaDrag[17]+f[2]*alphaDrag[16]+f[1]*alphaDrag[15]+f[0]*alphaDrag[14]); 

  return std::abs(0.0625*alphaDrag[0])+std::abs(0.0625*alphaDrag[7])+std::abs(0.0625*alphaDrag[14])+std::abs(0.4714045207910317*nuVtSqSum[0]*rdvxSq4)+std::abs(0.4714045207910317*nuVtSqSum[0]*rdvySq4)+std::abs(0.4714045207910317*nuVtSqSum[0]*rdvzSq4); 

} 
