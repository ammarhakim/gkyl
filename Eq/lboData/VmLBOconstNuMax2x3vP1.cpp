#include <VmLBOModDecl.h> 
double VmLBOconstNuVol2x3vMaxP1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[5]:      Cell-center coordinates. 
  // dxv[5]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[2]; 
  const double rdvxSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvy2 = 2.0/dxv[3]; 
  const double rdvySq4 = 4.0/(dxv[3]*dxv[3]); 
  const double rdvz2 = 2.0/dxv[4]; 
  const double rdvzSq4 = 4.0/(dxv[4]*dxv[4]); 

  double alphaDrag[18]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-5.656854249492382*w[2]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.828427124746191*nuUSum[2]*rdvx2; 
  alphaDrag[3] = -1.632993161855453*dxv[2]*nuSum*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[6] = (2.828427124746191*nuUSum[3]-5.656854249492382*w[3]*nuSum)*rdvy2; 
  alphaDrag[7] = 2.828427124746191*nuUSum[4]*rdvy2; 
  alphaDrag[8] = 2.828427124746191*nuUSum[5]*rdvy2; 
  alphaDrag[10] = -1.632993161855453*dxv[3]*nuSum*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[12] = (2.828427124746191*nuUSum[6]-5.656854249492382*w[4]*nuSum)*rdvz2; 
  alphaDrag[13] = 2.828427124746191*nuUSum[7]*rdvz2; 
  alphaDrag[14] = 2.828427124746191*nuUSum[8]*rdvz2; 
  alphaDrag[17] = -1.632993161855453*dxv[4]*nuSum*rdvz2; 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.3061862178478971*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[4]*alphaDrag[10]+f[2]*alphaDrag[8]+f[1]*alphaDrag[7]+f[0]*alphaDrag[6]); 
  out[5] += 0.3061862178478971*(f[5]*alphaDrag[17]+f[2]*alphaDrag[14]+f[1]*alphaDrag[13]+f[0]*alphaDrag[12]); 

  return std::abs(0.0883883476483184*alphaDrag[0])+std::abs(0.0883883476483184*alphaDrag[6])+std::abs(0.0883883476483184*alphaDrag[12])+std::abs(0.6666666666666666*nuVtSqSum[0]*rdvxSq4)+std::abs(0.6666666666666666*nuVtSqSum[0]*rdvySq4)+std::abs(0.6666666666666666*nuVtSqSum[0]*rdvzSq4); 

} 
