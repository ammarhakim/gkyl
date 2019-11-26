#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x2vSerP1(const double m_, const double *w, const double *dxv, double *cflRateByDir, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[3]:            Cell-center coordinates. 
  // dxv[3]:          Cell spacing. 
  // cflRatebyDir[3]: CFL rate in each direction. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // nuUSum:          sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:       sum of thermal speeds squared time their respective collisionalities. 
  // f[8]:            Input distribution function. 
  // out[8]:          Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 
  rdv2[1]   = 2.0/dxv[2]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 
  rdv2[1]   = rdv2[1]; 

  double alphaVpar[8]; 
  alphaVpar[0] = rdv2[0]*(2.0*nuUSum[0]-2.828427124746191*w[1]*nuSum); 
  alphaVpar[1] = 2.0*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = -1.632993161855453*nuSum; 

  double alphaMu[8]; 
  alphaMu[0] = (-5.656854249492382*rdv2[1]*w[2]*nuSum)+2.828427124746191*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_+2.828427124746191*BmagInv[0]*nuVtSqSum[0]*rdv2[1]*m_; 
  alphaMu[1] = 2.828427124746191*BmagInv[0]*nuVtSqSum[1]*rdv2[1]*m_+2.828427124746191*nuVtSqSum[0]*BmagInv[1]*rdv2[1]*m_; 
  alphaMu[3] = -3.265986323710906*nuSum; 

  cflRateByDir[0] = 0.; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 0.125*((-0.6123724356957944*alphaVpar[2])-0.3535533905932737*alphaVpar[1]+0.3535533905932737*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.125*(0.3535533905932737*(alphaVpar[1]+alphaVpar[0])-0.6123724356957944*alphaVpar[2]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.125*((-0.6123724356957944*alphaVpar[2])-0.3535533905932737*alphaVpar[1]+0.3535533905932737*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.125*(0.3535533905932737*(alphaVpar[1]+alphaVpar[0])-0.6123724356957944*alphaVpar[2]); 
  cflRateByDir[1] += std::abs(alphaL); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 0.125*(0.6123724356957944*alphaVpar[2]-0.3535533905932737*alphaVpar[1]+0.3535533905932737*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.125*(0.6123724356957944*alphaVpar[2]+0.3535533905932737*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.125*(0.6123724356957944*alphaVpar[2]-0.3535533905932737*alphaVpar[1]+0.3535533905932737*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.125*(0.6123724356957944*alphaVpar[2]+0.3535533905932737*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 0.125*((-0.6123724356957944*alphaMu[3])-0.3535533905932737*alphaMu[1]+0.3535533905932737*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.125*(0.3535533905932737*(alphaMu[1]+alphaMu[0])-0.6123724356957944*alphaMu[3]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.125*((-0.6123724356957944*alphaMu[3])-0.3535533905932737*alphaMu[1]+0.3535533905932737*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  alphaL = 0.125*(0.3535533905932737*(alphaMu[1]+alphaMu[0])-0.6123724356957944*alphaMu[3]); 
  cflRateByDir[2] += std::abs(alphaL); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 0.125*(0.6123724356957944*alphaMu[3]-0.3535533905932737*alphaMu[1]+0.3535533905932737*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.125*(0.6123724356957944*alphaMu[3]+0.3535533905932737*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.125*(0.6123724356957944*alphaMu[3]-0.3535533905932737*alphaMu[1]+0.3535533905932737*alphaMu[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  alphaR = 0.125*(0.6123724356957944*alphaMu[3]+0.3535533905932737*(alphaMu[1]+alphaMu[0])); 
  cflRateByDir[2] += std::abs(alphaR); 

  out[2] += 0.6123724356957944*(alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.6123724356957944*(alphaMu[3]*f[3]+alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphaVpar[2]*f[4]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 0.6123724356957944*(alphaMu[3]*f[5]+alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[6] += 0.6123724356957944*((alphaMu[3]+alphaVpar[2])*f[6]+alphaVpar[1]*f[5]+alphaMu[1]*f[4]+alphaVpar[0]*f[3]+alphaMu[0]*f[2]); 
  out[7] += 0.6123724356957944*((alphaMu[3]+alphaVpar[2])*f[7]+alphaVpar[0]*f[5]+alphaMu[0]*f[4]+alphaVpar[1]*f[3]+alphaMu[1]*f[2]); 

  return std::abs(0.3535533905932737*alphaVpar[0]) + 2.0*rdv2[1]*w[2]*nuSum+1.333333333333333*(BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*rdvSq4[1]*w[2]*m_+0.9428090415820625*nuVtSqSum[0]*rdvSq4[0]; 

} 
