#include <GkLBOModDecl.h> 
double GkLBOconstNuVol2x2vSerP1(const double m_, const double *w, const double *dxv, double *positivityWeightByDir, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[4]:            Cell-center coordinates. 
  // dxv[4]:          Cell spacing. 
  // cflRatebyDir[4]: CFL rate in each direction. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // nuUSum:          sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:       sum of thermal speeds squared time their respective collisionalities. 
  // f[16]:            Input distribution function. 
  // out[16]:          Incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[2]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 
  rdv2[1]   = 2.0/dxv[3]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 
  rdv2[1]   = rdv2[1]; 

  double alphaVpar[16]; 
  alphaVpar[0] = rdv2[0]*(2.0*nuUSum[0]-4.0*w[2]*nuSum); 
  alphaVpar[1] = 2.0*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = 2.0*rdv2[0]*nuUSum[2]; 
  alphaVpar[3] = -2.309401076758503*nuSum; 
  alphaVpar[5] = 2.0*rdv2[0]*nuUSum[3]; 

  double alphaMu[16]; 
  alphaMu[0] = (-8.0*rdv2[1]*w[3]*nuSum)+2.0*rdv2[1]*BmagInv[3]*nuVtSqSum[3]*m_+2.0*rdv2[1]*BmagInv[2]*nuVtSqSum[2]*m_+2.0*BmagInv[1]*nuVtSqSum[1]*rdv2[1]*m_+2.0*BmagInv[0]*nuVtSqSum[0]*rdv2[1]*m_; 
  alphaMu[1] = 2.0*rdv2[1]*BmagInv[2]*nuVtSqSum[3]*m_+2.0*rdv2[1]*nuVtSqSum[2]*BmagInv[3]*m_+2.0*BmagInv[0]*nuVtSqSum[1]*rdv2[1]*m_+2.0*nuVtSqSum[0]*BmagInv[1]*rdv2[1]*m_; 
  alphaMu[2] = 2.0*BmagInv[1]*rdv2[1]*nuVtSqSum[3]*m_+2.0*nuVtSqSum[1]*rdv2[1]*BmagInv[3]*m_+2.0*BmagInv[0]*rdv2[1]*nuVtSqSum[2]*m_+2.0*nuVtSqSum[0]*rdv2[1]*BmagInv[2]*m_; 
  alphaMu[4] = -4.618802153517007*nuSum; 
  alphaMu[5] = 2.0*BmagInv[0]*rdv2[1]*nuVtSqSum[3]*m_+2.0*nuVtSqSum[0]*rdv2[1]*BmagInv[3]*m_+2.0*BmagInv[1]*rdv2[1]*nuVtSqSum[2]*m_+2.0*nuVtSqSum[1]*rdv2[1]*BmagInv[2]*m_; 

  out[3] += 0.4330127018922193*(alphaVpar[5]*f[5]+alphaVpar[3]*f[3]+alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[4] += 0.4330127018922193*(alphaMu[5]*f[5]+alphaMu[4]*f[4]+alphaMu[2]*f[2]+alphaMu[1]*f[1]+alphaMu[0]*f[0]); 
  out[6] += 0.4330127018922193*(alphaVpar[3]*f[6]+alphaVpar[2]*f[5]+f[2]*alphaVpar[5]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[7] += 0.4330127018922193*(alphaVpar[3]*f[7]+alphaVpar[1]*f[5]+f[1]*alphaVpar[5]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2]); 
  out[8] += 0.4330127018922193*(alphaMu[4]*f[8]+alphaMu[2]*f[5]+f[2]*alphaMu[5]+alphaMu[0]*f[1]+f[0]*alphaMu[1]); 
  out[9] += 0.4330127018922193*(alphaMu[4]*f[9]+alphaMu[1]*f[5]+f[1]*alphaMu[5]+alphaMu[0]*f[2]+f[0]*alphaMu[2]); 
  out[10] += 0.4330127018922193*(alphaVpar[5]*f[12]+alphaMu[5]*f[11]+(alphaMu[4]+alphaVpar[3])*f[10]+alphaVpar[2]*f[9]+alphaVpar[1]*f[8]+alphaMu[2]*f[7]+alphaMu[1]*f[6]+alphaVpar[0]*f[4]+alphaMu[0]*f[3]); 
  out[11] += 0.4330127018922193*(alphaVpar[3]*f[11]+alphaVpar[0]*f[5]+f[0]*alphaVpar[5]+alphaVpar[1]*f[2]+f[1]*alphaVpar[2]); 
  out[12] += 0.4330127018922193*(alphaMu[4]*f[12]+alphaMu[0]*f[5]+f[0]*alphaMu[5]+alphaMu[1]*f[2]+f[1]*alphaMu[2]); 
  out[13] += 0.4330127018922193*((alphaMu[4]+alphaVpar[3])*f[13]+alphaVpar[2]*f[12]+alphaMu[2]*f[11]+alphaVpar[5]*f[9]+alphaVpar[0]*f[8]+alphaMu[5]*f[7]+alphaMu[0]*f[6]+alphaVpar[1]*f[4]+alphaMu[1]*f[3]); 
  out[14] += 0.4330127018922193*((alphaMu[4]+alphaVpar[3])*f[14]+alphaVpar[1]*f[12]+alphaMu[1]*f[11]+alphaVpar[0]*f[9]+alphaVpar[5]*f[8]+alphaMu[0]*f[7]+alphaMu[5]*f[6]+alphaVpar[2]*f[4]+alphaMu[2]*f[3]); 
  out[15] += 0.4330127018922193*((alphaMu[4]+alphaVpar[3])*f[15]+alphaVpar[0]*f[12]+alphaMu[0]*f[11]+alphaVpar[1]*f[9]+alphaVpar[2]*f[8]+alphaMu[1]*f[7]+alphaMu[2]*f[6]+f[4]*alphaVpar[5]+f[3]*alphaMu[5]); 

  positivityWeightByDir[1] = std::abs(0.25*alphaVpar[0]) + 0.6666666666666666*nuVtSqSum[0]*rdvSq4[0]; 
  positivityWeightByDir[2] = 2.0*rdv2[1]*w[3]*nuSum+rdvSq4[1]*(0.6666666666666666*BmagInv[3]*nuVtSqSum[3]+0.6666666666666666*BmagInv[2]*nuVtSqSum[2]+0.6666666666666666*BmagInv[1]*nuVtSqSum[1]+0.6666666666666666*BmagInv[0]*nuVtSqSum[0])*w[3]*m_; 
  return std::abs(0.25*alphaVpar[0]) + 2.0*rdv2[1]*w[3]*nuSum+rdvSq4[1]*(0.6666666666666666*BmagInv[3]*nuVtSqSum[3]+0.6666666666666666*BmagInv[2]*nuVtSqSum[2]+0.6666666666666666*BmagInv[1]*nuVtSqSum[1]+0.6666666666666666*BmagInv[0]*nuVtSqSum[0])*w[3]*m_+0.6666666666666666*nuVtSqSum[0]*rdvSq4[0]; 

} 
