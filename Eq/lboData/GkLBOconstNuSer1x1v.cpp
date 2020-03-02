#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f[5]:      Input distribution function. 
  // out[5]:    Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 

  double alphaVpar[5]; 
  alphaVpar[0] = rdv2[0]*(1.414213562373095*nuUSum[0]-2.0*w[1]*nuSum); 
  alphaVpar[1] = 1.414213562373095*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = -1.154700538379252*nuSum; 

  double facVpar[2]; 
  facVpar[0] = nuVtSqSum[0]*rdvSq4[0]; 
  facVpar[1] = rdvSq4[0]*nuVtSqSum[1]; 

  out[2] += 0.8660254037844386*(alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaVpar[2]*f[3]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[4] += 1.732050807568877*alphaVpar[2]*f[4]+1.936491673103709*(alphaVpar[1]*f[3]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2])+4.743416490252569*(f[1]*facVpar[1]+f[0]*facVpar[0]); 

  return std::abs(0.5*alphaVpar[0]) + 0.9428090415820625*nuVtSqSum[0]*rdvSq4[0]; 

} 
double GkLBOconstNuVol1x1vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f[8]:      Input distribution function. 
  // out[8]:    Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 

  double alphaVpar[8]; 
  alphaVpar[0] = rdv2[0]*(1.414213562373095*nuUSum[0]-2.0*w[1]*nuSum); 
  alphaVpar[1] = 1.414213562373095*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = -1.154700538379252*nuSum; 
  alphaVpar[4] = 1.414213562373095*rdv2[0]*nuUSum[2]; 

  double facVpar[3]; 
  facVpar[0] = nuVtSqSum[0]*rdvSq4[0]; 
  facVpar[1] = rdvSq4[0]*nuVtSqSum[1]; 
  facVpar[2] = rdvSq4[0]*nuVtSqSum[2]; 

  out[2] += 0.8660254037844386*(alphaVpar[4]*f[4]+alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.7745966692414833*(alphaVpar[1]*f[4]+f[1]*alphaVpar[4])+0.8660254037844386*(alphaVpar[2]*f[3]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 
  out[5] += 1.936491673103709*alphaVpar[4]*f[6]+1.732050807568877*alphaVpar[2]*f[5]+4.743416490252569*facVpar[2]*f[4]+1.936491673103709*(alphaVpar[1]*f[3]+alphaVpar[0]*f[2]+f[0]*alphaVpar[2])+4.743416490252569*(f[1]*facVpar[1]+f[0]*facVpar[0]); 
  out[6] += 0.8660254037844386*alphaVpar[2]*f[6]+0.5532833351724881*alphaVpar[4]*f[4]+0.8660254037844386*(alphaVpar[0]*f[4]+f[0]*alphaVpar[4])+0.7745966692414833*alphaVpar[1]*f[1]; 
  out[7] += 1.732050807568877*(alphaVpar[2]*f[7]+alphaVpar[1]*f[6])+4.242640687119286*facVpar[1]*f[4]+f[3]*(1.732050807568877*alphaVpar[4]+1.936491673103709*alphaVpar[0])+4.242640687119286*f[1]*facVpar[2]+1.936491673103709*(alphaVpar[1]*f[2]+f[1]*alphaVpar[2])+4.743416490252569*(f[0]*facVpar[1]+facVpar[0]*f[1]); 

  return std::abs(0.5*alphaVpar[0]-0.5590169943749475*alphaVpar[4]) + 1.272792206135784*facVpar[0]-1.42302494707577*facVpar[2]; 

} 
