#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, double *positivityWeightByDir, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[2]:            Cell-center coordinates. 
  // dxv[2]:          Cell spacing. 
  // cflRatebyDir[2]: CFL rate in each direction. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // nuUSum:          sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum:       sum of thermal speeds squared time their respective collisionalities. 
  // f[4]:            Input distribution function. 
  // out[4]:          Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[0]   = rdv2[0]; 

  double alphaVpar[4]; 
  alphaVpar[0] = rdv2[0]*(1.414213562373095*nuUSum[0]-2.0*w[1]*nuSum); 
  alphaVpar[1] = 1.414213562373095*rdv2[0]*nuUSum[1]; 
  alphaVpar[2] = -1.154700538379252*nuSum; 

  out[2] += 0.8660254037844386*(alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaVpar[2]*f[3]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 

  positivityWeightByDir[1] = std::abs(0.5*alphaVpar[0]) + 0.9428090415820625*nuVtSqSum[0]*rdvSq4[0]; 
  return std::abs(0.5*alphaVpar[0]) + 0.9428090415820625*nuVtSqSum[0]*rdvSq4[0]; 

} 
