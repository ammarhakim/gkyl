#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, double *cflRateByDir, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double *out) 
{ 
  // w[2]:            Cell-center coordinates. 
  // dxv[2]:          Cell spacing. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // cflRatebyDir[1]: collisionalities added (self and cross species collisionalities). 
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

  cflRateByDir[0] = 0.; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points. 
  alphaL = 0.25*((-0.8660254037844386*alphaVpar[2])-0.5*alphaVpar[1]+0.5*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  alphaL = 0.25*(0.5*(alphaVpar[1]+alphaVpar[0])-0.8660254037844386*alphaVpar[2]); 
  cflRateByDir[1] += std::abs(alphaL); 
  // Evaluate alpha at right surface quadrature points. 
  alphaR = 0.25*(0.8660254037844386*alphaVpar[2]-0.5*alphaVpar[1]+0.5*alphaVpar[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  alphaR = 0.25*(0.8660254037844386*alphaVpar[2]+0.5*(alphaVpar[1]+alphaVpar[0])); 
  cflRateByDir[1] += std::abs(alphaR); 

  out[2] += 0.8660254037844386*(alphaVpar[2]*f[2]+alphaVpar[1]*f[1]+alphaVpar[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaVpar[2]*f[3]+alphaVpar[0]*f[1]+f[0]*alphaVpar[1]); 

  return std::abs(0.5*alphaVpar[0]) + 0.9428090415820625*nuVtSqSum[0]*rdvSq4[0]; 

} 
