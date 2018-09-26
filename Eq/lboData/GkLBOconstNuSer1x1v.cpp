#include <GkLBOModDecl.h> 
double GkLBOconstNuVol1x1vSerP1(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 

  double alpha0[4]; 
  alpha0[0] = 1.414213562373095*rdv2[0]*nuU[0]-2.0*rdv2[0]*w[1]*nu; 
  alpha0[1] = 1.414213562373095*rdv2[0]*nuU[1]; 
  alpha0[2] = -1.154700538379252*nu; 

  double alpha1[4]; 
  alpha1[0] = 1.414213562373095*nuVtSq[0]*rdvSq4[0]; 
  alpha1[1] = 1.414213562373095*rdvSq4[0]*nuVtSq[1]; 

  out[2] += 0.8660254037844386*(alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[3] += 0.8660254037844386*(alpha0[2]*f[3]+alpha0[0]*f[1]+f[0]*alpha0[1]); 

  const double alpha0Mid = 0.5*alpha0[0]; 
  const double alpha1Mid = 0.6666666666666666*alpha1[0]; 
  return std::abs(alpha0Mid) + alpha1Mid; 

} 
double GkLBOconstNuVol1x1vSerP2(const double m_, const double *w, const double *dxv, const double *BmagInv, const double nu, const double *nuU, const double *nuVtSq, const double *f, double *out) 
{ 
  // w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu: diffusion coefficient (collisionality). nuU: bulk velocity times nu. nuVtSq: thermal speed squared times nu. f: Input distribution function. out: Incremented output 
  double rdv2[1]; 
  double rdvSq4[1]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 

  double alpha0[8]; 
  alpha0[0] = 1.414213562373095*rdv2[0]*nuU[0]-2.0*rdv2[0]*w[1]*nu; 
  alpha0[1] = 1.414213562373095*rdv2[0]*nuU[1]; 
  alpha0[2] = -1.154700538379252*nu; 
  alpha0[4] = 1.414213562373095*rdv2[0]*nuU[2]; 

  double alpha1[8]; 
  alpha1[0] = 1.414213562373095*nuVtSq[0]*rdvSq4[0]; 
  alpha1[1] = 1.414213562373095*rdvSq4[0]*nuVtSq[1]; 
  alpha1[4] = 1.414213562373095*rdvSq4[0]*nuVtSq[2]; 

  out[2] += 0.8660254037844386*(alpha0[4]*f[4]+alpha0[2]*f[2]+alpha0[1]*f[1]+alpha0[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha0[1]*f[4]+f[1]*alpha0[4])+0.8660254037844386*(alpha0[2]*f[3]+alpha0[0]*f[1]+f[0]*alpha0[1]); 
  out[5] += 1.936491673103709*alpha0[4]*f[6]+1.732050807568877*alpha0[2]*f[5]+3.354101966249685*alpha1[4]*f[4]+1.936491673103709*(alpha0[1]*f[3]+alpha0[0]*f[2]+f[0]*alpha0[2])+3.354101966249685*(alpha1[1]*f[1]+alpha1[0]*f[0]); 
  out[6] += 0.8660254037844386*alpha0[2]*f[6]+0.5532833351724881*alpha0[4]*f[4]+0.8660254037844386*(alpha0[0]*f[4]+f[0]*alpha0[4])+0.7745966692414833*alpha0[1]*f[1]; 
  out[7] += 1.732050807568877*(alpha0[2]*f[7]+alpha0[1]*f[6])+3.0*(alpha1[1]*f[4]+f[1]*alpha1[4])+1.732050807568877*f[3]*alpha0[4]+1.936491673103709*(alpha0[0]*f[3]+alpha0[1]*f[2]+f[1]*alpha0[2])+3.354101966249685*(alpha1[0]*f[1]+f[0]*alpha1[1]); 

  const double alpha0Mid = 0.5*alpha0[0]-0.5590169943749475*alpha0[4]; 
  const double alpha1Mid = (9*(0.5*alpha1[0]-0.5590169943749475*alpha1[4]))/5; 
  return std::abs(alpha0Mid) + alpha1Mid; 

} 
