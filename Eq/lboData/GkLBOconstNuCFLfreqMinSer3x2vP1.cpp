#include <GkLBOModDecl.h> 
void GkLBOconstNuCFLfreqMin3x2vSerP1(const double m_, const double *w, const double *dxv, const double *Lv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // m_:        species mass. 
  // w[5]:      Cell-center coordinates. 
  // dxv[5]:    Cell spacing. 
  // Lv[2]:    domain length in velocity space.
  // BmagInv[8]:reciprocal of the magnetic field magnitude, 1/B. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // cflFreqDrag: min CFL frequency supported by drag terms.
  // cflFreqDiff: min CFL frequency supported by diffusion terms.

  double kv[2]; 
  double kvSq[2]; 
  kv[0]   = 6.283185307179586/Lv[0]; 
  kvSq[0] = kv[0]*kv[0]; 
  kv[1]   = 6.283185307179586/Lv[1]; 
  kvSq[1] = kv[1]*kv[1]; 

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs(kv[0]*(0.3535533905932737*nuUSum[0]-1.0*w[3]*nuSum))+std::fabs(-2.0*kv[1]*w[4]*nuSum)); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], 1.8*(std::fabs(0.3535533905932737*kvSq[0]*nuVtSqSum[0])+std::fabs(0.25*kvSq[1]*w[4]*BmagInv[7]*nuVtSqSum[7]*m_+0.25*kvSq[1]*w[4]*BmagInv[6]*nuVtSqSum[6]*m_+0.25*kvSq[1]*w[4]*BmagInv[5]*nuVtSqSum[5]*m_+0.25*kvSq[1]*BmagInv[4]*nuVtSqSum[4]*w[4]*m_+0.25*kvSq[1]*BmagInv[3]*nuVtSqSum[3]*w[4]*m_+0.25*kvSq[1]*BmagInv[2]*nuVtSqSum[2]*w[4]*m_+0.25*BmagInv[1]*kvSq[1]*nuVtSqSum[1]*w[4]*m_+0.25*BmagInv[0]*nuVtSqSum[0]*kvSq[1]*w[4]*m_))); 

}
