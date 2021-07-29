#include <GkLBOModDecl.h> 
void GkLBOconstNuCFLfreqMin1x2vTensorP2(const double m_, const double *w, const double *dxv, const double *Lv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // m_:        species mass. 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // Lv[2]:    domain length in velocity space.
  // BmagInv[3]:reciprocal of the magnetic field magnitude, 1/B. 
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

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs(kv[0]*((-1.0*w[1]*nuSum)-0.7905694150420947*nuUSum[2]+0.7071067811865475*nuUSum[0]))+std::fabs(-2.0*kv[1]*w[2]*nuSum)); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], 1.8*(std::fabs(0.7071067811865475*kvSq[0]*nuVtSqSum[0]-0.7905694150420947*kvSq[0]*nuVtSqSum[2])+std::fabs(0.2857142857142857*kvSq[1]*BmagInv[2]*nuVtSqSum[2]*w[2]*m_-1.118033988749895*BmagInv[0]*kvSq[1]*nuVtSqSum[2]*w[2]*m_-1.118033988749895*nuVtSqSum[0]*kvSq[1]*BmagInv[2]*w[2]*m_+BmagInv[0]*nuVtSqSum[0]*kvSq[1]*w[2]*m_))); 

}
