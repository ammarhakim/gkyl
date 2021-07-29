#include <GkLBOModDecl.h> 
void GkLBOconstNuCFLfreqMin1x1vSerP2(const double m_, const double *w, const double *dxv, const double *Lv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // m_:        species mass. 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // Lv[1]:    domain length in velocity space.
  // BmagInv[3]:reciprocal of the magnetic field magnitude, 1/B. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // cflFreqDrag: min CFL frequency supported by drag terms.
  // cflFreqDiff: min CFL frequency supported by diffusion terms.

  double kv[1]; 
  double kvSq[1]; 
  kv[0]   = 6.283185307179586/Lv[0]; 
  kvSq[0] = kv[0]*kv[0]; 

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs(kv[0]*((-1.0*w[1]*nuSum)-0.7905694150420947*nuUSum[2]+0.7071067811865475*nuUSum[0]))); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], std::fabs((9*(0.7071067811865475*kvSq[0]*nuVtSqSum[0]-0.7905694150420947*kvSq[0]*nuVtSqSum[2]))/5)); 

}
