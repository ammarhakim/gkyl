#include <GkLBOModDecl.h> 
void GkLBOconstNuCFLfreqMin3x2vSerP2(const double m_, const double *w, const double *dxv, const double *Lv, const double *BmagInv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, double *cflFreqDrag, double *cflFreqDiff) 
{ 
  // m_:        species mass. 
  // w[5]:      Cell-center coordinates. 
  // dxv[5]:    Cell spacing. 
  // Lv[2]:    domain length in velocity space.
  // BmagInv[20]:reciprocal of the magnetic field magnitude, 1/B. 
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

  cflFreqDrag[0] = std::fmax(cflFreqDrag[0], std::fabs(kv[0]*((-1.0*w[3]*nuSum)-0.3952847075210473*nuUSum[9]-0.3952847075210473*nuUSum[8]-0.3952847075210473*nuUSum[7]+0.3535533905932737*nuUSum[0]))+std::fabs(-2.0*kv[1]*w[4]*nuSum)); 
  cflFreqDiff[0] = std::fmax(cflFreqDiff[0], 1.8*(std::fabs((-0.3952847075210473*kvSq[0]*nuVtSqSum[9])-0.3952847075210473*kvSq[0]*nuVtSqSum[8]-0.3952847075210473*kvSq[0]*nuVtSqSum[7]+0.3535533905932737*kvSq[0]*nuVtSqSum[0])+std::fabs((-0.4285714285714285*kvSq[1]*w[4]*BmagInv[19]*nuVtSqSum[19]*m_)-0.2795084971874737*kvSq[1]*BmagInv[4]*w[4]*nuVtSqSum[19]*m_-0.2795084971874737*kvSq[1]*nuVtSqSum[4]*w[4]*BmagInv[19]*m_-0.4285714285714285*kvSq[1]*w[4]*BmagInv[18]*nuVtSqSum[18]*m_-0.2795084971874737*kvSq[1]*w[4]*BmagInv[5]*nuVtSqSum[18]*m_-0.2795084971874737*kvSq[1]*w[4]*nuVtSqSum[5]*BmagInv[18]*m_-0.4285714285714285*kvSq[1]*w[4]*BmagInv[17]*nuVtSqSum[17]*m_-0.2795084971874737*kvSq[1]*w[4]*BmagInv[6]*nuVtSqSum[17]*m_-0.2795084971874737*kvSq[1]*w[4]*nuVtSqSum[6]*BmagInv[17]*m_-0.1785714285714286*kvSq[1]*w[4]*BmagInv[16]*nuVtSqSum[16]*m_-0.2795084971874738*kvSq[1]*BmagInv[2]*w[4]*nuVtSqSum[16]*m_-0.2795084971874738*kvSq[1]*nuVtSqSum[2]*w[4]*BmagInv[16]*m_-0.1785714285714286*kvSq[1]*w[4]*BmagInv[15]*nuVtSqSum[15]*m_-0.2795084971874738*BmagInv[1]*kvSq[1]*w[4]*nuVtSqSum[15]*m_-0.2795084971874738*kvSq[1]*nuVtSqSum[1]*w[4]*BmagInv[15]*m_-0.1785714285714286*kvSq[1]*w[4]*BmagInv[14]*nuVtSqSum[14]*m_-0.2795084971874738*kvSq[1]*BmagInv[3]*w[4]*nuVtSqSum[14]*m_-0.2795084971874738*kvSq[1]*nuVtSqSum[3]*w[4]*BmagInv[14]*m_-0.1785714285714286*kvSq[1]*w[4]*BmagInv[13]*nuVtSqSum[13]*m_-0.2795084971874738*kvSq[1]*BmagInv[3]*w[4]*nuVtSqSum[13]*m_-0.2795084971874738*kvSq[1]*nuVtSqSum[3]*w[4]*BmagInv[13]*m_-0.1785714285714286*kvSq[1]*w[4]*BmagInv[12]*nuVtSqSum[12]*m_-0.2795084971874738*BmagInv[1]*kvSq[1]*w[4]*nuVtSqSum[12]*m_-0.2795084971874738*kvSq[1]*nuVtSqSum[1]*w[4]*BmagInv[12]*m_-0.1785714285714286*kvSq[1]*w[4]*BmagInv[11]*nuVtSqSum[11]*m_-0.2795084971874738*kvSq[1]*BmagInv[2]*w[4]*nuVtSqSum[11]*m_-0.2795084971874738*kvSq[1]*nuVtSqSum[2]*w[4]*BmagInv[11]*m_-0.5*kvSq[1]*w[4]*BmagInv[10]*nuVtSqSum[10]*m_+0.07142857142857142*kvSq[1]*w[4]*BmagInv[9]*nuVtSqSum[9]*m_-0.2795084971874737*BmagInv[0]*kvSq[1]*w[4]*nuVtSqSum[9]*m_-0.2795084971874737*nuVtSqSum[0]*kvSq[1]*w[4]*BmagInv[9]*m_+0.07142857142857142*kvSq[1]*w[4]*BmagInv[8]*nuVtSqSum[8]*m_-0.2795084971874737*BmagInv[0]*kvSq[1]*w[4]*nuVtSqSum[8]*m_-0.2795084971874737*nuVtSqSum[0]*kvSq[1]*w[4]*BmagInv[8]*m_+0.07142857142857142*kvSq[1]*w[4]*BmagInv[7]*nuVtSqSum[7]*m_-0.2795084971874737*BmagInv[0]*kvSq[1]*w[4]*nuVtSqSum[7]*m_-0.2795084971874737*nuVtSqSum[0]*kvSq[1]*w[4]*BmagInv[7]*m_-0.25*kvSq[1]*w[4]*BmagInv[6]*nuVtSqSum[6]*m_-0.25*kvSq[1]*w[4]*BmagInv[5]*nuVtSqSum[5]*m_-0.25*kvSq[1]*BmagInv[4]*nuVtSqSum[4]*w[4]*m_+0.25*BmagInv[0]*nuVtSqSum[0]*kvSq[1]*w[4]*m_))); 

}
