#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf2x2vSer_Vpar_P1(const double m_, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:         Cell-center coordinates.
  // dxv[4]:       Cell spacing.
  // idx[4]:       current grid index.
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  if (idxr[2] == 1) {

    outr[3] += (nuVtSqSum[3]*(0.4330127018922193*fr[5]-0.75*fr[11])+nuVtSqSum[2]*(0.4330127018922193*fr[2]-0.75*fr[7])+nuVtSqSum[1]*(0.4330127018922193*fr[1]-0.75*fr[6])+nuVtSqSum[0]*(0.4330127018922193*fr[0]-0.75*fr[3]))*rdvSq4R; 
    outr[6] += (nuVtSqSum[2]*(0.4330127018922193*fr[5]-0.75*fr[11])+nuVtSqSum[3]*(0.4330127018922193*fr[2]-0.75*fr[7])+nuVtSqSum[0]*(0.4330127018922193*fr[1]-0.75*fr[6])+nuVtSqSum[1]*(0.4330127018922193*fr[0]-0.75*fr[3]))*rdvSq4R; 
    outr[7] += (nuVtSqSum[1]*(0.4330127018922193*fr[5]-0.75*fr[11])+nuVtSqSum[0]*(0.4330127018922193*fr[2]-0.75*fr[7])+nuVtSqSum[3]*(0.4330127018922193*fr[1]-0.75*fr[6])+nuVtSqSum[2]*(0.4330127018922193*fr[0]-0.75*fr[3]))*rdvSq4R; 
    outr[10] += (nuVtSqSum[3]*(0.4330127018922193*fr[12]-0.75*fr[15])+nuVtSqSum[2]*(0.4330127018922193*fr[9]-0.75*fr[14])+nuVtSqSum[1]*(0.4330127018922193*fr[8]-0.75*fr[13])+nuVtSqSum[0]*(0.4330127018922193*fr[4]-0.75*fr[10]))*rdvSq4R; 
    outr[11] += (nuVtSqSum[0]*(0.4330127018922193*fr[5]-0.75*fr[11])+nuVtSqSum[1]*(0.4330127018922193*fr[2]-0.75*fr[7])+nuVtSqSum[2]*(0.4330127018922193*fr[1]-0.75*fr[6])+(0.4330127018922193*fr[0]-0.75*fr[3])*nuVtSqSum[3])*rdvSq4R; 
    outr[13] += (nuVtSqSum[2]*(0.4330127018922193*fr[12]-0.75*fr[15])+nuVtSqSum[3]*(0.4330127018922193*fr[9]-0.75*fr[14])+nuVtSqSum[0]*(0.4330127018922193*fr[8]-0.75*fr[13])+nuVtSqSum[1]*(0.4330127018922193*fr[4]-0.75*fr[10]))*rdvSq4R; 
    outr[14] += (nuVtSqSum[1]*(0.4330127018922193*fr[12]-0.75*fr[15])+nuVtSqSum[0]*(0.4330127018922193*fr[9]-0.75*fr[14])+nuVtSqSum[3]*(0.4330127018922193*fr[8]-0.75*fr[13])+nuVtSqSum[2]*(0.4330127018922193*fr[4]-0.75*fr[10]))*rdvSq4R; 
    outr[15] += (nuVtSqSum[0]*(0.4330127018922193*fr[12]-0.75*fr[15])+nuVtSqSum[1]*(0.4330127018922193*fr[9]-0.75*fr[14])+nuVtSqSum[2]*(0.4330127018922193*fr[8]-0.75*fr[13])+nuVtSqSum[3]*(0.4330127018922193*fr[4]-0.75*fr[10]))*rdvSq4R; 

  } else {

    outl[3] += (nuVtSqSum[3]*((-0.75*fl[11])-0.4330127018922193*fl[5])+nuVtSqSum[2]*((-0.75*fl[7])-0.4330127018922193*fl[2])+nuVtSqSum[1]*((-0.75*fl[6])-0.4330127018922193*fl[1])+nuVtSqSum[0]*((-0.75*fl[3])-0.4330127018922193*fl[0]))*rdvSq4L; 
    outl[6] += (nuVtSqSum[2]*((-0.75*fl[11])-0.4330127018922193*fl[5])+nuVtSqSum[3]*((-0.75*fl[7])-0.4330127018922193*fl[2])+nuVtSqSum[0]*((-0.75*fl[6])-0.4330127018922193*fl[1])+nuVtSqSum[1]*((-0.75*fl[3])-0.4330127018922193*fl[0]))*rdvSq4L; 
    outl[7] += (nuVtSqSum[1]*((-0.75*fl[11])-0.4330127018922193*fl[5])+nuVtSqSum[0]*((-0.75*fl[7])-0.4330127018922193*fl[2])+nuVtSqSum[3]*((-0.75*fl[6])-0.4330127018922193*fl[1])+nuVtSqSum[2]*((-0.75*fl[3])-0.4330127018922193*fl[0]))*rdvSq4L; 
    outl[10] += (nuVtSqSum[3]*((-0.75*fl[15])-0.4330127018922193*fl[12])+nuVtSqSum[2]*((-0.75*fl[14])-0.4330127018922193*fl[9])+nuVtSqSum[1]*((-0.75*fl[13])-0.4330127018922193*fl[8])+nuVtSqSum[0]*((-0.75*fl[10])-0.4330127018922193*fl[4]))*rdvSq4L; 
    outl[11] += (nuVtSqSum[0]*((-0.75*fl[11])-0.4330127018922193*fl[5])+nuVtSqSum[1]*((-0.75*fl[7])-0.4330127018922193*fl[2])+nuVtSqSum[2]*((-0.75*fl[6])-0.4330127018922193*fl[1])+((-0.75*fl[3])-0.4330127018922193*fl[0])*nuVtSqSum[3])*rdvSq4L; 
    outl[13] += (nuVtSqSum[2]*((-0.75*fl[15])-0.4330127018922193*fl[12])+nuVtSqSum[3]*((-0.75*fl[14])-0.4330127018922193*fl[9])+nuVtSqSum[0]*((-0.75*fl[13])-0.4330127018922193*fl[8])+nuVtSqSum[1]*((-0.75*fl[10])-0.4330127018922193*fl[4]))*rdvSq4L; 
    outl[14] += (nuVtSqSum[1]*((-0.75*fl[15])-0.4330127018922193*fl[12])+nuVtSqSum[0]*((-0.75*fl[14])-0.4330127018922193*fl[9])+nuVtSqSum[3]*((-0.75*fl[13])-0.4330127018922193*fl[8])+nuVtSqSum[2]*((-0.75*fl[10])-0.4330127018922193*fl[4]))*rdvSq4L; 
    outl[15] += (nuVtSqSum[0]*((-0.75*fl[15])-0.4330127018922193*fl[12])+nuVtSqSum[1]*((-0.75*fl[14])-0.4330127018922193*fl[9])+nuVtSqSum[2]*((-0.75*fl[13])-0.4330127018922193*fl[8])+nuVtSqSum[3]*((-0.75*fl[10])-0.4330127018922193*fl[4]))*rdvSq4L; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf2x2vSer_Mu_P1(const double m_, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:         Cell-center coordinates.
  // dxv[4]:       Cell spacing.
  // idx[4]:       current grid index.
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  double diffFac[4]; 
  diffFac[0] = 0.5*(BmagInv[3]*nuVtSqSum[3]+BmagInv[2]*nuVtSqSum[2]+BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*m_; 
  diffFac[1] = 0.5*(BmagInv[2]*nuVtSqSum[3]+nuVtSqSum[2]*BmagInv[3]+BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*m_; 
  diffFac[2] = 0.5*(BmagInv[1]*nuVtSqSum[3]+nuVtSqSum[1]*BmagInv[3]+BmagInv[0]*nuVtSqSum[2]+nuVtSqSum[0]*BmagInv[2])*m_; 
  diffFac[3] = 0.5*(BmagInv[0]*nuVtSqSum[3]+nuVtSqSum[0]*BmagInv[3]+BmagInv[1]*nuVtSqSum[2]+nuVtSqSum[1]*BmagInv[2])*m_; 

  if (idxr[3] == 1) {

    outr[4] += (diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[12]+diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[9]+diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[8]+diffFac[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[5]+diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[4]+0.8660254037844386*(diffFac[2]*fr[2]+diffFac[1]*fr[1]+diffFac[0]*fr[0])*wr[3]-0.4330127018922193*(diffFac[2]*fr[2]+diffFac[1]*fr[1]+diffFac[0]*fr[0])*dxvr[3])*rdvSq4R; 
    outr[8] += (diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[12]+diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[9]+diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[8]+diffFac[2]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[5]+diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[4]+0.8660254037844386*(fr[2]*diffFac[3]+diffFac[0]*fr[1]+fr[0]*diffFac[1])*wr[3]-0.4330127018922193*(fr[2]*diffFac[3]+diffFac[0]*fr[1]+fr[0]*diffFac[1])*dxvr[3])*rdvSq4R; 
    outr[9] += (diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[12]+diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[9]+diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[8]+diffFac[1]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[5]+diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[4]+0.8660254037844386*(fr[1]*diffFac[3]+diffFac[0]*fr[2]+fr[0]*diffFac[2])*wr[3]-0.4330127018922193*(fr[1]*diffFac[3]+diffFac[0]*fr[2]+fr[0]*diffFac[2])*dxvr[3])*rdvSq4R; 
    outr[10] += (diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[15]+diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[14]+diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[13]+diffFac[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[11]+diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[10]+diffFac[2]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[7]+diffFac[1]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[6]+diffFac[0]*fr[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3]))*rdvSq4R; 
    outr[12] += (diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[12]+diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[9]+diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[8]+diffFac[0]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[5]+diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[4]+0.8660254037844386*(fr[0]*diffFac[3]+diffFac[1]*fr[2]+fr[1]*diffFac[2])*wr[3]-0.4330127018922193*(fr[0]*diffFac[3]+diffFac[1]*fr[2]+fr[1]*diffFac[2])*dxvr[3])*rdvSq4R; 
    outr[13] += (diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[15]+diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[14]+diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[13]+diffFac[2]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[11]+diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[10]+diffFac[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[7]+diffFac[0]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[6]+diffFac[1]*fr[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3]))*rdvSq4R; 
    outr[14] += (diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[15]+diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[14]+diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[13]+diffFac[1]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[11]+diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[10]+diffFac[0]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[7]+diffFac[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[6]+diffFac[2]*fr[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3]))*rdvSq4R; 
    outr[15] += (diffFac[0]*(0.75*dxvr[3]-1.5*wr[3])*fr[15]+diffFac[1]*(0.75*dxvr[3]-1.5*wr[3])*fr[14]+diffFac[2]*(0.75*dxvr[3]-1.5*wr[3])*fr[13]+diffFac[0]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[11]+diffFac[3]*(0.75*dxvr[3]-1.5*wr[3])*fr[10]+diffFac[1]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[7]+diffFac[2]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3])*fr[6]+diffFac[3]*fr[3]*(0.8660254037844386*wr[3]-0.4330127018922193*dxvr[3]))*rdvSq4R; 

  } else {

    outl[4] += (diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[12]+diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[9]+diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[8]+diffFac[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[5]+diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[4]-0.8660254037844386*(diffFac[2]*fl[2]+diffFac[1]*fl[1]+diffFac[0]*fl[0])*wl[3]-0.4330127018922193*(diffFac[2]*fl[2]+diffFac[1]*fl[1]+diffFac[0]*fl[0])*dxvl[3])*rdvSq4L; 
    outl[8] += (diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[12]+diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[9]+diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[8]+diffFac[2]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[5]+diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[4]-0.8660254037844386*(fl[2]*diffFac[3]+diffFac[0]*fl[1]+fl[0]*diffFac[1])*wl[3]-0.4330127018922193*(fl[2]*diffFac[3]+diffFac[0]*fl[1]+fl[0]*diffFac[1])*dxvl[3])*rdvSq4L; 
    outl[9] += (diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[12]+diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[9]+diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[8]+diffFac[1]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[5]+diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[4]-0.8660254037844386*(fl[1]*diffFac[3]+diffFac[0]*fl[2]+fl[0]*diffFac[2])*wl[3]-0.4330127018922193*(fl[1]*diffFac[3]+diffFac[0]*fl[2]+fl[0]*diffFac[2])*dxvl[3])*rdvSq4L; 
    outl[10] += (diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[15]+diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[14]+diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[13]+diffFac[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[11]+diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[10]+diffFac[2]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[7]+diffFac[1]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[6]+diffFac[0]*fl[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3]))*rdvSq4L; 
    outl[12] += (diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[12]+diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[9]+diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[8]+diffFac[0]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[5]+diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[4]-0.8660254037844386*(fl[0]*diffFac[3]+diffFac[1]*fl[2]+fl[1]*diffFac[2])*wl[3]-0.4330127018922193*(fl[0]*diffFac[3]+diffFac[1]*fl[2]+fl[1]*diffFac[2])*dxvl[3])*rdvSq4L; 
    outl[13] += (diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[15]+diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[14]+diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[13]+diffFac[2]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[11]+diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[10]+diffFac[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[7]+diffFac[0]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[6]+diffFac[1]*fl[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3]))*rdvSq4L; 
    outl[14] += (diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[15]+diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[14]+diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[13]+diffFac[1]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[11]+diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[10]+diffFac[0]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[7]+diffFac[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[6]+diffFac[2]*fl[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3]))*rdvSq4L; 
    outl[15] += (diffFac[0]*((-1.5*wl[3])-0.75*dxvl[3])*fl[15]+diffFac[1]*((-1.5*wl[3])-0.75*dxvl[3])*fl[14]+diffFac[2]*((-1.5*wl[3])-0.75*dxvl[3])*fl[13]+diffFac[0]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[11]+diffFac[3]*((-1.5*wl[3])-0.75*dxvl[3])*fl[10]+diffFac[1]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[7]+diffFac[2]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3])*fl[6]+diffFac[3]*fl[3]*((-0.8660254037844386*wl[3])-0.4330127018922193*dxvl[3]))*rdvSq4L; 

  }
  return 0.0; 
} 
