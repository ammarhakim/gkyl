#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf1x2vSer_Vpar_P1(const double m_, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates.
  // dxv[3]:       Cell spacing.
  // idx[3]:       current grid index.
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  if (idxr[1] == 1) {

    outr[2] += (nuVtSqSum[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[4])+nuVtSqSum[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[2]))*rdvSq4R; 
    outr[4] += (nuVtSqSum[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[4])+nuVtSqSum[1]*(0.6123724356957944*fr[0]-1.060660171779821*fr[2]))*rdvSq4R; 
    outr[6] += (nuVtSqSum[1]*(0.6123724356957944*fr[5]-1.060660171779821*fr[7])+nuVtSqSum[0]*(0.6123724356957944*fr[3]-1.060660171779821*fr[6]))*rdvSq4R; 
    outr[7] += (nuVtSqSum[0]*(0.6123724356957944*fr[5]-1.060660171779821*fr[7])+nuVtSqSum[1]*(0.6123724356957944*fr[3]-1.060660171779821*fr[6]))*rdvSq4R; 

  } else {

    outl[2] += (nuVtSqSum[1]*((-1.060660171779821*fl[4])-0.6123724356957944*fl[1])+nuVtSqSum[0]*((-1.060660171779821*fl[2])-0.6123724356957944*fl[0]))*rdvSq4L; 
    outl[4] += (nuVtSqSum[0]*((-1.060660171779821*fl[4])-0.6123724356957944*fl[1])+nuVtSqSum[1]*((-1.060660171779821*fl[2])-0.6123724356957944*fl[0]))*rdvSq4L; 
    outl[6] += (nuVtSqSum[1]*((-1.060660171779821*fl[7])-0.6123724356957944*fl[5])+nuVtSqSum[0]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[3]))*rdvSq4L; 
    outl[7] += (nuVtSqSum[0]*((-1.060660171779821*fl[7])-0.6123724356957944*fl[5])+nuVtSqSum[1]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[3]))*rdvSq4L; 

  }
  return 0.0; 
} 
double GkLBOconstNuBoundarySurf1x2vSer_Mu_P1(const double m_, const double *positivityWeightByDirL, const double *positivityWeightByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates.
  // dxv[3]:       Cell spacing.
  // idx[3]:       current grid index.
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  double diffFac[2]; 
  diffFac[0] = 0.7071067811865475*(BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*m_; 
  diffFac[1] = 0.7071067811865475*(BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*m_; 

  if (idxr[2] == 1) {

    outr[3] += (diffFac[1]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[5]+diffFac[0]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[3]+1.224744871391589*(diffFac[1]*fr[1]+diffFac[0]*fr[0])*wr[2]-0.6123724356957944*(diffFac[1]*fr[1]+diffFac[0]*fr[0])*dxvr[2])*rdvSq4R; 
    outr[5] += (diffFac[0]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[5]+diffFac[1]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[3]+1.224744871391589*(diffFac[0]*fr[1]+fr[0]*diffFac[1])*wr[2]-0.6123724356957944*(diffFac[0]*fr[1]+fr[0]*diffFac[1])*dxvr[2])*rdvSq4R; 
    outr[6] += (diffFac[1]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[7]+diffFac[0]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[6]+diffFac[1]*(1.224744871391589*wr[2]-0.6123724356957944*dxvr[2])*fr[4]+diffFac[0]*fr[2]*(1.224744871391589*wr[2]-0.6123724356957944*dxvr[2]))*rdvSq4R; 
    outr[7] += (diffFac[0]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[7]+diffFac[1]*(1.060660171779821*dxvr[2]-2.121320343559642*wr[2])*fr[6]+diffFac[0]*(1.224744871391589*wr[2]-0.6123724356957944*dxvr[2])*fr[4]+diffFac[1]*fr[2]*(1.224744871391589*wr[2]-0.6123724356957944*dxvr[2]))*rdvSq4R; 

  } else {

    outl[3] += (diffFac[1]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[5]+diffFac[0]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[3]-1.224744871391589*(diffFac[1]*fl[1]+diffFac[0]*fl[0])*wl[2]-0.6123724356957944*(diffFac[1]*fl[1]+diffFac[0]*fl[0])*dxvl[2])*rdvSq4L; 
    outl[5] += (diffFac[0]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[5]+diffFac[1]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[3]-1.224744871391589*(diffFac[0]*fl[1]+fl[0]*diffFac[1])*wl[2]-0.6123724356957944*(diffFac[0]*fl[1]+fl[0]*diffFac[1])*dxvl[2])*rdvSq4L; 
    outl[6] += (diffFac[1]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[7]+diffFac[0]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[6]+diffFac[1]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fl[4]+diffFac[0]*fl[2]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))*rdvSq4L; 
    outl[7] += (diffFac[0]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[7]+diffFac[1]*((-2.121320343559642*wl[2])-1.060660171779821*dxvl[2])*fl[6]+diffFac[0]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2])*fl[4]+diffFac[1]*fl[2]*((-1.224744871391589*wl[2])-0.6123724356957944*dxvl[2]))*rdvSq4L; 

  }
  return 0.0; 
} 
