#include <GkLBOModDecl.h> 
double GkLBOconstNuBoundarySurf1x1vSer_Vpar_P1(const double m_, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:         Cell-center coordinates.
  // dxv[2]:       Cell spacing.
  // idx[2]:       current grid index.
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  if (idxr[1] == 1) {

    outr[2] += (nuVtSqSum[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[3])+nuVtSqSum[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[2]))*rdvSq4R; 
    outr[3] += (nuVtSqSum[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[3])+nuVtSqSum[1]*(0.6123724356957944*fr[0]-1.060660171779821*fr[2]))*rdvSq4R; 

  } else {

    outl[2] += (nuVtSqSum[1]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[1])+nuVtSqSum[0]*((-1.060660171779821*fl[2])-0.6123724356957944*fl[0]))*rdvSq4L; 
    outl[3] += (nuVtSqSum[0]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[1])+nuVtSqSum[1]*((-1.060660171779821*fl[2])-0.6123724356957944*fl[0]))*rdvSq4L; 

  }
  return 0.0; 
} 
