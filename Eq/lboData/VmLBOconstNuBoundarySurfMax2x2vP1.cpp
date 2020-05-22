#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf2x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:       Cell-center coordinates.
  // dxv[4]:     Cell spacing.
  // idx[4]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  if (idxr[2] == 1) {

    outr[3] += (-0.75*nuVtSqSum[0]*fr[3]*rdvSq4R)+0.4330127018922193*fr[2]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[1]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[3] += (-0.75*nuVtSqSum[0]*fl[3]*rdvSq4L)-0.4330127018922193*fl[2]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[1]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf2x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:       Cell-center coordinates.
  // dxv[4]:     Cell spacing.
  // idx[4]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  if (idxr[3] == 1) {

    outr[4] += (-0.75*nuVtSqSum[0]*fr[4]*rdvSq4R)+0.4330127018922193*fr[2]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[1]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[4] += (-0.75*nuVtSqSum[0]*fl[4]*rdvSq4L)-0.4330127018922193*fl[2]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[1]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
