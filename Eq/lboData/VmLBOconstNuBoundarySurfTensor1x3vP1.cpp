#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf1x3vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:       Cell-center coordinates.
  // dxv[4]:     Cell spacing.
  // idx[4]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  if (idxr[1] == 1) {

    outr[2] += (-1.060660171779821*nuVtSqSum[1]*fr[5]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[2]*rdvSq4R+0.6123724356957944*fr[1]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[5] += (-1.060660171779821*nuVtSqSum[0]*fr[5]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[2]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[7] += (-1.060660171779821*nuVtSqSum[1]*fr[11]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[7]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[6]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[3]*rdvSq4R; 
    outr[9] += (-1.060660171779821*nuVtSqSum[1]*fr[12]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[9]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[8]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[4]*rdvSq4R; 
    outr[11] += (-1.060660171779821*nuVtSqSum[0]*fr[11]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[7]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[6]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[3]*rdvSq4R; 
    outr[12] += (-1.060660171779821*nuVtSqSum[0]*fr[12]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[9]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[8]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[4]*rdvSq4R; 
    outr[14] += (-1.060660171779821*nuVtSqSum[1]*fr[15]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[14]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[13]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[10]*rdvSq4R; 
    outr[15] += (-1.060660171779821*nuVtSqSum[0]*fr[15]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[14]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[13]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[10]*rdvSq4R; 

  } else {

    outl[2] += (-1.060660171779821*nuVtSqSum[1]*fl[5]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[2]*rdvSq4L-0.6123724356957944*fl[1]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[5] += (-1.060660171779821*nuVtSqSum[0]*fl[5]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[2]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[7] += (-1.060660171779821*nuVtSqSum[1]*fl[11]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[7]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[6]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[3]*rdvSq4L; 
    outl[9] += (-1.060660171779821*nuVtSqSum[1]*fl[12]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[9]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[8]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[4]*rdvSq4L; 
    outl[11] += (-1.060660171779821*nuVtSqSum[0]*fl[11]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[7]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[6]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[3]*rdvSq4L; 
    outl[12] += (-1.060660171779821*nuVtSqSum[0]*fl[12]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[9]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[8]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[4]*rdvSq4L; 
    outl[14] += (-1.060660171779821*nuVtSqSum[1]*fl[15]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[14]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[13]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[10]*rdvSq4L; 
    outl[15] += (-1.060660171779821*nuVtSqSum[0]*fl[15]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[14]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[13]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[10]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vTensor_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:       Cell-center coordinates.
  // dxv[4]:     Cell spacing.
  // idx[4]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  if (idxr[2] == 1) {

    outr[3] += (-1.060660171779821*nuVtSqSum[1]*fr[6]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[3]*rdvSq4R+0.6123724356957944*fr[1]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[6] += (-1.060660171779821*nuVtSqSum[0]*fr[6]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[3]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[7] += (-1.060660171779821*nuVtSqSum[1]*fr[11]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[7]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[5]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[2]*rdvSq4R; 
    outr[10] += (-1.060660171779821*nuVtSqSum[1]*fr[13]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[10]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[8]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[4]*rdvSq4R; 
    outr[11] += (-1.060660171779821*nuVtSqSum[0]*fr[11]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[7]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[5]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[2]*rdvSq4R; 
    outr[13] += (-1.060660171779821*nuVtSqSum[0]*fr[13]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[10]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[8]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[4]*rdvSq4R; 
    outr[14] += (-1.060660171779821*nuVtSqSum[1]*fr[15]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[14]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[12]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[9]*rdvSq4R; 
    outr[15] += (-1.060660171779821*nuVtSqSum[0]*fr[15]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[14]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[12]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[9]*rdvSq4R; 

  } else {

    outl[3] += (-1.060660171779821*nuVtSqSum[1]*fl[6]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[3]*rdvSq4L-0.6123724356957944*fl[1]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[6] += (-1.060660171779821*nuVtSqSum[0]*fl[6]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[3]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[7] += (-1.060660171779821*nuVtSqSum[1]*fl[11]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[7]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[5]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[2]*rdvSq4L; 
    outl[10] += (-1.060660171779821*nuVtSqSum[1]*fl[13]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[10]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[8]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[4]*rdvSq4L; 
    outl[11] += (-1.060660171779821*nuVtSqSum[0]*fl[11]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[7]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[5]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[2]*rdvSq4L; 
    outl[13] += (-1.060660171779821*nuVtSqSum[0]*fl[13]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[10]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[8]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[4]*rdvSq4L; 
    outl[14] += (-1.060660171779821*nuVtSqSum[1]*fl[15]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[14]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[12]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[9]*rdvSq4L; 
    outl[15] += (-1.060660171779821*nuVtSqSum[0]*fl[15]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[14]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[12]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[9]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x3vTensor_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:       Cell-center coordinates.
  // dxv[4]:     Cell spacing.
  // idx[4]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  if (idxr[3] == 1) {

    outr[4] += (-1.060660171779821*nuVtSqSum[1]*fr[8]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[4]*rdvSq4R+0.6123724356957944*fr[1]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[8] += (-1.060660171779821*nuVtSqSum[0]*fr[8]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[4]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[9] += (-1.060660171779821*nuVtSqSum[1]*fr[12]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[9]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[5]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[2]*rdvSq4R; 
    outr[10] += (-1.060660171779821*nuVtSqSum[1]*fr[13]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[10]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[6]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[3]*rdvSq4R; 
    outr[12] += (-1.060660171779821*nuVtSqSum[0]*fr[12]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[9]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[5]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[2]*rdvSq4R; 
    outr[13] += (-1.060660171779821*nuVtSqSum[0]*fr[13]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[10]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[6]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[3]*rdvSq4R; 
    outr[14] += (-1.060660171779821*nuVtSqSum[1]*fr[15]*rdvSq4R)-1.060660171779821*nuVtSqSum[0]*fr[14]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[11]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[7]*rdvSq4R; 
    outr[15] += (-1.060660171779821*nuVtSqSum[0]*fr[15]*rdvSq4R)-1.060660171779821*nuVtSqSum[1]*fr[14]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[11]*rdvSq4R+0.6123724356957944*nuVtSqSum[1]*fr[7]*rdvSq4R; 

  } else {

    outl[4] += (-1.060660171779821*nuVtSqSum[1]*fl[8]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[4]*rdvSq4L-0.6123724356957944*fl[1]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[8] += (-1.060660171779821*nuVtSqSum[0]*fl[8]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[4]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[9] += (-1.060660171779821*nuVtSqSum[1]*fl[12]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[9]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[5]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[2]*rdvSq4L; 
    outl[10] += (-1.060660171779821*nuVtSqSum[1]*fl[13]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[10]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[6]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[3]*rdvSq4L; 
    outl[12] += (-1.060660171779821*nuVtSqSum[0]*fl[12]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[9]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[5]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[2]*rdvSq4L; 
    outl[13] += (-1.060660171779821*nuVtSqSum[0]*fl[13]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[10]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[6]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[3]*rdvSq4L; 
    outl[14] += (-1.060660171779821*nuVtSqSum[1]*fl[15]*rdvSq4L)-1.060660171779821*nuVtSqSum[0]*fl[14]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[11]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[7]*rdvSq4L; 
    outl[15] += (-1.060660171779821*nuVtSqSum[0]*fl[15]*rdvSq4L)-1.060660171779821*nuVtSqSum[1]*fl[14]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[11]*rdvSq4L-0.6123724356957944*nuVtSqSum[1]*fl[7]*rdvSq4L; 

  }
  return 0.0; 
} 
