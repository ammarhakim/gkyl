#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf2x2vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:       Cell-center coordinates.
  // dxv[4]:     Cell spacing.
  // idx[4]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[8]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  if (idxr[2] == 1) {

    outr[3] += (-0.75*nuVtSqSum[3]*fr[11]*rdvSq4R)-0.75*nuVtSqSum[2]*fr[7]*rdvSq4R-0.75*nuVtSqSum[1]*fr[6]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[5]*rdvSq4R-0.75*nuVtSqSum[0]*fr[3]*rdvSq4R+0.4330127018922193*fr[2]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[1]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[6] += (-0.75*nuVtSqSum[2]*fr[11]*rdvSq4R)-0.75*nuVtSqSum[3]*fr[7]*rdvSq4R-0.75*nuVtSqSum[0]*fr[6]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[5]*rdvSq4R+0.4330127018922193*fr[2]*nuVtSqSum[3]*rdvSq4R-0.75*nuVtSqSum[1]*fr[3]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[1]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[7] += (-0.75*nuVtSqSum[1]*fr[11]*rdvSq4R)-0.75*nuVtSqSum[0]*fr[7]*rdvSq4R-0.75*nuVtSqSum[3]*fr[6]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[5]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[3]*rdvSq4R-0.75*nuVtSqSum[2]*fr[3]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[2]*rdvSq4R; 
    outr[10] += (-0.75*nuVtSqSum[3]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[2]*fr[14]*rdvSq4R-0.75*nuVtSqSum[1]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[12]*rdvSq4R-0.75*nuVtSqSum[0]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[9]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[4]*rdvSq4R; 
    outr[11] += (-0.75*nuVtSqSum[0]*fr[11]*rdvSq4R)-0.75*nuVtSqSum[1]*fr[7]*rdvSq4R-0.75*nuVtSqSum[2]*fr[6]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[5]*rdvSq4R-0.75*fr[3]*nuVtSqSum[3]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[3]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[2]*rdvSq4R; 
    outr[13] += (-0.75*nuVtSqSum[2]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[3]*fr[14]*rdvSq4R-0.75*nuVtSqSum[0]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[12]*rdvSq4R-0.75*nuVtSqSum[1]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[9]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[4]*rdvSq4R; 
    outr[14] += (-0.75*nuVtSqSum[1]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[0]*fr[14]*rdvSq4R-0.75*nuVtSqSum[3]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[12]*rdvSq4R-0.75*nuVtSqSum[2]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[9]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[4]*rdvSq4R; 
    outr[15] += (-0.75*nuVtSqSum[0]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[1]*fr[14]*rdvSq4R-0.75*nuVtSqSum[2]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[12]*rdvSq4R-0.75*nuVtSqSum[3]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[9]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[4]*rdvSq4R; 

  } else {

    outl[3] += (-0.75*nuVtSqSum[3]*fl[11]*rdvSq4L)-0.75*nuVtSqSum[2]*fl[7]*rdvSq4L-0.75*nuVtSqSum[1]*fl[6]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[5]*rdvSq4L-0.75*nuVtSqSum[0]*fl[3]*rdvSq4L-0.4330127018922193*fl[2]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[1]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[6] += (-0.75*nuVtSqSum[2]*fl[11]*rdvSq4L)-0.75*nuVtSqSum[3]*fl[7]*rdvSq4L-0.75*nuVtSqSum[0]*fl[6]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[5]*rdvSq4L-0.4330127018922193*fl[2]*nuVtSqSum[3]*rdvSq4L-0.75*nuVtSqSum[1]*fl[3]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[1]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[7] += (-0.75*nuVtSqSum[1]*fl[11]*rdvSq4L)-0.75*nuVtSqSum[0]*fl[7]*rdvSq4L-0.75*nuVtSqSum[3]*fl[6]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[5]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[3]*rdvSq4L-0.75*nuVtSqSum[2]*fl[3]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[2]*rdvSq4L; 
    outl[10] += (-0.75*nuVtSqSum[3]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[2]*fl[14]*rdvSq4L-0.75*nuVtSqSum[1]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[12]*rdvSq4L-0.75*nuVtSqSum[0]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[9]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[4]*rdvSq4L; 
    outl[11] += (-0.75*nuVtSqSum[0]*fl[11]*rdvSq4L)-0.75*nuVtSqSum[1]*fl[7]*rdvSq4L-0.75*nuVtSqSum[2]*fl[6]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[5]*rdvSq4L-0.75*fl[3]*nuVtSqSum[3]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[3]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[2]*rdvSq4L; 
    outl[13] += (-0.75*nuVtSqSum[2]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[3]*fl[14]*rdvSq4L-0.75*nuVtSqSum[0]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[12]*rdvSq4L-0.75*nuVtSqSum[1]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[9]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[4]*rdvSq4L; 
    outl[14] += (-0.75*nuVtSqSum[1]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[0]*fl[14]*rdvSq4L-0.75*nuVtSqSum[3]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[12]*rdvSq4L-0.75*nuVtSqSum[2]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[9]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[4]*rdvSq4L; 
    outl[15] += (-0.75*nuVtSqSum[0]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[1]*fl[14]*rdvSq4L-0.75*nuVtSqSum[2]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[12]*rdvSq4L-0.75*nuVtSqSum[3]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[9]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[4]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf2x2vTensor_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:       Cell-center coordinates.
  // dxv[4]:     Cell spacing.
  // idx[4]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[8]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  if (idxr[3] == 1) {

    outr[4] += (-0.75*nuVtSqSum[3]*fr[12]*rdvSq4R)-0.75*nuVtSqSum[2]*fr[9]*rdvSq4R-0.75*nuVtSqSum[1]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[5]*rdvSq4R-0.75*nuVtSqSum[0]*fr[4]*rdvSq4R+0.4330127018922193*fr[2]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[1]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[8] += (-0.75*nuVtSqSum[2]*fr[12]*rdvSq4R)-0.75*nuVtSqSum[3]*fr[9]*rdvSq4R-0.75*nuVtSqSum[0]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[5]*rdvSq4R-0.75*nuVtSqSum[1]*fr[4]*rdvSq4R+0.4330127018922193*fr[2]*nuVtSqSum[3]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[1]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[9] += (-0.75*nuVtSqSum[1]*fr[12]*rdvSq4R)-0.75*nuVtSqSum[0]*fr[9]*rdvSq4R-0.75*nuVtSqSum[3]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[5]*rdvSq4R-0.75*nuVtSqSum[2]*fr[4]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[3]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[2]*rdvSq4R; 
    outr[10] += (-0.75*nuVtSqSum[3]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[2]*fr[14]*rdvSq4R-0.75*nuVtSqSum[1]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[11]*rdvSq4R-0.75*nuVtSqSum[0]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[7]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[6]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[3]*rdvSq4R; 
    outr[12] += (-0.75*nuVtSqSum[0]*fr[12]*rdvSq4R)-0.75*nuVtSqSum[1]*fr[9]*rdvSq4R-0.75*nuVtSqSum[2]*fr[8]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[5]*rdvSq4R-0.75*nuVtSqSum[3]*fr[4]*rdvSq4R+0.4330127018922193*fr[0]*nuVtSqSum[3]*rdvSq4R+0.4330127018922193*fr[1]*nuVtSqSum[2]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[2]*rdvSq4R; 
    outr[13] += (-0.75*nuVtSqSum[2]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[3]*fr[14]*rdvSq4R-0.75*nuVtSqSum[0]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[11]*rdvSq4R-0.75*nuVtSqSum[1]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[7]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[6]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[3]*rdvSq4R; 
    outr[14] += (-0.75*nuVtSqSum[1]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[0]*fr[14]*rdvSq4R-0.75*nuVtSqSum[3]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[11]*rdvSq4R-0.75*nuVtSqSum[2]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[7]*rdvSq4R+0.4330127018922193*nuVtSqSum[3]*fr[6]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[3]*rdvSq4R; 
    outr[15] += (-0.75*nuVtSqSum[0]*fr[15]*rdvSq4R)-0.75*nuVtSqSum[1]*fr[14]*rdvSq4R-0.75*nuVtSqSum[2]*fr[13]*rdvSq4R+0.4330127018922193*nuVtSqSum[0]*fr[11]*rdvSq4R-0.75*nuVtSqSum[3]*fr[10]*rdvSq4R+0.4330127018922193*nuVtSqSum[1]*fr[7]*rdvSq4R+0.4330127018922193*nuVtSqSum[2]*fr[6]*rdvSq4R+0.4330127018922193*fr[3]*nuVtSqSum[3]*rdvSq4R; 

  } else {

    outl[4] += (-0.75*nuVtSqSum[3]*fl[12]*rdvSq4L)-0.75*nuVtSqSum[2]*fl[9]*rdvSq4L-0.75*nuVtSqSum[1]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[5]*rdvSq4L-0.75*nuVtSqSum[0]*fl[4]*rdvSq4L-0.4330127018922193*fl[2]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[1]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[8] += (-0.75*nuVtSqSum[2]*fl[12]*rdvSq4L)-0.75*nuVtSqSum[3]*fl[9]*rdvSq4L-0.75*nuVtSqSum[0]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[5]*rdvSq4L-0.75*nuVtSqSum[1]*fl[4]*rdvSq4L-0.4330127018922193*fl[2]*nuVtSqSum[3]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[1]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[9] += (-0.75*nuVtSqSum[1]*fl[12]*rdvSq4L)-0.75*nuVtSqSum[0]*fl[9]*rdvSq4L-0.75*nuVtSqSum[3]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[5]*rdvSq4L-0.75*nuVtSqSum[2]*fl[4]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[3]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[2]*rdvSq4L; 
    outl[10] += (-0.75*nuVtSqSum[3]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[2]*fl[14]*rdvSq4L-0.75*nuVtSqSum[1]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[11]*rdvSq4L-0.75*nuVtSqSum[0]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[7]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[6]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[3]*rdvSq4L; 
    outl[12] += (-0.75*nuVtSqSum[0]*fl[12]*rdvSq4L)-0.75*nuVtSqSum[1]*fl[9]*rdvSq4L-0.75*nuVtSqSum[2]*fl[8]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[5]*rdvSq4L-0.75*nuVtSqSum[3]*fl[4]*rdvSq4L-0.4330127018922193*fl[0]*nuVtSqSum[3]*rdvSq4L-0.4330127018922193*fl[1]*nuVtSqSum[2]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[2]*rdvSq4L; 
    outl[13] += (-0.75*nuVtSqSum[2]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[3]*fl[14]*rdvSq4L-0.75*nuVtSqSum[0]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[11]*rdvSq4L-0.75*nuVtSqSum[1]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[7]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[6]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[3]*rdvSq4L; 
    outl[14] += (-0.75*nuVtSqSum[1]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[0]*fl[14]*rdvSq4L-0.75*nuVtSqSum[3]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[11]*rdvSq4L-0.75*nuVtSqSum[2]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[7]*rdvSq4L-0.4330127018922193*nuVtSqSum[3]*fl[6]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[3]*rdvSq4L; 
    outl[15] += (-0.75*nuVtSqSum[0]*fl[15]*rdvSq4L)-0.75*nuVtSqSum[1]*fl[14]*rdvSq4L-0.75*nuVtSqSum[2]*fl[13]*rdvSq4L-0.4330127018922193*nuVtSqSum[0]*fl[11]*rdvSq4L-0.75*nuVtSqSum[3]*fl[10]*rdvSq4L-0.4330127018922193*nuVtSqSum[1]*fl[7]*rdvSq4L-0.4330127018922193*nuVtSqSum[2]*fl[6]*rdvSq4L-0.4330127018922193*fl[3]*nuVtSqSum[3]*rdvSq4L; 

  }
  return 0.0; 
} 
