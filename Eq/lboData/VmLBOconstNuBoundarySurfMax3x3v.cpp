#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf3x3vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:       Cell-center coordinates.
  // dxv[6]:     Cell spacing.
  // idx[6]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[12]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  if (idxr[3] == 1) {

    outr[4] += (-0.5303300858899105*nuVtSqSum[0]*fr[4]*rdvSq4R)+0.3061862178478971*fr[3]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[4] += (-0.5303300858899105*nuVtSqSum[0]*fl[4]*rdvSq4L)-0.3061862178478971*fl[3]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf3x3vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:       Cell-center coordinates.
  // dxv[6]:     Cell spacing.
  // idx[6]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[30]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[10]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  if (idxr[3] == 1) {

    outr[4] += 0.6846531968814573*nuVtSqSum[0]*fr[25]*rdvSq4R+0.3061862178478971*nuVtSqSum[9]*fr[24]*rdvSq4R+0.3061862178478971*nuVtSqSum[8]*fr[23]*rdvSq4R+0.3061862178478971*nuVtSqSum[7]*fr[22]*rdvSq4R-0.5303300858899105*nuVtSqSum[3]*fr[12]*rdvSq4R-0.5303300858899105*nuVtSqSum[2]*fr[11]*rdvSq4R-0.5303300858899105*nuVtSqSum[1]*fr[10]*rdvSq4R+0.3061862178478971*nuVtSqSum[6]*fr[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[5]*fr[8]*rdvSq4R+0.3061862178478971*nuVtSqSum[4]*fr[7]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[4]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[10] += 0.6846531968814573*nuVtSqSum[1]*fr[25]*rdvSq4R+0.273861278752583*nuVtSqSum[1]*fr[22]*rdvSq4R-0.5303300858899105*nuVtSqSum[5]*fr[12]*rdvSq4R-0.5303300858899105*nuVtSqSum[4]*fr[11]*rdvSq4R-0.4743416490252568*nuVtSqSum[7]*fr[10]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[10]*rdvSq4R+0.3061862178478971*nuVtSqSum[3]*fr[8]*rdvSq4R+0.273861278752583*fr[1]*nuVtSqSum[7]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[7]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[5]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[4]*rdvSq4R-0.5303300858899105*nuVtSqSum[1]*fr[4]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[11] += 0.6846531968814573*nuVtSqSum[2]*fr[25]*rdvSq4R+0.273861278752583*nuVtSqSum[2]*fr[23]*rdvSq4R-0.5303300858899105*nuVtSqSum[6]*fr[12]*rdvSq4R-0.4743416490252568*nuVtSqSum[8]*fr[11]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[11]*rdvSq4R-0.5303300858899105*nuVtSqSum[4]*fr[10]*rdvSq4R+0.3061862178478971*nuVtSqSum[3]*fr[9]*rdvSq4R+0.273861278752583*fr[2]*nuVtSqSum[8]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[7]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[6]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[4]*rdvSq4R-0.5303300858899105*nuVtSqSum[2]*fr[4]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[2]*rdvSq4R; 
    outr[12] += 0.6846531968814573*nuVtSqSum[3]*fr[25]*rdvSq4R+0.273861278752583*nuVtSqSum[3]*fr[24]*rdvSq4R-0.4743416490252568*nuVtSqSum[9]*fr[12]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[12]*rdvSq4R-0.5303300858899105*nuVtSqSum[6]*fr[11]*rdvSq4R-0.5303300858899105*nuVtSqSum[5]*fr[10]*rdvSq4R+0.273861278752583*fr[3]*nuVtSqSum[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[8]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[6]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[5]*rdvSq4R-0.5303300858899105*nuVtSqSum[3]*fr[4]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[3]*rdvSq4R; 
    outr[16] += (-0.5303300858899105*nuVtSqSum[0]*fr[16]*rdvSq4R)+0.3061862178478971*nuVtSqSum[3]*fr[15]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[14]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[13]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[5]*rdvSq4R; 
    outr[20] += (-0.5303300858899105*nuVtSqSum[0]*fr[20]*rdvSq4R)+0.3061862178478971*nuVtSqSum[3]*fr[19]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[18]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[17]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[6]*rdvSq4R; 
    outr[25] += (-2.651650429449552*nuVtSqSum[0]*fr[25]*rdvSq4R)-1.185854122563142*nuVtSqSum[9]*fr[24]*rdvSq4R-1.185854122563142*nuVtSqSum[8]*fr[23]*rdvSq4R-1.185854122563142*nuVtSqSum[7]*fr[22]*rdvSq4R+2.053959590644372*nuVtSqSum[3]*fr[12]*rdvSq4R+2.053959590644372*nuVtSqSum[2]*fr[11]*rdvSq4R+2.053959590644372*nuVtSqSum[1]*fr[10]*rdvSq4R-1.185854122563142*nuVtSqSum[6]*fr[9]*rdvSq4R-1.185854122563142*nuVtSqSum[5]*fr[8]*rdvSq4R-1.185854122563142*nuVtSqSum[4]*fr[7]*rdvSq4R+2.053959590644372*nuVtSqSum[0]*fr[4]*rdvSq4R-1.185854122563142*fr[3]*nuVtSqSum[3]*rdvSq4R-1.185854122563142*fr[2]*nuVtSqSum[2]*rdvSq4R-1.185854122563142*fr[1]*nuVtSqSum[1]*rdvSq4R-1.185854122563142*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[4] += (-0.6846531968814573*nuVtSqSum[0]*fl[25]*rdvSq4L)-0.3061862178478971*nuVtSqSum[9]*fl[24]*rdvSq4L-0.3061862178478971*nuVtSqSum[8]*fl[23]*rdvSq4L-0.3061862178478971*nuVtSqSum[7]*fl[22]*rdvSq4L-0.5303300858899105*nuVtSqSum[3]*fl[12]*rdvSq4L-0.5303300858899105*nuVtSqSum[2]*fl[11]*rdvSq4L-0.5303300858899105*nuVtSqSum[1]*fl[10]*rdvSq4L-0.3061862178478971*nuVtSqSum[6]*fl[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[5]*fl[8]*rdvSq4L-0.3061862178478971*nuVtSqSum[4]*fl[7]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[4]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[10] += (-0.6846531968814573*nuVtSqSum[1]*fl[25]*rdvSq4L)-0.273861278752583*nuVtSqSum[1]*fl[22]*rdvSq4L-0.5303300858899105*nuVtSqSum[5]*fl[12]*rdvSq4L-0.5303300858899105*nuVtSqSum[4]*fl[11]*rdvSq4L-0.4743416490252568*nuVtSqSum[7]*fl[10]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[10]*rdvSq4L-0.3061862178478971*nuVtSqSum[3]*fl[8]*rdvSq4L-0.273861278752583*fl[1]*nuVtSqSum[7]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[7]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[5]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[4]*rdvSq4L-0.5303300858899105*nuVtSqSum[1]*fl[4]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[11] += (-0.6846531968814573*nuVtSqSum[2]*fl[25]*rdvSq4L)-0.273861278752583*nuVtSqSum[2]*fl[23]*rdvSq4L-0.5303300858899105*nuVtSqSum[6]*fl[12]*rdvSq4L-0.4743416490252568*nuVtSqSum[8]*fl[11]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[11]*rdvSq4L-0.5303300858899105*nuVtSqSum[4]*fl[10]*rdvSq4L-0.3061862178478971*nuVtSqSum[3]*fl[9]*rdvSq4L-0.273861278752583*fl[2]*nuVtSqSum[8]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[7]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[6]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[4]*rdvSq4L-0.5303300858899105*nuVtSqSum[2]*fl[4]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[2]*rdvSq4L; 
    outl[12] += (-0.6846531968814573*nuVtSqSum[3]*fl[25]*rdvSq4L)-0.273861278752583*nuVtSqSum[3]*fl[24]*rdvSq4L-0.4743416490252568*nuVtSqSum[9]*fl[12]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[12]*rdvSq4L-0.5303300858899105*nuVtSqSum[6]*fl[11]*rdvSq4L-0.5303300858899105*nuVtSqSum[5]*fl[10]*rdvSq4L-0.273861278752583*fl[3]*nuVtSqSum[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[8]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[6]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[5]*rdvSq4L-0.5303300858899105*nuVtSqSum[3]*fl[4]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[3]*rdvSq4L; 
    outl[16] += (-0.5303300858899105*nuVtSqSum[0]*fl[16]*rdvSq4L)-0.3061862178478971*nuVtSqSum[3]*fl[15]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[14]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[13]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[5]*rdvSq4L; 
    outl[20] += (-0.5303300858899105*nuVtSqSum[0]*fl[20]*rdvSq4L)-0.3061862178478971*nuVtSqSum[3]*fl[19]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[18]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[17]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[6]*rdvSq4L; 
    outl[25] += (-2.651650429449552*nuVtSqSum[0]*fl[25]*rdvSq4L)-1.185854122563142*nuVtSqSum[9]*fl[24]*rdvSq4L-1.185854122563142*nuVtSqSum[8]*fl[23]*rdvSq4L-1.185854122563142*nuVtSqSum[7]*fl[22]*rdvSq4L-2.053959590644372*nuVtSqSum[3]*fl[12]*rdvSq4L-2.053959590644372*nuVtSqSum[2]*fl[11]*rdvSq4L-2.053959590644372*nuVtSqSum[1]*fl[10]*rdvSq4L-1.185854122563142*nuVtSqSum[6]*fl[9]*rdvSq4L-1.185854122563142*nuVtSqSum[5]*fl[8]*rdvSq4L-1.185854122563142*nuVtSqSum[4]*fl[7]*rdvSq4L-2.053959590644372*nuVtSqSum[0]*fl[4]*rdvSq4L-1.185854122563142*fl[3]*nuVtSqSum[3]*rdvSq4L-1.185854122563142*fl[2]*nuVtSqSum[2]*rdvSq4L-1.185854122563142*fl[1]*nuVtSqSum[1]*rdvSq4L-1.185854122563142*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf3x3vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:       Cell-center coordinates.
  // dxv[6]:     Cell spacing.
  // idx[6]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[12]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[4]*dxvl[4]); 
  double rdvSq4R = 4.0/(dxvr[4]*dxvr[4]); 

  if (idxr[4] == 1) {

    outr[5] += (-0.5303300858899105*nuVtSqSum[0]*fr[5]*rdvSq4R)+0.3061862178478971*fr[3]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[5] += (-0.5303300858899105*nuVtSqSum[0]*fl[5]*rdvSq4L)-0.3061862178478971*fl[3]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf3x3vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:       Cell-center coordinates.
  // dxv[6]:     Cell spacing.
  // idx[6]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[30]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[10]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[4]*dxvl[4]); 
  double rdvSq4R = 4.0/(dxvr[4]*dxvr[4]); 

  if (idxr[4] == 1) {

    outr[5] += 0.6846531968814573*nuVtSqSum[0]*fr[26]*rdvSq4R+0.3061862178478971*nuVtSqSum[9]*fr[24]*rdvSq4R+0.3061862178478971*nuVtSqSum[8]*fr[23]*rdvSq4R+0.3061862178478971*nuVtSqSum[7]*fr[22]*rdvSq4R-0.5303300858899105*nuVtSqSum[3]*fr[15]*rdvSq4R-0.5303300858899105*nuVtSqSum[2]*fr[14]*rdvSq4R-0.5303300858899105*nuVtSqSum[1]*fr[13]*rdvSq4R+0.3061862178478971*nuVtSqSum[6]*fr[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[5]*fr[8]*rdvSq4R+0.3061862178478971*nuVtSqSum[4]*fr[7]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[5]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[13] += 0.6846531968814573*nuVtSqSum[1]*fr[26]*rdvSq4R+0.273861278752583*nuVtSqSum[1]*fr[22]*rdvSq4R-0.5303300858899105*nuVtSqSum[5]*fr[15]*rdvSq4R-0.5303300858899105*nuVtSqSum[4]*fr[14]*rdvSq4R-0.4743416490252568*nuVtSqSum[7]*fr[13]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[13]*rdvSq4R+0.3061862178478971*nuVtSqSum[3]*fr[8]*rdvSq4R+0.273861278752583*fr[1]*nuVtSqSum[7]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[7]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[5]*rdvSq4R-0.5303300858899105*nuVtSqSum[1]*fr[5]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[4]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[14] += 0.6846531968814573*nuVtSqSum[2]*fr[26]*rdvSq4R+0.273861278752583*nuVtSqSum[2]*fr[23]*rdvSq4R-0.5303300858899105*nuVtSqSum[6]*fr[15]*rdvSq4R-0.4743416490252568*nuVtSqSum[8]*fr[14]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[14]*rdvSq4R-0.5303300858899105*nuVtSqSum[4]*fr[13]*rdvSq4R+0.3061862178478971*nuVtSqSum[3]*fr[9]*rdvSq4R+0.273861278752583*fr[2]*nuVtSqSum[8]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[7]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[6]*rdvSq4R-0.5303300858899105*nuVtSqSum[2]*fr[5]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[4]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[2]*rdvSq4R; 
    outr[15] += 0.6846531968814573*nuVtSqSum[3]*fr[26]*rdvSq4R+0.273861278752583*nuVtSqSum[3]*fr[24]*rdvSq4R-0.4743416490252568*nuVtSqSum[9]*fr[15]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[15]*rdvSq4R-0.5303300858899105*nuVtSqSum[6]*fr[14]*rdvSq4R-0.5303300858899105*nuVtSqSum[5]*fr[13]*rdvSq4R+0.273861278752583*fr[3]*nuVtSqSum[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[8]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[6]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[5]*rdvSq4R-0.5303300858899105*nuVtSqSum[3]*fr[5]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[3]*rdvSq4R; 
    outr[16] += (-0.5303300858899105*nuVtSqSum[0]*fr[16]*rdvSq4R)+0.3061862178478971*nuVtSqSum[3]*fr[12]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[11]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[10]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[4]*rdvSq4R; 
    outr[21] += (-0.5303300858899105*nuVtSqSum[0]*fr[21]*rdvSq4R)+0.3061862178478971*nuVtSqSum[3]*fr[19]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[18]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[17]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[6]*rdvSq4R; 
    outr[26] += (-2.651650429449552*nuVtSqSum[0]*fr[26]*rdvSq4R)-1.185854122563142*nuVtSqSum[9]*fr[24]*rdvSq4R-1.185854122563142*nuVtSqSum[8]*fr[23]*rdvSq4R-1.185854122563142*nuVtSqSum[7]*fr[22]*rdvSq4R+2.053959590644372*nuVtSqSum[3]*fr[15]*rdvSq4R+2.053959590644372*nuVtSqSum[2]*fr[14]*rdvSq4R+2.053959590644372*nuVtSqSum[1]*fr[13]*rdvSq4R-1.185854122563142*nuVtSqSum[6]*fr[9]*rdvSq4R-1.185854122563142*nuVtSqSum[5]*fr[8]*rdvSq4R-1.185854122563142*nuVtSqSum[4]*fr[7]*rdvSq4R+2.053959590644372*nuVtSqSum[0]*fr[5]*rdvSq4R-1.185854122563142*fr[3]*nuVtSqSum[3]*rdvSq4R-1.185854122563142*fr[2]*nuVtSqSum[2]*rdvSq4R-1.185854122563142*fr[1]*nuVtSqSum[1]*rdvSq4R-1.185854122563142*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[5] += (-0.6846531968814573*nuVtSqSum[0]*fl[26]*rdvSq4L)-0.3061862178478971*nuVtSqSum[9]*fl[24]*rdvSq4L-0.3061862178478971*nuVtSqSum[8]*fl[23]*rdvSq4L-0.3061862178478971*nuVtSqSum[7]*fl[22]*rdvSq4L-0.5303300858899105*nuVtSqSum[3]*fl[15]*rdvSq4L-0.5303300858899105*nuVtSqSum[2]*fl[14]*rdvSq4L-0.5303300858899105*nuVtSqSum[1]*fl[13]*rdvSq4L-0.3061862178478971*nuVtSqSum[6]*fl[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[5]*fl[8]*rdvSq4L-0.3061862178478971*nuVtSqSum[4]*fl[7]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[5]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[13] += (-0.6846531968814573*nuVtSqSum[1]*fl[26]*rdvSq4L)-0.273861278752583*nuVtSqSum[1]*fl[22]*rdvSq4L-0.5303300858899105*nuVtSqSum[5]*fl[15]*rdvSq4L-0.5303300858899105*nuVtSqSum[4]*fl[14]*rdvSq4L-0.4743416490252568*nuVtSqSum[7]*fl[13]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[13]*rdvSq4L-0.3061862178478971*nuVtSqSum[3]*fl[8]*rdvSq4L-0.273861278752583*fl[1]*nuVtSqSum[7]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[7]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[5]*rdvSq4L-0.5303300858899105*nuVtSqSum[1]*fl[5]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[4]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[14] += (-0.6846531968814573*nuVtSqSum[2]*fl[26]*rdvSq4L)-0.273861278752583*nuVtSqSum[2]*fl[23]*rdvSq4L-0.5303300858899105*nuVtSqSum[6]*fl[15]*rdvSq4L-0.4743416490252568*nuVtSqSum[8]*fl[14]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[14]*rdvSq4L-0.5303300858899105*nuVtSqSum[4]*fl[13]*rdvSq4L-0.3061862178478971*nuVtSqSum[3]*fl[9]*rdvSq4L-0.273861278752583*fl[2]*nuVtSqSum[8]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[7]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[6]*rdvSq4L-0.5303300858899105*nuVtSqSum[2]*fl[5]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[4]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[2]*rdvSq4L; 
    outl[15] += (-0.6846531968814573*nuVtSqSum[3]*fl[26]*rdvSq4L)-0.273861278752583*nuVtSqSum[3]*fl[24]*rdvSq4L-0.4743416490252568*nuVtSqSum[9]*fl[15]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[15]*rdvSq4L-0.5303300858899105*nuVtSqSum[6]*fl[14]*rdvSq4L-0.5303300858899105*nuVtSqSum[5]*fl[13]*rdvSq4L-0.273861278752583*fl[3]*nuVtSqSum[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[8]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[6]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[5]*rdvSq4L-0.5303300858899105*nuVtSqSum[3]*fl[5]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[3]*rdvSq4L; 
    outl[16] += (-0.5303300858899105*nuVtSqSum[0]*fl[16]*rdvSq4L)-0.3061862178478971*nuVtSqSum[3]*fl[12]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[11]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[10]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[4]*rdvSq4L; 
    outl[21] += (-0.5303300858899105*nuVtSqSum[0]*fl[21]*rdvSq4L)-0.3061862178478971*nuVtSqSum[3]*fl[19]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[18]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[17]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[6]*rdvSq4L; 
    outl[26] += (-2.651650429449552*nuVtSqSum[0]*fl[26]*rdvSq4L)-1.185854122563142*nuVtSqSum[9]*fl[24]*rdvSq4L-1.185854122563142*nuVtSqSum[8]*fl[23]*rdvSq4L-1.185854122563142*nuVtSqSum[7]*fl[22]*rdvSq4L-2.053959590644372*nuVtSqSum[3]*fl[15]*rdvSq4L-2.053959590644372*nuVtSqSum[2]*fl[14]*rdvSq4L-2.053959590644372*nuVtSqSum[1]*fl[13]*rdvSq4L-1.185854122563142*nuVtSqSum[6]*fl[9]*rdvSq4L-1.185854122563142*nuVtSqSum[5]*fl[8]*rdvSq4L-1.185854122563142*nuVtSqSum[4]*fl[7]*rdvSq4L-2.053959590644372*nuVtSqSum[0]*fl[5]*rdvSq4L-1.185854122563142*fl[3]*nuVtSqSum[3]*rdvSq4L-1.185854122563142*fl[2]*nuVtSqSum[2]*rdvSq4L-1.185854122563142*fl[1]*nuVtSqSum[1]*rdvSq4L-1.185854122563142*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf3x3vMax_VZ_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:       Cell-center coordinates.
  // dxv[6]:     Cell spacing.
  // idx[6]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[12]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[5]*dxvl[5]); 
  double rdvSq4R = 4.0/(dxvr[5]*dxvr[5]); 

  if (idxr[5] == 1) {

    outr[6] += (-0.5303300858899105*nuVtSqSum[0]*fr[6]*rdvSq4R)+0.3061862178478971*fr[3]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[6] += (-0.5303300858899105*nuVtSqSum[0]*fl[6]*rdvSq4L)-0.3061862178478971*fl[3]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf3x3vMax_VZ_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[6]:       Cell-center coordinates.
  // dxv[6]:     Cell spacing.
  // idx[6]:     current grid index.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[30]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[10]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[5]*dxvl[5]); 
  double rdvSq4R = 4.0/(dxvr[5]*dxvr[5]); 

  if (idxr[5] == 1) {

    outr[6] += 0.6846531968814573*nuVtSqSum[0]*fr[27]*rdvSq4R+0.3061862178478971*nuVtSqSum[9]*fr[24]*rdvSq4R+0.3061862178478971*nuVtSqSum[8]*fr[23]*rdvSq4R+0.3061862178478971*nuVtSqSum[7]*fr[22]*rdvSq4R-0.5303300858899105*nuVtSqSum[3]*fr[19]*rdvSq4R-0.5303300858899105*nuVtSqSum[2]*fr[18]*rdvSq4R-0.5303300858899105*nuVtSqSum[1]*fr[17]*rdvSq4R+0.3061862178478971*nuVtSqSum[6]*fr[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[5]*fr[8]*rdvSq4R+0.3061862178478971*nuVtSqSum[4]*fr[7]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[6]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[17] += 0.6846531968814573*nuVtSqSum[1]*fr[27]*rdvSq4R+0.273861278752583*nuVtSqSum[1]*fr[22]*rdvSq4R-0.5303300858899105*nuVtSqSum[5]*fr[19]*rdvSq4R-0.5303300858899105*nuVtSqSum[4]*fr[18]*rdvSq4R-0.4743416490252568*nuVtSqSum[7]*fr[17]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[17]*rdvSq4R+0.3061862178478971*nuVtSqSum[3]*fr[8]*rdvSq4R+0.273861278752583*fr[1]*nuVtSqSum[7]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[7]*rdvSq4R-0.5303300858899105*nuVtSqSum[1]*fr[6]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[5]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[4]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[1]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[18] += 0.6846531968814573*nuVtSqSum[2]*fr[27]*rdvSq4R+0.273861278752583*nuVtSqSum[2]*fr[23]*rdvSq4R-0.5303300858899105*nuVtSqSum[6]*fr[19]*rdvSq4R-0.4743416490252568*nuVtSqSum[8]*fr[18]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[18]*rdvSq4R-0.5303300858899105*nuVtSqSum[4]*fr[17]*rdvSq4R+0.3061862178478971*nuVtSqSum[3]*fr[9]*rdvSq4R+0.273861278752583*fr[2]*nuVtSqSum[8]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[7]*rdvSq4R+0.3061862178478971*fr[3]*nuVtSqSum[6]*rdvSq4R-0.5303300858899105*nuVtSqSum[2]*fr[6]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[4]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[2]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[2]*rdvSq4R; 
    outr[19] += 0.6846531968814573*nuVtSqSum[3]*fr[27]*rdvSq4R+0.273861278752583*nuVtSqSum[3]*fr[24]*rdvSq4R-0.4743416490252568*nuVtSqSum[9]*fr[19]*rdvSq4R-0.5303300858899105*nuVtSqSum[0]*fr[19]*rdvSq4R-0.5303300858899105*nuVtSqSum[6]*fr[18]*rdvSq4R-0.5303300858899105*nuVtSqSum[5]*fr[17]*rdvSq4R+0.273861278752583*fr[3]*nuVtSqSum[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[9]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[8]*rdvSq4R+0.3061862178478971*fr[2]*nuVtSqSum[6]*rdvSq4R-0.5303300858899105*nuVtSqSum[3]*fr[6]*rdvSq4R+0.3061862178478971*fr[1]*nuVtSqSum[5]*rdvSq4R+0.3061862178478971*fr[0]*nuVtSqSum[3]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[3]*rdvSq4R; 
    outr[20] += (-0.5303300858899105*nuVtSqSum[0]*fr[20]*rdvSq4R)+0.3061862178478971*nuVtSqSum[3]*fr[12]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[11]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[10]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[4]*rdvSq4R; 
    outr[21] += (-0.5303300858899105*nuVtSqSum[0]*fr[21]*rdvSq4R)+0.3061862178478971*nuVtSqSum[3]*fr[15]*rdvSq4R+0.3061862178478971*nuVtSqSum[2]*fr[14]*rdvSq4R+0.3061862178478971*nuVtSqSum[1]*fr[13]*rdvSq4R+0.3061862178478971*nuVtSqSum[0]*fr[5]*rdvSq4R; 
    outr[27] += (-2.651650429449552*nuVtSqSum[0]*fr[27]*rdvSq4R)-1.185854122563142*nuVtSqSum[9]*fr[24]*rdvSq4R-1.185854122563142*nuVtSqSum[8]*fr[23]*rdvSq4R-1.185854122563142*nuVtSqSum[7]*fr[22]*rdvSq4R+2.053959590644372*nuVtSqSum[3]*fr[19]*rdvSq4R+2.053959590644372*nuVtSqSum[2]*fr[18]*rdvSq4R+2.053959590644372*nuVtSqSum[1]*fr[17]*rdvSq4R-1.185854122563142*nuVtSqSum[6]*fr[9]*rdvSq4R-1.185854122563142*nuVtSqSum[5]*fr[8]*rdvSq4R-1.185854122563142*nuVtSqSum[4]*fr[7]*rdvSq4R+2.053959590644372*nuVtSqSum[0]*fr[6]*rdvSq4R-1.185854122563142*fr[3]*nuVtSqSum[3]*rdvSq4R-1.185854122563142*fr[2]*nuVtSqSum[2]*rdvSq4R-1.185854122563142*fr[1]*nuVtSqSum[1]*rdvSq4R-1.185854122563142*fr[0]*nuVtSqSum[0]*rdvSq4R; 

  } else {

    outl[6] += (-0.6846531968814573*nuVtSqSum[0]*fl[27]*rdvSq4L)-0.3061862178478971*nuVtSqSum[9]*fl[24]*rdvSq4L-0.3061862178478971*nuVtSqSum[8]*fl[23]*rdvSq4L-0.3061862178478971*nuVtSqSum[7]*fl[22]*rdvSq4L-0.5303300858899105*nuVtSqSum[3]*fl[19]*rdvSq4L-0.5303300858899105*nuVtSqSum[2]*fl[18]*rdvSq4L-0.5303300858899105*nuVtSqSum[1]*fl[17]*rdvSq4L-0.3061862178478971*nuVtSqSum[6]*fl[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[5]*fl[8]*rdvSq4L-0.3061862178478971*nuVtSqSum[4]*fl[7]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[6]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[17] += (-0.6846531968814573*nuVtSqSum[1]*fl[27]*rdvSq4L)-0.273861278752583*nuVtSqSum[1]*fl[22]*rdvSq4L-0.5303300858899105*nuVtSqSum[5]*fl[19]*rdvSq4L-0.5303300858899105*nuVtSqSum[4]*fl[18]*rdvSq4L-0.4743416490252568*nuVtSqSum[7]*fl[17]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[17]*rdvSq4L-0.3061862178478971*nuVtSqSum[3]*fl[8]*rdvSq4L-0.273861278752583*fl[1]*nuVtSqSum[7]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[7]*rdvSq4L-0.5303300858899105*nuVtSqSum[1]*fl[6]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[5]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[4]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[1]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[18] += (-0.6846531968814573*nuVtSqSum[2]*fl[27]*rdvSq4L)-0.273861278752583*nuVtSqSum[2]*fl[23]*rdvSq4L-0.5303300858899105*nuVtSqSum[6]*fl[19]*rdvSq4L-0.4743416490252568*nuVtSqSum[8]*fl[18]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[18]*rdvSq4L-0.5303300858899105*nuVtSqSum[4]*fl[17]*rdvSq4L-0.3061862178478971*nuVtSqSum[3]*fl[9]*rdvSq4L-0.273861278752583*fl[2]*nuVtSqSum[8]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[7]*rdvSq4L-0.3061862178478971*fl[3]*nuVtSqSum[6]*rdvSq4L-0.5303300858899105*nuVtSqSum[2]*fl[6]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[4]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[2]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[2]*rdvSq4L; 
    outl[19] += (-0.6846531968814573*nuVtSqSum[3]*fl[27]*rdvSq4L)-0.273861278752583*nuVtSqSum[3]*fl[24]*rdvSq4L-0.4743416490252568*nuVtSqSum[9]*fl[19]*rdvSq4L-0.5303300858899105*nuVtSqSum[0]*fl[19]*rdvSq4L-0.5303300858899105*nuVtSqSum[6]*fl[18]*rdvSq4L-0.5303300858899105*nuVtSqSum[5]*fl[17]*rdvSq4L-0.273861278752583*fl[3]*nuVtSqSum[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[9]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[8]*rdvSq4L-0.3061862178478971*fl[2]*nuVtSqSum[6]*rdvSq4L-0.5303300858899105*nuVtSqSum[3]*fl[6]*rdvSq4L-0.3061862178478971*fl[1]*nuVtSqSum[5]*rdvSq4L-0.3061862178478971*fl[0]*nuVtSqSum[3]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[3]*rdvSq4L; 
    outl[20] += (-0.5303300858899105*nuVtSqSum[0]*fl[20]*rdvSq4L)-0.3061862178478971*nuVtSqSum[3]*fl[12]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[11]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[10]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[4]*rdvSq4L; 
    outl[21] += (-0.5303300858899105*nuVtSqSum[0]*fl[21]*rdvSq4L)-0.3061862178478971*nuVtSqSum[3]*fl[15]*rdvSq4L-0.3061862178478971*nuVtSqSum[2]*fl[14]*rdvSq4L-0.3061862178478971*nuVtSqSum[1]*fl[13]*rdvSq4L-0.3061862178478971*nuVtSqSum[0]*fl[5]*rdvSq4L; 
    outl[27] += (-2.651650429449552*nuVtSqSum[0]*fl[27]*rdvSq4L)-1.185854122563142*nuVtSqSum[9]*fl[24]*rdvSq4L-1.185854122563142*nuVtSqSum[8]*fl[23]*rdvSq4L-1.185854122563142*nuVtSqSum[7]*fl[22]*rdvSq4L-2.053959590644372*nuVtSqSum[3]*fl[19]*rdvSq4L-2.053959590644372*nuVtSqSum[2]*fl[18]*rdvSq4L-2.053959590644372*nuVtSqSum[1]*fl[17]*rdvSq4L-1.185854122563142*nuVtSqSum[6]*fl[9]*rdvSq4L-1.185854122563142*nuVtSqSum[5]*fl[8]*rdvSq4L-1.185854122563142*nuVtSqSum[4]*fl[7]*rdvSq4L-2.053959590644372*nuVtSqSum[0]*fl[6]*rdvSq4L-1.185854122563142*fl[3]*nuVtSqSum[3]*rdvSq4L-1.185854122563142*fl[2]*nuVtSqSum[2]*rdvSq4L-1.185854122563142*fl[1]*nuVtSqSum[1]*rdvSq4L-1.185854122563142*fl[0]*nuVtSqSum[0]*rdvSq4L; 

  }
  return 0.0; 
} 
