#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf1x1vTensor_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const int edge, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:       Cell-center coordinates.
  // dxv[2]:     Cell spacing.
  // idx[2]:     current grid index.
  // edge:          =-1 for lower boundary, =1 for upper boundary.
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[3]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  if (edge < 0) {

    outr[2] += 1.369306393762915*nuVtSqSum[2]*fr[8]*rdvSq4R+1.369306393762915*nuVtSqSum[1]*fr[7]*rdvSq4R-1.060660171779821*nuVtSqSum[2]*fr[6]*rdvSq4R+1.369306393762915*nuVtSqSum[0]*fr[5]*rdvSq4R+0.6123724356957944*nuVtSqSum[2]*fr[4]*rdvSq4R-1.060660171779821*nuVtSqSum[1]*fr[3]*rdvSq4R-1.060660171779821*nuVtSqSum[0]*fr[2]*rdvSq4R+0.6123724356957944*fr[1]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[3] += 1.224744871391589*nuVtSqSum[1]*fr[8]*rdvSq4R+1.224744871391589*nuVtSqSum[2]*fr[7]*rdvSq4R+1.369306393762915*nuVtSqSum[0]*fr[7]*rdvSq4R-0.9486832980505138*nuVtSqSum[1]*fr[6]*rdvSq4R+1.369306393762915*nuVtSqSum[1]*fr[5]*rdvSq4R+0.5477225575051661*nuVtSqSum[1]*fr[4]*rdvSq4R-0.9486832980505137*nuVtSqSum[2]*fr[3]*rdvSq4R-1.060660171779821*nuVtSqSum[0]*fr[3]*rdvSq4R+0.5477225575051661*fr[1]*nuVtSqSum[2]*rdvSq4R-1.060660171779821*nuVtSqSum[1]*fr[2]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[1]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[5] += (-5.303300858899105*nuVtSqSum[2]*fr[8]*rdvSq4R)-5.303300858899106*nuVtSqSum[1]*fr[7]*rdvSq4R+4.107919181288745*nuVtSqSum[2]*fr[6]*rdvSq4R-5.303300858899105*nuVtSqSum[0]*fr[5]*rdvSq4R-2.371708245126284*nuVtSqSum[2]*fr[4]*rdvSq4R+4.107919181288745*nuVtSqSum[1]*fr[3]*rdvSq4R+4.107919181288745*nuVtSqSum[0]*fr[2]*rdvSq4R-2.371708245126284*fr[1]*nuVtSqSum[1]*rdvSq4R-2.371708245126284*fr[0]*nuVtSqSum[0]*rdvSq4R; 
    outr[6] += 0.8748177652797063*nuVtSqSum[2]*fr[8]*rdvSq4R+1.369306393762915*nuVtSqSum[0]*fr[8]*rdvSq4R+1.224744871391589*nuVtSqSum[1]*fr[7]*rdvSq4R-0.6776309271789384*nuVtSqSum[2]*fr[6]*rdvSq4R-1.060660171779821*nuVtSqSum[0]*fr[6]*rdvSq4R+1.369306393762915*nuVtSqSum[2]*fr[5]*rdvSq4R+0.3912303982179757*nuVtSqSum[2]*fr[4]*rdvSq4R+0.6123724356957944*nuVtSqSum[0]*fr[4]*rdvSq4R-0.9486832980505138*nuVtSqSum[1]*fr[3]*rdvSq4R-1.060660171779821*fr[2]*nuVtSqSum[2]*rdvSq4R+0.6123724356957944*fr[0]*nuVtSqSum[2]*rdvSq4R+0.5477225575051661*fr[1]*nuVtSqSum[1]*rdvSq4R; 
    outr[7] += (-4.743416490252569*nuVtSqSum[1]*fr[8]*rdvSq4R)-4.743416490252569*nuVtSqSum[2]*fr[7]*rdvSq4R-5.303300858899105*nuVtSqSum[0]*fr[7]*rdvSq4R+3.674234614174766*nuVtSqSum[1]*fr[6]*rdvSq4R-5.303300858899106*nuVtSqSum[1]*fr[5]*rdvSq4R-2.121320343559642*nuVtSqSum[1]*fr[4]*rdvSq4R+3.674234614174767*nuVtSqSum[2]*fr[3]*rdvSq4R+4.107919181288746*nuVtSqSum[0]*fr[3]*rdvSq4R-2.121320343559642*fr[1]*nuVtSqSum[2]*rdvSq4R+4.107919181288746*nuVtSqSum[1]*fr[2]*rdvSq4R-2.371708245126284*fr[0]*nuVtSqSum[1]*rdvSq4R-2.371708245126284*nuVtSqSum[0]*fr[1]*rdvSq4R; 
    outr[8] += (-3.388154635894693*nuVtSqSum[2]*fr[8]*rdvSq4R)-5.303300858899105*nuVtSqSum[0]*fr[8]*rdvSq4R-4.743416490252569*nuVtSqSum[1]*fr[7]*rdvSq4R+2.624453295839119*nuVtSqSum[2]*fr[6]*rdvSq4R+4.107919181288745*nuVtSqSum[0]*fr[6]*rdvSq4R-5.303300858899105*nuVtSqSum[2]*fr[5]*rdvSq4R-1.515228816828316*nuVtSqSum[2]*fr[4]*rdvSq4R-2.371708245126284*nuVtSqSum[0]*fr[4]*rdvSq4R+3.674234614174766*nuVtSqSum[1]*fr[3]*rdvSq4R+4.107919181288745*fr[2]*nuVtSqSum[2]*rdvSq4R-2.371708245126284*fr[0]*nuVtSqSum[2]*rdvSq4R-2.121320343559642*fr[1]*nuVtSqSum[1]*rdvSq4R; 

  } else {

    outl[2] += (-1.369306393762915*nuVtSqSum[2]*fl[8]*rdvSq4L)-1.369306393762915*nuVtSqSum[1]*fl[7]*rdvSq4L-1.060660171779821*nuVtSqSum[2]*fl[6]*rdvSq4L-1.369306393762915*nuVtSqSum[0]*fl[5]*rdvSq4L-0.6123724356957944*nuVtSqSum[2]*fl[4]*rdvSq4L-1.060660171779821*nuVtSqSum[1]*fl[3]*rdvSq4L-1.060660171779821*nuVtSqSum[0]*fl[2]*rdvSq4L-0.6123724356957944*fl[1]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[3] += (-1.224744871391589*nuVtSqSum[1]*fl[8]*rdvSq4L)-1.224744871391589*nuVtSqSum[2]*fl[7]*rdvSq4L-1.369306393762915*nuVtSqSum[0]*fl[7]*rdvSq4L-0.9486832980505138*nuVtSqSum[1]*fl[6]*rdvSq4L-1.369306393762915*nuVtSqSum[1]*fl[5]*rdvSq4L-0.5477225575051661*nuVtSqSum[1]*fl[4]*rdvSq4L-0.9486832980505137*nuVtSqSum[2]*fl[3]*rdvSq4L-1.060660171779821*nuVtSqSum[0]*fl[3]*rdvSq4L-0.5477225575051661*fl[1]*nuVtSqSum[2]*rdvSq4L-1.060660171779821*nuVtSqSum[1]*fl[2]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[1]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[5] += (-5.303300858899105*nuVtSqSum[2]*fl[8]*rdvSq4L)-5.303300858899106*nuVtSqSum[1]*fl[7]*rdvSq4L-4.107919181288745*nuVtSqSum[2]*fl[6]*rdvSq4L-5.303300858899105*nuVtSqSum[0]*fl[5]*rdvSq4L-2.371708245126284*nuVtSqSum[2]*fl[4]*rdvSq4L-4.107919181288745*nuVtSqSum[1]*fl[3]*rdvSq4L-4.107919181288745*nuVtSqSum[0]*fl[2]*rdvSq4L-2.371708245126284*fl[1]*nuVtSqSum[1]*rdvSq4L-2.371708245126284*fl[0]*nuVtSqSum[0]*rdvSq4L; 
    outl[6] += (-0.8748177652797063*nuVtSqSum[2]*fl[8]*rdvSq4L)-1.369306393762915*nuVtSqSum[0]*fl[8]*rdvSq4L-1.224744871391589*nuVtSqSum[1]*fl[7]*rdvSq4L-0.6776309271789384*nuVtSqSum[2]*fl[6]*rdvSq4L-1.060660171779821*nuVtSqSum[0]*fl[6]*rdvSq4L-1.369306393762915*nuVtSqSum[2]*fl[5]*rdvSq4L-0.3912303982179757*nuVtSqSum[2]*fl[4]*rdvSq4L-0.6123724356957944*nuVtSqSum[0]*fl[4]*rdvSq4L-0.9486832980505138*nuVtSqSum[1]*fl[3]*rdvSq4L-1.060660171779821*fl[2]*nuVtSqSum[2]*rdvSq4L-0.6123724356957944*fl[0]*nuVtSqSum[2]*rdvSq4L-0.5477225575051661*fl[1]*nuVtSqSum[1]*rdvSq4L; 
    outl[7] += (-4.743416490252569*nuVtSqSum[1]*fl[8]*rdvSq4L)-4.743416490252569*nuVtSqSum[2]*fl[7]*rdvSq4L-5.303300858899105*nuVtSqSum[0]*fl[7]*rdvSq4L-3.674234614174766*nuVtSqSum[1]*fl[6]*rdvSq4L-5.303300858899106*nuVtSqSum[1]*fl[5]*rdvSq4L-2.121320343559642*nuVtSqSum[1]*fl[4]*rdvSq4L-3.674234614174767*nuVtSqSum[2]*fl[3]*rdvSq4L-4.107919181288746*nuVtSqSum[0]*fl[3]*rdvSq4L-2.121320343559642*fl[1]*nuVtSqSum[2]*rdvSq4L-4.107919181288746*nuVtSqSum[1]*fl[2]*rdvSq4L-2.371708245126284*fl[0]*nuVtSqSum[1]*rdvSq4L-2.371708245126284*nuVtSqSum[0]*fl[1]*rdvSq4L; 
    outl[8] += (-3.388154635894693*nuVtSqSum[2]*fl[8]*rdvSq4L)-5.303300858899105*nuVtSqSum[0]*fl[8]*rdvSq4L-4.743416490252569*nuVtSqSum[1]*fl[7]*rdvSq4L-2.624453295839119*nuVtSqSum[2]*fl[6]*rdvSq4L-4.107919181288745*nuVtSqSum[0]*fl[6]*rdvSq4L-5.303300858899105*nuVtSqSum[2]*fl[5]*rdvSq4L-1.515228816828316*nuVtSqSum[2]*fl[4]*rdvSq4L-2.371708245126284*nuVtSqSum[0]*fl[4]*rdvSq4L-3.674234614174766*nuVtSqSum[1]*fl[3]*rdvSq4L-4.107919181288745*fl[2]*nuVtSqSum[2]*rdvSq4L-2.371708245126284*fl[0]*nuVtSqSum[2]*rdvSq4L-2.121320343559642*fl[1]*nuVtSqSum[1]*rdvSq4L; 

  }
  return 0.0; 
} 
