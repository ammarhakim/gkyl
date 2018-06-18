#include <VmLBOModDecl.h> 
double VmLBOconstNuBoundarySurf1x2vMax_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[4]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[2])+0.6123724356957944*fr[1]*nuVtSq[1]; 

    outr[2] += incr2[2]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[0]*((-1.060660171779821*fl[2])-0.6123724356957944*fl[0])-0.6123724356957944*fl[1]*nuVtSq[1]; 

    outl[2] += incr2[2]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x2vMax_VX_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[10]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[0]*(1.369306393762915*fr[8]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[2]*fr[7]+nuVtSq[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[4]); 
    incr2[4] = nuVtSq[1]*(1.369306393762915*fr[8]+0.5477225575051661*fr[7]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.5477225575051661*fr[1]-0.9486832980505137*fr[4])+nuVtSq[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[4]); 
    incr2[6] = nuVtSq[0]*(0.6123724356957944*fr[3]-1.060660171779821*fr[6])+0.6123724356957944*nuVtSq[1]*fr[5]; 
    incr2[8] = nuVtSq[0]*((-5.303300858899105*fr[8])+4.107919181288745*fr[2]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[2]*fr[7]+nuVtSq[1]*(4.107919181288745*fr[4]-2.371708245126284*fr[1]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[4] += incr2[4]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[8] += incr2[8]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[0]*((-1.369306393762915*fl[8])-1.060660171779821*fl[2]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[2]*fl[7]+nuVtSq[1]*((-1.060660171779821*fl[4])-0.6123724356957944*fl[1]); 
    incr2[4] = nuVtSq[1]*((-1.369306393762915*fl[8])-0.5477225575051661*fl[7]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.9486832980505137*fl[4])-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[4])-0.6123724356957944*fl[1]); 
    incr2[6] = nuVtSq[0]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[3])-0.6123724356957944*nuVtSq[1]*fl[5]; 
    incr2[8] = nuVtSq[0]*((-5.303300858899105*fl[8])-4.107919181288745*fl[2]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[2]*fl[7]+nuVtSq[1]*((-4.107919181288745*fl[4])-2.371708245126284*fl[1]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[4] += incr2[4]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[8] += incr2[8]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x2vMax_VX_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*4]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[4]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4r = 4.0/(dxvr[1]*dxvr[1]); 

  double incr2[20]; 
  if (idxr[1] == 1) {

    incr2[2] = nuVtSq[0]*((-1.620185174601965*fr[18])+1.369306393762915*fr[8]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[3]*fr[17]+nuVtSq[1]*(1.369306393762915*fr[12]-1.060660171779821*fr[4]+0.6123724356957944*fr[1])+nuVtSq[2]*(0.6123724356957944*fr[7]-1.060660171779821*fr[11]); 
    incr2[4] = nuVtSq[1]*((-1.620185174601965*fr[18])-0.9486832980505138*fr[11]+1.369306393762915*fr[8]+0.5477225575051661*fr[7]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.537852874200477*fr[17]+1.224744871391589*fr[12]-0.9486832980505137*fr[4]+0.5477225575051661*fr[1])+nuVtSq[0]*(1.369306393762915*fr[12]-1.060660171779821*fr[4]+0.6123724356957944*fr[1])+nuVtSq[3]*(0.537852874200477*fr[7]-0.931588505112178*fr[11]); 
    incr2[6] = nuVtSq[0]*(1.369306393762915*fr[14]-1.060660171779821*fr[6]+0.6123724356957944*fr[3])+0.6123724356957944*nuVtSq[2]*fr[13]+nuVtSq[1]*(0.6123724356957944*fr[5]-1.060660171779821*fr[10]); 
    incr2[8] = nuVtSq[0]*(6.274950199005565*fr[18]-5.303300858899105*fr[8]+4.107919181288745*fr[2]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[3]*fr[17]+nuVtSq[1]*((-5.303300858899106*fr[12])+4.107919181288745*fr[4]-2.371708245126284*fr[1])+nuVtSq[2]*(4.107919181288745*fr[11]-2.371708245126284*fr[7]); 
    incr2[10] = nuVtSq[1]*(1.369306393762915*fr[14]+0.5477225575051661*fr[13]-1.060660171779821*fr[6]+0.6123724356957944*fr[3])+0.5378528742004769*nuVtSq[3]*fr[13]+nuVtSq[2]*(0.5477225575051661*fr[5]-0.9486832980505137*fr[10])+nuVtSq[0]*(0.6123724356957944*fr[5]-1.060660171779821*fr[10]); 
    incr2[11] = nuVtSq[2]*((-1.620185174601965*fr[18])-0.6776309271789384*fr[11]+1.369306393762915*fr[8]+0.3912303982179757*fr[7]-1.060660171779821*fr[2]+0.6123724356957944*fr[0])+nuVtSq[1]*(0.5378528742004769*fr[17]+1.224744871391589*fr[12]-0.9486832980505138*fr[4]+0.5477225575051661*fr[1])+nuVtSq[3]*(0.3651483716701107*fr[17]+1.202675588605909*fr[12]-0.931588505112178*fr[4]+0.5378528742004769*fr[1])+nuVtSq[0]*(0.6123724356957944*fr[7]-1.060660171779821*fr[11]); 
    incr2[12] = nuVtSq[1]*(6.274950199005565*fr[18]+3.674234614174766*fr[11]-5.303300858899106*fr[8]-2.121320343559642*fr[7]+4.107919181288746*fr[2]-2.371708245126284*fr[0])+nuVtSq[2]*((-2.08309522448824*fr[17])-4.743416490252569*fr[12]+3.674234614174767*fr[4]-2.121320343559642*fr[1])+nuVtSq[0]*((-5.303300858899105*fr[12])+4.107919181288746*fr[4]-2.371708245126284*fr[1])+nuVtSq[3]*(3.608026765817728*fr[11]-2.08309522448824*fr[7]); 
    incr2[14] = nuVtSq[0]*((-5.303300858899105*fr[14])+4.107919181288745*fr[6]-2.371708245126284*fr[3])-2.371708245126284*nuVtSq[2]*fr[13]+nuVtSq[1]*(4.107919181288745*fr[10]-2.371708245126284*fr[5]); 
    incr2[16] = nuVtSq[0]*(0.6123724356957944*fr[9]-1.060660171779821*fr[16])+0.6123724356957944*nuVtSq[1]*fr[15]; 
    incr2[18] = nuVtSq[0]*((-14.8492424049175*fr[18])+12.54990039801113*fr[8]-9.721111047611789*fr[2]+5.612486080160912*fr[0])+5.612486080160912*nuVtSq[3]*fr[17]+nuVtSq[1]*(12.54990039801113*fr[12]-9.721111047611789*fr[4]+5.612486080160912*fr[1])+nuVtSq[2]*(5.612486080160912*fr[7]-9.721111047611789*fr[11]); 

    outr[2] += incr2[2]*rdvSq4r; 
    outr[4] += incr2[4]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[8] += incr2[8]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[11] += incr2[11]*rdvSq4r; 
    outr[12] += incr2[12]*rdvSq4r; 
    outr[14] += incr2[14]*rdvSq4r; 
    outr[16] += incr2[16]*rdvSq4r; 
    outr[18] += incr2[18]*rdvSq4r; 

  } else {

    incr2[2] = nuVtSq[0]*((-1.620185174601965*fl[18])-1.369306393762915*fl[8]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[3]*fl[17]+nuVtSq[1]*((-1.369306393762915*fl[12])-1.060660171779821*fl[4]-0.6123724356957944*fl[1])+nuVtSq[2]*((-1.060660171779821*fl[11])-0.6123724356957944*fl[7]); 
    incr2[4] = nuVtSq[1]*((-1.620185174601965*fl[18])-0.9486832980505138*fl[11]-1.369306393762915*fl[8]-0.5477225575051661*fl[7]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.537852874200477*fl[17])-1.224744871391589*fl[12]-0.9486832980505137*fl[4]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.369306393762915*fl[12])-1.060660171779821*fl[4]-0.6123724356957944*fl[1])+nuVtSq[3]*((-0.931588505112178*fl[11])-0.537852874200477*fl[7]); 
    incr2[6] = nuVtSq[0]*((-1.369306393762915*fl[14])-1.060660171779821*fl[6]-0.6123724356957944*fl[3])-0.6123724356957944*nuVtSq[2]*fl[13]+nuVtSq[1]*((-1.060660171779821*fl[10])-0.6123724356957944*fl[5]); 
    incr2[8] = nuVtSq[0]*((-6.274950199005565*fl[18])-5.303300858899105*fl[8]-4.107919181288745*fl[2]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[3]*fl[17]+nuVtSq[1]*((-5.303300858899106*fl[12])-4.107919181288745*fl[4]-2.371708245126284*fl[1])+nuVtSq[2]*((-4.107919181288745*fl[11])-2.371708245126284*fl[7]); 
    incr2[10] = nuVtSq[1]*((-1.369306393762915*fl[14])-0.5477225575051661*fl[13]-1.060660171779821*fl[6]-0.6123724356957944*fl[3])-0.5378528742004769*nuVtSq[3]*fl[13]+nuVtSq[2]*((-0.9486832980505137*fl[10])-0.5477225575051661*fl[5])+nuVtSq[0]*((-1.060660171779821*fl[10])-0.6123724356957944*fl[5]); 
    incr2[11] = nuVtSq[2]*((-1.620185174601965*fl[18])-0.6776309271789384*fl[11]-1.369306393762915*fl[8]-0.3912303982179757*fl[7]-1.060660171779821*fl[2]-0.6123724356957944*fl[0])+nuVtSq[3]*((-0.3651483716701107*fl[17])-1.202675588605909*fl[12]-0.931588505112178*fl[4]-0.5378528742004769*fl[1])+nuVtSq[1]*((-0.5378528742004769*fl[17])-1.224744871391589*fl[12]-0.9486832980505138*fl[4]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[11])-0.6123724356957944*fl[7]); 
    incr2[12] = nuVtSq[1]*((-6.274950199005565*fl[18])-3.674234614174766*fl[11]-5.303300858899106*fl[8]-2.121320343559642*fl[7]-4.107919181288746*fl[2]-2.371708245126284*fl[0])+nuVtSq[2]*((-2.08309522448824*fl[17])-4.743416490252569*fl[12]-3.674234614174767*fl[4]-2.121320343559642*fl[1])+nuVtSq[0]*((-5.303300858899105*fl[12])-4.107919181288746*fl[4]-2.371708245126284*fl[1])+nuVtSq[3]*((-3.608026765817728*fl[11])-2.08309522448824*fl[7]); 
    incr2[14] = nuVtSq[0]*((-5.303300858899105*fl[14])-4.107919181288745*fl[6]-2.371708245126284*fl[3])-2.371708245126284*nuVtSq[2]*fl[13]+nuVtSq[1]*((-4.107919181288745*fl[10])-2.371708245126284*fl[5]); 
    incr2[16] = nuVtSq[0]*((-1.060660171779821*fl[16])-0.6123724356957944*fl[9])-0.6123724356957944*nuVtSq[1]*fl[15]; 
    incr2[18] = nuVtSq[0]*((-14.8492424049175*fl[18])-12.54990039801113*fl[8]-9.721111047611789*fl[2]-5.612486080160912*fl[0])-5.612486080160912*nuVtSq[3]*fl[17]+nuVtSq[1]*((-12.54990039801113*fl[12])-9.721111047611789*fl[4]-5.612486080160912*fl[1])+nuVtSq[2]*((-9.721111047611789*fl[11])-5.612486080160912*fl[7]); 

    outl[2] += incr2[2]*rdvSq4l; 
    outl[4] += incr2[4]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[8] += incr2[8]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[11] += incr2[11]*rdvSq4l; 
    outl[12] += incr2[12]*rdvSq4l; 
    outl[14] += incr2[14]*rdvSq4l; 
    outl[16] += incr2[16]*rdvSq4l; 
    outl[18] += incr2[18]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x2vMax_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*2]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[2]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[4]; 
  if (idxr[2] == 1) {

    incr2[3] = nuVtSq[0]*(0.6123724356957944*fr[0]-1.060660171779821*fr[3])+0.6123724356957944*fr[1]*nuVtSq[1]; 

    outr[3] += incr2[3]*rdvSq4r; 

  } else {

    incr2[3] = nuVtSq[0]*((-1.060660171779821*fl[3])-0.6123724356957944*fl[0])-0.6123724356957944*fl[1]*nuVtSq[1]; 

    outl[3] += incr2[3]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x2vMax_VY_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*3]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[3]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[10]; 
  if (idxr[2] == 1) {

    incr2[3] = nuVtSq[0]*(1.369306393762915*fr[9]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[2]*fr[7]+nuVtSq[1]*(0.6123724356957944*fr[1]-1.060660171779821*fr[5]); 
    incr2[5] = nuVtSq[1]*(1.369306393762915*fr[9]+0.5477225575051661*fr[7]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.5477225575051661*fr[1]-0.9486832980505137*fr[5])+nuVtSq[0]*(0.6123724356957944*fr[1]-1.060660171779821*fr[5]); 
    incr2[6] = nuVtSq[0]*(0.6123724356957944*fr[2]-1.060660171779821*fr[6])+0.6123724356957944*nuVtSq[1]*fr[4]; 
    incr2[9] = nuVtSq[0]*((-5.303300858899105*fr[9])+4.107919181288745*fr[3]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[2]*fr[7]+nuVtSq[1]*(4.107919181288745*fr[5]-2.371708245126284*fr[1]); 

    outr[3] += incr2[3]*rdvSq4r; 
    outr[5] += incr2[5]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[9] += incr2[9]*rdvSq4r; 

  } else {

    incr2[3] = nuVtSq[0]*((-1.369306393762915*fl[9])-1.060660171779821*fl[3]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[2]*fl[7]+nuVtSq[1]*((-1.060660171779821*fl[5])-0.6123724356957944*fl[1]); 
    incr2[5] = nuVtSq[1]*((-1.369306393762915*fl[9])-0.5477225575051661*fl[7]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.9486832980505137*fl[5])-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[5])-0.6123724356957944*fl[1]); 
    incr2[6] = nuVtSq[0]*((-1.060660171779821*fl[6])-0.6123724356957944*fl[2])-0.6123724356957944*nuVtSq[1]*fl[4]; 
    incr2[9] = nuVtSq[0]*((-5.303300858899105*fl[9])-4.107919181288745*fl[3]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[2]*fl[7]+nuVtSq[1]*((-4.107919181288745*fl[5])-2.371708245126284*fl[1]); 

    outl[3] += incr2[3]*rdvSq4l; 
    outl[5] += incr2[5]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[9] += incr2[9]*rdvSq4l; 

  }
  return 0.0; 
} 
double VmLBOconstNuBoundarySurf1x2vMax_VY_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const int *idxl, const int *idxr, const double nu, const double vMuMidMax, const double *nuU, const double *nuVtSq, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:    Cell-center coordinates. dxv[NDIM]: Cell spacing. idx[NDIM]: current grid index.
  // nu:         constant collisionality. 
  // vMuMidMax:  maximum midpoint value of v-u. 
  // nuU[2*4]: bulk velocity (in 2 directions) times by nu. 
  // nuVtSq[4]: thermal speed squared times nu. 
  // fl/fr:      Distribution function in left/right cells 
  // outl/outr:  Incremented distribution function in left/right cells 
  double rdvSq4l = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4r = 4.0/(dxvr[2]*dxvr[2]); 

  double incr2[20]; 
  if (idxr[2] == 1) {

    incr2[3] = nuVtSq[0]*((-1.620185174601965*fr[19])+1.369306393762915*fr[9]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+0.6123724356957944*nuVtSq[3]*fr[17]+nuVtSq[1]*(1.369306393762915*fr[15]-1.060660171779821*fr[5]+0.6123724356957944*fr[1])+nuVtSq[2]*(0.6123724356957944*fr[7]-1.060660171779821*fr[13]); 
    incr2[5] = nuVtSq[1]*((-1.620185174601965*fr[19])-0.9486832980505138*fr[13]+1.369306393762915*fr[9]+0.5477225575051661*fr[7]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+nuVtSq[2]*(0.537852874200477*fr[17]+1.224744871391589*fr[15]-0.9486832980505137*fr[5]+0.5477225575051661*fr[1])+nuVtSq[0]*(1.369306393762915*fr[15]-1.060660171779821*fr[5]+0.6123724356957944*fr[1])+nuVtSq[3]*(0.537852874200477*fr[7]-0.931588505112178*fr[13]); 
    incr2[6] = nuVtSq[0]*(1.369306393762915*fr[16]-1.060660171779821*fr[6]+0.6123724356957944*fr[2])+0.6123724356957944*nuVtSq[2]*fr[11]+nuVtSq[1]*(0.6123724356957944*fr[4]-1.060660171779821*fr[10]); 
    incr2[9] = nuVtSq[0]*(6.274950199005565*fr[19]-5.303300858899105*fr[9]+4.107919181288745*fr[3]-2.371708245126284*fr[0])-2.371708245126284*nuVtSq[3]*fr[17]+nuVtSq[1]*((-5.303300858899106*fr[15])+4.107919181288745*fr[5]-2.371708245126284*fr[1])+nuVtSq[2]*(4.107919181288745*fr[13]-2.371708245126284*fr[7]); 
    incr2[10] = nuVtSq[1]*(1.369306393762915*fr[16]+0.5477225575051661*fr[11]-1.060660171779821*fr[6]+0.6123724356957944*fr[2])+0.5378528742004769*nuVtSq[3]*fr[11]+nuVtSq[2]*(0.5477225575051661*fr[4]-0.9486832980505137*fr[10])+nuVtSq[0]*(0.6123724356957944*fr[4]-1.060660171779821*fr[10]); 
    incr2[13] = nuVtSq[2]*((-1.620185174601965*fr[19])-0.6776309271789384*fr[13]+1.369306393762915*fr[9]+0.3912303982179757*fr[7]-1.060660171779821*fr[3]+0.6123724356957944*fr[0])+nuVtSq[1]*(0.5378528742004769*fr[17]+1.224744871391589*fr[15]-0.9486832980505138*fr[5]+0.5477225575051661*fr[1])+nuVtSq[3]*(0.3651483716701107*fr[17]+1.202675588605909*fr[15]-0.931588505112178*fr[5]+0.5378528742004769*fr[1])+nuVtSq[0]*(0.6123724356957944*fr[7]-1.060660171779821*fr[13]); 
    incr2[14] = nuVtSq[0]*(0.6123724356957944*fr[8]-1.060660171779821*fr[14])+0.6123724356957944*nuVtSq[1]*fr[12]; 
    incr2[15] = nuVtSq[1]*(6.274950199005565*fr[19]+3.674234614174766*fr[13]-5.303300858899106*fr[9]-2.121320343559642*fr[7]+4.107919181288746*fr[3]-2.371708245126284*fr[0])+nuVtSq[2]*((-2.08309522448824*fr[17])-4.743416490252569*fr[15]+3.674234614174767*fr[5]-2.121320343559642*fr[1])+nuVtSq[0]*((-5.303300858899105*fr[15])+4.107919181288746*fr[5]-2.371708245126284*fr[1])+nuVtSq[3]*(3.608026765817728*fr[13]-2.08309522448824*fr[7]); 
    incr2[16] = nuVtSq[0]*((-5.303300858899105*fr[16])+4.107919181288745*fr[6]-2.371708245126284*fr[2])-2.371708245126284*nuVtSq[2]*fr[11]+nuVtSq[1]*(4.107919181288745*fr[10]-2.371708245126284*fr[4]); 
    incr2[19] = nuVtSq[0]*((-14.8492424049175*fr[19])+12.54990039801113*fr[9]-9.721111047611789*fr[3]+5.612486080160912*fr[0])+5.612486080160912*nuVtSq[3]*fr[17]+nuVtSq[1]*(12.54990039801113*fr[15]-9.721111047611789*fr[5]+5.612486080160912*fr[1])+nuVtSq[2]*(5.612486080160912*fr[7]-9.721111047611789*fr[13]); 

    outr[3] += incr2[3]*rdvSq4r; 
    outr[5] += incr2[5]*rdvSq4r; 
    outr[6] += incr2[6]*rdvSq4r; 
    outr[9] += incr2[9]*rdvSq4r; 
    outr[10] += incr2[10]*rdvSq4r; 
    outr[13] += incr2[13]*rdvSq4r; 
    outr[14] += incr2[14]*rdvSq4r; 
    outr[15] += incr2[15]*rdvSq4r; 
    outr[16] += incr2[16]*rdvSq4r; 
    outr[19] += incr2[19]*rdvSq4r; 

  } else {

    incr2[3] = nuVtSq[0]*((-1.620185174601965*fl[19])-1.369306393762915*fl[9]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])-0.6123724356957944*nuVtSq[3]*fl[17]+nuVtSq[1]*((-1.369306393762915*fl[15])-1.060660171779821*fl[5]-0.6123724356957944*fl[1])+nuVtSq[2]*((-1.060660171779821*fl[13])-0.6123724356957944*fl[7]); 
    incr2[5] = nuVtSq[1]*((-1.620185174601965*fl[19])-0.9486832980505138*fl[13]-1.369306393762915*fl[9]-0.5477225575051661*fl[7]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])+nuVtSq[2]*((-0.537852874200477*fl[17])-1.224744871391589*fl[15]-0.9486832980505137*fl[5]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.369306393762915*fl[15])-1.060660171779821*fl[5]-0.6123724356957944*fl[1])+nuVtSq[3]*((-0.931588505112178*fl[13])-0.537852874200477*fl[7]); 
    incr2[6] = nuVtSq[0]*((-1.369306393762915*fl[16])-1.060660171779821*fl[6]-0.6123724356957944*fl[2])-0.6123724356957944*nuVtSq[2]*fl[11]+nuVtSq[1]*((-1.060660171779821*fl[10])-0.6123724356957944*fl[4]); 
    incr2[9] = nuVtSq[0]*((-6.274950199005565*fl[19])-5.303300858899105*fl[9]-4.107919181288745*fl[3]-2.371708245126284*fl[0])-2.371708245126284*nuVtSq[3]*fl[17]+nuVtSq[1]*((-5.303300858899106*fl[15])-4.107919181288745*fl[5]-2.371708245126284*fl[1])+nuVtSq[2]*((-4.107919181288745*fl[13])-2.371708245126284*fl[7]); 
    incr2[10] = nuVtSq[1]*((-1.369306393762915*fl[16])-0.5477225575051661*fl[11]-1.060660171779821*fl[6]-0.6123724356957944*fl[2])-0.5378528742004769*nuVtSq[3]*fl[11]+nuVtSq[2]*((-0.9486832980505137*fl[10])-0.5477225575051661*fl[4])+nuVtSq[0]*((-1.060660171779821*fl[10])-0.6123724356957944*fl[4]); 
    incr2[13] = nuVtSq[2]*((-1.620185174601965*fl[19])-0.6776309271789384*fl[13]-1.369306393762915*fl[9]-0.3912303982179757*fl[7]-1.060660171779821*fl[3]-0.6123724356957944*fl[0])+nuVtSq[3]*((-0.3651483716701107*fl[17])-1.202675588605909*fl[15]-0.931588505112178*fl[5]-0.5378528742004769*fl[1])+nuVtSq[1]*((-0.5378528742004769*fl[17])-1.224744871391589*fl[15]-0.9486832980505138*fl[5]-0.5477225575051661*fl[1])+nuVtSq[0]*((-1.060660171779821*fl[13])-0.6123724356957944*fl[7]); 
    incr2[14] = nuVtSq[0]*((-1.060660171779821*fl[14])-0.6123724356957944*fl[8])-0.6123724356957944*nuVtSq[1]*fl[12]; 
    incr2[15] = nuVtSq[1]*((-6.274950199005565*fl[19])-3.674234614174766*fl[13]-5.303300858899106*fl[9]-2.121320343559642*fl[7]-4.107919181288746*fl[3]-2.371708245126284*fl[0])+nuVtSq[2]*((-2.08309522448824*fl[17])-4.743416490252569*fl[15]-3.674234614174767*fl[5]-2.121320343559642*fl[1])+nuVtSq[0]*((-5.303300858899105*fl[15])-4.107919181288746*fl[5]-2.371708245126284*fl[1])+nuVtSq[3]*((-3.608026765817728*fl[13])-2.08309522448824*fl[7]); 
    incr2[16] = nuVtSq[0]*((-5.303300858899105*fl[16])-4.107919181288745*fl[6]-2.371708245126284*fl[2])-2.371708245126284*nuVtSq[2]*fl[11]+nuVtSq[1]*((-4.107919181288745*fl[10])-2.371708245126284*fl[4]); 
    incr2[19] = nuVtSq[0]*((-14.8492424049175*fl[19])-12.54990039801113*fl[9]-9.721111047611789*fl[3]-5.612486080160912*fl[0])-5.612486080160912*nuVtSq[3]*fl[17]+nuVtSq[1]*((-12.54990039801113*fl[15])-9.721111047611789*fl[5]-5.612486080160912*fl[1])+nuVtSq[2]*((-9.721111047611789*fl[13])-5.612486080160912*fl[7]); 

    outl[3] += incr2[3]*rdvSq4l; 
    outl[5] += incr2[5]*rdvSq4l; 
    outl[6] += incr2[6]*rdvSq4l; 
    outl[9] += incr2[9]*rdvSq4l; 
    outl[10] += incr2[10]*rdvSq4l; 
    outl[13] += incr2[13]*rdvSq4l; 
    outl[14] += incr2[14]*rdvSq4l; 
    outl[15] += incr2[15]*rdvSq4l; 
    outl[16] += incr2[16]*rdvSq4l; 
    outl[19] += incr2[19]*rdvSq4l; 

  }
  return 0.0; 
} 
